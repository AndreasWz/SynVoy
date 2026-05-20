#!/usr/bin/env python3
"""
plm_search.py — Protein Language Model embedding search for SynVoy

Uses ProtT5-XL-UniRef50 (encoder-only) to generate per-protein embeddings
and find remote homologs by cosine similarity in embedding space.  ProtT5
mean embeddings outperform ESM-2 and even MMseqs2-sensitive for remote
homology detection (Spearman ≥ 0.91 on structural similarity benchmarks).

This enables SynVoy to catch orthologs that have diverged beyond sequence-
based detection (<15 % identity), because protein language model embeddings
capture structural and functional similarity far better than substitution
matrices.

Two integration modes:
  1. ORF discovery  — predict ORFs (Prodigal) in a syntenic block, embed
     them, and find hits invisible to MMseqs2/tblastn/SW.
  2. Confidence boost — compute embedding similarity for annotated models
     and feed the score into the classification system.

Requirements (optional — pipeline works without them):
  pip install torch transformers sentencepiece
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    from sequence_utils import parse_fasta, write_fasta
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta

try:
    from gene_predictor import predict_orfs as _predict_orfs_unified
    _GENE_PREDICTOR_AVAILABLE = True
except ImportError:
    _GENE_PREDICTOR_AVAILABLE = False

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------
_PLM_AVAILABLE: Optional[bool] = None

def check_plm_available() -> bool:
    """Return True if torch + transformers are importable."""
    global _PLM_AVAILABLE
    if _PLM_AVAILABLE is None:
        try:
            import torch                          # noqa: F401
            from transformers import T5EncoderModel, T5Tokenizer  # noqa: F401
            _PLM_AVAILABLE = True
        except ImportError:
            _PLM_AVAILABLE = False
    return _PLM_AVAILABLE


PLM_MODEL_ID = "Rostlab/prot_t5_xl_half_uniref50-enc"

# ---------------------------------------------------------------------------
# Module-level model cache (singleton per process)
# ---------------------------------------------------------------------------
_model_cache: Dict[str, Any] = {
    "tokenizer": None,
    "model": None,
    "device": None,
}


def _chunk_sequence(seq: str, chunk_size: int, overlap: int) -> List[str]:
    """Split long protein sequences into overlapping chunks."""
    if chunk_size <= 0 or len(seq) <= chunk_size:
        return [seq]

    overlap = max(0, min(overlap, chunk_size - 1))
    step = max(1, chunk_size - overlap)
    chunks: List[str] = []
    for start in range(0, len(seq), step):
        chunk = seq[start : start + chunk_size]
        if chunk:
            chunks.append(chunk)
        if start + chunk_size >= len(seq):
            break
    return chunks


def _ensure_model_loaded(device: str = "cpu") -> None:
    """Load ProtT5 encoder into the process-level cache (once)."""
    if _model_cache["model"] is not None:
        return

    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    logger.info(f"Loading ProtT5 model ({PLM_MODEL_ID}) on {device} ...")
    tokenizer = T5Tokenizer.from_pretrained(PLM_MODEL_ID, do_lower_case=False)
    load_kwargs: Dict[str, Any] = {}
    if device != "cpu" and torch.cuda.is_available():
        load_kwargs["torch_dtype"] = torch.float16
    model = T5EncoderModel.from_pretrained(PLM_MODEL_ID, **load_kwargs)

    # Half-precision weights must be cast to float32 on CPU
    if device == "cpu":
        model = model.float()

    model = model.to(torch.device(device)).eval()

    _model_cache["tokenizer"] = tokenizer
    _model_cache["model"] = model
    _model_cache["device"] = device
    logger.info("ProtT5 model loaded successfully.")


# ---------------------------------------------------------------------------
# Embedding computation
# ---------------------------------------------------------------------------

def embed_proteins(
    sequences: List[Tuple[str, str]],
    device: str = "cpu",
    batch_size: int = 1,
    max_length: int = 1024,
    chunk_overlap: int = 128,
) -> Dict[str, np.ndarray]:
    """
    Compute ProtT5 mean-pool embeddings for protein sequences.

    Args:
        sequences: list of (id, amino_acid_sequence) pairs
        device:    'cpu' or 'cuda'
        batch_size: proteins per forward pass
        max_length: maximum residues per chunk
        chunk_overlap: overlap between chunks for long proteins

    Returns:
        dict  {sequence_id: np.ndarray(1024,)}
    """
    import torch

    _ensure_model_loaded(device)
    tokenizer = _model_cache["tokenizer"]
    model = _model_cache["model"]
    dev = torch.device(device)

    embeddings: Dict[str, np.ndarray] = {}

    for batch_start in range(0, len(sequences), batch_size):
        batch = sequences[batch_start : batch_start + batch_size]

        prepared: List[str] = []
        batch_ids: List[str] = []
        batch_chunk_map: List[Tuple[str, int]] = []
        for seq_id, seq in batch:
            if not seq or len(seq) < 5:
                continue
            # Replace rare / non-standard amino acids with X
            seq = re.sub(r"[UZOB]", "X", seq)
            for chunk_index, chunk in enumerate(_chunk_sequence(seq, max_length, chunk_overlap)):
                # ProtT5 expects space-separated single characters
                prepared.append(" ".join(list(chunk)))
                batch_ids.append(seq_id)
                batch_chunk_map.append((seq_id, chunk_index))

        if not prepared:
            continue

        encoded = tokenizer(
            prepared,
            add_special_tokens=True,
            padding="longest",
            truncation=True,
            max_length=max_length + 1,
            return_tensors="pt",
        )

        input_ids = encoded["input_ids"].to(dev)
        attention_mask = encoded["attention_mask"].to(dev)

        with torch.no_grad():
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)

        hidden = outputs.last_hidden_state  # (batch, seq_len, 1024)
        chunk_embeddings: Dict[str, List[Tuple[int, np.ndarray]]] = {}

        for i, seq_id in enumerate(batch_ids):
            mask = attention_mask[i].unsqueeze(-1).float()  # (seq_len, 1)
            mask_sum = mask.sum()
            if mask_sum.item() == 0:
                continue
            emb = (hidden[i] * mask).sum(dim=0) / mask_sum
            chunk_embeddings.setdefault(seq_id, []).append(
                (batch_chunk_map[i][1], emb.cpu().numpy().astype(np.float32))
            )

        for seq_id, seq_chunks in chunk_embeddings.items():
            seq_chunks.sort(key=lambda item: item[0])
            vectors = np.stack([emb for _, emb in seq_chunks])
            embeddings[seq_id] = vectors.mean(axis=0).astype(np.float32)

    return embeddings


# ---------------------------------------------------------------------------
# Similarity helpers
# ---------------------------------------------------------------------------

def cosine_similarity_matrix(
    query_embs: Dict[str, np.ndarray],
    target_embs: Dict[str, np.ndarray],
) -> Dict[Tuple[str, str], float]:
    """Pairwise cosine similarities between two sets of embeddings."""
    if not query_embs or not target_embs:
        return {}

    q_ids = list(query_embs.keys())
    t_ids = list(target_embs.keys())

    q_mat = np.stack([query_embs[k] for k in q_ids])   # (nq, d)
    t_mat = np.stack([target_embs[k] for k in t_ids])   # (nt, d)

    # L2-normalise rows
    q_norms = np.linalg.norm(q_mat, axis=1, keepdims=True)
    t_norms = np.linalg.norm(t_mat, axis=1, keepdims=True)
    q_norms[q_norms == 0] = 1.0
    t_norms[t_norms == 0] = 1.0

    sim = (q_mat / q_norms) @ (t_mat / t_norms).T  # (nq, nt)

    results: Dict[Tuple[str, str], float] = {}
    for qi, qid in enumerate(q_ids):
        for ti, tid in enumerate(t_ids):
            results[(qid, tid)] = float(sim[qi, ti])
    return results


def best_similarities(
    query_embs: Dict[str, np.ndarray],
    target_embs: Dict[str, np.ndarray],
) -> Dict[str, float]:
    """For each target, return maximum cosine similarity to any query."""
    sims = cosine_similarity_matrix(query_embs, target_embs)
    best: Dict[str, float] = {}
    for (_, tid), s in sims.items():
        if tid not in best or s > best[tid]:
            best[tid] = s
    return best


# ---------------------------------------------------------------------------
# Embedding I/O
# ---------------------------------------------------------------------------

def save_embeddings(embeddings: Dict[str, np.ndarray], path: str) -> None:
    """Save embeddings as  {ids: array, embeddings: matrix}  .npz file."""
    if not embeddings:
        return
    ids = list(embeddings.keys())
    matrix = np.stack([embeddings[k] for k in ids])
    np.savez_compressed(path, ids=np.array(ids, dtype=object), embeddings=matrix)


def load_embeddings(path: str) -> Dict[str, np.ndarray]:
    """Load embeddings written by save_embeddings()."""
    data = np.load(path, allow_pickle=True)
    ids = data["ids"]
    matrix = data["embeddings"]
    return {str(ids[i]): matrix[i] for i in range(len(ids))}


# ---------------------------------------------------------------------------
# ORF prediction (Prodigal)
# ---------------------------------------------------------------------------

def predict_orfs_prodigal(
    region_fasta: str,
    output_dir: str,
    meta_mode: bool = True,
    min_aa: int = 20,
    predictor: str = "auto",
    augustus_species: str = "fly",
) -> List[Dict[str, Any]]:
    """
    Predict ORFs/genes in a genomic region.

    Despite the legacy name, this now dispatches to the unified gene_predictor
    module which supports both Augustus (eukaryotes) and Prodigal (prokaryotes).
    Augustus is preferred for eukaryotic genomes because Prodigal cannot predict
    intron-containing genes.

    Args:
        region_fasta:    input region FASTA
        output_dir:      working directory
        meta_mode:       (legacy, ignored when using unified predictor)
        min_aa:          minimum protein length
        predictor:       'auto', 'augustus', or 'prodigal'
        augustus_species: Augustus species model (e.g. 'fly', 'human')

    Returns list of dicts:
        {id, seq, start (0-based), end (exclusive), strand}
    """
    # Use unified gene predictor if available
    if _GENE_PREDICTOR_AVAILABLE:
        return _predict_orfs_unified(
            region_fasta, output_dir,
            predictor=predictor,
            augustus_species=augustus_species,
            min_aa=min_aa,
        )

    # Fallback: direct Prodigal call (legacy path, only if gene_predictor
    # module is not importable)
    logger.warning("gene_predictor module not available; falling back to direct Prodigal")

    with tempfile.TemporaryDirectory(prefix="plm_search_prodigal_") as temp_dir:
        proteins_file = os.path.join(temp_dir, "prodigal_orfs.faa")
        gff_file = os.path.join(temp_dir, "prodigal_orfs.gff")

        cmd = ["prodigal", "-i", region_fasta, "-a", proteins_file, "-f", "gff", "-o", gff_file]
        if meta_mode:
            cmd.extend(["-p", "meta"])

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                logger.warning(f"Prodigal exited {result.returncode}: {result.stderr[:200]}")
                return []
        except subprocess.TimeoutExpired as exc:
            logger.warning(f"Prodigal timed out after 120 s: {exc}")
            return []
        except FileNotFoundError as exc:
            logger.warning(f"Prodigal unavailable: {exc}")
            return []

        if not os.path.exists(proteins_file) or os.path.getsize(proteins_file) == 0:
            return []

        orfs: List[Dict[str, Any]] = []
        current_id: Optional[str] = None
        current_seq: List[str] = []
        current_meta: Dict[str, Any] = {}

        with open(proteins_file) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_seq:
                        seq = "".join(current_seq).rstrip("*")
                        if len(seq) >= min_aa:
                            orfs.append({"id": current_id, "seq": seq, **current_meta})

                    # Prodigal header: >contig_N # start # end # strand # ID=...
                    parts = line[1:].split(" # ")
                    current_id = parts[0].strip()
                    current_seq = []
                    current_meta = {}
                    if len(parts) >= 4:
                        try:
                            current_meta["start"] = int(parts[1]) - 1   # 0-based
                            current_meta["end"] = int(parts[2])          # exclusive
                            current_meta["strand"] = "+" if parts[3] == "1" else "-"
                        except (ValueError, IndexError):
                            pass
                else:
                    current_seq.append(line)

        # last record
        if current_id and current_seq:
            seq = "".join(current_seq).rstrip("*")
            if len(seq) >= min_aa:
                orfs.append({"id": current_id, "seq": seq, **current_meta})

        return orfs


# ---------------------------------------------------------------------------
# Region-level PLM search  (called from iterative_search_runner)
# ---------------------------------------------------------------------------

def plm_search_region(
    goi_embeddings: Dict[str, np.ndarray],
    region_fasta: str,
    output_dir: str,
    similarity_threshold: float = 0.5,
    device: str = "cpu",
    predictor: str = "augustus",  # Changed: enforce Augustus for eukaryotic ORF prediction
    augustus_species: str = "fly",
) -> List[Dict[str, Any]]:
    """
    Search for GOI homologs in a genomic region via PLM embeddings.

    1. Predict ORFs/genes with Augustus (eukaryotic-capable; no Prodigal fallback to prevent silent multi-exon loss)
    2. Embed proteins with ProtT5
    3. Compare against pre-computed GOI embeddings
    4. Return hits above similarity threshold

    Returns list of hit dicts compatible with the augmented-search pipeline.
    """
    orfs = predict_orfs_prodigal(
        region_fasta, output_dir,
        predictor=predictor, augustus_species=augustus_species,
    )
    if not orfs:
        return []

    orf_sequences = [(orf["id"], orf["seq"]) for orf in orfs]
    orf_embeddings = embed_proteins(orf_sequences, device=device)
    if not orf_embeddings:
        return []

    similarities = cosine_similarity_matrix(goi_embeddings, orf_embeddings)

    orf_by_id = {orf["id"]: orf for orf in orfs}
    hits: List[Dict[str, Any]] = []

    for (goi_id, orf_id), sim in similarities.items():
        if sim < similarity_threshold:
            continue
        orf = orf_by_id.get(orf_id)
        if not orf:
            continue

        hits.append({
            "query": goi_id,
            "chrom": "region_seq",
            "start": orf.get("start", 0),
            "end": orf.get("end", 0),
            "strand": orf.get("strand", "+"),
            "pident": sim * 100.0,
            "identity": sim * 100.0,
            "evalue": max(1e-10, 1.0 - sim),
            "bits": sim * 200.0,
            "alnlen": len(orf.get("seq", "")),
            "qstart": 1,
            "qend": len(orf.get("seq", "")),
            "method": "plm_embedding",
            "embedding_similarity": sim,
        })

    hits.sort(key=lambda h: -h.get("embedding_similarity", 0))
    return hits


# ---------------------------------------------------------------------------
# Re-ranking helper  (compute embedding similarity for existing candidates)
# ---------------------------------------------------------------------------

def compute_candidate_similarities(
    candidate_sequences: List[Tuple[str, str]],
    goi_embeddings: Dict[str, np.ndarray],
    device: str = "cpu",
) -> Dict[str, float]:
    """
    Embed candidate protein models and return their max cosine similarity
    to any GOI query.  Used to feed embedding_similarity into classification.

    Args:
        candidate_sequences: [(model_id, protein_seq), ...]
        goi_embeddings:      pre-computed GOI embeddings

    Returns:
        {model_id: max_cosine_similarity}
    """
    if not candidate_sequences or not goi_embeddings:
        return {}
    cand_embs = embed_proteins(candidate_sequences, device=device)
    return best_similarities(goi_embeddings, cand_embs)


# ---------------------------------------------------------------------------
# CLI for pre-computing GOI embeddings
# ---------------------------------------------------------------------------

def precompute_goi_embeddings(
    db_fasta: str,
    output_path: str,
    device: str = "cpu",
    goi_prefix: str = "GOI_",
) -> Dict[str, np.ndarray]:
    """
    Read the initial database FASTA, extract GOI sequences, embed them,
    and save to an .npz file.
    """
    sequences: List[Tuple[str, str]] = []
    for header, clean_id, seq in parse_fasta(db_fasta):
        if clean_id.startswith(goi_prefix):
            sequences.append((clean_id, seq))

    if not sequences:
        logger.warning("No GOI sequences found in database for PLM embedding.")
        return {}

    logger.info(f"Embedding {len(sequences)} GOI sequence(s) with ProtT5 ...")
    embeddings = embed_proteins(sequences, device=device)
    save_embeddings(embeddings, output_path)
    logger.info(f"GOI embeddings saved to {output_path}")
    return embeddings


def main():
    parser = argparse.ArgumentParser(
        description="Pre-compute ProtT5 embeddings for SynVoy GOI sequences"
    )
    parser.add_argument("--db_fasta", required=True, help="Initial database FASTA")
    parser.add_argument("--output", required=True, help="Output .npz path")
    parser.add_argument("--device", default="cpu", help="cpu or cuda")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    if not check_plm_available():
        logger.error(
            "PLM search requires PyTorch and HuggingFace Transformers.\n"
            "  Install with: pip install torch transformers sentencepiece"
        )
        sys.exit(1)

    precompute_goi_embeddings(args.db_fasta, args.output, device=args.device)


if __name__ == "__main__":
    main()
