#!/usr/bin/env python3
"""
structural_search.py — Foldseek structural homology search for SynVoy

Uses ESMFold for on-the-fly protein structure prediction and Foldseek for
ultra-fast structural comparison via the 3Di structural alphabet.  This
catches orthologs that have diverged beyond sequence-based detection,
because protein structure diverges ~10x slower than sequence.

Integration modes:
  1. ORF discovery — predict ORFs (Prodigal) in a syntenic block, fold them
     with ESMFold, and compare against the pre-folded GOI structure via
     Foldseek.  Finds hits invisible to MMseqs2/tblastn/SW/PLM.
  2. Confidence boost — compute structural similarity (TM-score) for
     annotated models and feed the score into the classification system.

Requirements (optional — pipeline works without them):
  pip install torch
  conda install -c bioconda foldseek
  # ESMFold model is auto-downloaded from Meta/FAIR on first use
"""

import argparse
import logging
import os
import re
import shutil
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

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------------
_ESMFOLD_AVAILABLE: Optional[bool] = None
_FOLDSEEK_AVAILABLE: Optional[bool] = None


def check_esmfold_available() -> bool:
    """Return True if torch and ESMFold model are importable."""
    global _ESMFOLD_AVAILABLE
    if _ESMFOLD_AVAILABLE is None:
        try:
            import torch  # noqa: F401
            from transformers import EsmForProteinFolding, AutoTokenizer  # noqa: F401
            _ESMFOLD_AVAILABLE = True
        except ImportError:
            _ESMFOLD_AVAILABLE = False
    return _ESMFOLD_AVAILABLE


def check_foldseek_available() -> bool:
    """Return True if foldseek binary is on PATH."""
    global _FOLDSEEK_AVAILABLE
    if _FOLDSEEK_AVAILABLE is None:
        _FOLDSEEK_AVAILABLE = shutil.which("foldseek") is not None
    return _FOLDSEEK_AVAILABLE


def check_structural_search_available() -> bool:
    """Return True if both ESMFold and Foldseek are available."""
    return check_esmfold_available() and check_foldseek_available()


# ---------------------------------------------------------------------------
# Module-level model cache (singleton per process)
# ---------------------------------------------------------------------------
_esmfold_cache: Dict[str, Any] = {
    "tokenizer": None,
    "model": None,
    "device": None,
}


def _ensure_esmfold_loaded(device: str = "cpu") -> None:
    """Load ESMFold into the process-level cache (once)."""
    if _esmfold_cache["model"] is not None:
        return

    import torch
    from transformers import EsmForProteinFolding, AutoTokenizer

    logger.info(f"Loading ESMFold model on {device} ...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")

    # Always use float32 — half precision causes NaN in pTM score
    # computation for short sequences (transformers ≥5.x bug in
    # openfold_utils/loss.py compute_tm).
    model = model.float()

    model = model.to(torch.device(device)).eval()

    # ESMFold chunk_size trades throughput for peak memory. Auto-size based on
    # available VRAM so low-end GPUs (e.g. GTX 1650 / 4 GB) don't OOM mid-fold.
    chunk_size = _recommended_chunk_size(device, default=64)
    model.trunk.set_chunk_size(chunk_size)
    if chunk_size != 64:
        logger.info(f"ESMFold chunk_size set to {chunk_size} (VRAM-tuned)")

    _esmfold_cache["tokenizer"] = tokenizer
    _esmfold_cache["model"] = model
    _esmfold_cache["device"] = device
    logger.info("ESMFold model loaded successfully.")


# ---------------------------------------------------------------------------
# Structure prediction
# ---------------------------------------------------------------------------

def _probe_vram_gb() -> Optional[float]:
    """Return total VRAM (GB) on the first CUDA device, or None if unavailable."""
    try:
        import torch
        if not torch.cuda.is_available():
            return None
        total_bytes = torch.cuda.get_device_properties(0).total_memory
        return total_bytes / (1024 ** 3)
    except Exception as exc:
        logger.debug(f"VRAM probe failed: {exc}")
        return None


def _vram_tier_caps(vram_gb: float) -> Tuple[int, int]:
    """Return (max_length_cap, chunk_size) for a given VRAM.

    Tier table:
      <6 GB   → 150 aa, chunk_size 32   (GTX 1650 / student laptops)
      <10 GB  → 300 aa, chunk_size 48
      <20 GB  → 400 aa, chunk_size 64   (prior default for non-large GPUs)
      >=20 GB → no cap,  chunk_size 64
    """
    if vram_gb < 6:
        return 150, 32
    if vram_gb < 10:
        return 300, 48
    if vram_gb < 20:
        return 400, 64
    return 10_000, 64  # effectively uncapped


def _effective_max_length(device: str, requested: int) -> int:
    """Cap max_length based on available VRAM to avoid OOM at fold time.

    ESMFold memory scales quadratically in sequence length. The tiered caps
    in `_vram_tier_caps` match common consumer/student GPUs:
      - 4 GB  (GTX 1650)        → 150 aa
      - 8 GB  (RTX 3060/4060)   → 300 aa
      - 16 GB (RTX 3090 desktop chunks) → 400 aa
      - 24+ GB (A100/H100 slices)       → user's request honored.
    """
    if device != "cuda":
        return requested
    vram_gb = _probe_vram_gb()
    if vram_gb is None:
        return requested
    cap, _ = _vram_tier_caps(vram_gb)
    if requested > cap:
        logger.warning(
            f"GPU has {vram_gb:.1f} GB VRAM; capping structural_max_length "
            f"from {requested} to {cap} to avoid OOM. Folding will truncate "
            f"queries longer than {cap} aa."
        )
        return cap
    return requested


def _recommended_chunk_size(device: str, default: int = 64) -> int:
    """Return an ESMFold trunk chunk_size appropriate for the current GPU.

    Smaller chunk_size lowers peak memory at the cost of throughput. On tight
    VRAM (<6 GB) we drop to 32 which keeps the GTX 1650 class GPUs from OOM.
    On CPU / well-provisioned GPUs we keep the default 64.
    """
    if device != "cuda":
        return default
    vram_gb = _probe_vram_gb()
    if vram_gb is None:
        return default
    _, chunk = _vram_tier_caps(vram_gb)
    return chunk


def fold_protein(
    sequence: str,
    output_pdb: str,
    device: str = "cpu",
    max_length: int = 700,
) -> Optional[str]:
    """
    Predict 3D structure for a single protein sequence using ESMFold.

    Args:
        sequence:   amino acid sequence (single-letter, no spaces)
        output_pdb: path to write PDB output
        device:     'cpu' or 'cuda'
        max_length: truncate sequences longer than this (VRAM safety)

    Returns:
        path to PDB file, or None if folding failed
    """
    import torch

    _ensure_esmfold_loaded(device)
    tokenizer = _esmfold_cache["tokenizer"]
    model = _esmfold_cache["model"]

    max_length = _effective_max_length(device, max_length)
    seq = sequence[:max_length]
    # Replace non-standard amino acids
    seq = re.sub(r"[UZOB]", "X", seq)

    try:
        inputs = tokenizer(seq, return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.to(torch.device(device)) for k, v in inputs.items()}

        with torch.no_grad():
            output = model(**inputs)

        # Convert to PDB
        pdb_string = _output_to_pdb(output, seq)
        with open(output_pdb, "w") as f:
            f.write(pdb_string)

        return output_pdb

    except Exception as exc:
        logger.warning(f"ESMFold failed for sequence (len={len(seq)}): {exc}", exc_info=True)
        return None


def _output_to_pdb(output, sequence: str) -> str:
    """Convert ESMFold model output tensors to a PDB-format string."""
    import torch

    # Extract coordinates: (8, 1, L, 14, 3) — last recycling layer
    positions = output["positions"][-1]  # (1, L, 14, 3)
    positions = positions[0].cpu().numpy()  # (L, 14, 3)

    # atom14 layout: 0=N, 1=CA, 2=C, 3=O, 4=CB
    atom_names = ["N", "CA", "C", "O", "CB"]
    atom_indices = [0, 1, 2, 3, 4]

    aa_3letter = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
        'X': 'UNK',
    }

    lines = []
    atom_serial = 1
    for res_idx, aa in enumerate(sequence[:positions.shape[0]]):
        res_name = aa_3letter.get(aa.upper(), "UNK")
        res_num = res_idx + 1

        for atom_name, atom_idx in zip(atom_names, atom_indices):
            if atom_name == "CB" and aa.upper() == "G":
                continue  # Glycine has no CB
            if atom_idx >= positions.shape[1]:
                continue

            x, y, z = positions[res_idx, atom_idx]
            lines.append(
                f"ATOM  {atom_serial:5d}  {atom_name:<3s} {res_name:>3s} A{res_num:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom_name[0]:>2s}"
            )
            atom_serial += 1

    lines.append("END")
    return "\n".join(lines) + "\n"


def fold_proteins_batch(
    sequences: List[Tuple[str, str]],
    output_dir: str,
    device: str = "cpu",
    max_length: int = 700,
) -> Dict[str, str]:
    """
    Fold multiple proteins, saving each as a PDB file.

    Args:
        sequences: list of (id, amino_acid_sequence) pairs
        output_dir: directory for PDB files
        device:     'cpu' or 'cuda'
        max_length: truncate sequences longer than this

    Returns:
        dict {sequence_id: pdb_file_path} for successful folds
    """
    os.makedirs(output_dir, exist_ok=True)
    results: Dict[str, str] = {}

    for seq_id, seq in sequences:
        if not seq or len(seq) < 10:
            continue
        pdb_path = os.path.join(output_dir, f"{_safe_filename(seq_id)}.pdb")
        result = fold_protein(seq, pdb_path, device=device, max_length=max_length)
        if result:
            results[seq_id] = result

    return results


def _safe_filename(name: str) -> str:
    """Sanitise a sequence ID for use as a filename."""
    return re.sub(r"[^\w\-.]", "_", name)[:200]


# ---------------------------------------------------------------------------
# GOI structure caching
# ---------------------------------------------------------------------------

def prefold_goi_structures(
    db_fasta: str,
    output_dir: str,
    device: str = "cpu",
    goi_prefix: str = "GOI_",
    max_length: int = 700,
) -> Dict[str, str]:
    """
    Read the initial database FASTA, extract GOI sequences, fold them with
    ESMFold, and save PDB files.

    Returns:
        dict {goi_id: pdb_path}
    """
    sequences: List[Tuple[str, str]] = []
    for header, clean_id, seq in parse_fasta(db_fasta):
        if clean_id.startswith(goi_prefix):
            sequences.append((clean_id, seq))

    if not sequences:
        logger.warning("No GOI sequences found in database for structural folding.")
        return {}

    logger.info(f"Folding {len(sequences)} GOI structure(s) with ESMFold ...")
    goi_pdb_dir = os.path.join(output_dir, "goi_structures")
    results = fold_proteins_batch(sequences, goi_pdb_dir, device=device, max_length=max_length)
    logger.info(f"GOI structures folded: {len(results)}/{len(sequences)} successful")
    return results


def save_structure_index(index: Dict[str, str], path: str) -> None:
    """Save a {id: pdb_path} mapping as a TSV file."""
    with open(path, "w") as f:
        for seq_id, pdb_path in index.items():
            f.write(f"{seq_id}\t{pdb_path}\n")


def load_structure_index(path: str) -> Dict[str, str]:
    """Load a structure index written by save_structure_index()."""
    index: Dict[str, str] = {}
    if not os.path.exists(path):
        return index
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                index[parts[0]] = parts[1]
    return index


# ---------------------------------------------------------------------------
# Foldseek comparison
# ---------------------------------------------------------------------------

def foldseek_search(
    query_pdb: str,
    target_pdbs: Dict[str, str],
    output_dir: str,
    tm_threshold: float = 0.3,
    threads: int = 1,
) -> List[Dict[str, Any]]:
    """
    Compare a query structure against target structures using Foldseek.

    Creates a Foldseek database from target PDBs, runs easy-search, and
    parses the tabular output.

    Args:
        query_pdb:    path to the query PDB file
        target_pdbs:  dict {id: pdb_path} for target structures
        output_dir:   working directory for Foldseek
        tm_threshold: minimum TM-score to report
        threads:      CPU threads for Foldseek

    Returns:
        list of dicts with Foldseek hit information
    """
    if not target_pdbs:
        return []

    os.makedirs(output_dir, exist_ok=True)

    # Foldseek easy-search: query PDB against a directory of target PDBs
    target_dir = os.path.join(output_dir, "targets")
    os.makedirs(target_dir, exist_ok=True)

    # Symlink or copy target PDBs into target directory
    id_to_filename: Dict[str, str] = {}
    for seq_id, pdb_path in target_pdbs.items():
        safe_name = _safe_filename(seq_id) + ".pdb"
        target_path = os.path.join(target_dir, safe_name)
        id_to_filename[safe_name.replace(".pdb", "")] = seq_id
        if not os.path.exists(target_path):
            try:
                os.symlink(os.path.abspath(pdb_path), target_path)
            except OSError:
                shutil.copy2(pdb_path, target_path)

    result_file = os.path.join(output_dir, "foldseek_results.tsv")
    tmp_dir = os.path.join(output_dir, "foldseek_tmp")

    # Foldseek easy-search with tabular output
    # Output columns: query, target, fident, alnlen, mismatch, gapopen,
    #                 qstart, qend, tstart, tend, evalue, bits, alntmscore
    cmd = [
        "foldseek", "easy-search",
        query_pdb,
        target_dir,
        result_file,
        tmp_dir,
        "--threads", str(threads),
        "--format-output",
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore",
        "-e", "10",  # permissive e-value; filter by TM-score instead
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if proc.returncode != 0:
            logger.debug(f"Foldseek exited {proc.returncode}: {proc.stderr[:300]}")
            return []
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        logger.debug(f"Foldseek unavailable or timed out: {exc}")
        return []

    if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
        return []

    # Parse results
    hits: List[Dict[str, Any]] = []
    with open(result_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 13:
                continue

            try:
                tm_score = float(fields[12])
            except (ValueError, IndexError):
                tm_score = 0.0

            if tm_score < tm_threshold:
                continue

            target_name = fields[1].replace(".pdb", "")
            original_id = id_to_filename.get(target_name, target_name)

            try:
                hits.append({
                    "query_structure": fields[0],
                    "target_id": original_id,
                    "fident": float(fields[2]),
                    "alnlen": int(fields[3]),
                    "qstart": int(fields[6]),
                    "qend": int(fields[7]),
                    "tstart": int(fields[8]),
                    "tend": int(fields[9]),
                    "evalue": float(fields[10]),
                    "bits": float(fields[11]),
                    "tm_score": tm_score,
                })
            except (ValueError, IndexError):
                continue

    hits.sort(key=lambda h: -h.get("tm_score", 0))
    return hits


# ---------------------------------------------------------------------------
# Region-level structural search (called from iterative_search_runner)
# ---------------------------------------------------------------------------

def structural_search_region(
    goi_structures: Dict[str, str],
    region_fasta: str,
    output_dir: str,
    tm_threshold: float = 0.3,
    device: str = "cpu",
    max_length: int = 700,
    threads: int = 1,
    predictor: str = "auto",
    augustus_species: str = "fly",
) -> List[Dict[str, Any]]:
    """
    Search for GOI structural homologs in a genomic region.

    1. Predict ORFs/genes (Augustus for eukaryotes, Prodigal for prokaryotes)
    2. Fold ORFs with ESMFold
    3. Compare against pre-folded GOI structures via Foldseek
    4. Return hits above TM-score threshold

    Returns list of hit dicts compatible with the augmented-search pipeline.
    """
    # Import ORF prediction from plm_search (shared utility)
    try:
        from plm_search import predict_orfs_prodigal
    except ImportError:
        logger.debug("Cannot import predict_orfs_prodigal from plm_search")
        return []

    orfs = predict_orfs_prodigal(
        region_fasta, output_dir,
        predictor=predictor, augustus_species=augustus_species,
    )
    if not orfs:
        return []

    # Fold ORF proteins
    orf_sequences = [(orf["id"], orf["seq"]) for orf in orfs]
    orf_pdb_dir = os.path.join(output_dir, "orf_structures")
    orf_structures = fold_proteins_batch(
        orf_sequences, orf_pdb_dir, device=device, max_length=max_length,
    )
    if not orf_structures:
        return []

    orf_by_id = {orf["id"]: orf for orf in orfs}

    # Compare each GOI structure against all ORF structures
    all_hits: List[Dict[str, Any]] = []
    for goi_id, goi_pdb in goi_structures.items():
        if not os.path.exists(goi_pdb):
            continue

        fs_output = os.path.join(output_dir, f"foldseek_{_safe_filename(goi_id)}")
        fs_hits = foldseek_search(
            query_pdb=goi_pdb,
            target_pdbs=orf_structures,
            output_dir=fs_output,
            tm_threshold=tm_threshold,
            threads=threads,
        )

        for hit in fs_hits:
            orf_id = hit["target_id"]
            orf = orf_by_id.get(orf_id)
            if not orf:
                continue

            all_hits.append({
                "query": goi_id,
                "chrom": "region_seq",
                "start": orf.get("start", 0),
                "end": orf.get("end", 0),
                "strand": orf.get("strand", "+"),
                "pident": hit["fident"] * 100.0,
                "identity": hit["fident"] * 100.0,
                "evalue": hit["evalue"],
                "bits": hit["bits"],
                "alnlen": hit["alnlen"],
                "qstart": hit["qstart"],
                "qend": hit["qend"],
                "method": "foldseek_structural",
                "structural_similarity": hit["tm_score"],
            })

    all_hits.sort(key=lambda h: -h.get("structural_similarity", 0))
    return all_hits


# ---------------------------------------------------------------------------
# Re-ranking helper (compute structural similarity for existing candidates)
# ---------------------------------------------------------------------------

def compute_candidate_structural_similarities(
    candidate_sequences: List[Tuple[str, str]],
    goi_structures: Dict[str, str],
    output_dir: str,
    device: str = "cpu",
    max_length: int = 700,
    threads: int = 1,
) -> Dict[str, float]:
    """
    Fold candidate protein models, compare against GOI structures, and
    return the best TM-score for each candidate.

    Args:
        candidate_sequences: [(model_id, protein_seq), ...]
        goi_structures:      pre-folded GOI PDB paths {goi_id: pdb_path}
        output_dir:          working directory

    Returns:
        {model_id: best_tm_score}
    """
    if not candidate_sequences or not goi_structures:
        return {}

    cand_pdb_dir = os.path.join(output_dir, "candidate_structures")
    cand_structures = fold_proteins_batch(
        candidate_sequences, cand_pdb_dir, device=device, max_length=max_length,
    )
    if not cand_structures:
        return {}

    best_scores: Dict[str, float] = {}
    for goi_id, goi_pdb in goi_structures.items():
        if not os.path.exists(goi_pdb):
            continue

        fs_output = os.path.join(output_dir, f"rerank_{_safe_filename(goi_id)}")
        fs_hits = foldseek_search(
            query_pdb=goi_pdb,
            target_pdbs=cand_structures,
            output_dir=fs_output,
            tm_threshold=0.0,  # Return all for ranking
            threads=threads,
        )

        for hit in fs_hits:
            tid = hit["target_id"]
            tm = hit["tm_score"]
            if tid not in best_scores or tm > best_scores[tid]:
                best_scores[tid] = tm

    return best_scores


# ---------------------------------------------------------------------------
# CLI for pre-folding GOI structures
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Pre-fold GOI structures with ESMFold for SynVoy Foldseek search"
    )
    parser.add_argument("--db_fasta", required=True, help="Initial database FASTA")
    parser.add_argument("--output_dir", required=True, help="Output directory for PDB files")
    parser.add_argument("--device", default="cpu", help="cpu or cuda")
    parser.add_argument("--max_length", type=int, default=700,
                        help="Max sequence length for ESMFold (VRAM safety)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    if not check_esmfold_available():
        logger.error(
            "Structural search requires PyTorch and HuggingFace Transformers.\n"
            "  Install with: pip install torch transformers"
        )
        sys.exit(1)

    if not check_foldseek_available():
        logger.error(
            "Structural search requires Foldseek binary on PATH.\n"
            "  Install with: conda install -c bioconda foldseek"
        )
        sys.exit(1)

    structures = prefold_goi_structures(
        args.db_fasta, args.output_dir, device=args.device,
        max_length=args.max_length,
    )

    index_path = os.path.join(args.output_dir, "goi_structure_index.tsv")
    save_structure_index(structures, index_path)
    logger.info(f"Structure index saved to {index_path}")


if __name__ == "__main__":
    main()
