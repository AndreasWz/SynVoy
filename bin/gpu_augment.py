#!/usr/bin/env python3
"""
gpu_augment.py — main-process GPU augmentation phase for SynVoy (PLM + Foldseek).

WHY THIS EXISTS
---------------
The PLM (ProtT5) and structural (ESMFold+Foldseek) layers used to run *inside*
the per-genome ``ProcessPoolExecutor`` workers of iterative_search_runner. Those
workers are forked, and once the parent has initialised CUDA (it does, to fold /
embed the GOI templates) a forked child CANNOT re-initialise CUDA:

    RuntimeError: Cannot re-initialize CUDA in forked subprocess.
                  To use CUDA with multiprocessing, you must use the 'spawn' ...

So every target fold/embed died and was swallowed at ``logger.debug`` — the
structural/PLM layer was a *silent no-op on GPU* (observed on LRZ job 5699498:
1501/1501 folds failed, 0 structural scores). See CLAUDE.md "Known Bugs".

This module moves all GPU inference into the MAIN process (one CUDA context, no
fork), where it is inherently safe, and runs it in two model-sequential phases:

    embed everything with ProtT5   ->  unload ProtT5 (free ~5.5 GB)
    fold  everything with ESMFold  ->  unload ESMFold
    compare folds with Foldseek (CPU)

Loading the two models sequentially rather than co-resident is what reclaims the
VRAM that previously forced ESMFold to *truncate* the query to 400 aa (dropping
oskar's C-terminal OSK domain). Combined with domain-preserving windowed folding
(``structural_search.fold_protein_windows``), no residue/domain is dropped.

It is FAIL-LOUD: if a layer is *required* but its deps/GPU are unavailable, or if
every fold/embed of a non-empty input fails, it raises — rather than quietly
producing empty results (the old failure mode).

TWO PRODUCTS
------------
  * re-rank scores   {candidate_id: {embedding_similarity, structural_similarity}}
    fed back into _classify_goi_evidence by iterative_search_runner to
    boost/rescue confidence of already-modelled GOI candidates.
  * discovery hits   ORFs that are structurally/embedding-similar to the GOI but
    invisible to sequence search.

Usable both as an imported library (from iterative_search_runner.main, same
process) and as a standalone CLI (for testing / manual runs).
"""

import argparse
import json
import logging
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

try:
    from sequence_utils import parse_fasta, setup_logging, read_json, write_json
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, setup_logging, read_json, write_json

logger = logging.getLogger(__name__)

GOI_PREFIX = "GOI_"


# ---------------------------------------------------------------------------
# Availability (thin wrappers so callers don't need both modules)
# ---------------------------------------------------------------------------

def _plm_mod():
    import plm_search
    return plm_search


def _struct_mod():
    import structural_search
    return structural_search


def plm_available() -> bool:
    try:
        return _plm_mod().check_plm_available()
    except Exception:
        return False


def structural_available() -> bool:
    try:
        return _struct_mod().check_structural_search_available()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Work items
# ---------------------------------------------------------------------------

class WorkItem:
    """One protein to score: a re-rank candidate or a discovery ORF."""

    __slots__ = ("kind", "wid", "seq", "genome", "meta")

    def __init__(self, kind: str, wid: str, seq: str,
                 genome: str = "", meta: Optional[Dict[str, Any]] = None):
        self.kind = kind          # "candidate" | "orf"
        self.wid = wid            # unique id (candidate mRNA id, or orf id)
        self.seq = seq
        self.genome = genome
        self.meta = meta or {}


def load_manifest(path: str) -> List[WorkItem]:
    """Load a JSONL manifest of work items written by the Phase-1 workers.

    Each line: {"kind","id","seq",["genome","chrom","start","end","strand",
    "block"]}. Malformed / empty-seq lines are skipped (logged at debug).
    """
    items: List[WorkItem] = []
    if not path or not os.path.exists(path):
        return items
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                logger.debug("manifest %s: bad JSON at line %d", path, line_no)
                continue
            seq = (rec.get("seq") or "").strip()
            wid = rec.get("id")
            kind = rec.get("kind", "candidate")
            if not wid or not seq:
                continue
            meta = {k: rec[k] for k in
                    ("chrom", "start", "end", "strand", "block") if k in rec}
            items.append(WorkItem(kind, str(wid), seq, str(rec.get("genome", "")), meta))
    return items


def load_goi_templates_from_fasta(db_fasta: str,
                                  goi_prefix: str = GOI_PREFIX) -> List[Tuple[str, str]]:
    """Extract (id, seq) for the GOI query template(s) from the initial_db FASTA."""
    out: List[Tuple[str, str]] = []
    for _hdr, clean_id, seq in parse_fasta(db_fasta):
        if clean_id.startswith(goi_prefix) and seq:
            out.append((clean_id, seq))
    return out


# ---------------------------------------------------------------------------
# Phase 1 — embeddings (ProtT5), then free the model
# ---------------------------------------------------------------------------

def _embed_all(goi_templates: List[Tuple[str, str]],
               items: List[WorkItem],
               device: str) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """Embed GOI templates + all work items with ProtT5; return best cosine of
    each item to any GOI template, plus diagnostics. Frees ProtT5 afterwards."""
    plm = _plm_mod()
    diag: Dict[str, Any] = {"embedded": 0, "failed": 0}
    best: Dict[str, float] = {}

    to_embed = [(w.wid, w.seq) for w in items]
    try:
        goi_embs = plm.embed_proteins(goi_templates, device=device)
        if not goi_embs:
            raise RuntimeError("ProtT5 produced no GOI template embeddings")
        item_embs = plm.embed_proteins(to_embed, device=device) if to_embed else {}
    finally:
        # Always release ProtT5 so ESMFold has the VRAM next, even on error.
        plm.unload_model()

    diag["embedded"] = len(item_embs)
    diag["failed"] = len(to_embed) - len(item_embs)
    if item_embs:
        best = plm.best_similarities(goi_embs, item_embs)  # {item_id: max cosine}
    return best, diag


# ---------------------------------------------------------------------------
# Phase 2 — structures (ESMFold windowed), then Foldseek (CPU)
# ---------------------------------------------------------------------------

def _fold_all(goi_templates: List[Tuple[str, str]],
              items: List[WorkItem],
              work_dir: str,
              device: str,
              window: int,
              overlap: int) -> Tuple[Dict[str, str], Dict[str, List[str]], Dict[str, Any]]:
    """Fold GOI templates + all items with ESMFold (domain-preserving windows).

    Returns (goi_structures {win_id: pdb}, item_structures {item_id: [pdb,...]},
    diagnostics). Frees ESMFold afterwards. Raises if a non-empty fold set
    yields zero structures (the silent-no-op guard)."""
    ss = _struct_mod()
    diag: Dict[str, Any] = {"folded_items": 0, "failed_items": 0, "goi_windows": 0}

    goi_structures: Dict[str, str] = {}
    item_structures: Dict[str, List[str]] = {}
    try:
        goi_dir = os.path.join(work_dir, "goi_structures")
        for gid, gseq in goi_templates:
            for win_id, pdb in ss.fold_protein_windows(
                    gid, gseq, goi_dir, device=device, window=window, overlap=overlap):
                goi_structures[win_id] = pdb
        diag["goi_windows"] = len(goi_structures)
        if goi_templates and not goi_structures:
            raise RuntimeError(
                "ESMFold folded 0/%d GOI templates — GPU folding is broken "
                "(check CUDA / VRAM). Refusing to report a silent no-op."
                % len(goi_templates))

        item_dir = os.path.join(work_dir, "item_structures")
        for w in items:
            folds = ss.fold_protein_windows(
                w.wid, w.seq, item_dir, device=device, window=window, overlap=overlap)
            if folds:
                item_structures[w.wid] = [p for _wid, p in folds]
        diag["folded_items"] = len(item_structures)
        diag["failed_items"] = len(items) - len(item_structures)
        if items and not item_structures:
            raise RuntimeError(
                "ESMFold folded 0/%d target proteins — GPU folding is broken "
                "(the classic CUDA-in-fork silent no-op). Refusing to continue "
                "as if structural search succeeded." % len(items))
    finally:
        ss.unload_esmfold()

    return goi_structures, item_structures, diag


def _foldseek_best_tm(goi_structures: Dict[str, str],
                      item_structures: Dict[str, List[str]],
                      work_dir: str,
                      threads: int) -> Dict[str, float]:
    """Best TM-score of each item (over its windows) vs any GOI structure."""
    ss = _struct_mod()
    best_tm: Dict[str, float] = {}
    if not goi_structures or not item_structures:
        return best_tm

    # Flatten all item windows into one Foldseek target DB per GOI query, then
    # take the max TM across GOI queries and windows for each item.
    fs_dir = os.path.join(work_dir, "foldseek")
    # target_pdbs: window_pdb_id -> path ; map window id back to item id
    win_to_item: Dict[str, str] = {}
    target_pdbs: Dict[str, str] = {}
    for item_id, pdbs in item_structures.items():
        for k, p in enumerate(pdbs):
            win_id = item_id if len(pdbs) == 1 else f"{item_id}__w{k}"
            target_pdbs[win_id] = p
            win_to_item[win_id] = item_id

    for goi_id, goi_pdb in sorted(goi_structures.items()):
        if not os.path.exists(goi_pdb):
            continue
        hits = ss.foldseek_search(
            query_pdb=goi_pdb, target_pdbs=target_pdbs,
            output_dir=os.path.join(fs_dir, ss._safe_filename(goi_id)),
            tm_threshold=0.0, threads=threads)
        for h in hits:
            item_id = win_to_item.get(h["target_id"], h["target_id"])
            tm = float(h.get("tm_score", 0.0))
            if item_id not in best_tm or tm > best_tm[item_id]:
                best_tm[item_id] = tm
    return best_tm


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run_augmentation(
    db_fasta: str,
    items: List[WorkItem],
    work_dir: str,
    *,
    device: str = "cpu",
    do_plm: bool = False,
    do_structural: bool = False,
    require: bool = True,
    plm_threshold: float = 0.5,
    structural_tm_threshold: float = 0.3,
    fold_window: int = 380,
    fold_overlap: int = 60,
    threads: int = 1,
) -> Dict[str, Any]:
    """Run the GPU augmentation phase over `items`.

    Sequential by model: embed everything with ProtT5, free it, then fold
    everything with ESMFold, free it, then Foldseek-compare on CPU.

    `require`: if True (default) and a requested layer's deps/GPU are missing,
    raise instead of silently skipping — this is the fail-loud contract.

    Returns dict: {"scores": {id: {embedding_similarity, structural_similarity}},
    "discovery_hits": [...], "diagnostics": {...}}.
    """
    os.makedirs(work_dir, exist_ok=True)
    goi_templates = load_goi_templates_from_fasta(db_fasta)
    if not goi_templates:
        raise RuntimeError(
            "gpu_augment: no GOI templates (prefix %r) in %s" % (GOI_PREFIX, db_fasta))

    diagnostics: Dict[str, Any] = {
        "n_items": len(items),
        "n_goi_templates": len(goi_templates),
        "device": device,
        "plm": {"requested": do_plm},
        "structural": {"requested": do_structural},
    }
    emb_best: Dict[str, float] = {}
    tm_best: Dict[str, float] = {}

    # ---- PLM embedding phase --------------------------------------------
    if do_plm:
        if not plm_available():
            msg = ("gpu_augment: PLM search requested but ProtT5/torch is "
                   "unavailable in this environment.")
            if require:
                raise RuntimeError(msg)
            logger.warning("%s Skipping PLM.", msg)
            diagnostics["plm"]["skipped"] = "unavailable"
        else:
            emb_best, emb_diag = _embed_all(goi_templates, items, device)
            diagnostics["plm"].update(emb_diag)

    # ---- Structural folding + Foldseek phase ----------------------------
    if do_structural:
        if not structural_available():
            msg = ("gpu_augment: structural search requested but "
                   "ESMFold/Foldseek is unavailable in this environment.")
            if require:
                raise RuntimeError(msg)
            logger.warning("%s Skipping structural.", msg)
            diagnostics["structural"]["skipped"] = "unavailable"
        else:
            goi_structs, item_structs, fold_diag = _fold_all(
                goi_templates, items, work_dir, device, fold_window, fold_overlap)
            diagnostics["structural"].update(fold_diag)
            tm_best = _foldseek_best_tm(goi_structs, item_structs, work_dir, threads)
            diagnostics["structural"]["scored"] = len(tm_best)

    # ---- Assemble scores + discovery hits -------------------------------
    scores: Dict[str, Dict[str, Optional[float]]] = {}
    for w in items:
        emb = emb_best.get(w.wid)
        tm = tm_best.get(w.wid)
        if emb is None and tm is None:
            continue
        scores[w.wid] = {
            "embedding_similarity": emb,
            "structural_similarity": tm,
            "kind": w.kind,
            "genome": w.genome,
        }

    discovery_hits: List[Dict[str, Any]] = []
    for w in items:
        if w.kind != "orf":
            continue
        emb = emb_best.get(w.wid)
        tm = tm_best.get(w.wid)
        emb_ok = emb is not None and emb >= plm_threshold
        tm_ok = tm is not None and tm >= structural_tm_threshold
        if not (emb_ok or tm_ok):
            continue
        hit = {
            "id": w.wid, "genome": w.genome,
            "embedding_similarity": emb, "structural_similarity": tm,
            "method": ("foldseek_structural" if tm_ok else "plm_embedding"),
        }
        hit.update(w.meta)
        discovery_hits.append(hit)

    diagnostics["n_scored"] = len(scores)
    diagnostics["n_discovery_hits"] = len(discovery_hits)
    return {"scores": scores, "discovery_hits": discovery_hits,
            "diagnostics": diagnostics}


# ---------------------------------------------------------------------------
# CLI (standalone / testing)
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="SynVoy main-process GPU augmentation (PLM + Foldseek).")
    ap.add_argument("--db_fasta", required=True,
                    help="initial_db FASTA (GOI_-prefixed templates)")
    ap.add_argument("--manifest", required=True,
                    help="JSONL manifest of work items (candidates / ORFs)")
    ap.add_argument("--work_dir", required=True)
    ap.add_argument("--out", required=True, help="output JSON path")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--enable_plm_search", action="store_true")
    ap.add_argument("--enable_structural_search", action="store_true")
    ap.add_argument("--require", action="store_true",
                    help="fail loud if a requested layer is unavailable")
    ap.add_argument("--plm_threshold", type=float, default=0.5)
    ap.add_argument("--structural_tm_threshold", type=float, default=0.3)
    ap.add_argument("--fold_window", type=int, default=380)
    ap.add_argument("--fold_overlap", type=int, default=60)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    setup_logging()
    items = load_manifest(args.manifest)
    result = run_augmentation(
        args.db_fasta, items, args.work_dir,
        device=args.device,
        do_plm=args.enable_plm_search,
        do_structural=args.enable_structural_search,
        require=args.require,
        plm_threshold=args.plm_threshold,
        structural_tm_threshold=args.structural_tm_threshold,
        fold_window=args.fold_window,
        fold_overlap=args.fold_overlap,
        threads=args.threads,
    )
    write_json(args.out, result)
    d = result["diagnostics"]
    logger.info("gpu_augment: %d scored, %d discovery hit(s) over %d item(s).",
                d.get("n_scored", 0), d.get("n_discovery_hits", 0), d.get("n_items", 0))


if __name__ == "__main__":
    main()
