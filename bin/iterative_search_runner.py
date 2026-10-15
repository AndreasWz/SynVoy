#!/usr/bin/env python3

import argparse
import csv
import subprocess
import os
import shutil
import glob
import re
import concurrent.futures
import uuid
import sys
import json
import logging
import contextlib
import time
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from urllib.parse import unquote

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def has_parasail_available() -> bool:
    """
    Fast runtime probe for parasail availability in the current Python env.
    """
    try:
        import parasail  # type: ignore
        return True
    except Exception:
        return False


class SmithWatermanUnavailable(RuntimeError):
    """Raised when SW is requested (auto mode) but parasail is unavailable and the
    degraded path was not explicitly allowed. Callers turn this into a loud exit."""


def resolve_smith_waterman_enabled(enable: bool, sw_method: str,
                                   parasail_available: bool,
                                   allow_missing: bool) -> bool:
    """Pure decision for the auto-mode SW/parasail guard (testable in isolation).

    SW is the load-bearing tier for divergent GOI recovery, so an auto-mode request that
    can't be honoured (no parasail) must FAIL LOUD rather than silently drop to
    mmseqs/tblastn — unless the caller explicitly opts into the degraded path.
    Returns the resolved ``enable_smith_waterman`` flag, or raises
    ``SmithWatermanUnavailable`` when it must fail loud.
    """
    if enable and sw_method == "auto" and not parasail_available:
        if allow_missing:
            return False
        raise SmithWatermanUnavailable("parasail unavailable for requested Smith-Waterman")
    return enable


def _tail_lines(text: str, n: int = 20) -> str:
    if not text:
        return ""
    return "\n".join(text.strip().splitlines()[-n:])


def _format_mmseqs_failure(proc: subprocess.CompletedProcess) -> str:
    """
    Return compact failure details from an mmseqs subprocess result.
    """
    stderr_tail = _tail_lines(proc.stderr or "", n=20)
    stdout_tail = _tail_lines(proc.stdout or "", n=20)
    return stderr_tail or stdout_tail or "no tool output captured"


def _is_mmseqs_resource_failure(details: str) -> bool:
    """
    Detect common low-memory / prefilter-failure signatures from MMseqs.
    """
    txt = (details or "").lower()
    patterns = (
        "prefilter died",
        "search step died",
        "search died",
        "out of memory",
        "cannot allocate memory",
        "std::bad_alloc",
        "killed",
    )
    return any(p in txt for p in patterns)


# Shared helpers
# Use our own sequence utilities (no BioPython dependency)
try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate,
        parse_gff_attributes as _parse_gff_attributes,
        str2bool,
        longest_collinear_run as _shared_longest_collinear_run,
        sw_align, sw_align_span, sw_bits,
    )
    from annotate_goi_exons import (annotate_exons_from_hit_list, annotate_models_from_hit_list,
                                    MINIPROT_AVAILABLE)
    import synvoy_scale
except ImportError:
    # Fallback if not in path - add bin directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate,
        parse_gff_attributes as _parse_gff_attributes,
        str2bool,
        longest_collinear_run as _shared_longest_collinear_run,
        sw_align, sw_align_span, sw_bits,
    )
    from annotate_goi_exons import (annotate_exons_from_hit_list, annotate_models_from_hit_list,
                                    MINIPROT_AVAILABLE)
    import synvoy_scale

# Import fragment utilities if available
try:
    from fragment_query import generate_fragments, parse_fragment_id, merge_fragment_hits
    FRAGMENT_SUPPORT = True
except ImportError:
    FRAGMENT_SUPPORT = False

# Import PLM (Protein Language Model) search if available
try:
    from plm_search import (
        check_plm_available,
        plm_search_region,
        compute_candidate_similarities,
        precompute_goi_embeddings,
        load_embeddings,
        embed_proteins,
    )
    PLM_IMPORT_OK = True
except Exception as _plm_import_err:  # ImportError, or broken deps (e.g. AttributeError from a damaged numpy)
    PLM_IMPORT_OK = False
    sys.stderr.write(
        f"[plm] optional PLM search unavailable, continuing without it: "
        f"{type(_plm_import_err).__name__}: {_plm_import_err}\n"
    )

    def check_plm_available():
        return False

# Import structural search (Foldseek + ESMFold) if available
try:
    from structural_search import (
        check_structural_search_available,
        check_esmfold_available,
        check_foldseek_available,
        structural_search_region,
        compute_candidate_structural_similarities,
        prefold_goi_structures,
        save_structure_index,
        load_structure_index,
    )
    STRUCTURAL_IMPORT_OK = True
except Exception as _struct_import_err:  # ImportError, or broken deps (e.g. AttributeError from a damaged numpy)
    STRUCTURAL_IMPORT_OK = False
    sys.stderr.write(
        f"[struct] optional structural search unavailable, continuing without it: "
        f"{type(_struct_import_err).__name__}: {_struct_import_err}\n"
    )

    def check_structural_search_available():
        return False

    def check_esmfold_available():
        return False

    def check_foldseek_available():
        return False

# Import the main-process GPU augmentation phase (PLM + Foldseek). The GPU work
# used to run inside the forked ProcessPoolExecutor workers, but a fork cannot
# re-initialise CUDA once the parent has ("Cannot re-initialize CUDA in forked
# subprocess"), so every fold/embed died silently on GPU (LRZ job 5699498:
# 1501/1501 folds failed, 0 structural scores). It now runs once in the main
# process after the wave loop — see run_gpu_augment_phase() and gpu_augment.py.
try:
    import gpu_augment
    GPU_AUGMENT_OK = True
except Exception as _gpu_aug_err:  # optional deps; mirror the guards above
    GPU_AUGMENT_OK = False
    sys.stderr.write(
        f"[gpu_augment] optional GPU augmentation module unavailable: "
        f"{type(_gpu_aug_err).__name__}: {_gpu_aug_err}\n"
    )

# GPU inference (ProtT5 / ESMFold) must NOT run inside the forked workers (see
# above). This kill-switch keeps the legacy in-worker code paths present but
# disabled; the PLM/structural layers run via run_gpu_augment_phase() instead.
# Flip to True only for reference/experiments on a CUDA-free (CPU) device.
GPU_IN_WORKER = False

# Classification thresholds for GOI confidence assignment.  Populated from
# CLI args in main() so that they can be tuned per-run via nextflow.config.
CLASSIFY_THRESHOLDS = {
    # exon_annotation -> HIGH
    "high_min_exons": 2,
    "high_min_identity": 50.0,
    "high_min_flanking": 2,
    # §1 (2026-07-01): synteny ORDER, not just count. When the block has >= this many
    # flanking genes, HIGH additionally requires the longest collinear run of home-ordered
    # flanking to reach it — so a GOI in a SCRAMBLED paralog neighbourhood (same flanking
    # count, wrong order — the biglycan-as-decorin signature) is demoted to MEDIUM instead
    # of crowned HIGH on identity+count alone. Waived when there are fewer than this many
    # flanking (can't judge order) or when collinearity is unavailable (no home_rank).
    "high_min_collinear": 3,
    # exon_annotation -> MEDIUM
    "medium_min_identity": 40.0,
    "medium_min_flanking": 1,
    "medium_min_qcov": 0.65,
    # fallback_hit_span -> MEDIUM (flanking-supported)
    "fallback_med_min_flanking": 2,
    "fallback_med_min_qcov": 0.75,
    "fallback_med_min_identity": 60.0,
    # fallback_hit_span -> MEDIUM (strong flanking context)
    "fallback_strong_min_flanking": 5,
    "fallback_strong_min_qcov": 0.25,
    "fallback_strong_min_identity": 35.0,
    # §1v / F6 — ABSOLUTE identity floor for the strong-flanking fallback clause.
    # That clause is (flanking >= 5) AND (qcov >= 0.25 OR identity >= 35): the qcov
    # arm has no identity requirement at all, so a hit at ANY identity reached MEDIUM
    # as long as it sat in a flanking-rich block. Because flanking_support is a
    # property of the BLOCK, every junk hit anywhere inside that block inherited it —
    # this is the mechanism behind the yeast STE2 run's 32 MEDIUM calls at 21-27 %
    # identity against 3 real orthologs, and Drosophila Defensin's 170 GOI calls for 3.
    #
    # Calibrated, not guessed: the real Tetramorium melittin call that this clause
    # legitimately rescues sits at 34.7 % (fallback_span_with_strong_flanking_support,
    # BlockFlankingSupport=7, QueryCoverage=0.686), while the yeast false positives run
    # 21-27 %. A floor of 30 keeps the former and rejects the latter with margin on
    # both sides. Set 0 to restore the old un-floored behaviour.
    "fallback_strong_min_identity_floor": 30.0,
    # tandem_copy -> MEDIUM vs LOW (50 since 2026-10-09; the short-peptide preset: 40)
    "tandem_min_identity": 50.0,
    # A2 (TODO_JUN): a tandem copy needs real query coverage too, not just a
    # high local-window identity (QW2 surfaces the true coverage).
    "tandem_min_qcov": 0.35,
    # model_status thresholds
    "fragment_max_qcov": 0.4,
    "complete_min_qcov": 0.7,
    # Context promotion (2026-10-13, _promote_context_supported_goi): a genome with no
    # HIGH/MEDIUM call for a query gets its best weak model lifted to MEDIUM when its block
    # keeps >= this many anchors in home order. 0 = off (the default). It was on from
    # 2026-10-14 for one day: reversed and shuffled home proteins searched in the same Candida
    # windows meet the rule in 30 of 181 decoy cases with a kept alignment, so it calls chance
    # alignments. 4 is the value the measurements in docs/PARAMETERS.md used.
    "context_min_collinear": 0,
    # Sanity floor for that promotion; the search's E-value filter is the real homology test.
    "context_min_identity": 15.0,
    # Query coverage a promoted model needs (never below fragment_max_qcov). With the
    # fragment bound of 0.4 the rule promoted seven models in the family and melittin
    # benchmark runs, none on a curated gene, six of them hit chains at 0.42-0.47, two of
    # those in bumblebees that lost melittin; 0.5 removes those six and one of 46 found
    # pairs of the pilot.
    "context_min_qcov": 0.5,
}

# Family-consistency config (Change A in TP53_IMPLEMENTATION_FIX_PLAN.md).
# `tokens` is a set of uppercase-alphanum tokens derived from the query gene name.
# When populated, every GOI annotation receives a GoiFamilyConsistent attribute.
# When `strict` is True, family-inconsistent fallback/rescued_exon/raw_hit calls
# are downgraded to LOW/ambiguous so they do not masquerade as probable GOI.
FAMILY_CONFIG: Dict[str, Any] = {
    "strict": False,
    "tokens": set(),
    "strict_evidence_types": {"fallback_hit_span", "rescued_exon", "raw_hit"},
}

# Weak-GOI emission filter. Drops LOW-confidence fallback/rescued_exon/raw_hit
# calls whose (identity/100)*qcov is below `min_id_x_qcov`. Motivated by TP53
# runs where ~70% of LOW fallback rows had id*qcov < 0.10 (pure noise) while
# genuine MEDIUM/probable melittin fallback hits never fall below 0.12.
# Default 0.05 drops obvious garbage without touching any validated toxin call.
# Set to 0.0 to disable entirely.
EMISSION_CONFIG: Dict[str, Any] = {
    "min_id_x_qcov": 0.0,
    "weak_evidence_types": {"fallback_hit_span", "rescued_exon", "raw_hit"},
    "_dropped_counter": 0,
}


def _median_best_flanking_identity(flanking_hits: List[Dict[str, Any]]) -> Tuple[Optional[float], int]:
    """§1g distance proxy: median of the best %identity per flanking gene.

    `flanking_hits` are the non-GOI (flanking-anchor) MMseqs hits for one target
    genome. We collapse to one value per flanking gene (its strongest hit) so a
    multi-HSP gene doesn't dominate, then take the median. High median => the
    target is a close relative (flanking genes align well); low median => distant.
    Returns (median_identity_or_None, n_flanking_genes).

    Keyed on the base gene id: with exon-level search (the default) every exon is its own
    query record (``gene|exon_N``), and keying on the record made this a median over exons,
    weighting long genes by their exon count and letting one 3-exon gene pass the minimum of
    three "genes" (logs said "over 103 genes" for about 20 flanking genes).
    """
    by_gene: Dict[str, float] = {}
    for h in flanking_hits:
        q = extract_base_gene_id(h.get('query', '') or '')
        if not q:
            continue
        pid = float(h.get('pident', 0.0) or 0.0)
        if q not in by_gene or pid > by_gene[q]:
            by_gene[q] = pid
    if not by_gene:
        return None, 0
    vals = sorted(by_gene.values())
    n = len(vals)
    median = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2.0
    return median, n


def _apply_distance_adaptive_thresholds(
    args, flanking_hits: List[Dict[str, Any]], genome_name: str
) -> Dict[str, Any]:
    """§1g: relax the GOI classify identity bars for phylogenetically distant targets.

    Mutates the module-global CLASSIFY_THRESHOLDS for the current genome. This is
    safe under the ProcessPoolExecutor model: genomes run in separate worker
    processes (no shared mutable state across genomes), and even when a worker is
    reused for the next genome this function re-derives both adapted keys from the
    immutable `args` baseline first, so there is no cross-genome drift. The only
    keys touched are high_min_identity / medium_min_identity.

    By construction the 'close' tier is a pure no-op (thresholds == baseline), so
    well-conserved targets — and therefore the close-relative calls in existing
    benchmarks — are unaffected. Returns a per-genome decision dict (for logging).
    """
    base_high = float(args.classify_high_min_identity)
    base_med = float(args.classify_medium_min_identity)
    # Reset to baseline first (idempotent across worker reuse).
    CLASSIFY_THRESHOLDS["high_min_identity"] = base_high
    CLASSIFY_THRESHOLDS["medium_min_identity"] = base_med

    info: Dict[str, Any] = {
        "median_flank_identity": None,
        "n_flanking": 0,
        "tier": "disabled" if getattr(args, "disable_distance_autotune", False) else "baseline",
        "relax": 0.0,
        "high_min_identity": base_high,
        "medium_min_identity": base_med,
        "manual_review": False,
    }
    if getattr(args, "disable_distance_autotune", False):
        return info

    median, n = _median_best_flanking_identity(flanking_hits)
    info["median_flank_identity"] = median
    info["n_flanking"] = n

    min_flank = int(getattr(args, "distance_autotune_min_flanking", 3))
    if median is None or n < min_flank:
        logger.info(
            f"[{genome_name}] §1g distance auto-tune: only {n} flanking gene(s) "
            f"(< {min_flank}); keeping baseline thresholds "
            f"(HIGH {base_high:.0f}, MEDIUM {base_med:.0f})."
        )
        info["tier"] = "insufficient_evidence"
        return info

    close_p = float(getattr(args, "distance_autotune_close_pct", 70.0))
    far_p = float(getattr(args, "distance_autotune_far_pct", 40.0))
    max_relax = float(getattr(args, "distance_autotune_max_relax", 10.0))

    if median >= close_p:
        relax, tier = 0.0, "close"
    elif median <= far_p:
        relax, tier = max_relax, "far"
    else:
        relax = max_relax * (close_p - median) / max(close_p - far_p, 1e-6)
        tier = "mid"

    # Floors: the relaxation never takes a threshold below its floor, and a baseline that is
    # already at or below the floor (short peptides: 40 / 30) is not relaxed at all.
    high_floor = float(getattr(args, "distance_autotune_high_floor", 45.0))
    med_floor = float(getattr(args, "distance_autotune_medium_floor", 30.0))
    high_adj = max(base_high - relax, min(base_high, high_floor), 1.0)
    med_adj = max(base_med - relax, min(base_med, med_floor), 1.0)
    CLASSIFY_THRESHOLDS["high_min_identity"] = high_adj
    CLASSIFY_THRESHOLDS["medium_min_identity"] = med_adj
    info.update(
        tier=tier, relax=relax, high_min_identity=high_adj,
        medium_min_identity=med_adj, manual_review=(tier == "far"),
    )

    if relax > 0:
        _held = [name for name, base, adj in (("HIGH", base_high, high_adj), ("MEDIUM", base_med, med_adj))
                 if adj > base - relax + 1e-9]
        logger.info(
            f"[{genome_name}] §1g distance auto-tune: median flanking identity "
            f"{median:.1f}% over {n} genes -> tier={tier}, relax={relax:.1f} pts "
            f"(HIGH {base_high:.0f}->{high_adj:.0f}, MEDIUM {base_med:.0f}->{med_adj:.0f})"
            + (f"  [{' and '.join(_held)} held at the floor]" if _held else "")
            + ("  [manual_review]" if tier == "far" else "")
        )
    else:
        logger.info(
            f"[{genome_name}] §1g distance auto-tune: median flanking identity "
            f"{median:.1f}% over {n} genes -> tier=close; thresholds unchanged "
            f"(HIGH {base_high:.0f}, MEDIUM {base_med:.0f})."
        )
    return info


def _should_skip_weak_goi_emission(
    evidence_type: str,
    identity: float,
    query_cov: Optional[float],
    flanking_support: int = 0,
    embedding_similarity: Optional[float] = None,
    structural_similarity: Optional[float] = None,
) -> bool:
    """Return True if this GOI emission is weak enough to drop before GFF write.

    Only fires on LOW-confidence fallback/rescued_exon/raw_hit. A hit that
    would classify as MEDIUM/HIGH (strong flanking, PLM rescue, structural
    rescue, or adequate qcov+identity) is always kept.
    """
    threshold = float(EMISSION_CONFIG.get("min_id_x_qcov") or 0.0)
    if threshold <= 0.0:
        return False
    if evidence_type not in EMISSION_CONFIG.get("weak_evidence_types", set()):
        return False
    qcov = float(query_cov or 0.0)
    ident = float(identity or 0.0)
    # Safety: if qcov is missing (0), the emission lacks coverage info — keep.
    # Prevents accidentally dropping tandem-copy-like rows that slip through.
    if qcov <= 0.0 or ident <= 0.0:
        return False
    conf, _cls, _reason = _classify_goi_evidence(
        evidence_type=evidence_type,
        identity=ident,
        exon_count=1,
        query_cov=qcov,
        flanking_support=flanking_support,
        embedding_similarity=embedding_similarity,
        structural_similarity=structural_similarity,
    )
    if conf != "LOW":
        return False
    return (ident / 100.0) * qcov < threshold


def _normalize_family_token(value: Any) -> str:
    if not value:
        return ""
    return re.sub(r"[^A-Za-z0-9]", "", str(value)).upper()


def _check_family_consistency(target_gene: str, target_product: str) -> Tuple[bool, str]:
    tokens = FAMILY_CONFIG.get("tokens") or set()
    if not tokens:
        return True, "family_tokens_unset"
    gene_norm = _normalize_family_token(target_gene)
    product_norm = _normalize_family_token(target_product)
    if not gene_norm and not product_norm:
        return False, "no_target_annotation"
    for tok in sorted(tokens):   # set: sorted so the reported match is stable
        if not tok:
            continue
        if gene_norm and tok in gene_norm:
            return True, f"matched_gene:{tok}"
        if product_norm and tok in product_norm:
            return True, f"matched_product:{tok}"
    return False, "no_family_match"


@contextlib.contextmanager
def maybe_quiet_streams(quiet: bool):
    """
    Suppress noisy third-party stdout/stderr (e.g. miniprot diagnostics)
    when running in low-noise mode.
    """
    if not quiet:
        yield
        return
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            yield

def normalize_coordinates(start: int, end: int) -> Tuple[int, int]:
    return min(start, end), max(start, end)

def filter_exon_hits(hits: List[Dict[str, Any]], query_len: int) -> List[Dict[str, Any]]:
    """Hits of an exon-level query that carry query coordinates.

    A coverage and a length limit were applied here until March 2026. They removed the
    short exon hits that miniprot and the hit chain need, so they were switched off;
    the two settings that fed them are gone (`min_exon_query_cov`) or documented for
    what they still do (`min_exon_alnlen`: last-resort rows only).
    """
    if not hits or query_len <= 0:
        return hits
    return [h for h in hits if h.get('qstart') is not None and h.get('qend') is not None]


def parse_hits(
    hits_file: str,
    min_identity: float,
    min_length: int,
    evalue_thresh: float,
    query_lengths: Optional[Dict[str, int]] = None,
    short_query_frac: float = 0.60,
    short_query_min: int = 12,
) -> List[Dict[str, Any]]:
    """
    Parse MMseqs2 hits and return a list of hit dictionaries.
    Filters by basic quality metrics.
    Preserves qstart/qend (query protein positions) and strand for exon annotation.
    """
    hits = []
    if not os.path.exists(hits_file):
        return hits

    skipped_lines = 0
    total_lines = 0
    first_error = None
    with open(hits_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            total_lines += 1
            try:
                # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
                # 0      1       2       3       4         5        6       7     8       9     10      11
                if len(parts) < 11:
                    skipped_lines += 1
                    continue

                pident = float(parts[2])
                alnlen = int(parts[3])
                evalue = float(parts[10])

                q_id = parts[0]
                q_len = None
                if query_lengths:
                    q_len = query_lengths.get(q_id)
                    if q_len is None:
                        q_base = extract_base_gene_id(q_id)
                        q_len = query_lengths.get(q_base)

                effective_min_length = min_length
                if q_len and q_len > 0:
                    adaptive_min = max(short_query_min, int(round(q_len * short_query_frac)))
                    effective_min_length = min(min_length, adaptive_min)

                if (evalue <= evalue_thresh and
                    pident >= min_identity and
                    alnlen >= effective_min_length):

                    t_start = int(parts[8])
                    t_end = int(parts[9])
                    start, end = normalize_coordinates(t_start, t_end)
                    # Convert 1-based mmseqs/BLAST coordinates to 0-based half-open
                    # for Python slicing: start-1 becomes 0-based, end stays (exclusive)
                    start -= 1
                    strand = '+' if t_start <= t_end else '-'

                    q_start = int(parts[6])
                    q_end = int(parts[7])

                    bits = float(parts[11]) if len(parts) > 11 else 0.0

                    hits.append({
                        'query': q_id,
                        'target': parts[1], # Chromosome/Scaffold
                        'chrom': parts[1],
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'qstart': min(q_start, q_end),
                        'qend': max(q_start, q_end),
                        'evalue': evalue,
                        'pident': pident,
                        'alnlen': alnlen,
                        'bits': bits
                    })
            except (ValueError, IndexError) as e:
                skipped_lines += 1
                if first_error is None:
                    first_error = f"line {total_lines}: {e}"
                continue

    if skipped_lines > 0:
        logger.warning(
            f"parse_hits({os.path.basename(hits_file)}): skipped {skipped_lines}/{total_lines} "
            f"malformed lines (first: {first_error})"
        )
    # A0 reproducibility: mmseqs/tblastn .m8 row order varies run-to-run under
    # --threads>1, which flips downstream STABLE-sort tie-breaks. Impose a total
    # order here so every consumer sees a run-independent hit sequence
    # (TODO_JUN §A0 task 2b).
    hits.sort(key=lambda h: (
        h['query'], h['chrom'], h['start'], h['end'], h['strand'], -h['bits']
    ))
    return hits

def create_locus_object(query_id, hits):
    chrom = hits[0]['chrom']
    start = min(h['start'] for h in hits)
    end = max(h['end'] for h in hits)
    return {
        'query': query_id,
        'chrom': chrom,
        'start': start,
        'end': end,
        'hits': hits
    }

def extract_base_gene_id(query_id: str) -> str:
    """
    Extract the base gene ID from a query ID that may contain exon info.
    
    Handles formats:
    - "gene-LOC726866" -> "gene-LOC726866"
    - "gene-LOC726866|exon_1" -> "gene-LOC726866"
    - "gene-LOC726866|var1" -> "gene-LOC726866"
    - "gene-LOC726866|GCA_xxx_MP000001" -> "gene-LOC726866"
    """
    # Split on | and take the first part
    parts = query_id.split('|')
    base_id = parts[0]
    
    # Also handle cases where gene ID itself contains underscores but not exon info
    # The exon suffix is always "|exon_N"
    
    return base_id


def collapse_goi_queries_by_parent(goi_queries: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """
    One representative GOI query per parent ID: the longest sequence.

    Ties are broken explicitly and never by input order. On equal length the user's
    own query (id == parent) beats a wavefront expansion model, which shares its
    parent (`GOI_Melt|Apis_florea_fna_b0_l1_exon_ann` -> `GOI_Melt`); after that the
    smallest id wins. Before this rule the input order came from a set of strings, so
    the hash seed decided, independently for every search block, whether the home
    melittin or the equal-length A. florea model was used.
    """
    def _rank(q):
        qid = q.get('id', '')
        return (len(q.get('seq', '')), qid == extract_base_gene_id(qid))

    by_parent: Dict[str, Dict[str, Any]] = {}
    for q in sorted(goi_queries, key=lambda q: q.get('id', '')):
        parent = extract_base_gene_id(q.get('id', ''))
        if not parent:
            continue
        if parent not in by_parent or _rank(q) > _rank(by_parent[parent]):
            by_parent[parent] = q
    return by_parent


# The second GOI query of a window (the model of a relative searched earlier) is searched
# and modelled under its own parent id, the home parent plus this tag, so that its hits
# are never mixed with those of the home protein: the two sequences have different
# coordinates. The tag is removed again before a model is written.
RELATIVE_QUERY_TAG = '~rel'
_MODEL_GENOME_RE = re.compile(r'^(?P<genome>.+)_b\d+_[lr]\d+')
_RELATIVE_SCORE_CACHE: Dict[Tuple[str, str], float] = {}


def build_home_goi_ids(initial_db_path) -> set:
    """Ids of the GOI entries of the pristine query set: the home protein and its exons.

    Everything else that later carries a GOI id in the expanding query set is a gene model
    of a genome searched earlier. Returns an empty set on failure, which switches the
    query rule back to "the longest sequence of each GOI"."""
    ids = set()
    try:
        for _header, clean_id, _seq in parse_fasta(initial_db_path):
            if is_goi_query_id(clean_id):
                ids.add(clean_id)
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning(f"build_home_goi_ids failed ({exc}); GOI queries fall back to the longest sequence.")
        return set()
    return ids


def _model_genome_tag(model_id: str) -> str:
    """The genome a wave model came from, read from its id
    (``GOI_X|<genome>_b<block>_l<locus>_<kind>``); '' for anything else."""
    suffix = model_id.split('|', 1)[1] if '|' in model_id else ''
    m = _MODEL_GENOME_RE.match(suffix)
    return m.group('genome') if m else ''


def _relative_score(home: Dict[str, Any], model: Dict[str, Any]) -> float:
    key = (home.get('id', ''), model.get('id', ''))
    if key not in _RELATIVE_SCORE_CACHE:
        try:
            _RELATIVE_SCORE_CACHE[key] = sw_align(home.get('seq', ''), model.get('seq', ''))[0]
        except RuntimeError:
            _RELATIVE_SCORE_CACHE[key] = float(len(model.get('seq', '')))
    return _RELATIVE_SCORE_CACHE[key]


def select_goi_queries(db_sequences: Dict[str, Dict[str, Any]], home_goi_ids,
                       use_relative: bool = True):
    """The GOI queries of one search window.

    For every GOI the HOME protein is always a query, and what is reported about a gene
    model (identity, query coverage, confidence) always refers to it. In addition the
    gene model of the relative searched most recently is a second query: it can find a
    gene the home protein no longer reaches. Until 2026-10-08 a window was searched with
    ONE sequence per GOI, the longest of the home protein and all earlier models, so a
    model three residues longer than the home protein replaced it for the rest of the run
    and the identity of every later call referred to that model.

    The relative: of the models in the query set, those of the genome that was added
    last (the genome closest in search order to the current one); among them the model
    that aligns best to the home protein. A model identical to the home protein is not
    added.

    Returns ``(queries, by_parent, alias_of)``. ``queries`` are ``{'id', 'seq'}`` records;
    a relative's record carries the id ``<parent>~rel``. ``by_parent`` maps each of these
    ids to its record. ``alias_of`` maps the relative's own id in the query set to its
    alias. Without ``home_goi_ids`` the old rule applies and ``alias_of`` is empty.
    """
    records: Dict[str, Dict[str, Any]] = {}
    for rec in db_sequences.values():
        rid = rec.get('id', '')
        if rid and rid not in records and is_full_length_goi_query_id(rid):
            records[rid] = rec

    if not home_goi_ids:
        by_parent = collapse_goi_queries_by_parent(list(records.values()))
        return list(by_parent.values()), by_parent, {}

    order = {rid: i for i, rid in enumerate(records)}
    by_parent: Dict[str, Dict[str, Any]] = {}
    models_by_parent: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for rid in sorted(records):
        rec = records[rid]
        parent = extract_base_gene_id(rid)
        if rid in home_goi_ids:
            best = by_parent.get(parent)
            rank = (rid == parent, len(rec.get('seq', '')))
            if best is None or rank > (best['id'] == parent, len(best.get('seq', ''))):
                by_parent[parent] = rec
        else:
            models_by_parent[parent].append(rec)

    # A GOI known only from models (no home entry) keeps the old rule.
    for parent, models in models_by_parent.items():
        if parent not in by_parent:
            by_parent.update(collapse_goi_queries_by_parent(models))

    alias_of: Dict[str, str] = {}
    queries = [by_parent[p] for p in sorted(by_parent)]
    if use_relative:
        for parent in sorted(models_by_parent):
            home = by_parent.get(parent)
            if home is None or home.get('id') not in home_goi_ids:
                continue
            models = [m for m in models_by_parent[parent]
                      if m.get('seq') and m.get('seq') != home.get('seq')]
            if not models:
                continue
            last = max(models, key=lambda m: order[m['id']])
            tag = _model_genome_tag(last['id'])
            pool = [m for m in models if _model_genome_tag(m['id']) == tag]
            best = max(sorted(pool, key=lambda m: m['id'], reverse=True),
                       key=lambda m: _relative_score(home, m))
            alias = parent + RELATIVE_QUERY_TAG
            alias_rec = {'id': alias, 'seq': best['seq'],
                         'header': best.get('header', ''), 'model_id': best['id']}
            by_parent[alias] = alias_rec
            queries.append(alias_rec)
            alias_of[best['id']] = alias
    return queries, by_parent, alias_of


def _parse_gff_attr_string(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for pair in text.split(';'):
        if '=' in pair:
            k, v = pair.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def rebase_relative_candidates(raw_candidates, cand_parents, by_parent):
    """Report every gene model found with a relative's model against the HOME protein.

    ``cand_parents[i]`` is the parent id (home parent, or its ``~rel`` alias) under which
    ``raw_candidates[i]`` was built. A model built with the relative is dropped where the
    home protein produced a model of its own on the same DNA; otherwise its identity and
    query coverage are measured again by aligning its protein to the home protein, it is
    classified with those numbers, and it gets the home parent's id. The relative that
    found it and the identity to that relative stay on the record as ``SearchQuery`` and
    ``SearchQueryIdentity``.
    """
    real_models = {'exon_annotation', 'tandem_copy', 'fallback_hit_span'}
    home_spans = defaultdict(list)
    for cand, parent in zip(raw_candidates, cand_parents):
        if parent.endswith(RELATIVE_QUERY_TAG) or not cand.get('gff'):
            continue
        ev = _parse_gff_attr_string(cand['gff'][0].split('\t')[-1]).get('EvidenceType', '')
        if ev in real_models:
            home_spans[parent].append((cand['start'], cand['end']))

    kept = []
    for cand, parent in zip(raw_candidates, cand_parents):
        if not parent.endswith(RELATIVE_QUERY_TAG):
            kept.append(cand)
            continue
        home_parent = parent[:-len(RELATIVE_QUERY_TAG)]
        home = by_parent.get(home_parent)
        relative = by_parent.get(parent, {})
        if home is None or not cand.get('gff'):
            continue
        c_len = max(1, cand['end'] - cand['start'])
        if any(max(0, min(cand['end'], e) - max(cand['start'], b)) > 0.5 * min(c_len, max(1, e - b))
               for b, e in home_spans.get(home_parent, [])):
            continue                      # the home protein has its own model here

        protein = cand['record'].get('seq', '')
        try:
            aln = sw_align_span(home.get('seq', ''), protein)
        except RuntimeError:
            continue                      # nothing to report against without an aligner
        cols = cand['gff'][0].split('\t')
        if len(cols) < 9:
            continue
        old = _parse_gff_attr_string(cols[8])
        evidence = old.get('EvidenceType', '')

        def _rename(text):
            text = text.replace(parent, home_parent)
            return re.sub(r'(_b\d+)_l(\d+)', r'\1_r\2', text)

        base = {k: _rename(old[k]) for k in ('ID', 'Name', 'SynVoy_Parent') if k in old}
        for k in ('Type', 'Exons', 'StartCodon', 'StopCodon'):
            if k in old:
                base[k] = old[k]
        collinear = old.get('BlockCollinearSupport')
        attrs = _goi_feature_attrs(
            base,
            evidence_type=evidence,
            identity=aln['identity'],
            exon_count=int(float(old.get('Exons', 1) or 1)),
            query_cov=aln['query_cov'],
            flanking_support=int(float(old.get('BlockFlankingSupport', 0) or 0)),
            collinear_support=int(float(collinear)) if collinear not in (None, '') else None,
        )
        attrs['SearchQuery'] = relative.get('model_id', parent)
        attrs['SearchQueryIdentity'] = old.get('Identity', '')
        native = {'label': old.get('TargetGene', ''), 'product': old.get('TargetProduct', ''),
                  'feature_id': old.get('TargetID', '')}
        cols[5] = f"{aln['identity']:.1f}"
        cols[8] = _compose_gff_attrs(attrs, native if any(native.values()) else None)
        new_gff = ['\t'.join(cols)]
        for line in cand['gff'][1:]:
            lc = line.split('\t')
            if len(lc) >= 9:
                lc[8] = _rename(lc[8])
                line = '\t'.join(lc)
            new_gff.append(line)
        cand = dict(cand)
        cand['gff'] = new_gff
        cand['score'] = aln['identity']
        record = dict(cand['record'])
        record['id'] = _rename(record.get('id', ''))
        record['description'] = (_rename(record.get('description', ''))
                                 + f" search_query:{relative.get('model_id', parent)}")
        cand['record'] = record
        kept.append(cand)
    return kept


def score_goi_candidates(raw_candidates, by_parent) -> int:
    """Bit score of every GOI candidate against the HOME protein of its GOI.

    Smith-Waterman of the candidate's protein (the model's translation, or the translated
    hits of a chain, a copy, a rescued exon or a raw hit) against the home protein
    (sequence_utils.sw_bits: BLOSUM62, BLAST gaps 11/1, lambda 0.267, K 0.041). Kept as
    ``cand['bits']`` for the choice among overlapping candidates and written as ``Bits=``
    on the candidate's mRNA/gene row. Until 2026-10-14 only hit chains carried a bit score
    (that of their best single hit) and the choice went by identity; on the Candida pilot
    bits picked the model on the curated gene in 29 of 31 genome pairs with several
    models, identity in 16. Run after rebase_relative_candidates, so every candidate's
    parent is its home GOI. Returns how many candidates were scored.
    """
    n = 0
    for cand in raw_candidates:
        gff = cand.get('gff') or []
        cols = gff[0].split('\t') if gff else []
        if len(cols) < 9:
            continue
        parent = _parse_gff_attr_string(cols[8]).get('SynVoy_Parent', '')
        home = by_parent.get(parent) or {}
        protein = (cand.get('record') or {}).get('seq', '')
        if not home.get('seq') or not protein:
            continue
        try:
            bits = sw_bits(home['seq'], protein)
        except RuntimeError:
            continue
        cand['bits'] = bits
        cols[8] = f"{cols[8]};Bits={bits:.1f}"
        gff[0] = '\t'.join(cols)
        n += 1
    return n


def goi_candidate_order(cand) -> Tuple[float, int, int, int, str]:
    """Order of GOI candidates for the greedy choice among overlapping ones: bit score
    against the home protein (score_goi_candidates), then the SHORTER span (a chain over
    two copies must not win over the copy it contains at the same score), then block
    flanking support, start, id. Until 2026-10-14: identity first, longer span first."""
    return (-float(cand.get('bits', 0.0) or 0.0),
            int(cand['end']) - int(cand['start']),
            -int(cand.get('flanking_support', 0) or 0),
            int(cand.get('start', 0) or 0),
            (cand.get('record') or {}).get('id', ''))


def is_goi_query_id(query_id: str) -> bool:
    """
    Identify GOI-derived queries in the expanding DB.

    GOI IDs produced by annotate_goi_exons.py are prefixed with `GOI_`
    (full protein, exon entries, tandem copies).
    """
    if not query_id:
        return False
    base_id = extract_base_gene_id(query_id)
    # Backward compatibility: older runs could emit bare "exon_N" IDs
    # for GOI-derived exons (missing GOI_ prefix).
    bare_legacy_goi_exon = bool(re.fullmatch(r'exon_\d+', query_id))
    return (
        query_id.startswith('GOI_') or
        base_id.startswith('GOI_') or
        query_id.startswith('GOI_copy_') or
        bare_legacy_goi_exon
    )


def is_full_length_goi_query_id(query_id: str) -> bool:
    """
    Full-length GOI-like queries used for final gene model inference.
    Excludes exon and synthetic fragment entries.
    """
    if not is_goi_query_id(query_id):
        return False
    if '|exon_' in query_id:
        return False
    if '|frag_' in query_id:
        return False
    return True


def split_hits_into_loci(hits: List[Dict[str, Any]], max_gap: int) -> List[List[Dict[str, Any]]]:
    """
    Split hits into local loci so distant paralogs are not forced into one model.
    """
    if not hits:
        return []

    sorted_hits = sorted(hits, key=lambda h: (h.get('chrom', ''), h.get('gstart', h.get('start', 0))))
    loci: List[List[Dict[str, Any]]] = []
    current: List[Dict[str, Any]] = [sorted_hits[0]]

    for h in sorted_hits[1:]:
        prev = current[-1]
        same_chrom = h.get('chrom') == prev.get('chrom')
        prev_end = prev.get('gend', prev.get('end', 0))
        curr_start = h.get('gstart', h.get('start', 0))

        if same_chrom and (curr_start - prev_end) <= max_gap:
            current.append(h)
        else:
            loci.append(current)
            current = [h]

    loci.append(current)
    return loci


def build_flanking_query_by_parent(
    db_sequences: Dict[str, Dict[str, str]],
    parent_ids: set
) -> Dict[str, Dict[str, str]]:
    """
    Build one representative flanking query sequence per parent gene.

    Flanking DB entries are often exon-level (`gene|exon_N`). For these, we
    reconstruct a pseudo full-length protein by concatenating exon proteins in
    transcript order (reverse exon order on minus strand).

    Returns dict mapping parent_id -> {id, seq, header, exon_offsets}.
    exon_offsets maps exon_num -> (start_offset, length) in the reconstructed
    protein, enabling downstream qstart/qend remapping.
    """
    if not parent_ids:
        return {}

    # Deduplicate aliased entries (base IDs and explicit IDs can point to same record).
    unique_records: Dict[str, Dict[str, str]] = {}
    for rec in db_sequences.values():
        rid = rec.get('id', '')
        if rid and rid not in unique_records:
            unique_records[rid] = rec

    by_parent: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for rid, rec in unique_records.items():
        parent = extract_base_gene_id(rid)
        if parent in parent_ids and not is_goi_query_id(rid):
            by_parent[parent].append(rec)

    result: Dict[str, Dict[str, str]] = {}
    for parent, recs in by_parent.items():
        # Prefer true full-length record if present.
        full_records = [
            r for r in recs
            if r.get('id', '') == parent and '|exon_' not in r.get('id', '')
        ]
        if full_records:
            best = max(full_records, key=lambda r: len(r.get('seq', '')))
            if best.get('seq'):
                result[parent] = {
                    'id': parent,
                    'seq': best.get('seq', ''),
                    'header': best.get('header', parent),
                    'exon_offsets': {},
                }
                continue

        exon_parts = []
        strand_hint = '+'
        for r in recs:
            rid = r.get('id', '')
            m = re.search(r'\|exon_(\d+)$', rid)
            if not m:
                continue
            exon_num = int(m.group(1))
            seq = r.get('seq', '')
            if not seq:
                continue
            header = r.get('header', '')
            sm = re.search(r'\bstrand=([+-])\b', header)
            if sm:
                strand_hint = sm.group(1)
            exon_parts.append((exon_num, seq))

        if exon_parts:
            exon_parts.sort(key=lambda x: x[0], reverse=(strand_hint == '-'))
            # Build offset map: exon_num -> (start_offset_in_reconstructed, length)
            exon_offsets = {}
            offset = 0
            for exon_num, seq in exon_parts:
                exon_offsets[exon_num] = (offset, len(seq))
                offset += len(seq)
            recon_seq = ''.join(seq for _, seq in exon_parts)
            if recon_seq:
                result[parent] = {
                    'id': parent,
                    'seq': recon_seq,
                    'header': f"{parent}|reconstructed_exons={len(exon_parts)}|strand={strand_hint}",
                    'exon_offsets': exon_offsets,
                }
                continue

        # Last fallback: use longest available record.
        best = max(recs, key=lambda r: len(r.get('seq', '')))
        if best.get('seq'):
            result[parent] = {
                'id': parent,
                'seq': best.get('seq', ''),
                'header': best.get('header', parent),
                'exon_offsets': {},
            }

    return result


def _select_parent_id(attrs: Dict[str, str], model_id: str) -> str:
    """
    Select a stable parent identifier from heterogeneous GFF attributes.

    Supports SynVoy, NCBI/Ensembl-style, and custom tags.
    """
    for key in [
        "SynVoy_Parent",
        "ParentProtein",
        "protein_id",
        "gene_id",
        "gene",
        "locus_tag",
        "Name",
    ]:
        value = (attrs.get(key) or "").strip()
        if not value:
            continue
        # If multiple values are provided, keep the first token.
        value = value.split(",")[0].strip()
        if value:
            return value
    return extract_base_gene_id(model_id)


def _is_generic_label(label: str) -> bool:
    """Heuristic for non-informative locus-tag style labels."""
    if not label:
        return True
    txt = str(label).strip()
    if not txt:
        return True
    if txt.startswith("gene-"):
        txt = txt[5:]
    if re.match(r"^[A-Za-z]+\d*_\d+$", txt):
        return True
    if re.match(r"^[A-Z]{2,6}\d{0,3}_\d+$", txt):
        return True
    if re.match(r"^LOC\d+$", txt, re.IGNORECASE):
        return True
    return False


def _safe_gff_value(value: Any) -> str:
    """Sanitize attribute values so GFF parsing remains robust."""
    txt = str(value) if value is not None else ""
    txt = txt.replace(";", ",").replace("\t", " ").replace("\n", " ").strip()
    return txt


def _compose_gff_attrs(base_attrs: Dict[str, Any],
                       native_annot: Optional[Dict[str, str]] = None) -> str:
    """
    Compose a GFF attribute field from base attrs plus optional native-annotation attrs.
    """
    merged: Dict[str, Any] = {}
    for key, value in base_attrs.items():
        if value is None:
            continue
        sval = _safe_gff_value(value)
        if sval:
            merged[key] = sval

    label = ""
    product = ""
    if native_annot:
        if native_annot.get("label"):
            merged["TargetGene"] = _safe_gff_value(native_annot["label"])
            label = native_annot["label"]
        if native_annot.get("product"):
            merged["TargetProduct"] = _safe_gff_value(native_annot["product"])
            product = native_annot["product"]
        if native_annot.get("feature_id"):
            merged["TargetID"] = _safe_gff_value(native_annot["feature_id"])

    if base_attrs.get("SynVoyRole") == "goi" and FAMILY_CONFIG.get("tokens"):
        consistent, reason = _check_family_consistency(label, product)
        merged["GoiFamilyConsistent"] = "true" if consistent else "false"
        merged["GoiFamilyReason"] = reason
        # Strict downgrade only fires on informative mismatches (annotated locus
        # whose name disagrees with family tokens). Unannotated regions — common
        # for novel toxin/venom discovery where orthologs are de novo — are NOT
        # downgraded: synteny + sequence evidence must stand on their own there.
        if (FAMILY_CONFIG.get("strict")
                and not consistent
                and reason == "no_family_match"):
            ev = base_attrs.get("EvidenceType") or base_attrs.get("Type") or ""
            if ev in FAMILY_CONFIG.get("strict_evidence_types", set()):
                merged["Confidence"] = "LOW"
                merged["GOIClass"] = "ambiguous_goi_family_member"
                prev = merged.get("InferenceReason", "")
                merged["InferenceReason"] = (
                    f"{prev}|strict_family_downgrade" if prev else "strict_family_downgrade"
                )

    return ";".join(f"{k}={v}" for k, v in merged.items())


def _format_attr_float(value: Optional[float], digits: int = 3) -> Optional[str]:
    if value is None:
        return None
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return None


def _synteny_context_label(flanking_support: int) -> str:
    flanking_support = max(0, int(flanking_support or 0))
    if flanking_support >= 3:
        return "strong_flanking_support"
    if flanking_support >= 1:
        return "weak_flanking_support"
    return "no_flanking_support"


def _model_status(
    query_cov: float,
    exon_count: int = 1,
    evidence_type: str = "",
) -> str:
    """
    Label a gene model's completeness independent of confidence.

    Returns one of:
      - ``complete``  — query coverage >= complete_min_qcov and multi-exon
        (or single-exon gene that is expected to be single-exon, e.g.
        tandem_copy)
      - ``partial``   — between fragment and complete
      - ``fragment``  — query coverage < fragment_max_qcov or rescued_exon /
        raw_hit
    """
    ct = CLASSIFY_THRESHOLDS
    qcov = float(query_cov or 0.0)
    exons = max(1, int(exon_count or 1))

    if evidence_type in {"rescued_exon", "raw_hit"}:
        return "fragment"
    if qcov < ct["fragment_max_qcov"]:
        return "fragment"
    if qcov >= ct["complete_min_qcov"] and (exons >= 2 or evidence_type == "tandem_copy"):
        return "complete"
    return "partial"


def _classify_goi_evidence(
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    flanking_support: int = 0,
    embedding_similarity: Optional[float] = None,
    structural_similarity: Optional[float] = None,
    collinear_support: Optional[int] = None,
    ambiguous_tier_enabled: bool = True,
) -> Tuple[str, str, str]:
    """
    Assign a conservative confidence/class label to GOI-derived candidates.

    The goal is not to prove orthology here, but to prevent fallback-heavy
    output from masquerading as confident GOI evidence downstream.

    embedding_similarity (0-1 cosine) from ProtT5 can BOOST confidence when
    the structural/functional signal is strong even if sequence identity is low.
    structural_similarity (0-1 TM-score) from ESMFold+Foldseek can BOOST
    confidence when 3D structure matches despite extreme sequence divergence.
    Neither signal ever reduces confidence — sequence methods already confirmed
    a match.

    ambiguous_tier_enabled adds a fourth verdict, AMBIGUOUS, for candidates promoted
    only by flanking support (see the fallback_hit_span branch). Pass False to restore
    the pre-2026-08-31 behaviour where those were labelled MEDIUM/probable_goi.
    """
    identity = float(identity or 0.0)
    exon_count = max(1, int(exon_count or 1))
    qcov = float(query_cov or 0.0)
    context = _synteny_context_label(flanking_support)

    ct = CLASSIFY_THRESHOLDS
    if ambiguous_tier_enabled is None:
        # Routed through CLASSIFY_THRESHOLDS rather than a call argument: that dict is
        # already the mechanism carrying classifier settings into every
        # ProcessPoolExecutor worker, so all five call sites pick it up unchanged.
        ambiguous_tier_enabled = bool(ct.get("ambiguous_tier_enabled", True))

    # §1 synteny-ORDER gate: the flanking are collinear enough to trust the neighbourhood as
    # this GOI's true (orthologous) locus, not a scrambled paralog neighbourhood. Waived when
    # collinearity is unavailable (collinear_support is None → no home_rank / feature off) or
    # when the block has too few flanking to judge order (< high_min_collinear).
    collinear_ok = (
        collinear_support is None
        or flanking_support < ct["high_min_collinear"]
        or collinear_support >= ct["high_min_collinear"]
    )

    if evidence_type == "tandem_copy":
        # A2 (TODO_JUN): MEDIUM requires BOTH identity AND real query coverage.
        # A short high-identity local window (e.g. a 7-aa micro-window or a 23-aa
        # hit on a 70-aa query — Vollenhovia/Polistes "79%"/"71%" overcalls) now
        # carries its true low qcov (QW2) and is demoted to LOW. Missing/zero
        # coverage is treated as failing the floor.
        if identity < ct["tandem_min_identity"]:
            return "LOW", "tandem_goi_copy", "goi_tandem_copy_low_identity"
        if qcov < ct.get("tandem_min_qcov", 0.0):
            return "LOW", "tandem_goi_copy", "goi_tandem_copy_low_coverage"
        # Context, like every other evidence type: MEDIUM needs the sequence to carry
        # the call (exon_annotation's qcov arm) or a strongly supported block (the
        # strong-flanking fallback bar). Tandem detection runs on raw window hits, and
        # a short query's best chance alignment in a window (~25-35 aa at 40-50 %) clears
        # 40 % / 0.35 by itself. Family benchmark 2026-09-11: 77 such MEDIUM calls in
        # the SECP and SPIN runs, all off-truth, all 16-33 bits, qcov <= 0.59, block
        # flanking 1-4 in 71 of them; the real tandem copies (KAZA, melittin bees) sit
        # at flanking 7-13 or qcov >= 0.6.
        if (flanking_support >= ct["fallback_strong_min_flanking"]
                or qcov >= ct["medium_min_qcov"]):
            return "MEDIUM", "tandem_goi_copy", "goi_tandem_copy_detected"
        return "LOW", "tandem_goi_copy", "goi_tandem_copy_weak_context"

    if evidence_type == "exon_annotation":
        # A model of one exon can be HIGH only where the home gene itself has one exon
        # (high_min_exons is then 1), and then it has to cover the query like a complete
        # gene: one exon is also what a single chance alignment looks like.
        one_exon_ok = exon_count >= 2 or qcov >= ct.get("high_single_exon_min_qcov", 0.0)
        # The flanking arm does not lift a model this pipeline itself calls a FRAGMENT
        # (ModelStatus, qcov < fragment_max_qcov): family benchmark 2026-09-11, single-exon
        # 6-9 %-coverage models in 2-flanking blocks were the MEDIUM calls left off-truth,
        # while every on-truth call had qcov >= 0.5. Unknown coverage is not low coverage:
        # the floor only applies when it was measured. The same holds for HIGH: a two-exon
        # model over 8 % of the query at 60 % was HIGH and seeded the next wave while the
        # same model at 45 % was LOW (bug hunt 2026-10-13, A4; once in about 380 HIGH rows).
        fragment = query_cov is not None and qcov < ct["fragment_max_qcov"]
        if (exon_count >= ct["high_min_exons"]
                and one_exon_ok
                and not fragment
                and identity >= ct["high_min_identity"]
                and flanking_support >= ct["high_min_flanking"]):
            if collinear_ok:
                if exon_count < 2:
                    return "HIGH", "confident_goi", "single_exon_model_as_at_home_with_flanking_support"
                return "HIGH", "confident_goi", "multi_exon_model_with_flanking_support"
            # §1: HIGH-quality model + enough flanking, but the flanking are NOT in conserved
            # home order → likely a paralog neighbourhood, not the true ortholog locus. Demote
            # to MEDIUM (still a probable GOI) rather than assert HIGH on identity+count alone.
            return "MEDIUM", "probable_goi", "high_identity_but_flanking_not_collinear"
        if identity >= ct["medium_min_identity"] and (
                (flanking_support >= ct["medium_min_flanking"] and not fragment)
                or qcov >= ct["medium_min_qcov"]):
            confidence = "MEDIUM"
            reason = "modeled_goi_with_partial_support"
            # PLM boost: MEDIUM → HIGH when embedding strongly agrees
            if (embedding_similarity is not None
                    and embedding_similarity >= ct.get("plm_high_threshold", 0.85)
                    and flanking_support >= ct["high_min_flanking"]):
                confidence = "HIGH"
                reason += "_plm_boosted"
            # Structural boost: MEDIUM → HIGH when TM-score strongly agrees
            if (confidence == "MEDIUM"
                    and structural_similarity is not None
                    and structural_similarity >= ct.get("structural_high_threshold", 0.7)
                    and flanking_support >= ct["high_min_flanking"]):
                confidence = "HIGH"
                reason += "_structural_boosted"
            return confidence, "probable_goi", reason
        # LOW from sequence — but PLM may rescue
        confidence = "LOW"
        goi_class = "ambiguous_goi_family_member"
        reason = ("fragment_model_with_flanking_only"
                  if fragment and identity >= ct["medium_min_identity"]
                  and flanking_support >= ct["medium_min_flanking"]
                  else "modeled_goi_but_family_context_is_weak")
        if (embedding_similarity is not None
                and embedding_similarity >= ct.get("plm_medium_threshold", 0.7)
                and flanking_support >= 1):
            confidence = "MEDIUM"
            goi_class = "probable_goi"
            reason = "modeled_goi_plm_rescued"
        # Structural rescue: LOW → MEDIUM when TM-score confirms structural match
        if (confidence == "LOW"
                and structural_similarity is not None
                and structural_similarity >= ct.get("structural_medium_threshold", 0.5)
                and flanking_support >= 1):
            confidence = "MEDIUM"
            goi_class = "probable_goi"
            reason = "modeled_goi_structural_rescued"
        return confidence, goi_class, reason

    if evidence_type == "fallback_hit_span":
        if (flanking_support >= ct["fallback_med_min_flanking"]
                and qcov >= ct["fallback_med_min_qcov"]
                and identity >= ct["fallback_med_min_identity"]):
            return "MEDIUM", "probable_goi", "fallback_span_supported_by_flanking_context"
        # Strong syntenic evidence upgrades even low-qcov hits.
        # Short secreted peptides (e.g. melittin, ~70 aa) can only reach ~45% qcov
        # because the signal peptide exon is too diverged to align; the flanking
        # gene context provides orthology evidence beyond per-alignment statistics.
        if (flanking_support >= ct["fallback_strong_min_flanking"]
                and identity >= ct.get("fallback_strong_min_identity_floor", 0.0)
                and (qcov >= ct["fallback_strong_min_qcov"]
                     or identity >= ct["fallback_strong_min_identity"])):
            # AMBIGUOUS, not MEDIUM, unless the sequence itself carries the call.
            #
            # This branch promotes on FLANKING SUPPORT, and flanking support cannot
            # distinguish "the gene is here" from "the gene was LOST from a locus whose
            # neighbourhood is still conserved" -- the second is exactly what gene loss
            # looks like. Measured on the 2026-08-31 benchmark (19 species, per-call
            # evidence dumped from the run):
            #
            #   gene PRESENT   identity 84.3 73.5 63.5 59.8 54.5 30.2   flanking 6-10
            #   gene LOST      identity 93.8 57.7 37.0 31.9 30.4 30.2   flanking 7-8
            #
            # The distributions overlap completely, the single highest-identity call in
            # the benchmark (93.8%, Bombus terrestris) is a FALSE POSITIVE, and flanking
            # support is 6-10 on BOTH sides. No threshold on this evidence separates them,
            # so asserting "probable_goi" here is asserting more than the data supports.
            #
            # AMBIGUOUS says the defensible thing instead: a candidate sits in the right
            # neighbourhood, and the sequence evidence does not establish that it is the
            # gene. It stays in the output and in the plots; it is excluded from the
            # ortholog counts. Independent confirmation (RBH against the home paralog
            # panel, PLM/structural rescue below, or curation) is what promotes it.
            if ambiguous_tier_enabled:
                return ("AMBIGUOUS", "syntenic_candidate_unconfirmed",
                        "fallback_span_with_strong_flanking_support")
            return "MEDIUM", "probable_goi", "fallback_span_with_strong_flanking_support"
        # PLM rescue for fallback hits
        if (embedding_similarity is not None
                and embedding_similarity >= ct.get("plm_medium_threshold", 0.7)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "fallback_span_plm_rescued"
        # Structural rescue for fallback hits
        if (structural_similarity is not None
                and structural_similarity >= ct.get("structural_medium_threshold", 0.5)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "fallback_span_structural_rescued"
        return "LOW", "ambiguous_goi_family_member", "fallback_span_only"

    if evidence_type == "rescued_exon":
        # PLM can rescue isolated exons when embedding strongly matches
        if (embedding_similarity is not None
                and embedding_similarity >= ct.get("plm_medium_threshold", 0.7)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "rescued_exon_plm_boosted"
        # Structural rescue for isolated exons
        if (structural_similarity is not None
                and structural_similarity >= ct.get("structural_medium_threshold", 0.5)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "rescued_exon_structural_boosted"
        return "LOW", "ambiguous_goi_family_member", "isolated_rescued_exon"

    if evidence_type == "raw_hit":
        return "LOW", "ambiguous_goi_family_member", "single_raw_hit_only"

    return "LOW", "ambiguous_goi_family_member", f"unclassified_{evidence_type or 'goi'}"


# Methods whose per-hit pident / bits / qcov are SYNTHETIC proxies (ProtT5 cosine
# or Foldseek fident/TM-score scaled to look like alignment statistics) rather
# than real sequence-alignment values. Feeding those into the identity/coverage
# branches of _classify_goi_evidence conflates embedding/structural similarity
# with sequence identity, which can manufacture a confident GOI call for a
# same-fold non-homolog. Such hits are routed through the dedicated embedding/
# structural rescue branches instead (which additionally require flanking support).
ML_DISCOVERY_METHODS = frozenset({"plm_embedding", "foldseek_structural"})


def _pooled_identity(hits: List[Dict[str, Any]]) -> float:
    """Percent identity of a model assembled from several local alignments.

    Identity is identical columns over aligned columns, so the per-hit values are
    weighted by alignment length: a 12-aa hit at 80 % and a 100-aa hit at 30 % are one
    35 % model, not a 55 % one (the unweighted mean this replaced let a short
    high-identity chance window lift a model over the identity bars). Falls back to
    the plain mean when any hit lacks an alignment length.
    """
    if not hits:
        return 0.0
    lens = [int(h.get("alnlen", 0) or 0) for h in hits]
    pids = [float(h.get("pident", 0) or 0.0) for h in hits]
    if all(n > 0 for n in lens):
        return sum(p * n for p, n in zip(pids, lens)) / sum(lens)
    return sum(pids) / len(pids)


def _resolve_fallback_signals(
    ordered_hits: List[Dict[str, Any]],
) -> Tuple[float, Optional[float], Optional[float], bool]:
    """Separate real-sequence identity from ML-discovery similarity for a
    fallback locus.

    Returns ``(seq_identity, embedding_similarity, structural_similarity,
    ml_only)``:

      * ``seq_identity`` — pooled (alignment-length-weighted, see _pooled_identity)
        pident of the NON-ML (real alignment) hits, or
        0.0 if the locus is supported only by ML-discovery seeds. This is the
        value that should be written as the model Identity and passed to the
        identity-based classifier — never the synthetic cosine/TM proxy.
      * ``embedding_similarity`` / ``structural_similarity`` — the best ML signal
        among contributing ML hits (or ``None``), surfaced so classification can
        use the principled rescue path.
      * ``ml_only`` — ``True`` when the locus has no real-sequence support, in
        which case the caller should also suppress the (synthetic) query
        coverage so MEDIUM can only be reached via the ML rescue branch.

    With no ML hits (PLM/structural search off), ``seq_identity`` is the pooled
    identity of all hits and ``ml_only`` is False.
    """
    seq_hits = [h for h in ordered_hits if h.get("method") not in ML_DISCOVERY_METHODS]
    ml_hits = [h for h in ordered_hits if h.get("method") in ML_DISCOVERY_METHODS]

    seq_identity = _pooled_identity(seq_hits) if seq_hits else 0.0

    embedding_similarity: Optional[float] = None
    structural_similarity: Optional[float] = None
    for h in ml_hits:
        es = h.get("embedding_similarity")
        if es is not None:
            embedding_similarity = float(es) if embedding_similarity is None else max(embedding_similarity, float(es))
        ss = h.get("structural_similarity")
        if ss is not None:
            structural_similarity = float(ss) if structural_similarity is None else max(structural_similarity, float(ss))

    return seq_identity, embedding_similarity, structural_similarity, (not seq_hits and bool(ml_hits))


CONTEXT_PROMOTABLE_EVIDENCE = frozenset({"exon_annotation", "fallback_hit_span"})


def _promote_context_supported_goi(gff_lines: List[str]) -> Tuple[List[str], int]:
    """Lift one weak GOI model per query to MEDIUM when gene ORDER, not identity, carries it.

    The identity bars (MEDIUM 40, relaxed to 30 at most by section 1g) assume the ortholog
    keeps a share of its sequence that scales with the distance of its neighbours. Fast
    genes do not: on the Candida pilot (2026-10-13) the curated orthologs that a protein
    search misses were modelled on the right gene at 18-38 % identity, their neighbours at
    58-85 %, and labelled LOW or AMBIGUOUS, so the report dropped 27 of 34 of them. Their
    blocks kept 4-15 anchors in home order. The models beside them were fragments, or sat
    in genomes that already had a confident call.

    Rule, per query (SynVoy_Parent) and genome, applied after every block is written:
      * no HIGH or MEDIUM call for this query in the genome, a tandem copy included
        (tandem copies are `gene` rows, every other model an `mRNA` row);
      * candidates: exon models and hit chains (not rescued exons, raw hits, tandem copies)
        labelled LOW or AMBIGUOUS, coverage >= context_min_qcov (and never a fragment),
        identity >= context_min_identity, BlockCollinearSupport >= context_min_collinear;
      * the candidate with the longest collinear run wins (then identity x coverage) and
        becomes MEDIUM / probable_goi with ContextPromoted=<old label>.
    The promoted call does not seed the next wave (_classify_goi_for_seed_and_tree), and it
    still goes through the report's paralog-ownership check, which has to catch a family
    relative sitting in the same conserved neighbourhood.

    Off by default (context_min_collinear 0; it was on, at 4, for one day from 2026-10-14:
    decoys -- reversed and shuffled home proteins searched in the same windows -- meet the
    rule in 30 of 181 Candida decoy cases with a kept alignment). Gene order cannot tell "the
    gene is here" from "the gene was lost and a relative or a chance match sits in its
    place", which is why the AMBIGUOUS tier exists and why the report names these calls
    apart (gene_order_only_goi). Measured: with a coverage floor of 0.4 the rule promoted
    seven models in the apyrase/DPP4/serine-protease/melittin runs, none on a curated gene,
    two of them the AMBIGUOUS calls in the bumblebees that lost melittin; 0.5 leaves one (a
    relative, which the ownership check then demotes). On 12 family-benchmark seeds not
    used to choose either value: 44 promoted, 2 demoted by the ownership check, 29 of the
    other 42 on a curated family gene (9 of the 13 wrong in one seed whose home locus is
    not where the curated genes are). SECP from no call at all to locus recall 0.50.
    """
    ct = CLASSIFY_THRESHOLDS
    min_coll = int(ct.get("context_min_collinear", 0) or 0)
    if min_coll <= 0 or not gff_lines:
        return gff_lines, 0
    min_ident = float(ct.get("context_min_identity", 0.0) or 0.0)
    min_qcov = max(float(ct.get("context_min_qcov", 0.0) or 0.0),
                   float(ct.get("fragment_max_qcov", 0.4)))

    confident: set = set()
    best: Dict[str, Tuple[Tuple[int, float], int]] = {}
    for i, line in enumerate(gff_lines):
        parts = line.split("\t")
        if len(parts) < 9 or parts[2] not in ("mRNA", "gene"):
            continue
        attrs = _parse_gff_attributes(parts[8])
        if attrs.get("SynVoyRole") != "goi":
            continue
        query = attrs.get("SynVoy_Parent") or attrs.get("Name", "")
        conf = (attrs.get("Confidence") or "").upper()
        if conf in ("HIGH", "MEDIUM"):
            confident.add(query)
            continue
        if (parts[2] != "mRNA" or conf not in ("LOW", "AMBIGUOUS")
                or attrs.get("EvidenceType") not in CONTEXT_PROMOTABLE_EVIDENCE):
            continue
        try:
            ident = float(attrs.get("Identity", 0) or 0)
            qcov = float(attrs.get("QueryCoverage", 0) or 0)
            coll = int(attrs.get("BlockCollinearSupport", 0) or 0)
        except ValueError:
            continue
        if qcov < min_qcov or ident < min_ident or coll < min_coll:
            continue
        key = (coll, ident * qcov)
        if query not in best or key > best[query][0]:
            best[query] = (key, i)

    out = list(gff_lines)
    n = 0
    for query, (_key, i) in best.items():
        if query in confident:
            continue
        parts = out[i].split("\t")
        old = _parse_gff_attributes(parts[8]).get("Confidence", "")
        new_attrs = []
        for pair in parts[8].split(";"):
            k = pair.split("=", 1)[0]
            if k == "Confidence":
                pair = "Confidence=MEDIUM"
            elif k == "GOIClass":
                pair = "GOIClass=probable_goi"
            elif k == "InferenceReason":
                pair = "InferenceReason=low_identity_in_conserved_gene_order"
            new_attrs.append(pair)
        new_attrs.append(f"ContextPromoted={old}")
        parts[8] = ";".join(new_attrs)
        out[i] = "\t".join(parts)
        n += 1
    return out, n


def _classify_flanking_evidence(
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
) -> Tuple[str, str]:
    identity = float(identity or 0.0)
    exon_count = max(1, int(exon_count or 1))
    qcov = float(query_cov or 0.0)

    if evidence_type in {"flanking_miniprot", "rearranged_flanking"}:
        if exon_count >= 2 or identity >= 55.0:
            return "HIGH", "modeled_flanking_gene"
        return "MEDIUM", "single_exon_flanking_model"

    if evidence_type in {"flanking_hit_span", "rearranged_flanking_fallback"}:
        if qcov >= 0.65 and identity >= 55.0:
            return "MEDIUM", "coarse_flanking_span_with_good_support"
        return "LOW", "coarse_flanking_span_only"

    return "MEDIUM", f"unclassified_{evidence_type or 'flanking'}"


def _goi_feature_attrs(
    base_attrs: Dict[str, Any],
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    flanking_support: int = 0,
    embedding_similarity: Optional[float] = None,
    structural_similarity: Optional[float] = None,
    collinear_support: Optional[int] = None,
    exons: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    confidence, goi_class, reason = _classify_goi_evidence(
        evidence_type=evidence_type,
        identity=identity,
        exon_count=exon_count,
        query_cov=query_cov,
        flanking_support=flanking_support,
        embedding_similarity=embedding_similarity,
        structural_similarity=structural_similarity,
        collinear_support=collinear_support,
    )
    attrs = dict(base_attrs)
    attrs.setdefault("Identity", f"{float(identity or 0.0):.1f}")
    attrs["SynVoyRole"] = "goi"
    attrs["EvidenceType"] = evidence_type
    attrs["Confidence"] = confidence
    attrs["GOIClass"] = goi_class
    attrs["ModelStatus"] = _model_status(query_cov, exon_count, evidence_type)
    # Whether the model carries the gene's own termini. The per-exon flags existed but only
    # ever reached CDS lines, so a model starting three codons inside the gene, with no
    # initiator Met, was indistinguishable from a complete one (family benchmark 2026-09-15:
    # 62.6 % of curated exons exact, the rest almost all terminal).
    if exons:
        attrs["StartCodon"] = "yes" if exons[0].get("has_start_codon") else "no"
        attrs["StopCodon"] = "yes" if exons[-1].get("has_stop_codon") else "no"
    attrs["SyntenyContext"] = _synteny_context_label(flanking_support)
    attrs["BlockFlankingSupport"] = str(max(0, int(flanking_support or 0)))
    if collinear_support is not None:
        attrs["BlockCollinearSupport"] = str(max(0, int(collinear_support)))
    attrs["InferenceReason"] = reason
    if query_cov is not None:
        attrs["QueryCoverage"] = _format_attr_float(query_cov)
    if exon_count:
        attrs.setdefault("Exons", str(int(exon_count)))
    if embedding_similarity is not None:
        attrs["EmbeddingSimilarity"] = f"{embedding_similarity:.3f}"
    if structural_similarity is not None:
        attrs["StructuralSimilarity"] = f"{structural_similarity:.3f}"
    return attrs


def _flanking_feature_attrs(
    base_attrs: Dict[str, Any],
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    context: Optional[str] = None,
) -> Dict[str, Any]:
    confidence, reason = _classify_flanking_evidence(
        evidence_type=evidence_type,
        identity=identity,
        exon_count=exon_count,
        query_cov=query_cov,
    )
    attrs = dict(base_attrs)
    attrs.setdefault("Identity", f"{float(identity or 0.0):.1f}")
    attrs["SynVoyRole"] = "flanking"
    attrs["EvidenceType"] = evidence_type
    attrs["Confidence"] = confidence
    attrs["ModelStatus"] = _model_status(query_cov, exon_count, evidence_type)
    attrs["InferenceReason"] = reason
    if context:
        attrs["SyntenyContext"] = context
    if query_cov is not None:
        attrs["QueryCoverage"] = _format_attr_float(query_cov)
    if exon_count:
        attrs.setdefault("Exons", str(int(exon_count)))
    return attrs


def find_native_annotation_path(genome_path: str) -> Optional[str]:
    """
    Locate the native annotation file next to a genome FASTA.
    Supports .gff/.gff3 (+ gz variants).
    """
    p = Path(genome_path)
    candidates = [
        p.with_suffix(".gff"),
        p.with_suffix(".gff3"),
        Path(str(p) + ".gff"),
        Path(str(p) + ".gff3"),
        p.with_suffix(".gff.gz"),
        p.with_suffix(".gff3.gz"),
        Path(str(p) + ".gff.gz"),
        Path(str(p) + ".gff3.gz"),
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


def load_native_annotation_index(gff_path: Optional[str]) -> Dict[str, Dict[str, Any]]:
    """
    Build a per-chromosome interval index from native target annotations.
    """
    if not gff_path or not os.path.exists(gff_path):
        return {}

    entries_by_chrom: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    accepted_types = {"gene", "mRNA", "mrna", "transcript"}

    for feat in parse_gff(gff_path):
        ftype = feat.get("type")
        if ftype not in accepted_types:
            continue

        attrs = feat.get("attributes", {}) or {}
        label = (
            attrs.get("gene")
            or attrs.get("Name")
            or attrs.get("locus_tag")
            or attrs.get("gene_id")
            or attrs.get("ID")
            or ""
        )
        product = attrs.get("product", "") or ""
        feature_id = attrs.get("ID") or attrs.get("locus_tag") or ""

        label = unquote(str(label)).strip()
        product = unquote(str(product)).strip()
        feature_id = unquote(str(feature_id)).strip()

        if not label and not product and not feature_id:
            continue

        try:
            s = int(feat.get("start", 0))
            e = int(feat.get("end", 0))
        except (ValueError, TypeError):
            continue
        if e < s:
            s, e = e, s

        entries_by_chrom[str(feat.get("seqid", ""))].append({
            "start": s,
            "end": e,
            "strand": feat.get("strand", "."),
            "label": label,
            "product": product,
            "feature_id": feature_id,
        })

    index: Dict[str, Dict[str, Any]] = {}
    for chrom, entries in entries_by_chrom.items():
        if not chrom:
            continue
        entries.sort(key=lambda x: (x["start"], x["end"]))
        index[chrom] = {
            "entries": entries,
            "starts": [e["start"] for e in entries],
        }
    return index


def lookup_native_annotation(native_index: Dict[str, Dict[str, Any]],
                             chrom: str,
                             start: int,
                             end: int,
                             strand: Optional[str] = None,
                             max_distance: int = 5000) -> Optional[Dict[str, str]]:
    """
    Find the best native annotation overlapping (or very near) an interval.
    Coordinates are expected as 1-based inclusive.
    """
    if not native_index or not chrom:
        return None
    chrom_data = native_index.get(chrom)
    if not chrom_data:
        return None

    s, e = int(start), int(end)
    if e < s:
        s, e = e, s

    entries = chrom_data["entries"]
    starts = chrom_data["starts"]
    i = bisect_right(starts, e + max_distance)
    if i <= 0:
        return None

    best = None
    best_score = None
    j = i - 1
    while j >= 0:
        ent = entries[j]
        if ent["end"] < s - max_distance:
            break
        ov = max(0, min(e, ent["end"]) - max(s, ent["start"]) + 1)
        if ov > 0:
            dist = 0
        else:
            dist = min(abs(s - ent["end"]), abs(e - ent["start"]))
            if dist > max_distance:
                j -= 1
                continue

        strand_score = 1 if strand and strand in {"+", "-"} and ent.get("strand") == strand else 0
        label_quality = 0 if _is_generic_label(ent.get("label", "")) else 1
        score = (
            1 if ov > 0 else 0,
            ov,
            -dist,
            strand_score,
            label_quality,
            (ent["end"] - ent["start"] + 1),
        )
        if best_score is None or score > best_score:
            best_score = score
            best = ent
        j -= 1

    if not best:
        return None

    return {
        "label": best.get("label", ""),
        "product": best.get("product", ""),
        "feature_id": best.get("feature_id", ""),
    }


def _flanking_chain_consistent(cds_intervals: List[Tuple[int, int]], strand: str) -> bool:
    """
    Validate exon chain direction from emitted CDS order for one transcript.

    '+' strand expects left-to-right exon centers.
    '-' strand expects right-to-left exon centers.
    """
    if len(cds_intervals) <= 1:
        return True

    centers = [0.5 * (s + e) for s, e in cds_intervals]
    if strand == "-":
        return all(centers[i] >= centers[i + 1] for i in range(len(centers) - 1))
    return all(centers[i] <= centers[i + 1] for i in range(len(centers) - 1))


def _split_models_into_loci(
    model_infos: List[Dict[str, Any]],
    locus_gap_bp: int = 50000
) -> List[List[Dict[str, Any]]]:
    """
    Split models into genomic loci for a single parent ID.

    Locus definition:
    - same chromosome
    - same strand
    - model start within `locus_gap_bp` of current locus end
    """
    if not model_infos:
        return []

    ordered = sorted(
        model_infos,
        key=lambda m: (
            m.get("chrom", ""),
            m.get("strand", "+"),
            m.get("start", 0),
            m.get("end", 0),
        ),
    )

    loci: List[List[Dict[str, Any]]] = []
    current: List[Dict[str, Any]] = []
    current_chrom = ""
    current_strand = "+"
    current_end = -1

    for info in ordered:
        chrom = info.get("chrom", "")
        strand = info.get("strand", "+")
        start = info.get("start", 0)
        end = info.get("end", 0)

        if not current:
            current = [info]
            current_chrom = chrom
            current_strand = strand
            current_end = end
            continue

        same_locus = (
            chrom == current_chrom and
            strand == current_strand and
            start <= (current_end + locus_gap_bp)
        )
        if same_locus:
            current.append(info)
            current_end = max(current_end, end)
        else:
            loci.append(current)
            current = [info]
            current_chrom = chrom
            current_strand = strand
            current_end = end

    if current:
        loci.append(current)
    return loci


# Consecutive exons of one chained model may overlap by at most this much genome
# (10 codons: local alignments bleed a few residues past a splice site).
CHAIN_MAX_EXON_OVERLAP_BP = 30


def _query_union_coverage(hits: List[Dict[str, Any]], query_len: int) -> float:
    """Fraction of the query inside the UNION of the hits' aligned intervals.

    The span from the lowest to the highest aligned query position is not coverage:
    family benchmark 2026-09-11, two 18- and 29-aa hits from opposite ends of the
    736-aa HYAL query, 33 kb apart, spanned 0.90 of it while aligning 0.06 -- and the
    fallback classifier's qcov >= 0.75 arm called that MEDIUM.
    """
    intervals = sorted(
        (min(int(h.get("qstart", 0)), int(h.get("qend", 0))),
         max(int(h.get("qstart", 0)), int(h.get("qend", 0))))
        for h in hits
    )
    covered, cur_s, cur_e = 0, None, None
    for s, e in intervals:
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                covered += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        covered += cur_e - cur_s + 1
    return min(1.0, covered / max(1, int(query_len)))


def _longest_monotonic_query_chain(
    ordered_hits: List[Dict[str, Any]],
    strand: str,
    max_gap: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Keep the longest genomic-order/query-order-consistent hit chain.

    `ordered_hits` must already be sorted by genomic start.
    For '+' strands, query centers should increase with genomic position.
    For '-' strands, query centers should decrease with genomic position.

    When ``max_gap`` is set (the max-intron distance), two hits are only chained
    as consecutive exons of one gene if the genomic gap between them is
    <= ``max_gap``. This stops two far-apart spurious hits — e.g. two ~15-aa
    melittin matches sitting 60 kb apart on *different* genes — from being
    stitched into a single multi-exon mRNA whose span is almost entirely empty
    intron (the "two seqs far apart called one gene" failure). If the gap bound
    leaves only single-hit chains, the single best hit (highest bitscore, then
    aligned length) is returned so a genuine compact locus is still emitted.
    """
    n = len(ordered_hits)
    if n <= 1:
        return ordered_hits

    qcenters = []
    for h in ordered_hits:
        qs = h.get("qstart", 0)
        qe = h.get("qend", 0)
        qcenters.append(0.5 * (min(qs, qe) + max(qs, qe)))

    def _quality(idx: int):
        h = ordered_hits[idx]
        return (float(h.get("bits", 0.0)), int(h.get("alnlen", 0)))

    dp = [1] * n
    prev = [-1] * n
    best_idx = 0

    for i in range(n):
        for j in range(i):
            ok = qcenters[j] <= qcenters[i] if strand != "-" else qcenters[j] >= qcenters[i]
            if ok and max_gap is not None:
                # Hits are sorted by genomic start, so a j->i transition makes j
                # the immediate upstream exon of i; bound the intron between them.
                gap = int(ordered_hits[i].get("gstart", 0)) - int(ordered_hits[j].get("gend", 0))
                if gap > max_gap:
                    ok = False
            if ok:
                # Two exons of one gene never share genome. Hits overlapping by more
                # than an alignment end's bleed are two alignments of ONE target
                # segment to different parts of the query, and chaining them fakes
                # coverage: family benchmark 2026-09-11, a 40-aa HYAL window aligned
                # to both ends of the 736-aa query read as qcov 0.92 -> MEDIUM.
                overlap = int(ordered_hits[j].get("gend", 0)) - int(ordered_hits[i].get("gstart", 0))
                if overlap > CHAIN_MAX_EXON_OVERLAP_BP:
                    ok = False
            if ok and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j
        # Prefer the longest chain; break ties toward the higher-quality
        # endpoint so a broken (all length-1) chain returns the best single hit.
        if dp[i] > dp[best_idx] or (dp[i] == dp[best_idx] and _quality(i) > _quality(best_idx)):
            best_idx = i

    chain_idx = []
    cur = best_idx
    while cur != -1:
        chain_idx.append(cur)
        cur = prev[cur]
    chain_idx.reverse()
    return [ordered_hits[i] for i in chain_idx]


def deduplicate_flanking_models(
    gene_records: List[Dict[str, Any]],
    gff_lines: List[str],
    genome_name: str = "",
    locus_gap_bp: int = 50000,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Keep one best flanking model per parent ID per genomic locus.

    GOI annotations are untouched. Flanking best model is chosen by:
    1) strand-consistent exon chain
    2) miniprot model over hit-span fallback
    3) more exons
    4) longer total CDS
    5) longer transcript span
    6) higher identity
    """
    if not gff_lines:
        return gene_records, gff_lines

    flanking_sources = {"flanking_annotation", "flanking_hits"}
    models: Dict[str, Dict[str, Any]] = {}
    cds_by_model: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    for line in gff_lines:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        source = parts[1]
        ftype = parts[2]
        if source not in flanking_sources:
            continue
        attrs = _parse_gff_attributes(parts[8])

        if ftype == "mRNA":
            model_id = attrs.get("ID")
            if not model_id:
                continue
            parent_id = _select_parent_id(attrs, model_id)
            try:
                identity = float(parts[5]) if parts[5] != "." else 0.0
            except (ValueError, TypeError):
                identity = 0.0
            try:
                start = int(parts[3])
                end = int(parts[4])
            except (ValueError, TypeError):
                continue
            models[model_id] = {
                "id": model_id,
                "parent": parent_id,
                "source": source,
                "chrom": parts[0],
                "strand": parts[6],
                "start": start,
                "end": end,
                "identity": identity,
            }
        elif ftype == "CDS":
            parent_model = attrs.get("Parent")
            if not parent_model:
                continue
            try:
                cds_by_model[parent_model].append((int(parts[3]), int(parts[4])))
            except (ValueError, TypeError):
                continue

    if not models:
        return gene_records, gff_lines

    models_by_parent: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for info in models.values():
        models_by_parent[info["parent"]].append(info)

    keep_model_ids = set()
    for parent_id, parent_models in models_by_parent.items():
        loci = _split_models_into_loci(parent_models, locus_gap_bp=max(0, int(locus_gap_bp)))
        for locus_models in loci:
            best = None
            best_score = None
            for info in locus_models:
                model_id = info["id"]
                cds_intervals = cds_by_model.get(model_id, [])
                exon_count = len(cds_intervals)
                cds_total_nt = sum(max(0, e - s + 1) for s, e in cds_intervals)
                span_nt = max(0, info["end"] - info["start"] + 1)
                consistent = _flanking_chain_consistent(cds_intervals, info.get("strand", "+"))
                source_rank = 1 if info.get("source") == "flanking_annotation" else 0

                score = (
                    1 if consistent else 0,
                    source_rank,
                    exon_count,
                    cds_total_nt,
                    span_nt,
                    info.get("identity", 0.0),
                )
                if best_score is None or score > best_score:
                    best_score = score
                    best = model_id

            if best:
                keep_model_ids.add(best)

    if len(keep_model_ids) == len(models):
        return gene_records, gff_lines

    filtered_gff: List[str] = []
    for line in gff_lines:
        parts = line.split("\t")
        if len(parts) < 9:
            filtered_gff.append(line)
            continue
        source = parts[1]
        ftype = parts[2]
        if source not in flanking_sources:
            filtered_gff.append(line)
            continue
        attrs = _parse_gff_attributes(parts[8])
        if ftype == "mRNA":
            model_id = attrs.get("ID", "")
            if model_id in keep_model_ids:
                filtered_gff.append(line)
        elif ftype == "CDS":
            model_id = attrs.get("Parent", "")
            if model_id in keep_model_ids:
                filtered_gff.append(line)
        else:
            filtered_gff.append(line)

    filtered_gene_records = []
    for rec in gene_records:
        rec_id = rec.get("id", "")
        if rec_id in models and rec_id not in keep_model_ids:
            continue
        filtered_gene_records.append(rec)

    dropped = len(models) - len(keep_model_ids)
    if dropped > 0:
        gname = genome_name or "genome"
        logger.info(
            f"[{gname}] Flanking model dedup (per-locus): kept {len(keep_model_ids)}/{len(models)} "
            f"models across {len(models_by_parent)} parent IDs (dropped {dropped} duplicates)."
        )

    return filtered_gene_records, filtered_gff


# ---------------------------------------------------------------------------
# Collinearity-aware seed placement (docs/TODO.md §1m — proper fix for the
# "fumble", replacing the §17 hull rescue's after-the-fact patch).
#
# The home flanking genes have a fixed order on the home chromosome. A genuine
# syntenic block reproduces that order (forward or reversed) in the target; a
# coincidental cluster of paralog cross-hits does not. Two consequences:
#   (#1) gap-bridging: two same-chromosome clusters separated by a large gap are
#        kept in ONE block when the flanking continue the same home-order run
#        across the gap — so a GOI whose immediate neighbours are rearranged
#        (and therefore sits in a flanking gap, e.g. cow decorin at chr5:21 Mb)
#        is inside the search window from the start, no rescue pass needed.
#   (#3) scoring: blocks are ranked by the length of their collinear run, not by
#        raw gene count, so the true ordered block outranks a scrambled
#        paralog-neighbourhood cluster (the source of the biglycan false call).
#
# `home_rank` maps a flanking gene's base ID -> its ordinal position in home
# genomic order (built once from the pristine initial_db; see build_home_rank).
# When `home_rank` is falsy every helper degrades to the legacy behaviour, so
# importers/tests that call identify_synteny_blocks without it are unaffected.
# ---------------------------------------------------------------------------
def _longest_collinear_run(ranks):
    """Thin alias — the implementation lives in sequence_utils so the region scorer
    (cluster_grs.py) measures collinearity the same way this search core does. See
    F5 in docs/STATE_OF_THE_PROJECT.md."""
    return _shared_longest_collinear_run(ranks)


def _block_anchor_ranks(block, home_rank):
    """Home ranks of the anchored (non-GOI, known) loci in a block, in the block's
    existing target order."""
    out = []
    for locus in block:
        r = home_rank.get(extract_base_gene_id(locus.get('query', '')))
        if r is not None:
            out.append(r)
    return out


def _can_bridge(block, candidate, home_rank, max_rank_gap, min_anchors):
    """True if `candidate` (the next locus after a large genomic gap) collinearly
    continues `block`'s home-order run — i.e. the gap is an insertion/rearrangement
    inside one conserved neighbourhood, not the boundary to an unrelated cluster.

    Conservative by design: requires the block to already carry `min_anchors`
    distinct anchored flanking genes (so a lone spurious hit can't seed a giant
    bridged block) and only bridges a candidate that extends BEYOND the furthest
    rank the block has reached, in the block's own direction, by at most
    `max_rank_gap` positions (a few un-mapped flanking genes tolerated).

    NB the direction is the BLOCK's own, taken from `_longest_collinear_run` over
    its anchors — not a fixed ascending order. A uniformly inverted neighbourhood
    therefore bridges normally; what is refused is a direction REVERSAL across the
    gap. Verified 2026-07-21 (see tests/test_synteny_collinearity.py):

        A B C … D E F  (forward, continues forward)   -> bridged
        F E D … C B A  (inverted, continues inverted) -> bridged
        D E F … C B A  (direction reverses at the gap)-> NOT bridged
        A B   … D E F  (2 anchors: bridged at the default min_anchors 2 since
                        2026-10-14, NOT bridged at 3, the default before)
        A     … B C D  (1 anchor: never bridged)

    With 2 anchors the direction rests on those two genes alone; the other limits
    (same scaffold, rank step <= max_rank_gap, at most synteny_bridge_max_per_block
    bridges) are unchanged.
    """
    cand_rank = home_rank.get(extract_base_gene_id(candidate.get('query', '')))
    if cand_rank is None:
        return False
    anchors = _block_anchor_ranks(block, home_rank)
    if len(set(anchors)) < min_anchors:
        return False
    _, direction = _longest_collinear_run(anchors)
    if direction == '+':
        frontier = max(anchors)
        return frontier < cand_rank <= frontier + max_rank_gap
    frontier = min(anchors)
    return frontier - max_rank_gap <= cand_rank < frontier


def _right_cluster(all_loci, i, cluster_distance):
    """The run of loci starting at index `i` that ordinary proximity clustering would
    join into one block (same chromosome, successive gaps < cluster_distance).

    This is the half of the evidence the one-sided rule never sees.
    """
    out = [all_loci[i]]
    j = i + 1
    while j < len(all_loci):
        prev, cur = out[-1], all_loci[j]
        if (cur['chrom'] == prev['chrom']
                and cur['start'] - prev['end'] < cluster_distance):
            out.append(cur)
            j += 1
        else:
            break
    return out


def _can_bridge_two_sided(block, right, home_rank, max_rank_gap, min_anchors):
    """True if the loci on BOTH sides of the gap together form one collinear run.

    Why this exists (measured 2026-09-09, docs/CAN_BRIDGE_ANALYSIS.md): `_can_bridge`
    counts anchors on the LEFT only and inspects exactly ONE locus on the right, because
    it is called from a left-to-right walk that has not read the rest yet. The result is
    that identical evidence decides differently depending on where the gap happens to
    fall --

        A B   ...gap... C D E   -> refused   (2 anchors left)
        A B C ...gap... D E     -> bridged   (3 anchors left)

    -- and that the first locus after the gap acts as a gatekeeper for everything behind
    it, so `A B C ...gap... D A A` bridges (D passes; the duplicates then join by ordinary
    proximity) while `A B C ...gap... A A D` does not.

    Here the whole right cluster is weighed with the block: the union must carry
    `min_anchors` distinct anchors AND be collinear over at least that many, and the rank
    step across the gap must continue the union's own direction by at most `max_rank_gap`.
    Inversions still bridge (the direction is derived, not assumed); a scrambled or
    duplicated continuation still does not.
    """
    left_ranks = [r for r in _block_anchor_ranks(block, home_rank) if r is not None]
    right_ranks = [home_rank.get(extract_base_gene_id(l.get('query', '')))
                   for l in right]
    right_ranks = [r for r in right_ranks if r is not None]
    if not left_ranks or not right_ranks:
        return False
    if len(set(left_ranks) | set(right_ranks)) < min_anchors:
        return False

    run, direction = _longest_collinear_run(left_ranks + right_ranks)
    if run < min_anchors:
        return False

    if direction == '+':
        frontier = max(left_ranks)
        return frontier < min(right_ranks) <= frontier + max_rank_gap
    frontier = min(left_ranks)
    return frontier - max_rank_gap <= max(right_ranks) < frontier


def build_home_rank(initial_db_path):
    """Map each flanking gene's base ID -> ordinal rank in home genomic order.

    PREPARE_INITIAL_DB writes the flanking genes first, in the order produced by
    extract_flanking_genes (sorted by chrom,start), then appends GOI proxy
    sequences. We rank flanking (non-GOI) base IDs by first appearance and skip
    GOI entries (they are not synteny anchors). Built from the pristine
    initial_db — NOT the expanding per-wave DB, whose order is not stable.
    Returns {} on any failure (collinearity then disabled, legacy behaviour).
    """
    home_rank = {}
    try:
        rank = 0
        for _header, clean_id, _seq in parse_fasta(initial_db_path):
            if is_goi_query_id(clean_id):
                continue
            base = extract_base_gene_id(clean_id)
            if base not in home_rank:
                home_rank[base] = rank
                rank += 1
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning(f"build_home_rank failed ({exc}); collinearity scoring disabled.")
        return {}
    return home_rank


def derive_goi_fallback_max_gap(goi_info_path, margin, floor, default_max_gap):
    """Data-driven max-intron for the GOI fallback hit-chain.

    The fallback chain (``_longest_monotonic_query_chain``) joins hits into one
    multi-exon gene; without a gap bound two far-apart spurious hits (e.g. two
    ~15-aa melittin matches 60 kb apart on unrelated genes) get stitched into one
    mRNA. Rather than a hand-set ``--max_intron``, scale the home GOI's OWN largest
    intron (from ANNOTATE_GOI's ``goi_info.json``) by a generous ``margin`` — the
    home annotation tells us this gene's natural intron sizes, so a big margin
    covers intron expansion in divergent targets while still rejecting the
    pathological joins. ``floor`` keeps a sane minimum for compact / single-exon
    genes (whose home ``max_intron_bp`` may be ~0). Falls back to
    ``default_max_gap`` (the static ``--max_intron``) when the home GOI structure
    is unavailable (no GFF / alignment-only with no coords).

    Returns ``(max_gap_bp, reason)``.
    """
    if not goi_info_path or not os.path.exists(goi_info_path):
        return default_max_gap, f"no goi_info; using --max_intron={default_max_gap}"
    try:
        with open(goi_info_path) as fh:
            info = json.load(fh)
    except Exception as exc:
        return default_max_gap, f"goi_info unreadable ({exc}); using --max_intron={default_max_gap}"
    home_max_intron = info.get("max_intron_bp")
    if home_max_intron is None:
        return default_max_gap, f"goi_info has no max_intron_bp; using --max_intron={default_max_gap}"
    try:
        derived = int(round(float(home_max_intron) * float(margin)))
    except (TypeError, ValueError):
        return default_max_gap, f"bad max_intron_bp={home_max_intron!r}; using --max_intron={default_max_gap}"
    max_gap = max(int(floor), derived)
    return max_gap, (
        f"home_max_intron={home_max_intron} bp x margin={margin} = {derived} bp, "
        f"floored at {floor} -> {max_gap} bp"
    )


def resolve_anchor_hits(hits, min_overlap=0.5):
    """One anchor per place: where hits of DIFFERENT flanking genes cover the same target
    DNA, only the best-scoring one stays.

    A target gene is the ortholog of at most one home gene. Two flanking genes that are
    relatives of each other (two kinases, two members of a tandem pair) both hit every
    member of their family, so one target gene used to count as two anchors: a place
    with a single gene passed the "at least two flanking genes" test for a block, and
    the number of flanking genes around a GOI call was inflated the same way. After this
    step a tandem pair A, B is counted as A on its own ortholog and B on its own.

    Two hits compete when they overlap by at least ``min_overlap`` of the shorter one,
    on either strand. The higher bit score wins; ties go to the lower E-value, then to
    the smaller query id. GOI hits are returned untouched.
    """
    import bisect

    keep = []
    by_chrom = defaultdict(list)
    for h in hits:
        if is_goi_query_id(h.get('query', '')):
            keep.append(h)
        else:
            by_chrom[h.get('chrom')].append(h)

    for chrom in sorted(by_chrom, key=str):
        # Strongest first; a hit stays unless a stronger hit of another gene already
        # holds the same place. (Greedy by strength, so a hit that lost cannot push out
        # a third one.)
        ranked = sorted(by_chrom[chrom],
                        key=lambda h: (-float(h.get('bits', 0) or 0), float(h.get('evalue', 0) or 0),
                                       h['query'], h['start'], h['end']))
        starts = []          # starts of the accepted hits, sorted
        accepted = []        # the accepted hits in the same order
        longest = 0
        for h in ranked:
            gene = extract_base_gene_id(h['query'])
            length = h['end'] - h['start']
            lost = False
            i = bisect.bisect_left(starts, h['end']) - 1
            while i >= 0 and starts[i] > h['start'] - longest:
                o = accepted[i]
                overlap = min(h['end'], o['end']) - max(h['start'], o['start'])
                shorter = min(length, o['end'] - o['start'])
                if (overlap > 0 and shorter > 0 and overlap >= min_overlap * shorter
                        and extract_base_gene_id(o['query']) != gene):
                    lost = True
                    break
                i -= 1
            if lost:
                continue
            pos = bisect.bisect_left(starts, h['start'])
            starts.insert(pos, h['start'])
            accepted.insert(pos, h)
            longest = max(longest, length)
        keep.extend(accepted)
    return keep


def limit_anchor_places(hits, max_places, gap_limit=None, default_gap=20000):
    """Flanking hits restricted to the ``max_places`` best places of each flanking gene.

    A place is a run of hits of one flanking gene on one target sequence whose neighbours
    are at most that gene's gap limit apart (``gap_limit(gene)``, the limit that separates
    two occurrences of the gene); its weight is the sum of its bit scores. The ortholog
    of a flanking gene is one place; its further places are relatives and shared domains.
    They are what makes two anchors lie near each other by chance: on the saved hits of 19
    genome searches (yeast, fly, bee) the blocks with two or more anchors fell from 258 to
    53 with five places per gene, and every block that held a HIGH, MEDIUM or ambiguous
    call was still there. Off by default (``max_places`` <= 0 returns the hits unchanged):
    not yet run on the family benchmark.

    Only the seeding of blocks uses the result. The other hits still count inside a
    window (flanking support, flanking gene models).
    """
    if not max_places or int(max_places) <= 0:
        return list(hits)
    by_gene_seq = defaultdict(list)
    for h in hits:
        by_gene_seq[(extract_base_gene_id(h.get('query', '')), h.get('chrom'))].append(h)
    places = defaultdict(list)
    for (gene, _seq), gene_hits in by_gene_seq.items():
        limit = gap_limit(gene) if gap_limit else default_gap
        gene_hits = sorted(gene_hits, key=lambda h: (min(h['start'], h['end']), max(h['start'], h['end'])))
        current = [gene_hits[0]]
        current_end = max(gene_hits[0]['start'], gene_hits[0]['end'])
        for h in gene_hits[1:]:
            start, end = min(h['start'], h['end']), max(h['start'], h['end'])
            if start - current_end <= limit:
                current.append(h)
                current_end = max(current_end, end)
            else:
                places[gene].append(current)
                current, current_end = [h], end
        places[gene].append(current)
    keep = set()
    for gene, gene_places in places.items():
        ranked = sorted(gene_places, key=lambda pl: (
            -sum(float(x.get('bits', 0) or 0) for x in pl),
            str(pl[0].get('chrom')), min(pl[0]['start'], pl[0]['end'])))
        for place in ranked[:int(max_places)]:
            keep.update(id(x) for x in place)
    return [h for h in hits if id(h) in keep]


def identify_synteny_blocks(hits, max_intron=20000, cluster_distance=50000,
                            home_rank=None, bridge_max_gap=0,
                            bridge_max_rank_gap=5, bridge_min_anchors=2,
                            bridge_two_sided=False, bridge_max_per_block=2,
                            gene_gap_limit=None):
    """
    Identify all synteny blocks from hits.
    
    Key Logic:
    - Group hits by query gene (handling multi-exon genes AND exon-level queries)
    - Cluster loci that are close together (likely same gene/region)
    - Score blocks by number of unique flanking genes found
    
    Updated to handle exon-level queries like "gene-LOC726866|exon_1"

    Collinearity (§1m): when `home_rank` is supplied, two refinements activate
    (otherwise behaviour is identical to the legacy gene-count version):
      - gap-bridging (#1): a large genomic gap between two same-chromosome
        clusters is NOT split into separate blocks when the flanking continue the
        same home-order run across it (see _can_bridge), so a GOI sitting in such
        a gap is inside the search window.
      - collinear scoring (#3): blocks are ranked by the length of their longest
        collinear run of flanking genes, then by gene count, so a genuine ordered
        block outranks a scrambled paralog-neighbourhood cluster.

    Args:
        hits: List of hit dictionaries
        max_intron: Maximum distance between exons of same gene (bp)
        cluster_distance: Maximum distance to cluster genes into synteny block (bp)
        home_rank: {flanking_base_id: ordinal home-order rank} or None/{} (legacy)
        bridge_max_gap: max bp gap to bridge when collinear (0 disables bridging)
        bridge_max_rank_gap: max home-rank jump still treated as a continuation
        bridge_min_anchors: min distinct anchored flanking genes before bridging
        bridge_two_sided: weigh the whole cluster on the far side of a gap instead of
            only its first locus (docs/CAN_BRIDGE_ANALYSIS.md). False = legacy rule.
        bridge_max_per_block: max bridges one block may chain (0 = unlimited). Bounds a
            runaway block: each bridge can span bridge_max_gap, so without a cap a chain
            of them could swallow a chromosome.
        gene_gap_limit: optional ``f(gene_id) -> bp``, the largest gap between two hits of
            that gene that are still one occurrence of it (its own span at home, scaled).
            Without it ``max_intron`` applies to every gene.

    Returns:
        List of dictionaries (blocks), sorted by score (descending). Each block
        carries 'collinear_chain_len'/'collinear_direction'/'bridged'.
    """
    use_collinearity = bool(home_rank)
    if not hits:
        return []
        
    # --- Step 1: Group hits by Query Gene (Base ID) ---
    # Multiple hits from same query = different exons or duplicates
    # With exon_mode, query may be "gene-LOC726866|exon_1" -> base is "gene-LOC726866"
    hits_by_query = defaultdict(list)
    for h in hits:
        # Extract base query ID (handles both |var and |exon_ suffixes)
        base_query = extract_base_gene_id(h['query'])
        hits_by_query[base_query].append(h)
        
    # --- Step 2: Define Gene Loci per Query ---
    # Each query gene may have multiple loci (paralogs/duplications)
    # Hits close together (<max_intron) = same locus (multi-exon gene)
    all_loci = []
    for query_id, q_hits in hits_by_query.items():
        # Sort hits by genomic position
        q_hits.sort(key=lambda x: (x['chrom'], x['start']))
        gene_limit = gene_gap_limit(query_id) if gene_gap_limit else max_intron

        current_locus_hits = []
        for h in q_hits:
            if not current_locus_hits:
                current_locus_hits.append(h)
                continue
            
            last_hit = current_locus_hits[-1]
            
            # Same chromosome and close enough = same locus (exons)
            if (h['chrom'] == last_hit['chrom'] and 
                h['start'] - last_hit['end'] < gene_limit):
                current_locus_hits.append(h)
            else:
                # Start new locus
                all_loci.append(create_locus_object(query_id, current_locus_hits))
                current_locus_hits = [h]
        
        # Don't forget last locus
        if current_locus_hits:
            all_loci.append(create_locus_object(query_id, current_locus_hits))

    # --- Step 3: Cluster Loci into Synteny Blocks ---
    # Loci from different genes that are close = synteny block
    all_loci.sort(key=lambda x: (x['chrom'], x['start']))
    if not all_loci: 
        return []
        
    bridging_on = use_collinearity and bridge_max_gap > 0
    synteny_blocks = []
    current_block = [all_loci[0]]
    bridges_in_block = 0

    i = 1
    while i < len(all_loci):
        locus = all_loci[i]
        last_locus = current_block[-1]
        same_chrom = locus['chrom'] == last_locus['chrom']
        gap = locus['start'] - last_locus['end'] if same_chrom else None

        if same_chrom and gap < cluster_distance:
            # Normal clustering: genes close together = same neighbourhood.
            current_block.append(locus)
            i += 1
            continue

        bridge_allowed = (bridging_on and same_chrom and gap < bridge_max_gap
                          and (bridge_max_per_block <= 0
                               or bridges_in_block < bridge_max_per_block))

        if bridge_allowed and bridge_two_sided:
            # #1b two-sided gap-bridging: weigh the whole cluster on the far side of the
            # gap, not just its first locus, so the decision does not depend on where the
            # rearrangement happened to fall. See _can_bridge_two_sided.
            right = _right_cluster(all_loci, i, cluster_distance)
            if _can_bridge_two_sided(current_block, right, home_rank,
                                     bridge_max_rank_gap, bridge_min_anchors):
                current_block.extend(right)
                bridges_in_block += 1
                i += len(right)
                continue
        elif bridge_allowed and _can_bridge(current_block, locus, home_rank,
                                            bridge_max_rank_gap, bridge_min_anchors):
            # #1 gap-bridging: the flanking continue the home-order run across a
            # large gap (a local rearrangement that the GOI may sit inside), so
            # keep one block spanning the gap instead of dropping the interior.
            current_block.append(locus)
            bridges_in_block += 1
            i += 1
            continue

        synteny_blocks.append(current_block)
        current_block = [locus]
        bridges_in_block = 0
        i += 1
    synteny_blocks.append(current_block)

    # --- Step 4: Format, Score (collinearity) and Sort Blocks ---
    final_blocks = []
    for block in synteny_blocks:
        chrom = block[0]['chrom']
        start = min(l['start'] for l in block)
        end = max(l['end'] for l in block)
        genes_list = sorted(set(extract_base_gene_id(l['query']) for l in block))

        if use_collinearity:
            # One representative locus per flanking GENE, in target order — the same
            # basis process_region_block uses for the classifier's collinear_support.
            # Before this, every locus contributed a rank, and _longest_collinear_run's
            # LIS is NON-decreasing, so repeated ranks counted as a run: a block that is
            # ONE gene hit four times reported collinear_chain_len=4 (observed on real
            # data, oskar NC_045757.1 / gene-Dmel_CG31100) and outranked a block of three
            # genuinely ordered genes in the sort below and in the per-genome cap.
            first_pos = {}
            for l in block:
                gid = extract_base_gene_id(l.get('query', ''))
                pos = int(l.get('start', 0) or 0)
                if gid not in first_pos or pos < first_pos[gid][0]:
                    first_pos[gid] = (pos, home_rank.get(gid))
            ranks = [rank for _gid, (_pos, rank)
                     in sorted(first_pos.items(), key=lambda kv: kv[1][0])]
            chain_len, direction = _longest_collinear_run(ranks)
            # A block "bridged" if any adjacent loci span >= cluster_distance.
            bridged = any(
                block[i + 1]['chrom'] == block[i]['chrom'] and
                (block[i + 1]['start'] - block[i]['end']) >= cluster_distance
                for i in range(len(block) - 1)
            )
        else:
            chain_len, direction, bridged = 0, '+', False

        final_blocks.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'genes_count': len(genes_list),
            'loci_count': len(block),
            'genes': genes_list,
            'collinear_chain_len': chain_len,
            'collinear_direction': direction,
            'bridged': bridged,
        })

    if use_collinearity:
        # #3: rank by collinear-run length first, gene count as tiebreak.
        final_blocks.sort(key=lambda x: (x['collinear_chain_len'], x['genes_count']),
                          reverse=True)
    else:
        # Legacy: sort by raw gene count (unchanged for callers without home_rank).
        final_blocks.sort(key=lambda x: x['genes_count'], reverse=True)

    return final_blocks


def identify_best_synteny_block(hits, max_intron=20000, cluster_distance=50000):
    """
    Backwards-compatible wrapper kept for tests/importers.
    """
    blocks = identify_synteny_blocks(hits, max_intron=max_intron, cluster_distance=cluster_distance)
    return blocks[0] if blocks else None

def calculate_adaptive_padding(hits: List[Dict[str, Any]], best_region: Dict[str, Any],
                               default: int = 100000, min_pad: int = 50000, max_pad: int = 200000) -> int:
    """
    Calculate region padding based on gene spacing in hits.
    Returns padding distance in base pairs.
    """
    # Filter hits to the best region's chromosome
    region_hits = [h for h in hits if h['chrom'] == best_region['chrom']]
    
    if len(region_hits) < 2:
        return default
    
    # Sort by position
    sorted_hits = sorted(region_hits, key=lambda h: h['start'])
    
    # Calculate inter-gene gaps
    gaps = []
    for i in range(len(sorted_hits) - 1):
        gap = sorted_hits[i+1]['start'] - sorted_hits[i]['end']
        if gap > 0:  # Only positive gaps
            gaps.append(gap)
    
    if not gaps:
        return default
    
    # Average gap * 2 (to cover one gene on each side)
    avg_gap = sum(gaps) / len(gaps)
    adaptive_padding = int(avg_gap * 2)
    
    # Clamp to reasonable range
    final_padding = max(min_pad, min(max_pad, adaptive_padding))
    
    return final_padding


def bp_or_auto(value):
    """argparse type of a distance that is either derived ('auto') or fixed in bp."""
    text = str(value).strip().lower()
    if text in ('auto', ''):
        return None
    return int(float(text))


def _expansion_note(source, run_anchors):
    """How the expansion of one genome came about, for the log line."""
    if source == 'anchors':
        return f"measured over {run_anchors} anchors in a row"
    longest = f"longest run of anchors in home order: {run_anchors}"
    if source == 'genome size':
        return f"from the genome sizes; {longest}"
    return f"not measured; {longest}"


def genome_scale_args(args, flanking_hits, genome_path, genome_name=''):
    """A copy of ``args`` whose distances are in base pairs for THIS target genome.

    Distances that describe a neighbourhood are given in gene pitches and converted with
    the pitch of the home neighbourhood and the expansion measured in this genome
    (synvoy_scale). A distance the user set in bp is used as it is. Without home
    coordinates in the query set, every distance keeps its former fixed value.

    The copy also carries ``distances`` (what was used and where each value came from,
    written into the result GFF) and ``expansion``.
    """
    import copy as _copy
    geometry = getattr(args, 'home_geometry', None)
    size_ratio = None
    home_bp = (geometry or {}).get('home_genome_bp')
    if home_bp:
        try:
            target_bp = os.path.getsize(genome_path) * (4 if str(genome_path).endswith(('.gz', '.gzip')) else 1)
            size_ratio = target_bp / float(home_bp)
        except OSError:
            size_ratio = None
    expansion, run_anchors, source = synvoy_scale.target_expansion(geometry, flanking_hits, size_ratio=size_ratio)
    distances = synvoy_scale.derive_distances(
        geometry, expansion,
        genes={
            'block_gap_genes': getattr(args, 'block_gap_genes', None),
            'window_margin_genes': getattr(args, 'window_margin_genes', None),
            'goi_margin_genes': getattr(args, 'goi_margin_genes', None),
            'bridge_max_gap_genes': getattr(args, 'synteny_bridge_max_gap_genes', None),
            'envelope_margin_genes': getattr(args, 'dispersed_envelope_margin_genes', None),
        },
        fixed_bp={
            'cluster_distance': getattr(args, 'cluster_distance', None),
            'region_padding': getattr(args, 'region_padding', None),
            'gap_search_window': getattr(args, 'goi_margin', None),
            'synteny_bridge_max_gap': getattr(args, 'synteny_bridge_max_gap', None),
            'dispersed_envelope_margin': getattr(args, 'dispersed_envelope_margin', None),
        },
    )
    distances['expansion_source'] = source
    distances['expansion_anchors'] = run_anchors
    out = _copy.copy(args)
    out.cluster_distance = distances['cluster_distance']
    # A window margin given in bp keeps the old behaviour (adapted to the spacing of the
    # hits, between padding_min and padding_max); a derived one is used as it is.
    out.adaptive_padding = getattr(args, 'region_padding', None) is not None
    out.region_padding = distances['region_padding']
    out.goi_margin = distances['goi_margin']
    # DNA cut out around a block. Only hits that touch the core (block +- GOI margin) are
    # used, and a gene model may use the DNA up to the model margin beyond its hits, so
    # anything further out would be searched for nothing. The window margin stays the
    # extent of the neighbourhood (which blocks merge, which anchors count as support);
    # the DNA searched is the smaller of the two. With a window margin given in bp the
    # whole margin is searched, as before. (At a human locus the difference is megabases:
    # 7 pitches of 335 kb on each side against 2 pitches plus the room for one gene.)
    out.search_margin = None
    if not out.adaptive_padding and geometry:
        model_margin = max(int(getattr(args, 'gap_search_window', 0) or 0),
                           int(float(getattr(args, 'goi_home_span', 0) or 0) * float(expansion or 1.0)))
        out.search_margin = int(min(distances['region_padding'], distances['goi_margin'] + model_margin))
        distances['search_margin'] = out.search_margin
    # 0 switches bridging off, in either unit.
    bridge_off = (getattr(args, 'synteny_bridge_max_gap', None) == 0
                  or (getattr(args, 'synteny_bridge_max_gap', None) is None
                      and float(getattr(args, 'synteny_bridge_max_gap_genes', 1) or 0) <= 0))
    out.synteny_bridge_max_gap = 0 if bridge_off else distances['synteny_bridge_max_gap']
    distances['synteny_bridge_max_gap'] = out.synteny_bridge_max_gap
    out.dispersed_envelope_margin = distances['dispersed_envelope_margin']
    out.expansion = expansion
    out.distances = distances
    if geometry:
        logger.info(
            f"[{genome_name}] distances: home pitch {geometry['pitch'] / 1000:.1f} kb x expansion "
            f"{expansion:.2f} ({_expansion_note(source, run_anchors)}) -> block gap "
            f"{distances['cluster_distance'] / 1000:.0f} kb, window margin "
            f"{distances['region_padding'] / 1000:.0f} kb"
            + (f" (DNA searched: {out.search_margin / 1000:.0f} kb per side)" if out.search_margin else "")
            + f", GOI margin {distances['goi_margin'] / 1000:.0f} kb, "
            f"bridge up to {out.synteny_bridge_max_gap / 1000:.0f} kb.")
    else:
        logger.info(f"[{genome_name}] distances: no home coordinates in the query set; fixed values "
                    f"(block gap {distances['cluster_distance'] / 1000:.0f} kb).")
    return out


def _gene_limit(args, gene_id, floor_bp):
    """Largest gap inside one occurrence of a flanking gene in this target (bp)."""
    return synvoy_scale.gene_gap_limit(getattr(args, 'home_geometry', None),
                                       extract_base_gene_id(gene_id),
                                       float(getattr(args, 'expansion', 1.0) or 1.0), floor_bp)


def run_augmented_mmseqs_with_retries(
    variants_fasta: str,
    region_fasta: str,
    aug_hits_file: str,
    aug_tmp_dir: str,
    args,
    threads: int,
    genome_name: str,
) -> Tuple[bool, str]:
    """
    Run augmented MMseqs search with low-memory retries.

    Unlike the primary whole-genome search, this runs on per-block regions, but
    still needs retry logic because strict split-memory limits (e.g. 1G) can
    fail on translated ORF DB construction for certain blocks.
    """
    base_attempts = [
        {
            "label": "aug_primary",
            "threads": max(1, int(threads)),
            "sens": float(args.mmseqs_sens),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "aug_lowmem_auto",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 7.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "aug_lowmem_2g",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 7.0),
            "split": "2G",
        },
        {
            "label": "aug_lowmem_verysafe",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 6.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
    ]

    attempts = []
    seen = set()
    for a in base_attempts:
        key = (a["threads"], a["sens"], a["split"])
        if key in seen:
            continue
        seen.add(key)
        attempts.append(a)

    os.makedirs(aug_tmp_dir, exist_ok=True)
    last_details = ""

    # Create query and target databases ONCE.
    query_db = os.path.join(aug_tmp_dir, "queryDB")
    target_db = os.path.join(aug_tmp_dir, "targetDB")
    fmt_output = "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    proc_qdb = subprocess.run(
        ["mmseqs", "createdb", variants_fasta, query_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_qdb.returncode != 0:
        return False, _format_mmseqs_failure(proc_qdb)

    proc_tdb = subprocess.run(
        ["mmseqs", "createdb", region_fasta, target_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_tdb.returncode != 0:
        return False, _format_mmseqs_failure(proc_tdb)

    for i, attempt in enumerate(attempts, start=1):
        if os.path.exists(aug_hits_file):
            try:
                os.remove(aug_hits_file)
            except OSError:
                pass

        attempt_tmp = os.path.join(aug_tmp_dir, f"attempt_{i}")
        os.makedirs(attempt_tmp, exist_ok=True)
        result_db = os.path.join(attempt_tmp, "resultDB")

        search_cmd = [
            "mmseqs", "search",
            query_db, target_db, result_db, attempt_tmp,
            "--search-type", "2",
            "--threads", str(attempt["threads"]),
            "-s", str(attempt["sens"]),
            "-e", str(min(args.aug_relaxed_evalue_cap, args.evalue * args.aug_relaxed_evalue_mult)),
            "--min-seq-id", "0.0",
            "--split-memory-limit", str(attempt["split"]),
            "-v", str(args.mmseqs_verbosity),
        ]

        proc = subprocess.run(
            search_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

        if proc.returncode != 0:
            details = _format_mmseqs_failure(proc)
            last_details = details
            resource_fail = _is_mmseqs_resource_failure(details)
            logger.warning(
                f"[{genome_name}] Augmented MMseqs search attempt {i}/{len(attempts)} failed "
                f"({attempt['label']}): {details}"
            )
            if not resource_fail:
                break
            continue

        # convertalis
        conv_cmd = [
            "mmseqs", "convertalis",
            query_db, target_db, result_db, aug_hits_file,
            "--format-output", fmt_output,
            "-v", str(args.mmseqs_verbosity),
        ]
        conv_proc = subprocess.run(
            conv_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        if conv_proc.returncode != 0:
            last_details = _format_mmseqs_failure(conv_proc)
            break

        if i > 1:
            logger.warning(
                f"[{genome_name}] Augmented MMseqs succeeded on retry {i}/{len(attempts)} "
                f"({attempt['label']}: threads={attempt['threads']}, "
                f"s={attempt['sens']}, split={attempt['split']})."
            )
        return True, ""

    return False, last_details or "unknown augmented MMseqs failure"

# E-value handed to tblastn in a search window, whatever the window search's own settings.
# Not only an output filter: tblastn links HSPs by sum statistics only among those it keeps,
# so its cutoff also changes the E-value of what it reports. With the old cutoff of 10
# (min(aug_relaxed_evalue_cap, search_evalue x aug_relaxed_evalue_mult)) the curated
# ortholog of one hidden pair of the Candida pilot (2026-10-14) was absent from tblastn's
# output in SynVoy's own 43-kb window and reported at E 0.005 with 1000. The kept threshold
# (E <= 0.1, identity >= 15 %, length >= 15) is applied when the output is read, as before.
TBLASTN_REPORT_EVALUE = 1000.0


def run_augmented_search(region_fasta: str, goi_queries: List[Dict[str, str]], 
                        genome_name: str, args, unique_id: str, threads: int) -> List[Dict[str, Any]]:
    """
    Run augmented search (MMseqs2 + tblastn + Smith-Waterman) for GOI queries.
    
    Uses both methods for maximum sensitivity:
    1. MMseqs2 with query fragments (fast, good for similar sequences)
    2. tblastn over the region (sensitive translated search)
    3. Smith-Waterman via parasail/ssearch36 (slower, better for divergent sequences)
    
    Args:
        region_fasta: Path to extracted region FASTA
        goi_queries: List of GOI query dicts with 'id' and 'seq'
        genome_name: Name of current genome
        args: Command line arguments
        unique_id: Unique ID for temp files
        threads: Number of threads to use
        
    Returns:
        List of hit dictionaries combining MMseqs2 and Smith-Waterman results
    """
    all_hits = []
    
    try:
        # Generate variants for each GOI query
        variants_fasta = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_variants.faa"
        query_fasta = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_goi_query.faa"
        all_variants = []
        
        # Write full query sequences for Smith-Waterman
        write_fasta([(q['id'], q['seq']) for q in goi_queries], query_fasta)
        goi_query_lengths = {
            q['id']: len(q['seq'])
            for q in goi_queries
            if q.get('id') and q.get('seq')
        }
        
        if not FRAGMENT_SUPPORT:
            print(f"[{genome_name}] Warning: fragment_query module not available, using full sequences only", flush=True)
            # Just use the original sequences
            all_variants = [(q['id'], q['seq']) for q in goi_queries]
        else:
            for query in goi_queries:
                # Generate fragments (halves, thirds, quarters)
                fragments = generate_fragments(query['seq'], query['id'], min_size=15)
                all_variants.extend([(f[0], f[1]) for f in fragments])
        
        write_fasta(all_variants, variants_fasta)
        variant_query_lengths = {
            vid: len(vseq)
            for vid, vseq in all_variants
            if vid and vseq
        }
        
        # Shared relaxed thresholds for augmented search.
        relaxed_identity = max(args.aug_relaxed_identity_min, args.min_identity * args.aug_relaxed_identity_factor)
        relaxed_length = max(args.aug_relaxed_length_min, int(args.min_length // args.aug_relaxed_length_div))

        # ========== 1. MMseqs2 Search ==========
        aug_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_aug_hits.m8"
        aug_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_aug_mmseqs"
        
        os.makedirs(aug_tmp_dir, exist_ok=True)
        
        # A0: the augmented GOI search runs over a tiny per-block region, so
        # forcing it single-threaded costs little but removes mmseqs' multithread
        # marginal-hit jitter — the root of the GOI confidence flips that cascade
        # through the wave seed set (TODO_JUN §A0, pinned source #2).
        aug_threads = 1 if getattr(args, 'deterministic_goi_search', True) else threads
        mmseqs_ok, mmseqs_details = run_augmented_mmseqs_with_retries(
            variants_fasta=variants_fasta,
            region_fasta=region_fasta,
            aug_hits_file=aug_hits_file,
            aug_tmp_dir=aug_tmp_dir,
            args=args,
            threads=aug_threads,
            genome_name=genome_name,
        )

        if mmseqs_ok:
            # Parse hits - use relaxed thresholds for augmented search
            relaxed_identity = max(args.aug_relaxed_identity_min, args.min_identity * args.aug_relaxed_identity_factor)
            relaxed_length = max(args.aug_relaxed_length_min, int(args.min_length // args.aug_relaxed_length_div))
            
            mmseqs_hits = parse_hits(
                aug_hits_file,
                relaxed_identity,
                relaxed_length,
                args.evalue * args.aug_relaxed_parse_evalue_mult,
                query_lengths=variant_query_lengths
            )
            if mmseqs_hits:
                print(f"[{genome_name}] MMseqs2 augmented search found {len(mmseqs_hits)} hits.", flush=True)
                all_hits.extend(mmseqs_hits)
        else:
            print(f"[{genome_name}] Augmented MMseqs failed after retries: {mmseqs_details}", flush=True)

        # ========== 2. tblastn Search ==========
        tblastn_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_tblastn_hits.m8"
        blast_db_prefix = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_tblastn_db"
        try:
            makedb_result = subprocess.run(
                ["makeblastdb", "-in", region_fasta, "-dbtype", "nucl", "-out", blast_db_prefix],
                check=True,
                capture_output=True, text=True
            )
            tblastn_result = subprocess.run(
                [
                    "tblastn",
                    "-query", query_fasta,
                    "-db", blast_db_prefix,
                    "-evalue", str(TBLASTN_REPORT_EVALUE),
                    "-seg", "no",
                    "-max_target_seqs", "5000",
                    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                    "-out", tblastn_hits_file
                ],
                check=True,
                capture_output=True, text=True
            )

            tblastn_hits = parse_hits(
                tblastn_hits_file,
                relaxed_identity,
                relaxed_length,
                args.evalue * args.aug_relaxed_parse_evalue_mult,
                query_lengths=goi_query_lengths
            )
            if tblastn_hits:
                print(f"[{genome_name}] tblastn augmented search found {len(tblastn_hits)} hits.", flush=True)
                for hit in tblastn_hits:
                    hit['method'] = 'tblastn'
                all_hits.extend(tblastn_hits)
        except FileNotFoundError:
            logger.warning(f"[{genome_name}] tblastn/makeblastdb not found on PATH. Install BLAST+ to enable tblastn search.")
        except subprocess.CalledProcessError as blast_err:
            stderr_snippet = (blast_err.stderr or '')[:300].strip()
            logger.warning(
                f"[{genome_name}] tblastn failed (exit {blast_err.returncode}), continuing with MMseqs/SW only."
                + (f" stderr: {stderr_snippet}" if stderr_snippet else "")
            )
        
        # ========== 3. Smith-Waterman Search ==========
        # Use Smith-Waterman for very divergent sequences (more sensitive than MMseqs2)
        if args.enable_smith_waterman:
            sw_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_sw_hits.m8"
            
            try:
                # Use smith_waterman_search.py script
                sw_cmd = [
                    "python3", os.path.join(os.path.dirname(__file__), "smith_waterman_search.py"),
                    "--query", query_fasta,
                    "--target", region_fasta,
                    "--output", sw_hits_file,
                    "--min_score", str(int(round(args.sw_min_score))),
                    "--min_identity", str(args.sw_min_identity),
                    "--threads", str(threads),
                    "--method", str(args.sw_method)
                ]
                
                result = subprocess.run(sw_cmd, capture_output=True, text=True, timeout=args.sw_timeout_seconds)
                
                if result.returncode == 0 and os.path.exists(sw_hits_file):
                    # Parse Smith-Waterman hits (BLAST m8 format)
                    sw_hits = parse_hits(
                        sw_hits_file,
                        args.sw_min_identity,
                        10,    # min 10 aa (was 2 – way too permissive)
                        # E <= 1 over the window's six-frame search space. Only a real
                        # filter since 2026-09-11: SW used to write E=0.001 for every hit.
                        # On shuffled queries 2 % of per-frame best hits pass (~0.1 per
                        # window); the weak melittin-family members SW exists for sit at
                        # E 0.002-0.25 (smith_waterman_search.KARLIN_ALTSCHUL_BLOSUM62).
                        1.0,
                        query_lengths=goi_query_lengths
                    )
                    if sw_hits:
                        print(f"[{genome_name}] Smith-Waterman found {len(sw_hits)} additional hits.", flush=True)
                        # Mark hits as from Smith-Waterman
                        for hit in sw_hits:
                            hit['method'] = 'smith_waterman'
                        all_hits.extend(sw_hits)
                elif result.stderr:
                    print(f"[{genome_name}] Smith-Waterman warning: {result.stderr[:200]}", flush=True)
                    
            except subprocess.TimeoutExpired:
                print(f"[{genome_name}] Smith-Waterman timed out after {args.sw_timeout_seconds} seconds, using MMseqs2 only.", flush=True)
            except FileNotFoundError:
                print(f"[{genome_name}] Smith-Waterman script not found, using MMseqs2 only.", flush=True)
            except Exception as sw_err:
                print(f"[{genome_name}] Smith-Waterman failed: {sw_err}, using MMseqs2 only.", flush=True)
        
        # ========== 4. PLM Embedding Search ==========
        # Use ProtT5 protein language model to find ORFs that are functionally
        # similar to the GOI even when sequence identity is undetectable (<15%).
        if (GPU_IN_WORKER and getattr(args, 'enable_plm_search', False)
                and PLM_IMPORT_OK and check_plm_available()):
            try:
                goi_emb_path = os.path.join(args.output_dir, "goi_embeddings.npz")
                if os.path.exists(goi_emb_path):
                    goi_embeddings = load_embeddings(goi_emb_path)
                    if goi_embeddings:
                        plm_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_plm"
                        os.makedirs(plm_tmp_dir, exist_ok=True)

                        # Short queries produce noisy mean-pooled embeddings —
                        # raise the similarity threshold to suppress false hits.
                        max_goi_len = max(goi_query_lengths.values()) if goi_query_lengths else 0
                        plm_thresh = getattr(args, 'plm_similarity_threshold', 0.5)
                        if 0 < max_goi_len < 100:
                            plm_thresh = max(plm_thresh, 0.75)

                        plm_hits = plm_search_region(
                            goi_embeddings=goi_embeddings,
                            region_fasta=region_fasta,
                            output_dir=plm_tmp_dir,
                            similarity_threshold=plm_thresh,
                            device=getattr(args, 'plm_device', 'cpu'),
                            predictor='augustus',  # Enforce Augustus for eukaryotic ORF prediction; no Prodigal fallback (loses multi-exon genes)
                            augustus_species=getattr(args, 'augustus_species', 'fly'),
                        )

                        if plm_hits:
                            print(
                                f"[{genome_name}] PLM embedding search found "
                                f"{len(plm_hits)} candidate ORF(s).",
                                flush=True,
                            )
                            all_hits.extend(plm_hits)

                        if os.path.exists(plm_tmp_dir):
                            shutil.rmtree(plm_tmp_dir, ignore_errors=True)
            except Exception as plm_err:
                logger.debug(f"[{genome_name}] PLM search failed: {plm_err}")

        # ========== 5. Foldseek Structural Search ==========
        # Use ESMFold to predict ORF structures and Foldseek 3Di alphabet to
        # find structural homologs that have diverged beyond all sequence methods.
        # Skip structural search for very short queries — ESMFold produces
        # degenerate structures below ~50 aa and TM-score is meaningless for
        # single short helices (e.g. melittin, defensins).
        _max_goi_len = max(goi_query_lengths.values()) if goi_query_lengths else 0
        _struct_too_short = 0 < _max_goi_len < 50
        if (GPU_IN_WORKER and getattr(args, 'enable_structural_search', False)
                and STRUCTURAL_IMPORT_OK and check_structural_search_available()
                and _struct_too_short):
            print(
                f"[{genome_name}] Skipping structural search: "
                f"longest GOI query is {_max_goi_len} aa (<50 aa threshold).",
                flush=True,
            )

        if (GPU_IN_WORKER and getattr(args, 'enable_structural_search', False)
                and STRUCTURAL_IMPORT_OK and check_structural_search_available()
                and not _struct_too_short):
            try:
                goi_struct_index = os.path.join(args.output_dir, "goi_structure_index.tsv")
                if os.path.exists(goi_struct_index):
                    goi_structures = load_structure_index(goi_struct_index)
                    if goi_structures:
                        struct_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_struct"
                        os.makedirs(struct_tmp_dir, exist_ok=True)

                        struct_hits = structural_search_region(
                            goi_structures=goi_structures,
                            region_fasta=region_fasta,
                            output_dir=struct_tmp_dir,
                            tm_threshold=getattr(
                                args, 'structural_tm_threshold', 0.3
                            ),
                            device=getattr(args, 'structural_device', 'cpu'),
                            max_length=getattr(args, 'structural_max_length', 700),
                            threads=getattr(args, 'threads', 1),
                            predictor='augustus',  # Enforce Augustus for eukaryotic ORF prediction; no Prodigal fallback (loses multi-exon genes)
                            augustus_species=getattr(args, 'augustus_species', 'fly'),
                        )

                        if struct_hits:
                            print(
                                f"[{genome_name}] Foldseek structural search found "
                                f"{len(struct_hits)} candidate ORF(s).",
                                flush=True,
                            )
                            all_hits.extend(struct_hits)

                        if os.path.exists(struct_tmp_dir):
                            shutil.rmtree(struct_tmp_dir, ignore_errors=True)
            except Exception as struct_err:
                logger.debug(f"[{genome_name}] Structural search failed: {struct_err}")

        # Clean up temp files
        files_to_remove = [variants_fasta, query_fasta, aug_hits_file]
        if 'sw_hits_file' in locals():
             files_to_remove.append(sw_hits_file)
        files_to_remove.append(tblastn_hits_file)
            
        for f in files_to_remove:
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(aug_tmp_dir):
            shutil.rmtree(aug_tmp_dir, ignore_errors=True)
        for f in glob.glob(f"{blast_db_prefix}*"):
            if os.path.exists(f):
                os.remove(f)
        
        # Deduplicate hits by position (Greedy Spatial Clustering)
        if all_hits:
            # Sort by bitscore desc, then evalue asc
            # Ensure bits is float
            for h in all_hits:
                if 'bits' not in h: h['bits'] = 0.0
                h['bits'] = float(h['bits'])
            
            all_hits.sort(key=lambda x: (-x['bits'], x['evalue'],
                                         x['chrom'], x['start'], x['end'], x['query']))
            
            kept_hits = []
            for hit in all_hits:
                # Check overlap with kept hits
                is_overlapped = False
                h_start = hit['start']
                h_end = hit['end']
                h_len = h_end - h_start
                
                for kept in kept_hits:
                    if hit['chrom'] != kept['chrom']: continue
                    
                    k_start = kept['start']
                    k_end = kept['end']
                    
                    # Calculate overlap
                    o_start = max(h_start, k_start)
                    o_end = min(h_end, k_end)
                    overlap = max(0, o_end - o_start)
                    
                    # If overlap covers > 50% of THIS hit, drop it (we already have a better one)
                    if overlap > 0.5 * h_len:
                        is_overlapped = True
                        break
                
                if not is_overlapped:
                    kept_hits.append(hit)
            
            all_hits = kept_hits
            print(f"[{genome_name}] Combined augmented search: {len(all_hits)} unique hits (spatially filtered).", flush=True)
        
        return all_hits
        
    except Exception as e:
        import traceback
        logger.warning(
            f"[{genome_name}] Augmented search failed: {e}\n"
            f"  Traceback: {traceback.format_exc().strip().splitlines()[-3:]}"
        )
        return []

def batch_rbh_check(
    candidates,
    home_db,
    cand_map,
    threads=1,
    evalue=1e-5,
    min_coverage=0.5,
    min_identity=25.0,
    split_memory_limit="8G",
    mmseqs_verbosity=1,
):
    """
    Perform Reciprocal Best Hit check for multiple candidates at once.
    
    Enhanced validation:
    1. RBH to home genome (traditional)
    2. Coverage check: alignment must cover >50% of both query and target
    3. Identity must be reasonable (>25%)
    
    candidates: list of dicts with 'id' and 'seq' keys (or BioPython SeqRecord)
    """
    if not candidates:
        return set()
    
    unique_id = uuid.uuid4().hex
    query_fasta = f"batch_candidates_{unique_id}.fasta"
    rbh_out = f"batch_rbh_{unique_id}.m8"
    tmp_subdir = f"tmp_rbh_batch_{unique_id}"
    
    valid_ids = set()

    try:
        def _ids_match(parent_name: str, target_name: str) -> bool:
            parent_base = extract_base_gene_id(parent_name).strip()
            target_base = extract_base_gene_id(target_name).strip()
            return (
                parent_base == target_base or
                parent_name == target_name or
                target_name == parent_name or
                parent_base in target_base or
                target_base in parent_base
            )

        def _best_ungapped_metrics(query_seq: str, target_seq: str) -> Tuple[float, float, float]:
            """
            Best ungapped overlap identity between two protein sequences.
            Returns (identity_pct, qcov_pct, tcov_pct).
            """
            q = query_seq or ""
            t = target_seq or ""
            if not q or not t:
                return 0.0, 0.0, 0.0

            best_identity = 0.0
            best_qcov = 0.0
            best_tcov = 0.0
            best_overlap = 0

            # offset < 0 means target starts before query
            for offset in range(-len(t) + 1, len(q)):
                q_start = max(0, offset)
                t_start = max(0, -offset)
                overlap = min(len(q) - q_start, len(t) - t_start)
                if overlap <= 0:
                    continue

                matches = 0
                for idx in range(overlap):
                    if q[q_start + idx] == t[t_start + idx]:
                        matches += 1

                identity = (matches / overlap) * 100.0
                qcov = (overlap / len(q)) * 100.0
                tcov = (overlap / len(t)) * 100.0

                if (
                    identity > best_identity or
                    (identity == best_identity and overlap > best_overlap)
                ):
                    best_identity = identity
                    best_qcov = qcov
                    best_tcov = tcov
                    best_overlap = overlap

            return best_identity, best_qcov, best_tcov

        def _fallback_rbh_local(records_local, db_for_convert, mapping) -> set:
            """
            Fallback RBH path for short proteins where mmseqs prefilter cannot run.
            """
            fallback_valid = set()
            tmp_home_faa = f"home_export_{uuid.uuid4().hex}.faa"
            try:
                subprocess.run(
                    ["mmseqs", "convert2fasta", db_for_convert, tmp_home_faa],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                home_records = [(hid, hseq) for _, hid, hseq in parse_fasta(tmp_home_faa)]
                if not home_records:
                    return fallback_valid

                for cand_id, cand_seq in records_local:
                    if cand_id not in mapping:
                        continue

                    best_target = None
                    best_identity = -1.0
                    best_qcov = 0.0
                    best_tcov = 0.0

                    for target_id, target_seq in home_records:
                        ident, qcov, tcov = _best_ungapped_metrics(cand_seq, target_seq)
                        if (
                            ident > best_identity or
                            (ident == best_identity and (qcov + tcov) > (best_qcov + best_tcov))
                        ):
                            best_identity = ident
                            best_target = target_id
                            best_qcov = qcov
                            best_tcov = tcov

                    if not best_target:
                        continue

                    parent = mapping[cand_id]
                    coverage_ok = (
                        best_qcov >= min_coverage * 100 and
                        best_tcov >= min_coverage * 100
                    )
                    identity_ok = best_identity >= min_identity
                    if _ids_match(parent, best_target) and coverage_ok and identity_ok:
                        fallback_valid.add(cand_id)

            finally:
                if os.path.exists(tmp_home_faa):
                    os.remove(tmp_home_faa)

            return fallback_valid

        def _candidate_to_pair(candidate):
            # Dict form used by pipeline internals
            if isinstance(candidate, dict):
                cid = candidate.get('id')
                cseq = candidate.get('seq')
                return cid, str(cseq) if cseq is not None else None

            # SeqRecord-like form used by unit tests
            cid = getattr(candidate, 'id', None)
            cseq_obj = getattr(candidate, 'seq', None)
            cseq = str(cseq_obj) if cseq_obj is not None else None
            return cid, cseq

        # Write FASTA using our utility
        records = []
        for c in candidates:
            cid, cseq = _candidate_to_pair(c)
            if cid and cseq:
                records.append((cid, cseq))

        if not records:
            return set()

        write_fasta(records, query_fasta)
            
        db_path = home_db
        if os.path.isdir(home_db):
            db_path = os.path.join(home_db, "db")
            
        cmd = [
            "mmseqs", "easy-search",
            query_fasta, db_path, rbh_out, tmp_subdir,
            "-e", str(evalue),
            "--format-output", "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen",
            "--max-seqs", "1", # Top hit only
            # Parameters, not the module-global `args` main() creates: this function has
            # no caller in the pipeline, and the global does not exist when it is imported.
            "--split-memory-limit", str(split_memory_limit),
            "-v", str(mmseqs_verbosity),
            "--threads", str(threads)
        ]

        mmseqs_ok = True
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            mmseqs_ok = False
            logger.warning(
                "RBH mmseqs easy-search failed; using local fallback "
                "(expected for very short peptides)."
            )
            valid_ids.update(_fallback_rbh_local(records, db_path, cand_map))

        seen_candidates = set()
        if mmseqs_ok and os.path.exists(rbh_out):
            with open(rbh_out) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 9: continue
                    
                    cand_id = parts[0]
                    seen_candidates.add(cand_id)
                    target_id = parts[1]
                    pident = float(parts[2])
                    qcov = float(parts[3]) if len(parts) > 3 else 100
                    tcov = float(parts[4]) if len(parts) > 4 else 100

                    # MMseqs coverage can be reported as fraction [0,1] or percent [0,100]
                    if qcov <= 1.0:
                        qcov *= 100.0
                    if tcov <= 1.0:
                        tcov *= 100.0
                    
                    # Enhanced validation
                    if cand_id not in cand_map:
                        continue
                    
                    parent = cand_map[cand_id]
                    
                    # Check 1: ID matching (exact or close)
                    ids_match = _ids_match(parent, target_id)
                    
                    # Check 2: Coverage (both query and target must be well-covered)
                    coverage_ok = (qcov >= min_coverage * 100 and 
                                  tcov >= min_coverage * 100)
                    
                    # Check 3: Identity must be reasonable
                    identity_ok = pident >= min_identity
                    
                    if ids_match and coverage_ok and identity_ok:
                        valid_ids.add(cand_id)
                    elif ids_match and not coverage_ok:
                        print(f"RBH: {cand_id} matches {target_id} but low coverage "
                              f"(qcov={qcov:.0f}%, tcov={tcov:.0f}%). Likely fragment/paralog.",
                              file=sys.stderr)
                    elif ids_match and not identity_ok:
                        logger.debug(f"RBH: {cand_id} matches {target_id} but very low identity "
                              f"({pident:.1f}%). Possible pseudogene.")

        # Peptides below the k-mer prefilter's reach get NO row rather than an error
        # (current MMseqs2 exits 0 with an empty table), so the fallback above never
        # ran for them and every short candidate was silently rejected.
        if mmseqs_ok:
            short_missed = [(cid, cseq) for cid, cseq in records
                            if cid not in seen_candidates and len(cseq) < 30]
            if short_missed:
                valid_ids.update(_fallback_rbh_local(short_missed, db_path, cand_map))

    except Exception as e:
        logger.error(f"RBH check failed: {e}")
        return set()
    finally:
        # Cleanup
        if os.path.exists(query_fasta): os.remove(query_fasta)
        if os.path.exists(rbh_out): os.remove(rbh_out)
        if os.path.exists(tmp_subdir):
             shutil.rmtree(tmp_subdir, ignore_errors=True)
                          
    return valid_ids


def is_valid_fallback(
    ordered_hit_count: int,
    qcov: float,
    single_aln: int,
    aln_total: int,
    best_bits: float,
    query_len: int,
    short_query_len: int,
    short_min_aln_aa: int,
    short_min_bits: float,
) -> bool:
    """A1 "the fumble" (TODO_JUN §A1, maps §1.2): decide if a GOI fallback model is emitted.

    The legacy gate vetoed on query COVERAGE, which is structurally low for a
    short precursor whose conserved mature peptide is only a fraction of the
    query (melittin: a 26-aa mature region in a 70-aa preprotein). A real
    high-bitscore divergent hit (Xylocopa b0, 50%/bits≈91) was rejected at
    qcov≈0.18 while a spurious high-coverage/low-bitscore other-block artifact
    (b1, 28%) passed at qcov=0.40 — so the GOI was emitted in the WRONG block.

    For short queries (< ``short_query_len``) we instead gate on absolute aligned
    length + bitscore (no coverage veto), so the real model is emitted in its own
    flanking-seeded block. Long queries keep the coverage floor (LY6/TP53
    unaffected); Tetramorium (76% cov) clears either path. ``short_query_len <= 0``
    reverts to the legacy coverage gate for every query.
    """
    if 0 < short_query_len and query_len < short_query_len:
        aln_for_gate = single_aln if ordered_hit_count == 1 else aln_total
        return aln_for_gate >= short_min_aln_aa and best_bits >= short_min_bits
    if ordered_hit_count == 1:
        return not (qcov < 0.25 or single_aln < 25)
    return not (qcov < 0.25 and aln_total < 35 and best_bits < 60.0)


def define_waves(genome_entries, args):
    """Partition genomes into execution waves (the expanding-database schedule).

    Three modes, in precedence order: ``--disable_wavefront`` (F4 control: one
    parallel wave, no iteration at all), rank/quantile binning (§A4, the default),
    and the legacy absolute-distance binning. Extracted from ``main()`` so the
    schedule is unit-testable.
    """
    if getattr(args, "disable_wavefront", False):
        # F4 control: one parallel wave => the query DB is never augmented between
        # genomes, so every target is searched with the initial DB only.
        waves = [list(genome_entries)]
        logger.info(
            "F4 wavefront DISABLED: all %d genome(s) in ONE parallel wave — "
            "no expanding-database iteration (control arm).", len(genome_entries)
        )
    elif args.rank_wave_binning:
        # §A4 (TODO_JUN): bin by phylo-distance RANK so a saturated metric
        # (tight clade, all dist≈0.99) still grades close→far.
        waves = assign_waves_by_rank(genome_entries)
        logger.info(
            "§A4 rank/quantile wave binning: ENABLED — "
            f"{len(waves)} graded waves (closest serial → farthest parallel), "
            "robust to saturated phylo-distance."
        )
    else:
        waves = []

        # IMPROVED WAVEFRONT STRATEGY (legacy absolute-distance binning):
        # - Closest genomes (dist < 0.05): Process strictly serially for maximum sensitivity
        # - Medium distance (0.05 - 0.15): Small waves (2-3 genomes)
        # - Distant genomes (> 0.15): Larger waves (can parallelize more)

        i = 0
        while i < len(genome_entries):
            curr = genome_entries[i]

            if curr['dist'] < 0.05:
                # Very close: Serial processing (wave of 1)
                waves.append([curr])
                i += 1
            elif curr['dist'] < 0.15:
                # Medium distance: Small waves of 2-3 genomes with similar distance
                wave = [curr]
                i += 1
                while i < len(genome_entries) and abs(genome_entries[i]['dist'] - curr['dist']) < 0.01:
                    wave.append(genome_entries[i])
                    i += 1
                    if len(wave) >= 3:  # Max 3 per wave for medium distance
                        break
                waves.append(wave)
            else:
                # Distant: Can parallelize more (waves of up to 5)
                wave = [curr]
                i += 1
                while i < len(genome_entries) and abs(genome_entries[i]['dist'] - curr['dist']) < 0.02:
                    wave.append(genome_entries[i])
                    i += 1
                    if len(wave) >= 5:  # Max 5 per wave for distant genomes
                        break
                waves.append(wave)
    return waves


def assign_waves_by_rank(genome_entries):
    """§A4 (TODO_JUN): rank/quantile wave binning.

    ``genome_entries`` must already be sorted closest-first (PHYLO_SORT order).
    The legacy absolute-distance binning collapses for a tight clade whose phylo
    distances saturate — e.g. all 19 bees land at dist≈0.99 after max-normalization
    — producing ~4 undifferentiated waves with no close→far gradient. Binning by
    RANK instead always grades: the closest genomes go in small (serial) waves so
    their GOI seeds reach ``expanded_db.faa`` before the divergent genomes are
    searched (the §A3 precondition for closest-first seeding), while distant
    genomes go in larger waves for parallelism. Pure function of the sorted order
    → deterministic (A0-safe). Returns a list of waves (each a list of entries).
    """
    n = len(genome_entries)
    waves = []
    i = 0
    while i < n:
        q = i / max(1, n - 1)   # rank-quantile in [0, 1]; 0.0 = closest
        if q < 0.10:
            size = 1            # closest tier: strictly serial (max seeding)
        elif q < 0.35:
            size = 2
        elif q < 0.70:
            size = 3
        else:
            size = 5            # farthest tier: parallelize
        waves.append(genome_entries[i:i + size])
        i += size
    return waves


def process_region_block(block_idx, block, hits, genome_seqs, db_sequences, genome_name, args, unique_id,
                         threads_per_job, native_annot_index=None):
    """
    Process a single synteny block: extract region, identify queries, run augmented search, and annotate.
    Returns: (list_of_new_genes, list_of_gff_lines)
    """
    chrom = block['chrom']
    if chrom not in genome_seqs:
        return [], []

    slen = len(genome_seqs[chrom])
    if getattr(args, 'adaptive_padding', True):
        padding = calculate_adaptive_padding(
            hits,
            block,
            default=args.region_padding,
            min_pad=args.padding_min,
            max_pad=args.padding_max
        )
    else:
        padding = args.region_padding

    # The neighbourhood of the block (anchors counted as support) and the DNA that is cut
    # out and searched; the second is the smaller one when distances are derived
    # (genome_scale_args: search_margin).
    hood_start = max(0, block['start'] - padding)
    hood_end = min(slen, block['end'] + padding)
    search_pad = padding
    if not getattr(args, 'adaptive_padding', True) and getattr(args, 'search_margin', None):
        search_pad = min(padding, int(args.search_margin))
    w_start = max(0, block['start'] - search_pad)
    w_end = min(slen, block['end'] + search_pad)
    subseq = genome_seqs[chrom][w_start:w_end]
    # Keep final GOI calls anchored near the flanking-defined block.
    core_margin = getattr(args, 'goi_margin', None) or args.gap_search_window
    core_start = max(0, block['start'] - core_margin)
    core_end = min(slen, block['end'] + core_margin)
    # DNA around a candidate's hits that the gene model may use: never less than the
    # gene's own span at home.
    goi_model_margin = max(args.gap_search_window,
                           int(getattr(args, 'goi_home_span', 0) * float(getattr(args, 'expansion', 1.0) or 1.0)))

    temp_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_block{block_idx}.fasta"
    query_mini_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_block{block_idx}_query.faa"
    write_fasta([("region_seq", subseq)], temp_fa)

    neighbourhood_hits = [
        h for h in hits
        if h['chrom'] == chrom and h['end'] > hood_start and h['start'] < hood_end
    ]
    relevant_hits = [h for h in neighbourhood_hits if h['end'] > w_start and h['start'] < w_end]
    block_flanking_support = len({
        extract_base_gene_id(h.get('query', ''))
        for h in neighbourhood_hits
        if h.get('query') and not is_goi_query_id(h.get('query', ''))
    })
    # §1 synteny ORDER: longest collinear run of the block's flanking genes in home order
    # (one representative position per flanking gene, sorted by TARGET position → home ranks
    # → LIS). Feeds _classify_goi_evidence so a scrambled paralog neighbourhood can't reach
    # HIGH on flanking COUNT alone. None when collinearity is unavailable (no home_rank /
    # --disable_synteny_collinearity) → the classifier waives the order gate (legacy behaviour).
    home_rank = getattr(args, 'home_rank', None)
    if home_rank:
        flank_first_pos: Dict[str, int] = {}
        for h in neighbourhood_hits:
            q = h.get('query', '')
            if q and not is_goi_query_id(q):
                gid = extract_base_gene_id(q)
                pos = int(h.get('start', 0) or 0)
                if gid not in flank_first_pos or pos < flank_first_pos[gid]:
                    flank_first_pos[gid] = pos
        ordered_ranks = [home_rank.get(gid) for gid, _pos
                         in sorted(flank_first_pos.items(), key=lambda kv: kv[1])]
        block_collinear_support, _dir = _longest_collinear_run(ordered_ranks)
    else:
        block_collinear_support = None
    unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)

    found_queries = []

    # Always include GOI-derived queries so blocks remain anchored to the target family.
    goi_ids = [k for k in db_sequences.keys() if is_goi_query_id(k)]
    full_length_goi = {gid for gid in goi_ids if '|frag_' not in gid}
    if full_length_goi:
        unique_queries.update(full_length_goi)
    else:
        unique_queries.update(goi_ids)

    # sorted(): `unique_queries` is a set of strings, whose iteration order changes with
    # the per-process PYTHONHASHSEED. That order used to decide which of two equal-length
    # GOI representatives survived the parent collapse below, so every block searched after
    # the first wavefront expansion flipped a coin between the home melittin and the
    # A. florea model (LRZ job 5778368: 17/20 GFFs differed, exactly the genomes searched
    # after wave 2).
    for query_id in sorted(unique_queries):
        if query_id in db_sequences:
            found_queries.append(db_sequences[query_id])
            continue
        base = extract_base_gene_id(query_id)
        if base in db_sequences:
            found_queries.append(db_sequences[base])

    if not found_queries:
        if os.path.exists(temp_fa):
            os.remove(temp_fa)
        return [], []

    # Region annotation should be driven by full-length GOI queries.
    # Flanking genes are used to define candidate regions, not to produce GOI models.
    full_goi_queries = [q for q in found_queries if is_full_length_goi_query_id(q.get('id', ''))]
    goi_queries = full_goi_queries if full_goi_queries else [q for q in found_queries if is_goi_query_id(q.get('id', ''))]
    if not goi_queries:
        if os.path.exists(temp_fa):
            os.remove(temp_fa)
        return [], []

    # One query per GOI is the home protein; the model of the relative searched most
    # recently is a second one, kept apart under its own parent id (select_goi_queries).
    home_goi_ids = getattr(args, 'home_goi_ids', None)
    relative_alias: Dict[str, str] = {}
    cand_marks: List[Tuple[int, str]] = []
    if home_goi_ids:
        goi_queries, goi_query_by_parent, relative_alias = select_goi_queries(
            db_sequences, home_goi_ids,
            use_relative=not getattr(args, 'disable_relative_query', False))
    else:
        goi_query_by_parent = collapse_goi_queries_by_parent(goi_queries)
        goi_queries = list(goi_query_by_parent.values())

    write_fasta([(q['id'], q['seq']) for q in goi_queries], query_mini_fa)

    annotated_records_raw = []
    valid_gff_lines = []
    raw_candidates = []
    flanking_candidates = []
    clean_gname = genome_name.replace('.', '_').replace('-', '_').replace(' ', '_')

    def _mRNA_attrs(base_attrs: Dict[str, Any], g_start: int, g_end: int, g_strand: str) -> str:
        native = lookup_native_annotation(
            native_annot_index or {},
            chrom=chrom,
            start=g_start,
            end=g_end,
            strand=g_strand,
        )
        return _compose_gff_attrs(base_attrs, native)

    def _emit_exon_model(parent_id, parent_query_seq, exons, model_protein, tag):
        """One miniprot gene model as a GOI candidate. ``tag`` names it inside the block
        (``l1`` for the first model of locus 1, ``l1m2`` for a second gene found in the
        same locus). Returns the model's exon intervals in window coordinates."""
        exons.sort(key=lambda e: e.get('qstart', 0))
        # miniprot's own translation of the model (in frame across
        # every intron phase and modelled frameshift); see
        # annotate_using_miniprot. Per-exon concatenation is the fallback.
        exon_protein = model_protein or ''.join(e['seq'] for e in exons)
        strand = exons[0].get('strand', '+')
        avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
        model_qcov = None
        if parent_query_seq:
            qmin_model = min(e.get('qstart', 1) for e in exons)
            qmax_model = max(e.get('qend', 1) for e in exons)
            model_qcov = (qmax_model - qmin_model + 1) / max(1, len(parent_query_seq))

        global_start = w_start + min(e['gstart'] for e in exons) + 1
        global_end = w_start + max(e['gend'] for e in exons)
        new_id = f"{parent_id}|{clean_gname}_b{block_idx}_{tag}_exon_ann"

        cand_gff = [
            (
                f"{chrom}\texon_annotation\tmRNA\t{global_start}\t{global_end}\t"
                f"{avg_pident:.1f}\t{strand}\t.\t"
                f"{_mRNA_attrs(_goi_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id}, evidence_type='exon_annotation', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, flanking_support=block_flanking_support, collinear_support=block_collinear_support, exons=exons), global_start, global_end, strand)}"
            )
        ]
        for eidx, exon in enumerate(exons, 1):
            exon_gs = w_start + exon['gstart'] + 1
            exon_ge = w_start + exon['gend']
            exon_phase = int(exon.get('phase', 0) or 0)
            attrs = f"ID={new_id}_CDS{eidx};Parent={new_id}"
            if exon.get('splice_acceptor'):
                attrs += f";SpliceAcceptor={exon['splice_acceptor']}"
            if exon.get('splice_donor'):
                attrs += f";SpliceDonor={exon['splice_donor']}"
            if exon.get('has_start_codon'):
                attrs += ";StartCodon=ATG"
            if exon.get('has_stop_codon'):
                attrs += ";StopCodon=yes"
            cand_gff.append(
                f"{chrom}\texon_annotation\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t{exon_phase}\t{attrs}"
            )

        raw_candidates.append({
            'start': global_start,
            'end': global_end,
            'score': avg_pident,
            'record': {
                'id': new_id,
                'seq': exon_protein,
                'description': (
                    f"coords:{global_start}-{global_end} "
                    f"parent:{parent_id} exons:{len(exons)} "
                    f"identity:{avg_pident:.1f}"
                )
            },
            'gff': cand_gff
        })
        return [(e['gstart'], e['gend']) for e in exons]

    def _emit_tandem_copy(parent_id, parent_query_seq, copy, locus_idx):
        """One single-segment copy (see detect_tandem_duplications) as a GOI candidate."""
        global_start = w_start + copy['gstart'] + 1
        global_end = w_start + copy['gend']
        strand = copy.get('strand', '+')
        copy_id = f"{copy['id']}|{clean_gname}_b{block_idx}_l{locus_idx}"
        # QW2 (docs/TODO_JUN.md): record the copy's TRUE query
        # coverage so a short high-identity local window (e.g. a
        # 23-aa hit on a 70-aa query) can't masquerade as a full
        # tandem copy.
        _qlen = len(parent_query_seq) if parent_query_seq else 0
        _qs, _qe = copy.get('qstart'), copy.get('qend')
        copy_qcov = None
        if _qlen and _qs is not None and _qe is not None:
            copy_qcov = max(0.0, min(1.0, (int(_qe) - int(_qs) + 1) / _qlen))
        raw_candidates.append({
            'start': global_start,
            'end': global_end,
            'score': copy.get('pident', 0),
            'record': {
                'id': copy_id,
                'seq': copy['seq'],
                'description': (
                    f"coords:{global_start}-{global_end} "
                    f"parent:{parent_id} tandem_copy "
                    f"identity:{copy.get('pident', 0):.1f}"
                )
            },
            'gff': [
                f"{chrom}\ttandem_copy\tgene\t{global_start}\t{global_end}\t"
                f"{copy.get('pident', 0):.1f}\t{strand}\t.\t"
                f"{_mRNA_attrs(_goi_feature_attrs({'ID': copy_id, 'Name': copy['id'], 'SynVoy_Parent': parent_id, 'Type': 'tandem_copy'}, evidence_type='tandem_copy', identity=copy.get('pident', 0), exon_count=1, query_cov=copy_qcov, flanking_support=block_flanking_support, collinear_support=block_collinear_support), global_start, global_end, strand)}"
            ]
        })
        return [(copy['gstart'], copy['gend'])]

    try:
        augmented_hits_mmseqs = run_augmented_search(
            temp_fa,
            goi_queries,
            f"{genome_name}_b{block_idx}",
            args,
            f"{unique_id}_b{block_idx}",
            threads_per_job
        )

        all_search_hits = list(augmented_hits_mmseqs) if augmented_hits_mmseqs else []

        # Include original block-supporting hits (converted to local coords).
        for hit in relevant_hits:
            if not is_goi_query_id(hit.get('query', '')):
                continue
            region_hit = dict(hit)
            if home_goi_ids:
                hit_query = hit.get('query', '')
                if hit_query in relative_alias:
                    region_hit['query'] = relative_alias[hit_query]
                elif hit_query not in home_goi_ids:
                    # A model of another genome: its positions are not positions in the
                    # home protein, and it is not a query of this window.
                    continue
            local_start = max(0, hit['start'] - w_start)
            local_end = min(w_end - w_start, hit['end'] - w_start)
            if local_end <= local_start:
                continue
            region_hit['start'] = local_start
            region_hit['end'] = local_end
            all_search_hits.append(region_hit)

        if all_search_hits:
            hits_by_gene = defaultdict(list)
            for hit in all_search_hits:
                local_start = max(0, int(hit.get('start', 0)))
                local_end = min(len(subseq), int(hit.get('end', 0)))
                if local_end <= local_start:
                    continue

                global_start = w_start + local_start
                global_end = w_start + local_end
                if global_end <= core_start or global_start >= core_end:
                    continue

                query_id = hit.get('query', '')
                if not query_id:
                    continue

                parent = query_id.split('|frag_')[0] if '|frag_' in query_id else query_id
                parent = extract_base_gene_id(parent)

                # Calculate offset for fragment or exon queries so qstart/qend
                # are expressed in full-protein coordinates.
                q_offset = 0
                if '|frag_' in query_id and FRAGMENT_SUPPORT:
                    try:
                        frag_info = parse_fragment_id(query_id)
                        q_offset = frag_info['start'] - 1
                    except Exception:
                        q_offset = 0
                elif '|exon_' in query_id:
                    # Build exon offset map from DB sequences on first encounter.
                    if not hasattr(run_augmented_search, '_goi_exon_offsets'):
                        run_augmented_search._goi_exon_offsets = {}
                    if parent not in run_augmented_search._goi_exon_offsets:
                        exon_seqs = []
                        for rid, rec in db_sequences.items():
                            if extract_base_gene_id(rid) == parent:
                                m = re.search(r'\|exon_(\d+)$', rid)
                                if m and rec.get('seq'):
                                    exon_seqs.append((int(m.group(1)), len(rec['seq'])))
                        exon_seqs.sort()
                        offmap = {}
                        off = 0
                        for enum, elen in exon_seqs:
                            offmap[enum] = off
                            off += elen
                        run_augmented_search._goi_exon_offsets[parent] = offmap
                    offmap = run_augmented_search._goi_exon_offsets.get(parent, {})
                    m = re.search(r'\|exon_(\d+)$', query_id)
                    if m:
                        enum = int(m.group(1))
                        q_offset = offmap.get(enum, 0)

                hits_by_gene[parent].append({
                    'qstart': hit.get('qstart', 1) + q_offset,
                    'qend': hit.get('qend', 100) + q_offset,
                    'gstart': local_start,
                    'gend': local_end,
                    'evalue': hit.get('evalue', 1),
                    'pident': hit.get('pident', 0),
                    'alnlen': hit.get('alnlen', 0),
                    'bits': hit.get('bits', 0),
                    'strand': hit.get('strand', '+'),
                    'chrom': chrom,
                    # Preserve ML-discovery provenance: PLM/structural hits carry a
                    # SYNTHETIC pident (cosine*100 / fident*100) and bits — these must
                    # NOT be read as sequence identity downstream (see
                    # _resolve_fallback_signals). The real embedding/TM-score signal
                    # rides along separately so classification can use the principled
                    # rescue path instead.
                    'method': hit.get('method'),
                    'embedding_similarity': hit.get('embedding_similarity'),
                    'structural_similarity': hit.get('structural_similarity'),
                })

            exon_parents = set()
            for q in found_queries:
                qid = q.get('id', '')
                if '|exon_' in qid:
                    exon_parents.add(extract_base_gene_id(qid))

            locus_gap = max(5000, min(100000, args.max_intron * 2),
                            int(getattr(args, 'goi_fallback_max_gap', 0) or 0))

            for parent_id, parent_hits in hits_by_gene.items():
                is_goi_parent = is_goi_query_id(parent_id)
                if not is_goi_parent:
                    continue
                cand_marks.append((len(raw_candidates), parent_id))

                parent_query = goi_query_by_parent.get(parent_id)
                if not parent_query:
                    continue
                parent_query_seq = parent_query.get('seq', '')
                if not parent_query_seq:
                    continue

                parent_loci = split_hits_into_loci(parent_hits, max_gap=locus_gap)
                logger.debug(f"[PARENT] parent_id={parent_id} is_goi={is_goi_parent} num_loci={len(parent_loci)} exon_parents={exon_parents}")
                if not parent_loci:
                    continue

                def _locus_rank(locus_hits):
                    # Aligned query residues (union), not the query SPAN: two short
                    # chance hits from opposite ends of the query spanned ~0.9 and
                    # outranked the real locus for one of the two kept slots.
                    qcov = _query_union_coverage(locus_hits, len(parent_query_seq))
                    best_bits = max(h.get('bits', 0) for h in locus_hits)
                    return (qcov, best_bits, len(locus_hits))

                # Keep only the strongest loci per parent to avoid noisy paralog floods.
                parent_loci = sorted(parent_loci, key=_locus_rank, reverse=True)[:2]

                for locus_idx, gene_hits in enumerate(parent_loci, start=1):
                    work_hits = gene_hits
                    if parent_id in exon_parents:
                        logger.debug(f"[EXON] {parent_id} starting with {len(work_hits)} hits")
                        work_hits = filter_exon_hits(work_hits, len(parent_query_seq))
                        logger.debug(f"[EXON] {parent_id} after filter: {len(work_hits)} hits")
                        if not work_hits:
                            logger.debug(f"[EXON] {parent_id} dropped entirely by filter_exon_hits!")
                            continue

                    exons = []
                    model_protein = None
                    extra_models = []       # further genes miniprot found in this locus
                    extra_copies = []       # single-segment copies where miniprot has no model
                    try:
                        is_tandem = False
                        tandem_copies = []
                        if is_goi_parent:
                            from annotate_goi_exons import detect_tandem_duplications
                            is_tandem, tandem_copies = detect_tandem_duplications(
                                work_hits, parent_query_seq, subseq, chrom
                            )

                        # Every gene model of the locus, not only the best: in a tandem
                        # array miniprot aligns the query to each copy it can reach, and the
                        # copies score within a per cent of each other, so keeping the first
                        # alone reported one copy picked almost at random. Single-segment
                        # copies from the hit pattern are kept where no model lies.
                        with maybe_quiet_streams(args.quiet_subtools):
                            models = annotate_models_from_hit_list(
                                work_hits,
                                parent_query_seq,
                                subseq,
                                chrom,
                                gap_search_window=goi_model_margin,
                                sensitive=is_goi_parent,
                                max_models=(1 if getattr(args, 'disable_multi_model', False)
                                            else int(getattr(args, 'goi_max_models_per_locus', 10) or 1)),
                            )
                        if models:
                            exons, model_protein = models[0]
                            extra_models = models[1:]
                            if is_tandem and tandem_copies:
                                spans = [(min(e['gstart'] for e in m_ex), max(e['gend'] for e in m_ex))
                                         for m_ex, _prot in models]
                                extra_copies = [
                                    c for c in tandem_copies
                                    if not any(min(c['gend'], e) - max(c['gstart'], b) > 0 for b, e in spans)
                                ]
                        elif is_tandem and tandem_copies:
                            exons = tandem_copies
                    except Exception as e:
                        import traceback
                        print(f"[ERROR] miniprot annotation failed for {parent_id} at {chrom}:{w_start}-{w_end}: {e}")
                        traceback.print_exc()
                        exons = []
                        extra_models = []
                        extra_copies = []

                    if not exons:
                        if work_hits:
                            strand_votes = {'+': 0.0, '-': 0.0}
                            for h in work_hits:
                                strand_votes[h.get('strand', '+')] += float(h.get('bits', 1))
                            strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

                            ordered_hits = sorted(
                                [h for h in work_hits if h.get('strand', '+') == strand],
                                key=lambda h: h.get('gstart', 0)
                            )
                            if ordered_hits:
                                # GOI-only branch (is_goi_parent gate above): bound the
                                # inter-exon gap by the home-GOI-derived max-intron so two
                                # far-apart spurious hits aren't stitched into one gene.
                                _goi_gap = getattr(args, 'goi_fallback_max_gap', args.max_intron)
                                ordered_hits = _longest_monotonic_query_chain(ordered_hits, strand, max_gap=_goi_gap)
                                if ordered_hits:
                                    query_len = max(1, len(parent_query_seq))
                                    qcov = _query_union_coverage(ordered_hits, query_len)
                                    aln_total = sum(max(0, int(h.get('alnlen', 0))) for h in ordered_hits)
                                    best_bits = max(float(h.get('bits', 0)) for h in ordered_hits)
                                    logger.debug(f"[FALLBACK] GOI={parent_id} qcov={qcov:.2f} aln_total={aln_total} bits={best_bits:.1f} len(hits)={len(ordered_hits)}")

                                    # A1 "the fumble": short queries gate on aligned length +
                                    # bitscore, not coverage (see is_valid_fallback).
                                    valid_fallback = is_valid_fallback(
                                        ordered_hit_count=len(ordered_hits),
                                        qcov=qcov,
                                        single_aln=int(ordered_hits[0].get('alnlen', 0)),
                                        aln_total=aln_total,
                                        best_bits=best_bits,
                                        query_len=query_len,
                                        short_query_len=args.fallback_short_query_len,
                                        short_min_aln_aa=args.fallback_short_min_aln_aa,
                                        short_min_bits=args.fallback_short_min_bits,
                                    )

                                    # Span sanity check: reject fallback hits whose
                                    # genomic span is wildly disproportionate to the
                                    # query protein length. For short peptides
                                    # (<150 aa) the intron:exon ratio is dominated
                                    # by introns, so a flat 30x multiplier rejects
                                    # valid multi-exon loci (e.g. melittin at 70 aa
                                    # can span tens of kb in bee genomes).
                                    if valid_fallback:
                                        # Short queries (<150 aa, e.g. melittin/toxin peptides):
                                        #   genomic span is dominated by introns; allow up to
                                        #   600x the CDS length so real multi-exon loci pass.
                                        # Long queries (>=150 aa, e.g. TP53/TP63/TP73):
                                        #   raised from 30x -> 50x after TP53 runs showed real
                                        #   TP63/TP73 loci spanning >100 kb being rejected.
                                        if query_len < 150:
                                            max_span_nt = max(100_000, query_len * 3 * 200)
                                        else:
                                            max_span_nt = max(5000, query_len * 3 * 50)
                                        gspan_min = min(h.get('gstart', 0) for h in ordered_hits)
                                        gspan_max = max(h.get('gend', 0) for h in ordered_hits)
                                        actual_span = gspan_max - gspan_min
                                        if actual_span > max_span_nt:
                                            logger.debug(
                                                f"[FALLBACK] REJECTED span={actual_span} "
                                                f"> max_span={max_span_nt} for {parent_id}"
                                            )
                                            valid_fallback = False

                                    logger.debug(f"[FALLBACK] valid_fallback={valid_fallback}")

                                    if valid_fallback:
                                        coding_frags = []
                                        cds_intervals = []
                                        for h in ordered_hits:
                                            hs = h.get('gstart', 0)
                                            he = h.get('gend', 0)
                                            if he <= hs: continue
                                            dna = subseq[hs:he]
                                            if strand == '-': dna = reverse_complement(dna)
                                            dna = dna[:len(dna) - (len(dna) % 3)]
                                            if len(dna) < 9: continue
                                            aa = translate(dna).replace('*', '')
                                            if not aa: continue
                                            coding_frags.append(aa)
                                            cds_intervals.append((hs, he))

                                        if coding_frags and cds_intervals:
                                            if strand == '-': coding_frags = list(reversed(coding_frags))
                                            flank_protein = ''.join(coding_frags)
                                            global_start = w_start + min(s for s, _ in cds_intervals) + 1
                                            global_end = w_start + max(e for _, e in cds_intervals)
                                            # Keep the model Identity derived ONLY from real
                                            # sequence hits; route ML-discovery (PLM/structural)
                                            # signal through the dedicated rescue branch so a
                                            # synthetic cosine/TM proxy can't masquerade as
                                            # sequence identity (see _resolve_fallback_signals).
                                            avg_pident, fb_embed_sim, fb_struct_sim, fb_ml_only = \
                                                _resolve_fallback_signals(ordered_hits)
                                            # For an ML-only locus the per-hit query coverage is
                                            # also synthetic (PLM hits report full coverage), so
                                            # suppress it — MEDIUM must come from the ML rescue
                                            # branch (similarity + flanking), not fake qcov.
                                            classify_qcov = None if fb_ml_only else qcov

                                            copy_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_fallback"

                                            cand_gff = [
                                                (
                                                    f"{chrom}\tfallback_hits\tmRNA\t{global_start}\t{global_end}\t"
                                                    f"{avg_pident:.1f}\t{strand}\t.\t"
                                                    f"{_mRNA_attrs(_goi_feature_attrs({'ID': copy_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'fallback_hit_span'}, evidence_type='fallback_hit_span', identity=avg_pident, exon_count=len(cds_intervals), query_cov=classify_qcov, flanking_support=block_flanking_support, embedding_similarity=fb_embed_sim, structural_similarity=fb_struct_sim, collinear_support=block_collinear_support), global_start, global_end, strand)}"
                                                )
                                            ]
                                            for eidx, (hs, he) in enumerate(cds_intervals, 1):
                                                exon_gs = w_start + hs + 1
                                                exon_ge = w_start + he
                                                cand_gff.append(
                                                    f"{chrom}\tfallback_hits\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\tID={copy_id}_CDS{eidx};Parent={copy_id}"
                                                )
                                            
                                            raw_candidates.append({
                                                'start': global_start,
                                                'end': global_end,
                                                'score': avg_pident,
                                                # A1: carry bitscore + flanking support so GOI
                                                # candidate ranking honours flanking-supported
                                                # bitscore, not raw coverage/identity.
                                                'bits': best_bits,
                                                'flanking_support': block_flanking_support,
                                                'record': {
                                                    'id': copy_id,
                                                    'seq': flank_protein,
                                                    'description': f"coords:{global_start}-{global_end} parent:{parent_id} type:fallback identity:{avg_pident:.1f}"
                                                },
                                                'gff': cand_gff
                                            })

                    if exons:
                        is_tandem_result = any(e.get('id', '').startswith('GOI_copy_') for e in exons)
                        if is_tandem_result:
                            for copy in exons:
                                _emit_tandem_copy(parent_id, parent_query_seq, copy, locus_idx)
                        else:
                            covered_intervals = _emit_exon_model(
                                parent_id, parent_query_seq, exons, model_protein, f"l{locus_idx}")
                            for model_no, (m_exons, m_protein) in enumerate(extra_models, start=2):
                                covered_intervals += _emit_exon_model(
                                    parent_id, parent_query_seq, m_exons, m_protein,
                                    f"l{locus_idx}m{model_no}")
                            for copy in extra_copies:
                                covered_intervals += _emit_tandem_copy(
                                    parent_id, parent_query_seq, copy, locus_idx)

                            # Rescue only for GOI-derived parents to prevent fallback explosions.
                            if is_goi_parent:
                                sorted_raw = sorted(work_hits, key=lambda h: h.get('gstart', h.get('start', 0)))
                                extra_count = 0
                                for hit in sorted_raw:
                                    h_start = hit.get('gstart', hit.get('start', 0))
                                    h_end = hit.get('gend', hit.get('end', 0))
                                    h_len = h_end - h_start
                                    if h_len < 1:
                                        continue

                                    is_covered = False
                                    for es, ee in covered_intervals:
                                        overlap = max(0, min(h_end, ee) - max(h_start, es))
                                        if overlap > 0.5 * h_len:
                                            is_covered = True
                                            break
                                    if is_covered:
                                        continue

                                    extra_count += 1
                                    region_dna = subseq[h_start:h_end]
                                    if hit.get('strand', '+') == '-':
                                        region_dna = reverse_complement(region_dna)
                                    region_dna = region_dna[:len(region_dna) - len(region_dna) % 3]
                                    if len(region_dna) < 6:
                                        continue

                                    hit_prot = translate(region_dna).replace('*', '')
                                    if not hit_prot:
                                        continue

                                    gs = w_start + h_start + 1
                                    ge = w_start + h_end
                                    hit_qspan = abs(hit.get('qend', 0) - hit.get('qstart', 0)) + 1
                                    hit_qcov = (hit_qspan / len(parent_query_seq)) if parent_query_seq else None
                                    raw_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_extra{extra_count}"
                                    raw_candidates.append({
                                        'start': gs,
                                        'end': ge,
                                        'score': hit.get('pident', 0),
                                        'record': {
                                            'id': raw_id,
                                            'seq': hit_prot,
                                            'description': f"coords:{gs}-{ge} parent:{parent_id} type:rescued_exon"
                                        },
                                        'gff': [
                                            f"{chrom}\trescued_exon\tmRNA\t{gs}\t{ge}\t{hit.get('pident',0):.1f}\t{hit.get('strand','+')}\t.\t{_mRNA_attrs(_goi_feature_attrs({'ID': raw_id, 'Name': parent_id, 'SynVoy_Parent': parent_id}, evidence_type='rescued_exon', identity=hit.get('pident', 0), exon_count=1, query_cov=hit_qcov, flanking_support=block_flanking_support, collinear_support=block_collinear_support), gs, ge, hit.get('strand','+'))}",
                                            f"{chrom}\trescued_exon\tCDS\t{gs}\t{ge}\t.\t{hit.get('strand','+')}\t0\tID={raw_id}_CDS1;Parent={raw_id}"
                                        ]
                                    })

                    # Structure inference failed: keep only a very strict GOI fallback set.
                    elif is_goi_parent:
                        locus_qcov = _query_union_coverage(work_hits, len(parent_query_seq))
                        best_bits = max(h.get('bits', 0) for h in work_hits) if work_hits else 0.0

                        # If miniprot is available, only allow fallback on very strong evidence.
                        if MINIPROT_AVAILABLE and (len(work_hits) < 2 or locus_qcov < 0.75 or best_bits < 60):
                            continue

                        # Fallback should be conservative; otherwise short/noisy hits explode.
                        min_fallback_identity = max(args.gap_min_identity, 90.0)
                        min_fallback_alnlen = max(
                            args.min_exon_alnlen,
                            int(max(15, len(parent_query_seq) * 0.6))
                        )
                        # Bit score first (was identity, then length, then bits).
                        sorted_hits = sorted(
                            work_hits,
                            key=lambda h: (h.get('bits', 0), h.get('pident', 0), h.get('alnlen', 0)),
                            reverse=True
                        )[:max(1, args.gap_max_hits)]

                        accepted_fallback = 0
                        for h_idx, hit in enumerate(sorted_hits, 1):
                            if hit.get('pident', 0) < min_fallback_identity:
                                continue
                            if hit.get('alnlen', 0) < min_fallback_alnlen:
                                continue
                            qspan = abs(hit.get('qend', 0) - hit.get('qstart', 0)) + 1
                            qcov = (qspan / len(parent_query_seq)) if len(parent_query_seq) > 0 else 0.0
                            if qcov < 0.6:
                                continue

                            g_s, g_e = hit['gstart'], hit['gend']
                            region_dna = subseq[g_s:g_e]
                            if hit.get('strand', '+') == '-':
                                region_dna = reverse_complement(region_dna)
                            region_dna = region_dna[:len(region_dna) - len(region_dna) % 3]
                            if len(region_dna) < 9:
                                continue

                            hit_protein = translate(region_dna).replace('*', '')
                            if not hit_protein:
                                continue

                            nt_s = w_start + g_s + 1
                            nt_e = w_start + g_e
                            qcov = (qspan / len(parent_query_seq)) if len(parent_query_seq) > 0 else 0.0
                            new_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_fallback_{h_idx}"
                            raw_candidates.append({
                                'start': nt_s,
                                'end': nt_e,
                                'score': hit.get('pident', 0),
                                'record': {
                                    'id': new_id,
                                    'seq': hit_protein,
                                    'description': f"coords:{nt_s}-{nt_e} parent:{parent_id} identity:{hit.get('pident', 0):.1f}"
                                },
                                'gff': [
                                    f"{chrom}\traw_hit\tmRNA\t{nt_s}\t{nt_e}\t{hit.get('pident', 0):.1f}\t{hit.get('strand', '+')}\t.\t{_mRNA_attrs(_goi_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id}, evidence_type='raw_hit', identity=hit.get('pident', 0), exon_count=1, query_cov=qcov, flanking_support=block_flanking_support, collinear_support=block_collinear_support), nt_s, nt_e, hit.get('strand', '+'))}",
                                    f"{chrom}\traw_hit\tCDS\t{nt_s}\t{nt_e}\t.\t{hit.get('strand', '+')}\t0\tID={new_id}_CDS1;Parent={new_id}"
                                ]
                            })
                            accepted_fallback += 1
                            if accepted_fallback >= 1:
                                break

        if relative_alias and raw_candidates:
            cand_parents = []
            mark_idx = 0
            for cand_idx in range(len(raw_candidates)):
                while mark_idx + 1 < len(cand_marks) and cand_marks[mark_idx + 1][0] <= cand_idx:
                    mark_idx += 1
                cand_parents.append(cand_marks[mark_idx][1] if cand_marks else '')
            raw_candidates[:] = rebase_relative_candidates(
                raw_candidates, cand_parents, goi_query_by_parent)

        # Flanking gene annotation pass (kept separate from GOI logic).
        flanking_parent_ids = {
            extract_base_gene_id(h.get('query', ''))
            for h in relevant_hits
            if h.get('query') and not is_goi_query_id(h.get('query', ''))
        }
        flanking_query_by_parent = build_flanking_query_by_parent(db_sequences, flanking_parent_ids)
        flanking_hits_by_parent = defaultdict(list)
        flanking_locus_gap = max(5000, min(100000, args.max_intron * 2))

        for hit in relevant_hits:
            query_id = hit.get('query', '')
            if not query_id or is_goi_query_id(query_id):
                continue

            parent_id = extract_base_gene_id(query_id)
            if parent_id not in flanking_query_by_parent:
                continue

            local_start = max(0, int(hit.get('start', 0)) - w_start)
            local_end = min(len(subseq), int(hit.get('end', 0)) - w_start)
            if local_end <= local_start:
                continue

            global_start = w_start + local_start
            global_end = w_start + local_end
            if global_end <= core_start or global_start >= core_end:
                continue

            flanking_hits_by_parent[parent_id].append({
                'query': query_id,
                'qstart': hit.get('qstart', 1),
                'qend': hit.get('qend', 1),
                'gstart': local_start,
                'gend': local_end,
                'evalue': hit.get('evalue', 1),
                'pident': hit.get('pident', 0),
                'alnlen': hit.get('alnlen', 0),
                'bits': hit.get('bits', 0),
                'strand': hit.get('strand', '+'),
                'chrom': chrom
            })

        # Remap per-exon qstart/qend to reconstructed full-protein coordinates.
        # Without this, a hit from gene-TOP1MT|exon_5 (qstart=1,qend=37) would
        # be compared against the 601aa full protein → qcov=6% → wrongly fails
        # coverage thresholds.  After remapping, it becomes qstart=170,qend=207
        # which correctly reflects its position in the reconstructed protein.
        for parent_id, parent_hits in flanking_hits_by_parent.items():
            query_rec = flanking_query_by_parent.get(parent_id)
            if not query_rec:
                continue
            exon_offsets = query_rec.get('exon_offsets', {})
            if not exon_offsets:
                continue
            for h in parent_hits:
                qid = h.get('query', '')
                m = re.search(r'\|exon_(\d+)$', qid)
                if not m:
                    continue
                exon_num = int(m.group(1))
                if exon_num not in exon_offsets:
                    continue
                offset, _exon_len = exon_offsets[exon_num]
                # Shift from exon-local to full-protein coordinates (1-based)
                h['qstart'] = h.get('qstart', 1) + offset
                h['qend'] = h.get('qend', 1) + offset

        for parent_id, parent_hits in flanking_hits_by_parent.items():
            query_rec = flanking_query_by_parent.get(parent_id)
            if not query_rec:
                continue
            parent_query_seq = query_rec.get('seq', '')
            if not parent_query_seq:
                continue

            parent_loci = split_hits_into_loci(
                parent_hits, max_gap=_gene_limit(args, parent_id, flanking_locus_gap))
            if not parent_loci:
                continue

            # Keep strongest loci first; flanking is anchor support, not exhaustive paralog mining.
            parent_loci = sorted(
                parent_loci,
                key=lambda hs: (max(h.get('bits', 0) for h in hs), len(hs)),
                reverse=True
            )[:2]

            for locus_idx, work_hits in enumerate(parent_loci, start=1):
                model_protein = None
                try:
                    with maybe_quiet_streams(args.quiet_subtools):
                        exons, model_protein = annotate_exons_from_hit_list(
                            work_hits,
                            parent_query_seq,
                            subseq,
                            chrom,
                            gap_search_window=max(args.gap_search_window,
                                                  _gene_limit(args, parent_id, 0)),
                            sensitive=False
                        )
                except Exception as e:
                    logger.warning(f"[{genome_name}] Flanking exon annotation failed for {parent_id} locus {locus_idx}: {e}")
                    exons = []

                if not exons:
                    # Conservative fallback: use hit-defined exon blocks from first to last hit.
                    if not work_hits:
                        continue

                    strand_votes = {'+': 0, '-': 0}
                    for h in work_hits:
                        strand_votes[h.get('strand', '+')] += 1
                    strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

                    ordered_hits = sorted(
                        [h for h in work_hits if h.get('strand', '+') == strand],
                        key=lambda h: h.get('gstart', 0)
                    )
                    if not ordered_hits:
                        continue

                    # Keep only the longest chain consistent with query-order direction.
                    ordered_hits = _longest_monotonic_query_chain(
                        ordered_hits, strand, max_gap=_gene_limit(args, parent_id, args.max_intron))
                    if not ordered_hits:
                        continue

                    # Require minimum support before emitting coarse flanking fallback models.
                    qmin = min(min(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                    qmax = max(max(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                    query_len = max(1, len(parent_query_seq))
                    qcov = (qmax - qmin + 1) / query_len
                    aln_total = sum(max(0, int(h.get('alnlen', 0))) for h in ordered_hits)
                    best_bits = max(float(h.get('bits', 0)) for h in ordered_hits)

                    # Single-hit fallback must look like substantial coverage.
                    if len(ordered_hits) == 1:
                        min_single_qcov = 0.45
                        min_single_aln = max(30, int(0.35 * query_len))
                        if qcov < min_single_qcov or int(ordered_hits[0].get('alnlen', 0)) < min_single_aln:
                            continue
                    else:
                        # Multi-hit chain fallback: require meaningful overall support.
                        if qcov < 0.35 and aln_total < max(45, int(0.45 * query_len)) and best_bits < 80.0:
                            continue

                    coding_frags = []
                    cds_intervals = []
                    for h in ordered_hits:
                        hs = h.get('gstart', 0)
                        he = h.get('gend', 0)
                        if he <= hs:
                            continue
                        dna = subseq[hs:he]
                        if strand == '-':
                            dna = reverse_complement(dna)
                        dna = dna[:len(dna) - (len(dna) % 3)]
                        if len(dna) < 9:
                            continue
                        aa = translate(dna).replace('*', '')
                        if not aa:
                            continue
                        coding_frags.append(aa)
                        cds_intervals.append((hs, he))

                    if not coding_frags or not cds_intervals:
                        continue

                    if strand == '-':
                        coding_frags = list(reversed(coding_frags))

                    flank_protein = ''.join(coding_frags)
                    global_start = w_start + min(s for s, _ in cds_intervals) + 1
                    global_end = w_start + max(e for _, e in cds_intervals)
                    avg_pident = _pooled_identity(ordered_hits)
                    new_id = f"{parent_id}|{clean_gname}_b{block_idx}_fl{locus_idx}_flank_hits"

                    cand_gff = [
                        (
                            f"{chrom}\tflanking_hits\tmRNA\t{global_start}\t{global_end}\t"
                            f"{avg_pident:.1f}\t{strand}\t.\t"
                            f"{_mRNA_attrs(_flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'flanking_hit_span'}, evidence_type='flanking_hit_span', identity=avg_pident, exon_count=len(cds_intervals), query_cov=qcov, context='candidate_region_anchor'), global_start, global_end, strand)}"
                        )
                    ]
                    for eidx, (hs, he) in enumerate(cds_intervals, 1):
                        exon_gs = w_start + hs + 1
                        exon_ge = w_start + he
                        cand_gff.append(
                            f"{chrom}\tflanking_hits\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\tID={new_id}_CDS{eidx};Parent={new_id}"
                        )

                    flanking_candidates.append({
                        'start': global_start,
                        'end': global_end,
                        'score': avg_pident,
                        'record': {
                            'id': new_id,
                            'seq': flank_protein,
                            'description': (
                                f"coords:{global_start}-{global_end} parent:{parent_id} "
                                f"type:flanking_hit_span identity:{avg_pident:.1f}"
                            )
                        },
                        'gff': cand_gff
                    })
                    continue

                exons.sort(key=lambda e: e.get('qstart', 0))
                # miniprot's in-frame model translation (see annotate_using_miniprot).
                flank_protein = model_protein or ''.join(e.get('seq', '') for e in exons if e.get('seq'))
                if not flank_protein:
                    continue

                strand = exons[0].get('strand', '+')
                avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
                model_qcov = None
                if parent_query_seq:
                    qmin_model = min(e.get('qstart', 1) for e in exons)
                    qmax_model = max(e.get('qend', 1) for e in exons)
                    model_qcov = (qmax_model - qmin_model + 1) / max(1, len(parent_query_seq))
                global_start = w_start + min(e['gstart'] for e in exons) + 1
                global_end = w_start + max(e['gend'] for e in exons)
                new_id = f"{parent_id}|{clean_gname}_b{block_idx}_fl{locus_idx}_flank_ann"

                cand_gff = [
                    (
                        f"{chrom}\tflanking_annotation\tmRNA\t{global_start}\t{global_end}\t"
                        f"{avg_pident:.1f}\t{strand}\t.\t"
                        f"{_mRNA_attrs(_flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'flanking_miniprot'}, evidence_type='flanking_miniprot', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, context='candidate_region_anchor'), global_start, global_end, strand)}"
                    )
                ]
                for eidx, exon in enumerate(exons, 1):
                    exon_gs = w_start + exon['gstart'] + 1
                    exon_ge = w_start + exon['gend']
                    exon_phase = int(exon.get('phase', 0) or 0)
                    cand_gff.append(
                        f"{chrom}\tflanking_annotation\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t{exon_phase}\tID={new_id}_CDS{eidx};Parent={new_id}"
                    )

                flanking_candidates.append({
                    'start': global_start,
                    'end': global_end,
                    'score': avg_pident,
                    'record': {
                        'id': new_id,
                        'seq': flank_protein,
                        'description': (
                            f"coords:{global_start}-{global_end} parent:{parent_id} "
                            f"exons:{len(exons)} identity:{avg_pident:.1f}"
                        )
                    },
                    'gff': cand_gff
                })

    except Exception as e:
        import traceback
        logger.warning(
            f"[{genome_name}] Block {block_idx} processing failed "
            f"({chrom}:{block.get('start','?')}-{block.get('end','?')}): {e}\n"
            f"  {traceback.format_exc().strip().splitlines()[-1]}"
        )

    if os.path.exists(temp_fa):
        os.remove(temp_fa)
    if os.path.exists(query_mini_fa):
        os.remove(query_mini_fa)

    # ── PLM re-ranking: batch-embed GOI model proteins, compute embedding
    #    similarity, re-run classification where PLM boost applies, and
    #    update GFF attributes accordingly. ──
    if (GPU_IN_WORKER and getattr(args, 'enable_plm_search', False)
            and PLM_IMPORT_OK and check_plm_available()
            and raw_candidates):
        try:
            goi_emb_path = os.path.join(args.output_dir, "goi_embeddings.npz")
            if os.path.exists(goi_emb_path):
                goi_embeddings = load_embeddings(goi_emb_path)
                goi_cand_seqs = [
                    (cand['record']['id'], cand['record']['seq'])
                    for cand in raw_candidates
                    if cand['record'].get('seq') and len(cand['record']['seq']) >= 5
                ]
                if goi_cand_seqs and goi_embeddings:
                    cand_sims = compute_candidate_similarities(
                        goi_cand_seqs, goi_embeddings,
                        device=getattr(args, 'plm_device', 'cpu'),
                    )
                    boosted = 0
                    for cand in raw_candidates:
                        cid = cand['record']['id']
                        sim = cand_sims.get(cid)
                        if sim is None:
                            continue
                        cand['embedding_similarity'] = sim

                        if not cand.get('gff'):
                            continue
                        mRNA_line = cand['gff'][0]

                        # Parse existing classification from GFF attrs to re-run
                        # with embedding_similarity.
                        attrs_str = mRNA_line.split('\t')[-1] if '\t' in mRNA_line else ""
                        attr_kv = {}
                        for pair in attrs_str.split(';'):
                            if '=' in pair:
                                k, v = pair.split('=', 1)
                                attr_kv[k.strip()] = v.strip()

                        old_conf = attr_kv.get('Confidence', '')
                        ev_type = attr_kv.get('EvidenceType', '')
                        try:
                            ident = float(attr_kv.get('Identity', 0))
                        except (ValueError, TypeError):
                            ident = 0.0
                        try:
                            exons = int(attr_kv.get('Exons', 1))
                        except (ValueError, TypeError):
                            exons = 1
                        try:
                            qcov = float(attr_kv.get('QueryCoverage', 0))
                        except (ValueError, TypeError):
                            qcov = None
                        try:
                            fl_sup = int(attr_kv.get('BlockFlankingSupport', 0))
                        except (ValueError, TypeError):
                            fl_sup = 0

                        # Re-classify with embedding similarity
                        new_conf, new_class, new_reason = _classify_goi_evidence(
                            evidence_type=ev_type,
                            identity=ident,
                            exon_count=exons,
                            query_cov=qcov,
                            flanking_support=fl_sup,
                            embedding_similarity=sim,
                        )

                        # Update GFF attrs if confidence was boosted
                        updated_attrs = attrs_str + f";EmbeddingSimilarity={sim:.3f}"
                        if new_conf != old_conf and old_conf in ('LOW', 'MEDIUM'):
                            updated_attrs = updated_attrs.replace(
                                f"Confidence={old_conf}", f"Confidence={new_conf}"
                            )
                            old_reason = attr_kv.get('InferenceReason', '')
                            if old_reason:
                                updated_attrs = updated_attrs.replace(
                                    f"InferenceReason={old_reason}",
                                    f"InferenceReason={new_reason}",
                                )
                            old_class = attr_kv.get('GOIClass', '')
                            if old_class:
                                updated_attrs = updated_attrs.replace(
                                    f"GOIClass={old_class}",
                                    f"GOIClass={new_class}",
                                )
                            boosted += 1

                        # Rebuild the mRNA line with updated attributes
                        cols = mRNA_line.split('\t')
                        if len(cols) >= 9:
                            cols[8] = updated_attrs
                            cand['gff'][0] = '\t'.join(cols)

                    if boosted:
                        print(
                            f"[{genome_name}] PLM re-ranking: boosted {boosted}/{len(raw_candidates)} "
                            f"candidate(s) via embedding similarity.",
                            flush=True,
                        )
        except Exception as plm_rerank_err:
            logger.debug(
                f"[{genome_name}] PLM re-ranking failed: {plm_rerank_err}"
            )

    # ── Structural re-ranking: fold GOI candidate proteins with ESMFold,
    #    compare against GOI structures via Foldseek, re-run classification
    #    where structural boost applies, and update GFF attributes. ──
    if (GPU_IN_WORKER and getattr(args, 'enable_structural_search', False)
            and STRUCTURAL_IMPORT_OK and check_structural_search_available()
            and raw_candidates):
        try:
            goi_struct_index = os.path.join(args.output_dir, "goi_structure_index.tsv")
            if os.path.exists(goi_struct_index):
                goi_structures = load_structure_index(goi_struct_index)
                goi_cand_seqs = [
                    (cand['record']['id'], cand['record']['seq'])
                    for cand in raw_candidates
                    if cand['record'].get('seq') and len(cand['record']['seq']) >= 10
                ]
                if goi_cand_seqs and goi_structures:
                    struct_rerank_dir = os.path.join(
                        args.output_dir, f"tmp_struct_rerank_{genome_name}"
                    )
                    os.makedirs(struct_rerank_dir, exist_ok=True)

                    cand_tm_scores = compute_candidate_structural_similarities(
                        goi_cand_seqs, goi_structures,
                        output_dir=struct_rerank_dir,
                        device=getattr(args, 'structural_device', 'cpu'),
                        max_length=getattr(args, 'structural_max_length', 700),
                        threads=getattr(args, 'threads', 1),
                    )

                    boosted = 0
                    for cand in raw_candidates:
                        cid = cand['record']['id']
                        tm = cand_tm_scores.get(cid)
                        if tm is None:
                            continue
                        cand['structural_similarity'] = tm

                        if not cand.get('gff'):
                            continue
                        mRNA_line = cand['gff'][0]

                        attrs_str = mRNA_line.split('\t')[-1] if '\t' in mRNA_line else ""
                        attr_kv = {}
                        for pair in attrs_str.split(';'):
                            if '=' in pair:
                                k, v = pair.split('=', 1)
                                attr_kv[k.strip()] = v.strip()

                        old_conf = attr_kv.get('Confidence', '')
                        ev_type = attr_kv.get('EvidenceType', '')
                        try:
                            ident = float(attr_kv.get('Identity', 0))
                        except (ValueError, TypeError):
                            ident = 0.0
                        try:
                            exons = int(attr_kv.get('Exons', 1))
                        except (ValueError, TypeError):
                            exons = 1
                        try:
                            qcov = float(attr_kv.get('QueryCoverage', 0))
                        except (ValueError, TypeError):
                            qcov = None
                        try:
                            fl_sup = int(attr_kv.get('BlockFlankingSupport', 0))
                        except (ValueError, TypeError):
                            fl_sup = 0
                        try:
                            emb_sim = float(attr_kv.get('EmbeddingSimilarity', 0))
                        except (ValueError, TypeError):
                            emb_sim = None

                        new_conf, new_class, new_reason = _classify_goi_evidence(
                            evidence_type=ev_type,
                            identity=ident,
                            exon_count=exons,
                            query_cov=qcov,
                            flanking_support=fl_sup,
                            embedding_similarity=emb_sim if emb_sim else None,
                            structural_similarity=tm,
                        )

                        updated_attrs = attrs_str + f";StructuralSimilarity={tm:.3f}"
                        if new_conf != old_conf and old_conf in ('LOW', 'MEDIUM'):
                            updated_attrs = updated_attrs.replace(
                                f"Confidence={old_conf}", f"Confidence={new_conf}"
                            )
                            old_reason = attr_kv.get('InferenceReason', '')
                            if old_reason:
                                updated_attrs = updated_attrs.replace(
                                    f"InferenceReason={old_reason}",
                                    f"InferenceReason={new_reason}",
                                )
                            old_class = attr_kv.get('GOIClass', '')
                            if old_class:
                                updated_attrs = updated_attrs.replace(
                                    f"GOIClass={old_class}",
                                    f"GOIClass={new_class}",
                                )
                            boosted += 1

                        cols = mRNA_line.split('\t')
                        if len(cols) >= 9:
                            cols[8] = updated_attrs
                            cand['gff'][0] = '\t'.join(cols)

                    if boosted:
                        print(
                            f"[{genome_name}] Structural re-ranking: boosted {boosted}/{len(raw_candidates)} "
                            f"candidate(s) via TM-score.",
                            flush=True,
                        )

                    if os.path.exists(struct_rerank_dir):
                        shutil.rmtree(struct_rerank_dir, ignore_errors=True)
        except Exception as struct_rerank_err:
            logger.debug(
                f"[{genome_name}] Structural re-ranking failed: {struct_rerank_err}"
            )

    found_ids = set()

    # Greedy non-overlapping selection for GOI/tandem candidates, by bit score against
    # the home protein (score_goi_candidates, goi_candidate_order).
    if raw_candidates:
        score_goi_candidates(raw_candidates, goi_query_by_parent)
        raw_candidates.sort(key=goi_candidate_order)

        final_candidates = []
        accepted_intervals = []

        for cand in raw_candidates:
            c_start, c_end = cand['start'], cand['end']
            c_len = c_end - c_start
            if c_len <= 0:
                continue

            is_overlapping = False
            for a_start, a_end in accepted_intervals:
                overlap = max(0, min(c_end, a_end) - max(c_start, a_start))
                if overlap > 0.5 * c_len:
                    is_overlapping = True
                    break

            if not is_overlapping:
                final_candidates.append(cand)
                accepted_intervals.append((c_start, c_end))

        for cand in final_candidates:
            rec = cand['record']
            if rec['id'] in found_ids:
                continue
            found_ids.add(rec['id'])
            annotated_records_raw.append(rec)
            valid_gff_lines.extend(cand['gff'])

    # Add flanking candidates without interfering with GOI selection.
    if flanking_candidates:
        flanking_candidates.sort(
            key=lambda x: (x.get('score', 0), x['end'] - x['start'],
                           x.get('start', 0), x.get('record', {}).get('id', '')),
            reverse=True
        )
        flanking_intervals = []
        for cand in flanking_candidates:
            rec = cand['record']
            if rec['id'] in found_ids:
                continue

            c_start, c_end = cand['start'], cand['end']
            c_len = c_end - c_start
            if c_len <= 0:
                continue

            is_overlapping = False
            for a_start, a_end in flanking_intervals:
                overlap = max(0, min(c_end, a_end) - max(c_start, a_start))
                if overlap > 0.5 * c_len:
                    is_overlapping = True
                    break
            if is_overlapping:
                continue

            flanking_intervals.append((c_start, c_end))
            found_ids.add(rec['id'])
            annotated_records_raw.append(rec)
            valid_gff_lines.extend(cand['gff'])

    return annotated_records_raw, valid_gff_lines


def merge_synteny_blocks(blocks, padding):
    """
    Merge overlapping synteny blocks to prevent redundant searches.
    
    If the search windows (block +/- padding) overlap, the blocks are merged.
    """
    if not blocks:
        return []

    # Sort by chrom, then start
    # Ensure all blocks have 'chrom' and 'start'
    valid_blocks = [b for b in blocks if 'chrom' in b and 'start' in b]
    if len(valid_blocks) < len(blocks):
        print(f"Warning: Dropped {len(blocks) - len(valid_blocks)} blocks missing coordinates during merge.")
    
    valid_blocks.sort(key=lambda b: (b['chrom'], b['start']))

    merged = []
    if not valid_blocks:
        return []
        
    current_block = valid_blocks[0]
    current_chrom = current_block['chrom']
    # Effective end of the search region for the current block
    current_search_end = current_block['end'] + padding

    for i in range(1, len(valid_blocks)):
        next_block = valid_blocks[i]
        next_search_start = next_block['start'] - padding
        
        # Check for overlap
        if next_block['chrom'] == current_chrom and next_search_start <= current_search_end:
            # Merge!
            # 1. Update core coordinates to cover both blocks
            current_block['start'] = min(current_block['start'], next_block['start'])
            current_block['end'] = max(current_block['end'], next_block['end'])
            
            # 2. Update stats: the merged block holds the UNION of both parts' genes.
            # The counts used to be summed and `genes` kept the first part's list, so
            # two loci of ONE flanking gene less than 2 x padding apart read as a
            # 2-gene block and passed min_block_genes (local luciferase + yeast runs,
            # 2026-10-06: 89 of 544 retained blocks held a single distinct gene), and
            # the proxy / cap checks never saw the genes of the later parts.
            merged_genes = set(current_block.get('genes') or []) | set(next_block.get('genes') or [])
            if merged_genes:
                current_block['genes'] = sorted(merged_genes)
                current_block['genes_count'] = len(merged_genes)
            else:
                # Blocks without a gene list (older callers): keep the summed count.
                current_block['genes_count'] = (current_block.get('genes_count', 0)
                                                + next_block.get('genes_count', 0))
            # Max score (best local evidence)
            current_block['score'] = max(current_block.get('score', 0), next_block.get('score', 0))
            # Preserve the stronger collinear-run length so the §1m cap re-sort
            # (process_single_genome) still ranks merged blocks correctly.
            current_block['collinear_chain_len'] = max(
                current_block.get('collinear_chain_len', 0),
                next_block.get('collinear_chain_len', 0))
            current_block['bridged'] = bool(
                current_block.get('bridged') or next_block.get('bridged'))
            
            # 3. Combine hits if present
            if 'hits' in next_block:
                if 'hits' not in current_block:
                    current_block['hits'] = []
                current_block['hits'].extend(next_block['hits'])
            
            # 4. Extend the search horizon
            current_search_end = max(current_search_end, next_block['end'] + padding)
            
            # Note: We keep modifying 'current_block' in place
        else:
            # No overlap, finalize current and move to next
            merged.append(current_block)
            current_block = next_block
            current_chrom = current_block['chrom']
            current_search_end = current_block['end'] + padding
            
    merged.append(current_block)
    return merged


def run_mmseqs_easy_search_with_retries(
    db_path: str,
    genome_path: str,
    hits_file: str,
    tmp_dir: str,
    args,
    threads_per_job: int,
    genome_name: str,
) -> Tuple[bool, str, bool]:
    """
    Run mmseqs easy-search with conservative retries for resource failures.

    Returns:
      (success, details, resource_related_failure)
    """
    base_attempts = [
        {
            "label": "primary",
            "threads": max(1, int(threads_per_job)),
            "sens": float(args.mmseqs_sens),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "lowmem_retry_1",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 7.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "lowmem_retry_2",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 6.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
    ]

    attempts = []
    seen = set()
    for a in base_attempts:
        key = (a["threads"], a["sens"], a["split"])
        if key in seen:
            continue
        seen.add(key)
        attempts.append(a)

    last_details = ""
    last_resource_fail = False
    os.makedirs(tmp_dir, exist_ok=True)

    # Create query and target databases ONCE (outside retry loop).
    query_db = os.path.join(tmp_dir, "queryDB")
    target_db = os.path.join(tmp_dir, "targetDB")
    fmt_output = "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    # -- createdb for query (protein FASTA) --
    proc_qdb = subprocess.run(
        ["mmseqs", "createdb", db_path, query_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_qdb.returncode != 0:
        details = _format_mmseqs_failure(proc_qdb)
        logger.warning(f"[{genome_name}] mmseqs createdb (query) failed: {details}")
        return False, details, False

    # -- createdb for target (nucleotide genome) --
    proc_tdb = subprocess.run(
        ["mmseqs", "createdb", genome_path, target_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_tdb.returncode != 0:
        details = _format_mmseqs_failure(proc_tdb)
        logger.warning(f"[{genome_name}] mmseqs createdb (target) failed: {details}")
        return False, details, False

    for i, attempt in enumerate(attempts, start=1):
        if os.path.exists(hits_file):
            try:
                os.remove(hits_file)
            except OSError:
                pass

        attempt_tmp = os.path.join(tmp_dir, f"attempt_{i}")
        os.makedirs(attempt_tmp, exist_ok=True)
        result_db = os.path.join(attempt_tmp, "resultDB")

        # -- mmseqs search (protein query vs translated nucleotide target) --
        search_cmd = [
            "mmseqs", "search",
            query_db, target_db, result_db, attempt_tmp,
            "--search-type", "2",
            "--threads", str(attempt["threads"]),
            "-s", str(attempt["sens"]),
            "-e", str(args.evalue),
            "--split-memory-limit", str(attempt["split"]),
            "-v", str(args.mmseqs_verbosity),
        ]

        search_proc = subprocess.run(
            search_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

        if search_proc.returncode != 0:
            details = _format_mmseqs_failure(search_proc)
            resource_fail = _is_mmseqs_resource_failure(details)
            last_details = details
            last_resource_fail = resource_fail
            logger.warning(
                f"[{genome_name}] MMseqs search attempt {i}/{len(attempts)} failed "
                f"({attempt['label']}: threads={attempt['threads']}, "
                f"s={attempt['sens']}, split={attempt['split']}): {details}"
            )
            if not resource_fail:
                break
            continue

        # -- mmseqs convertalis → m8 output --
        conv_cmd = [
            "mmseqs", "convertalis",
            query_db, target_db, result_db, hits_file,
            "--format-output", fmt_output,
            "-v", str(args.mmseqs_verbosity),
        ]

        conv_proc = subprocess.run(
            conv_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

        if conv_proc.returncode != 0:
            details = _format_mmseqs_failure(conv_proc)
            last_details = details
            logger.warning(f"[{genome_name}] mmseqs convertalis failed: {details}")
            break

        if i > 1:
            logger.warning(
                f"[{genome_name}] MMseqs succeeded on retry {i}/{len(attempts)} "
                f"({attempt['label']}: threads={attempt['threads']}, "
                f"s={attempt['sens']}, split={attempt['split']})."
            )
        return True, "", False

    return False, last_details, last_resource_fail


_SEED_QUALIFYING_CLASSES = {"confident_goi", "probable_goi"}
_TREE_QUALIFYING_CONFIDENCE = {"HIGH", "MEDIUM"}


def _safe_float(value, default=None):
    """Parse a GFF-attr string ("" / "0.37" / "5") to float; default on failure."""
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _classify_goi_for_seed_and_tree(
    all_genes, feature_meta, is_goi_query_id_fn,
    seed_on_flanking_support=False,
    seed_min_flanking=2,
    seed_min_qcov=0.5,
):
    """Split GOI-role genes into (seed, tree_extra, suppressed_count).

    Seed genes feed `expanded_db.faa` (used as the next wave's MMseqs2 query DB
    AND historically as the tree input). Strict filter: HIGH/MEDIUM confidence
    AND goi_class in {confident_goi, probable_goi}.

    Tree-extra genes feed `goi_for_tree.faa` *in addition to* the seed list:
    HIGH/MEDIUM-confidence GOI hits whose goi_class disqualified them from
    seeding (typically `tandem_goi_copy`). They appear in MAFFT/IQ-TREE so
    species like Apis cerana — which only carry tandem-duplicate evidence —
    are not silently dropped from the phylogeny.

    A3 (TODO_JUN, opt-in via ``seed_on_flanking_support``): the strict class
    filter "neutered" the wave — divergent/tandem orthologs (e.g. Apis cerana's
    100%-identical melittin, classed `tandem_goi_copy`) never seeded the next
    wave, so a close relative's perfect copy couldn't bridge to a more divergent
    bee. When enabled, a HIGH/MEDIUM GOI feature ALSO seeds — regardless of
    goi_class — on a generic synteny signal: strong flanking support
    (``block_flanking_support >= seed_min_flanking``) AND real query coverage
    (``query_coverage >= seed_min_qcov``). It keys on coverage+flanking, not
    identity, so it doesn't preferentially admit paralogs; §1m collinearity and
    §1m ownership remain the drift backstops.

    A context-promoted call (``_promote_context_supported_goi``) rests on gene order,
    not on sequence. It feeds the tree but never seeds the next wave, so the search
    itself runs as without the promotion and a wrong promotion cannot carry into
    later waves.

    Suppressed_count is the total of GOI-role genes that did NOT qualify for
    seeding (used for logging).
    """
    seed = []
    tree_extra = []
    suppressed = 0
    for g in all_genes:
        gid = g.get("id", "")
        meta = feature_meta.get(gid, {})
        role = meta.get("role", "goi" if is_goi_query_id_fn(gid) else "flanking")
        if role != "goi":
            continue
        confidence = meta.get("confidence", "")
        goi_class = meta.get("goi_class", "")
        promoted = bool(meta.get("context_promoted"))
        qualifies = (not promoted
                     and confidence in _TREE_QUALIFYING_CONFIDENCE
                     and goi_class in _SEED_QUALIFYING_CLASSES)
        if (not qualifies and not promoted and seed_on_flanking_support
                and confidence in _TREE_QUALIFYING_CONFIDENCE):
            fs = _safe_float(meta.get("block_flanking_support"))
            qc = _safe_float(meta.get("query_coverage"))
            if (fs is not None and qc is not None
                    and fs >= seed_min_flanking and qc >= seed_min_qcov):
                qualifies = True
        if qualifies:
            seed.append(g)
        else:
            suppressed += 1
            if confidence in _TREE_QUALIFYING_CONFIDENCE:
                tree_extra.append(g)
    return seed, tree_extra, suppressed


# ---------------------------------------------------------------------------
# §1o — dispersed-GOI rescue.
#
# Block seeding uses FLANKING hits only (see process_single_genome), then the GOI
# is searched inside each flanking-seeded block. That silently discards a strong
# GOI hit that dispersed away from its flanking anchors. Motivating case: Aedes
# oskar hit at 45.3% id / e-2.8e-29 sat 62.8 Mb from its flanking anchors — which
# are on the SAME chromosome (NC_035109.1) — so no block formed there and 0 GOI
# models were emitted, even though the search had already found the real oskar.
#
# The neighbourhood being SHATTERED is what licenses trusting a distant same-
# chromosome GOI hit: if the flanking order is destroyed, "the GOI must be next to
# its flanking" is no longer a valid constraint. In a collinear (well-conserved)
# genome an isolated GOI hit far from flanking is still treated as a paralog
# (strict), so close-relative benchmarks are untouched. These two pure helpers
# encode that decision; process_single_genome adds the qualifying hits to the
# block seed set so they get modelled by the normal machinery (as LOW fallback).
# ---------------------------------------------------------------------------

def _flanking_anchor_positions(flanking_hits):
    """Collapse the noisy raw hits to ONE canonical location per flanking gene
    (its highest-scoring hit), grouped by chromosome. The raw m8 has many
    partial/spurious flanking hits scattered genome-wide; the best hit per gene is
    its real anchor. Returns {chrom: [(gene, mid, start, end), ...]}."""
    best = {}
    for h in flanking_hits:
        q = h.get('query', '')
        if not q or is_goi_query_id(q):
            continue
        g = extract_base_gene_id(q)
        if g not in best or h.get('bits', 0) > best[g].get('bits', 0):
            best[g] = h
    by_chrom = defaultdict(list)
    for g, h in best.items():
        s, e = h.get('start', 0), h.get('end', 0)
        by_chrom[h.get('chrom')].append((g, (s + e) // 2, s, e))
    return by_chrom


def compute_rearrangement_score(flanking_hits, min_genes=3):
    """Diagnostic only: fraction of distinct flanking genes NOT on the dominant
    chromosome (cross-chromosome scatter, [0,1]). Logged for context; the rescue
    gate below is the per-chromosome flanking ENVELOPE, which validation showed to
    be the correct, self-scaling signal (a stretched neighbourhood on ONE
    chromosome — the real Aedes case — has LOW scatter but a WIDE envelope)."""
    by_chrom = _flanking_anchor_positions(flanking_hits)
    n = sum(len(v) for v in by_chrom.values())
    if n < min_genes:
        return 0.0
    best = max((len(v) for v in by_chrom.values()), default=0)
    return 1.0 - best / float(n)


def _nearest_anchor(chrom, hs, he, by_chrom):
    """(distance_bp, gene_id) to the closest flanking anchor anywhere in the genome.

    Distance is 0 when the hit overlaps an anchor. Same-chromosome anchors always win
    over any off-chromosome one, so the number answers "how far is this hit from the
    neighbourhood *on its own scaffold*" and falls back to the global nearest only when
    the chromosome carries no anchor at all.
    """
    best = (None, "", False)  # dist, gene, same_chrom
    for c, genes in by_chrom.items():
        same = (c == chrom)
        for gene, _m, s, e in genes:
            d = 0 if not (he < s or hs > e) else min(abs(s - he), abs(hs - e))
            if best[0] is None or (same, -d) > (best[2], -best[0]):
                best = (d, gene, same)
    return best[0], best[1], best[2]


NONSYNTENIC_COLUMNS = [
    'genome', 'query', 'chrom', 'start', 'end', 'identity', 'evalue', 'alnlen',
    'rejected_by', 'nearest_anchor_gene', 'nearest_anchor_bp',
    'nearest_anchor_same_chrom',
]


NOT_SEARCHED_COLUMNS = ['genome', 'reason', 'details']


def _write_not_searched(args, prefix, genome_name, reason, details=''):
    """Record a genome this search could not search, next to the raw hits (the report
    stages hits/ like the §1x candidates). Without it the genome produced no result file
    and vanished from the report: neither a call nor "absent" (bug hunt 2026-10-13, A10)."""
    path = os.path.join(args.output_dir, 'hits', f"{prefix}{genome_name}.not_searched.tsv")
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=NOT_SEARCHED_COLUMNS, delimiter='\t')
            w.writeheader()
            w.writerow({'genome': genome_name, 'reason': reason,
                        'details': ' '.join(str(details).split())[:500]})
    except OSError as exc:
        logger.warning(f"[{genome_name}] could not record the skipped search: {exc}")


def _write_nonsyntenic_candidates(args, prefix, genome_name, records):
    """§1x: write the refused-but-strong GOI candidates next to the raw hits.

    Written INTO `hits/`, which `modules/generate_report.nf` already stages wholesale
    (`cp -r "$d"/* staged_results/hits/`). That means the report picks these up with no
    new Nextflow channel and no edit to `main.nf`'s `workflow {}` — which matters,
    because that method body is at the JVM 64 kB limit (CLAUDE.md §17).
    """
    if not records:
        return
    path = os.path.join(args.output_dir, 'hits',
                        f"{prefix}{genome_name}.nonsyntenic.tsv")
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=NONSYNTENIC_COLUMNS, delimiter='\t')
            w.writeheader()
            for r in records:
                w.writerow({k: r.get(k, '') for k in NONSYNTENIC_COLUMNS})
        logger.info(f"[{genome_name}] §1x wrote {len(records)} non-syntenic GOI "
                    f"candidate(s) to {os.path.basename(path)}")
    except OSError as exc:
        logger.warning(f"[{genome_name}] §1x could not write {path}: {exc}")


def select_dispersed_goi_seeds(hits, flanking_hits, *,
                               min_identity=40.0, max_evalue=1e-10,
                               min_alnlen=50, min_chrom_anchors=3,
                               region_padding=50000, envelope_margin=200000,
                               rejected_out=None, genome_name=""):
    """Rescue a strong GOI hit that dispersed WITHIN the flanking ENVELOPE of its
    chromosome, so it can seed its own block.

    The envelope — [min .. max] anchor position on a chromosome carrying
    >= `min_chrom_anchors` distinct flanking genes — IS the data-driven,
    per-target 'how far to trust a GOI hit' distance the flanking spacing implies:
      * conserved target  -> tight envelope; a real GOI is right next to its
        flanking (within region_padding), so the normal path already models it and
        this is a NO-OP; a paralog Mb away is OUTSIDE the tight envelope -> rejected.
      * rearranged target -> wide envelope (the Aedes oskar case: 11 anchors span
        15-243 Mb on NC_035109.1); a dispersed-but-syntenic GOI inside it is
        recovered.
    A strong hit on a chromosome with too few anchors, or outside the envelope, is
    rejected as a probable paralog / fold-analogue (the biglycan-on-chrX and
    structdisc cases have no such envelope). Hits already within region_padding of
    an anchor are skipped (handled by the normal block path).

    Returns a list of GOI hit dicts to add to the block seed set.

    §1x — `rejected_out`: when a list is passed, every hit that CLEARS the quality bar
    (identity / e-value / alignment length) but is then refused by a SYNTENY gate is
    appended to it as a record. Refusing such a hit is correct — at long divergences a
    translocated ortholog is genuinely indistinguishable from a paralog by synteny alone
    — but refusing it *silently* is not: the run reports "no ortholog" while holding a
    54-72 % hit inside the real gene (the Anopheles defensin case). The caller writes
    these to `hits/<genome>.nonsyntenic.tsv` so the report can say "found, not syntenic"
    instead of "absent".
    """
    by_chrom = _flanking_anchor_positions(flanking_hits)

    def _reject(h, gate):
        if rejected_out is None:
            return
        hs, he = h.get('start', 0), h.get('end', 0)
        dist, anchor, same = _nearest_anchor(h.get('chrom'), hs, he, by_chrom)
        rejected_out.append({
            'genome': genome_name,
            'query': h.get('query', ''),
            'chrom': h.get('chrom', ''),
            'start': hs,
            'end': he,
            'identity': round(float(h.get('pident', 0.0) or 0.0), 2),
            'evalue': h.get('evalue', ''),
            'alnlen': h.get('alnlen', 0),
            'rejected_by': gate,
            'nearest_anchor_gene': anchor,
            'nearest_anchor_bp': dist if dist is not None else '',
            'nearest_anchor_same_chrom': bool(same),
        })
        logger.info(
            f"[{genome_name}] §1x non-syntenic GOI candidate REJECTED ({gate}): "
            f"{h.get('chrom')}:{hs}-{he} id={h.get('pident', 0):.1f}% "
            f"e={h.get('evalue')} aln={h.get('alnlen')} — nearest flanking anchor "
            + (f"{anchor} at {dist:,} bp"
               + ("" if same else " on a DIFFERENT chromosome") if dist is not None
               else "none in this genome"))

    def _clears_quality(h):
        return (is_goi_query_id(h.get('query', ''))
                and h.get('pident', 0.0) >= min_identity
                and h.get('evalue', 1.0) <= max_evalue
                and h.get('alnlen', 0) >= min_alnlen)

    envelopes = {}
    for c, genes in by_chrom.items():
        if len(genes) >= min_chrom_anchors:
            lo = min(g[2] for g in genes) - envelope_margin
            hi = max(g[3] for g in genes) + envelope_margin
            envelopes[c] = (lo, hi, genes)
    if not envelopes:
        # No chromosome carries enough anchors to define an envelope anywhere. Strong
        # GOI hits still deserve reporting — this is the "we could not judge synteny at
        # all" case, which is different again from "judged and refused".
        for h in hits:
            if _clears_quality(h):
                _reject(h, 'no_flanking_envelope_in_genome')
        return []

    seeds, seen = [], set()
    for h in hits:
        if not is_goi_query_id(h.get('query', '')):
            continue
        if h.get('pident', 0.0) < min_identity or h.get('evalue', 1.0) > max_evalue:
            continue
        if h.get('alnlen', 0) < min_alnlen:
            continue
        env = envelopes.get(h.get('chrom'))
        if not env:                              # chromosome doesn't carry the neighbourhood
            _reject(h, 'chromosome_lacks_min_anchors')
            continue
        lo, hi, genes = env
        hs, he = h.get('start', 0), h.get('end', 0)
        if he < lo or hs > hi:                   # outside the flanking envelope
            _reject(h, 'outside_flanking_envelope')
            continue
        near_anchor = any(
            not (he < s - region_padding or hs > e + region_padding)
            for _g, _m, s, e in genes)
        if near_anchor:                          # already covered by a flanking block
            continue
        key = (h.get('chrom'), hs // max(1, region_padding))
        if key in seen:
            continue
        seen.add(key)
        seeds.append(h)
    return seeds


def process_single_genome(genome_path, db_path, args, home_db_dir, prefix, threads_per_job):
    """
    Worker function to search a single genome.
    Returns: (genome_name, list_of_new_genes, list_of_tree_extra_genes,
             error_message_or_none)
    """
    genome_name = os.path.basename(genome_path)
    if not os.path.exists(genome_path):
        msg = "Genome file not found"
        logger.error(f"[{genome_name}] {msg}.")
        return genome_name, [], [], msg
    
    unique_id = uuid.uuid4().hex
    hits_file = f"{args.output_dir}/hits/{prefix}{genome_name}.m8"
    tmp_dir = f"{args.output_dir}/tmp_mmseqs_{unique_id}_{genome_name}"
    
    new_genes = []
    tree_extra_genes = []
    error_message = None

    try:
        # Parse DB once so we can adapt hit-length filtering for short proteins.
        base_db_sequences = {}
        query_lengths = {}
        for header, clean_id, seq in parse_fasta(db_path):
            base_db_sequences[clean_id] = {'id': clean_id, 'seq': seq, 'header': header}
            query_lengths[clean_id] = len(seq)
            base = extract_base_gene_id(clean_id)
            if base != clean_id and base not in base_db_sequences:
                base_db_sequences[base] = {'id': clean_id, 'seq': seq, 'header': header}
            query_lengths.setdefault(base, len(seq))

        # Identify flanking proteins that are proxies for the GOI query
        # (e.g., gene-LY6E is a proxy for GOI_LY6E — same gene, different prefix).
        # Proxy blocks seeded by these genes will bypass the min_block_genes filter.
        _goi_bare_names: set = set()
        for _rid in list(base_db_sequences.keys()):
            if is_goi_query_id(_rid) and '|' not in _rid:
                _bare = _rid
                if _bare.startswith('GOI_copy_'):
                    _bare = _bare[9:]
                elif _bare.startswith('GOI_'):
                    _bare = _bare[4:]
                _goi_bare_names.add(_bare)
        goi_proxy_flanking_parents: set = set()
        for _rid in list(base_db_sequences.keys()):
            if is_goi_query_id(_rid) or '|' in _rid:
                continue
            _bare = _rid[5:] if _rid.startswith('gene-') else _rid
            if _bare in _goi_bare_names:
                goi_proxy_flanking_parents.add(_rid)
        if goi_proxy_flanking_parents:
            logger.info(
                f"[{genome_name}] GOI proxy flanking parents: "
                f"{sorted(goi_proxy_flanking_parents)}"
            )

        # §1m collinearity-aware seed placement. home_rank is built once in main()
        # from the pristine initial_db and rides on `args`; disable_synteny_collinearity
        # falls fully back to the legacy gene-count clustering/scoring.
        _collinearity_on = not getattr(args, 'disable_synteny_collinearity', False)
        home_rank = getattr(args, 'home_rank', None) if _collinearity_on else None
        home_rank = home_rank or None  # treat empty dict as "off"
        _bridge_max_rank_gap = getattr(args, 'synteny_bridge_max_rank_gap', 5)
        _bridge_min_anchors = getattr(args, 'synteny_bridge_min_anchors', 2)
        _bridge_two_sided = bool(getattr(args, 'synteny_bridge_two_sided', False))
        _bridge_max_per_block = getattr(args, 'synteny_bridge_max_per_block', 2)

        # 1. Search (MMseqs) with low-memory retries for fragmented genomes.
        mmseqs_ok, mmseqs_details, resource_fail = run_mmseqs_easy_search_with_retries(
            db_path=db_path,
            genome_path=genome_path,
            hits_file=hits_file,
            tmp_dir=tmp_dir,
            args=args,
            threads_per_job=threads_per_job,
            genome_name=genome_name,
        )
        if not mmseqs_ok:
            if resource_fail:
                # Skipped, not failed: one fragmented genome must not end the whole search.
                # But it is not "searched, nothing found" either: the report lists it.
                logger.error(
                    f"[{genome_name}] NOT SEARCHED: MMseqs resource failures persisted across "
                    "retries (prefilter/search crash). The report lists this genome under "
                    "genomes_not_searched."
                )
                _write_not_searched(args, prefix, genome_name, 'mmseqs_resource_failure',
                                    mmseqs_details)
                return genome_name, [], [], None
            raise RuntimeError(
                f"MMseqs easy-search failed for {genome_name}: {mmseqs_details}"
            )

        hits = parse_hits(
            hits_file,
            args.min_identity,
            args.min_length,
            args.evalue,
            query_lengths=query_lengths,
        )
        if not hits:
            logger.info(f"[{genome_name}] No hits found in MMseqs output.")
            return genome_name, [], [], None
            
        logger.info(f"[{genome_name}] Parsed {len(hits)} hits.")
        if not getattr(args, 'disable_anchor_resolution', False):
            _n_before = len(hits)
            hits = resolve_anchor_hits(hits)
            if len(hits) < _n_before:
                logger.info(
                    f"[{genome_name}] One anchor per place: {_n_before - len(hits)} flanking hit(s) "
                    f"removed where a stronger hit of another flanking gene covers the same DNA.")

        # 2. Identify Synteny (Multiple Blocks)
        # Anchor blocks with non-GOI (flanking) hits first, then search GOI inside them.
        # This follows the intended flow: flanking -> region -> GOI modeling.
        flanking_seed_hits = [h for h in hits if not is_goi_query_id(h.get('query', ''))]
        seed_hits = flanking_seed_hits if flanking_seed_hits else hits
        if flanking_seed_hits:
            logger.info(f"[{genome_name}] Using {len(flanking_seed_hits)} flanking-anchor hits for block seeding.")
        else:
            logger.info(f"[{genome_name}] No flanking anchors found; falling back to GOI/mixed hits for block seeding.")

        # From here on `args` holds this genome's distances in base pairs.
        args = genome_scale_args(args, flanking_seed_hits, genome_path, genome_name)
        c_dist = args.cluster_distance
        _bridge_max_gap = args.synteny_bridge_max_gap
        _gene_gap = (lambda _g: _gene_limit(args, _g, args.max_intron))
        _max_places = int(getattr(args, 'anchor_max_places', 0) or 0)
        if _max_places > 0 and flanking_seed_hits:
            seed_hits = limit_anchor_places(flanking_seed_hits, _max_places, _gene_gap)
            logger.info(
                f"[{genome_name}] Anchor places: blocks are seeded from {len(seed_hits)} of "
                f"{len(flanking_seed_hits)} flanking hits (the {_max_places} best places of each gene).")

        # §1o: rescue a strong GOI hit that dispersed away from its flanking
        # anchors (same-chromosome macro-synteny) by adding it as a block seed, so
        # the normal machinery models it. No-op unless the neighbourhood is
        # shattered (rearrangement gate) — protects well-conserved benchmarks.
        rescue_goi_genes = set()
        # §1x: collect the strong-but-non-syntenic candidates the gate refuses, so the
        # run can distinguish "no ortholog here" from "found one, could not place it".
        _nonsyntenic = [] if getattr(args, 'report_nonsyntenic_candidates', True) else None
        if getattr(args, 'enable_dispersed_goi_rescue', True) and flanking_seed_hits:
            _rscore = compute_rearrangement_score(flanking_seed_hits)
            _disp_seeds = select_dispersed_goi_seeds(
                hits, flanking_seed_hits,
                min_identity=getattr(args, 'dispersed_goi_min_identity', 40.0),
                max_evalue=getattr(args, 'dispersed_goi_max_evalue', 1e-10),
                min_alnlen=getattr(args, 'dispersed_goi_min_alnlen', 50),
                min_chrom_anchors=getattr(args, 'dispersed_min_chrom_anchors', 3),
                region_padding=args.region_padding,
                envelope_margin=getattr(args, 'dispersed_envelope_margin', 200000),
                rejected_out=_nonsyntenic, genome_name=genome_name)
            _write_nonsyntenic_candidates(args, prefix, genome_name, _nonsyntenic)
            if _disp_seeds:
                seed_hits = seed_hits + _disp_seeds
                rescue_goi_genes = {extract_base_gene_id(h.get('query', ''))
                                    for h in _disp_seeds}
                logger.info(
                    f"[{genome_name}] §1o dispersed-GOI rescue: added "
                    f"{len(_disp_seeds)} strong GOI hit(s) inside the flanking "
                    f"envelope as rescue seed(s) (scatter={_rscore:.2f}, "
                    f"id>={getattr(args, 'dispersed_goi_min_identity', 40.0)}).")

        # §1g: estimate this target's phylogenetic distance from the median identity
        # of its flanking-gene matches and relax the GOI classify identity bars for
        # distant targets (no-op for close relatives). Must run before any GOI
        # classification for this genome (process_region_block below).
        _dist_info = _apply_distance_adaptive_thresholds(args, flanking_seed_hits, genome_name)

        synteny_blocks = identify_synteny_blocks(
            seed_hits,
            max_intron=args.max_intron,
            cluster_distance=c_dist,
            home_rank=home_rank,
            bridge_max_gap=_bridge_max_gap,
            bridge_max_rank_gap=_bridge_max_rank_gap,
            bridge_min_anchors=_bridge_min_anchors,
            bridge_two_sided=_bridge_two_sided,
            bridge_max_per_block=_bridge_max_per_block,
            gene_gap_limit=_gene_gap,
        )
        if home_rank:
            _n_bridged = sum(1 for _b in synteny_blocks if _b.get('bridged'))
            if _n_bridged:
                logger.info(
                    f"[{genome_name}] §1m collinearity: bridged a flanking gap in "
                    f"{_n_bridged} block(s) (home-order continuation across rearrangement)."
                )

        if not synteny_blocks:
            logger.info(f"[{genome_name}] No valid syntenic region found.")
            return genome_name, [], [], None
            
        # Optimization: Merge overlapping search regions
        pre_merge_count = len(synteny_blocks)
        synteny_blocks = merge_synteny_blocks(synteny_blocks, args.region_padding)
        if len(synteny_blocks) < pre_merge_count:
            logger.info(f"[{genome_name}] Merged {pre_merge_count} blocks into {len(synteny_blocks)} discrete search regions.")
        else:
            logger.info(f"[{genome_name}] Found {len(synteny_blocks)} discrete syntenic blocks.")

        # Keep only blocks with enough anchor genes to avoid spending hours
        # on singleton/noise loci in fragmented genomes.
        # Exception: single-gene blocks seeded by a GOI proxy (e.g., gene-LY6E)
        # are retained so the full GOI search runs on the correct scaffold.
        if args.min_block_genes > 1:
            pre_filter_count = len(synteny_blocks)
            filtered_blocks = []
            proxy_blocks = []
            for _b in synteny_blocks:
                if _b.get('genes_count', 0) >= args.min_block_genes:
                    filtered_blocks.append(_b)
                elif goi_proxy_flanking_parents and any(
                    _g in goi_proxy_flanking_parents for _g in _b.get('genes', [])
                ):
                    proxy_blocks.append(_b)
                elif any(is_goi_query_id(_g) or _g in rescue_goi_genes
                         for _g in _b.get('genes', [])):
                    # §1o: a GOI-seeded (dispersed-rescue) block is single-gene by
                    # construction — keep it despite min_block_genes.
                    proxy_blocks.append(_b)
            synteny_blocks = filtered_blocks + proxy_blocks
            logger.info(
                f"[{genome_name}] Block filter (min_block_genes={args.min_block_genes}): "
                f"{len(filtered_blocks)} standard + {len(proxy_blocks)} GOI-proxy "
                f"= {len(synteny_blocks)}/{pre_filter_count} retained."
            )

        # Fallback: if filtering removed everything, still keep top blocks by score.
        if not synteny_blocks:
            logger.info(
                f"[{genome_name}] No blocks after filtering. "
                f"Falling back to top {args.max_blocks_per_genome} unfiltered blocks."
            )
            synteny_blocks = identify_synteny_blocks(
                seed_hits,
                max_intron=args.max_intron,
                cluster_distance=c_dist,
                home_rank=home_rank,
                bridge_max_gap=_bridge_max_gap,
                bridge_max_rank_gap=_bridge_max_rank_gap,
                bridge_min_anchors=_bridge_min_anchors,
                bridge_two_sided=_bridge_two_sided,
                bridge_max_per_block=_bridge_max_per_block,
                gene_gap_limit=_gene_gap,
            )
            synteny_blocks = merge_synteny_blocks(synteny_blocks, args.region_padding)

        # §1m: when capping bites, keep the most collinear blocks (merge re-sorted
        # by position). GOI-proxy blocks (e.g. gene-LY6E) are protected first so a
        # low-collinearity proxy seed is never dropped by the cap.
        if home_rank and len(synteny_blocks) > 1:
            def _cap_key(b):
                genes = set(b.get('genes', []))
                is_proxy = (bool(genes & goi_proxy_flanking_parents)
                            or bool(genes & rescue_goi_genes)
                            or any(is_goi_query_id(g) for g in genes))
                return (1 if is_proxy else 0,
                        b.get('collinear_chain_len', 0),
                        b.get('genes_count', 0))
            synteny_blocks.sort(key=_cap_key, reverse=True)

        # Hard cap per genome to prevent pathological runtimes.
        if args.max_blocks_per_genome > 0 and len(synteny_blocks) > args.max_blocks_per_genome:
            logger.info(
                f"[{genome_name}] Capping blocks: {len(synteny_blocks)} -> "
                f"{args.max_blocks_per_genome} (max_blocks_per_genome)."
            )
            synteny_blocks = synteny_blocks[:args.max_blocks_per_genome]
        
        all_genes = []
        all_gff_lines = []
        
        genome_seqs = load_genome(genome_path)
        native_gff_path = find_native_annotation_path(genome_path)
        native_annot_index = load_native_annotation_index(native_gff_path)
        if native_gff_path and native_annot_index:
            logger.info(
                f"[{genome_name}] Native annotation index loaded from "
                f"{os.path.basename(native_gff_path)} ({len(native_annot_index)} contigs)."
            )
        elif native_gff_path:
            logger.info(
                f"[{genome_name}] Native annotation file found ({os.path.basename(native_gff_path)}) "
                f"but no usable gene/transcript features parsed."
            )
        else:
            logger.info(f"[{genome_name}] No native annotation GFF found next to genome FASTA.")

        def _target_mrna_attrs(
            chrom_name: str,
            g_start: int,
            g_end: int,
            g_strand: str,
            base_attrs: Dict[str, Any],
        ) -> str:
            native = lookup_native_annotation(
                native_annot_index or {},
                chrom=chrom_name,
                start=g_start,
                end=g_end,
                strand=g_strand,
            )
            return _compose_gff_attrs(base_attrs, native)

        for i, block in enumerate(synteny_blocks):
            block_genes, block_gff = process_region_block(
                i,
                block,
                hits,
                genome_seqs,
                base_db_sequences,
                genome_name,
                args,
                unique_id,
                threads_per_job,
                native_annot_index=native_annot_index,
            )
            all_genes.extend(block_genes)
            all_gff_lines.extend(block_gff)

            if (i + 1) % 10 == 0:
                logger.info(f"[{genome_name}] Block progress: {i + 1}/{len(synteny_blocks)}.")

        # Flanking-only post-pass: keep one best chain per flanking parent per locus across blocks.
        all_genes, all_gff_lines = deduplicate_flanking_models(
            all_genes,
            all_gff_lines,
            genome_name=genome_name,
            locus_gap_bp=max(5000, int(args.cluster_distance)),
        )
        # Flanking models keep one CDS row per exon. They used to be collapsed to a
        # single gene-span CDS "for plotting", which threw the exon positions away
        # while the mRNA kept Exons=N -- so the plot drew N evenly spaced exons that
        # do not exist. (Rearranged flanking models, written below, always kept theirs.)

        # --- Cross-chromosome flanking recovery ---
        # Some flanking genes may have translocated to a different chromosome.
        # After standard block processing, check which flanking parents are still
        # missing, look for strong hits on other chromosomes, and annotate them
        # with a "rearranged" flag so the plot can show them.
        all_db_flanking_parents = {
            extract_base_gene_id(rid)
            for rid in base_db_sequences
            if not is_goi_query_id(rid) and '|' not in rid
        }
        found_flanking_parents = set()
        for g in all_genes:
            gid = g.get('id', '')
            if not is_goi_query_id(gid):
                found_flanking_parents.add(extract_base_gene_id(gid))

        missing_flanking = all_db_flanking_parents - found_flanking_parents
        block_chroms = {b['chrom'] for b in synteny_blocks}

        if missing_flanking and hits:
            # Gather off-block-chrom hits for missing parents
            off_chrom_hits_by_parent = defaultdict(list)
            for h in hits:
                parent = extract_base_gene_id(h.get('query', ''))
                if parent not in missing_flanking:
                    continue
                if h['chrom'] in block_chroms:
                    continue  # Already considered during block processing
                off_chrom_hits_by_parent[parent].append(h)

            if off_chrom_hits_by_parent:
                logger.info(
                    f"[{genome_name}] Cross-chromosome flanking recovery: "
                    f"{len(off_chrom_hits_by_parent)} missing parents have off-block hits."
                )
                flanking_query_map = build_flanking_query_by_parent(
                    base_db_sequences, set(off_chrom_hits_by_parent.keys())
                )
                clean_gname = genome_name.replace('.', '_').replace('-', '_').replace(' ', '_')
                flanking_locus_gap = max(5000, min(100000, args.max_intron * 2))

                for parent_id, parent_hits in off_chrom_hits_by_parent.items():
                    query_rec = flanking_query_map.get(parent_id)
                    if not query_rec:
                        continue
                    parent_query_seq = query_rec.get('seq', '')
                    if not parent_query_seq:
                        continue

                    # Remap exon-level qstart/qend
                    exon_offsets = query_rec.get('exon_offsets', {})
                    if exon_offsets:
                        for h in parent_hits:
                            m = re.search(r'\|exon_(\d+)$', h.get('query', ''))
                            if m:
                                exon_num = int(m.group(1))
                                if exon_num in exon_offsets:
                                    offset, _ = exon_offsets[exon_num]
                                    h['qstart'] = h.get('qstart', 1) + offset
                                    h['qend'] = h.get('qend', 1) + offset

                    parent_loci = split_hits_into_loci(
                        parent_hits, max_gap=_gene_limit(args, parent_id, flanking_locus_gap))
                    if not parent_loci:
                        continue

                    # Only process the single strongest off-chrom locus
                    best_locus = max(
                        parent_loci,
                        key=lambda hs: (max(h.get('bits', 0) for h in hs), len(hs))
                    )
                    best_bits = max(h.get('bits', 0) for h in best_locus)
                    if best_bits < 40:
                        continue  # Too weak to be a real ortholog

                    off_chrom = best_locus[0].get('chrom', '')
                    if off_chrom not in genome_seqs:
                        continue

                    # Build a local region around the off-chrom hits
                    off_start = min(h.get('start', 0) for h in best_locus)
                    off_end = max(h.get('end', 0) for h in best_locus)

                    # If this off-chrom hit is a GOI proxy (e.g., gene-LY6E),
                    # run a full GOI block search instead of just annotating it
                    # as rearranged_flanking.
                    if parent_id in goi_proxy_flanking_parents:
                        _proxy_block = {
                            'chrom': off_chrom,
                            'start': off_start,
                            'end': off_end,
                            'genes_count': 1,
                            'loci_count': len(best_locus),
                            'genes': [parent_id],
                        }
                        logger.info(
                            f"[{genome_name}] Cross-chrom GOI proxy {parent_id} on "
                            f"{off_chrom}:{off_start}-{off_end} "
                            f"({best_bits:.0f} bits) — running full GOI search."
                        )
                        _proxy_idx = len(synteny_blocks) + 1000 + len(all_genes)
                        _p_genes, _p_gff = process_region_block(
                            _proxy_idx,
                            _proxy_block,
                            hits,
                            genome_seqs,
                            base_db_sequences,
                            genome_name,
                            args,
                            unique_id,
                            threads_per_job,
                            native_annot_index=native_annot_index,
                        )
                        all_genes.extend(_p_genes)
                        all_gff_lines.extend(_p_gff)
                        found_flanking_parents.add(parent_id)
                        continue

                    off_slen = len(genome_seqs[off_chrom])
                    off_pad = min(args.region_padding, 50000)
                    off_w_start = max(0, off_start - off_pad)
                    off_w_end = min(off_slen, off_end + off_pad)
                    off_subseq = genome_seqs[off_chrom][off_w_start:off_w_end]

                    # Remap hits to local coordinates for annotation
                    local_hits = []
                    for h in best_locus:
                        local_start = max(0, int(h.get('start', 0)) - off_w_start)
                        local_end = min(len(off_subseq), int(h.get('end', 0)) - off_w_start)
                        if local_end <= local_start:
                            continue
                        local_hits.append({
                            'query': h.get('query', ''),
                            'qstart': h.get('qstart', 1),
                            'qend': h.get('qend', 1),
                            'gstart': local_start,
                            'gend': local_end,
                            'evalue': h.get('evalue', 1),
                            'pident': h.get('pident', 0),
                            'alnlen': h.get('alnlen', 0),
                            'bits': h.get('bits', 0),
                            'strand': h.get('strand', '+'),
                            'chrom': off_chrom,
                        })

                    if not local_hits:
                        continue

                    # Attempt miniprot-style exon annotation
                    exons = []
                    model_protein = None
                    try:
                        with maybe_quiet_streams(args.quiet_subtools):
                            exons, model_protein = annotate_exons_from_hit_list(
                                local_hits,
                                parent_query_seq,
                                off_subseq,
                                off_chrom,
                                gap_search_window=max(args.gap_search_window,
                                                      _gene_limit(args, parent_id, 0)),
                                sensitive=False
                            )
                    except Exception as e:
                        logger.warning(f"[{genome_name}] Off-block exon annotation failed at {off_chrom}:{off_w_start}-{off_w_end}: {e}")
                        exons = []

                    if exons:
                        exons.sort(key=lambda e: e.get('qstart', 0))
                        # miniprot's in-frame model translation (see annotate_using_miniprot).
                        exon_protein = model_protein or ''.join(e['seq'] for e in exons)
                        strand = exons[0].get('strand', '+')
                        avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
                        model_qcov = None
                        if parent_query_seq:
                            qmin_model = min(e.get('qstart', 1) for e in exons)
                            qmax_model = max(e.get('qend', 1) for e in exons)
                            model_qcov = (qmax_model - qmin_model + 1) / max(1, len(parent_query_seq))
                        global_start = off_w_start + min(e['gstart'] for e in exons) + 1
                        global_end = off_w_start + max(e['gend'] for e in exons)
                        new_id = f"{parent_id}|{clean_gname}_rearranged"

                        rearr_gff = [
                            f"{off_chrom}\trearranged_flanking\tmRNA\t{global_start}\t{global_end}\t"
                            f"{avg_pident:.1f}\t{strand}\t.\t"
                            f"{_target_mrna_attrs(off_chrom, global_start, global_end, strand, _flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'rearranged_flanking', 'Rearranged_from': ','.join(sorted(block_chroms))}, evidence_type='rearranged_flanking', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, context='cross_chromosome_rearranged'))}"
                        ]
                        for eidx, e in enumerate(exons, 1):
                            exon_gs = off_w_start + e['gstart'] + 1
                            exon_ge = off_w_start + e['gend']
                            exon_phase = int(e.get('phase', 0) or 0)
                            rearr_gff.append(
                                f"{off_chrom}\trearranged_flanking\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t{exon_phase}\t"
                                f"ID={new_id}_CDS{eidx};Parent={new_id}"
                            )

                        all_genes.append({
                            'id': new_id,
                            'seq': exon_protein,
                            'description': (
                                f"coords:{off_chrom}:{global_start}-{global_end} "
                                f"parent:{parent_id} type:rearranged_flanking "
                                f"identity:{avg_pident:.1f}"
                            )
                        })
                        all_gff_lines.extend(rearr_gff)
                        found_flanking_parents.add(parent_id)
                        logger.info(
                            f"[{genome_name}] Recovered rearranged flanking gene "
                            f"{parent_id} on {off_chrom} ({avg_pident:.1f}% identity, "
                            f"{len(exons)} exons)."
                        )
                    else:
                        # Fallback: use concatenated hit-derived protein
                        strand_votes = {'+': 0, '-': 0}
                        for h in local_hits:
                            strand_votes[h.get('strand', '+')] += 1
                        strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

                        ordered_hits = sorted(
                            [h for h in local_hits if h.get('strand', '+') == strand],
                            key=lambda h: h.get('gstart', 0)
                        )
                        ordered_hits = _longest_monotonic_query_chain(
                            ordered_hits, strand,
                            max_gap=_gene_limit(args, parent_id, args.max_intron)) if ordered_hits else []

                        if ordered_hits:
                            coding_frags = []
                            cds_intervals = []
                            for h in ordered_hits:
                                hs = h.get('gstart', 0)
                                he = h.get('gend', 0)
                                if he <= hs:
                                    continue
                                dna = off_subseq[hs:he]
                                if strand == '-':
                                    dna = reverse_complement(dna)
                                dna = dna[:len(dna) - (len(dna) % 3)]
                                if len(dna) < 9:
                                    continue
                                aa = translate(dna).replace('*', '')
                                if not aa:
                                    continue
                                coding_frags.append(aa)
                                cds_intervals.append((hs, he))

                            if coding_frags and cds_intervals:
                                if strand == '-':
                                    coding_frags = list(reversed(coding_frags))
                                flank_protein = ''.join(coding_frags)
                                global_start = off_w_start + min(s for s, _ in cds_intervals) + 1
                                global_end = off_w_start + max(e for _, e in cds_intervals)
                                avg_pident = _pooled_identity(ordered_hits)
                                qmin = min(min(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                                qmax = max(max(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                                model_qcov = ((qmax - qmin + 1) / len(parent_query_seq)) if parent_query_seq else None
                                new_id = f"{parent_id}|{clean_gname}_rearranged_fallback"

                                rearr_gff = [
                                    f"{off_chrom}\trearranged_flanking\tmRNA\t{global_start}\t{global_end}\t"
                                    f"{avg_pident:.1f}\t{strand}\t.\t"
                                    f"{_target_mrna_attrs(off_chrom, global_start, global_end, strand, _flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'rearranged_flanking_fallback', 'Rearranged_from': ','.join(sorted(block_chroms))}, evidence_type='rearranged_flanking_fallback', identity=avg_pident, exon_count=len(cds_intervals), query_cov=model_qcov, context='cross_chromosome_rearranged'))}"
                                ]
                                for cidx, (hs, he) in enumerate(cds_intervals, 1):
                                    exon_gs = off_w_start + hs + 1
                                    exon_ge = off_w_start + he
                                    rearr_gff.append(
                                        f"{off_chrom}\trearranged_flanking\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\t"
                                        f"ID={new_id}_CDS{cidx};Parent={new_id}"
                                    )

                                all_genes.append({
                                    'id': new_id,
                                    'seq': flank_protein,
                                    'description': (
                                        f"coords:{off_chrom}:{global_start}-{global_end} "
                                        f"parent:{parent_id} type:rearranged_flanking_fallback "
                                        f"identity:{avg_pident:.1f}"
                                    )
                                })
                                all_gff_lines.extend(rearr_gff)
                                found_flanking_parents.add(parent_id)
                                logger.info(
                                    f"[{genome_name}] Recovered rearranged flanking gene "
                                    f"{parent_id} on {off_chrom} (fallback, {avg_pident:.1f}% identity)."
                                )

        # --- GOI proxy sweep on block chromosomes (Step 3b) ---
        # Some GOI proxies (e.g., gene-LY6E) may have hits on block chromosomes
        # that fall OUTSIDE the padded regions already processed.  This happens
        # when the proxy hit is on the same large scaffold as another flanking
        # gene (e.g., CYP11B on the same NW_ contig) but far from that anchor.
        # We sweep all such loci and run a full GOI block search on each one.
        if goi_proxy_flanking_parents:
            _proxy_gap = max(5000, min(100000, args.max_intron * 2))
            # Build rough coverage set from already-processed blocks
            _covered: list = []
            for _b in synteny_blocks:
                _pad = args.region_padding
                _covered.append((_b['chrom'], _b['start'] - _pad, _b['end'] + _pad))

            # sorted(): each sweep marks its locus covered, so iteration order decides
            # which proxy's search runs where two loci overlap (and shifts _sw_idx).
            for _proxy_id in sorted(goi_proxy_flanking_parents):
                if _proxy_id in found_flanking_parents:
                    continue
                _proxy_on_block = [
                    h for h in hits
                    if extract_base_gene_id(h.get('query', '')) == _proxy_id
                    and h['chrom'] in block_chroms
                ]
                if not _proxy_on_block:
                    continue
                _proxy_loci = split_hits_into_loci(
                    _proxy_on_block, max_gap=_gene_limit(args, _proxy_id, _proxy_gap))
                for _lh in _proxy_loci:
                    _lbits = max(h.get('bits', 0) for h in _lh)
                    if _lbits < 40:
                        continue
                    _lchrom = _lh[0]['chrom']
                    _lstart = min(h.get('start', 0) for h in _lh)
                    _lend = max(h.get('end', 0) for h in _lh)
                    # Skip if already covered by a processed block
                    if any(
                        _cc == _lchrom and _lstart < _ce and _lend > _cs
                        for _cc, _cs, _ce in _covered
                    ):
                        continue
                    _sweep_block = {
                        'chrom': _lchrom,
                        'start': _lstart,
                        'end': _lend,
                        'genes_count': 1,
                        'loci_count': len(_lh),
                        'genes': [_proxy_id],
                    }
                    logger.info(
                        f"[{genome_name}] On-block GOI proxy sweep: {_proxy_id} on "
                        f"{_lchrom}:{_lstart}-{_lend} ({_lbits:.0f} bits) outside "
                        f"processed regions — running full GOI search."
                    )
                    _sw_idx = len(synteny_blocks) + 2000 + len(all_genes)
                    _sw_genes, _sw_gff = process_region_block(
                        _sw_idx,
                        _sweep_block,
                        hits,
                        genome_seqs,
                        base_db_sequences,
                        genome_name,
                        args,
                        unique_id,
                        threads_per_job,
                        native_annot_index=native_annot_index,
                    )
                    all_genes.extend(_sw_genes)
                    all_gff_lines.extend(_sw_gff)
                    if any(is_goi_query_id(_g.get('id', '')) for _g in _sw_genes):
                        found_flanking_parents.add(_proxy_id)
                    # Add to covered so we don't re-run overlapping loci
                    _covered.append((_lchrom, _lstart - args.region_padding,
                                     _lend + args.region_padding))

        # --- Flanking gene tracking summary (Step 4) ---
        still_missing = all_db_flanking_parents - found_flanking_parents
        if all_db_flanking_parents:
            found_count = len(found_flanking_parents & all_db_flanking_parents)
            total_count = len(all_db_flanking_parents)
            logger.info(
                f"[{genome_name}] Flanking gene summary: "
                f"{found_count}/{total_count} parents found."
            )
            if still_missing:
                # For each missing parent, report best hit info from raw hits
                missing_details = []
                for parent_id in sorted(still_missing):
                    parent_hits = [
                        h for h in hits
                        if extract_base_gene_id(h.get('query', '')) == parent_id
                    ]
                    if parent_hits:
                        best = max(parent_hits, key=lambda h: h.get('bits', 0))
                        missing_details.append(
                            f"{parent_id}(best_hit:{best['chrom']},"
                            f"bits={best.get('bits', 0):.1f},"
                            f"id={best.get('pident', 0):.1f}%)"
                        )
                    else:
                        missing_details.append(f"{parent_id}(no_hits)")
                logger.info(
                    f"[{genome_name}] Missing flanking parents: "
                    f"{'; '.join(missing_details)}"
                )

        feature_meta = {}

        all_gff_lines, n_promoted = _promote_context_supported_goi(all_gff_lines)
        if n_promoted:
            logger.info(
                f"[{genome_name}] Context promotion: {n_promoted} GOI model(s) LOW/AMBIGUOUS -> "
                f"MEDIUM on >= {CLASSIFY_THRESHOLDS['context_min_collinear']} anchors in home order."
            )

        # Write aggregated GFF
        if all_gff_lines:
             gff_out = f"{args.output_dir}/regions/{genome_name}.gff"
             faa_out = f"{args.output_dir}/regions/{genome_name}.faa"
             tsv_out = f"{args.output_dir}/regions/{genome_name}.homology.tsv"
             
             with open(gff_out, 'w') as gf:
                 gf.write("##gff-version 3\n")
                 if getattr(args, 'distances', None):
                     gf.write(synvoy_scale.format_gff_header(args.distances) + "\n")
                 for gl in all_gff_lines:
                     gf.write(gl + "\n")

             write_fasta([(g['id'], g['seq']) for g in all_genes], faa_out)

             for gl in all_gff_lines:
                 parts = gl.split("\t")
                 if len(parts) < 9 or parts[2] not in {"mRNA", "gene"}:
                     continue
                 attrs = _parse_gff_attributes(parts[8])
                 model_id = attrs.get("ID")
                 if not model_id:
                     continue
                 feature_meta[model_id] = {
                     "parent": _select_parent_id(attrs, model_id),
                     "role": attrs.get("SynVoyRole", "goi" if is_goi_query_id(model_id) else "flanking"),
                     "confidence": attrs.get("Confidence", ""),
                     "goi_class": attrs.get("GOIClass", ""),
                     "model_status": attrs.get("ModelStatus", ""),
                     "evidence_type": attrs.get("EvidenceType", attrs.get("Type", "")),
                     "identity": attrs.get("Identity", ""),
                     "n_exons": attrs.get("Exons", ""),
                     "synteny_context": attrs.get("SyntenyContext", ""),
                     "block_flanking_support": attrs.get("BlockFlankingSupport", ""),
                     "query_coverage": attrs.get("QueryCoverage", ""),
                     "target_gene": attrs.get("TargetGene", ""),
                     "target_product": attrs.get("TargetProduct", ""),
                     "embedding_similarity": attrs.get("EmbeddingSimilarity", ""),
                     "structural_similarity": attrs.get("StructuralSimilarity", ""),
                     "context_promoted": attrs.get("ContextPromoted", ""),
                 }

             with open(tsv_out, 'w') as tf:
                 tf.write(
                     "target_id\thome_id\trole\tconfidence\tgoi_class\tmodel_status\t"
                     "evidence_type\tidentity\tn_exons\tsynteny_context\t"
                     "block_flanking_support\tquery_coverage\ttarget_gene\ttarget_product\t"
                     "embedding_similarity\tstructural_similarity\n"
                 )
                 for rec in all_genes:
                     meta = feature_meta.get(rec['id'], {})
                     parent = meta.get("parent") or extract_base_gene_id(rec['id'])
                     tf.write(
                         "\t".join([
                             rec['id'],
                             parent,
                             meta.get("role", "goi" if is_goi_query_id(rec['id']) else "flanking"),
                             meta.get("confidence", ""),
                             meta.get("goi_class", ""),
                             meta.get("model_status", ""),
                             meta.get("evidence_type", ""),
                             meta.get("identity", ""),
                             meta.get("n_exons", ""),
                             meta.get("synteny_context", ""),
                             meta.get("block_flanking_support", ""),
                             meta.get("query_coverage", ""),
                             meta.get("target_gene", ""),
                             meta.get("target_product", ""),
                             meta.get("embedding_similarity", ""),
                             meta.get("structural_similarity", ""),
                         ]) + "\n"
                     )
        
        # Expansion DB should only be augmented with GOI-derived models that
        # survived confidence/ambiguity triage. Low-confidence ambiguous/tandem
        # calls are still reported in GFF/plots but should not recursively seed
        # later waves.
        new_genes, tree_extra_genes, suppressed_seed_count = (
            _classify_goi_for_seed_and_tree(
                all_genes, feature_meta, is_goi_query_id,
                seed_on_flanking_support=getattr(args, 'seed_on_flanking_support', False),
                seed_min_flanking=getattr(args, 'seed_flanking_min_count', 2),
                seed_min_qcov=getattr(args, 'seed_flanking_min_qcov', 0.5),
            )
        )
        if all_genes:
            logger.info(
                f"[{genome_name}] Expansion payload: {len(new_genes)} GOI-derived / "
                f"{len(all_genes)} total annotations"
                + (
                    f" ({suppressed_seed_count} GOI-like annotations withheld from seeding"
                    + (f", {len(tree_extra_genes)} of which still feed the tree"
                       if tree_extra_genes else "")
                    + ")."
                    if suppressed_seed_count
                    else "."
                )
            )

    except Exception as e:
        error_message = str(e)
        logger.error(f"[{genome_name}] Error processing: {error_message}")
    finally:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)

    return genome_name, new_genes, tree_extra_genes, error_message


# ---------------------------------------------------------------------------
# Main-process GPU augmentation phase (PLM + Foldseek), run AFTER the wave loop.
#
# GPU inference can't live in the forked workers (CUDA-in-fork), so it runs here
# once, in a single CUDA context, over the GOI candidate models the workers
# already wrote to regions/<genome>.gff + regions/<genome>.faa. gpu_augment
# loads ProtT5 then ESMFold sequentially (freeing each) and folds with
# domain-preserving windows, then Foldseek-compares — producing an
# embedding_similarity / structural_similarity per candidate, which we feed back
# through _classify_goi_evidence to boost/rescue confidence in the GFF.
# ---------------------------------------------------------------------------

def _select_goi_candidates_from_region(gff_path, faa_path):
    """Return [(mrna_id, protein_seq)] for GOI candidate models in a region GFF."""
    id_to_seq = {clean_id: seq for _h, clean_id, seq in parse_fasta(faa_path)}
    out = []
    seen = set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or '\t' not in line:
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue
            attrs = _parse_gff_attributes(parts[8])
            mid = attrs.get('ID')
            if not mid or mid in seen:
                continue
            is_goi = (attrs.get('SynVoyRole') == 'goi'
                      or is_goi_query_id(mid)
                      or bool(attrs.get('GOIClass')))
            if not is_goi:
                continue
            seq = id_to_seq.get(mid)
            if seq and len(seq) >= 10:
                out.append((mid, seq))
                seen.add(mid)
    return out


def _apply_ml_scores_to_mrna_line(mRNA_line, embedding_similarity, structural_similarity):
    """Attach ML scores to a candidate mRNA GFF line and re-classify confidence.

    Mirrors the legacy in-worker re-rank (a strict superset of it): ML signals
    only ever *raise* confidence (LOW/MEDIUM → higher). Idempotent — will not
    double-append EmbeddingSimilarity/StructuralSimilarity. Returns
    (new_line, boosted:int)."""
    cols = mRNA_line.split('\t')
    if len(cols) < 9:
        return mRNA_line, 0
    attrs_str = cols[8]
    attr_kv = {}
    for pair in attrs_str.split(';'):
        if '=' in pair:
            k, v = pair.split('=', 1)
            attr_kv[k.strip()] = v.strip()

    old_conf = attr_kv.get('Confidence', '')
    ev_type = attr_kv.get('EvidenceType', '')
    try:
        ident = float(attr_kv.get('Identity', 0))
    except (ValueError, TypeError):
        ident = 0.0
    try:
        exons = int(attr_kv.get('Exons', 1))
    except (ValueError, TypeError):
        exons = 1
    try:
        qcov = float(attr_kv.get('QueryCoverage', 0))
    except (ValueError, TypeError):
        qcov = None
    try:
        fl_sup = int(attr_kv.get('BlockFlankingSupport', 0))
    except (ValueError, TypeError):
        fl_sup = 0

    new_conf, new_class, new_reason = _classify_goi_evidence(
        evidence_type=ev_type,
        identity=ident,
        exon_count=exons,
        query_cov=qcov,
        flanking_support=fl_sup,
        embedding_similarity=embedding_similarity,
        structural_similarity=structural_similarity,
    )

    updated = attrs_str
    if embedding_similarity is not None and 'EmbeddingSimilarity=' not in attrs_str:
        updated += f";EmbeddingSimilarity={embedding_similarity:.3f}"
    if structural_similarity is not None and 'StructuralSimilarity=' not in attrs_str:
        updated += f";StructuralSimilarity={structural_similarity:.3f}"

    boosted = 0
    # A context-promoted MEDIUM was not reached by this classifier; re-running it here
    # without the block's gene order must not take it back down.
    demotes_promoted = 'ContextPromoted' in attr_kv and new_conf != 'HIGH'
    if new_conf != old_conf and old_conf in ('LOW', 'MEDIUM') and not demotes_promoted:
        updated = updated.replace(f"Confidence={old_conf}", f"Confidence={new_conf}")
        old_reason = attr_kv.get('InferenceReason', '')
        if old_reason:
            updated = updated.replace(
                f"InferenceReason={old_reason}", f"InferenceReason={new_reason}")
        old_class = attr_kv.get('GOIClass', '')
        if old_class:
            updated = updated.replace(f"GOIClass={old_class}", f"GOIClass={new_class}")
        boosted = 1

    cols[8] = updated
    return '\t'.join(cols), boosted


def _apply_scores_to_region_gff(gff_path, scores):
    """Rewrite a region GFF, attaching ML scores + re-classifying. Returns
    (n_attached, n_boosted)."""
    with open(gff_path) as fh:
        lines = fh.readlines()
    out = []
    attached = 0
    boosted = 0
    for line in lines:
        if line.startswith('#') or '\t' not in line:
            out.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'mRNA':
            out.append(line)
            continue
        attrs = _parse_gff_attributes(parts[8])
        mid = attrs.get('ID')
        sc = scores.get(mid) if mid else None
        if not sc:
            out.append(line)
            continue
        new_line, did_boost = _apply_ml_scores_to_mrna_line(
            line.rstrip('\n'),
            sc.get('embedding_similarity'),
            sc.get('structural_similarity'),
        )
        out.append(new_line + '\n')
        attached += 1
        boosted += did_boost
    if attached:
        with open(gff_path, 'w') as fh:
            fh.writelines(out)
    return attached, boosted


def run_gpu_augment_phase(args, genome_map=None):
    """Score GOI candidate models with ProtT5 + ESMFold/Foldseek in the main
    process (no fork), and fold the scores back into the region GFFs.

    Opt-in via --enable_plm_search / --enable_structural_search. Fail-loud by
    default (--gpu_augment_fail_loud) so a requested-but-broken GPU layer errors
    rather than silently producing no ML evidence (the old fork-CUDA no-op).

    With --enable_structural_discovery, also predicts ORFs in each syntenic
    window (from the flanking genes), folds/embeds the promising ones, and
    materialises any that structurally match the GOI as new gene models —
    catching orthologs sequence search missed entirely. `genome_map` maps a
    region-file stem (basename of the genome FASTA) to that FASTA path; required
    for discovery (to read the target sequence), unused for re-rank."""
    if not (getattr(args, 'enable_plm_search', False)
            or getattr(args, 'enable_structural_search', False)):
        return

    fail_loud = getattr(args, 'gpu_augment_fail_loud', True)
    if not GPU_AUGMENT_OK:
        msg = ("[gpu_augment] GPU augmentation requested but the gpu_augment "
               "module / its deps are unavailable.")
        if fail_loud:
            raise RuntimeError(msg)
        logger.warning("%s Continuing without ML scores.", msg)
        return

    regions_dir = os.path.join(args.output_dir, "regions")
    if not os.path.isdir(regions_dir):
        logger.info("[gpu_augment] no regions/ dir; nothing to score.")
        return

    items = []
    genome_gffs = {}
    for fn in sorted(os.listdir(regions_dir)):
        if not fn.endswith(".gff"):
            continue
        genome = fn[:-4]
        gff_path = os.path.join(regions_dir, fn)
        faa_path = os.path.join(regions_dir, genome + ".faa")
        if not os.path.exists(faa_path):
            continue
        genome_gffs[genome] = gff_path
        for cid, seq in _select_goi_candidates_from_region(gff_path, faa_path):
            items.append(gpu_augment.WorkItem("candidate", cid, seq, genome))

    # ── Discovery: predict ORFs in each syntenic window and queue them so the
    #    GPU phase can find orthologs sequence search missed. Non-fatal per
    #    genome (additive; a discovery bug must not sink the validated re-rank).
    orf_records = {}   # wid -> (genome, gff_path, orf dict)
    genome_synteny = {}  # genome -> {"flanking": [...], "seq_goi_chroms": set()}
    if getattr(args, 'enable_structural_discovery', False) and genome_map:
        try:
            import discover_regions as dr
        except Exception as _dr_err:
            dr = None
            logger.warning("[gpu_augment] discovery module unavailable: %s", _dr_err)
        if dr is not None:
            for genome, gff_path in sorted(genome_gffs.items()):
                gfasta = genome_map.get(genome)
                if not gfasta or not os.path.exists(gfasta):
                    continue
                try:
                    flanking, existing = dr.read_region_features(gff_path)
                    if not flanking:
                        continue
                    genome_synteny[genome] = {
                        "flanking": flanking,
                        "seq_goi_chroms": dr.seq_goi_chroms(gff_path, min_confidence="MEDIUM"),
                    }
                    genome_seqs = load_genome(gfasta)
                    chrom_lengths = {c: len(s) for c, s in genome_seqs.items()}
                    windows = dr.cluster_intervals(
                        flanking,
                        max_gap=getattr(args, 'discovery_window_gap', 100000),
                        pad=getattr(args, 'discovery_window_pad', 20000),
                        chrom_lengths=chrom_lengths)
                    disc_work = os.path.join(args.output_dir, "gpu_augment_work",
                                             "discovery", genome)
                    for wi, window in enumerate(windows):
                        orfs = dr.predict_region_orfs(
                            genome_seqs, window,
                            os.path.join(disc_work, f"w{wi}"),
                            predictor='augustus',
                            augustus_species=getattr(args, 'augustus_species', 'fly'),
                            min_aa=getattr(args, 'discovery_min_aa', 30))
                        for orf in orfs:
                            if not orf.get('seq') or len(orf['seq']) < 20:
                                continue
                            if dr.orf_overlaps_any(orf, existing, min_frac=0.3):
                                continue  # already modelled (re-rank owns it)
                            wid = ("%s__disc__%s_%d_%d_%s" %
                                   (genome, orf['chrom'], orf['start'],
                                    orf['end'], orf['strand']))
                            if wid in orf_records:
                                continue
                            items.append(gpu_augment.WorkItem(
                                'orf', wid, orf['seq'], genome,
                                {'chrom': orf['chrom'], 'start': orf['start'],
                                 'end': orf['end'], 'strand': orf['strand']}))
                            orf_records[wid] = (genome, gff_path, orf)
                    del genome_seqs
                except Exception as exc:
                    logger.error("[gpu_augment] discovery failed for %s "
                                 "(continuing): %s", genome, exc)
            if orf_records:
                logger.info("[gpu_augment] discovery: %d ORF candidate(s) queued "
                            "across %d genome(s).", len(orf_records), len(genome_gffs))

    if not items:
        logger.info("[gpu_augment] no candidate models or ORFs to score; skipping.")
        return

    # Ensure nothing is resident so gpu_augment loads models sequentially.
    try:
        import plm_search as _p
        import structural_search as _s
        _p.unload_model()
        _s.unload_esmfold()
    except Exception:
        pass

    device = (args.structural_device if getattr(args, 'enable_structural_search', False)
              else getattr(args, 'plm_device', 'cpu'))
    work_dir = os.path.join(args.output_dir, "gpu_augment_work")
    logger.info(
        "[gpu_augment] scoring %d GOI candidate model(s) across %d genome(s) "
        "on device=%s (PLM=%s, structural=%s) ...",
        len(items), len(genome_gffs), device,
        getattr(args, 'enable_plm_search', False),
        getattr(args, 'enable_structural_search', False),
    )
    try:
        result = gpu_augment.run_augmentation(
            db_fasta=args.initial_db,
            items=items,
            work_dir=work_dir,
            device=device,
            do_plm=getattr(args, 'enable_plm_search', False),
            do_structural=getattr(args, 'enable_structural_search', False),
            require=fail_loud,
            plm_threshold=getattr(args, 'plm_similarity_threshold', 0.5),
            structural_tm_threshold=getattr(args, 'structural_tm_threshold', 0.3),
            fold_window=getattr(args, 'structural_fold_window', 380),
            fold_overlap=getattr(args, 'structural_fold_overlap', 60),
            threads=getattr(args, 'threads', 1),
            orf_fold_min_embedding=getattr(args, 'discovery_orf_fold_min_embedding', 0.4),
            orf_fold_max=getattr(args, 'discovery_orf_fold_max', 40),
        )
    except Exception as exc:
        if fail_loud:
            raise
        logger.warning("[gpu_augment] phase failed (continuing): %s", exc)
        return

    scores = result.get("scores", {})
    try:
        with open(os.path.join(args.output_dir, "gpu_augment_report.json"), 'w') as fh:
            json.dump(result.get("diagnostics", {}), fh, indent=2)
    except Exception:
        pass

    total_attached = 0
    total_boosted = 0
    for genome, gff_path in sorted(genome_gffs.items()):
        attached, boosted = _apply_scores_to_region_gff(gff_path, scores)
        total_attached += attached
        total_boosted += boosted

    # ── Materialise structural discovery hits as new GOI gene models ──
    # Conservative: only a strong structural match (TM >= discovery_min_tm,
    # default = structural_high_threshold) becomes a model, so discovery cannot
    # pollute a clean close-relative run with false positives. Appends AFTER the
    # re-rank rewrite so the new models survive.
    disc_hits = result.get("discovery_hits", [])
    materialised = 0
    if disc_hits and orf_records:
        try:
            import discover_regions as dr
        except Exception:
            dr = None
        if dr is not None:
            from collections import defaultdict
            tm_bar = getattr(args, 'discovery_min_tm',
                             getattr(args, 'structural_high_threshold', 0.7))
            # GOI templates for the sequence guards (coverage/length/paralog).
            # Use only the FULL-length GOI template(s); the exon-split fragments
            # (`GOI_osk|exon_N`, as short as 33 aa) would let a shared-domain FP
            # score ~100% coverage against a tiny fragment and slip guard 2.
            try:
                _all_goi = gpu_augment.load_goi_templates_from_fasta(args.initial_db)
            except Exception:
                _all_goi = []
            goi_templates = [(cid, seq) for cid, seq in _all_goi
                             if "exon" not in cid.lower()]
            if not goi_templates:            # all were fragments → fall back
                goi_templates = _all_goi
            # Optional home-paralog panel (guard 3). Absent by default → guard 3
            # is a no-op; guards 1/2/4 still fire.
            panel = None
            panel_path = getattr(args, 'home_paralog_panel', None)
            if panel_path and os.path.exists(panel_path):
                try:
                    panel = [(cid, seq) for _h, cid, seq in parse_fasta(panel_path) if seq]
                except Exception:
                    panel = None
            guard_rejects = defaultdict(int)   # reason -> count
            new_by_genome = defaultdict(list)   # genome -> [(gff_lines, (id,seq))]
            gff_for_genome = {}
            for hit in disc_hits:
                rec = orf_records.get(hit.get("id"))
                if not rec:
                    continue
                tm = hit.get("structural_similarity")
                if tm is None or tm < tm_bar:
                    continue
                genome, gff_path, orf = rec
                # ── Discovery guards: reject promiscuous-fold false positives ──
                syn = genome_synteny.get(genome, {})
                flank = syn.get("flanking", [])
                sgc = syn.get("seq_goi_chroms", set())
                flank_support = dr.count_flanking_support(
                    orf, flank,
                    window_bp=getattr(args, 'discovery_window_gap', 100000))
                ok, reason, gmetrics = dr.validate_discovery_hit(
                    orf, goi_templates,
                    panel=panel,
                    flanking_support=flank_support,
                    genome_has_seq_goi=bool(sgc),
                    seq_goi_chroms=sgc,
                    min_len_ratio=getattr(args, 'discovery_min_len_ratio', 0.5),
                    max_len_ratio=getattr(args, 'discovery_max_len_ratio', 2.0),
                    min_goi_coverage=getattr(args, 'discovery_min_goi_coverage', 0.5),
                    paralog_margin=getattr(args, 'discovery_paralog_margin', 1.0),
                    min_block_flanking=getattr(args, 'discovery_min_block_flanking', 2),
                    require_same_chrom_as_seq_goi=not getattr(
                        args, 'discovery_allow_offblock', False))
                if not ok:
                    guard_rejects[reason] += 1
                    logger.info("[gpu_augment] discovery hit %s rejected (%s): "
                                "%s", hit.get("id"), reason,
                                {k: gmetrics[k] for k in
                                 ("orf_len", "goi_len", "goi_coverage",
                                  "flanking_support") if k in gmetrics})
                    continue
                conf = dr.confidence_from_scores(
                    tm, hit.get("embedding_similarity"),
                    tm_high=getattr(args, 'structural_high_threshold', 0.7),
                    tm_medium=getattr(args, 'structural_medium_threshold', 0.5))
                idx = len(new_by_genome[genome])
                gff_lines, faa = dr.build_discovery_gff(
                    orf, genome, idx, tm, hit.get("embedding_similarity"), conf)
                new_by_genome[genome].append((gff_lines, faa))
                gff_for_genome[genome] = gff_path
                materialised += 1
            if guard_rejects:
                logger.info("[gpu_augment] discovery guards rejected %d hit(s): %s",
                            sum(guard_rejects.values()), dict(guard_rejects))
                result.setdefault("diagnostics", {})["discovery_guard_rejects"] = \
                    dict(guard_rejects)
            for genome, entries in new_by_genome.items():
                gpath = gff_for_genome[genome]
                fpath = gpath[:-4] + ".faa"
                with open(gpath, 'a') as gf:
                    for gff_lines, _faa in entries:
                        for gl in gff_lines:
                            gf.write(gl + "\n")
                with open(fpath, 'a') as ff:
                    for _gl, (mid, seq) in entries:
                        ff.write(f">{mid}\n")
                        for i in range(0, len(seq), 80):
                            ff.write(seq[i:i + 80] + "\n")
            if materialised:
                logger.info("[gpu_augment] discovery: materialised %d structural "
                            "model(s) (TM>=%.2f).", materialised, tm_bar)
    result.setdefault("diagnostics", {})["discovery_materialised"] = materialised

    logger.info(
        "[gpu_augment] scored %d candidate(s); attached ML scores to %d model(s), "
        "confidence-boosted %d.", len(scores), total_attached, total_boosted)


def main():
    parser = argparse.ArgumentParser(description="Iterative Genome Search Runner (Wavefront Parallel)")
    parser.add_argument("--initial_db", required=True)
    parser.add_argument("--sorted_genomes", required=True, 
                        help="Tab-separated file: genome_path\\tdistance. Genomes sorted by distance.")
    parser.add_argument("--genomes_dir", help="Directory containing genome files (if paths in sorted_genomes are relative)")
    parser.add_argument("--home_db_dir", help="Home Proteome MMseqs DB for RBH")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min_identity", type=float, default=40.0)
    parser.add_argument("--min_length", type=int, default=50)
    parser.add_argument("--max_intron", type=int, default=20000)
    parser.add_argument("--goi_info", default=None,
                        help="ANNOTATE_GOI goi_info.json — supplies the home GOI's "
                             "genomic footprint so the GOI fallback hit-chain max-gap "
                             "is derived from the gene's own introns instead of a fixed "
                             "--max_intron. Optional; falls back to --max_intron when absent.")
    parser.add_argument("--goi_fallback_intron_margin", type=float, default=5.0,
                        help="Multiplier applied to the home GOI's largest intron to "
                             "bound the GOI fallback hit-chain (big = lenient).")
    parser.add_argument("--goi_fallback_intron_floor", type=int, default=10000,
                        help="Minimum GOI fallback hit-chain max-gap (bp); covers compact "
                             "or single-exon home genes whose largest intron is ~0.")
    parser.add_argument("--threads", type=int, default=4, help="Total threads available for parallel processing")
    # Distances that describe a neighbourhood are given in gene pitches (the mean step
    # between two anchors at home, times the expansion measured in each target genome;
    # see synvoy_scale.py). Each has a base-pair option that fixes it: 'auto' = derive.
    parser.add_argument("--block_gap_genes", type=float, default=7.0,
                        help="Two flanking genes further apart than this many gene pitches "
                             "start a new block (default: 7)")
    parser.add_argument("--cluster_distance", type=bp_or_auto, default=None,
                        help="The same limit in bp, for every genome. 'auto' (default) = "
                             "derive it from --block_gap_genes")
    parser.add_argument("--window_margin_genes", type=float, default=7.0,
                        help="DNA searched on each side of a block, in gene pitches (default: 7)")
    parser.add_argument("--goi_margin_genes", type=float, default=2.0,
                        help="A GOI hit must lie in the block or at most this many gene "
                             "pitches outside it (default: 2)")
    parser.add_argument("--goi_margin", type=bp_or_auto, default=None,
                        help="The same limit in bp. 'auto' (default) = derive it")
    parser.add_argument("--synteny_bridge_max_gap_genes", type=float, default=20.0,
                        help="Largest gap a collinear block may bridge, in gene pitches; "
                             "0 switches bridging off (default: 20)")
    parser.add_argument("--dispersed_envelope_margin_genes", type=float, default=9.0,
                        help="Margin of the flanking envelope for a GOI hit away from its "
                             "anchors, in gene pitches (default: 9)")
    parser.add_argument("--disable_relative_query", type=str2bool, default=False,
                        help="Search each window with the home protein only; do not add the "
                             "gene model of the relative searched most recently")
    parser.add_argument("--disable_anchor_resolution", type=str2bool, default=False,
                        help="Count every flanking gene with a hit, also where several "
                             "flanking genes hit the same target DNA (behaviour before 2026-10-08)")
    parser.add_argument("--anchor_max_places", type=int, default=0,
                        help="Each flanking gene seeds blocks only at its N best places in a "
                             "target genome (by summed bit score); 0 = at every place (default)")
    parser.add_argument("--goi_max_models_per_locus", type=int, default=10,
                        help="Gene models kept per GOI locus of a window: every miniprot model "
                             "that does not lie on a better one, up to this number; 1 = the "
                             "best model only (default: 10)")
    parser.add_argument("--mmseqs_sens", type=float, default=7.5, help="MMseqs2 sensitivity (higher = more sensitive but slower)")
    parser.add_argument("--mmseqs_split_memory_limit", default="0", help="MMseqs split memory limit (e.g. 3G, 8000M, 0=auto)")
    parser.add_argument("--mmseqs_verbosity", type=int, default=1, help="MMseqs verbosity (0-3)")
    parser.add_argument(
        "--deterministic_goi_search", type=str2bool, default=True,
        help="A0 reproducibility: force the per-region augmented GOI MMseqs search to "
             "--threads 1 so marginal/divergent GOI hits don't jitter run-to-run "
             "(tblastn is already single-threaded; regions are tiny so the speed cost "
             "is small). Set false to restore the multithreaded (non-deterministic) path."
    )
    parser.add_argument("--region_padding", type=bp_or_auto, default=None,
                        help="Window margin in bp. 'auto' (default) = derive it from "
                             "--window_margin_genes. With a number the margin of each block is "
                             "adapted to the spacing of its hits between --padding_min and "
                             "--padding_max, as before 2026-10-08")
    parser.add_argument("--padding_min", type=int, default=50000,
                        help="Smallest margin, only with a fixed --region_padding")
    parser.add_argument("--padding_max", type=int, default=200000,
                        help="Largest margin, only with a fixed --region_padding")
    parser.add_argument("--enable_smith_waterman", type=str2bool, default=True, help="Enable Smith-Waterman search")
    parser.add_argument("--sw_method", type=str, default="auto", help="Smith-Waterman backend (auto, parasail, ssearch36)")
    parser.add_argument("--allow_missing_smith_waterman", type=str2bool, default=False,
                        help="Proceed (with only mmseqs/tblastn) when SW is requested but parasail is "
                             "unavailable, instead of failing loud. SW is load-bearing for divergent "
                             "GOI recovery, so this degrades results silently — opt in deliberately.")
    parser.add_argument("--sw_min_score", type=float, default=50.0)
    parser.add_argument("--sw_min_identity", type=float, default=20.0)
    parser.add_argument("--sw_timeout_seconds", type=int, default=300)
    parser.add_argument("--aug_relaxed_evalue_mult", type=float, default=1000.0)
    parser.add_argument("--aug_relaxed_evalue_cap", type=float, default=10.0)
    parser.add_argument("--aug_relaxed_parse_evalue_mult", type=float, default=10.0)
    parser.add_argument("--aug_relaxed_identity_factor", type=float, default=0.6)
    parser.add_argument("--aug_relaxed_identity_min", type=float, default=25.0)
    parser.add_argument("--aug_relaxed_length_div", type=float, default=2.0)
    parser.add_argument("--aug_relaxed_length_min", type=int, default=15)
    parser.add_argument("--gap_search_window", type=int, default=50000,
                        help="DNA added around a candidate's hits before the gene model is "
                             "built; never less than the gene's own span at home (bp)")
    parser.add_argument("--gap_min_identity", type=float, default=25.0)
    parser.add_argument("--gap_max_hits", type=int, default=5)
    parser.add_argument("--min_exon_alnlen", type=int, default=30,
                        help="Shortest alignment (aa) accepted for a last-resort raw_hit row "
                             "of the GOI; never less than 60 %% of the query")
    parser.add_argument(
        "--max_blocks_per_genome",
        type=int,
        default=80,
        help="Maximum synteny blocks to process per genome (0 = no cap)"
    )
    parser.add_argument(
        "--min_block_genes",
        type=int,
        default=2,
        help="Minimum unique anchor genes required for a synteny block"
    )
    # §1m collinearity-aware seed placement (proper fix for the flanking-gap fumble).
    parser.add_argument(
        "--disable_synteny_collinearity", type=str2bool, default=False,
        help="Disable collinearity-aware block scoring + gap-bridging (legacy gene-count clustering)."
    )
    parser.add_argument(
        "--synteny_bridge_max_gap", type=bp_or_auto, default=None,
        help="Largest gap a collinear block may bridge, in bp, for every genome. 'auto' "
             "(default) = derive it from --synteny_bridge_max_gap_genes; 0 disables bridging."
    )
    parser.add_argument(
        "--synteny_bridge_max_rank_gap", type=int, default=5,
        help="Max home-rank jump still treated as a collinear continuation when bridging."
    )
    parser.add_argument(
        "--synteny_bridge_min_anchors", type=int, default=2,
        help="Min distinct anchored flanking genes a block needs before it may bridge a gap "
             "(2 since 2026-10-14, was 3). Also the minimum collinear run of the two-sided "
             "rule (--synteny_bridge_two_sided, off)."
    )
    parser.add_argument(
        "--synteny_bridge_two_sided", type=str2bool, nargs='?', const=True, default=False,
        help="Weigh the whole cluster on the FAR side of a gap when deciding to bridge, "
             "not just its first locus. Fixes the asymmetry where identical collinear "
             "evidence bridges or not depending on where the gap falls "
             "(docs/CAN_BRIDGE_ANALYSIS.md). Default off pending the A/B measurement."
    )
    parser.add_argument(
        "--synteny_bridge_max_per_block", type=int, default=2,
        help="Max gaps a single block may bridge (0 = unlimited). Each bridge can span "
             "synteny_bridge_max_gap, so this bounds a runaway chained block."
    )
    # §A1 (TODO_JUN) — "the fumble": short-query fallback gate keyed on aligned
    # length + bitscore instead of query coverage (a short mature peptide in a
    # longer precursor is structurally low-coverage even on a perfect hit).
    parser.add_argument(
        "--fallback_short_query_len", type=int, default=150,
        help="§A1: queries shorter than this (aa) gate their GOI fallback on aligned "
             "length + bitscore, NOT query coverage (the mature peptide may be a small "
             "fraction of a precursor). 0 = legacy coverage gate for all queries."
    )
    parser.add_argument(
        "--fallback_short_min_aln_aa", type=int, default=15,
        help="§A1: minimum aligned length (aa) for a short-query GOI fallback; rejects "
             "micro-window overcalls while admitting short mature peptides (e.g. 26-aa melittin)."
    )
    parser.add_argument(
        "--fallback_short_min_bits", type=float, default=30.0,
        help="§A1: minimum bitscore for a short-query GOI fallback (significance floor that "
             "lets a high-bitscore divergent hit pass where coverage would have vetoed it)."
    )
    parser.add_argument(
        "--quiet_subtools",
        type=str2bool,
        default=True,
        help="Suppress noisy third-party tool output (miniprot/tblastn diagnostics)"
    )
    parser.add_argument("--prefix", default="", help="Prefix for output files (e.g. locus ID)")
    parser.add_argument("--resume", action="store_true", help="Resume from previous checkpoint if available")

    # Classification thresholds (GOI confidence / model status)
    parser.add_argument("--classify_high_min_identity", type=float, default=50.0,
                        help="Min identity for HIGH confidence exon_annotation "
                             "(lowered from 60 to accommodate cross-vertebrate orthologs)")
    parser.add_argument("--classify_medium_min_identity", type=float, default=40.0,
                        help="Min identity for MEDIUM confidence exon_annotation (40: in five "
                             "benchmark families no gene model below 40 %% identity lay on a "
                             "curated gene)")
    parser.add_argument("--classify_tandem_min_identity", type=float, default=50.0,
                        help="Min identity for MEDIUM confidence tandem_copy (below = LOW)")
    parser.add_argument("--classify_high_min_collinear", type=int, default=3,
                        help="§1: HIGH also requires the longest collinear run of home-ordered "
                             "flanking to reach this, once a block has >= this many flanking "
                             "genes — so a scrambled paralog neighbourhood is demoted to MEDIUM. "
                             "Set 0 to disable the order gate (count-only, legacy).")
    parser.add_argument("--classify_context_min_collinear", type=lambda v: int(float(v)), default=0,
                        help="A genome with no HIGH/MEDIUM call for a query gets its best "
                             "LOW/AMBIGUOUS exon model or hit chain lifted to MEDIUM when its "
                             "block keeps at least this many anchors in home order "
                             "(InferenceReason low_identity_in_conserved_gene_order). "
                             "0 = off (default).")
    parser.add_argument("--classify_context_min_identity", type=float, default=15.0,
                        help="Identity floor (%%) for that promotion.")
    parser.add_argument("--classify_context_min_qcov", type=float, default=0.5,
                        help="Query coverage a model needs for that promotion "
                             "(never below --classify_fragment_max_qcov).")
    parser.add_argument("--report_nonsyntenic_candidates", type=str2bool, default=True,
                        help="§1x: record GOI hits that clear the quality bar but are "
                             "refused by a synteny gate to hits/<genome>.nonsyntenic.tsv, "
                             "so the report can distinguish 'no ortholog found' from "
                             "'found one, could not place it syntenically'. Default on.")
    parser.add_argument("--classify_fallback_strong_min_identity_floor", type=float,
                        default=30.0,
                        help="§1v: absolute identity floor for the strong-flanking "
                             "fallback->MEDIUM clause. Its qcov arm has no identity "
                             "requirement, so without this a hit at any identity reached "
                             "MEDIUM inside a flanking-rich block (yeast STE2: 32 MEDIUM "
                             "at 21-27%%). Default 30 keeps the 34.7%% ant-melittin call "
                             "and rejects the 21-27%% false positives. 0 = legacy.")
    parser.add_argument("--classify_tandem_min_qcov", type=float, default=0.35,
                        help="§A2: min query coverage for MEDIUM tandem_copy; a high-identity "
                             "but short local window (true low qcov) is demoted to LOW.")
    # §A3 (TODO_JUN) — un-neuter the wave seed (opt-in). A HIGH/MEDIUM GOI feature
    # seeds the expanding DB on strong flanking + real coverage regardless of
    # goi_class, so a close relative's tandem/divergent ortholog can bridge to a
    # more divergent species in the next wave.
    parser.add_argument("--seed_on_flanking_support", type=str2bool, default=False,
                        help="§A3: also seed the next wave with a HIGH/MEDIUM GOI feature "
                             "(any goi_class, incl. tandem_goi_copy) when flanking support and "
                             "query coverage clear the floors below. Default off until validated.")
    parser.add_argument("--seed_flanking_min_count", type=int, default=2,
                        help="§A3: min flanking-gene support for flanking-supported seeding.")
    parser.add_argument("--seed_flanking_min_qcov", type=float, default=0.5,
                        help="§A3: min query coverage for flanking-supported seeding.")
    # §A4 (TODO_JUN) — rank/quantile wave binning (opt-in). The phylo-distance
    # metric saturates for tight clades (all bees → dist≈0.99 after max-norm),
    # collapsing absolute-threshold binning into ~4 undifferentiated waves with
    # no close→far gradient. When on, waves are assigned by RANK (genome_entries
    # is pre-sorted closest-first): closest genomes in small serial waves for
    # maximum seed propagation (the §A3 precondition), distant genomes in larger
    # waves for parallelism. Deterministic (pure function of the sorted order).
    parser.add_argument("--rank_wave_binning", type=str2bool, default=True,
                        help="§A4 / F4: bin waves by phylo-distance RANK/quantile instead "
                             "of absolute-distance thresholds. Default ON since 2026-07-26: "
                             "the legacy binning compares ABSOLUTE cut-offs (0.05/0.15) "
                             "against a divide-by-max NORMALISED distance, which is "
                             "incoherent by construction — the farthest genome is always "
                             "1.0, so the closest reaches the serial tier only if it is "
                             ">=20x nearer. Measured over 44 runs on disk: legacy gave the "
                             "closest genome a serial wave in 7, rank binning in 44. Pass "
                             "false to restore the legacy behaviour.")
    # F4 control arm: the wavefront is the mechanism the method is named for, but no
    # run has ever demonstrated that it changes the result. Comparing rank-vs-legacy
    # binning cannot answer that (both still iterate); the control is NO iteration at
    # all — every genome searched with the initial DB only. This flag provides it.
    parser.add_argument("--disable_ambiguous_tier", type=str2bool, default=False,
                        help="Restore the pre-2026-08-31 behaviour where a candidate "
                             "promoted only by flanking support was labelled "
                             "MEDIUM/probable_goi. By default such a call is AMBIGUOUS: "
                             "measured on the melittin benchmark, identity and flanking "
                             "support are indistinguishable between species that HAVE the "
                             "gene and species that LOST it, so 'probable ortholog' "
                             "asserts more than the evidence supports.")
    parser.add_argument("--disable_wavefront", type=str2bool, default=False,
                        help="F4 control arm: place ALL genomes in ONE parallel wave, so "
                             "no genome is ever searched with a nearer relative's recovered "
                             "ortholog (the expanding database never expands). This is the "
                             "measurement baseline for 'does the wavefront earn its place'; "
                             "it is not a performance option. Overrides --rank_wave_binning.")
    parser.add_argument("--classify_fragment_max_qcov", type=float, default=0.4,
                        help="Query coverage below this marks model as fragment")
    parser.add_argument("--classify_complete_min_qcov", type=float, default=0.7,
                        help="Query coverage above this (with multi-exon) marks model as complete")
    parser.add_argument("--strict_goi_family", type=str2bool, default=False,
                        help="Downgrade fallback/rescued_exon/raw_hit GOI calls whose annotated "
                             "TargetGene/TargetProduct does not match the family tokens")
    parser.add_argument("--goi_family_tokens", type=str, default="",
                        help="Comma-separated family name tokens for strict mode; "
                             "without them the family-name check is off")

    # PLM (Protein Language Model) embedding search
    parser.add_argument("--enable_plm_search", type=str2bool, default=False,
                        help="Enable ProtT5 embedding search for remote homolog detection")
    parser.add_argument("--plm_device", type=str, default="cpu",
                        help="Device for PLM inference (cpu or cuda)")
    parser.add_argument("--plm_similarity_threshold", type=float, default=0.5,
                        help="Minimum cosine similarity for PLM-discovered ORFs (0-1)")
    parser.add_argument("--plm_medium_threshold", type=float, default=0.7,
                        help="Embedding similarity above this can boost LOW → MEDIUM")
    parser.add_argument("--plm_high_threshold", type=float, default=0.85,
                        help="Embedding similarity above this can boost MEDIUM → HIGH")

    # Gene predictor selection (Augustus for eukaryotes, Prodigal for prokaryotes)
    parser.add_argument("--gene_predictor", type=str, default="auto",
                        help="Gene predictor: auto (prefer Augustus), augustus, or prodigal")
    parser.add_argument("--augustus_species", type=str, default="fly",
                        help="Augustus species model (e.g. fly, human, honeybee1)")

    # Structural search (ESMFold + Foldseek)
    parser.add_argument("--enable_structural_search", type=str2bool, default=False,
                        help="Enable ESMFold + Foldseek structural search for remote homologs")
    parser.add_argument("--structural_device", type=str, default="cpu",
                        help="Device for ESMFold inference (cpu or cuda)")
    parser.add_argument("--structural_tm_threshold", type=float, default=0.3,
                        help="Minimum TM-score for Foldseek-discovered ORFs (0-1)")
    parser.add_argument("--structural_medium_threshold", type=float, default=0.5,
                        help="TM-score above this can boost LOW → MEDIUM")
    parser.add_argument("--structural_high_threshold", type=float, default=0.7,
                        help="TM-score above this can boost MEDIUM → HIGH")
    parser.add_argument("--structural_max_length", type=int, default=700,
                        help="Max sequence length for ESMFold (VRAM safety, default 700)")
    parser.add_argument("--structural_fold_window", type=int, default=380,
                        help="Windowed-fold size for the main-process GPU phase; long "
                             "proteins are tiled into overlapping windows (each <= VRAM "
                             "cap) so no domain is truncated (default 380)")
    parser.add_argument("--structural_fold_overlap", type=int, default=60,
                        help="Overlap between fold windows (default 60)")
    parser.add_argument("--gpu_augment_fail_loud", type=str2bool, default=True,
                        help="Fail the run if PLM/structural search is requested but its "
                             "deps/GPU are unavailable or every fold/embed fails, instead "
                             "of silently producing no ML scores (default True)")
    # Structural/PLM DISCOVERY: predict ORFs in each syntenic window, fold the
    # promising ones, and materialise strong structural matches as new GOI models
    # (finds orthologs sequence search missed — Run-2 divergence ladder). Opt-in.
    parser.add_argument("--enable_structural_discovery", type=str2bool, default=False,
                        help="Predict + fold region ORFs to discover sequence-invisible "
                             "GOI orthologs (needs --enable_structural_search; default off)")
    parser.add_argument("--discovery_window_gap", type=int, default=100000,
                        help="Merge flanking genes within this bp into one discovery window")
    parser.add_argument("--discovery_window_pad", type=int, default=20000,
                        help="Pad each discovery window by this many bp")
    parser.add_argument("--discovery_min_aa", type=int, default=30,
                        help="Minimum ORF length (aa) for discovery")
    parser.add_argument("--discovery_min_tm", type=float, default=0.7,
                        help="Minimum Foldseek TM-score to materialise a discovery model "
                             "(conservative default = structural_high_threshold)")
    parser.add_argument("--discovery_orf_fold_min_embedding", type=float, default=0.4,
                        help="Only fold discovery ORFs whose ProtT5 cosine >= this (cost guard)")
    parser.add_argument("--discovery_orf_fold_max", type=int, default=40,
                        help="Cap on discovery ORFs folded per run (top by embedding)")
    # Discovery GUARDS — a structural TM match alone cannot tell an ortholog from
    # a fold analogue (oskar's OSK domain is a lipase fold; a 2187-aa lipase folded
    # to TM 0.80 and was materialised as a HIGH 'oskar' at 5.8% query coverage).
    # These add the sequence + synteny sanity the fold score lacks.
    parser.add_argument("--discovery_min_len_ratio", type=float, default=0.5,
                        help="Guard 1: reject a discovery ORF shorter than this x the GOI length")
    parser.add_argument("--discovery_max_len_ratio", type=float, default=2.0,
                        help="Guard 1: reject a discovery ORF longer than this x the GOI length")
    parser.add_argument("--discovery_min_goi_coverage", type=float, default=0.5,
                        help="Guard 2: SW alignment must span >= this fraction of the GOI "
                             "(separates a diverged ortholog from a shared-domain fold analogue)")
    parser.add_argument("--discovery_paralog_margin", type=float, default=1.0,
                        help="Guard 3: GOI SW score must beat the best home-paralog score by "
                             "this factor (needs --home_paralog_panel; else no-op)")
    parser.add_argument("--discovery_min_block_flanking", type=int, default=2,
                        help="Guard 4a: require >= this many flanking anchors near the ORF's window")
    parser.add_argument("--discovery_allow_offblock", type=str2bool, default=False,
                        help="Guard 4b: if True, allow discovery on a chromosome without a "
                             "sequence-found GOI even when the genome has one elsewhere (default off)")
    parser.add_argument("--home_paralog_panel", default=None,
                        help="Optional home-paralog panel FASTA for discovery guard 3 (RBH). "
                             "When absent, guard 3 is skipped and guards 1/2/4 still apply.")
    # §1o — dispersed-GOI rescue: recover a strong GOI hit that drifted away from
    # its flanking anchors (same-chromosome macro-synteny) in a shattered
    # neighbourhood. Gated by rearrangement score so it is a no-op on well-
    # conserved targets (the melittin/LY6/TP53 benchmarks).
    parser.add_argument("--enable_dispersed_goi_rescue", type=str2bool, default=True,
                        help="§1o: seed a block from a strong GOI hit that dispersed away "
                             "from flanking but shares a chromosome with an anchor, when the "
                             "neighbourhood is rearranged (default on)")
    parser.add_argument("--dispersed_goi_min_identity", type=float, default=40.0,
                        help="§1o: minimum %%identity for a dispersed GOI hit to be rescued")
    parser.add_argument("--dispersed_goi_max_evalue", type=float, default=1e-10,
                        help="§1o: maximum e-value for a dispersed GOI hit to be rescued")
    parser.add_argument("--dispersed_goi_min_alnlen", type=int, default=50,
                        help="§1o: minimum alignment length (aa/bp) for a dispersed GOI hit")
    parser.add_argument("--dispersed_min_chrom_anchors", type=int, default=3,
                        help="§1o: a chromosome must carry >= this many distinct flanking "
                             "genes to define a rescue envelope (macro-synteny evidence)")
    parser.add_argument("--dispersed_envelope_margin", type=bp_or_auto, default=None,
                        help="§1o: pad the flanking envelope by this many bp on each side. "
                             "'auto' (default) = derive it from --dispersed_envelope_margin_genes")

    # §1g: per-target phylogenetic-distance auto-tune of the GOI classify bars.
    parser.add_argument("--disable_distance_autotune", type=str2bool, default=False,
                        help="§1g: disable per-target stringency auto-tune (keep one global stringency)")
    parser.add_argument("--distance_autotune_close_pct", type=float, default=70.0,
                        help="§1g: median flanking identity (pct) >= this -> no relaxation (close relative)")
    parser.add_argument("--distance_autotune_far_pct", type=float, default=40.0,
                        help="§1g: median flanking identity (pct) <= this -> full relaxation + manual_review")
    parser.add_argument("--distance_autotune_max_relax", type=float, default=10.0,
                        help="§1g: max identity-points subtracted from the HIGH/MEDIUM bars for distant targets")
    parser.add_argument("--distance_autotune_min_flanking", type=int, default=3,
                        help="§1g: minimum flanking genes needed to trust the distance estimate")
    parser.add_argument("--distance_autotune_high_floor", type=float, default=45.0,
                        help="§1g: the relaxation never takes the HIGH identity threshold below this "
                             "(a baseline at or below it is not relaxed); 0 = no floor")
    parser.add_argument("--distance_autotune_medium_floor", type=float, default=30.0,
                        help="§1g: the same for the MEDIUM identity threshold; 0 = no floor")

    args = parser.parse_args()

    # Populate classification thresholds from CLI args
    CLASSIFY_THRESHOLDS["ambiguous_tier_enabled"] = not args.disable_ambiguous_tier
    # "At least two exons" for HIGH follows the home gene: a gene without introns at home
    # is expected without introns in the target, and could otherwise never be HIGH.
    args.goi_home_exons = None
    if getattr(args, "goi_info", None):
        try:
            with open(args.goi_info) as _fh:
                args.goi_home_exons = int(json.load(_fh).get("num_exons") or 0) or None
        except (OSError, TypeError, ValueError):
            args.goi_home_exons = None
    if args.goi_home_exons == 1:
        CLASSIFY_THRESHOLDS["high_min_exons"] = 1
        CLASSIFY_THRESHOLDS["high_single_exon_min_qcov"] = args.classify_complete_min_qcov
    CLASSIFY_THRESHOLDS["high_min_identity"] = args.classify_high_min_identity
    CLASSIFY_THRESHOLDS["medium_min_identity"] = args.classify_medium_min_identity
    CLASSIFY_THRESHOLDS["tandem_min_identity"] = args.classify_tandem_min_identity
    CLASSIFY_THRESHOLDS["tandem_min_qcov"] = args.classify_tandem_min_qcov
    # §1: 0 disables the order gate. Use a value larger than any real block so the
    # `flanking_support < high_min_collinear` waiver always fires (→ order never gates).
    CLASSIFY_THRESHOLDS["high_min_collinear"] = (
        args.classify_high_min_collinear if args.classify_high_min_collinear > 0 else 10 ** 9)
    CLASSIFY_THRESHOLDS["fallback_strong_min_identity_floor"] = (
        args.classify_fallback_strong_min_identity_floor)
    CLASSIFY_THRESHOLDS["fragment_max_qcov"] = args.classify_fragment_max_qcov
    CLASSIFY_THRESHOLDS["context_min_collinear"] = max(0, args.classify_context_min_collinear)
    CLASSIFY_THRESHOLDS["context_min_identity"] = args.classify_context_min_identity
    CLASSIFY_THRESHOLDS["context_min_qcov"] = args.classify_context_min_qcov
    CLASSIFY_THRESHOLDS["complete_min_qcov"] = args.classify_complete_min_qcov
    CLASSIFY_THRESHOLDS["plm_medium_threshold"] = args.plm_medium_threshold
    CLASSIFY_THRESHOLDS["plm_high_threshold"] = args.plm_high_threshold
    CLASSIFY_THRESHOLDS["structural_medium_threshold"] = args.structural_medium_threshold
    CLASSIFY_THRESHOLDS["structural_high_threshold"] = args.structural_high_threshold

    # Family-consistency config (strict GOI mode)
    FAMILY_CONFIG["strict"] = bool(args.strict_goi_family)
    override_tokens = {
        _normalize_family_token(t) for t in (args.goi_family_tokens or "").split(",")
    }
    override_tokens.discard("")
    # The names come from --goi_family_tokens only. (Until 2026-10-09 they were also
    # "derived" from the first header of the query set, which is a flanking exon record
    # written by the pipeline itself, so nothing was ever derived.) Three presets set
    # strict_goi_family; without names that is a plain no-op, not a problem to warn about.
    FAMILY_CONFIG["tokens"] = override_tokens
    if FAMILY_CONFIG["strict"] and not FAMILY_CONFIG["tokens"]:
        logger.info("Family-name check: off (no --goi_family_tokens given).")
        FAMILY_CONFIG["strict"] = False
    if FAMILY_CONFIG["tokens"]:
        logger.info(
            f"GOI family tokens: {sorted(FAMILY_CONFIG['tokens'])} "
            f"(strict={FAMILY_CONFIG['strict']})"
        )

    # INPUT VALIDATION
    # 1. Validate required files exist
    if not os.path.exists(args.initial_db):
        logger.error(
            f"Initial database file not found: {args.initial_db}. "
            f"This FASTA is produced by EXTRACT_FLANKING_GENES and contains "
            f"the flanking-gene proteins used as synteny parents. "
            f"Try: (1) rerun without -resume to regenerate it; "
            f"(2) check the Nextflow work dir for a failed EXTRACT_FLANKING_GENES "
            f"process; (3) ensure the home GFF has protein-coding genes "
            f"flanking the GOI."
        )
        sys.exit(1)

    if not os.path.exists(args.sorted_genomes):
        logger.error(
            f"Sorted genomes file not found: {args.sorted_genomes}. "
            f"This is the ranked target-genome list produced by "
            f"RANK_BY_SIMILARITY. "
            f"Try: rerun without -resume; if it still fails, check the "
            f"upstream RANK_BY_SIMILARITY process logs in Nextflow's work/ dir."
        )
        sys.exit(1)

    # 2. Validate initial_db is not empty
    if os.path.getsize(args.initial_db) == 0:
        logger.error(
            f"Initial database file is empty: {args.initial_db}. "
            f"EXTRACT_FLANKING_GENES produced a zero-byte FASTA, which means "
            f"no flanking genes were found. "
            f"Try: (1) verify the home GFF has 'gene' and 'CDS' features "
            f"around the GOI coordinates; "
            f"(2) increase --n_flanking_genes; "
            f"(3) check for GOI-vs-flanking similarity filtering being too "
            f"aggressive (see EXTRACT_FLANKING_GENES logs)."
        )
        sys.exit(1)

    # §1m: build the home-order rank map ONCE from the pristine initial_db (flanking
    # genes are written first, in home genomic order). Attached to `args` so it
    # pickles to every ProcessPoolExecutor worker (process_single_genome). When
    # disabled or empty, process_single_genome falls back to legacy clustering.
    if args.disable_synteny_collinearity:
        args.home_rank = {}
        logger.info("§1m collinearity-aware seed placement: DISABLED (legacy gene-count clustering).")
    else:
        args.home_rank = build_home_rank(args.initial_db)
        logger.info(
            f"§1m collinearity-aware seed placement: ENABLED "
            f"({len(args.home_rank)} flanking anchors ranked)."
        )

    # The GOI entries of the pristine query set are the home protein; everything that
    # joins them later is a gene model of a genome searched earlier (select_goi_queries).
    args.home_goi_ids = build_home_goi_ids(args.initial_db)

    # Home neighbourhood geometry for the distances in gene pitches (synvoy_scale).
    args.home_geometry = synvoy_scale.home_geometry(args.initial_db)
    if args.home_geometry:
        logger.info(
            f"Home neighbourhood: {args.home_geometry['n']} anchors over "
            f"{args.home_geometry['span'] / 1000:.0f} kb, pitch "
            f"{args.home_geometry['pitch'] / 1000:.1f} kb, largest gap "
            f"{args.home_geometry['max_gap'] / 1000:.0f} kb.")
    else:
        logger.info("Home neighbourhood: no anchor coordinates in the query set; "
                    "distances fall back to fixed base pairs.")
    args.goi_home_span = 0
    if getattr(args, "goi_info", None):
        try:
            with open(args.goi_info) as _fh:
                args.goi_home_span = int(json.load(_fh).get("genomic_span_bp") or 0)
        except (OSError, TypeError, ValueError):
            pass

    # GOI fallback hit-chain max-gap: derived from the home GOI's own largest
    # intron (goi_info.json) x margin, not a hand-set --max_intron. Attached to
    # `args` so it pickles to every ProcessPoolExecutor worker, like home_rank.
    args.goi_fallback_max_gap, _gap_reason = derive_goi_fallback_max_gap(
        getattr(args, "goi_info", None),
        args.goi_fallback_intron_margin,
        args.goi_fallback_intron_floor,
        args.max_intron,
    )
    logger.info(f"GOI fallback chain max-gap (data-driven max-intron): {_gap_reason}")

    # 3. Validate parameters are in valid ranges
    if args.min_identity < 0 or args.min_identity > 100:
        logger.error(f"Invalid min_identity: {args.min_identity}. Must be between 0 and 100")
        sys.exit(1)
    
    if args.min_length < 1:
        logger.error(f"Invalid min_length: {args.min_length}. Must be >= 1")
        sys.exit(1)
    
    if args.evalue <= 0:
        logger.error(f"Invalid evalue: {args.evalue}. Must be > 0")
        sys.exit(1)
    
    if args.threads < 1:
        logger.error(f"Invalid threads: {args.threads}. Must be >= 1")
        sys.exit(1)
    
    if args.mmseqs_sens < 1 or args.mmseqs_sens > 9:
        logger.warning(f"MMseqs sensitivity {args.mmseqs_sens} outside typical range (1-9)")
    if args.mmseqs_verbosity < 0 or args.mmseqs_verbosity > 3:
        logger.error("mmseqs_verbosity must be between 0 and 3")
        sys.exit(1)
    if args.padding_min < 0 or args.padding_max < 0 or (args.region_padding or 0) < 0:
        logger.error("Padding values must be non-negative")
        sys.exit(1)
    if args.padding_max < args.padding_min:
        logger.error("padding_max must be >= padding_min")
        sys.exit(1)
    if args.max_intron < 0:
        logger.error("max_intron must be >= 0")
        sys.exit(1)
    if args.gap_max_hits < 1:
        logger.error("gap_max_hits must be >= 1")
        sys.exit(1)
    if args.min_exon_alnlen < 1:
        logger.error("min_exon_alnlen must be >= 1")
        sys.exit(1)
    if args.max_blocks_per_genome < 0:
        logger.error("max_blocks_per_genome must be >= 0")
        sys.exit(1)
    if args.anchor_max_places < 0:
        logger.error("anchor_max_places must be >= 0 (0 = every place)")
        sys.exit(1)
    if args.min_block_genes < 1:
        logger.error("min_block_genes must be >= 1")
        sys.exit(1)
    if args.aug_relaxed_length_div <= 0:
        logger.error("aug_relaxed_length_div must be > 0")
        sys.exit(1)

    # Classification threshold validation
    for name in ['classify_high_min_identity', 'classify_medium_min_identity', 'classify_tandem_min_identity']:
        val = getattr(args, name)
        if val < 0 or val > 100:
            logger.error(f"{name} must be between 0 and 100 (got {val})")
            sys.exit(1)
    for name in ['classify_fragment_max_qcov', 'classify_complete_min_qcov']:
        val = getattr(args, name)
        if val < 0 or val > 1:
            logger.error(f"{name} must be between 0 and 1 (got {val})")
            sys.exit(1)
    if args.classify_fragment_max_qcov >= args.classify_complete_min_qcov:
        logger.warning(
            f"classify_fragment_max_qcov ({args.classify_fragment_max_qcov}) >= "
            f"classify_complete_min_qcov ({args.classify_complete_min_qcov}); "
            f"fragment and complete model status ranges overlap"
        )

    # Smith-Waterman validation
    if args.sw_min_score < 0:
        logger.error(f"sw_min_score must be >= 0 (got {args.sw_min_score})")
        sys.exit(1)
    if args.sw_min_identity < 0 or args.sw_min_identity > 100:
        logger.error(f"sw_min_identity must be between 0 and 100 (got {args.sw_min_identity})")
        sys.exit(1)
    if args.sw_timeout_seconds < 1:
        logger.error(f"sw_timeout_seconds must be >= 1 (got {args.sw_timeout_seconds})")
        sys.exit(1)

    # SW is the load-bearing tier for divergent GOI recovery (e.g. re-finding a 26-40 %
    # myrmicitoxin from a melittin query). Silently dropping it on a missing-parasail env —
    # the classic case is a `.venv` shadowing the conda env on PATH — changes results on
    # exactly the hard cases SynVoy exists for, with no loud signal. So in auto mode, a
    # requested-but-unavailable SW is a HARD ERROR unless the user explicitly opts into the
    # degraded path.
    _sw_before = args.enable_smith_waterman
    try:
        args.enable_smith_waterman = resolve_smith_waterman_enabled(
            _sw_before, args.sw_method,
            has_parasail_available(), args.allow_missing_smith_waterman,
        )
        if _sw_before and not args.enable_smith_waterman:  # degraded path was opted into
            logger.warning(
                "parasail not available in auto mode; Smith-Waterman DISABLED for this run "
                "(--allow_missing_smith_waterman set). Divergent GOI recovery is degraded to "
                "mmseqs/tblastn only."
            )
    except SmithWatermanUnavailable:
        venv_hint = os.environ.get("VIRTUAL_ENV")
        env_note = (
            f"\n  NOTE: a Python virtualenv is active (VIRTUAL_ENV={venv_hint}); it is likely "
            f"shadowing the conda env that has parasail. Deactivate it / strip it from PATH."
            if venv_hint else ""
        )
        logger.error(
            "Smith-Waterman was requested (--enable_smith_waterman true, --sw_method auto) "
            "but parasail is not importable in this environment. SW is load-bearing for "
            "divergent GOI recovery, so refusing to silently continue without it." + env_note +
            "\n  Fix one of:\n"
            "    - install parasail into the run's env  (conda install -c bioconda parasail-python)\n"
            "    - --sw_method ssearch36                 (use the slower FASTA/ssearch36 backend)\n"
            "    - --enable_smith_waterman false         (deliberately search without SW)\n"
            "    - --allow_missing_smith_waterman true   (proceed degraded, keep SW auto-off)"
        )
        sys.exit(1)

    # Gene predictor validation
    if args.gene_predictor not in ("auto", "augustus", "prodigal"):
        logger.error(f"gene_predictor must be 'auto', 'augustus', or 'prodigal' (got '{args.gene_predictor}')")
        sys.exit(1)

    # PLM search validation & availability check
    if args.enable_plm_search:
        if not PLM_IMPORT_OK:
            logger.error(
                "PLM search enabled but plm_search module could not be imported.\n"
                "  Ensure bin/plm_search.py exists and torch/transformers are installed:\n"
                "  pip install torch transformers sentencepiece"
            )
            sys.exit(1)
        if not check_plm_available():
            logger.error(
                "PLM search enabled but PyTorch/Transformers not installed.\n"
                "  Install with: pip install torch transformers sentencepiece"
            )
            sys.exit(1)
        if args.plm_similarity_threshold < 0 or args.plm_similarity_threshold > 1:
            logger.error("plm_similarity_threshold must be between 0 and 1")
            sys.exit(1)
        if args.plm_device not in ("cpu", "cuda"):
            logger.error("plm_device must be 'cpu' or 'cuda'")
            sys.exit(1)

    # Structural search (ESMFold + Foldseek) validation & availability check
    if args.enable_structural_search:
        if not STRUCTURAL_IMPORT_OK:
            logger.error(
                "Structural search enabled but structural_search module could not be imported.\n"
                "  Ensure bin/structural_search.py exists and dependencies are installed:\n"
                "  pip install torch transformers\n"
                "  conda install -c bioconda foldseek"
            )
            sys.exit(1)
        if not check_esmfold_available():
            logger.error(
                "Structural search enabled but ESMFold not available.\n"
                "  Install with: pip install torch transformers"
            )
            sys.exit(1)
        if not check_foldseek_available():
            logger.error(
                "Structural search enabled but Foldseek binary not found.\n"
                "  Install with: conda install -c bioconda foldseek"
            )
            sys.exit(1)
        if args.structural_tm_threshold < 0 or args.structural_tm_threshold > 1:
            logger.error("structural_tm_threshold must be between 0 and 1")
            sys.exit(1)
        if args.structural_device not in ("cpu", "cuda"):
            logger.error("structural_device must be 'cpu' or 'cuda'")
            sys.exit(1)
        if args.structural_max_length < 10:
            logger.error("structural_max_length must be >= 10")
            sys.exit(1)

    logger.info("")
    logger.info("═══ Iterative Search Configuration ═══")
    logger.info(f"  Threads:          {args.threads}")
    logger.info(f"  Min identity:     {args.min_identity}%")
    logger.info(f"  Min length:       {args.min_length}")
    logger.info(f"  E-value cutoff:   {args.evalue}")
    logger.info(f"  MMseqs sens:      {args.mmseqs_sens}")
    logger.info(f"  MMseqs mem limit: {args.mmseqs_split_memory_limit}")
    logger.info(f"  Gene predictor:   {args.gene_predictor}" + (f" (species={args.augustus_species})" if args.gene_predictor != "prodigal" else ""))
    logger.info(f"  Smith-Waterman:   {'enabled (' + args.sw_method + ')' if args.enable_smith_waterman else 'disabled'}")
    if args.enable_plm_search:
        logger.info(f"  PLM search:       enabled (ProtT5, device={args.plm_device})")
        logger.info(f"  PLM threshold:    {args.plm_similarity_threshold}")
        logger.info(f"  PLM boost:        MEDIUM≥{args.plm_medium_threshold}, HIGH≥{args.plm_high_threshold}")
    else:
        logger.info(f"  PLM search:       disabled")
    if args.enable_structural_search:
        logger.info(f"  Structural search: enabled (ESMFold+Foldseek, device={args.structural_device})")
        logger.info(f"  Struct TM thresh: {args.structural_tm_threshold}")
        logger.info(f"  Struct boost:     MEDIUM≥{args.structural_medium_threshold}, HIGH≥{args.structural_high_threshold}")
        logger.info(f"  Struct max len:   {args.structural_max_length}")
    else:
        logger.info(f"  Structural search: disabled")
    logger.info(f"  Max blocks/genome:{args.max_blocks_per_genome}")
    logger.info(f"  Region padding:   {args.padding_min}-{args.padding_max} bp")
    if args.prefix:
        logger.info(f"  Locus prefix:     {args.prefix}")
    logger.info("══════════════════════════════════════")
    
    prefix = f"{args.prefix}_" if args.prefix else ""
    
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/hits", exist_ok=True)
    os.makedirs(f"{args.output_dir}/regions", exist_ok=True)
    
    # CHECKPOINTING: Check for resume
    checkpoint_file = f"{args.output_dir}/.checkpoint"
    start_wave = 0
    current_db = f"{args.output_dir}/current_db.faa"
    
    if args.resume and os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as cf:
            checkpoint_data = json.loads(cf.read())
            start_wave = checkpoint_data.get('completed_waves', 0)
            total_saved_waves = checkpoint_data.get('total_waves', '?')
            last_db = checkpoint_data.get('last_db', None)

            if last_db and os.path.exists(last_db):
                logger.info(
                    f"Resuming from checkpoint: {start_wave}/{total_saved_waves} waves already completed. "
                    f"Continuing from wave {start_wave + 1}."
                )
                current_db = last_db
            else:
                logger.warning(
                    f"Checkpoint found ({start_wave} waves completed) but database file is missing "
                    f"(expected: {last_db}). Starting from beginning."
                )
                start_wave = 0
                shutil.copyfile(args.initial_db, current_db)
    else:
        if args.resume:
            logger.info("No checkpoint found; starting fresh.")
        shutil.copyfile(args.initial_db, current_db)
    
    # Parse Genomes and Distances
    genome_entries = []
    parse_warnings = []
    with open(args.sorted_genomes, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            gname = parts[0]
            if len(parts) < 2:
                parse_warnings.append(
                    f"  Line {line_num}: missing distance column (using 0.0): {line!r}"
                )
            try:
                dist = float(parts[1]) if len(parts) > 1 else 0.0
            except ValueError:
                logger.error(
                    f"sorted_genomes line {line_num}: cannot parse distance "
                    f"from {parts[1]!r}. Expected tab-separated: genome_path<TAB>distance"
                )
                sys.exit(1)

            gpath = gname
            if args.genomes_dir:
                if not os.path.isabs(gname):
                    gpath = os.path.join(args.genomes_dir, os.path.basename(gname))

            genome_entries.append({'name': gname, 'path': gpath, 'dist': dist})

    if parse_warnings:
        for w in parse_warnings:
            logger.warning(w)

    # Validate we loaded genomes
    if not genome_entries:
        logger.error(
            f"No genomes found in sorted_genomes file: {args.sorted_genomes}\n"
            f"  The file should contain one genome per line: genome_path<TAB>distance\n"
            f"  Check that the upstream PHYLO_SORT step completed successfully."
        )
        sys.exit(1)

    # Pre-flight: verify all genome files are accessible before starting the search
    missing_genomes = []
    for entry in genome_entries:
        if not os.path.exists(entry['path']):
            missing_genomes.append(f"  {entry['name']} -> {entry['path']}")
    if missing_genomes:
        logger.error(
            f"{len(missing_genomes)} of {len(genome_entries)} genome file(s) not found:\n"
            + "\n".join(missing_genomes[:10])
            + (f"\n  ... and {len(missing_genomes) - 10} more" if len(missing_genomes) > 10 else "")
            + f"\n  genomes_dir={args.genomes_dir or '(not set)'}"
            + "\n  Check that target genomes were downloaded/staged correctly."
        )
        sys.exit(1)

    logger.info(f"Loaded {len(genome_entries)} genomes (all files accessible).")

    # Normalize phylogenetic distances if they are not in [0,1]
    finite_dists = [g['dist'] for g in genome_entries if g.get('dist') not in [None, float('inf')]]
    if finite_dists:
        max_dist = max(finite_dists)
        if max_dist > 1.0:
            for g in genome_entries:
                if g['dist'] == float('inf'):
                    g['dist'] = 1.0
                else:
                    g['dist'] = g['dist'] / max_dist
            logger.info(f"Normalized phylogenetic distances by max {max_dist:.3f}")
    
    # Define Waves (see define_waves)
    waves = define_waves(genome_entries, args)

    logger.info(f"Defined {len(waves)} waves of execution.")

    # ── PLM: Pre-compute GOI embeddings before waves start ──
    # (Legacy in-worker path only. The main-process phase — run_gpu_augment_phase
    #  — prepares templates itself, sequentially, so we do NOT pre-load models
    #  here: co-resident ProtT5+ESMFold is exactly what forced the 400 aa
    #  truncation on a 16 GB card.)
    if GPU_IN_WORKER and args.enable_plm_search:
        goi_emb_path = os.path.join(args.output_dir, "goi_embeddings.npz")
        if not os.path.exists(goi_emb_path):
            logger.info("Pre-computing GOI embeddings with ProtT5 ...")
            goi_embs = precompute_goi_embeddings(
                db_fasta=args.initial_db,
                output_path=goi_emb_path,
                device=args.plm_device,
            )
            if not goi_embs:
                logger.warning(
                    "No GOI sequences found for PLM embedding. "
                    "PLM search will be skipped."
                )
                args.enable_plm_search = False
            else:
                logger.info(
                    f"GOI embeddings ready: {len(goi_embs)} sequence(s) → {goi_emb_path}"
                )
        else:
            logger.info(f"Using cached GOI embeddings: {goi_emb_path}")

    # ── Structural: Pre-fold GOI structures before waves start ──
    # (Legacy in-worker path only — see the PLM note above.)
    if GPU_IN_WORKER and args.enable_structural_search:
        goi_struct_index = os.path.join(args.output_dir, "goi_structure_index.tsv")
        if not os.path.exists(goi_struct_index):
            logger.info("Pre-folding GOI structures with ESMFold ...")
            goi_structs = prefold_goi_structures(
                db_fasta=args.initial_db,
                output_dir=args.output_dir,
                device=args.structural_device,
                max_length=args.structural_max_length,
            )
            if not goi_structs:
                logger.warning(
                    "No GOI sequences found for structural folding. "
                    "Structural search will be skipped."
                )
                args.enable_structural_search = False
            else:
                save_structure_index(goi_structs, goi_struct_index)
                logger.info(
                    f"GOI structures ready: {len(goi_structs)} structure(s) → {goi_struct_index}"
                )
        else:
            logger.info(f"Using cached GOI structures: {goi_struct_index}")

    # Compute cumulative genome index for progress tracking
    total_genomes = len(genome_entries)
    cumulative_genome_idx = 0  # genomes completed so far

    latest_db = current_db
    total_new_genes = 0
    genomes_with_hits = 0
    genomes_without_hits = 0
    # Tree-only GOI hits (e.g. tandem_goi_copy) — accumulated across waves
    # and written to goi_for_tree.faa so that COMPUTE_TREE sees species like
    # Apis cerana whose only GOI evidence is a tandem duplicate.
    tree_extras_total = []
    tree_extras_seen_ids = set()

    for i, wave in enumerate(waves):
        # Skip already completed waves
        if i < start_wave:
            cumulative_genome_idx += len(wave)
            logger.info(f"Skipping wave {i+1}/{len(waves)} (already completed, resuming)")
            continue

        wave_genome_names = [e['name'] for e in wave]
        rank_q = cumulative_genome_idx / max(1, total_genomes - 1)
        logger.info(
            f"═══ Wave {i+1}/{len(waves)} "
            f"[genomes {cumulative_genome_idx+1}-{cumulative_genome_idx+len(wave)}/{total_genomes}] "
            f"({len(wave)} genome(s), dist≈{wave[0]['dist']:.3f}, rank≈{rank_q:.2f}) ═══"
        )
        for entry in wave:
            gbase = os.path.basename(entry['name'])
            logger.info(f"  • {gbase} (dist={entry['dist']:.4f})")

        # Parallel Execution
        max_workers = min(len(wave), args.threads)
        threads_per_job = max(1, args.threads // max_workers)

        wave_results = []
        wave_errors = []
        # A0 reproducibility: collect each genome's contribution keyed by name
        # and assemble in the wave's fixed order AFTER the pool drains, instead
        # of in as_completed() order. A run-dependent expanded_db.faa byte-order
        # otherwise flips marginal-hit tie-breaks in the next wave's mmseqs
        # query DB → seed-set cascade (TODO_JUN §A0, pinned source #1).
        # Keyed by wave-position index (not name) so duplicate genome names in a
        # single wave can't collide and silently drop a genome's genes.
        wave_new_genes_by_idx = {}
        wave_tree_extras_by_idx = {}
        wave_start_time = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for idx, entry in enumerate(wave):
                future = executor.submit(
                    process_single_genome,
                    entry['path'], latest_db, args, args.home_db_dir, prefix, threads_per_job
                )
                futures[future] = (idx, entry)

            completed = 0
            for future in concurrent.futures.as_completed(futures):
                idx, entry = futures[future]
                completed += 1
                cumulative_genome_idx += 1
                gbase = os.path.basename(entry['name'])
                try:
                    gname, new_genes, tree_extras, genome_error = future.result()
                    if genome_error:
                        wave_errors.append((gname, genome_error))
                        logger.error(
                            f"  ✗ {gbase} FAILED ({completed}/{len(wave)} in wave, "
                            f"{cumulative_genome_idx}/{total_genomes} overall): {genome_error}"
                        )
                        continue
                    gene_count = len(new_genes) if new_genes else 0
                    extra_count = len(tree_extras) if tree_extras else 0
                    if gene_count > 0 or extra_count > 0:
                        genomes_with_hits += 1
                        msg = f"  ✓ {gbase}: {gene_count} new gene(s)"
                        if extra_count:
                            msg += f" + {extra_count} tree-only hit(s)"
                        logger.info(
                            f"{msg} ({cumulative_genome_idx}/{total_genomes} overall)"
                        )
                    else:
                        genomes_without_hits += 1
                        logger.info(
                            f"  – {gbase}: no new genes "
                            f"({cumulative_genome_idx}/{total_genomes} overall)"
                        )
                    if new_genes:
                        wave_new_genes_by_idx[idx] = new_genes
                    if tree_extras:
                        wave_tree_extras_by_idx[idx] = tree_extras
                except Exception as exc:
                    wave_errors.append((entry.get('name', 'unknown'), str(exc)))
                    logger.error(
                        f"  ✗ {gbase} CRASHED ({cumulative_genome_idx}/{total_genomes}): {exc}"
                    )

        # A0: deterministic assembly — walk the wave in its fixed genome order
        # and sort each genome's genes by id, so wave_results (→ expanded_db.faa
        # → next wave's mmseqs query) is byte-identical run-to-run.
        for idx in range(len(wave)):
            ng = wave_new_genes_by_idx.get(idx)
            if ng:
                wave_results.extend(sorted(ng, key=lambda g: g.get('id', '')))
            te = wave_tree_extras_by_idx.get(idx)
            if te:
                for g in sorted(te, key=lambda g: g.get('id', '')):
                    gid = g.get('id', '')
                    if gid and gid not in tree_extras_seen_ids:
                        tree_extras_seen_ids.add(gid)
                        tree_extras_total.append(g)

        wave_elapsed = time.time() - wave_start_time

        if wave_errors:
            logger.error(f"Wave {i+1} had {len(wave_errors)} error(s):")
            for gname, msg in wave_errors:
                logger.error(f"  • {os.path.basename(gname)}: {msg}")
            # Classify errors to give actionable advice
            oom_errors = [m for _, m in wave_errors if 'memory' in m.lower() or 'killed' in m.lower() or 'oom' in m.lower()]
            if oom_errors:
                logger.error(
                    "  Hint: some failures look like out-of-memory (OOM). Try:\n"
                    "    --mmseqs_split_memory_limit <lower value>  or\n"
                    "    increase the process memory in nextflow.config"
                )
            logger.error("Aborting iterative search due to per-genome processing errors.")
            sys.exit(1)

        # Update DB after Wave
        total_new_genes += len(wave_results)
        if wave_results:
            logger.info(
                f"Wave {i+1} complete: {len(wave_results)} new gene(s) in {wave_elapsed:.0f}s. "
                f"Running total: {total_new_genes} gene(s) from {genomes_with_hits} genome(s)."
            )

            new_genes_fasta = f"{args.output_dir}/iter_{i+1}_new_genes.faa"
            write_fasta([(g['id'], g['seq']) for g in wave_results], new_genes_fasta)

            next_db = f"{args.output_dir}/db_iter_{i+1}.faa"
            with open(next_db, 'w') as ndb:
                with open(latest_db, 'r') as old_db:
                    shutil.copyfileobj(old_db, ndb)
                with open(new_genes_fasta, 'r') as new_g:
                    shutil.copyfileobj(new_g, ndb)

            # Clean up previous DB if it's not the initial one
            if i > 0 and latest_db != current_db:
                try:
                    os.remove(latest_db)
                except OSError as e:
                    logger.warning(f"Could not remove old DB file {latest_db}: {e}")

            latest_db = next_db
        else:
            logger.info(f"Wave {i+1} complete: no new genes ({wave_elapsed:.0f}s).")

        # CHECKPOINT: Save progress after each wave
        with open(checkpoint_file, 'w') as cf:
            json.dump({
                'completed_waves': i + 1,
                'last_db': latest_db,
                'total_waves': len(waves)
            }, cf)

    # ── GPU augmentation (PLM + Foldseek), main process, single CUDA context ──
    # Runs after all waves so every genome's regions/<g>.gff + .faa exist. No-op
    # unless --enable_plm_search / --enable_structural_search. genome_map lets the
    # discovery path read each target FASTA (keyed by the region-file stem =
    # basename of the genome path).
    genome_map = {os.path.basename(e['path']): e['path'] for e in genome_entries}
    run_gpu_augment_phase(args, genome_map)

    expanded_db = f"{args.output_dir}/expanded_db.faa"
    if os.path.exists(latest_db):
        shutil.move(latest_db, expanded_db)

    # Build a tree-input FASTA = expanded_db + tree-only GOI hits (tandem_copy
    # and other MEDIUM/HIGH GOI hits that were withheld from wave seeding).
    # This is what COMPUTE_TREE consumes — keeping seeding strict but the tree
    # complete.
    goi_for_tree = f"{args.output_dir}/goi_for_tree.faa"
    if os.path.exists(expanded_db):
        shutil.copyfile(expanded_db, goi_for_tree)
    else:
        # Defensive: emit an empty file so downstream wiring doesn't trip.
        open(goi_for_tree, 'w').close()
    if tree_extras_total:
        existing_ids = set()
        if os.path.exists(goi_for_tree):
            with open(goi_for_tree) as fh:
                for line in fh:
                    if line.startswith('>'):
                        existing_ids.add(line[1:].split()[0].strip())
        appended = 0
        with open(goi_for_tree, 'a') as fh:
            for g in tree_extras_total:
                gid = g.get('id', '')
                seq = g.get('seq', '')
                if not gid or not seq or gid in existing_ids:
                    continue
                fh.write(f">{gid}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i+80] + "\n")
                existing_ids.add(gid)
                appended += 1
        logger.info(
            f"  Tree-input FASTA: {goi_for_tree} (added {appended} tree-only hit(s) "
            f"on top of expanded_db.faa)"
        )

    # Final summary
    logger.info("")
    logger.info("═══ Iterative Search Summary ═══")
    logger.info(f"  Genomes searched:  {total_genomes}")
    logger.info(f"  Genomes with hits: {genomes_with_hits}")
    logger.info(f"  Genomes no hits:   {genomes_without_hits}")
    logger.info(f"  Total new genes:   {total_new_genes}")
    logger.info(f"  Final database:    {expanded_db}")
    logger.info("════════════════════════════════")

if __name__ == "__main__":
    main()
