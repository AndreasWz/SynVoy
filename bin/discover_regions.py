#!/usr/bin/env python3
"""
discover_regions.py — structural/PLM ORF *discovery* helpers for SynVoy.

Counterpart to the re-rank path: re-ranking only scores GOI candidate models
that sequence search already found. DISCOVERY predicts every ORF in a syntenic
neighbourhood, folds/embeds them, and Foldseek/cosine-matches against the GOI
template — so it can surface an ortholog that diverged past all sequence
detection (<~25% id, "twilight zone") but kept its 3D fold. This is what Run 2's
divergence ladder (mosquito -> Nasonia -> Gryllus) needs: baseline BLAST misses,
structure rescues.

This module holds only the pure, testable orchestration around the GPU engine
(gpu_augment.py does the actual embed/fold/compare):
  * derive syntenic search windows from the flanking genes in a region GFF,
  * predict ORFs in each window (Augustus via gene_predictor) and map their
    coordinates back to the chromosome,
  * drop ORFs overlapping already-modelled genes (re-rank owns those),
  * materialise a surviving structural hit as a GOI gene-model GFF + protein.

Conservative by design: discovery models default to a HIGH structural-TM bar so
they cannot pollute a clean close-relative run with false positives.
"""

import logging
import os
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

try:
    from sequence_utils import parse_gff_attributes as _parse_gff_attributes
    from sequence_utils import load_genome
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_gff_attributes as _parse_gff_attributes
    from sequence_utils import load_genome

try:
    import parasail  # optional; only needed for the discovery sequence guards
except Exception:  # pragma: no cover - parasail present in synvoy_env
    parasail = None

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Region GFF -> intervals
# ---------------------------------------------------------------------------

def _is_goi_attrs(attrs: Dict[str, str]) -> bool:
    return attrs.get("SynVoyRole") == "goi" or bool(attrs.get("GOIClass"))


def read_region_features(gff_path: str) -> Tuple[List[Tuple[str, int, int]],
                                                 List[Tuple[str, int, int]]]:
    """Return (flanking_intervals, existing_model_intervals) as (chrom, start,
    end) 0-based half-open, parsed from a region GFF's mRNA/gene features.

    flanking = non-GOI models (the syntenic anchors that define where to look);
    existing = ALL modelled spans (GOI + flanking) — discovery skips ORFs that
    overlap these, since they are already annotated.
    """
    flanking: List[Tuple[str, int, int]] = []
    existing: List[Tuple[str, int, int]] = []
    if not gff_path or not os.path.exists(gff_path):
        return flanking, existing
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("mRNA", "gene"):
                continue
            try:
                start = int(parts[3]) - 1   # GFF 1-based -> 0-based
                end = int(parts[4])
            except (ValueError, IndexError):
                continue
            chrom = parts[0]
            attrs = _parse_gff_attributes(parts[8])
            existing.append((chrom, start, end))
            if not _is_goi_attrs(attrs):
                flanking.append((chrom, start, end))
    return flanking, existing


_CONF_RANK = {"LOW": 1, "MEDIUM": 2, "HIGH": 3}


def seq_goi_chroms(gff_path: str, min_confidence: str = "MEDIUM") -> set:
    """Chromosomes that already carry a sequence-found GOI model at >=
    min_confidence. Used by guard 4b so discovery only rescues a genome on a
    chromosome where sequence search found *nothing* confident (or, when it did
    find the GOI, on that same chromosome)."""
    want = _CONF_RANK.get(min_confidence.upper(), 2)
    chroms: set = set()
    if not gff_path or not os.path.exists(gff_path):
        return chroms
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("mRNA", "gene"):
                continue
            attrs = _parse_gff_attributes(parts[8])
            if not _is_goi_attrs(attrs):
                continue
            conf = (attrs.get("Confidence") or "LOW").upper()
            if _CONF_RANK.get(conf, 0) >= want:
                chroms.add(parts[0])
    return chroms


def count_flanking_support(orf: Dict[str, Any],
                           flanking: List[Tuple[str, int, int]],
                           window_bp: int = 100_000) -> int:
    """Number of flanking anchors on the ORF's chromosome within window_bp of its
    span — a proxy for 'how syntenic is the block this ORF sits in'."""
    chrom = orf.get("chrom")
    s, e = orf.get("start", 0), orf.get("end", 0)
    n = 0
    for fc, fs, fe in flanking:
        if fc != chrom:
            continue
        if fe >= s - window_bp and fs <= e + window_bp:
            n += 1
    return n


def cluster_intervals(intervals: List[Tuple[str, int, int]],
                      max_gap: int = 100_000,
                      pad: int = 20_000,
                      chrom_lengths: Optional[Dict[str, int]] = None
                      ) -> List[Tuple[str, int, int]]:
    """Merge per-chromosome intervals within `max_gap`, pad each window by
    `pad`, and clamp to [0, chrom_len]. Returns windows (chrom, start, end)."""
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, s, e in intervals:
        by_chrom.setdefault(chrom, []).append((s, e))

    windows: List[Tuple[str, int, int]] = []
    for chrom in sorted(by_chrom):
        spans = sorted(by_chrom[chrom])
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e + max_gap:
                cur_e = max(cur_e, e)
            else:
                windows.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        windows.append((chrom, cur_s, cur_e))

    padded: List[Tuple[str, int, int]] = []
    for chrom, s, e in windows:
        ps = max(0, s - pad)
        pe = e + pad
        if chrom_lengths and chrom in chrom_lengths:
            pe = min(pe, chrom_lengths[chrom])
        padded.append((chrom, ps, pe))
    return padded


# ---------------------------------------------------------------------------
# ORF prediction in a window (coords mapped back to the chromosome)
# ---------------------------------------------------------------------------

def predict_region_orfs(genome_seqs: Dict[str, str],
                        window: Tuple[str, int, int],
                        work_dir: str,
                        predictor: str = "augustus",
                        augustus_species: str = "fly",
                        min_aa: int = 30) -> List[Dict[str, Any]]:
    """Predict ORFs in a chromosome window and return them with CHROMOSOME
    coordinates. Each ORF: {id, seq, chrom, start(0-based), end(excl), strand}.

    Window sequence is sliced from the already-loaded genome, so no re-read.
    """
    chrom, ws, we = window
    seq = genome_seqs.get(chrom)
    if not seq:
        return []
    sub = seq[ws:we]
    if len(sub) < 60:
        return []

    os.makedirs(work_dir, exist_ok=True)
    win_fa = os.path.join(work_dir, "window.fa")
    with open(win_fa, "w") as fh:
        fh.write(f">{chrom}\n")
        for i in range(0, len(sub), 80):
            fh.write(sub[i:i + 80] + "\n")

    try:
        from gene_predictor import predict_orfs
    except ImportError:
        logger.warning("gene_predictor unavailable; discovery ORF prediction skipped.")
        return []

    try:
        raw = predict_orfs(win_fa, work_dir, predictor=predictor,
                           augustus_species=augustus_species, min_aa=min_aa)
    except Exception as exc:  # predictor failures are non-fatal for discovery
        logger.warning("ORF prediction failed for %s:%d-%d: %s", chrom, ws, we, exc)
        return []

    out: List[Dict[str, Any]] = []
    for orf in raw:
        rel_s = int(orf.get("start", 0))
        rel_e = int(orf.get("end", len(orf.get("seq", "")) * 3))
        out.append({
            "id": orf.get("id", ""),
            "seq": orf.get("seq", ""),
            "chrom": chrom,
            "start": ws + rel_s,        # chromosome, 0-based
            "end": ws + rel_e,          # chromosome, exclusive
            "strand": orf.get("strand", "+"),
        })
    return out


def _overlap_frac(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    """Fraction of the SHORTER span covered by the overlap of a and b."""
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    if e <= s:
        return 0.0
    ov = e - s
    shorter = max(1, min(a[1] - a[0], b[1] - b[0]))
    return ov / shorter


def orf_overlaps_any(orf: Dict[str, Any],
                     intervals: List[Tuple[str, int, int]],
                     min_frac: float = 0.3) -> bool:
    """True if the ORF overlaps any interval on the same chromosome by >= min_frac."""
    span = (orf["start"], orf["end"])
    for chrom, s, e in intervals:
        if chrom != orf["chrom"]:
            continue
        if _overlap_frac(span, (s, e)) >= min_frac:
            return True
    return False


# ---------------------------------------------------------------------------
# Materialise a structural discovery hit as a GOI gene model
# ---------------------------------------------------------------------------

def confidence_from_scores(tm: Optional[float],
                           emb: Optional[float],
                           tm_high: float = 0.7,
                           tm_medium: float = 0.5,
                           emb_high: float = 0.85,
                           emb_medium: float = 0.7) -> str:
    """Confidence for a discovery model from its structural TM + embedding sim."""
    tm = tm or 0.0
    emb = emb or 0.0
    if tm >= tm_high or emb >= emb_high:
        return "HIGH"
    if tm >= tm_medium or emb >= emb_medium:
        return "MEDIUM"
    return "LOW"


def build_discovery_gff(orf: Dict[str, Any],
                        genome_name: str,
                        model_index: int,
                        tm: Optional[float],
                        emb: Optional[float],
                        confidence: str,
                        source_query: str = "GOI") -> Tuple[List[str], Tuple[str, str]]:
    """Build a GOI gene-model GFF (mRNA + exon lines) + (id, protein) for a
    structural discovery hit. Coordinates are chromosome 0-based half-open in the
    ORF dict; GFF is 1-based inclusive."""
    chrom = orf["chrom"]
    gstart = orf["start"] + 1
    gend = orf["end"]
    strand = orf.get("strand", "+")
    mid = f"{genome_name}_structdisc_{model_index}"

    parts = [f"StructuralSimilarity={tm:.3f}"] if tm is not None else []
    if emb is not None:
        parts.append(f"EmbeddingSimilarity={emb:.3f}")
    score_attrs = (";" + ";".join(parts)) if parts else ""

    mrna_attrs = (
        f"ID={mid};SynVoyRole=goi;GOIClass=probable_goi;Confidence={confidence};"
        f"EvidenceType=structural_discovery;"
        f"InferenceReason=structural_discovery_no_sequence_hit;"
        f"SourceQuery={source_query}{score_attrs}"
    )
    mrna = "\t".join([chrom, "SynVoy", "mRNA", str(gstart), str(gend), ".",
                      strand, ".", mrna_attrs])
    exon = "\t".join([chrom, "SynVoy", "exon", str(gstart), str(gend), ".",
                      strand, ".", f"ID={mid}.exon1;Parent={mid}"])
    return [mrna, exon], (mid, orf.get("seq", ""))


# ---------------------------------------------------------------------------
# Discovery guards — reject promiscuous-fold false positives before a
# structural hit is ever materialised as a GOI model.
#
# Motivation: the oskar OSK domain is a GDSL/SGNH-lipase fold; genomes are full
# of lipase-domain proteins that fold like it. The first discovery-enabled run
# (LRZ 5702387) materialised a 2187-aa lipase (query oskar = 606 aa) as a
# HIGH-confidence "oskar" purely because ONE sub-domain window folded to TM 0.80
# — 5.8% query coverage, ~zero sequence homology, on a different scaffold than
# the real syntenic oskar. A TM>=0.7 gate alone cannot tell an ortholog from a
# fold analogue. These guards add the sequence + synteny sanity that the fold
# score lacks. All are pure/testable; parasail (already a synvoy_env dep) powers
# the alignment guards and degrades to a no-op only if it is unavailable.
# ---------------------------------------------------------------------------

def passes_length_ratio(orf_len: int, goi_len: int,
                        min_ratio: float = 0.5,
                        max_ratio: float = 2.0) -> bool:
    """Guard 1 — length sanity. True iff orf_len is within
    [min_ratio, max_ratio] x goi_len. A real ortholog is length-comparable to
    the query; a fold analogue (a domain inside a much larger protein) is not.
    2187 vs 606 -> ratio 3.6 -> rejected."""
    if goi_len <= 0:
        return True
    ratio = orf_len / float(goi_len)
    return min_ratio <= ratio <= max_ratio


def _sw_align(a: str, b: str) -> Tuple[float, int, int]:
    """Smith-Waterman of a (query) vs b. Returns (score, aligned_query_residues,
    matches). aligned_query_residues counts query-consuming aligned columns
    (M/=/X); matches counts identical columns when the cigar exposes them."""
    if parasail is None:
        raise RuntimeError("parasail is required for the discovery sequence guards")
    res = parasail.sw_trace_striped_16(a, b, 11, 1, parasail.blosum62)
    score = float(res.score)
    cigar = getattr(res, "cigar", None)
    ops = getattr(cigar, "seq", None) if cigar is not None else None
    if ops is None or len(ops) == 0:
        return score, 0, 0
    q_aligned = 0
    matches = 0
    for op in ops:
        length = op >> 4
        code = op & 0xF
        if code in (0, 7, 8):     # M / = / X all consume the query
            q_aligned += length
            if code == 7:
                matches += length
    return score, q_aligned, matches


def best_goi_alignment(orf_seq: str,
                       goi_templates: List[Tuple[str, str]]
                       ) -> Tuple[float, float, str, float]:
    """Guard 2 support — align orf_seq against every GOI template and return the
    best (coverage, identity_pct, goi_id, score). Coverage is the fraction of the
    GOI template spanned by the alignment (aligned GOI residues / len(GOI)); it
    separates a diverged-but-real ortholog (high coverage, low identity — the
    twilight zone) from a shared-domain fold analogue (tiny coverage)."""
    best_cov, best_id, best_gid, best_score = 0.0, 0.0, "", 0.0
    if not orf_seq:
        return best_cov, best_id, best_gid, best_score
    for gid, gseq in goi_templates:
        if not gseq:
            continue
        # GOI as the query so aligned columns measure coverage OF the GOI.
        score, q_aligned, matches = _sw_align(gseq, orf_seq)
        cov = q_aligned / float(len(gseq)) if gseq else 0.0
        ident = (100.0 * matches / q_aligned) if q_aligned else 0.0
        if cov > best_cov or (cov == best_cov and score > best_score):
            best_cov, best_id, best_gid, best_score = cov, ident, gid, score
    return best_cov, best_id, best_gid, best_score


def best_panel_score(orf_seq: str,
                     panel: Optional[List[Tuple[str, str]]]) -> Tuple[float, str]:
    """Guard 3 support — best SW score of orf_seq against a home-paralog panel
    (id, seq). Returns (score, name); (0.0, '') when the panel is empty/None."""
    best_score, best_name = 0.0, ""
    if not panel:
        return best_score, best_name
    for name, seq in panel:
        if not seq:
            continue
        score, _q, _m = _sw_align(orf_seq, seq)
        if score > best_score:
            best_score, best_name = score, name
    return best_score, best_name


def validate_discovery_hit(
    orf: Dict[str, Any],
    goi_templates: List[Tuple[str, str]],
    *,
    panel: Optional[List[Tuple[str, str]]] = None,
    flanking_support: int = 0,
    genome_has_seq_goi: bool = False,
    seq_goi_chroms: Optional[set] = None,
    min_len_ratio: float = 0.5,
    max_len_ratio: float = 2.0,
    min_goi_coverage: float = 0.5,
    paralog_margin: float = 1.0,
    min_block_flanking: int = 2,
    require_same_chrom_as_seq_goi: bool = True,
) -> Tuple[bool, str, Dict[str, Any]]:
    """Apply the four discovery guards to a candidate structural hit.

    Returns (ok, reason, metrics). `ok=False` means the hit is a likely
    false positive and must not be materialised; `reason` is a short slug for
    logging/diagnostics. `metrics` carries the computed coverage/identity/score
    so they can be recorded on the surviving model.

    Guards, in order (cheapest first):
      1. length ratio — orf length within [min,max] x GOI length.
      2. GOI coverage — SW alignment spans >= min_goi_coverage of the GOI.
      3. paralog discrimination — GOI SW score must beat the best home-paralog
         score by `paralog_margin` (no-op when no panel is supplied).
      4. syntenic primacy — the window must carry >= min_block_flanking flanking
         anchors, and (when the genome already has a sequence-found GOI) the hit
         must sit on a chromosome that carries one (a rescue is only needed where
         sequence search found nothing).
    """
    metrics: Dict[str, Any] = {}
    orf_seq = orf.get("seq", "") or ""
    orf_len = len(orf_seq)
    goi_len = max((len(s) for _i, s in goi_templates), default=0)
    metrics["orf_len"] = orf_len
    metrics["goi_len"] = goi_len

    # Guard 4a — synteny anchor floor (cheap, no alignment).
    metrics["flanking_support"] = flanking_support
    if flanking_support < min_block_flanking:
        return False, "weak_synteny", metrics

    # Guard 4b — prefer the primary block: if sequence search already placed the
    # GOI in this genome, only rescue on a chromosome that carries it.
    seq_goi_chroms = seq_goi_chroms or set()
    if (require_same_chrom_as_seq_goi and genome_has_seq_goi
            and orf.get("chrom") not in seq_goi_chroms):
        return False, "off_primary_block", metrics

    # Guard 1 — length ratio.
    if not passes_length_ratio(orf_len, goi_len, min_len_ratio, max_len_ratio):
        return False, "length_ratio", metrics

    # Guards 2/3 need parasail; if it's missing, keep the fold-only decision but
    # flag it so the caller can decide (default: allow, but record).
    if parasail is None:
        metrics["seq_guards"] = "parasail_unavailable"
        return True, "ok_no_seq_guards", metrics

    cov, ident, gid, goi_score = best_goi_alignment(orf_seq, goi_templates)
    metrics.update({"goi_coverage": round(cov, 4),
                    "goi_identity": round(ident, 2),
                    "goi_ref": gid, "goi_sw_score": goi_score})

    # Guard 2 — coverage floor.
    if cov < min_goi_coverage:
        return False, "low_goi_coverage", metrics

    # Guard 3 — paralog discrimination (only when a panel is supplied).
    panel_score, panel_name = best_panel_score(orf_seq, panel)
    metrics["panel_score"] = panel_score
    metrics["panel_best"] = panel_name
    if panel and panel_score > 0 and goi_score < paralog_margin * panel_score:
        return False, "beaten_by_paralog", metrics

    return True, "ok", metrics
