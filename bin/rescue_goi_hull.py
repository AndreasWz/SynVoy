#!/usr/bin/env python3
"""GOI synteny-hull rescue (docs/TODO.md §1m fumble fix).

The fumble: SynVoy only models the GOI inside flanking-SEEDED sub-blocks. When the
GOI sits in a *gap* between flanking blocks because its immediate neighbours are
rearranged in the target, no search window ever covers it and the true ortholog is
missed — even though it is a near-perfect match. Real case: human decorin's flanking
brackets cow chr5 17.8-28.4 Mb, but cow decorin at 21.0 Mb falls between the seeded
sub-blocks, so the run recovered only paralog fragments (max 55 %) while miniprot finds
the true 96 %-identity decorin there in 2 s.

This pass takes the flanking neighbourhood (the dominant cluster of HIGH-confidence
flanking genes per target chromosome — the syntenic HULL) and runs miniprot of the GOI
query across that whole window, so a GOI sitting in a flanking gap is still modelled.
It only fires when the hull has strong flanking support and no HIGH-confidence GOI model
already exists inside it (so it never duplicates a good call). Emits standard SynVoy GFF
rows tagged ``EvidenceType=synteny_hull_rescue`` with confidence set by miniprot identity.

Reuses the miniprot/parse helpers from rescue_strong_synteny.py.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from sequence_utils import parse_gff, parse_fasta, translate, reverse_complement  # noqa: E402
from rescue_strong_synteny import (  # noqa: E402
    _read_query,
    _run_miniprot_relaxed,
    _parse_miniprot_gff,
    _get_attr,
    _identity_from_mrna_attrs,
)


def _interval_overlaps(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 <= b1 and b0 <= a1


def _identity_pct(mrna_attrs: str) -> float:
    """miniprot --gff emits Identity as a FRACTION (0.962); return it as a percent."""
    v = _identity_from_mrna_attrs(mrna_attrs)
    return v * 100.0 if v <= 1.0 else v


def collect_features(gff_path: str):
    """Return (flanking_by_chrom, high_goi_by_chrom) from a SynVoy region GFF.

    flanking: gene features with SynVoyRole=flanking AND Confidence=HIGH.
    high_goi: gene features with SynVoyRole=goi AND Confidence=HIGH (already-good models).
    """
    flanking: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    high_goi: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    # SynVoy emits flanking/GOI as either `gene` or `mRNA` depending on EvidenceType
    # (flanking_miniprot → mRNA only; some GOI → gene). Read both and dedup by base ID
    # so a gene+its mRNA aren't double-counted.
    seen = set()
    for feat in parse_gff(gff_path, feature_types=["gene", "mRNA"]):
        attrs = feat["attributes"]
        role = attrs.get("SynVoyRole", "")
        conf = attrs.get("Confidence", "")
        fid = (attrs.get("ID") or "").replace(".mRNA", "")
        key = (feat["seqid"], fid, role)
        if key in seen:
            continue
        seen.add(key)
        s, e = feat["start"], feat["end"]
        if role == "flanking" and conf == "HIGH":
            flanking[feat["seqid"]].append((s, e))
        elif role == "goi" and conf == "HIGH":
            high_goi[feat["seqid"]].append((s, e))
    return flanking, high_goi


def dominant_cluster(intervals: List[Tuple[int, int]], max_gap: int) -> List[Tuple[int, int]]:
    """Cluster intervals by genomic gap; return the cluster with the most members."""
    if not intervals:
        return []
    ivs = sorted(intervals)
    clusters: List[List[Tuple[int, int]]] = [[ivs[0]]]
    for s, e in ivs[1:]:
        if s - clusters[-1][-1][1] <= max_gap:
            clusters[-1].append((s, e))
        else:
            clusters.append([(s, e)])
    return max(clusters, key=len)


def compute_hull(intervals: List[Tuple[int, int]], cluster_max_gap: int,
                 max_window: int) -> Tuple[int, int, int]:
    """Syntenic hull bracketing the GOI: the span of all flanking genes within
    ``max_window`` of the dominant flanking cluster's centre.

    This is the key to catching a GOI in a flanking GAP: the dominant cluster alone
    (e.g. cow chr5 22.5-24 Mb) can sit to one side of the GOI (decorin at 21.0 Mb),
    so we widen to the whole neighbourhood (17.8-28.4 Mb) while still excluding a far
    rearranged/spurious flanking hit (the 105 Mb outlier, ~80 Mb away).

    Returns (hull_start, hull_end, n_flanking_in_hull).
    """
    dom = dominant_cluster(intervals, cluster_max_gap)
    if not dom:
        return 0, 0, 0
    centre = (min(s for s, _e in dom) + max(e for _s, e in dom)) // 2
    lo, hi = centre - max_window // 2, centre + max_window // 2
    in_hull = [(s, e) for s, e in intervals if lo <= (s + e) // 2 <= hi]
    if not in_hull:
        in_hull = dom
    return min(s for s, _e in in_hull), max(e for _s, e in in_hull), len(in_hull)


def _extract_window(genome_path: str, chrom: str, start: int, end: int,
                    samtools_bin: str = "samtools") -> str:
    """genome[chrom:start..end] (1-based inclusive). samtools faidx if available,
    else a streaming FASTA fallback (no external dependency, works in tests)."""
    import subprocess
    try:
        res = subprocess.run([samtools_bin, "faidx", genome_path, f"{chrom}:{start}-{end}"],
                             check=True, capture_output=True, text=True)
        return "".join(l for l in res.stdout.splitlines() if not l.startswith(">"))
    except (FileNotFoundError, subprocess.CalledProcessError):
        for _h, cid, seq in parse_fasta(genome_path):
            if cid == chrom:
                return seq[max(0, start - 1):end]
    return ""


def _confidence_for(identity: float, high: float, medium: float) -> str:
    if identity >= high:
        return "HIGH"
    if identity >= medium:
        return "MEDIUM"
    return "LOW"


def _build_hull_gff_rows(mrna, cds_rows, chrom, window_start, query_id, parent_id,
                         confidence, query_coverage=None) -> List[str]:
    """Map miniprot window-relative coords to genome coords; SynVoy GFF block tagged
    EvidenceType=synteny_hull_rescue, confidence by identity."""
    mp_start, mp_end, strand = int(mrna[3]), int(mrna[4]), mrna[6]
    identity = _identity_pct(mrna[8])
    gstart = window_start + mp_start - 1
    gend = window_start + mp_end - 1
    # QueryCoverage is emitted so the report can tell "coverage is low" apart from
    # "coverage was never recorded". Without it every rescue model looked like the
    # latter, and a coverage-based demotion could not fire without also demoting
    # good full-length rescues.
    cov_attr = "" if query_coverage is None else f";QueryCoverage={query_coverage:.3f}"
    common = (f"Name={query_id};SynVoyRole=goi;Confidence={confidence};"
              f"EvidenceType=synteny_hull_rescue;GOIClass=synteny_hull_rescue;"
              f"ModelStatus=rescue;Identity={identity:.1f};TargetGene={query_id}"
              f"{cov_attr}")
    rows = [
        "\t".join([chrom, "SynVoy_hull", "gene", str(gstart), str(gend), ".", strand, ".",
                   f"ID={parent_id};{common}"]),
        "\t".join([chrom, "SynVoy_hull", "mRNA", str(gstart), str(gend), ".", strand, ".",
                   f"ID={parent_id}.mRNA;Parent={parent_id};{common}"]),
    ]
    for i, cds in enumerate(cds_rows, start=1):
        cs = window_start + int(cds[3]) - 1
        ce = window_start + int(cds[4]) - 1
        rows.append("\t".join([chrom, "SynVoy_hull", "CDS", str(cs), str(ce), ".",
                               strand, cds[7],
                               f"ID={parent_id}.CDS.{i};Parent={parent_id}.mRNA;"
                               f"SynVoyRole=goi;Confidence={confidence};"
                               f"EvidenceType=synteny_hull_rescue"]))
    return rows


def _translate_model(window_seq: str, cds_rows) -> str:
    """Protein for a miniprot model, from CDS rows in WINDOW coordinates.

    Needed because the rescue passes emit only a GFF, but the §1m ownership check
    RBH-aligns a protein FASTA. Without this a rescue model can never be checked
    against the home-paralog panel — the guard that exists to catch a mislabelled
    GOI never sees the models most likely to carry one.

    CDS rows are concatenated in transcript order (reversed on the minus strand,
    where miniprot still reports ascending coordinates), the first row's phase is
    trimmed, and any trailing stop is dropped.
    """
    if not cds_rows:
        return ""
    strand = cds_rows[0][6]
    rows = sorted(cds_rows, key=lambda c: int(c[3]), reverse=(strand == "-"))
    pieces = []
    for c in rows:
        s, e = int(c[3]), int(c[4])
        seg = window_seq[max(0, s - 1):e]
        if strand == "-":
            seg = reverse_complement(seg)
        pieces.append(seg)
    dna = "".join(pieces)
    try:
        phase = int(rows[0][7])
    except (ValueError, TypeError, IndexError):
        phase = 0
    dna = dna[phase:]
    dna = dna[:len(dna) - (len(dna) % 3)]
    return translate(dna).rstrip("*").replace("*", "X")


def _cds_coverage(cds_rows, query_len: int) -> float:
    """Fraction of the query covered by the model's CDS (aa)."""
    if query_len <= 0:
        return 0.0
    aa = sum((int(c[4]) - int(c[3]) + 1) for c in cds_rows) / 3.0
    return min(1.0, aa / query_len)


def main() -> int:
    ap = argparse.ArgumentParser(description="GOI synteny-hull rescue (fumble fix)")
    ap.add_argument("--region_gff", required=True, help="Per-genome SynVoy region GFF (flanking + goi)")
    ap.add_argument("--target_genome", required=True, help="Target genome FASTA")
    ap.add_argument("--query", required=True, help="GOI query protein FASTA")
    ap.add_argument("--genome_name", required=True, help="Logical genome name for ID stems")
    ap.add_argument("--output", required=True, help="Output GFF")
    ap.add_argument("--output_faa",
                    help="Optional protein FASTA of the rescued model(s), so the "
                         "§1m locus-ownership check can RBH them against the home "
                         "paralog panel. Written even when empty.")
    ap.add_argument("--samtools_bin", default="samtools")
    ap.add_argument("--min_flanking", type=int, default=4,
                    help="Min HIGH flanking genes in the dominant cluster to attempt a hull rescue")
    ap.add_argument("--cluster_max_gap", type=int, default=3000000,
                    help="bp gap that splits flanking into clusters (the hull = densest cluster)")
    ap.add_argument("--window_pad", type=int, default=100000, help="bp pad around the hull window")
    ap.add_argument("--max_window", type=int, default=20000000,
                    help="Skip hull windows larger than this (bp) — guards against runaway spans")
    ap.add_argument("--min_identity", type=float, default=40.0,
                    help="Min miniprot %identity to emit a rescued model")
    ap.add_argument("--min_coverage", type=float, default=0.5,
                    help="Min query coverage (CDS aa / query len) to emit a model")
    ap.add_argument("--classify_high_min_identity", type=float, default=70.0)
    ap.add_argument("--classify_medium_min_identity", type=float, default=45.0)
    ap.add_argument("--miniprot_outc", type=float, default=0.3)
    ap.add_argument("--miniprot_timeout", type=int, default=600)
    args = ap.parse_args()

    import re
    flanking, high_goi = collect_features(args.region_gff)
    query_id, query_seq = _read_query(args.query)
    # Sanitise for GFF IDs (UniProt headers carry pipes, e.g. sp|P07585|PGS2_HUMAN).
    query_id = re.sub(r"[^A-Za-z0-9_.-]", "_", query_id) or "GOI"
    qlen = len(query_seq)
    out_lines = ["##gff-version 3",
                 f"# synteny_hull_rescue for {args.genome_name} vs {query_id}"]
    faa_records = []  # (header, protein) for --output_faa

    if not query_seq:
        with open(args.output, "w") as fh:
            fh.write("\n".join(out_lines) + "\n")
        if args.output_faa:
            open(args.output_faa, "w").close()
        print(f"[hull] empty query {args.query}", file=sys.stderr)
        return 0

    import tempfile
    n_emitted = 0
    with tempfile.TemporaryDirectory(prefix="synvoy_hull_") as tmp:
        query_fa = os.path.join(tmp, "query.faa")
        with open(query_fa, "w") as fh:
            fh.write(f">{query_id}\n{query_seq}\n")

        for chrom, ivs in flanking.items():
            hull_s, hull_e, n_in_hull = compute_hull(ivs, args.cluster_max_gap, args.max_window)
            if n_in_hull < args.min_flanking:
                continue
            # Skip if a HIGH GOI model already sits inside the hull (already recovered).
            if any(_interval_overlaps(hull_s, hull_e, gs, ge)
                   for gs, ge in high_goi.get(chrom, [])):
                continue
            ws = max(1, hull_s - args.window_pad)
            we = hull_e + args.window_pad
            if we - ws > args.max_window:
                print(f"[hull] {chrom}: window {(we-ws)/1e6:.1f} Mb > max_window — skipped",
                      file=sys.stderr)
                continue
            window = _extract_window(args.target_genome, chrom, ws, we, args.samtools_bin)
            if not window:
                continue
            win_fa = os.path.join(tmp, "window.fa")
            with open(win_fa, "w") as fh:
                fh.write(f">{chrom}\n{window}\n")
            gff_text = _run_miniprot_relaxed(win_fa, query_fa, args.miniprot_outc,
                                             args.miniprot_timeout)
            # Best model by identity that clears coverage + identity gates.
            best = None
            for mrna, cds_rows in _parse_miniprot_gff(gff_text):
                ident = _identity_pct(mrna[8])
                cov = _cds_coverage(cds_rows, qlen)
                if ident >= args.min_identity and cov >= args.min_coverage:
                    if best is None or ident > best[0]:
                        best = (ident, mrna, cds_rows)
            if best is None:
                continue
            ident, mrna, cds_rows = best
            conf = _confidence_for(ident, args.classify_high_min_identity,
                                   args.classify_medium_min_identity)
            parent_id = f"GOI_{query_id}|{args.genome_name}_{chrom}_hull_rescue"
            best_cov = _cds_coverage(cds_rows, qlen)
            out_lines.extend(_build_hull_gff_rows(mrna, cds_rows, chrom, ws,
                                                  query_id, parent_id, conf,
                                                  query_coverage=best_cov))
            prot = _translate_model(window, cds_rows)
            if prot:
                faa_records.append((f"{parent_id}.mRNA", prot))
            n_emitted += 1
            print(f"[hull] {args.genome_name} {chrom}:{hull_s}-{hull_e} "
                  f"({n_in_hull} HIGH flanking, no HIGH GOI) -> rescued model "
                  f"id={ident:.1f}% cov={best_cov:.2f} conf={conf}", file=sys.stderr)

    with open(args.output, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    if args.output_faa:
        # Always written (possibly empty) so the Nextflow output glob never misses.
        with open(args.output_faa, "w") as fh:
            for header, prot in faa_records:
                fh.write(f">{header}\n{prot}\n")
    print(f"[hull] {args.genome_name}: {n_emitted} model(s) rescued.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
