#!/usr/bin/env python3
"""
Validate SynVoy GOI predictions against a ground-truth GFF/GFF3.

This script is intentionally lightweight (no pandas) and focuses on:
  - CDS tracking % (ground-truth CDS bases covered by predicted GOI intervals)
  - False positive rate (predicted GOI intervals outside any GT gene region)
  - Optional locus boundary deltas (GT anchors vs predicted locus)

Example:
  python scripts/validate_ly6_results.py \
    --ground-truth ground_truth/ly6_3ftx/Ly6/Hosa_NC_000008.gff \
    --pred-gff results/ly6e_definitive_pro/plot_inputs_synteny_block_locus_1/GCF_000001635.27.fna.gff \
    --gt-goi-regex 'LY6E|PSCA|LYPD2' \
    --pred-locus-bed results/ly6e_definitive_pro/plot_inputs_synteny_block_locus_1/GCF_000001635.27.fna.regions.bed
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Feature:
    seqid: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    ftype: str
    attrs: Dict[str, str]


def parse_attributes(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_str.split(";"):
        if not part:
            continue
        if "=" in part:
            key, val = part.split("=", 1)
            attrs[key] = val
    return attrs


def _strip_version(seqid: str) -> str:
    """Strip NCBI version suffix (e.g. NW_018151385.1 -> NW_018151385)."""
    dot = seqid.rfind(".")
    if dot >= 0 and seqid[dot + 1:].isdigit():
        return seqid[:dot]
    return seqid


def parse_gff(path: Path, normalize_seqids: bool = False) -> List[Feature]:
    feats: List[Feature] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _source, ftype, start_s, end_s, _score, _strand, _phase, attrs_s = fields
            if normalize_seqids:
                seqid = _strip_version(seqid)
            try:
                start = int(start_s) - 1
                end = int(end_s)
            except ValueError:
                continue
            attrs = parse_attributes(attrs_s)
            feats.append(Feature(seqid=seqid, start=start, end=end, ftype=ftype, attrs=attrs))
    return feats


def parse_bed(path: Path) -> List[Tuple[str, int, int]]:
    intervals: List[Tuple[str, int, int]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            seqid = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            intervals.append((seqid, start, end))
    return intervals


def collect_intervals(features: Iterable[Feature]) -> Dict[str, List[Tuple[int, int]]]:
    by_seq: Dict[str, List[Tuple[int, int]]] = {}
    for feat in features:
        by_seq.setdefault(feat.seqid, []).append((feat.start, feat.end))
    return by_seq


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def total_bases(intervals: Dict[str, List[Tuple[int, int]]]) -> int:
    total = 0
    for seqid, spans in intervals.items():
        for start, end in spans:
            total += max(0, end - start)
    return total


def overlap_bases(a: Dict[str, List[Tuple[int, int]]],
                  b: Dict[str, List[Tuple[int, int]]]) -> int:
    total = 0
    for seqid, a_spans in a.items():
        if seqid not in b:
            continue
        b_spans = b[seqid]
        i = j = 0
        while i < len(a_spans) and j < len(b_spans):
            a_start, a_end = a_spans[i]
            b_start, b_end = b_spans[j]
            start = max(a_start, b_start)
            end = min(a_end, b_end)
            if end > start:
                total += end - start
            if a_end < b_end:
                i += 1
            else:
                j += 1
    return total


def interval_overlaps_any(interval: Tuple[int, int],
                          spans: List[Tuple[int, int]]) -> bool:
    start, end = interval
    for s, e in spans:
        if e <= start:
            continue
        if s >= end:
            return False
        return True
    return False


def bounds_from_intervals(intervals: Dict[str, List[Tuple[int, int]]]) -> Optional[Tuple[str, int, int]]:
    seqids = [s for s, spans in intervals.items() if spans]
    if not seqids:
        return None
    if len(seqids) > 1:
        # Multi-contig. Caller can decide how to handle.
        seqid = "MULTI"
        all_spans = [span for spans in intervals.values() for span in spans]
        start = min(s for s, _ in all_spans)
        end = max(e for _, e in all_spans)
        return (seqid, start, end)
    seqid = seqids[0]
    spans = intervals[seqid]
    start = min(s for s, _ in spans)
    end = max(e for _, e in spans)
    return (seqid, start, end)


def attr_matches(attrs: Dict[str, str], regex: re.Pattern) -> bool:
    for key in ("Name", "ID", "Parent", "SynTerra_Parent", "gene", "product"):
        val = attrs.get(key)
        if val and regex.search(val):
            return True
    return False


def is_pred_goi(feat: Feature, regex: Optional[re.Pattern]) -> bool:
    if feat.attrs.get("SynTerraRole", "").lower() == "goi":
        return True
    if regex and attr_matches(feat.attrs, regex):
        return True
    for key in ("Name", "ID", "Parent", "SynTerra_Parent"):
        val = feat.attrs.get(key, "")
        if val.startswith("GOI_") or val.startswith("GOI_copy"):
            return True
    return False


def is_gt_goi(feat: Feature, regex: Optional[re.Pattern]) -> bool:
    if regex is None:
        return True
    return attr_matches(feat.attrs, regex)


def parse_types(types_s: str) -> List[str]:
    return [t.strip() for t in types_s.split(",") if t.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate SynVoy GOI predictions against ground-truth GFF.")
    parser.add_argument("--ground-truth", required=True, help="Ground-truth GFF/GFF3 path.")
    parser.add_argument("--pred-gff", required=True, help="Predicted SynVoy GFF path.")
    parser.add_argument("--pred-locus-bed", help="Optional predicted locus BED to compute locus deltas.")
    parser.add_argument("--gt-goi-regex", help="Regex for GT GOI selection (matches Name/ID/Parent).")
    parser.add_argument("--pred-goi-regex", help="Regex for predicted GOI selection (matches Name/ID/Parent).")
    parser.add_argument("--gt-goi-types", default="CDS", help="Comma list of GT feature types to evaluate.")
    parser.add_argument("--pred-goi-types", default="mRNA,gene", help="Comma list of predicted GOI feature types.")
    parser.add_argument("--gt-gene-types", default="gene,mRNA,CDS",
                        help="Comma list of GT feature types considered 'annotated gene space'.")
    parser.add_argument("--gt-anchor-regex",
                        help="Regex for GT anchor genes to define locus bounds (else uses all GT genes).")
    parser.add_argument("--normalize-seqids", action="store_true",
                        help="Strip NCBI version suffixes (.1, .2, …) from seqids in both GFFs "
                             "so GT 'NW_018151385' matches pred 'NW_018151385.1'.")
    parser.add_argument("--output", help="Optional JSON output path.")
    args = parser.parse_args()

    gt_path = Path(args.ground_truth)
    pred_path = Path(args.pred_gff)
    pred_bed_path = Path(args.pred_locus_bed) if args.pred_locus_bed else None

    norm = getattr(args, "normalize_seqids", False)
    gt_feats = parse_gff(gt_path, normalize_seqids=norm)
    pred_feats = parse_gff(pred_path, normalize_seqids=norm)

    gt_goi_re = re.compile(args.gt_goi_regex) if args.gt_goi_regex else None
    pred_goi_re = re.compile(args.pred_goi_regex) if args.pred_goi_regex else None
    gt_anchor_re = re.compile(args.gt_anchor_regex) if args.gt_anchor_regex else None

    gt_goi_types = set(parse_types(args.gt_goi_types))
    pred_goi_types = set(parse_types(args.pred_goi_types))
    gt_gene_types = set(parse_types(args.gt_gene_types))

    gt_goi_feats = [f for f in gt_feats if f.ftype in gt_goi_types and is_gt_goi(f, gt_goi_re)]
    pred_goi_feats = [f for f in pred_feats if f.ftype in pred_goi_types and is_pred_goi(f, pred_goi_re)]
    gt_gene_feats = [f for f in gt_feats if f.ftype in gt_gene_types]

    gt_goi_intervals = {k: merge_intervals(v) for k, v in collect_intervals(gt_goi_feats).items()}
    pred_goi_intervals = {k: merge_intervals(v) for k, v in collect_intervals(pred_goi_feats).items()}
    gt_gene_intervals = {k: merge_intervals(v) for k, v in collect_intervals(gt_gene_feats).items()}

    gt_goi_bases = total_bases(gt_goi_intervals)
    covered_bases = overlap_bases(gt_goi_intervals, pred_goi_intervals)
    cds_tracking_pct = (covered_bases / gt_goi_bases * 100.0) if gt_goi_bases else 0.0

    # False positives: predicted GOI intervals not overlapping any GT gene region
    pred_total = 0
    pred_fp = 0
    for seqid, spans in pred_goi_intervals.items():
        pred_total += len(spans)
        gt_spans = gt_gene_intervals.get(seqid, [])
        for span in spans:
            if not interval_overlaps_any(span, gt_spans):
                pred_fp += 1
    fpr = (pred_fp / pred_total * 100.0) if pred_total else 0.0

    # Locus bounds
    if gt_anchor_re:
        gt_anchor_feats = [f for f in gt_feats if f.ftype == "gene" and attr_matches(f.attrs, gt_anchor_re)]
        gt_anchor_intervals = {k: merge_intervals(v) for k, v in collect_intervals(gt_anchor_feats).items()}
    else:
        gt_anchor_intervals = gt_gene_intervals

    gt_bounds = bounds_from_intervals(gt_anchor_intervals)

    if pred_bed_path and pred_bed_path.exists():
        bed_intervals = parse_bed(pred_bed_path)
        pred_locus_intervals: Dict[str, List[Tuple[int, int]]] = {}
        for seqid, start, end in bed_intervals:
            pred_locus_intervals.setdefault(seqid, []).append((start, end))
        pred_locus_intervals = {k: merge_intervals(v) for k, v in pred_locus_intervals.items()}
    else:
        pred_locus_intervals = pred_goi_intervals

    pred_bounds = bounds_from_intervals(pred_locus_intervals)

    locus_deltas = {}
    if gt_bounds and pred_bounds and gt_bounds[0] == pred_bounds[0]:
        locus_deltas = {
            "locus_seqid": gt_bounds[0],
            "gt_start": gt_bounds[1],
            "gt_end": gt_bounds[2],
            "pred_start": pred_bounds[1],
            "pred_end": pred_bounds[2],
            "locus_start_delta_bp": pred_bounds[1] - gt_bounds[1],
            "locus_end_delta_bp": pred_bounds[2] - gt_bounds[2],
            "locus_span_delta_bp": (pred_bounds[2] - pred_bounds[1]) - (gt_bounds[2] - gt_bounds[1]),
        }
    else:
        locus_deltas = {
            "locus_seqid": None,
            "gt_start": gt_bounds[1] if gt_bounds else None,
            "gt_end": gt_bounds[2] if gt_bounds else None,
            "pred_start": pred_bounds[1] if pred_bounds else None,
            "pred_end": pred_bounds[2] if pred_bounds else None,
            "note": "Bounds are on different contigs or unavailable. Provide --pred-locus-bed for clarity.",
        }

    report = {
        "ground_truth_gff": str(gt_path),
        "pred_gff": str(pred_path),
        "gt_goi_regex": args.gt_goi_regex,
        "pred_goi_regex": args.pred_goi_regex,
        "gt_goi_types": sorted(gt_goi_types),
        "pred_goi_types": sorted(pred_goi_types),
        "gt_gene_types": sorted(gt_gene_types),
        "gt_goi_cds_bases": gt_goi_bases,
        "covered_gt_goi_cds_bases": covered_bases,
        "cds_tracking_pct": round(cds_tracking_pct, 3),
        "pred_goi_intervals": pred_total,
        "pred_goi_false_positives": pred_fp,
        "pred_goi_false_positive_rate_pct": round(fpr, 3),
        "locus_bounds": locus_deltas,
    }

    out_json = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        Path(args.output).write_text(out_json + "\n")
    else:
        print(out_json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
