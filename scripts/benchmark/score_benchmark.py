#!/usr/bin/env python3
"""Score competitor tool outputs against the curated melittin ortholog truth set.

Each tool wrapper script must emit a normalized TSV at:
    benchmark_results/<tool>/calls.tsv

with columns (one row per target species):

    species               # Snake_case species name matching melittin_orthologs.tsv
    accession             # NCBI accession (or 'local' for the 5 Koludarov genomes)
    called_status         # one of: PRESENT, ABSENT, AMBIGUOUS
    locus_chrom           # chromosome / scaffold ID where ortholog was called (or '-')
    locus_start           # 1-based start (or '-')
    locus_end             # 1-based end (or '-')
    strand                # +, -, or '-'
    confidence            # tool-native confidence (free text, e.g. HIGH/MEDIUM/LOW or e-value)
    extra                 # tool-specific notes (free text)

This script joins each tool's calls.tsv against
tests/benchmark_truth/melittin_orthologs.tsv and emits:
    benchmark_results/master_table.tsv     # rows = species x tool
    benchmark_results/confusion_per_tool.tsv  # TP/FP/FN/TN/precision/recall/F1 per tool
    benchmark_results/confusion_per_tier.tsv  # same, broken out by tier

Scoring rules:
- truth=PRESENT and called=PRESENT  -> TP
- truth=PRESENT and called=ABSENT   -> FN
- truth=LOST    and called=ABSENT   -> TN
- truth=LOST    and called=PRESENT  -> FP
- truth=ABSENT  and called=ABSENT   -> TN
- truth=ABSENT  and called=PRESENT  -> FP
- truth=PSEUDOGENE -> partial credit: locus found at syntenic position with low
  confidence is the correct answer; counted in PSEUDOGENE_correct/PSEUDOGENE_wrong
  buckets, not in TP/FP/FN/TN.
- truth=PRESENT_NONSYNTENIC -> the gene exists per cDNA literature but synteny
  may not be conserved. Reported separately; not counted in confusion matrix.
- truth=UNCERTAIN -> not scored. Reported as method-agreement matrix only.
- called=AMBIGUOUS -> counted as ABSENT for confusion-matrix purposes; flagged
  in the master table for inspection.

Usage:
    python scripts/benchmark/score_benchmark.py \\
        --truth tests/benchmark_truth/melittin_orthologs.tsv \\
        --calls-glob 'benchmark_results/*/calls.tsv' \\
        --outdir benchmark_results/

Tool wrappers (synvoy, orthofinder, sonicparanoid, mmseqs2_rbh, genespace,
optionally toga2) live alongside this file as run_<tool>.sh and produce the
normalized calls.tsv described above.
"""

from __future__ import annotations

import argparse
import csv
import glob as _glob
import sys
from collections import defaultdict
from pathlib import Path


TRUTH_COLUMNS = ["species", "accession", "tier", "status", "confidence", "source",
                 "syntenic_locus", "notes"]
CALL_COLUMNS = ["species", "accession", "called_status", "locus_chrom",
                "locus_start", "locus_end", "strand", "confidence", "extra"]


def read_truth(path: Path) -> dict[str, dict]:
    rows = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(TRUTH_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"truth file {path} missing columns: {sorted(missing)}")
        for row in reader:
            rows[row["species"]] = row
    return rows


def read_calls(path: Path) -> dict[str, dict]:
    rows = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(CALL_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"calls file {path} missing columns: {sorted(missing)}")
        for row in reader:
            rows[row["species"]] = row
    return rows


def classify(truth_status: str, called_status: str) -> str:
    """Return one of: TP, FP, FN, TN, PSEUDOGENE_correct, PSEUDOGENE_wrong,
    NONSYNTENIC, UNCERTAIN, INVALID."""
    if called_status == "AMBIGUOUS":
        called_status = "ABSENT"
    if truth_status == "UNCERTAIN":
        return "UNCERTAIN"
    if truth_status == "PRESENT_NONSYNTENIC":
        return "NONSYNTENIC"
    if truth_status == "PSEUDOGENE":
        return "PSEUDOGENE_correct" if called_status == "ABSENT" else "PSEUDOGENE_wrong"
    if truth_status == "PRESENT":
        return "TP" if called_status == "PRESENT" else "FN"
    if truth_status in ("LOST", "ABSENT"):
        return "FP" if called_status == "PRESENT" else "TN"
    return "INVALID"


def confusion(rows: list[dict]) -> dict:
    """Aggregate TP/FP/FN/TN counts and precision/recall/F1."""
    counts = defaultdict(int)
    for r in rows:
        counts[r["bucket"]] += 1
    tp = counts["TP"]
    fp = counts["FP"]
    fn = counts["FN"]
    tn = counts["TN"]
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision + recall else float("nan"))
    return {
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "PSEUDOGENE_correct": counts["PSEUDOGENE_correct"],
        "PSEUDOGENE_wrong": counts["PSEUDOGENE_wrong"],
        "NONSYNTENIC": counts["NONSYNTENIC"],
        "UNCERTAIN": counts["UNCERTAIN"],
        "precision": f"{precision:.3f}",
        "recall": f"{recall:.3f}",
        "F1": f"{f1:.3f}",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--calls-glob", required=True,
                    help="Glob matching <toolname>/calls.tsv files")
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    truth = read_truth(args.truth)

    call_files = sorted(Path(p) for p in _glob.glob(args.calls_glob))
    if not call_files:
        sys.exit(f"no calls.tsv files matched glob: {args.calls_glob}")

    master_rows = []
    per_tool_rows = defaultdict(list)
    per_tool_tier_rows = defaultdict(lambda: defaultdict(list))

    for cf in call_files:
        tool = cf.parent.name
        calls = read_calls(cf)
        for species, t_row in truth.items():
            if t_row["tier"] == "home":
                continue

            c_row = calls.get(species)
            called = c_row["called_status"] if c_row else "ABSENT"
            bucket = classify(t_row["status"], called)
            entry = {
                "tool": tool,
                "species": species,
                "tier": t_row["tier"],
                "truth_status": t_row["status"],
                "truth_confidence": t_row["confidence"],
                "called_status": called,
                "bucket": bucket,
                "locus_chrom": c_row["locus_chrom"] if c_row else "-",
                "locus_start": c_row["locus_start"] if c_row else "-",
                "locus_end":   c_row["locus_end"]   if c_row else "-",
                "tool_confidence": c_row["confidence"] if c_row else "-",
                "tool_extra":      c_row["extra"]      if c_row else "-",
            }
            master_rows.append(entry)
            per_tool_rows[tool].append(entry)
            per_tool_tier_rows[tool][t_row["tier"]].append(entry)

    master_path = args.outdir / "master_table.tsv"
    with master_path.open("w") as fh:
        cols = ["tool", "species", "tier", "truth_status", "truth_confidence",
                "called_status", "bucket", "locus_chrom", "locus_start",
                "locus_end", "tool_confidence", "tool_extra"]
        fh.write("\t".join(cols) + "\n")
        for r in master_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    confusion_path = args.outdir / "confusion_per_tool.tsv"
    with confusion_path.open("w") as fh:
        cols = ["tool", "TP", "FP", "FN", "TN", "PSEUDOGENE_correct",
                "PSEUDOGENE_wrong", "NONSYNTENIC", "UNCERTAIN",
                "precision", "recall", "F1"]
        fh.write("\t".join(cols) + "\n")
        for tool, rows in per_tool_rows.items():
            cm = confusion(rows)
            fh.write(tool + "\t" + "\t".join(str(cm[c]) for c in cols[1:]) + "\n")

    tier_path = args.outdir / "confusion_per_tier.tsv"
    with tier_path.open("w") as fh:
        cols = ["tool", "tier", "TP", "FP", "FN", "TN",
                "precision", "recall", "F1"]
        fh.write("\t".join(cols) + "\n")
        for tool, by_tier in per_tool_tier_rows.items():
            for tier, rows in by_tier.items():
                cm = confusion(rows)
                fh.write(f"{tool}\t{tier}\t" +
                         "\t".join(str(cm[c]) for c in cols[2:]) + "\n")

    print(f"wrote {master_path}")
    print(f"wrote {confusion_path}")
    print(f"wrote {tier_path}")


if __name__ == "__main__":
    main()
