#!/usr/bin/env python3
"""
Run Ly6 ground-truth checks against available SynVoy outputs.

Behavior:
  - Skips *_extraction.gff ground-truth files
  - Matches GT files to predicted GFFs by contig token in the filename
  - Uses CDS features when available for GT GOI; falls back to exons otherwise
  - Produces JSON + Markdown summary in ground_truth/ly6_3ftx
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple


GT_DIR = Path("ground_truth/ly6_3ftx/Ly6")
OUT_JSON = Path("ground_truth/ly6_3ftx/ly6_validation_report.json")
OUT_MD = Path("ground_truth/ly6_3ftx/ly6_validation_report.md")

GOI_REGEX = r"LY6|LYPD|LYNX|SLURP|PSCA|3FTx"


def parse_pred_seqids(path: Path) -> List[str]:
    seqids: List[str] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            seqid = line.split("\t", 1)[0]
            if seqid not in seqids:
                seqids.append(seqid)
    return seqids


def extract_tokens(gt_file: Path) -> List[str]:
    base = gt_file.stem
    tokens = [base]
    if "_" in base:
        prefix, rest = base.split("_", 1)
        if prefix.isalpha() and len(prefix) == 4:
            tokens.append(rest)
    return tokens


def choose_match(tokens: List[str], pred_seqids: Dict[Path, List[str]]) -> Optional[Tuple[Path, str, str]]:
    best = None
    best_score = -1
    for pred_file, seqids in pred_seqids.items():
        for token in tokens:
            for seqid in seqids:
                score = 0
                if seqid == token:
                    score = 3
                elif seqid.startswith(token):
                    score = 2
                elif token in seqid:
                    score = 1
                if score > best_score:
                    best_score = score
                    best = (pred_file, seqid, token)
    return best if best_score > 0 else None


def gt_goi_type(gt_file: Path, regex: re.Pattern) -> str:
    cds_hits = 0
    exon_hits = 0
    with gt_file.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = parts[8]
            if not regex.search(attrs):
                continue
            if ftype == "CDS":
                cds_hits += 1
            elif ftype == "exon":
                exon_hits += 1
    return "CDS" if cds_hits > 0 else "exon" if exon_hits > 0 else "CDS"


def run_validator(gt_file: Path, pred_file: Path, bed_file: Optional[Path], goi_type: str) -> dict:
    cmd = [
        "python",
        "scripts/validate_ly6_results.py",
        "--ground-truth",
        str(gt_file),
        "--pred-gff",
        str(pred_file),
        "--gt-goi-regex",
        GOI_REGEX,
        "--gt-goi-types",
        goi_type,
        "--normalize-seqids",
    ]
    if bed_file and bed_file.exists():
        cmd.extend(["--pred-locus-bed", str(bed_file)])
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Ly6 ground-truth validation")
    parser.add_argument(
        "--outdir",
        default="results/ly6_gt_v3",
        help="SynVoy output directory (default: results/ly6_gt_v3)",
    )
    args = parser.parse_args()
    pred_dir = Path(args.outdir) / "plot_inputs_synteny_block_locus_1"

    gt_files = sorted(GT_DIR.glob("*.gff"))
    pred_files = [p for p in sorted(pred_dir.glob("*.gff")) if p.name != "home_genome.gff"]
    pred_seqids = {p: parse_pred_seqids(p) for p in pred_files}

    results = []
    skipped = []
    unmatched = []

    regex = re.compile(GOI_REGEX)

    for gt_file in gt_files:
        if gt_file.name.endswith("_extraction.gff"):
            skipped.append(gt_file.name)
            continue
        tokens = extract_tokens(gt_file)
        match = choose_match(tokens, pred_seqids)
        if not match:
            unmatched.append(gt_file.name)
            continue
        pred_file, seqid, token = match
        bed_file = pred_file.with_suffix(".regions.bed")
        goi_type = gt_goi_type(gt_file, regex)
        report = run_validator(gt_file, pred_file, bed_file, goi_type)
        report["match"] = {
            "gt_file": gt_file.name,
            "pred_file": pred_file.name,
            "matched_token": token,
            "matched_seqid": seqid,
            "gt_goi_type_used": goi_type,
        }
        results.append(report)

    OUT_JSON.write_text(json.dumps({
        "goi_regex": GOI_REGEX,
        "skipped_extractions": skipped,
        "unmatched_gt_files": unmatched,
        "results": results,
    }, indent=2, sort_keys=True) + "\n")

    # Markdown summary
    lines = []
    lines.append("# Ly6 Ground-Truth Validation Summary")
    lines.append("")
    lines.append(f"GOI regex: `{GOI_REGEX}`")
    lines.append("")
    lines.append("Matches")
    lines.append("")
    lines.append("GT GFF\tPred GFF\tGOI Type\tCDS Tracking %\tFPR %\tLocus Start Δ\tLocus End Δ")
    lines.append("---\t---\t---\t---:\t---:\t---:\t---:")
    for r in results:
        m = r["match"]
        lb = r.get("locus_bounds", {})
        lines.append(
            f"{m['gt_file']}\t{m['pred_file']}\t{m['gt_goi_type_used']}\t"
            f"{r.get('cds_tracking_pct', 0):.3f}\t{r.get('pred_goi_false_positive_rate_pct', 0):.3f}\t"
            f"{lb.get('locus_start_delta_bp', 'NA')}\t{lb.get('locus_end_delta_bp', 'NA')}"
        )
    lines.append("")
    lines.append(f"Skipped extractions ({len(skipped)}): " + (", ".join(skipped) if skipped else "none"))
    lines.append(f"Unmatched GT files ({len(unmatched)}): " + (", ".join(unmatched) if unmatched else "none"))
    OUT_MD.write_text("\n".join(lines) + "\n")

    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
