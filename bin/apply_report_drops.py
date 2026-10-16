#!/usr/bin/env python3
"""Apply the report's decisions to the per-genome GFFs before the figures are drawn.

The report (generate_report.py) offers a weak match in the slot where the gene was expected
(Basis=gene_order_only) for curation only in a genome without a counted call of the gene;
elsewhere it drops the record and lists it (``goi_dedup.gene_order_only_dropped_records``).
This script copies the per-genome GFFs of one home locus into the plot-input folder without
those rows: a GOI model with Basis=gene_order_only that overlaps a dropped record of the same
genome, and every row below it (Parent=). The published GFFs (plot_inputs_*/) and every
figure drawn from them then show what the report shows (review fix RF-3). Without a report
(the NO_REPORT sentinel, or a report from before this rule) the GFFs are copied unchanged.

    apply_report_drops.py --report synvoy_report.json --out_dir plot_inputs_X a.fna.gff b.fna.gff
"""
import argparse
import json
import os
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_report import canonical_genome_id  # noqa: E402
from sequence_utils import parse_gff_attributes, setup_logging  # noqa: E402

logger = setup_logging(name=__name__)


def dropped_spans(report_path):
    """{genome: [(chrom, start, end)]} of the gene_order_only records the report dropped."""
    if not report_path or not os.path.isfile(report_path):
        return {}
    with open(report_path) as fh:
        try:
            report = json.load(fh)
        except ValueError:
            logger.warning(f"{report_path}: not a report (unparsable JSON); GFFs copied unchanged")
            return {}
    spans = defaultdict(list)
    for r in (report.get("goi_dedup") or {}).get("gene_order_only_dropped_records") or []:
        spans[r["genome"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    return dict(spans)


def _parents(attrs):
    return set(filter(None, (attrs.get("Parent") or "").split(",")))


def filter_gff(path, spans, out_path):
    """Write ``path`` to ``out_path`` without the GOI models with Basis=gene_order_only that
    overlap one of ``spans`` ([(chrom, start, end)]) and without the rows below them.
    Returns the ids of the models left out."""
    with open(path) as fh:
        lines = fh.readlines()
    rows = []
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        rows.append(parse_gff_attributes(f[8]) if not ln.startswith("#") and len(f) >= 9 else None)
    drop = set()
    for ln, a in zip(lines, rows):
        if a is None or a.get("Basis") != "gene_order_only" or (a.get("SynVoyRole") or "").lower() != "goi":
            continue
        f = ln.split("\t")
        s, e = int(f[3]), int(f[4])
        if any(c == f[0] and s <= ce and cs <= e for c, cs, ce in spans):
            drop.add(a.get("ID", ""))
    drop.discard("")
    models = sorted(drop)
    changed = True
    while changed:                       # rows below a dropped row, at any depth
        changed = False
        for a in rows:
            if a and a.get("ID") and a["ID"] not in drop and _parents(a) & drop:
                drop.add(a["ID"])
                changed = True
    with open(out_path, "w") as out:
        for ln, a in zip(lines, rows):
            if a is not None and (a.get("ID") in drop or _parents(a) & drop):
                continue
            out.write(ln)
    return models


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--report", default="", help="synvoy_report.json (or the NO_REPORT sentinel)")
    ap.add_argument("--out_dir", required=True, help="folder the GFFs are written to")
    ap.add_argument("gffs", nargs="*", help="per-genome GFFs of one home locus")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    spans = dropped_spans(a.report)
    for path in a.gffs:
        out_path = os.path.join(a.out_dir, os.path.basename(path))
        genome_spans = spans.get(canonical_genome_id(path), [])
        if not genome_spans:
            shutil.copyfile(path, out_path)
            continue
        left_out = filter_gff(path, genome_spans, out_path)
        if left_out:
            print(f"[report] {os.path.basename(path)}: left out {len(left_out)} gene-order-only "
                  f"model(s) the report dropped (the genome has a counted call): {', '.join(left_out)}")


if __name__ == "__main__":
    main()
