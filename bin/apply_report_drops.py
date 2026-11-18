#!/usr/bin/env python3
"""Apply the report's decisions to the per-genome GFFs before the figures are drawn.

The report (generate_report.py) offers a weak match in the slot where the gene was expected
(Basis=gene_order_only) for curation only in a genome without a counted call of the gene;
elsewhere it drops the record and lists it (``goi_dedup.gene_order_only_dropped_records``). It
keeps an answer of the genome-wide path (GOIClass=genome_wide_strong_hit) only as the best one of
a genome without a counted call and without another record at its place, and lists the others
(``goi_dedup.genome_wide_strong_hits_dropped_records``). This script copies the per-genome GFFs
of one home locus into the plot-input folder without those rows: a GOI model of the rule's kind
(Basis=gene_order_only, or GOIClass=genome_wide_strong_hit) that overlaps a dropped record of the
same genome, and every row below it (Parent=). The published GFFs (plot_inputs_*/) and every
figure drawn from them then show what the report shows (review fix RF-3). Without a report
(the NO_REPORT sentinel, or a report from before these rules) the GFFs are copied unchanged.

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
from synvoy_label import GENOME_WIDE_CLASS  # noqa: E402

logger = setup_logging(name=__name__)


# The lists of dropped records in the report, and the GOI rows each one leaves out.
DROP_LISTS = (("gene_order_only_dropped_records", "gene_order_only"),
              ("genome_wide_strong_hits_dropped_records", GENOME_WIDE_CLASS))


def _is_dropped_kind(attrs, rule):
    """Whether a GOI row is of the kind a rule drops: Basis=gene_order_only, or an answer of the
    genome-wide path (GOIClass=genome_wide_strong_hit)."""
    if rule == GENOME_WIDE_CLASS:
        return attrs.get("GOIClass") == GENOME_WIDE_CLASS
    return attrs.get("Basis") == "gene_order_only"


def dropped_spans(report_path):
    """{genome: [(chrom, start, end, rule)]} of the records the report dropped: the
    gene_order_only records and the answers of the genome-wide path (DROP_LISTS)."""
    if not report_path or not os.path.isfile(report_path):
        return {}
    with open(report_path) as fh:
        try:
            report = json.load(fh)
        except ValueError:
            logger.warning(f"{report_path}: not a report (unparsable JSON); GFFs copied unchanged")
            return {}
    spans = defaultdict(list)
    dedup = report.get("goi_dedup") or {}
    for key, rule in DROP_LISTS:
        for r in dedup.get(key) or []:
            spans[r["genome"]].append((r["chrom"], int(r["start"]), int(r["end"]), rule))
    return dict(spans)


def _parents(attrs):
    return set(filter(None, (attrs.get("Parent") or "").split(",")))


def filter_gff(path, spans, out_path):
    """Write ``path`` to ``out_path`` without the GOI models that overlap one of ``spans``
    ([(chrom, start, end, rule)]; a span without its rule is a gene_order_only record) and are of
    that rule's kind (_is_dropped_kind), and without the rows below them. Returns the ids of the
    models left out."""
    with open(path) as fh:
        lines = fh.readlines()
    rows = []
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        rows.append(parse_gff_attributes(f[8]) if not ln.startswith("#") and len(f) >= 9 else None)
    drop = set()
    spans = [tuple(sp) if len(sp) > 3 else tuple(sp) + ("gene_order_only",) for sp in spans]
    for ln, a in zip(lines, rows):
        if a is None or (a.get("SynVoyRole") or "").lower() != "goi":
            continue
        f = ln.split("\t")
        s, e = int(f[3]), int(f[4])
        if any(c == f[0] and s <= ce and cs <= e and _is_dropped_kind(a, rule) for c, cs, ce, rule in spans):
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
            print(f"[report] {os.path.basename(path)}: left out {len(left_out)} model(s) the report "
                  f"dropped (gene-order-only beside a counted call, or a strong hit without synteny "
                  f"the report did not keep): {', '.join(left_out)}")


if __name__ == "__main__":
    main()
