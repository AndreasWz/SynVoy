#!/usr/bin/env python3
"""Apply the report's decisions to the per-genome GFFs before the figures are drawn.

The report (generate_report.py) offers a weak match in the slot where the gene was expected
(Basis=gene_order_only) for curation only in a genome without a counted call of the gene;
elsewhere it drops the record and lists it (``goi_dedup.gene_order_only_dropped_records``). It
keeps an answer of the genome-wide path (GOIClass=genome_wide_strong_hit) only as the best one of
a genome without a counted call and without another record at its place, and lists the others
(``goi_dedup.genome_wide_strong_hits_dropped_records``). This script copies the per-genome GFFs
of one home locus into the plot-input folder without those rows, and every row below them
(Parent=). The published GFFs (plot_inputs_*/) and every figure drawn from them then show what
the report shows (review fix RF-3). Without a report (the NO_REPORT sentinel, or a report from
before these rules) the GFFs are copied unchanged.

With ``--locus`` and the report's ``figures.<locus>`` block (generate_report.figure_marks), the
rows left out are the models that block marks ``dropped``, by their ids: the members of the
dropped records that this locus's search modelled. Without the block, they are the GOI models of
a rule's kind (Basis=gene_order_only, or GOIClass=genome_wide_strong_hit) that overlap a dropped
record of the same genome; that rule also takes a kept record's model of the same kind, such as a
kept answer that overlaps an answer the report dropped at another home locus.

With ``--locus``, it also marks the GOI models of that home locus with the report's outcome:
``ReportOutcome`` (merged, overruled_piece, paralog_not_goi, identity_coverage_decoupled,
genome_wide_strong_hit, not_in_report) and ``ReportRecord`` (the span of the model's record) on
the model's row, ``ReportDrawnBy`` (the model that draws the record) on a merged model,
``AlsoFoundFrom`` (the other home loci) on the model that draws a record found from several, and
the header line ``##synvoy-report-marks <locus>`` after ``##gff-version``. No row is removed for
these marks: the figures leave out the models marked merged, so each record is drawn once, and
the header tells them that the report's marks are applied (a GFF without it, from a report
without the block, is drawn by the figures' own overlap rule).

    apply_report_drops.py --report synvoy_report.json --locus synteny_block_locus_1 \\
        --out_dir plot_inputs_X a.fna.gff b.fna.gff
"""
import argparse
import json
import os
import shutil
import sys
from collections import Counter, defaultdict
from urllib.parse import quote

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_report import canonical_genome_id  # noqa: E402
from sequence_utils import parse_gff_attributes, setup_logging  # noqa: E402
from synvoy_label import GENOME_WIDE_CLASS  # noqa: E402

logger = setup_logging(name=__name__)


# The lists of dropped records in the report, and the GOI rows each one leaves out.
DROP_LISTS = (("gene_order_only_dropped_records", "gene_order_only"),
              ("genome_wide_strong_hits_dropped_records", GENOME_WIDE_CLASS))
# The header line that says the report's marks are applied to a GFF, and the attributes of a mark.
MARKS_HEADER = "##synvoy-report-marks"
MARK_KEYS = ("ReportOutcome", "ReportRecord", "ReportDrawnBy", "AlsoFoundFrom")


def _load_report(report_path):
    """The report as a dict, or None (no file, the NO_REPORT sentinel, unparsable JSON)."""
    if not report_path or not os.path.isfile(report_path):
        return None
    with open(report_path) as fh:
        try:
            report = json.load(fh)
        except ValueError:
            logger.warning(f"{report_path}: not a report (unparsable JSON); GFFs copied unchanged")
            return None
    return report if isinstance(report, dict) else None


def _is_dropped_kind(attrs, rule):
    """Whether a GOI row is of the kind a rule drops: Basis=gene_order_only, or an answer of the
    genome-wide path (GOIClass=genome_wide_strong_hit)."""
    if rule == GENOME_WIDE_CLASS:
        return attrs.get("GOIClass") == GENOME_WIDE_CLASS
    return attrs.get("Basis") == "gene_order_only"


def dropped_spans(report_path):
    """{genome: [(chrom, start, end, rule)]} of the records the report dropped: the
    gene_order_only records and the answers of the genome-wide path (DROP_LISTS)."""
    return _dropped_spans_of(_load_report(report_path))


def _dropped_spans_of(report):
    spans = defaultdict(list)
    dedup = (report or {}).get("goi_dedup") or {}
    for key, rule in DROP_LISTS:
        for r in dedup.get(key) or []:
            spans[r["genome"]].append((r["chrom"], int(r["start"]), int(r["end"]), rule))
    return dict(spans)


def locus_of(name):
    """The home locus id of a plot-input folder's home BED name: synteny_block_locus_1(.bed) ->
    locus_1; a bare id is returned as it is."""
    name = os.path.basename(name or "")
    name = name[:-len(".bed")] if name.endswith(".bed") else name
    return name[len("synteny_block_"):] if name.startswith("synteny_block_") else name


def report_marks(report, locus):
    """{genome: {model id: {attribute: value}}} of the report's figure marks for ``locus``
    (``figures.<locus>``), or None when the report has no such block (no report, a report from
    before the block, no --locus): the GFFs then carry no marks."""
    block = ((report or {}).get("figures") or {}).get(locus) if locus else None
    if not isinstance(block, dict):
        return None
    marks = defaultdict(dict)
    for e in block.get("models") or ():
        marks[e["genome"]].setdefault(e["model"], {})["ReportOutcome"] = e["outcome"]
        if e.get("record"):
            marks[e["genome"]][e["model"]]["ReportRecord"] = e["record"]
        if e.get("drawn_by"):
            marks[e["genome"]][e["model"]]["ReportDrawnBy"] = e["drawn_by"]
    for e in block.get("also_found_from") or ():
        if e.get("loci"):
            marks[e["genome"]].setdefault(e["model"], {})["AlsoFoundFrom"] = ",".join(e["loci"])
    return dict(marks)


def _parents(attrs):
    return set(filter(None, (attrs.get("Parent") or "").split(",")))


def _gff_value(text):
    """An attribute value with the characters GFF3 reserves (; = & , tab, newline) escaped; the
    commas between the loci of AlsoFoundFrom stay (a list of values)."""
    return ",".join(quote(part, safe=" :|._-~/()") for part in str(text).split(","))


def _rows(lines):
    return [parse_gff_attributes(ln.rstrip("\n").split("\t")[8])
            if not ln.startswith("#") and len(ln.rstrip("\n").split("\t")) >= 9 else None
            for ln in lines]


def _drop_rows(lines, spans):
    """``lines`` without the GOI models that overlap one of ``spans`` ([(chrom, start, end, rule)];
    a span without its rule is a gene_order_only record) and are of that rule's kind
    (_is_dropped_kind), and without the rows below them; and the ids of the models left out."""
    rows = _rows(lines)
    drop = set()
    spans = [tuple(sp) if len(sp) > 3 else tuple(sp) + ("gene_order_only",) for sp in spans]
    for ln, a in zip(lines, rows):
        if a is None or (a.get("SynVoyRole") or "").lower() != "goi":
            continue
        f = ln.split("\t")
        s, e = int(f[3]), int(f[4])
        if any(c == f[0] and s <= ce and cs <= e and _is_dropped_kind(a, rule) for c, cs, ce, rule in spans):
            drop.add(a.get("ID", ""))
    return _without_models(lines, rows, drop)


def drop_models(lines, model_ids):
    """``lines`` without the rows of the models ``model_ids`` and the rows below them; and the ids
    of the models left out (those of ``model_ids`` the lines hold)."""
    rows, wanted = _rows(lines), set(model_ids)
    return _without_models(lines, rows, {a.get("ID", "") for a in rows
                                         if a is not None and a.get("ID") in wanted})


def _without_models(lines, rows, drop):
    drop = set(drop)
    drop.discard("")
    models = sorted(drop)
    changed = True
    while changed:                       # rows below a dropped row, at any depth
        changed = False
        for a in rows:
            if a and a.get("ID") and a["ID"] not in drop and _parents(a) & drop:
                drop.add(a["ID"])
                changed = True
    kept = [ln for ln, a in zip(lines, rows)
            if not (a is not None and (a.get("ID") in drop or _parents(a) & drop))]
    return kept, models


def filter_gff(path, spans, out_path):
    """Write ``path`` to ``out_path`` without the GOI models that overlap one of ``spans``
    ([(chrom, start, end, rule)]; a span without its rule is a gene_order_only record) and are of
    that rule's kind (_is_dropped_kind), and without the rows below them. Returns the ids of the
    models left out."""
    with open(path) as fh:
        lines = fh.readlines()
    kept, models = _drop_rows(lines, spans)
    with open(out_path, "w") as out:
        out.writelines(kept)
    return models


def mark_lines(lines, marks, locus):
    """``lines`` with the report's marks of one genome (``marks``: {model id: {attribute: value}})
    appended to the rows (mRNA or gene) of those models, and the header line
    ``##synvoy-report-marks <locus>`` after ``##gff-version`` (first when there is none). Marks a
    GFF already carries are replaced, so a GFF marked twice reads as marked once; every other row is
    kept as it is. Returns (lines, {model id: its marks} of the models marked)."""
    out, marked = [], {}
    for ln in lines:
        if ln.startswith(MARKS_HEADER):
            continue
        f = ln.rstrip("\n").split("\t")
        if ln.startswith("#") or len(f) < 9 or f[2] not in ("mRNA", "gene"):
            out.append(ln)
            continue
        attrs = parse_gff_attributes(f[8])
        new = marks.get(attrs.get("ID", ""))
        if not new and not any(k in attrs for k in MARK_KEYS):
            out.append(ln)
            continue
        parts = [p for p in f[8].split(";") if p.strip() and p.split("=", 1)[0].strip() not in MARK_KEYS]
        if new:
            parts += [f"{k}={_gff_value(new[k])}" for k in MARK_KEYS if new.get(k)]
            marked[attrs.get("ID", "")] = new
        f[8] = ";".join(parts)
        out.append("\t".join(f) + ("\n" if ln.endswith("\n") else ""))
    at = 1 if out and out[0].startswith("##gff-version") else 0
    out.insert(at, f"{MARKS_HEADER} {locus}\n")
    return out, marked


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--report", default="", help="synvoy_report.json (or the NO_REPORT sentinel)")
    ap.add_argument("--locus", default="",
                    help="the home locus of the GFFs (locus_1, or the home BED's name "
                         "synteny_block_locus_1): marks their GOI models with the report's outcome "
                         "(figures.<locus>) and leaves out the models that block marks dropped; "
                         "without it no model is marked, and the dropped records' spans decide "
                         "what is left out")
    ap.add_argument("--out_dir", required=True, help="folder the GFFs are written to")
    ap.add_argument("gffs", nargs="*", help="per-genome GFFs of one home locus")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    report = _load_report(args.report)
    spans = _dropped_spans_of(report)
    locus = locus_of(args.locus)
    marks = report_marks(report, locus)
    for path in args.gffs:
        out_path = os.path.join(args.out_dir, os.path.basename(path))
        genome = canonical_genome_id(path)
        genome_spans = spans.get(genome, [])
        if not genome_spans and marks is None:
            shutil.copyfile(path, out_path)
            continue
        with open(path) as fh:
            lines = fh.readlines()
        marked = {}
        if marks is not None:
            # The block names the models of the dropped records: those, not every model of the
            # dropped kind that overlaps a dropped record (a kept answer beside an answer dropped
            # at another home locus, say).
            genome_marks = marks.get(genome, {})
            lines, left_out = drop_models(lines, [mid for mid, m in genome_marks.items()
                                                  if m.get("ReportOutcome") == "dropped"])
            lines, marked = mark_lines(lines, genome_marks, locus)
        else:
            lines, left_out = _drop_rows(lines, genome_spans)
        with open(out_path, "w") as out:
            out.writelines(lines)
        if left_out:
            print(f"[report] {os.path.basename(path)}: left out {len(left_out)} model(s) the report "
                  f"dropped (gene-order-only beside a counted call, or a strong hit without synteny "
                  f"the report did not keep): {', '.join(left_out)}")
        outcomes = Counter(m["ReportOutcome"] for m in marked.values() if m.get("ReportOutcome"))
        also = sum(1 for m in marked.values() if m.get("AlsoFoundFrom"))
        if outcomes or also:
            what = [f"{n} {o}" for o, n in sorted(outcomes.items())]
            what += [f"{also} also found from another home locus"] if also else []
            print(f"[report] {os.path.basename(path)}: marked the report's outcome on "
                  f"{len(marked)} GOI model(s) of {locus}: {', '.join(what)}")


if __name__ == "__main__":
    main()
