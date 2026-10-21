#!/usr/bin/env python3
"""Compare what the figures show by the candidate regions with the choice by the label's own
scoring that the figures make with block tables (plot selection, steps 2-3, 2026-10-19/21).
Read-only: it changes nothing in the pipeline and writes one table.

    scripts/compare_plot_selection.py <run dir> [--max_regions 6] [--out <tsv>]
        [--call_labels HIGH,MEDIUM,AMBIGUOUS] [--context_max_esyn E] [--hide_fragments]

Per home locus and target genome it reads what the figures read:
plot_inputs_synteny_block_locus_<N>/<genome>.gff (the models, after the report's drops) and the
candidate regions of CLUSTER_REGIONS (plot_inputs_*/*.regions.bed; their scores in
regions/locus_<N>/<genome>.scores.tsv), and the search's block table
blocks/locus_<N>/locus_<N>_<genome>.blocks.tsv.

today   plot_synteny.select_target_genes, the code the figures run without block tables,
        with the pipeline's settings (fragment GOI models hidden, at most 10 GOI models steer
        the choice, --flank_fallback_bp 1000000).
blocks  plot_synteny.choose_blocks + models_in_choice, the code the figures run with block
        tables (rule "Calls+significant", the defaults):
        1. the GOI models of the row (plot_synteny.choice_models: the models with fragments
           hidden plus every fragment GOI row), in the judge's order
           (plot_synteny._goi_priority_key);
        2. the search block (the N of _b<N>_ in its id) of each HIGH, MEDIUM and AMBIGUOUS one
           is drawn, whatever its E_syn;
        3. a significant block (E_syn <= 0.1) without a HIGH/MEDIUM model is a predicted
           place without a counted call;
        4. the other significant blocks with a model, by E_syn ascending, until --max_regions
           blocks are drawn (6, call blocks included);
        5. a model without a block number (cross-chromosome flanking recovery) is drawn when it
           lies inside the neighbourhood (hood_start-hood_end) of a drawn block.
        Every model of a drawn block is drawn, fragments hidden as in the default figures.
        Options try other rules: --call_labels, the labels whose GOI models draw their block
        in step 2 (a weaker model is drawn only when its block is drawn otherwise);
        --context_max_esyn, the largest E_syn a block may have to be drawn in step 4 instead
        of the search's `significant`; --hide_fragments, fragments left out of the choice.
        The rule the step-2 comparison ran: --call_labels HIGH,MEDIUM,AMBIGUOUS,LOW
        --context_max_esyn inf --hide_fragments.

Output (default <run>/plot_selection_compare.tsv): one row per (locus, genome, model) drawn by
exactly one of the two choices, with
  drawn_by  today or blocks
  why       call_block         step 2: the block holds a GOI model of the row (a HIGH/MEDIUM
                               one when the block is significant)
            predicted_no_call  E_syn <= 0.1 and no HIGH/MEDIUM model in the block (drawn by
                               step 2 for a weaker call, or by step 4)
            context            step 4, E_syn > 0.1 (only with --context_max_esyn above 0.1)
            in_drawn_span      step 5
            region_only        drawn today only; `region` names today's region that holds it
and one summary line per genome on stdout. A run from before the block table has no E_syn for
blocks without a GOI model: only steps 1-2 run there, a block's span is the span of its own
models, and the summary says so.
"""
import argparse
import contextlib
import csv
import glob
import io
import os
import re
import sys
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import plot_synteny as ps  # noqa: E402
from sequence_utils import block_index_of  # noqa: E402
from synvoy_evalue import SIGNIFICANT_BLOCK_ESYN  # noqa: E402

TARGET_GFF = re.compile(r"\.(fna|fa|fasta)(\.gz)?\.gff$")    # a target's models, as replot.sh
COUNTED = ps.COUNTED_LABELS
CALL_LABELS = ps.BLOCK_RULE_CALL_LABELS
FLANK_FALLBACK_BP = 1_000_000                                 # params.flank_fallback_bp
COLUMNS = ["locus", "genome", "model", "role", "label", "goi_class", "block", "block_esyn",
           "drawn_by", "why", "region"]


read_block_table = ps.read_block_table


def block_table_path(run, locus, genome):
    """The published table, else the copy a local gate run archived from the work folder."""
    for p in (os.path.join(run, "blocks", locus, f"{locus}_{genome}.blocks.tsv"),
              os.path.join(run, "archive", locus, "hits", f"{locus}_{genome}.blocks.tsv")):
        if os.path.exists(p):
            return p
    return None


def read_regions(run, locus, genome):
    """Today's regions of this genome: [(chrom, start0, end, name, class, confidence)]."""
    path = os.path.join(run, "regions", locus, f"{genome}.scores.tsv")
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return [(r["chrom"], int(r["start"]), int(r["end"]), r["region_name"].split("|")[-1],
                 r.get("region_class", ""), r.get("confidence", ""))
                for r in csv.DictReader(fh, delimiter="\t")]


def today_choice(gff, genome_id, regions_by_genome):
    """{model id: gene} the figure draws today, and the plot's notes on how it chose."""
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        genes, _regions = ps.select_target_genes(gff, genome_id, regions_by_genome, FLANK_FALLBACK_BP)
    notes = [line.split(": ", 1)[1].split(" (")[0].rstrip(".")
             for line in log.getvalue().splitlines() if line.startswith(f"[plot] {genome_id}: ")]
    return {g["gff_id"]: g for g in genes}, notes


def blocks_choice(genes_all, tag, table, max_regions, call_labels=CALL_LABELS,
                  context_max_esyn=None, choice_genes=None):
    """The block-based choice (module docstring), made by the figures' own code: {model id:
    (why, block)} of the models drawn (from genes_all, the default figures' models), {block:
    why}, the predicted places without a counted call, E_syn per block and {model id: block}.
    choice_genes: the models the choice reads (default genes_all)."""
    choice = ps.choose_blocks(genes_all if choice_genes is None else choice_genes, tag, table,
                              max_regions, call_labels, context_max_esyn)
    chosen = ps.models_in_choice(genes_all, tag, choice, call_labels)
    block_of = {g["gff_id"]: block_index_of(g["gff_id"], tag) for g in genes_all}
    return chosen, choice.drawn, choice.predicted, choice.esyn.get, block_of


def compare_genome(run, locus, gff, regions_by_genome, max_regions, call_labels=CALL_LABELS,
                   context_max_esyn=None, hide_fragments=False):
    """(rows, summary line) of one target genome."""
    genome = os.path.basename(gff)[:-len(".gff")]
    genome_id = ps.clean_genome_name(genome)
    tag = genome.replace(".", "_").replace("-", "_").replace(" ", "_")
    table = read_block_table(block_table_path(run, locus, genome))
    with contextlib.redirect_stdout(io.StringIO()):
        genes_all = ps.parse_target_gff(gff, cap=False)
    today, notes = today_choice(gff, genome_id, regions_by_genome)
    choice_genes = None if hide_fragments else ps.choice_models(gff)
    chosen, drawn, predicted, esyn, block_of = blocks_choice(genes_all, tag, table, max_regions,
                                                             call_labels, context_max_esyn,
                                                             choice_genes)
    regions = read_regions(run, locus, genome)
    genes = {g["gff_id"]: g for g in genes_all}
    genes.update({k: v for k, v in today.items() if k not in genes})

    def region_of(g):
        names = [f"{n}:{cls}" for c, s, e, n, cls, _conf in regions
                 if c == g["chrom"] and g["start"] - 1 < e and g["end"] > s]
        return ",".join(names) or ("goi_centred_fallback" if any("fallback" in n for n in notes) else "")

    rows = []
    for mid in sorted(set(today) ^ set(chosen), key=lambda m: (genes[m]["chrom"], genes[m]["start"], m)):
        g = genes[mid]
        b = block_of.get(mid, block_index_of(mid, tag))
        e = esyn(b) if b is not None else None
        rows.append({"locus": locus, "genome": genome, "model": mid, "role": g["role"],
                     "label": g["confidence"], "goi_class": g.get("goi_class", ""),
                     "block": "" if b is None else b, "block_esyn": "" if e is None else f"{e:.3g}",
                     "drawn_by": "today" if mid in today else "blocks",
                     "why": "region_only" if mid in today else chosen[mid][0],
                     "region": region_of(g) if mid in today else ""})
    calls = [g for g in genes_all if ps._is_goi_target_gene(g) and g["confidence"] in COUNTED]
    missing = [g for g in calls if g["gff_id"] not in today]
    whys = defaultdict(int)
    for w in drawn.values():
        whys[w] += 1
    line = (f"{locus} {genome}: both {len(set(today) & set(chosen))}, today only "
            f"{len(set(today) - set(chosen))}, blocks only {len(set(chosen) - set(today))}; "
            f"counted calls {len(calls)}, missing today {len(missing)}")
    if missing:
        line += " (" + ", ".join(f"{g['confidence']} b{block_of[g['gff_id']]}" for g in missing) + ")"
    missing_blocks = [g for g in calls if g["gff_id"] not in chosen]
    if missing_blocks:
        line += f", missing blocks {len(missing_blocks)}"
    line += (f"; blocks drawn {len(drawn)} (call_block {whys['call_block']}, predicted_no_call "
             f"{whys['predicted_no_call']}, context {whys['context']})")
    if table is not None:
        line += (f" of {len(table)} searched, {sum(1 for r in table.values() if r['esyn'] <= SIGNIFICANT_BLOCK_ESYN)}"
                 f" significant, predicted places {len(predicted)}")
        not_drawn = sorted(b for b in predicted if b not in drawn)
        if not_drawn:
            line += f" ({len(not_drawn)} not drawn)"
    else:
        line += "; no blocks.tsv: steps 1-2 only, block spans from the models"
    line += f"; today regions {len(regions)}"
    if notes:
        line += " [" + "; ".join(notes) + "]"
    return rows, line


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", help="a finished run's --outdir")
    ap.add_argument("--max_regions", type=int, default=ps.BLOCK_RULE_MAX_BLOCKS,
                    help="blocks drawn per genome at most, call blocks included (default 6, "
                         "the figures' cap)")
    ap.add_argument("--out", help="output TSV (default <run>/plot_selection_compare.tsv)")
    ap.add_argument("--call_labels", default=",".join(CALL_LABELS),
                    help="labels whose GOI models draw their block in step 2 (default HIGH, "
                         "MEDIUM, AMBIGUOUS)")
    ap.add_argument("--context_max_esyn", type=float, default=None,
                    help="largest E_syn of a block drawn in step 4 (default: the blocks the "
                         "search calls significant, E_syn <= 0.1; inf: any)")
    ap.add_argument("--hide_fragments", action="store_true",
                    help="leave fragment GOI models out of the choice, as the step-2 "
                         "comparison did (default: they take part, hidden from view)")
    args = ap.parse_args()
    run = os.path.abspath(args.run)
    ps._SHOW_GOI_FRAGMENTS = False
    ps._MAX_GOI_PER_GENOME = 10
    plot_dirs = sorted(glob.glob(os.path.join(run, "plot_inputs_synteny_block_locus_*")),
                       key=lambda p: int(p.rsplit("_", 1)[1]))
    if not plot_dirs:
        sys.exit(f"no plot_inputs_synteny_block_locus_* folder in {run}")
    out_rows = []
    for plot_dir in plot_dirs:
        locus = "locus_" + plot_dir.rsplit("_", 1)[1]
        with contextlib.redirect_stdout(io.StringIO()):
            regions_by_genome = ps.parse_candidate_regions(sorted(glob.glob(os.path.join(plot_dir, "*.regions.bed"))))
        for gff in sorted(p for p in glob.glob(os.path.join(plot_dir, "*.gff")) if TARGET_GFF.search(p)):
            rows, line = compare_genome(run, locus, gff, regions_by_genome, args.max_regions,
                                        tuple(args.call_labels.split(",")), args.context_max_esyn,
                                        args.hide_fragments)
            out_rows.extend(rows)
            print(line)
    out = args.out or os.path.join(run, "plot_selection_compare.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"{len(out_rows)} rows -> {out}")


if __name__ == "__main__":
    main()
