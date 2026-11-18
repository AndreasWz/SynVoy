"""Item 36 (bug list of 16 Nov 2026): the anchor grids drew what the report merged.

The figures deduped the GOI models with a rule of their own (reciprocal overlap >= 0.5), while
the report makes one record of two models that overlap by more than 80 %, that lie one at least
90 % inside the other, or that are the two halves of one reading frame, and keeps a compact model
over a wide one that adds barely any query coverage (generate_report.dedupe_goi_annotations). So
the grid drew a superseded wide model and its compact model as two copies, the two halves of one
reading frame as two copies, and two genes the report keeps apart (60 % reciprocal overlap) as one.

Now the report writes, per home locus, the outcome of each GOI model of the search (figures.<locus>:
merged, dropped, or why it does not count a call drawn, and the other home loci that found a record
drawn there); apply_report_drops.py --locus marks the plot-input GFFs with it (ReportOutcome,
ReportRecord, AlsoFoundFrom, header ##synvoy-report-marks); and the figures leave out the models the
report merged or dropped, keep the calls of two records apart, and name the outcome in the tooltips
and in the source data. A GFF without the marks is drawn as before.
"""
import contextlib
import csv
import html
import io
import json
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

import apply_report_drops as ard  # noqa: E402
import generate_report as gr  # noqa: E402
import plot_synteny as ps  # noqa: E402
import synvoy_grid as sg  # noqa: E402

CLASS = {"HIGH": "confident_goi", "MEDIUM": "probable_goi", "AMBIGUOUS": "syntenic_candidate_unconfirmed",
         "LOW": "ambiguous_goi_family_member"}


def _model(mid, s, e, conf, qcov, cds, status="complete", chrom="chr1", ident="80.0"):
    """A GOI gene model of the search: its mRNA row and CDS rows (phase 0, + strand)."""
    rows = [f"{chrom}\tminiprot\tmRNA\t{s}\t{e}\t{ident}\t+\t.\tID={mid};Name=GOI_X;SynVoy_Parent=GOI_X;"
            f"Identity={ident};Confidence={conf};GOIClass={CLASS[conf]};SynVoyRole=goi;"
            f"EvidenceType=exon_annotation;ModelStatus={status};QueryCoverage={qcov};Bits=300.0;"
            f"Eeff=1e-20;BlockEsyn=1e-5;Exons={len(cds)}"]
    rows += [f"{chrom}\tminiprot\tCDS\t{cs}\t{ce}\t.\t+\t0\tID={mid}_CDS{i};Parent={mid}"
             for i, (cs, ce) in enumerate(cds, 1)]
    return rows


def _flank(chrom, s, e, gid):
    return (f"{chrom}\tflanking_miniprot\tmRNA\t{s}\t{e}\t80.0\t+\t.\tID={gid}|G_fna_{s}_flank_ann;"
            f"Name={gid};SynVoy_Parent={gid};Identity=80.0;Confidence=HIGH;SynVoyRole=flanking;"
            f"EvidenceType=flanking_miniprot")


ID1, ID2 = "GOI_X|G_fna_b1_l1_exon_ann", "GOI_X|G_fna_b1_l1m2_exon_ann"
# (a) an 11.4-kb HIGH model swallowing the 2.6-kb MEDIUM model of the same gene: the report keeps
# the compact one (superseded_wide_model).
WIDE = _model(ID1, 1000, 12400, "HIGH", 0.988, [(1000, 1500), (5000, 5600), (12000, 12400)])
COMPACT = _model(ID2, 5000, 7600, "MEDIUM", 0.904, [(5000, 5600), (7000, 7600)])
# (b) two halves of one reading frame, meeting in 10 bp of one frame (same_orf_merged).
UP = _model(ID1, 1000, 1800, "HIGH", 0.55, [(1000, 1200), (1500, 1800)], status="partial")
DOWN = _model(ID2, 1791, 2600, "MEDIUM", 0.5, [(1791, 2100), (2400, 2600)], status="partial")
# (c) two models overlapping 60 % reciprocally, neither inside the other, in other frames: two genes.
P = _model(ID1, 1000, 2000, "HIGH", 0.9, [(1000, 2000)])
Q = _model(ID2, 1400, 2400, "MEDIUM", 0.9, [(1400, 2400)])


def _report(tmp_path, staged):
    """build_report on staged region GFFs ({staged name: rows}), written to report.json."""
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir(exist_ok=True)
    for name, rows in staged.items():
        (tmp_path / "regions" / name).write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    report = gr.build_report(str(tmp_path))
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    return report, path


def _search_gff(tmp_path, rows, name="G.fna.gff"):
    (tmp_path / "search").mkdir(exist_ok=True)
    path = tmp_path / "search" / name
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return path


def _apply(tmp_path, report_path, gffs, locus="synteny_block_locus_1", out="plot_inputs"):
    proc = subprocess.run([sys.executable, os.path.join(BIN, "apply_report_drops.py"), "--report",
                           str(report_path), "--locus", locus, "--out_dir", str(tmp_path / out)]
                          + [str(g) for g in gffs], capture_output=True, text=True, check=True)
    return proc.stdout


def _copies(gff):
    """The figure's GOI calls of one GFF as [(model id, is a copy)], and its console."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        genes = ps.parse_target_gff(str(gff), cap=False)
    goi = [g for g in genes if ps._is_goi_target_gene(g)]
    return [(g["gff_id"], c) for g, c in ps._goi_copies_and_weak(goi)], out.getvalue()


# ------------------------------------------------------------- the report and the marks ----

def test_a_superseded_wide_model_is_not_drawn_beside_its_compact_model(tmp_path):
    report, path = _report(tmp_path, {"locus_1__G.fna.gff": WIDE + COMPACT})
    (rec,) = report["goi_dedup"]["records"]
    assert rec["mrna_id"] == ID2 and rec["superseded_wide_model"]
    assert report["figures"] == {"locus_1": {
        "models": [{"genome": "G", "model": ID1, "outcome": "merged", "record": "chr1:5000-7600"}],
        "also_found_from": []}}
    assert "trees" not in report                      # written without a tree
    gff = _search_gff(tmp_path, WIDE + COMPACT)
    assert _copies(gff)[0] == [(ID1, True), (ID2, True)]      # before: two copies
    log = _apply(tmp_path, path, [gff])
    assert "marked the report's outcome on 1 GOI model(s) of locus_1: 1 merged" in log
    copies, console = _copies(tmp_path / "plot_inputs" / "G.fna.gff")
    assert copies == [(ID2, True)]
    assert "left out 1 GOI model(s) the report merged into the call drawn for their gene" in console


def test_the_two_halves_of_one_reading_frame_are_one_call(tmp_path):
    report, path = _report(tmp_path, {"locus_1__G.fna.gff": UP + DOWN})
    (rec,) = report["goi_dedup"]["records"]
    assert rec["mrna_id"] == ID1 and rec["same_orf_merged"]
    assert [(m["model"], m["outcome"]) for m in report["figures"]["locus_1"]["models"]] == [(ID2, "merged")]
    gff = _search_gff(tmp_path, UP + DOWN)
    assert _copies(gff)[0] == [(ID1, True), (ID2, True)]
    _apply(tmp_path, path, [gff])
    assert _copies(tmp_path / "plot_inputs" / "G.fna.gff")[0] == [(ID1, True)]


def test_two_genes_the_report_keeps_apart_are_both_drawn(tmp_path):
    report, path = _report(tmp_path, {"locus_1__G.fna.gff": P + Q})
    assert len(report["goi_dedup"]["records"]) == 2
    assert report["figures"]["locus_1"] == {"models": [], "also_found_from": []}
    gff = _search_gff(tmp_path, P + Q)
    assert _copies(gff)[0] == [(ID1, True)]          # the figures' own dedup dropped Q
    _apply(tmp_path, path, [gff])
    marked = (tmp_path / "plot_inputs" / "G.fna.gff").read_text().splitlines()
    assert marked[:2] == ["##gff-version 3", "##synvoy-report-marks locus_1"]
    assert marked[2:] == P + Q                       # no row marked, none removed
    assert _copies(tmp_path / "plot_inputs" / "G.fna.gff")[0] == [(ID1, True), (ID2, True)]


def test_a_gene_found_from_two_home_loci_names_the_other_locus(tmp_path):
    """One record (cross_locus_duplicate), drawn once in each locus's figures, as in its trees;
    each drawing names the other locus."""
    b2, b7 = "GOI_X|G_fna_b2_l1_exon_ann", "GOI_X|G_fna_b7_l1_exon_ann"
    cds = [(5000, 5600), (7000, 7600)]
    report, path = _report(tmp_path, {"locus_1__G.fna.gff": _model(b2, 5000, 7600, "HIGH", 0.95, cds),
                                      "locus_2__G.fna.gff": _model(b7, 5000, 7600, "HIGH", 0.95, cds)})
    (rec,) = report["goi_dedup"]["records"]
    assert rec["cross_locus_duplicate"]
    assert report["figures"] == {
        "locus_1": {"models": [], "also_found_from": [{"genome": "G", "model": b2, "loci": ["locus_2"]}]},
        "locus_2": {"models": [], "also_found_from": [{"genome": "G", "model": b7, "loci": ["locus_1"]}]}}
    gff = _search_gff(tmp_path, _model(b2, 5000, 7600, "HIGH", 0.95, cds))
    log = _apply(tmp_path, path, [gff])
    assert "1 also found from another home locus" in log
    assert f"ID={b2};" in (row := next(ln for ln in (tmp_path / "plot_inputs" / "G.fna.gff").read_text()
                                       .splitlines() if "\tmRNA\t" in ln))
    assert row.endswith(";AlsoFoundFrom=locus_2")
    with contextlib.redirect_stdout(io.StringIO()):
        (g,) = ps.parse_target_gff(str(tmp_path / "plot_inputs" / "G.fna.gff"), cap=False)
    assert g["also_found_from"] == ["locus_2"] and g["report_outcome"] == ""
    assert ps._report_note(g) == "also found from locus_2: one record of the report"


def _row(mid, start, end, conf="HIGH", genome="G", source="locus_1", rescue_pass="", status="complete"):
    return {"genome": genome, "source": source, "chrom": "c1", "start": start, "end": end, "mrna_id": mid,
            "confidence": conf, "goi_class": CLASS[conf], "rescue_pass": rescue_pass, "model_status": status}


def _rec(rep, start, end, conf="HIGH", cls=None, members=(), overruled=(), prov=("locus_1",), **kw):
    return dict({"genome": "G", "chrom": "c1", "start": start, "end": end, "mrna_id": rep[1],
                 "confidence": conf, "goi_class": cls or CLASS[conf], "provenance": list(prov),
                 "_rep": rep, "_members": list(members) or [rep + (start, end)],
                 "_overruled": list(overruled)}, **kw)


DEMOTED = {"paralog_not_goi", "identity_coverage_decoupled"}


def test_figure_marks_name_why_a_call_drawn_is_not_counted():
    rows = [_row("A", 100, 900), _row("A2", 120, 880, conf="MEDIUM"),          # paralog_not_goi, 2 models
            _row("B", 2000, 2900, conf="AMBIGUOUS"), _row("B2", 2050, 2850, conf="AMBIGUOUS"),
            _row("C", 5000, 5900, conf="MEDIUM"),                              # a kept answer
            _row("D", 8000, 8900),                                             # counted, its record a rescue's
            _row("D.mRNA", 8000, 8950, rescue_pass="hull_rescue"),
            _row("E", 9500, 9600, conf="MEDIUM"),                              # in no record
            _row("L", 300, 600, conf="LOW")]                                   # LOW: no figure mark
    records = [
        _rec(("locus_1", "A"), 100, 900, cls="paralog_not_goi",
             members=[("locus_1", "A", 100, 900), ("locus_1", "A2", 120, 880)]),
        _rec(("locus_1", "B"), 2000, 2900, conf="AMBIGUOUS",
             members=[("locus_1", "B", 2000, 2900), ("locus_1", "B2", 2050, 2850)]),
        _rec(("locus_1", "C"), 5000, 5900, conf="MEDIUM", cls="genome_wide_strong_hit", genome_wide_strong_hit=True),
        _rec(("locus_1", "D.mRNA"), 8000, 8950,
             members=[("locus_1", "D.mRNA", 8000, 8950), ("locus_1", "D", 8000, 8900)])]
    out = gr.figure_marks(rows, records, "locus_1", DEMOTED)
    assert [(m["model"], m["outcome"], m["record"]) for m in out["models"]] == [
        ("A", "paralog_not_goi", "c1:100-900"), ("A2", "merged", "c1:100-900"),
        ("B2", "merged", "c1:2000-2900"),
        ("C", "genome_wide_strong_hit", "c1:5000-5900"),
        ("E", "not_in_report", "")]
    assert out["also_found_from"] == []
    assert {m["outcome"] for m in out["models"]} <= set(ps.REPORT_OUTCOME_TEXT)   # the figures name each
    # another home locus sees none of them
    assert gr.figure_marks(rows, records, "locus_2", DEMOTED) == {"models": [], "also_found_from": []}


def test_figure_marks_of_an_overruled_piece_and_of_a_dropped_record():
    rows = [_row("W", 1000, 9000, conf="AMBIGUOUS"), _row("P", 2000, 2600, status="fragment"),
            _row("X", 20000, 20900, conf="AMBIGUOUS"), _row("X2", 20010, 20890, conf="AMBIGUOUS"),
            _row("Y", 30000, 31000, conf="MEDIUM")]
    whole = _rec(("locus_1", "W"), 1000, 9000, conf="AMBIGUOUS", basis="not_specific",
                 members=[("locus_1", "W", 1000, 9000), ("locus_1", "P", 2000, 2600)],
                 overruled=[("locus_1", "P")])
    # the whole gene found from another locus only: this locus draws the record by its piece
    elsewhere = _rec(("locus_2", "W9"), 1000, 9000, conf="AMBIGUOUS", prov=("locus_1", "locus_2"),
                     members=[("locus_2", "W9", 1000, 9000), ("locus_1", "P", 2000, 2600)],
                     overruled=[("locus_1", "P")])
    dropped = _rec(("locus_1", "X"), 20000, 20900, conf="AMBIGUOUS", basis="gene_order_only",
                   members=[("locus_1", "X", 20000, 20900), ("locus_1", "X2", 20010, 20890)])
    stand = _rec(("locus_1", "Y"), 30000, 31000, conf="MEDIUM")
    out = gr.figure_marks(rows[:2] + rows[2:], [whole, stand], "locus_1", DEMOTED,
                          dropped=[(dropped, "gene_order_only")])
    assert [(m["model"], m["outcome"]) for m in out["models"]] == [
        ("P", "merged"), ("X", "dropped"), ("X2", "dropped")]
    out = gr.figure_marks([rows[1]], [elsewhere], "locus_1", DEMOTED)
    assert [(m["model"], m["outcome"]) for m in out["models"]] == [("P", "overruled_piece")]
    assert {"merged", "dropped", "overruled_piece"} <= set(ps.REPORT_OUTCOME_TEXT)
    assert set(ps.REPORT_LEFT_OUT) == {"merged", "dropped"}
    assert out["also_found_from"] == [{"genome": "G", "model": "P", "loci": ["locus_2"]}]


def test_the_trees_keep_their_rule():
    """The figures choose a record's drawing by the trees' rule (_shown_member); the trees' list is
    unchanged."""
    rows = [_row("A", 1000, 3000), _row("A2", 1100, 2900, conf="MEDIUM")]
    rec = _rec(("locus_1", "A"), 1000, 3000, members=[("locus_1", "A", 1000, 3000), ("locus_1", "A2", 1100, 2900)])
    assert gr.tree_leaves_not_counted(rows, [rec], "locus_1", ["A", "A2"]) == [
        {"leaf": "A2", "reason": "merged", "record": "G c1:1000-3000"}]
    assert gr._shown_member(rec, "locus_1", {"A2"}) == "A2"
    assert gr._shown_member(rec, "locus_1", set()) is None


def test_rescue_gffs_are_no_figure_input(tmp_path):
    assert [gr.rescue_pass_of(n) for n in ("locus_1__G.hull_rescue.gff", "locus_2__G.rescue.gff",
                                           "locus_1__G.fna.gff", "G.fa.gff")] == ["hull_rescue", "rescue", "", ""]
    hull = ["chr1\tminiprot\tgene\t5000\t7600\t.\t+\t.\tID=hull1;SynVoyRole=goi;Confidence=HIGH;"
            "GOIClass=synteny_hull_rescue;EvidenceType=relaxed_miniprot;Identity=80.0;QueryCoverage=0.95",
            "chr1\tminiprot\tmRNA\t5000\t7600\t.\t+\t.\tID=hull1.mRNA;Parent=hull1;SynVoyRole=goi;"
            "Confidence=HIGH;GOIClass=synteny_hull_rescue;EvidenceType=relaxed_miniprot;Identity=80.0;"
            "QueryCoverage=0.95;Bits=900"]
    report, _ = _report(tmp_path, {"locus_1__G.fna.gff": COMPACT, "locus_1__G.hull_rescue.gff": hull})
    assert sorted(report["figures"]) == ["locus_1"]
    # the hull model represents the record; the search's model draws it, unmarked
    assert report["figures"]["locus_1"] == {"models": [], "also_found_from": []}


# ---------------------------------------------------------------- apply_report_drops.py ----

def test_marks_are_appended_and_rows_kept(tmp_path):
    report = {"figures": {"locus_1": {
        "models": [{"genome": "G", "model": ID1, "outcome": "merged", "record": "chr1:5000-7600"}],
        "also_found_from": [{"genome": "G", "model": ID2, "loci": ["locus_2", "locus_3"]}]}}}
    lines = ["##gff-version 3\n", "##synvoy-distances {}\n"] + [r + "\n" for r in WIDE + COMPACT]
    marks = ard.report_marks(report, ard.locus_of("synteny_block_locus_1"))
    out, marked = ard.mark_lines(lines, marks["G"], "locus_1")
    assert out[:3] == ["##gff-version 3\n", "##synvoy-report-marks locus_1\n", "##synvoy-distances {}\n"]
    assert len(out) == len(lines) + 1
    assert out[3] == WIDE[0] + ";ReportOutcome=merged;ReportRecord=chr1:5000-7600\n"
    assert out[7] == COMPACT[0] + ";AlsoFoundFrom=locus_2,locus_3\n"
    assert [ln for ln in out if "\tCDS\t" in ln] == [ln for ln in lines if "\tCDS\t" in ln]
    assert sorted(marked) == [ID1, ID2]
    again, _ = ard.mark_lines(out, marks["G"], "locus_1")      # marked twice reads as once
    assert again == out
    cleared, none = ard.mark_lines(out, {}, "locus_1")        # a mark the report no longer has goes
    assert cleared[3] == WIDE[0] + "\n" and none == {}


def test_locus_of_the_home_bed_name():
    assert ard.locus_of("synteny_block_locus_12") == "locus_12"
    assert ard.locus_of("synteny_block_locus_1.bed") == "locus_1"
    assert ard.locus_of("locus_2") == "locus_2"


@pytest.mark.parametrize("report,locus", [
    ({"goi_dedup": {}}, "locus_1"),                              # a report from before the block
    ({"figures": {"locus_2": {"models": [], "also_found_from": []}}}, "locus_1"),   # not this locus
    ({"figures": {"locus_1": {"models": [], "also_found_from": []}}}, ""),          # no --locus
])
def test_without_the_block_the_gffs_are_copied_byte_for_byte(tmp_path, report, locus):
    (tmp_path / "r.json").write_text(json.dumps(report))
    gff = _search_gff(tmp_path, WIDE + COMPACT)
    cmd = [sys.executable, os.path.join(BIN, "apply_report_drops.py"), "--report", str(tmp_path / "r.json"),
           "--out_dir", str(tmp_path / "o"), str(gff)] + (["--locus", locus] if locus else [])
    subprocess.run(cmd, check=True, capture_output=True)
    assert (tmp_path / "o" / "G.fna.gff").read_bytes() == gff.read_bytes()
    assert not ps._has_report_marks(str(tmp_path / "o" / "G.fna.gff"))


def test_the_module_passes_the_locus():
    nf = open(os.path.join(ROOT, "modules", "plot_synteny.nf")).read()
    assert "apply_report_drops.py --report $report --locus ${home_bed.baseName} --out_dir" in nf


# ------------------------------------------------------------------------- the figures ----

def test_without_marks_the_figures_dedup_as_before(tmp_path):
    for rows, expected in ((WIDE + COMPACT, [(ID1, True), (ID2, True)]), (P + Q, [(ID1, True)])):
        assert _copies(_search_gff(tmp_path, rows))[0] == expected


def test_a_low_call_under_a_record_stays_hidden(tmp_path):
    """The report merges no LOW model, so the figures' overlap rule still hides one under a call,
    also under a model the report merged."""
    low = _model("GOI_X|G_fna_b1_l1_fallback", 1100, 12300, "LOW", 0.3, [(1100, 12300)])
    gff = tmp_path / "G.fna.gff"
    gff.write_text("##gff-version 3\n##synvoy-report-marks locus_1\n" + "\n".join(
        [WIDE[0] + ";ReportOutcome=merged;ReportRecord=chr1:5000-7600"] + WIDE[1:] + COMPACT + low) + "\n")
    assert _copies(gff)[0] == [(ID2, True)]


def _tracks():
    names = [f"gene-F{k}" for k in range(4)]
    home = [{"chrom": "h", "start": 1000 + 2000 * k, "end": 2000 + 2000 * k, "start_plot": 1000 + 2000 * k,
             "end_plot": 2000 + 2000 * k, "name": n, "strand": "+", "exon_coords": []} for k, n in enumerate(names)]
    home.insert(2, {"chrom": "h", "start": 4600, "end": 4900, "start_plot": 4600, "end_plot": 4900,
                    "name": "GOI_T", "strand": "+", "exon_coords": []})

    def goi(start, conf="HIGH", **kw):
        g = {"chrom": "c1", "start": start, "end": start + 300, "start_plot": start, "end_plot": start + 300,
             "name": f"GOI_T|m{start}", "gff_id": f"GOI_T|m{start}", "strand": "+", "role": "goi",
             "home_gene_id": "GOI_T", "confidence": conf, "identity": 70.0, "query_coverage": 0.95,
             "evidence_type": "exon_annotation", "exon_coords": [(start, start + 100), (start + 200, start + 300)],
             "report_outcome": "", "report_record": "", "also_found_from": []}
        g.update(kw)
        return g
    genes = [{"chrom": "c1", "start": 1000 + 2000 * k, "end": 2000 + 2000 * k, "start_plot": 1000 + 2000 * k,
              "end_plot": 2000 + 2000 * k, "name": f"t_{k}", "gff_id": f"t_{k}", "strand": "+",
              "exon_coords": [], "home_gene_id": n, "identity": 80.0, "confidence": "HIGH",
              "query_coverage": 0.95} for k, n in enumerate(names)]
    genes += [goi(4700, also_found_from=["locus_2"]),
              goi(5200, "MEDIUM", report_outcome="paralog_not_goi", report_record="c1:5200-5500")]
    left_out = [goi(4650, "HIGH", report_outcome="merged", report_record="c1:4700-5000")]
    target = {"label": "T0", "species": "Bombus species0", "is_home": False, "genome_id": "g0", "offset": 0,
              "breaks": [], "genes": genes, "goi_outside_view": [], "goi_report_left_out": left_out}
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home",
               "offset": 0, "breaks": [], "genes": home}, target]
    return tracks, {n: ps.GENE_PALETTE[k] for k, n in enumerate(names)}


def test_the_grid_lists_what_the_report_did_and_draws_each_record_once():
    tracks, colours = _tracks()
    _svg, legend, source, info = ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace())
    rows = [r for r in source if r["row"] == 2 and r["column"] == "T"]
    assert [(r["model_id"], r["drawn_as"], r["report"], r["report_record"], r["also_found_from"])
            for r in rows] == [
        ("GOI_T|m4700", "model", "", "", "locus_2"),
        ("GOI_T|m5200", "model", "paralog_not_goi", "c1:5200-5500", ""),
        ("GOI_T|m4650", "not drawn", "merged", "c1:4700-5000", "")]
    assert rows[-1]["state"] == "placed" and rows[-1]["strand_drawn"] == "" and rows[-1]["position"] == ""
    assert set(source[0]) == set(sg.SOURCE_COLUMNS)
    assert info["goi_copies"] == 2 and info["goi_models_left_out_by_report"] == 1
    assert info["goi_calls_not_counted_by_report"] == {"paralog_not_goi": 1}
    assert info["goi_calls_also_found_from"] == 1
    assert ("1 further T model is not drawn: the report merged it into the call drawn for its gene"
            in legend)
    assert ("1 T call in the T column is not counted by the report: another home gene's call "
            "(paralog_not_goi, 1) (source data, column report).") in legend
    assert ("1 T call was also found from another home locus; the report holds it as one record "
            "(source data, column also_found_from).") in legend
    full = html.unescape(ps.render_anchor_grid(tracks, colours, {}, {}, SimpleNamespace()))
    assert "not counted by the report: another home gene's call (paralog_not_goi)" in full
    assert "also found from locus_2: one record of the report" in full
    positions = html.unescape(ps.render_anchor_grid_positional(tracks, colours, {}, {}, SimpleNamespace()))
    assert "2 copies; also found from locus_2: one record of the report" in positions


def test_a_grid_without_report_outcomes_says_nothing_of_them():
    tracks, colours = _tracks()
    for g in tracks[1]["genes"]:
        g.update(report_outcome="", report_record="", also_found_from=[])
    tracks[1]["goi_report_left_out"] = []
    _svg, legend, source, info = ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace())
    assert info["goi_models_left_out_by_report"] == 0 and info["goi_calls_not_counted_by_report"] == {}
    assert "report" not in legend.replace("reported", "")
    assert all(r["report"] == "" for r in source)


def test_the_ribbon_tooltip_names_the_outcome():
    tracks, _ = _tracks()
    paralog, also = tracks[1]["genes"][-1], tracks[1]["genes"][-2]
    data = json.loads(html.unescape(ps._build_tooltip_json(paralog, tracks[1], {})))
    assert data["report"] == "not counted: another home gene's call (paralog_not_goi)"
    data = json.loads(html.unescape(ps._build_tooltip_json(also, tracks[1], {})))
    assert data["alsoFound"] == "locus_2 (one record of the report)" and "report" not in data


# ------------------------------------------------------------------------- end to end ----

HOME = [("FLA", 1000, 2000), ("FLB", 3000, 4000), ("QRY1", 5000, 7600), ("FLC", 13000, 14000),
        ("FLD", 15000, 16000)]


def test_report_marks_reach_the_figures(tmp_path):
    """build_report, apply_report_drops.py and plot_synteny.py on one synthetic target: the wide
    model the report merged into the compact one is not drawn, and has its row of the source data."""
    target = [_flank("tchr5", s + 500_000, e + 500_000, f"gene-{g}") for g, s, e in HOME if g != "QRY1"]
    shift = 500_000

    def moved(rows):
        out = []
        for r in rows:
            f = r.split("\t")
            f[0], f[3], f[4] = "tchr5", str(int(f[3]) + shift), str(int(f[4]) + shift)
            out.append("\t".join(f).replace("GOI_X", "GOI_QRY1").replace("|G_fna", "|GCA_000000009_1_fna"))
        return out
    goi_rows = moved(WIDE) + moved(COMPACT)
    _report_dir = tmp_path / "rep"
    _report_dir.mkdir()
    report, report_path = _report(_report_dir, {"locus_1__GCA_000000009.1.fna.gff": target + goi_rows})
    wide_id = ID1.replace("GOI_X", "GOI_QRY1").replace("|G_fna", "|GCA_000000009_1_fna")
    assert [(m["model"], m["outcome"]) for m in report["figures"]["locus_1"]["models"]] == [(wide_id, "merged")]
    (tmp_path / "synteny_block_locus_1.bed").write_text("".join(
        f"hchr1\t{s - 1}\t{e}\tgene-{g}\t.\t+\t{g}\n" for g, s, e in HOME if g != "QRY1"))
    (tmp_path / "locus_1.bed").write_text("hchr1\t5100\t7500\tgene_loc\t1e-40\t+\n")
    home = ["##gff-version 3"]
    for g, s, e in HOME:
        home += [f"hchr1\tRefSeq\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{g};Name={g}",
                 f"hchr1\tRefSeq\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{g};Parent=gene-{g}",
                 f"hchr1\tRefSeq\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{g};Parent=rna-{g}"]
    (tmp_path / "home.gff").write_text("\n".join(home) + "\n")
    gff = _search_gff(tmp_path, target + goi_rows, "GCA_000000009.1.fna.gff")
    _apply(tmp_path, report_path, [gff], out="plot_inputs_synteny_block_locus_1")
    marked = tmp_path / "plot_inputs_synteny_block_locus_1" / "GCA_000000009.1.fna.gff"
    out = tmp_path / "synteny_block_locus_1_synteny_plot.html"
    proc = subprocess.run(
        [sys.executable, os.path.join(BIN, "plot_synteny.py"),
         "--home_bed", str(tmp_path / "synteny_block_locus_1.bed"),
         "--query_bed", str(tmp_path / "locus_1.bed"), "--home_gff", str(tmp_path / "home.gff"),
         "--target_gffs", str(marked), "--report", str(report_path),
         "--output", str(out), "--no_network", "--common_names", "off", "--no_fragment_variant"],
        capture_output=True, text=True, cwd=tmp_path,
        env=dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1"))
    assert proc.returncode == 0, proc.stderr
    assert ("GCA_000000009.1.fna.gff: left out 1 GOI model(s) the report merged into the call drawn "
            "for their gene, or dropped") in proc.stdout
    assert "carry no report marks" not in proc.stdout
    with open(tmp_path / "synteny_block_locus_1_anchor_grid_print.source.tsv", newline="") as fh:
        goi = [r for r in csv.DictReader(fh, delimiter="\t")
               if r["genome_id"] == "GCA_000000009.1" and r["column"] == "QRY1"]
    assert [(r["confidence"], r["drawn_as"], r["report"]) for r in goi] == [
        ("MEDIUM", "model", ""), ("HIGH", "not drawn", "merged")]
    assert goi[1]["model_id"] == wide_id and goi[1]["report_record"] == "tchr5:505000-507600"
    assert "1 further QRY1 model is not drawn" in (
        tmp_path / "synteny_block_locus_1_anchor_grid_print.legend.txt").read_text()
    # the search's own GFF, without the marks: both models are drawn, and the log says why
    proc = subprocess.run(
        [sys.executable, os.path.join(BIN, "plot_synteny.py"),
         "--home_bed", str(tmp_path / "synteny_block_locus_1.bed"),
         "--query_bed", str(tmp_path / "locus_1.bed"), "--home_gff", str(tmp_path / "home.gff"),
         "--target_gffs", str(gff), "--output", str(tmp_path / "plain_synteny_plot.html"),
         "--no_network", "--common_names", "off", "--no_fragment_variant"],
        capture_output=True, text=True, cwd=tmp_path,
        env=dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1"))
    assert proc.returncode == 0, proc.stderr
    assert "1 of 1 target GFF(s) carry no report marks (##synvoy-report-marks)" in proc.stdout
    with open(tmp_path / "plain_anchor_grid_print.source.tsv", newline="") as fh:
        goi = [r for r in csv.DictReader(fh, delimiter="\t")
               if r["genome_id"] == "GCA_000000009.1" and r["column"] == "QRY1"]
    assert sorted(r["drawn_as"] for r in goi) == ["model", "model"]
