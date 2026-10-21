"""scripts/compare_plot_selection.py (plot selection, steps 2-3, 2026-10-19/21): what the
figures show by the candidate regions (plot_synteny.select_target_genes) against the choice by
the label's scoring (the call's own block first, then the blocks by E_syn) that the figures make
with block tables (plot_synteny.choose_blocks), on a small run folder. Its default is the
figures' rule (Calls+significant); FIRST_RULE is the rule the step-2 comparison ran."""
import csv
import importlib.util
import os
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import iterative_search_runner as isr  # noqa: E402

FIRST_RULE = ("--call_labels", "HIGH,MEDIUM,AMBIGUOUS,LOW", "--context_max_esyn", "inf")

_spec = importlib.util.spec_from_file_location(
    "compare_plot_selection", os.path.join(ROOT, "scripts", "compare_plot_selection.py"))
cps = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cps)

# block: (span, neighbourhood, E_syn); block 4 has no model in the plot GFF
BLOCKS = {0: ((100_000, 200_000), (90_000, 210_000), 1e-6),
          1: ((500_000, 520_000), (480_000, 540_000), 0.5),
          2: ((800_000, 900_000), (780_000, 920_000), 1e-3),
          3: ((1_500_000, 1_600_000), (1_480_000, 1_620_000), 3.0),
          4: ((2_000_000, 2_100_000), (1_980_000, 2_120_000), 0.05)}


def _model(mid, start, end, role, conf, extra=""):
    return (f"c1\tsynvoy\tmRNA\t{start}\t{end}\t.\t+\t.\tID={mid};Name={mid.split('|')[0]};"
            f"SynVoy_Parent={mid.split('|')[0]};SynVoyRole={role};Confidence={conf};Identity=80{extra}")


MODELS = [
    _model("gene-F1|G_fna_b0_fl1_flank_ann", 101_001, 102_000, "flanking", "HIGH"),
    _model("GOI_q|G_fna_b0_l1_exon_ann", 120_001, 121_500, "goi", "HIGH",
           ";ModelStatus=complete;Bits=200;Eeff=1e-30;BlockEsyn=1e-06"),
    _model("gene-F2|G_fna_b0_fl1_flank_ann", 150_001, 151_000, "flanking", "HIGH"),
    _model("GOI_q|G_fna_b1_l1_exon_ann", 505_001, 506_000, "goi", "MEDIUM",
           ";ModelStatus=complete;Bits=60;Eeff=1e-5;BlockEsyn=0.5"),
    _model("gene-F3|G_fna_b2_fl1_flank_ann", 810_001, 811_000, "flanking", "HIGH"),
    _model("gene-F5|G_fna_rearranged", 815_001, 816_000, "flanking", "MEDIUM"),
    _model("gene-F4|G_fna_b3_fl1_flank_ann", 1_510_001, 1_511_000, "flanking", "HIGH"),
    _model("gene-F6|G_fna_rearranged", 3_000_001, 3_001_000, "flanking", "MEDIUM"),
]


def _run(tmp_path, with_table=True, extra=()):
    plot = tmp_path / "plot_inputs_synteny_block_locus_1"
    plot.mkdir(parents=True)
    (plot / "G.fna.gff").write_text("##gff-version 3\n" + "\n".join(MODELS + list(extra)) + "\n")
    (plot / "G.fna.regions.bed").write_text("c1\t90000\t210000\tG|Reg1_G2_CLOW_S0.20\t0.20\t+\n")
    (tmp_path / "regions" / "locus_1").mkdir(parents=True)
    (tmp_path / "regions" / "locus_1" / "G.fna.scores.tsv").write_text(
        "region_rank\tregion_name\tchrom\tstart\tend\tregion_class\tconfidence\n"
        "1\tG|Reg1_G2_CLOW_S0.20\tc1\t90000\t210000\tgoi_overlap\tLOW\n")
    if with_table:
        (tmp_path / "blocks" / "locus_1").mkdir(parents=True)
        with open(tmp_path / "blocks" / "locus_1" / "locus_1_G.fna.blocks.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=isr.BLOCK_TABLE_COLUMNS, delimiter="\t")
            w.writeheader()
            for b, ((s, e), (h0, h1), esyn) in BLOCKS.items():
                w.writerow({"genome": "G.fna", "block": b, "chrom": "c1", "start": s, "end": e,
                            "hood_start": h0, "hood_end": h1, "esyn": repr(esyn),
                            "significant": int(esyn <= 0.1)})
    return tmp_path


def _compare(tmp_path, monkeypatch, capsys, *extra, with_table=True, models=()):
    run = _run(tmp_path, with_table, models)
    monkeypatch.setattr(sys, "argv", ["compare_plot_selection.py", str(run), *extra])
    cps.main()
    lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.startswith("locus_1 ")]
    with open(run / "plot_selection_compare.tsv", newline="") as fh:
        rows = {r["model"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return lines, rows


def test_the_calls_block_first_then_the_blocks_by_esyn(tmp_path, monkeypatch, capsys):
    (line,), rows = _compare(tmp_path, monkeypatch, capsys, *FIRST_RULE)
    assert {m: (r["drawn_by"], r["why"]) for m, r in rows.items()} == {
        "GOI_q|G_fna_b1_l1_exon_ann": ("blocks", "call_block"),            # missed today
        "gene-F3|G_fna_b2_fl1_flank_ann": ("blocks", "predicted_no_call"),
        "gene-F5|G_fna_rearranged": ("blocks", "in_drawn_span"),           # inside block 2
        "gene-F4|G_fna_b3_fl1_flank_ann": ("blocks", "context")}
    assert rows["GOI_q|G_fna_b1_l1_exon_ann"]["block_esyn"] == "0.5"
    assert "both 3, today only 0, blocks only 4" in line
    assert "counted calls 2, missing today 1 (MEDIUM b1)" in line
    assert ("blocks drawn 4 (call_block 2, predicted_no_call 1, context 1) of 5 searched, "
            "3 significant, predicted places 2 (1 not drawn)") in line


def test_the_cap_counts_the_call_blocks(tmp_path, monkeypatch, capsys):
    (line,), rows = _compare(tmp_path, monkeypatch, capsys, *FIRST_RULE, "--max_regions", "3")
    assert "gene-F4|G_fna_b3_fl1_flank_ann" not in rows
    assert rows["gene-F3|G_fna_b2_fl1_flank_ann"]["why"] == "predicted_no_call"
    assert "blocks drawn 3" in line


LOW_CALL = _model("GOI_q|G_fna_b3_l1_fallback", 1_520_001, 1_520_600, "goi", "LOW",
                  ";ModelStatus=partial;Bits=20;Eeff=40;BlockEsyn=3")


def test_a_low_call_draws_its_block_unless_the_labels_say_otherwise(tmp_path, monkeypatch, capsys):
    (_line,), rows = _compare(tmp_path / "a", monkeypatch, capsys, *FIRST_RULE, models=[LOW_CALL])
    assert rows["gene-F4|G_fna_b3_fl1_flank_ann"]["why"] == "call_block"
    assert rows["GOI_q|G_fna_b3_l1_fallback"]["why"] == "call_block"
    # the default, the figures' rule: no LOW call block, no non-significant fill-up
    (line,), rows = _compare(tmp_path / "b", monkeypatch, capsys, models=[LOW_CALL])
    assert {m: r["why"] for m, r in rows.items()} == {
        "GOI_q|G_fna_b1_l1_exon_ann": "call_block",
        "gene-F3|G_fna_b2_fl1_flank_ann": "predicted_no_call",
        "gene-F5|G_fna_rearranged": "in_drawn_span"}
    assert "blocks drawn 3 (call_block 2, predicted_no_call 1, context 0)" in line


FRAGMENT = _model("GOI_q|G_fna_b3_r1_fallback", 1_520_001, 1_520_300, "goi", "MEDIUM",
                  ";ModelStatus=fragment;Bits=40;Eeff=1e-4;BlockEsyn=3")


def test_a_fragment_call_draws_its_block_and_stays_hidden(tmp_path, monkeypatch, capsys):
    (line,), rows = _compare(tmp_path / "a", monkeypatch, capsys, models=[FRAGMENT])
    assert rows["gene-F4|G_fna_b3_fl1_flank_ann"]["why"] == "call_block"   # E_syn 3, its block
    assert "GOI_q|G_fna_b3_r1_fallback" not in rows                       # hidden in both
    assert "blocks drawn 4 (call_block 3, predicted_no_call 1, context 0)" in line
    (line,), rows = _compare(tmp_path / "b", monkeypatch, capsys, "--hide_fragments", models=[FRAGMENT])
    assert "gene-F4|G_fna_b3_fl1_flank_ann" not in rows
    assert "blocks drawn 3 (call_block 2, predicted_no_call 1, context 0)" in line


def test_without_a_block_table_only_the_calls_blocks(tmp_path, monkeypatch, capsys):
    (line,), rows = _compare(tmp_path, monkeypatch, capsys, with_table=False)
    assert set(rows) == {"GOI_q|G_fna_b1_l1_exon_ann"}
    assert "no blocks.tsv: steps 1-2 only" in line


def test_today_is_the_figures_own_choice(tmp_path, monkeypatch):
    run = _run(tmp_path)
    cps.ps._SHOW_GOI_FRAGMENTS = False
    regions = cps.ps.parse_candidate_regions([str(run / "plot_inputs_synteny_block_locus_1" / "G.fna.regions.bed")])
    today, _notes = cps.today_choice(str(run / "plot_inputs_synteny_block_locus_1" / "G.fna.gff"), "G", regions)
    assert sorted(today) == ["GOI_q|G_fna_b0_l1_exon_ann", "gene-F1|G_fna_b0_fl1_flank_ann",
                             "gene-F2|G_fna_b0_fl1_flank_ann"]


def test_the_figures_call_the_same_selection():
    src = open(os.path.join(ROOT, "bin", "plot_synteny.py")).read()
    assert src.count("_cap_goi_genes(parse_target_gff(gff_file, cap=False)") == 1
    assert ("genes, candidate_regions = select_target_genes(\n"
            "                gff_file, genome_id, candidate_regions_by_genome, args.flank_fallback_bp)") in src
    assert "genes, candidate_regions = select_target_genes_by_blocks(gff_file, genome_id, block_choice)" in src
    cps_src = open(os.path.join(ROOT, "scripts", "compare_plot_selection.py")).read()
    assert "ps.choose_blocks(" in cps_src and "ps.models_in_choice(" in cps_src
