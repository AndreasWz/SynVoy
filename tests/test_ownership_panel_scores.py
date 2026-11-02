"""Fix-phase item 10, the ownership half: the check writes every panel score, the report reads the GOI's best entry."""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import generate_report as gr  # noqa: E402

pytest.importorskip("parasail")


def test_the_table_holds_every_entry(tmp_path):
    q = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTL"
    panel = tmp_path / "panel.faa"
    panel.write_text(f">locus_1|gene-Q\n{q}\n>home_gene-P|gene-P\n{q[:40]}\n>locus_1|gene-N\nMLLLLLLLLLLLLLLLLLLLLLLL\n")
    faa = tmp_path / "t.faa"
    faa.write_text(f">GOI_x|m1\n{q}\n")
    gff = tmp_path / "t.gff"
    gff.write_text("c1\texon_annotation\tmRNA\t100\t400\t90\t+\t.\tID=GOI_x|m1;SynVoyRole=goi;Confidence=HIGH\n")
    out = tmp_path / "out.tsv"
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "reciprocal_best_paralog_check.py"), "--home_query", str(panel),
                    "--target_faa", str(faa), "--target_gff", str(gff), "--locus_id", "locus_1", "--genome_name", "G",
                    "--output", str(out)], check=True)
    head, row = [ln.split("\t") for ln in out.read_text().splitlines()]
    rec = dict(zip(head, row))
    scores = gr._panel_scores(rec["panel_scores"])
    assert set(scores) == {"locus_1|gene-Q", "home_gene-P|gene-P", "locus_1|gene-N"}
    assert scores["locus_1|gene-Q"] == float(rec["best_bit"]) and rec["best_paralog"] == "locus_1|gene-Q"


@pytest.mark.parametrize("field,want", [
    ("", {}), (None, {}),
    ("locus_1|gene-Q=430.0;home_gene-P|gene-P=212.5", {"locus_1|gene-Q": 430.0, "home_gene-P|gene-P": 212.5}),
    ("locus_1|a=b=7.0", {"locus_1|a=b": 7.0}),
])
def test_panel_scores_read_back(field, want):
    assert gr._panel_scores(field) == want


def test_a_table_without_the_column_reads_its_two_best():
    row = {"best_paralog": "locus_4|gene-P", "best_bit_f": 205.0, "second_paralog": "locus_1|gene-Q", "second_bit_f": 200.0}
    assert gr._ownership_scores(row) == {"locus_4|gene-P": 205.0, "locus_1|gene-Q": 200.0}
    assert gr._ownership_scores(dict(row, second_paralog="")) == {"locus_4|gene-P": 205.0}
    assert gr._ownership_scores(dict(row, panel_scores_d={"a": 1.0})) == {"a": 1.0}


META = {"locus_1|gene-Q": {"is_goi": True, "goi_found_by": "goi_info"},
        "home_gene-C|gene-C": {"is_goi": True, "goi_found_by": "own_copy"},
        "locus_4|gene-P": {"is_goi": False, "goi_found_by": ""}}


@pytest.mark.parametrize("copy_score,goi", [
    (500.0, True),      # a tie
    (525.0, True),      # 25 raw = 9.6 bits over the GOI's own entry: not told apart
    (526.0, False),     # 26 raw = 10.0 bits
    (700.0, False),
    (300.0, True),
])
def test_an_own_copy_counts_as_the_goi_unless_the_call_favours_it(copy_score, goi):
    scores = {"locus_1|gene-Q": 500.0, "home_gene-C|gene-C": copy_score, "locus_4|gene-P": 100.0}
    assert gr._is_goi_for_call("home_gene-C|gene-C", scores, META) is goi
    assert gr._is_goi_for_call("locus_1|gene-Q", scores, META) is True
    assert gr._is_goi_for_call("locus_4|gene-P", scores, META) is False


def test_an_own_copy_without_the_goi_entrys_score_counts_as_the_goi():
    assert gr._is_goi_for_call("home_gene-C|gene-C", {"home_gene-C|gene-C": 700.0}, META) is True


def test_a_name_with_separators_does_not_lose_the_row():
    assert gr._panel_scores("locus_1|GOI_Q;1=430.0;x=2.5") == {"locus_1|GOI_Q;1": 430.0, "x": 2.5}
    assert gr._panel_scores("locus_1|GOI_Q=a;b=7.0") == {"locus_1|GOI_Q=a;b": 7.0}
    assert gr._panel_scores("garbage") == {}


def test_panel_names_hold_no_separator():
    import build_home_paralog_panel as bhp
    assert bhp.panel_name("Q;1=x y|z") == "Q_1_x_y_z"


def test_the_goi_entry_summary_names_the_gene_not_a_copy():
    meta = {"home_gene-A|gene-A": {"gene": "gene-A", "is_goi": True, "goi_found_by": "own_copy"},
            "locus_1|gene-Q": {"gene": "gene-Q", "is_goi": True, "goi_found_by": "goi_info"},
            "locus_1|gene-N": {"gene": "gene-N", "is_goi": False, "goi_found_by": ""}}
    assert gr.goi_panel_entry(meta) == {"gene": "gene-Q", "panel_ids": ["home_gene-A|gene-A", "locus_1|gene-Q"], "found_by": "goi_info"}
