"""Fix-phase item 3: "N need curation" leaves out another home gene's calls.

needs_curation_goi and ambiguous_goi counted the AMBIGUOUS records that ownership relabelled paralog_not_goi. Those
are another home gene's calls, so the report names them apart (ambiguous_other_home_gene_goi).
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import generate_report as gr  # noqa: E402

PANEL = ("panel_id\tlocus_id\tgene\tchrom\tis_goi_gene\tgoi_found_by\n"
         "locus_1|gene-G\tlocus_1\tgene-G\tc1\t1\tgoi_info\n"
         "home_gene-P|gene-P\thome_gene-P\tgene-P\tc9\t0\t\n")
OWN_HEAD = "genome\tlocus\tchrom\tstart\tend\tmrna_id\tbest_paralog\tbest_bit\tsecond_paralog\tsecond_bit\tbitscore_gap\n"


def _gff_row(start, end, mid, conf, basis, cls="syntenic_candidate_unconfirmed"):
    attrs = (f"ID={mid};Name=GOI_G;SynVoy_Parent=GOI_G;SynVoyRole=goi;Confidence={conf};GOIClass={cls};"
             f"EvidenceType=exon_annotation;Identity=60.0;QueryCoverage=0.900;ModelStatus=complete;Basis={basis}")
    return "\t".join(["c1", "exon_annotation", "mRNA", str(start), str(end), "60", "+", ".", attrs]) + "\n"


def _report(tmp_path):
    """Genome Gx: a HIGH call, an AMBIGUOUS call of the GOI and one that ownership gives to gene-P. Genome Gy, with
    no counted call: an AMBIGUOUS gene-order call, which stays."""
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "locus_ownership").mkdir()
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text(
        "##gff-version 3\n"
        + _gff_row(1000, 2000, "GOI_G|m1", "HIGH", "significant_specific", "confident_goi")
        + _gff_row(5000, 6000, "GOI_G|m2", "AMBIGUOUS", "not_specific")
        + _gff_row(20000, 21000, "GOI_G|m4", "AMBIGUOUS", "not_specific"))
    (stage / "regions" / "locus_1__Gy.fna.gff").write_text(
        "##gff-version 3\n" + _gff_row(9000, 10000, "GOI_G|m3", "AMBIGUOUS", "gene_order_only"))
    (stage / "locus_ownership" / "home_paralog_panel.faa.meta.tsv").write_text(PANEL)
    (stage / "locus_ownership" / "locus_1__Gx.locus_ownership.tsv").write_text(
        OWN_HEAD
        + "Gx\tlocus_1\tc1\t1000\t2000\tGOI_G|m1\tlocus_1|gene-G\t500\thome_gene-P|gene-P\t100\t400\n"
        + "Gx\tlocus_1\tc1\t5000\t6000\tGOI_G|m2\tlocus_1|gene-G\t300\thome_gene-P|gene-P\t250\t50\n"
        + "Gx\tlocus_1\tc1\t20000\t21000\tGOI_G|m4\thome_gene-P|gene-P\t400\tlocus_1|gene-G\t100\t300\n")
    out = tmp_path / "report.json"
    proc = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "generate_report.py"), "--results_dir", str(stage),
                           "--output", str(out), "--allow-empty"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return json.loads(out.read_text())


def test_another_home_genes_call_is_named_apart(tmp_path):
    summary = _report(tmp_path)["summary"]
    assert summary["needs_curation_goi"] == 2
    assert summary["ambiguous_goi"] == 2                         # the old name of the same number
    assert summary["needs_curation_by_basis"] == {"not_specific": 1, "gene_order_only": 1}
    assert summary["ambiguous_other_home_gene_goi"] == 1
    assert summary["high_confidence_goi"] == 1
    assert "2 need curation" in summary["headline"] and "1 of another home gene" in summary["headline"]


def test_the_record_stays_in_the_report(tmp_path):
    """Only the count moves: the record is there, relabelled, with its owner."""
    rec = [r for r in _report(tmp_path)["goi_dedup"]["records"] if r["start"] == 20000][0]
    assert (rec["confidence"], rec["goi_class"], rec["inferred_paralog"]) == ("AMBIGUOUS", "paralog_not_goi", "gene-P")


def test_a_run_without_ownership_counts_as_before(tmp_path):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text(
        "##gff-version 3\n" + _gff_row(5000, 6000, "GOI_G|m2", "AMBIGUOUS", "not_specific"))
    out = tmp_path / "report.json"
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "generate_report.py"), "--results_dir", str(stage),
                    "--output", str(out), "--allow-empty"], check=True, capture_output=True)
    summary = json.loads(out.read_text())["summary"]
    assert summary["needs_curation_goi"] == 1 and summary["ambiguous_other_home_gene_goi"] == 0
    assert "another home gene" not in summary["headline"].replace("another home gene fits as well", "")


def test_the_headline():
    h = gr._format_goi_headline(1, 0, 5, 3, ambiguous=2, curation_by_basis={"not_specific": 1, "gene_order_only": 1},
                                other_home_gene=4)
    assert h == ("1 ortholog(s) (1 high confidence); 2 need curation (1 gene order only, 1 another home gene fits as "
                 "well); 4 of another home gene; 5 no significant match; across 3 genome(s).")
    assert gr._format_goi_headline(0, 0, 0, 1, other_home_gene=2) == "NO ortholog; 2 of another home gene; across 1 genome(s)."
    assert gr._format_goi_headline(0, 0, 0, 0) == "No GOI ortholog annotations were produced."
