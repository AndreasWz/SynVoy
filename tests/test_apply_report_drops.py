"""Review fix RF-3 (report): the per-genome GFFs and the figures leave out the weak matches in
the gene's slot (Basis=gene_order_only) that the report drops because the genome has a
counted call (bin/apply_report_drops.py, run by PLOT_SYNTENY after GENERATE_REPORT)."""
import json
import os
import subprocess
import sys

import pytest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)

import apply_report_drops as ard  # noqa: E402
import generate_report as gr  # noqa: E402

GOO = "SynVoyRole=goi;EvidenceType=fallback_hit_span;Confidence=AMBIGUOUS;Basis=gene_order_only"


def _gff(path, rows):
    path.write_text("##gff-version 3\n##synvoy-distances {}\n" + "".join(r + "\n" for r in rows))


def _rows(chrom="c1"):
    return [
        f"{chrom}\tfallback_hits\tmRNA\t1000\t1900\t30.0\t+\t.\tID=GOI_q|G_b0_l1_fallback;Name=GOI_q;{GOO}",
        f"{chrom}\tfallback_hits\tCDS\t1000\t1900\t.\t+\t0\tID=GOI_q|G_b0_l1_fallback_CDS1;"
        f"Parent=GOI_q|G_b0_l1_fallback",
        f"{chrom}\tminiprot\tmRNA\t1200\t2400\t80.0\t+\t.\tID=GOI_q|G_b1_l1;Name=GOI_q;SynVoyRole=goi;"
        f"EvidenceType=exon_annotation;Confidence=HIGH;Basis=significant_specific",
        f"{chrom}\tminiprot\tCDS\t1200\t2400\t.\t+\t0\tID=GOI_q|G_b1_l1_CDS1;Parent=GOI_q|G_b1_l1",
        f"{chrom}\tfallback_hits\tmRNA\t90000\t90300\t28.0\t-\t.\tID=GOI_q|G_b3_l1_fallback;Name=GOI_q;{GOO}",
        f"{chrom}\tflanking_annotation\tmRNA\t5000\t6000\t70.0\t+\t.\tID=gene-A|G_fl1;SynVoyRole=flanking;"
        f"Confidence=HIGH",
    ]


def _report(path, records):
    path.write_text(json.dumps({"goi_dedup": {"gene_order_only_dropped_records": records}}))


def test_the_dropped_model_and_its_rows_go_nothing_else(tmp_path):
    _gff(tmp_path / "Gx.fna.gff", _rows())
    _gff(tmp_path / "Gy.fna.gff", _rows())        # another genome: the report kept its rows
    _report(tmp_path / "r.json", [{"genome": "Gx", "chrom": "c1", "start": 1000, "end": 1900,
                                   "mrna_id": "GOI_q|G_b0_l1_fallback", "provenance": ["locus_1"]}])
    out = tmp_path / "plot_inputs_x"
    subprocess.run([sys.executable, os.path.join(BIN, "apply_report_drops.py"), "--report",
                    str(tmp_path / "r.json"), "--out_dir", str(out),
                    str(tmp_path / "Gx.fna.gff"), str(tmp_path / "Gy.fna.gff")], check=True)
    kept = (out / "Gx.fna.gff").read_text().splitlines()
    assert kept[:2] == ["##gff-version 3", "##synvoy-distances {}"]
    assert kept[2:] == _rows()[2:]                 # the HIGH model at the same place stays
    assert (out / "Gy.fna.gff").read_text() == (tmp_path / "Gy.fna.gff").read_text()


def test_genome_named_by_accession_as_the_report_names_it(tmp_path):
    name = "GCA_000123456.1_Tgenus_species.fna.gff"
    assert gr.canonical_genome_id("locus_2__" + name) == gr.canonical_genome_id(name) == "GCA_000123456.1"
    _gff(tmp_path / name, _rows())
    _report(tmp_path / "r.json", [{"genome": "GCA_000123456.1", "chrom": "c1", "start": 1500,
                                   "end": 1600, "mrna_id": "x", "provenance": ["locus_2"]}])
    spans = ard.dropped_spans(str(tmp_path / "r.json"))
    left_out = ard.filter_gff(str(tmp_path / name), spans["GCA_000123456.1"], str(tmp_path / "o.gff"))
    assert left_out == ["GOI_q|G_b0_l1_fallback"]


@pytest.mark.parametrize("report", ["NO_REPORT", "old.json", "broken.json"])
def test_without_the_list_the_gffs_are_copied_unchanged(tmp_path, report):
    (tmp_path / "NO_REPORT").write_text("sentinel: no synvoy_report.json\n")
    (tmp_path / "old.json").write_text(json.dumps({"goi_dedup": {"gene_order_only_dropped": 2}}))
    (tmp_path / "broken.json").write_text("{")
    _gff(tmp_path / "Gx.fna.gff", _rows())
    out = tmp_path / "o"
    subprocess.run([sys.executable, os.path.join(BIN, "apply_report_drops.py"), "--report",
                    str(tmp_path / report), "--out_dir", str(out), str(tmp_path / "Gx.fna.gff")],
                   check=True)
    assert (out / "Gx.fna.gff").read_text() == (tmp_path / "Gx.fna.gff").read_text()


def test_plot_module_draws_from_the_filtered_copies():
    nf = open(os.path.join(BIN, "..", "modules", "plot_synteny.nf")).read()
    assert "path report" in nf and "apply_report_drops.py --report $report" in nf
    assert '"--target_gffs \\$target_gffs_in"' in nf
    main = open(os.path.join(BIN, "..", "main.nf")).read()
    assert main.index("ADJUDICATE_AND_REPORT(") < main.index("PLOT_SYNTENY(")
    assert "ADJUDICATE_AND_REPORT.out.report_for_plot" in main
    sub = open(os.path.join(BIN, "..", "subworkflows", "adjudicate_and_report.nf")).read()
    assert "report_for_plot = GENERATE_REPORT.out.report" in sub and ".first()" in sub
