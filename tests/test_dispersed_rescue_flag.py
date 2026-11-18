"""Calls found through the dispersed-GOI rescue carry a flag (bug list items 13 and 8, 2026-11-16).

A GOI model of a block the rescue seeded that lies outside every GOI zone of the blocks built
without the rescue is a model only the rescue made the search find: its feature row carries
DispersedRescue=1, RescueRule (envelope, near_neighbour, or no_flanking_hits in a genome without
any flanking hit), RescueAnchors and RescueAnchorBp (iterative_search_runner.mark_rescue_blocks,
mark_dispersed_rescue). The block table names the rule, the decoy table carries the flag, and the
report counts the counted calls with it (summary.dispersed_rescue_goi, dispersed_rescue_by_rule)
and lists the genomes without any flanking hit (summary.genomes_without_flanking_hits).
"""
import csv
import os
import re
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402

RESCUE = {"rule": "near_neighbour", "anchors": 1, "anchor_bp": 200}


def _cand(mid, start, end, role="goi"):
    attrs = f"ID={mid};SynVoy_Parent=GOI_q;SynVoyRole={role};Confidence=LOW;EvidenceType=exon_annotation"
    return {"start": start, "end": end,
            "gff": [f"ctg7\texon_annotation\tmRNA\t{start}\t{end}\t.\t+\t.\t{attrs}",
                    f"ctg7\texon_annotation\tCDS\t{start}\t{end}\t.\t+\t0\tID={mid}.cds1;Parent={mid}"]}


def test_the_flag_goes_on_the_feature_row_of_models_outside_the_first_pass_zones():
    outside = _cand("GOI_q|g_b1_l1_exon_ann", 51_201, 51_800)
    inside = _cand("GOI_q|g_b1_l2_exon_ann", 40_001, 40_600)
    flanking = _cand("gene-E|g_b1_fl1_flank_ann", 50_001, 51_000, role="flanking")
    zones = [("ctg7", 30_000, 45_000), ("chr1", 0, 10**9)]
    assert isr.mark_dispersed_rescue([outside, inside, flanking], "ctg7", RESCUE, zones) == 1
    assert outside["gff"][0].endswith(";DispersedRescue=1;RescueRule=near_neighbour;RescueAnchors=1;"
                                      "RescueAnchorBp=200")
    assert "DispersedRescue" not in outside["gff"][1]                    # the CDS row is left as it is
    assert "DispersedRescue" not in inside["gff"][0] and "DispersedRescue" not in flanking["gff"][0]
    # a block the rescue did not seed: nothing is flagged
    other = _cand("GOI_q|g_b0_l1_exon_ann", 51_201, 51_800)
    assert isr.mark_dispersed_rescue([other], "ctg7", None, []) == 0 and "Dispersed" not in other["gff"][0]
    # a genome without any flanking hit: every model, no distance
    nf = _cand("GOI_q|g_b0_l1_exon_ann", 51_201, 51_800)
    isr.mark_dispersed_rescue([nf], "ctg7", dict(isr.NO_FLANKING_RESCUE), [])
    assert nf["gff"][0].endswith(";DispersedRescue=1;RescueRule=no_flanking_hits;RescueAnchors=0")


def test_decoy_candidates_are_flagged_and_their_table_carries_it(tmp_path):
    decoy = _cand("GOI_q~rev|g_b1_l1_exon_ann", 51_201, 51_800)
    decoy["gff"][0] = decoy["gff"][0].replace("SynVoy_Parent=GOI_q", "SynVoy_Parent=GOI_q~rev")
    isr.mark_dispersed_rescue([decoy], "ctg7", RESCUE, [])
    args = types.SimpleNamespace(output_dir=str(tmp_path))
    assert isr.write_decoy_outputs(args, "", "g", decoy["gff"]) == 1
    with open(tmp_path / "decoys" / "g.decoy.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert (rows[0]["DispersedRescue"], rows[0]["RescueRule"]) == ("1", "near_neighbour")


def test_the_relabel_keeps_the_flag():
    line = ("ctg7\tx\tmRNA\t51201\t51800\t.\t+\t.\tID=GOI_q|g_b1_l1;SynVoy_Parent=GOI_q;SynVoyRole=goi;"
            "Confidence=LOW;GOIClass=x;InferenceReason=awaiting_label;Place=genome;EPlace=1e-30;"
            "DispersedRescue=1;RescueRule=near_neighbour;RescueAnchors=1;RescueAnchorBp=200;Eeff=1e-30;"
            "QueryCoverage=0.9;SpecMargin=50;RealCompetitor=0")
    out, _ = isr.relabel_goi_lines([line])
    a = isr._parse_gff_attr_string(out[0].split("\t")[8])
    assert (a["DispersedRescue"], a["RescueRule"], a["RescueAnchors"], a["RescueAnchorBp"]) == (
        "1", "near_neighbour", "1", "200")


def test_a_block_takes_the_rule_of_its_first_seed_by_position():
    seeds = [dict(chrom="c1", start=5_000, end=5_300, query="GOI_q", rescue_rule="envelope",
                  rescue_anchors=4, rescue_anchor_bp=3_000),
             dict(chrom="c1", start=1_000, end=1_300, query="GOI_q", rescue_rule="near_neighbour",
                  rescue_anchors=1, rescue_anchor_bp=100)]
    blocks = [{"chrom": "c1", "start": 900, "end": 6_000, "genes": ["GOI_q", "gene-E"]},
              {"chrom": "c2", "start": 900, "end": 6_000, "genes": ["gene-A", "gene-B"]}]
    assert isr.mark_rescue_blocks(blocks, seeds) == 1
    assert blocks[0]["dispersed_rescue"] == {"rule": "near_neighbour", "anchors": 1, "anchor_bp": 100}
    assert "dispersed_rescue" not in blocks[1]


def test_the_block_table_names_the_rule(tmp_path):
    args = types.SimpleNamespace(output_dir=str(tmp_path), genome_ctx={"block_esyn": {0: 1.0, 1: 1.0}})
    isr.record_block_row(args, 0, {"chrom": "c1", "start": 0, "end": 10, "genes": ["gene-A"],
                                   "dispersed_rescue": RESCUE}, 0, 10, 1, None)
    isr.record_block_row(args, 1, {"chrom": "c2", "start": 0, "end": 10, "genes": ["gene-B"]}, 0, 10, 1, None)
    isr._write_block_table(args, "", "g.fna", [], args.genome_ctx)
    with open(tmp_path / "blocks" / "g.fna.blocks.tsv", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == isr.BLOCK_TABLE_COLUMNS and reader.fieldnames[-1] == "dispersed_rescue"
    assert [r["dispersed_rescue"] for r in rows] == ["near_neighbour", ""]


def test_the_search_flags_the_candidates_of_a_rescued_block():
    """process_region_block flags the real and the decoy candidates of a rescued block (not a
    window of the genome-wide path), and a model the leftover join rebuilds."""
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    body = src.split("def process_region_block(", 1)[1].split("\ndef ", 1)[0]
    assert "rescue_info = None if genome_wide else block.get('dispersed_rescue')" in body
    for target in ("raw_candidates", "decoy_candidates", "[joined]"):
        assert f"mark_dispersed_rescue({target}, chrom, rescue_info, pre_rescue_zones)" in body
    genome = src.split("def process_single_genome(", 1)[1].split("\ndef ", 1)[0]
    assert "_b['dispersed_rescue'] = dict(NO_FLANKING_RESCUE)" in genome
    assert re.search(r"logger\.warning\(\s*f\"\[\{genome_name\}\] No flanking gene was found", genome)


def _report_dir(tmp_path, rows, policy=None):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir(exist_ok=True)
    lines = []
    for mid, chrom, start, conf, cls, ident, qcov, extra in rows:
        lines.append(f"{chrom}\tx\tmRNA\t{start}\t{start + 899}\t.\t+\t.\tID={mid};SynVoy_Parent=GOI_q;"
                     f"SynVoyRole=goi;Confidence={conf};GOIClass={cls};EvidenceType=exon_annotation;"
                     f"Identity={ident};QueryCoverage={qcov}{extra}")
    (tmp_path / "regions" / "locus_1__G.fna.gff").write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
    for name, row in (policy or {}).items():
        cols = list(row)
        (tmp_path / "hits" / f"locus_1_{name}.target_policy.tsv").write_text(
            "\t".join(cols) + "\n" + "\t".join(str(row[c]) for c in cols) + "\n")
    return gr.build_report(str(tmp_path))


FLAG = ";DispersedRescue=1;RescueRule=near_neighbour;RescueAnchors=1;RescueAnchorBp=200"


def test_the_report_counts_the_counted_calls_found_by_the_rescue(tmp_path):
    report = _report_dir(tmp_path, [
        ("GOI_q|G_fna_b1_l1", "c1", 1_001, "HIGH", "confident_goi", 80, 0.9, FLAG),
        # demoted by the coverage rule (high identity over a fifth of the query): not counted
        ("GOI_q|G_fna_b2_l1", "c2", 1_001, "MEDIUM", "probable_goi", 90, 0.2,
         FLAG.replace("near_neighbour", "envelope")),
        ("GOI_q|G_fna_b0_l1", "c3", 1_001, "HIGH", "confident_goi", 80, 0.9, "")])
    records = {r["chrom"]: r for r in report["goi_dedup"]["records"]}
    assert (records["c1"]["dispersed_rescue"], records["c1"]["rescue_rule"], records["c1"]["rescue_anchors"],
            records["c1"]["rescue_anchor_bp"]) == (True, "near_neighbour", "1", "200")
    assert records["c2"]["goi_class"] == "identity_coverage_decoupled"
    assert "dispersed_rescue" not in records["c3"]
    assert report["summary"]["dispersed_rescue_goi"] == 1
    assert report["summary"]["dispersed_rescue_by_rule"] == {"near_neighbour": 1}


def test_the_report_names_the_genomes_without_flanking_hits(tmp_path):
    base = {"genome": "", "median_flank_identity": "NA", "n_flanking_genes": 0, "blocks_searched": 3}
    report = _report_dir(tmp_path, [("GOI_q|G_fna_b0_l1", "c1", 1_001, "HIGH", "confident_goi", 80, 0.9,
                                     ";DispersedRescue=1;RescueRule=no_flanking_hits;RescueAnchors=0")],
                         policy={"G.fna": dict(base, genome="G.fna"),
                                 "H.fna": dict(base, genome="H.fna", blocks_searched=0),
                                 "K.fna": dict(base, genome="K.fna", n_flanking_genes=7)})
    assert report["summary"]["genomes_without_flanking_hits"] == ["G"]
    assert any("no hit of any flanking gene" in w and "G" in w for w in report["warnings"])
    assert report["summary"]["dispersed_rescue_by_rule"] == {"no_flanking_hits": 1}
    other = tmp_path / "clean"
    other.mkdir()
    clean = _report_dir(other, [("GOI_q|G_fna_b0_l1", "c1", 1_001, "HIGH", "confident_goi", 80, 0.9, "")])
    assert clean["summary"]["genomes_without_flanking_hits"] == []
    assert not any("flanking gene" in w for w in clean["warnings"])


def _two_locus_report(tmp_path, flags):
    """One gene of genome G modelled at the same coordinates from two home loci; ``flags`` says
    which locus's model carries DispersedRescue=1."""
    for sub in ("regions", "hits", "scores"):
        (tmp_path / sub).mkdir(exist_ok=True)
    for locus, flagged in flags:
        (tmp_path / "regions" / f"{locus}__G.fna.gff").write_text(
            "##gff-version 3\n"
            "c1\tx\tmRNA\t901\t1500\t.\t+\t.\tID=GOI_Def|G_fna_b0_l1;SynVoy_Parent=GOI_Def;"
            "SynVoyRole=goi;Confidence=HIGH;GOIClass=probable_goi;EvidenceType=exon_annotation;"
            f"Identity=70;QueryCoverage=0.9{FLAG if flagged else ''}\n")
    return gr.build_report(str(tmp_path))


def test_a_gene_one_locus_found_in_its_own_block_is_no_rescue_born_call(tmp_path):
    """Review fix (2026-11-17): the dedup record took the flag from its representative alone, so a
    gene found through the rescue by one home locus and in a normal block by another was counted as
    rescue-born, or not, depending on which locus came first. Now the flag needs every member."""
    for i, flags in enumerate(((("locus_1", True), ("locus_2", False)),
                               (("locus_1", False), ("locus_2", True)))):
        d = tmp_path / str(i)
        d.mkdir()
        report = _two_locus_report(d, flags)
        (rec,) = report["goi_dedup"]["records"]
        assert rec["provenance"] == ["locus_1", "locus_2"]
        assert "dispersed_rescue" not in rec and report["summary"]["dispersed_rescue_goi"] == 0
    d = tmp_path / "both"
    d.mkdir()
    report = _two_locus_report(d, (("locus_1", True), ("locus_2", True)))
    (rec,) = report["goi_dedup"]["records"]
    assert rec["dispersed_rescue"] is True and report["summary"]["dispersed_rescue_goi"] == 1
