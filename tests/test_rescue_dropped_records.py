"""The GOI hits the dispersed-GOI rescue dropped leave a record (bug list item 10, 2026-11-16).

Before, a GOI hit that failed the rescue's quality bar was only counted, and in a genome without any
flanking envelope not even that: the outputs could not say which test a hit failed, or whether it was
seen. Now every such hit outside every searched GOI zone is a row of hits/<locus>_<genome>.nonsyntenic.tsv
(outcome dropped, rejected_by the test), the strongest NONSYNTENIC_DROPPED_MAX_ROWS per genome search,
with the total in target_policy.tsv (rescue_dropped_outside_zones); a hit the rescue left to a GOI zone
whose block was not searched is a dropped row too (dead_zone).
"""
import csv
import os
import sys
import types

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402


def _goi(chrom, start, pident=60.0, evalue=1e-20, alnlen=80, bits=100.0, qend=104):
    return {"query": "GOI_q", "chrom": chrom, "start": start, "end": start + 300, "strand": "-",
            "qstart": 5, "qend": qend, "pident": pident, "evalue": evalue, "alnlen": alnlen, "bits": bits}


def _anchor(gene, chrom, start):
    return {"query": gene, "chrom": chrom, "start": start, "end": start + 1000, "strand": "+",
            "pident": 90.0, "evalue": 1e-50, "alnlen": 300, "bits": 500.0}


# The third: 20 aligned query residues (5-24), below the rescue's 50; MMseqs2's alnlen of such a
# hit is its target span in nucleotides (60), which the length test does not read.
WEAK = [_goi("c7", 9_000_000, pident=30.0), _goi("c7", 9_100_000, evalue=1e-3),
        _goi("c7", 9_200_000, alnlen=60, qend=24)]


def test_the_branch_without_any_envelope_counts_its_drops():
    """One anchor in the genome: no sequence has an envelope. The three weak hits used to be
    neither counted nor recorded."""
    drops, dropped = {}, []
    seeds = isr.select_dispersed_goi_seeds(WEAK + [_anchor("gene-A", "c1", 1_000)], [_anchor("gene-A", "c1", 1_000)],
                                           drops=drops, dropped_out=dropped)
    assert seeds == []
    assert drops == {"dispersed_identity": 1, "dispersed_evalue": 1, "dispersed_alnlen": 1}
    assert [(r["rejected_by"], r["outcome"], r["envelope"]) for r in dropped] == [
        ("dispersed_identity", "dropped", "none_in_genome"), ("dispersed_evalue", "dropped", "none_in_genome"),
        ("dispersed_alnlen", "dropped", "none_in_genome")]


def test_a_dropped_record_carries_the_hit_and_its_place():
    anchors = [_anchor(f"gene-a{i}", "c1", 1_000_000 + i * 100_000) for i in range(4)]
    zone = [("c1", 990_000, 1_011_000)]
    hits = [_goi("c1", 1_150_000, pident=30.0), _goi("c1", 1_005_000)] + anchors
    dropped, in_zone = [], []
    isr.select_dispersed_goi_seeds(hits, anchors, goi_zones=zone, region_padding=10_000,
                                   envelope_margin=500_000, dropped_out=dropped, in_zone_out=in_zone)
    weak, near = dropped
    assert (weak["rejected_by"], weak["envelope"], weak["strand"], weak["qstart"], weak["qend"], weak["bits"]) == (
        "dispersed_identity", "inside", "-", 5, 104, 100.0)
    assert (weak["nearest_anchor_gene"], weak["nearest_anchor_bp"], weak["anchors_on_sequence"]) == (
        "gene-a1", 49_000, 4)
    assert near["rejected_by"] == "dispersed_in_goi_zone" and in_zone[0]["start"] == 1_005_000
    assert set(isr.NONSYNTENIC_COLUMNS) >= set(weak) - {"locus", "modelled_here", "modelled_by"}


def test_the_table_lists_the_drops_outside_every_searched_zone_strongest_first():
    ctx = {"rescue_refused": [{"chrom": "c9", "start": 5, "end": 10, "outcome": "refused",
                               "rejected_by": "chromosome_lacks_min_anchors"}],
           "rescue_dropped": [
               {"chrom": "c1", "start": 100, "end": 400, "bits": 50.0, "evalue": 1e-5, "query": "GOI_q",
                "rejected_by": "dispersed_identity"},                                   # in a searched zone
               {"chrom": "c2", "start": 100, "end": 400, "bits": 40.0, "evalue": 1e-5, "query": "GOI_q",
                "rejected_by": "dispersed_identity"},
               {"chrom": "c2", "start": 900, "end": 1200, "bits": 90.0, "evalue": 1e-12, "query": "GOI_q",
                "rejected_by": "dispersed_in_goi_zone"},                                # its block not searched
               {"chrom": "c3", "start": 100, "end": 400, "bits": 40.0, "evalue": 1e-9, "query": "GOI_q",
                "rejected_by": "dispersed_alnlen"}]}
    cores = [("c1", 0, 1_000)]
    rows, n = isr.rescue_table_rows(ctx, cores)
    assert n == 3
    assert [(r["chrom"], r["rejected_by"]) for r in rows] == [
        ("c9", "chromosome_lacks_min_anchors"), ("c2", "dead_zone"), ("c3", "dispersed_alnlen"),
        ("c2", "dispersed_identity")]
    rows, n = isr.rescue_table_rows(ctx, cores, max_dropped=1)
    assert n == 3 and [r["chrom"] for r in rows] == ["c9", "c2"]
    assert ctx["rescue_dropped"][2]["rejected_by"] == "dispersed_in_goi_zone"     # not changed in place


def _args(tmp_path, prefix="locus_2"):
    return types.SimpleNamespace(output_dir=str(tmp_path), prefix=prefix, decoy_mode=False)


ROW = {"genome": "G.fna", "query": "GOI_q", "chrom": "c1", "start": 1_000, "end": 1_300, "strand": "+",
       "qstart": 1, "qend": 100, "bits": 120.0, "identity": 70.0, "evalue": 1e-30, "alnlen": 300,
       "rejected_by": "outside_flanking_envelope", "nearest_anchor_gene": "gene-A", "nearest_anchor_bp": 5_000,
       "nearest_anchor_same_chrom": True, "outcome": "refused", "envelope": "outside",
       "anchors_on_sequence": 3, "genes_on_sequence": 4}


def test_the_writer_fills_the_locus_and_the_model_and_is_reproducible(tmp_path):
    gff = ["c1\tx\tmRNA\t1201\t2000\t.\t+\t.\tID=GOI_q|G_fna_b0_l1;SynVoyRole=goi;Confidence=MEDIUM",
           "c1\tx\tmRNA\t901\t1100\t.\t+\t.\tID=GOI_q|G_fna_b0_l2;SynVoyRole=goi;Confidence=LOW",
           "c1\tx\tmRNA\t1001\t1300\t.\t+\t.\tID=gene-A|G_fna_b0_fl1;SynVoyRole=flanking;Confidence=HIGH"]
    far = dict(ROW, start=50_000, end=50_300)
    path = tmp_path / "hits" / "locus_2_G.fna.nonsyntenic.tsv"
    isr._write_nonsyntenic_candidates(_args(tmp_path), "locus_2_", "G.fna", [ROW, far], gff_lines=gff)
    first = path.read_bytes()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == isr.NONSYNTENIC_COLUMNS
    assert [(r["locus"], r["modelled_here"], r["modelled_by"]) for r in rows] == [
        ("locus_2", "MEDIUM", "GOI_q|G_fna_b0_l1"), ("locus_2", "", "")]
    isr._write_nonsyntenic_candidates(_args(tmp_path), "locus_2_", "G.fna", [ROW, far], gff_lines=gff)
    assert path.read_bytes() == first
    # before the models (the first write): no model columns; a table left without rows goes
    isr._write_nonsyntenic_candidates(_args(tmp_path), "locus_2_", "G.fna", [ROW])
    with open(path, newline="") as fh:
        assert [r["modelled_here"] for r in csv.DictReader(fh, delimiter="\t")] == [""]
    isr._write_nonsyntenic_candidates(_args(tmp_path), "locus_2_", "G.fna", [])
    assert not path.exists()


def test_the_target_policy_records_the_total(tmp_path):
    args = _args(tmp_path)
    isr._write_target_policy(args, "locus_2_", "G.fna", 60.0, 9, {"rescue_dropped_outside_zones": 345})
    with open(tmp_path / "hits" / "locus_2_G.fna.target_policy.tsv", newline="") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert row["rescue_dropped_outside_zones"] == "345"
    path = str(tmp_path / "hits" / "locus_2_G.fna.target_policy.tsv")
    assert gr.sum_target_policy_column([path, path], "rescue_dropped_outside_zones") == 690
    assert gr.sum_target_policy_column([path], "no_such_column") is None
