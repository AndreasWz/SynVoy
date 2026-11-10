"""Fix item 27b (9 Nov 2026): the GOI's name no longer breaks ties.

The judge's sort key (synvoy_label.rank_key, which orders greedy_select) ended in the candidate's id, and the dispersed
seeds took the first hit per bin in query order. So a GOI name that sorts before 'GOI_copy_' (digits, capitals, a, b)
and one that sorts after it (orf19.*) chose differently between a tandem copy and a hit-built candidate at one place
(measured on one Candida gene, 28 rename runs: the copy or the rescued exon, both LOW, same bits). The key now ends in a
tail without the name: the kind of model (exon_annotation, tandem copy, hit chain, rescued exon, raw hit), strand, CDS
rows and the id behind its '|'. The kinds' order keeps the measured choice of names after 'GOI_copy_' (the copy).
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import iterative_search_runner as isr  # noqa: E402

NAMES = ["orf19.6528", "zbc-x", "ABC-x", "12-x", "abc-x", "P00334", "LOC1"]


def _cand(rid, parent, ev, start=764882, end=765106, strand="-"):
    attrs = (f"ID={rid};Name={parent};SynVoy_Parent={parent};Identity=30.0;SynVoyRole=goi;EvidenceType={ev};"
             "Confidence=LOW;Place=slot;BlockEsyn=0.126;EPlace=inf;SpecMargin=5;RealCompetitor=1;"
             "BestOtherInNeighbourhood=0;QueryCoverage=0.4")
    mrna = f"trop_s_7\tsrc\tmRNA\t{start}\t{end}\t30.0\t{strand}\t.\t{attrs}"
    cds = f"trop_s_7\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={rid}.cds1;Parent={rid}"
    return {"start": start, "end": end, "bits": 26.2, "record": {"id": rid, "seq": "M"}, "gff": [mrna, cds]}


@pytest.mark.parametrize("name", NAMES)
def test_a_copy_and_a_rescued_exon_at_one_place_choose_the_same_for_every_name(name):
    goi = f"GOI_{name}"
    copy = _cand(isr.tandem_copy_record_id("GOI_copy_1", goi, "Candida_tropicalis_fna", 1, 1).replace("_l1_", "_r1_"),
                 "GOI_copy_1", "tandem_copy")
    exon = _cand(f"{goi}|Candida_tropicalis_fna_b1_l1_extra1", goi, "rescued_exon")
    isr.judge_candidates([copy, exon])
    kept = isr.greedy_select([exon, copy])
    assert [isr._parse_gff_attr_string(c["gff"][0].split("\t")[8])["EvidenceType"] for c in kept] == ["tandem_copy"]


@pytest.mark.parametrize("name", NAMES)
def test_the_tail_holds_no_name(name):
    goi = f"GOI_{name}"
    copy = _cand(isr.tandem_copy_record_id("GOI_copy_1", goi, "G_fna", 1, 1), "GOI_copy_1", "tandem_copy")
    tail = isr.name_free_tail(copy["gff"], copy["record"]["id"])
    assert name not in repr(tail) and tail[0] == isr.TIE_KIND["tandem_copy"] and tail[3] == "G_fna_b1_l1"


def _hit(query, bits, start=600000, end=600300):
    return {"query": query, "chrom": "chr2", "start": start, "end": end, "strand": "+", "pident": 70.0,
            "evalue": 1e-30, "alnlen": 100, "bits": bits}


@pytest.mark.parametrize("name", NAMES)
def test_a_dispersed_seed_is_the_home_proteins_hit_for_every_name(name):
    """One bin with a 300-bit hit of the home protein and a 150-bit hit of a seeded tandem copy: in query
    order an orf19 name seeded the copy's hit and P00334 the home protein's."""
    goi = f"GOI_{name}"
    flanking = [{"query": f"gene-F{i}", "chrom": "chr2", "start": 100000 + i * 200000, "end": 102000 + i * 200000,
                 "strand": "+", "pident": 80.0, "evalue": 1e-50, "alnlen": 300, "bits": 400.0} for i in range(5)]
    hits = sorted([_hit(goi, 300.0), _hit(f"GOI_copy_1|G_fna_b0_l1_by_{goi}", 150.0, start=600010)],
                  key=lambda h: h["query"])
    seeds = isr.select_dispersed_goi_seeds(hits + flanking, flanking, min_chrom_anchors=3)
    assert [s["query"] for s in seeds] == [goi]


@pytest.mark.parametrize("name", NAMES)
def test_the_block_tables_best_call_does_not_depend_on_the_name(name, tmp_path):
    """The block table's best call (blocks/<genome>.blocks.tsv) is the copy for every name: its key ends in
    the name-free tail of the model's mRNA row, then the model id."""
    import csv
    import types
    goi = f"GOI_{name}"
    copy = _cand(isr.tandem_copy_record_id("GOI_copy_1", goi, "G_fna", 0, 1), "GOI_copy_1", "tandem_copy")["gff"]
    exon = _cand(f"{goi}|G_fna_b0_l1_extra1", goi, "rescued_exon")["gff"]
    ctx = {"block_rows": {0: {"chrom": "trop_s_7", "start": 700000, "end": 800000, "hood_start": 650000,
                              "hood_end": 850000, "flanking_support": 3, "collinear_support": None}},
           "block_esyn": {0: 0.126}}
    isr._write_block_table(types.SimpleNamespace(output_dir=str(tmp_path)), "", "G.fna", gff_lines=exon + copy,
                           ctx=ctx)
    with open(tmp_path / "blocks" / "G.fna.blocks.tsv") as fh:
        (row,) = list(csv.DictReader(fh, delimiter="\t"))
    assert row["best_call_id"] == isr.tandem_copy_record_id("GOI_copy_1", goi, "G_fna", 0, 1)
    assert (row["n_low"], row["best_call"]) == ("2", "LOW")
