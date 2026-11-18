"""The strong hit without synteny (bug list items 2 and 3, 2026-11-16): the genome-wide path
(--goi_genome_wide_path) gives a GOI hit outside every searched GOI zone one answer when the
judge's HIGH terms hold in the genome place -- E <= 1e-6 over the genome, >= 0.7 of the query
covered by the hits' HSPs, no HSP across a home exon junction, no real competitor within 10 bits
-- and the model there is judged HIGH, is no retrocopy and not frame-broken without an intron. It
is written as MEDIUM with its own class (genome_wide_strong_hit), never seeds, never enters the
tree, and the report keeps it out of the counted calls: one per genome, only where nothing is
counted. Its decoys are searched in the genome place in decoy mode. The switch
goi_genome_wide_without_blocks is gone: the path runs in every searched genome. All synthetic."""
import json
import logging
import os
import random
import shutil
import subprocess
import sys
import types

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.join(ROOT, "scripts", "benchmark"))

import annotate_goi_exons as age  # noqa: E402
import apply_report_drops as ard  # noqa: E402
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402
import synvoy_label  # noqa: E402
import synvoy_specificity  # noqa: E402

CLASS = synvoy_label.GENOME_WIDE_CLASS
Q = "M" + ("ACDEFGHIKLMNPQRSTVWY" * 15)[:299]                   # a 300-aa home protein


# ====================================================================== the label
def test_the_answer_is_the_judges_high_written_as_medium():
    assert synvoy_label.genome_wide_answer("HIGH", "significant_specific") == (
        "MEDIUM", "significant_specific", CLASS)
    for conf, basis in (("MEDIUM", "significant_specific"), ("AMBIGUOUS", "not_specific"),
                        ("LOW", "not_significant"), (None, None)):
        assert synvoy_label.genome_wide_answer(conf, basis) is None


# ====================================================================== place decisions
def _place(parent="GOI_q", chrom="c1", start=0, e=1e-30, cov=0.9, retro="", real=1, margin=20.0, bits=200.0):
    return {"parent": parent, "chrom": chrom, "start": start, "end": start + 300, "evalue": e,
            "coverage_value": cov, "retrocopy_hits": retro, "real_competitor": real,
            "margin_value": margin, "bits_value": bits, "hits": 1}


def test_a_place_is_tried_only_with_the_judges_high_terms():
    places = [_place(e=1e-5), _place(chrom="c2", real=1, margin=5.0), _place(chrom="c3", real=1, margin=12.0),
              _place(chrom="c4", real=0, margin=-40.0), _place(chrom="c5", cov=0.67),
              _place(chrom="c6", retro=1), _place(chrom="c7", real="", margin=None)]
    tried = isr.choose_genome_wide_places(places)
    assert [p["chrom"] for p in tried] == ["c3", "c4"] or [p["chrom"] for p in tried] == ["c4", "c3"]
    assert {p["chrom"]: p.get("decision") for p in places if p not in tried} == {
        "c1": "evalue", "c2": "specificity", "c5": "coverage", "c6": "retrocopy_hits", "c7": "specificity"}
    assert sorted(p["rank"] for p in tried) == [1, 2]


def test_places_are_ranked_by_e_then_bits_not_by_the_sequence_name():
    """Two places at E 0 went by chromosome name; below 1e-10 E ties, and bits decide."""
    a, b = _place(chrom="a", e=0.0, bits=300.0), _place(chrom="b", e=0.0, bits=400.0)
    c, d = _place(chrom="c", e=1e-12, bits=250.0), _place(chrom="d", e=1e-11, bits=260.0)
    tried = isr.choose_genome_wide_places([a, b, c, d])
    assert [p["chrom"] for p in tried] == ["b", "a", "d", "c"]
    assert [p["rank"] for p in tried] == [1, 2, 3, 4]
    x, y = _place(chrom="x", e=1e-40, bits=300.0, cov=0.8), _place(chrom="y", e=1e-40, bits=300.0, cov=0.95)
    assert [p["chrom"] for p in isr.choose_genome_wide_places([x, y])] == ["y", "x"]


def _evaluate(monkeypatch, hsps, junctions=(), q=Q):
    called = []
    monkeypatch.setattr(isr, "_hit_pieces", lambda genome_seqs, h: ["PIECE" + str(h["start"])])
    monkeypatch.setattr(synvoy_evalue, "blastp_hsps", lambda query, subjects, sp: list(hsps))
    monkeypatch.setattr(isr, "first_best_subject", lambda query, subjects: (321.0, subjects[0]))

    def spec(*a, **k):
        called.append(k.get("own_gene_id"))
        return {"margin": 45.0, "real_competitor": True, "best_other": "gene-p", "in_neighbourhood": False}
    monkeypatch.setattr(synvoy_specificity, "specificity", spec)
    p = {"parent": "GOI_q", "chrom": "c1", "start": 100, "end": 1000,
         "hit_list": [{"start": 100, "end": 1000, "pident": 81.25}], "hits": 1}
    isr.evaluate_genome_wide_place(p, q, {}, 10_000_000, "proteome", "hood", 10, own_gene_id="gene-q",
                                   junctions=list(junctions))
    return p, called


def test_coverage_is_the_union_of_the_significant_hsps(monkeypatch):
    p, called = _evaluate(monkeypatch, [(1e-40, 200.0, 1, 100, 0), (1e-30, 150.0, 200, 300, 0)])
    assert p["coverage"] == "0.670" and isr.choose_genome_wide_places([p]) == []
    assert p["decision"] == "coverage" and called == []           # no specificity below the bar
    p, called = _evaluate(monkeypatch, [(1e-40, 200.0, 1, 150, 0), (1e-30, 150.0, 140, 300, 0)])
    assert p["coverage"] == "1.000" and p["evalue_genome"] == "1e-40" and called == ["gene-q"]
    assert (p["bits"], p["margin"], p["real_competitor"], p["best_other"]) == ("321.0", "45.0", 1, "gene-p")
    assert [x["chrom"] for x in isr.choose_genome_wide_places([p])] == ["c1"]
    # an HSP that chance alone gives once in the genome place does not cover anything
    p, _ = _evaluate(monkeypatch, [(1e-40, 200.0, 1, 150, 0), (5.0, 20.0, 151, 300, 0)])
    assert p["coverage"] == "0.500"
    # the reversed control of the bug list: strong E over a sliver of the query
    p, _ = _evaluate(monkeypatch, [(2e-94, 300.0, 10, 40, 0)])
    assert float(p["coverage"]) < 0.2 and isr.choose_genome_wide_places([p]) == []


def test_with_core_completeness_the_core_coverage_counts(monkeypatch):
    monkeypatch.setattr(isr, "CORE_COMPLETENESS", True)
    monkeypatch.setattr(isr, "core_coverage", lambda spans, seq: 0.95)
    p, _ = _evaluate(monkeypatch, [(1e-40, 200.0, 30, 210, 0)])
    assert p["coverage"] == "0.950"
    monkeypatch.setattr(isr, "CORE_COMPLETENESS", False)
    p, _ = _evaluate(monkeypatch, [(1e-40, 200.0, 30, 210, 0)])
    assert p["coverage"] == "0.603"


def test_home_exon_junctions_are_the_last_residues_of_all_exons_but_the_last():
    assert isr.home_exon_junctions([17, 39, 73]) == [17, 56]
    assert isr.home_exon_junctions([120]) == [] and isr.home_exon_junctions([]) == []


def test_a_piece_spans_a_junction_with_ten_residues_on_each_side():
    junctions = [17, 56]
    assert isr.spans_a_junction([(5, 30)], junctions)                   # 8..17 and 18..27 aligned
    assert not isr.spans_a_junction([(9, 30)], junctions)               # only 9 residues before
    assert not isr.spans_a_junction([(5, 26)], junctions)               # only 9 residues after
    assert not isr.spans_a_junction([(1, 17), (18, 56), (57, 129)], junctions)   # one per exon
    assert isr.spans_a_junction([(1, 17), (40, 80)], junctions)         # exons 2-3 contiguous
    assert not isr.spans_a_junction([(5, 30)], [])                      # single-exon home gene


def test_an_hsp_across_a_home_junction_is_a_retrocopy_place(monkeypatch):
    junctions = isr.home_exon_junctions([100, 80, 120])           # at 100 and 180
    p, called = _evaluate(monkeypatch, [(1e-80, 500.0, 5, 295, 0)], junctions)
    assert p["retrocopy_hits"] == 1 and called == []
    isr.choose_genome_wide_places([p])
    assert p["decision"] == "retrocopy_hits"
    p, _ = _evaluate(monkeypatch, [(1e-30, 180.0, 3, 100, 0), (1e-25, 150.0, 101, 180, 1),
                                   (1e-30, 190.0, 181, 298, 2)], junctions)
    assert p["retrocopy_hits"] == 0 and isr.choose_genome_wide_places([p]) == [p]
    p, _ = _evaluate(monkeypatch, [(1e-80, 500.0, 5, 295, 0)], ())
    assert p["retrocopy_hits"] == "" and isr.choose_genome_wide_places([p]) == [p]


# ====================================================================== the tries
def _pass_args(tmp_path, **over):
    a = types.SimpleNamespace(genome_ctx={"genome_len": 10_000_000, "block_esyn": {}, "goi_cores": [],
                                          "flank_median": 72.25},
                              gap_search_window=1000, goi_home_span=0, expansion=1.0,
                              goi_fallback_max_gap=0, max_intron=5000, home_goi_ids={"GOI_q"},
                              home_goi_gene_id="gene-q", n_flanking_genes=10, output_dir=str(tmp_path),
                              home_goi_junctions=[], max_blocks_per_genome=0)
    for k, v in over.items():
        setattr(a, k, v)
    return a


def _stub_pass(monkeypatch, places, outcomes, bits=None):
    """genome_wide_goi_pass with its steps stubbed: ``places`` as genome_wide_places would give
    them (every one meets the HIGH terms, ``bits`` by start), the window of each place given by
    ``outcomes`` (start -> its window's outcome). Returns (windows searched, the table's rows)."""
    searched, table = [], []
    monkeypatch.setattr(isr, "_home_proteome", lambda args: ("proteome", "hood"))
    monkeypatch.setattr(isr, "select_goi_queries", lambda *a, **k: (
        None, {p["parent"]: {"id": p["parent"], "seq": Q} for p in places}, {}))
    monkeypatch.setattr(isr, "home_goi_parents", lambda args: sorted({p["parent"] for p in places}))
    monkeypatch.setattr(isr, "genome_wide_places",
                        lambda hits, zones, parents, gap, records=(), decoy_parents=None:
                        [dict(p) for p in places if (decoy_parents is None) != isr.is_decoy_parent(p["parent"])])

    def evaluate(p, q, *a, **k):
        p.update(evalue=1e-50, evalue_genome="1e-50", coverage_value=0.95, coverage="0.950",
                 retrocopy_hits="", real_competitor=0, margin_value=30.0, margin="30.0",
                 bits_value=(bits or {}).get(p["start"], 300.0), bits="300.0", best_hit_identity=80.0)
    monkeypatch.setattr(isr, "evaluate_genome_wide_place", evaluate)

    def fake_block(number, block, hits, genome_seqs, db_sequences, genome_name, args, *a, **k):
        searched.append((number, block["start"], block.get("decoy_case")))
        args.genome_ctx.setdefault("genome_wide_outcome", {})[number] = outcomes.get(block["start"], "no_candidate")
        return [], []
    monkeypatch.setattr(isr, "process_region_block", fake_block)
    monkeypatch.setattr(isr, "_write_genome_wide_table", lambda args, prefix, g, ps: table.extend(ps))
    return searched, table


def _gw_place(start, parent="GOI_q", chrom="c1"):
    return {"parent": parent, "chrom": chrom, "start": start, "end": start + 900, "hits": 1,
            "hit_list": [{"pident": 80.0}]}


def test_the_first_answer_ends_the_tries_of_its_home_gene(monkeypatch, tmp_path):
    places = [_gw_place(1000), _gw_place(50_000), _gw_place(90_000)]
    searched, table = _stub_pass(monkeypatch, places, {1000: "judge_not_high", 50_000: "answer"},
                                 bits={1000: 500.0, 50_000: 400.0, 90_000: 300.0})
    args = _pass_args(tmp_path)
    isr.genome_wide_goi_pass(args, [], {"c1": "A" * 200_000}, {}, "G.fna", "u", 1, None, "", 5000)
    assert [s[1] for s in searched] == [1000, 50_000]            # a failed model moves on to the next
    by_start = {p["start"]: p for p in table}
    assert (by_start[1000]["decision"], by_start[1000]["model_decision"]) == ("modelled", "judge_not_high")
    assert (by_start[50_000]["decision"], by_start[50_000]["model_decision"]) == ("modelled", "answer")
    assert by_start[90_000]["decision"] == "not_tried"
    assert {p["decision"] for p in table} <= set(isr.GENOME_WIDE_DECISIONS)
    assert {p["model_decision"] for p in table if p.get("model_decision")} <= set(isr.GENOME_WIDE_MODEL_DECISIONS)
    assert [by_start[s]["block"] for s in (1000, 50_000)] == [5000, 5001]
    assert by_start[1000]["significant_blocks"] == 0 and by_start[1000]["median_flank_identity"] == "72.2"


def test_each_home_gene_gets_its_own_tries(monkeypatch, tmp_path):
    places = [_gw_place(1000, "GOI_a"), _gw_place(60_000, "GOI_a"), _gw_place(1000, "GOI_b", "c2")]
    searched, table = _stub_pass(monkeypatch, places, {1000: "answer"})
    isr.genome_wide_goi_pass(_pass_args(tmp_path), [], {"c1": "A" * 200_000, "c2": "A" * 9_000}, {},
                             "G.fna", "u", 1, None, "", 5000)
    assert sorted((p["parent"], p["start"], p["decision"]) for p in table) == [
        ("GOI_a", 1000, "modelled"), ("GOI_a", 60_000, "not_tried"), ("GOI_b", 1000, "modelled")]


def test_the_tries_stop_at_the_block_cap_and_are_counted(monkeypatch, tmp_path):
    places = [_gw_place(1000), _gw_place(50_000), _gw_place(90_000)]
    searched, table = _stub_pass(monkeypatch, places, {}, bits={1000: 500.0, 50_000: 400.0, 90_000: 300.0})
    args = _pass_args(tmp_path, max_blocks_per_genome=1)
    isr.genome_wide_goi_pass(args, [], {"c1": "A" * 200_000}, {}, "G.fna", "u", 1, None, "", 5000)
    assert [s[1] for s in searched] == [1000]
    assert [p["decision"] for p in sorted(table, key=lambda p: p["start"])] == ["modelled", "not_tried", "not_tried"]
    assert args.genome_ctx["limits"] == {"genome_wide_tries_cap": 2}
    assert "genome_wide_tries_cap" in isr.LIMIT_KEYS


def test_a_place_inside_a_block_window_is_near_that_block(monkeypatch, tmp_path):
    places = [_gw_place(40_000), _gw_place(150_000, chrom="c2")]
    searched, table = _stub_pass(monkeypatch, places, {})
    args = _pass_args(tmp_path)
    args.genome_ctx["block_windows"] = [("c1", 5000, 45_000, 3), ("c1", 30_000, 60_000, 1)]
    isr.genome_wide_goi_pass(args, [], {"c1": "A" * 200_000, "c2": "A" * 200_000}, {}, "G.fna", "u", 1,
                             None, "", 5000)
    assert {p["chrom"]: p["near_block"] for p in table} == {"c1": 1, "c2": ""}
    place = {"chrom": "c1", "start": 100, "end": 200}
    assert isr.near_block_of(place, [("c1", 150, 900, 7), ("c2", 0, 999, 2)]) == 7
    assert isr.near_block_of(place, [("c1", 200, 900, 7)]) == ""


def test_the_pass_runs_without_a_significant_block_and_without_any_block(monkeypatch, tmp_path):
    for ctx in ({"genome_len": 10_000_000, "block_esyn": {0: 0.5}}, {"genome_len": 10_000_000}):
        searched, _table = _stub_pass(monkeypatch, [_gw_place(1000)], {1000: "answer"})
        args = _pass_args(tmp_path)
        args.genome_ctx = dict(ctx)
        isr.genome_wide_goi_pass(args, [], {"c1": "A" * 9000}, {}, "G.fna", "u", 1, None, "", 5000)
        assert [s[1] for s in searched] == [1000]


def test_a_genome_without_blocks_is_searched_on_with_the_path():
    """build_genome_blocks gives no block as [] (searched on, the genome-wide path runs) with the
    path and as None (the search of the genome ends) without it."""
    args = types.SimpleNamespace(max_intron=5000, cluster_distance=10000, region_padding=1000,
                                 min_block_genes=2, max_blocks_per_genome=80)
    opts = dict(home_rank=None, bridge={}, gene_gap=None, goi_seeded=lambda b: False, proxy_parents=set())
    assert isr.build_genome_blocks([], args, "G", gw_any_genome=False, **opts) == (None, 0)
    assert isr.build_genome_blocks([], args, "G", gw_any_genome=True, **opts) == ([], 0)
    src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    genome = src.split("def process_single_genome(", 1)[1].split("\ndef ", 1)[0]
    assert "_gw_any_genome = bool(getattr(args, 'goi_genome_wide_path', False))" in genome
    assert "goi_genome_wide_without_blocks" not in src


def test_rows_written_before_the_pass_hold_no_place():
    rows = ["c1\trearranged_flanking\tmRNA\t7001\t8000\t80.0\t+\t.\tID=gene-A|G_rearranged;SynVoyRole=flanking",
            "c1\trearranged_flanking\tCDS\t7001\t8000\t.\t+\t0\tID=gene-A|G_rearranged_CDS1;Parent=gene-A|G_rearranged",
            "c2\tx\tgene\t1\t300\t.\t+\t.\tID=GOI_q|G_b0_l1_copy"]
    assert isr.written_rows(rows) == [("c1", 7000, 8000), ("c2", 0, 300)]
    hit = {"query": "GOI_q", "chrom": "c1", "start": 7500, "end": 7800, "pident": 90.0}
    assert isr.genome_wide_places([hit], [], {"GOI_q"}, 1000, records=isr.written_rows(rows)) == []
    assert len(isr.genome_wide_places([hit], [], {"GOI_q"}, 1000)) == 1


# ====================================================================== the window's answer
def _cand(start, end, label="HIGH", basis="significant_specific", extra="", cds=None, rank=None,
          evidence="exon_annotation"):
    rows = [f"c1\t{evidence}\tmRNA\t{start}\t{end}\t80.0\t+\t.\tID=GOI_q|G_b5000_l1_{start};"
            f"SynVoy_Parent=GOI_q;SynVoyRole=goi;EvidenceType={evidence};GenomeWide=1{extra}"]
    for i, (s, e) in enumerate(cds or [(start, end)], 1):
        rows.append(f"c1\t{evidence}\tCDS\t{s}\t{e}\t.\t+\t0\tID=x_CDS{i};Parent=GOI_q|G_b5000_l1_{start}")
    return {"start": start, "end": end, "gff": rows, "window_label": label, "window_basis": basis,
            "rank_key": rank if rank is not None else (0, start)}


def test_the_window_keeps_its_best_answer_only():
    best, second = _cand(1001, 2000, rank=(0,)), _cand(5001, 6000, rank=(1,))
    kept, dropped = isr.genome_wide_keep([second, best])
    assert kept == [best] and dropped == [(second, "second_in_window")]


def test_what_is_no_answer():
    medium = _cand(1001, 2000, label="MEDIUM")
    retro = _cand(3001, 4000, extra=";RetrocopyCandidate=1")
    broken_no_intron = _cand(5001, 6000, extra=";FrameBroken=1", cds=[(5001, 5500), (5511, 6000)])
    broken_intron = _cand(7001, 8000, extra=";FrameBroken=1", cds=[(7001, 7300), (7501, 8000)])
    for c, reason in ((medium, "judge_not_high"), (retro, "retrocopy_cds"),
                      (broken_no_intron, "frame_broken_no_intron")):
        assert isr.genome_wide_drop_reason(c) == reason
    assert isr.genome_wide_drop_reason(broken_intron) is None
    assert isr.genome_wide_drop_reason(_cand(1001, 2000, extra=";RetrocopyCandidate=0")) is None
    index = isr.intervals_by_chrom([("c1", 1500, 1600)])
    assert isr.genome_wide_drop_reason(_cand(1001, 2000), index) == "overlaps_record"
    assert isr.genome_wide_drop_reason(_cand(1601, 2000), index) is None
    kept, dropped = isr.genome_wide_keep([medium, retro, broken_no_intron, broken_intron])
    assert kept == [broken_intron] and [r for _c, r in dropped] == [
        "judge_not_high", "retrocopy_cds", "frame_broken_no_intron"]
    assert {r for _c, r in dropped} <= set(isr.GENOME_WIDE_MODEL_DECISIONS)
    assert isr.genome_wide_keep([]) == ([], [])


def test_the_marks_of_a_genome_wide_candidate():
    junctions = [100, 180]
    hits_by_gene = {"GOI_q": [{"qstart": 5, "qend": 295, "gstart": 0, "gend": 870}]}
    c = {"start": 1001, "end": 1870, "gff": ["c1\tx\tmRNA\t1001\t1870\t.\t+\t.\tID=m;SynVoy_Parent=GOI_q"]}
    isr.mark_genome_wide_candidates([c], hits_by_gene, 1000, junctions,
                                    {"coverage": 0.9123, "neighbour_median": 68.04, "near_block": 4})
    a = isr._parse_gff_attr_string(c["gff"][0].split("\t")[8])
    assert (a["GenomeWide"], a["RetrocopyCandidate"], a["GenomeWideCoverage"],
            a["NeighbourMedianIdentity"], a["NearBlock"]) == ("1", "1", "0.912", "68.0", "4")
    exon_by_exon = {"GOI_q": [{"qstart": 3, "qend": 100, "gstart": 0, "gend": 294},
                              {"qstart": 101, "qend": 180, "gstart": 400, "gend": 640}]}
    c = {"start": 1001, "end": 1640, "gff": ["c1\tx\tmRNA\t1001\t1640\t.\t+\t.\tID=m;SynVoy_Parent=GOI_q"]}
    isr.mark_genome_wide_candidates([c], exon_by_exon, 1000, junctions,
                                    {"coverage": 0.6, "neighbour_median": None, "near_block": ""})
    a = isr._parse_gff_attr_string(c["gff"][0].split("\t")[8])
    assert a["RetrocopyCandidate"] == "0" and "NeighbourMedianIdentity" not in a and "NearBlock" not in a
    c = {"start": 1001, "end": 1640, "gff": ["c1\tx\tmRNA\t1001\t1640\t.\t+\t.\tID=m;SynVoy_Parent=GOI_q"]}
    isr.mark_genome_wide_candidates([c], hits_by_gene, 1000, [], {})
    assert c["gff"][0].endswith(";GenomeWide=1")               # untestable without junctions


# ====================================================================== the window itself
class _Args(types.SimpleNamespace):
    def __getattr__(self, name):
        return None


def _window_args(tmp_path, **over):
    a = _Args(output_dir=str(tmp_path), region_padding=5_000, adaptive_padding=False, max_intron=1_000,
              goi_fallback_max_gap=1_000, goi_margin=5_000, gap_search_window=500, goi_home_span=0,
              expansion=1.0, fallback_short_query_len=150, fallback_short_min_aln_aa=15,
              fallback_short_min_bits=30.0, home_rank=None, goi_ranks={}, home_geometry={},
              home_goi_ids={"GOI_x"},
              genome_ctx={"genome_len": 100_000_000, "places": {}, "block_esyn": {}, "decoy_lines": [],
                          "goi_cores": []})
    for k, v in over.items():
        setattr(a, k, v)
    return a


def _stub_window(monkeypatch, window_hits, label="HIGH", seen=None):
    def search(region, queries, *a, **k):
        if seen is not None:
            seen.append([q["id"] for q in queries])
        return [dict(h, query=queries[0]["id"]) for h in window_hits]
    monkeypatch.setattr(isr, "run_augmented_search", search)
    monkeypatch.setattr(isr, "annotate_models_from_hit_list", lambda *a, **k: [])
    monkeypatch.setattr(age, "detect_tandem_duplications", lambda *a, **k: (False, []))
    for name in ("prepare_label_proteins", "prepare_label_subjects", "score_goi_candidates",
                 "annotate_place_evalues", "annotate_model_extras", "annotate_goi_specificity"):
        monkeypatch.setattr(isr, name, lambda *a, **k: None)

    def judge(cands):
        for i, c in enumerate(cands):
            c["window_label"], c["window_basis"], c["rank_key"] = label, "significant_specific", (i,)
    monkeypatch.setattr(isr, "judge_candidates", judge)


def _hit(qstart, qend, start=1_000):
    return {"query": "GOI_x", "chrom": "region_seq", "start": start, "end": start + 3 * (qend - qstart + 1),
            "strand": "+", "qstart": qstart, "qend": qend, "evalue": 1e-60, "pident": 85.0,
            "alnlen": 3 * (qend - qstart + 1), "bits": 400.0}


def _gw_block(**over):
    b = {"chrom": "c1", "start": 41_000, "end": 41_840, "genes": ["GOI_x"], "genes_count": 1, "loci_count": 1,
         "genome_wide": True, "window": (40_000, 43_000), "margin": 1_000, "records_index": {},
         "place_attrs": {"coverage": 0.93, "neighbour_median": 70.0, "near_block": 2}}
    b.update(over)
    return b


DB = {"GOI_x": {"id": "GOI_x", "seq": "M" * 300, "header": "GOI_x"}}


def test_a_genome_wide_window_writes_its_answer_with_the_marks(tmp_path, monkeypatch):
    _stub_window(monkeypatch, [_hit(1, 280)])
    args = _window_args(tmp_path)
    _genes, gff = isr.process_region_block(5000, _gw_block(), [], {"c1": "ATG" * 30_000}, DB, "G.fna",
                                           args, "u", 1)
    rows = [ln for ln in gff if ln.split("\t")[2] in ("mRNA", "gene") and "SynVoyRole=goi" in ln]
    assert len(rows) == 1
    a = isr._parse_gff_attr_string(rows[0].split("\t")[8])
    assert (a["GenomeWide"], a["GenomeWideCoverage"], a["NeighbourMedianIdentity"], a["NearBlock"]) == (
        "1", "0.930", "70.0", "2")
    assert args.genome_ctx["genome_wide_outcome"] == {5000: "answer"}
    assert args.genome_ctx["goi_cores"] == [("c1", 40_000, 43_000)]       # a searched zone as a whole
    assert "block_windows" not in args.genome_ctx and 5000 not in args.genome_ctx.get("block_context", {})


def test_a_genome_wide_window_without_an_answer_writes_nothing(tmp_path, monkeypatch):
    _stub_window(monkeypatch, [_hit(1, 280)], label="MEDIUM")
    args = _window_args(tmp_path)
    _genes, gff = isr.process_region_block(5000, _gw_block(), [], {"c1": "ATG" * 30_000}, DB, "G.fna",
                                           args, "u", 1)
    assert gff == []
    assert args.genome_ctx["genome_wide_outcome"] == {5000: "judge_not_high"}
    assert args.genome_ctx["genome_wide_dropped"] == {"judge_not_high": 1}


def test_a_window_of_the_decoy_pass_searches_its_decoy_alone(tmp_path, monkeypatch):
    seen = []
    _stub_window(monkeypatch, [_hit(1, 280)], seen=seen)
    args = _window_args(tmp_path, decoy_mode=True)
    decoy = {"id": "GOI_x~rev", "seq": ("M" * 300)[::-1]}
    block = _gw_block(decoy_case="GOI_x~rev", decoy_query=decoy)
    _genes, gff = isr.process_region_block(5001, block, [], {"c1": "ATG" * 30_000}, DB, "G.fna", args, "u", 1)
    assert gff == [] and seen == [["GOI_x~rev"]]
    lines = args.genome_ctx["decoy_lines"]
    assert len([ln for ln in lines if ln.split("\t")[2] in ("mRNA", "gene")]) == 1
    assert "SynVoy_Parent=GOI_x~rev" in lines[0] and ";GenomeWide=1" in lines[0]
    assert args.genome_ctx["goi_cores"] == []             # the real search's zones are not its
    assert args.genome_ctx["genome_wide_outcome"] == {5001: "answer"}
    assert "genome_wide_dropped" not in args.genome_ctx


# ====================================================================== the genome-end label
def _row(extra, eeff="1e-40", evidence="exon_annotation"):
    return (f"c1\t{evidence}\tmRNA\t1\t900\t80.0\t+\t.\tID=GOI_q|G_b5000_l1;SynVoy_Parent=GOI_q;"
            f"SynVoyRole=goi;EvidenceType={evidence};Confidence=LOW;GOIClass=ambiguous_goi_family_member;"
            f"InferenceReason=awaiting_label;QueryCoverage=0.950;Place=genome;BlockEsyn=inf;"
            f"SpecMargin=40.0;RealCompetitor=1;EPlace={eeff};Eeff={eeff}{extra}")


def test_the_genome_end_label_of_a_genome_wide_row(caplog):
    out, counts = isr.relabel_goi_lines([_row(";GenomeWide=1"), _row(""),
                                         _row(";GenomeWide=1", evidence="tandem_copy")])
    attrs = [isr._parse_gff_attr_string(ln.split("\t")[8]) for ln in out]
    assert [(a["Confidence"], a["GOIClass"], a["Basis"]) for a in attrs] == [
        ("MEDIUM", CLASS, "significant_specific"), ("HIGH", "confident_goi", "significant_specific"),
        ("MEDIUM", CLASS, "significant_specific")]
    assert counts == {"MEDIUM/genome_wide": 2, "HIGH/significant_specific": 1}
    with caplog.at_level(logging.WARNING, logger=isr.logger.name):
        out, counts = isr.relabel_goi_lines([_row(";GenomeWide=1", eeff="1e-4")])
    a = isr._parse_gff_attr_string(out[0].split("\t")[8])
    assert (a["Confidence"], a["GOIClass"]) == ("LOW", "ambiguous_goi_family_member")
    assert "not the HIGH its window kept it for" in caplog.text


# ====================================================================== seeds and tree
def test_an_answer_never_seeds_and_never_feeds_the_tree():
    genes = [{"id": "GOI_q|G_b5000_l1", "seq": "M" * 300}, {"id": "GOI_q|G_b1_l1", "seq": "M" * 280}]
    answer = {"role": "goi", "parent": "GOI_q", "confidence": "MEDIUM", "goi_class": CLASS,
              "genome_wide": True, "real_competitor": "0"}
    main = {"role": "goi", "parent": "GOI_q", "confidence": "MEDIUM", "goi_class": "probable_goi",
            "real_competitor": "0"}
    seed, tree, suppressed = isr._classify_goi_for_seed_and_tree(
        genes[:1], {genes[0]["id"]: answer}, isr.is_goi_query_id)
    assert (seed, tree, suppressed) == ([], [], 1)
    seed, tree, suppressed = isr._classify_goi_for_seed_and_tree(
        genes, {genes[0]["id"]: answer, genes[1]["id"]: main}, isr.is_goi_query_id)
    assert [g["id"] for g in seed] == ["GOI_q|G_b1_l1"] and tree == [] and suppressed == 1
    # the class alone is enough (a row whose flag a reader lost)
    seed, _tree, _n = isr._classify_goi_for_seed_and_tree(
        genes[:1], {genes[0]["id"]: dict(answer, genome_wide=False)}, isr.is_goi_query_id)
    assert seed == []


# ====================================================================== decoys in the genome place
def _m8(rows):
    return "".join("\t".join(map(str, r)) + "\n" for r in rows)


def test_decoy_rows_leave_the_genome_table_and_the_rest_stays_byte_for_byte(tmp_path):
    real = [("GOI_q", "c1", 88.7, 753, 28, 0, 1, 251, 5604, 6356, "4.061E-152", 466),
            ("gene-A", "c1", 100.0, 603, 0, 0, 1, 201, 3001, 3603, "3.570E-149", 454),
            ("GOI_q|exon_2", "c2", 74.6, 300, 63, 0, 3, 103, 5747, 5448, "7.245E-22", 78)]
    decoy = [("GOI_q~rev", "c1", 27.3, 249, 57, 0, 66, 144, 5802, 6050, "7.857E-05", 36),
             ("GOI_q~shuf2", "c2", 30.0, 90, 20, 0, 5, 35, 100, 10, "0.004", 31)]
    m8 = tmp_path / "G.m8"
    m8.write_text(_m8([real[0], decoy[0], real[1], decoy[1], real[2]]))
    assert isr.split_decoy_rows(str(m8), str(tmp_path / "decoy.m8")) == 2
    assert m8.read_bytes() == _m8(real).encode()
    assert (tmp_path / "decoy.m8").read_text() == _m8(decoy)


def test_the_query_set_with_decoys_is_the_query_set_and_the_decoys(tmp_path):
    db = tmp_path / "q.faa"
    db.write_text(">GOI_q some header\nMKTAYIAK\n>gene-A\nMSSSS")           # no newline at the end
    decoys = [{"id": "GOI_q~rev", "seq": "KAIYATKM"}]
    out = isr.query_db_with_decoys(str(db), decoys, str(tmp_path / "qd.faa"))
    assert open(out).read() == ">GOI_q some header\nMKTAYIAK\n>gene-A\nMSSSS\n>GOI_q~rev\nKAIYATKM\n"
    db.write_text(">GOI_q\nMKTAYIAK\n")
    assert open(isr.query_db_with_decoys(str(db), decoys, out)).read() == ">GOI_q\nMKTAYIAK\n>GOI_q~rev\nKAIYATKM\n"


def test_the_genome_search_decoys_are_those_of_the_home_protein():
    db = {"GOI_q": {"id": "GOI_q", "seq": "MKTAYIAKQRQISFVKSHFSRQ"},
          "GOI_q|exon_1": {"id": "GOI_q|exon_1", "seq": "MKTAYIAK"},
          "GOI_q|Sp_fna_b0_l1_exon_ann": {"id": "GOI_q|Sp_fna_b0_l1_exon_ann", "seq": "MKTAYIAKQRQISFVKSHFSRA"},
          "gene-A": {"id": "gene-A", "seq": "MSSSSSSS"}}
    decoys = isr.genome_search_decoys(db, {"GOI_q", "GOI_q|exon_1"})
    assert [d["id"] for d in decoys] == ["GOI_q~rev", "GOI_q~shuf1", "GOI_q~shuf2"]
    assert decoys[0]["seq"] == db["GOI_q"]["seq"][::-1]


def _codon_table():
    import itertools
    from sequence_utils import translate
    out = {}
    for c in ("".join(p) for p in itertools.product("ACGT", repeat=3)):
        aa = translate(c)
        if aa != "*":
            out.setdefault(aa, c)
    return out


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")
def test_with_mmseqs_the_real_rows_are_the_same_with_and_without_the_decoys(tmp_path):
    """The decoys ride in the genome search's MMseqs2 call; the real rows of hits/<genome>.m8 must
    be the rows of the same search without them, byte for byte (MMseqs2 searches each query on its
    own). A small synthetic genome: three neighbours, two diverged copies of the GOI."""
    rng = random.Random(11)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    codon = _codon_table()

    def protein(n):
        return "M" + "".join(rng.choice(aa) for _ in range(n))

    def back(p):
        return "".join(codon[a] for a in p)

    def dna(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def mutate(p, rate):
        return "".join(rng.choice(aa) if rng.random() < rate else a for a in p)
    goi, flanks = protein(250), [protein(200) for _ in range(3)]
    chrom1 = dna(3000) + back(flanks[0]) + dna(2000) + back(mutate(goi, 0.2)) + dna(2000) + back(flanks[1]) + dna(3000)
    chrom2 = (dna(5000) + isr.reverse_complement(back(mutate(goi, 0.35))) + dna(4000)
              + back(mutate(flanks[2], 0.1)) + dna(2000))
    genome = tmp_path / "genome.fna"
    genome.write_text(f">chr1\n{chrom1}\n>chr2\n{chrom2}\n")
    db = tmp_path / "q.faa"
    db.write_text("".join(f">{i}\n{s}\n" for i, s in [("GOI_q", goi), ("gene-A", flanks[0]),
                                                       ("gene-B", flanks[1]), ("gene-C", flanks[2])]))
    decoys = isr.genome_search_decoys({"GOI_q": {"id": "GOI_q", "seq": goi}}, {"GOI_q"})
    with_decoys = isr.query_db_with_decoys(str(db), decoys, str(tmp_path / "qd.faa"))
    args = types.SimpleNamespace(mmseqs_sens=9.5, mmseqs_split_memory_limit="2G", evalue=0.01,
                                 mmseqs_verbosity=1)
    plain, mixed = tmp_path / "plain.m8", tmp_path / "mixed.m8"
    assert isr.run_mmseqs_easy_search_with_retries(str(db), str(genome), str(plain), str(tmp_path / "t1"),
                                                   args, 1, "G")[0]
    assert isr.run_mmseqs_easy_search_with_retries(with_decoys, str(genome), str(mixed),
                                                   str(tmp_path / "t2"), args, 1, "G")[0]
    isr.split_decoy_rows(str(mixed), str(tmp_path / "decoy.m8"))
    assert plain.read_bytes()
    assert mixed.read_bytes() == plain.read_bytes()
    assert all(isr.is_decoy_parent(ln.split("\t")[0]) for ln in (tmp_path / "decoy.m8").read_text().splitlines())


def _decoy_row(parent, conf, genome_wide=False, cls="probable_goi"):
    return (f"c1\tx\tmRNA\t1\t900\t.\t+\t.\tID={parent}|G_b1_l1;SynVoy_Parent={parent};SynVoyRole=goi;"
            f"Confidence={conf};GOIClass={cls}" + (";GenomeWide=1" if genome_wide else ""))


def test_genome_wide_decoys_are_counted_apart():
    lines = [_decoy_row("GOI_q~rev", "MEDIUM", True, CLASS),                     # an answer of the pass
             _decoy_row("GOI_q~shuf1", "MEDIUM", True, CLASS),
             _decoy_row("GOI_q~shuf1", "MEDIUM"),                                 # counted in a window
             _decoy_row("GOI_q~shuf2", "LOW")]
    assert isr.decoy_gw_counts(lines, ["GOI_q"]) == {"decoy_gw_cases": 3, "decoy_gw_answers": 1}
    assert isr.decoy_gw_counts(lines, []) == {"decoy_gw_cases": 0, "decoy_gw_answers": 0}
    # the window decoy rate does not see the genome-wide rows
    assert isr.decoy_case_counts(lines, ["GOI_q"]) == {"decoy_cases": 3, "decoy_high": 0,
                                                       "decoy_medium_plus": 1, "decoy_ambiguous_plus": 1}


def test_the_decoy_table_and_the_target_policy_columns(tmp_path):
    args = types.SimpleNamespace(output_dir=str(tmp_path), decoy_mode=True, goi_genome_wide_path=True)
    n = isr.write_decoy_outputs(args, "locus_1_", "G.fna", [_decoy_row("GOI_q~rev", "MEDIUM", True, CLASS)])
    rows = (tmp_path / "decoys" / "locus_1_G.fna.decoy.tsv").read_text().splitlines()
    head, row = rows[0].split("\t"), rows[1].split("\t")
    assert n == 1 and dict(zip(head, row))["GenomeWide"] == "1" and head[-2:] == ["GenomeWide", "NearBlock"]
    ctx = {"decoy_gw_counts": {"decoy_gw_cases": 3, "decoy_gw_answers": 0}}
    isr._write_target_policy(args, "locus_1_", "G.fna", 70.0, 5, ctx)
    text = (tmp_path / "hits" / "locus_1_G.fna.target_policy.tsv").read_text().splitlines()
    cols = dict(zip(text[0].split("\t"), text[1].split("\t")))
    assert (cols["decoy_gw_cases"], cols["decoy_gw_answers"]) == ("3", "0")
    args.goi_genome_wide_path = False
    isr._write_target_policy(args, "locus_1_", "G.fna", 70.0, 5, ctx)
    assert "decoy_gw_cases" not in (tmp_path / "hits" / "locus_1_G.fna.target_policy.tsv").read_text()


def test_the_decoy_pass_tries_its_cases_like_the_real_pass(monkeypatch, tmp_path):
    places = [_gw_place(1000), _gw_place(9000, "GOI_q~rev"), _gw_place(70_000, "GOI_q~rev"),
              _gw_place(3000, "GOI_q~shuf1", "c2")]
    searched, table = _stub_pass(monkeypatch, places, {9000: "answer", 1000: "no_candidate"})
    args = _pass_args(tmp_path)
    decoys = [{"id": "GOI_q~rev", "seq": "R" * 50}, {"id": "GOI_q~shuf1", "seq": "S" * 50},
              {"id": "GOI_q~shuf2", "seq": "T" * 50}]
    isr.genome_wide_goi_pass(args, [], {"c1": "A" * 200_000, "c2": "A" * 20_000}, {}, "G.fna", "u", 1, None,
                             "", 5000, decoy_hits=[], decoy_query_list=decoys)
    assert [(s[1], s[2]) for s in searched] == [(1000, None), (9000, "GOI_q~rev"), (3000, "GOI_q~shuf1")]
    rows = sorted((p["decoy"], p["parent"], p["start"], p["decision"]) for p in table)
    assert rows == [("", "GOI_q", 1000, "modelled"), ("rev", "GOI_q", 9000, "modelled"),
                    ("rev", "GOI_q", 70_000, "not_tried"), ("shuf1", "GOI_q", 3000, "modelled")]
    assert args.genome_ctx["decoy_gw_parents"] == ["GOI_q"]


# ====================================================================== the report
DEMOTED = {"paralog_not_goi", "identity_coverage_decoupled"}


def _record(genome, chrom, start, conf="MEDIUM", cls="probable_goi", answer=False, near="", eeff="1e-40",
            source="locus_1"):
    r = {"genome": genome, "source": source, "chrom": chrom, "start": start, "end": start + 899,
         "strand": "+", "evidence": "exon_annotation", "cds": [], "cds_frames": [], "confidence": conf,
         "goi_class": cls, "identity": "80.0", "query_cov": "0.95", "core_cov": "", "bits": "300.0",
         "own_copy_set_aside": "", "own_copy_bits": "", "target_gene": "", "model_status": "complete",
         "refusal_lifted": "", "mrna_id": f"GOI_q|{genome}_{chrom}_{start}", "basis": "significant_specific",
         "eeff": eeff, "place": "genome", "block_esyn": "inf", "copies_in_slot": "", "slot_rank": "",
         "dispersed_rescue": False, "rescue_rule": "", "rescue_anchors": "", "rescue_anchor_bp": "",
         "frame_broken": False, "genome_wide": answer, "genome_wide_coverage": "0.950" if answer else "",
         "retrocopy_candidate": "0" if answer else "", "neighbour_median_identity": "71.0" if answer else "",
         "near_block": near}
    if answer:
        r["goi_class"] = CLASS
    return r


def test_answers_are_deduped_apart_and_the_other_records_are_unchanged():
    main = [_record("G1", "c1", 1000, "HIGH", "confident_goi"), _record("G2", "c1", 5000)]
    answers = [_record("G3", "c9", 100, answer=True), _record("G3", "c9", 100, answer=True, source="locus_2")]
    without = gr.dedupe_goi_records([dict(r) for r in main])
    with_answers = gr.dedupe_goi_records([dict(r) for r in main] + [dict(r) for r in answers])
    n = len(without["records"])
    assert with_answers["records"][:n] == without["records"]
    assert {k: v for k, v in with_answers.items() if k != "records"} == {
        k: v for k, v in without.items() if k != "records"}
    extra = with_answers["records"][n:]
    assert len(extra) == 1 and extra[0]["genome_wide_strong_hit"] is True
    # the same answer from two loci at equal identity: relabelled, still an answer, still not counted
    assert extra[0]["goi_class"] == "cross_locus_duplicate" and extra[0]["provenance"] == ["locus_1", "locus_2"]
    assert not gr._is_counted(extra[0], DEMOTED)
    assert extra[0]["genome_wide"] is True and extra[0]["neighbour_median_identity"] == "71.0"


def test_the_report_keeps_one_answer_per_genome_only_where_nothing_is_counted():
    recs = [_record("G1", "c1", 1000, "HIGH", "confident_goi"), _record("G1", "c5", 100, answer=True),
            _record("G2", "c1", 1000, "HIGH", "paralog_not_goi"), _record("G2", "c7", 100, answer=True),
            _record("G3", "c1", 100, answer=True, eeff="1e-8"), _record("G3", "c2", 100, answer=True, eeff="1e-60"),
            _record("G4", "c1", 1000, "AMBIGUOUS", "syntenic_candidate_unconfirmed", source="locus_2"),
            _record("G4", "c1", 1500, answer=True),
            _record("G5", "c1", 100, answer=True, near="3")]
    for r in recs:
        if r.get("genome_wide"):
            r["genome_wide_strong_hit"] = True
    dedup = {"records": recs}
    dropped = gr.drop_genome_wide_strong_hits(dedup, DEMOTED)
    reasons = {(d["genome"], d["chrom"]): d["reason"] for d in dedup["genome_wide_strong_hits_dropped_records"]}
    assert reasons == {("G1", "c5"): "genome_has_counted_call", ("G3", "c1"): "not_best_in_genome",
                       ("G4", "c1"): "overlaps_record"}
    assert len(dropped) == 3 and dedup["genome_wide_strong_hits_dropped"] == 3
    kept = gr.genome_wide_answer_records(dedup["records"], DEMOTED)
    assert sorted((r["genome"], r["chrom"]) for r in kept) == [("G2", "c7"), ("G3", "c2"), ("G5", "c1")]
    absence = gr.genome_absence_categories(dedup["records"], {"G6": {"best_block_context": "complete"}}, DEMOTED)
    assert absence["strong_hit_without_synteny"] == ["G2", "G3"]
    assert absence["strong_hit_near_block"] == ["G5"]
    assert absence["absent_context_complete"] == ["G6"] and "G2" not in absence["absent_no_block"]


def test_an_answer_demoted_by_ownership_is_no_answer():
    r = _record("G1", "c1", 100, answer=True)
    r.update(genome_wide_strong_hit=True, goi_class="paralog_not_goi")
    dedup = {"records": [r]}
    assert gr.drop_genome_wide_strong_hits(dedup, DEMOTED) == []
    assert dedup["records"] == [r] and gr.genome_wide_answer_records(dedup["records"], DEMOTED) == []


def test_an_answer_does_not_drop_the_gene_order_candidate_of_its_genome():
    goo = _record("G1", "c1", 1000, "AMBIGUOUS", "syntenic_candidate_unconfirmed")
    goo["basis"] = "gene_order_only"
    answer = dict(_record("G1", "c9", 100, answer=True), genome_wide_strong_hit=True)
    dedup = {"records": [goo, answer]}
    assert gr.drop_gene_order_only(dedup, DEMOTED) == []
    assert dedup["post_dedup_records"] == 1 and dedup["ambiguous_goi_deduped"] == 1


def test_the_headline_names_the_answers_apart():
    text = gr._format_goi_headline(2, 1, 4, 5, genome_wide=1, near_block=2)
    assert text == ("3 ortholog(s) (2 high, 1 medium confidence); 1 strong hit(s) without synteny (not "
                    "counted); 2 strong hit(s) near a block, outside its GOI zone (not counted); 4 no "
                    "significant match; across 5 genome(s).")
    assert gr._format_goi_headline(0, 0, 0, 0, genome_wide=1).startswith("NO ortholog; 1 strong hit(s)")
    assert gr._format_goi_headline(0, 0, 0, 0) == "No GOI ortholog annotations were produced."


def _report_dir(tmp_path, genomes, policy=None, decoy_gw=None):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir(exist_ok=True)
    for genome, rows in genomes.items():
        lines = []
        for mid, chrom, start, conf, cls, extra in rows:
            lines.append(f"{chrom}\tx\tmRNA\t{start}\t{start + 899}\t.\t+\t.\tID={mid};SynVoy_Parent=GOI_q;"
                         f"SynVoyRole=goi;Confidence={conf};GOIClass={cls};EvidenceType=exon_annotation;"
                         f"Identity=80.0;QueryCoverage=0.95;Basis=significant_specific;Eeff=1e-40;"
                         f"Place=genome;BlockEsyn=inf{extra}")
        (tmp_path / "regions" / f"locus_1__{genome}.fna.gff").write_text(
            "##gff-version 3\n" + "\n".join(lines) + "\n")
    for name, row in (policy or {}).items():
        cols = list(row)
        (tmp_path / "hits" / f"locus_1_{name}.target_policy.tsv").write_text(
            "\t".join(cols) + "\n" + "\t".join(str(row[c]) for c in cols) + "\n")
    return gr.build_report(str(tmp_path))


ANSWER = ";GenomeWide=1;GenomeWideCoverage=0.950;NeighbourMedianIdentity=71.0;RetrocopyCandidate=0"


def test_the_report_counts_no_answer_and_lists_the_kept_ones(tmp_path):
    main = {"G1": [("GOI_q|G1_b0_l1", "c1", 1001, "HIGH", "confident_goi", "")],
            "G2": [("GOI_q|G2_b0_l1", "c1", 1001, "MEDIUM", "probable_goi", "")]}
    (tmp_path / "plain").mkdir()
    plain = _report_dir(tmp_path / "plain", main)
    genomes = {g: list(rows) for g, rows in main.items()}
    genomes["G1"].append(("GOI_q|G1_b5000_l1", "c8", 1001, "MEDIUM", CLASS, ANSWER))       # dropped
    genomes["G3"] = [("GOI_q|G3_b5000_l1", "c8", 1001, "MEDIUM", CLASS, ANSWER)]            # kept
    genomes["G4"] = [("GOI_q|G4_b5001_l1", "c2", 2001, "MEDIUM", CLASS, ANSWER + ";NearBlock=0")]
    (tmp_path / "with").mkdir()
    report = _report_dir(tmp_path / "with", genomes)
    s, p = report["summary"], plain["summary"]
    for key in ("high_confidence_goi", "medium_confidence_goi", "headline_metric",
                "high_confidence_goi_pre_dedup", "medium_confidence_goi_pre_dedup", "total_goi_annotations",
                "genomes_with_goi_annotations", "low_confidence_goi", "needs_curation_goi"):
        assert s[key] == p[key], key
    assert (s["genome_wide_strong_hits"], s["genome_wide_strong_hit_genomes"]) == (1, ["G3"])
    assert (s["near_block_strong_hits"], s["near_block_strong_hit_genomes"]) == (1, ["G4"])
    assert s["genome_wide_strong_hits_dropped"] == 1
    assert s["strong_hit_without_synteny"] == ["G3"] and s["strong_hit_near_block"] == ["G4"]
    assert "G3" not in s["goi_absent_genomes"] and "G4" not in s["goi_absent_genomes"]
    assert s["headline"] == ("2 ortholog(s) (1 high, 1 medium confidence); 1 strong hit(s) without synteny "
                             "(not counted); 1 strong hit(s) near a block, outside its GOI zone (not "
                             "counted); across 2 genome(s).")
    assert report["annotations"]["genome_wide_answer_rows"] == 3
    assert any("hold a strong hit without synteny" in w for w in report["warnings"])
    assert any("near a block, outside its GOI zone" in w for w in report["warnings"])
    assert not any("strong hit" in w for w in plain["warnings"])
    answers = [r for r in report["goi_dedup"]["records"] if r.get("genome_wide_strong_hit")]
    assert sorted(r["genome"] for r in answers) == ["G3", "G4"]
    main_records = [r for r in report["goi_dedup"]["records"] if not r.get("genome_wide_strong_hit")]
    assert main_records == plain["goi_dedup"]["records"]
    assert report["goi_dedup"]["genome_wide_strong_hits_dropped_records"][0]["reason"] == "genome_has_counted_call"


def test_the_decoy_rate_of_the_genome_place(tmp_path):
    row = {"genome": "G.fna", "decoy_cases": 3, "decoy_high": 0, "decoy_medium_plus": 0,
           "decoy_ambiguous_plus": 0, "decoy_gw_cases": 3, "decoy_gw_answers": 0}
    old = {k: v for k, v in row.items() if not k.startswith("decoy_gw")}
    for name, r in (("G.fna", row), ("H.fna", dict(row, genome="H.fna")), ("K.fna", dict(old, genome="K.fna"))):
        (tmp_path / f"locus_1_{name}.target_policy.tsv").write_text(
            "\t".join(r) + "\n" + "\t".join(str(v) for v in r.values()) + "\n")
    rate = gr.collect_decoy_rate(sorted(str(p) for p in tmp_path.glob("*.tsv")))
    assert rate["cases"] == 9
    assert rate["genome_wide"]["cases"] == 6 and rate["genome_wide"]["genome_searches"] == 2
    assert rate["genome_wide"]["answers"]["calls"] == 0 and 0 < rate["genome_wide"]["answers"]["ci95_upper"] < 0.5
    (tmp_path / "locus_1_G.fna.target_policy.tsv").unlink()
    (tmp_path / "locus_1_H.fna.target_policy.tsv").unlink()
    assert "genome_wide" not in gr.collect_decoy_rate([str(p) for p in tmp_path.glob("*.tsv")])


def test_the_figures_leave_out_a_dropped_answer(tmp_path):
    rows = ["c8\tx\tmRNA\t1001\t1900\t.\t+\t.\tID=GOI_q|G1_b5000_l1;SynVoyRole=goi;Confidence=MEDIUM;"
            f"GOIClass={CLASS};GenomeWide=1",
            "c8\tx\tCDS\t1001\t1900\t.\t+\t0\tID=GOI_q|G1_b5000_l1_CDS1;Parent=GOI_q|G1_b5000_l1",
            "c8\tx\tmRNA\t1001\t1900\t.\t+\t.\tID=GOI_q|G1_b0_l1;SynVoyRole=goi;Confidence=HIGH;"
            "GOIClass=confident_goi"]
    gff = tmp_path / "G1.fna.gff"
    gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    (tmp_path / "r.json").write_text(json.dumps({"goi_dedup": {
        "gene_order_only_dropped_records": [],
        "genome_wide_strong_hits_dropped_records": [{"genome": "G1", "chrom": "c8", "start": 1001, "end": 1900,
                                                     "mrna_id": "GOI_q|G1_b5000_l1", "provenance": ["locus_1"],
                                                     "reason": "genome_has_counted_call"}]}}))
    spans = ard.dropped_spans(str(tmp_path / "r.json"))
    assert spans == {"G1": [("c8", 1001, 1900, CLASS)]}
    left_out = ard.filter_gff(str(gff), spans["G1"], str(tmp_path / "o.gff"))
    assert left_out == ["GOI_q|G1_b5000_l1"]
    assert (tmp_path / "o.gff").read_text().splitlines()[1:] == rows[2:]


def test_the_family_scorer_never_counts_an_answer(tmp_path):
    import family_recovery as fr
    recs = [dict(_record("Gx", "c1", 100), mrna_id="m1"),
            dict(_record("Gx", "c2", 100, answer=True), mrna_id="m2", genome_wide_strong_hit=True)]
    (tmp_path / "synvoy_report.json").write_text(json.dumps({"goi_dedup": {"records": recs}}))
    loaded = fr.load_report_records(str(tmp_path))
    assert [(r["id"], r["final"], r["genome_wide_strong_hit"]) for r in loaded] == [
        ("m1", True, False), ("m2", False, True)]


# ====================================================================== blastp_hsps
@pytest.mark.skipif(synvoy_evalue._blastp() is None, reason="blastp not on PATH")
def test_blastp_evalue_is_the_smallest_of_blastp_hsps():
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    q = "".join(rng.choice(aa) for _ in range(300))
    subjects = ["".join(rng.choice(aa) for _ in range(40)), q[150:300], q[1:120]]
    space = synvoy_evalue.search_space(len(q), 12_000_000)
    hsps = synvoy_evalue.blastp_hsps(q, subjects, space)
    assert {h[4] for h in hsps} == {1, 2}
    assert synvoy_evalue.blastp_evalue(q, subjects, space) == min(
        ((e, b) for e, b, *_rest in hsps), key=lambda x: (x[0], -x[1]))
    by_subject = {h[4]: h for h in hsps}
    assert (by_subject[1][2], by_subject[1][3]) == (151, 300) and (by_subject[2][2], by_subject[2][3]) == (2, 120)
    cov = isr._query_union_coverage([{"qstart": h[2], "qend": h[3]} for h in hsps], len(q))
    assert round(cov, 3) == 0.897
