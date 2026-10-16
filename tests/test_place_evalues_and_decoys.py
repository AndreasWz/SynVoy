"""Places, E-values in their place and the decoy mode of the search (label rework, step 3,
2026-10-14; iterative_search_runner + bin/synvoy_evalue.py).

A kept GOI candidate gets Place (slot / core / window / genome), BlockEsyn, EPlace and
EPlaceTblastn; at the end of a genome Eeff and EeffTblastn = those times the number of
significant blocks (SigBlocks). In decoy mode the reversed home protein and two shuffles
of it are searched in every window; their models are chosen and judged apart and written
to decoys/, never to the result GFF.
"""
import math
import os
import random
import sys
import types

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue as ev  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {a: c for a, c in zip(AA, ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT",
                                   "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"])}

def _tabular(cmd):
    """A tblastn or blast_formatter call that writes tabular output (-outfmt 6 ...)."""
    return ("-outfmt" in cmd and str(cmd[cmd.index("-outfmt") + 1]).startswith("6")
            and os.path.basename(str(cmd[0])) in ("tblastn", "blast_formatter"))



def _protein(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice(AA) for _ in range(n))


def _cand(mid, start, end, seq, parent="GOI_q", score=50.0, role="goi"):
    attrs = (f"ID={parent}|g_b0_{mid};Name={parent};SynVoy_Parent={parent};SynVoyRole={role};"
             f"Confidence=MEDIUM;EvidenceType=exon_annotation")
    return {"start": start, "end": end, "score": score, "record": {"id": f"{parent}|g_b0_{mid}", "seq": seq},
            "gff": [f"chr1\texon_annotation\tmRNA\t{start}\t{end}\t{score}\t+\t.\t{attrs}"]}


# ---------------------------------------------------------------- decoys
def test_decoys_are_reversed_and_shuffled_and_repeatable():
    q = _protein(80, 1)
    d = isr.decoy_queries([{"id": "GOI_q", "seq": q}, {"id": "GOI_q~rel", "seq": q}])
    # the relative's query gets its own three (RF-3, 2026-10-15), filed under ~<kind>~rel
    assert [x["id"] for x in d] == ["GOI_q~rev", "GOI_q~shuf1", "GOI_q~shuf2",
                                    "GOI_q~rev~rel", "GOI_q~shuf1~rel", "GOI_q~shuf2~rel"]
    assert d[0]["seq"] == q[::-1]
    assert sorted(d[1]["seq"]) == sorted(q) and d[1]["seq"] != q and d[1]["seq"] != d[2]["seq"]
    assert isr.decoy_queries([{"id": "GOI_q", "seq": q}]) == d[:3]   # same shuffles on a rerun
    assert all(isr.is_decoy_parent(x["id"]) for x in d) and not isr.is_decoy_parent("GOI_q~rel")


def test_greedy_select_keeps_the_better_of_two_overlapping():
    a = dict(_cand("a", 100, 900, "M"), bits=80.0)
    b = dict(_cand("b", 150, 850, "M"), bits=95.0)
    c = dict(_cand("c", 2000, 2600, "M"), bits=10.0)
    assert [x["record"]["id"] for x in isr.greedy_select([a, b, c])] == ["GOI_q|g_b0_b", "GOI_q|g_b0_c"]


# ---------------------------------------------------------------- block terms
def test_block_terms_record_esyn_and_find_the_slot():
    args = types.SimpleNamespace(
        genome_ctx={"places": {f"gene-{i}": 1 for i in range(6)}, "genome_len": 10_000_000,
                    "block_esyn": {}},
        home_rank={f"gene-{i}": i for i in range(6)},
        goi_ranks={"gene-0": -3, "gene-1": -2, "gene-2": -1, "gene-3": 1, "gene-4": 2, "gene-5": 3},
        home_geometry={"pitch": 1000.0}, expansion=1.0)
    flank = [_cand("f2", 20_000, 21_000, "M", parent="gene-2", role="flanking"),
             _cand("f3", 30_000, 31_000, "M", parent="gene-3", role="flanking"),
             _cand("f0", 5_000, 6_000, "M", parent="gene-0", role="flanking")]
    t = isr.block_place_terms(args, 7, 0, 60_000, 4_000, 40_000, 6, flank)
    assert args.genome_ctx["block_esyn"][7] == pytest.approx(t["esyn"])
    assert t["esyn"] == pytest.approx(ev.e_syn(6, [1] * 6, 60_000, 10_000_000))
    assert t["significant"] and t["slot"] == (20_000, 31_000)


def test_without_order_the_block_is_not_significant():
    args = types.SimpleNamespace(genome_ctx={"places": {}, "genome_len": 1_000_000, "block_esyn": {}},
                                 home_rank={}, goi_ranks={}, home_geometry={}, expansion=1.0)
    t = isr.block_place_terms(args, 0, 0, 50_000, 0, 50_000, None, [])
    assert not t["significant"] and t["slot"] is None


# ---------------------------------------------------------------- E-values in the place
def _window_with_gene(q, gene_at=3000, window_len=9000):
    rng = random.Random(9)
    dna = "".join(rng.choice("ACGT") for _ in range(window_len))
    cds = "".join(CODON[a] for a in q)
    return dna[:gene_at] + cds + dna[gene_at + len(cds):], gene_at, gene_at + len(cds)


def test_place_evalues_and_their_multiplier():
    if not ev._blastp():
        pytest.skip("blastp not installed")
    q = _protein(120, 4)
    window, g0, g1 = _window_with_gene(q)
    w_start = 100_000
    cand = _cand("l1", w_start + g0 + 1, w_start + g1, q)
    hits = {"GOI_q": [{"gstart": g0, "gend": g1, "strand": "+", "method": "mmseqs", "evalue": 1e-30}]}
    # tblastn's own rows, as it reported them for the window (before keep floor and dedup)
    tblastn_raw = [{"query": "GOI_q", "start": g0, "end": g1, "evalue": 1e-20},
                   {"query": "GOI_q~rev", "start": g0, "end": g1, "evalue": 1e-40},
                   {"query": "GOI_q", "start": g0 - 3000, "end": g0 - 2000, "evalue": 1e-50}]
    terms = {"esyn": 1e-4, "significant": True, "slot": (w_start + g0 - 500, w_start + g1 + 500),
             "core": (w_start, w_start + 9000), "window": (w_start, w_start + 9000),
             "genome_len": 12_000_000, "window_len": 9000}
    isr.annotate_place_evalues([cand], {"GOI_q": {"seq": q}}, hits, terms, window, w_start, tblastn_raw)
    a = isr._parse_gff_attr_string(cand["gff"][0].split("\t")[8])
    assert a["Place"] == "slot" and float(a["EPlace"]) < 1e-30
    assert float(a["EPlaceTblastn"]) == pytest.approx(1e-20 * (g1 - g0 + 1000) / 9000, rel=0.01)
    ctx = {"block_esyn": {0: 1e-4, 1: 0.05, 2: 3.0}}
    out = isr.finalize_place_evalues(cand["gff"], ctx)
    b = isr._parse_gff_attr_string(out[0].split("\t")[8])
    assert b["SigBlocks"] == "2"
    assert float(b["Eeff"]) == pytest.approx(2 * float(a["EPlace"]), rel=0.01)


def test_genome_place_is_not_multiplied():
    line = ("chr1\tx\tmRNA\t1\t300\t.\t+\t.\tID=GOI_q|m;SynVoy_Parent=GOI_q;SynVoyRole=goi;"
            "Place=genome;EPlace=0.004;EPlaceTblastn=inf")
    out = isr.finalize_place_evalues([line], {"block_esyn": {0: 1e-5, 1: 1e-6}})
    a = isr._parse_gff_attr_string(out[0].split("\t")[8])
    assert a["Eeff"] == "0.004" and a["EeffTblastn"] == "inf" and a["SigBlocks"] == "2"


def test_decoy_outputs(tmp_path):
    args = types.SimpleNamespace(output_dir=str(tmp_path))
    line = ("chr1\tx\tmRNA\t10\t300\t.\t+\t.\tID=GOI_q~shuf1|g_b0_l1;SynVoy_Parent=GOI_q~shuf1;"
            "SynVoyRole=goi;Confidence=LOW;EvidenceType=fallback_hit_span;Place=core;Eeff=3.2")
    assert isr.write_decoy_outputs(args, "", "g1", [line]) == 1
    rows = (tmp_path / "decoys" / "g1.decoy.tsv").read_text().splitlines()
    assert rows[1].split("\t")[:4] == ["g1", "GOI_q", "shuf1", "GOI_q~shuf1|g_b0_l1"]
    assert "Place" in rows[0] and "\tcore\t" in rows[1]
    assert (tmp_path / "decoys" / "g1.decoy.gff").read_text().count("\n") == 2


def test_decoy_mode_is_wired_and_off_by_default():
    config = open(os.path.join(REPO, "nextflow.config")).read()
    assert "decoy_mode = false" in config
    module = open(os.path.join(REPO, "modules", "iterative_search.nf")).read()
    assert "--decoy_mode ${params.decoy_mode}" in module


def test_anchor_ranks_around_the_gene():
    genes = {"a": {"start": 100, "end": 200}, "b": {"start": 300, "end": 400},
             "c": {"start": 900, "end": 1000}, "d": {"start": 1100, "end": 1200},
             "host": {"start": 450, "end": 800}}
    assert ev.ranks_around_goi(genes, (500, 600)) == {"b": -1, "a": -2, "c": 1, "d": 2}
    assert ev.ranks_around_goi(genes, None) == {}


def test_window_search_hands_over_every_tblastn_row(tmp_path, monkeypatch):
    """The keep floor (E <= 0.1) still applies to the search; the sink gets every row."""
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")

    def fake_run(cmd, *a, **k):
        # tabular output, by tblastn or by blast_formatter from its archive (RF-3, F4)
        if _tabular(cmd):
            with open(cmd[cmd.index("-out") + 1], "w") as fh:
                fh.write("GOI_q\tregion_seq\t40.0\t50\t30\t0\t1\t50\t100\t249\t1e-05\t60.0\n")
                fh.write("GOI_q\tregion_seq\t12.0\t20\t18\t0\t1\t20\t900\t959\t3.5\t20.0\n")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(isr.subprocess, "run", fake_run)
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", lambda **k: (False, "skipped"))
    args = types.SimpleNamespace(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_evalue_mult=1000.0,
                                 aug_relaxed_evalue_cap=10.0, aug_relaxed_parse_evalue_mult=10.0,
                                 aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6,
                                 min_identity=10.0, aug_relaxed_length_min=15, aug_relaxed_length_div=2,
                                 min_length=10, enable_smith_waterman=False, deterministic_goi_search=True)
    sink = []
    hits = isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": _protein(60, 2)}], "g1", args,
                                    "t1", 1, tblastn_sink=sink)
    assert [h["evalue"] for h in hits] == [1e-05]
    assert sorted(h["evalue"] for h in sink) == [1e-05, 3.5]


# ---------------------------------------------------------------- hit regions (2026-10-15)
def test_a_hit_is_clipped_to_the_model_in_its_own_frame():
    assert isr.clip_hit_in_frame(100, 400, '+', 150, 1000) == (151, 400)    # 51 bp = 17 codons in
    assert isr.clip_hit_in_frame(100, 400, '-', 0, 350) == (100, 349)       # read from 400 down
    assert isr.clip_hit_in_frame(100, 400, '+', 120, 300) == (121, 300)
    assert isr.clip_hit_in_frame(100, 400, '+', 390, 900) is None           # < 15 bp left
    assert isr.clip_hit_in_frame(100, 400, '+', 500, 900) is None


def test_a_broken_model_is_judged_on_hit_regions_that_stick_out():
    """miniprot's span is the aligned part; the window's hits often reach a few codons
    beyond it. A model whose CDS has a broken frame has no protein of its own, so these
    hits are all it can be judged on (the one Candida case lost its E-value without them)."""
    if not ev._blastp():
        pytest.skip("blastp not installed")
    q = _protein(150, 11)
    window, g0, g1 = _window_with_gene(q)
    w_start = 50_000
    cand = _cand("l1", w_start + g0 + 31, w_start + g1 - 30, q)       # 10 codons in from each end
    cand["label_protein"] = ""                                         # FrameBroken
    hits = {"GOI_q": [{"gstart": g0, "gend": g1, "strand": "+", "method": "mmseqs", "evalue": 1e-30}]}
    terms = {"esyn": 1e-6, "significant": True, "slot": (w_start, w_start + 9000),
             "core": (w_start, w_start + 9000), "window": (w_start, w_start + 9000),
             "genome_len": 12_000_000, "window_len": 9000}
    isr.annotate_place_evalues([cand], {"GOI_q": {"seq": q}}, hits, terms, window, w_start, [])
    a = isr._parse_gff_attr_string(cand["gff"][0].split("\t")[8])
    assert float(a["EPlace"]) < 1e-30


def test_decoys_only_of_the_home_protein_named_by_its_parent():
    """A seeded model that is a query of its own (a GOI known only from models) gets no
    decoy: named from its record id, ``GOI_copy_1|g_b0_l1~shuf1`` was filed under the real
    parent ``GOI_copy_1`` and its models counted as real (one Candida pilot gene, 2026-10-15)."""
    q = _protein(60, 3)
    queries = [{"id": "GOI_q", "seq": q},
               {"id": "GOI_copy_1|g_b0_l1", "seq": _protein(25, 4)},
               {"id": "GOI_q~rel", "seq": q}]
    ids = [d["id"] for d in isr.decoy_queries(queries)]
    rel = ["GOI_q~rev~rel", "GOI_q~shuf1~rel", "GOI_q~shuf2~rel"]     # the relative's query
    assert ids == ["GOI_q~rev", "GOI_q~shuf1", "GOI_q~shuf2"] + rel
    assert all(isr.extract_base_gene_id(i) == i for i in ids)
    # with the home ids of the run: a home record whose id carries more than its parent
    with_ids = isr.decoy_queries([{"id": "GOI_q|isoform_2", "seq": q}] + queries[1:],
                                 home_ids={"GOI_q|isoform_2"})
    assert [d["id"] for d in with_ids] == ["GOI_q~rev", "GOI_q~shuf1", "GOI_q~shuf2"] + rel
