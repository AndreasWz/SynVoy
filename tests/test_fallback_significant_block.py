"""A GOI hit chain that the fallback gate refuses is written in a significant block (bug list
item 7, 2026-11-16).

Where miniprot finds no model, the window's GOI hits are chained, and the chain was refused when
it was short or covered little of the query (fallback_refusal): a small or hidden fragment at the
right place got no model, while a leftover hit beside a model (rescued_exon) had no such gate. In
a block with E_syn <= 0.1, which predicted the gene's place, the refusal is now lifted: the chain
is written with RefusalLifted=<the tests it failed> and the judge's E-value in its place decides
its label (a fragment is never HIGH: HIGH needs 0.7 of the query). Elsewhere the refusal stands.
A lifted call neither seeds the next wave nor enters the tree, and the report counts it.
"""
import csv
import os
import sys
import types

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import annotate_goi_exons as age  # noqa: E402
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402

RANKS = {f"gene-n{i}": i for i in range(1, 5)}


def _ctx_args(home_rank=None, **over):
    a = types.SimpleNamespace(
        home_rank=home_rank, goi_ranks={}, home_geometry={}, expansion=1.0,
        genome_ctx={"genome_len": 100_000_000, "places": {g: 1 for g in (home_rank or {})},
                    "block_esyn": {}})
    for k, v in over.items():
        setattr(a, k, v)
    return a


# ------------------------------------------------------------------ the block's E_syn, early
@pytest.mark.parametrize("home_rank,k,significant", [(RANKS, 4, True), (RANKS, 1, False), (None, None, False)])
def test_block_esyn_is_the_value_block_place_terms_records(home_rank, k, significant):
    a = _ctx_args(home_rank)
    early = isr.block_esyn(a, k, 14_000, 5_000, 19_000)
    terms = isr.block_place_terms(a, 0, 5_000, 19_000, 5_000, 19_000, k, [], hood_len=14_000)
    assert early == terms["esyn"] == a.genome_ctx["block_esyn"][0]
    assert terms["significant"] is significant
    assert (early <= synvoy_evalue.SIGNIFICANT_BLOCK_ESYN) is significant


def test_block_esyn_records_nothing():
    a = _ctx_args(RANKS)
    isr.block_esyn(a, 4, 14_000, 5_000, 19_000)
    assert a.genome_ctx["block_esyn"] == {}


def test_the_decision():
    assert isr.fallback_decision(["long_qcov"], True) == ([], ["long_qcov"])
    assert isr.fallback_decision(["long_qcov", "long_aln"], False) == (["long_qcov", "long_aln"], [])
    assert isr.fallback_decision([], True) == ([], [])
    assert isr.fallback_decision([], False) == ([], [])


def test_the_lift_is_counted_apart():
    assert "fallback_lifted" in isr.LIMIT_KEYS and "limit_fallback_lifted" in isr.LIMIT_COLUMNS
    assert isr.LIMIT_KEYS.index("fallback_lifted") == isr.LIMIT_KEYS.index("fallback_long_bits") + 1


# ------------------------------------------------------------------ the window
class _Args(types.SimpleNamespace):
    def __getattr__(self, name):        # every setting the window does not need here is unset
        return None


def _window_args(tmp_path, home_rank):
    return _Args(output_dir=str(tmp_path), region_padding=5_000, adaptive_padding=False, max_intron=1_000,
                 goi_fallback_max_gap=1_000, goi_margin=5_000, gap_search_window=500, goi_home_span=0,
                 expansion=1.0, fallback_short_query_len=150, fallback_short_min_aln_aa=15,
                 fallback_short_min_bits=30.0, home_rank=home_rank, goi_ranks={}, home_geometry={},
                 genome_ctx={"genome_len": 100_000_000, "places": {g: 1 for g in (home_rank or {})},
                             "block_esyn": {}, "decoy_lines": []})


def _anchor(gene, start):
    return {"query": gene, "chrom": "c1", "start": start, "end": start + 500, "strand": "+", "qstart": 1,
            "qend": 160, "evalue": 1e-50, "pident": 80.0, "alnlen": 480, "bits": 300.0}


def _window_hit(qstart, qend, start=6_000, query="GOI_x"):
    # a window hit, window coordinates (the window starts at 5,000); alnlen as MMseqs2 writes it
    return {"query": query, "chrom": "region_seq", "start": start, "end": start + 3 * (qend - qstart + 1),
            "strand": "+", "qstart": qstart, "qend": qend, "evalue": 1e-6, "pident": 70.0,
            "alnlen": 3 * (qend - qstart + 1), "bits": 45.0}


def _stub_window(monkeypatch, window_hits):
    """The window's search, miniprot and the judge's blastp/Smith-Waterman steps stubbed: miniprot
    finds no model, so the GOI hits go to the hit chain."""
    monkeypatch.setattr(isr, "run_augmented_search", lambda *a, **k: [dict(h) for h in window_hits])
    monkeypatch.setattr(isr, "annotate_models_from_hit_list", lambda *a, **k: [])
    monkeypatch.setattr(age, "detect_tandem_duplications", lambda *a, **k: (False, []))
    for name in ("prepare_label_proteins", "prepare_label_subjects", "score_goi_candidates",
                 "annotate_place_evalues", "annotate_model_extras", "annotate_goi_specificity",
                 "judge_candidates"):
        monkeypatch.setattr(isr, name, lambda *a, **k: None)


def _search_window(tmp_path, monkeypatch, window_hits, home_rank, query_len=300):
    """process_region_block over one block of four anchors (in home order with ``home_rank``),
    stubbed as _stub_window says. Returns (GOI mRNA rows, limit counts, args)."""
    _stub_window(monkeypatch, window_hits)
    args = _window_args(tmp_path, home_rank)
    anchors = [_anchor(g, 10_000 + 1_000 * i) for i, g in enumerate(RANKS)]
    block = {"chrom": "c1", "start": 10_000, "end": 14_000, "genes": list(RANKS), "genes_count": 4}
    db = {"GOI_x": {"id": "GOI_x", "seq": "M" * query_len, "header": "GOI_x"}}
    _genes, gff = isr.process_region_block(0, block, anchors, {"c1": "ATG" * 10_000}, db, "G.fna", args,
                                           "u", 1)
    rows = [ln for ln in gff if ln.split("\t")[2] == "mRNA" and "SynVoyRole=goi" in ln]
    return rows, args.genome_ctx.get("limits", {}), args


def test_a_short_chain_in_a_significant_block_is_written_and_flagged(tmp_path, monkeypatch):
    """20 of 300 query residues in one hit: refused by coverage and length, lifted here."""
    rows, limits, args = _search_window(tmp_path, monkeypatch, [_window_hit(101, 120)], RANKS)
    assert len(rows) == 1
    a = isr._parse_gff_attr_string(rows[0].split("\t")[8])
    assert a["RefusalLifted"] == "long_qcov,long_aln"
    assert (a["EvidenceType"], a["ModelStatus"], a["QueryCoverage"]) == ("fallback_hit_span", "fragment", "0.067")
    assert limits == {"fallback_lifted": 1}
    assert args.genome_ctx["block_esyn"][0] <= synvoy_evalue.SIGNIFICANT_BLOCK_ESYN


def test_the_same_chain_in_a_block_without_home_order_is_refused(tmp_path, monkeypatch):
    rows, limits, args = _search_window(tmp_path, monkeypatch, [_window_hit(101, 120)], None)
    assert rows == []
    assert limits == {"fallback_refused": 1, "fallback_long_qcov": 1, "fallback_long_aln": 1}
    assert args.genome_ctx["block_esyn"][0] > synvoy_evalue.SIGNIFICANT_BLOCK_ESYN


def test_a_chain_the_gate_accepts_carries_no_flag(tmp_path, monkeypatch):
    rows, limits, _args = _search_window(tmp_path, monkeypatch, [_window_hit(1, 200)], RANKS)
    assert len(rows) == 1 and "RefusalLifted" not in rows[0]
    assert limits == {}


def test_a_decoys_lift_is_flagged_but_not_counted(tmp_path, monkeypatch):
    """A decoy query's chain takes the same path (decoy.tsv shows the false-call rate of lifted
    calls); the limit counts leave decoys out."""
    monkeypatch.setattr(isr, "decoy_queries", lambda queries, home_ids=None: [
        {"id": "GOI_x~rev", "seq": "M" * 300}])
    _stub_window(monkeypatch, [_window_hit(101, 120, query="GOI_x~rev")])
    args = _window_args(tmp_path, RANKS)
    args.decoy_mode = True
    anchors = [_anchor(g, 10_000 + 1_000 * i) for i, g in enumerate(RANKS)]
    block = {"chrom": "c1", "start": 10_000, "end": 14_000, "genes": list(RANKS), "genes_count": 4}
    db = {"GOI_x": {"id": "GOI_x", "seq": "M" * 300, "header": "GOI_x"}}
    isr.process_region_block(0, block, anchors, {"c1": "ATG" * 10_000}, db, "G.fna", args, "u", 1)
    lines = args.genome_ctx["decoy_lines"]
    mrna = [ln for ln in lines if ln.split("\t")[2] == "mRNA"]
    assert len(mrna) == 1 and "RefusalLifted=long_qcov,long_aln" in mrna[0]
    assert "fallback_lifted" not in args.genome_ctx.get("limits", {})
    out = types.SimpleNamespace(output_dir=str(tmp_path))
    assert isr.write_decoy_outputs(out, "", "G.fna", lines) == 1
    with open(tmp_path / "decoys" / "G.fna.decoy.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["RefusalLifted"] == "long_qcov,long_aln"


def test_a_relatives_lifted_chain_the_rebase_drops_is_not_counted(tmp_path, monkeypatch):
    """Review fix (2026-11-17): the lift was counted when the chain was built, before the rebase
    drops a relative's (~rel) chain beside a model of the home protein; such a chain never reached
    the judge, so fallback_lifted counted more than the RefusalLifted rows. It is counted after the
    rebase now."""
    home_seq, rel_seq = "M" * 300, "MW" * 145
    hits = [_window_hit(1, 280), _window_hit(101, 120, start=6_300, query="GOI_x~rel")]
    _stub_window(monkeypatch, hits)
    monkeypatch.setattr(isr, "_relative_score", lambda home, m: 1.0)

    def miniprot(work_hits, query_seq, *a, **k):
        if query_seq != home_seq:
            return []                                   # the relative's hits: no model, a chain
        return [([{"qstart": 1, "qend": 280, "gstart": 6_000, "gend": 6_840, "strand": "+",
                   "pident": 70.0, "seq": "M" * 280, "phase": 0}], "M" * 280)]
    monkeypatch.setattr(isr, "annotate_models_from_hit_list", miniprot)
    args = _window_args(tmp_path, RANKS)
    args.home_goi_ids = {"GOI_x"}
    anchors = [_anchor(g, 10_000 + 1_000 * i) for i, g in enumerate(RANKS)]
    block = {"chrom": "c1", "start": 10_000, "end": 14_000, "genes": list(RANKS), "genes_count": 4}
    db = {"GOI_x": {"id": "GOI_x", "seq": home_seq, "header": "GOI_x"},
          "GOI_x|H_fna_b0_l1_exon_ann": {"id": "GOI_x|H_fna_b0_l1_exon_ann", "seq": rel_seq,
                                         "header": "GOI_x|H_fna_b0_l1_exon_ann"}}
    _genes, gff = isr.process_region_block(0, block, anchors, {"c1": "ATG" * 10_000}, db, "G.fna", args,
                                           "u", 1)
    rows = [ln for ln in gff if ln.split("\t")[2] == "mRNA" and "SynVoyRole=goi" in ln]
    assert rows and not any("RefusalLifted" in ln for ln in rows)
    assert "fallback_lifted" not in args.genome_ctx.get("limits", {})


def test_the_lifted_chains_are_counted_where_they_reach_the_judge():
    lifted = _cand("GOI_x", 11_001, 11_060, "fallback_hit_span", ";RefusalLifted=long_qcov")
    decoy = _cand("GOI_x~rev", 11_001, 11_060, "fallback_hit_span", ";RefusalLifted=long_qcov")
    plain = _cand("GOI_x", 12_001, 12_060, "fallback_hit_span")
    assert isr.count_lifted_chains([lifted, decoy, plain, {"gff": []}]) == 1
    body = open(os.path.join(ROOT, "bin", "iterative_search_runner.py")).read().split(
        "def process_region_block(", 1)[1].split("\ndef ", 1)[0]
    assert body.count("'fallback_lifted'") == 1
    assert body.index("rebase_relative_candidates(") < body.index("'fallback_lifted'")


# ------------------------------------------------------------------ the relative's query
def _cand(parent, start, end, evidence, extra=""):
    attrs = (f"ID={parent}|G_fna_b0_l1_x;Name={parent};SynVoy_Parent={parent};SynVoyRole=goi;"
             f"EvidenceType={evidence};Identity=60.0;Exons=1;BlockFlankingSupport=4{extra}")
    return {"start": start, "end": end, "score": 60.0,
            "record": {"id": f"{parent}|G_fna_b0_l1_x", "seq": "M" * 40, "description": ""},
            "gff": [f"c1\tx\tmRNA\t{start}\t{end}\t60.0\t+\t.\t{attrs}"]}


def test_a_lifted_home_chain_does_not_drop_the_relatives_model(monkeypatch):
    """A model found with the relative's query is dropped where the home protein made its own model
    on the same DNA. A lifted chain is no such model: both reach the window's judge."""
    monkeypatch.setattr(isr, "sw_align_span",
                        lambda home, prot: {"identity": 55.0, "query_cov": 0.8, "q_start": 1, "q_end": 240})
    by_parent = {"GOI_x": {"seq": "M" * 300}, "GOI_x~rel": {"seq": "M" * 290, "model_id": "GOI_x|H_b0_l1"}}
    lifted = _cand("GOI_x", 11_001, 11_060, "fallback_hit_span", ";RefusalLifted=long_qcov,long_aln")
    rel_model = _cand("GOI_x~rel", 10_801, 11_900, "exon_annotation")
    kept = isr.rebase_relative_candidates([lifted, rel_model], ["GOI_x", "GOI_x~rel"], by_parent,
                                          genome_tag="G_fna")
    assert len(kept) == 2
    plain = _cand("GOI_x", 11_001, 11_060, "fallback_hit_span")
    kept = isr.rebase_relative_candidates([plain, rel_model], ["GOI_x", "GOI_x~rel"], by_parent,
                                          genome_tag="G_fna")
    assert len(kept) == 1                                    # as before: the home chain's DNA


def test_a_rebased_lifted_chain_keeps_its_flag(monkeypatch):
    monkeypatch.setattr(isr, "sw_align_span",
                        lambda home, prot: {"identity": 50.0, "query_cov": 0.1, "q_start": 1, "q_end": 30})
    by_parent = {"GOI_x": {"seq": "M" * 300}, "GOI_x~rel": {"seq": "M" * 290, "model_id": "GOI_x|H_b0_l1"}}
    rel_chain = _cand("GOI_x~rel", 11_001, 11_060, "fallback_hit_span", ";RefusalLifted=long_qcov")
    (kept,) = isr.rebase_relative_candidates([rel_chain], ["GOI_x~rel"], by_parent, genome_tag="G_fna")
    a = isr._parse_gff_attr_string(kept["gff"][0].split("\t")[8])
    assert a["RefusalLifted"] == "long_qcov" and a["SynVoy_Parent"] == "GOI_x"


# ------------------------------------------------------------------ seeds and the tree
def _meta(conf, lifted=""):
    return {"role": "goi", "parent": "GOI_x", "confidence": conf, "refusal_lifted": lifted,
            "real_competitor": "0"}


def _split(metas):
    genes = [{"id": f"GOI_x|G_fna_b0_l{i}_m", "seq": "M"} for i in range(len(metas))]
    meta = {g["id"]: m for g, m in zip(genes, metas)}
    seed, tree, suppressed = isr._classify_goi_for_seed_and_tree(genes, meta, isr.is_goi_query_id)
    ids = lambda gs: [int(g["id"].split("_l")[1].split("_")[0]) for g in gs]  # noqa: E731
    return ids(seed), ids(tree), suppressed


def test_a_lifted_call_neither_seeds_nor_feeds_the_tree():
    assert _split([_meta("MEDIUM", "long_qcov")]) == ([], [], 1)
    assert _split([_meta("MEDIUM")]) == ([0], [], 0)                     # regression


def test_a_lifted_call_takes_no_part_in_the_pool():
    """A lifted HIGH would make the genome's MEDIUM calls tree-only; it is left out first, so the
    MEDIUM call still seeds."""
    assert _split([_meta("HIGH", "short_bits"), _meta("MEDIUM")]) == ([1], [], 1)
    assert _split([_meta("HIGH"), _meta("MEDIUM", "long_qcov"), _meta("MEDIUM")]) == ([0], [2], 2)


def test_the_genome_end_reads_the_flag_into_the_seed_table():
    src = open(os.path.join(ROOT, "bin", "iterative_search_runner.py")).read()
    genome = src.split("def process_single_genome(", 1)[1].split("\ndef ", 1)[0]
    assert '"refusal_lifted": attrs.get("RefusalLifted", "")' in genome


# ------------------------------------------------------------------ the report
def _report_dir(tmp_path, rows):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir(exist_ok=True)
    lines = []
    for mid, chrom, conf, ident, qcov, status, extra in rows:
        lines.append(f"{chrom}\tx\tmRNA\t1001\t1900\t.\t+\t.\tID={mid};SynVoy_Parent=GOI_q;"
                     f"SynVoyRole=goi;Confidence={conf};GOIClass=probable_goi;EvidenceType=fallback_hit_span;"
                     f"Identity={ident};QueryCoverage={qcov};ModelStatus={status}{extra}")
    (tmp_path / "regions" / "locus_1__G.fna.gff").write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
    return gr.build_report(str(tmp_path))


def test_the_report_counts_counted_fragments_and_lifted_calls(tmp_path):
    report = _report_dir(tmp_path, [
        ("GOI_q|G_fna_b0_l1_fallback", "c1", "MEDIUM", 42, 0.2, "fragment", ";RefusalLifted=long_qcov"),
        ("GOI_q|G_fna_b0_l2_fallback", "c2", "MEDIUM", 42, 0.25, "fragment", ""),
        # demoted by the coverage rule (high identity over a fifth of the query): not counted
        ("GOI_q|G_fna_b0_l3_fallback", "c3", "MEDIUM", 90, 0.2, "fragment", ";RefusalLifted=long_qcov"),
        ("GOI_q|G_fna_b0_l4_fallback", "c4", "MEDIUM", 60, 0.9, "partial", ""),
        ("GOI_q|G_fna_b0_l5_fallback", "c5", "LOW", 30, 0.1, "fragment", ";RefusalLifted=short_bits")])
    records = {r["chrom"]: r for r in report["goi_dedup"]["records"]}
    assert records["c1"]["refusal_lifted"] == "long_qcov"
    assert "refusal_lifted" not in records["c2"] and "refusal_lifted" not in records["c4"]
    assert records["c3"]["goi_class"] == "identity_coverage_decoupled"
    assert report["summary"]["counted_fragment_calls"] == 2
    assert report["summary"]["refusal_lifted_calls"] == 1
