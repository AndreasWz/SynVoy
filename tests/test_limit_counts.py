"""Stage 2 of the profile handoff: every limit that can drop a result is counted, and the
counts change nothing (iterative_search_runner.LIMIT_KEYS, the limit_* columns of
hits/<genome>.target_policy.tsv, the report's `limits` block; annotate_goi_exons.LIMIT_COUNTS;
split_loci's --max_loci WARNING; rescue_strong_synteny's timeout line).

Each counter: one case that trips it and one that does not.
"""
import csv
import logging
import os
import subprocess
import sys
import types

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import annotate_goi_exons as age  # noqa: E402
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402


def _ctx_args(**over):
    a = types.SimpleNamespace(genome_ctx={}, **over)
    return a


def _limits(args):
    return args.genome_ctx.get("limits", {})


@pytest.fixture(autouse=True)
def _clean_counts():
    age.LIMIT_COUNTS.clear()
    yield
    age.LIMIT_COUNTS.clear()


# ============================================================ helpers
def test_count_limit_needs_a_genome_and_a_count():
    isr.count_limit(types.SimpleNamespace(), "sw_timeout")          # no genome: no-op
    args = _ctx_args()
    isr.count_limit(args, "sw_timeout", 0)
    assert _limits(args) == {}
    isr.count_limit(args, "sw_timeout")
    isr.count_limit(args, "sw_timeout", 2)
    assert _limits(args) == {"sw_timeout": 3}


def test_limit_query_kind():
    assert isr.limit_query_kind("GOI_q|frag_1-40") == "goi"
    assert isr.limit_query_kind("GOI_q~rel") == "rel"
    assert isr.limit_query_kind("GOI_q~rev|frag_1-40") == "decoy"
    assert isr.limit_query_kind("GOI_q~shuf2~rel") == "decoy"
    assert isr.limit_query_kind("GOI_copy_1") == "goi"


def test_every_key_has_a_column():
    assert len(set(isr.LIMIT_KEYS)) == len(isr.LIMIT_KEYS)
    assert isr.LIMIT_COLUMNS == [f"limit_{k}" for k in isr.LIMIT_KEYS]


# ============================================================ miniprot calls
class _Proc:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode, self.stdout, self.stderr = returncode, stdout, stderr


@pytest.fixture
def miniprot_on(monkeypatch):
    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    return monkeypatch


@pytest.mark.parametrize("sensitive,kind", [(True, "goi"), (False, "flanking")])
def test_miniprot_timeout_is_counted_per_kind(miniprot_on, sensitive, kind):
    def boom(cmd, **k):
        raise subprocess.TimeoutExpired(cmd, 120)
    miniprot_on.setattr(age.subprocess, "run", boom)
    assert age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=sensitive) == []
    assert age.LIMIT_COUNTS == {f"miniprot_timeout_{kind}": 1}


def test_miniprot_failure_and_no_alignment_are_counted_apart(miniprot_on):
    miniprot_on.setattr(age.subprocess, "run", lambda cmd, **k: _Proc(returncode=1, stderr="bad"))
    age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=True)
    miniprot_on.setattr(age.subprocess, "run", lambda cmd, **k: _Proc(stdout="  \n"))
    age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=True)
    miniprot_on.setattr(age.subprocess, "run", lambda cmd, **k: _Proc(stdout="##gff-version 3\n"))
    miniprot_on.setattr(age, "_parse_miniprot_gff", lambda text: [])
    age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=False)
    assert age.LIMIT_COUNTS == {"miniprot_failed_goi": 1, "miniprot_no_alignment_goi": 1,
                                "miniprot_no_alignment_flanking": 1}


def test_a_decoy_querys_miniprot_calls_are_not_counted(miniprot_on):
    def boom(cmd, **k):
        raise subprocess.TimeoutExpired(cmd, 120)
    miniprot_on.setattr(age.subprocess, "run", boom)
    with isr.miniprot_counting(False):
        age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=True)
    assert age.LIMIT_COUNTS == {} and age.LIMIT_COUNTING is True
    with isr.miniprot_counting(True):
        age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=True)
    assert age.LIMIT_COUNTS == {"miniprot_timeout_goi": 1}


def test_a_miniprot_model_counts_nothing(miniprot_on):
    miniprot_on.setattr(age.subprocess, "run", lambda cmd, **k: _Proc(stdout="x\n"))
    miniprot_on.setattr(age, "_parse_miniprot_gff", lambda text: [{"score": 10}])
    assert age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=True) == [{"score": 10}]
    assert age.LIMIT_COUNTS == {}


# ============================================================ models per GOI locus
def _model(start, end, score):
    return {"score": score, "cds_list": [{"gstart": start, "gend": end}]}


def _fake_models(monkeypatch, rounds):
    """_run_miniprot returns rounds[i] in round i (a callable may count a timeout)."""
    calls = {"n": 0}

    def run(query, target, chrom, max_intron, sensitive):
        i = calls["n"]
        calls["n"] += 1
        r = rounds[i] if i < len(rounds) else []
        return r() if callable(r) else list(r)
    monkeypatch.setattr(age, "_run_miniprot", run)
    monkeypatch.setattr(age, "_build_model",
                        lambda m, *a, **k: ([{"gstart": m["cds_list"][0]["gstart"]}], "M"))
    return calls


def test_model_cap_stop_is_counted(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100), _model(1000, 1300, 90)]])
    kept = age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=1)
    assert len(kept) == 1
    assert age.LIMIT_COUNTS["goi_model_cap_stopped"] == 1
    assert age.LIMIT_COUNTS["goi_model_loci"] == 1
    assert age.LIMIT_COUNTS["goi_models_found"] == 2
    assert age.LIMIT_COUNTS["goi_models_kept"] == 1


def test_cap_is_not_counted_when_the_rest_is_the_same_gene(monkeypatch):
    # The second alignment lies on the kept model: it would be dropped anyway.
    _fake_models(monkeypatch, [[_model(0, 300, 100), _model(10, 290, 90)]])
    age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=1)
    assert "goi_model_cap_stopped" not in age.LIMIT_COUNTS


def test_cap_is_not_counted_below_the_cap(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100), _model(1000, 1300, 90)]])
    kept = age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=3)
    assert len(kept) == 2
    assert "goi_model_cap_stopped" not in age.LIMIT_COUNTS


def test_cap_counted_when_a_hit_is_still_unexplained(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100)]])
    age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=1,
                             evidence=[(1500, 1600)])
    assert age.LIMIT_COUNTS["goi_model_cap_stopped"] == 1


def test_masked_round_below_the_score_bound_is_counted(monkeypatch):
    # Round 2 (masked) returns a model at 0.2 x the best: below the 0.3 bound.
    _fake_models(monkeypatch, [[_model(0, 300, 100)], [_model(1500, 1600, 20)]])
    age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=5,
                             evidence=[(1500, 1600)])
    assert age.LIMIT_COUNTS["goi_models_below_outs"] == 1


def test_masked_round_above_the_bound_is_kept_and_not_counted(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100)], [_model(1500, 1600, 50)]])
    kept = age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=5,
                                    evidence=[(1500, 1600)])
    assert len(kept) == 2
    assert "goi_models_below_outs" not in age.LIMIT_COUNTS
    assert "goi_models_near_outs" not in age.LIMIT_COUNTS or age.LIMIT_COUNTS["goi_models_near_outs"] == 0


def test_kept_model_near_the_bound_is_counted(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100), _model(1000, 1300, 35)]])
    age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=5)
    assert age.LIMIT_COUNTS["goi_models_near_outs"] == 1


def test_timeout_that_ends_the_masked_rounds_is_counted(monkeypatch):
    def timeout():
        age.LIMIT_COUNTS["miniprot_timeout_goi"] = age.LIMIT_COUNTS.get("miniprot_timeout_goi", 0) + 1
        return []
    calls = _fake_models(monkeypatch, [[_model(0, 300, 100)], timeout, [_model(1500, 1600, 90)]])
    kept = age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=5,
                                    evidence=[(1500, 1600)])
    # The timeout ends the copy search: the third round never runs.
    assert len(kept) == 1 and calls["n"] == 2
    assert age.LIMIT_COUNTS["goi_model_rounds_ended_by_timeout"] == 1


def test_masked_round_without_alignment_is_no_timeout(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100)], []])
    age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=True, max_models=5,
                             evidence=[(1500, 1600)])
    assert "goi_model_rounds_ended_by_timeout" not in age.LIMIT_COUNTS


def test_flanking_models_are_not_counted(monkeypatch):
    _fake_models(monkeypatch, [[_model(0, 300, 100), _model(1000, 1300, 99.5)]])
    age.miniprot_gene_models("M" * 100, "A" * 2000, "c1", sensitive=False, max_models=1)
    assert age.LIMIT_COUNTS == {}


# ============================================================ window hits
def _m8(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def _row(q, alnlen, pident=40.0, evalue=1e-5):
    # query target pident alnlen mismatch gapopen qstart qend tstart tend evalue bits
    return [q, "region_seq", pident, alnlen, 0, 0, 1, alnlen, 100, 100 + 3 * alnlen, evalue, 40.0]


def test_parse_hits_counts_the_length_floor(tmp_path):
    p = tmp_path / "h.m8"
    _m8(p, [_row("GOI_q", 40), _row("GOI_q", 12), _row("GOI_q", 12, evalue=5.0),
            _row("GOI_q~rev", 12)])
    drops = {}
    hits = isr.parse_hits(str(p), 15.0, 15, 0.1, query_lengths={"GOI_q": 200, "GOI_q~rev": 200},
                          drops=drops)
    assert len(hits) == 1
    # Only the 12-aa real row: the one failing E is not a length drop, the decoy is left out.
    assert drops == {"len_floor": 1, "len_floor_lt15_long_query": 1}


def test_parse_hits_without_drops_or_short_query(tmp_path):
    p = tmp_path / "h.m8"
    _m8(p, [_row("GOI_q", 40), _row("GOI_s", 13)])
    drops = {}
    isr.parse_hits(str(p), 15.0, 15, 0.1, query_lengths={"GOI_q": 200, "GOI_s": 80}, drops=drops)
    # GOI_s is 80 aa: its 13-aa row fails the floor (15) but is no long-query drop.
    assert drops == {"len_floor": 1}
    assert isr.parse_hits(str(p), 15.0, 15, 0.1, query_lengths={"GOI_q": 200}) == \
        isr.parse_hits(str(p), 15.0, 15, 0.1, query_lengths={"GOI_q": 200}, drops={})


def test_queries_at_max_seqs(tmp_path):
    p = tmp_path / "h.m8"
    _m8(p, [_row("GOI_a", 40)] * 300 + [_row("GOI_b", 40)] * 299 + [_row("GOI_a~rev", 40)] * 300)
    assert isr.queries_at_max_seqs(str(p)) == 1
    assert isr.queries_at_max_seqs(str(tmp_path / "missing.m8")) == 0


def _aug_args(tmp_path, **over):
    a = dict(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_evalue_mult=1000.0,
             aug_relaxed_evalue_cap=10.0, aug_relaxed_parse_evalue_mult=10.0,
             aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6, min_identity=10.0,
             aug_relaxed_length_min=15, aug_relaxed_length_div=2, min_length=10,
             enable_smith_waterman=True, deterministic_goi_search=True,
             sw_min_score=20, sw_min_identity=15.0, sw_method="parasail",
             sw_max_alignments_per_piece=1, sw_timeout_seconds=300, genome_ctx={})
    a.update(over)
    return types.SimpleNamespace(**a)


def _run_window(tmp_path, monkeypatch, args, sw="ok"):
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")
    query = [{"id": "GOI_q", "seq": "M" + "A" * 119}]

    def fake_mmseqs(**k):
        _m8(k["aug_hits_file"], [_row("GOI_q", 40), _row("GOI_q", 12)] + [_row("GOI_x", 40)] * 300)
        return True, ""

    def fake_run(cmd, *a, **k):
        name = os.path.basename(str(cmd[0]))
        if name == "tblastn" and "-out" in cmd:
            _m8(cmd[cmd.index("-out") + 1], [_row("GOI_q", 30), _row("GOI_q", 11)])
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    def fake_guarded(cmd, guard, what, param=None, **k):      # the Smith-Waterman child
        assert os.path.basename(str(cmd[0])) == "python3" and param == "sw_timeout_seconds"
        if sw == "timeout":
            raise isr.StepFailed(f"{what} did not finish within {guard} s (--{param})")
        _m8(cmd[cmd.index("--output") + 1], [_row("GOI_q", 9, evalue=0.01), _row("GOI_q", 30, evalue=0.01)])
        return types.SimpleNamespace(
            returncode=0, stdout="",
            stderr="Warning: Smith-Waterman skipped for GOI_long: 12000 aa is longer than 10000 aa\n")

    monkeypatch.setattr(isr.subprocess, "run", fake_run)
    monkeypatch.setattr(isr, "run_guarded", fake_guarded)
    monkeypatch.setattr(isr.shutil, "which", lambda name: None)
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", fake_mmseqs)
    return isr.run_augmented_search(str(region), query, "g1", args, "t1", 1)


def test_window_search_counts_each_engines_length_floor_and_sw(tmp_path, monkeypatch):
    args = _aug_args(tmp_path)
    _run_window(tmp_path, monkeypatch, args)
    lim = _limits(args)
    assert lim["window_len_floor_mmseqs"] >= 1
    assert lim["window_len_floor_tblastn"] == 1
    assert lim["window_len_floor_sw"] == 1           # 9 aa < SW's floor of 10
    assert lim["window_len_floor_lt15_long_query"] >= 2
    assert lim["mmseqs_maxseqs_window"] == 1         # GOI_x: 300 rows
    assert lim["sw_long_query_skipped"] == 1
    assert "sw_timeout" not in lim


def test_a_sw_guard_failure_stops_the_window_search(tmp_path, monkeypatch):
    """Fix-phase item 18: a Smith-Waterman search that does not finish within the guard stops
    the genome's search (StepFailed, through run_augmented_search's catch-all); until
    2026-10-28 it was counted as sw_timeout and the window kept its other hits."""
    args = _aug_args(tmp_path)
    with pytest.raises(isr.StepFailed, match="Smith-Waterman"):
        _run_window(tmp_path, monkeypatch, args, sw="timeout")
    assert "sw_timeout" not in _limits(args)


def test_window_search_counts_nothing_without_a_genome(tmp_path, monkeypatch):
    args = _aug_args(tmp_path)
    del args.genome_ctx
    hits = _run_window(tmp_path, monkeypatch, args)
    assert hits and not hasattr(args, "genome_ctx")


# ============================================================ GOI candidates of a window
def test_core_filter():
    args = _ctx_args()
    assert isr.outside_goi_core(args, "GOI_q", 0, 100, 500, 1000)
    assert isr.outside_goi_core(args, "GOI_q~rev", 0, 100, 500, 1000)      # decoy: not counted
    assert not isr.outside_goi_core(args, "GOI_q", 400, 600, 500, 1000)
    assert _limits(args) == {"goi_core_filter": 1}


def test_regions_cap_per_query_kind():
    args = _ctx_args(goi_candidate_regions_per_query=2)
    loci = [[{"bits": b}] for b in (10, 30, 20, 5)]
    kept = isr.cap_goi_parent_loci(args, "GOI_q", loci, lambda l: l[0]["bits"])
    assert [l[0]["bits"] for l in kept] == [30, 20]
    isr.cap_goi_parent_loci(args, "GOI_q~rel", loci[:3], lambda l: l[0]["bits"])
    isr.cap_goi_parent_loci(args, "GOI_q~shuf1", loci, lambda l: l[0]["bits"])
    assert _limits(args) == {"regions_cap": 2, "regions_cap_rel": 1, "regions_cap_decoy": 2}


def test_regions_cap_not_counted_within_the_cap():
    args = _ctx_args(goi_candidate_regions_per_query=4)
    isr.cap_goi_parent_loci(args, "GOI_q", [[{"bits": 1}], [{"bits": 2}]], lambda l: l[0]["bits"])
    assert _limits(args) == {}


def test_flanking_loci_cap():
    args = _ctx_args()
    loci = [[{"bits": b}] for b in (5, 9, 7)]
    assert [l[0]["bits"] for l in isr.cap_flanking_loci(args, loci)] == [9, 7]
    assert _limits(args) == {"flanking_loci_cap": 1}
    args2 = _ctx_args()
    isr.cap_flanking_loci(args2, loci[:2])
    assert _limits(args2) == {}


@pytest.mark.parametrize("case,expected", [
    (dict(ordered_hit_count=1, qcov=0.5, single_aln=10, aln_total=10, best_bits=20.0, query_len=60), ["short_aln", "short_bits"]),
    (dict(ordered_hit_count=1, qcov=0.5, single_aln=40, aln_total=40, best_bits=20.0, query_len=60), ["short_bits"]),
    (dict(ordered_hit_count=1, qcov=0.5, single_aln=40, aln_total=40, best_bits=60.0, query_len=60), []),
    (dict(ordered_hit_count=1, qcov=0.1, single_aln=40, aln_total=40, best_bits=60.0, query_len=300), ["long_qcov"]),
    (dict(ordered_hit_count=1, qcov=0.5, single_aln=20, aln_total=20, best_bits=60.0, query_len=300), ["long_aln"]),
    (dict(ordered_hit_count=3, qcov=0.1, single_aln=10, aln_total=30, best_bits=50.0, query_len=300), ["long_qcov", "long_aln", "long_bits"]),
    (dict(ordered_hit_count=3, qcov=0.1, single_aln=10, aln_total=30, best_bits=70.0, query_len=300), []),
])
def test_fallback_refusal_reasons_agree_with_the_gate(case, expected):
    gate = dict(short_query_len=100, short_min_aln_aa=30, short_min_bits=35.0)
    assert isr.fallback_refusal(**case, **gate) == expected
    assert isr.is_valid_fallback(**case, **gate) == (not expected)


def test_fallback_refusal_is_counted_per_test():
    args = _ctx_args()
    isr.count_fallback_refusal(args, "GOI_q", ["long_qcov", "long_aln", "long_bits"])
    isr.count_fallback_refusal(args, "GOI_q~rev", ["short_aln"])            # decoy
    isr.count_fallback_refusal(args, "GOI_q", [])                          # emitted
    assert _limits(args) == {"fallback_refused": 1, "fallback_long_qcov": 1,
                             "fallback_long_aln": 1, "fallback_long_bits": 1}


def test_fallback_span_cap():
    args = _ctx_args()
    near = [{"gstart": 0, "gend": 100}, {"gstart": 90_000, "gend": 90_100}]
    far = [{"gstart": 0, "gend": 100}, {"gstart": 150_000, "gend": 150_100}]
    assert not isr.fallback_span_too_wide(args, "GOI_q", near, 70)      # < 100 kb for 70 aa
    assert isr.fallback_span_too_wide(args, "GOI_q", far, 70)
    assert isr.fallback_span_too_wide(args, "GOI_q", far, 300)         # 150 x 300 = 45 kb
    assert _limits(args) == {"fallback_span_cap": 2}


# ============================================================ blocks and the dispersed rescue
def test_capped_blocks_are_logged_and_counted(caplog):
    args = _ctx_args()
    blocks = [{"chrom": "c1", "start": 10, "end": 20, "genes": ["g1", "g2"], "genes_count": 2},
              {"chrom": "c2", "start": 30, "end": 40, "genes": ["gene-GOI"], "genes_count": 1}]
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.log_capped_blocks(args, "G.fna", blocks, 80, lambda b: "gene-GOI" in b["genes"],
                              {"gene-GOI"}, False)
    assert _limits(args) == {"blocks_capped": 2, "blocks_capped_goi_seeded": 1,
                             "blocks_capped_without_home_order": 2}
    lines = [r.getMessage() for r in caplog.records if "block cap dropped" in r.getMessage()]
    assert len(lines) == 2 and "block 81 c2:30-40" in lines[1] and "GOI proxy yes" in lines[1]


def test_no_capped_block_counts_nothing():
    args = _ctx_args()
    isr.log_capped_blocks(args, "G.fna", [], 80, lambda b: False, set(), True)
    assert _limits(args) == {}


def _goi(chrom, start, pident=60.0, evalue=1e-20, alnlen=80):
    return {"query": "GOI_q", "chrom": chrom, "start": start, "end": start + 300,
            "pident": pident, "evalue": evalue, "alnlen": alnlen}


def _anchor(gene, chrom, start):
    return {"query": gene, "chrom": chrom, "start": start, "end": start + 1000,
            "pident": 90.0, "evalue": 1e-50, "alnlen": 300, "bits": 500.0}


def test_dispersed_rescue_counts_each_test_in_order():
    anchors = [_anchor(f"gene-a{i}", "c1", 1_000_000 + i * 100_000) for i in range(4)]
    hits = [_goi("c1", 1_150_000, pident=30.0),                 # identity
            _goi("c1", 1_150_000, pident=30.0, evalue=1.0),     # identity first
            _goi("c1", 1_150_000, evalue=1e-3),                 # E
            _goi("c1", 1_150_000, alnlen=40),                   # length (the 50-aa floor)
            _goi("c1", 1_005_000),                              # near an anchor
            _goi("c1", 1_600_000)]                              # seed (inside the envelope)
    drops, near = {}, []
    seeds = isr.select_dispersed_goi_seeds(hits, anchors, region_padding=10_000,
                                           envelope_margin=500_000, drops=drops, near_out=near)
    assert len(seeds) == 1
    assert drops == {"dispersed_identity": 2, "dispersed_evalue": 1, "dispersed_alnlen": 1,
                     "dispersed_near_anchor": 1}
    assert [h["start"] for h in near] == [1_005_000]


def test_dispersed_rescue_without_counters_is_unchanged():
    anchors = [_anchor(f"gene-a{i}", "c1", 1_000_000 + i * 100_000) for i in range(4)]
    hits = [_goi("c1", 1_600_000), _goi("c1", 1_005_000)]
    a = isr.select_dispersed_goi_seeds(hits, anchors, region_padding=10_000, envelope_margin=500_000)
    b = isr.select_dispersed_goi_seeds(hits, anchors, region_padding=10_000, envelope_margin=500_000,
                                       drops={}, near_out=[])
    assert a == b and len(a) == 1


def test_dead_zone_and_genome_line(caplog):
    args = _ctx_args()
    args.genome_ctx.update({
        "goi_cores": [("c1", 1000, 5000)],
        "near_anchor_hits": [{"chrom": "c1", "start": 2000, "end": 2300},      # in a core
                             {"chrom": "c1", "start": 9000, "end": 9300},      # dead zone
                             {"chrom": "c2", "start": 2000, "end": 2300}]})    # dead zone
    age.LIMIT_COUNTS["miniprot_timeout_goi"] = 2
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.log_genome_limits(args, "G.fna")
    assert _limits(args) == {"miniprot_timeout_goi": 2, "dead_zone": 2}
    assert age.LIMIT_COUNTS == {}
    msg = [r.getMessage() for r in caplog.records if "Limit counts" in r.getMessage()]
    assert msg == ["[G.fna] Limit counts (what each limit dropped; docs/OUTPUT.md): "
                   "miniprot_timeout_goi 2, dead_zone 2"]


def test_genome_line_says_none(caplog):
    args = _ctx_args()
    args.genome_ctx["goi_cores"] = [("c1", 1000, 5000)]
    args.genome_ctx["near_anchor_hits"] = [{"chrom": "c1", "start": 2000, "end": 2300}]
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.log_genome_limits(args, "G.fna")
    assert _limits(args) == {}
    assert any(r.getMessage().endswith("none") for r in caplog.records)


# ============================================================ target_policy.tsv and the report
def _policy(tmp_path, ctx=None, name="G.fna"):
    args = types.SimpleNamespace(output_dir=str(tmp_path), decoy_mode=False)
    isr._write_target_policy(args, "locus_1_", name, 60.0, 10, ctx)
    path = tmp_path / "hits" / f"locus_1_{name}.target_policy.tsv"
    with open(path, newline="") as fh:
        return path, list(csv.DictReader(fh, delimiter="\t"))


def test_target_policy_carries_every_limit_column(tmp_path):
    _path, rows = _policy(tmp_path, {"limits": {"sw_timeout": 3}})
    cols = [k for k in rows[0] if k.startswith("limit_")]
    assert cols == isr.LIMIT_COLUMNS
    assert rows[0]["limit_sw_timeout"] == "3" and rows[0]["limit_dead_zone"] == "0"


def test_target_policy_before_the_genome_end_has_zeros(tmp_path):
    _path, rows = _policy(tmp_path)
    assert all(rows[0][c] == "0" for c in isr.LIMIT_COLUMNS)


def test_report_sums_the_limits_and_keeps_them_out_of_target_policy(tmp_path):
    p1, _ = _policy(tmp_path, {"limits": {"sw_timeout": 2, "regions_cap": 1}}, name="A.fna")
    p2, _ = _policy(tmp_path, {"limits": {"sw_timeout": 1}}, name="B.fna")
    out = gr.collect_limits([str(p1), str(p2)])
    assert out["genome_searches"] == 2
    assert out["counts"]["sw_timeout"] == {"total": 3, "genomes": ["A", "B"]}
    assert out["counts"]["regions_cap"] == {"total": 1, "genomes": ["A"]}
    assert out["counts"]["dead_zone"] == {"total": 0, "genomes": []}
    assert list(out["counts"]) == isr.LIMIT_KEYS
    policy = gr.collect_target_policy([str(p1), str(p2)])
    assert not any(k.startswith("limit_") for row in policy.values() for k in row)


def test_report_without_limit_columns_has_no_block(tmp_path):
    p = tmp_path / "old.target_policy.tsv"
    p.write_text("genome\tmedian_flank_identity\nG.fna\t50\n")
    assert gr.collect_limits([str(p)]) is None


# ============================================================ split_loci and the rescue
def test_max_loci_warning_whenever_the_cap_drops_a_locus(tmp_path, capsys):
    import split_loci
    bed = tmp_path / "hits.bed"
    # Four loci 1 Mb apart, bit scores 100, 90, 80, 70: all pass the ratio test (0.5).
    rows = [f"chr1\t{i * 1_000_000}\t{i * 1_000_000 + 300}\tGOI_q\t1e-30\t+\t{b}"
            for i, b in enumerate((100, 90, 80, 70))]
    bed.write_text("\n".join(rows) + "\n")
    sys.argv = ["split_loci.py", "--bed", str(bed), "--output_prefix", str(tmp_path / "locus"),
                "--max_loci", "2", "--locus_gap", "50000"]
    split_loci.main()
    err = capsys.readouterr().err
    assert "WARNING: --max_loci 2 dropped 2 home locus/loci" in err
    assert "candidate loci before filtering" not in err        # below the family threshold

    sys.argv = ["split_loci.py", "--bed", str(bed), "--output_prefix", str(tmp_path / "loc2"),
                "--max_loci", "5", "--locus_gap", "50000"]
    split_loci.main()
    assert "--max_loci" not in capsys.readouterr().err


def test_rescue_logs_its_miniprot_timeouts(tmp_path, monkeypatch, capsys):
    import rescue_strong_synteny as rss
    monkeypatch.setattr(rss, "TIMEOUTS", {"windows": 0, "runs": 0})

    def boom(cmd, **k):
        raise subprocess.TimeoutExpired(cmd, 5)
    monkeypatch.setattr(rss.subprocess, "run", boom)
    assert rss._run_miniprot_relaxed("w.fa", "q.fa", 0.05, 5) == ""
    monkeypatch.setattr(rss.subprocess, "run", lambda cmd, **k: _Proc(stdout="##gff\n"))
    assert rss._run_miniprot_relaxed("w.fa", "q.fa", 0.05, 5) == "##gff\n"
    assert rss.TIMEOUTS == {"windows": 1, "runs": 2}
