"""A leftover hit joined into a confident gene model (batch of 2026-10-15, G3B;
--goi_join_leftover_hits, off by default). A hit no model of its locus covers is offered to a
HIGH or MEDIUM model of its own query (home query or decoy, never the relative's) when it lies
outside the model on its strand, continues the model's query from its first or last CDS by G1's
rule and has no other model in the gap; the re-modelled gene replaces the model only when it
covers the hit and ranks strictly first. All synthetic."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402


def _cds(qs, qe, gs, ge):
    return {"qstart": qs, "qend": qe, "gstart": gs, "gend": ge}


def _hit(qs, qe, gs, ge, strand="+", bits=40.0):
    return dict(_cds(qs, qe, gs, ge), strand=strand, bits=bits)


def _model(exons, rank, label="MEDIUM", parent="GOI_q", strand="+", ident="m"):
    return {"model_exons": exons, "model_parent": parent, "model_strand": strand, "model_tag": "l1",
            "window_label": label, "rank_key": rank, "record": {"id": ident}}


def _leftover(hit, parent="GOI_q", ident="r"):
    return {"leftover_hit": hit, "leftover_parent": parent, "record": {"id": ident}}


MODEL = [_cds(41, 120, 5000, 5240)]


def test_the_side_of_a_hit_that_continues_the_query():
    assert isr.leftover_join_side(MODEL, _hit(1, 40, 2000, 2120), "+", 50000) == "before"
    assert isr.leftover_join_side([_cds(1, 80, 1000, 1240)], _hit(81, 120, 3000, 3120), "+", 50000) == "after"
    # minus strand: the genomically first CDS holds the query's end
    minus = [_cds(41, 120, 1000, 1240)]
    assert isr.leftover_join_side(minus, _hit(121, 160, 200, 320, "-"), "-", 50000) == "before"
    assert isr.leftover_join_side(minus, _hit(1, 40, 3000, 3120, "-"), "-", 50000) == "after"


def test_refused_restart_strand_gap_fit_and_inside():
    assert isr.leftover_join_side(MODEL, _hit(50, 90, 2000, 2120), "+", 50000) is None      # restart
    assert isr.leftover_join_side(MODEL, _hit(1, 40, 2000, 2120, "-"), "+", 50000) is None  # strand
    assert isr.leftover_join_side(MODEL, _hit(1, 40, 2000, 2120), "+", 1000) is None        # max gap
    far = [_cds(200, 300, 5000, 5300)]                                  # 179 residues skipped, 40 bp
    assert isr.leftover_join_side(far, _hit(1, 20, 4900, 4960), "+", 50000) is None
    two = [_cds(1, 40, 1000, 1120), _cds(81, 120, 5000, 5120)]
    assert isr.leftover_join_side(two, _hit(41, 80, 3000, 3120), "+", 50000) is None        # inside
    minus = [_cds(41, 120, 1000, 1240)]
    assert isr.leftover_join_side(minus, _hit(1, 40, 200, 320, "-"), "-", 50000) is None    # restart


def test_who_may_be_joined():
    home = {"GOI_q"}
    assert isr.leftover_join_eligible(_model(MODEL, (), "HIGH"), home)
    assert isr.leftover_join_eligible(_model(MODEL, (), "MEDIUM", parent="GOI_q~rev"), home)   # decoy
    assert not isr.leftover_join_eligible(_model(MODEL, (), "AMBIGUOUS"), home)
    assert not isr.leftover_join_eligible(_model(MODEL, (), "LOW"), home)
    assert not isr.leftover_join_eligible(_model(MODEL, (), "HIGH", parent="GOI_q~rel"), home)
    assert not isr.leftover_join_eligible(_model(MODEL, (), "HIGH", parent="GOI_copy_1|g_b0_l1"), home)


def _rebuilt(rank, exons=None):
    calls = []

    def rebuild(model, hit):
        calls.append((model["record"]["id"], hit["qstart"]))
        joined = _model(exons or [_cds(1, 40, 2000, 2120)] + model["model_exons"], rank, ident="m")
        joined["record"] = {"id": "m"}
        return joined
    return rebuild, calls


def _run(cands, rebuild, max_gap=50000):
    return isr.join_leftover_hits(cands, max_gap, lambda c: isr.leftover_join_eligible(c, {"GOI_q"}), rebuild)


def test_a_join_that_ranks_first_replaces_the_model_and_its_leftover():
    model, left = _model(MODEL, (-3, 1e-5)), _leftover(_hit(1, 40, 2000, 2120))
    rebuild, calls = _rebuilt((-4, 1e-9))
    out, n = _run([model, left], rebuild)
    assert n == 1 and calls == [("m", 1)]
    assert len(out) == 1 and out[0]["rank_key"] == (-4, 1e-9) and out[0] is not model


def test_no_join_unless_strictly_first_and_covering_the_hit():
    for rank in ((-3, 1e-5), (-3, 1e-3)):                               # a tie, or worse
        model, left = _model(MODEL, (-3, 1e-5)), _leftover(_hit(1, 40, 2000, 2120))
        out, n = _run([model, left], _rebuilt(rank)[0])
        assert n == 0 and out[0] is model and out[1] is left
    model, left = _model(MODEL, (-3, 1e-5)), _leftover(_hit(1, 40, 2000, 2120))
    out, n = _run([model, left], _rebuilt((-4, 0.0), exons=[_cds(41, 120, 5000, 5240)])[0])
    assert n == 0 and out[0] is model                                     # the hit is not on it
    out, n = _run([model, left], lambda m, h: None)
    assert n == 0 and out == [model, left]


def test_a_model_or_copy_in_between_and_other_parents_are_never_offered():
    left = _leftover(_hit(1, 40, 2000, 2120))
    rebuild, calls = _rebuilt((-4, 0.0))
    other = _model([_cds(1, 30, 3000, 3090)], (-2, 1.0), "LOW", ident="o")
    _run([_model(MODEL, (-3, 1e-5)), left, other], rebuild)
    copy = {"copy_intervals": [(3000, 3200)], "record": {"id": "c"}}
    _run([_model(MODEL, (-3, 1e-5)), left, copy], rebuild)
    _run([_model(MODEL, (-3, 1e-5), "LOW"), left], rebuild)                       # not confident
    _run([_model(MODEL, (-3, 1e-5)), _leftover(_hit(1, 40, 2000, 2120), parent="GOI_q~rel")], rebuild)
    assert calls == []


def test_the_strongest_leftover_first_and_every_covered_leftover_dropped():
    model = _model(MODEL, (-3, 1e-5))
    weak, strong = _leftover(_hit(1, 40, 2000, 2120, bits=30), ident="w"), \
        _leftover(_hit(2, 40, 2003, 2120, bits=60), ident="s")
    rebuild, calls = _rebuilt((-4, 0.0))
    out, n = _run([model, weak, strong], rebuild)
    assert calls == [("m", 2)] and n == 1 and [c["record"]["id"] for c in out] == ["m"]


def test_the_window_label_is_kept_by_the_judge():
    cand = {"start": 1, "end": 90, "record": {"id": "x"},
            "gff": ["c\tsrc\tmRNA\t1\t90\t.\t+\t.\tID=x;SynVoy_Parent=GOI_q"]}
    isr.judge_candidates([cand])
    assert cand["window_label"] == "LOW"
