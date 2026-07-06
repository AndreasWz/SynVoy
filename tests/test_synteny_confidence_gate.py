"""Tests for §1 synteny-ORDER confidence gating in _classify_goi_evidence.

Synteny already entered confidence as flanking_support (a COUNT). §1 adds the ORDER
dimension: a GOI in a SCRAMBLED paralog neighbourhood (enough flanking, but not in
conserved home order — the biglycan-as-decorin signature) is demoted from HIGH to
MEDIUM, while small ordered blocks (melittin/LY6) and the no-collinearity path are
untouched.
"""
import importlib.util
import os

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")
spec = importlib.util.spec_from_file_location(
    "iterative_search_runner", os.path.join(BIN, "iterative_search_runner.py"))
isr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(isr)


def classify(**kw):
    base = dict(evidence_type="exon_annotation", identity=80.0, exon_count=3,
                query_cov=0.8, flanking_support=6)
    base.update(kw)
    return isr._classify_goi_evidence(**base)


def test_collinear_none_is_legacy_high():
    # No collinearity info (home_rank unavailable) → order gate waived → HIGH on count.
    conf, cls, _ = classify(collinear_support=None)
    assert conf == "HIGH" and cls == "confident_goi"


def test_ordered_neighbourhood_stays_high():
    # 6 flanking, all collinear (run 6) → true ortholog locus → HIGH.
    conf, cls, _ = classify(flanking_support=6, collinear_support=6)
    assert conf == "HIGH" and cls == "confident_goi"


def test_scrambled_neighbourhood_demoted_to_medium():
    # 6 flanking but longest collinear run is only 2 (< high_min_collinear=3) → paralog
    # neighbourhood → demote to MEDIUM, don't crown HIGH on identity+count alone.
    conf, cls, reason = classify(flanking_support=6, collinear_support=2)
    assert conf == "MEDIUM" and cls == "probable_goi"
    assert reason == "high_identity_but_flanking_not_collinear"


def test_small_ordered_block_never_gated():
    # 2 flanking (< high_min_collinear) → can't judge order → waived → HIGH preserved.
    # This is the melittin/LY6 single-copy safety case.
    conf, cls, _ = classify(flanking_support=2, collinear_support=1)
    assert conf == "HIGH" and cls == "confident_goi"


def test_gate_only_touches_high_not_medium_path():
    # Below HIGH identity: the MEDIUM path is unchanged regardless of collinearity.
    conf, cls, _ = classify(identity=40.0, flanking_support=6, collinear_support=1)
    assert conf == "MEDIUM"  # medium_min_identity path, not the collinear demotion


def test_disable_via_huge_threshold_is_count_only(monkeypatch):
    # classify_high_min_collinear=0 maps to a huge threshold → order gate never fires.
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "high_min_collinear", 10 ** 9)
    conf, cls, _ = classify(flanking_support=6, collinear_support=1)
    assert conf == "HIGH"


def test_longest_collinear_run_helper():
    assert isr._longest_collinear_run([1, 2, 3, 4])[0] == 4        # ordered
    assert isr._longest_collinear_run([4, 3, 2, 1])[0] == 4        # reverse-ordered (still collinear)
    assert isr._longest_collinear_run([1, 9, 2, 8, 3])[0] == 3     # scrambled → run of 3 (1,2,3)
    assert isr._longest_collinear_run([5, 1, 4, 2, 3])[0] <= 3     # scrambled
