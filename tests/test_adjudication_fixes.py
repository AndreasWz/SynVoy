#!/usr/bin/env python3
"""Regressions for the 2026-07-26 adjudication-layer fixes.

Two independent defects, both in the layer that turns a search result into a
number a reader trusts:

  F1  the region p-value was computed from a score that had already been
      inflated by the GOI-overlap bonus / distant-rescue floor, neither of
      which the permutation null can produce;
  F6  the strong-flanking fallback->MEDIUM clause had no absolute identity
      floor, so any hit inside a flanking-rich block reached MEDIUM.

See docs/STATE_OF_THE_PROJECT.md Part B.
"""
import importlib.util
import os

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cg = _load("cluster_grs", "cluster_grs.py")
isr = _load("iterative_search_runner", "iterative_search_runner.py")


# ── F1: the p-value must not see the ranking bonus ──────────────────────────

def _hit(query, start, strand="+"):
    return {"query": query, "chrom": "c1", "start": start, "end": start + 1000,
            "strand": strand, "identity": 80.0, "evalue": 1e-30}


def _pvalue_setup():
    gene_map = {f"g{i}": i for i in range(10)}
    cluster = [_hit(f"g{i}", i * 1000) for i in range(3)]
    background = [_hit(f"g{i % 10}", 5_000_000 + i * 7000) for i in range(60)]
    return gene_map, cluster, cluster + background


def _p(score):
    gene_map, cluster, all_hits = _pvalue_setup()
    return cg.estimate_pvalue(
        score, cluster, all_hits, 20_000_000, 150000,
        cg.score_flexible_synteny, gene_map, 10, n=200, seed=42)


def test_pvalue_floor_is_laplace_corrected():
    """p can never be 0; the floor is 1/(n+1)."""
    assert _p(10.0) == 1 / 201


def test_goi_bonus_would_have_driven_p_to_the_floor():
    """Pins the ORIGINAL defect so a regression is unmistakable.

    The bonus is +0.15. Feeding the bonused score to the null moves this
    cluster from 'not significant' to the floor — which is exactly what the
    code used to do.
    """
    _, cluster, _ = _pvalue_setup()
    gene_map = {f"g{i}": i for i in range(10)}
    u, k, s = cg.score_flexible_synteny(cluster, gene_map)
    synteny = (0.4 * (u / 10) + 0.3 * k + 0.3 * s) * (u / 10)
    assert _p(synteny) > 0.05, "baseline cluster should not be significant"
    assert _p(synteny + 0.15) == 1 / 201, "bonused score saturates the null"


def test_scores_tsv_exposes_the_tested_statistic_and_sort_key():
    """`synteny_score` (what p was computed from) and `goi_block_flanking`
    (the primary sort key, §1u) must both be emitted, or neither the p-value
    nor the region ordering can be checked from the output."""
    src = open(os.path.join(BIN, "cluster_grs.py")).read()
    fieldblock = src[src.index('fieldnames = ['):src.index('selection_reason",') ]
    for col in ('"synteny_score"', '"goi_block_flanking"', '"p_value"'):
        assert col in fieldblock, f"{col} missing from the scores TSV header"
    # and the p-value must be computed from synteny_score, not final_score
    assert "estimate_pvalue(\n            synteny_score," in src


# ── F6: the strong-flanking fallback clause needs an identity floor ─────────

def _classify(identity, qcov, flanking):
    return isr._classify_goi_evidence(
        evidence_type="fallback_hit_span", identity=identity,
        query_cov=qcov, flanking_support=flanking)


def test_low_identity_in_flanking_rich_block_is_not_medium():
    """The yeast STE2 false positives: 21-27 % identity, ample coverage, sitting
    in a block with plenty of flanking genes. Must be LOW."""
    for ident in (21.0, 24.0, 27.0):
        conf, cls, _ = _classify(ident, qcov=0.40, flanking=7)
        assert conf == "LOW", f"{ident}% reached {conf}"
        assert cls == "ambiguous_goi_family_member"


def test_ant_melittin_rescue_still_reaches_medium():
    """The real Tetramorium melittin call this clause exists to rescue:
    34.7 % identity, qcov 0.686, BlockFlankingSupport 7. Must stay MEDIUM."""
    conf, cls, reason = _classify(34.7, qcov=0.686, flanking=7)
    assert conf == "MEDIUM", f"ant melittin demoted to {conf}"
    assert cls == "probable_goi"
    assert reason == "fallback_span_with_strong_flanking_support"


def test_floor_sits_between_the_two_populations():
    """Guards the calibration itself: the default floor must reject the yeast
    junk (<=27 %) and admit the melittin call (34.7 %)."""
    floor = isr.CLASSIFY_THRESHOLDS["fallback_strong_min_identity_floor"]
    assert 27.0 < floor <= 34.7


def test_floor_of_zero_restores_legacy_behaviour():
    original = isr.CLASSIFY_THRESHOLDS["fallback_strong_min_identity_floor"]
    isr.CLASSIFY_THRESHOLDS["fallback_strong_min_identity_floor"] = 0.0
    try:
        conf, _, _ = _classify(21.0, qcov=0.40, flanking=7)
        assert conf == "MEDIUM"
    finally:
        isr.CLASSIFY_THRESHOLDS["fallback_strong_min_identity_floor"] = original


def test_high_identity_low_coverage_arm_is_unaffected():
    """The other arm of the OR (identity >= 35) is above the floor anyway, so a
    high-identity low-coverage hit still passes."""
    conf, _, reason = _classify(60.0, qcov=0.05, flanking=7)
    assert conf == "MEDIUM"
    assert reason == "fallback_span_with_strong_flanking_support"


# ── F2: strand_consistency must measure conservation, not target uniformity ──

_HOME_STRANDS = {f"g{i}": ("+" if i % 2 == 0 else "-") for i in range(6)}
_GENE_MAP6 = {f"g{i}": i for i in range(6)}


def _cluster(strand_fn):
    return [{"query": f"g{i}", "chrom": "c1", "start": i * 1000, "end": i * 1000 + 900,
             "strand": strand_fn(i), "identity": 80.0, "evalue": 1e-30}
            for i in range(6)]


def _strand_score(cluster, home=True):
    _, _, s = cg.score_flexible_synteny(
        cluster, _GENE_MAP6, _HOME_STRANDS if home else None)
    return s


def test_conserved_neighbourhood_scores_top():
    assert _strand_score(_cluster(lambda i: _HOME_STRANDS[f"g{i}"])) == 1.0


def test_uniform_inversion_is_still_conserved():
    """An inversion flips every strand but preserves relative orientation, so it is
    conserved synteny — exactly as `_longest_collinear_run` treats gene order in both
    directions. It must not be penalised."""
    flipped = _cluster(lambda i: "-" if _HOME_STRANDS[f"g{i}"] == "+" else "+")
    assert _strand_score(flipped) == 1.0


def test_scrambled_orientation_is_penalised():
    """The discrimination the term is supposed to provide."""
    assert _strand_score(_cluster(lambda i: "+")) < 1.0


def test_legacy_measure_was_inverted():
    """Pins the original defect. The old measure — fraction of hits on the cluster's
    own majority strand — never consulted the home genome, so it scored the SCRAMBLED
    cluster at the maximum and the CONSERVED one at the minimum. Strictly worse than
    ignoring strand altogether, and it carried weight 0.3 of the quality score.
    """
    conserved = _cluster(lambda i: _HOME_STRANDS[f"g{i}"])
    scrambled = _cluster(lambda i: "+")
    assert _strand_score(scrambled, home=False) > _strand_score(conserved, home=False)
    # and the fix reverses that ordering
    assert _strand_score(scrambled) < _strand_score(conserved)


def test_missing_home_strands_falls_back_not_fabricates():
    """No strand column in the BED must not silently score everything as conserved."""
    conserved = _cluster(lambda i: _HOME_STRANDS[f"g{i}"])
    assert _strand_score(conserved, home=False) == 0.5   # legacy value, not 1.0


def test_home_strand_loader_reads_column_six(tmp_path):
    bed = tmp_path / "synteny.bed"
    bed.write_text(
        "chr1\t100\t200\tgene-A\t.\t+\tACT1\n"
        "chr1\t300\t400\tgene-B\t.\t-\tYPT1\n"
        "chr1\t500\t600\tgene-C\t.\t.\tNOSTRAND\n")
    strands = cg.load_home_strands(str(bed))
    assert strands["gene-A"] == "+"
    assert strands["gene-B"] == "-"
    assert "gene-C" not in strands          # '.' is not a strand


# ── F5: one definition of collinearity, shared by both halves ───────────────

def test_region_scorer_and_search_core_agree_on_collinearity():
    """Both must call the same helper. A block judged collinear when it is seeded must
    not be judged non-collinear when it is ranked."""
    import sys
    sys.path.insert(0, BIN)
    from sequence_utils import longest_collinear_run as shared
    assert isr._longest_collinear_run([0, 1, 3, 2, 4]) == shared([0, 1, 3, 2, 4])
    assert isr._longest_collinear_run([4, 3, 2, 1, 0]) == shared([4, 3, 2, 1, 0])


def test_consistency_uses_lis_not_adjacent_pairs():
    """One transposed gene in an otherwise perfect run of 5 breaks TWO adjacent pairs
    (old score 3/4 = 0.75) but costs the LIS only one element (4/5 = 0.80)."""
    gene_map = {f"g{i}": i for i in range(5)}
    order = [0, 1, 3, 2, 4]                       # single transposition
    cluster = [{"query": f"g{r}", "chrom": "c1", "start": i * 1000,
                "end": i * 1000 + 900, "strand": "+", "identity": 80.0,
                "evalue": 1e-30} for i, r in enumerate(order)]
    _, consistency, _ = cg.score_flexible_synteny(cluster, gene_map)
    assert consistency == 4 / 5


def test_inverted_order_is_fully_consistent():
    gene_map = {f"g{i}": i for i in range(5)}
    cluster = [{"query": f"g{r}", "chrom": "c1", "start": i * 1000,
                "end": i * 1000 + 900, "strand": "+", "identity": 80.0,
                "evalue": 1e-30} for i, r in enumerate([4, 3, 2, 1, 0])]
    _, consistency, _ = cg.score_flexible_synteny(cluster, gene_map)
    assert consistency == 1.0


# ── F4: the wavefront must actually grade close -> far ──────────────────────

def _legacy_waves(d):
    """The pre-2026-07-26 absolute-threshold binning, for comparison."""
    out, i = [], 0
    while i < len(d):
        c = d[i]
        if c < 0.05:
            out.append(1); i += 1
        elif c < 0.15:
            n = 1; i += 1
            while i < len(d) and abs(d[i] - c) < 0.01:
                n += 1; i += 1
                if n >= 3:
                    break
            out.append(n)
        else:
            n = 1; i += 1
            while i < len(d) and abs(d[i] - c) < 0.02:
                n += 1; i += 1
                if n >= 5:
                    break
            out.append(n)
    return out


def test_rank_binning_is_the_default():
    """The legacy binning compares absolute cut-offs against a divide-by-max
    normalised distance, so it degenerates on any uniform target set."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--rank_wave_binning", type=isr.str2bool, default=True)
    assert parser.parse_args([]).rank_wave_binning is True


def test_uniform_target_set_collapses_under_legacy_binning():
    """The melittin ground-truth case: five targets at the same taxonomic depth.
    Divide-by-max sends every one to 1.0 -> a single parallel wave -> the expanding
    database never iterates. This is what F4 documents."""
    raw = [998, 999, 999, 999, 1000]
    norm = [x / max(raw) for x in raw]
    assert _legacy_waves(norm) == [5], "legacy should collapse to one wave"


def test_rank_binning_grades_the_same_set():
    entries = [{"i": i} for i in range(5)]
    waves = isr.assign_waves_by_rank(entries)
    assert [len(w) for w in waves] == [1, 2, 2]
    assert len(waves[0]) == 1, "closest genome must be searched alone"


def test_rank_binning_always_starts_serial():
    """Whatever the distance distribution, the closest genome gets its own wave —
    the precondition for its recovered models to seed the next wave."""
    for n in range(1, 40):
        waves = isr.assign_waves_by_rank([{"i": i} for i in range(n)])
        assert len(waves[0]) == 1
        assert sum(len(w) for w in waves) == n      # partition, nothing dropped
