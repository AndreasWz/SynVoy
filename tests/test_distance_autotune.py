#!/usr/bin/env python3
"""Unit tests for §1g per-target phylogenetic-distance auto-tune.

`iterative_search_runner.py` uses one global GOI identity stringency for every
target genome, so a divergent ortholog in a distant relative is dropped at the
same bar a close relative clears. §1g estimates each target's distance from the
median identity of its flanking-gene matches (data-driven, no tree/network) and
relaxes the HIGH/MEDIUM classify bars for distant targets only.

Pins:
  - the median-best-per-gene distance proxy,
  - the close (no-op) / mid (partial) / far (full + manual_review) tiers,
  - the no-op guarantees (disabled, too-few-flanking),
  - idempotency of the module-global mutation across worker reuse.
"""

import os
import sys
import types
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

import iterative_search_runner as isr  # noqa: E402
from iterative_search_runner import (  # noqa: E402
    CLASSIFY_THRESHOLDS,
    _median_best_flanking_identity,
    _apply_distance_adaptive_thresholds,
)


def _args(**overrides):
    base = dict(
        classify_high_min_identity=50.0,
        classify_medium_min_identity=35.0,
        disable_distance_autotune=False,
        distance_autotune_close_pct=70.0,
        distance_autotune_far_pct=40.0,
        distance_autotune_max_relax=10.0,
        distance_autotune_min_flanking=3,
    )
    base.update(overrides)
    return types.SimpleNamespace(**base)


def _hits(pairs):
    """pairs = list of (query, pident)."""
    return [{"query": q, "pident": p} for q, p in pairs]


class TestMedianBestFlankingIdentity(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(_median_best_flanking_identity([]), (None, 0))

    def test_collapses_to_best_per_gene(self):
        # gene A has two HSPs; only its best (90) should count, not drag the median.
        hits = _hits([("A", 50.0), ("A", 90.0), ("B", 80.0), ("C", 70.0)])
        median, n = _median_best_flanking_identity(hits)
        self.assertEqual(n, 3)            # 3 distinct genes
        self.assertEqual(median, 80.0)    # median of [70, 80, 90]

    def test_even_count_averages(self):
        median, n = _median_best_flanking_identity(_hits([("A", 40.0), ("B", 60.0)]))
        self.assertEqual(n, 2)
        self.assertEqual(median, 50.0)


class TestApplyDistanceAdaptiveThresholds(unittest.TestCase):
    def setUp(self):
        # Restore baseline globals before each test (mirrors main()'s population).
        CLASSIFY_THRESHOLDS["high_min_identity"] = 50.0
        CLASSIFY_THRESHOLDS["medium_min_identity"] = 35.0

    def test_close_tier_is_noop(self):
        hits = _hits([("A", 88.0), ("B", 85.0), ("C", 90.0)])
        info = _apply_distance_adaptive_thresholds(_args(), hits, "Close_sp")
        self.assertEqual(info["tier"], "close")
        self.assertEqual(info["relax"], 0.0)
        self.assertFalse(info["manual_review"])
        self.assertEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 50.0)
        self.assertEqual(CLASSIFY_THRESHOLDS["medium_min_identity"], 35.0)

    def test_far_tier_full_relax_and_manual_review(self):
        hits = _hits([("A", 30.0), ("B", 35.0), ("C", 38.0), ("D", 33.0)])
        info = _apply_distance_adaptive_thresholds(_args(), hits, "Far_sp")
        self.assertEqual(info["tier"], "far")
        self.assertEqual(info["relax"], 10.0)
        self.assertTrue(info["manual_review"])
        self.assertEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 40.0)
        self.assertEqual(CLASSIFY_THRESHOLDS["medium_min_identity"], 25.0)

    def test_mid_tier_partial_relax(self):
        # median 55 -> halfway between 40 and 70 -> 5pt relax
        hits = _hits([("A", 55.0), ("B", 55.0), ("C", 55.0)])
        info = _apply_distance_adaptive_thresholds(_args(), hits, "Mid_sp")
        self.assertEqual(info["tier"], "mid")
        self.assertAlmostEqual(info["relax"], 5.0, places=6)
        self.assertFalse(info["manual_review"])
        self.assertAlmostEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 45.0, places=6)
        self.assertAlmostEqual(CLASSIFY_THRESHOLDS["medium_min_identity"], 30.0, places=6)

    def test_disabled_is_noop(self):
        hits = _hits([("A", 30.0), ("B", 30.0), ("C", 30.0)])
        info = _apply_distance_adaptive_thresholds(_args(disable_distance_autotune=True), hits, "X")
        self.assertEqual(info["tier"], "disabled")
        self.assertEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 50.0)
        self.assertEqual(CLASSIFY_THRESHOLDS["medium_min_identity"], 35.0)

    def test_insufficient_flanking_keeps_baseline(self):
        hits = _hits([("A", 30.0), ("B", 30.0)])  # only 2 genes, min is 3
        info = _apply_distance_adaptive_thresholds(_args(), hits, "X")
        self.assertEqual(info["tier"], "insufficient_evidence")
        self.assertEqual(info["n_flanking"], 2)
        self.assertEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 50.0)
        self.assertEqual(CLASSIFY_THRESHOLDS["medium_min_identity"], 35.0)

    def test_idempotent_across_worker_reuse(self):
        # Simulate a pool worker handling a far genome then a close one: the
        # second call must reset from the args baseline, not the far-relaxed state.
        far = _hits([("A", 30.0), ("B", 32.0), ("C", 31.0)])
        close = _hits([("A", 85.0), ("B", 88.0), ("C", 90.0)])
        _apply_distance_adaptive_thresholds(_args(), far, "Far")
        self.assertEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 40.0)
        _apply_distance_adaptive_thresholds(_args(), close, "Close")
        self.assertEqual(CLASSIFY_THRESHOLDS["high_min_identity"], 50.0)
        self.assertEqual(CLASSIFY_THRESHOLDS["medium_min_identity"], 35.0)

    def test_respects_preset_adjusted_baseline(self):
        # §1f may have already lowered the bars (e.g. short-peptide preset).
        # §1g relaxation should subtract from THAT baseline, not the hardcoded 50.
        args = _args(classify_high_min_identity=40.0, classify_medium_min_identity=25.0)
        far = _hits([("A", 30.0), ("B", 32.0), ("C", 31.0)])
        info = _apply_distance_adaptive_thresholds(args, far, "Far")
        self.assertEqual(info["high_min_identity"], 30.0)   # 40 - 10
        self.assertEqual(info["medium_min_identity"], 15.0)  # 25 - 10


if __name__ == "__main__":
    unittest.main()
