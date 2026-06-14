#!/usr/bin/env python3
"""
A2 + A3 regression tests (TODO_JUN).

A2 — coverage-aware tandem classification: a tandem_copy reaches MEDIUM only when
identity AND real query coverage clear their floors; a short high-identity local
window (the Vollenhovia/Polistes "79%"/"71%" overcalls, now carrying their true
low qcov via QW2) is demoted to LOW.

A3 — un-neuter the wave seed (opt-in): when `seed_on_flanking_support` is on, a
HIGH/MEDIUM GOI feature seeds the expanding DB regardless of goi_class (incl.
tandem_goi_copy) provided flanking support + query coverage clear the floors — so
Apis cerana's 100%-identical melittin can bridge to the next wave. Default off ⇒
behaviour identical to before.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

import iterative_search_runner as isr
from iterative_search_runner import (
    _classify_goi_evidence,
    _classify_goi_for_seed_and_tree,
    CLASSIFY_THRESHOLDS,
)


def _gene(gid, seq="MKT"):
    return {"id": gid, "seq": seq}


def _is_goi(gid):
    return gid.startswith("GOI_")


class TestA2TandemCoverage(unittest.TestCase):
    """A2: tandem_copy MEDIUM needs identity AND coverage."""

    def setUp(self):
        # Pin thresholds to defaults in case another test mutated the global dict.
        self._saved = dict(CLASSIFY_THRESHOLDS)
        CLASSIFY_THRESHOLDS["tandem_min_identity"] = 40.0
        CLASSIFY_THRESHOLDS["tandem_min_qcov"] = 0.35

    def tearDown(self):
        CLASSIFY_THRESHOLDS.clear()
        CLASSIFY_THRESHOLDS.update(self._saved)

    def _conf(self, identity, qcov):
        return _classify_goi_evidence("tandem_copy", identity=identity, query_cov=qcov)

    def test_high_identity_low_coverage_demoted(self):
        """Vollenhovia: 79% identity but qcov 0.30 → LOW (low coverage)."""
        conf, cls, reason = self._conf(79.0, 0.30)
        self.assertEqual(conf, "LOW")
        self.assertEqual(cls, "tandem_goi_copy")
        self.assertEqual(reason, "goi_tandem_copy_low_coverage")

    def test_missing_coverage_demoted(self):
        """query_cov None (unknown) → LOW, not MEDIUM."""
        conf, _cls, reason = self._conf(100.0, None)
        self.assertEqual(conf, "LOW")
        self.assertEqual(reason, "goi_tandem_copy_low_coverage")

    def test_identity_and_coverage_pass(self):
        """Real full-length tandem copy: 79% / 0.5 → MEDIUM."""
        conf, cls, reason = self._conf(79.0, 0.5)
        self.assertEqual(conf, "MEDIUM")
        self.assertEqual(cls, "tandem_goi_copy")
        self.assertEqual(reason, "goi_tandem_copy_detected")

    def test_low_identity_still_low_identity_reason(self):
        """Low identity keeps the identity reason even if coverage is fine."""
        conf, _cls, reason = self._conf(30.0, 0.9)
        self.assertEqual(conf, "LOW")
        self.assertEqual(reason, "goi_tandem_copy_low_identity")

    def test_exactly_on_floors_passes(self):
        conf, _cls, _r = self._conf(40.0, 0.35)
        self.assertEqual(conf, "MEDIUM")


class TestA3FlankingSeed(unittest.TestCase):
    """A3: opt-in flanking-supported seeding of tandem/divergent orthologs."""

    def _apis_cerana_meta(self, flanking="9", qcov="1.0"):
        return {
            "GOI_copy_1|Apis_cerana_fna_b0_l1": {
                "role": "goi", "confidence": "MEDIUM", "goi_class": "tandem_goi_copy",
                "block_flanking_support": flanking, "query_coverage": qcov,
            },
        }

    def test_default_off_tandem_does_not_seed(self):
        """Default (flag off): tandem copy goes to tree, not seed — unchanged."""
        genes = [_gene("GOI_copy_1|Apis_cerana_fna_b0_l1")]
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            genes, self._apis_cerana_meta(), _is_goi
        )
        self.assertEqual(seed, [])
        self.assertEqual(len(tree_extra), 1)
        self.assertEqual(suppressed, 1)

    def test_flag_on_strong_flanking_seeds(self):
        """Flag on + strong flanking + high coverage → tandem copy seeds."""
        genes = [_gene("GOI_copy_1|Apis_cerana_fna_b0_l1")]
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            genes, self._apis_cerana_meta(), _is_goi,
            seed_on_flanking_support=True, seed_min_flanking=2, seed_min_qcov=0.5,
        )
        self.assertEqual(len(seed), 1)
        self.assertEqual(tree_extra, [])
        self.assertEqual(suppressed, 0)

    def test_flag_on_weak_flanking_does_not_seed(self):
        """Flag on but flanking below floor → still withheld (tree only)."""
        genes = [_gene("GOI_copy_1|Apis_cerana_fna_b0_l1")]
        seed, tree_extra, _ = _classify_goi_for_seed_and_tree(
            genes, self._apis_cerana_meta(flanking="1"), _is_goi,
            seed_on_flanking_support=True, seed_min_flanking=2, seed_min_qcov=0.5,
        )
        self.assertEqual(seed, [])
        self.assertEqual(len(tree_extra), 1)

    def test_flag_on_low_coverage_does_not_seed(self):
        """Flag on, strong flanking, but coverage below floor → withheld."""
        genes = [_gene("GOI_copy_1|Apis_cerana_fna_b0_l1")]
        seed, tree_extra, _ = _classify_goi_for_seed_and_tree(
            genes, self._apis_cerana_meta(qcov="0.3"), _is_goi,
            seed_on_flanking_support=True, seed_min_flanking=2, seed_min_qcov=0.5,
        )
        self.assertEqual(seed, [])
        self.assertEqual(len(tree_extra), 1)

    def test_flag_on_low_confidence_never_seeds(self):
        """A LOW-confidence feature never seeds even with strong flanking+coverage."""
        genes = [_gene("GOI_copy_1|noise")]
        meta = {
            "GOI_copy_1|noise": {
                "role": "goi", "confidence": "LOW", "goi_class": "ambiguous_goi_family_member",
                "block_flanking_support": "20", "query_coverage": "0.9",
            },
        }
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            genes, meta, _is_goi, seed_on_flanking_support=True,
        )
        self.assertEqual(seed, [])
        self.assertEqual(tree_extra, [])  # LOW never reaches the tree either
        self.assertEqual(suppressed, 1)

    def test_flag_on_missing_coverage_attr_does_not_seed(self):
        """Empty QueryCoverage attr (None) can't satisfy the coverage floor."""
        genes = [_gene("GOI_copy_1|Apis_cerana_fna_b0_l1")]
        seed, _tree, _s = _classify_goi_for_seed_and_tree(
            genes, self._apis_cerana_meta(qcov=""), _is_goi,
            seed_on_flanking_support=True,
        )
        self.assertEqual(seed, [])

    def test_confident_goi_seeds_regardless_of_flag(self):
        """Normal confident_goi path is untouched by A3."""
        genes = [_gene("GOI_Melt|Apis_florea_fna_b0_l1_exon_ann")]
        meta = {
            "GOI_Melt|Apis_florea_fna_b0_l1_exon_ann": {
                "role": "goi", "confidence": "HIGH", "goi_class": "confident_goi",
            },
        }
        for flag in (False, True):
            seed, _t, _s = _classify_goi_for_seed_and_tree(
                genes, meta, _is_goi, seed_on_flanking_support=flag
            )
            self.assertEqual(len(seed), 1)


if __name__ == '__main__':
    unittest.main()
