#!/usr/bin/env python3
"""
A2 + A3 regression tests (TODO_JUN).

A2 — coverage-aware tandem classification: a tandem_copy reaches MEDIUM only when
identity AND real query coverage clear their floors; a short high-identity local
window (the Vollenhovia/Polistes "79%"/"71%" overcalls, now carrying their true
low qcov via QW2) is demoted to LOW.

A3 — the opt-in flanking-supported seeding was replaced on 2026-10-14 by one seed
rule (label rework, step 2): per query and genome the HIGH calls, else the MEDIUM calls,
each only when specific against the home proteome. Its tests are TestSeedRule below.
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

    def _conf(self, identity, qcov, flanking=8):
        # Default: a strongly supported block, so these tests isolate identity/coverage.
        return _classify_goi_evidence("tandem_copy", identity=identity, query_cov=qcov,
                                      flanking_support=flanking)

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


class TestTandemContext(unittest.TestCase):
    """Family benchmark 2026-09-11: a tandem copy needs context, like every other
    evidence type -- strong block flanking, or the sequence covering the query."""

    def setUp(self):
        self._saved = dict(CLASSIFY_THRESHOLDS)
        CLASSIFY_THRESHOLDS.update(tandem_min_identity=40.0, tandem_min_qcov=0.35,
                                   fallback_strong_min_flanking=5, medium_min_qcov=0.65)

    def tearDown(self):
        CLASSIFY_THRESHOLDS.clear()
        CLASSIFY_THRESHOLDS.update(self._saved)

    def _conf(self, identity, qcov, flanking):
        return _classify_goi_evidence("tandem_copy", identity=identity, query_cov=qcov,
                                      flanking_support=flanking)

    def test_chance_window_hit_in_weak_block_is_low(self):
        """SECP junk: 50 % over 0.39 of a 77-aa query, block flanking 4 -> LOW."""
        conf, cls, reason = self._conf(50.0, 0.39, 4)
        self.assertEqual((conf, cls, reason),
                         ("LOW", "tandem_goi_copy", "goi_tandem_copy_weak_context"))

    def test_strong_block_keeps_partial_copy(self):
        """Melittin bee copy: 57 % over ~0.5, block flanking 8 -> MEDIUM."""
        self.assertEqual(self._conf(57.0, 0.5, 8)[0], "MEDIUM")

    def test_full_length_copy_needs_no_flanking(self):
        """KAZA copy in a 2-flanking block: 87 % over 0.87 -> MEDIUM on sequence alone."""
        self.assertEqual(self._conf(87.4, 0.869, 2)[0], "MEDIUM")

    def test_context_bars_are_inclusive(self):
        self.assertEqual(self._conf(40.0, 0.35, 5)[0], "MEDIUM")
        self.assertEqual(self._conf(40.0, 0.65, 0)[0], "MEDIUM")
        self.assertEqual(self._conf(40.0, 0.64, 4)[0], "LOW")

    def test_identity_and_coverage_floors_still_come_first(self):
        self.assertEqual(self._conf(39.9, 0.9, 10)[2], "goi_tandem_copy_low_identity")
        self.assertEqual(self._conf(90.0, 0.30, 10)[2], "goi_tandem_copy_low_coverage")


class TestSeedRule(unittest.TestCase):
    """Seeds of the next wave (label rework 2026-10-14, step 2; replaced the A3 opt-in).

    Per query and genome: the HIGH calls if there is one, else the MEDIUM calls, each only
    when it is specific (no real competitor in the home proteome, or a margin of >= 10 bits
    over it). The class filter is gone: a tandem copy seeds like any call."""

    def _seed(self, meta):
        genes = [_gene(gid) for gid in meta]
        return _classify_goi_for_seed_and_tree(genes, meta, _is_goi)

    @staticmethod
    def _call(conf, parent="GOI_Melt", goi_class="probable_goi", margin="", real="", promoted=""):
        return {"role": "goi", "confidence": conf, "goi_class": goi_class, "parent": parent,
                "spec_margin": margin, "real_competitor": real, "context_promoted": promoted}

    def test_tandem_copy_seeds_now(self):
        """Apis cerana's melittin copy (tandem_goi_copy) seeds: the class filter is gone."""
        seed, tree_extra, suppressed = self._seed({
            "GOI_copy_1|Apis_cerana_fna_b0_l1": self._call("MEDIUM", goi_class="tandem_goi_copy",
                                                           margin="40.0", real="1")})
        self.assertEqual(len(seed), 1)
        self.assertEqual((tree_extra, suppressed), ([], 0))

    def test_high_takes_the_seat_from_a_medium_relative(self):
        """A genome with a HIGH call seeds only its HIGH calls; the paralog MEDIUM next to it
        (apyrase: three genomes, margins -339 to -363 bits) feeds the tree only."""
        seed, tree_extra, _ = self._seed({
            "GOI_A|g_b0_l1_exon_ann": self._call("HIGH", parent="GOI_A", margin="350.0", real="1"),
            "GOI_A|g_b0_l1m2_exon_ann": self._call("MEDIUM", parent="GOI_A", margin="-350.0", real="1"),
        })
        self.assertEqual([g["id"] for g in seed], ["GOI_A|g_b0_l1_exon_ann"])
        self.assertEqual([g["id"] for g in tree_extra], ["GOI_A|g_b0_l1m2_exon_ann"])

    def test_medium_seeds_where_there_is_no_high(self):
        seed, _t, _s = self._seed({"GOI_A|g_b0_l1": self._call("MEDIUM", parent="GOI_A",
                                                                margin="12.0", real="1")})
        self.assertEqual(len(seed), 1)

    def test_specificity_bar_is_10_bits_against_a_real_competitor(self):
        for margin, seeds in (("10.0", 1), ("9.9", 0), ("-1.2", 0)):
            seed, tree_extra, _ = self._seed({"GOI_A|g_b0_l1": self._call(
                "HIGH", parent="GOI_A", margin=margin, real="1")})
            self.assertEqual(len(seed), seeds, margin)
            self.assertEqual(len(tree_extra), 1 - seeds, margin)

    def test_no_real_competitor_seeds_whatever_the_margin(self):
        seed, _t, _s = self._seed({"GOI_A|g_b0_l1": self._call("MEDIUM", parent="GOI_A",
                                                                margin="-5.0", real="0")})
        self.assertEqual(len(seed), 1)

    def test_without_specificity_terms_the_label_decides(self):
        """No home proteome in the run: no terms, nothing held back."""
        seed, _t, _s = self._seed({"GOI_A|g_b0_l1": self._call("MEDIUM", parent="GOI_A")})
        self.assertEqual(len(seed), 1)

    def test_queries_are_pooled_apart(self):
        seed, tree_extra, _ = self._seed({
            "GOI_A|g_b0_l1": self._call("HIGH", parent="GOI_A", margin="50", real="1"),
            "GOI_B|g_b1_l1": self._call("MEDIUM", parent="GOI_B", margin="50", real="1"),
        })
        self.assertEqual(len(seed), 2)
        self.assertEqual(tree_extra, [])

    def test_low_and_ambiguous_never_seed(self):
        """Neither label seeds, whatever its basis (gene_order_only or not_specific), and
        neither reaches the tree."""
        seed, tree_extra, suppressed = self._seed({
            "GOI_A|low": self._call("LOW", parent="GOI_A"),
            "GOI_A|amb": self._call("AMBIGUOUS", parent="GOI_A", goi_class="syntenic_candidate_unconfirmed"),
        })
        self.assertEqual((seed, tree_extra, suppressed), ([], [], 2))

if __name__ == '__main__':
    unittest.main()
