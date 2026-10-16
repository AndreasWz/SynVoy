#!/usr/bin/env python3
"""
A3 regression tests (TODO_JUN): the seed rule.

The opt-in flanking-supported seeding was replaced on 2026-10-14 by one seed rule (label
rework, step 2): per query and genome the HIGH calls, else the MEDIUM calls, each only when
specific against the home proteome. Its tests are TestSeedRule below.

The A2 tests (a tandem copy reached MEDIUM only with identity AND coverage, and with
context) were removed on 2026-10-15 with the identity-based label: a tandem copy is
labelled by synvoy_label.judge like every GOI call.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

import iterative_search_runner as isr
from iterative_search_runner import _classify_goi_for_seed_and_tree


def _gene(gid, seq="MKT"):
    return {"id": gid, "seq": seq}


def _is_goi(gid):
    return gid.startswith("GOI_")


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
