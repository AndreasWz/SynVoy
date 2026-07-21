#!/usr/bin/env python3
"""Tests for §1o dispersed-GOI rescue (bin/iterative_search_runner.py).

Recovers a strong GOI hit that dispersed WITHIN the flanking envelope of its
chromosome (a stretched-but-real neighbourhood), while staying a no-op on
well-conserved targets and rejecting anchorless paralog / fold-analogue hits.

Load-bearing REAL case: oskar_ladder run 5702852, Aedes aegypti. The GOI query
hit the annotated oskar (LOC5577338, NC_035109.1:184,139,916-184,154,866) at
45.3% id / e-2.77e-29; on that chromosome 11 flanking anchors span 15-243 Mb and
oskar sits inside that envelope, but no single anchor is within 62 Mb, so no
block formed there and 0 GOI models were emitted.
"""

import os
import sys
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

import iterative_search_runner as isr  # noqa: E402


def _flank(gene, chrom, start, end):
    return {"query": f"gene-Dmel_{gene}", "chrom": chrom,
            "start": start, "end": end, "pident": 80.0,
            "evalue": 1e-40, "alnlen": 300, "strand": "+", "bits": 400.0}


# Mirrors the real Aedes anchors: a stretched envelope on the oskar chromosome
# (>=3 anchors spanning tens of Mb) plus a couple on other chromosomes.
AEDES_FLANKING = [
    _flank("CG11967", "NC_035109.1", 15_100_000, 15_106_000),
    _flank("CG11968", "NC_035109.1", 121_344_000, 121_351_000),
    _flank("CG11963", "NC_035109.1", 121_366_000, 121_405_000),
    _flank("CG8379",  "NC_035109.1", 203_400_000, 203_460_000),
    _flank("CG9836",  "NC_035109.1", 243_500_000, 243_506_000),
    _flank("CG11971", "NC_035107.1", 194_442_000, 194_442_400),
    _flank("CG11975", "NC_035108.1", 335_570_000, 335_571_000),
]

# The real dispersed oskar hit — inside the NC_035109.1 envelope, 62 Mb from
# the nearest anchor.
AEDES_OSKAR_HIT = {"query": "GOI_osk", "chrom": "NC_035109.1",
                   "start": 184_153_826, "end": 184_154_285, "pident": 45.3,
                   "evalue": 2.769e-29, "alnlen": 459, "strand": "+", "bits": 130.0}


class TestRearrangementScoreDiagnostic(unittest.TestCase):
    def test_in_range(self):
        s = isr.compute_rearrangement_score(AEDES_FLANKING)
        self.assertGreaterEqual(s, 0.0)
        self.assertLessEqual(s, 1.0)

    def test_too_few_flanking_is_zero(self):
        self.assertEqual(isr.compute_rearrangement_score(AEDES_FLANKING[:2]), 0.0)

    def test_best_hit_per_gene_collapses_noise(self):
        # two hits of the same gene, different scores -> one canonical anchor
        noisy = AEDES_FLANKING + [
            dict(_flank("CG11967", "NC_035199.9", 5, 305), bits=10.0)]  # weaker, spurious
        pos = isr._flanking_anchor_positions(noisy)
        # CG11967's best hit (bits 400) wins -> it stays on NC_035109.1, not the junk chrom
        genes_on_junk = [g for g, *_ in pos.get("NC_035199.9", [])]
        self.assertNotIn("gene-Dmel_CG11967", genes_on_junk)


class TestSelectDispersedGoiSeeds(unittest.TestCase):
    def test_recovers_the_real_aedes_oskar(self):
        seeds = isr.select_dispersed_goi_seeds(
            AEDES_FLANKING + [AEDES_OSKAR_HIT], AEDES_FLANKING)
        self.assertEqual(len(seeds), 1)
        self.assertEqual(seeds[0]["query"], "GOI_osk")
        self.assertEqual(seeds[0]["chrom"], "NC_035109.1")

    def test_rejects_chromosome_without_enough_anchors(self):
        # Only 2 anchors on the oskar chromosome -> no envelope (biglycan-on-chrX
        # / structdisc analogue: a strong hit with no macro-synteny envelope).
        few = [h for h in AEDES_FLANKING if h["chrom"] == "NC_035109.1"][:2]
        seeds = isr.select_dispersed_goi_seeds(few + [AEDES_OSKAR_HIT], few)
        self.assertEqual(seeds, [])

    def test_rejects_hit_outside_envelope(self):
        far = dict(AEDES_OSKAR_HIT, start=500_000_000, end=500_000_459)
        seeds = isr.select_dispersed_goi_seeds(
            AEDES_FLANKING + [far], AEDES_FLANKING)
        self.assertEqual(seeds, [])

    def test_rejects_weak_identity(self):
        weak = dict(AEDES_OSKAR_HIT, pident=28.0)   # fold-analogue level
        seeds = isr.select_dispersed_goi_seeds(
            AEDES_FLANKING + [weak], AEDES_FLANKING)
        self.assertEqual(seeds, [])

    def test_rejects_insignificant_evalue(self):
        weak = dict(AEDES_OSKAR_HIT, evalue=1e-3)
        seeds = isr.select_dispersed_goi_seeds(
            AEDES_FLANKING + [weak], AEDES_FLANKING)
        self.assertEqual(seeds, [])

    def test_conserved_genome_is_a_noop(self):
        # Tight neighbourhood (all anchors within ~40 kb) + the real GOI right
        # next to them: the only strong hit is within region_padding of an anchor,
        # so nothing new is seeded (the normal block path already models it).
        tight = [_flank(f"CG{i}", "chr1", 1_000_000 + i * 8000, 1_000_000 + i * 8000 + 500)
                 for i in range(5)]
        goi_near = {"query": "GOI_x", "chrom": "chr1", "start": 1_020_000,
                    "end": 1_020_459, "pident": 90.0, "evalue": 1e-50,
                    "alnlen": 400, "strand": "+", "bits": 800.0}
        seeds = isr.select_dispersed_goi_seeds(
            tight + [goi_near], tight, region_padding=50000)
        self.assertEqual(seeds, [])

    def test_skips_hit_already_next_to_an_anchor(self):
        near = dict(AEDES_OSKAR_HIT, start=121_360_000, end=121_360_459)
        seeds = isr.select_dispersed_goi_seeds(
            AEDES_FLANKING + [near], AEDES_FLANKING, region_padding=50000)
        self.assertEqual(seeds, [])

    def test_dedupes_multiple_hits_in_same_window(self):
        h2 = dict(AEDES_OSKAR_HIT, start=184_153_900, end=184_154_200)
        seeds = isr.select_dispersed_goi_seeds(
            AEDES_FLANKING + [AEDES_OSKAR_HIT, h2], AEDES_FLANKING)
        self.assertEqual(len(seeds), 1)


if __name__ == "__main__":
    unittest.main()
