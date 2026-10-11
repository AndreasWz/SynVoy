#!/usr/bin/env python3
"""Tests for the hit-profile / auto-preset selection (docs/TODO.md §1f).

profile_hits.py summarises the LOCATE_GENE hit distribution; auto_select_preset.py
maps it to a query-type preset. The committed fixture tests/fixtures/luciferase_hits_blast.txt
is the real 185-hit Photinus pyralis luciferase BLAST output that motivated §1f — it
must land on preset_paralog_discrimination + WARN.

Since 2026-10-08 the profile counts places, not alignment rows: a locus is scored by one
copy of the gene (all its exons), and the number of copies at the strongest locus is the
depth of query coverage there. Before, more than 15 rows meant "large paralog family",
which every single-copy gene with eight exons or more has.
"""

import json
import os
import sys
import tempfile
import unittest

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

import profile_hits  # noqa: E402
from auto_select_preset import select_preset  # noqa: E402

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "luciferase_hits_blast.txt")


class TestProfileHits(unittest.TestCase):
    def test_luciferase_profile_matches_known_numbers(self):
        hits = profile_hits.parse_hits([FIXTURE])
        self.assertEqual(len(hits), 185)
        prof = profile_hits.build_profile(hits, query_len=None, window=100000, strong_frac=0.5)
        self.assertEqual(prof["n_hits"], 185)
        self.assertEqual(prof["n_hits_raw"], 185)
        # A locus scores one copy of the gene: the best set of hits that covers the query
        # once (826 bits), not its best single exon (415 bits, the value before).
        self.assertEqual(prof["top_bit"], 826.3)
        self.assertEqual(prof["second_bit"], 618.0)
        self.assertAlmostEqual(prof["bit_ratio"], 826.3 / 618.0, places=2)
        # §1d kept 4 loci on this data; the strong-locus count should agree.
        self.assertEqual(prof["n_strong_loci"], 4)
        self.assertGreater(prof["n_independent_loci"], 15)
        # Several luciferase-like genes lie next to each other at the strongest locus.
        self.assertGreater(prof["primary_copies"], 4)
        # query_len falls back to max qend when no FASTA is given.
        self.assertEqual(prof["query_len"], 550)
        # Gradual identity decay: wide spread from a near-identical top to divergent paralogs.
        self.assertLess(prof["identity_median"], 50)
        self.assertGreater(prof["identity_max"], 90)

    def test_query_len_from_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "q.faa")
            with open(fa, "w") as fh:
                fh.write(">q\n" + "M" * 42 + "\n")
            self.assertEqual(profile_hits._read_query_len(fa), 42)

    def test_empty_hits(self):
        prof = profile_hits.build_profile([], query_len=300, window=100000, strong_frac=0.5)
        self.assertEqual(prof["n_hits"], 0)
        self.assertIsNone(prof["top_bit"])

    @staticmethod
    def _hit(chrom, tstart, qstart, qend, bits, evalue=1e-20, pident=90.0):
        return {"chrom": chrom, "tstart": tstart, "tend": tstart + 3 * (qend - qstart + 1),
                "qstart": qstart, "qend": qend, "bits": bits, "evalue": evalue, "pident": pident}

    def _exons(self, chrom, start, n=11, exon_aa=50, bits=60.0):
        """One copy of a gene with n exons, each reported by both search programs."""
        hits = []
        for i in range(n):
            q0 = 1 + i * exon_aa
            for _program in range(2):
                hits.append(self._hit(chrom, start + i * 1000, q0, q0 + exon_aa - 1, bits))
        return hits

    def test_many_exons_are_one_copy_at_one_locus(self):
        """An 11-exon single-copy gene: 22 rows, one locus, covered once."""
        prof = profile_hits.build_profile(self._exons("chr1", 10_000), 550, 100000, 0.5)
        self.assertEqual(prof["n_hits"], 22)
        self.assertEqual(prof["n_independent_loci"], 1)
        self.assertEqual(prof["n_strong_loci"], 1)
        self.assertAlmostEqual(prof["primary_copies"], 1.0, places=2)
        self.assertEqual(prof["top_bit"], 11 * 60.0)
        self.assertEqual(select_preset(prof)[0], "preset_single_copy")

    def test_two_copies_next_to_each_other_are_a_tandem_array(self):
        hits = self._exons("chr1", 10_000) + self._exons("chr1", 40_000, bits=55.0)
        prof = profile_hits.build_profile(hits, 550, 100000, 0.5)
        self.assertEqual(prof["n_independent_loci"], 1)
        self.assertAlmostEqual(prof["primary_copies"], 2.0, places=2)
        # the score of the locus is ONE copy, not the sum of both
        self.assertEqual(prof["top_bit"], 11 * 60.0)
        self.assertEqual(select_preset(prof)[0], "preset_tandem_family")

    def test_a_relative_elsewhere_is_scored_as_a_whole_gene(self):
        """One well-conserved exon elsewhere is not a strong locus; a full relative is."""
        home = self._exons("chr1", 10_000)
        one_exon = [self._hit("chr5", 500_000, 101, 150, 58.0)]
        prof = profile_hits.build_profile(home + one_exon, 550, 100000, 0.5)
        self.assertEqual(prof["n_independent_loci"], 2)
        self.assertEqual(prof["n_strong_loci"], 1)
        self.assertEqual(select_preset(prof)[0], "preset_single_copy")
        relative = self._exons("chr7", 900_000, bits=40.0)
        prof = profile_hits.build_profile(home + relative, 550, 100000, 0.5)
        self.assertEqual(prof["n_strong_loci"], 2)
        self.assertEqual(select_preset(prof)[0], "preset_paralog_discrimination")

    def test_hits_above_the_evalue_limit_are_not_profiled(self):
        hits = self._exons("chr1", 10_000) + [
            self._hit(f"chr{i}", 1000, 1, 30, 25.0, evalue=3.0) for i in range(2, 30)]
        prof = profile_hits.build_profile(hits, 550, 100000, 0.5, max_evalue=0.01)
        self.assertEqual(prof["n_hits_raw"], 22 + 28)
        self.assertEqual(prof["n_hits"], 22)
        self.assertEqual(prof["n_independent_loci"], 1)


class TestSelectPreset(unittest.TestCase):
    def test_luciferase_is_paralog_discrimination_with_warning(self):
        hits = profile_hits.parse_hits([FIXTURE])
        prof = profile_hits.build_profile(hits, None, 100000, 0.5)
        preset, warn, reasons = select_preset(prof)
        self.assertEqual(preset, "preset_paralog_discrimination")
        self.assertTrue(warn)
        self.assertTrue(any("4 strong loci among 50" in r for r in reasons))

    def test_the_number_of_rows_decides_nothing(self):
        """26 rows of one 11-exon gene used to mean 'large paralog family'."""
        prof = {"query_len": 576, "n_hits": 26, "n_independent_loci": 1,
                "n_strong_loci": 1, "bit_ratio": None, "primary_copies": 1.0}
        self.assertEqual(select_preset(prof), ("preset_single_copy", False, select_preset(prof)[2]))

    def test_short_peptide_overrides(self):
        # A short query wins even with a paralog-like hit count.
        prof = {"query_len": 26, "n_hits": 40, "n_independent_loci": 12,
                "n_strong_loci": 3, "bit_ratio": 1.1}
        preset, warn, _ = select_preset(prof)
        self.assertEqual(preset, "preset_short_peptide")
        self.assertFalse(warn)

    def test_single_copy(self):
        prof = {"query_len": 300, "n_hits": 2, "n_independent_loci": 1,
                "n_strong_loci": 1, "bit_ratio": 5.0}
        self.assertEqual(select_preset(prof)[0], "preset_single_copy")

    def test_tandem_family(self):
        """Several copies at the strongest locus, no comparable locus elsewhere."""
        prof = {"query_len": 120, "n_hits": 8, "n_independent_loci": 2,
                "n_strong_loci": 1, "bit_ratio": 6.0, "primary_copies": 2.3}
        self.assertEqual(select_preset(prof)[0], "preset_tandem_family")

    def test_relatives_elsewhere_are_paralog_discrimination(self):
        prof = {"query_len": 120, "n_hits": 8, "n_independent_loci": 5,
                "n_strong_loci": 3, "bit_ratio": 1.25, "primary_copies": 1.0}
        preset, warn, _ = select_preset(prof)
        self.assertEqual(preset, "preset_paralog_discrimination")
        self.assertFalse(warn)          # five places are not a large family
        prof["n_independent_loci"] = 40
        self.assertTrue(select_preset(prof)[1])

    def test_no_hits_defaults_single_copy(self):
        prof = {"query_len": 300, "n_hits": 0}
        self.assertEqual(select_preset(prof)[0], "preset_single_copy")


if __name__ == "__main__":
    unittest.main()
