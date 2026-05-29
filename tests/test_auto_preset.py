#!/usr/bin/env python3
"""Tests for the hit-profile / auto-preset selection (docs/TODO.md §1f).

profile_hits.py summarises the LOCATE_GENE hit distribution; auto_select_preset.py
maps it to a query-type preset. The committed fixture tests/fixtures/luciferase_hits_blast.txt
is the real 185-hit Photinus pyralis luciferase BLAST output that motivated §1f — it
must land on preset_paralog_discrimination + WARN.
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
        self.assertEqual(prof["top_bit"], 415.0)
        self.assertEqual(prof["second_bit"], 354.0)
        self.assertAlmostEqual(prof["bit_ratio"], 415.0 / 354.0, places=2)
        # §1d kept 4 loci on this data; the strong-locus count should agree.
        self.assertEqual(prof["n_strong_loci"], 4)
        self.assertGreater(prof["n_independent_loci"], 15)
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


class TestSelectPreset(unittest.TestCase):
    def test_luciferase_is_paralog_discrimination_with_warning(self):
        hits = profile_hits.parse_hits([FIXTURE])
        prof = profile_hits.build_profile(hits, None, 100000, 0.5)
        preset, warn, reasons = select_preset(prof)
        self.assertEqual(preset, "preset_paralog_discrimination")
        self.assertTrue(warn)
        self.assertTrue(any("185" in r for r in reasons))

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
        prof = {"query_len": 120, "n_hits": 8, "n_independent_loci": 5,
                "n_strong_loci": 3, "bit_ratio": 1.25}
        self.assertEqual(select_preset(prof)[0], "preset_tandem_family")

    def test_no_hits_defaults_single_copy(self):
        prof = {"query_len": 300, "n_hits": 0}
        self.assertEqual(select_preset(prof)[0], "preset_single_copy")


if __name__ == "__main__":
    unittest.main()
