#!/usr/bin/env python3
"""Regression tests for the home-genome quality gate (docs/TODO.md §1l).

Bug: get_assembly_quality() parsed NCBI assembly docsums via positional
`xtract -element ...`. NCBI docsums carry no top-level <ScaffoldCount>/<ContigCount>
(those live in the Meta/Stats block), so xtract collapsed the missing fields and the
present <ScaffoldN50> slid into the scaffold_count column. For Photinus pyralis
GCF_008802855.1 that meant a 47 Mb scaffold N50 was read as "scaffold_count=47017841",
which exceeds bad_max_scaffolds=500000, so the default `drop` policy rejected an
excellent reference assembly and the easy-mode run crashed at FETCH_HOME_GENOME.

Fix: `xtract -def NA` emits a placeholder for absent elements, keeping columns
aligned. The counts then arrive as NA (-> None) and the gate relies on N50 +
assembly level. These tests exercise the pure parsing/decision logic (no network,
no entrez-direct), so they run in CI.
"""
import os
import sys
import types
import unittest

BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))
sys.path.insert(0, BIN_DIR)

from fetch_home_genome import parse_quality_line, is_bad_quality  # noqa: E402


def _args(max_contigs=500000, max_scaffolds=500000, min_n50=5000):
    return types.SimpleNamespace(
        bad_max_contigs=max_contigs,
        bad_max_scaffolds=max_scaffolds,
        bad_min_n50=min_n50,
    )


# Real `xtract ... -def NA` output for GCF_008802855.1 (Photinus pyralis):
# no top-level ScaffoldCount/ContigCount -> NA; ScaffoldN50=47 Mb, ContigN50=170 kb.
PHOTINUS_LINE = "GCF_008802855.1\tPhotinus pyralis\treference genome\tScaffold\tNA\tNA\t47017841\t170308\tNA\tNA"


class TestParseQualityLine(unittest.TestCase):
    def test_scaffold_n50_not_misread_as_scaffold_count(self):
        e = parse_quality_line(PHOTINUS_LINE)
        # The whole point: 47017841 is the 47 Mb N50, NOT a scaffold count.
        self.assertIsNone(e["scaffold_count"])
        self.assertIsNone(e["contig_count"])
        self.assertEqual(e["scaffold_n50"], 47017841)
        self.assertEqual(e["contig_n50"], 170308)
        self.assertEqual(e["assembly_status"], "Scaffold")

    def test_real_counts_parse_when_present(self):
        # If a record DID expose counts, they must still parse as integers.
        line = "GCA_x\tFoo bar\tna\tScaffold\t2160\t7909\t47017841\t170308\tNA\tNA"
        e = parse_quality_line(line)
        self.assertEqual(e["scaffold_count"], 2160)
        self.assertEqual(e["contig_count"], 7909)

    def test_short_line_is_padded(self):
        e = parse_quality_line("GCF_1\tFoo bar")  # only 2 columns
        self.assertEqual(e["accession"], "GCF_1")
        self.assertIsNone(e["scaffold_n50"])


class TestQualityGate(unittest.TestCase):
    def test_good_scaffold_assembly_not_rejected(self):
        e = parse_quality_line(PHOTINUS_LINE)
        bad, reasons = is_bad_quality(e, _args())
        self.assertFalse(bad, f"excellent reference assembly wrongly flagged: {reasons}")

    def test_chromosome_level_always_accepted(self):
        line = "GCF_z\tFoo bar\treference genome\tChromosome\tNA\tNA\tNA\tNA\tNA\tNA"
        e = parse_quality_line(line)
        bad, _ = is_bad_quality(e, _args())
        self.assertFalse(bad)

    def test_fragmented_assembly_still_flagged_by_low_n50(self):
        # A genuinely poor contig-level assembly (N50 800 bp) must still be caught
        # via the N50 check, since the count check is now N50-driven.
        line = "GCA_y\tFoo bar\tna\tContig\tNA\tNA\tNA\t800\tNA\tNA"
        e = parse_quality_line(line)
        bad, reasons = is_bad_quality(e, _args())
        self.assertTrue(bad)
        self.assertTrue(any("N50" in r for r in reasons), reasons)


if __name__ == "__main__":
    unittest.main()
