#!/usr/bin/env python3
"""Tests for the LOW-confidence GOI filter in load_goi_intervals_from_gff
(roadmap Item E — Melipona contig_140 false-positive regression).

Background:
    iterative_search_runner.py classifies each GOI hit as HIGH/MEDIUM/LOW
    based on identity, query coverage, and flanking support. Prior to this
    fix, cluster_grs.py:load_goi_intervals_from_gff() loaded ALL GOI rows
    regardless of confidence, treating them as equally weighted synteny
    anchors. The Melipona case showed this produced HIGH/S1.00 anchor
    regions on a contig where the underlying GOI hits were all LOW.

This test covers:
    1. Rows with Confidence=LOW are excluded from goi_intervals.
    2. Rows with Confidence=MEDIUM/HIGH are kept.
    3. Rows without a Confidence attribute (legacy GFF) are kept (back-compat).
"""

import os
import sys
import tempfile
import unittest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

from cluster_grs import load_goi_intervals_from_gff  # noqa: E402


def _gff_row(chrom, ftype, start, end, attrs):
    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())
    return "\t".join([chrom, "test_src", ftype, str(start), str(end),
                      ".", "+", ".", attr_str]) + "\n"


class TestLowConfidenceGoiFilter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gff = os.path.join(self.tmp, "iterative.gff")

    def _write(self, lines):
        with open(self.gff, "w") as fh:
            fh.write("##gff-version 3\n")
            for ln in lines:
                fh.write(ln)

    def test_low_confidence_rows_are_skipped(self):
        # Melipona contig_140 pattern: 3 LOW GOI hits clustered closely.
        self._write([
            _gff_row("contig_140", "gene", 1505499, 1505591, {
                "ID": "GOI_copy_2|...", "Name": "GOI_copy_2",
                "SynVoy_Parent": "GOI_Melt", "Confidence": "LOW",
                "EvidenceType": "tandem_copy",
            }),
            _gff_row("contig_140", "mRNA", 1587903, 1588001, {
                "ID": "GOI_Melt|fb", "Name": "GOI_Melt",
                "SynVoy_Parent": "GOI_Melt", "Confidence": "LOW",
                "EvidenceType": "fallback_hit_span",
            }),
            _gff_row("contig_140", "gene", 1498789, 1498908, {
                "ID": "GOI_copy_1|...", "Name": "GOI_copy_1",
                "SynVoy_Parent": "GOI_Melt", "Confidence": "LOW",
                "EvidenceType": "tandem_copy",
            }),
        ])
        intervals = load_goi_intervals_from_gff(self.gff, padding_bp=20000)
        self.assertEqual(intervals, [],
                         msg=f"Expected no intervals from all-LOW GFF, got {intervals}")

    def test_medium_and_high_rows_are_kept(self):
        self._write([
            _gff_row("chr1", "mRNA", 1000, 1100, {
                "ID": "GOI_a|x", "Name": "GOI_a",
                "SynVoy_Parent": "GOI_Q", "Confidence": "MEDIUM",
            }),
            _gff_row("chr1", "gene", 5000, 5100, {
                "ID": "GOI_b|x", "Name": "GOI_b",
                "SynVoy_Parent": "GOI_Q", "Confidence": "HIGH",
            }),
        ])
        intervals = load_goi_intervals_from_gff(self.gff, padding_bp=0)
        self.assertEqual(len(intervals), 2,
                         msg=f"Expected 2 intervals from MEDIUM+HIGH, got {intervals}")

    def test_missing_confidence_attribute_is_kept_for_back_compat(self):
        # Legacy GFF without Confidence — should still be loaded.
        self._write([
            _gff_row("chr1", "mRNA", 1000, 1100, {
                "ID": "GOI_legacy|x", "Name": "GOI_legacy",
                "SynVoy_Parent": "GOI_Q",
            }),
        ])
        intervals = load_goi_intervals_from_gff(self.gff, padding_bp=0)
        self.assertEqual(len(intervals), 1)

    def test_mixed_confidence_only_keeps_non_low(self):
        # Real-world: some rows LOW, some MEDIUM. Only MEDIUM survives.
        self._write([
            _gff_row("chr1", "gene", 100, 200, {
                "ID": "GOI_low|x", "Name": "GOI_low",
                "SynVoy_Parent": "GOI_Q", "Confidence": "LOW",
            }),
            _gff_row("chr1", "mRNA", 5000, 5100, {
                "ID": "GOI_med|x", "Name": "GOI_med",
                "SynVoy_Parent": "GOI_Q", "Confidence": "MEDIUM",
            }),
            _gff_row("chr1", "gene", 10000, 10200, {
                "ID": "GOI_low2|x", "Name": "GOI_low2",
                "SynVoy_Parent": "GOI_Q", "Confidence": "low",  # case-insensitive
            }),
        ])
        intervals = load_goi_intervals_from_gff(self.gff, padding_bp=0)
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0]["start"], 5000)
        self.assertEqual(intervals[0]["end"], 5100)
        # Confidence should be carried through for anchor-score scaling
        self.assertEqual(intervals[0]["confidence"], "MEDIUM")

    def test_confidence_carried_through_merge(self):
        # Two MEDIUM hits within padding distance merge → max conf preserved.
        # A HIGH hit + MEDIUM hit overlapping → merged interval gets HIGH.
        from cluster_grs import merge_intervals  # noqa: E402
        merged = merge_intervals([
            {"chrom": "c1", "start": 100, "end": 200, "confidence": "MEDIUM"},
            {"chrom": "c1", "start": 150, "end": 300, "confidence": "HIGH"},
        ])
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["confidence"], "HIGH")


class TestAnchorScoreScaling(unittest.TestCase):
    def test_anchor_score_for_confidence_table(self):
        from cluster_grs import _anchor_score_for_confidence
        self.assertEqual(_anchor_score_for_confidence("HIGH"), 0.60)
        self.assertEqual(_anchor_score_for_confidence("MEDIUM"), 0.40)
        self.assertEqual(_anchor_score_for_confidence("LOW"), 0.20)
        # Missing or unknown defaults to MEDIUM (legacy GFF compatibility)
        self.assertEqual(_anchor_score_for_confidence(""), 0.40)
        self.assertEqual(_anchor_score_for_confidence(None), 0.40)
        self.assertEqual(_anchor_score_for_confidence("weird"), 0.40)
        # Case-insensitive
        self.assertEqual(_anchor_score_for_confidence("medium"), 0.40)

    def test_anchor_cluster_inherits_score(self):
        # Regression: Colletes Reg1 case — a single MEDIUM hit with no
        # overlapping synteny cluster should produce an anchor at score 0.40,
        # NOT 1.00 (which would outrank genuine multi-gene synteny clusters).
        from cluster_grs import build_goi_anchor_clusters
        anchors = build_goi_anchor_clusters(
            goi_intervals=[
                {"chrom": "c1", "start": 1000, "end": 2000, "confidence": "MEDIUM"},
                {"chrom": "c1", "start": 5000, "end": 6000, "confidence": "HIGH"},
            ],
            existing_clusters=[],  # nothing overlaps → both anchors injected
        )
        self.assertEqual(len(anchors), 2)
        self.assertEqual(anchors[0]["score"], 0.40)
        self.assertTrue(anchors[0]["is_goi_anchor"])
        self.assertEqual(anchors[1]["score"], 0.60)

    def test_anchor_skipped_when_synteny_cluster_overlaps(self):
        # If an existing scored cluster already covers the GOI interval, no
        # anchor should be injected — the cluster's own score wins.
        from cluster_grs import build_goi_anchor_clusters
        anchors = build_goi_anchor_clusters(
            goi_intervals=[{"chrom": "c1", "start": 1000, "end": 2000, "confidence": "HIGH"}],
            existing_clusters=[{"chrom": "c1", "start": 1500, "end": 1700}],
        )
        self.assertEqual(anchors, [])


if __name__ == "__main__":
    unittest.main()
