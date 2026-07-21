#!/usr/bin/env python3
"""Tests for bin/discover_regions.py — structural ORF discovery helpers.

Pure logic (region clustering, coord mapping, overlap filtering, GFF
materialisation); the ORF predictor is mocked so these run without Augustus.
"""

import os
import sys
import unittest
from unittest import mock

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

import discover_regions as dr  # noqa: E402


def _mrna(chrom, s, e, attrs):
    return "\t".join([chrom, "SynVoy", "mRNA", str(s), str(e), ".", "+", ".", attrs])


class TestReadRegionFeatures(unittest.TestCase):
    def test_splits_goi_from_flanking(self):
        import tempfile
        d = tempfile.mkdtemp()
        gff = os.path.join(d, "g.gff")
        with open(gff, "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write(_mrna("chr1", 100, 200, "ID=GOI_x;SynVoyRole=goi;GOIClass=probable_goi") + "\n")
            fh.write(_mrna("chr1", 500, 600, "ID=flank_1;SynVoyRole=flanking") + "\n")
            fh.write(_mrna("chr1", 800, 900, "ID=flank_2") + "\n")  # no role => flanking
        flanking, existing = dr.read_region_features(gff)
        # flanking excludes the GOI; existing has all three
        self.assertEqual(len(flanking), 2)
        self.assertEqual(len(existing), 3)
        # 0-based half-open
        self.assertIn(("chr1", 99, 200), existing)

    def test_missing_file(self):
        self.assertEqual(dr.read_region_features("/no/such"), ([], []))


class TestClusterIntervals(unittest.TestCase):
    def test_merges_within_gap_splits_beyond(self):
        iv = [("chr1", 1000, 1200), ("chr1", 1300, 1500), ("chr1", 900000, 900100)]
        w = dr.cluster_intervals(iv, max_gap=100000, pad=0)
        self.assertEqual(w, [("chr1", 1000, 1500), ("chr1", 900000, 900100)])

    def test_pad_and_clamp(self):
        w = dr.cluster_intervals([("chr1", 1000, 1200)], max_gap=1000, pad=20000,
                                 chrom_lengths={"chr1": 5000})
        self.assertEqual(w, [("chr1", 0, 5000)])   # clamped low + high

    def test_per_chromosome(self):
        iv = [("chr1", 100, 200), ("chr2", 100, 200)]
        w = dr.cluster_intervals(iv, max_gap=10, pad=0)
        self.assertEqual(sorted(w), [("chr1", 100, 200), ("chr2", 100, 200)])


class TestOverlap(unittest.TestCase):
    def test_overlap_frac_and_chrom_guard(self):
        orf = {"chrom": "chr1", "start": 1050, "end": 1180}
        self.assertTrue(dr.orf_overlaps_any(orf, [("chr1", 1000, 1200)], 0.3))
        self.assertFalse(dr.orf_overlaps_any(orf, [("chr2", 1000, 1200)], 0.3))
        # tiny overlap below threshold
        self.assertFalse(dr.orf_overlaps_any(
            {"chrom": "chr1", "start": 1190, "end": 2000}, [("chr1", 1000, 1200)], 0.3))


class TestConfidence(unittest.TestCase):
    def test_tm_and_embedding_tiers(self):
        self.assertEqual(dr.confidence_from_scores(0.9, None), "HIGH")
        self.assertEqual(dr.confidence_from_scores(0.55, None), "MEDIUM")
        self.assertEqual(dr.confidence_from_scores(0.2, None), "LOW")
        self.assertEqual(dr.confidence_from_scores(0.2, 0.9), "HIGH")   # embedding rescues
        self.assertEqual(dr.confidence_from_scores(None, 0.72), "MEDIUM")


class TestBuildDiscoveryGff(unittest.TestCase):
    def test_coords_1based_and_attrs(self):
        orf = {"chrom": "chr5", "start": 21015209, "end": 21017000,
               "strand": "+", "seq": "MKAYIL"}
        gff, faa = dr.build_discovery_gff(orf, "cow.fna", 0, 0.90, 0.5, "HIGH")
        mrna = gff[0].split("\t")
        self.assertEqual(mrna[2], "mRNA")
        self.assertEqual(mrna[3], "21015210")   # 0-based 21015209 -> GFF 21015210
        self.assertEqual(mrna[4], "21017000")
        self.assertIn("EvidenceType=structural_discovery", mrna[8])
        self.assertIn("Confidence=HIGH", mrna[8])
        self.assertIn("StructuralSimilarity=0.900", mrna[8])
        self.assertIn("SynVoyRole=goi", mrna[8])
        self.assertEqual(faa, ("cow.fna_structdisc_0", "MKAYIL"))
        # exon child references the mRNA
        self.assertIn("Parent=cow.fna_structdisc_0", gff[1])

    def test_no_tm_omits_structural_attr(self):
        orf = {"chrom": "c", "start": 0, "end": 30, "strand": "-", "seq": "MK"}
        gff, _ = dr.build_discovery_gff(orf, "g", 3, None, 0.8, "MEDIUM")
        self.assertNotIn("StructuralSimilarity=", gff[0])
        self.assertIn("EmbeddingSimilarity=0.800", gff[0])


class TestPredictRegionOrfsCoordMapping(unittest.TestCase):
    def test_rel_coords_mapped_to_chromosome(self):
        genome = {"chr1": "A" * 100000}
        window = ("chr1", 40000, 45000)   # slice offset 40000
        fake_orfs = [
            {"id": "g1", "seq": "M" * 100, "start": 100, "end": 400, "strand": "+"},
            {"id": "g2", "seq": "M" * 50, "start": 1000, "end": 1150, "strand": "-"},
        ]
        with mock.patch("gene_predictor.predict_orfs", return_value=fake_orfs):
            import tempfile
            out = dr.predict_region_orfs(genome, window, tempfile.mkdtemp())
        self.assertEqual(out[0]["chrom"], "chr1")
        self.assertEqual(out[0]["start"], 40100)   # 40000 + 100
        self.assertEqual(out[0]["end"], 40400)
        self.assertEqual(out[1]["start"], 41000)
        self.assertEqual(out[1]["strand"], "-")

    def test_missing_chrom_returns_empty(self):
        self.assertEqual(dr.predict_region_orfs({}, ("chrX", 0, 100), "/tmp/x"), [])


try:
    import parasail  # noqa: F401
    _HAVE_PARASAIL = True
except Exception:
    _HAVE_PARASAIL = False


# A deterministic ~120 aa "GOI" and derived ORFs for the sequence guards.
_AA = "ACDEFGHIKLMNPQRSTVWY"
_GOI = "".join(_AA[(i * 7 + 3) % 20] for i in range(120))
_UNREL = "".join(_AA[(i * 11 + 5) % 20] for i in range(120))  # different phase


class TestLengthRatioGuard(unittest.TestCase):
    def test_band(self):
        self.assertTrue(dr.passes_length_ratio(606, 606))
        self.assertTrue(dr.passes_length_ratio(900, 606))     # 1.49x ok
        self.assertFalse(dr.passes_length_ratio(2187, 606))   # 3.6x -> reject (the FP)
        self.assertFalse(dr.passes_length_ratio(200, 606))    # 0.33x -> reject

    def test_zero_goi_len_is_permissive(self):
        self.assertTrue(dr.passes_length_ratio(100, 0))


class TestSeqGoiChroms(unittest.TestCase):
    def test_confidence_gate(self):
        import tempfile
        d = tempfile.mkdtemp()
        gff = os.path.join(d, "g.gff")
        with open(gff, "w") as fh:
            fh.write(_mrna("chr1", 100, 200,
                           "ID=GOI_a;SynVoyRole=goi;Confidence=HIGH") + "\n")
            fh.write(_mrna("chr2", 100, 200,
                           "ID=GOI_b;SynVoyRole=goi;Confidence=LOW") + "\n")
            fh.write(_mrna("chr3", 100, 200, "ID=flank;SynVoyRole=flanking") + "\n")
        chroms = dr.seq_goi_chroms(gff, min_confidence="MEDIUM")
        self.assertEqual(chroms, {"chr1"})   # HIGH kept, LOW GOI + flanking excluded


class TestFlankingSupport(unittest.TestCase):
    def test_counts_same_chrom_within_window(self):
        orf = {"chrom": "chr1", "start": 1_000_000, "end": 1_001_000}
        flank = [("chr1", 950_000, 960_000), ("chr1", 1_050_000, 1_060_000),
                 ("chr1", 5_000_000, 5_001_000), ("chr2", 1_000_000, 1_001_000)]
        # window_bp=100k -> first two count, far one + other chrom don't
        self.assertEqual(dr.count_flanking_support(orf, flank, window_bp=100_000), 2)


@unittest.skipUnless(_HAVE_PARASAIL, "parasail required for sequence guards")
class TestGoiAlignment(unittest.TestCase):
    def test_full_length_high_coverage(self):
        cov, ident, gid, score = dr.best_goi_alignment(_GOI, [("GOI_x", _GOI)])
        self.assertGreater(cov, 0.95)
        self.assertGreater(ident, 95.0)
        self.assertEqual(gid, "GOI_x")

    def test_shared_domain_low_coverage(self):
        # ORF shares only a 15-aa patch of the GOI, rest unrelated -> low coverage.
        orf = _UNREL[:60] + _GOI[40:55] + _UNREL[60:]
        cov, ident, gid, score = dr.best_goi_alignment(orf, [("GOI_x", _GOI)])
        self.assertLess(cov, 0.4)


@unittest.skipUnless(_HAVE_PARASAIL, "parasail required for sequence guards")
class TestValidateDiscoveryHit(unittest.TestCase):
    def _orf(self, seq, chrom="chr1", start=1_000_000):
        return {"chrom": chrom, "start": start, "end": start + len(seq) * 3,
                "strand": "+", "seq": seq}

    def _gois(self):
        return [("GOI_x", _GOI)]

    def test_passes_clean_ortholog(self):
        ok, reason, m = dr.validate_discovery_hit(
            self._orf(_GOI), self._gois(), flanking_support=5)
        self.assertTrue(ok, reason)
        self.assertEqual(reason, "ok")
        self.assertGreater(m["goi_coverage"], 0.9)

    def test_guard4a_weak_synteny(self):
        ok, reason, _m = dr.validate_discovery_hit(
            self._orf(_GOI), self._gois(), flanking_support=1, min_block_flanking=2)
        self.assertFalse(ok)
        self.assertEqual(reason, "weak_synteny")

    def test_guard4b_off_primary_block(self):
        ok, reason, _m = dr.validate_discovery_hit(
            self._orf(_GOI, chrom="chrOTHER"), self._gois(), flanking_support=5,
            genome_has_seq_goi=True, seq_goi_chroms={"chrMAIN"})
        self.assertFalse(ok)
        self.assertEqual(reason, "off_primary_block")

    def test_guard1_length_ratio(self):
        # 3x the GOI length, in-frame concat so coverage of the GOI is fine —
        # only the length guard should catch it.
        long_orf = self._orf(_GOI * 3)
        ok, reason, _m = dr.validate_discovery_hit(
            long_orf, self._gois(), flanking_support=5)
        self.assertFalse(ok)
        self.assertEqual(reason, "length_ratio")

    def test_guard2_low_coverage_reproduces_the_FP(self):
        # A protein whose only similarity to the GOI is one short patch (the
        # promiscuous-fold false positive at the sequence level), kept in the
        # length band so guard 1 passes and guard 2 is what fires.
        fp = _UNREL[:100] + _GOI[30:45] + _UNREL[:5]   # 120 aa, one 15-aa GOI patch
        ok, reason, m = dr.validate_discovery_hit(
            self._orf(fp), self._gois(), flanking_support=5, min_goi_coverage=0.5)
        self.assertFalse(ok)
        self.assertEqual(reason, "low_goi_coverage")
        self.assertLess(m["goi_coverage"], 0.5)

    def test_guard3_beaten_by_paralog(self):
        # ORF is actually a paralog: aligns better to the panel than to the GOI.
        orf = self._orf(_UNREL)   # identical to the panel entry, unrelated to GOI
        ok, reason, _m = dr.validate_discovery_hit(
            orf, self._gois(), panel=[("home_paralog", _UNREL)],
            flanking_support=5, min_goi_coverage=0.0, paralog_margin=1.0)
        self.assertFalse(ok)
        self.assertEqual(reason, "beaten_by_paralog")


if __name__ == "__main__":
    unittest.main()
