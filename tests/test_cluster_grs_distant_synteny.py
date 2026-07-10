#!/usr/bin/env python3
"""Tests for the distant-macrosynteny GOI rescue in cluster_grs.py.

Motivation (oskar divergence-ladder run 5702852): in Aedes aegypti the oskar
query hit the real ortholog at NC_035109.1:184.15M with 45.3% identity /
e-value 2.77e-29 — but the neighbourhood had dispersed (the nearest flanking
anchors sit 62.8 Mb away on the same chromosome), so the region scored ~0.15
for lack of LOCAL flanking and was dropped. The rescue promotes such a strong,
same-chromosome, dispersed hit in a rearranged genome, while still rejecting
paralogs (biglycan, on a different chromosome from decorin's flanking) and
fold-analogues (the structdisc lipase, low identity).
"""

import os
import subprocess
import sys
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

import cluster_grs as cg  # noqa: E402


class TestRearrangementScore(unittest.TestCase):
    def _hit(self, query, chrom, start, ev=1e-20):
        return {"query": query, "chrom": chrom, "start": start,
                "end": start + 300, "evalue": ev, "identity": 80.0}

    def test_collinear_scores_zero(self):
        # flanking in home order, one scaffold, tightly spaced → conserved.
        gene_map = {f"f{i}": i for i in range(5)}
        hits = [self._hit(f"f{i}", "chr1", 1000 + i * 5000) for i in range(5)]
        self.assertLess(cg.target_rearrangement_score(hits, gene_map), 0.1)

    def test_shattered_scores_high(self):
        # Aedes topology: home-consecutive flanking land on 5 scaffolds, adjacency
        # destroyed; only two stay together on one chromosome (far from the rest).
        gene_map = {"f0": 0, "f1": 1, "f2": 2, "f3": 3, "f4": 4, "f5": 5}
        hits = [
            self._hit("f0", "s5", 100),
            self._hit("f1", "chr1", 121_000_000),
            self._hit("f2", "chr1", 121_050_000),   # adjacent to f1
            self._hit("f3", "s7", 100),
            self._hit("f4", "s8", 100),
            self._hit("f5", "s9", 100),
        ]
        self.assertGreater(cg.target_rearrangement_score(hits, gene_map), 0.5)

    def test_too_few_flanking_returns_zero(self):
        gene_map = {"f0": 0, "f1": 1}
        self.assertEqual(
            cg.target_rearrangement_score([self._hit("f0", "chr1", 1)], gene_map), 0.0)

    def test_ignores_goi_and_unknown_queries(self):
        gene_map = {"f0": 0, "f1": 1}
        hits = [self._hit("GOI_osk", "chr1", 1), self._hit("stranger", "chr9", 1)]
        self.assertEqual(cg.target_rearrangement_score(hits, gene_map), 0.0)


class TestRescueDecision(unittest.TestCase):
    """The exact cases the user asked to separate."""

    FLANK = {"NC_035109.1", "NC_035107.1", "NC_035108.1"}   # Aedes anchor scaffolds

    def test_aedes_oskar_rescued(self):
        self.assertTrue(cg.qualifies_distant_macrosynteny_rescue(
            45.3, 2.77e-29, "NC_035109.1", self.FLANK, 0.73))

    def test_structdisc_lipase_rejected_low_identity(self):
        # fold analogue: weak identity AND on a scaffold with no anchor
        self.assertFalse(cg.qualifies_distant_macrosynteny_rescue(
            28.6, 1e-3, "STHB01000011.1", self.FLANK, 0.73))

    def test_biglycan_rejected_different_chromosome(self):
        # decorin paralog: strong identity, but chrX carries no decorin flanking
        # (they are on chr5) → no shared-chromosome macrosynteny.
        self.assertFalse(cg.qualifies_distant_macrosynteny_rescue(
            55.0, 1e-20, "chrX", {"chr5", "chr9"}, 0.6))

    def test_collinear_genome_rejected(self):
        # strong isolated hit but the genome is NOT rearranged → likely a paralog,
        # stay strict.
        self.assertFalse(cg.qualifies_distant_macrosynteny_rescue(
            50.0, 1e-15, "chr1", {"chr1"}, 0.1))

    def test_borderline_identity_rejected(self):
        self.assertFalse(cg.qualifies_distant_macrosynteny_rescue(
            39.9, 1e-29, "NC_035109.1", self.FLANK, 0.73))


class TestEndToEndRescue(unittest.TestCase):
    """Drive cluster_grs.py on an Aedes-topology m8 and confirm the dispersed
    oskar region is emitted as a MEDIUM goi_dispersed_macrosynteny call."""

    def _run(self, tmp, extra=None):
        # home order: f0 f1 GOI_osk f2 f3 f4 f5
        bed = os.path.join(tmp, "syn.bed")
        with open(bed, "w") as fh:
            for i, name in enumerate(["f0", "f1", "GOI_osk", "f2", "f3", "f4", "f5"]):
                fh.write(f"home\t{i*1000}\t{i*1000+500}\t{name}\n")
        # m8: GOI hits its real (dispersed) locus on chr1 @184M at 45%/e-29;
        # flanking f1,f2 land on chr1 @121M (62.8 Mb away), rest shattered.
        m8 = os.path.join(tmp, "hits.m8")
        rows = [
            ("GOI_osk", "chr1", 45.3, 184_000_000, 184_000_459, 2.77e-29),
            ("f1", "chr1", 90.0, 121_000_000, 121_001_000, 1e-40),
            ("f2", "chr1", 84.0, 121_050_000, 121_060_000, 1e-40),
            ("f0", "s5", 40.0, 100, 400, 1e-9),
            ("f3", "s7", 40.0, 100, 400, 1e-9),
            ("f4", "s8", 40.0, 100, 400, 1e-9),
            ("f5", "s9", 40.0, 100, 400, 1e-9),
        ]
        with open(m8, "w") as fh:
            for q, t, pid, ts, te, ev in rows:
                fh.write(f"{q}\t{t}\t{pid}\t100\t0\t0\t1\t100\t{ts}\t{te}\t{ev}\t120\n")
        out = os.path.join(tmp, "regions.bed")
        scores = os.path.join(tmp, "scores.tsv")
        cmd = [sys.executable, os.path.join(REPO_ROOT, "bin", "cluster_grs.py"),
               "--hits", m8, "--synteny_bed", bed,
               "--output", out, "--scores_output", scores,
               "--flanking_count", "6"] + (extra or [])
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        with open(scores) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows_out = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh]
        return rows_out

    def test_dispersed_oskar_is_rescued(self):
        import tempfile
        rows = self._run(tempfile.mkdtemp())
        disp = [r for r in rows if r["region_class"] == "goi_dispersed_macrosynteny"]
        self.assertEqual(len(disp), 1, f"expected one rescued region, got {rows}")
        r = disp[0]
        self.assertEqual(r["chrom"], "chr1")
        self.assertEqual(r["confidence"], "MEDIUM")
        self.assertEqual(r["selection_reason"], "distant_macrosynteny_rescue")

    def test_kill_switch_disables_rescue(self):
        import tempfile
        rows = self._run(tempfile.mkdtemp(), extra=["--disable_distant_synteny_rescue"])
        disp = [r for r in rows if r["region_class"] == "goi_dispersed_macrosynteny"]
        self.assertEqual(disp, [])


if __name__ == "__main__":
    unittest.main()
