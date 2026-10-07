#!/usr/bin/env python3
"""Regression tests for the goi_missing_but_strong_synteny region class (docs/TODO.md §1e).

Background: in the luciferase rerun, Aquatica's real Luc1 neighbourhood (block b7,
~8 HIGH-confidence flanking genes, GOI itself too divergent to model) was dropped from
the output BED — all 6 emitted regions were weak GOI-fallback score_floor picks that get
a +0.15 overlap bonus and crowd out the genuine flanking-only neighbourhood.

cluster_grs.py now:
  * counts HIGH-confidence flanking genes per cluster from the target GFF,
  * surfaces a block with >= strong_synteny_min_flanking such genes and no GOI as
    region_class=goi_missing_but_strong_synteny / selection_reason=flanking_only_no_goi_call,
    promoted to >= MEDIUM, guaranteed a slot even when the adaptive cap is full of weak
    GOI-fallback regions.
"""

import csv
import os
import subprocess
import sys
import tempfile
import unittest

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

import cluster_grs  # noqa: E402

CLUSTER_SCRIPT = os.path.abspath(os.path.join(BIN_DIR, "cluster_grs.py"))


def _flank_mrna(chrom, start, end, gene, confidence):
    attrs = (
        f"ID={gene}|{chrom}_{start}_flank;Name={gene};SynVoy_Parent={gene};"
        f"SynVoyRole=flanking;Confidence={confidence};EvidenceType=flanking_miniprot;"
        f"ModelStatus=complete;Identity=70.0"
    )
    return "\t".join([chrom, "flanking_annotation", "mRNA", str(start), str(end),
                      "70.0", "+", ".", attrs])


def _goi_mrna(chrom, start, end, identity, confidence, goi_class):
    attrs = (
        f"ID=GOI_LOC1|{chrom}_fallback;Name=GOI_LOC1;SynVoy_Parent=GOI_LOC1;"
        f"SynVoyRole=goi;Confidence={confidence};GOIClass={goi_class};"
        f"EvidenceType=fallback_hit_span;Identity={identity}"
    )
    return "\t".join([chrom, "fallback_hits", "mRNA", str(start), str(end),
                      str(identity), "+", ".", attrs])


class TestStrongSyntenyUnit(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def _gff(self, lines):
        path = os.path.join(self.dir, "target.gff")
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n" + "\n".join(lines) + "\n")
        return path

    def test_load_flanking_intervals_only_high(self):
        gff = self._gff([
            _flank_mrna("c1", 100, 200, "geneA", "HIGH"),
            _flank_mrna("c1", 300, 400, "geneB", "MEDIUM"),  # excluded
            _flank_mrna("c1", 500, 600, "geneC", "HIGH"),
            _goi_mrna("c1", 700, 800, 30.0, "LOW", "ambiguous_goi_family_member"),  # not flanking
        ])
        ivs = cluster_grs.load_flanking_intervals_from_gff(gff, confidences=("HIGH",))
        genes = sorted(iv["gene"] for iv in ivs)
        self.assertEqual(genes, ["geneA", "geneC"])

    def test_count_high_flanking_dedups_by_gene_and_respects_span(self):
        ivs = [
            {"chrom": "c1", "start": 100, "end": 200, "gene": "geneA"},
            {"chrom": "c1", "start": 100, "end": 250, "gene": "geneA"},  # same gene
            {"chrom": "c1", "start": 5000, "end": 5100, "gene": "geneB"},
            {"chrom": "c2", "start": 100, "end": 200, "gene": "geneC"},  # other chrom
        ]
        self.assertEqual(cluster_grs.count_high_flanking("c1", 0, 6000, ivs), 2)
        self.assertEqual(cluster_grs.count_high_flanking("c1", 0, 300, ivs), 1)
        self.assertEqual(cluster_grs.count_high_flanking("c2", 0, 300, ivs), 1)


class TestStrongSyntenyIntegration(unittest.TestCase):
    """Drive cluster_grs.py end-to-end via its CLI on synthetic inputs."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        # 6 weak GOI-fallback neighbourhoods (each its own chrom) + 1 strong flanking-only
        # block. The weak ones get the GOI overlap bonus and would fill the cap of 6.
        self.weak_chroms = [f"cW{i}" for i in range(1, 7)]
        self.strong_chrom = "cStrong"

        m8_lines = []
        gff_lines = []
        bed_lines = []
        rank = 0

        # Strong block: 6 HIGH flanking genes clustered within ~50 kb, NO GOI.
        for j in range(6):
            gene = f"sg{j}"
            pos = 200000 + j * 8000
            m8_lines.append(self._m8(gene, self.strong_chrom, pos, pos + 1500))
            gff_lines.append(_flank_mrna(self.strong_chrom, pos, pos + 1500, gene, "HIGH"))
            bed_lines.append(f"chrHome\t{rank*1000}\t{rank*1000+900}\t{gene}\t.\t+")
            rank += 1

        # Weak GOI clusters: 2 flanking hits + 1 MEDIUM probable_goi each.
        for ci, chrom in enumerate(self.weak_chroms):
            for j in range(2):
                gene = f"wg{ci}_{j}"
                pos = 100000 + j * 5000
                m8_lines.append(self._m8(gene, chrom, pos, pos + 1200))
                gff_lines.append(_flank_mrna(chrom, pos, pos + 1200, gene, "MEDIUM"))
                bed_lines.append(f"chrHome\t{rank*1000}\t{rank*1000+900}\t{gene}\t.\t+")
                rank += 1
            gff_lines.append(_goi_mrna(chrom, 100500, 101500, 24.0, "MEDIUM", "probable_goi"))

        self.hits = self._write("hits.m8", "\n".join(m8_lines) + "\n")
        self.gff = self._write("target.gff", "##gff-version 3\n" + "\n".join(gff_lines) + "\n")
        self.bed = self._write("synteny.bed", "\n".join(bed_lines) + "\n")

    @staticmethod
    def _m8(query, chrom, tstart, tend):
        # query target pid aln mm gap qs qe ts te evalue bits
        return "\t".join([query, chrom, "70.0", "400", "1", "0", "1", "400",
                          str(tstart), str(tend), "1e-30", "200"])

    def _write(self, name, content):
        path = os.path.join(self.dir, name)
        with open(path, "w") as fh:
            fh.write(content)
        return path

    def _run(self, extra_args=None):
        out_bed = os.path.join(self.dir, "regions.bed")
        out_scores = os.path.join(self.dir, "regions.scores.tsv")
        cmd = [
            sys.executable, CLUSTER_SCRIPT,
            "--hits", self.hits,
            "--target_gff", self.gff,
            "--synteny_bed", self.bed,
            "--output", out_bed,
            "--scores_output", out_scores,
            "--cluster_distance", "50000",
        ] + (extra_args or [])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, msg=f"stderr={proc.stderr}")
        with open(out_scores) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows, proc.stderr

    def test_strong_block_surfaced_despite_full_cap(self):
        rows, _ = self._run()
        by_chrom = {r["chrom"]: r for r in rows}
        # The flanking-only strong block must be present even though 6 weak GOI clusters
        # (with the overlap bonus) would otherwise fill the adaptive cap of 6.
        self.assertIn(self.strong_chrom, by_chrom,
                      msg=f"strong block dropped; emitted chroms={sorted(by_chrom)}")
        strong = by_chrom[self.strong_chrom]
        self.assertEqual(strong["region_class"], "goi_missing_but_strong_synteny")
        self.assertEqual(strong["selection_reason"], "flanking_only_no_goi_call")
        self.assertEqual(strong["goi_missing"], "True")
        self.assertIn(strong["confidence"], {"MEDIUM", "HIGH"})
        self.assertGreaterEqual(int(strong["high_flanking_count"]), 5)
        self.assertEqual(strong["goi_overlap"], "False")

    def test_feature_disabled_does_not_classify(self):
        rows, _ = self._run(extra_args=["--strong_synteny_min_flanking", "0"])
        classes = {r["region_class"] for r in rows}
        self.assertNotIn("goi_missing_but_strong_synteny", classes)

    def test_high_threshold_not_met_not_classified(self):
        # Require 99 HIGH flanking genes -> the 6-gene block no longer qualifies.
        rows, _ = self._run(extra_args=["--strong_synteny_min_flanking", "99"])
        classes = {r["region_class"] for r in rows}
        self.assertNotIn("goi_missing_but_strong_synteny", classes)


class TestStrongSyntenyUnboundRegression(unittest.TestCase):
    """§1q: `strong_synteny` was READ by the distant-macrosynteny gate before it was
    assigned, in the same per-cluster loop.

    Reproducing it needs the FIRST scored cluster to overlap a GOI — otherwise an
    earlier iteration has already bound the name and the bug degrades from a crash
    into a silent stale read. This fixture therefore has exactly ONE cluster, and it
    overlaps a GOI. Against the pre-fix ordering this raises
    ``UnboundLocalError: cannot access local variable 'strong_synteny'``; it is how a
    real run died (yeast STE2 demo, 2026-07-21).
    """

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        chrom = "cA"
        # one flanking anchor + a GOI model overlapping it => goi_overlap on cluster #1
        m8 = "\t".join([ "fg0", chrom, "70.0", "400", "1", "0", "1", "400",
                          "100000", "101200", "1e-30", "200"])
        gff = [_flank_mrna(chrom, 100000, 101200, "fg0", "HIGH"),
               _goi_mrna(chrom, 100400, 101000, 24.0, "MEDIUM", "probable_goi")]
        self.hits = self._write("hits.m8", m8 + "\n")
        self.gff = self._write("target.gff", "##gff-version 3\n" + "\n".join(gff) + "\n")
        self.bed = self._write("synteny.bed", "chrHome\t0\t900\tfg0\t.\t+\n")

    def _write(self, name, content):
        path = os.path.join(self.dir, name)
        with open(path, "w") as fh:
            fh.write(content)
        return path

    def test_single_goi_overlapping_cluster_does_not_crash(self):
        out_bed = os.path.join(self.dir, "regions.bed")
        out_scores = os.path.join(self.dir, "regions.scores.tsv")
        proc = subprocess.run(
            [sys.executable, CLUSTER_SCRIPT,
             "--hits", self.hits, "--target_gff", self.gff, "--synteny_bed", self.bed,
             "--output", out_bed, "--scores_output", out_scores,
             "--cluster_distance", "50000"],
            capture_output=True, text=True)
        self.assertNotIn("UnboundLocalError", proc.stderr)
        self.assertEqual(proc.returncode, 0, msg=f"stderr={proc.stderr[-2000:]}")
        with open(out_scores) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertTrue(any(r["goi_overlap"] == "True" for r in rows),
                        msg="fixture must produce a GOI-overlapping cluster")


if __name__ == "__main__":
    unittest.main()
