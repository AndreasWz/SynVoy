#!/usr/bin/env python3
"""Regression tests for the self_consistency report block (docs/TODO.md §1j).

Covers two flag types currently emitted:
  * strong_synteny_no_goi  — (genome, locus, chrom) with >= strong_flanking_min HIGH
    flanking genes but zero HIGH-confidence GOI calls. Advises the user to lower the
    classify thresholds or inspect the block manually (the orthologous neighbourhood
    is there but the GOI itself is too divergent for the current model).
  * cross_locus_duplicate  — same target gene identically modeled from >=2 home loci
    (already collapsed by §1h dedup; surfaced here for visibility).

The reciprocal-best paralog check (TODO §1j first bullet) is deferred; that's
documented in build_self_consistency.deferred_checks.
"""

import os
import shutil
import sys
import tempfile
import unittest

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

from generate_report import (  # noqa: E402
    build_report,
    build_self_consistency,
    collect_goi_annotations,
    collect_high_flanking_per_locus,
    dedupe_goi_annotations,
)


def _flank_mrna(chrom, start, end, gene, conf="HIGH"):
    attrs = (
        f"ID={gene}|{chrom}_{start};Name={gene};SynVoy_Parent={gene};"
        f"SynVoyRole=flanking;Confidence={conf};EvidenceType=flanking_miniprot;"
        f"ModelStatus=complete;Identity=70.0"
    )
    return "\t".join([chrom, "flanking_annotation", "mRNA", str(start), str(end),
                      "70.0", "+", ".", attrs])


def _goi_mrna(chrom, start, end, identity, conf="HIGH", goi_class="confident_goi",
              target_gene=""):
    attrs = (
        f"ID=GOI_LOC1|{chrom}_{start};Name=GOI_LOC1;SynVoy_Parent=GOI_LOC1;"
        f"SynVoyRole=goi;Confidence={conf};GOIClass={goi_class};"
        f"EvidenceType=exon_annotation;Identity={identity};ModelStatus=complete"
    )
    if target_gene:
        attrs += f";TargetGene={target_gene}"
    return "\t".join([chrom, "exon_annotation", "mRNA", str(start), str(end),
                      str(identity), "+", ".", attrs])


def _gff(lines):
    return "##gff-version 3\n" + "\n".join(lines) + "\n"


class TestStrongSyntenyNoGoiFlag(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "regions"))
        os.makedirs(os.path.join(self.dir, "scores"))

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _write(self, rel, content):
        with open(os.path.join(self.dir, rel), "w") as fh:
            fh.write(content)

    def test_flag_fires_when_strong_flanking_no_high_goi(self):
        # locus_1 / GCA_x has 25 HIGH flanking genes on chrA and NO HIGH GOI -> flag.
        flanks = [_flank_mrna("chrA", 1000 + i*5000, 1500 + i*5000, f"flA{i}")
                  for i in range(25)]
        self._write("regions/locus_1__GCA_111.1.fna.gff", _gff(flanks))

        report = build_report(self.dir)
        sc = report["self_consistency"]
        flags = [f for f in sc["flags"] if f["type"] == "strong_synteny_no_goi"]
        self.assertEqual(len(flags), 1)
        self.assertEqual(flags[0]["genome"], "GCA_111.1")
        self.assertEqual(flags[0]["locus"], "locus_1")
        self.assertEqual(flags[0]["chrom"], "chrA")
        self.assertEqual(flags[0]["high_flanking_count"], 25)
        self.assertIn("classify", flags[0]["advice"].lower())
        self.assertEqual(report["summary"]["strong_synteny_no_goi_flags"], 1)

    def test_no_flag_when_high_goi_present(self):
        # Same strong flanking, but a HIGH GOI call on the same chrom -> no flag.
        lines = [_goi_mrna("chrA", 50000, 51000, "85.0")]
        lines += [_flank_mrna("chrA", 1000 + i*5000, 1500 + i*5000, f"flA{i}")
                  for i in range(25)]
        self._write("regions/locus_1__GCA_222.1.fna.gff", _gff(lines))

        report = build_report(self.dir)
        flags = [f for f in report["self_consistency"]["flags"]
                 if f["type"] == "strong_synteny_no_goi"]
        self.assertEqual(flags, [])

    def test_no_flag_below_threshold(self):
        # Only 10 HIGH flanking, default threshold 20 -> no flag.
        flanks = [_flank_mrna("chrA", 1000 + i*5000, 1500 + i*5000, f"flA{i}")
                  for i in range(10)]
        self._write("regions/locus_1__GCA_333.1.fna.gff", _gff(flanks))

        report = build_report(self.dir)
        flags = [f for f in report["self_consistency"]["flags"]
                 if f["type"] == "strong_synteny_no_goi"]
        self.assertEqual(flags, [])

    def test_dedup_of_high_flank_by_gene_name(self):
        # Same gene appearing twice (e.g. annotation + hit-span) counts once.
        flanks = [_flank_mrna("chrA", 1000, 1500, "geneA"),
                  _flank_mrna("chrA", 1000, 2500, "geneA"),
                  _flank_mrna("chrA", 5000, 5500, "geneB")]
        self._write("regions/locus_1__GCA_444.1.fna.gff", _gff(flanks))
        counts = collect_high_flanking_per_locus(
            [os.path.join(self.dir, "regions/locus_1__GCA_444.1.fna.gff")]
        )
        self.assertEqual(counts[("GCA_444.1", "locus_1", "chrA")], 2)


class TestCrossLocusDuplicateSurfaced(unittest.TestCase):
    def test_dup_records_appear_in_self_consistency(self):
        # Two GFFs from different loci with the same target chrom/coords/identity
        # -> §1h dedup flags it cross_locus_duplicate; §1j surfaces it here.
        dir_ = tempfile.mkdtemp()
        try:
            os.makedirs(os.path.join(dir_, "regions"))
            for loc in ("locus_3", "locus_4"):
                with open(os.path.join(dir_, f"regions/{loc}__GCA_555.1.fna.gff"), "w") as fh:
                    fh.write(_gff([_goi_mrna("chrZ", 100, 200, "69.0",
                                             target_gene="RN001_000713")]))
            report = build_report(dir_)
            dup_flags = [f for f in report["self_consistency"]["flags"]
                         if f["type"] == "cross_locus_duplicate"]
            self.assertEqual(len(dup_flags), 1)
            self.assertEqual(dup_flags[0]["loci"], ["locus_3", "locus_4"])
            self.assertEqual(dup_flags[0]["identity"], "69.0")
            self.assertEqual(dup_flags[0]["chrom"], "chrZ")
        finally:
            shutil.rmtree(dir_, ignore_errors=True)


class TestUnit(unittest.TestCase):
    def test_build_self_consistency_threshold_arg(self):
        # Threshold pass-through: 10 flanking + threshold=5 fires; threshold=15 doesn't.
        flanking = {("G", "locus_1", "c1"): 10}
        annos = []
        dedup = {"records": []}
        out_lo = build_self_consistency(annos, dedup, flanking, strong_flanking_min=5)
        out_hi = build_self_consistency(annos, dedup, flanking, strong_flanking_min=15)
        self.assertEqual(out_lo["summary"]["n_strong_synteny_no_goi"], 1)
        self.assertEqual(out_hi["summary"]["n_strong_synteny_no_goi"], 0)
        # Reciprocal-best is explicitly listed as deferred (transparency for users).
        self.assertIn("reciprocal_best_paralog", out_lo["deferred_checks"])


if __name__ == "__main__":
    unittest.main()
