#!/usr/bin/env python3
"""Regression tests for cross-locus GOI dedup in bin/generate_report.py (docs/TODO.md §1h).

Background: a target ortholog is found once per home locus whose flanking neighbourhood
overlaps it. After the §1d locus cap the luciferase rerun still reported the same Aquatica
gene identically from locus_3 AND locus_4. The report headline should count distinct
ortholog genes, not per-locus seeds, and flag bit-for-bit-identical multi-locus hits as
``cross_locus_duplicate``.

These fixtures mirror the real luciferase_rerun_auto pattern with tiny synthetic GFFs:
  * Pyrocoelia-like genome  : two genuine HIGH genes on different chroms (83.5 %, 56.5 %)
  * Aquatica-like genome    : one HIGH gene reported identically from locus_3 and locus_4
  * Luciola-like genome     : one HIGH gene from a single locus
Expected post-dedup: 4 distinct HIGH genes (was 5), 1 cross_locus_duplicate.
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
    canonical_genome_id,
    collect_goi_annotations,
    dedupe_goi_annotations,
    staging_source_label,
)


def _goi_mrna(chrom, start, end, identity, confidence="HIGH",
              goi_class="confident_goi", target_gene="", parent="GOI_LOC1"):
    attrs = (
        f"ID={parent}|blk_{chrom}_{start};Name={parent};SynVoy_Parent={parent};"
        f"SynVoyRole=goi;Confidence={confidence};GOIClass={goi_class};"
        f"EvidenceType=exon_annotation;Identity={identity};ModelStatus=complete"
    )
    if target_gene:
        attrs += f";TargetGene={target_gene}"
    return "\t".join([chrom, "exon_annotation", "mRNA", str(start), str(end),
                      str(identity), "+", ".", attrs])


def _gff(*mrna_lines):
    return "##gff-version 3\n" + "\n".join(mrna_lines) + "\n"


class TestStagingPrefixParsing(unittest.TestCase):
    def test_canonical_genome_id_strips_locus_prefix(self):
        # GCA names already survive via the embedded-accession regex...
        self.assertEqual(canonical_genome_id("locus_3__GCA_036250285.1.fna.gff"),
                         "GCA_036250285.1")
        # ...but for Pro-mode names the invariant that matters is that the prefixed and
        # non-prefixed forms collapse to the SAME id, so one genome's annotations don't
        # get split across loci. (Suffix .fa/.fasta is left intact either way.)
        self.assertEqual(canonical_genome_id("locus_3__Colletes_gigas.fa.gff"),
                         canonical_genome_id("Colletes_gigas.fa.gff"))
        self.assertEqual(canonical_genome_id("locus_4__Colletes_gigas.fa.gff"),
                         canonical_genome_id("locus_3__Colletes_gigas.fa.gff"))
        # Numeric staging fallback (modules without locus threading) is also stripped.
        self.assertEqual(canonical_genome_id("7__home.fasta.gff"),
                         canonical_genome_id("home.fasta.gff"))

    def test_canonical_genome_id_collapses_rescue_passes(self):
        # §1m/§17 hull rescue and §1e strong-synteny rescue emit
        # "<genome>.hull_rescue.gff" / "<genome>.rescue.gff". These MUST collapse to
        # the real genome so the rescue model is deduped + owned alongside the main
        # calls — otherwise a redundant paralog rescue leaks as a phantom genome and
        # escapes ownership demotion (the dcn_18 mouse.hull_rescue biglycan leak).
        self.assertEqual(canonical_genome_id("mouse.hull_rescue.gff"), "mouse")
        self.assertEqual(canonical_genome_id("mouse.hull_rescue.gff"),
                         canonical_genome_id("mouse.fna.gff"))
        self.assertEqual(canonical_genome_id("zebrafish.fna.rescue.gff"), "zebrafish")
        self.assertEqual(canonical_genome_id("cow.hull_rescue.gff"),
                         canonical_genome_id("cow.fna"))

    def test_canonical_genome_id_strips_fa_fasta(self):
        # docs/TODO_JUN.md QW1: Pro-mode ".fa"/".fasta" genomes must collapse like
        # ".fna"/".faa" — otherwise "Colletes_gigas.fa.gff" -> "Colletes_gigas.fa"
        # while the bare-stem hull-rescue "Colletes_gigas.hull_rescue.gff" ->
        # "Colletes_gigas" split one genome into a real record + an empty phantom
        # that falsely populated summary.goi_absent_genomes (melittin_full: all 5
        # solitary/stingless bees read as "melittin absent").
        self.assertEqual(canonical_genome_id("Colletes_gigas.fa.gff"), "Colletes_gigas")
        self.assertEqual(canonical_genome_id("Xylocopa_violacea.fasta.gff"), "Xylocopa_violacea")
        # The bug: the hull-rescue (bare stem) and the main .fa GFF MUST agree now.
        self.assertEqual(canonical_genome_id("Colletes_gigas.hull_rescue.gff"),
                         canonical_genome_id("Colletes_gigas.fa.gff"))
        # ".fa" must NOT over-strip ".faa"/".fasta"/".fna" (longer suffixes win).
        self.assertEqual(canonical_genome_id("Apis_mellifera.faa"), "Apis_mellifera")
        self.assertEqual(canonical_genome_id("GCF_123.1.fa"), "GCF_123.1")

    def test_staging_source_label(self):
        self.assertEqual(staging_source_label("locus_12__GCA_1.1.fna.gff"), "locus_12")
        self.assertEqual(staging_source_label("3__GCA_1.1.fna.gff"), "3")
        # No prefix -> the bare filename is the source label (always non-empty).
        self.assertEqual(staging_source_label("GCA_1.1.fna.gff"), "GCA_1.1.fna.gff")


class TestDedupeUnit(unittest.TestCase):
    def test_identical_across_loci_is_cross_locus_duplicate(self):
        annos = [
            {"genome": "G", "source": "locus_3", "chrom": "c1", "start": 500, "end": 700,
             "confidence": "HIGH", "goi_class": "confident_goi", "identity": "69.0",
             "target_gene": "RN001_000713", "model_status": "complete"},
            {"genome": "G", "source": "locus_4", "chrom": "c1", "start": 500, "end": 700,
             "confidence": "HIGH", "goi_class": "confident_goi", "identity": "69.0",
             "target_gene": "RN001_000713", "model_status": "complete"},
        ]
        out = dedupe_goi_annotations(annos)
        self.assertEqual(out["post_dedup_records"], 1)
        self.assertEqual(out["high_confidence_goi_deduped"], 1)
        self.assertEqual(out["cross_locus_duplicates"], 1)
        self.assertEqual(out["hits_collapsed_by_dedup"], 1)
        rec = out["records"][0]
        self.assertEqual(rec["goi_class"], "cross_locus_duplicate")
        self.assertEqual(rec["provenance"], ["locus_3", "locus_4"])
        self.assertEqual(rec["n_merged_hits"], 2)

    def test_overlap_but_different_identity_merges_without_dup_class(self):
        # Same locus neighbourhood, different model identities -> still one gene,
        # but NOT the bit-for-bit "one gene, many seeds" signal.
        annos = [
            {"genome": "G", "source": "locus_3", "chrom": "c1", "start": 500, "end": 700,
             "confidence": "MEDIUM", "goi_class": "probable_goi", "identity": "48.0",
             "target_gene": "", "model_status": "partial"},
            {"genome": "G", "source": "locus_4", "chrom": "c1", "start": 505, "end": 705,
             "confidence": "HIGH", "goi_class": "confident_goi", "identity": "69.0",
             "target_gene": "", "model_status": "complete"},
        ]
        out = dedupe_goi_annotations(annos)
        self.assertEqual(out["post_dedup_records"], 1)
        rec = out["records"][0]
        # Representative is the HIGH/higher-identity member.
        self.assertEqual(rec["confidence"], "HIGH")
        self.assertEqual(rec["identity"], "69.0")
        self.assertFalse(rec["cross_locus_duplicate"])
        self.assertEqual(rec["goi_class"], "confident_goi")
        self.assertEqual(rec["provenance"], ["locus_3", "locus_4"])

    def test_non_overlapping_same_chrom_stay_distinct(self):
        annos = [
            {"genome": "G", "source": "locus_1", "chrom": "c1", "start": 100, "end": 200,
             "confidence": "HIGH", "goi_class": "confident_goi", "identity": "83.5",
             "target_gene": "", "model_status": "complete"},
            {"genome": "G", "source": "locus_1", "chrom": "c1", "start": 9000, "end": 9100,
             "confidence": "HIGH", "goi_class": "confident_goi", "identity": "55.0",
             "target_gene": "", "model_status": "complete"},
        ]
        out = dedupe_goi_annotations(annos)
        self.assertEqual(out["post_dedup_records"], 2)
        self.assertEqual(out["cross_locus_duplicates"], 0)

    def test_low_confidence_is_ignored(self):
        annos = [
            {"genome": "G", "source": "locus_1", "chrom": "c1", "start": 100, "end": 200,
             "confidence": "LOW", "goi_class": "ambiguous_goi_family_member",
             "identity": "22.0", "target_gene": "", "model_status": "fragment"},
        ]
        out = dedupe_goi_annotations(annos)
        self.assertEqual(out["post_dedup_records"], 0)
        self.assertEqual(out["pre_dedup_high_medium"], 0)


class TestDedupeIntegration(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.test_dir, "regions"))
        os.makedirs(os.path.join(self.test_dir, "scores"))

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def _write(self, rel, content):
        with open(os.path.join(self.test_dir, rel), "w") as fh:
            fh.write(content)

    def test_luciferase_pattern_collapses_cross_locus_duplicates(self):
        pyro = "GCA_036250285.1.fna.gff"
        aqua = "GCA_035610365.1.fna.gff"
        luci = "GCA_054859625.1.fna.gff"

        # Pyrocoelia: two genuine HIGH genes on DIFFERENT chroms (must stay distinct).
        self._write(f"regions/locus_1__{pyro}",
                    _gff(_goi_mrna("chrP1", 63185487, 63187465, "83.5", target_gene="RI129_003463")))
        self._write(f"regions/locus_4__{pyro}",
                    _gff(_goi_mrna("chrP2", 11800794, 11802726, "56.5", target_gene="RI129_009591")))

        # Aquatica: SAME gene reported identically from locus_3 and locus_4.
        for loc in ("locus_3", "locus_4"):
            self._write(f"regions/{loc}__{aqua}",
                        _gff(_goi_mrna("chrA1", 49613237, 49615185, "69.0", target_gene="RN001_000713")))

        # Luciola: one HIGH gene from a single locus.
        self._write(f"regions/locus_4__{luci}",
                    _gff(_goi_mrna("chrL1", 72775747, 72777668, "59.5")))

        report = build_report(self.test_dir)
        summary = report["summary"]

        # 5 raw HIGH seeds -> 4 distinct genes.
        self.assertEqual(summary["high_confidence_goi_pre_dedup"], 5)
        self.assertEqual(summary["high_confidence_goi"], 4)
        self.assertEqual(summary["headline_metric"], 4)
        self.assertEqual(summary["cross_locus_duplicate_goi"], 1)
        self.assertEqual(summary["goi_hits_collapsed_by_dedup"], 1)
        self.assertIn("4 high-confidence", summary["headline"])

        records = report["goi_dedup"]["records"]
        high = [r for r in records if r["confidence"] == "HIGH"]
        self.assertEqual(len(high), 4)

        # The Aquatica gene is the single cross-locus duplicate, with both loci as provenance.
        dups = [r for r in high if r["cross_locus_duplicate"]]
        self.assertEqual(len(dups), 1)
        self.assertEqual(dups[0]["genome"], "GCA_035610365.1")
        self.assertEqual(dups[0]["provenance"], ["locus_3", "locus_4"])
        self.assertEqual(dups[0]["goi_class"], "cross_locus_duplicate")

        # Pyrocoelia's two genes remain distinct (different chroms, not collapsed).
        pyro_high = sorted(r["identity"] for r in high if r["genome"] == "GCA_036250285.1")
        self.assertEqual(pyro_high, ["56.5", "83.5"])


if __name__ == "__main__":
    unittest.main()
