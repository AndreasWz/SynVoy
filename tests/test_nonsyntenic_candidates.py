#!/usr/bin/env python3
"""§1x — the silent discard.

The most damaging failure mode SynVoy had: the search finds the right gene, the
synteny gate correctly declines to *call* it, and the pipeline then reports the
genome as having no ortholog at all. Refusing the hit is defensible; refusing it
silently is not — "nothing found" and "found, could not place it" are opposite
scientific claims.

Motivating case, traced in full (results/demo_Def, 2026-07-21): D. melanogaster
Defensin vs Anopheles gambiae (~250 My). The raw m8 holds

    GOI_Def   id=54.2  aln=192  e=5.845e-16  NC_064602.1:78,296,974-78,297,165

which lands inside RefSeq LOC1270637 (product=defensin). Every quality gate passes;
both synteny gates reject, because the flanking neighbourhood is recovered on a
DIFFERENT chromosome (16 anchors on NC_064601.1) and Anopheles has translocated the
gene off its ancestral neighbourhood. SynVoy emitted no model, no report entry and
no warning.

See docs/STATE_OF_THE_PROJECT.md Part C §1x.
"""
import importlib.util
import json
import os
import shutil
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


isr = _load("iterative_search_runner", "iterative_search_runner.py")
gr = _load("generate_report", "generate_report.py")


def _hit(query, chrom, start, end, pident=54.2, evalue=5.845e-16, alnlen=192):
    return {"query": query, "chrom": chrom, "start": start, "end": end,
            "pident": pident, "evalue": evalue, "alnlen": alnlen, "strand": "+"}


def _anopheles_case():
    """16 flanking anchors on NC_064601.1; the defensin hit on NC_064602.1."""
    flanking = [_hit(f"gene-FL{i}", "NC_064601.1",
                     13_300_000 + i * 3_000_000, 13_300_000 + i * 3_000_000 + 2000,
                     pident=70.0)
                for i in range(16)]
    # two stray anchors on the defensin's chromosome — below min_chrom_anchors=3
    flanking += [_hit("gene-STRAY1", "NC_064602.1", 10_000_000, 10_002_000, pident=65.0),
                 _hit("gene-STRAY2", "NC_064602.1", 12_000_000, 12_002_000, pident=64.0)]
    goi = _hit("GOI_Def", "NC_064602.1", 78_296_974, 78_297_165)
    return flanking + [goi], flanking


class TestRejectionIsRecorded(unittest.TestCase):
    def test_strong_hit_is_still_refused(self):
        """The gate's behaviour must not change — this is a reporting fix."""
        hits, flanking = _anopheles_case()
        seeds = isr.select_dispersed_goi_seeds(hits, flanking, min_chrom_anchors=3)
        self.assertEqual(seeds, [], "the synteny gate must still reject this hit")

    def test_rejection_is_recorded_with_its_gate(self):
        hits, flanking = _anopheles_case()
        rejected = []
        isr.select_dispersed_goi_seeds(hits, flanking, min_chrom_anchors=3,
                                       rejected_out=rejected,
                                       genome_name="Anopheles_gambiae.fna")
        self.assertEqual(len(rejected), 1)
        r = rejected[0]
        self.assertEqual(r["rejected_by"], "chromosome_lacks_min_anchors")
        self.assertEqual(r["chrom"], "NC_064602.1")
        self.assertEqual(r["start"], 78_296_974)
        self.assertAlmostEqual(r["identity"], 54.2, places=1)
        self.assertEqual(r["genome"], "Anopheles_gambiae.fna")

    def test_nearest_anchor_distance_is_reported(self):
        """The number that tells a reader whether 'not syntenic' means 'a bit off' or
        'a different chromosome entirely'."""
        hits, flanking = _anopheles_case()
        rejected = []
        isr.select_dispersed_goi_seeds(hits, flanking, min_chrom_anchors=3,
                                       rejected_out=rejected, genome_name="ag")
        r = rejected[0]
        # nearest same-chromosome anchor is gene-STRAY2 at 12.0 Mb; hit at 78.3 Mb
        self.assertTrue(r["nearest_anchor_same_chrom"])
        self.assertEqual(r["nearest_anchor_gene"], "gene-STRAY2")
        self.assertGreater(r["nearest_anchor_bp"], 66_000_000)

    def test_weak_hits_are_not_reported(self):
        """Only hits that CLEAR the quality bar are worth surfacing; sub-threshold
        noise would drown the signal."""
        flanking = [_hit(f"gene-FL{i}", "chr1", i * 100_000, i * 100_000 + 2000, pident=70.0)
                    for i in range(5)]
        weak = _hit("GOI_x", "chr9", 5_000_000, 5_000_100, pident=12.0, evalue=0.5, alnlen=20)
        rejected = []
        isr.select_dispersed_goi_seeds(flanking + [weak], flanking,
                                       rejected_out=rejected, genome_name="g")
        self.assertEqual(rejected, [])

    def test_no_envelope_anywhere_still_reports(self):
        """When no chromosome carries enough anchors the function returns early. A
        strong hit must still be surfaced — that is 'could not judge', not 'absent'."""
        flanking = [_hit("gene-A", "chr1", 1000, 2000, pident=70.0)]  # 1 < min_chrom_anchors
        strong = _hit("GOI_x", "chr7", 9_000_000, 9_000_500)
        rejected = []
        seeds = isr.select_dispersed_goi_seeds(flanking + [strong], flanking,
                                               min_chrom_anchors=3,
                                               rejected_out=rejected, genome_name="g")
        self.assertEqual(seeds, [])
        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0]["rejected_by"], "no_flanking_envelope_in_genome")

    def test_opt_out_records_nothing(self):
        hits, flanking = _anopheles_case()
        seeds = isr.select_dispersed_goi_seeds(hits, flanking, min_chrom_anchors=3,
                                               rejected_out=None)
        self.assertEqual(seeds, [])  # and no crash with the channel disabled


class TestReportSurfacesRejections(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        for sub in ("regions", "hits", "scores"):
            os.makedirs(os.path.join(self.dir, sub), exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _write_nonsyntenic(self, genome, identity=54.2):
        path = os.path.join(self.dir, "hits", f"{genome}{gr.NONSYNTENIC_SUFFIX}")
        cols = isr.NONSYNTENIC_COLUMNS
        row = {"genome": genome, "query": "GOI_Def", "chrom": "NC_064602.1",
               "start": 78296974, "end": 78297165, "identity": identity,
               "evalue": "5.845e-16", "alnlen": 192,
               "rejected_by": "chromosome_lacks_min_anchors",
               "nearest_anchor_gene": "gene-STRAY2", "nearest_anchor_bp": 66296974,
               "nearest_anchor_same_chrom": "True"}
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
        return path

    def test_loader_reads_the_tsv(self):
        self._write_nonsyntenic("Anopheles_gambiae.fna")
        recs, per_genome = gr.collect_nonsyntenic_candidates(
            [os.path.join(self.dir, "hits", f"Anopheles_gambiae.fna{gr.NONSYNTENIC_SUFFIX}")])
        self.assertEqual(len(recs), 1)
        self.assertIn("Anopheles_gambiae", per_genome)
        self.assertAlmostEqual(recs[0]["identity"], 54.2, places=1)

    def test_summary_ranks_best_hit_per_genome(self):
        recs = [
            {"genome": "g1", "identity": 54.2, "chrom": "c1", "start": 1, "end": 2,
             "rejected_by": "chromosome_lacks_min_anchors", "nearest_anchor_bp": 10,
             "nearest_anchor_same_chrom": True},
            {"genome": "g1", "identity": 71.7, "chrom": "c1", "start": 3, "end": 4,
             "rejected_by": "outside_flanking_envelope", "nearest_anchor_bp": 20,
             "nearest_anchor_same_chrom": True},
        ]
        block = gr.summarize_nonsyntenic_candidates(recs, {"g1": recs})
        self.assertEqual(block["total"], 2)
        self.assertAlmostEqual(block["best_per_genome"]["g1"]["identity"], 71.7, places=1)
        self.assertEqual(block["by_gate"]["outside_flanking_envelope"], 1)

    def test_malformed_file_does_not_break_the_report(self):
        bad = os.path.join(self.dir, "hits", f"broken{gr.NONSYNTENIC_SUFFIX}")
        with open(bad, "w") as fh:
            fh.write("not\ta\tvalid\theader\n\x00\x00\n")
        recs, per_genome = gr.collect_nonsyntenic_candidates([bad])
        self.assertEqual(recs, [])
        self.assertEqual(per_genome, {})

    def test_genome_moves_out_of_goi_absent(self):
        """The headline behaviour change: a genome holding a strong unplaceable hit
        must no longer be counted as simply 'absent'."""
        self._write_nonsyntenic("Anopheles_gambiae.fna")
        report = gr.build_report(self.dir)
        summary = report["summary"]
        self.assertIn("goi_found_but_not_syntenic", summary)
        self.assertIn("Anopheles_gambiae", summary["goi_found_but_not_syntenic"])
        self.assertNotIn("Anopheles_gambiae", summary["goi_absent_genomes"])
        self.assertEqual(summary["nonsyntenic_candidate_count"], 1)
        block = report["rejected_candidates"]
        self.assertEqual(block["total"], 1)
        self.assertEqual(block["by_gate"]["chromosome_lacks_min_anchors"], 1)
        # and it must survive a JSON round-trip
        json.loads(json.dumps(report))


if __name__ == "__main__":
    unittest.main()
