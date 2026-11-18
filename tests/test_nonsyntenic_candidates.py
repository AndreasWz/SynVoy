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

    def test_an_anchor_on_another_chromosome_gives_no_distance(self):
        """Fix item 6: with no anchor on the hit's chromosome the record carried the coordinate
        difference to an anchor on another chromosome (an arbitrary gene), and 0 bp when that
        anchor's coordinates overlapped the hit's. Now the distance and the gene are empty."""
        flanking = [_hit(f"gene-FL{i}", "chr1", 1_000 + i * 50_000, 3_000 + i * 50_000, pident=70.0)
                    for i in range(4)]
        for goi_start in (9_000_000, 1_500):   # far off, and overlapping an anchor's coordinates
            rejected = []
            isr.select_dispersed_goi_seeds([_hit("GOI_Def", "chr7", goi_start, goi_start + 300)] + flanking,
                                           flanking, min_chrom_anchors=3, rejected_out=rejected,
                                           genome_name="g")
            self.assertEqual(len(rejected), 1)
            r = rejected[0]
            self.assertFalse(r["nearest_anchor_same_chrom"])
            self.assertEqual((r["nearest_anchor_bp"], r["nearest_anchor_gene"]), ("", ""))

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
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
        return path

    def test_loader_reads_the_tsv(self):
        self._write_nonsyntenic("Anopheles_gambiae.fna")
        recs, per_genome = gr.collect_nonsyntenic_candidates(
            [os.path.join(self.dir, "hits", f"Anopheles_gambiae.fna{gr.NONSYNTENIC_SUFFIX}")])
        self.assertEqual(len(recs), 1)
        self.assertIn("Anopheles_gambiae", per_genome)
        self.assertAlmostEqual(recs[0]["identity"], 54.2, places=1)

    def test_an_empty_distance_is_read_as_null(self):
        """Fix item 6: the loader read an empty distance as 0, i.e. 'overlaps an anchor'."""
        path = self._write_nonsyntenic("Anopheles_gambiae.fna")
        lines = open(path).read().splitlines()
        cols = lines[0].split("\t")
        row = lines[1].split("\t")
        row[cols.index("nearest_anchor_bp")] = ""
        row[cols.index("nearest_anchor_same_chrom")] = "False"
        with open(path, "w") as fh:
            fh.write(lines[0] + "\n" + "\t".join(row) + "\n")
        recs, _ = gr.collect_nonsyntenic_candidates([path])
        self.assertIsNone(recs[0]["nearest_anchor_bp"])
        self.assertFalse(recs[0]["nearest_anchor_same_chrom"])

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


def _row(genome, chrom, start, query="GOI_Def", identity=60.0, bits=200.0, locus="locus_1", **kw):
    row = {"genome": genome, "query": query, "chrom": chrom, "start": start, "end": start + 300,
           "strand": "+", "qstart": 1, "qend": 100, "bits": bits, "identity": identity, "evalue": 1e-30,
           "alnlen": 300, "rejected_by": "chromosome_lacks_min_anchors", "nearest_anchor_gene": "gene-A",
           "nearest_anchor_bp": 2_000_000, "nearest_anchor_same_chrom": True, "outcome": "refused",
           "envelope": "chromosome_without", "anchors_on_sequence": 1, "genes_on_sequence": 2,
           "locus": locus, "modelled_here": "", "modelled_by": ""}
    row.update(kw)
    return row


def _write_table(directory, name, rows, cols=None):
    cols = cols or isr.NONSYNTENIC_COLUMNS
    path = os.path.join(directory, "hits", f"{name}{gr.NONSYNTENIC_SUFFIX}")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    return path


class TestRefusedRowsAndPlaces(unittest.TestCase):
    """Bug list item 11 (2026-11-16): each row carries strand, query span and bits; the report
    sums rows to places, and a place that a counted call of the run overlaps (another home
    locus's, a rescue pass's) is modelled, not "found, not syntenic"."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        for sub in ("regions", "hits", "scores"):
            os.makedirs(os.path.join(self.dir, sub), exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_a_record_carries_the_hit(self):
        hits, flanking = _anopheles_case()
        hits[-1].update(qstart=12, qend=75, bits=81.5, strand="-")
        rejected = []
        isr.select_dispersed_goi_seeds(hits, flanking, rejected_out=rejected, genome_name="ag")
        r = rejected[0]
        self.assertEqual((r["strand"], r["qstart"], r["qend"], r["bits"]), ("-", 12, 75, 81.5))
        self.assertEqual((r["outcome"], r["envelope"], r["anchors_on_sequence"], r["genes_on_sequence"]),
                         ("refused", "chromosome_without", 2, 2))
        self.assertEqual(isr.NONSYNTENIC_COLUMNS[:9], ["genome", "query", "chrom", "start", "end", "strand",
                                                       "qstart", "qend", "bits"])

    def test_rows_of_one_place_are_one_place(self):
        rows = [_row("G.fna", "c1", 1_000), _row("G.fna", "c1", 1_100, query="GOI_Def|exon_1", bits=90.0),
                _row("G.fna", "c1", 4_000, query="GOI_Def|G2_fna_b0_l1", bits=150.0),
                _row("G.fna", "c1", 900_000)]
        path = _write_table(self.dir, "locus_1_G.fna", rows)
        recs, per_genome = gr.collect_nonsyntenic_candidates([path])
        places = gr.nonsyntenic_places(recs)
        self.assertEqual([(p["start"], p["end"], p["n_rows"]) for p in places], [(1_000, 4_300, 3), (900_000, 900_300, 1)])
        self.assertEqual(places[0]["queries"], ["GOI_Def", "GOI_Def|G2_fna_b0_l1", "GOI_Def|exon_1"])
        self.assertEqual((places[0]["bits"], places[0]["strand"]), (200.0, "+"))
        # another home locus is another place
        self.assertEqual(len(gr.nonsyntenic_places(recs + [dict(recs[0], locus="locus_2")])), 3)

    def test_a_place_with_a_counted_call_of_another_locus_is_modelled(self):
        _write_table(self.dir, "locus_1_G.fna", [_row("G.fna", "c1", 1_000)])
        _write_table(self.dir, "locus_1_H.fna", [_row("H.fna", "c1", 1_000)])
        gff = ("##gff-version 3\n"
               "c1\tx\tmRNA\t901\t1500\t.\t+\t.\tID=GOI_Def|G_fna_b0_l1;SynVoy_Parent=GOI_Def;SynVoyRole=goi;"
               "Confidence=MEDIUM;GOIClass=probable_goi;EvidenceType=exon_annotation;Identity=70;QueryCoverage=0.9\n")
        with open(os.path.join(self.dir, "regions", "locus_2__G.fna.gff"), "w") as fh:
            fh.write(gff)
        report = gr.build_report(self.dir)
        block, summary = report["rejected_candidates"], report["summary"]
        places = {p["genome"]: p for p in block["places"]}
        self.assertEqual(places["G"]["modelled_in_run"]["provenance"], ["locus_2"])
        self.assertIs(places["G"]["modelled_in_run"]["other_locus"], True)
        self.assertIsNone(places["H"]["modelled_in_run"])
        self.assertEqual((summary["nonsyntenic_place_count"], summary["nonsyntenic_unmodelled_place_count"]), (2, 1))
        self.assertEqual(summary["goi_found_but_not_syntenic"], ["H"])
        self.assertEqual(summary["nonsyntenic_candidate_count"], 2)        # rows, as before
        self.assertEqual(sorted(block["best_per_genome"]), ["H"])

    def test_a_demoted_call_does_not_model_the_place(self):
        """Only a counted call (HIGH or MEDIUM, not demoted) models a place: a MEDIUM call that the
        coverage rule demotes (identity_coverage_decoupled) leaves it unmodelled, as a call that
        ownership gives to another home gene (paralog_not_goi) does."""
        _write_table(self.dir, "locus_1_G.fna", [_row("G.fna", "c1", 1_000)])
        gff = ("##gff-version 3\n"
               "c1\tx\tmRNA\t901\t1500\t.\t+\t.\tID=GOI_Def|G_fna_b0_l1;SynVoy_Parent=GOI_Def;SynVoyRole=goi;"
               "Confidence=MEDIUM;GOIClass=probable_goi;EvidenceType=exon_annotation;Identity=90;QueryCoverage=0.2\n")
        with open(os.path.join(self.dir, "regions", "locus_2__G.fna.gff"), "w") as fh:
            fh.write(gff)
        report = gr.build_report(self.dir)
        self.assertEqual(report["goi_dedup"]["records"][0]["goi_class"], "identity_coverage_decoupled")
        place = report["rejected_candidates"]["places"][0]
        self.assertIsNone(place["modelled_in_run"])
        self.assertEqual(report["summary"]["nonsyntenic_unmodelled_place_count"], 1)
        self.assertEqual(sorted(report["rejected_candidates"]["best_per_genome"]), ["G"])
        counted = dict(report["goi_dedup"]["records"][0], goi_class="probable_goi")
        recs, per_genome = gr.collect_nonsyntenic_candidates(
            [os.path.join(self.dir, "hits", f"locus_1_G.fna{gr.NONSYNTENIC_SUFFIX}")])
        block = gr.summarize_nonsyntenic_candidates(recs, per_genome, counted_records=[counted])
        self.assertEqual(block["unmodelled_place_count"], 0)
        self.assertIs(block["places"][0]["modelled_in_run"]["other_locus"], True)

    def test_a_table_of_the_old_format_still_loads(self):
        old = ["genome", "query", "chrom", "start", "end", "identity", "evalue", "alnlen", "rejected_by",
               "nearest_anchor_gene", "nearest_anchor_bp", "nearest_anchor_same_chrom"]
        path = _write_table(self.dir, "locus_3_G.fna", [_row("G.fna", "c1", 1_000)], cols=old)
        recs, per_genome = gr.collect_nonsyntenic_candidates([path])
        self.assertEqual(len(recs), 1)
        r = recs[0]
        self.assertEqual((r["locus"], r["strand"], r["qstart"], r["bits"], r["modelled_here"]),
                         ("locus_3", "", None, None, ""))
        self.assertEqual(gr.nonsyntenic_places(recs)[0]["n_rows"], 1)

    def test_dropped_rows_are_no_candidates(self):
        rows = [_row("G.fna", "c1", 1_000),
                _row("K.fna", "c1", 5_000, outcome="dropped", rejected_by="dispersed_identity", identity=30.0),
                _row("K.fna", "c2", 5_000, outcome="dropped", rejected_by="dead_zone")]
        _write_table(self.dir, "locus_1_mixed.fna", rows)
        with open(os.path.join(self.dir, "hits", "locus_1_K.fna.target_policy.tsv"), "w") as fh:
            fh.write("genome\tn_flanking_genes\tblocks_searched\trescue_dropped_outside_zones\n"
                     "K.fna\t9\t2\t417\n")
        report = gr.build_report(self.dir)
        block, summary = report["rejected_candidates"], report["summary"]
        self.assertEqual(block["total"], 1)
        self.assertEqual(summary["goi_found_but_not_syntenic"], ["G"])
        self.assertNotIn("K", summary["goi_absent_genomes"] + summary["goi_found_but_not_syntenic"])
        dropped = block["dropped"]
        self.assertEqual((dropped["total"], dropped["hits"], dropped["genomes"]), (2, 417, ["K"]))
        self.assertEqual(dropped["by_test"], {"dead_zone": 1, "dispersed_identity": 1})



def _goi_gff_row(chrom, start, end, conf, goi_class, identity=70, qcov=0.9, basis="", mid="GOI_Def|G_fna_b0_l1_m1"):
    return ("\t".join([chrom, "x", "mRNA", str(start), str(end), ".", "+", ".",
                       f"ID={mid};SynVoy_Parent=GOI_Def;SynVoyRole=goi;Confidence={conf};"
                       f"GOIClass={goi_class};EvidenceType=exon_annotation;Identity={identity};"
                       f"QueryCoverage={qcov}" + (f";Basis={basis}" if basis else "")]) + "\n")


class TestFoundButNotSyntenicCountsCalls(unittest.TestCase):
    """Bug list item 12: the list subtracted every genome with any GOI row, a LOW one too, so a
    genome with a refused hit and only uncounted models was in neither goi_absent_genomes nor
    goi_found_but_not_syntenic. Only a counted call (HIGH or MEDIUM, not demoted) takes a genome
    out of it now."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        for sub in ("regions", "hits", "scores"):
            os.makedirs(os.path.join(self.dir, sub), exist_ok=True)
        # The refused hit lies on c1 at 1,000 bp; the models below lie elsewhere on c9.
        _write_table(self.dir, "locus_1_G.fna", [_row("G.fna", "c1", 1_000, identity=70.0)])

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _gff(self, *rows):
        with open(os.path.join(self.dir, "regions", "locus_1__G.fna.gff"), "w") as fh:
            fh.write("##gff-version 3\n" + "".join(rows))
        return gr.build_report(self.dir)["summary"]

    def test_only_a_low_model(self):
        summary = self._gff(_goi_gff_row("c9", 50_001, 50_600, "LOW", "ambiguous_goi_family_member",
                                         identity=30, basis="not_significant"))
        self.assertEqual(summary["goi_found_but_not_syntenic"], ["G"])
        self.assertNotIn("G", summary["goi_absent_genomes"])

    def test_only_an_ambiguous_record(self):
        summary = self._gff(_goi_gff_row("c9", 50_001, 50_600, "AMBIGUOUS", "syntenic_candidate_unconfirmed",
                                         basis="not_specific"))
        self.assertEqual(summary["goi_found_but_not_syntenic"], ["G"])

    def test_a_demoted_call(self):
        summary = self._gff(_goi_gff_row("c9", 50_001, 50_600, "HIGH", "confident_goi", identity=60, qcov=0.2))
        self.assertEqual(summary["coverage_demoted_calls"], 1)
        self.assertEqual(summary["goi_found_but_not_syntenic"], ["G"])

    def test_a_counted_call_elsewhere(self):
        summary = self._gff(_goi_gff_row("c9", 50_001, 50_600, "MEDIUM", "probable_goi"))
        self.assertEqual(summary["goi_found_but_not_syntenic"], [])
        self.assertNotIn("G", summary["goi_absent_genomes"])


if __name__ == "__main__":
    unittest.main()
