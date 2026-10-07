#!/usr/bin/env python3
"""§1y — freeze the melittin benchmark by COORDINATE, not by species presence.

Why this file exists
--------------------
SynVoy hit the *T. bicarinatum* melittin (`A0A6M3Z554_1`) in 8 runs and then
silently lost it: the 2026-06-30 run's nearest call is 2,543 bp from the nearest
truth model, and nobody noticed, because the benchmark scored SPECIES PRESENCE and
the species flag stayed up. These tests make that regression detectable.

The two call coordinates below are frozen from real runs on disk and re-verified
2026-08-30 against `local_runs/`:

  mel_det_rep2               OV788322.1:15,634,953-15,635,084  34.7%  MEDIUM  (HIT)
  melittin_val_sw_20260630   OV788322.1:15,638,808-15,638,999  31.9%  MEDIUM  (LOST)

They are literals rather than fixtures because `local_runs/` is not committed; the
point is to freeze the *numbers*, so a future change that stops recovering the
first, or stops flagging the second, fails here.

DATA AVAILABILITY -- read before trusting a green CI run
--------------------------------------------------------
`tests/benchmark_truth/` is gitignored (`.gitignore:66`): the curated coordinates
are unpublished annotation data. So the tests that need the real truth table CANNOT
run on a fresh CI checkout, and they SKIP there. What still runs in public CI is the
scoring logic on synthetic truth -- which is worth having, but it is not the freeze.

To make the freeze binding where the data does exist (a local run, or a private
runner with the truth table mounted), set:

    SYNVOY_REQUIRE_TRUTH=1 pytest tests/test_coordinate_benchmark.py

which turns "truth table missing" from a skip into a failure. Do that in any job you
intend to actually gate on §1y.
"""
import importlib.util
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRUTH = os.path.join(ROOT, "tests", "benchmark_truth", "melittin_loci.tsv")
HAVE_TRUTH = os.path.isfile(TRUTH)
REQUIRE_TRUTH = os.environ.get("SYNVOY_REQUIRE_TRUTH") == "1"


def truth_required(cls):
    """Skip when the gitignored truth table is absent -- unless explicitly required."""
    if HAVE_TRUTH:
        return cls
    if REQUIRE_TRUTH:
        return unittest.skipIf(False, "")(cls)   # let it run and fail loudly
    return unittest.skip(
        f"curated truth table not present ({TRUTH} is gitignored); "
        "set SYNVOY_REQUIRE_TRUTH=1 to make this a failure")(cls)


def synthetic_truth():
    """A stand-in truth set covering every truth_kind, for the logic tests.

    Deliberately fake coordinates: these tests must run in public CI, where the real
    curated table is absent.
    """
    return [
        {"species": "Spec_a", "assembly": "GCA_1", "chrom": "chr1",
         "start": 1000, "end": 1854, "strand": "+", "model_id": "gene_a",
         "truth_kind": "gene", "confidence": "gold", "source": "synthetic", "notes": ""},
        {"species": "Spec_a", "assembly": "GCA_1", "chrom": "chr1",
         "start": 9000, "end": 9100, "strand": "-", "model_id": "frag_a",
         "truth_kind": "fragment", "confidence": "putative", "source": "synthetic",
         "notes": ""},
        {"species": "Spec_b", "assembly": "GCA_2", "chrom": "chr9",
         "start": 500, "end": 90000, "strand": "+", "model_id": "syntenic_locus",
         "truth_kind": "locus_gene_lost", "confidence": "gold", "source": "synthetic",
         "notes": ""},
        {"species": "Spec_c", "assembly": "GCA_3", "chrom": "chrZ",
         "start": 100, "end": 900, "strand": "+", "model_id": "gene_c",
         "truth_kind": "gene", "confidence": "gold", "source": "synthetic", "notes": ""},
    ]


def _load():
    path = os.path.join(ROOT, "scripts", "benchmark", "score_coordinates.py")
    spec = importlib.util.spec_from_file_location("score_coordinates", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sc = _load()

# The ant melittin truth model, from the curated manual annotation
# (Ant_Venoms/data/annotations/toxins_in_Tbic_OV788322_annotations.gff).
TBIC_CHROM = "OV788322.1"
TBIC_START, TBIC_END = 15634953, 15635807

HIT_CALL = {"chrom": TBIC_CHROM, "start": 15634953, "end": 15635084, "strand": "+",
            "identity": "34.7", "confidence": "MEDIUM", "goi_class": "probable_goi",
            "query_cov": "0.686", "exons": "1", "model_status": "partial",
            "id": "GOI_Melt|mel_det_rep2", "src": "frozen"}
LOST_CALL = {"chrom": TBIC_CHROM, "start": 15638808, "end": 15638999, "strand": "+",
             "identity": "31.9", "confidence": "MEDIUM", "goi_class": "probable_goi",
             "query_cov": "", "exons": "1", "model_status": "partial",
             "id": "GOI_Melt|melittin_val_sw", "src": "frozen"}


class TestOverlapMath(unittest.TestCase):
    def test_overlap_is_inclusive(self):
        self.assertEqual(sc.overlap_bp(10, 20, 10, 20), 11)
        self.assertEqual(sc.overlap_bp(10, 20, 20, 30), 1)
        self.assertEqual(sc.overlap_bp(10, 20, 21, 30), 0)

    def test_distance_zero_when_overlapping(self):
        self.assertEqual(sc.distance_bp(10, 20, 15, 25), 0)

    def test_distance_symmetric_gap(self):
        # 21..29 exclusive => 9 bases between them, both directions.
        self.assertEqual(sc.distance_bp(10, 20, 30, 40), 10)
        self.assertEqual(sc.distance_bp(30, 40, 10, 20), 10)


@truth_required
class TestTruthSet(unittest.TestCase):
    def setUp(self):
        self.truth = sc.load_truth(TRUTH)

    def test_ant_melittin_coordinates_are_frozen(self):
        """The curated A0A6M3Z554_1 span must not drift."""
        m = [t for t in self.truth if t["model_id"] == "A0A6M3Z554_1"]
        self.assertEqual(len(m), 1)
        self.assertEqual((m[0]["chrom"], m[0]["start"], m[0]["end"]),
                         (TBIC_CHROM, TBIC_START, TBIC_END))
        self.assertEqual(m[0]["truth_kind"], "gene")

    def test_lost_loci_are_marked(self):
        """Melipona and Tetragonula are the deliberate false-positive tests."""
        lost = {t["species"] for t in self.truth if t["truth_kind"] == "locus_gene_lost"}
        self.assertEqual(lost, {"Melipona_beecheii", "Tetragonula_carbonaria"})

    def test_region_spans_are_not_typed_as_curated_genes(self):
        """A 110 kb syntenic region must not be scored for overlap FRACTION."""
        for t in self.truth:
            if t["model_id"] == "syntenic_locus":
                self.assertIn(t["truth_kind"], ("locus_gene_present", "locus_gene_lost"))

    def test_every_model_has_a_source(self):
        for t in self.truth:
            self.assertTrue(t["source"].strip(), f"{t['model_id']} has no provenance")


@truth_required
class TestAntMelittinRegression(unittest.TestCase):
    """The §1y freeze: the hit must stay a hit and the loss must stay detectable."""

    def setUp(self):
        self.truth = sc.load_truth(TRUTH)
        self.kw = dict(searched_species={"Tetramorium_bicarinatum"},
                       searched_chroms={TBIC_CHROM},
                       exclude_species=["Apis_mellifera"])

    def _row(self, rows, model_id):
        return next(r for r in rows if r["model_id"] == model_id)

    def test_mel_det_rep2_recovers_the_ant_melittin(self):
        rows, _calls, summary = sc.score(self.truth, [HIT_CALL], **self.kw)
        r = self._row(rows, "A0A6M3Z554_1")
        self.assertEqual(r["hit"], "yes")
        self.assertEqual(r["overlap_bp"], 132)
        self.assertEqual(r["start_offset"], 0, "start was exact to the base")
        self.assertAlmostEqual(r["overlap_frac_of_truth"], 132 / 855, places=4)
        self.assertEqual(summary["curated_models_hit"], 1)

    def test_the_hit_is_only_15_percent_of_the_gene(self):
        """Presence scoring called this a clean TP. It is one exon of 855 bp."""
        rows, _c, summary = sc.score(self.truth, [HIT_CALL], **self.kw)
        r = self._row(rows, "A0A6M3Z554_1")
        self.assertLess(r["overlap_frac_of_truth"], 0.20)
        self.assertEqual(summary["models_hit_over_50pct"], 0,
                         "nothing here recovered half the gene")

    def test_val_sw_run_lost_the_gene(self):
        rows, _c, summary = sc.score(self.truth, [LOST_CALL], **self.kw)
        r = self._row(rows, "A0A6M3Z554_1")
        self.assertEqual(r["hit"], "no")
        self.assertEqual(r["overlap_bp"], 0)
        self.assertEqual(summary["curated_models_hit"], 0)

    def test_val_sw_nearest_truth_model_is_2543_bp(self):
        """The number in §1y: nearest model is the fragment, not A0A6M3Z554_1."""
        rows, _c, _s = sc.score(self.truth, [LOST_CALL], **self.kw)
        self.assertEqual(self._row(rows, "A0A6M3Z541_fragment")["nearest_call_distance"],
                         2543)
        self.assertEqual(self._row(rows, "A0A6M3Z554_1")["nearest_call_distance"], 3001)

    def test_region_recall_can_rise_while_the_gene_is_lost(self):
        """The trap §1y is about: more regions hit, actual recovery zero."""
        _r, _c, hit = sc.score(self.truth, [HIT_CALL], **self.kw)
        _r2, _c2, lost = sc.score(self.truth, [LOST_CALL], **self.kw)
        self.assertEqual(hit["curated_models_hit"], 1)
        self.assertEqual(lost["curated_models_hit"], 0)


class TestScoreability(unittest.TestCase):
    """Runs on SYNTHETIC truth so it works in public CI without the curated table."""

    def setUp(self):
        self.truth = synthetic_truth()

    def test_unsearched_species_is_not_a_miss(self):
        rows, _c, summary = sc.score(self.truth, [], searched_species={"Spec_a"},
                                     searched_chroms={"chr1"})
        other = [r for r in rows if r["species"] != "Spec_a"]
        self.assertTrue(all(r["hit"] == "not_searched" for r in other))
        self.assertEqual(summary["truth_models_real"], 2)   # gene_a + frag_a

    def test_assembly_mismatch_is_not_a_miss(self):
        """Species WAS a target, but none of its curated scaffolds appear => the run
        used a different assembly. Scoring that as a miss invents a failure."""
        rows, _c, summary = sc.score(self.truth, [], searched_species={"Spec_a"},
                                     searched_chroms={"chr_other"})
        sa = [r for r in rows if r["species"] == "Spec_a"]
        self.assertTrue(all(r["hit"] == "assembly_mismatch" for r in sa))
        self.assertEqual(summary["truth_models_assembly_mismatch"], 2)
        self.assertEqual(summary["truth_models_real"], 0)

    def test_genuine_miss_when_another_scaffold_of_that_species_was_seen(self):
        """Same assembly confirmed by a sibling scaffold => a real miss."""
        truth = synthetic_truth() + [
            {"species": "Spec_a", "assembly": "GCA_1", "chrom": "chr2", "start": 1,
             "end": 100, "strand": "+", "model_id": "gene_a2", "truth_kind": "gene",
             "confidence": "gold", "source": "synthetic", "notes": ""}]
        rows, _c, _s = sc.score(truth, [], searched_species={"Spec_a"},
                                searched_chroms={"chr2"})
        self.assertEqual(next(r for r in rows if r["model_id"] == "gene_a")["hit"], "no")

    def test_call_on_a_lost_locus_is_a_false_positive(self):
        call = dict(HIT_CALL, chrom="chr9", start=1000, end=1200)
        _r, call_rows, summary = sc.score(self.truth, [call],
                                          searched_species={"Spec_b"},
                                          searched_chroms={"chr9"})
        self.assertEqual(summary["false_positive_calls_on_lost_loci"], 1)
        self.assertEqual(call_rows[0]["on_lost_locus"], "yes")

    def test_lost_locus_is_excluded_from_recall(self):
        """A lost locus is a false-positive trap, not a model to recover."""
        _r, _c, summary = sc.score(self.truth, [], searched_species={"Spec_b"},
                                   searched_chroms={"chr9"})
        self.assertEqual(summary["truth_models_real"], 0)
        self.assertEqual(summary["lost_loci_total"], 1)

    def test_region_span_excluded_from_fraction_stats(self):
        """A hit on a region counts for recall but not for overlap fraction."""
        truth = [dict(synthetic_truth()[2], truth_kind="locus_gene_present")]
        call = dict(HIT_CALL, chrom="chr9", start=1000, end=1200)
        _r, _c, summary = sc.score(truth, [call], searched_species={"Spec_b"},
                                   searched_chroms={"chr9"})
        self.assertEqual(summary["truth_models_hit"], 1)
        self.assertEqual(summary["curated_models_searched"], 0)
        self.assertEqual(summary["median_overlap_frac"], 0.0)

    def test_confidence_floor_filters_calls(self):
        low = dict(HIT_CALL, chrom="chr1", start=1000, end=1131, confidence="LOW")
        _r, _c, summary = sc.score(self.truth, [low], min_confidence="MEDIUM",
                                   searched_species={"Spec_a"},
                                   searched_chroms={"chr1"})
        self.assertEqual(summary["calls_total"], 0)

    def test_overlap_fraction_is_reported_not_a_boolean(self):
        """The §1y lesson, on synthetic data: 132 bp of an 855 bp gene is 15%."""
        call = dict(HIT_CALL, chrom="chr1", start=1000, end=1131)
        rows, _c, summary = sc.score(self.truth, [call], searched_species={"Spec_a"},
                                     searched_chroms={"chr1"})
        r = next(x for x in rows if x["model_id"] == "gene_a")
        self.assertEqual(r["hit"], "yes")
        self.assertEqual(r["overlap_bp"], 132)
        self.assertAlmostEqual(r["overlap_frac_of_truth"], 132 / 855, places=4)
        self.assertEqual(summary["models_hit_over_50pct"], 0)


class TestGffParsing(unittest.TestCase):
    def test_only_goi_mrna_features_are_scored(self):
        import tempfile
        gff = ("##gff-version 3\n"
               "OV788322.1\tx\tmRNA\t100\t200\t50\t+\t.\tID=a;SynVoyRole=goi;Confidence=HIGH\n"
               "OV788322.1\tx\tmRNA\t300\t400\t50\t+\t.\tID=b;SynVoyRole=flanking\n"
               "OV788322.1\tx\tCDS\t100\t200\t.\t+\t0\tID=a_CDS;Parent=a\n")
        with tempfile.NamedTemporaryFile("w", suffix=".gff", delete=False) as fh:
            fh.write(gff)
            path = fh.name
        try:
            calls = sc.load_calls_from_gff([path])
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0]["start"], 100)
            self.assertEqual(calls[0]["confidence"], "HIGH")
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
