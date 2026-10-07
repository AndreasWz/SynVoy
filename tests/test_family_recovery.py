#!/usr/bin/env python3
"""Family-recovery benchmark harness (scripts/benchmark/family_recovery.py).

The logic runs on SYNTHETIC truth so it works in public CI. The curated Ant_Venoms
dataset is unpublished and lives outside this repo; `TestRealDataset` pins the numbers
the design doc (docs/NEXT_SESSION_FAMILY_BENCHMARK.md) relies on and SKIPS when the
data is absent -- set SYNVOY_REQUIRE_FAMILY_DATA=1 to turn that skip into a failure.
"""
import importlib.util
import json
import os
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name):
    path = os.path.join(ROOT, "scripts", "benchmark", f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fr = _load("family_recovery")
sc = _load("score_coordinates")

HOME = "Home_sp"
FAM = "FAM_X"


def gene(gm, species, locus, chrom, start, end, family=FAM):
    return {"gm": gm, "species": species, "family": family, "locus": locus,
            "chrom": chrom, "start": start, "end": end, "strand": "+",
            "status": "complete-functional", "prot_len": 300}


def synthetic_genes():
    """Home has X_loc1 (seed) and X_loc2; X_tandem shares loc1's address; X_ant has no
    home member at all. Sp_c lacks the family entirely (the absence test)."""
    return [
        gene("H1", HOME, "X_loc1", "hchr1", 1000, 2000),
        gene("H2", HOME, "X_loc2", "hchr2", 5000, 6000),
        gene("A1", "Sp_a", "X_loc1", "achr1", 10000, 11000),
        gene("A2", "Sp_a", "X_tandem", "achr1", 14000, 15000),
        gene("A3", "Sp_a", "X_loc2", "achr7", 1000, 2000),
        gene("A4", "Sp_a", "X_ant", "achr9", 1000, 2000),
        gene("B1", "Sp_b", "X_loc1", "bchr1", 20000, 21000),
        gene("O1", "Sp_a", "Y_loc", "achr1", 30000, 31000, family="FAM_Y"),
    ]


ADDRESSES = {
    "X_loc1": ("FAM_X", "upA", "dnB"),
    "X_tandem": ("FAM_X", "upA", "dnB"),     # same address as X_loc1
    "X_loc2": ("FAM_X", "upC", "dnD"),
    "X_ant": None,                           # uncurated flanks
    "Y_loc": ("FAM_Y", "upA", "dnB"),        # same flanks, OTHER family: not equivalent
}


def call(species, chrom, start, end, conf="MEDIUM", cls="probable_goi", loci=("locus_1",),
         final=None, source="report"):
    c = {"species": species, "chrom": chrom, "start": start, "end": end,
         "confidence": conf, "goi_class": cls, "identity": "50", "query_cov": "0.9",
         "id": f"GOI|{species}_{start}", "home_loci": set(loci), "source": source}
    if source == "report":
        c["final"] = (conf in ("HIGH", "MEDIUM") and cls not in fr.DEMOTED_CLASSES
                      if final is None else final)
    return c


class TestTruthHelpers(unittest.TestCase):
    def test_species_key_matches_accessions_spelling(self):
        self.assertEqual(fr.species_key("Probolomyrmex sp. GAGA-0580"),
                         "Probolomyrmex_sp_GAGA-0580")
        self.assertEqual(fr.species_key("Apis mellifera"), "Apis_mellifera")

    def test_addresses_need_both_flanks(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("stable_id\tfamily_stable_id\tanchor_up\tanchor_dn\n"
                     "L1\tF\tup\tdn\nL2\tF\tup\t\nL3\tF\t\t\n")
            path = fh.name
        try:
            a = fr.load_addresses(path)
        finally:
            os.unlink(path)
        self.assertEqual(a["L1"], ("F", "up", "dn"))
        self.assertIsNone(a["L2"], "a one-sided anchor is not an address")
        self.assertIsNone(a["L3"])

    def test_tiers(self):
        t = {g["gm"]: g["tier"] for g in
             fr.assign_tiers(synthetic_genes(), FAM, HOME, {"X_loc1"}, ADDRESSES)}
        self.assertEqual(t["A1"], "home_locus")
        self.assertEqual(t["B1"], "home_locus")
        self.assertEqual(t["A2"], "home_address")
        self.assertEqual(t["A3"], "unsearched_home", "home has loc2, the run did not search it")
        self.assertEqual(t["A4"], "no_home_anchor")
        self.assertNotIn("H1", t, "home-genome genes are never targets")
        self.assertNotIn("O1", t, "other families are not tiered")

    def test_address_equivalence_is_within_family(self):
        """Y_loc shares X_loc1's flanks but is FAM_Y: must not become reachable."""
        genes = synthetic_genes() + [gene("A9", "Sp_a", "Y_loc", "achr1", 1, 2)]
        t = {g["gm"]: g["tier"] for g in
             fr.assign_tiers(genes, FAM, HOME, {"X_loc1"}, ADDRESSES)}
        self.assertEqual(t["A9"], "no_home_anchor")

    def test_label_spread(self):
        genes = [gene("p1", HOME, "L", "c1", 1, 10), gene("p2", HOME, "L", "c2", 1, 10),
                 gene("p3", HOME, "M", "c1", 1, 10), gene("p4", HOME, "M", "c1", 2_000_000,
                                                          2_000_010),
                 gene("p5", HOME, "N", "c1", 1, 10), gene("p6", HOME, "N", "c1", 500, 900)]
        self.assertTrue(fr.label_is_spread(genes, HOME, "L"), "two scaffolds")
        self.assertTrue(fr.label_is_spread(genes, HOME, "M"), "over 1 Mb")
        self.assertFalse(fr.label_is_spread(genes, HOME, "N"))


class TestPlan(unittest.TestCase):
    def test_one_row_per_home_seed_with_ceiling(self):
        rows = fr.build_plan(synthetic_genes(), ADDRESSES, HOME)
        x = {r["seed_gm"]: r for r in rows if r["family"] == FAM}
        self.assertEqual(set(x), {"H1", "H2"})
        self.assertEqual((x["H1"]["home_locus"], x["H1"]["home_address"],
                          x["H1"]["unsearched_home"], x["H1"]["no_home_anchor"]),
                         (2, 1, 1, 1))
        self.assertIn("shared_address_unresolvable", x["H1"]["warnings"])
        # Seeding from loc2 reaches only A3; loc1's genes become unsearched_home.
        self.assertEqual((x["H2"]["home_locus"], x["H2"]["unsearched_home"]), (1, 3))

    def test_family_without_home_member_is_flagged(self):
        rows = fr.build_plan(synthetic_genes(), ADDRESSES, HOME)
        y = [r for r in rows if r["family"] == "FAM_Y"]
        self.assertEqual(y[0]["warnings"], "no_home_member")

    def test_split_is_preregistered(self):
        """The hold-out families must not drift into tuning."""
        self.assertEqual({f for f, s in fr.SPLIT.items() if s == "holdout"},
                         {"FAM_VA", "FAM_SP", "FAM_APH"})
        self.assertEqual({f for f, s in fr.SPLIT.items() if s == "tune"},
                         {"FAM_SECP", "FAM_PLA2"})


class TestHomeLoci(unittest.TestCase):
    def test_home_locus_maps_to_curated_label(self):
        labels = fr.label_home_loci(
            {"locus_1": ("hchr1", 1200, 1300), "locus_2": ("hchr9", 1, 100)},
            synthetic_genes(), FAM, HOME)
        self.assertEqual(labels["locus_1"], "X_loc1")
        self.assertIsNone(labels["locus_2"], "an uncurated home region stays None")

    def test_slack_reaches_a_nearby_gene(self):
        labels = fr.label_home_loci({"locus_1": ("hchr1", 2500, 2600)},
                                    synthetic_genes(), FAM, HOME, slack=1000)
        self.assertEqual(labels["locus_1"], "X_loc1")


class RunDirMixin:
    """A minimal SynVoy output tree: one home locus, two target GFFs, a report."""

    def make_run(self, gff_by_species, records, bed=("hchr1", 1199, 1300)):
        self.tmp = tempfile.TemporaryDirectory()
        run = self.tmp.name
        pdir = os.path.join(run, "plot_inputs_synteny_block_locus_1")
        os.makedirs(pdir)
        os.makedirs(os.path.join(run, "regions"))
        with open(os.path.join(pdir, "locus_1.bed"), "w") as fh:
            fh.write("\t".join(map(str, bed)) + "\tgene_loc\t1e-10\t+\t60\n")
        with open(os.path.join(pdir, "home_genome.gff"), "w") as fh:
            fh.write("hchr1\tx\tmRNA\t1\t9\t.\t+\t.\tID=h;SynVoyRole=goi\n")
        for sp, lines in gff_by_species.items():
            with open(os.path.join(pdir, f"{sp}.fna.gff"), "w") as fh:
                fh.write("##gff-version 3\n" + "".join(l + "\n" for l in lines))
            open(os.path.join(run, "regions", f"{sp}.fna.regions.bed"), "w").close()
        with open(os.path.join(run, "synvoy_report.json"), "w") as fh:
            json.dump({"goi_dedup": {"records": records}}, fh)
        # Real runs publish per-task log DIRECTORIES named like GFFs.
        os.makedirs(os.path.join(run, "logs", "X_RESCUE_GOI_HULL__Sp_a.hull_rescue.gff"))
        return run

    def tearDown(self):
        if hasattr(self, "tmp"):
            self.tmp.cleanup()


def gff(chrom, ftype, start, end, attrs):
    return f"{chrom}\tsrc\t{ftype}\t{start}\t{end}\t50\t+\t.\t{attrs}"


class TestRunLoading(RunDirMixin, unittest.TestCase):
    def test_models_blocks_and_home_locus(self):
        run = self.make_run({"Sp_a": [
            gff("achr1", "mRNA", 10000, 10500,
                "ID=GOI_X|Sp_a_fna_b0_l1_exon_ann;SynVoyRole=goi;Confidence=MEDIUM"),
            gff("achr1", "gene", 14000, 14200,
                "ID=GOI_copy_1|Sp_a_fna_b0_l1;SynVoyRole=goi;Confidence=LOW"),
            gff("achr1", "mRNA", 60000, 61000,
                "ID=gene-F1|Sp_a_fna_b0_fl1_flank_ann;SynVoyRole=flanking"),
            gff("achr1", "gene", 70000, 71000, "ID=parentgene;SynVoyRole=goi"),
            gff("achr1", "mRNA", 70000, 71000,
                "ID=child;Parent=parentgene;SynVoyRole=goi;Confidence=LOW"),
        ]}, [])
        calls, spans = fr.load_run_models(fr.gather_run_gffs(run), HOME)
        kinds = sorted((c["start"], c["confidence"]) for c in calls)
        self.assertEqual(kinds, [(10000, "MEDIUM"), (14000, "LOW"), (70000, "LOW")],
                         "tandem gene counted; gene that parents an mRNA is not")
        self.assertTrue(all(c["home_loci"] == {"locus_1"} for c in calls))
        self.assertEqual(spans[("locus_1", "Sp_a", "achr1", 0)], [10000, 61000])
        self.assertNotIn(HOME, {c["species"] for c in calls})
        self.assertNotIn("home_genome", {c["species"] for c in calls})
        self.assertEqual(fr.load_home_loci(run), {"locus_1": ("hchr1", 1200, 1300)})
        self.assertEqual(fr.load_searched_species(run), {"Sp_a"})

    def test_report_records_final_rule(self):
        recs = [
            {"genome": "Sp_a", "chrom": "c", "start": 1, "end": 9, "confidence": "HIGH",
             "goi_class": "confident_goi", "provenance": ["locus_1"]},
            {"genome": "Sp_a", "chrom": "c", "start": 1, "end": 9, "confidence": "HIGH",
             "goi_class": "paralog_not_goi", "provenance": ["locus_1"]},
            {"genome": "Sp_a", "chrom": "c", "start": 1, "end": 9, "confidence": "MEDIUM",
             "goi_class": "identity_coverage_decoupled"},
            {"genome": "Sp_a", "chrom": "c", "start": 1, "end": 9,
             "confidence": "AMBIGUOUS", "goi_class": "syntenic_candidate_unconfirmed"},
        ]
        run = self.make_run({}, recs)
        self.assertEqual([r["final"] for r in fr.load_report_records(run)],
                         [True, False, False, False])

    def test_owning_locus_is_the_final_attribution(self):
        """§1m ownership may reassign a call; that, not the search provenance, is
        the ortholog assignment the report makes."""
        run = self.make_run({}, [
            {"genome": "Sp_a", "chrom": "c", "start": 1, "end": 9, "confidence": "HIGH",
             "goi_class": "confident_goi", "provenance": ["locus_1"],
             "owning_locus": "locus_2"},
            {"genome": "Sp_a", "chrom": "c", "start": 1, "end": 9, "confidence": "HIGH",
             "goi_class": "confident_goi", "provenance": ["locus_1", "locus_3"]},
        ])
        self.assertEqual([r["home_loci"] for r in fr.load_report_records(run)],
                         [{"locus_2"}, {"locus_1", "locus_3"}])

    def test_genome_name_forms(self):
        self.assertEqual(fr.genome_of("locus_3__Sp_a.fna.gff"), "Sp_a")
        self.assertEqual(fr.genome_of("Sp_a.hull_rescue.gff"), "Sp_a")
        self.assertEqual(fr.home_locus_of("/r/plot_inputs_synteny_block_locus_2/Sp_a.fna.gff"),
                         "locus_2")
        self.assertEqual(fr.home_locus_of("/r/x/locus_3__Sp_a.fna.gff"), "locus_3")
        self.assertIsNone(fr.home_locus_of("/r/x/Sp_a.fna.gff"))


class TestScoring(unittest.TestCase):
    """score_run on synthetic truth: verdicts, layers, call classes."""

    LABELS = {"locus_1": "X_loc1"}
    SEARCHED = {"Sp_a", "Sp_b", "Sp_c"}

    def run_score(self, report, gff_calls=(), spans=None, labels=None, **kw):
        return fr.score_run(synthetic_genes(), FAM, HOME, labels or self.LABELS, ADDRESSES,
                            list(gff_calls), spans or {}, list(report),
                            searched_species=kw.pop("searched", self.SEARCHED), **kw)

    @staticmethod
    def by_gm(rows):
        return {r["gm"]: r for r in rows}

    def test_recovered_and_unresolvable(self):
        rows, calls, s = self.run_score([
            call("Sp_a", "achr1", 10100, 10900),     # on A1 (X_loc1)
            call("Sp_a", "achr1", 14100, 14900),     # on A2 (X_tandem, same address)
        ])
        g = self.by_gm(rows)
        self.assertEqual(g["A1"]["verdict"], "recovered")
        self.assertEqual(g["A2"]["verdict"], "unresolvable_by_synteny",
                         "a shared-address hit is never scored as an error")
        self.assertEqual([c["call_class"] for c in calls], ["on_locus", "on_address"])
        self.assertEqual(s["locus_precision"], 1.0)

    def test_family_only_when_found_from_the_wrong_home_locus(self):
        labels = {"locus_1": "X_loc1", "locus_2": "X_loc2"}
        rows, calls, _s = self.run_score(
            [call("Sp_a", "achr1", 10100, 10900, loci=("locus_2",))], labels=labels)
        self.assertEqual(self.by_gm(rows)["A1"]["verdict"], "family_only")
        self.assertEqual(calls[0]["call_class"], "on_other_family_locus")

    def test_layer_a_no_block_nearby(self):
        spans = {("locus_1", "Sp_a", "achr1", 0): [500000, 600000]}
        rows, _c, s = self.run_score([], spans=spans)
        a1 = self.by_gm(rows)["A1"]
        self.assertEqual((a1["verdict"], a1["layer"]), ("missed", "A_neighbourhood"))
        self.assertEqual(a1["nearest_block_bp"], 489000)

    def test_layer_b_block_but_no_model(self):
        spans = {("locus_1", "Sp_a", "achr1", 0): [30000, 60000]}
        rows, _c, _s = self.run_score([], spans=spans)
        self.assertEqual(self.by_gm(rows)["A1"]["layer"], "B_modelling")

    def test_layer_c_model_below_floor(self):
        spans = {("locus_1", "Sp_a", "achr1", 0): [10000, 60000]}
        low = call("Sp_a", "achr1", 10100, 10900, conf="LOW", source="gff")
        rows, _c, _s = self.run_score([], gff_calls=[low], spans=spans)
        self.assertEqual(self.by_gm(rows)["A1"]["layer"], "C_classification")

    def test_layer_d_demoted_by_adjudication(self):
        spans = {("locus_1", "Sp_a", "achr1", 0): [10000, 60000]}
        med = call("Sp_a", "achr1", 10100, 10900, conf="MEDIUM", source="gff")
        demoted = call("Sp_a", "achr1", 10100, 10900, conf="MEDIUM", cls="paralog_not_goi")
        rows, calls, _s = self.run_score([demoted], gff_calls=[med], spans=spans)
        self.assertEqual(self.by_gm(rows)["A1"]["layer"], "D_adjudication")
        self.assertEqual(calls, [], "a demoted record is not a final call")

    def test_unreachable_tiers_are_not_misses(self):
        rows, _c, s = self.run_score([])
        g = self.by_gm(rows)
        self.assertEqual(g["A4"]["verdict"], "unreachable")
        self.assertEqual(g["A3"]["verdict"], "unreachable")
        self.assertEqual(g["A4"]["layer"], "")
        self.assertEqual(s["reachable_genes"], 3, "A1, A2, B1")
        self.assertEqual(sum(s["layer_census"].values()), 3)

    def test_found_without_anchor(self):
        rows, _c, _s = self.run_score([call("Sp_a", "achr9", 1100, 1900)])
        self.assertEqual(self.by_gm(rows)["A4"]["verdict"], "found_without_anchor")

    def test_not_searched_species_excluded(self):
        rows, _c, s = self.run_score([], searched={"Sp_a"})
        self.assertEqual(self.by_gm(rows)["B1"]["verdict"], "not_searched")
        self.assertEqual(s["reachable_genes"], 2)

    def test_call_classes_off_truth_and_other_family(self):
        _r, calls, s = self.run_score([
            call("Sp_a", "achr1", 30100, 30900),   # on the FAM_Y gene
            call("Sp_a", "achr5", 1, 900),         # nothing curated there
            call("Sp_c", "cchr1", 1, 900),         # species lacks the family
        ])
        self.assertEqual([c["call_class"] for c in calls],
                         ["other_family", "off_truth", "off_truth"])
        self.assertEqual(calls[0]["other_families"], "FAM_Y")
        self.assertEqual(s["final_calls_in_species_lacking_family"], 1)
        self.assertEqual(s["species_lacking_family"], ["Sp_c"])

    def test_single_locus_run_attributes_unprovenanced_records(self):
        rows, _c, _s = self.run_score([call("Sp_a", "achr1", 10100, 10900, loci=())])
        self.assertEqual(self.by_gm(rows)["A1"]["verdict"], "recovered")

    def test_summary_recalls(self):
        _r, _c, s = self.run_score([call("Sp_a", "achr1", 10100, 10900),
                                    call("Sp_a", "achr1", 14100, 14900)])
        self.assertEqual(s["locus_recall"], 0.5, "A1 of {A1, B1}")
        self.assertEqual(s["address_recall"], round(2 / 3, 4), "A1+A2 of {A1, A2, B1}")
        self.assertEqual(s["family_species_recall"], 0.5, "Sp_a of {Sp_a, Sp_b}")


class TestAmendment1(unittest.TestCase):
    """Unassigned labels and seed-unrelated truth genes (added after the first run)."""

    def genes(self):
        return synthetic_genes() + [
            gene("T1", "Sp_t", "TETLOC_0001", "tchr1", 1000, 2000),
            gene("V1", "Sp_v", "scaf:vchr3", "vchr3", 1000, 2000),
        ]

    def test_unassigned_labels_get_their_own_tier(self):
        t = {g["gm"]: g["tier"] for g in
             fr.assign_tiers(self.genes(), FAM, HOME, {"X_loc1"}, ADDRESSES)}
        self.assertEqual((t["T1"], t["V1"]), ("unassigned_label", "unassigned_label"))
        self.assertEqual(t["A4"], "no_home_anchor", "ordinary labels unchanged")

    def test_call_on_unassigned_gene_leaves_the_locus_precision_denominator(self):
        rows, calls, s = fr.score_run(
            self.genes(), FAM, HOME, {"locus_1": "X_loc1"}, ADDRESSES, [], {},
            [call("Sp_t", "tchr1", 1100, 1900), call("Sp_a", "achr1", 10100, 10900)])
        self.assertEqual({c["call_class"] for c in calls}, {"on_unassigned_label", "on_locus"})
        self.assertEqual(s["locus_precision"], 1.0, "1 on_locus of 1 judged call")
        self.assertEqual(s["family_precision"], 1.0)
        self.assertEqual({r["gm"]: r["verdict"] for r in rows}["T1"], "unassigned_found")

    def test_seed_unrelated_truth_is_flagged_and_reported_both_ways(self):
        evalues = {"A1": 1e-80, "B1": 0.7}
        rows, _c, s = fr.score_run(
            synthetic_genes(), FAM, HOME, {"locus_1": "X_loc1"}, ADDRESSES, [], {},
            [call("Sp_a", "achr1", 10100, 10900)], seed_evalues=evalues)
        g = {r["gm"]: r for r in rows}
        self.assertEqual((g["A1"]["seed_unrelated"], g["B1"]["seed_unrelated"]), ("", "yes"))
        self.assertEqual(s["seed_unrelated_truth"], ["B1"])
        self.assertEqual(s["locus_recall"], 0.5, "raw: A1 of {A1, B1}")
        self.assertEqual(s["locus_recall_excl_seed_unrelated"], 1.0)

    def test_assembly_limited_truth_is_flagged_and_reported_both_ways(self):
        """A gene on a contig too short to carry a flanking gene is not a search miss."""
        lengths = {("Sp_a", "achr1"): 2_000_000, ("Sp_b", "bchr1"): 3_510}
        rows, _c, s = fr.score_run(
            synthetic_genes(), FAM, HOME, {"locus_1": "X_loc1"}, ADDRESSES, [], {},
            [call("Sp_a", "achr1", 10100, 10900)], contig_lengths=lengths)
        g = {r["gm"]: r for r in rows}
        self.assertEqual((g["A1"]["assembly_limited"], g["B1"]["assembly_limited"]),
                         ("", "yes"))
        self.assertEqual(g["B1"]["contig_bp"], 3510)
        self.assertEqual(g["A3"]["contig_bp"], "", "no index for that contig: unflagged")
        self.assertEqual(s["assembly_limited_truth"], ["B1"])
        self.assertEqual((s["locus_recall"], s["locus_recall_excl_assembly_limited"]),
                         (0.5, 1.0))

    def test_record_census_counts_every_tier_against_the_truth(self):
        """Deliverable 3: AMBIGUOUS / LOW / demoted records are judged too, per phylo verdict."""
        recs = [dict(call("Sp_a", "achr1", 10100, 10900, conf="AMBIGUOUS",
                          cls="syntenic_candidate_unconfirmed"), phylo_verdict="phylo_concordant"),
                dict(call("Sp_a", "achr5", 1, 900, conf="AMBIGUOUS",
                          cls="syntenic_candidate_unconfirmed"), phylo_verdict="phylo_discordant"),
                dict(call("Sp_b", "bchr1", 20100, 20900, conf="MEDIUM", cls="paralog_not_goi"),
                     phylo_verdict="phylo_concordant")]
        _r, calls, s = fr.score_run(synthetic_genes(), FAM, HOME, {"locus_1": "X_loc1"},
                                    ADDRESSES, [], {}, recs)
        self.assertEqual(calls, [], "none of these is a final call")
        bc = s["record_census"]["by_confidence"]
        self.assertEqual(bc["AMBIGUOUS"], {"on_locus": 1, "off_truth": 1})
        self.assertEqual(bc["MEDIUM->paralog_not_goi"], {"on_locus": 1})
        bp = s["record_census"]["by_phylo_verdict"]
        self.assertEqual(bp["phylo_discordant"], {"off_truth": 1})
        self.assertEqual(bp["phylo_concordant"], {"on_locus": 2})

    @unittest.skipUnless(importlib.util.find_spec("parasail"), "parasail not installed")
    def test_seed_similarity_separates_homologs_from_noise(self):
        seed = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ" * 4
        mutated = seed[:100] + "A" * 20 + seed[120:]
        unrelated = "GGSGGSWPPPHHHHCCCCNNNNQQQQ" * 10
        ev = fr.seed_similarity(seed, {"h": mutated, "u": unrelated, "x": ""},
                                ["h", "u", "x", "missing"])
        self.assertLess(ev["h"], 1e-50)
        self.assertGreater(ev["u"], fr.MAX_SEED_EVALUE)
        self.assertNotIn("x", ev)
        self.assertNotIn("missing", ev)

    def test_seed_evalue_is_length_aware(self):
        """Amendment 2: the raw floor of 100 called real 77-aa secapin orthologs
        (scores 53-99) unrelated. The same score means more on a shorter comparison."""
        self.assertLess(fr.sw_evalue(53, 77, 110), fr.MAX_SEED_EVALUE)
        self.assertGreater(fr.sw_evalue(38, 560, 600), fr.MAX_SEED_EVALUE,
                           "the Vespa APYR curation error stays flagged")
        self.assertLess(fr.sw_evalue(60, 80, 100), fr.sw_evalue(60, 800, 1000))


class TestStage(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        fd = os.path.join(d, "final_dataset")
        os.makedirs(os.path.join(fd, "proteins"))
        with open(os.path.join(fd, "EXPORT_MANIFEST.tsv"), "w") as fh:
            fh.write("gm\tspecies\trole\tfamily\tlocus\tscaffold\tstart\tend\tstrand\t"
                     "n_exons\tprot_len\tis_venom\tstructural_status\n"
                     "H1\tHome sp\tout\tFAM_X\tX_loc1\thchr1\t1\t9\t+\t1\t5\t1\tcf\n"
                     "A1\tSp a\tant\tFAM_X\tX_loc1\tachr1\t1\t9\t+\t1\t5\t1\tcf\n"
                     "B1\tSp b.\tant\tFAM_X\tX_loc1\tbchr1\t1\t9\t+\t1\t5\t1\tcf\n"
                     "A2\tSp a\tant\tFAM_X\tX_loc2\tachr2\t1\t9\t+\t1\t5\t1\tcf\n")
        with open(os.path.join(fd, "proteins", "Home_sp.faa"), "w") as fh:
            fh.write(">H0 other\nMMMM\n>H1 family=FAM_X\nMKV\nLLA\n>H2\nMA\n")
        g = os.path.join(d, "genomes")
        for acc, files in {"GH": ("GH_genomic.fna", "genomic.gff"),
                           "GA": ("GA_genomic.fna",)}.items():
            os.makedirs(os.path.join(g, acc))
            for f in files:
                open(os.path.join(g, acc, f), "w").close()
        with open(os.path.join(g, "ACCESSIONS.tsv"), "w") as fh:
            fh.write("species\taccession\nHome_sp\tGH\nSp_a\tGA\nSp_b\tGB\n")
        self.data = d

    def tearDown(self):
        self.tmp.cleanup()

    def test_stage_writes_query_and_species_named_targets(self):
        out = os.path.join(self.data, "run")
        res = fr.stage(self.data, "H1", out, home_species="Home_sp")
        with open(os.path.join(out, "query.faa")) as fh:
            self.assertEqual(fh.read(), ">H1 family=FAM_X\nMKVLLA\n")
        self.assertEqual(res["staged"], ["Sp_a"])
        self.assertEqual(res["missing"], ["Sp_b"], "no genome dir for GB")
        self.assertTrue(os.path.islink(os.path.join(out, "targets", "Sp_a.fna")))
        self.assertTrue(os.path.islink(os.path.join(out, "home", "Home_sp.gff")))

    def test_seed_must_be_from_the_home_genome(self):
        with self.assertRaises(SystemExit):
            fr.stage(self.data, "A1", os.path.join(self.data, "run"), home_species="Home_sp")

    def test_curated_seed_writes_seed_json_and_checks_family(self):
        out = os.path.join(self.data, "run")
        fr.stage(self.data, "H1", out, home_species="Home_sp", family="FAM_X")
        with open(os.path.join(out, "seed.json")) as fh:
            info = json.load(fh)
        self.assertEqual((info["kind"], info["seed"], info["locus"]),
                         ("curated", "H1", "X_loc1"))
        with self.assertRaises(SystemExit):
            fr.stage(self.data, "H1", out, home_species="Home_sp", family="FAM_Y")

    def write_home_annotation(self):
        gh = os.path.join(self.data, "genomes", "GH")
        with open(os.path.join(gh, "genomic.gff"), "w") as fh:
            fh.write("##gff-version 3\n"
                     + gff("hchr5", "gene", 500, 900, "ID=gene-LOCR;Name=LOCR;gene=LOCR") + "\n"
                     + gff("hchr5", "CDS", 500, 600, "ID=c1;gene=LOCR;protein_id=XP_1.1") + "\n"
                     + gff("hchr5", "CDS", 700, 900, "ID=c2;gene=LOCR;protein_id=XP_2.1") + "\n"
                     + gff("hchr5", "CDS", 800, 900, "ID=c3;gene=LOCR;protein_id=XP_2.1") + "\n"
                     + gff("hchr6", "CDS", 1, 90, "ID=c4;gene=OTHER;protein_id=XP_9.1") + "\n")
        with open(os.path.join(gh, "protein.faa"), "w") as fh:
            fh.write(">XP_1.1 short\nMKV\n>XP_2.1 long\nMKVLL\n>XP_9.1\nMAAAAAAA\n")

    def test_home_gene_seed_is_the_longest_isoform_with_the_family_label(self):
        """Amendment 2 RBH arm: a home gene that is NOT a curated member."""
        self.write_home_annotation()
        out = os.path.join(self.data, "rbh")
        res = fr.stage(self.data, "LOCR", out, home_species="Home_sp", family="FAM_X")
        with open(os.path.join(out, "query.faa")) as fh:
            self.assertEqual(fh.read().split("\n")[1], "MKVLL", "longest isoform")
        with open(os.path.join(out, "seed.json")) as fh:
            info = json.load(fh)
        self.assertEqual((info["kind"], info["protein_id"], info["locus"], info["chrom"],
                          info["start"], info["end"]),
                         ("home_gene", "XP_2.1", "X_loc1", "hchr5", 500, 900))
        self.assertEqual(res["staged"], ["Sp_a"])

    def test_home_gene_seed_needs_a_family_and_a_real_label(self):
        self.write_home_annotation()
        out = os.path.join(self.data, "rbh")
        with self.assertRaises(SystemExit):
            fr.stage(self.data, "LOCR", out, home_species="Home_sp")
        with self.assertRaises(SystemExit):
            fr.stage(self.data, "LOCR", out, home_species="Home_sp", family="FAM_X",
                     label="not_a_home_label")
        with self.assertRaises(SystemExit):
            fr.stage(self.data, "NOSUCH", out, home_species="Home_sp", family="FAM_X")

    def test_home_gene_seed_may_stand_for_a_label_home_lacks(self):
        """The address audit can put a home gene at a label no home member carries."""
        self.write_home_annotation()
        out = os.path.join(self.data, "rbh")
        fr.stage(self.data, "LOCR", out, home_species="Home_sp", family="FAM_X",
                 label="X_loc2")
        with open(os.path.join(out, "seed.json")) as fh:
            self.assertEqual(json.load(fh)["locus"], "X_loc2")

    def test_curated_seed_relabelled_by_the_address_audit(self):
        out = os.path.join(self.data, "run")
        fr.stage(self.data, "H1", out, home_species="Home_sp", label="X_loc2")
        with open(os.path.join(out, "seed.json")) as fh:
            info = json.load(fh)
        self.assertEqual((info["kind"], info["locus"], info["manifest_locus"]),
                         ("curated", "X_loc2", "X_loc1"))
        with self.assertRaises(SystemExit):
            fr.stage(self.data, "H1", out, home_species="Home_sp", label="Y_loc9")

    def test_orf_seed(self):
        orf = os.path.join(self.data, "orf.faa")
        with open(orf, "w") as fh:
            fh.write(">whatever\nMRFQ\nLYIL\n")
        out = os.path.join(self.data, "orf")
        fr.stage(self.data, "ORF_X", out, home_species="Home_sp", family="FAM_X",
                 label="X_loc2", seed_faa=orf, seed_coords="hchr9:100-450:-")
        with open(os.path.join(out, "query.faa")) as fh:
            self.assertEqual(fh.read(), ">ORF_X home ORF seed for FAM_X\nMRFQLYIL\n")
        with open(os.path.join(out, "seed.json")) as fh:
            info = json.load(fh)
        self.assertEqual((info["kind"], info["locus"], info["chrom"], info["start"],
                          info["end"], info["strand"], info["prot_len"]),
                         ("orf", "X_loc2", "hchr9", 100, 450, "-", 8))
        for bad in ({"family": None}, {"seed_coords": None}):
            kw = dict(family="FAM_X", seed_faa=orf, seed_coords="hchr9:100-450")
            kw.update(bad)
            with self.assertRaises(SystemExit):
                fr.stage(self.data, "ORF_X", out, home_species="Home_sp", **kw)
        with self.assertRaises(SystemExit, msg="an ORF seed may not reuse a GM id"):
            fr.stage(self.data, "H1", out, home_species="Home_sp", family="FAM_X",
                     seed_faa=orf, seed_coords="hchr9:100-450")

    def test_parse_coords(self):
        self.assertEqual(fr.parse_coords("NC_1.1:5-9:-"), ("NC_1.1", 5, 9, "-"))
        self.assertEqual(fr.parse_coords("scaf:x:5-9"), ("scaf:x", 5, 9, "."))
        with self.assertRaises(SystemExit):
            fr.parse_coords("NC_1.1:5")


class TestApplySeed(unittest.TestCase):
    def test_relabelled_curated_seed_moves_the_searched_label(self):
        genes = synthetic_genes()
        home_gm = next(g["gm"] for g in genes if g["species"] == HOME)
        info = {"kind": "curated", "seed": home_gm, "locus": "X_loc2",
                "manifest_locus": "X_loc1"}
        out = fr.apply_seed(genes, info, HOME)
        self.assertEqual(next(g for g in out if g["gm"] == home_gm)["locus"], "X_loc2")
        self.assertEqual(next(g for g in genes if g["gm"] == home_gm)["locus"], "X_loc1",
                         "input untouched")

    def test_plain_curated_seed_and_no_seed_change_nothing(self):
        genes = synthetic_genes()
        self.assertIs(fr.apply_seed(genes, {"kind": "curated", "seed": "H"}, HOME), genes)
        self.assertIs(fr.apply_seed(genes, None, HOME), genes)

    def test_orf_seed_joins_as_a_home_member(self):
        info = {"kind": "orf", "seed": "ORF_X", "family": FAM, "locus": "X_loc1",
                "chrom": "hchr9", "start": 100, "end": 450, "strand": "-"}
        out = fr.apply_seed(synthetic_genes(), info, HOME)
        self.assertEqual(out[-1]["gm"], "SEED:ORF_X")
        self.assertEqual(out[-1]["status"], "orf_seed")


class TestSeedHomeMember(unittest.TestCase):
    def test_home_gene_seed_labels_the_run_home_locus(self):
        info = {"kind": "home_gene", "seed": "LOCR", "family": FAM, "locus": "X_loc1",
                "chrom": "hchr5", "start": 500, "end": 900, "strand": "+"}
        member = fr.seed_home_member(info, HOME)
        genes = synthetic_genes()
        home = {"locus_1": ("hchr5", 520, 880)}
        self.assertIsNone(fr.label_home_loci(home, genes, FAM, HOME)["locus_1"],
                          "without the seed member the RBH locus is uncurated")
        self.assertEqual(fr.label_home_loci(home, genes + [member], FAM, HOME)["locus_1"],
                         "X_loc1")
        tiers = fr.assign_tiers(genes + [member], FAM, HOME, {"X_loc1"}, ADDRESSES)
        self.assertNotIn("SEED:LOCR", {t["gm"] for t in tiers}, "home genes are not scored")

    def test_curated_seed_adds_nothing(self):
        self.assertIsNone(fr.seed_home_member({"kind": "curated", "seed": "H1"}))
        self.assertIsNone(fr.seed_home_member(None))


class TestSeedcheck(unittest.TestCase):
    def test_rbh_summary(self):
        top = {"a": "SEED", "b": "SEED", "c": "ALT", "d": "ALT", "e": "ALT", "f": None}
        self.assertEqual(fr.rbh_summary(list(top), top, "SEED"), (0.4, "ALT", 3))
        self.assertEqual(fr.rbh_summary(["a", "b"], top, "SEED"), (1.0, None, 0))
        self.assertEqual(fr.rbh_summary(["f"], top, "SEED"), (None, None, 0),
                         "no target with a home hit: undefined, not 0")

    def test_protein_gene_map_first_cds_and_trimmed_product(self):
        with tempfile.NamedTemporaryFile("w", suffix=".gff", delete=False) as fh:
            fh.write("#x\n"
                     + gff("c1", "CDS", 700, 800, "gene=G;protein_id=P1;product=venom X%2C isoform 2") + "\n"
                     + gff("c1", "CDS", 100, 200, "gene=G;protein_id=P1;product=venom X") + "\n"
                     + gff("c1", "mRNA", 1, 900, "gene=G;protein_id=P9") + "\n")
            path = fh.name
        try:
            m = fr.protein_gene_map(path)
        finally:
            os.unlink(path)
        self.assertEqual(m, {"P1": ("G", "c1", 700, "venom X")})


class TestCallcheck(unittest.TestCase):
    def test_extract_windows_one_pass_multiline(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fna", delete=False) as fh:
            fh.write(">c1 desc\nACGTA\nCGTAC\n>c2\nTTTT\n>c3\nGGGG\n")
            path = fh.name
        try:
            w = fr.extract_windows(path, {"c1": [(2, 7, "k1"), (9, 10, "k2")],
                                          "c3": [(1, 2, "k3")], "absent": [(1, 5, "k4")]})
        finally:
            os.unlink(path)
        self.assertEqual(w, {"k1": "CGTACG", "k2": "AC", "k3": "GG"})

    def test_best_home_genes_collapses_isoforms(self):
        pmap = {"P1": ("G1", "c", 1, ""), "P1b": ("G1", "c", 1, ""), "P2": ("G2", "c", 9, "")}
        hits = [("P1", 90.0, 0.5), ("P1b", 120.0, 0.5), ("P2", 100.0, 0.4),
                ("Pnone", 30.0, 0.3)]
        self.assertEqual(fr.best_home_genes(hits, pmap), [("G1", 120.0), ("G2", 100.0)])

    def test_call_origin(self):
        self.assertEqual(fr.call_origin([("S", 150.0), ("P", 110.0)], "S"), "query_gene")
        self.assertEqual(fr.call_origin([("P", 150.0), ("S", 110.0)], "S"),
                         "other_home_gene")
        self.assertEqual(fr.call_origin([("S", 40.0)], "S"), "no_home_hit",
                         "below MIN_CALL_BITS nothing explains the window")
        self.assertEqual(fr.call_origin([], "S"), "no_home_hit")


class TestScoreCoordinatesTandemCopies(unittest.TestCase):
    """The §1y scorer missed `gene`-only tandem-copy GOI models (LRZ det_rep R2)."""

    def test_tandem_gene_is_a_call_and_pairs_are_not_doubled(self):
        text = ("##gff-version 3\n"
                + gff("OV1", "gene", 15634947, 15635084,
                      "ID=GOI_copy_2|x_b0_l1;SynVoyRole=goi;Confidence=LOW") + "\n"
                + gff("OV1", "gene", 100, 200, "ID=g;SynVoyRole=goi") + "\n"
                + gff("OV1", "mRNA", 100, 200, "ID=m;Parent=g;SynVoyRole=goi") + "\n")
        with tempfile.NamedTemporaryFile("w", suffix=".gff", delete=False) as fh:
            fh.write(text)
            path = fh.name
        try:
            calls = sc.load_calls_from_gff([path])
        finally:
            os.unlink(path)
        self.assertEqual(sorted(c["start"] for c in calls), [100, 15634947])


DATA = fr.DEFAULT_DATA
HAVE_DATA = os.path.isfile(os.path.join(DATA, "final_dataset", "EXPORT_MANIFEST.tsv"))


@unittest.skipUnless(HAVE_DATA or os.environ.get("SYNVOY_REQUIRE_FAMILY_DATA") == "1",
                     f"curated family dataset not present at {DATA}")
class TestRealDataset(unittest.TestCase):
    """Pins the facts the design doc was built on (verified 2026-09-10).

    FAMILY-LEVEL AGGREGATES ONLY: the curated set is unpublished and this repo is
    public, so no locus label or gene-model id may appear here (they live in the
    gitignored tests/benchmark_truth/family_benchmark_private.md).
    """

    @classmethod
    def setUpClass(cls):
        fd = os.path.join(DATA, "final_dataset")
        cls.all = fr.load_manifest(os.path.join(fd, "EXPORT_MANIFEST.tsv"),
                                   exclude_families=())
        cls.genes = fr.load_manifest(os.path.join(fd, "EXPORT_MANIFEST.tsv"))
        cls.addr = fr.load_addresses(os.path.join(fd, "LOCUS_ANCHORS.tsv"))
        cls.by_family = {}
        for r in fr.build_plan(cls.genes, cls.addr):
            if r.get("seed_gm"):
                cls.by_family.setdefault(r["family"], []).append(r)

    def best(self, family):
        """The seed with the most reachable genes -- the one the benchmark runs."""
        return max(self.by_family[family],
                   key=lambda r: (r["home_locus"] + r["home_address"], r["seed_gm"]))

    def test_counts(self):
        self.assertEqual(len(self.all), 939)
        no_gr = [g for g in self.all
                 if g["family"] not in ("FAM_GR1", "FAM_GR2", "FAM_GR3", "FAM_GR7")]
        self.assertEqual(len(no_gr), 747)

    def test_the_three_shared_addresses(self):
        groups = {}
        for loc, a in self.addr.items():
            if a is not None:
                groups.setdefault(a, set()).add(loc)
        shared = {k[0]: len(v) for k, v in groups.items() if len(v) > 1
                  and k[0] not in fr.DEFAULT_EXCLUDE}
        self.assertEqual(shared, {"FAM_APH": 3, "FAM_HYAL": 2, "FAM_SP": 2})

    def test_sp_is_reachable_only_through_the_shared_address(self):
        """No ant SP gene carries an Apis SP label: the best seed reaches 30 genes, all
        through the shared address; the other Apis SP seed reaches nothing."""
        b = self.best("FAM_SP")
        self.assertEqual((b["home_locus"], b["home_address"]), (0, 30))
        self.assertEqual(min(r["reachable_species"] for r in self.by_family["FAM_SP"]), 0)

    def test_single_copy_ladder_is_reachable(self):
        for fam in ("FAM_APYR", "FAM_SPIN", "FAM_DPP4"):
            b = self.best(fam)
            self.assertGreaterEqual(b["home_locus"], 28, fam)
            self.assertLessEqual(b["no_home_anchor"], 1, fam)

    def test_pla2_home_label_is_not_an_address(self):
        self.assertTrue(all("home_label_not_an_address" in r["warnings"]
                            for r in self.by_family["FAM_PLA2"]))

    def test_no_home_anchor_total(self):
        """Of the 724 non-Apis genes, even with every Apis locus of the family searched:
        412 sit at loci with no Apis LABEL; 82 of those share an address with an Apis
        locus and stay reachable; 43 carry no locus assignment at all (amendment 1);
        the remaining 287 are unreachable by construction."""
        total = unanchored = unassigned = 0
        for fam in self.by_family:
            home_labels = {g["locus"] for g in self.genes
                           if g["family"] == fam and g["species"] == fr.DEFAULT_HOME}
            tiers = fr.assign_tiers(self.genes, fam, fr.DEFAULT_HOME, home_labels, self.addr)
            total += len(tiers)
            unanchored += sum(t["tier"] == "no_home_anchor" for t in tiers)
            unassigned += sum(t["tier"] == "unassigned_label" for t in tiers)
        self.assertEqual((unanchored, unassigned, total), (287, 43, 724))


if __name__ == "__main__":
    unittest.main()
