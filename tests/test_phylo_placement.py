#!/usr/bin/env python3
"""Regressions for the phylogenetic placement check (docs/STATE_OF_THE_PROJECT.md F9).

The check exists because three per-call approaches failed: identity, coverage, flanking
support and home-paralog RBH cannot tell an ortholog from a lineage-specific paralog,
and the measured proof is that the highest-identity call in the melittin benchmark
(93.8 %, *Bombus terrestris*, melittin known lost) outranks every true positive.

The question this asks instead — does sequence divergence track SPECIES divergence? —
is only answerable across the whole set of calls.
"""
import importlib.util
import os
import unittest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin")


def _load():
    path = os.path.join(BIN, "phylo_placement_check.py")
    spec = importlib.util.spec_from_file_location("phylo_placement_check", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pp = _load()


def _rec(species, sp_dist, identity):
    return {"model_id": f"GOI_X|{species}_fna_b0_l1_fallback", "species": species,
            "species_distance": sp_dist, "identity": identity,
            "sw_score": identity * 3, "target_len_aa": 70,
            "gene_divergence": round(1.0 - identity / 100.0, 4)}


class TestIdParsing(unittest.TestCase):
    def test_genome_stem_recovered_from_model_id(self):
        for mid, want in (
            ("GOI_Melt|Colletes_gigas_fa_b0_l1_fallback", "Colletes_gigas"),
            ("GOI_Melt|Apis_cerana_fna_b2_l3_fallback", "Apis_cerana"),
            ("gene-LOC1|Bombus_impatiens_fna_b0_fl1_flank_ann", "Bombus_impatiens"),
        ):
            self.assertEqual(pp.species_from_model_id(mid), want)

    def test_id_without_pipe_is_tolerated(self):
        self.assertEqual(pp.species_from_model_id("Apis_mellifera_fna_b0_l1_x"),
                         "Apis_mellifera")

    def test_bare_block_locus_suffix_resolves(self):
        # The tandem-copy path emits "GOI_copy_1|<stem>_b0_l1" with NOTHING after the
        # locus id. The first regex required a trailing token, so these resolved to a
        # species absent from sorted_genomes.txt -> no distance -> every verdict came
        # back insufficient_data. That silently killed 7 of 8 calls in the first F9
        # cluster run (job 5777904) while the run still reported success.
        self.assertEqual(pp.species_from_model_id("GOI_copy_1|Apis_cerana_fna_b0_l1"),
                         "Apis_cerana")
        self.assertEqual(
            pp.species_from_model_id("GOI_copy_2|Tetragonula_carbonaria_fa_b0_l1"),
            "Tetragonula_carbonaria")

    def test_known_stems_resolve_every_emitting_path(self):
        # Enumerated from real run GFFs; the rescue paths were never handled at all,
        # and fig4 attributes 6 of 7 HIGH calls to the hull rescue.
        stems = ["Apis_cerana", "Colletes_gigas", "cow", "GCF_002204515.2"]
        for mid, want in (
            ("GOI_copy_1|Apis_cerana_fna_b0_l1", "Apis_cerana"),
            ("GOI_Melt|Colletes_gigas_WUUM01000001.1_hull_rescue", "Colletes_gigas"),
            ("GOI_rescue_cow_locus_Reg1_0", "cow"),
            ("GOI_copy_1|GCF_002204515_2_fna_b0_l1", "GCF_002204515.2"),
        ):
            self.assertEqual(pp.species_from_model_id(mid, known_stems=stems), want)

    def test_longest_stem_wins_so_prefixes_do_not_collide(self):
        stems = ["Apis_cerana", "Apis_cerana_2"]
        self.assertEqual(
            pp.species_from_model_id("GOI_copy_1|Apis_cerana_2_fna_b0_l1", known_stems=stems),
            "Apis_cerana_2")
        self.assertEqual(
            pp.species_from_model_id("GOI_copy_1|Apis_cerana_fna_b0_l1", known_stems=stems),
            "Apis_cerana")


class TestRecordBuilding(unittest.TestCase):
    def _faa(self, entries):
        import tempfile
        fh = tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False)
        for name, seq in entries:
            fh.write(f">{name}\n{seq}\n")
        fh.close()
        return fh.name

    def test_short_fragments_are_skipped(self):
        # The first F9 run fitted on a 7-aa "protein" at 40 % identity. Over 7 residues
        # that identity is noise, and one such point can set the slope for a whole run.
        path = self._faa([("GOI_X|Apis_cerana_fna_b0_l1", "MKV"),
                          ("GOI_X|Colletes_gigas_fa_b0_l1", "MKVLIAALFVAVMTLSHGAWQ")])
        try:
            recs = pp.build_records("MKVLIAALFVAVMTLSHGAWQ",
                                    path, {"Apis_cerana": 2.0, "Colletes_gigas": 6.0},
                                    min_target_len=20)
            self.assertEqual([r["species"] for r in recs], ["Colletes_gigas"])
        finally:
            os.unlink(path)

    def test_species_resolution_uses_the_known_stems(self):
        path = self._faa([("GOI_copy_1|Apis_cerana_fna_b0_l1", "MKVLIAALFVAVMTLSHGAWQ")])
        try:
            recs = pp.build_records("MKVLIAALFVAVMTLSHGAWQ", path, {"Apis_cerana": 2.0})
            self.assertEqual(recs[0]["species"], "Apis_cerana")
            self.assertEqual(recs[0]["species_distance"], 2.0)
        finally:
            os.unlink(path)


class TestSpeciesDistances(unittest.TestCase):
    def test_parses_and_strips_extensions(self):
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write("Apis_cerana.fna\t2\nEuglossa_dilemma.fa\t6\nbad_line\n")
            p = fh.name
        try:
            d = pp.load_species_distances(p)
            self.assertEqual(d, {"Apis_cerana": 2.0, "Euglossa_dilemma": 6.0})
        finally:
            os.unlink(p)

    def test_missing_file_is_empty_not_an_error(self):
        self.assertEqual(pp.load_species_distances("/nonexistent/x.txt"), {})


class TestTheilSen(unittest.TestCase):
    def test_recovers_a_clean_line(self):
        xs = [1, 2, 3, 4, 5]
        ys = [2 * x + 1 for x in xs]
        slope, intercept = pp.theil_sen(xs, ys)
        self.assertAlmostEqual(slope, 2.0)
        self.assertAlmostEqual(intercept, 1.0)

    def test_is_not_dragged_by_one_outlier(self):
        """A least-squares fit would tilt toward the outlier and hide it. That is the
        whole reason for a median-of-slopes estimator here."""
        xs = [1, 2, 3, 4, 5, 6]
        ys = [2 * x for x in xs]
        ys[-1] = 100.0
        slope, _ = pp.theil_sen(xs, ys)
        self.assertAlmostEqual(slope, 2.0, places=6)

    def test_no_spread_returns_none(self):
        self.assertEqual(pp.theil_sen([3, 3, 3], [1, 2, 3])[0], None)


class TestEvaluate(unittest.TestCase):
    def test_orthologs_tracking_species_distance_are_concordant(self):
        recs = [_rec("A", 1, 95), _rec("B", 2, 90), _rec("C", 3, 85),
                _rec("D", 4, 80), _rec("E", 5, 75)]
        out, summary = pp.evaluate(recs)
        self.assertEqual(summary["n_discordant"], 0)
        self.assertTrue(all(r["verdict"] == "phylo_concordant" for r in out))

    def test_the_bombus_case_is_flagged(self):
        """The scenario this check exists for: a CLOSE species whose recovered
        sequence is FAR more diverged than its species distance predicts, because it
        is a paralog that duplicated before the speciation."""
        recs = [_rec("Apis_cerana", 2, 97), _rec("Apis_florea", 2, 84),
                _rec("Euglossa", 6, 58), _rec("Colletes", 7, 37),
                _rec("Xylocopa", 7, 41),
                _rec("Bombus_terrestris", 3, 53)]     # close species, far sequence
        out, summary = pp.evaluate(recs, min_calls=4, z_threshold=2.0)
        bad = [r for r in out if r["verdict"] == "phylo_discordant"]
        self.assertEqual([r["species"] for r in bad], ["Bombus_terrestris"])
        self.assertEqual(summary["n_discordant"], 1)

    def test_better_than_expected_is_never_flagged(self):
        """One-sided by design: a call CLOSER to the query than its species distance
        predicts is a well-conserved ortholog, not a paralog."""
        recs = [_rec("A", 1, 90), _rec("B", 2, 80), _rec("C", 3, 70),
                _rec("D", 4, 60), _rec("E", 5, 99)]   # far species, very close sequence
        out, _ = pp.evaluate(recs)
        self.assertEqual([r["verdict"] for r in out if r["species"] == "E"],
                         ["phylo_concordant"])

    def test_too_few_calls_reports_insufficient_not_a_guess(self):
        recs = [_rec("A", 1, 90), _rec("B", 2, 80)]
        out, summary = pp.evaluate(recs, min_calls=4)
        self.assertTrue(all(r["verdict"] == "insufficient_data" for r in out))
        self.assertEqual(summary["n_insufficient"], 2)
        self.assertEqual(summary["n_fitted"], 0)

    def test_no_species_distance_spread_is_insufficient(self):
        recs = [_rec(s, 5, i) for s, i in zip("ABCDE", (90, 80, 70, 60, 50))]
        out, summary = pp.evaluate(recs, min_calls=4)
        self.assertTrue(all(r["verdict"] == "insufficient_data" for r in out))

    def test_near_constant_fallback_distances_are_refused(self):
        # PHYLO_SORT emits near-constant sentinel distances (993-1000, or a flat 999)
        # when the taxonomy walk cannot separate the targets — seen in local_runs/
        # melittin_full_qw (14 of 19 genomes at exactly 999). Fitting a slope through
        # that is fitting noise, so a fit needs real levels, not merely >= 2 values.
        recs = [_rec(s, d, i) for s, d, i in
                zip("ABCDE", (999, 999, 999, 999, 1000), (90, 80, 70, 60, 50))]
        out, summary = pp.evaluate(recs, min_calls=4, min_distance_levels=3)
        self.assertTrue(all(r["verdict"] == "insufficient_data" for r in out))
        self.assertEqual(summary["distance_levels"], 2)
        self.assertEqual(summary["distance_range"], "999-1000")

    def test_real_distance_spread_is_fitted(self):
        recs = [_rec(s, d, i) for s, d, i in
                zip("ABCDE", (2, 6, 9, 12, 13), (95, 85, 78, 70, 68))]
        _out, summary = pp.evaluate(recs, min_calls=4, min_distance_levels=3)
        self.assertEqual(summary["n_fitted"], 5)
        self.assertEqual(summary["distance_levels"], 5)

    def test_call_with_unknown_species_distance_is_not_fitted(self):
        recs = [_rec("A", 1, 95), _rec("B", 2, 90), _rec("C", 3, 85),
                _rec("D", 4, 80), _rec("U", None, 30)]
        out, summary = pp.evaluate(recs)
        u = next(r for r in out if r["species"] == "U")
        self.assertEqual(u["verdict"], "insufficient_data")
        self.assertEqual(summary["n_fitted"], 4)

    def test_identical_residuals_do_not_divide_by_zero(self):
        recs = [_rec("A", 1, 90), _rec("B", 2, 80), _rec("C", 3, 70), _rec("D", 4, 60)]
        out, _ = pp.evaluate(recs)
        self.assertTrue(all(r["verdict"] == "phylo_concordant" for r in out))
        self.assertTrue(all(isinstance(r["z_score"], float) for r in out))


if __name__ == "__main__":
    unittest.main()


# ── report-side integration ────────────────────────────────────────────────

def _load_gr():
    path = os.path.join(BIN, "generate_report.py")
    spec = importlib.util.spec_from_file_location("generate_report_pp", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gr = _load_gr()


def _dedup(confidence, mrna_id, goi_class="syntenic_candidate_unconfirmed"):
    return {"genome": "G", "chrom": "c1", "start": 100, "end": 400,
            "confidence": confidence, "identity": "40.0", "goi_class": goi_class,
            "query_cov": "0.5", "mrna_id": mrna_id}


def _row(mrna_id, verdict, z="0.4"):
    return {"model_id": mrna_id, "species": "G", "species_distance": "3",
            "identity": "40.0", "gene_divergence": "0.60",
            "expected_divergence": "0.58", "z_score": z, "verdict": verdict}


class TestReportIntegration(unittest.TestCase):
    def test_discordant_call_raises_a_flag(self):
        rec = _dedup("MEDIUM", "m1", "probable_goi")
        flags, summary = gr.apply_phylo_placement(
            {"records": [rec]}, {"m1": _row("m1", "phylo_discordant", "3.1")})
        self.assertEqual(len(flags), 1)
        self.assertEqual(flags[0]["type"], "phylo_discordant")
        self.assertEqual(summary["discordant"], 1)
        self.assertEqual(rec["phylo_verdict"], "phylo_discordant")

    def test_promotion_is_off_by_default_but_counted(self):
        """The measurement-first pattern: report what it WOULD do before doing it."""
        rec = _dedup("AMBIGUOUS", "m1")
        _flags, summary = gr.apply_phylo_placement(
            {"records": [rec]}, {"m1": _row("m1", "phylo_concordant")})
        self.assertEqual(summary["would_promote"], 1)
        self.assertEqual(summary["promoted"], 0)
        self.assertEqual(rec["confidence"], "AMBIGUOUS", "must not move by default")

    def test_promotion_when_enabled_restores_medium(self):
        rec = _dedup("AMBIGUOUS", "m1")
        _flags, summary = gr.apply_phylo_placement(
            {"records": [rec]}, {"m1": _row("m1", "phylo_concordant")}, promote=True)
        self.assertEqual(summary["promoted"], 1)
        self.assertEqual(rec["confidence"], "MEDIUM")
        self.assertEqual(rec["goi_class"], "probable_goi")
        self.assertTrue(rec["phylo_promoted"])

    def test_discordant_ambiguous_is_never_promoted(self):
        rec = _dedup("AMBIGUOUS", "m1")
        _flags, summary = gr.apply_phylo_placement(
            {"records": [rec]}, {"m1": _row("m1", "phylo_discordant", "4.0")},
            promote=True)
        self.assertEqual(rec["confidence"], "AMBIGUOUS")
        self.assertEqual(summary["promoted"], 0)

    def test_insufficient_data_changes_nothing(self):
        rec = _dedup("AMBIGUOUS", "m1")
        flags, summary = gr.apply_phylo_placement(
            {"records": [rec]}, {"m1": _row("m1", "insufficient_data", "")},
            promote=True)
        self.assertEqual(flags, [])
        self.assertEqual(rec["confidence"], "AMBIGUOUS")
        self.assertEqual(summary["insufficient_data"], 1)

    def test_unmatched_record_is_untouched(self):
        rec = _dedup("AMBIGUOUS", "no_such_model")
        _flags, summary = gr.apply_phylo_placement(
            {"records": [rec]}, {"m1": _row("m1", "phylo_concordant")}, promote=True)
        self.assertEqual(summary["evaluated"], 0)
        self.assertNotIn("phylo_verdict", rec)

    def test_no_rows_is_a_clean_noop(self):
        rec = _dedup("AMBIGUOUS", "m1")
        flags, summary = gr.apply_phylo_placement({"records": [rec]}, {})
        self.assertEqual((flags, summary["evaluated"]), ([], 0))
