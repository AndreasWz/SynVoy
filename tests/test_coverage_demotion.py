#!/usr/bin/env python3
"""Regressions for the 2026-08-22 coverage demotion (F8 follow-up).

A GOI call with high %identity over a tiny slice of the query is a short local
window, not an ortholog. ``build_self_consistency`` has flagged these as
``identity_coverage_decoupled`` since QW3, but only advisorily -- the headline
counted them anyway. ``apply_coverage_demotion`` now relabels them and drops them
from the headline, exactly as §1m does for ``paralog_not_goi``.

The load-bearing case is test_unknown_coverage_is_not_demoted: before
``rescue_goi_hull.py`` learned to emit ``QueryCoverage``, EVERY rescue model had
no recorded coverage, so demoting on absence would have removed the one true call
(Apis cerana melittin, 97.1%) along with the artefacts.

See docs/STATE_OF_THE_PROJECT.md F8 and A.8.
"""
import importlib.util
import os

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gr = _load("generate_report", "generate_report.py")
hull = _load("rescue_goi_hull", "rescue_goi_hull.py")


def _rec(identity, qcov, confidence="HIGH", goi_class="synteny_hull_rescue"):
    return {
        "genome": "G1", "chrom": "c1", "start": 100, "end": 132,
        "confidence": confidence, "identity": str(identity),
        "goi_class": goi_class, "query_cov": ("" if qcov is None else str(qcov)),
    }


# ── the rule fires on known-low coverage ────────────────────────────────────

def test_high_identity_low_coverage_is_demoted():
    dedup = {"records": [_rec(100.0, 0.14)]}
    summary = gr.apply_coverage_demotion(dedup)
    assert summary["demoted"] == 1
    rec = dedup["records"][0]
    assert rec["goi_class"] == "identity_coverage_decoupled"
    assert rec["original_goi_class"] == "synteny_hull_rescue"
    assert rec["coverage_demoted"] is True


def test_medium_calls_are_demoted_too():
    dedup = {"records": [_rec(95.0, 0.10, confidence="MEDIUM")]}
    assert gr.apply_coverage_demotion(dedup)["demoted"] == 1


# ── and must NOT fire otherwise ─────────────────────────────────────────────

def test_unknown_coverage_is_not_demoted():
    """The melittin guard: absence of coverage is not evidence of low coverage."""
    dedup = {"records": [_rec(97.1, None)]}
    summary = gr.apply_coverage_demotion(dedup)
    assert summary["demoted"] == 0
    assert summary["coverage_not_recorded"] == 1
    assert dedup["records"][0]["goi_class"] == "synteny_hull_rescue"


def test_full_length_match_is_kept():
    dedup = {"records": [_rec(97.1, 0.98)]}
    assert gr.apply_coverage_demotion(dedup)["demoted"] == 0
    assert dedup["records"][0]["goi_class"] == "synteny_hull_rescue"


def test_low_identity_is_out_of_scope():
    """A 30% call over 10% of the query is already handled by the identity floors."""
    dedup = {"records": [_rec(30.0, 0.10)]}
    assert gr.apply_coverage_demotion(dedup)["demoted"] == 0


def test_low_confidence_calls_are_untouched():
    dedup = {"records": [_rec(99.0, 0.05, confidence="LOW")]}
    assert gr.apply_coverage_demotion(dedup)["demoted"] == 0


def test_thresholds_are_honoured():
    dedup = {"records": [_rec(60.0, 0.30)]}
    assert gr.apply_coverage_demotion(dedup, min_identity=70.0)["demoted"] == 0
    dedup = {"records": [_rec(60.0, 0.30)]}
    assert gr.apply_coverage_demotion(dedup, max_qcov=0.20)["demoted"] == 0


# ── the real melittin case, reconstructed ───────────────────────────────────

def test_melittin_lrz_case():
    """Job 5757307's three artefacts go; job 5757853's true ortholog stays.

    The artefacts were 32 bp (~10 aa) models of a flanking gene against a 70 aa
    query -> coverage ~0.15. The real call is 135 aa -> coverage ~0.97.
    """
    dedup = {"records": [
        _rec(100.0, 0.15), _rec(100.0, 0.15), _rec(100.0, 0.15),
        _rec(97.1, 0.97, goi_class="synteny_hull_rescue"),
    ]}
    summary = gr.apply_coverage_demotion(dedup)
    assert summary["demoted"] == 3
    survivors = [r for r in dedup["records"]
                 if r["goi_class"] != "identity_coverage_decoupled"]
    assert len(survivors) == 1
    assert survivors[0]["identity"] == "97.1"


# ── the rescue must emit the coverage the rule depends on ───────────────────

def test_hull_rescue_emits_query_coverage():
    mrna = ["c1", "miniprot", "mRNA", "1", "406", ".", "+", ".", "Identity=0.971"]
    cds = [["c1", "miniprot", "CDS", "1", "406", ".", "+", "0"]]
    rows = hull._build_hull_gff_rows(mrna, cds, "c1", 1000, "P01501", "GOI_x",
                                     "HIGH", query_coverage=0.97)
    assert any("QueryCoverage=0.970" in r for r in rows)


def test_hull_rescue_omits_attribute_when_unknown():
    mrna = ["c1", "miniprot", "mRNA", "1", "406", ".", "+", ".", "Identity=0.971"]
    cds = [["c1", "miniprot", "CDS", "1", "406", ".", "+", "0"]]
    rows = hull._build_hull_gff_rows(mrna, cds, "c1", 1000, "P01501", "GOI_x", "HIGH")
    assert not any("QueryCoverage" in r for r in rows)


# ── F4 control arm: --disable_wavefront ─────────────────────────────────────

isr = _load("iterative_search_runner", "iterative_search_runner.py")


class _Args:
    def __init__(self, disable_wavefront=False, rank_wave_binning=True):
        self.disable_wavefront = disable_wavefront
        self.rank_wave_binning = rank_wave_binning


def _entries(dists):
    return [{"name": "g%d" % i, "dist": d} for i, d in enumerate(dists)]


def test_disable_wavefront_gives_one_parallel_wave():
    """The F4 control: no genome is ever searched with a relative's ortholog."""
    entries = _entries([0.2, 0.7, 0.7, 1.0])
    waves = isr.define_waves(entries, _Args(disable_wavefront=True))
    assert len(waves) == 1
    assert len(waves[0]) == 4


def test_disable_wavefront_overrides_rank_binning():
    entries = _entries([0.2, 0.7, 0.7, 1.0])
    both = isr.define_waves(entries, _Args(disable_wavefront=True, rank_wave_binning=True))
    assert len(both) == 1


def test_default_still_grades_close_to_far():
    """The observed LRZ melittin set: 4 targets at three taxonomic ranks."""
    entries = _entries([0.2, 0.7, 0.7, 1.0])
    waves = isr.define_waves(entries, _Args())
    assert len(waves) > 1
    assert len(waves[0]) == 1          # closest genome searched serially
    assert sum(len(w) for w in waves) == 4


def test_every_genome_appears_exactly_once_in_all_modes():
    entries = _entries([0.0, 0.1, 0.3, 0.6, 0.9, 1.0])
    for args in (_Args(), _Args(rank_wave_binning=False), _Args(disable_wavefront=True)):
        waves = isr.define_waves(entries, args)
        names = [g["name"] for w in waves for g in w]
        assert sorted(names) == sorted(g["name"] for g in entries)


# ── rescue models must reach the §1m ownership check (they need a protein) ──

def _cds(start, end, strand="+", phase="0"):
    return ["c1", "miniprot", "CDS", str(start), str(end), ".", strand, phase, ""]


def test_translate_plus_strand():
    assert hull._translate_model("ATGGCTTAA", [_cds(1, 9)]) == "MA"


def test_translate_minus_strand():
    # reverse complement of ATGGCTTAA
    assert hull._translate_model("TTAAGCCAT", [_cds(1, 9, strand="-")]) == "MA"


def test_translate_honours_phase():
    assert hull._translate_model("GGATGGCTTAA", [_cds(1, 11, phase="2")]) == "MA"


def test_translate_joins_exons_in_transcript_order():
    assert hull._translate_model("ATGGCTTAA", [_cds(1, 3), _cds(4, 9)]) == "MA"


def test_translate_minus_strand_orders_exons_descending():
    """miniprot reports ascending coords even on the minus strand."""
    got = hull._translate_model("TTAAGCCAT", [_cds(1, 6, strand="-"), _cds(7, 9, strand="-")])
    assert got == "MA"


def test_translate_trims_partial_codon_and_internal_stops():
    out = hull._translate_model("ATGTAAGCTA", [_cds(1, 10)])
    assert out.startswith("M")
    assert "*" not in out


def test_translate_empty_input():
    assert hull._translate_model("ACGT", []) == ""


# ── easy mode must survive UniProt's strain-qualified organism names ────────

fhg = _load("fetch_home_genome", "fetch_home_genome.py")


def test_strain_qualifier_is_stripped():
    """D6VTK4 (yeast STE2) — the exact string that broke case_ste2 on LRZ."""
    got = fhg.normalize_species_name("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)")
    assert got == "Saccharomyces cerevisiae"


def test_isolate_and_serotype_qualifiers():
    assert fhg.normalize_species_name("Foo bar (isolate XYZ)") == "Foo bar"
    assert fhg.normalize_species_name("Foo bar (serotype 2)") == "Foo bar"


def test_plain_binomial_is_untouched():
    for name in ("Apis mellifera", "Homo sapiens", "Bombus terrestris"):
        assert fhg.normalize_species_name(name) == name


def test_unrelated_parenthetical_is_kept():
    """Only known qualifier keywords are stripped, so nothing else is mangled."""
    assert fhg.normalize_species_name("Foo bar (nickname)") == "Foo bar (nickname)"


def test_empty_input_is_safe():
    assert fhg.normalize_species_name("") == ""
    assert fhg.normalize_species_name(None) is None
