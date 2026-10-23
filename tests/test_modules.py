"""Preset modules (2026-10-15, PR-1 of the preset plan; one rule per key since the profile
handoff of 2026-10-20, stage 1): features, profile, rules, constraints, resolution.

* query_features.py   features of the query sequence alone
* profile_hits.py     per-locus union coverage, loci, secondary_full_length, shared_domain_only
* home_context.py     family members next to the home gene (the label's family test)
* select_modules.py   triggers, the design profile, one rule per key (owners, formula, cap),
                      the derivation record, the constraint checks, the --modules value
* resolve_effective_params.py  --modules legacy reproduces the old preset exactly; auto,
                      lists and deltas apply the module rules below the preset layer, and a
                      value the user set always wins
"""
from __future__ import annotations

import itertools
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import auto_select_preset as asp  # noqa: E402
import generate_report  # noqa: E402
import home_context as hc  # noqa: E402
import profile_hits  # noqa: E402
import query_features as qf  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
import select_modules as sm  # noqa: E402
import synvoy_specificity  # noqa: E402

# UniProt P01501 (melittin precursor) and Q16553 (LY6E), public sequences.
MELITTIN = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"
LY6E = ("MKIFLPVLLAALLGVERASSLMCFSCLNQKSNLYCLKPTICSDQDNYCVTVSASAGIGNLVTFGHSLSKTCSPACPIPEGVNVG"
        "VASMGISCCQSFLCNFSAADGGLRASVTLLGAGLLLSLLPALLRFGP")

DEFAULTS = {
    "sw_min_score": 20, "sw_min_identity": 10.0, "min_hit_identity": 10, "min_hit_length": 10,
    "classify_fragment_max_qcov": 0.4, "classify_complete_min_qcov": 0.7,
    "expand_goi_similar": False, "goi_family_tokens": "", "goi_relative_evalue": 1e-5,
    "adaptive_score_floor_frac": 0.3, "adaptive_score_floor_abs": 0.03,
    "adaptive_max_regions": 6, "adaptive_unique_gene_floor": 3,
}


def _profile(query_len=300, n_hits=4, n_strong=1, n_loci=1, ratio=None, copies=1.0,
             secondary_qcov=None, seq=None):
    loci = [{"chrom": "c1", "start": 1, "end": 900, "bit": 400.0, "qcov": 1.0}]
    if secondary_qcov is not None:
        loci.append({"chrom": "c2", "start": 1, "end": 900, "bit": 300.0, "qcov": secondary_qcov})
    return {
        "query_len": query_len, "n_hits": n_hits, "n_hits_raw": n_hits,
        "n_independent_loci": n_loci, "n_strong_loci": n_strong, "top_bit": 400.0,
        "second_bit": 300.0 if n_strong >= 2 else None, "bit_ratio": ratio,
        "primary_copies": copies, "loci": loci,
        "secondary_full_length": bool(n_strong >= 2 and (secondary_qcov or 0) >= 0.5),
        "shared_domain_only": bool(n_strong >= 2 and (secondary_qcov or 0) < 0.5),
        "query_features": qf.query_features(seq or ("MKT" + "AILVG" * (query_len // 5))[:query_len]),
    }


# ============================================================ query features
def test_ly6e_has_a_signal_peptide_a_gpi_tail_and_ten_cysteines():
    f = qf.query_features(LY6E)
    assert (f["sp"], f["sp_len"]) == (True, 20)            # UniProt: signal 1-20
    assert f["gpi"] and 25 <= f["gpi_len"] <= 35            # UniProt: omega at 101 of 131
    assert f["n_cys"] == 10 and f["cys_max90"] == 10
    assert LY6E.endswith("GGLRASVTLLGAGLLLSLLPALLRFGP")


def test_melittin_precursor_has_its_signal_peptide_and_low_complexity():
    f = qf.query_features(MELITTIN)
    assert (f["sp"], f["sp_len"]) == (True, 21)            # cleaved after ...YIYA
    assert not f["gpi"] and f["n_cys"] == 0
    assert f["seg_frac"] >= 0.10                            # the PEPEP pro-region


def test_low_complexity_and_repeats_are_measured():
    polyq = "MKT" + "Q" * 30 + "LVNVALVFMVVYISYIYAWERTYIPASDFGHKLCVNM"
    assert qf.query_features(polyq)["seg_frac"] >= 0.10
    assert qf.query_features("MKTAPEPAPEPAPEPAPEPAPEPLVNW")["repeat_frac"] >= 0.5
    plain = qf.query_features("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    assert not plain["sp"] and not plain["gpi"]


def test_query_features_are_deterministic_and_rounded():
    assert qf.query_features(LY6E) == qf.query_features(LY6E)
    assert json.dumps(qf.query_features(LY6E), sort_keys=True) == json.dumps(
        qf.query_features(LY6E.lower()), sort_keys=True)


# ============================================================ profile_hits
def _hit(chrom, qs, qe, ts, te, bits, ev=1e-30, pid=90.0):
    return {"chrom": chrom, "qstart": qs, "qend": qe, "tstart": ts, "tend": te,
            "pident": pid, "bits": bits, "evalue": ev}


def test_locus_coverage_is_a_union_not_a_sum():
    hits = [_hit("c1", 1, 60, 1000, 1180, 100.0), _hit("c1", 40, 100, 2000, 2183, 90.0)]
    p = profile_hits.build_profile(hits, 100, 50_000, 0.5)
    assert p["loci"][0]["qcov"] == 1.0                     # 1-100, not 60 + 61
    assert p["primary_copies"] == 1.21                     # the copy number keeps the sum


def test_profile_keeps_every_old_key_and_adds_the_new_ones():
    hits = [_hit("c1", 1, 100, 1000, 1300, 200.0), _hit("c2", 1, 40, 5000, 5120, 150.0)]
    p = profile_hits.build_profile(hits, 100, 50_000, 0.5, query_seq=MELITTIN)
    old = {"query_len", "n_hits", "n_hits_raw", "n_independent_loci", "n_strong_loci", "top_bit",
           "second_bit", "bit_ratio", "primary_copies", "identity_min", "identity_q1",
           "identity_median", "identity_q3", "identity_max", "identity_iqr"}
    assert old <= set(p)
    assert (p["n_strong_loci"], p["secondary_full_length"], p["shared_domain_only"]) == (2, False, True)
    assert p["query_features"]["sp_len"] == 21
    empty = profile_hits.build_profile([], 70, 50_000, 0.5, query_seq=MELITTIN)
    assert empty["loci"] == [] and empty["secondary_full_length"] is False
    assert empty["query_features"]["query_len"] == 70


# ============================================================ the rules
def _preset_default_keys():
    return set(re.findall(r"'([a-z_]+)'", re.search(
        r"def presetDefaultKeys\(\) \{\s*return \[(.*?)\]", (ROOT / "main.nf").read_text(), re.S).group(1)))


def test_every_rule_has_owners_a_default_a_settings_reader_and_no_params_reader():
    """A rule reaches the search only through the settings map: main.nf feeds its key
    (presetDefaultKeys) and a process reads ${settings.x}. A ${params.x} reader in a process
    would win over the rule without a word (main.nf's start-up check of the value a user
    passes is not a reader)."""
    keys = _preset_default_keys()
    readers = "".join(p.read_text() for d in ("modules", "subworkflows") for p in (ROOT / d).glob("*.nf"))
    ids = {m.id for m in sm.MODULES}
    for r in sm.RULES:
        assert r.owners and set(r.owners) <= ids, r.key
        assert r.key in keys, r.key
        assert f"settings.{r.key}" in readers, r.key
        assert not re.search(rf"params\.{r.key}\b", readers), r.key
        assert r.key not in sm.FORBIDDEN_FOR_MODULES, r.key
        assert r.evidence in sm.EVIDENCE and r.text and r.doc, r.key
    assert list(sm.RULES_BY_KEY) == [r.key for r in sm.RULES]
    shipped = rep.load_shipped_defaults(str(ROOT / "nextflow.config"), sm.RULES_BY_KEY)
    assert set(shipped) == set(sm.RULES_BY_KEY)                    # each has a shipped default


# Features just inside and just outside each trigger that owns a rule (the trigger constants
# of select_modules.py).
EDGES = {
    "A1": ({"primary_copies": 1.5}, {"primary_copies": 1.4}),
    "A2": ({"n_strong_loci": 2, "secondary_full_length": True}, {"n_strong_loci": 1, "secondary_full_length": True}),
    "B1": ({"query_len": 99}, {"query_len": 100}),
    "B2": ({"seg_frac": 0.10, "repeat_frac": 0.0}, {"seg_frac": 0.099, "repeat_frac": 0.249}),
    "B3": ({"sp": True, "sp_len": 20}, {"sp": False, "sp_len": 0, "gpi": False}),
}


@pytest.mark.parametrize("key", [r.key for r in sm.RULES])
def test_each_rule_at_the_edge_of_its_trigger(key):
    rule = sm.RULES_BY_KEY[key]
    for mid in rule.owners:
        inside, outside = EDGES[mid]
        res = rep.resolve(FULL, auto_preset=None, modules_spec="auto", features=inside)
        want = rule.formula(sm.design_profile(inside), frozenset({mid}))
        assert res["settings"][key] == want and res["source"][key] == f"module:{mid}", (key, mid)
        res = rep.resolve(FULL, auto_preset=None, modules_spec="auto", features=outside)
        assert res["source"][key] == "default" and res["settings"][key] == FULL[key], (key, mid)


@pytest.mark.parametrize("copies,want,capped", [(9.0, 10, False), (10.0, 10, True), (31.0, 10, True)])
def test_the_regions_rule_at_its_cap(copies, want, capped):
    settings, _source, derivation = sm.apply_rules({"A1"}, sm.design_profile({"primary_copies": copies}))
    d = derivation["goi_candidate_regions_per_query"]
    assert settings["goi_candidate_regions_per_query"] == want
    assert ("uncapped" in d) is capped and d.get("uncapped", want) == int(copies) + 1


# Settings keys that no rule writes: no preset or module has set them since 2026-10-09/15, so
# they are params read through the settings map. The profile handoff (§7.4) wants the map to
# hold the rule keys only; moving these readers to ${params.x} changes the report's
# run_settings, so it waits for stage 4. This set may only shrink.
SETTINGS_KEYS_WITHOUT_A_RULE = {"sw_min_score", "sw_min_identity", "min_hit_identity", "min_hit_length",
                                "expand_goi_similar", "goi_family_tokens", "adaptive_unique_gene_floor"}


def test_the_settings_map_holds_the_rule_keys_and_known_leftovers_only():
    assert _preset_default_keys() == set(sm.RULES_BY_KEY) | SETTINGS_KEYS_WITHOUT_A_RULE
    assert not SETTINGS_KEYS_WITHOUT_A_RULE & set(sm.RULES_BY_KEY)


def test_the_rules_of_stage_1_are_the_fifteen_module_keys_of_2026_10_15():
    """Stage 1 of the profile handoff moved the keys into rules without changing a value."""
    assert {m.id: sm.rule_keys(m.id) for m in sm.MODULES} == {
        "A1": ["goi_relative_evalue", "goi_candidate_regions_per_query", "report_copies_in_slot"],
        "A2": ["adaptive_score_floor_frac", "adaptive_score_floor_abs", "adaptive_max_regions"],
        "A3": [],
        "B1": ["classify_complete_min_qcov", "classify_fragment_max_qcov", "window_tblastn_word_size",
               "window_tblastn_comp_stats", "window_mmseqs_min_length", "miniprot_short_orf_max_len"],
        "B2": ["goi_relative_statistic", "report_low_complexity"],
        "B3": ["core_completeness"],
        "B4": [], "C1": [], "C2": [],
    }
    assert {k: r.cap for k, r in sm.RULES_BY_KEY.items() if r.cap is not None} == {
        "goi_candidate_regions_per_query": 10}


def test_a_rule_applies_when_any_owner_is_applied_and_its_formula_alone_decides():
    rule = sm.Rule("goi_candidate_regions_per_query", ("C1", "A1"),
                   lambda p, on: p["c"] + (1 if "C1" in on else 0), "c, + 1 with C1",
                   inputs=("c",), cap=5, doc="a test rule")
    assert sm.apply_rules(set(), {"c": 3}, [rule]) == ({}, {}, {})
    assert sm.apply_rules({"A1", "B2"}, {"c": 3}, [rule])[0] == {"goi_candidate_regions_per_query": 3}
    settings, source, derivation = sm.apply_rules({"C1", "A1"}, {"c": 3}, [rule])
    assert settings == {"goi_candidate_regions_per_query": 4}
    assert source == {"goi_candidate_regions_per_query": "module:A1+C1"}
    d = derivation["goi_candidate_regions_per_query"]
    assert d["line"] == "goi_candidate_regions_per_query 4 = c, + 1 with C1, c = 3 (A1+C1)"
    assert (d["owners"], d["inputs"], d["evidence"]) == (["A1", "C1"], {"c": 3}, "inferred")
    assert "cap" not in d
    settings, _source, derivation = sm.apply_rules({"A1"}, {"c": 9}, [rule])
    d = derivation["goi_candidate_regions_per_query"]
    assert settings["goi_candidate_regions_per_query"] == 5 and (d["cap"], d["uncapped"]) == (5, 9)
    assert d["line"].endswith(", capped at 5 (wanted 9) (A1)")


def _rule(key, owners, **kw):
    return sm.Rule(key, owners, lambda p, on: 1, "1", **dict({"doc": "test"}, **kw))


@pytest.mark.parametrize("rules,why", [
    ([_rule("goi_relative_evalue", ("A1",)), _rule("goi_relative_evalue", ("B1",))], "two rules"),
    ([_rule("decoy_mode", ("B2",))], "forbidden"),
    ([_rule("n_flanking_genes", ("A1",))], "forbidden"),
    ([_rule("tandem_copy_parent", ("A1",))], "forbidden"),
    ([_rule("goi_relative_evalue", ("Z9",))], "not module ids"),
    ([_rule("goi_relative_evalue", ())], "not module ids"),
    ([_rule("goi_relative_evalue", ("A1",), inputs=("Lmystery",))], "does not have"),
    ([_rule("goi_relative_evalue", ("A1",), evidence="hunch")], "evidence"),
    ([_rule("goi_relative_evalue", ("A1",), doc="")], "doc"),
])
def test_the_registry_refuses_a_bad_rule(rules, why):
    with pytest.raises(ValueError, match=why):
        sm.validate_registry(rules, sm.MODULES)


class _Reads(dict):
    """A profile that records which quantities a formula reads."""

    def __init__(self, values):
        super().__init__(values)
        self.read = set()

    def __getitem__(self, k):
        self.read.add(k)
        return super().__getitem__(k)

    def get(self, k, default=None):
        self.read.add(k)
        return super().get(k, default)


@pytest.mark.parametrize("c", [None, 1, 3, 40])
def test_every_formula_reads_only_its_declared_profile_inputs(c):
    values = {"c": c, "R": 2, "Fam": 3, "Lq": 60, "Lx": 30, "Lpiece": 30, "Lcore": 40, "nx": 2, "n": 10}
    assert set(values) == set(sm.PROFILE_SYMBOLS)
    for r in sm.RULES:
        p = _Reads(values)
        r.formula(p, frozenset(r.owners))
        assert p.read <= set(r.inputs), (r.key, p.read)


# ============================================================ the profile
def test_the_design_profile_of_a_secreted_gene_in_an_array():
    hp = _profile(query_len=len(LY6E), copies=2.6, n_loci=3, seq=LY6E)
    home = {"relatives_nearby": [
        {"gene_id": "gene-a", "rank": -4, "evalue_seg_no": 1e-9},
        {"gene_id": "gene-b", "rank": 2, "evalue_seg_no": 1e-6}]}
    goi = {"method": "gff_annotation", "exons": [
        {"id": "GOI_q|exon_1", "length_aa": 35}, {"id": "GOI_q|exon_2", "length_aa": 96},
        {"id": "GOI_copy_1", "length_aa": 300}]}
    f = sm.flat_features(hp, goi, home)
    p = sm.design_profile(f, goi, "10")
    assert f["sp_len"] == 20 and f["gpi"]
    assert p == {"c": 3, "R": 4, "Fam": 3, "Lq": 131, "Lx": 96, "Lpiece": 96,
                 "Lcore": 131 - 20 - f["gpi_len"], "nx": 2, "n": 10}
    # a short query: the longest piece is the query
    short = sm.design_profile(dict(f, query_len=70, sp=False, gpi=False), goi, 10)
    assert (short["Lpiece"], short["Lcore"]) == (70, 70)


def test_the_design_profile_without_its_inputs():
    assert sm.design_profile({}, None, None) == {
        "c": None, "R": 0, "Fam": None, "Lq": None, "Lx": None, "Lpiece": None, "Lcore": None,
        "nx": None, "n": None}
    f = sm.flat_features(_profile(query_len=80), None, None)
    p = sm.design_profile(f, {"method": "tandem_duplication", "exons": [{"id": "e", "length_aa": 40}]})
    assert (p["Lx"], p["Lpiece"], p["nx"], p["c"]) == (None, 80, None, 1)
    assert sm.design_profile(f, {"method": "gff_annotation", "exons": [{"id": "GOI_copy_1",
                                                                       "length_aa": 40}]})["Lx"] is None


# ============================================================ derivation and constraints
FULL = dict(DEFAULTS, goi_candidate_regions_per_query=2, report_copies_in_slot=False,
            window_tblastn_word_size=3, window_tblastn_comp_stats=2, window_mmseqs_min_length=30,
            miniprot_short_orf_max_len=60, goi_relative_statistic="sw", report_low_complexity=False,
            core_completeness=False)


def test_the_derivation_says_what_set_each_key_and_from_what():
    f = sm.flat_features(_profile(query_len=50, copies=3.0), None, None)
    res = rep.resolve(dict(FULL, goi_relative_evalue="1e-6"), auto_preset=None, modules_spec="auto",
                      features=f, user_set={"goi_relative_evalue": "1e-6"}, shipped={"goi_relative_evalue": "1e-5"})
    d = res["derivation"]
    assert d["goi_candidate_regions_per_query"]["line"] == "goi_candidate_regions_per_query 4 = max(2, c + 1), c = 3 (A1)"
    assert (d["goi_candidate_regions_per_query"]["default"], d["goi_candidate_regions_per_query"]["inputs"]) == (2, {"c": 3})
    assert d["classify_complete_min_qcov"]["line"] == "classify_complete_min_qcov 0.5 (B1)"
    assert d["goi_relative_evalue"] == {"value": "1e-6", "default": 1e-5, "rule": "user_set", "wanted": 1e-3,
                                        "line": "goi_relative_evalue 1e-6 (set by the user)"}
    assert set(d) == {k for k, s in res["source"].items() if s != "default"}
    assert res["profile"]["c"] == 3 and res["profile"]["Lq"] == 50
    legacy = rep.resolve(FULL, auto_preset={"preset": "preset_short_peptide"}, modules_spec="legacy", features=f)
    assert {k: v["rule"] for k, v in legacy["derivation"].items()} == {
        "classify_complete_min_qcov": "preset:preset_short_peptide",
        "classify_fragment_max_qcov": "preset:preset_short_peptide"}


def _constraint(res, cid):
    return {c["id"]: c for c in res["constraints"]}[cid]


def test_a_binding_cap_is_named_and_is_no_violation():
    f = sm.flat_features(_profile(copies=40.0), None, None)
    res = rep.resolve(FULL, auto_preset=None, modules_spec="auto", features=f)
    assert res["settings"]["goi_candidate_regions_per_query"] == 10
    assert res["derivation"]["goi_candidate_regions_per_query"]["uncapped"] == 41
    c3 = _constraint(res, "C3")
    assert c3["ok"] and "the cap of 10 binds" in c3["detail"]
    assert _constraint(res, "C11")["ok"]
    legacy = rep.resolve(FULL, auto_preset={"preset": "preset_tandem_family"}, modules_spec="legacy", features=f)
    assert not _constraint(legacy, "C3")["applies"] and not _constraint(legacy, "C11")["applies"]


def test_a_value_the_user_set_against_a_rule_is_a_warning():
    f = sm.flat_features(_profile(query_len=50, copies=3.0), None, None)
    res = rep.resolve(dict(FULL, classify_fragment_max_qcov=0.6, goi_candidate_regions_per_query=2),
                      auto_preset=None, modules_spec="auto", features=f,
                      user_set={"classify_fragment_max_qcov": 0.6, "goi_candidate_regions_per_query": 2})
    for cid in ("C3", "C12"):
        c = _constraint(res, cid)
        assert c["applies"] and not c["ok"] and c["on_violation"] == "warn", cid
    assert res["settings"]["classify_fragment_max_qcov"] == 0.6
    odd = rep.resolve(dict(FULL, classify_fragment_max_qcov="low"), auto_preset=None, modules_spec="auto",
                      features=f, user_set={"classify_fragment_max_qcov": "low"})
    c12 = _constraint(odd, "C12")
    assert not c12["ok"] and c12["on_violation"] == "warn" and c12["detail"].startswith("not checked")


def test_a_violation_the_rules_make_alone_stops_the_run(monkeypatch, tmp_path, capsys):
    bug = sm.Rule("classify_fragment_max_qcov", ("B1",), lambda p, on: 0.9, "0.9", doc="a bug")
    monkeypatch.setattr(sm, "RULES", [bug if r.key == bug.key else r for r in sm.RULES])
    f = sm.flat_features(_profile(query_len=50), None, None)
    c12 = _constraint(rep.resolve(FULL, auto_preset=None, modules_spec="auto", features=f), "C12")
    assert not c12["ok"] and c12["on_violation"] == "stop"
    over = sm.check_constraints({"goi_candidate_regions_per_query": 12},
                                {"goi_candidate_regions_per_query": "module:A1"}, {"c": 3}, {"A1"})
    assert {c["id"]: (c.get("ok"), c.get("on_violation")) for c in over}["C11"] == (False, "stop")
    (tmp_path / "defaults.json").write_text(json.dumps(FULL))
    (tmp_path / "hp.json").write_text(json.dumps(_profile(query_len=50)))
    monkeypatch.setattr(sys, "argv", ["resolve_effective_params.py", "--hit_profile", str(tmp_path / "hp.json"),
                                      "--modules", "auto", "--defaults", str(tmp_path / "defaults.json"),
                                      "--output", str(tmp_path / "eff.json")])
    with pytest.raises(SystemExit) as stop:
        rep.main()
    assert stop.value.code == 2
    assert "constraint C12 violated" in capsys.readouterr().err


# ============================================================ forbidden keys
def test_switches_n_and_the_mechanism_are_forbidden_to_the_rules():
    """C10 of the profile handoff: no rule writes a switch, n_flanking_genes, decoy_mode,
    expand_goi_similar or the module mechanism; every forbidden key is a real param."""
    assert len(sm.SWITCHES) == 27
    declared = rep.load_shipped_defaults(str(ROOT / "nextflow.config"), sm.FORBIDDEN_FOR_MODULES)
    assert set(declared) == set(sm.FORBIDDEN_FOR_MODULES)
    assert {"n_flanking_genes", "decoy_mode", "expand_goi_similar", "modules", "tandem_copy_parent",
            "goi_genome_wide_path", "sw_method"} <= sm.FORBIDDEN_FOR_MODULES
    assert not set(sm.RULES_BY_KEY) & sm.FORBIDDEN_FOR_MODULES


def test_merge_does_not_depend_on_the_order_of_the_modules():
    f = sm.flat_features(_profile(query_len=50, copies=3.0, n_strong=3, secondary_qcov=0.9), None, None)
    chosen = [m for m in sm.MODULES if m.id in ("A1", "A2", "B1", "B2", "C1")]
    first = sm.merge(chosen, f)
    for perm in itertools.permutations(chosen):
        assert sm.merge(list(perm), f) == first


# Keys the old presets never had (added by later PRs of the plan); --modules legacy keeps
# applying PRESET_OVERRIDES alone, so a legacy run never sees them.
KEYS_AFTER_V1 = {"goi_candidate_regions_per_query", "report_copies_in_slot",
                 "window_tblastn_word_size", "window_tblastn_comp_stats",
                 "window_mmseqs_min_length", "miniprot_short_orf_max_len"}


def test_the_old_presets_are_module_sets_with_the_same_values():
    by_id = {m.id: m for m in sm.MODULES}
    f = sm.flat_features(_profile(), None, None)
    for preset, ids in sm.LEGACY_ALIAS.items():
        settings, _src = sm.merge([by_id[i] for i in ids], f)
        old = rep.PRESET_OVERRIDES[preset]
        assert {k: settings[k] for k in old} == old, preset
        assert set(settings) - set(old) <= KEYS_AFTER_V1, preset
    assert set(sm.LEGACY_ALIAS) == set(rep.PRESET_OVERRIDES)


def test_groovy_module_list_matches_the_registry():
    text = (ROOT / "main.nf").read_text()
    block = re.search(r"def presetModules\(\) \{\s*return \[(.*?)\]", text, re.S).group(1)
    groovy = dict(re.findall(r"([A-C]\d): '([a-z_]+)'", block))
    assert groovy == {m.id: m.name for m in sm.MODULES}


@pytest.mark.parametrize("spec,expected", [
    ("legacy", ("legacy", set(), set())),
    ("", ("legacy", set(), set())),
    ("auto", ("auto", set(), set())),
    ("none", ("none", set(), set())),
    ("short_query,A1", ("list", {"B1", "A1"}, set())),
    ("auto,-dispersed_family,+cys_framework", ("auto", {"B4"}, {"A2"})),
    ("none,+short_query", ("list", {"B1"}, set())),
])
def test_modules_spec(spec, expected):
    assert sm.parse_spec(spec) == expected


@pytest.mark.parametrize("spec", [True, "true", "bogus", "legacy,short_query", "none,-B1", "-B1"])
def test_bad_modules_specs_are_refused(spec):
    with pytest.raises(ValueError):
        sm.parse_spec(spec)


# ============================================================ triggers
def test_a_family_member_next_to_the_gene_fires_local_family_under_auto_only():
    """LY6E-like: one copy at home, SLURP1-like relative four genes to the left."""
    hp = _profile(query_len=131, seq=LY6E)
    home = {"relatives_nearby": [{"gene_id": "gene-X", "rank": -4, "distance_bp": 276000,
                                  "evalue_seg_no": 2.8e-7, "evalue_seg_yes": 3e-7}],
            "relative_adjacent": False, "relative_proximal": True, "statistic": "blastp"}
    f = sm.flat_features(hp, None, home)
    fired = {d["id"] for d in sm.evaluate(f)}
    assert {"A1", "B3", "B4"} <= fired and "B1" not in fired
    auto = rep.resolve(DEFAULTS, auto_preset={"preset": asp.select_preset(hp)[0]},
                       modules_spec="auto", features=f)
    assert auto["settings"]["goi_relative_evalue"] == 1e-3
    assert auto["source"]["goi_relative_evalue"] == "module:A1"
    legacy = rep.resolve(DEFAULTS, auto_preset={"preset": asp.select_preset(hp)[0]},
                         modules_spec="legacy", features=f)
    assert legacy["preset_applied"] == "preset_single_copy"
    assert legacy["settings"] == DEFAULTS
    assert "A1" in {d["id"] for d in legacy["modules_would_fire"]}


def test_shared_domain_and_dispersed_family_are_told_apart():
    luciferase = sm.flat_features(_profile(n_strong=4, n_loci=50, ratio=1.3, secondary_qcov=0.3), None, None)
    fired = {d["id"] for d in sm.evaluate(luciferase)}
    assert "A3" in fired and "A2" not in fired
    tp53 = sm.flat_features(_profile(n_strong=3, n_loci=3, ratio=1.2, secondary_qcov=0.8), None, None)
    fired = {d["id"] for d in sm.evaluate(tp53)}
    assert "A2" in fired and "A3" not in fired


def test_exon_structure_counts_real_exons_only():
    exons = [{"id": "GOI_q|exon_1"}, {"id": "GOI_q|exon_2"}, {"id": "GOI_copy_1"}]
    f = sm.flat_features(None, {"method": "gff_annotation", "exons": exons}, None)
    assert f["goi_real_exons"] == 2 and sm._t_exon_structure(f)[0]
    f = sm.flat_features(None, {"method": "tandem_duplication", "exons": exons}, None)
    assert f["goi_real_exons"] is None and not sm._t_exon_structure(f)[0]


def test_contig_end_needs_a_long_scaffold():
    near = {"dist_to_scaffold_end": 12_000, "scaffold_len": 3_000_000}
    assert sm._t_home_context(sm.flat_features(None, None, near))[0]
    assert not sm._t_home_context(sm.flat_features(None, None, dict(near, scaffold_len=200_000)))[0]
    assert not sm._t_home_context(sm.flat_features(None, None, dict(near, dist_to_scaffold_end=90_000)))[0]


# ============================================================ resolution
GRID = [
    dict(n_hits=0, query_len=300),
    dict(query_len=50), dict(query_len=99), dict(query_len=100), dict(query_len=101),
    dict(n_strong=2, ratio=1.2, secondary_qcov=0.9), dict(n_strong=2, ratio=2.0, secondary_qcov=0.9),
    dict(copies=1.4), dict(copies=1.5), dict(copies=2.0), dict(),
]


@pytest.mark.parametrize("kw", GRID)
def test_legacy_reproduces_the_old_preset_exactly(kw):
    """Same settings, preset and per-key source as the resolver before the modules."""
    hp = _profile(**kw)
    preset = asp.select_preset(hp)[0]
    res = rep.resolve(DEFAULTS, auto_preset={"preset": preset}, modules_spec="legacy",
                      features=sm.flat_features(hp, None, None))
    expected = dict(DEFAULTS, **rep.PRESET_OVERRIDES[preset])
    assert res["preset_applied"] == preset and res["settings"] == expected
    for k in DEFAULTS:
        want = f"preset:{preset}" if k in rep.PRESET_OVERRIDES[preset] else "default"
        assert res["source"][k] == want, (k, res["source"][k])
    assert [d["id"] for d in res["modules_applied"]] == list(sm.LEGACY_ALIAS[preset])


def test_none_applies_nothing_and_a_list_forces_its_modules():
    f = sm.flat_features(_profile(), None, None)
    none = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"}, modules_spec="none",
                       features=f)
    assert none["settings"] == DEFAULTS and none["preset_applied"] is None
    forced = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="short_query", features=f)
    assert forced["settings"]["classify_complete_min_qcov"] == 0.5
    assert forced["modules_applied"][0]["reason"].startswith("named in --modules")


def test_preset_override_sits_above_the_modules():
    f = sm.flat_features(_profile(copies=2.0), None, None)
    res = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="auto", features=f,
                      preset_override="preset_short_peptide")
    assert res["settings"]["goi_relative_evalue"] == 1e-3                 # A1 fired
    assert res["source"]["classify_fragment_max_qcov"] == "preset:preset_short_peptide"
    assert res["preset_applied"] == "preset_short_peptide"


def test_a_value_the_user_set_survives_a_module_and_is_reported():
    f = sm.flat_features(_profile(copies=2.0), None, None)
    res = rep.resolve(dict(DEFAULTS, goi_relative_evalue="1e-6"), auto_preset=None,
                      modules_spec="auto", features=f, user_set={"goi_relative_evalue": "1e-6"})
    assert res["settings"]["goi_relative_evalue"] == "1e-6"
    assert res["source"]["goi_relative_evalue"] == "user_set"
    assert res["user_set_kept"] == ["goi_relative_evalue"]


def test_a_profile_preset_enters_as_user_set():
    """-profile preset_X sets params at launch: the resolver sees them as user-set values."""
    f = sm.flat_features(_profile(), None, None)
    preset_vals = rep.PRESET_OVERRIDES["preset_paralog_discrimination"]
    res = rep.resolve(dict(DEFAULTS, **preset_vals), auto_preset={"preset": "preset_single_copy"},
                      modules_spec="legacy", features=f, user_set=dict(preset_vals))
    for k, v in preset_vals.items():
        assert res["settings"][k] == v and res["source"][k] == "user_set"


def test_switching_auto_apply_off_applies_nothing_but_records_the_triggers():
    f = sm.flat_features(_profile(query_len=50), None, None)
    res = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"}, modules_spec="auto",
                      features=f, disable_auto_apply=True)
    assert res["settings"] == DEFAULTS and res["modules_applied"] == []
    assert "B1" in {d["id"] for d in res["modules_would_fire"]}


def test_resolver_cli_with_sentinels_and_modules(tmp_path):
    (tmp_path / "defaults.json").write_text(json.dumps(DEFAULTS))
    (tmp_path / "hp.json").write_text(json.dumps(_profile(query_len=50)))
    sentinel = tmp_path / "NO_HOME_CONTEXT"
    sentinel.write_text("sentinel: none\n")
    out = tmp_path / "eff.json"
    (tmp_path / "goi.json").write_text(json.dumps({"method": "gff_annotation", "exons": [
        {"id": "GOI_q|exon_1", "length_aa": 18}, {"id": "GOI_q|exon_2", "length_aa": 32}]}))
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--auto_preset", str(tmp_path / "NO_AUTO_PRESET"), "--hit_profile", str(tmp_path / "hp.json"),
         "--goi_info", str(tmp_path / "goi.json"), "--home_context", str(sentinel),
         "--modules", "auto", "--n_flanking_genes", "10",
         "--defaults", str(tmp_path / "defaults.json"), "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    res = json.loads(out.read_text())
    assert res["modules_mode"] == "auto" and res["settings"]["classify_fragment_max_qcov"] == 0.3
    assert res["features"]["relatives_nearby"] is None
    assert "MODULES: mode auto; applied: B1 short_query" in proc.stdout
    assert "RULE: classify_fragment_max_qcov 0.3 (B1)" in proc.stdout
    assert (res["profile"]["n"], res["profile"]["Lx"], res["profile"]["Lpiece"], res["profile"]["nx"]) == (10, 32, 32, 2)
    assert res["derivation"]["classify_fragment_max_qcov"]["default"] == 0.4
    assert all(c["ok"] for c in res["constraints"] if c["applies"])
    bad = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--modules", "bogus", "--defaults", str(tmp_path / "defaults.json"), "--output", str(out)],
        capture_output=True, text=True)
    assert bad.returncode == 2 and "unknown preset module" in bad.stderr


def test_report_reads_the_run_settings(tmp_path):
    eff = {"modules_mode": "legacy", "modules_applied": [], "modules_would_fire": [{"id": "B3"}],
           "preset_applied": "preset_single_copy", "apply_reason": "x", "settings": {"a": 1},
           "features": {"sp": True}}
    p = tmp_path / "effective_params.json"
    p.write_text(json.dumps(eff))
    rs = generate_report.load_run_settings(str(p))
    assert rs["modules_would_fire"] == [{"id": "B3"}] and rs["features"] == {"sp": True}
    (tmp_path / "NO_EFFECTIVE_PARAMS").write_text("sentinel: none\n")
    assert generate_report.load_run_settings(str(tmp_path / "NO_EFFECTIVE_PARAMS")) is None


# ============================================================ home context
CODON = {"M": "ATG", "K": "AAA", "I": "ATT", "F": "TTT", "L": "CTG", "P": "CCG", "V": "GTG",
         "A": "GCG", "G": "GGC", "E": "GAA", "R": "CGT", "S": "AGC", "C": "TGC", "N": "AAC",
         "Q": "CAG", "Y": "TAT", "T": "ACC", "D": "GAT", "H": "CAT", "W": "TGG"}


def _gene_dna(protein):
    return "".join(CODON[a] for a in protein) + "TAA"


@pytest.mark.skipif(shutil.which("blastp") is None and synvoy_specificity._blastp() is None,
                    reason="blastp not available")
def test_home_context_finds_a_family_member_among_the_neighbours(tmp_path):
    import random
    rnd = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    relative = "".join(c if rnd.random() < 0.6 else rnd.choice(aa) for c in LY6E)
    others = ["M" + "".join(rnd.choice(aa) for _ in range(150)) for _ in range(5)]
    proteins = others[:2] + [relative] + [LY6E] + others[2:]     # relative at rank -1
    spacer = "".join(rnd.choice("ACGT") for _ in range(3000))
    seq, gff, pos = spacer, ["##gff-version 3"], len(spacer)
    for i, prot in enumerate(proteins):
        dna = _gene_dna(prot)
        s, e = pos + 1, pos + len(dna)
        gff += [f"chr1\tt\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-g{i}",
                f"chr1\tt\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-g{i};Parent=gene-g{i}",
                f"chr1\tt\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-g{i};Parent=rna-g{i}"]
        seq += dna + spacer
        pos = len(seq)
    (tmp_path / "g.fna").write_text(">chr1\n" + seq + "\n")
    (tmp_path / "g.gff").write_text("\n".join(gff) + "\n")
    goi_line = [l for l in gff if l.startswith("chr1\tt\tCDS") and "cds-g3" in l][0].split("\t")
    (tmp_path / "goi.bed").write_text(f"chr1\t{int(goi_line[3]) - 1}\t{goi_line[4]}\tGOI_x|exon_1\t100\t+\n")
    shipped = synvoy_specificity.FAMILY_TEST_SEG
    ctx = hc.home_context(LY6E, str(tmp_path / "g.fna"), str(tmp_path / "g.gff"),
                          str(tmp_path / "goi.bed"), 10)
    assert ctx["goi_gene"] == "gene-g3" and ctx["n_neighbours_tested"] == 6
    assert [r["gene_id"] for r in ctx["relatives_nearby"]] == ["gene-g2"]
    assert ctx["relatives_nearby"][0]["rank"] == -1 and ctx["relative_adjacent"] is True
    assert ctx["scaffold_len"] == len(seq) and ctx["home_genome_bp"] == len(seq)
    assert synvoy_specificity.FAMILY_TEST_SEG == shipped       # restored after the test calls
    no_gff = hc.home_context(LY6E, str(tmp_path / "g.fna"), "NO_GFF", str(tmp_path / "goi.bed"), 10)
    assert no_gff["relatives_nearby"] is None and no_gff["dist_to_scaffold_end"] is not None


def test_home_context_without_blastp_reports_no_relatives(tmp_path, monkeypatch):
    (tmp_path / "g.fna").write_text(">chr1\n" + "ACGT" * 3000 + "\n")
    (tmp_path / "g.gff").write_text("chr1\tt\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=r1\n")
    (tmp_path / "goi.bed").write_text("chr1\t99\t300\tGOI_x|exon_1\t100\t+\n")
    monkeypatch.setattr(synvoy_specificity, "_blastp", lambda: None)
    ctx = hc.home_context("MKT", str(tmp_path / "g.fna"), str(tmp_path / "g.gff"),
                          str(tmp_path / "goi.bed"), 10)
    assert ctx["statistic"] == "fallback" and ctx["relatives_nearby"] is None
