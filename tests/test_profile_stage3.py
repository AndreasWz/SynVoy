"""Stage 3 of the profile handoff: the search's rules and constraints (bin/select_modules.py).

Each rule at its trigger's edge and at its cap, each constraint it brings (C1, C3 for the
models, C7), C8 taken out again, the resolver's W, and the GOI miniprot timeout (a new param
with no rule).
`--modules legacy` never applies any of them (test_modules.KEYS_AFTER_V1).
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import annotate_goi_exons as age  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
import select_modules as sm  # noqa: E402

DEFAULTS = {
    "goi_candidate_regions_per_query": 2, "goi_max_models_per_locus": 10, "goi_margin_genes": 2,
    "sw_max_alignments_per_piece": 1, "aug_relaxed_length_min": 15, "fallback_short_min_aln_aa": 15,
    "dispersed_goi_min_alnlen": 50, "window_tblastn_comp_stats": 2, "window_tblastn_word_size": 3,
    "classify_complete_min_qcov": 0.7, "classify_fragment_max_qcov": 0.4,
}


def _rel(*ranks):
    return [{"gene_id": f"gene-{r}", "rank": r, "evalue_seg_no": 1e-6} for r in ranks]


def _goi(*lengths):
    return {"method": "gff_annotation",
            "exons": [{"id": f"GOI_q|exon_{i + 1}", "length_aa": n} for i, n in enumerate(lengths)]}


def _rules(on, features, goi_info=None, w=None):
    """settings and derivation of the rules of these modules on these features."""
    p = sm.design_profile(features, goi_info, 10, w)
    settings, _source, derivation = sm.apply_rules(on, p)
    return settings, derivation, p


def _constraint(res, cid):
    return {c["id"]: c for c in res["constraints"]}[cid]


# ============================================================ A1: models per GOI locus
@pytest.mark.parametrize("copies,want,uncapped", [
    (1.5, 10, None), (9.0, 10, None), (10.0, 11, None), (31.0, 32, None), (59.0, 60, None), (60.0, 60, 61)])
def test_a1_models_per_locus_are_10_up_to_9_copies_then_c_plus_1_at_most_60(copies, want, uncapped):
    settings, derivation, _p = _rules({"A1"}, {"primary_copies": copies})
    assert settings["goi_max_models_per_locus"] == want
    assert derivation["goi_max_models_per_locus"].get("uncapped") == uncapped


# ============================================================ A1: the GOI margin
@pytest.mark.parametrize("ranks,w,want", [
    ((-4, 2), 7, 5),        # the farthest family member at rank 4: 4 + 1
    ((9,), 7, 7),           # beyond the window margin: W
    ((), 7, 2),             # A1 by coverage alone, no member in the walk: the default 2
    ((1,), 7, 2),           # a member next door: max(2, 2)
    ((-4,), 3, 3),          # a narrower window margin caps it
    ((12,), None, 7),       # W unknown: the shipped 7
])
def test_a1_goi_margin_reaches_the_home_array_within_the_window(ranks, w, want):
    f = {"primary_copies": 2.0, "relatives_nearby": _rel(*ranks)}
    settings, derivation, p = _rules({"A1"}, f, w=w)
    assert settings["goi_margin_genes"] == want
    assert derivation["goi_margin_genes"]["line"] == (
        f"goi_margin_genes {want} = min(W, max(2, R + 1)), R = {p['R']}, "
        f"W = {'unknown' if w is None else w} (A1)")


def test_c1_the_goi_margin_never_passes_the_window_margin():
    f = {"primary_copies": 2.0, "relatives_nearby": _rel(9)}
    res = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="auto", features=f,
                      design=sm.design_profile(f, None, 10, 5))
    c1 = _constraint(res, "C1")
    assert res["settings"]["goi_margin_genes"] == 5 and c1["applies"] and c1["ok"]
    # a value the user set beyond the window margin wins, with a warning
    user = rep.resolve(dict(DEFAULTS, goi_margin_genes=9), auto_preset=None, modules_spec="auto",
                       features=f, user_set={"goi_margin_genes": 9},
                       design=sm.design_profile(f, None, 10, 5))
    c1 = _constraint(user, "C1")
    assert user["settings"]["goi_margin_genes"] == 9
    assert not c1["ok"] and c1["on_violation"] == "warn" and "window margin 5" in c1["detail"]
    # legacy: the default 2 lies within any shipped window margin
    legacy = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="legacy", features=f)
    assert _constraint(legacy, "C1")["ok"]


# ============================================================ A1 / C1: SW alignments per piece
@pytest.mark.parametrize("on,features,goi,want", [
    ({"A1"}, {"primary_copies": 1.5}, None, 3),                  # two home copies
    ({"A1"}, {"primary_copies": None, "relatives_nearby": _rel(3)}, None, 3),
    ({"C1"}, {"goi_real_exons": 2}, _goi(30, 40), 3),           # two exons
    ({"C1"}, {"primary_copies": 1.0, "goi_real_exons": 1}, _goi(80), 1),   # forced, one exon, one copy
])
def test_sw_alignments_per_piece_are_3_for_copies_or_exons(on, features, goi, want):
    settings, _d, _p = _rules(on, features, goi)
    assert settings["sw_max_alignments_per_piece"] == want


def test_under_auto_a_single_exon_single_copy_gene_keeps_one_sw_alignment():
    f = {"query_len": 300, "primary_copies": 1.0, "goi_real_exons": 1}
    res = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="auto", features=f)
    assert res["settings"]["sw_max_alignments_per_piece"] == 1 and res["source"]["sw_max_alignments_per_piece"] == "default"
    f = dict(f, goi_real_exons=2)
    res = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="auto", features=f)
    assert res["settings"]["sw_max_alignments_per_piece"] == 3 and res["source"]["sw_max_alignments_per_piece"] == "module:C1"


# ============================================================ B1 / C1: the short-piece floors
@pytest.mark.parametrize("query_len,exons,want", [
    (70, (35, 35), 21),     # a 2-exon 70-aa precursor: hits of about 35 aa
    (50, (50,), 30),        # a 50-aa single-exon query
    (300, (200, 100), 50),  # B1 does not fire, C1 does: 0.6 x 200 is above 50
    (12, (12,), 10),        # never below 10
    (None, (), 50),         # piece length unknown: the default
])
def test_the_dispersed_rescue_floor_follows_the_longest_piece(query_len, exons, want):
    goi = _goi(*exons) if exons else None
    on = {"C1"} if query_len and query_len >= 100 else {"B1"}
    settings, derivation, p = _rules(on, {"query_len": query_len}, goi)
    assert settings["dispersed_goi_min_alnlen"] == want
    assert derivation["dispersed_goi_min_alnlen"]["inputs"] == {"Lpiece": p["Lpiece"]}


def test_b1_lowers_the_window_and_fallback_floors_to_10_and_c1_alone_does_not():
    settings, _d, _p = _rules({"B1"}, {"query_len": 60})
    assert (settings["aug_relaxed_length_min"], settings["fallback_short_min_aln_aa"]) == (10, 10)
    settings, _d, _p = _rules({"C1"}, {"query_len": 300, "goi_real_exons": 3}, _goi(10, 12, 15))
    assert "aug_relaxed_length_min" not in settings and "fallback_short_min_aln_aa" not in settings


def test_c7_a_length_floor_from_a_rule_lies_between_10_and_the_piece_bound(monkeypatch):
    f = {"query_len": 70}
    res = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="auto", features=f,
                      design=sm.design_profile(f, _goi(35, 35), 10, 7))
    c7 = _constraint(res, "C7")
    assert c7["applies"] and c7["ok"] and "bound max(10, round(0.6 x Lpiece)) = 21" in c7["detail"]
    bug = sm.Rule("dispersed_goi_min_alnlen", ("B1", "C1"), lambda p, on: 40, "40", doc="a bug")
    monkeypatch.setattr(sm, "RULES", [bug if r.key == bug.key else r for r in sm.RULES])
    res = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="auto", features=f,
                      design=sm.design_profile(f, _goi(35, 35), 10, 7))
    c7 = _constraint(res, "C7")
    assert not c7["ok"] and c7["on_violation"] == "stop"
    # a floor the user set is theirs (warn), and no rule-set floor means C7 does not apply
    plain = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="legacy", features=f)
    assert not _constraint(plain, "C7")["applies"]


# ============================================================ C8 taken out again
def test_b1_turns_composition_statistics_off_also_for_a_low_complexity_query():
    """C8 (keep tblastn's composition statistics, 2, for a short low-complexity query) was tried
    in stage 3 and removed: on the family pilot it cost real tblastn window hits on curated SPIN
    copies (Harpegnathos 33 -> 19 hits in one window) and shortened a counted call, with decoys
    at 0 MEDIUM+ either way. B1 gives 0 for every short query again, as in stage 1."""
    import query_features as qf
    mel = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"
    f = {"query_len": len(mel), **{k: v for k, v in qf.query_features(mel).items() if k in ("seg_frac", "repeat_frac")}}
    res = rep.resolve(dict(DEFAULTS), auto_preset=None, modules_spec="auto", features=f)
    assert {"B1", "B2"} <= {m["id"] for m in res["modules_applied"]}      # short and low-complexity
    assert res["settings"]["window_tblastn_comp_stats"] == 0
    assert res["source"]["window_tblastn_comp_stats"] == "module:B1"
    assert "C8" not in {c["id"] for c in res["constraints"]}
    assert "lowc" not in sm.PROFILE_SYMBOLS and "lowc" not in res["profile"]
    rule = sm.RULES_BY_KEY["window_tblastn_comp_stats"]
    assert rule.inputs == () and rule.formula({}, frozenset({"B1"})) == 0


# ============================================================ the resolver passes W
def test_the_resolver_cli_records_w_and_caps_the_margin_with_it(tmp_path, monkeypatch):
    hp = {"query_len": 300, "primary_copies": 2.0, "loci": [], "query_features": {"query_len": 300}}
    home = {"relatives_nearby": _rel(-6, 2)}
    (tmp_path / "hp.json").write_text(json.dumps(hp))
    (tmp_path / "home.json").write_text(json.dumps(home))
    (tmp_path / "defaults.json").write_text(json.dumps(DEFAULTS))
    monkeypatch.setattr(sys, "argv", [
        "resolve_effective_params.py", "--hit_profile", str(tmp_path / "hp.json"),
        "--home_context", str(tmp_path / "home.json"), "--modules", "auto",
        "--n_flanking_genes", "10", "--window_margin_genes", "5",
        "--defaults", str(tmp_path / "defaults.json"), "--output", str(tmp_path / "eff.json")])
    rep.main()
    out = json.loads((tmp_path / "eff.json").read_text())
    assert out["profile"]["W"] == 5 and out["profile"]["R"] == 6
    assert out["settings"]["goi_margin_genes"] == 5
    assert out["derivation"]["goi_margin_genes"]["line"] == "goi_margin_genes 5 = min(W, max(2, R + 1)), R = 6, W = 5 (A1)"


def test_the_module_passes_the_window_margin_to_the_resolver():
    text = (ROOT / "modules" / "resolve_effective_params.nf").read_text()
    assert "--window_margin_genes ${params.window_margin_genes}" in text


# ============================================================ goi_miniprot_timeout
class _Proc:
    returncode, stdout, stderr = 0, "", ""


@pytest.mark.parametrize("sensitive,want", [(True, 300), (False, 900)])
def test_the_goi_miniprot_timeout_reaches_the_sensitive_runs_only(monkeypatch, sensitive, want):
    seen = {}

    def run(cmd, **kw):
        seen["timeout"] = kw.get("timeout")
        return _Proc()
    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(age, "GOI_MINIPROT_TIMEOUT", 300)
    monkeypatch.setattr(age, "FLANKING_MINIPROT_TIMEOUT", 900)
    monkeypatch.setattr(age.subprocess, "run", run)
    age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=sensitive)
    assert seen["timeout"] == want


def test_a_goi_run_beyond_its_guard_names_its_limit_and_stops(monkeypatch):
    """Fix-phase item 18: the guard stops the task (StepFailed) and names its parameter; until
    2026-10-28 a timeout gave no model and was counted as miniprot_timeout_goi."""
    def boom(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout"))
    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(age, "GOI_MINIPROT_TIMEOUT", 45)
    monkeypatch.setattr(age.subprocess, "run", boom)
    with pytest.raises(age.StepFailed, match=r"within 45 s \(--goi_miniprot_timeout\)"):
        age._run_miniprot("MAAA", "ACGT" * 30, "c1", sensitive=True)
    assert "miniprot_timeout_goi" not in age.LIMIT_COUNTS


def test_the_search_takes_the_goi_timeout_and_sets_it_per_genome_task():
    help_text = subprocess.run([sys.executable, str(ROOT / "bin" / "iterative_search_runner.py"), "--help"],
                               capture_output=True, text=True).stdout
    assert "--goi_miniprot_timeout" in help_text
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    assert "_age.GOI_MINIPROT_TIMEOUT = int(getattr(args, 'goi_miniprot_timeout', 3600))" in src
    nf = (ROOT / "modules" / "iterative_search.nf").read_text()
    assert "--goi_miniprot_timeout ${params.goi_miniprot_timeout}" in nf
    assert age.GOI_MINIPROT_TIMEOUT == 3600
