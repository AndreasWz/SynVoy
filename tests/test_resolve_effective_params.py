"""Tests for docs/TODO.md §1f Phase B — preset auto-apply via dynamic injection.

Three layers:

1. ``resolve()`` unit logic — defaults < preset < user_cli, with edge cases.
2. The CLI script as a black box — sentinel handling, JSON round-trip.
3. **Sync test** — ``PRESET_OVERRIDES`` (Python) must match
   ``conf/presets/*.config`` (Groovy). Single source of truth lives in
   ``bin/resolve_effective_params.py`` for testability + Nextflow channel use;
   the .config files remain for the static `-profile preset_X` startup path
   and must be kept in sync by hand. This test catches drift.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import resolve_effective_params as rep  # noqa: E402


# ──────────────────────────────────────────────────────────────────────
# 1) Unit logic
# ──────────────────────────────────────────────────────────────────────

# A representative slice of the defaults the workflow ships with. Picked to
# overlap every preset so we can verify each override path lands.
DEFAULTS = {
    "sw_min_score": 50,
    "sw_min_identity": 30.0,
    "min_hit_identity": 25,
    "min_hit_length": 20,
    "classify_fragment_max_qcov": 0.4,
    "classify_complete_min_qcov": 0.7,
    "expand_goi_similar": True,
    "goi_family_tokens": "",
    "goi_relative_evalue": 1e-5,
    "adaptive_score_floor_frac": 0.3,
    "adaptive_score_floor_abs": 0.05,
    "adaptive_max_regions": 8,
    "adaptive_unique_gene_floor": 2,
}


def test_no_auto_preset_and_no_override_is_passthrough():
    res = rep.resolve(DEFAULTS, auto_preset=None, preset_override="")
    assert res["preset_applied"] is None
    # Settings equal defaults, source is all 'default'.
    assert res["settings"] == DEFAULTS
    assert set(res["source"].values()) == {"default"}
    assert "no auto_preset" in res["apply_reason"]


def test_auto_preset_applies_overrides_and_traces_source():
    auto = {"preset": "preset_paralog_discrimination"}
    res = rep.resolve(DEFAULTS, auto_preset=auto)
    assert res["preset_applied"] == "preset_paralog_discrimination"
    # Spot-check the key overrides land.
    assert res["settings"]["adaptive_score_floor_frac"] == 0.45
    assert res["settings"]["adaptive_max_regions"] == 4
    # Unaffected keys keep defaults.
    assert res["settings"]["sw_min_score"] == DEFAULTS["sw_min_score"]
    assert res["settings"]["classify_complete_min_qcov"] == DEFAULTS["classify_complete_min_qcov"]
    # Source trace marks the touched keys.
    assert res["source"]["adaptive_max_regions"].startswith("preset:")
    assert res["source"]["goi_relative_evalue"] == "default"
    assert res["source"]["sw_min_score"] == "default"


def test_explicit_override_beats_auto_preset():
    """User-pinned preset wins even if auto_preset.json suggested another."""
    auto = {"preset": "preset_short_peptide"}
    res = rep.resolve(DEFAULTS, auto_preset=auto, preset_override="preset_single_copy")
    assert res["preset_applied"] == "preset_single_copy"
    # single_copy sets nothing since 2026-10-15; the short-peptide values stay out.
    assert res["settings"] == DEFAULTS
    assert "user pinned" in res["apply_reason"]


def test_disable_auto_apply_bypasses_everything():
    auto = {"preset": "preset_paralog_discrimination"}
    res = rep.resolve(DEFAULTS, auto_preset=auto, preset_override="preset_single_copy",
                      disable_auto_apply=True)
    assert res["preset_applied"] is None
    assert res["settings"] == DEFAULTS
    assert "auto_apply_preset=false" in res["apply_reason"]


def test_user_cli_overrides_beat_preset():
    """Explicit user CLI values must take precedence over preset overrides."""
    auto = {"preset": "preset_short_peptide"}
    user = {"classify_complete_min_qcov": 0.99}
    res = rep.resolve(DEFAULTS, auto_preset=auto, user_overrides=user)
    assert res["settings"]["classify_complete_min_qcov"] == 0.99
    assert res["source"]["classify_complete_min_qcov"] == "user_cli"
    # Other preset overrides still apply.
    assert res["settings"]["classify_fragment_max_qcov"] == 0.3


def test_unknown_preset_in_auto_preset_falls_back_to_defaults():
    """A typo (e.g. `preset_typo`) in auto_preset.json must not crash; we
    keep defaults and report the failure mode."""
    auto = {"preset": "preset_typo"}
    res = rep.resolve(DEFAULTS, auto_preset=auto)
    assert res["preset_applied"] is None
    assert res["settings"] == DEFAULTS
    assert "unknown" in res["apply_reason"].lower()


def test_unknown_explicit_override_does_not_crash():
    res = rep.resolve(DEFAULTS, auto_preset=None, preset_override="preset_typo")
    assert res["preset_applied"] is None
    assert res["settings"] == DEFAULTS
    assert "not a known preset" in res["apply_reason"]


# ──────────────────────────────────────────────────────────────────────
# 2) CLI black-box tests
# ──────────────────────────────────────────────────────────────────────

SCRIPT = ROOT / "bin" / "resolve_effective_params.py"


def _run_cli(*argv, cwd: Path) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        check=True, cwd=cwd, capture_output=True,
    )
    with open(cwd / "out.json") as fh:
        return json.load(fh)


def test_cli_sentinel_auto_preset_is_noop(tmp_path: Path):
    """A sentinel auto_preset file (starts with 'sentinel:') must be treated as
    'no recommendation' rather than crashing the resolver."""
    sentinel = tmp_path / "NO_AUTO_PRESET"
    sentinel.write_text("sentinel: nothing to apply\n")
    defaults_path = tmp_path / "defaults.json"
    defaults_path.write_text(json.dumps(DEFAULTS))
    out = tmp_path / "out.json"

    result = _run_cli(
        "--auto_preset", str(sentinel),
        "--defaults", str(defaults_path),
        "--output", str(out),
        cwd=tmp_path,
    )
    assert result["preset_applied"] is None
    assert result["settings"] == DEFAULTS


def test_cli_missing_auto_preset_is_noop(tmp_path: Path):
    defaults_path = tmp_path / "defaults.json"
    defaults_path.write_text(json.dumps(DEFAULTS))
    result = _run_cli(
        "--auto_preset", str(tmp_path / "nope.json"),  # does not exist
        "--defaults", str(defaults_path),
        "--output", str(tmp_path / "out.json"),
        cwd=tmp_path,
    )
    assert result["preset_applied"] is None


def test_cli_typed_user_override_via_defaults(tmp_path: Path):
    """--user_override coerces VAL to the type of the defaults entry."""
    defaults_path = tmp_path / "defaults.json"
    defaults_path.write_text(json.dumps(DEFAULTS))
    auto = tmp_path / "auto.json"
    auto.write_text(json.dumps({"preset": "preset_single_copy"}))
    result = _run_cli(
        "--auto_preset", str(auto),
        "--defaults", str(defaults_path),
        "--user_override", "classify_complete_min_qcov=0.65",
        "--user_override", "adaptive_max_regions=12",
        "--user_override", "expand_goi_similar=false",
        "--output", str(tmp_path / "out.json"),
        cwd=tmp_path,
    )
    assert result["settings"]["classify_complete_min_qcov"] == 0.65  # float
    assert result["settings"]["adaptive_max_regions"] == 12          # int
    assert result["settings"]["expand_goi_similar"] is False         # bool


# ──────────────────────────────────────────────────────────────────────
# 3) Sync test — PRESET_OVERRIDES must match conf/presets/*.config.
# ──────────────────────────────────────────────────────────────────────

PRESET_CONFIG_MAP = {
    "preset_short_peptide": ROOT / "conf" / "presets" / "short_divergent_peptide.config",
    "preset_tandem_family": ROOT / "conf" / "presets" / "tandem_family.config",
    "preset_single_copy": ROOT / "conf" / "presets" / "single_copy_conserved.config",
    "preset_paralog_discrimination": ROOT / "conf" / "presets" / "paralog_discrimination.config",
}


def _parse_groovy_params_block(text: str) -> dict:
    """Pull `key = value` lines out of a Groovy `params { ... }` block.

    Strips comments and trailing whitespace. Coerces the RHS into Python types
    (true/false → bool, numeric → int/float, quoted → str). Lines with
    interpolation or comments-only are skipped."""
    m = re.search(r"params\s*\{([\s\S]*?)\n\}", text)
    if not m:
        return {}
    body = m.group(1)
    out: dict[str, object] = {}
    for raw_line in body.splitlines():
        line = raw_line.split("//", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip().rstrip(",")
        if val in ("true", "false"):
            out[key] = (val == "true")
        elif val.startswith("'") and val.endswith("'"):
            out[key] = val[1:-1]
        elif val.startswith('"') and val.endswith('"'):
            out[key] = val[1:-1]
        else:
            try:
                out[key] = int(val)
            except ValueError:
                try:
                    out[key] = float(val)
                except ValueError:
                    out[key] = val  # leave as-is; sync test will catch
    return out


def test_preset_overrides_match_groovy_configs():
    """Python's ``PRESET_OVERRIDES`` must mirror the Groovy ``.config`` files
    exactly. If you intentionally diverge, update both files together."""
    drifts: list[str] = []
    for preset, cfg_path in PRESET_CONFIG_MAP.items():
        assert cfg_path.exists(), f"missing config: {cfg_path}"
        groovy_vals = _parse_groovy_params_block(cfg_path.read_text())
        python_vals = rep.PRESET_OVERRIDES[preset]

        # Keys: Python must declare every key the Groovy file declares.
        py_keys = set(python_vals.keys())
        gr_keys = set(groovy_vals.keys())
        if py_keys != gr_keys:
            missing_py = gr_keys - py_keys
            missing_gr = py_keys - gr_keys
            if missing_py:
                drifts.append(f"{preset}: Python missing keys {sorted(missing_py)}")
            if missing_gr:
                drifts.append(f"{preset}: Groovy missing keys {sorted(missing_gr)}")
            continue

        # Values: numeric equality (50 == 50.0 OK), bool/str exact.
        for k, gv in groovy_vals.items():
            pv = python_vals[k]
            ok = (isinstance(gv, (int, float)) and isinstance(pv, (int, float)) and
                  float(gv) == float(pv)) or gv == pv
            if not ok:
                drifts.append(f"{preset}.{k}: groovy={gv!r} python={pv!r}")
    assert not drifts, "PRESET_OVERRIDES drift:\n  " + "\n  ".join(drifts)


# ──────────────────────────────────────────────────────────────────────
# 4) User-set values win over a preset
#
# The --defaults map is the launch-time `params`, so it already holds whatever
# the user set. Before --shipped_config existed a preset silently overwrote
# those values. A value counts as user-set when it differs from the default
# shipped in nextflow.config.
# ──────────────────────────────────────────────────────────────────────

# Must list exactly the keys main.nf feeds the resolver (presetDefaultKeys()).
def _preset_default_keys_from_main_nf() -> list[str]:
    text = (ROOT / "main.nf").read_text()
    m = re.search(r"def presetDefaultKeys\(\) \{\s*return \[(.*?)\]", text, re.S)
    assert m, "main.nf no longer defines presetDefaultKeys()"
    return re.findall(r"'([a-z_]+)'", m.group(1))


def test_shipped_defaults_are_found_for_every_preset_key():
    """If nextflow.config changes shape and a key stops being parsed, user-set
    detection silently stops working for it. Fail here instead."""
    keys = _preset_default_keys_from_main_nf()
    assert len(keys) == 19
    shipped = rep.load_shipped_defaults(str(ROOT / "nextflow.config"), keys)
    assert sorted(shipped) == sorted(keys), sorted(set(keys) - set(shipped))
    # Spot-check literals and that trailing comments/quotes are stripped.
    assert rep._norm(shipped["classify_fragment_max_qcov"]) == ("num", 0.4)
    assert rep._norm(shipped["expand_goi_similar"]) == ("bool", False)
    assert shipped["goi_family_tokens"] == ""
    assert rep._norm(shipped["classify_complete_min_qcov"]) == ("num", 0.7)
    assert rep._norm(shipped["goi_relative_evalue"]) == ("num", 1e-5)
    assert rep._norm(shipped["goi_candidate_regions_per_query"]) == ("num", 2)
    assert rep._norm(shipped["report_copies_in_slot"]) == ("bool", False)
    assert rep._norm(shipped["window_tblastn_comp_stats"]) == ("num", 2)
    assert rep._norm(shipped["miniprot_short_orf_max_len"]) == ("num", 60)


def test_every_preset_key_is_fed_by_main_nf():
    """A preset key missing from preset_default_keys never reaches `settings`
    with a default, and modules reading settings.<key> would get null."""
    keys = set(_preset_default_keys_from_main_nf())
    for preset, overrides in rep.PRESET_OVERRIDES.items():
        assert set(overrides) <= keys, (preset, sorted(set(overrides) - keys))


def test_detect_user_set_is_type_tolerant():
    shipped = {"a": "50.0", "b": "false", "c": "", "d": "20", "e": "true"}
    # Nextflow hands CLI values over as strings: none of these is a change.
    same = {"a": 50, "b": False, "c": "", "d": "20", "e": "True"}
    assert rep.detect_user_set(same, shipped) == {}
    changed = {"a": "35", "b": "true", "c": "TP53", "d": 10, "e": False}
    assert rep.detect_user_set(changed, shipped) == changed
    # A key the shipped config does not declare is never called user-set.
    assert rep.detect_user_set({"zzz": 1}, shipped) == {}


def test_user_set_value_beats_auto_preset_and_is_reported():
    defaults = dict(DEFAULTS, classify_complete_min_qcov="0.6")   # user typed 0.6
    res = rep.resolve(defaults, auto_preset={"preset": "preset_short_peptide"},
                      user_set={"classify_complete_min_qcov": "0.6"})
    assert res["settings"]["classify_complete_min_qcov"] == "0.6"
    assert res["source"]["classify_complete_min_qcov"] == "user_set"
    assert res["user_set_kept"] == ["classify_complete_min_qcov"]
    # Keys the user did not touch still take the preset's value.
    assert res["settings"]["classify_fragment_max_qcov"] == 0.3
    assert res["source"]["classify_fragment_max_qcov"] == "preset:preset_short_peptide"


def test_user_set_equal_to_preset_value_is_not_reported_as_kept():
    res = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"},
                      user_set={"classify_complete_min_qcov": 0.5})
    assert res["settings"]["classify_complete_min_qcov"] == 0.5
    assert res["user_set_kept"] == []


def test_explicit_user_override_still_beats_user_set():
    res = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"},
                      user_set={"sw_min_score": 15}, user_overrides={"sw_min_score": 12})
    assert res["settings"]["sw_min_score"] == 12
    assert res["source"]["sw_min_score"] == "user_cli"


def test_cli_keeps_a_value_the_user_set(tmp_path: Path):
    """End to end against the real nextflow.config: a launch-time value that
    differs from the shipped default survives an auto-applied preset."""
    keys = _preset_default_keys_from_main_nf()
    shipped = rep.load_shipped_defaults(str(ROOT / "nextflow.config"), keys)
    defaults = dict(shipped)
    defaults["classify_fragment_max_qcov"] = "0.2"        # as typed on the CLI
    (tmp_path / "defaults.json").write_text(json.dumps(defaults))
    (tmp_path / "auto.json").write_text(json.dumps({"preset": "preset_short_peptide"}))
    out = tmp_path / "eff.json"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--auto_preset", str(tmp_path / "auto.json"),
         "--defaults", str(tmp_path / "defaults.json"),
         "--shipped_config", str(ROOT / "nextflow.config"),
         "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    res = json.loads(out.read_text())
    assert res["preset_applied"] == "preset_short_peptide"
    assert res["settings"]["classify_fragment_max_qcov"] == "0.2"
    assert res["user_set_kept"] == ["classify_fragment_max_qcov"]
    assert "kept 1 user-set value(s)" in proc.stdout
    # An untouched key takes the preset.
    assert res["settings"]["classify_complete_min_qcov"] == 0.5


def test_cli_without_shipped_config_keeps_the_old_behaviour(tmp_path: Path):
    defaults = dict(DEFAULTS, classify_fragment_max_qcov="0.2")
    (tmp_path / "defaults.json").write_text(json.dumps(defaults))
    (tmp_path / "auto.json").write_text(json.dumps({"preset": "preset_short_peptide"}))
    out = tmp_path / "eff.json"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--auto_preset", str(tmp_path / "auto.json"),
         "--defaults", str(tmp_path / "defaults.json"), "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert json.loads(out.read_text())["settings"]["classify_fragment_max_qcov"] == 0.3


def test_unreadable_shipped_config_warns_and_falls_back(tmp_path: Path):
    assert rep.load_shipped_defaults(str(tmp_path / "missing.config"), ["a"]) == {}
    (tmp_path / "defaults.json").write_text(json.dumps(DEFAULTS))
    out = tmp_path / "eff.json"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--defaults", str(tmp_path / "defaults.json"),
         "--shipped_config", str(tmp_path / "missing.config"), "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0
    assert "no defaults could be read" in proc.stderr


def test_strip_config_comment_respects_quotes():
    assert rep._strip_config_comment("a = 5 // note") == "a = 5 "
    assert rep._strip_config_comment("a = 'x//y'  // note") == "a = 'x//y'  "
    assert rep._strip_config_comment('url = "http://h/p"') == 'url = "http://h/p"'


def test_no_preset_carries_gene_specific_family_tokens():
    """Presets are auto-applied to arbitrary queries. The paralog preset used to
    ship the TP53 family names, so an Adh query was searched with strict TP53
    tokens and its hit-chain calls on annotated genes were demoted to LOW."""
    for preset, overrides in rep.PRESET_OVERRIDES.items():
        assert "goi_family_tokens" not in overrides, preset
    for cfg in (ROOT / "conf" / "presets").glob("*.config"):
        code = "\n".join(l for l in cfg.read_text().splitlines()
                         if not l.strip().startswith("//"))
        assert "goi_family_tokens" not in code, cfg.name


# The parameters of the identity-based GOI label, removed on 2026-10-15 (main.nf removedParams).
REMOVED_LABEL_PARAMS = (
    "classify_high_min_identity", "classify_medium_min_identity", "classify_tandem_min_identity",
    "classify_tandem_min_qcov", "classify_high_min_collinear",
    "classify_fallback_strong_min_identity_floor", "disable_ambiguous_tier", "strict_goi_family",
    "disable_distance_autotune", "distance_autotune_close_pct", "distance_autotune_far_pct",
    "distance_autotune_max_relax", "distance_autotune_min_flanking",
    "distance_autotune_high_floor", "distance_autotune_medium_floor",
)


def test_presets_hold_only_what_still_acts():
    """2026-10-15: one judge labels every GOI call (bin/synvoy_label.py), so a preset may
    not move a label threshold. What is left: the model-status coverage pair for short
    peptides, the anchor filter for tandem families, the region floors for paralog families."""
    assert rep.PRESET_OVERRIDES == {
        "preset_short_peptide": {"classify_complete_min_qcov": 0.5, "classify_fragment_max_qcov": 0.3},
        "preset_tandem_family": {"goi_relative_evalue": 1e-3},
        "preset_single_copy": {},
        "preset_paralog_discrimination": {"adaptive_score_floor_frac": 0.45,
                                          "adaptive_score_floor_abs": 0.08,
                                          "adaptive_max_regions": 4},
    }


def test_removed_label_parameters_are_set_nowhere():
    keys = set(_preset_default_keys_from_main_nf())
    configs = sorted((ROOT / "conf").rglob("*.config"))
    for name in REMOVED_LABEL_PARAMS:
        assert name not in keys, name
        for overrides in rep.PRESET_OVERRIDES.values():
            assert name not in overrides, name
        for cfg in configs:
            code = "\n".join(l.split("//", 1)[0] for l in cfg.read_text().splitlines())
            assert not re.search(rf"\b{name}\s*=", code), (cfg.name, name)


def test_no_preset_sets_a_parameter_that_does_nothing():
    """`min_query_length` is checked before a preset can be chosen, `auto_params` cannot
    change after launch, `multi_profile` was never implemented."""
    for name, overrides in rep.PRESET_OVERRIDES.items():
        for dead in ("min_query_length", "auto_params", "multi_profile",
                     "max_flanking_goi_similarity") + REMOVED_LABEL_PARAMS:
            assert dead not in overrides, (name, dead)


def test_short_peptide_preset_no_longer_lowers_floors_that_do_not_bind():
    """Until 2026-10-09 the short-peptide preset lowered `sw_min_score`, `sw_min_identity`,
    `min_hit_identity` and `min_hit_length`. That changed nothing: the E-value filters are
    stricter than these floors. The numbers behind the statement are pinned here, so the
    entries come back only together with a reason."""
    import smith_waterman_search as sw

    short = rep.PRESET_OVERRIDES["preset_short_peptide"]
    for key in ("sw_min_score", "sw_min_identity", "min_hit_identity", "min_hit_length"):
        assert key not in short, key

    shipped = rep.load_shipped_defaults(
        str(ROOT / "nextflow.config"),
        ["sw_min_score", "min_hit_identity", "min_hit_length", "aug_relaxed_identity_min",
         "aug_relaxed_identity_factor", "aug_relaxed_length_min", "aug_relaxed_length_div",
         "min_query_length"])
    num = {k: rep._norm(v)[1] for k, v in shipped.items()}

    # Smith-Waterman: the lowest raw score with E <= 1 for the shortest query the pipeline
    # accepts, in a window of only 2 kb (six frames = 2 x the window in residues).
    lowest = next(s for s in range(1, 300)
                  if sw.karlin_altschul(s, int(num["min_query_length"]), 2 * 2000, 11, 1)[1] <= 1.0)
    assert lowest == 31
    assert num["sw_min_score"] < lowest

    # Window hits: the floors in use are the same for the default and the former preset value.
    for identity, length in ((num["min_hit_identity"], num["min_hit_length"]), (8, 8)):
        assert max(num["aug_relaxed_identity_min"], identity * num["aug_relaxed_identity_factor"]) == 15
        assert max(num["aug_relaxed_length_min"], length // num["aug_relaxed_length_div"]) == 15


def test_paralog_preset_leaves_the_fragment_limit_alone():
    """classify_fragment_max_qcov 0.35 left the paralog preset on 2026-10-09: in 2,651 GOI
    rows of three benchmark families no gene model lay between 0.35 and 0.40 coverage."""
    assert "classify_fragment_max_qcov" not in rep.PRESET_OVERRIDES["preset_paralog_discrimination"]
