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
    "min_query_length": 30,
    "classify_high_min_identity": 50.0,
    "classify_medium_min_identity": 35.0,
    "classify_tandem_min_identity": 40.0,
    "classify_fragment_max_qcov": 0.4,
    "classify_complete_min_qcov": 0.7,
    "expand_goi_similar": True,
    "strict_goi_family": False,
    "goi_family_tokens": "",
    "max_flanking_goi_similarity": 35.0,
    "adaptive_score_floor_frac": 0.3,
    "adaptive_score_floor_abs": 0.05,
    "adaptive_max_regions": 8,
    "adaptive_unique_gene_floor": 2,
    "auto_params": True,
    "multi_profile": True,
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
    assert res["settings"]["classify_high_min_identity"] == 55.0
    assert res["settings"]["adaptive_max_regions"] == 4
    assert res["settings"]["strict_goi_family"] is True
    # auto_params/multi_profile flipped off (§paralog_discrimination preset).
    assert res["settings"]["auto_params"] is False
    assert res["settings"]["multi_profile"] is False
    # Unaffected keys keep defaults.
    assert res["settings"]["sw_min_score"] == DEFAULTS["sw_min_score"]
    # Source trace marks the touched keys.
    assert res["source"]["classify_high_min_identity"].startswith("preset:")
    assert res["source"]["sw_min_score"] == "default"


def test_explicit_override_beats_auto_preset():
    """User-pinned preset wins even if auto_preset.json suggested another."""
    auto = {"preset": "preset_short_peptide"}
    res = rep.resolve(DEFAULTS, auto_preset=auto, preset_override="preset_single_copy")
    assert res["preset_applied"] == "preset_single_copy"
    assert res["settings"]["classify_high_min_identity"] == 55.0
    # Short-peptide-only key remains at the default (single_copy doesn't touch it).
    assert res["settings"]["sw_min_score"] == DEFAULTS["sw_min_score"]
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
    auto = {"preset": "preset_single_copy"}
    user = {"classify_high_min_identity": 99.9}
    res = rep.resolve(DEFAULTS, auto_preset=auto, user_overrides=user)
    assert res["settings"]["classify_high_min_identity"] == 99.9
    assert res["source"]["classify_high_min_identity"] == "user_cli"
    # Other preset overrides still apply.
    assert res["settings"]["strict_goi_family"] is True


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
        "--user_override", "classify_high_min_identity=42.5",
        "--user_override", "adaptive_max_regions=12",
        "--user_override", "strict_goi_family=false",
        "--output", str(tmp_path / "out.json"),
        cwd=tmp_path,
    )
    assert result["settings"]["classify_high_min_identity"] == 42.5  # float
    assert result["settings"]["adaptive_max_regions"] == 12         # int
    assert result["settings"]["strict_goi_family"] is False         # bool


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
