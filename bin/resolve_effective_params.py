#!/usr/bin/env python3
"""Resolve effective preset-affected params at workflow runtime (docs/TODO.md §1f).

Background: `bin/auto_select_preset.py` produces an *advisory* `auto_preset.json`
after LOCATE_GENE. Nextflow params are immutable after launch, so historically
the recommendation could not be applied within the same run — users had to
re-launch with `-profile <preset>`. §1f turns that advisory into an actual
runtime override by emitting a flat JSON map of resolved settings, which the
downstream processes (ITERATIVE_SEARCH, CLUSTER_REGIONS, EXTRACT_FLANKING)
read via Nextflow value-channel staging instead of `${params.X}`.

This file is the *single source of truth* for what each preset overrides
(``PRESET_OVERRIDES`` below). The Groovy ``conf/presets/*.config`` files
remain for the static `-profile preset_X` startup path and must be kept in
sync — ``tests/test_resolve_effective_params.py::test_presets_match_groovy_configs``
diffs them and fails CI on drift.

Resolution order (later wins):
  1. defaults from --defaults <json>          (passed in from nextflow.config params)
  2. preset overrides from PRESET_OVERRIDES[preset]
  3. --user_override KEY=VAL  (explicit user CLI flags collected by the caller)

If --disable_auto_apply is set, skip step 2 entirely — defaults flow through
unchanged. If both --preset_override and an auto_preset.json are provided, the
explicit override wins. An unparseable / sentinel auto_preset is a no-op.
"""
from __future__ import annotations

import argparse
import json
import os
import sys


# ──────────────────────────────────────────────────────────────────────
# Preset overrides — single source of truth, mirror of conf/presets/*.config.
# Update both files together. The sync test compares this dict against the
# Groovy configs by regex-extracting `param = value` lines.
# ──────────────────────────────────────────────────────────────────────

PRESET_OVERRIDES: dict[str, dict[str, object]] = {
    "preset_short_peptide": {
        "sw_min_score": 10,
        "sw_min_identity": 8.0,
        "min_hit_identity": 8,
        "min_hit_length": 8,
        "min_query_length": 20,
        "classify_complete_min_qcov": 0.5,
        "classify_fragment_max_qcov": 0.3,
        "classify_high_min_identity": 40.0,
        "classify_medium_min_identity": 25.0,
        "expand_goi_similar": False,
        "strict_goi_family": False,
    },
    "preset_tandem_family": {
        "expand_goi_similar": False,
        "max_flanking_goi_similarity": 25.0,
        "strict_goi_family": True,
        "classify_tandem_min_identity": 35.0,
        "classify_medium_min_identity": 30.0,
    },
    "preset_single_copy": {
        "strict_goi_family": True,
        "expand_goi_similar": False,
        "classify_high_min_identity": 55.0,
        "classify_medium_min_identity": 40.0,
        "classify_complete_min_qcov": 0.75,
    },
    "preset_paralog_discrimination": {
        "strict_goi_family": True,
        "goi_family_tokens": "TP53,TP63,TP73,TRP53,TRP63,TRP73,P53,P63,P73",
        "classify_high_min_identity": 55.0,
        "classify_medium_min_identity": 40.0,
        "classify_tandem_min_identity": 50.0,
        "classify_fragment_max_qcov": 0.35,
        "classify_complete_min_qcov": 0.75,
        "adaptive_score_floor_frac": 0.45,
        "adaptive_score_floor_abs": 0.08,
        "adaptive_max_regions": 4,
        "adaptive_unique_gene_floor": 3,
        "auto_params": False,
        "multi_profile": False,
    },
}

VALID_PRESETS = set(PRESET_OVERRIDES.keys())


def resolve(defaults: dict, *, auto_preset: dict | None, preset_override: str = "",
            disable_auto_apply: bool = False,
            user_overrides: dict | None = None) -> dict:
    """Apply the resolution order (defaults < preset < user CLI) and return a
    flat ``{settings: ..., preset_applied: ..., source: ...}`` result.

    ``source`` is a per-key trace: ``default`` | ``preset`` | ``user_cli``."""
    settings = dict(defaults)
    source = {k: "default" for k in defaults}
    preset_applied: str | None = None
    apply_reason = ""

    if disable_auto_apply:
        apply_reason = "auto_apply_preset=false → defaults flow through"
    else:
        if preset_override:
            if preset_override not in VALID_PRESETS:
                # Loud failure mode: a user-pinned override that doesn't exist is
                # likely a typo. Fall back to defaults but flag it in the JSON.
                apply_reason = (
                    f"preset_override='{preset_override}' is not a known preset — "
                    f"valid: {sorted(VALID_PRESETS)} — defaults applied"
                )
            else:
                preset_applied = preset_override
                apply_reason = f"user pinned preset_override={preset_override}"
        elif auto_preset:
            chosen = (auto_preset.get("preset") or "").strip()
            if chosen in VALID_PRESETS:
                preset_applied = chosen
                apply_reason = f"auto_preset.json chose {chosen}"
            elif chosen:
                apply_reason = f"auto_preset.json's preset='{chosen}' is unknown — defaults applied"
            else:
                apply_reason = "auto_preset.json had no preset field — defaults applied"
        else:
            apply_reason = "no auto_preset.json and no override — defaults applied"

    if preset_applied:
        for k, v in PRESET_OVERRIDES[preset_applied].items():
            if k not in settings:
                # Preset declares a key we didn't pass in defaults for. That means
                # the upstream params dict missed something — log it but still
                # apply so the preset can actually take effect.
                source[k] = f"preset:{preset_applied}(new)"
            else:
                source[k] = f"preset:{preset_applied}"
            settings[k] = v

    if user_overrides:
        for k, v in user_overrides.items():
            settings[k] = v
            source[k] = "user_cli"

    return {
        "preset_applied": preset_applied,
        "apply_reason": apply_reason,
        "settings": settings,
        "source": source,
    }


def _read_json_or_empty(path: str) -> dict | None:
    """Return the JSON, or None for a sentinel/missing/unparseable file."""
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            text = fh.read().strip()
        if not text:
            return None
        # Sentinel files (NO_AUTO_PRESET etc) start with "sentinel:" by
        # SynVoy convention, but they may also just be empty.
        if text.lower().startswith("sentinel"):
            return None
        return json.loads(text)
    except (OSError, json.JSONDecodeError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--auto_preset", default="",
                    help="auto_preset.json from LOCATE_GENE (may be a sentinel or missing)")
    ap.add_argument("--preset_override", default="",
                    help="Pin a specific preset, bypassing auto-detection")
    ap.add_argument("--disable_auto_apply", action="store_true",
                    help="Pass defaults through without applying any preset")
    ap.add_argument("--defaults", required=True,
                    help="JSON string OR path to a JSON file with the default param values")
    ap.add_argument("--user_override", action="append", default=[], metavar="KEY=VAL",
                    help="Override a single param; repeat for multiple. CLI-style"
                         " explicit user overrides take precedence over presets.")
    ap.add_argument("--output", required=True, help="Where to write the resolved JSON")
    args = ap.parse_args()

    # Defaults: JSON string OR file path
    defaults_arg = args.defaults
    if os.path.exists(defaults_arg):
        with open(defaults_arg) as fh:
            defaults = json.load(fh)
    else:
        try:
            defaults = json.loads(defaults_arg)
        except json.JSONDecodeError as exc:
            print(f"ERROR: --defaults is neither a path nor valid JSON: {exc}", file=sys.stderr)
            sys.exit(2)

    user_overrides: dict[str, object] = {}
    for kv in args.user_override:
        if "=" not in kv:
            print(f"WARN: ignoring malformed --user_override '{kv}' (expected KEY=VAL)",
                  file=sys.stderr)
            continue
        k, v = kv.split("=", 1)
        # Preserve type of the defaults entry when possible (so "55.0" → float).
        if k in defaults:
            ref = defaults[k]
            try:
                if isinstance(ref, bool):
                    v = v.strip().lower() in ("true", "1", "yes", "on")
                elif isinstance(ref, int):
                    v = int(v)
                elif isinstance(ref, float):
                    v = float(v)
            except ValueError:
                pass  # leave as string
        user_overrides[k] = v

    auto_preset = _read_json_or_empty(args.auto_preset)

    result = resolve(
        defaults,
        auto_preset=auto_preset,
        preset_override=args.preset_override,
        disable_auto_apply=args.disable_auto_apply,
        user_overrides=user_overrides,
    )

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)

    print(f"EFFECTIVE-PRESET: {result['preset_applied']} ({result['apply_reason']})")


if __name__ == "__main__":
    main()
