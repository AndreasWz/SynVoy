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
  2. preset modules (--modules auto or a list; select_modules.py): one rule per key
  3. preset overrides from PRESET_OVERRIDES[preset]: the chooser's preset under
     --modules legacy (the default), or --preset_override in any mode
  4. values the user set themselves           (detected via --shipped_config)
  5. --user_override KEY=VAL                  (explicit, for callers that know)

Under --modules legacy the module triggers are evaluated and recorded only
(modules_would_fire); the preset is applied exactly as before the modules existed, and
its nominal module set is listed in modules_applied (reason "legacy: <preset>").

Besides the settings, the output records (any mode): ``profile``, the design quantities
the module rules read (select_modules.design_profile); ``derivation``, per key that a
module rule, the preset or the user set, its value, the default and what set it (for a
module rule also the formula and its inputs: "goi_candidate_regions_per_query 4 =
max(2, c + 1), c = 3 (A1)"); and ``constraints``, the checks of how the settings must
relate to each other. A violation that only the module rules made is a bug: the script
exits 2 and the run stops before the search.

Step 3 exists because the --defaults map is the launch-time ``params``: it
already contains whatever the user passed on the command line, in a params
file, a ``-c`` config or a ``-profile``. Without step 3 a preset silently
overwrote those values -- ``--goi_relative_evalue 1e-4`` did nothing whenever
the tandem-family preset fired. A value counts as user-set when it differs from
the default shipped in ``nextflow.config``; the one case this cannot see is a user
explicitly passing the shipped default, for which ``--auto_apply_preset false``
remains the way to pin it.

If --disable_auto_apply is set, skip step 2 entirely — defaults flow through
unchanged. If both --preset_override and an auto_preset.json are provided, the
explicit override wins. An unparseable / sentinel auto_preset is a no-op.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

try:
    import select_modules
except ImportError:  # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import select_modules


# ──────────────────────────────────────────────────────────────────────
# Preset overrides — single source of truth, mirror of conf/presets/*.config.
# Update both files together. The sync test compares this dict against the
# Groovy configs by regex-extracting `param = value` lines.
# ──────────────────────────────────────────────────────────────────────

PRESET_OVERRIDES: dict[str, dict[str, object]] = {
    # A preset holds only what differs from the defaults FOR A REASON.
    # 2026-10-08: the identity thresholds of the presets for queries of 100 aa or more went
    # (two raised HIGH/MEDIUM to 55/40 whatever the chooser decided).
    # 2026-10-09: four hit floors left the short-peptide preset (sw_min_score 10,
    # sw_min_identity 8, min_hit_identity 8, min_hit_length 8): the E-value filters of the
    # search are stricter than these floors, so lowering them changed nothing.
    # 2026-10-15: every identity threshold and strict_goi_family went with the
    # identity-based GOI label: one judge (bin/synvoy_label.py) labels every GOI call from
    # the E-value in the expected place, specificity and completeness, and no preset moves it.
    "preset_short_peptide": {
        # Model status only (fragment / partial / complete): one missed exon is a large share
        # of a short protein.
        "classify_complete_min_qcov": 0.5,
        "classify_fragment_max_qcov": 0.3,
    },
    "preset_tandem_family": {
        # Several copies next to each other: also weaker relatives of the GOI are kept
        # out of the anchors.
        "goi_relative_evalue": 1e-3,
    },
    "preset_single_copy": {},
    # No goi_family_tokens here on purpose: this preset is auto-applied to ANY large
    # gene family, and a hard-coded TP53 list made the strict family gate (removed
    # 2026-10-15) demote hit-chain calls of unrelated queries (Adh, luciferase) for not
    # being named TP53.
    # classify_fragment_max_qcov 0.35 left this preset on 2026-10-09: in the three benchmark
    # families run with it (2,651 GOI rows) no gene model lay between 0.35 and 0.40 coverage,
    # so it decided nothing, and a LOWER fragment limit is the lenient direction.
    "preset_paralog_discrimination": {
        "adaptive_score_floor_frac": 0.45,
        "adaptive_score_floor_abs": 0.08,
        "adaptive_max_regions": 4,
    },
}

VALID_PRESETS = set(PRESET_OVERRIDES.keys())


def _norm(value):
    """Type-tolerant comparison key: Nextflow hands CLI values over as strings, so
    "35" must equal 35.0 and "false" must equal False."""
    if isinstance(value, bool):
        return ("bool", value)
    if value is None:
        return ("none", None)
    text = str(value).strip()
    if text.lower() in ("true", "false"):
        return ("bool", text.lower() == "true")
    try:
        return ("num", float(text))
    except ValueError:
        return ("str", text)


def _strip_config_comment(line: str) -> str:
    """Drop a trailing ``//`` comment that is not inside a quoted string."""
    quote = ""
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = ""
        elif ch in "'\"":
            quote = ch
        elif ch == "/" and line[i:i + 2] == "//":
            return line[:i]
    return line


def load_shipped_defaults(config_path: str, keys) -> dict:
    """Literal defaults of ``keys`` from the top-level ``params {}`` block of the
    shipped nextflow.config. Missing/unreadable file or key -> simply absent, which
    makes the caller fall back to "not user-set" for that key (the old behaviour)."""
    wanted = set(keys)
    found: dict = {}
    try:
        with open(config_path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return found
    depth = 0
    in_params = False
    for raw in lines:
        line = _strip_config_comment(raw)
        if not in_params:
            if re.match(r"^params\s*\{", line):
                in_params = True
                depth = line.count("{") - line.count("}")
            continue
        if depth == 1:
            m = re.match(r"^\s*([A-Za-z_]\w*)\s*=\s*(.*?)\s*$", line)
            if m and m.group(1) in wanted and m.group(1) not in found:
                val = m.group(2)
                if len(val) >= 2 and val[0] == val[-1] and val[0] in "'\"":
                    val = val[1:-1]
                found[m.group(1)] = val
        depth += line.count("{") - line.count("}")
        if depth <= 0:
            break
    return found


def _typed(text):
    """A shipped default as written in nextflow.config, typed for the JSON record."""
    kind, value = _norm(text)
    if kind == "num":
        return int(value) if value.is_integer() and re.fullmatch(r"-?\d+", str(text).strip()) else value
    return value


def detect_user_set(defaults: dict, shipped: dict) -> dict:
    """Keys whose launch-time value differs from the shipped default, i.e. values
    the user chose (CLI flag, params file, -c config, or a profile)."""
    return {k: v for k, v in defaults.items()
            if k in shipped and _norm(v) != _norm(shipped[k])}


def resolve(defaults: dict, *, auto_preset: dict | None, preset_override: str = "",
            disable_auto_apply: bool = False,
            user_overrides: dict | None = None,
            user_set: dict | None = None,
            modules_spec: object = "legacy",
            features: dict | None = None,
            design: dict | None = None,
            shipped: dict | None = None) -> dict:
    """Apply the resolution order (defaults < modules < preset < user-set < user CLI) and
    return a flat ``{settings: ..., preset_applied: ..., source: ...}`` result.

    ``source`` is a per-key trace: ``default`` | ``module:<ids>`` | ``preset:<name>`` |
    ``user_set`` | ``user_cli``. ``user_set_kept`` lists the keys where a value the user
    set was kept although a module or the applied preset wanted a different one.
    ``modules_spec`` is the --modules value (ValueError when it is not one); ``features``
    the flat dict of select_modules.flat_features; ``design`` the profile of
    select_modules.design_profile (default: the profile of the features alone);
    ``shipped`` the defaults of nextflow.config, for the derivation of user-set keys."""
    settings = dict(defaults)
    source = {k: "default" for k in defaults}
    preset_applied: str | None = None
    apply_reason = ""
    feats = features or {}
    prof = design if design is not None else select_modules.design_profile(feats)
    mode, chosen, fired = select_modules.choose(modules_spec, feats)
    modules_applied: list = []
    wanted: dict = {}
    derivation: dict = {}
    rules_on: list = []

    if disable_auto_apply:
        apply_reason = "auto_apply_preset=false → defaults flow through"
    else:
        if mode in ("auto", "list") and chosen:
            rules_on = [d["id"] for d in chosen]
            mod_settings, mod_source, mod_derivation = select_modules.apply_rules(rules_on, prof)
            for k, v in mod_settings.items():
                source[k] = mod_source[k] if k in settings else mod_source[k] + "(new)"
                settings[k] = v
                wanted[k] = v
                derivation[k] = dict(mod_derivation[k], default=defaults.get(k))
            modules_applied = chosen
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
        elif mode != "legacy":
            apply_reason = (f"modules {mode}: "
                            + (", ".join(f"{d['id']} {d['name']}" for d in chosen) or "none"))
        elif auto_preset:
            chosen_preset = (auto_preset.get("preset") or "").strip()
            if chosen_preset in VALID_PRESETS:
                preset_applied = chosen_preset
                apply_reason = f"auto_preset.json chose {chosen_preset}"
            elif chosen_preset:
                apply_reason = f"auto_preset.json's preset='{chosen_preset}' is unknown — defaults applied"
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
            wanted[k] = v
            derivation[k] = {"value": v, "default": defaults.get(k), "rule": f"preset:{preset_applied}",
                             "line": f"{k} {v} (preset {preset_applied})"}
        if mode == "legacy":
            by_id = {m.id: m for m in select_modules.MODULES}
            modules_applied = [{"id": i, "name": by_id[i].name, "reason": f"legacy: {preset_applied}"}
                               for i in select_modules.LEGACY_ALIAS.get(preset_applied, ())]

    user_set_kept = []
    if user_set:
        for k in sorted(user_set):
            v = user_set[k]
            if k in wanted and _norm(wanted[k]) != _norm(v):
                user_set_kept.append(k)
            settings[k] = v
            source[k] = "user_set"
            derivation[k] = {"value": v, "default": _typed((shipped or {}).get(k)), "rule": "user_set",
                             "line": f"{k} {v} (set by the user)"}
            if k in wanted:
                derivation[k]["wanted"] = wanted[k]

    if user_overrides:
        for k, v in user_overrides.items():
            settings[k] = v
            source[k] = "user_cli"
            derivation[k] = {"value": v,
                             "default": _typed(shipped[k]) if k in (shipped or {}) else defaults.get(k),
                             "rule": "user_cli", "line": f"{k} {v} (--user_override)"}
            if k in wanted:
                derivation[k]["wanted"] = wanted[k]

    legacy_pick = (auto_preset or {}).get("preset") if isinstance(auto_preset, dict) else None
    return {
        "preset_applied": preset_applied,
        "apply_reason": apply_reason,
        "settings": settings,
        "source": source,
        "user_set_kept": user_set_kept,
        "modules_mode": mode,
        "modules_applied": modules_applied,
        "modules_would_fire": fired,
        "legacy_preset": legacy_pick,
        "features": feats,
        "home_copies": select_modules.home_copies(feats),
        "profile": prof,
        "derivation": derivation,
        "constraints": select_modules.check_constraints(settings, source, prof, rules_on),
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
    ap.add_argument("--shipped_config", default="",
                    help="Path to the shipped nextflow.config. A --defaults value that "
                         "differs from the default declared there was set by the user "
                         "and is kept over the preset. Omit for the old behaviour "
                         "(a preset overrides every key it names).")
    ap.add_argument("--hit_profile", default="",
                    help="hit_profile.json of LOCATE_GENE (with query_features); a sentinel or "
                         "missing file = no hit features")
    ap.add_argument("--goi_info", default="", help="goi_info.json of ANNOTATE_GOI")
    ap.add_argument("--home_context", default="",
                    help="home_context.json of HOME_CONTEXT (a sentinel or missing file = none)")
    ap.add_argument("--modules", default="legacy",
                    help="Preset modules: legacy (default: the old preset, triggers recorded "
                         "only), auto, none, a list (short_query,local_family) or deltas "
                         "(auto,-dispersed_family). See select_modules.py")
    ap.add_argument("--n_flanking_genes", type=int, default=None,
                    help="n_flanking_genes of the run: the profile's n (recorded; no rule reads "
                         "it yet)")
    ap.add_argument("--window_margin_genes", type=float, default=None,
                    help="window_margin_genes of the run: the profile's W, the cap of the GOI "
                         "margin rule (constraint C1)")
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

    user_set: dict[str, object] = {}
    shipped: dict = {}
    if args.shipped_config:
        shipped = load_shipped_defaults(args.shipped_config, defaults.keys())
        if not shipped:
            print(f"WARN: no defaults could be read from --shipped_config "
                  f"'{args.shipped_config}'; user-set values cannot be told apart "
                  f"from defaults, so a preset may override them", file=sys.stderr)
        user_set = detect_user_set(defaults, shipped)

    goi_info = _read_json_or_empty(args.goi_info)
    features = select_modules.flat_features(_read_json_or_empty(args.hit_profile),
                                            goi_info,
                                            _read_json_or_empty(args.home_context))
    try:
        result = resolve(
            defaults,
            auto_preset=auto_preset,
            preset_override=args.preset_override,
            disable_auto_apply=args.disable_auto_apply,
            user_overrides=user_overrides,
            user_set=user_set,
            modules_spec=args.modules,
            features=features,
            design=select_modules.design_profile(features, goi_info, args.n_flanking_genes,
                                                 args.window_margin_genes),
            shipped=shipped,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)

    print(f"EFFECTIVE-PRESET: {result['preset_applied']} ({result['apply_reason']})")
    applied = ", ".join(f"{d['id']} {d['name']} ({d['reason']})" for d in result["modules_applied"])
    fired = ", ".join(f"{d['id']} {d['name']} ({d['reason']})" for d in result["modules_would_fire"])
    print(f"MODULES: mode {result['modules_mode']}; applied: {applied or 'none'}; "
          f"triggers firing: {fired or 'none'}")
    for k in sorted(result["derivation"]):
        if str(result["source"].get(k, "")).startswith("module:"):
            print(f"RULE: {result['derivation'][k]['line']}")
    if result["user_set_kept"]:
        kept = ", ".join(f"{k}={result['settings'][k]}" for k in result["user_set_kept"])
        print(f"EFFECTIVE-PRESET: kept {len(result['user_set_kept'])} user-set value(s) "
              f"over {result['preset_applied'] or 'the modules'}: {kept}")
    errors = []
    for c in result["constraints"]:
        if c["applies"] and not c["ok"]:
            msg = f"constraint {c['id']} violated ({c['doc']}): {c['detail']}"
            if c["on_violation"] == "stop":
                errors.append(msg)
            else:
                print(f"WARN: {msg}; kept, because a value the user set (or a default) is involved",
                      file=sys.stderr)
    if errors:
        for msg in errors:
            print(f"ERROR: {msg}. The module rules made it: a bug, the run stops before the search.",
                  file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
