#!/usr/bin/env python3
"""The run's settings, read-only (import-only helper; profile handoff, stage 4).

RESOLVE_EFFECTIVE_PARAMS writes effective_params.json once per run. Its ``settings`` map
holds every key of main.nf's presetDefaultKeys() at the value the run uses: the shipped
default, the old preset or a preset module's rule (bin/select_modules.py), or a value the
user set, which always wins. EXTRACT_FLANKING, ITERATIVE_SEARCH and CLUSTER_REGIONS get the
map as ``${settings.x}``; a step after the search gets the file itself and reads its keys
here:

    settings = synvoy_settings.load_settings(args.effective_params)
    cap, _src = synvoy_settings.pick(args.max_goi_per_genome, settings,
                                     "max_goi_per_genome", 10, synvoy_settings.whole)

A value given on the script's own command line wins (a direct call, a test, an option of
--plot_extra_args); else the run's settings value; else the script's own default, which is
what a call without the file gets (a re-draw of a run older than stage 4). The mapping is
read-only, and nothing here writes into an argparse namespace: the value a step uses is the
one the run decided, or the one its caller gave.
"""
from __future__ import annotations

import json
import os
from types import MappingProxyType
from typing import Any, Callable, Mapping, Optional, Tuple

# key -> the scripts that read it through load_settings. tests/test_modules.py checks that
# each rule key has a reader (this map, or ${settings.x} in a module) and that every script
# named here reads its key through this helper.
READERS: Mapping[str, Tuple[str, ...]] = MappingProxyType({
    "locus_ownership_max_genes_per_locus": ("build_home_paralog_panel.py",),
    "paralog_panel_max_family": ("build_home_paralog_panel.py",),
    "goi_hull_min_coverage": ("rescue_goi_hull.py",),
    "max_goi_per_genome": ("plot_synteny.py",),
})

_EMPTY: Mapping[str, Any] = MappingProxyType({})


def load_settings(json_path: Optional[str]) -> Mapping[str, Any]:
    """The ``settings`` map of effective_params.json as a read-only mapping. Empty for no
    path, a missing or empty file, or a sentinel (NO_EFFECTIVE_PARAMS, "sentinel: ...").
    ValueError for a file that is no such JSON: the step would otherwise run on its defaults
    without a word."""
    if not json_path or not os.path.exists(json_path) or os.path.getsize(json_path) == 0:
        return _EMPTY
    with open(json_path) as fh:
        text = fh.read().strip()
    if not text or text.lower().startswith("sentinel"):
        return _EMPTY
    try:
        data = json.loads(text)
    except ValueError as exc:
        raise ValueError(f"{json_path} is not a JSON file: {exc}") from None
    settings = data.get("settings") if isinstance(data, dict) else None
    if not isinstance(settings, dict):
        raise ValueError(f"{json_path} has no settings map (is it effective_params.json?)")
    return MappingProxyType(dict(settings))


def whole(value: Any) -> int:
    """A count: 8, 8.0, "8" and "8.0" give 8. ValueError for a fraction, a word or a
    boolean (Nextflow hands a value from the command line over as text)."""
    if isinstance(value, bool):
        raise ValueError(f"not a count: {value!r}")
    number = float(str(value).strip())
    if not number.is_integer():
        raise ValueError(f"not a whole number: {value!r}")
    return int(number)


def pick(given: Any, settings: Mapping[str, Any], key: str, default: Any,
         kind: Callable[[Any], Any]) -> Tuple[Any, str]:
    """(value, where it came from) of one key: ``given`` (the script's own option, None when
    the call did not pass it), else the settings value, else ``default``; converted with
    ``kind``. Where: "option", "settings" or "default"."""
    if given is not None:
        return kind(given), "option"
    if key in settings and settings[key] is not None:
        try:
            return kind(settings[key]), "settings"
        except (TypeError, ValueError) as exc:
            raise ValueError(f"settings value of {key} ({settings[key]!r}): {exc}") from None
    return kind(default), "default"
