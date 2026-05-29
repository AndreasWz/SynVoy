#!/usr/bin/env python3
"""Pick a query-type preset from a hit_profile.json (docs/TODO.md §1f).

Maps the LOCATE_GENE hit distribution onto one of the bundled `-profile preset_*`
configurations so a user who doesn't know which preset to choose still gets sensible
parameters. Decision policy (first match wins):

  query_len < short_peptide_max          -> preset_short_peptide
  n_hits > paralog_hits  OR  n_strong_loci > tandem_max_strong
                                          -> preset_paralog_discrimination  (+WARN)
  n_strong_loci >= 2  and  bit_ratio < single_copy_min_ratio
                                          -> preset_tandem_family
  otherwise                               -> preset_single_copy

Xaver's luciferase run (185 hits, 50 loci, bit_ratio 1.17) lands on
preset_paralog_discrimination + WARN — which with the §1d locus cap would have cut
the 55 h run to a few hours and surfaced the real ortholog calls.

Prints the chosen preset + reasoning; the workflow can echo this and let the user
override with --preset_override <profile>.
"""

import argparse
import json
import sys

VALID_PRESETS = {
    "preset_short_peptide",
    "preset_tandem_family",
    "preset_single_copy",
    "preset_paralog_discrimination",
}


def select_preset(profile, *, short_peptide_max=100, paralog_hits=15,
                  tandem_max_strong=5, single_copy_min_ratio=2.0):
    """Return (preset_name, warn: bool, reasons: list[str])."""
    n_hits = profile.get("n_hits") or 0
    query_len = profile.get("query_len")
    n_strong = profile.get("n_strong_loci") or 0
    n_loci = profile.get("n_independent_loci") or 0
    ratio = profile.get("bit_ratio")
    reasons = []

    if n_hits == 0:
        reasons.append("no home-genome hits; defaulting to single-copy parameters")
        return "preset_single_copy", False, reasons

    if query_len is not None and query_len < short_peptide_max:
        reasons.append(
            f"query is {query_len} aa (< {short_peptide_max}) — short peptide, where "
            f"flanking synteny matters more than the (tiny) GOI signal"
        )
        return "preset_short_peptide", False, reasons

    if n_hits > paralog_hits or n_strong > tandem_max_strong:
        reasons.append(
            f"{n_hits} hits across {n_loci} independent loci ({n_strong} strong); "
            f"bit_ratio={ratio} (primary barely beats the next locus) — large paralog "
            f"family. Capping loci and tightening family discrimination avoids a runaway "
            f"multi-locus search."
        )
        return "preset_paralog_discrimination", True, reasons

    if n_strong >= 2 and (ratio is None or ratio < single_copy_min_ratio):
        reasons.append(
            f"{n_strong} strong loci with bit_ratio={ratio} (plateau, no single dominant "
            f"hit) — tandem / multi-copy family"
        )
        return "preset_tandem_family", False, reasons

    reasons.append(
        f"{n_strong} strong locus/loci, bit_ratio={ratio} (a single dominant primary) "
        f"over {n_hits} hit(s) — single-copy / conserved"
    )
    return "preset_single_copy", False, reasons


def main():
    ap = argparse.ArgumentParser(description="Auto-select a query-type preset from hit_profile.json (TODO §1f)")
    ap.add_argument("--hit_profile", required=True, help="hit_profile.json from profile_hits.py")
    ap.add_argument("--output", help="Optional JSON with the decision")
    ap.add_argument("--preset_override", default=None,
                    help="Force a preset (one of preset_short_peptide/tandem_family/single_copy/"
                         "paralog_discrimination); bypasses auto-selection but still records the profile")
    ap.add_argument("--short_peptide_max", type=int, default=100)
    ap.add_argument("--paralog_hits", type=int, default=15)
    ap.add_argument("--tandem_max_strong", type=int, default=5)
    ap.add_argument("--single_copy_min_ratio", type=float, default=2.0)
    args = ap.parse_args()

    with open(args.hit_profile) as fh:
        profile = json.load(fh)

    if args.preset_override:
        if args.preset_override not in VALID_PRESETS:
            print(f"ERROR: --preset_override '{args.preset_override}' is not one of {sorted(VALID_PRESETS)}",
                  file=sys.stderr)
            sys.exit(2)
        preset, warn, reasons = args.preset_override, False, [f"user override: {args.preset_override}"]
    else:
        preset, warn, reasons = select_preset(
            profile,
            short_peptide_max=args.short_peptide_max,
            paralog_hits=args.paralog_hits,
            tandem_max_strong=args.tandem_max_strong,
            single_copy_min_ratio=args.single_copy_min_ratio,
        )

    print(f"AUTO-PRESET: {preset}")
    for r in reasons:
        print(f"  reason: {r}")
    if warn:
        print(f"  WARNING: large paralog family detected — review regions.bed and consider "
              f"a more specific query or stricter --search_evalue / --max_loci.", file=sys.stderr)
    print(f"  override with: --preset_override <name>  (or -profile {preset} to confirm)")

    if args.output:
        with open(args.output, "w") as fh:
            json.dump({"preset": preset, "warn": warn, "reasons": reasons, "profile": profile}, fh, indent=2)


if __name__ == "__main__":
    main()
