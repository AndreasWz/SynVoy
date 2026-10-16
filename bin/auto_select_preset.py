#!/usr/bin/env python3
"""Pick a query-type preset from a hit_profile.json (docs/TODO.md §1f).

Maps the LOCATE_GENE hit distribution onto one of the bundled `-profile preset_*`
configurations so a user who doesn't know which preset to choose still gets sensible
parameters. Decision policy (first match wins):

  query_len < short_peptide_max          -> preset_short_peptide
  n_strong_loci >= 2  and  bit_ratio < single_copy_min_ratio
                                          -> preset_paralog_discrimination
                                             (relatives at other places of the genome;
                                             a strong locus scores >= half of the top one,
                                             so the ratio is <= 2 and the rule is in
                                             effect n_strong_loci >= 2)
  primary_copies >= tandem_min_copies    -> preset_tandem_family
                                             (several copies next to each other)
  otherwise                               -> preset_single_copy

A WARNING is added when the query aligns at many places (more than `many_loci`, or
more than `tandem_max_strong` strong loci): a large family.

The decision counts places (loci, copies), not alignment rows. A row is one exon in one
search program; until 2026-10-08 more than 15 rows meant "large paralog family", which
was true of every single-copy gene with eight exons or more, and a tandem array was
never recognised because all its copies fall into one locus. Xaver's luciferase run
(50 loci, four of them strong) still lands on preset_paralog_discrimination + WARN.

Prints the chosen preset + reasoning; the workflow can echo this and let the user
override with --preset_override <profile>. Under --modules legacy (the default) this
choice is the run's preset; under --modules auto the preset modules
(select_modules.py) decide instead and this choice is only recorded.
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


def select_preset(profile, *, short_peptide_max=100, tandem_max_strong=5,
                  single_copy_min_ratio=2.0, tandem_min_copies=1.5, many_loci=10):
    """Return (preset_name, warn: bool, reasons: list[str])."""
    n_hits = profile.get("n_hits") or 0
    query_len = profile.get("query_len")
    n_strong = profile.get("n_strong_loci") or 0
    n_loci = profile.get("n_independent_loci") or 0
    ratio = profile.get("bit_ratio")
    copies = profile.get("primary_copies")
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

    large_family = n_strong > tandem_max_strong or n_loci > many_loci

    if n_strong >= 2 and (ratio is None or ratio < single_copy_min_ratio):
        reasons.append(
            f"{n_strong} strong loci among {n_loci} (bit_ratio={ratio}: the next locus scores "
            f"more than half of the primary) — relatives of the query at other places of "
            f"the genome, which have to be told apart from it"
        )
        return "preset_paralog_discrimination", large_family, reasons

    if copies is not None and copies >= tandem_min_copies:
        reasons.append(
            f"the query is covered {copies} times at its strongest locus — several copies "
            f"next to each other (tandem array)"
        )
        return "preset_tandem_family", large_family, reasons

    reasons.append(
        f"{n_strong} strong locus/loci among {n_loci}, bit_ratio={ratio}, query covered "
        f"{copies} time(s) at the strongest — single-copy / conserved"
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
    ap.add_argument("--tandem_max_strong", type=int, default=5)
    ap.add_argument("--single_copy_min_ratio", type=float, default=2.0)
    ap.add_argument("--tandem_min_copies", type=float, default=1.5,
                    help="Query coverage depth at the strongest locus from which it counts "
                         "as a tandem array (default 1.5)")
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
            tandem_max_strong=args.tandem_max_strong,
            single_copy_min_ratio=args.single_copy_min_ratio,
            tandem_min_copies=args.tandem_min_copies,
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
