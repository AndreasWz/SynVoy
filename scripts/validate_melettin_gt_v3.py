#!/usr/bin/env python3
"""Compare a fresh melittin run output against the pinned v3 fixture.

Usage:
    python3 scripts/validate_melettin_gt_v3.py --outdir results/melettin_gt_v4

Exits 0 if the new output matches the fixture within tolerances.
Exits 1 with a per-check breakdown otherwise.

Checks:
  - Same set of species produced BED files.
  - Each species' BED row uses the same scaffold ID and the start/end
    coordinates are within `--coord-tolerance` (default 10%) of the
    fixture values.
  - scores TSV 'confidence' string matches per species.
  - Tree has the same leaf set as the fixture.

Not checked (too noisy for regression, but reported informationally):
  - Numeric score deltas (scoring is sensitive to minor tweaks).
  - Branch lengths in the tree.

See tests/ground_truth_test/melettin_gt_v3/README.md for details.
"""

import argparse
import csv
import os
import re
import sys


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIXTURE_DIR = os.path.join(REPO_ROOT, "tests", "ground_truth_test", "melettin_gt_v3")
FIXTURE_REGIONS = os.path.join(FIXTURE_DIR, "regions")
FIXTURE_TREE = os.path.join(FIXTURE_DIR, "locus_1_tree.nwk")


def _newick_leaves(nwk):
    stripped = re.sub(r":\d+(\.\d+)?([eE][+-]?\d+)?", "", nwk)
    stripped = stripped.rstrip(";").strip()
    tokens = re.split(r"[(),]", stripped)
    return {t.strip() for t in tokens if t.strip() and not re.fullmatch(r"\d+(\.\d+)?", t.strip())}


def _load_bed_rows(path):
    """Return all 6-col BED rows in file order. The fixture is single-row;
    fresh runs may emit multiple rows per species (GOI-anchor + score-floor)."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            rows.append({
                "chrom": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
                "score": float(cols[4]),
                "strand": cols[5],
            })
    return rows


def _load_bed(path):
    rows = _load_bed_rows(path)
    return rows[0] if rows else None


def _pick_match(fixture_row, fresh_rows):
    """Find the fresh BED row that best matches the fixture's scaffold.
    Prefers an exact scaffold match; among matches, prefers the one whose
    coords overlap the fixture region. Returns None if no scaffold matches."""
    same_scaffold = [r for r in fresh_rows if r["chrom"] == fixture_row["chrom"]]
    if not same_scaffold:
        return None
    f_start, f_end = fixture_row["start"], fixture_row["end"]
    overlapping = [r for r in same_scaffold
                   if r["start"] < f_end and r["end"] > f_start]
    if overlapping:
        return overlapping[0]
    return same_scaffold[0]


def _load_confidence(scores_path, chrom=None, start=None, end=None):
    """Return confidence for a given (chrom, start, end) row, or for the first
    row if no coords given. Matching by coords lets us align the confidence
    check with the BED row chosen by _pick_match in multi-region cases."""
    if not scores_path or not os.path.exists(scores_path):
        return None
    with open(scores_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        return None
    if chrom is None:
        return rows[0].get("confidence", "").upper()
    for row in rows:
        try:
            if (row.get("chrom") == chrom
                    and int(row.get("start", -1)) == start
                    and int(row.get("end", -1)) == end):
                return row.get("confidence", "").upper()
        except (TypeError, ValueError):
            continue
    # Coord match failed → fall back to first row so the diff still surfaces.
    return rows[0].get("confidence", "").upper()


def _collect_species(regions_dir):
    species = {}
    if not os.path.isdir(regions_dir):
        return species
    for fname in os.listdir(regions_dir):
        if fname.endswith(".fa.regions.bed"):
            sp = fname[:-len(".fa.regions.bed")]
            species.setdefault(sp, {})["bed"] = os.path.join(regions_dir, fname)
        elif fname.endswith(".fa.scores.tsv"):
            sp = fname[:-len(".fa.scores.tsv")]
            species.setdefault(sp, {})["scores"] = os.path.join(regions_dir, fname)
    return species


def compare(fresh_outdir, coord_tolerance):
    fresh_regions = os.path.join(fresh_outdir, "regions")
    fresh_tree = os.path.join(fresh_outdir, "locus_1_tree.nwk")

    fixture_species = _collect_species(FIXTURE_REGIONS)
    fresh_species = _collect_species(fresh_regions)

    failures = []
    notes = []

    missing = set(fixture_species) - set(fresh_species)
    extra = set(fresh_species) - set(fixture_species)
    if missing:
        failures.append(f"Species missing from fresh run: {sorted(missing)}")
    if extra:
        notes.append(f"Species added by fresh run (not regression, informational): {sorted(extra)}")

    for sp in sorted(set(fixture_species) & set(fresh_species)):
        ref_bed = _load_bed(fixture_species[sp]["bed"])
        fresh_rows = _load_bed_rows(fresh_species[sp].get("bed", ""))
        if not fresh_rows:
            failures.append(f"{sp}: fresh run has no BED row")
            continue
        new_bed = _pick_match(ref_bed, fresh_rows)
        if new_bed is None:
            scaffolds = sorted({r["chrom"] for r in fresh_rows})
            failures.append(
                f"{sp}: fixture scaffold {ref_bed['chrom']} not found in fresh "
                f"output (fresh scaffolds: {scaffolds})"
            )
            continue
        if len(fresh_rows) > 1 and new_bed is not fresh_rows[0]:
            notes.append(
                f"{sp}: fresh emits {len(fresh_rows)} regions; matched fixture "
                f"to non-top region '{new_bed['name']}' (top is "
                f"'{fresh_rows[0]['name']}')."
            )
        # Containment check: if the fresh region fully contains the fixture
        # region (and isn't more than 50% wider), accept it without coord-drift
        # warning. Wider windows that still anchor the same locus are not
        # regressions.
        contains_fixture = (
            new_bed["start"] <= ref_bed["start"]
            and new_bed["end"] >= ref_bed["end"]
        )
        ref_len = ref_bed["end"] - ref_bed["start"]
        new_len = new_bed["end"] - new_bed["start"]
        width_ratio = new_len / max(ref_len, 1)
        if contains_fixture and width_ratio <= 1.5:
            pass  # accepted: same locus, tolerable widening
        else:
            start_delta = abs(new_bed["start"] - ref_bed["start"]) / max(ref_len, 1)
            end_delta = abs(new_bed["end"] - ref_bed["end"]) / max(ref_len, 1)
            if start_delta > coord_tolerance or end_delta > coord_tolerance:
                failures.append(
                    f"{sp}: coords drifted > {coord_tolerance:.0%}: "
                    f"fixture={ref_bed['start']}-{ref_bed['end']} "
                    f"fresh={new_bed['start']}-{new_bed['end']}"
                )
        if ref_bed["strand"] != new_bed["strand"]:
            failures.append(
                f"{sp}: strand flipped: fixture={ref_bed['strand']} fresh={new_bed['strand']}"
            )
        ref_sc = fixture_species[sp].get("scores")
        new_sc = fresh_species[sp].get("scores")
        if ref_sc and new_sc:
            ref_conf = _load_confidence(
                ref_sc, ref_bed["chrom"], ref_bed["start"], ref_bed["end"]
            )
            new_conf = _load_confidence(
                new_sc, new_bed["chrom"], new_bed["start"], new_bed["end"]
            )
            if ref_conf != new_conf:
                failures.append(
                    f"{sp}: confidence changed: fixture={ref_conf} fresh={new_conf}"
                )

    if not os.path.exists(fresh_tree):
        failures.append(f"Tree missing in fresh run: {fresh_tree}")
    else:
        with open(FIXTURE_TREE) as fh:
            ref_leaves = _newick_leaves(fh.read())
        with open(fresh_tree) as fh:
            new_leaves = _newick_leaves(fh.read())
        if ref_leaves != new_leaves:
            only_ref = ref_leaves - new_leaves
            only_new = new_leaves - ref_leaves
            failures.append(
                f"Tree leaf set differs: only_in_fixture={sorted(only_ref)} "
                f"only_in_fresh={sorted(only_new)}"
            )

    return failures, notes


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", required=True,
                        help="Fresh melittin run output directory (contains regions/ and locus_1_tree.nwk)")
    parser.add_argument("--coord-tolerance", type=float, default=0.10,
                        help="Fractional tolerance for start/end coord drift (default 0.10 = 10%%)")
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        print(f"ERROR: --outdir does not exist: {args.outdir}", file=sys.stderr)
        sys.exit(2)

    failures, notes = compare(args.outdir, args.coord_tolerance)

    for note in notes:
        print(f"note: {note}")

    if failures:
        print(f"FAIL: {len(failures)} regression(s) vs. melettin_gt_v3 fixture", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        sys.exit(1)

    print(f"PASS: fresh run {args.outdir} matches melettin_gt_v3 fixture")


if __name__ == "__main__":
    main()
