#!/usr/bin/env python3
"""Score SynVoy calls against the melittin truth set by COORDINATE OVERLAP (§1y).

Why this exists
---------------
`score_benchmark.py` scores **species presence**: did the tool say PRESENT for a
species whose truth is PRESENT. That is the number in `docs/STATE_OF_THE_PROJECT.md`
D.1 (SynVoy F1 0.700, recall 0.875), and it is the number that let two overclaims
through:

  * The *T. bicarinatum* melittin "recovery" was **132 bp of an 855 bp gene** -- one
    exon, `ModelStatus=partial`. Presence scoring records that as a clean TP.
  * The 2026-06-30 run **lost that hit entirely** (its nearest call is 2,543 bp away
    from the nearest truth model) and still scored PRESENT for the species, because
    some other call on the same scaffold kept the species flag up.

So a call is scored here by *where it landed*, and the overlap FRACTION is reported
rather than a boolean, because the fraction is what distinguishes "recovered the
gene" from "landed on its first exon".

Truth kinds
-----------
`gene` / `fragment`     a CURATED model with real exon coordinates. Hit/miss AND
                        overlap fraction are both meaningful.
`locus_gene_present`    the gene is known present but only the syntenic REGION is
                        curated (Colletes' span is 110 kb). Hit/miss is meaningful;
                        an overlap fraction against a region is not, so these are
                        excluded from the fraction statistics.
`locus_gene_lost`       the syntenic locus is present but the GENE IS LOST (Melipona,
                        Tetragonula). A call here is a FALSE POSITIVE -- the deliberate
                        false-positive test, which presence scoring cannot express.

A truth model whose chromosome never appears anywhere in the run is reported
`not_searched` and excluded from the denominators: a species that was not a target
is not a miss, and silently counting it as one understates recall.

Usage
-----
    score_coordinates.py --truth tests/benchmark_truth/melittin_loci.tsv \\
                         --run   local_runs/mel_det_rep2 \\
                         [--min-confidence MEDIUM] [--out report.tsv]

`--run` may be a SynVoy output directory (every `*.gff` under it is scanned for
`SynVoyRole=goi` mRNA features) or a single GFF file.
"""
import argparse
import csv
import glob
import os
import re
import sys

CONF_RANK = {"LOW": 0, "MEDIUM": 1, "HIGH": 2}


# ---------------------------------------------------------------- loading ---
def load_truth(path):
    """Read the coordinate truth set. Returns a list of model dicts.

    `chrom_aliases` holds equivalent accessions for the SAME physical sequence -- a
    GenBank/RefSeq mirror pair such as CM079762.1 == NC_091864.1, verified by sequence
    length. Without it, curating on one accession and running on the other made every
    Cardiocondyla model unscoreable. `chrom_keys` is the set the scorer matches against.
    """
    models = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["start"] = int(r["start"])
            r["end"] = int(r["end"])
            aliases = [a for a in (r.get("chrom_aliases") or "").split(";") if a]
            r["chrom_keys"] = {r["chrom"], *aliases}
            models.append(r)
    if not models:
        sys.exit(f"no truth models in {path}")
    return models


def _attrs(field):
    out = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def load_searched_chroms(paths):
    """Every seqid appearing in the run's GFFs, GOI or not.

    NOTE this is NOT proof of what was searched: SynVoy only emits features for
    scaffolds where it found something, so a target scaffold with no hits leaves no
    trace here. Used only to tell an assembly mismatch from a genuine miss -- see
    classify_scoreability.
    """
    chroms = set()
    for p in paths:
        try:
            fh = open(p)
        except OSError:
            continue
        with fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t", 1)
                if f[0]:
                    chroms.add(f[0])
    return chroms


def load_searched_species(run_dir):
    """Target species of a run, from `regions/<species>.<ext>.regions.bed`.

    This is the authoritative list -- every target gets a regions file whether or
    not anything was found in it, so a species here with no hit is a REAL miss.
    Returns None when the layout is absent (single-GFF input), which disables the
    not_searched/assembly_mismatch distinction rather than guessing.
    """
    regions = os.path.join(run_dir, "regions") if os.path.isdir(run_dir) else None
    if not regions or not os.path.isdir(regions):
        return None
    species = set()
    for name in os.listdir(regions):
        m = re.match(r"^(.+?)\.(?:fna|fa|fasta)\.(?:regions\.bed|scores\.tsv)$", name)
        if m:
            species.add(m.group(1))
    return species or None


def classify_scoreability(model, searched_species, truth_by_species, searched_chroms):
    """Can this truth model be scored against this run?

    'not_searched'      the species was not a target -- not a miss.
    'assembly_mismatch' the species WAS a target but none of its curated scaffolds
                        appear in the run, i.e. the run used a different assembly
                        (Cardiocondyla is curated on GCA_019399895.2 while the truth
                        table names GCF_019399895.1). Coordinates are not comparable;
                        scoring it as a miss would be a fabricated failure.
    'scoreable'         otherwise -- a miss here is a real miss.
    """
    if searched_species is None:
        return "scoreable"
    if model["species"] not in searched_species:
        return "not_searched"
    sp_chroms = truth_by_species.get(model["species"], set())
    if sp_chroms and not (sp_chroms & searched_chroms):
        return "assembly_mismatch"
    return "scoreable"


def load_calls_from_gff(paths):
    """Extract GOI model features from SynVoy GFFs.

    Only `SynVoyRole=goi` models count: flanking-gene models and CDS children are not
    ortholog calls and must not be scored as if they were. A model is an `mRNA` line
    OR a `gene` line -- tandem copies are written as `gene` only
    (iterative_search_runner.py, `tandem_copy`), and reading mRNA alone scored LRZ
    det_rep R2 as having LOST the T. bicarinatum melittin when it had modelled the
    same locus as a tandem copy. A `gene` that is the Parent of an mRNA is skipped so a
    gene+mRNA pair is never counted twice.
    """
    calls, seen = [], set()
    for p in paths:
        try:
            fh = open(p)
        except OSError:
            continue
        with fh:
            rows = [line.rstrip("\n").split("\t") for line in fh
                    if not line.startswith("#") and line.strip()]
        parents = {_attrs(f[8]).get("Parent") for f in rows
                   if len(f) >= 9 and f[2] == "mRNA"}
        for f in rows:
            if len(f) < 9 or f[2] not in ("mRNA", "gene"):
                continue
            a = _attrs(f[8])
            if a.get("SynVoyRole") != "goi":
                continue
            if f[2] == "gene" and a.get("ID") in parents:
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            key = (f[0], start, end, a.get("ID", ""))
            if key in seen:          # the same model is staged under several dirs
                continue
            seen.add(key)
            calls.append({
                "chrom": f[0], "start": start, "end": end, "strand": f[6],
                "identity": a.get("Identity", f[5]),
                "confidence": a.get("Confidence", "").upper(),
                "goi_class": a.get("GOIClass", ""),
                "query_cov": a.get("QueryCoverage", ""),
                "exons": a.get("Exons", ""),
                "model_status": a.get("ModelStatus", ""),
                "id": a.get("ID", ""),
                "src": os.path.basename(p),
            })
    return calls


def gather_gffs(target):
    if os.path.isfile(target):
        return [target]
    if not os.path.isdir(target):
        sys.exit(f"--run not found: {target}")
    return sorted(glob.glob(os.path.join(target, "**", "*.gff*"), recursive=True))


# ---------------------------------------------------------------- scoring ---
def overlap_bp(a_start, a_end, b_start, b_end):
    """Inclusive 1-based overlap length; 0 when disjoint."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def distance_bp(a_start, a_end, b_start, b_end):
    """Gap between two disjoint intervals; 0 if they overlap."""
    if overlap_bp(a_start, a_end, b_start, b_end):
        return 0
    return b_start - a_end if b_start > a_end else a_start - b_end


def score(truth, calls, min_confidence=None, searched_chroms=None, exclude_species=(),
          searched_species=None):
    """Match calls to truth models by coordinate overlap.

    Returns (model_rows, call_rows, summary). Every truth model gets a row even if
    nothing hit it -- a silent miss is the failure mode §1y is about.
    """
    if min_confidence:
        floor = CONF_RANK.get(min_confidence.upper(), 0)
        calls = [c for c in calls
                 if CONF_RANK.get(c["confidence"], -1) >= floor]
    if exclude_species:
        truth = [m for m in truth if m["species"] not in set(exclude_species)]

    truth_by_species = {}
    for m in truth:
        truth_by_species.setdefault(m["species"], set()).update(
            m.get("chrom_keys") or {m["chrom"]})
    searched_chroms = searched_chroms or set()

    model_rows = []
    hit_call_ids = set()
    for m in truth:
        cands = []
        keys = m.get("chrom_keys") or {m["chrom"]}
        for i, c in enumerate(calls):
            if c["chrom"] not in keys:
                continue
            ov = overlap_bp(m["start"], m["end"], c["start"], c["end"])
            if ov > 0:
                cands.append((ov, i, c))
        cands.sort(key=lambda t: -t[0])
        span = m["end"] - m["start"] + 1
        state = classify_scoreability(m, searched_species, truth_by_species,
                                      searched_chroms)
        if cands:
            ov, i, c = cands[0]
            hit_call_ids.update(j for _, j, _ in cands)
            model_rows.append({
                **{k: m[k] for k in ("species", "chrom", "start", "end", "model_id",
                                     "truth_kind", "confidence")},
                "truth_span": span, "hit": "yes", "n_calls_overlapping": len(cands),
                "best_call": f"{c['chrom']}:{c['start']}-{c['end']}",
                "call_confidence": c["confidence"], "call_identity": c["identity"],
                "overlap_bp": ov, "overlap_frac_of_truth": round(ov / span, 4),
                "start_offset": c["start"] - m["start"],
                "nearest_call_distance": 0,
            })
        else:
            same = [c for c in calls if c["chrom"] in keys]
            nd = min((distance_bp(m["start"], m["end"], c["start"], c["end"])
                      for c in same), default="")
            model_rows.append({
                **{k: m[k] for k in ("species", "chrom", "start", "end", "model_id",
                                     "truth_kind", "confidence")},
                "truth_span": span,
                "hit": ("no" if state == "scoreable" else state),
                "n_calls_overlapping": 0,
                "best_call": "", "call_confidence": "", "call_identity": "",
                "overlap_bp": 0, "overlap_frac_of_truth": 0.0, "start_offset": "",
                "nearest_call_distance": nd,
            })

    call_rows = []
    for i, c in enumerate(calls):
        on_truth = i in hit_call_ids
        lost = any(m["truth_kind"] == "locus_gene_lost"
                   and c["chrom"] in (m.get("chrom_keys") or {m["chrom"]})
                   and overlap_bp(m["start"], m["end"], c["start"], c["end"]) > 0
                   for m in truth)
        call_rows.append({
            "chrom": c["chrom"], "start": c["start"], "end": c["end"],
            "confidence": c["confidence"], "identity": c["identity"],
            "goi_class": c["goi_class"], "query_cov": c["query_cov"],
            "on_truth_model": "yes" if on_truth else "no",
            "on_lost_locus": "yes" if lost else "no",
        })

    UNSCOREABLE = ("not_searched", "assembly_mismatch")
    searchable = [m for m in model_rows if m["hit"] not in UNSCOREABLE]
    # Recall counts every truth model that a call COULD have hit; the fraction stats
    # count only curated gene/fragment models, where a fraction has a meaning.
    real = [m for m in searchable
            if m["truth_kind"] in ("gene", "fragment", "locus_gene_present")]
    curated = [m for m in searchable if m["truth_kind"] in ("gene", "fragment")]
    lost_models = [m for m in searchable if m["truth_kind"] == "locus_gene_lost"]
    hits = [m for m in real if m["hit"] == "yes"]
    curated_hits = [m for m in curated if m["hit"] == "yes"]
    fracs = sorted(m["overlap_frac_of_truth"] for m in curated_hits)
    summary = {
        "truth_models_not_searched": sum(1 for m in model_rows
                                         if m["hit"] == "not_searched"),
        "truth_models_assembly_mismatch": sum(1 for m in model_rows
                                              if m["hit"] == "assembly_mismatch"),
        "truth_models_real": len(real),
        "truth_models_hit": len(hits),
        "locus_recall": round(len(hits) / len(real), 4) if real else 0.0,
        "curated_models_searched": len(curated),
        "curated_models_hit": len(curated_hits),
        "median_overlap_frac": (fracs[len(fracs) // 2] if fracs else 0.0),
        "models_hit_over_50pct": sum(1 for m in curated_hits
                                     if m["overlap_frac_of_truth"] >= 0.5),
        "calls_total": len(call_rows),
        "calls_on_truth": sum(1 for c in call_rows if c["on_truth_model"] == "yes"),
        "locus_precision": (round(sum(1 for c in call_rows if c["on_truth_model"] == "yes")
                                  / len(call_rows), 4) if call_rows else 0.0),
        "lost_loci_total": len(lost_models),
        "false_positive_calls_on_lost_loci": sum(1 for c in call_rows
                                                 if c["on_lost_locus"] == "yes"),
    }
    return model_rows, call_rows, summary


# -------------------------------------------------------------------- cli ---
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default="tests/benchmark_truth/melittin_loci.tsv")
    ap.add_argument("--run", required=True,
                    help="SynVoy output directory, or a single GFF file")
    ap.add_argument("--min-confidence", choices=["LOW", "MEDIUM", "HIGH"],
                    help="only score calls at or above this confidence")
    ap.add_argument("--exclude-species", nargs="*", default=["Apis_mellifera"],
                    help="species to drop from the truth set (default: the home "
                         "genome, which is never a target)")
    ap.add_argument("--out", help="write the per-model TSV here")
    ap.add_argument("--out-calls", help="write the per-call TSV here")
    args = ap.parse_args()

    truth = load_truth(args.truth)
    gffs = gather_gffs(args.run)
    calls = load_calls_from_gff(gffs)
    searched = load_searched_chroms(gffs)
    sp = load_searched_species(args.run)
    model_rows, call_rows, summary = score(truth, calls, args.min_confidence,
                                           searched_chroms=searched,
                                           exclude_species=args.exclude_species,
                                           searched_species=sp)

    print(f"truth models : {len(truth)}   GOI calls parsed: {len(calls)}"
          + (f"   (confidence >= {args.min_confidence})" if args.min_confidence else ""))
    print("\n=== TRUTH MODELS")
    print(f"  {'species':<26}{'model':<26}{'kind':<16}{'hit':<5}{'overlap':>9}"
          f"{'frac':>8}{'offset':>9}{'nearest':>10}")
    for m in model_rows:
        print(f"  {m['species']:<26}{m['model_id']:<26}{m['truth_kind']:<16}"
              f"{m['hit']:<5}{m['overlap_bp']:>9}{m['overlap_frac_of_truth']:>8.3f}"
              f"{str(m['start_offset']):>9}{str(m['nearest_call_distance']):>10}")

    print("\n=== SUMMARY")
    for k, v in summary.items():
        print(f"  {k:<38}{v}")

    print("\n=== READ THIS NUMBER, NOT THE PRESENCE NUMBER")
    print(f"  locus recall      {summary['truth_models_hit']}/{summary['truth_models_real']} "
          f"searched truth models overlapped by a call")
    print(f"  substantial       {summary['models_hit_over_50pct']}/"
          f"{summary['curated_models_searched']} curated models recovered over "
          f">=50% of their span")
    print(f"  median overlap    {summary['median_overlap_frac']:.1%} of the curated model")
    if summary["truth_models_not_searched"]:
        print(f"  ({summary['truth_models_not_searched']} truth model(s) whose species was "
              f"not a target -- excluded, not counted as misses)")
    if summary["truth_models_assembly_mismatch"]:
        print(f"  ** {summary['truth_models_assembly_mismatch']} truth model(s) are curated "
              f"on a DIFFERENT ASSEMBLY than the run used -- excluded as unscoreable, "
              f"not as misses. Fix the truth table's accession to compare them.")
    if summary["false_positive_calls_on_lost_loci"]:
        print(f"  ** {summary['false_positive_calls_on_lost_loci']} call(s) land on a locus "
              f"where the gene is LOST -- these are false positives")

    for path, rows in ((args.out, model_rows), (args.out_calls, call_rows)):
        if path and rows:
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                   lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
