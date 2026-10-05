#!/usr/bin/env python3
"""Phylogenetic placement check — is a recovered GOI call where orthology predicts?

docs/STATE_OF_THE_PROJECT.md F9. Three earlier attempts to curb the over-calling all
failed, and they failed for the same reason: identity, coverage, flanking support and
home-paralog RBH are all evaluated **one call at a time**, and one call in isolation
carries no information about whether it is the ortholog or a lineage-specific paralog.
The measured proof is that the highest-identity call in the melittin benchmark (93.8 %,
*Bombus terrestris*, where melittin is known lost) outranks every true positive.

This check asks a question no single call can answer:

    does sequence divergence track SPECIES divergence?

Under orthology it must — the gene and the species diverged together. A gene that
duplicated *before* the speciation (a paralog retained while the true ortholog was lost,
which is exactly the bombolitin/melittin case in *Bombus*) sits anomalously FAR from the
query for how CLOSE its species is. That discordance is visible only across the whole set
of calls, which is why it is computed here rather than in the per-call classifier.

Method
------
For every recovered GOI protein:

    gene_divergence   d_g = 1 - identity(query, protein)
    species_distance  d_s = taxonomic rank distance (PHYLO_SORT's sorted_genomes.txt)

Fit d_g ~ d_s with a **Theil-Sen** estimator (median of pairwise slopes) — robust to the
very outliers we are hunting, and dependency-free. Residuals are scaled by their MAD to
give a z-score; a call whose divergence is far above the fit is `phylo_discordant`.

Deliberately NOT a tree. With a 70 aa peptide at ~30 % identity a gene tree has no
resolvable topology and its bootstrap support would be noise. A rate/distance concordance
test extracts the orthology signal that IS present without pretending to a topology that
is not. `compute_tree.py` still builds the display tree.

Output is advisory by default: it emits verdicts, and `generate_report.py` decides what to
do with them.
"""
import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import parse_fasta, setup_logging, str2bool, sw_align  # noqa: E402

logger = setup_logging(name="phylo_placement_check")

# GOI model ids carry the genome stem plus a suffix that differs per emitting path.
# Observed forms (enumerated from real run GFFs, 2026-09-09):
#     GOI_Melt|Colletes_gigas_fa_b0_l1_fallback      seeded block, fallback span
#     GOI_Melt|Apis_florea_fna_b0_l1_exon_ann        seeded block, exon model
#     GOI_copy_1|Apis_cerana_fna_b0_l1               tandem copy  <-- NO trailing token
#     GOI_Melt|Colletes_gigas_WUUM01000001.1_hull_rescue   hull rescue (chrom, not block)
#     GOI_rescue_<genome>_locus_<name>_<idx>         strong-synteny rescue (no pipe)
#
# The first regex here required a token AFTER the locus id, so the bare `_b0_l1` form
# never matched and those calls resolved to a species that does not exist in
# sorted_genomes.txt -> no distance -> `insufficient_data`. In the first F9 cluster run
# that silently killed 7 of 8 calls and the check reported nothing at all.
#
# Parsing arbitrary suffixes is the wrong shape of solution: the genome stems are KNOWN
# (they are the keys of sorted_genomes.txt), so resolve by longest known-stem prefix and
# keep the regex only as a fallback for callers that pass no stem list.
_ID_SUFFIX = re.compile(r"_b\d+_(?:l\d+|fl\d+)(?:_.*)?$")
_EXT = re.compile(r"[._](?:fna|fa|fasta)$")
_RESCUE_PREFIX = re.compile(r"^GOI_rescue_")


def _normalize_stem(text):
    """Model ids sanitize '.' to '_' (GCF_002204515.2.fna -> GCF_002204515_2_fna), so
    compare stems in a single normalized space."""
    return text.replace(".", "_")


def species_from_model_id(model_id, known_stems=None):
    """Recover the genome stem from a GOI model id, or '' if it cannot be resolved.

    `known_stems` (the genome stems from sorted_genomes.txt) makes this exact for every
    emitting path; without it the regex fallback handles the seeded-block forms only.
    """
    core = model_id.split("|", 1)[-1].strip()
    core = _RESCUE_PREFIX.sub("", core)

    if known_stems:
        norm_core = _normalize_stem(core)
        best = ""
        for stem in known_stems:
            norm_stem = _normalize_stem(stem)
            if not norm_stem:
                continue
            # A stem matches when the id IS it, or continues with a separator — so
            # "Apis_cerana" never swallows an id belonging to "Apis_cerana_2".
            if norm_core == norm_stem or norm_core.startswith(norm_stem + "_"):
                if len(norm_stem) > len(best):
                    best = norm_stem
        if best:
            # Return the caller's own spelling of the stem so distances.get() hits.
            for stem in known_stems:
                if _normalize_stem(stem) == best:
                    return stem

    core = _ID_SUFFIX.sub("", core)
    core = _EXT.sub("", core)
    return core


def load_species_distances(path):
    """`sorted_genomes.txt` is '<genome file>\\t<taxonomic distance>' per line."""
    dist = {}
    if not path or not os.path.isfile(path):
        return dist
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            stem = _EXT.sub("", re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", parts[0].strip()))
            try:
                dist[stem] = float(parts[1])
            except ValueError:
                continue
    return dist


def theil_sen(xs, ys):
    """Median of pairwise slopes, plus the median intercept it implies.

    Robust to the outliers this check exists to find: a least-squares fit would be
    dragged toward a discordant call and hide it.
    """
    slopes = []
    n = len(xs)
    for i in range(n):
        for j in range(i + 1, n):
            dx = xs[j] - xs[i]
            if dx != 0:
                slopes.append((ys[j] - ys[i]) / dx)
    if not slopes:
        return None, None
    slope = _median(slopes)
    intercept = _median([y - slope * x for x, y in zip(xs, ys)])
    return slope, intercept


def _median(vals):
    v = sorted(vals)
    n = len(v)
    if not n:
        return 0.0
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def _mad(vals, centre):
    """Median absolute deviation, scaled to be a consistent sigma estimator."""
    if not vals:
        return 0.0
    return 1.4826 * _median([abs(v - centre) for v in vals])


def evaluate(records, min_calls=4, z_threshold=2.0, min_distance_levels=3):
    """Score every call for phylogenetic concordance.

    ``records`` are dicts with species / species_distance / gene_divergence. Returns the
    same list with fit columns added, plus a summary. With fewer than ``min_calls``
    species, or no spread in species distance, no fit is attempted and everything is
    reported ``insufficient_data`` — silence is the correct output when the data cannot
    support the test.
    """
    usable = [r for r in records if r["species_distance"] is not None]
    xs = [r["species_distance"] for r in usable]
    ys = [r["gene_divergence"] for r in usable]
    summary = {"n_calls": len(records), "n_fitted": 0, "n_concordant": 0,
               "n_discordant": 0, "n_insufficient": 0,
               "slope": "", "intercept": "", "resid_sigma": ""}

    # Report the spread so a degenerate distance set is visible rather than inferred.
    # PHYLO_SORT emits near-constant fallback distances (e.g. 993-1000, or a flat 999)
    # when the taxonomy walk cannot separate the targets; fitting a slope through that
    # is fitting noise, so require real levels before attempting it.
    summary["distance_levels"] = len(set(xs))
    summary["distance_range"] = (f"{min(xs):g}-{max(xs):g}" if xs else "")

    if len(usable) < min_calls or len(set(xs)) < max(2, min_distance_levels):
        for r in records:
            r.update(expected_divergence="", residual="", z_score="",
                     verdict="insufficient_data")
        summary["n_insufficient"] = len(records)
        return records, summary

    slope, intercept = theil_sen(xs, ys)
    if slope is None:
        for r in records:
            r.update(expected_divergence="", residual="", z_score="",
                     verdict="insufficient_data")
        summary["n_insufficient"] = len(records)
        return records, summary

    resids = [y - (slope * x + intercept) for x, y in zip(xs, ys)]
    centre = _median(resids)
    sigma = _mad(resids, centre)
    summary.update(slope=round(slope, 6), intercept=round(intercept, 6),
                   resid_sigma=round(sigma, 6), n_fitted=len(usable))

    for r in records:
        if r["species_distance"] is None:
            r.update(expected_divergence="", residual="", z_score="",
                     verdict="insufficient_data")
            summary["n_insufficient"] += 1
            continue
        exp = slope * r["species_distance"] + intercept
        resid = r["gene_divergence"] - exp
        # A flat residual distribution (sigma 0) means every call fits the trend
        # exactly; treat that as concordant rather than dividing by zero.
        z = 0.0 if sigma <= 0 else (resid - centre) / sigma
        # ONE-SIDED on purpose: only "more diverged than its species predicts" is
        # evidence of a paralog. A call closer to the query than expected is a
        # well-conserved ortholog, not a problem.
        verdict = "phylo_discordant" if z >= z_threshold else "phylo_concordant"
        r.update(expected_divergence=round(exp, 4), residual=round(resid, 4),
                 z_score=round(z, 3), verdict=verdict)
        summary["n_discordant" if verdict == "phylo_discordant" else "n_concordant"] += 1
    return records, summary


def build_records(query_seq, goi_faa, distances, min_target_len=20):
    """One record per recovered GOI protein, with its species distance attached.

    `min_target_len` drops fragments too short to carry a divergence estimate. The first
    F9 run fitted on a 7-aa "protein" at 40 % identity — over 7 residues that number is
    noise, and one such point can set the slope for the whole genome set.
    """
    records = []
    stems = list(distances.keys())
    unresolved, too_short = [], []
    for _header, model_id, seq in parse_fasta(goi_faa):
        if not seq:
            continue
        if len(seq) < min_target_len:
            too_short.append((model_id, len(seq)))
            continue
        sp = species_from_model_id(model_id, known_stems=stems)
        dist = distances.get(sp)
        if dist is None:
            unresolved.append(model_id)
        score, ident = sw_align(query_seq, seq)
        records.append({
            "model_id": model_id, "species": sp,
            "species_distance": dist,
            "identity": round(ident, 2), "sw_score": round(score, 1),
            "target_len_aa": len(seq),
            "gene_divergence": round(1.0 - (ident / 100.0), 4),
        })

    if too_short:
        logger.info("skipped %d call(s) shorter than %d aa (e.g. %s)",
                    len(too_short), min_target_len,
                    ", ".join(f"{m} ({n} aa)" for m, n in too_short[:3]))
    # Fail LOUD on species resolution: a silent miss here degrades every verdict to
    # `insufficient_data` while the run still reports success (the 2026-09-09 F9 bug).
    if unresolved:
        logger.warning(
            "%d/%d call(s) could not be matched to a genome in sorted_genomes.txt and "
            "will have NO species distance — verdicts will be weaker or absent. "
            "Unmatched ids: %s | known stems: %s",
            len(unresolved), len(records), ", ".join(unresolved[:5]),
            ", ".join(sorted(stems)[:8]))
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--goi_faa", required=True,
                    help="goi_for_tree.faa from ITERATIVE_SEARCH (recovered GOI proteins)")
    ap.add_argument("--query", required=True, help="normalized GOI query protein FASTA")
    ap.add_argument("--sorted_genomes", required=True,
                    help="PHYLO_SORT sorted_genomes.txt (<genome>\\t<distance>)")
    ap.add_argument("--output", required=True, help="per-call TSV")
    ap.add_argument("--min_calls", type=int, default=4,
                    help="minimum species with calls before a fit is attempted")
    ap.add_argument("--min_target_len", type=int, default=20,
                    help="skip recovered proteins shorter than this (aa); a 7-aa "
                         "fragment carries no usable divergence estimate")
    ap.add_argument("--min_distance_levels", type=int, default=3,
                    help="minimum distinct species-distance values before a fit is "
                         "attempted (guards PHYLO_SORT's near-constant fallbacks)")
    ap.add_argument("--z_threshold", type=float, default=2.0,
                    help="residual z at/above which a call is phylo_discordant")
    ap.add_argument("--disable_phylo_placement", type=str2bool, default=False)
    args = ap.parse_args()

    cols = ["model_id", "species", "species_distance", "identity", "sw_score",
            "target_len_aa", "gene_divergence", "expected_divergence", "residual",
            "z_score", "verdict"]

    def write(rows):
        with open(args.output, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    if args.disable_phylo_placement:
        logger.info("phylo placement check disabled; emitting header-only TSV")
        write([])
        return 0

    queries = [s for _h, _i, s in parse_fasta(args.query) if s]
    if not queries:
        logger.warning("empty query FASTA %s — nothing to do", args.query)
        write([])
        return 0
    # Multi-FASTA queries are paralog panels; the first record is the representative,
    # matching how the rest of the pipeline treats them.
    query_seq = queries[0]

    if not os.path.isfile(args.goi_faa) or os.path.getsize(args.goi_faa) == 0:
        logger.info("no recovered GOI proteins (%s) — nothing to place", args.goi_faa)
        write([])
        return 0

    distances = load_species_distances(args.sorted_genomes)
    records = build_records(query_seq, args.goi_faa, distances,
                           min_target_len=args.min_target_len)
    if not records:
        logger.info("no parsable GOI proteins in %s", args.goi_faa)
        write([])
        return 0

    records, summary = evaluate(records, args.min_calls, args.z_threshold,
                                min_distance_levels=args.min_distance_levels)
    write(records)

    logger.info(
        "phylo placement: %d call(s), %d fitted — %d concordant, %d discordant, "
        "%d insufficient (slope=%s, sigma=%s)",
        summary["n_calls"], summary["n_fitted"], summary["n_concordant"],
        summary["n_discordant"], summary["n_insufficient"],
        summary["slope"], summary["resid_sigma"])
    for r in records:
        if r.get("verdict") == "phylo_discordant":
            logger.info(
                "  DISCORDANT %s (%s): divergence %.3f vs %.3f expected at species "
                "distance %s (z=%.2f) — more diverged than orthology predicts",
                r["species"], r["model_id"], r["gene_divergence"],
                r["expected_divergence"], r["species_distance"], r["z_score"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
