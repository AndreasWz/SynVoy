#!/usr/bin/env python3
"""Summarise the LOCATE_GENE hit distribution into a compact hit_profile.json.

The shape of the home-genome BLAST/MMseqs hit distribution tells you what kind of
query this is — single-copy, a tandem family, or a sprawling paralog superfamily —
before any expensive iterative search runs. auto_select_preset.py consumes this to
pick a query-type preset, so a user who doesn't know which `-profile preset_*` to
choose still gets sensible parameters (docs/TODO.md §1f).

Metrics emitted:
  query_len            query length in aa (from --query FASTA, else max qend across hits)
  n_hits               total raw HSP rows (BLAST + MMseqs combined)
  n_independent_loci   non-overlapping target loci (HSPs merged within --locus_window)
  n_strong_loci        independent loci whose best bit >= --strong_frac * top_bit
  top_bit, second_bit  best bit of the strongest / next-strongest independent locus
  bit_ratio            top_bit / second_bit  (a primary that barely beats the next
                       locus -> paralog-rich; one that dominates -> single-copy)
  identity_{min,q1,median,q3,max}, identity_iqr   percent-identity spread of all HSPs

Reads standard `-outfmt 6` / MMseqs `--format-output query,target,pident,alnlen,
mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits` (12 cols).
"""

import argparse
import json
import os
import sys


def _read_query_len(query_fasta):
    """Length (aa) of the first record in a FASTA, or None."""
    if not query_fasta or not os.path.exists(query_fasta):
        return None
    length = 0
    started = False
    with open(query_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            length += len(line.strip())
    return length or None


def parse_hits(paths):
    """Parse 12-column BLAST/MMseqs hit rows. Returns list of dicts."""
    hits = []
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#") or line.startswith("query\t"):
                    continue
                p = line.split("\t")
                if len(p) < 12:
                    continue
                try:
                    pident = float(p[2])
                    qend = int(p[7])
                    tstart, tend = int(p[8]), int(p[9])
                    bits = float(p[11])
                except ValueError:
                    continue
                if tend < tstart:
                    tstart, tend = tend, tstart
                hits.append({
                    "chrom": p[1], "qend": qend,
                    "tstart": tstart, "tend": tend,
                    "pident": pident, "bits": bits,
                })
    return hits


def _cluster_loci(hits, window):
    """Merge HSPs into independent target loci (same chrom, within `window` bp).
    Returns the best bit score per locus, sorted descending."""
    loci = []
    for h in sorted(hits, key=lambda x: (x["chrom"], x["tstart"])):
        if loci and loci[-1]["chrom"] == h["chrom"] and h["tstart"] <= loci[-1]["end"] + window:
            loci[-1]["end"] = max(loci[-1]["end"], h["tend"])
            loci[-1]["bit"] = max(loci[-1]["bit"], h["bits"])
        else:
            loci.append({"chrom": h["chrom"], "end": h["tend"], "bit": h["bits"]})
    return sorted((l["bit"] for l in loci), reverse=True)


def _quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * frac


def build_profile(hits, query_len, window, strong_frac):
    if not hits:
        return {
            "query_len": query_len, "n_hits": 0, "n_independent_loci": 0,
            "n_strong_loci": 0, "top_bit": None, "second_bit": None,
            "bit_ratio": None, "identity_min": None, "identity_q1": None,
            "identity_median": None, "identity_q3": None, "identity_max": None,
            "identity_iqr": None,
        }

    if query_len is None:
        query_len = max(h["qend"] for h in hits)

    locus_bits = _cluster_loci(hits, window)
    top_bit = locus_bits[0]
    second_bit = locus_bits[1] if len(locus_bits) > 1 else None
    bit_ratio = (top_bit / second_bit) if second_bit else None
    n_strong = sum(1 for b in locus_bits if b >= strong_frac * top_bit)

    idents = sorted(h["pident"] for h in hits)
    q1 = _quantile(idents, 0.25)
    q3 = _quantile(idents, 0.75)

    return {
        "query_len": query_len,
        "n_hits": len(hits),
        "n_independent_loci": len(locus_bits),
        "n_strong_loci": n_strong,
        "top_bit": round(top_bit, 1),
        "second_bit": round(second_bit, 1) if second_bit is not None else None,
        "bit_ratio": round(bit_ratio, 3) if bit_ratio is not None else None,
        "identity_min": round(idents[0], 1),
        "identity_q1": round(q1, 1),
        "identity_median": round(_quantile(idents, 0.5), 1),
        "identity_q3": round(q3, 1),
        "identity_max": round(idents[-1], 1),
        "identity_iqr": round(q3 - q1, 1),
    }


def main():
    ap = argparse.ArgumentParser(description="Profile LOCATE_GENE hits into hit_profile.json (TODO §1f)")
    ap.add_argument("--blast", help="BLAST -outfmt 6 hits (12 cols)")
    ap.add_argument("--mmseqs", help="MMseqs m8 hits (12 cols)")
    ap.add_argument("--query", help="Query protein FASTA (for query_len)")
    ap.add_argument("--locus_window", type=int, default=100000,
                    help="Merge HSPs within this many bp into one target locus (default 100000)")
    ap.add_argument("--strong_frac", type=float, default=0.5,
                    help="A locus counts as 'strong' if its best bit >= strong_frac * top_bit (default 0.5)")
    ap.add_argument("--output", required=True, help="Output hit_profile.json")
    args = ap.parse_args()

    hits = parse_hits([args.blast, args.mmseqs])
    profile = build_profile(hits, _read_query_len(args.query), args.locus_window, args.strong_frac)
    with open(args.output, "w") as fh:
        json.dump(profile, fh, indent=2)
    print(
        f"hit_profile: n_hits={profile['n_hits']} loci={profile['n_independent_loci']} "
        f"strong_loci={profile['n_strong_loci']} top_bit={profile['top_bit']} "
        f"bit_ratio={profile['bit_ratio']} query_len={profile['query_len']}"
    )


if __name__ == "__main__":
    main()
