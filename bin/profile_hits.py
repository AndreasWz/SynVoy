#!/usr/bin/env python3
"""Summarise the LOCATE_GENE hit distribution into a compact hit_profile.json.

The shape of the home-genome BLAST/MMseqs hit distribution tells you what kind of
query this is — single-copy, a tandem family, or a sprawling paralog superfamily —
before any expensive iterative search runs. auto_select_preset.py consumes this to
pick a query-type preset, so a user who doesn't know which `-profile preset_*` to
choose still gets sensible parameters (docs/TODO.md §1f).

Metrics emitted:
  query_len            query length in aa (from --query FASTA, else max qend across hits)
  n_hits               alignment rows with E <= --max_evalue (BLAST + MMseqs combined);
                       n_hits_raw counts every row. A row is one exon in one search
                       program, so neither says how many genes there are.
  n_independent_loci   places where the query aligns: hits with at most --locus_gap_genes
                       annotated genes between them are one place (the definition of a
                       home locus in split_loci.py; --locus_gap bp without an annotation)
  n_strong_loci        independent loci whose score >= --strong_frac * top score
  top_bit, second_bit  score of the strongest / next-strongest independent locus
  bit_ratio            top_bit / second_bit  (a primary that barely beats the next
                       locus -> paralog-rich; one that dominates -> single-copy)
  primary_copies       how many times the query is covered at the strongest locus:
                       aligned query residues of its hits / query length. About 1 for
                       one gene, whatever its number of exons; about 2 for a tandem pair.
  identity_{min,q1,median,q3,max}, identity_iqr   percent-identity spread of the HSPs

The score of a locus is the score of one copy of the gene there: the best total bit
score of hits that cover the query at most once (sequence_utils.one_copy_score), after
hits that cover the same DNA are reduced to the best one (the two search programs report
every exon twice). Until 2026-10-08 it was the best single hit, i.e. the best exon, and
the decision downstream counted rows: an 11-exon single-copy gene had 26 rows and was
called a "large paralog family". tblastn rows up to E = 10 were counted as well.

Reads standard `-outfmt 6` / MMseqs `--format-output query,target,pident,alnlen,
mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits` (12 cols).
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import one_copy_score  # noqa: E402
from synvoy_scale import FALLBACK_LOCUS_GAP_BP, genes_between, load_gene_spans  # noqa: E402


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
                    qstart, qend = int(p[6]), int(p[7])
                    tstart, tend = int(p[8]), int(p[9])
                    evalue = float(p[10])
                    bits = float(p[11])
                except ValueError:
                    continue
                if tend < tstart:
                    tstart, tend = tend, tstart
                if qend < qstart:
                    qstart, qend = qend, qstart
                hits.append({
                    "chrom": p[1], "qstart": qstart, "qend": qend,
                    "tstart": tstart, "tend": tend,
                    "pident": pident, "bits": bits, "evalue": evalue,
                })
    return hits


def nonredundant(hits):
    """One hit per stretch of DNA: of hits that overlap by more than half of the shorter
    one, the best-scoring stays. Removes the second search program's copy of every exon."""
    kept = []
    for h in sorted(hits, key=lambda x: (-x["bits"], x["chrom"], x["tstart"])):
        length = h["tend"] - h["tstart"]
        clash = False
        for k in kept:
            if k["chrom"] != h["chrom"]:
                continue
            overlap = min(h["tend"], k["tend"]) - max(h["tstart"], k["tstart"])
            if overlap > 0.5 * min(length, k["tend"] - k["tstart"]):
                clash = True
                break
        if not clash:
            kept.append(h)
    return kept


def _cluster_loci(hits, window, gene_spans=None, gap_genes=5):
    """Merge HSPs into independent loci of the home genome.

    Two hits on one sequence belong to one locus when at most ``gap_genes`` annotated
    genes lie between them (``gene_spans``: {sequence: sorted (start, end)}), and on a
    sequence without annotated genes when they are at most ``window`` bp apart. With
    ``window`` given and ``gene_spans`` None the distance decides everywhere.
    Returns the loci, strongest first: ``{'bit': score of one copy of the gene there,
    'aligned': aligned query residues of all hits}`` each. Expects non-redundant hits."""
    def same_locus(locus, h):
        if locus["chrom"] != h["chrom"]:
            return False
        spans = (gene_spans or {}).get(h["chrom"])
        if spans:
            return genes_between(spans, locus["end"], h["tstart"]) <= gap_genes
        return h["tstart"] <= locus["end"] + window

    loci = []
    for h in sorted(hits, key=lambda x: (x["chrom"], x["tstart"])):
        aligned = h["qend"] - h["qstart"] + 1
        if loci and same_locus(loci[-1], h):
            loci[-1]["end"] = max(loci[-1]["end"], h["tend"])
            loci[-1]["hits"].append(h)
            loci[-1]["aligned"] += aligned
        else:
            loci.append({"chrom": h["chrom"], "end": h["tend"], "hits": [h], "aligned": aligned})
    for locus in loci:
        locus["bit"] = one_copy_score(locus.pop("hits"))
    return sorted(loci, key=lambda l: -l["bit"])


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


def build_profile(hits, query_len, window, strong_frac, max_evalue=None, gene_spans=None, gap_genes=5):
    n_raw = len(hits)
    if max_evalue is not None:
        hits = [h for h in hits if h.get("evalue", 0.0) <= max_evalue]
    if not hits:
        return {
            "query_len": query_len, "n_hits": 0, "n_hits_raw": n_raw, "n_independent_loci": 0,
            "n_strong_loci": 0, "top_bit": None, "second_bit": None,
            "bit_ratio": None, "primary_copies": None, "identity_min": None,
            "identity_q1": None, "identity_median": None, "identity_q3": None,
            "identity_max": None, "identity_iqr": None,
        }

    if query_len is None:
        query_len = max(h["qend"] for h in hits)

    loci = _cluster_loci(nonredundant(hits), window, gene_spans, gap_genes)
    top_bit = loci[0]["bit"]
    second_bit = loci[1]["bit"] if len(loci) > 1 else None
    bit_ratio = (top_bit / second_bit) if second_bit else None
    n_strong = sum(1 for l in loci if l["bit"] >= strong_frac * top_bit)
    primary_copies = loci[0]["aligned"] / query_len if query_len else None

    idents = sorted(h["pident"] for h in hits)
    q1 = _quantile(idents, 0.25)
    q3 = _quantile(idents, 0.75)

    return {
        "query_len": query_len,
        "n_hits": len(hits),
        "n_hits_raw": n_raw,
        "n_independent_loci": len(loci),
        "n_strong_loci": n_strong,
        "top_bit": round(top_bit, 1),
        "second_bit": round(second_bit, 1) if second_bit is not None else None,
        "bit_ratio": round(bit_ratio, 3) if bit_ratio is not None else None,
        "primary_copies": round(primary_copies, 2) if primary_copies is not None else None,
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
    ap.add_argument("--gff", default=None,
                    help="Home genome annotation (GFF3). With it, loci are separated by counting "
                         "annotated genes; 'NO_GFF' or a missing file = none")
    ap.add_argument("--locus_gap_genes", type=int, default=5,
                    help="Hits with at most this many annotated genes between them are one "
                         "locus (default 5; the rule split_loci.py uses for the home loci)")
    ap.add_argument("--locus_gap", default="auto",
                    help="The gap in base pairs instead. 'auto' (default): count genes when an "
                         f"annotation is given, else {FALLBACK_LOCUS_GAP_BP} bp")
    ap.add_argument("--strong_frac", type=float, default=0.5,
                    help="A locus counts as 'strong' if its best bit >= strong_frac * top_bit (default 0.5)")
    ap.add_argument("--max_evalue", type=float, default=0.01,
                    help="Only hits at or below this E-value are profiled (default 0.01, the "
                         "limit the home step applies to the hits it uses)")
    ap.add_argument("--output", required=True, help="Output hit_profile.json")
    args = ap.parse_args()

    hits = parse_hits([args.blast, args.mmseqs])
    gap_text = str(args.locus_gap).strip().lower()
    gap_bp = None if gap_text in ("auto", "") else int(float(gap_text))
    gene_spans = None if gap_bp is not None else load_gene_spans(args.gff, {h["chrom"] for h in hits})
    profile = build_profile(hits, _read_query_len(args.query),
                            gap_bp if gap_bp is not None else FALLBACK_LOCUS_GAP_BP, args.strong_frac,
                            max_evalue=args.max_evalue, gene_spans=gene_spans,
                            gap_genes=args.locus_gap_genes)
    profile["locus_unit"] = ("bp" if gap_bp is not None else
                             "genes" if gene_spans else "bp (no annotation)")
    with open(args.output, "w") as fh:
        json.dump(profile, fh, indent=2)
    print(
        f"hit_profile: n_hits={profile['n_hits']} (raw {profile['n_hits_raw']}) "
        f"loci={profile['n_independent_loci']} strong_loci={profile['n_strong_loci']} "
        f"primary_copies={profile['primary_copies']} top_bit={profile['top_bit']} "
        f"bit_ratio={profile['bit_ratio']} query_len={profile['query_len']}"
    )


if __name__ == "__main__":
    main()
