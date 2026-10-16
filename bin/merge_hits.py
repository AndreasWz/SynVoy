#!/usr/bin/env python3
"""The places where the query aligns in the home genome, one row per gene copy.

The home search (MMseqs2 and tblastn, 12-column tables) reports each exon of each copy of the
query's gene. Hits with E <= --max_evalue are kept, one per stretch of DNA
(sequence_utils.load_home_hits: both programs report every exon), and chained per scaffold and
strand in coding order. A hit extends a chain when

  * it continues the query after the chain's last hit (sequence_utils.continues_query: it
    starts at most HOME_CHAIN_QUERY_TOL = 20 residues before that hit ends, and ends after
    it), so a second copy of the gene starts a chain of its own however close it lies;
  * the genomic gap to that hit is at most --max_intron (an intron of the gene); and
  * the two do not lie in two different annotated genes, and no annotated gene lies entirely
    between them (with --gff; genes nested in an intron of the gene stay possible).

When several chains could take a hit, it goes to the one whose last query end is closest to
the hit's query start. Each chain is one row of the 7-column BED that split_loci.py reads:
scaffold, start, end, ``gene_loc``, the best E-value, the strand and the best bit score, over the
chain's hits and the duplicates they replaced (as before; split_loci.py scores a locus from
the hit tables themselves, as one copy of the gene, so it does not need the copy's score here).

--max_intron 'auto' (the default) is the largest intron of the annotated genes that the hits
overlap, times --intron_margin, at least --intron_floor and at most MAX_INTRON_CAP (1.1 Mb:
the longest introns of human genes are about 1 Mb); the floor without an annotation. The
gene model step's goi_info.json, from which the search derives its own intron bound, does not
exist yet when this runs.

Until 2026-10-15 (batch HG-1A) hits on one strand were joined when the next began less than
1 kb after the current end, whatever their query coordinates: two copies 500 bp apart became
one interval (fly Adh and Adhr), and the exons of one gene with introns over 1 kb stayed apart
(DPP4 in the honeybee: 15 intervals, firefly luciferase 99).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import (cds_models_by_gene, continues_query, in_coding_order,  # noqa: E402
                            load_home_hits)
from synvoy_scale import genes_between, load_gene_spans  # noqa: E402

MAX_INTRON_CAP = 1_100_000      # bp; the bound of 'auto' (Andreas, 2026-10-15)
# How far a hit may start in the query before the chain's last hit ends and still continue it
# (Andreas, 2026-10-16). The alignment of an exon often runs into the intron before it, whose
# translation matches the end of the previous exon's part of the query: on 140 links between
# hits inside one annotated gene (family seeds, firefly luciferase, fly Adh) the overlap was
# at most 19 residues (6 links over 10, among them the home genes of DPP4, APYR and
# luciferase), while a second copy of the gene overlapped by 43 or more. S0's shared
# QUERY_CONTINUATION_TOL (10) stays for its other callers (the flank walk, the GOI hit
# chain); whether they need 20 as well is an open question.
HOME_CHAIN_QUERY_TOL = 20


def gene_models_of(gff, chroms):
    """Coding models of the annotated genes on ``chroms`` (cds_models_by_gene), {} without a
    usable GFF."""
    if not gff or os.path.basename(str(gff)) == 'NO_GFF' or not os.path.exists(gff):
        return {}
    return cds_models_by_gene(gff, set(chroms))


def largest_intron(hits, models):
    """The largest intron (gap between two CDS rows of the longest coding model) of the
    annotated genes that a hit overlaps; 0 when none has one."""
    by_chrom = {}
    for h in hits:
        by_chrom.setdefault(h['chrom'], []).append((h['start'], h['end']))
    best = 0
    for info in models.values():
        ivs = by_chrom.get(info['chrom'])
        if not ivs:
            continue
        exons = sorted((s - 1, e) for s, e, _p in info['exons'])
        lo, hi = exons[0][0], max(e for _s, e in exons)
        if not any(s < hi and e > lo for s, e in ivs):
            continue
        for (_s1, e1), (s2, _e2) in zip(exons, exons[1:]):
            best = max(best, s2 - e1)
    return best


def resolve_max_intron(value, hits, models, margin, floor):
    """--max_intron: a number as given, or 'auto' (see the module docstring)."""
    text = str(value).strip().lower()
    if text not in ('auto', ''):
        return int(float(text)), 'set by the user'
    largest = largest_intron(hits, models)
    bound = int(min(MAX_INTRON_CAP, max(floor, largest * margin)))
    if not models:
        return bound, f'no annotation: the floor of {floor} bp'
    return bound, (f'{largest} bp, the largest intron of the annotated genes the hits overlap, '
                   f'x {margin:g}, within {floor}..{MAX_INTRON_CAP} bp')


def _genes_at(h, spans):
    return {(s, e) for s, e in spans if s < h['end'] and e > h['start']}


def may_join(a, b, spans):
    """Whether hits ``a`` and ``b`` (one scaffold and strand) can belong to one gene copy as
    far as the annotation tells: not two different annotated genes, and no annotated gene
    entirely between them. A gene both overlap (nested genes in its introns) allows it."""
    if not spans:
        return True
    ga, gb = _genes_at(a, spans), _genes_at(b, spans)
    if ga & gb:
        return True
    if ga and gb:
        return False
    left_end, right_start = min(a['end'], b['end']), max(a['start'], b['start'])
    return right_start <= left_end or genes_between(spans, left_end, right_start) == 0


def chain_hits(hits, max_intron, gene_spans):
    """Chains of hits, one per gene copy (see the module docstring), each a list of hits in
    coding order; the chains in the order of their first hit along the scaffold."""
    groups = {}
    for h in hits:
        groups.setdefault((h['chrom'], h['strand']), []).append(h)
    chains = []
    for (chrom, strand), group in sorted(groups.items()):
        spans = gene_spans.get(chrom, [])
        open_chains = []
        for h in in_coding_order(group, strand):
            best = None
            for k, chain in enumerate(open_chains):
                a = chain[-1]
                gap = h['start'] - a['end'] if strand == '+' else a['start'] - h['end']
                if (gap > max_intron or not continues_query(a, h, tol=HOME_CHAIN_QUERY_TOL)
                        or not may_join(a, h, spans)):
                    continue
                key = (abs(h['qstart'] - a['qend']), max(gap, 0), k)
                if best is None or key < best[0]:
                    best = (key, chain)
            if best is None:
                open_chains.append([h])
            else:
                best[1].append(h)
        chains.extend(open_chains)
    chains.sort(key=lambda c: (c[0]['chrom'], min(h['start'] for h in c), max(h['end'] for h in c)))
    return chains


def chain_row(chain):
    """The BED row of a chain: the stretch, best E-value and best bit score of its hits and of
    the duplicates they replaced (the other program's hits of the same exons), as the 1-kb
    merge reported them."""
    every = [x for h in chain for x in [h] + h.get('duplicates', [])]
    start = min(h['start'] for h in every)
    end = max(h['end'] for h in every)
    evalue = min(h['evalue'] for h in every)
    bits = max(h['bits'] for h in every)
    return (f"{chain[0]['chrom']}\t{start}\t{end}\tgene_loc\t{evalue}\t{chain[0]['strand']}\t"
            f"{bits:g}")


def main():
    parser = argparse.ArgumentParser(
        description="Chain the home hits of the query into one BED row per gene copy")
    parser.add_argument("--hits", action="append", default=[],
                        help="Raw hit table of the home search (BLAST outfmt 6 / MMseqs2 m8, 12 "
                             "columns); repeat for several")
    parser.add_argument("--output", required=True, help="Output BED file")
    parser.add_argument("--max_evalue", type=float, default=0.01,
                        help="Maximum E-value of a hit (default: 0.01; the pipeline passes "
                             "search_evalue)")
    parser.add_argument("--gff", default=None,
                        help="Home genome annotation: the intron bound of 'auto' and the rule "
                             "that a chain never joins two annotated genes; 'NO_GFF' = none")
    parser.add_argument("--max_intron", default="auto",
                        help="Largest genomic gap between two hits of one gene copy, in bp, or "
                             "'auto' (default; see the module docstring)")
    parser.add_argument("--intron_margin", type=float, default=5.0,
                        help="'auto': factor on the largest intron of the overlapped genes "
                             "(default 5.0; the pipeline passes goi_fallback_intron_margin)")
    parser.add_argument("--intron_floor", type=int, default=10000,
                        help="'auto': the smallest bound, also the bound without an annotation "
                             "(default 10000; the pipeline passes goi_fallback_intron_floor)")
    args = parser.parse_args()

    hits = load_home_hits(args.hits, args.max_evalue, with_duplicates=True)
    if not hits:
        print(f"No hit with E <= {args.max_evalue}.")
        open(args.output, 'w').close()
        return
    chroms = {h['chrom'] for h in hits}
    models = gene_models_of(args.gff, chroms)
    gene_spans = load_gene_spans(args.gff, chroms)
    max_intron, why = resolve_max_intron(args.max_intron, hits, models, args.intron_margin,
                                         args.intron_floor)
    chains = chain_hits(hits, max_intron, gene_spans)
    print(f"{len(hits)} hit(s) with E <= {args.max_evalue} -> {len(chains)} gene cop(y/ies); "
          f"largest gap within a copy {max_intron} bp ({why}).")
    with open(args.output, 'w') as out:
        for chain in chains:
            out.write(chain_row(chain) + "\n")


if __name__ == "__main__":
    main()
