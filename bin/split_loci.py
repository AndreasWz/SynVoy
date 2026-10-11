#!/usr/bin/env python3
"""Group the places where the query aligns in the home genome into loci.

Two neighbouring hit intervals belong to one locus when few annotated genes lie between
them (--locus_gap_genes, default 5): the flanking genes of two such places would largely
be the same genes, so searching them separately would search one neighbourhood twice.
Counting genes instead of base pairs makes the rule mean the same in every genome: 50 kb,
the fixed value used until 2026-10-08, are about 30 genes in yeast and one in human.
Without an annotation, or with --locus_gap set to a number, the gap is taken in base pairs.
"""
import argparse
import bisect
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import one_copy_score  # noqa: E402
from synvoy_scale import FALLBACK_LOCUS_GAP_BP, load_gene_spans  # noqa: E402,F401

def parse_args():
    parser = argparse.ArgumentParser(description="Split BED file into multiple files based on locus clustering.")
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output BED files")
    parser.add_argument("--gff", default=None,
                        help="Home genome annotation (GFF3). With it, loci are separated by "
                             "counting annotated genes; 'NO_GFF' or a missing file = none")
    parser.add_argument("--locus_gap_genes", type=int, default=5,
                        help="Two hit intervals with at most this many annotated genes between "
                             "them are one locus (default: 5). Needs --gff")
    parser.add_argument("--locus_gap", default="auto",
                        help="The gap in base pairs instead. 'auto' (default): count genes "
                             "when an annotation is given, else 50000 bp")
    parser.add_argument("--home_gene_bed", default=None,
                        help="goi_annotation.bed of the home gene model step: the coding exons "
                             "of the annotated home gene that IS the query")
    parser.add_argument("--home_gene_info", default=None,
                        help="goi_info.json of the same step; the gene is used only when its "
                             "method is 'gff_annotation' (an annotated gene matched the query)")
    parser.add_argument("--hits", action="append", default=[],
                        help="Raw hit table of the home search (BLAST outfmt 6 / MMseqs m8, 12 "
                             "columns); repeat for several. With them a locus is scored as "
                             "ONE copy of the gene: the best set of its hits that covers the "
                             "query at most once, i.e. the whole gene, not its best exon and "
                             "not the sum over several copies")
    parser.add_argument("--max_evalue", type=float, default=0.01,
                        help="E-value limit for the hits of --hits (default: 0.01)")
    parser.add_argument("--home_locus_mode", default="auto",
                        choices=["auto", "gene_only", "alignment"],
                        help="auto (default): the annotated home gene is the primary locus, "
                             "further loci come from the alignments. gene_only: only the "
                             "annotated gene's locus. alignment: loci from the alignments "
                             "alone, as before 2026-10-08")
    parser.add_argument("--max_loci", type=int, default=5,
                        help="Hard cap on the number of distinct loci to keep (default: 5). "
                             "Protects against large paralog families (e.g. luciferase / "
                             "AMP-binding) producing dozens of spurious loci.")
    parser.add_argument("--min_bit_ratio", type=float, default=0.5,
                        help="Keep a secondary locus only if its best bit score is at least this "
                             "fraction of the primary locus's best bit score (default: 0.5). "
                             "Ignored when the input BED carries no bit-score column.")
    parser.add_argument("--family_warning_threshold", type=int, default=10,
                        help="Emit a large-paralog-family WARNING when the pre-filter locus "
                             "count exceeds this (default: 10).")
    return parser.parse_args()


def _cluster_best_evalue(cl):
    return min(e['evalue'] for e in cl)


def _cluster_best_bit(cl):
    """Score of a locus. With raw hits (see score_clusters_by_hits): the score of one
    copy of the gene there, stored on the first entry. Otherwise the best bit score of
    its intervals, which is the best exon."""
    if cl and cl[0].get('locus_score') is not None:
        return cl[0]['locus_score']
    bits = [e['bit'] for e in cl if e.get('bit') is not None]
    return max(bits) if bits else None


def load_raw_hits(paths, max_evalue):
    """Hits of the home search as dicts (chrom, start, end, bits), E <= max_evalue, one
    per stretch of DNA: of hits overlapping by more than half of the shorter one the
    best-scoring stays (both search programs report every exon)."""
    hits = []
    for path in paths or []:
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) < 12:
                    continue
                try:
                    a, b = int(p[8]), int(p[9])
                    evalue, bits = float(p[10]), float(p[11])
                except ValueError:
                    continue
                if evalue > max_evalue:
                    continue
                try:
                    qa, qb = int(p[6]), int(p[7])
                except ValueError:
                    continue
                hits.append({'chrom': p[1], 'start': min(a, b) - 1, 'end': max(a, b), 'bits': bits,
                             'qstart': min(qa, qb), 'qend': max(qa, qb)})
    kept = []
    by_chrom = {}
    for h in sorted(hits, key=lambda x: (-x['bits'], x['chrom'], x['start'])):
        length = h['end'] - h['start']
        clash = False
        for k in by_chrom.get(h['chrom'], []):
            overlap = min(h['end'], k['end']) - max(h['start'], k['start'])
            if overlap > 0.5 * min(length, k['end'] - k['start']):
                clash = True
                break
        if not clash:
            by_chrom.setdefault(h['chrom'], []).append(h)
            kept.append(h)
    return kept


def score_clusters_by_hits(clusters, raw_hits):
    """Give every locus the score of one copy of the gene there (see one_copy_score).

    The interval table that defines the loci keeps one bit score per interval, the best
    one, so a locus used to be ranked by its best exon. A relative that shares one
    well-conserved domain then scores like the gene itself (human decorin: four loci of
    other leucine-rich-repeat genes passed the "half of the primary" test and each was
    searched in every genome, while the locus of its real relatives ranked sixth and was
    cut). Summed over the exons of one copy, the score measures the whole gene.
    """
    if not raw_hits:
        return
    by_chrom = {}
    for h in raw_hits:
        by_chrom.setdefault(h['chrom'], []).append(h)
    for cl in clusters:
        chrom, start, end = _cluster_span(cl)
        total = one_copy_score([h for h in by_chrom.get(chrom, [])
                                if h['start'] < end and h['end'] > start])
        cl[0]['locus_score'] = total if total > 0 else None


def _cluster_span(cl):
    return cl[0]['chrom'], min(e['start'] for e in cl), max(e['end'] for e in cl)


FALLBACK_GAP_BP = FALLBACK_LOCUS_GAP_BP


def genes_between(sorted_spans, left_end, right_start):
    """Annotated genes that lie entirely between two positions."""
    count = 0
    for i in range(bisect.bisect_left(sorted_spans, (left_end, -1)), len(sorted_spans)):
        start, end = sorted_spans[i]
        if start >= right_start:
            break
        if end <= right_start:
            count += 1
    return count


def load_home_gene(bed_path, info_path):
    """``(chrom, start, end, strand)`` of the annotated home gene that matched the query,
    or None: no file, no rows, or the home gene model was not taken from the annotation."""
    import json
    if not bed_path or not os.path.exists(bed_path):
        return None
    if info_path:
        try:
            with open(info_path) as fh:
                if json.load(fh).get('method') != 'gff_annotation':
                    return None
        except (OSError, ValueError):
            return None
    rows = []
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3 and parts[1].isdigit() and parts[2].isdigit():
                rows.append((parts[0], int(parts[1]), int(parts[2]), parts[5] if len(parts) > 5 else '+'))
    if not rows:
        return None
    chrom = rows[0][0]
    rows = [r for r in rows if r[0] == chrom]
    return chrom, min(r[1] for r in rows), max(r[2] for r in rows), rows[0][3]


def apply_home_gene(clusters, gene):
    """Make the annotated home gene the primary locus.

    Every cluster of hits that overlaps the gene is replaced by ONE entry, the gene's own
    span, which carries the best E-value and bit score of the hits it replaces. The gene
    is what the query is; the hits around it are that gene's exons, and often also hits
    in a relative next door, which made the locus several intervals with a walk for
    flanking genes from each. Returns ``(primary_cluster, other_clusters)``.
    """
    chrom, start, end, strand = gene
    inside, others = [], []
    for cl in clusters:
        if any(e['chrom'] == chrom and e['start'] < end and e['end'] > start for e in cl):
            inside.extend(cl)
        else:
            others.append(cl)
    evalue = min((e['evalue'] for e in inside), default=0.0)
    bits = [e['bit'] for e in inside if e.get('bit') is not None]
    bit = max(bits) if bits else None
    bit_text = f"{bit:g}" if bit is not None else "."
    primary = [{
        'chrom': chrom, 'start': start, 'end': end, 'evalue': evalue, 'bit': bit,
        'line': f"{chrom}\t{start}\t{end}\tgene_loc\t{evalue}\t{strand}\t{bit_text}",
    }]
    return primary, others


def split_into_loci(entries, gene_spans, gap_genes, gap_bp):
    """Cluster hit intervals (dicts with chrom, start, end), sorted by position.

    ``gap_bp`` (a number) separates loci by distance. With ``gap_bp`` None, loci on a
    sequence with annotated genes are separated when more than ``gap_genes`` genes lie
    between them, and by FALLBACK_GAP_BP elsewhere."""
    clusters, current, current_end = [], [], 0
    for entry in sorted(entries, key=lambda x: (x['chrom'], x['start'])):
        if current and entry['chrom'] == current[-1]['chrom']:
            spans = gene_spans.get(entry['chrom'])
            if gap_bp is not None:
                same = entry['start'] - current_end <= gap_bp
            elif spans:
                same = genes_between(spans, current_end, entry['start']) <= gap_genes
            else:
                same = entry['start'] - current_end <= FALLBACK_GAP_BP
            if same:
                current.append(entry)
                current_end = max(current_end, entry['end'])
                continue
        if current:
            clusters.append(current)
        current, current_end = [entry], entry['end']
    if current:
        clusters.append(current)
    return clusters

def main():
    args = parse_args()
    
    entries = []
    try:
        with open(args.bed) as f:
            for line in f:
                if not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                
                evalue = float(parts[4]) if len(parts) > 4 and parts[4] not in ('.', '') else float('inf')
                bit = float(parts[6]) if len(parts) > 6 and parts[6] not in ('.', '') else None
                entries.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'evalue': evalue,
                    'bit': bit,
                    'line': line.strip()
                })
    except FileNotFoundError:
        print(f"ERROR: BED file {args.bed} not found. Expected output from LOCATE_GENE.", file=sys.stderr)
        sys.exit(1)

    if not entries:
        print(
            "ERROR: Gene not found in home genome. Both BLAST and MMseqs2 returned no hits "
            f"above the e-value threshold.\n"
            "Why: The query sequence may be too divergent from the home genome, the gene may be "
            "absent, or the e-value cutoff may be too strict.\n"
            "Try: Lower --search_evalue (e.g. 1e-2), verify the query FASTA is the correct "
            "protein/CDS sequence, or check that --home_genome points to the right assembly.",
            file=sys.stderr
        )
        sys.exit(1)

    gap_text = str(args.locus_gap).strip().lower()
    gap_bp = None if gap_text in ('auto', '') else int(float(gap_text))
    gene_spans = {} if gap_bp is not None else load_gene_spans(
        args.gff, {e['chrom'] for e in entries})
    if gap_bp is not None:
        print(f"Locus gap: {gap_bp} bp (set by the user).")
    elif gene_spans:
        print(f"Locus gap: more than {args.locus_gap_genes} annotated gene(s) between two hits "
              f"start a new locus.")
    else:
        print(f"Locus gap: {FALLBACK_GAP_BP} bp (no annotation to count genes in).")
    clusters = split_into_loci(entries, gene_spans, args.locus_gap_genes, gap_bp)

    # Rank clusters by significance (best e-value in cluster)
    for cl in clusters:
        cl_evalues = [e['evalue'] for e in cl]
        cl.sort(key=lambda x: x['evalue'])  # sort within cluster for reference
    
    clusters.sort(key=lambda cl: min(e['evalue'] for e in cl))

    # The annotated home gene, when the query is one, is the primary locus whatever the
    # ranking of the hits says (a repeat domain can make a relative's locus score higher).
    raw_hits = load_raw_hits(args.hits, args.max_evalue)
    home_gene = None if args.home_locus_mode == 'alignment' else load_home_gene(
        args.home_gene_bed, args.home_gene_info)
    if home_gene:
        primary, others = apply_home_gene(clusters, home_gene)
        print(f"Primary locus = the annotated home gene {home_gene[0]}:{home_gene[1]}-{home_gene[2]} "
              f"({home_gene[3]}); {len(others)} further place(s) where the query aligns.")
        clusters = [primary] + ([] if args.home_locus_mode == 'gene_only' else others)
    elif args.home_locus_mode == 'gene_only':
        print("WARNING: --home_locus_mode gene_only, but no annotated home gene matched the "
              "query; the loci are taken from the alignments.", file=sys.stderr)
    score_clusters_by_hits(clusters, raw_hits)
    if raw_hits and not home_gene:
        # Without an annotated gene the primary locus is the one that explains most of
        # the query.
        clusters.sort(key=lambda cl: -(_cluster_best_bit(cl) or 0.0))
    
    # Filter: keep the primary locus + only secondary loci of comparable strength.
    #
    # The old filter OR'd "has >=2 HSPs" with an e-value-ratio test. Both fail for
    # large paralog families: most paralogs have multiple HSPs, and when the primary
    # hit's e-value is BLAST's saturated 0.0 the ratio test degenerates to "keep
    # anything <= 1e-10". Firefly luciferase (an AMP-binding family member) leaked
    # all 49 loci through this way and the run took ~55 h (docs/TODO.md §1d).
    #
    # Bit scores do not saturate, so rank by best_bit / primary_bit instead. The
    # max_loci cap is a hard backstop regardless of the ratio test, and applies even
    # when the input carries no bit column (older 6-column BEDs).
    n_prefilter = len(clusters)
    if n_prefilter > args.family_warning_threshold:
        print(
            f"WARNING: {n_prefilter} candidate loci before filtering. The query likely belongs "
            f"to a large paralog family (e.g. AMP-binding/acyl-CoA synthetase, LY6/3FTx, MRJP). "
            f"Only the strongest {args.max_loci} will be kept. To adjust: raise/lower --max_loci, "
            f"tighten --search_evalue, or use a more specific query.",
            file=sys.stderr,
        )

    if len(clusters) > 1:
        primary = clusters[0]                       # best e-value cluster (clusters already sorted)
        primary_bit = _cluster_best_bit(primary)
        has_bits = all(_cluster_best_bit(cl) is not None for cl in clusters)

        kept = [primary]
        dropped = []  # list of (cluster, reason)
        for cl in clusters[1:]:
            if has_bits and primary_bit and primary_bit > 0:
                ratio = _cluster_best_bit(cl) / primary_bit
                if ratio >= args.min_bit_ratio:
                    kept.append(cl)
                else:
                    dropped.append((cl, f"bit {_cluster_best_bit(cl):.0f} = {ratio:.2f}x primary "
                                        f"{primary_bit:.0f}, below --min_bit_ratio {args.min_bit_ratio}"))
            else:
                # No usable bit scores: defer entirely to the max_loci cap below.
                kept.append(cl)

        # Hard cap: never keep more than max_loci, whatever the ratio test allowed. The
        # primary locus always stays: it is the query's own gene (or its best hit), and a
        # relative's locus with one better-scoring exon must not push it out.
        if len(kept) > args.max_loci:
            rank = (lambda cl: -(_cluster_best_bit(cl) or 0.0)) if has_bits \
                else (lambda cl: _cluster_best_evalue(cl))
            survivors = {id(primary)} | set(
                map(id, sorted(kept[1:], key=rank)[:max(0, args.max_loci - 1)]))
            dropped.extend((cl, f"exceeds --max_loci {args.max_loci}")
                           for cl in kept if id(cl) not in survivors)
            kept = [cl for cl in kept if id(cl) in survivors]

        for cl, reason in dropped:
            chrom, s, e = _cluster_span(cl)
            print(f"Filtered out locus {chrom}:{s}-{e} "
                  f"(best e-value={_cluster_best_evalue(cl):.2e}; {reason})")
        clusters = kept

    print(f"Found {len(clusters)} distinct loci.")
    
    for i, cl in enumerate(clusters):
        # We name the file so Nextflow can pick it up.
        # Use simple indexing.
        out_name = f"{args.output_prefix}_{i+1}.bed"
        with open(out_name, 'w') as f_out:
            for entry in cl:
                f_out.write(entry['line'] + "\n")
        best_ev = min(e['evalue'] for e in cl)
        print(f"Wrote locus {i+1} to {out_name} ({len(cl)} hits, best e-value={best_ev:.2e})")

if __name__ == "__main__":
    main()
