#!/usr/bin/env python3
import argparse
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Split BED file into multiple files based on locus clustering.")
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output BED files")
    parser.add_argument("--dist_threshold", type=int, default=50000, help="Distance threshold for clustering (bp)")
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
    bits = [e['bit'] for e in cl if e.get('bit') is not None]
    return max(bits) if bits else None


def _cluster_span(cl):
    return cl[0]['chrom'], min(e['start'] for e in cl), max(e['end'] for e in cl)

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

    # Sort
    entries.sort(key=lambda x: (x['chrom'], x['start']))
    
    clusters = []
    current_cluster = []
    
    for entry in entries:
        if not current_cluster:
            current_cluster.append(entry)
            continue
            
        last = current_cluster[-1]
        
        # Check same chrom
        if entry['chrom'] != last['chrom']:
            clusters.append(current_cluster)
            current_cluster = [entry]
            continue
            
        # Check distance
        dist = entry['start'] - last['end']
        if dist > args.dist_threshold:
            clusters.append(current_cluster)
            current_cluster = [entry]
            continue
            
        current_cluster.append(entry)
        
    if current_cluster:
        clusters.append(current_cluster)
    
    # Rank clusters by significance (best e-value in cluster)
    for cl in clusters:
        cl_evalues = [e['evalue'] for e in cl]
        cl.sort(key=lambda x: x['evalue'])  # sort within cluster for reference
    
    clusters.sort(key=lambda cl: min(e['evalue'] for e in cl))
    
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

        # Hard cap: never keep more than max_loci, whatever the ratio test allowed.
        if len(kept) > args.max_loci:
            rank = (lambda cl: -(_cluster_best_bit(cl) or 0.0)) if has_bits \
                else (lambda cl: _cluster_best_evalue(cl))
            survivors = set(map(id, sorted(kept, key=rank)[:args.max_loci]))
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
