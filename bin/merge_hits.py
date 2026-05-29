#!/usr/bin/env python3

import sys
import argparse


def _parse_bit(parts):
    """Bit score lives in the optional 7th BED column (emitted by locate_gene.nf).
    Returns None when absent so older 6-column inputs still work."""
    if len(parts) > 6 and parts[6] not in ('.', ''):
        try:
            return float(parts[6])
        except ValueError:
            return None
    return None


def _max_bit(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return max(a, b)


def main():
    parser = argparse.ArgumentParser(description="Merge BLAST/MMseqs hits to BED")
    parser.add_argument("--mmseqs", help="MMseqs2 m8 output (converted to bed-like)")
    parser.add_argument("--blast", help="BLAST output (converted to bed-like)")
    parser.add_argument("--output", required=True, help="Output BED file")
    parser.add_argument("--max_evalue", type=float, default=1e-3,
                        help="Maximum e-value to keep a hit (default: 1e-3)")
    
    args = parser.parse_args()
    
    hits = []
    
    # Read MMseqs (if provided/exists)
    if args.mmseqs:
        try:
            with open(args.mmseqs) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        hits.append({
                            'chrom': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'source': 'mmseqs',
                            'score': float(parts[4]), # evalue/score
                            'strand': parts[5],
                            'bit': _parse_bit(parts)
                        })
        except FileNotFoundError:
            pass

    # Read BLAST (if provided/exists)
    if args.blast:
        try:
            with open(args.blast) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        hits.append({
                            'chrom': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'source': 'blast',
                            'score': float(parts[4]),
                            'strand': parts[5],
                            'bit': _parse_bit(parts)
                        })
        except FileNotFoundError:
            pass

    if not hits:
        open(args.output, 'w').close()
        return

    # Filter by e-value threshold
    before = len(hits)
    hits = [h for h in hits if h['score'] <= args.max_evalue]
    if before != len(hits):
        print(f"E-value filter ({args.max_evalue}): kept {len(hits)}/{before} hits")

    if not hits:
        print("No hits passed e-value filter.")
        open(args.output, 'w').close()
        return

    # Sort by chrom, start
    hits.sort(key=lambda x: (x['chrom'], x['start']))

    merged = []
    if hits:
        curr = hits[0]
        for next_hit in hits[1:]:
            # Check overlap
            # Same chrom and start within range of end (allow gaps)
            if (next_hit['chrom'] == curr['chrom'] and
                    next_hit['start'] < curr['end'] + 1000 and
                    next_hit['strand'] == curr['strand']):
                # Merge
                curr['end'] = max(curr['end'], next_hit['end'])
                # Keep best e-value (min) and best bit score (max) across merged HSPs.
                curr['score'] = min(curr['score'], next_hit['score'])
                curr['bit'] = _max_bit(curr.get('bit'), next_hit.get('bit'))
                # Strand voting? 
                if curr['strand'] != next_hit['strand']:
                    pass # Ambiguous, keep curr
            else:
                merged.append(curr)
                curr = next_hit
        merged.append(curr)

    with open(args.output, 'w') as out:
        for m in merged:
            # BED: chrom, start, end, name, evalue, strand, bitscore
            # The 7th (bit) column lets split_loci.py rank loci robustly where
            # BLAST e-values saturate at 0.0 (docs/TODO.md §1d). '.' when unknown.
            bit = m.get('bit')
            bit_str = f"{bit:g}" if bit is not None else "."
            out.write(f"{m['chrom']}\t{m['start']}\t{m['end']}\tgene_loc\t{m['score']}\t{m['strand']}\t{bit_str}\n")

if __name__ == "__main__":
    main()
