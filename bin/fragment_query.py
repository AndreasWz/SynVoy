#!/usr/bin/env python3
"""
fragment_query.py - Generate sequence fragments for improved gene discovery

This script generates overlapping fragments of a query protein sequence
at different granularities (halves, thirds, quarters) to improve detection
of divergent or partial gene matches in target genomes.

Usage:
    python fragment_query.py --query input.faa --output fragments.faa [--min_size 20]

Output format:
    Each fragment has an ID like: {original_id}|frag_{start}_{end}_{total}
    Description includes: fragment_type=half|third|quarter position=1|2|3|4
"""

import argparse
import os
import sys

# Use our own sequence utilities (no BioPython)
try:
    from sequence_utils import parse_fasta, write_fasta
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta


def generate_fragments(seq, seq_id, min_size=20):
    """
    Generate sequence fragments at different granularities.
    
    Args:
        seq: Sequence string
        seq_id: Original sequence ID
        min_size: Minimum fragment size in amino acids
    
    Returns:
        List of (id, seq, description) tuples for each fragment
    """
    fragments = []
    length = len(seq)
    
    # 1. Full sequence (always include)
    fragments.append((
        f"{seq_id}|frag_1_{length}_{length}",
        seq,
        f"fragment_type=full position=1/1 length={length}"
    ))
    
    # 2. Halves
    if length // 2 >= min_size:
        half = length // 2
        
        # First half
        frag1 = seq[:half]
        fragments.append((
            f"{seq_id}|frag_1_{half}_{length}",
            frag1,
            f"fragment_type=half position=1/2 length={len(frag1)}"
        ))
        
        # Second half
        frag2 = seq[half:]
        fragments.append((
            f"{seq_id}|frag_{half+1}_{length}_{length}",
            frag2,
            f"fragment_type=half position=2/2 length={len(frag2)}"
        ))
    
    # 3. Thirds
    if length // 3 >= min_size:
        third = length // 3
        
        # First third
        frag1 = seq[:third]
        fragments.append((
            f"{seq_id}|frag_1_{third}_{length}",
            frag1,
            f"fragment_type=third position=1/3 length={len(frag1)}"
        ))
        
        # Second third
        frag2 = seq[third:2*third]
        fragments.append((
            f"{seq_id}|frag_{third+1}_{2*third}_{length}",
            frag2,
            f"fragment_type=third position=2/3 length={len(frag2)}"
        ))
        
        # Third third (includes remainder)
        frag3 = seq[2*third:]
        fragments.append((
            f"{seq_id}|frag_{2*third+1}_{length}_{length}",
            frag3,
            f"fragment_type=third position=3/3 length={len(frag3)}"
        ))
    
    # 4. Quarters
    if length // 4 >= min_size:
        quarter = length // 4
        
        for i in range(4):
            start = i * quarter
            end = (i + 1) * quarter if i < 3 else length  # Last quarter gets remainder
            
            frag = seq[start:end]
            fragments.append((
                f"{seq_id}|frag_{start+1}_{end}_{length}",
                frag,
                f"fragment_type=quarter position={i+1}/4 length={len(frag)}"
            ))
    
    # 5. Sliding window (optional, for very long sequences)
    # This creates overlapping windows to catch edge cases
    SLIDING_WINDOW_SIZE = 50
    SLIDING_STEP = 25
    
    if length > SLIDING_WINDOW_SIZE * 2:  # Only for longer sequences
        pos = 0
        window_num = 1
        while pos + SLIDING_WINDOW_SIZE <= length:
            frag = seq[pos:pos + SLIDING_WINDOW_SIZE]
            fragments.append((
                f"{seq_id}|slide_{pos+1}_{pos+SLIDING_WINDOW_SIZE}_{length}",
                frag,
                f"fragment_type=sliding position=window_{window_num} length={len(frag)}"
            ))
            pos += SLIDING_STEP
            window_num += 1
    
    return fragments


def parse_fragment_id(frag_id):
    """
    Parse a fragment ID to extract original gene ID and fragment info.
    
    Args:
        frag_id: Fragment ID like "gene-LOC726866|frag_1_50_100"
    
    Returns:
        dict with keys: gene_id, start, end, total, fragment_type
    """
    parts = frag_id.split('|')
    gene_id = parts[0]
    
    result = {
        'gene_id': gene_id,
        'start': 1,
        'end': None,
        'total': None,
        'fragment_type': 'full'
    }
    
    if len(parts) > 1:
        frag_part = parts[-1]  # Last part contains fragment info
        
        if frag_part.startswith('frag_'):
            # Format: frag_START_END_TOTAL
            try:
                coords = frag_part.replace('frag_', '').split('_')
                result['start'] = int(coords[0])
                result['end'] = int(coords[1])
                result['total'] = int(coords[2])
                
                # Determine type based on positions
                if result['end'] == result['total'] and result['start'] == 1:
                    result['fragment_type'] = 'full'
                elif result['total'] and result['end'] - result['start'] + 1 > result['total'] // 2:
                    result['fragment_type'] = 'half'
                elif result['total'] and result['end'] - result['start'] + 1 > result['total'] // 4:
                    result['fragment_type'] = 'third'
                else:
                    result['fragment_type'] = 'quarter'
                    
            except (ValueError, IndexError):
                pass
                
        elif frag_part.startswith('slide_'):
            result['fragment_type'] = 'sliding'
            try:
                coords = frag_part.replace('slide_', '').split('_')
                result['start'] = int(coords[0])
                result['end'] = int(coords[1])
                result['total'] = int(coords[2])
            except (ValueError, IndexError):
                pass
                
        elif frag_part.startswith('exon_'):
            result['fragment_type'] = 'exon'
    
    return result


def merge_fragment_hits(hits, max_gap=50):
    """
    Merge hits from different fragments of the same gene.
    
    Args:
        hits: List of hit dictionaries with fragment IDs
        max_gap: Maximum gap between fragments to merge (in target coordinates)
    
    Returns:
        List of merged hit dictionaries
    """
    from collections import defaultdict
    
    # Group hits by base gene ID
    hits_by_gene = defaultdict(list)
    for h in hits:
        parsed = parse_fragment_id(h['query'])
        gene_id = parsed['gene_id']
        h['_parsed_fragment'] = parsed
        hits_by_gene[gene_id].append(h)
    
    merged_hits = []
    
    for gene_id, gene_hits in hits_by_gene.items():
        # Sort by target position
        gene_hits.sort(key=lambda x: (x['chrom'], x['start']))
        
        # Merge overlapping/adjacent hits
        current_group = [gene_hits[0]]
        
        for hit in gene_hits[1:]:
            last = current_group[-1]
            
            # Check if same chromosome and close enough
            if hit['chrom'] == last['chrom'] and hit['start'] - last['end'] <= max_gap:
                current_group.append(hit)
            else:
                # Output merged group
                merged_hits.append(_merge_hit_group(gene_id, current_group))
                current_group = [hit]
        
        # Don't forget last group
        if current_group:
            merged_hits.append(_merge_hit_group(gene_id, current_group))
    
    return merged_hits


def _merge_hit_group(gene_id, hits):
    """Merge a group of hits into a single hit."""
    if len(hits) == 1:
        return hits[0]
    
    # Combine metadata
    merged = {
        'query': gene_id,  # Use base gene ID
        'chrom': hits[0]['chrom'],
        'start': min(h['start'] for h in hits),
        'end': max(h['end'] for h in hits),
        'evalue': min(h.get('evalue', 1) for h in hits),
        'pident': max(h.get('pident', 0) for h in hits),
        'alnlen': sum(h.get('alnlen', 0) for h in hits),
        'fragments_merged': len(hits),
        'fragment_types': list(set(h['_parsed_fragment']['fragment_type'] for h in hits))
    }
    
    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Generate sequence fragments for improved gene discovery"
    )
    parser.add_argument("--query", required=True, 
                        help="Input FASTA file with query sequence(s)")
    parser.add_argument("--output", required=True,
                        help="Output FASTA file with fragments")
    parser.add_argument("--min_size", type=int, default=20,
                        help="Minimum fragment size in amino acids (default: 20)")
    parser.add_argument("--no_sliding", action="store_true",
                        help="Disable sliding window fragments")
    parser.add_argument("--summary", help="Output summary TSV file")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.query):
        print(f"Error: Query file not found: {args.query}", file=sys.stderr)
        sys.exit(1)
    
    all_fragments = []
    summary_data = []
    
    # Process each sequence in input
    for _, rec_id, rec_seq in parse_fasta(args.query):
        print(f"Processing {rec_id} (length: {len(rec_seq)})")
        
        fragments = generate_fragments(rec_seq, rec_id, args.min_size)
        
        # Optionally filter out sliding windows
        if args.no_sliding:
            fragments = [f for f in fragments if 'slide_' not in f[0]]
        
        all_fragments.extend(fragments)
        
        # Summary stats
        frag_types = {}
        for f in fragments:
            # f is (id, seq, description)
            ftype = f[2].split('fragment_type=')[1].split()[0]
            frag_types[ftype] = frag_types.get(ftype, 0) + 1
        
        summary_data.append({
            'gene': rec_id,
            'length': len(rec_seq),
            'total_fragments': len(fragments),
            **frag_types
        })
    
    # Write fragments - convert to (id, seq) tuples for write_fasta
    output_records = [(f[0], f[1]) for f in all_fragments]
    write_fasta(output_records, args.output)
    print(f"Wrote {len(all_fragments)} fragments to {args.output}")
    
    # Write summary if requested
    if args.summary:
        with open(args.summary, 'w') as f:
            f.write("gene\tlength\ttotal_fragments\tfull\thalf\tthird\tquarter\tsliding\n")
            for row in summary_data:
                f.write(f"{row['gene']}\t{row['length']}\t{row['total_fragments']}\t")
                f.write(f"{row.get('full', 0)}\t{row.get('half', 0)}\t")
                f.write(f"{row.get('third', 0)}\t{row.get('quarter', 0)}\t")
                f.write(f"{row.get('sliding', 0)}\n")
        print(f"Summary written to {args.summary}")


if __name__ == "__main__":
    main()
