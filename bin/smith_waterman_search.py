#!/usr/bin/env python3
"""
smith_waterman_search.py - Smith-Waterman local alignment for GOI search

Uses parasail library for vectorized Smith-Waterman alignment.
This provides more sensitive local alignment than MMseqs2 for divergent sequences.

Usage:
    python smith_waterman_search.py --query goi.faa --target region.fna \\
        --output hits.tsv --matrix BLOSUM62 --gap_open 10 --gap_extend 1
"""

import argparse
import sys
import os

try:
    import parasail
    PARASAIL_AVAILABLE = True
except ImportError:
    PARASAIL_AVAILABLE = False
    print("WARNING: parasail not installed. Falling back to ssearch36.", file=sys.stderr)

# Use our own sequence utilities
try:
    from sequence_utils import parse_fasta, translate, reverse_complement
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, translate, reverse_complement


def translate_in_six_frames(dna_seq):
    """
    Translate DNA sequence in all 6 reading frames.
    Returns list of (frame_id, protein_seq, start_offset) tuples.
    """
    frames = []
    
    # Forward strand (frames +1, +2, +3)
    for frame in range(3):
        prot = translate(dna_seq[frame:])
        frames.append((f"+{frame+1}", prot, frame))
    
    # Reverse strand (frames -1, -2, -3)
    rev_seq = reverse_complement(dna_seq)
    for frame in range(3):
        prot = translate(rev_seq[frame:])
        frames.append((f"-{frame+1}", prot, frame))
    
    return frames


def smith_waterman_parasail(query_seq, target_seq, matrix_name='BLOSUM62', 
                            gap_open=10, gap_extend=1):
    """
    Perform Smith-Waterman alignment using parasail (fast vectorized implementation).
    
    Returns dict with alignment details.
    """
    if not PARASAIL_AVAILABLE:
        raise ImportError("parasail library required for Smith-Waterman search")
    
    # Get substitution matrix
    matrix = parasail.blosum62 if matrix_name == 'BLOSUM62' else parasail.pam100
    
    # Run Smith-Waterman with traceback
    result = parasail.sw_trace_striped_32(
        query_seq, target_seq, gap_open, gap_extend, matrix
    )
    
    # Calculate percent identity if traceback available
    if hasattr(result, 'traceback'):
        traceback = result.traceback
        identity = (traceback.comp.count('|') / len(traceback.comp)) * 100 if traceback.comp else 0
    else:
        # Estimate from score
        identity = (result.score / (len(query_seq) * 5)) * 100  # Rough estimate
    
    return {
        'score': result.score,
        'end_query': result.end_query,
        'end_ref': result.end_ref,
        'identity': identity,
        'length': result.end_query - result.end_ref if hasattr(result, 'end_ref') else len(query_seq)
    }


def smith_waterman_ssearch(query_faa, target_fna, output_tsv, threads=1):
    """
    Fallback: Use FASTA package's ssearch36 for Smith-Waterman.
    
    ssearch36 performs rigorous Smith-Waterman with no heuristics.
    """
    import subprocess
    
    # Check if ssearch36 is available
    try:
        subprocess.run(['ssearch36', '-h'], stdout=subprocess.DEVNULL, 
                      stderr=subprocess.DEVNULL, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: ssearch36 not found. Install FASTA package:", file=sys.stderr)
        print("  conda install -c bioconda fasta3", file=sys.stderr)
        sys.exit(1)
    
    # Run ssearch36: protein query vs translated DNA target
    cmd = [
        'ssearch36',
        '-m', '8',  # Tabular output (BLAST-like)
        '-T', str(threads),
        '-E', '10',  # Relaxed E-value
        '-3',  # Query is protein, target is DNA (translate target)
        query_faa,
        target_fna
    ]
    
    with open(output_tsv, 'w') as out:
        subprocess.run(cmd, stdout=out, check=True)
    
    print(f"Smith-Waterman search complete: {output_tsv}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Smith-Waterman search for GOI in genomic regions")
    parser.add_argument("--query", required=True, help="Query protein FASTA (GOI)")
    parser.add_argument("--target", required=True, help="Target DNA FASTA (genomic region)")
    parser.add_argument("--output", required=True, help="Output TSV (BLAST m8 format)")
    parser.add_argument("--matrix", default="BLOSUM62", choices=["BLOSUM62", "PAM100"],
                       help="Substitution matrix")
    parser.add_argument("--gap_open", type=int, default=10, help="Gap opening penalty")
    parser.add_argument("--gap_extend", type=int, default=1, help="Gap extension penalty")
    parser.add_argument("--min_score", type=int, default=50, help="Minimum alignment score")
    parser.add_argument("--min_identity", type=float, default=20.0, help="Minimum percent identity")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--method", default="auto", choices=["auto", "parasail", "ssearch36"],
                       help="Smith-Waterman implementation")
    
    args = parser.parse_args()
    
    # Decide method
    method = args.method
    if method == "auto":
        method = "parasail" if PARASAIL_AVAILABLE else "ssearch36"
    
    print(f"Using Smith-Waterman implementation: {method}", file=sys.stderr)
    
    if method == "ssearch36":
        # Use external ssearch36 program
        smith_waterman_ssearch(args.query, args.target, args.output, args.threads)
    
    elif method == "parasail":
        # Use Python parasail library (faster, in-memory)
        if not PARASAIL_AVAILABLE:
            print("ERROR: parasail not installed but requested. Install with:", file=sys.stderr)
            print("  pip install parasail", file=sys.stderr)
            sys.exit(1)
        
        # Load query proteins
        queries = list(parse_fasta(args.query))
        
        # Load target DNA and translate in 6 frames
        targets = list(parse_fasta(args.target))
        
        hits = []
        
        for q_header, q_id, q_seq in queries:
            print(f"Searching with query: {q_id}", file=sys.stderr)
            
            for t_header, t_id, t_seq in targets:
                # Translate target in all 6 frames
                frames = translate_in_six_frames(t_seq)
                
                for frame_id, prot_seq, offset in frames:
                    # Run Smith-Waterman
                    try:
                        result = smith_waterman_parasail(
                            q_seq, prot_seq, 
                            args.matrix, args.gap_open, args.gap_extend
                        )
                        
                        # Filter by thresholds
                        if result['score'] >= args.min_score and result['identity'] >= args.min_identity:
                            # Convert protein coordinates back to DNA coordinates
                            dna_start = offset + (result['end_ref'] - result['length']) * 3
                            dna_end = offset + result['end_ref'] * 3
                            
                            # BLAST m8 format output
                            hits.append({
                                'query': q_id,
                                'target': t_id,
                                'pident': result['identity'],
                                'alnlen': result['length'],
                                'mismatch': 0,  # Not calculated
                                'gapopen': 0,   # Not calculated
                                'qstart': 1,
                                'qend': result['end_query'],
                                'tstart': dna_start,
                                'tend': dna_end,
                                'evalue': 0.001,  # Placeholder
                                'bits': result['score']
                            })
                    except Exception as e:
                        print(f"Warning: SW alignment failed for {q_id} vs {t_id} frame {frame_id}: {e}",
                              file=sys.stderr)
                        continue
        
        # Write output
        with open(args.output, 'w') as f:
            for hit in hits:
                f.write(f"{hit['query']}\\t{hit['target']}\\t{hit['pident']:.1f}\\t"
                       f"{hit['alnlen']}\\t{hit['mismatch']}\\t{hit['gapopen']}\\t"
                       f"{hit['qstart']}\\t{hit['qend']}\\t{hit['tstart']}\\t{hit['tend']}\\t"
                       f"{hit['evalue']:.2e}\\t{hit['bits']:.1f}\\n")
        
        print(f"Found {len(hits)} Smith-Waterman alignments above threshold", file=sys.stderr)


if __name__ == "__main__":
    main()
