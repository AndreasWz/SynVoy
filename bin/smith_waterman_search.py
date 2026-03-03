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
from sequence_utils import parse_fasta, translate, reverse_complement

try:
    import parasail
    HAS_PARASAIL = True
except ImportError:
    HAS_PARASAIL = False
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
    if not HAS_PARASAIL:
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
    Perform iterative Smith-Waterman search using Parasail (if available) or ssearch36.
    
    If Parasail is available:
        - Translates target DNA in 6 frames.
        - Runs striped SW (protein-protein).
        - Iteratively masks hits to find secondary alignments (tandem dupes).
        - Maps coordinates back to genomic DNA.
        
    If Parasail is missing:
        - Falls back to ssearch36 (external binary) with iterative masking loop.
    """
    import subprocess
    import shutil
    
    if HAS_PARASAIL:
        run_parasail_sw(query_faa, target_fna, output_tsv)
        return

    # Check if ssearch36 is available
    try:
        subprocess.run(['ssearch36', '-h'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print("ERROR: ssearch36 not found and parasail-python not installed.", file=sys.stderr)
        print("  Please install one of them:", file=sys.stderr)
        print("  pip install parasail", file=sys.stderr)
        print("  conda install -c bioconda fasta3", file=sys.stderr)
        sys.exit(1)
    
    # Read original target sequence
    targets = list(parse_fasta(target_fna))
    if not targets:
        return
        
    # We assume one target sequence for simplicity in this context
    # (Checking one region against one query set)
    t_header, t_id, t_seq_orig = targets[0]
    t_seq_mutable = list(t_seq_orig) # Mutable list of chars
    
    all_hits_lines = []
    
    # Iteration loop
    max_iter = 20
    tmp_target = f"{output_tsv}.tmp.fna"
    
    for i in range(max_iter):
        # Write current masked target
        with open(tmp_target, 'w') as f:
            f.write(f">{t_id}\n{''.join(t_seq_mutable)}\n")
            
        # Run ssearch36
        cmd = [
            'ssearch36',
            '-m', '8',  # Tabular output
            '-T', str(threads),
            '-E', '20000',  # Relaxed E-value
            '-3',  # Query is protein, target is DNA
            query_faa,
            tmp_target
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout
        except subprocess.CalledProcessError as e:
            print(f"ssearch36 failed at iteration {i}: {e}", file=sys.stderr)
            break
            
        if not output.strip():
            break # No more hits
            
        # Parse hits from this iteration
        new_hits_found = False
        lines = output.strip().split('\n')
        for line in lines:
            if not line.strip() or line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 12: continue
            
            # Hit found!
            new_hits_found = True
            all_hits_lines.append(line)
            
            # Mask region
            try:
                ts = int(parts[8])
                te = int(parts[9])
                start = min(ts, te)
                end = max(ts, te)
                
                # 1-based to 0-based
                start_0 = max(0, start - 1)
                end_0 = min(len(t_seq_mutable), end)
                
                # Mask with N
                for k in range(start_0, end_0):
                    t_seq_mutable[k] = 'N'
                    
            except ValueError:
                continue
        
        if not new_hits_found:
            break
            
    # Cleanup
    if os.path.exists(tmp_target):
        os.remove(tmp_target)
        
    # Write aggregated hits
    with open(output_tsv, 'w') as out:
        for line in all_hits_lines:
            out.write(line + '\n')
    
    print(f"Smith-Waterman search complete (ssearch36): {output_tsv} ({len(all_hits_lines)} hits)", file=sys.stderr)


def run_parasail_sw(query_faa, target_fna, output_tsv):
    """
    In-memory iterative Smith-Waterman using Parasail.
    Uses 'Best Hit per Iteration' strategy:
    1. Iterate searching all frames.
    2. Pick BEST hit globally.
    3. Mask that region in ALL frames.
    4. Repeat.
    """
    # 1. Load Sequences
    queries = list(parse_fasta(query_faa))
    targets = list(parse_fasta(target_fna))
    
    if not queries or not targets:
        return

    # Assume single target region
    t_header, t_id, t_seq_dna = targets[0]
    
    # 2. Translate Target (6 frames)
    # Store: {'seq': mutable_list_aa, 'frame': int, 'strand': str, 'offset': int}
    frames = []
    
    # Forward (+1, +2, +3)
    for i in range(3):
        seq_trans = list(translate(t_seq_dna[i:]))
        frames.append({'seq': seq_trans, 'frame': i+1, 'strand': '+', 'offset': i})
        
    # Reverse (-1, -2, -3)
    # RevComp DNA first
    rc_dna = reverse_complement(t_seq_dna)
    for i in range(3):
        seq_trans = list(translate(rc_dna[i:]))
        frames.append({'seq': seq_trans, 'frame': -(i+1), 'strand': '-', 'offset': i})
        
    all_hits = []
    
    # 3. For each query
    for q_head, q_id, q_seq in queries:
        try:
            # Use 32-bit profile
            profile = parasail.profile_create_32(q_seq, parasail.blosum62)
        except Exception as e:
            print(f"Parasail profile creation failed for {q_id}: {e}", file=sys.stderr)
            continue
            
        # Iteration Loop (Find multiple non-overlapping hits)
        for iteration in range(20): 
            best_score = 0
            best_result = None
            best_frame_idx = -1
            
            # Search all frames
            for f_idx, frame in enumerate(frames):
                seq_str = "".join(frame['seq'])
                # sw_trace_striped_profile_32
                result = parasail.sw_trace_striped_profile_32(profile, seq_str, 11, 1)
                
                if result.score > best_score:
                    best_score = result.score
                    best_result = result
                    best_frame_idx = f_idx
            
            # Check threshold
            if best_score < 40: # Raw score threshold
                break
                
            # Process Best Hit
            result = best_result
            frame = frames[best_frame_idx]
            
            # CIGAR Parsing for coordinates
            try: 
                cigar_decoded = result.cigar.decode.decode()
            except:
                break
                
            # Parse CIGAR to find start and alignment length
            # We must handle leading 'D's which shift the logical start
            cigar_ops = []
            curr_num = ""
            for char in cigar_decoded:
                if char.isdigit():
                    curr_num += char
                else:
                    n = int(curr_num) if curr_num else 1
                    curr_num = ""
                    cigar_ops.append((char, n))
            
            # Calculate total trace length
            t_len_aln = 0
            q_len_aln = 0
            aln_len = 0
            
            for op, n in cigar_ops:
                aln_len += n
                if op in ['M', '=', 'X', 'D']:
                    t_len_aln += n
                if op in ['M', '=', 'X', 'I']:
                    q_len_aln += n

            t_end_aa = result.end_ref
            q_end_aa = result.end_query
            
            t_start_aa = t_end_aa - t_len_aln + 1
            q_start_aa = q_end_aa - q_len_aln + 1
            
            # Adjust start for leading gaps (D or I)
            # Leading D means we skipped Ref bases -> t_start shifts forward
            # Leading I means we skipped Query bases -> q_start shifts forward
            
            # Trim leading D (Deletion from Ref)
            # Usually Parasail S-W trace starts with D if it skips beginning of ref?
            # We iterate ops to adjust starts
            for op, n in cigar_ops:
                if op == 'D':
                    t_start_aa += n # Advance match start on Ref
                elif op == 'I':
                    q_start_aa += n # Advance match start on Query
                else:
                    break # First match/sub stops the trimming

            
            # DEBUG
            # print(f"DEBUG_HIT: Iter={iteration} Frame={frame['frame']} Score={result.score}")
            # print(f"DEBUG_HIT: Cigar={cigar_decoded} T_Range=[{t_start_aa}, {t_end_aa}]")
            
            # DNA Coordinates of the Hit
            if frame['strand'] == '+':
                t_start_dna = t_start_aa * 3 + frame['offset']
                t_end_dna = t_end_aa * 3 + frame['offset'] + 2
                ts = t_start_dna + 1
                te = t_end_dna + 1
            else:
                rc_start = t_start_aa * 3 + frame['offset']
                rc_end = t_end_aa * 3 + frame['offset'] + 2
                L = len(t_seq_dna)
                te_fwd = L - rc_start
                ts_fwd = L - rc_end
                ts = ts_fwd
                te = te_fwd
                
                # For masking, we need the internal "frame coordinates" for every frame.
                # Simplest: Define Genomic Interval [min_dna, max_dna]
                t_start_dna = min(ts, te) - 1
                t_end_dna = max(ts, te) - 1 + 1 # exclusive upper bound?
            
            # print(f"DEBUG_DNA: [{t_start_dna}, {t_end_dna}]")

            # Report
            hit_line = f"{q_id}\t{t_id}\t99.9\t{aln_len}\t0\t0\t{q_start_aa+1}\t{q_end_aa+1}\t{ts}\t{te}\t1e-10\t{result.score}"
            all_hits.append(hit_line)
            
            # MASKING: Mask this genomic region in ALL frames
            # Genomic Interval: (ts-1) to (te) (0-based)
            g_start = min(ts, te) - 1
            g_end = max(ts, te)
            
            # print(f"DEBUG_MASKING: Interval [{g_start}, {g_end}]")
            
            for f_i, f in enumerate(frames):
                # Map genomic [g_start, g_end] to frame AA coords
                # Frame Fwd: AA_idx = (DNA_idx - offset) / 3
                # We mask conservatively: any AA that overlaps the genomic region
                
                if f['strand'] == '+':
                    start_i = int((g_start - f['offset'] - 2)/3)
                    end_i = int((g_end - f['offset'])/3) + 1 # Exclusive
                    
                else:
                    rc_g_start = len(t_seq_dna) - g_end
                    rc_g_end = len(t_seq_dna) - g_start
                    
                    start_i = int((rc_g_start - f['offset'] - 2)/3)
                    end_i = int((rc_g_end - f['offset'])/3) + 1
                
                # Clamp
                start_i = max(0, start_i)
                end_i = min(len(f['seq']), end_i)
                
                # DEBUG BEFORE
                # if f_i == 1: # Check Frame 2 specifically
                #     sample = "".join(f['seq'][30:40])
                #     print(f"DEBUG_F2_BEFORE: ...{sample}...")

                # Apply Mask
                changed = False
                for k in range(start_i, end_i):
                    f['seq'][k] = 'X'
                    changed = True
                
                # if f_i == 1:
                #      sample = "".join(f['seq'][30:40])
                #      print(f"DEBUG_F2_AFTER: ...{sample}... Changed={changed} Range=[{start_i}, {end_i}]")
                     
    # Write output
                    
    # Write output
    with open(output_tsv, 'w') as f:
        for line in all_hits:
            f.write(line + "\n")
    
    print(f"Smith-Waterman search complete (Parasail): {output_tsv} ({len(all_hits)} hits)", file=sys.stderr)


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
        method = "parasail" if HAS_PARASAIL else "ssearch36"
    
    print(f"Using Smith-Waterman implementation: {method}", file=sys.stderr)
    
    if method == "ssearch36":
        # Use external ssearch36 program
        smith_waterman_ssearch(args.query, args.target, args.output, args.threads)
    
    elif method == "parasail":
        # Use Python parasail library (faster, in-memory)
        if not HAS_PARASAIL:
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
                f.write(f"{hit['query']}\t{hit['target']}\t{hit['pident']:.1f}\t"
                       f"{hit['alnlen']}\t{hit['mismatch']}\t{hit['gapopen']}\t"
                       f"{hit['qstart']}\t{hit['qend']}\t{hit['tstart']}\t{hit['tend']}\t"
                       f"{hit['evalue']:.2e}\t{hit['bits']:.1f}\n")
        
        print(f"Found {len(hits)} Smith-Waterman alignments above threshold", file=sys.stderr)


if __name__ == "__main__":
    main()
