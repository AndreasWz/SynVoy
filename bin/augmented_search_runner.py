#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys

# Use our own sequence utilities (no BioPython)
try:
    from sequence_utils import parse_fasta, write_fasta, load_genome, reverse_complement
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta, load_genome, reverse_complement

def run_command(cmd):
    print(f"CMD: {' '.join(cmd)}")
    subprocess.check_call(cmd)

def parse_args():
    parser = argparse.ArgumentParser(description="Run augmented search in Genomic Regions")
    parser.add_argument("--regions_bed", required=True, help="Input BED file of regions")
    parser.add_argument("--target_genome", required=True, help="Target genome FASTA")
    parser.add_argument("--query_gene", required=True, help="Query gene FASTA (protein or DNA)")
    parser.add_argument("--output_base", required=True, help="Output basename for .bed and .faa")
    parser.add_argument("--padding", type=int, default=10000, help="Padding around regions")
    parser.add_argument("--mmseqs_sens", type=str, default="8.5")
    return parser.parse_args()

def extract_regions(bed_file, genome_file, padding):
    genome_seqs = load_genome(genome_file)
    regions_seqs = []
    
    with open(bed_file) as f:
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else f"{chrom}_{start}_{end}"
            
            if chrom not in genome_seqs: continue
            
            # Pad
            slen = len(genome_seqs[chrom])
            p_start = max(0, start - padding)
            p_end = min(slen, end + padding)
            
            # Extract
            seq = genome_seqs[chrom][p_start:p_end]
            
            # ID contains coordinate info for mapping back
            # ID: name|chrom|p_start
            rid = f"{name}|{chrom}|{p_start}"
            regions_seqs.append((rid, seq))
            
    return regions_seqs
            
    return regions_seqs

def main():
    args = parse_args()
    
    # Input Validation
    if not os.path.exists(args.regions_bed):
        sys.stderr.write(f"ERROR: regions_bed file not found: {args.regions_bed}\n")
        sys.exit(1)
        
    if os.path.getsize(args.regions_bed) == 0:
        sys.stderr.write(f"INFO: regions_bed is empty. No regions to search.\n")
        # Create empty outputs for workflow consistency
        open(f"{args.output_base}.bed", 'w').close()
        open(f"{args.output_base}.fna", 'w').close()
        open(f"{args.output_base}_regions.fna", 'w').close()
        open(f"{args.output_base}_variants.faa", 'w').close()
        open(f"{args.output_base}_hits.m8", 'w').close()
        sys.exit(0)
    
    # 1. Extract Regions
    regions_fna = f"{args.output_base}_regions.fna"
    regions = extract_regions(args.regions_bed, args.target_genome, args.padding)
    write_fasta(regions, regions_fna)
    
    if not regions:
        print("No regions extracted. Exiting.")
        # Create empty outputs
        open(f"{args.output_base}.bed", 'w').close()
        open(f"{args.output_base}.fna", 'w').close()
        return

    # 2. Generate Variants
    variants_faa = f"{args.output_base}_variants.faa"
    # Call generate_variants.py
    # Assuming it's in PATH or same dir
    script_dir = os.path.dirname(os.path.realpath(__file__))
    gen_var_script = os.path.join(script_dir, "generate_variants.py")
    
    run_command([
        "python3", gen_var_script,
        "--query", args.query_gene,
        "--output", variants_faa,
        "--mutation_rate", "0.05",
        "--num_variants", "10"
    ])
    
    # 3. Search (MMseqs2 tblastn)
    hits_m8 = f"{args.output_base}_hits.m8"
    tmp_dir = f"{args.output_base}_tmp"
    
    # Check if query is DNA or Protein?
    # Usually assume Protein for "gene.fasta" in these workflows, or DNA.
    # If DNA, we should ideally translate or use search-type 3 (nucleotide->nucleotide) or 2 (translated).
    # Augmented search implies finding orthologs. Protein is best.
    # If input is DNA, we should translate it first?
    # generate_variants.py handles protein seqs mostly (AA alphabet).
    # If input is DNA, generate_variants might fail or produce DNA variants.
    # Let's assume input is Protein for now, or check?
    # The instructions say "gene sequence (DNA or amino)".
    # If DNA, we need to translate.
    # Simple check: look for T vs U? or just try to translate.
    # For now, let's trust MMseqs handles it or user provided relevant file.
    # Using search-type 2 (Protein query -> Translated Target)
    
    run_command([
        "mmseqs", "easy-search",
        variants_faa, regions_fna, hits_m8, tmp_dir,
        "--search-type", "2",
        "-s", args.mmseqs_sens,
        "--min-seq-id", "0.2", # Very relaxed
        "-e", "10", # Relaxed
        "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    ])
    
    # 4. Parse hits and map coordinates
    candidates = []
    
    if os.path.exists(hits_m8):
        with open(hits_m8) as f:
            for line in f:
                parts = line.strip().split('\t')
                # target is regions_fna ID: name|chrom|p_start
                tid_parts = parts[1].split('|')
                if len(tid_parts) < 3: continue
                
                chrom = tid_parts[1]
                p_start = int(tid_parts[2])
                
                # Hit coords in region
                h_start = int(parts[8])
                h_end = int(parts[9])
                
                # Global coords
                g_start = p_start + h_start
                g_end = p_start + h_end
                
                candidates.append({
                    'chrom': chrom,
                    'start': min(g_start, g_end), # Standardize
                    'end': max(g_start, g_end),
                    'score': float(parts[11]), # bitscore
                    'name': f"cand_{chrom}_{g_start}",
                    'strand': '+' if h_start < h_end else '-' # Crude strand guess from coords? 
                    # Actually MMseqs tstart/tend usually denotes strand. if tstart > tend, minus.
                })
                # Re-check strand logic for MMseqs2:
                # If tstart > tend, it's minus strand?
                # MMseqs2 output format: tstart always < tend?
                # Need to check documentation.
                # Usually standard BLAST output 6 puts start > end for minus.
                
                if h_start > h_end:
                     candidates[-1]['strand'] = '-'
                     candidates[-1]['start'], candidates[-1]['end'] = g_end, g_start
                
    # 5. Write candidates BED
    with open(f"{args.output_base}.bed", 'w') as f:
        for c in candidates:
            f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{c['name']}\t{c['score']}\t{c['strand']}\n")
            
    # 6. Extract Candidates FASTA
    # We can extract from genome using global coords
    # output.faa
    genome_seqs = load_genome(args.target_genome)
    cand_recs = []
    
    for c in candidates:
        if c['chrom'] in genome_seqs:
            seq = genome_seqs[c['chrom']][c['start']:c['end']]
            # Reverse complement if strand -
            if c['strand'] == '-':
                seq = reverse_complement(seq)
            
            cand_recs.append((c['name'], seq))
            
    write_fasta(cand_recs, f"{args.output_base}.fna")
    print(f"Found {len(candidates)} candidate genes.")

if __name__ == "__main__":
    main()
