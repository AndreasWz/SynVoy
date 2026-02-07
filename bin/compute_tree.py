#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

# No BioPython needed - we count FASTA headers directly

def run_mafft(input_fasta, output_aln, threads=1):
    """
    Run MAFFT alignment.
    """
    cmd = ["mafft", "--amino", "--auto", "--thread", str(threads), input_fasta]
    with open(output_aln, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)

def run_fasttree(input_aln, output_nwk):
    """
    Run FastTree to infer phylogeny.
    """
    cmd = ["FastTree", "-lg", "-gamma"] # standard protein model
    # FastTree reads from stdin or file, prints to stdout
    
    with open(output_nwk, "w") as out_f:
        # Pass input_aln as argument or redirect? FastTree accepts filename.
        # fasttree alignment.fasta > tree.nwk
        subprocess.run(cmd + [input_aln], stdout=out_f, check=True)

def main():
    parser = argparse.ArgumentParser(description="Compute Phylogenetic Tree from Fasta")
    parser.add_argument("--input", required=True, help="Input protein fasta (candidates)")
    parser.add_argument("--output", required=True, help="Output Newick tree file")
    parser.add_argument("--threads", type=int, default=1)
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Input file {args.input} not found.")
        sys.exit(1)
        
    # Check if we have enough sequences
    count = 0
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"): count += 1
            
    if count < 3:
        print("Not enough sequences to build a tree (<3).")
        # Write a stub tree so downstream pipeline steps don't fail
        with open(args.output, 'w') as f:
            f.write("(placeholder:0.0);\n")
        sys.exit(0)

    # 1. Align
    aln_file = args.input + ".aln"
    print(f"Aligning {count} sequences with MAFFT...")
    try:
        run_mafft(args.input, aln_file, args.threads)
    except Exception as e:
        print(f"MAFFT failed: {e}")
        sys.exit(1)
        
    # 2. Tree
    print("Inferring Tree with FastTree...")
    try:
        run_fasttree(aln_file, args.output)
    except Exception as e:
        print(f"FastTree failed: {e}")
        sys.exit(1)
        
    print(f"Tree written to {args.output}")

if __name__ == "__main__":
    main()
