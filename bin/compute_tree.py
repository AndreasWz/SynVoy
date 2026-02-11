#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import shutil

# No BioPython needed - we count FASTA headers directly

def run_mafft(input_fasta, output_aln, threads=1):
    """
    Run MAFFT alignment.
    """
    cmd = ["mafft", "--amino", "--auto", "--thread", str(threads), input_fasta]
    with open(output_aln, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)

def run_iqtree(input_aln, output_nwk, threads=1):
    """
    Run IQ-TREE with automatic model selection and ultrafast bootstrap.
    
    Uses ModelFinder Plus (-m MFP) for automatic substitution model selection
    and ultrafast bootstrap (-B 1000) for branch support values.
    This is vastly superior to FastTree's fixed LG+Gamma model.
    """
    prefix = output_nwk.replace('.nwk', '')
    cmd = [
        "iqtree2",
        "-s", input_aln,
        "--prefix", prefix,
        "-m", "MFP",        # ModelFinder Plus (auto model selection)
        "-B", "1000",        # Ultrafast bootstrap
        "-T", str(threads),
        "--quiet",
        "--redo"             # Overwrite previous runs
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        # Try 'iqtree' (v1) if 'iqtree2' not found
        cmd[0] = "iqtree"
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    
    # IQ-TREE writes to {prefix}.treefile, rename to expected output
    treefile = prefix + ".treefile"
    if os.path.exists(treefile):
        shutil.move(treefile, output_nwk)
    
    # Cleanup IQ-TREE auxiliary files
    for ext in ['.iqtree', '.log', '.mldist', '.model.gz', '.bionj',
                '.ckp.gz', '.contree', '.splits.nex', '.uniqueseq.phy']:
        aux = prefix + ext
        if os.path.exists(aux):
            os.remove(aux)

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
    print("Inferring tree with IQ-TREE (auto model selection + ultrafast bootstrap)...")
    try:
        run_iqtree(aln_file, args.output, args.threads)
    except Exception as e:
        print(f"IQ-TREE failed: {e}")
        sys.exit(1)
        
    print(f"Tree written to {args.output}")

if __name__ == "__main__":
    main()

