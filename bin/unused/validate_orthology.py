#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Validate orthology via Reciprocal Best Hit (RBH)")
    parser.add_argument("--candidate", required=True, help="Candidate gene FASTA")
    parser.add_argument("--home_proteins", required=True, help="Home organism proteins FASTA (BLAST/MMseqs db)")
    parser.add_argument("--query_id", required=True, help="Original Query Gene ID")
    parser.add_argument("--output", required=True, help="Validation result TSV")
    parser.add_argument("--method", default="mmseqs", help="Search method: mmseqs or blast")
    return parser.parse_args()

def check_rbh(candidate_fasta, home_db, query_id):
    """
    Run search of candidate against home_db.
    Return True if best hit is query_id.
    """
    tmp_out = f"rbh_check_{os.getpid()}.m8"
    tmp_dir = f"tmp_rbh_{os.getpid()}"
    cmd = [
        "mmseqs", "easy-search", 
        candidate_fasta, home_db, tmp_out, tmp_dir,
        "--search-type", "1", # protein vs protein
        "--format-output", "query,target,pident,evalue,bits",
        "-e", "1e-3",
        "--max-seqs", "1" # We only care about the top hit
    ]
    
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("MMseqs2 search failed.")
        return False, None
        
    best_hit = None
    if os.path.exists(tmp_out):
        with open(tmp_out) as f:
            for line in f:
                parts = line.strip().split('\t')
                best_hit = parts[1] # target is col 2
                break # only top hit
                
    if best_hit:
        # Check match (allow partial match if IDs modified)
        if query_id in best_hit or best_hit in query_id:
            return True, best_hit
            
    return False, best_hit

def main():
    args = parse_args()
    
    # Run Reciprocal Search
    is_rbh, best_hit = check_rbh(args.candidate, args.home_proteins, args.query_id)
    
    result = "ORTHOLOG" if is_rbh else "PARALOG_OR_NOISE"
    
    print(f"Validation Result: {result}")
    print(f"  Candidate: {args.candidate}")
    print(f"  Expected: {args.query_id}")
    print(f"  Found Best Hit: {best_hit}")

    with open(args.output, 'w') as f:
        f.write(f"candidate_file\texpected_id\tbest_hit_id\tstatus\n")
        f.write(f"{args.candidate}\t{args.query_id}\t{best_hit}\t{result}\n")

if __name__ == "__main__":
    main()
