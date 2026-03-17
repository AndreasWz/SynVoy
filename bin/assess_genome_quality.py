#!/usr/bin/env python3
import argparse
import sys
import json
import os

# Use our own sequence utilities (no BioPython)
try:
    from sequence_utils import parse_fasta
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta

def calculate_n50(lengths):
    sorted_lens = sorted(lengths, reverse=True)
    total = sum(lengths)
    half_total = total / 2
    
    current_sum = 0
    n50 = 0
    l50 = 0
    
    for i, l in enumerate(sorted_lens):
        current_sum += l
        if current_sum >= half_total:
            n50 = l
            l50 = i + 1
            break
            
    return n50, l50, total

def main():
    parser = argparse.ArgumentParser(description="Assess Genome Assembly Quality")
    parser.add_argument("--genome", required=True)
    parser.add_argument("--output", required=True, help="JSON output file")
    # Match nextflow.config easy-mode defaults unless explicitly overridden.
    parser.add_argument("--min_n50", type=int, default=5000,
                        help="N50 below this → bad quality (default: 5000)")
    parser.add_argument("--max_contigs", type=int, default=500000,
                        help="Contig count above this → bad quality (default: 500000)")
    
    args = parser.parse_args()
    
    lengths = []
    try:
        for _, _, rec_seq in parse_fasta(args.genome):
            lengths.append(len(rec_seq))
    except Exception as e:
        print(f"Error reading genome: {e}")
        sys.exit(1)
        
    if not lengths:
        stats = {
            "n50": 0, "l50": 0, "total_len": 0, "num_contigs": 0,
            "status": "FAIL", "msg": "Empty or invalid FASTA"
        }
    else:
        n50, l50, total = calculate_n50(lengths)
        num_contigs = len(lengths)
        
        status = "PASS"
        reasons = []
        if n50 < args.min_n50:
            reasons.append(f"N50 ({n50}) below threshold ({args.min_n50})")
        if num_contigs > args.max_contigs:
            reasons.append(f"contig count ({num_contigs}) above threshold ({args.max_contigs})")
        
        if reasons:
            status = "FAIL"
            msg = "; ".join(reasons)
        else:
            msg = "OK"
        
        stats = {
            "genome": os.path.basename(args.genome),
            "n50": n50,
            "l50": l50,
            "total_len": total,
            "num_contigs": num_contigs,
            "status": status,
            "msg": msg,
            "thresholds": {
                "min_n50": args.min_n50,
                "max_contigs": args.max_contigs,
            },
        }
        
    with open(args.output, 'w') as f:
        json.dump(stats, f, indent=2)
        
    print(f"QC Complete. Status: {stats['status']}")

if __name__ == "__main__":
    main()
