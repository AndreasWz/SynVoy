#!/usr/bin/env python3
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Split BED file into multiple files based on locus clustering.")
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output BED files")
    parser.add_argument("--dist_threshold", type=int, default=50000, help="Distance threshold for clustering (bp)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    entries = []
    try:
        with open(args.bed) as f:
            for line in f:
                if not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                
                entries.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'line': line.strip()
                })
    except FileNotFoundError:
        print(f"Error: File {args.bed} not found.")
        return

    if not entries:
        print("No entries found.")
        return

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
        
    print(f"Found {len(clusters)} distinct loci.")
    
    for i, cl in enumerate(clusters):
        # We name the file so Nextflow can pick it up.
        # Use simple indexing.
        out_name = f"{args.output_prefix}_{i+1}.bed"
        with open(out_name, 'w') as f_out:
            for entry in cl:
                f_out.write(entry['line'] + "\n")
        print(f"Wrote locus {i+1} to {out_name} ({len(cl)} hits)")

if __name__ == "__main__":
    main()
