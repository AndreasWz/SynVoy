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
                
                evalue = float(parts[4]) if len(parts) > 4 and parts[4] not in ('.', '') else float('inf')
                entries.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'evalue': evalue,
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
    
    # Rank clusters by significance (best e-value in cluster)
    for cl in clusters:
        cl_evalues = [e['evalue'] for e in cl]
        cl.sort(key=lambda x: x['evalue'])  # sort within cluster for reference
    
    clusters.sort(key=lambda cl: min(e['evalue'] for e in cl))
    
    # Filter: keep only clusters with sufficient evidence
    # Primary locus = best cluster; secondary loci must have meaningful hits
    if len(clusters) > 1:
        primary_best = min(e['evalue'] for e in clusters[0])
        filtered = [clusters[0]]
        for cl in clusters[1:]:
            best_ev = min(e['evalue'] for e in cl)
            n_hits = len(cl)
            # Keep secondary locus if:
            # 1. It has at least 2 significant hits, OR
            # 2. It has 1 hit with e-value within 1e6 of primary (likely real homolog)
            if n_hits >= 2:
                filtered.append(cl)
            elif best_ev <= 1e-10 and (primary_best == 0 or best_ev / max(primary_best, 1e-200) < 1e6):
                filtered.append(cl)
            else:
                print(f"Filtered out weak locus on {cl[0]['chrom']}: "
                      f"{n_hits} hit(s), best e-value={best_ev:.2e}")
        clusters = filtered
        
    print(f"Found {len(clusters)} distinct loci.")
    
    for i, cl in enumerate(clusters):
        # We name the file so Nextflow can pick it up.
        # Use simple indexing.
        out_name = f"{args.output_prefix}_{i+1}.bed"
        with open(out_name, 'w') as f_out:
            for entry in cl:
                f_out.write(entry['line'] + "\n")
        best_ev = min(e['evalue'] for e in cl)
        print(f"Wrote locus {i+1} to {out_name} ({len(cl)} hits, best e-value={best_ev:.2e})")

if __name__ == "__main__":
    main()
