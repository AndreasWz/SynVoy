#!/usr/bin/env python3
import sys
import json
import argparse
import os

def parse_bed(filepath):
    genes = []
    if not os.path.exists(filepath):
        return genes
    with open(filepath) as fh:
        for line in fh:
            p = line.strip().split('\t')
            if len(p) >= 4:
                genes.append({
                    "chrom": p[0],
                    "start": int(p[1]),
                    "end": int(p[2]),
                    "name": p[3],
                })
    return genes

def main():
    parser = argparse.ArgumentParser("Analyze home locus and derive dynamic parameters for iterative search")
    parser.add_argument("--goi_info", required=True)
    parser.add_argument("--flanking_bed", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    goi_info = {}
    if os.path.exists(args.goi_info):
        with open(args.goi_info) as f:
            goi_info = json.load(f)
    
    flanking = parse_bed(args.flanking_bed)

    # Base parameters
    params = {
        "mmseqs_sensitivity": 6.0,
        "min_synteny_score": 0.2, # default baseline
        "flank_search_radius": 1000000,
        "locus_notes": []
    }

    # Evaluate GOI characteristics
    exons = goi_info.get("exons", [])
    num_exons = len(exons)

    if num_exons > 15:
        params["locus_notes"].append(f"Highly fragmented GOI ({num_exons} exons), boosting sensitivity")
        params["mmseqs_sensitivity"] = 7.5
        params["min_synteny_score"] = 0.15
        
    if num_exons > 0 and num_exons <= 2:
        params["locus_notes"].append(f"Small/compact GOI ({num_exons} exons), relying more on flanking synteny")
        params["mmseqs_sensitivity"] = 5.5
        params["min_synteny_score"] = 0.25

    # Check for clusters in flanking genes (naive name-based)
    q_id = goi_info.get("query_id", "GOI").lower()
    cluster_count = sum(1 for g in flanking if q_id in g["name"].lower() or "like" in g["name"].lower())
    if cluster_count > 1:
        params["locus_notes"].append(f"Potential tandem cluster detected ({cluster_count} related genes), expanding search radius")
        params["flank_search_radius"] = 2000000
        params["min_synteny_score"] = 0.25 # Require stronger synteny to avoid paralog bleeding

    print("\n".join(params["locus_notes"]))

    with open(args.output, "w") as out:
        json.dump(params, out, indent=2)

if __name__ == "__main__":
    main()
