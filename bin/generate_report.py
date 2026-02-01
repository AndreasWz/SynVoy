#!/usr/bin/env python3
import argparse
import json
import os
import glob

def main():
    parser = argparse.ArgumentParser(description="Generate SynTerra Final Report")
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--qc_json", help="Path to QC summary JSON")
    parser.add_argument("--output", required=True, help="Report JSON")
    
    args = parser.parse_args()
    
    report = {
        "genome_qc": [],
        "synteny_results": {
             "genes_discovered": {},
             "synteny_hits_count": {}
        },
        "summary": {}
    }
    
    # Load QC
    if args.qc_json and os.path.exists(args.qc_json):
        try:
            with open(args.qc_json) as f:
                content = f.read().strip()
                # Fix common Nextflow aggregation issues (trailing commas)
                if content.endswith(',]'):
                    content = content[:-2] + ']'
                elif content.endswith(','):
                    content = content[:-1]
                
                report["genome_qc"] = json.loads(content)
        except Exception as e:
            print(f"Warning: Failed to load QC JSON: {e}")
        
    # Scan Results - Regions (New Genes)
    genes_added_per_genome = {}
    regions_dir = os.path.join(args.results_dir, "regions")
    if os.path.exists(regions_dir):
        for f in glob.glob(os.path.join(regions_dir, "*_new_genes.faa")):
            gname = os.path.basename(f).replace("_new_genes.faa", "")
            # Strip locus prefix if present (e.g. locus_1.bed_)
            if "_GCA_" in gname or "_GCF_" in gname:
                gname = gname.split('_GCA_')[1] if '_GCA_' in gname else gname.split('_GCF_')[1]
                gname = ("GCA_" if "GCA_" in os.path.basename(f) else "GCF_") + gname
            
            count = 0
            with open(f) as fa:
                count = sum(1 for line in fa if line.startswith('>'))
            genes_added_per_genome[gname] = genes_added_per_genome.get(gname, 0) + count
            
            # Scan Results - Hits (Found via synteny)
    hits_per_genome = {}
    hits_dir = os.path.join(args.results_dir, "hits")
    if os.path.exists(hits_dir):
        for f in glob.glob(os.path.join(hits_dir, "*.m8")):
            gname = os.path.basename(f).replace(".m8", "")
            # Robust stripping: many genomes start with GCA_ or GCF_
            # If we find GCA_ or GCF_, take everything from there to the end
            if "GCA_" in gname:
                gname = "GCA_" + gname.split("GCA_")[-1]
            elif "GCF_" in gname:
                gname = "GCF_" + gname.split("GCF_")[-1]
            
            count = 0
            with open(f) as hit_f:
                count = sum(1 for line in hit_f if line.strip())
            hits_per_genome[gname] = hits_per_genome.get(gname, 0) + count

    report["synteny_results"]["genes_discovered"] = genes_added_per_genome
    report["synteny_results"]["synteny_hits_count"] = hits_per_genome
    
    report["summary"]["total_new_genes"] = sum(genes_added_per_genome.values())
    report["summary"]["genomes_with_hits"] = len(hits_per_genome)
    report["summary"]["total_hits"] = sum(hits_per_genome.values())
    
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
        
    print(f"Report generated: {args.output}")

if __name__ == "__main__":
    main()
