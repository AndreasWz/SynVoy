#!/usr/bin/env python3
import argparse
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Generate interactive synteny plot")
    parser.add_argument("--home_bed", required=True, help="Synteny block BED for home genome")
    parser.add_argument("--target_gffs", nargs='+', required=True, help="List of annotated GFF files for targets")
    parser.add_argument("--target_names", nargs='+', required=True, help="List of target genome names corresponding to GFFs")
    parser.add_argument("--candidate_beds", nargs='+', help="List of candidate BED files (optional, order matching targets)")
    parser.add_argument("--homology_tsvs", nargs='+', help="List of homology TSV files (optional)")
    parser.add_argument("--output", required=True, help="Output HTML file")
    return parser.parse_args()

def parse_bed(bed_file):
    genes = []
    try:
        with open(bed_file) as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) < 6: continue
                genes.append({
                    'chrom': p[0],
                    'start': int(p[1]),
                    'end': int(p[2]),
                    'name': p[3],
                    'strand': p[5],
                    'type': 'gene'
                })
    except Exception as e:
        print(f"Error reading BED {bed_file}: {e}")
    return genes

def parse_gff(gff_file):
    genes = []
    try:
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.strip().split('\t')
                if len(p) < 9: continue
                if p[2] != 'gene' and p[2] != 'CDS': continue # Augustus uses gene/CDS
                
                # Extract simple name (e.g. from ID=... or Name=...)
                attr = p[8]
                name = "gene"
                if "ID=" in attr:
                    name = attr.split("ID=")[1].split(";")[0]
                elif "Name=" in attr:
                    name = attr.split("Name=")[1].split(";")[0]
                elif "Parent=" in attr:
                    name = attr.split("Parent=")[1].split(";")[0]
                    
                genes.append({
                    'chrom': p[0],
                    'start': int(p[3]),
                    'end': int(p[4]),
                    'name': name,
                    'strand': p[6],
                    'type': p[2]
                })
    except Exception as e:
        print(f"Error reading GFF {gff_file}: {e}")
    return genes

def parse_homology(tsv_files):
    """
    Parse homology TSV files. 
    Format: TargetGene \t HomeGene
    Returns: Dict[TargetGene, HomeGene]
    """
    mapping = {}
    if not tsv_files:
        return mapping
        
    for tsv in tsv_files:
        if tsv == "NO_HOMOLOGY": continue
        try:
            with open(tsv) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        mapping[parts[0]] = parts[1]
        except Exception as e:
            print(f"Error reading homology {tsv}: {e}")
    return mapping

def main():
    args = parse_args()
    
    # 1. Parse Homology
    homology_map = parse_homology(args.homology_tsvs)
    
    # Data structure: list of tracks
    tracks = []
    
    # 1. Home Genome
    home_genes = parse_bed(args.home_bed)
    tracks.append({'name': 'Home Genome', 'genes': home_genes})
    
    # Generate Color Map for Home Genes
    # Assign a unique color to each Home Gene
    # Cycle through a palette
    palette = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5'
    ]
    
    home_gene_colors = {}
    for i, g in enumerate(home_genes):
        home_gene_colors[g['name']] = palette[i % len(palette)]
    
    # 2. Target Genomes
    if len(args.target_gffs) != len(args.target_names):
        print("Error: Number of GFFs and Names must match.")
        length = min(len(args.target_gffs), len(args.target_names))
    else:
        length = len(args.target_gffs)

    for i in range(length):
        gff_file = args.target_gffs[i]
        name = args.target_names[i]
        genes = parse_gff(gff_file)
        
        # Filter GFF to only keep genes for plotting tracks (ignore CDS/exons for now or map them)
        # Augustus outputs 'gene' features.
        plot_genes = [g for g in genes if g['type'] == 'gene']
        if not plot_genes:
             # Fallback if no gene features, use CDS but grouping might be needed
             # For now, let's assume 'gene' exists or take CDS as genes
             plot_genes = genes
        
        cands = []
        if args.candidate_beds and i < len(args.candidate_beds):
            cand_file = args.candidate_beds[i]
            # Handle Nextflow's empty list passing as needed, or file existence
            if cand_file != "NO_CANDIDATES" and os.path.exists(cand_file):
                 cands = parse_bed(cand_file)
        
        tracks.append({'name': name, 'genes': plot_genes, 'cands': cands})
        
    # Plotting
    fig = make_subplots(rows=len(tracks), cols=1, shared_xaxes=False, vertical_spacing=0.1, subplot_titles=[t['name'] for t in tracks])
    
    for i, track in enumerate(tracks):
        genes = track['genes']
        row = i + 1
        
        if not genes:
            continue
            
        min_start = min(g['start'] for g in genes)
        
        for g in genes:
            # Determine Color
            color = '#cccccc' # Default Grey
            
            # Identify gene name for coloring
            # If Home Genome (i==0), use its own name
            if i == 0:
                color = home_gene_colors.get(g['name'], '#cccccc')
                hover_text = f"{g['name']} (Home)"
            else:
                # If Target Genome, check homology
                target_name = g['name']
                if target_name in homology_map:
                    homolog = homology_map[target_name]
                    color = home_gene_colors.get(homolog, '#cccccc')
                    hover_text = f"{target_name}<br>Homolog: {homolog}"
                else:
                    hover_text = f"{target_name}<br>No Homolog Found"
            
            # Coordinates
            start = g['start']
            end = g['end']
            
            # Add Trace
            fig.add_trace(go.Scatter(
                x=[start, end, end, start, start],
                y=[0, 0, 1, 1, 0],
                fill="toself",
                mode='lines',
                name=g['name'],
                text=f"{hover_text}<br>{g['chrom']}:{start}-{end} ({g['strand']})",
                line=dict(color=color, width=0),
                fillcolor=color,
                showlegend=False
            ), row=row, col=1)
            
            # Add Label
            fig.add_annotation(
                x=(start + end)/2,
                y=0.5,
                text=g['name'] if i==0 else (homology_map.get(g['name'], g['name'])), # Label with Home Name if possible
                showarrow=False,
                font=dict(size=10, color="black"),
                row=row, col=1
            )

    fig.update_layout(height=300 * len(tracks), title_text="Synteny Plot (Homology Colored)")
    fig.write_html(args.output)
    print(f"Plot saved to {args.output}")

if __name__ == "__main__":
    main()
