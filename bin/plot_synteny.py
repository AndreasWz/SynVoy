#!/usr/bin/env python3
import argparse
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import os
import colorsys

try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False

def parse_args():
    parser = argparse.ArgumentParser(description="Generate interactive synteny plot")
    parser.add_argument("--home_bed", required=True, help="Synteny block BED for home genome")
    parser.add_argument("--home_gff", help="Annotated GFF for home genome (optional, for product names)")
    parser.add_argument("--target_gffs", nargs='*', default=[], help="List of annotated GFF files for targets")
    parser.add_argument("--target_names", nargs='*', default=[], help="List of target genome names corresponding to GFFs")
    parser.add_argument("--candidate_beds", nargs='*', default=[], help="List of candidate BED files (optional, order matching targets)")
    parser.add_argument("--homology_tsvs", nargs='*', default=[], help="List of homology TSV files (optional)")
    parser.add_argument("--tree", help="Newick tree file for coloring by clade")
    parser.add_argument("--query_bed", help="BED file with query location for highlighting")
    parser.add_argument("--output", required=True, help="Output HTML file")
    return parser.parse_args()

# ... (existing functions) ...

def parse_product_map(gff_file):
    """
    Parse GFF to extract ID -> Product Name mapping.
    Handles 'Parent' in mRNA linking to 'ID' in gene/CDS, or direct on gene.
    """
    product_map = {}
    if not gff_file or not os.path.exists(gff_file) or gff_file == "NO_GFF":
        return product_map
        
    try:
        from urllib.parse import unquote
        # First pass: Link Transcript Parent -> Product
        transcript_products = {}
        
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.strip().split('\t')
                if len(p) < 9: continue
                
                feat = p[2]
                attr = p[8]
                attrs = {}
                for x in attr.split(';'):
                    if '=' in x:
                        k, v = x.split('=', 1)
                        attrs[k] = unquote(v)
                
                # Check for product description
                product = attrs.get('product')
                
                if feat == 'mRNA' and 'Parent' in attrs and product:
                    transcript_products[attrs['Parent']] = product
                elif feat == 'gene' and product:
                    # sometimes gene has product
                    if 'ID' in attrs:
                        product_map[attrs['ID']] = product
                    if 'Name' in attrs:
                         product_map[attrs['Name']] = product

        # Second Pass: If needed, or just allow gene ID lookup from transcript map
        # Actually 'Parent' of mRNA is usually the Gene ID.
        product_map.update(transcript_products)
        
    except Exception as e:
        print(f"Error parsing Home GFF for products: {e}")
        
    return product_map

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
    """
    Parse GFF to extract genes for synteny plotting.
    Handles both standard NCBI GFF and Miniprot GFF (with SynTerra_* attributes).
    
    For Miniprot GFF:
    - Uses SynTerra_Parent for home gene ID mapping (color assignment)
    - Uses SynTerra_ID for unique identification
    - Extracts Identity for quality info
    
    For standard GFF:
    - Uses ID/Name attributes as before
    """
    gene_products = {}
    genes = []
    
    try:
        from urllib.parse import unquote
        
        # First Pass: Collect Products from mRNA (standard GFF)
        # And collect genes/mRNA features
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.strip().split('\t')
                if len(p) < 9: continue
                
                feat_type = p[2]
                attr = p[8]
                source = p[1]  # e.g., "miniprot", "NCBI", etc.
                
                # Helper to parse attributes
                attrs = {}
                for x in attr.split(';'):
                    if '=' in x:
                        k, v = x.split('=', 1)
                        attrs[k] = unquote(v)
                
                if feat_type == 'mRNA' and 'Parent' in attrs and 'product' in attrs:
                    gene_products[attrs['Parent']] = attrs['product']
                
                # Handle SynTerra annotated mRNA (miniprot, augmented_search, or mmseqs2)
                if source in ('miniprot', 'augmented_search', 'mmseqs2') and feat_type == 'mRNA':
                    # mRNA is the main feature for target gene annotation
                    gene_id = attrs.get('ID', '')
                    
                    # SynTerra attributes for home gene mapping
                    synterra_parent = attrs.get('SynTerra_Parent', '')
                    synterra_id = attrs.get('SynTerra_ID', '')
                    identity = attrs.get('Identity', '0')
                    
                    # Extract base home gene name for consistent coloring
                    # SynTerra_Parent format: gene-LOC726866 or GOI_P01501
                    home_gene_base = synterra_parent
                    if home_gene_base:
                        # Remove prefix like 'gene-' for cleaner matching
                        if home_gene_base.startswith('gene-'):
                            home_gene_base = home_gene_base  # Keep as-is for lookup
                        elif home_gene_base.startswith('GOI_'):
                            home_gene_base = home_gene_base  # Keep GOI prefix
                    
                    # For display name, use SynTerra_Parent (the home gene name)
                    name = synterra_parent if synterra_parent else gene_id
                    
                    genes.append({
                        'chrom': p[0],
                        'start': int(p[3]),
                        'end': int(p[4]),
                        'name': name,
                        'id': gene_id,
                        'home_gene_id': synterra_parent,  # Key for color mapping!
                        'synterra_id': synterra_id,
                        'identity': float(identity) if identity else 0,
                        'strand': p[6],
                        'type': 'gene'  # Treat mRNA as gene for plotting
                    })
                
                elif feat_type == 'gene' or feat_type == 'CDS':
                    # Standard GFF handling
                    gene_id = attrs.get('ID', '')
                    name = attrs.get('Name', '')
                    if not name: name = gene_id
                    
                    genes.append({
                        'chrom': p[0],
                        'start': int(p[3]),
                        'end': int(p[4]),
                        'name': name,
                        'id': gene_id,
                        'home_gene_id': gene_id,  # Same as ID for standard
                        'strand': p[6],
                        'type': feat_type
                    })
                    
        # Update Names with Products (for standard GFF)
        for g in genes:
            # Try to find product using ID
            if g['id'] in gene_products:
                g['name'] = gene_products[g['id']]
            # Fallback: try removing 'gene-' prefix matches
            elif g['id'].startswith('gene-') and g['id'] in gene_products:
                g['name'] = gene_products[g['id']]
            
            # Clean up long names?
            if len(g['name']) > 30:
                g['name'] = g['name'][:27] + "..."

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

def generate_tree_colors(tree_file, color_mode='clade'):
    """
    Load tree, assign colors based on phylogenetic relationships.
    
    Args:
        tree_file: Path to Newick tree file
        color_mode: 
            'clade' - Group by monophyletic clades (similar sequences same color)
            'order' - Color by tree traversal order (gradient across tree)
    
    Returns: 
        Dict[leaf_name, hex_color]
    """
    if not ETE3_AVAILABLE:
        print("Warning: ETE3 not installed. Cannot parse tree for coloring.")
        return {}
        
    try:
        t = Tree(tree_file)
        leaves = [leaf for leaf in t.iter_leaves()]
        n_leaves = len(leaves)
        color_map = {}
        
        if n_leaves == 0:
            return {}
        
        if color_mode == 'clade' and n_leaves > 2:
            # Clade-based coloring: Group leaves by common ancestors
            # This gives similar colors to closely related sequences
            
            # Calculate pairwise distances and identify clusters
            # Using tree topology: leaves from same subtree get similar colors
            
            # Get midpoint or arbitrary internal node for reference
            midpoint = t.get_tree_root()
            
            # Get all children of root (major clades)
            root_children = midpoint.children
            
            if len(root_children) >= 2:
                # Color each major clade with a different base hue
                hue_step = 1.0 / len(root_children)
                
                for clade_idx, child in enumerate(root_children):
                    base_hue = clade_idx * hue_step
                    clade_leaves = [l for l in child.iter_leaves()]
                    
                    # Within clade, vary saturation/value
                    for leaf_idx, leaf in enumerate(clade_leaves):
                        # Slight hue variation within clade
                        if len(clade_leaves) > 1:
                            hue_var = 0.05 * (leaf_idx / (len(clade_leaves) - 1) - 0.5)
                        else:
                            hue_var = 0
                        
                        hue = (base_hue + hue_var) % 1.0
                        saturation = 0.7 + 0.2 * (leaf_idx / max(1, len(clade_leaves) - 1))
                        value = 0.9
                        
                        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                        r, g, b = [int(x * 255) for x in rgb]
                        hex_col = f"#{r:02x}{g:02x}{b:02x}"
                        clean_name = leaf.name.replace("'", "").replace('"', "")
                        color_map[clean_name] = hex_col
            else:
                # Fallback to order-based
                color_mode = 'order'
        
        if color_mode == 'order' or not color_map:
            # Simple order-based coloring (gradient across tree)
            for i, leaf in enumerate(leaves):
                hue = i / n_leaves
                rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
                r, g, b = [int(x * 255) for x in rgb]
                hex_col = f"#{r:02x}{g:02x}{b:02x}"
                clean_name = leaf.name.replace("'", "").replace('"', "")
                color_map[clean_name] = hex_col
        
        print(f"Generated colors for {len(color_map)} tree leaves (mode: {color_mode})")
        return color_map
        
    except Exception as e:
        print(f"Error parsing tree {tree_file}: {e}")
        return {}

def main():
    args = parse_args()
    
    # 0. Parse Home Products (Link ID -> Real Name)
    home_product_map = parse_product_map(args.home_gff) if args.home_gff else {}
    
    # 1. Parse Homology
    # Maps TargetGene -> HomeGeneID
    homology_map = parse_homology(args.homology_tsvs)
    
    # 1.5 Parse Query Location (for highlighting)
    query_intervals = []
    if args.query_bed and os.path.exists(args.query_bed):
        try:
             with open(args.query_bed) as f:
                 for line in f:
                     p = line.strip().split('\t')
                     if len(p) >= 3:
                         query_intervals.append({
                             'chrom': p[0],
                             'start': int(p[1]),
                             'end': int(p[2])
                         })
        except: pass
    
    # NEW: Parse Tree Colors
    tree_colors = {}
    if args.tree:
        tree_colors = generate_tree_colors(args.tree)
        if tree_colors:
            print(f"Loaded {len(tree_colors)} colors from Phylogenetic Tree.")
    
    # Data structure: list of tracks
    tracks = []
    
    # 1. Home Genome
    home_genes = parse_bed(args.home_bed)
    
    # Separate ID from Display Name for proper color mapping
    for g in home_genes:
        # Store original ID for color/tree lookup
        if 'id' not in g:
            g['id'] = g['name']
        
        # Now update display name with product description
        gid = g['id']
        if gid in home_product_map:
            g['display_name'] = home_product_map[gid]
        elif gid.replace('gene-', '') in home_product_map:
            g['display_name'] = home_product_map[gid.replace('gene-', '')]
        elif gid.startswith('gene-'):
            short_id = gid.split('gene-')[1]
            if short_id in home_product_map:
                g['display_name'] = home_product_map[short_id]
            else:
                g['display_name'] = g['name']
        else:
            g['display_name'] = g['name']
        
        # Truncate long names
        if len(g['display_name']) > 30:
            g['display_name'] = g['display_name'][:27] + "..."
             
    tracks.append({'name': 'Home Genome', 'genes': home_genes})

    # Target Genomes Parsing...
    # ...
    # COLOR LOGIC UPDATE NEEDED?
    # get_color(gene_name, homolog_name)
    # gene_name for home was ID. Now it is Product.
    # Tree likely has IDs.
    # So we should pass ID to get_color.
    
    # Need to update get_color function signature or usuage.
    
    palette = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5'
    ]
    
    # Default Palette Map (map original IDs)
    default_home_colors = {}
    for i, g in enumerate(home_genes):
        default_home_colors[g['id']] = palette[i % len(palette)]
    
    def get_color(gene_id, home_gene_id=None):
        """
        Get color for gene based on its home gene ID.
        
        Args:
            gene_id: The local gene ID (e.g., MP000001)
            home_gene_id: The corresponding home genome gene ID (e.g., gene-LOC726866)
                         This is the key for consistent coloring across genomes.
        """
        # Use home_gene_id as primary key for consistent coloring
        lookup_id = home_gene_id if home_gene_id else gene_id
        
        # 1. Try Tree Color (by ID)
        if lookup_id in tree_colors:
            return tree_colors[lookup_id]
        
        # Try without gene- prefix
        clean_id = lookup_id.replace('gene-', '') if lookup_id else ''
        if clean_id and clean_id in tree_colors:
            return tree_colors[clean_id]
            
        # 2. Try Default Palette (by Home Gene ID)
        if lookup_id in default_home_colors:
            return default_home_colors[lookup_id]
        
        # Try with gene- prefix if not present
        if lookup_id and not lookup_id.startswith('gene-'):
            prefixed_id = f"gene-{lookup_id}"
            if prefixed_id in default_home_colors:
                return default_home_colors[prefixed_id]
             
        # 3. Fallback gray
        return '#cccccc'

    # ...
    
    # 2. Target Genomes Loop
    if len(args.target_gffs) != len(args.target_names):
        print("Error: Number of GFFs and Names must match.", file=sys.stderr)
        length = min(len(args.target_gffs), len(args.target_names))
    else:
        length = len(args.target_gffs)
    
    # Validate other arrays
    if args.candidate_beds and len(args.candidate_beds) != length:
        print(f"Warning: Expected {length} candidate BEDs but got {len(args.candidate_beds)}. Some may be missing.", file=sys.stderr)
    
    if args.homology_tsvs and len(args.homology_tsvs) != length:
        print(f"Warning: Expected {length} homology TSVs but got {len(args.homology_tsvs)}. Some may be missing.", file=sys.stderr)

    for i in range(length):
        gff_file = args.target_gffs[i]
        name = args.target_names[i]
        # parse_gff now extracts home_gene_id from SynTerra_Parent attribute
        genes = parse_gff(gff_file) 
        
        plot_genes = [g for g in genes if g['type'] == 'gene']
        if not plot_genes:
             plot_genes = genes
        
        # Determine Chromosome(s)
        chroms = sorted(list(set(g['chrom'] for g in plot_genes)))
        chrom_str = ", ".join(chroms[:2]) # Limit to 2
        if len(chroms) > 2: chrom_str += "..."
        
        full_title = f"{name} ({chrom_str})" if chrom_str else name
        
        tracks.append({'name': full_title, 'genes': plot_genes, 'cands': []})
        
        # Update cands logic safely
        if args.candidate_beds and i < len(args.candidate_beds):
            cand_file = args.candidate_beds[i]
            # Handle Nextflow's empty list passing as needed, or file existence
            if cand_file != "NO_CANDIDATES" and os.path.exists(cand_file):
                 tracks[-1]['cands'] = parse_bed(cand_file)

    # Dynamic vertical spacing based on number of tracks
    n_tracks = len(tracks)
    vertical_spacing = min(0.05, 0.8 / max(1, n_tracks - 1)) if n_tracks > 2 else 0.1
    
    fig = make_subplots(rows=n_tracks, cols=1, shared_xaxes=False, vertical_spacing=vertical_spacing, subplot_titles=[t['name'] for t in tracks])
    
    for i, track in enumerate(tracks):
        genes = track['genes']
        row = i + 1
        
        if not genes:
            continue
            
        min_start = min(g['start'] for g in genes)
        
        for g in genes:
            # Determine Color using home_gene_id (from SynTerra_Parent attribute)
            gid = g.get('id', g['name'])
            
            # For Miniprot GFF: home_gene_id is set directly from SynTerra_Parent
            # For standard GFF or home: home_gene_id equals gid
            home_gene_id = g.get('home_gene_id', gid)
            
            # Also check homology_map as fallback (for old TSV format)
            if not home_gene_id or home_gene_id == gid:
                mapped_id = homology_map.get(gid)
                if mapped_id:
                    home_gene_id = mapped_id
            
            # Resolve Homolog Display Name
            homolog_display = home_gene_id
            if home_gene_id:
                 # Clean ID if needed
                 clean_hid = home_gene_id.replace('gene-', '')
                 if home_gene_id in home_product_map:
                     homolog_display = home_product_map[home_gene_id]
                 elif clean_hid in home_product_map:
                     homolog_display = home_product_map[clean_hid]
                 
                 # Truncate
                 if len(homolog_display) > 30:
                     homolog_display = homolog_display[:27] + "..."

            # Get color using home_gene_id for consistent cross-genome coloring
            color = get_color(gid, home_gene_id)

            # Hover Text - use display_name if available, otherwise name
            display_name = g.get('display_name', g['name'])
            if i == 0:
                 hover_text = f"{display_name} (Home)"
            else:
                if home_gene_id and home_gene_id != gid:
                    hover_text = f"{display_name}<br>Homolog: {homolog_display}"
                else:
                    hover_text = f"{display_name}<br>No Homolog Found"
            
            # Coordinates
            start = g['start']
            end = g['end']
            
            # Check for Query Overlap (Home Track Only)
            is_query = False
            if i == 0 and query_intervals:
                for q in query_intervals:
                    if g['chrom'] == q['chrom']:
                        # Overlap logic
                        if not (g['end'] < q['start'] or g['start'] > q['end']):
                            is_query = True
                            break

            # Add Trace
            line_dict = dict(color=color, width=0)
            if is_query:
                line_dict = dict(color="black", width=2.5) # Thick black border for Query
                
            fig.add_trace(go.Scatter(
                x=[start, end, end, start, start],
                y=[0, 0, 1, 1, 0],
                fill="toself",
                mode='lines',
                name=g['name'],
                text=f"{hover_text}<br>TYPE: {'QUERY' if is_query else 'Gene'}<br>{g['chrom']}:{start}-{end} ({g['strand']})",
                line=line_dict,
                fillcolor=color,
                showlegend=False
            ), row=row, col=1)
            
            # Add Label
            display_name = g.get('display_name', g['name'])
            label_text = display_name
            font_dict = dict(size=10, color="black")
            
            if i > 0 and home_gene_id:
                label_text = homolog_display
            
            if is_query:
                label_text = "★ " + label_text
                font_dict = dict(size=12, color="red") # Highlight label
                
            fig.add_annotation(
                x=(start + end)/2,
                y=0.5,
                text=label_text,
                showarrow=False,
                font=font_dict,
                row=row, col=1
            )

        # Plot Candidates (if any)
        if 'cands' in track and track['cands']:
             for cand in track['cands']:
                  start = cand['start']
                  end = cand['end']
                  
                  # Candidates Style: Red, distinct
                  color = "#FF0000"
                  
                  fig.add_trace(go.Scatter(
                    x=[start, end, end, start, start],
                    y=[0, 0, 1, 1, 0],
                    fill="toself",
                    mode='lines',
                    name="Candidate",
                    text=f"CANDIDATE: {cand['name']}<br>{track['name']}:{start}-{end}",
                    line=dict(color="red", width=2, dash="dot"),
                    fillcolor="rgba(255, 0, 0, 0.5)",
                    showlegend=False
                ), row=row, col=1)
                
                  fig.add_annotation(
                    x=(start + end)/2,
                    y=-0.2, # Below the gene track
                    text="★ CANDIDATE",
                    showarrow=False,
                    font=dict(size=9, color="red"),
                    row=row, col=1
                )

    # Dynamic height: minimum 200px per track, maximum 400px per track
    track_height = max(150, min(300, 2000 // max(1, n_tracks)))
    fig_height = track_height * n_tracks + 100  # +100 for title/margins
    
    fig.update_layout(
        height=fig_height, 
        title_text="Synteny Plot (Phylo-Colored)",
        showlegend=False
    )
    
    # Hide y-axes for cleaner look (gene bars use y=0-1 range)
    for i in range(1, n_tracks + 1):
        yaxis_name = f'yaxis{i}' if i > 1 else 'yaxis'
        fig.update_layout(**{yaxis_name: dict(showticklabels=False, showgrid=False, zeroline=False)})
    
    fig.write_html(args.output)
    print(f"Plot saved to {args.output}")

if __name__ == "__main__":
    main()
