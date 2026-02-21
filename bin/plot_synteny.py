#!/usr/bin/env python3
"""
plot_synteny.py  –  Interactive synteny visualization for SynTerra

Layout
──────
  •  Home genome at top, target genomes below (ordered by phylogenetic distance)
  •  Gene arrows (pentagons) coloured by homology group
  •  Connecting ribbons between homologous genes in adjacent tracks
  •  GOI highlighted with warm/red clade colours from the phylogenetic tree
  •  Flanking genes share a consistent colour derived from the home-genome name

Inputs
──────
  --home_bed        Synteny-block BED for the home genome
  --home_gff        NCBI GFF for the home genome (product-name lookup)
  --query_bed       BED file with query-gene location (GOI identification)
  --target_gffs     Target-genome GFFs (SynTerra exon_annotation format)
  --target_names    Display names (optional – derived from GFF filename if absent)
  --candidate_beds  Cluster-region BED files (used to filter target genes to candidate loci)
  --homology_tsvs   Homology TSV files (target -> home mapping, fallback)
  --tree            Newick tree for GOI clade colouring + target ordering

Output
──────
  --output          Interactive HTML file (Plotly)
"""

import argparse
import colorsys
import json
import os
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

import plotly.graph_objects as go

try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False


# ======================================================================
# Colour palettes
# ======================================================================

# Tableau-20 style qualitative palette for flanking genes
GENE_PALETTE = [
    "#4e79a7", "#f28e2b", "#59a14f", "#b07aa1", "#76b7b2",
    "#edc948", "#ff9da7", "#9c755f", "#86bcb6", "#e15759",
    "#8cd17d", "#499894", "#d4a6c8", "#a0cbe8", "#ffbe7d",
    "#d37295", "#fabfd2", "#b6992d", "#7b848f", "#f1ce63",
]

GOI_COLOUR    = "#e31a1c"   # bright red (default for GOI)
GOI_BORDER    = "#8b0000"   # dark red
UNMATCHED_CLR = "#d9d9d9"   # light gray
TRACK_BG_CLR  = "#f5f5f5"   # very light gray track background


# ======================================================================
# Parsing helpers
# ======================================================================

def parse_bed(bed_file):
    """Parse a BED file -> list of dicts with chrom/start/end/name/strand."""
    genes = []
    if not bed_file or not os.path.exists(bed_file):
        return genes
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 4:
                continue
            genes.append({
                "chrom":  p[0],
                "start":  int(p[1]),
                "end":    int(p[2]),
                "name":   p[3],
                "strand": p[5] if len(p) > 5 else "+",
                # Optional display label (BED col7 from extract_flanking_genes.py)
                "display_name": p[6] if len(p) > 6 else "",
            })
    return genes


def parse_candidate_regions(candidate_beds):
    """
    Parse candidate region BEDs grouped by genome ID inferred from filename.
    """
    regions_by_genome = defaultdict(list)
    for bed in candidate_beds or []:
        if not bed or not os.path.exists(bed):
            continue
        genome_id = clean_genome_name(
            os.path.basename(bed).replace(".regions.bed", "").replace(".bed", "")
        )
        with open(bed) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                p = line.split("\t")
                if len(p) < 3:
                    continue
                try:
                    chrom = p[0]
                    start = int(p[1])
                    end = int(p[2])
                except ValueError:
                    continue
                if end > start:
                    regions_by_genome[genome_id].append((chrom, start, end))
    return regions_by_genome


def _match_regions_for_genome(regions_by_genome, genome_id):
    """
    Resolve candidate regions for a target genome ID with tolerant matching.
    """
    if genome_id in regions_by_genome:
        return regions_by_genome[genome_id]
    for rid, regs in regions_by_genome.items():
        if rid in genome_id or genome_id in rid:
            return regs
    return []


def filter_genes_to_candidate_regions(genes, candidate_regions):
    """
    Keep only genes overlapping at least one candidate region.
    Candidate BED is 0-based half-open; parsed GFF genes are treated 1-based.
    """
    if not candidate_regions:
        return genes
    kept = []
    for g in genes:
        g_start0 = max(0, int(g["start"]) - 1)
        g_end0 = int(g["end"])
        for chrom, rs, re in candidate_regions:
            if g["chrom"] != chrom:
                continue
            ov = max(0, min(g_end0, re) - max(g_start0, rs))
            if ov > 0:
                kept.append(g)
                break
    return kept


def _is_goi_target_gene(gene):
    name = gene.get("name", "") or ""
    home_id = gene.get("home_gene_id", "") or ""
    return name.startswith("GOI_") or home_id.startswith("GOI_")


def _region_overlaps_gene(region, gene):
    chrom, rs, re = region
    if gene.get("chrom") != chrom:
        return False
    # Candidate BED is 0-based half-open; GFF genes are 1-based closed.
    gs0 = max(0, int(gene["start"]) - 1)
    ge0 = int(gene["end"])
    ov = max(0, min(ge0, re) - max(gs0, rs))
    return ov > 0


def _candidate_regions_with_goi(candidate_regions, genes):
    goi_genes = [g for g in genes if _is_goi_target_gene(g)]
    if not goi_genes:
        return []
    return [
        reg for reg in candidate_regions
        if any(_region_overlaps_gene(reg, gg) for gg in goi_genes)
    ]


def _select_goi_context_genes(genes, flank_bp=200000):
    """
    Fallback when candidate regions miss GOI loci:
    keep genes on the dominant GOI chromosome around GOI coordinates.
    """
    goi_genes = [g for g in genes if _is_goi_target_gene(g)]
    if not goi_genes:
        return []

    per_chrom = defaultdict(list)
    for g in goi_genes:
        per_chrom[g["chrom"]].append(g)

    def _chrom_key(item):
        chrom, glist = item
        best_identity = max((x.get("identity", 0.0) for x in glist), default=0.0)
        return (len(glist), best_identity)

    goi_chrom, goi_list = max(per_chrom.items(), key=_chrom_key)
    goi_min = min(g["start"] for g in goi_list)
    goi_max = max(g["end"] for g in goi_list)
    win_s = max(1, goi_min - max(0, int(flank_bp)))
    win_e = goi_max + max(0, int(flank_bp))

    selected = [
        g for g in genes
        if g["chrom"] == goi_chrom and g["end"] >= win_s and g["start"] <= win_e
    ]
    return selected


def _parse_gff_attrs(attr_field):
    attrs = {}
    for kv in (attr_field or "").split(";"):
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        attrs[k] = unquote(v)
    return attrs


def _is_generic_gene_label(name):
    """Return True for non-informative locus-tag style labels."""
    if not name:
        return True
    txt = clean_gene_label(str(name).strip())
    if not txt:
        return True
    if re.match(r"^[A-Za-z]{1,8}\d*_\d+$", txt):
        return True
    if re.match(r"^LOC\d+$", txt, re.IGNORECASE):
        return True
    return False


def _is_noninformative_product(product):
    """Return True for generic/placeholder product annotations."""
    if not product:
        return True
    txt = str(product).strip().lower()
    generic = (
        "hypothetical protein",
        "uncharacterized protein",
        "unknown protein",
        "predicted protein",
    )
    return any(tok in txt for tok in generic)


def _format_product_label(product, max_words=5, max_chars=42):
    """
    Convert long product descriptions into compact labels suitable for plotting.
    """
    if not product:
        return ""
    txt = str(product).strip()
    txt = re.sub(r"^(putative|probable|predicted)\s+", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r"\s+", " ", txt)
    words = txt.split(" ")
    if len(words) > max_words:
        txt = " ".join(words[:max_words])
    if len(txt) > max_chars:
        txt = txt[: max_chars - 3].rstrip() + "..."
    return txt


def _preferred_target_label(gene):
    """
    Prefer native target annotation labels when informative.
    Fallback order: target_gene -> target_product -> name -> home_gene_id.
    """
    target_gene = gene.get("target_gene", "")
    if target_gene and not _is_generic_gene_label(target_gene):
        return target_gene

    target_product = gene.get("target_product", "")
    if target_product and not _is_noninformative_product(target_product):
        pretty = _format_product_label(target_product)
        if pretty:
            return pretty

    name = gene.get("name", "")
    if name and not _is_generic_gene_label(name):
        return name

    return target_gene or name or gene.get("home_gene_id", "")


def parse_target_gff(gff_file):
    """
    Parse a SynTerra target-genome GFF.

    Extracts mRNA plus gene-level features for tandem copies.
    Returns list of gene dicts with 'home_gene_id' from SynTerra_Parent.
    Deduplicates overlapping entries (same region annotated by different queries).
    """
    genes = []
    if not gff_file or not os.path.exists(gff_file):
        return genes
    with open(gff_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            ftype = p[2]
            if ftype not in ("mRNA", "gene", "tandem_copy"):
                continue
            attrs = _parse_gff_attrs(p[8])
            raw_name = attrs.get("Name", attrs.get("ID", ""))
            target_gene = attrs.get("TargetGene", "")
            target_product = attrs.get("TargetProduct", "")
            target_id = attrs.get("TargetID", "")

            try:
                identity = float(attrs.get("Identity", "0"))
            except Exception:
                identity = 0.0

            genes.append({
                "chrom":        p[0],
                "start":        int(p[3]),
                "end":          int(p[4]),
                "name":         raw_name,
                "target_gene":  target_gene,
                "target_product": target_product,
                "target_id":    target_id,
                "strand":       p[6],
                "identity":     identity,
                "home_gene_id": attrs.get("SynTerra_Parent", attrs.get("Parent", "")),
                "n_exons":      int(attrs.get("Exons", "1")),
            })

    # Deduplicate overlapping entries (same genomic region from different queries)
    if len(genes) > 1:
        genes.sort(key=lambda g: -g["identity"])  # best identity first
        kept = []
        for g in genes:
            is_dup = False
            for k in kept:
                if g["chrom"] != k["chrom"]:
                    continue
                ov = max(0, min(g["end"], k["end"]) - max(g["start"], k["start"]))
                len_g = max(1, g["end"] - g["start"])
                len_k = max(1, k["end"] - k["start"])
                if min(ov / len_g, ov / len_k) >= 0.50:
                    is_dup = True
                    break
            if not is_dup:
                kept.append(g)
        genes = kept

    return genes


def parse_homology_tsvs(tsv_files):
    """Parse homology TSVs -> dict mapping target_gene -> home_gene."""
    mapping = {}
    if not tsv_files:
        return mapping
    for tsv in tsv_files:
        if not tsv or tsv == "NO_HOMOLOGY" or not os.path.exists(tsv):
            continue
        with open(tsv) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]
    return mapping


def parse_home_gff_products(gff_file):
    """Parse home GFF -> dict mapping gene ID/Name -> product description."""
    products = {}
    if not gff_file or not os.path.exists(gff_file) or gff_file == "NO_GFF":
        return products
    try:
        from urllib.parse import unquote
        with open(gff_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.strip().split("\t")
                if len(p) < 9:
                    continue
                attrs = {}
                for kv in p[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k] = unquote(v)
                product = attrs.get("product", "")
                if product:
                    for key in ("ID", "Name", "Parent", "gene", "locus_tag"):
                        if key in attrs:
                            products[attrs[key]] = product
    except Exception as exc:
        print(f"Warning: could not parse home GFF products: {exc}")
    return products


# ======================================================================
# Tree helpers
# ======================================================================

def _genome_id_from_leaf(leaf_name):
    """
    Extract a GCF/GCA genome accession from a tree leaf name.

    Leaf format examples:
      GOI_P01501|GCF_029169275_1_fna_exon_ann  ->  GCF_029169275.1
      GOI_P01501                                ->  None (home)
    """
    if "|" not in leaf_name:
        return None
    for part in leaf_name.split("|"):
        if part.startswith("GCF_") or part.startswith("GCA_"):
            # GCF_029169275_1_fna_exon_ann -> GCF_029169275.1
            pieces = part.replace("_fna_exon_ann", "").replace("_fna", "").split("_")
            if len(pieces) >= 3:
                return f"{pieces[0]}_{pieces[1]}.{pieces[2]}"
            return "_".join(pieces)
    return None


def parse_tree_clade_colours(tree_file):
    """
    Assign warm-palette colours to GOI leaves based on phylogenetic tree.

    Returns
    -------
    goi_genome_colours : dict   genome_id|'home' -> hex colour
    target_order       : list   genome_ids sorted by distance to home (closest first)
    """
    goi_colours = {}
    target_order = []

    if not tree_file or not os.path.exists(tree_file) or not ETE3_AVAILABLE:
        return goi_colours, target_order

    try:
        t = Tree(tree_file)
        leaves = list(t.iter_leaves())
        n = len(leaves)
        if n == 0:
            return goi_colours, target_order

        # Identify home-genome leaf (no genome ID in name)
        home_leaves   = [l for l in leaves if _genome_id_from_leaf(l.name) is None]
        target_leaves = [l for l in leaves if _genome_id_from_leaf(l.name) is not None]

        # Assign warm colours along tree-traversal order (red -> amber)
        for i, leaf in enumerate(leaves):
            hue = 0.0 + (i / max(1, n - 1)) * 0.20
            sat = 0.90 - (i / max(1, n - 1)) * 0.20
            r, g, b = colorsys.hsv_to_rgb(hue, sat, 0.90)
            colour = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
            gid = _genome_id_from_leaf(leaf.name)
            if gid:
                goi_colours[gid] = colour
            else:
                goi_colours["home"] = colour

        # Order targets by phylogenetic distance to home
        if home_leaves:
            ref = home_leaves[0]
            dist_map = {}
            for tl in target_leaves:
                gid = _genome_id_from_leaf(tl.name)
                if gid:
                    d = t.get_distance(ref, tl)
                    if gid not in dist_map or d < dist_map[gid]:
                        dist_map[gid] = d
            target_order = sorted(dist_map, key=dist_map.get)

        print(f"Tree: assigned {len(goi_colours)} GOI colours, "
              f"target order = {target_order}")
    except Exception as exc:
        print(f"Warning: could not parse tree: {exc}")
    return goi_colours, target_order


# ======================================================================
# Colour assignment
# ======================================================================

# Module-level set populated during main() with names of GOI genes
_GOI_NAMES = set()


def is_goi(name):
    """Return True if *name* represents the Gene of Interest."""
    if not name:
        return False
    if name.startswith("GOI_") or "|exon_" in name:
        return True
    return name in _GOI_NAMES


def _overlaps_any(gene, intervals):
    for q in intervals:
        if gene["chrom"] == q["chrom"] and gene["start"] < q["end"] and gene["end"] > q["start"]:
            return True
    return False


def identify_goi_names(home_genes, query_intervals):
    """
    Identify which home-gene names are GOI by overlapping with query_bed.
    Populates the module-level _GOI_NAMES set.
    Only marks genes that are *small* relative to the query span as GOI.
    Large container genes (e.g. LOC726866 spanning 17kb) are excluded.
    """
    _GOI_NAMES.clear()
    if not query_intervals:
        return
    # Query span
    q_span = sum(q["end"] - q["start"] for q in query_intervals)
    max_goi_size = max(q_span * 20, 5000)  # generous but bounded
    best_overlap = 0
    best_genes = []
    for gene in home_genes:
        gsize = gene["end"] - gene["start"]
        if gsize > max_goi_size:
            continue
        overlap = 0
        for q in query_intervals:
            if gene["chrom"] != q["chrom"]:
                continue
            ov = min(gene["end"], q["end"]) - max(gene["start"], q["start"])
            if ov > 0:
                overlap += ov
        if overlap <= 0:
            continue
        # TODO: Review GOI selection logic for multi-locus runs (audit if this is too strict)
        if overlap > best_overlap:
            best_overlap = overlap
            best_genes = [gene]
        elif overlap == best_overlap:
            best_genes.append(gene)
    for gene in best_genes:
        _GOI_NAMES.add(gene["name"])
    # Always include GOI_ prefixed names
    for gene in home_genes:
        if gene["name"].startswith("GOI_"):
            _GOI_NAMES.add(gene["name"])
    if _GOI_NAMES:
        print(f"GOI genes identified: {_GOI_NAMES}")


def assign_gene_colours(home_genes, query_intervals=None):
    """
    Map each home-gene name -> hex colour.

    GOI genes -> GOI_COLOUR (will be overridden per-genome with tree colours).
    Flanking genes -> GENE_PALETTE (deterministic order).
    Uses name-based GOI check only (not coordinate overlap) to avoid
    marking large container-loci as GOI.
    """
    cmap = {}
    idx = 0
    for gene in home_genes:
        name = gene["name"]
        if name in cmap:
            continue
        if is_goi(name):
            cmap[name] = GOI_COLOUR
        else:
            cmap[name] = GENE_PALETTE[idx % len(GENE_PALETTE)]
            idx += 1
    return cmap


# ======================================================================
# Genome / gene-name helpers
# ======================================================================

def clean_genome_name(name):
    """GCF_029169275.1.fna -> GCF_029169275.1"""
    name = os.path.basename(name)
    for sfx in (".fna", ".fa", ".fasta", ".gz"):
        if name.endswith(sfx):
            name = name[: -len(sfx)]
    return name


def clean_gene_label(name):
    """gene-LOC412898 -> LOC412898 ;  GOI_P01501 -> P01501"""
    if name is None:
        return ""
    name = str(name)
    if name.startswith("gene-"):
        return name[5:]
    if name.startswith("GOI_"):
        return name[4:]
    return name


# ======================================================================
# Drawing primitives
# ======================================================================


def _arrow_xy(x0, x1, y_base, height, strand):
    """Pentagon vertices for a gene arrow."""
    w = x1 - x0
    aw = min(w * 0.25, height * 2.5)       # arrow-head width (capped)
    if aw < 1:
        aw = min(w * 0.5, 1)
    ym = y_base + height / 2
    yt = y_base + height
    if strand == "+":
        xs = [x0, x1 - aw, x1, x1 - aw, x0, x0]
        ys = [y_base, y_base, ym, yt, yt, y_base]
    else:
        xs = [x0, x0 + aw, x1, x1, x0 + aw, x0]
        ys = [ym, y_base, y_base, yt, yt, ym]
    return xs, ys


def _hex_to_rgba(hexc, alpha):
    hexc = hexc.lstrip("#")
    r, g, b = int(hexc[0:2], 16), int(hexc[2:4], 16), int(hexc[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _get_coords(gene):
    return gene.get("start_plot", gene["start"]), gene.get("end_plot", gene["end"])


def add_gene(fig, gene, x_off, y_base, h, colour, border_clr, border_w,
             hover, show_legend, legend_group):
    g_start, g_end = _get_coords(gene)
    xs, ys = _arrow_xy(g_start - x_off, g_end - x_off, y_base, h, gene["strand"])
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        fill="toself",
        fillcolor=colour,
        line=dict(color=border_clr, width=border_w),
        mode="lines",
        hoverinfo="text",
        text=hover,
        showlegend=show_legend,
        legendgroup=legend_group,
        name=clean_gene_label(legend_group),
    ))


def add_ribbon(fig, g_upper, g_lower, off_u, off_l, y_u_bot, y_l_top, colour, alpha=0.18):
    u_start, u_end = _get_coords(g_upper)
    l_start, l_end = _get_coords(g_lower)
    
    u0 = u_start - off_u
    u1 = u_end   - off_u
    l0 = l_start - off_l
    l1 = l_end   - off_l
    fill = _hex_to_rgba(colour, alpha)
    edge = _hex_to_rgba(colour, alpha * 1.8)
    fig.add_trace(go.Scatter(
        x=[u0, u1, l1, l0, u0],
        y=[y_u_bot, y_u_bot, y_l_top, y_l_top, y_u_bot],
        fill="toself", fillcolor=fill,
        line=dict(color=edge, width=0.5),
        mode="lines", hoverinfo="skip", showlegend=False,
    ))


def add_label(fig, gene, x_off, y_base, h, text, fsize=8, fcolour="black",
              is_goi_flag=False):
    g_start, g_end = _get_coords(gene)
    xc = (g_start + g_end) / 2 - x_off
    gw = g_end - g_start
    if is_goi_flag:
        text = "* " + text
        fcolour = GOI_BORDER
        fsize = max(fsize, 10)
    fig.add_annotation(
        x=xc, y=y_base + h + h * 0.35,
        text=text, showarrow=False,
        font=dict(size=fsize, color=fcolour),
        textangle=-35 if gw < 5000 else 0,
        xanchor="center", yanchor="bottom",
    )



def draw_gap_break(fig, x_pos, y_pos, height, text):
    """Draw a visual break mark (zigzag) and text label for a compressed gap."""
    # Zigzag path
    w = 200  # visual width of break (in genomic coordinates, scaled by Layout)
    # Actually, x_pos is the center.
    # We draw two parallel lines with a slash? Or distinct "break" symbol?
    # Simple: Text annotation with a small vertical tick or " // "
    
    fig.add_annotation(
        x=x_pos, y=y_pos + height/2,
        text=f"<b>//</b><br>{text}",
        showarrow=False,
        font=dict(size=9, color="black"),
        xanchor="center", yanchor="middle",
        yshift=0
    )


def compress_track_coordinates(genes, threshold=50000, visual_gap=2000):
    """
    Compress large gaps between genes.
    
    Returns
    -------
    compressed_genes : list of dicts (with added 'start_plot', 'end_plot')
    breaks           : list of dicts (x, gap_size, text)
    """
    if not genes:
        return [], []
    
    # Sort by start position
    # Note: genes might overlap, but generally are sequential.
    sorted_genes = sorted(genes, key=lambda g: g["start"])
    compressed = []
    breaks = []
    
    current_shift = 0
    
    # Initialize first gene
    # We preserve the absolute coordinate of the first gene *minus 0 shift* initially
    # effectively keeping the first gene at its real coordinate relative to start of cluster?
    # Yes, but we will shift everything by anchor later anyway.
    
    for i, g in enumerate(sorted_genes):
        new_g = g.copy()
        
        # Check gap from previous gene
        if i > 0:
            prev = sorted_genes[i-1]
            # Use raw coordinates for gap calculation
            gap = g["start"] - prev["end"]
            
            if gap > threshold:
                remove = gap - visual_gap
                current_shift += remove
                
                # The visual center of the gap in PLOT coordinates
                # prev_end_plot = prev['end'] - (current_shift - remove) 
                #               = prev['end_plot']
                # gap_start_plot = prev_end_plot
                # gap_end_plot   = gap_start_plot + visual_gap
                # center = gap_start_plot + visual_gap/2
                
                prev_end_plot = prev["end"] - (current_shift - remove) 
                break_x = prev_end_plot + visual_gap / 2
                
                breaks.append({
                    "x": break_x,
                    "gap_size": gap,
                    "text": f"{gap/1e6:.2f} Mb" if gap >= 1e6 else f"{gap/1e3:.0f} kb"
                })

        new_g["start_plot"] = g["start"] - current_shift
        new_g["end_plot"]   = g["end"]   - current_shift
        compressed.append(new_g)
        
    return compressed, breaks


def get_anchor_center(genes):
    """
    Find the center coordinate (plot) of the GOI.
    If multiple GOIs, average them. If none, usage center of the range.
    """
    goi_centers = []
    for g in genes:
        # Check 'name' and 'home_gene_id'
        if is_goi(g.get("name")) or is_goi(g.get("home_gene_id")):
            goi_centers.append((g["start_plot"] + g["end_plot"]) / 2)
            
    if goi_centers:
        return sum(goi_centers) / len(goi_centers)
        
    # Fallback: center of the entire cluster
    if not genes:
        return 0
    start = min(g["start_plot"] for g in genes)
    end   = max(g["end_plot"]   for g in genes)
    return (start + end) / 2

def _goi_colour_for_genome(genome_id, goi_genome_colours):
    """Look up GOI colour for a specific genome, with fuzzy matching."""
    if not goi_genome_colours:
        return GOI_COLOUR
    # Exact match
    if genome_id in goi_genome_colours:
        return goi_genome_colours[genome_id]
    # Prefix match  (e.g. "GCF_029169275.1" in "GCF_029169275.1.fna")
    for key, clr in goi_genome_colours.items():
        if key in genome_id or genome_id in key:
            return clr
    # Home fallback
    return goi_genome_colours.get("home", GOI_COLOUR)


def _lookup_product(gene_name, products):
    """Fuzzy product-name lookup."""
    for candidate in (gene_name, gene_name.replace("gene-", ""),
                      "gene-" + gene_name if not gene_name.startswith("gene-") else ""):
        if candidate in products:
            return products[candidate]
    return ""


def _preferred_home_label(gene, home_products):
    """
    Prefer informative home labels:
    gene symbol/name first, then product for generic locus-tag IDs.
    """
    display_name = gene.get("display_name", "")
    cleaned_display = clean_gene_label(display_name)
    if cleaned_display and not _is_generic_gene_label(cleaned_display):
        return cleaned_display

    name = gene.get("name", "")
    cleaned = clean_gene_label(name)
    if cleaned and not _is_generic_gene_label(cleaned):
        return cleaned

    product = _lookup_product(name, home_products)
    if product and not _is_noninformative_product(product):
        pretty = _format_product_label(product)
        if pretty:
            return pretty

    return cleaned or name


# ======================================================================
# Tree visualization
# ======================================================================

def _render_tree_html(tree_file, goi_genome_colours, output_path,
                      species_map=None):
    """
    Render a horizontal dendrogram of the GOI phylogenetic tree as an
    interactive Plotly HTML file.  Leaf nodes are coloured with the same
    clade palette used in the synteny plot.
    """
    if not tree_file or not os.path.exists(tree_file) or not ETE3_AVAILABLE:
        return

    try:
        t = Tree(tree_file)
    except Exception as exc:
        print(f"Warning: could not parse tree for rendering: {exc}")
        return

    leaves = list(t.iter_leaves())
    if len(leaves) < 2:
        return

    # --- 1. Assign (x, y) coordinates via recursive DFS ----------------
    # x = branch length (horizontal), y = leaf index (vertical)
    node_coords = {}           # node -> (x, y)
    leaf_counter = [0]         # mutable counter

    def _layout(node, x_offset):
        if node.is_leaf():
            y = leaf_counter[0]
            leaf_counter[0] += 1
            node_coords[node] = (x_offset + node.dist, y)
        else:
            child_ys = []
            for child in node.children:
                _layout(child, x_offset + node.dist)
                child_ys.append(node_coords[child][1])
            node_coords[node] = (x_offset + node.dist, sum(child_ys) / len(child_ys))

    _layout(t, 0)

    # --- 2. Build Plotly traces ----------------------------------------
    fig = go.Figure()

    # Branch lines (parent -> child: horizontal then vertical)
    for node in t.traverse():
        if node.is_root():
            continue
        parent = node.up
        px, py = node_coords[parent]
        cx, cy = node_coords[node]
        # Horizontal line from parent x to child x, at child y
        fig.add_trace(go.Scatter(
            x=[px, px, cx], y=[py, cy, cy],
            mode="lines",
            line=dict(color="black", width=1.5),
            hoverinfo="skip", showlegend=False,
        ))

    # Leaf dots + labels
    for leaf in leaves:
        lx, ly = node_coords[leaf]
        gid = _genome_id_from_leaf(leaf.name)
        key = gid if gid else "home"

        # Colour: try exact, then fuzzy
        colour = GOI_COLOUR
        if goi_genome_colours:
            if key in goi_genome_colours:
                colour = goi_genome_colours[key]
            else:
                for k, c in goi_genome_colours.items():
                    if k in key or key in k:
                        colour = c
                        break

        # Clean label
        label = leaf.name
        if "|" in label:
            parts = label.split("|")
            goi_part = parts[0]
            genome_part = parts[1] if len(parts) > 1 else ""
            # Prettify genome ID — use species name when available
            genome_pretty = genome_part.replace("_fna_exon_ann", "").replace("_fna", "")
            if species_map:
                for acc, sp_name in species_map.items():
                    if acc in genome_pretty:
                        genome_pretty = sp_name
                        break
            label = f"{goi_part} | {genome_pretty}"
        else:
            label = f"{label} (home)"

        fig.add_trace(go.Scatter(
            x=[lx], y=[ly],
            mode="markers+text",
            marker=dict(size=14, color=colour, line=dict(color="black", width=1)),
            text=[label],
            textposition="middle right",
            textfont=dict(size=11),
            hovertext=f"<b>{leaf.name}</b><br>Branch length: {leaf.dist:.6f}",
            hoverinfo="text",
            showlegend=False,
        ))

    # --- 3. Layout -----------------------------------------------------
    n_leaves = len(leaves)
    fig.update_layout(
        title=dict(
            text="<b>SynTerra GOI Phylogenetic Tree</b>",
            x=0.5, font=dict(size=15),
        ),
        height=max(300, n_leaves * 60 + 100),
        width=900,
        xaxis=dict(
            title="Evolutionary distance",
            showgrid=True, gridcolor="rgba(200,200,200,0.3)",
            zeroline=True,
        ),
        yaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.5, n_leaves - 0.5],
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=40, r=300, t=60, b=50),
    )

    fig.write_html(output_path)
    print(f"Tree plot saved to {output_path}")


# ======================================================================
# Main
# ======================================================================

def main():
    ap = argparse.ArgumentParser(description="SynTerra synteny plot")
    ap.add_argument("--home_bed",       required=True)
    ap.add_argument("--home_gff",       default=None)
    ap.add_argument("--query_bed",      default=None)
    ap.add_argument("--target_gffs",    nargs="*", default=[])
    ap.add_argument("--target_names",   nargs="*", default=[])
    ap.add_argument("--candidate_beds", nargs="*", default=[])
    ap.add_argument("--homology_tsvs",  nargs="*", default=[])
    ap.add_argument("--tree",           default=None)
    ap.add_argument("--species_map",    default=None,
                    help="TSV mapping accession → species name")
    ap.add_argument("--output",         required=True)
    args = ap.parse_args()

    # -- 0. Load species mapping -----------------------------------------
    species_map = {}  # accession -> species name
    if args.species_map and os.path.exists(args.species_map):
        with open(args.species_map) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    # Format can be:
                    # 2 columns: accession<TAB>species
                    # 3 columns: accession<TAB>species<TAB>tax_level
                    species_map[parts[0]] = parts[1]
        print(f"[plot] Loaded species mapping for {len(species_map)} genomes")

    # -- 1. Parse inputs -------------------------------------------------

    home_genes = parse_bed(args.home_bed)
    if not home_genes:
        msg = f"ERROR: empty home BED: {args.home_bed}"
        print(msg, file=sys.stderr)
        fig = go.Figure()
        fig.add_annotation(
            text=msg,
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False,
            font=dict(size=14, color="crimson"),
        )
        fig.update_layout(
            title="SynTerra Synteny Plot (Failed: empty home BED)",
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        fig.write_html(args.output)
        sys.exit(2)
    home_genes.sort(key=lambda g: g["start"])

    query_intervals = []
    if args.query_bed and os.path.exists(args.query_bed):
        for g in parse_bed(args.query_bed):
            query_intervals.append({"chrom": g["chrom"],
                                    "start": g["start"], "end": g["end"]})

    # Identify GOI gene names dynamically from query_bed overlap
    identify_goi_names(home_genes, query_intervals)

    home_products = parse_home_gff_products(args.home_gff) if args.home_gff else {}
    homology_map  = parse_homology_tsvs(args.homology_tsvs)

    goi_genome_colours, tree_target_order = parse_tree_clade_colours(args.tree)

    # -- 2. Build target tracks (matched by filename, not positional index)
    candidate_regions_by_genome = parse_candidate_regions(args.candidate_beds)

    target_tracks = []
    for gff_file in args.target_gffs:
        genome_id = clean_genome_name(
            os.path.basename(gff_file).replace(".gff", ""))
        genes_all = parse_target_gff(gff_file)
        candidate_regions = _match_regions_for_genome(candidate_regions_by_genome, genome_id)
        genes = filter_genes_to_candidate_regions(genes_all, candidate_regions)

        # Prefer candidate regions that actually overlap GOI models.
        goi_candidate_regions = _candidate_regions_with_goi(candidate_regions, genes_all)
        if goi_candidate_regions:
            genes = filter_genes_to_candidate_regions(genes_all, goi_candidate_regions)
            candidate_regions = goi_candidate_regions

        # If candidate regions exist but miss GOI, recover a GOI-centered context.
        if candidate_regions and not any(_is_goi_target_gene(g) for g in genes):
            fallback_genes = _select_goi_context_genes(genes_all, flank_bp=200000)
            if fallback_genes:
                genes = fallback_genes
                print(
                    f"[plot] {genome_id}: candidate regions missed GOI; "
                    f"using GOI-centered fallback ({len(genes)} genes)."
                )

        # Richest-block fallback: if candidate regions yielded very few genes
        # (≤2), they likely cover low-quality scoring windows rather than the
        # block with the best synteny evidence.  Find the densest genomic
        # cluster of genes from the full GFF and show that instead.
        if len(genes) <= 2 and len(genes_all) > len(genes):
            from collections import Counter
            block_counter = Counter()
            block_genes = defaultdict(list)
            for g in genes_all:
                # Group by block ID extracted from gene name (e.g. "...b15_fl1...")
                gname = g.get("name", "")
                import re as _re
                m = _re.search(r'_b(\d+)_', gname)
                if m:
                    bid = m.group(1)
                    block_counter[bid] += 1
                    block_genes[bid].append(g)
            if block_counter:
                best_block = block_counter.most_common(1)[0][0]
                best_genes = block_genes[best_block]
                if len(best_genes) > len(genes):
                    genes = best_genes
                    print(
                        f"[plot] {genome_id}: candidate regions had ≤2 genes; "
                        f"using richest block b{best_block} ({len(genes)} genes)."
                    )

        # Reduce clutter: if GOI is present, keep the GOI chromosome only.
        if any(_is_goi_target_gene(g) for g in genes):
            goi_chroms = {g["chrom"] for g in genes if _is_goi_target_gene(g)}
            if len(goi_chroms) == 1:
                goi_chrom = next(iter(goi_chroms))
                genes = [g for g in genes if g["chrom"] == goi_chrom]

        # Don't skip target track if there are no genes found; we want to show it's empty
        genes.sort(key=lambda g: g["start"])
        # Use species name from mapping if available
        display = genome_id
        for acc, sp_name in species_map.items():
            if acc in genome_id:
                display = sp_name
                break
        target_chrom = genes[0]["chrom"] if genes else (candidate_regions[0][0] if candidate_regions else "unknown")
        
        target_tracks.append({
            "genome_id":    genome_id,
            "display_name": display,
            "genes":        genes,
            "chrom":        target_chrom,
        })

    # Order targets by phylogenetic distance (if tree available), else alphabetically by species name
    if tree_target_order:
        def _tree_key(t):
            for i, gid in enumerate(tree_target_order):
                if gid in t["genome_id"]:
                    return i
            return 999
        target_tracks.sort(key=_tree_key)
    else:
        # Fallback: order by display name (species name) alphabetically
        target_tracks.sort(key=lambda t: t["display_name"])

    # -- 3. Colour map ---------------------------------------------------

    gene_colours = assign_gene_colours(home_genes, query_intervals)

    # -- 4. Assemble track list & Compress -------------------------------

    home_chrom = home_genes[0]["chrom"]
    
    # Initial list with raw genes
    raw_tracks = [{
        "label":     f"Home genome ({home_chrom})",
        "genes":     home_genes,
        "is_home":   True,
        "genome_id": "home",
        "chrom":     home_chrom
    }]
    for tt in target_tracks:
        raw_tracks.append({
            "label":     f"{tt['display_name']} ({tt['chrom']})",
            "genes":     tt["genes"],
            "is_home":   False,
            "genome_id": tt["genome_id"],
            "chrom":     tt["chrom"]
        })

    all_tracks = []
    for track in raw_tracks:
        # 1. Compress
        c_genes, breaks = compress_track_coordinates(track["genes"], threshold=50000, visual_gap=3000)
        track["genes"]  = c_genes
        track["breaks"] = breaks
        
        # 2. Find Anchor (GOI center) to align at x=0
        anchor = get_anchor_center(c_genes)
        track["offset"] = anchor  # This effectively centers the plot on the GOI
        
        all_tracks.append(track)

    n_tracks = len(all_tracks)

    # -- 5. Layout geometry -----------------------------------------------

    GENE_H      = 0.35          # gene arrow height
    TRACK_SPACE = 1.4           # vertical pitch between tracks
    RIBBON_GAP  = 0.12          # gap between gene arrow and ribbon edge

    fig = go.Figure()

    # -- 5a. Track background bands --------------------------------------
    for ti, track in enumerate(all_tracks):
        yb    = (n_tracks - 1 - ti) * TRACK_SPACE
        x_off = track["offset"]
        if not track["genes"]:
            continue
        # Use plotted coordinates for background extent
        x_min = min(g["start_plot"] for g in track["genes"]) - x_off - 1000
        x_max = max(g["end_plot"]   for g in track["genes"]) - x_off + 1000
        fig.add_shape(
            type="rect",
            x0=x_min, x1=x_max, y0=yb - 0.02, y1=yb + GENE_H + 0.02,
            fillcolor=TRACK_BG_CLR, line=dict(width=0), layer="below",
        )

    # -- 5b. Ribbons (draw first so they sit behind genes) ---------------
    for ti in range(n_tracks - 1):
        upper = all_tracks[ti]
        lower = all_tracks[ti + 1]
        y_u = (n_tracks - 1 - ti)       * TRACK_SPACE
        y_l = (n_tracks - 1 - (ti + 1)) * TRACK_SPACE
        y_ribbon_top = y_u - RIBBON_GAP
        y_ribbon_bot = y_l + GENE_H + RIBBON_GAP

        for lg in lower["genes"]:
            home_id = lg.get("home_gene_id", "")
            if not home_id:
                continue
            for ug in upper["genes"]:
                u_name = ug["name"]
                u_home = ug.get("home_gene_id", u_name)
                match = (u_home == home_id or u_name == home_id
                         or (is_goi(u_home) and is_goi(home_id))
                         or (is_goi(u_name) and is_goi(home_id)))
                if match:
                    # Determine ribbon colour
                    colour = gene_colours.get(home_id,
                             gene_colours.get(u_name, UNMATCHED_CLR))
                    if is_goi(home_id):
                        colour = _goi_colour_for_genome(
                            lower["genome_id"], goi_genome_colours)
                    add_ribbon(fig, ug, lg,
                               upper["offset"], lower["offset"],
                               y_ribbon_top, y_ribbon_bot,
                               colour, alpha=0.20)
                    break

    # -- 5c. Gene arrows -------------------------------------------------
    legend_shown = set()

    for ti, track in enumerate(all_tracks):
        yb    = (n_tracks - 1 - ti) * TRACK_SPACE
        x_off = track["offset"]

        # Draw large genes first so small genes render on top
        # Use plot coordinates for size sorting? Yes.
        sorted_genes = sorted(track["genes"],
                               key=lambda g: g["end_plot"] - g["start_plot"],
                               reverse=True)

        for gene in sorted_genes:
            name = gene["name"]
            home_id = gene.get("home_gene_id", name)
            target_label = _preferred_target_label(gene)
            goi_f   = is_goi(name) or is_goi(home_id)

            # --- colour ---
            if goi_f:
                colour = _goi_colour_for_genome(
                    track["genome_id"], goi_genome_colours)
                bclr, bw = GOI_BORDER, 2.5
            elif home_id in gene_colours:
                colour = gene_colours[home_id]
                bclr, bw = "rgba(0,0,0,0.35)", 1
            elif name in gene_colours:
                colour = gene_colours[name]
                bclr, bw = "rgba(0,0,0,0.35)", 1
            else:
                colour = UNMATCHED_CLR
                bclr, bw = "rgba(0,0,0,0.15)", 0.5

            # --- hover text ---
            if track["is_home"]:
                cn = _preferred_home_label(gene, home_products)
            else:
                cn = clean_gene_label(target_label)
            if track["is_home"]:
                product = _lookup_product(name, home_products)
                hover = f"<b>{cn}</b>"
                if product:
                    hover += f"<br><i>{product}</i>"
                # Use raw coords for hover
                hover += (f"<br>{gene['chrom']}:{gene['start']:,}-{gene['end']:,}"
                          f"<br>Strand: {gene['strand']}")
                if goi_f:
                    hover += "<br><b>GENE OF INTEREST</b>"
            else:
                hover = f"<b>{cn}</b>"
                if gene.get("target_product"):
                    hover += f"<br><i>{gene['target_product']}</i>"
                if home_id:
                    hover += f"<br>Homolog: {clean_gene_label(home_id)}"
                if "identity" in gene:
                    hover += f"<br>Identity: {gene['identity']:.1f}%"
                if "n_exons" in gene:
                    hover += f"<br>Exons: {gene['n_exons']}"
                hover += (f"<br>{gene['chrom']}:{gene['start']:,}-{gene['end']:,}"
                          f"<br>Strand: {gene['strand']}")
                if goi_f:
                    hover += "<br><b>GENE OF INTEREST</b>"

            # --- legend (one entry per home-gene name) ---
            lg_key = home_id if home_id else name
            show_leg = lg_key not in legend_shown
            if show_leg:
                legend_shown.add(lg_key)

            add_gene(fig, gene, x_off, yb, GENE_H, colour, bclr, bw,
                     hover, show_leg, lg_key)

    # -- 5d. Gene labels -------------------------------------------------
    for ti, track in enumerate(all_tracks):
        yb    = (n_tracks - 1 - ti) * TRACK_SPACE
        x_off = track["offset"]
        for gene in track["genes"]:
            name    = gene["name"]
            home_id = gene.get("home_gene_id", name)
            goi_f   = is_goi(name) or is_goi(home_id)
            if not track["is_home"]:
                label = clean_gene_label(_preferred_target_label(gene))
                if not label and home_id:
                    label = clean_gene_label(home_id)
            else:
                label = _preferred_home_label(gene, home_products)
            add_label(fig, gene, x_off, yb, GENE_H, label,
                      fsize=8, is_goi_flag=goi_f)

    # -- 5e. Track labels (left margin) ----------------------------------
    for ti, track in enumerate(all_tracks):
        yb = (n_tracks - 1 - ti) * TRACK_SPACE
        fig.add_annotation(
            x=-0.01, y=yb + GENE_H / 2,
            text=f"<b>{track['label']}</b>",
            showarrow=False,
            font=dict(size=11, color="black"),
            xref="paper", yref="y",
            xanchor="right", yanchor="middle",
        )

    # -- 5f. Gap breaks --------------------------------------------------
    for ti, track in enumerate(all_tracks):
        yb    = (n_tracks - 1 - ti) * TRACK_SPACE
        x_off = track["offset"]
        for brk in track.get("breaks", []):
            draw_gap_break(fig, brk["x"] - x_off, yb, GENE_H, brk["text"])

    # -- 6. Figure styling -----------------------------------------------

    # Compute plotted range
    # Collect all plotted X coordinates
    all_x = []
    for track in all_tracks:
        x_off = track["offset"]
        for g in track["genes"]:
            all_x.append(g["start_plot"] - x_off)
            all_x.append(g["end_plot"]   - x_off)

    if not all_x:
        x_min, x_max = -1000, 1000
    else:
        x_min, x_max = min(all_x), max(all_x)
        
    pad = (x_max - x_min) * 0.05 + 5000
    x_range = [x_min - pad, x_max + pad]

    fig_height = max(500, n_tracks * 200 + 80)

    fig.update_layout(
        title=dict(
            text=("<b>SynTerra Synteny Plot</b>"
                  "<br><sup>Genes coloured by homology group | "
                  "* = Gene of Interest | Ribbons connect orthologs | // = Compressed Gaps</sup>"),
            x=0.5, font=dict(size=15),
        ),
        height=fig_height,
        width=1500,
        xaxis=dict(
            title="", # No title since numbers are relative/discontinuous
            showgrid=False, 
            zeroline=False,
            showticklabels=False, # Hide ticks as they are discontinuous
            range=x_range,
        ),
        yaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.6,
                   (n_tracks - 1) * TRACK_SPACE + GENE_H + 1.0],
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=260, r=60, t=85, b=55), # Increased left margin
        legend=dict(
            title="<b>Gene (home ID)</b>",
            orientation="v", x=1.01, y=1.0,
            font=dict(size=9),
            tracegroupgap=2,
        ),
        hovermode="closest",
    )
    
    # Add Scale Bar (10 kb)
    # Place it in bottom right? Or top left?
    # Let's put it at bottom right
    scale_len = 10000
    sb_x1 = x_max
    sb_x0 = x_max - scale_len
    sb_y  = -0.4
    
    fig.add_shape(
        type="line",
        x0=sb_x0, x1=sb_x1, y0=sb_y, y1=sb_y,
        line=dict(color="black", width=3),
    )
    fig.add_annotation(
        x=(sb_x0 + sb_x1)/2, y=sb_y - 0.1,
        text="<b>10 kb</b>",
        showarrow=False,
        font=dict(size=10, color="black"),
        yanchor="top"
    )

    fig.write_html(args.output)
    print(f"Synteny plot saved to {args.output}")
    print(f"  Tracks: {n_tracks} (1 home + {len(target_tracks)} targets)")
    print(f"  Gap compression: active (>50kb -> 3kb visual)")

    # -- 7. Tree plot (separate HTML) ------------------------------------
    tree_output = args.output.replace("_synteny_plot.html", "_tree.html")
    if tree_output == args.output:
        tree_output = args.output.replace(".html", "_tree.html")
    _render_tree_html(args.tree, goi_genome_colours, tree_output,
                      species_map=species_map)


if __name__ == "__main__":
    main()
