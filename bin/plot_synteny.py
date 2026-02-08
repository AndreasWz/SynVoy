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
  --candidate_beds  Cluster-region BED files (optional, not currently drawn)
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
import sys

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
            })
    return genes


def parse_target_gff(gff_file):
    """
    Parse a SynTerra target-genome GFF.

    Extracts **mRNA** features only (avoids CDS double-counting).
    Returns list of gene dicts with 'home_gene_id' from SynTerra_Parent.
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
            if len(p) < 9 or p[2] != "mRNA":
                continue
            attrs = {}
            for kv in p[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            genes.append({
                "chrom":        p[0],
                "start":        int(p[3]),
                "end":          int(p[4]),
                "name":         attrs.get("Name", attrs.get("ID", "")),
                "strand":       p[6],
                "identity":     float(attrs.get("Identity", "0")),
                "home_gene_id": attrs.get("SynTerra_Parent", ""),
                "n_exons":      int(attrs.get("Exons", "1")),
            })
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
                    for key in ("ID", "Name", "Parent"):
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

def is_goi(name):
    """Return True if *name* represents the Gene of Interest."""
    if not name:
        return False
    return name.startswith("GOI_") or name == "gene-Melt" or "|exon_" in name


def _overlaps_any(gene, intervals):
    for q in intervals:
        if gene["chrom"] == q["chrom"] and gene["start"] < q["end"] and gene["end"] > q["start"]:
            return True
    return False


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


def add_gene(fig, gene, x_off, y_base, h, colour, border_clr, border_w,
             hover, show_legend, legend_group):
    xs, ys = _arrow_xy(gene["start"] - x_off, gene["end"] - x_off, y_base, h, gene["strand"])
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
    u0 = g_upper["start"] - off_u
    u1 = g_upper["end"]   - off_u
    l0 = g_lower["start"] - off_l
    l1 = g_lower["end"]   - off_l
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
    xc = (gene["start"] + gene["end"]) / 2 - x_off
    gw = gene["end"] - gene["start"]
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


# ======================================================================
# Internal helpers
# ======================================================================

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
    ap.add_argument("--output",         required=True)
    args = ap.parse_args()

    # -- 1. Parse inputs -------------------------------------------------

    home_genes = parse_bed(args.home_bed)
    if not home_genes:
        print("ERROR: empty home BED", file=sys.stderr)
        go.Figure().write_html(args.output)
        return
    home_genes.sort(key=lambda g: g["start"])

    query_intervals = []
    if args.query_bed and os.path.exists(args.query_bed):
        for g in parse_bed(args.query_bed):
            query_intervals.append({"chrom": g["chrom"],
                                    "start": g["start"], "end": g["end"]})

    home_products = parse_home_gff_products(args.home_gff) if args.home_gff else {}
    homology_map  = parse_homology_tsvs(args.homology_tsvs)

    goi_genome_colours, tree_target_order = parse_tree_clade_colours(args.tree)

    # -- 2. Build target tracks (matched by filename, not positional index)

    target_tracks = []
    for gff_file in args.target_gffs:
        genome_id = clean_genome_name(
            os.path.basename(gff_file).replace(".gff", ""))
        genes = parse_target_gff(gff_file)
        if not genes:
            continue
        genes.sort(key=lambda g: g["start"])
        target_tracks.append({
            "genome_id":    genome_id,
            "display_name": genome_id,
            "genes":        genes,
            "chrom":        genes[0]["chrom"],
        })

    # Order targets by phylogenetic distance (if tree available)
    if tree_target_order:
        def _tree_key(t):
            for i, gid in enumerate(tree_target_order):
                if gid in t["genome_id"]:
                    return i
            return 999
        target_tracks.sort(key=_tree_key)

    # -- 3. Colour map ---------------------------------------------------

    gene_colours = assign_gene_colours(home_genes, query_intervals)

    # -- 4. Assemble track list ------------------------------------------

    home_chrom  = home_genes[0]["chrom"]
    home_offset = min(g["start"] for g in home_genes) - 2000

    all_tracks = [{
        "label":     f"Home genome ({home_chrom})",
        "genes":     home_genes,
        "offset":    home_offset,
        "is_home":   True,
        "genome_id": "home",
    }]
    for tt in target_tracks:
        off = min(g["start"] for g in tt["genes"]) - 2000
        all_tracks.append({
            "label":     f"{tt['display_name']} ({tt['chrom']})",
            "genes":     tt["genes"],
            "offset":    off,
            "is_home":   False,
            "genome_id": tt["genome_id"],
        })

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
        x_min = min(g["start"] for g in track["genes"]) - x_off - 1000
        x_max = max(g["end"]   for g in track["genes"]) - x_off + 1000
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
        sorted_genes = sorted(track["genes"],
                               key=lambda g: g["end"] - g["start"],
                               reverse=True)

        for gene in sorted_genes:
            name    = gene["name"]
            home_id = gene.get("home_gene_id", name)
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
            cn = clean_gene_label(name)
            if track["is_home"]:
                product = _lookup_product(name, home_products)
                hover = f"<b>{cn}</b>"
                if product:
                    hover += f"<br><i>{product}</i>"
                hover += (f"<br>{gene['chrom']}:{gene['start']:,}-{gene['end']:,}"
                          f"<br>Strand: {gene['strand']}")
                if goi_f:
                    hover += "<br><b>GENE OF INTEREST</b>"
            else:
                hover = f"<b>{cn}</b>"
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
            if not track["is_home"] and home_id:
                label = clean_gene_label(home_id)
            else:
                label = clean_gene_label(name)
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

    # -- 6. Figure styling -----------------------------------------------

    # Compute sensible x range across all tracks
    all_x_max = max(
        g["end"] - track["offset"]
        for track in all_tracks for g in track["genes"]
    )
    fig_height = max(500, n_tracks * 200 + 80)

    fig.update_layout(
        title=dict(
            text=("<b>SynTerra Synteny Plot</b>"
                  "<br><sup>Genes coloured by homology group | "
                  "* = Gene of Interest | Ribbons connect orthologs</sup>"),
            x=0.5, font=dict(size=15),
        ),
        height=fig_height,
        width=1500,
        xaxis=dict(
            title="Position (bp, relative to region start)",
            showgrid=True, gridcolor="rgba(200,200,200,0.3)",
            zeroline=False,
            range=[-all_x_max * 0.01, all_x_max * 1.04],
        ),
        yaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.6,
                   (n_tracks - 1) * TRACK_SPACE + GENE_H + 1.0],
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=220, r=60, t=85, b=55),
        legend=dict(
            title="<b>Gene (home ID)</b>",
            orientation="v", x=1.01, y=1.0,
            font=dict(size=9),
            tracegroupgap=2,
        ),
        hovermode="closest",
    )

    fig.write_html(args.output)
    print(f"Synteny plot saved to {args.output}")
    print(f"  Tracks: {n_tracks} (1 home + {len(target_tracks)} targets)")
    print(f"  Home genes: {len(home_genes)}")
    for tt in target_tracks:
        print(f"  {tt['display_name']}: {len(tt['genes'])} genes")


if __name__ == "__main__":
    main()
