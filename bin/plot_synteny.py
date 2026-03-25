#!/usr/bin/env python3
"""
plot_synteny.py  –  Interactive synteny visualization for SynVoy

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
  --target_gffs     Target-genome GFFs (SynVoy exon_annotation format)
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


# ---- Lightweight Newick parser (fallback when ete3 is unavailable) ----

class _SimpleNode:
    """Minimal tree node for Newick parsing when ete3 is broken."""

    __slots__ = ("name", "dist", "children", "up")

    def __init__(self, name="", dist=0.0):
        self.name = name
        self.dist = dist
        self.children = []
        self.up = None

    def is_leaf(self):
        return len(self.children) == 0

    def is_root(self):
        return self.up is None

    def iter_leaves(self):
        if self.is_leaf():
            yield self
        else:
            for child in self.children:
                yield from child.iter_leaves()

    def traverse(self):
        yield self
        for child in self.children:
            yield from child.traverse()

    def get_distance(self, other):
        """Compute patristic distance via LCA (simple BFS approach)."""
        # Build path-to-root for both nodes
        def _path_to_root(node):
            path = {}
            d = 0.0
            n = node
            while n is not None:
                path[id(n)] = d
                d += n.dist
                n = n.up
            return path

        path_self = _path_to_root(self)
        n = other
        d_other = 0.0
        while n is not None:
            if id(n) in path_self:
                return path_self[id(n)] + d_other
            d_other += n.dist
            n = n.up
        return float("inf")


def _parse_newick(newick_str):
    """Parse a Newick string into a _SimpleNode tree."""
    s = newick_str.strip().rstrip(";").strip()
    if not s:
        return _SimpleNode()

    pos = [0]

    def _parse():
        node = _SimpleNode()
        if s[pos[0]] == "(":
            pos[0] += 1  # skip '('
            while True:
                child = _parse()
                child.up = node
                node.children.append(child)
                if pos[0] < len(s) and s[pos[0]] == ",":
                    pos[0] += 1
                elif pos[0] < len(s) and s[pos[0]] == ")":
                    pos[0] += 1
                    break
                else:
                    break

        # Read label and/or distance
        label_chars = []
        while pos[0] < len(s) and s[pos[0]] not in (",", ")", ";", "("):
            label_chars.append(s[pos[0]])
            pos[0] += 1
        label = "".join(label_chars).strip()
        if ":" in label:
            parts = label.rsplit(":", 1)
            node.name = parts[0].strip()
            try:
                node.dist = float(parts[1])
            except ValueError:
                node.dist = 0.0
        else:
            node.name = label
        return node

    return _parse()


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


def _confidence_rank(value):
    return {"LOW": 0, "MEDIUM": 1, "HIGH": 2}.get((value or "").upper(), -1)


def _is_goi_target_gene(gene):
    role = (gene.get("role") or "").strip().lower()
    if role:
        return role == "goi"
    name = gene.get("name", "") or ""
    home_id = gene.get("home_gene_id", "") or ""
    return name.startswith("GOI_") or home_id.startswith("GOI_")


def _is_resolved_goi_target_gene(gene):
    if not _is_goi_target_gene(gene):
        return False
    goi_class = (gene.get("goi_class") or "").strip().lower()
    confidence = (gene.get("confidence") or "").strip().upper()
    if goi_class in {"ambiguous_goi_family_member", "tandem_goi_copy"}:
        return False
    return confidence != "LOW"


def _track_goi_status(track):
    if any(_is_resolved_goi_target_gene(g) for g in track.get("genes", [])):
        return "resolved"
    if any(_is_goi_target_gene(g) for g in track.get("genes", [])):
        return "ambiguous"
    return "absent"


def _format_bp_label(length_bp):
    value = max(0, int(length_bp or 0))
    if value >= 1_000_000:
        if value % 1_000_000 == 0:
            return f"{value // 1_000_000} Mb"
        return f"{value / 1_000_000:.1f} Mb"
    if value >= 1_000:
        if value % 1_000 == 0:
            return f"{value // 1_000} kb"
        return f"{value / 1_000:.1f} kb"
    return f"{value} bp"


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


def _goi_priority_key(gene):
    goi_like = 1 if _is_goi_target_gene(gene) else 0
    resolved = 1 if _is_resolved_goi_target_gene(gene) else 0
    goi_class = (gene.get("goi_class") or "").strip().lower()
    class_rank = {
        "confident_goi": 3,
        "probable_goi": 2,
        "tandem_goi_copy": 1,
        "ambiguous_goi_family_member": 0,
    }.get(goi_class, 1 if goi_like else -1)
    return (
        goi_like,
        resolved,
        _confidence_rank(gene.get("confidence")),
        class_rank,
        float(gene.get("identity", 0.0)),
    )


def parse_target_gff(gff_file):
    """
    Parse a SynVoy target-genome GFF.

    Extracts mRNA plus gene-level features for tandem copies.
    Also collects CDS sub-features to build exon coordinate lists.
    Returns list of gene dicts with 'home_gene_id' from SynVoy_Parent.
    Deduplicates overlapping entries (same region annotated by different queries).
    """
    genes = []
    cds_by_parent = defaultdict(list)  # mRNA_ID -> [(start, end), ...]
    if not gff_file or not os.path.exists(gff_file):
        return genes
    # First pass: collect CDS sub-features
    with open(gff_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            if p[2] == "CDS":
                attrs = _parse_gff_attrs(p[8])
                parent = attrs.get("Parent", "")
                if parent:
                    try:
                        cds_by_parent[parent].append((int(p[3]), int(p[4])))
                    except ValueError:
                        pass
    # Deduplicate and sort CDS intervals per parent
    for parent in cds_by_parent:
        coords = sorted(set(cds_by_parent[parent]))
        # Merge overlapping CDS intervals
        merged = [coords[0]]
        for s, e in coords[1:]:
            if s <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        cds_by_parent[parent] = merged
    # Second pass: collect gene/mRNA features
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

            try:
                query_coverage = float(attrs.get("QueryCoverage", "nan"))
            except Exception:
                query_coverage = None
            if query_coverage is not None and query_coverage != query_coverage:
                query_coverage = None

            role = (attrs.get("SynVoyRole") or "").strip().lower()
            if not role:
                role = "goi" if raw_name.startswith("GOI_") or attrs.get("SynVoy_Parent", "").startswith("GOI_") else "flanking"

            gene_id = attrs.get("ID", "")
            exon_coords = cds_by_parent.get(gene_id, [])
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
                "home_gene_id": attrs.get("SynVoy_Parent", attrs.get("Parent", "")),
                "n_exons":      int(attrs.get("Exons", "1")),
                "exon_coords":  exon_coords,
                "role":         role,
                "confidence":   (attrs.get("Confidence", "") or "").upper(),
                "goi_class":    attrs.get("GOIClass", ""),
                "evidence_type": attrs.get("EvidenceType", attrs.get("Type", "")),
                "synteny_context": attrs.get("SyntenyContext", ""),
                "query_coverage": query_coverage,
                "inference_reason": attrs.get("InferenceReason", ""),
            })

    # Deduplicate overlapping entries (same genomic region from different queries)
    # GOI entries are ALWAYS preferred over non-GOI entries at the same locus.
    if len(genes) > 1:
        # Sort: GOI first (always kept), then by descending identity
        def _dedup_sort_key(g):
            return (
                -_goi_priority_key(g)[0],
                -_goi_priority_key(g)[1],
                -_goi_priority_key(g)[2],
                -_goi_priority_key(g)[3],
                -g["identity"],
            )
        genes.sort(key=_dedup_sort_key)

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

    # Cap GOI entries per genome: keep only the N best by identity.
    # Iterative search can produce hundreds of low-quality fallback GOI
    # annotations (especially without target GFFs) scattered across many
    # chromosomes.  Keeping all of them clutters the plot with noisy
    # connections.  Retain the best MAX_GOI_PER_GENOME entries.
    MAX_GOI_PER_GENOME = 10
    goi_genes = [g for g in genes if _is_goi_target_gene(g)]
    if len(goi_genes) > MAX_GOI_PER_GENOME:
        goi_genes.sort(key=_goi_priority_key, reverse=True)
        goi_to_drop = set(id(g) for g in goi_genes[MAX_GOI_PER_GENOME:])
        genes = [g for g in genes if id(g) not in goi_to_drop]
        print(
            f"[plot] GOI cap: kept {MAX_GOI_PER_GENOME}/{len(goi_genes)} GOI entries "
            f"(dropped {len(goi_to_drop)} lower-priority GOI-like entries)"
        )

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
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[0] != "target_id":
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


def parse_home_gff_exons(gff_file, gene_names):
    """Parse home GFF to extract exon/CDS coordinates for genes in the plot.

    The NCBI GFF hierarchy is: gene -> mRNA -> exon/CDS.
    We build a mapping: gene_name -> [(start, end), ...] merged CDS intervals.

    Parameters
    ----------
    gff_file : str
        Path to the home genome GFF.
    gene_names : set
        Gene names (e.g. 'gene-LOC412108') present in the home BED.

    Returns
    -------
    dict : gene_name -> list of (start, end) tuples (sorted, merged CDS coords)
    """
    exons_by_gene = {}
    if not gff_file or not os.path.exists(gff_file) or gff_file == "NO_GFF":
        return exons_by_gene
    try:
        from urllib.parse import unquote

        # Step 1: map gene-ID -> gene-name, mRNA-ID -> gene-name
        gene_id_to_name = {}  # gene ID -> gene name from BED
        mrna_to_gene = {}     # mRNA ID -> gene name

        with open(gff_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.strip().split("\t")
                if len(p) < 9:
                    continue
                ftype = p[2]
                attrs = {}
                for kv in p[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k] = unquote(v)
                if ftype == "gene":
                    gid = attrs.get("ID", "")
                    # Check if this gene appears in our plot
                    if gid in gene_names:
                        gene_id_to_name[gid] = gid
                elif ftype == "mRNA":
                    parent = attrs.get("Parent", "")
                    if parent in gene_id_to_name:
                        mrna_id = attrs.get("ID", "")
                        if mrna_id:
                            mrna_to_gene[mrna_id] = gene_id_to_name[parent]

        # Step 2: collect CDS coordinates keyed by gene name
        cds_by_gene = defaultdict(list)
        with open(gff_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.strip().split("\t")
                if len(p) < 9:
                    continue
                if p[2] not in ("CDS", "exon"):
                    continue
                attrs = {}
                for kv in p[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k] = unquote(v)
                parent = attrs.get("Parent", "")
                gene_name = mrna_to_gene.get(parent)
                if gene_name:
                    try:
                        cds_by_gene[gene_name].append((int(p[3]), int(p[4])))
                    except ValueError:
                        pass

        # Merge overlapping intervals per gene
        for gene_name, coords in cds_by_gene.items():
            coords = sorted(set(coords))
            merged = [coords[0]]
            for s, e in coords[1:]:
                if s <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            exons_by_gene[gene_name] = merged
    except Exception as exc:
        print(f"Warning: could not parse home GFF exons: {exc}")
    return exons_by_gene


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

    if not tree_file or not os.path.exists(tree_file):
        return goi_colours, target_order

    try:
        if ETE3_AVAILABLE:
            t = Tree(tree_file)
            leaves = list(t.iter_leaves())
        else:
            with open(tree_file) as fh:
                newick_str = fh.read().strip()
            t = _parse_newick(newick_str)
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
                    if ETE3_AVAILABLE:
                        d = t.get_distance(ref, tl)
                    else:
                        d = ref.get_distance(tl)
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


def clean_gene_label(name, keep_goi_prefix=False):
    """gene-LOC412898 -> LOC412898 ;  GOI_P01501 -> P01501"""
    if name is None:
        return ""
    name = str(name)
    if name.startswith("gene-"):
        return name[5:]
    if name.startswith("GOI_"):
        if keep_goi_prefix:
            # For target tracks: translate 'GOI_copy_3' -> 'GOI #3'
            suffix = name[4:]
            import re as _re
            m = _re.match(r'copy_(\d+)', suffix)
            if m:
                return f"GOI #{m.group(1)}"
            m = _re.match(r'(.*?)_copy_(\d+)', suffix)
            if m:
                return f"GOI #{m.group(2)}"
            return f"GOI {suffix}" if suffix else "GOI"
        return name[4:]
    return name


# ======================================================================
# Drawing primitives
# ======================================================================


def _arrow_xy(x0, x1, y_base, height, strand, junctions=None):
    """Pentagon vertices for a gene arrow, optionally with V-notch
    indentations at exon-boundary *junctions* (list of x-positions)."""
    w = x1 - x0
    aw = min(w * 0.18, height * 3.5)       # arrow-head width (prominent)
    if aw < 1:
        aw = min(w * 0.5, 1)
    ym = y_base + height / 2
    yt = y_base + height
    yb = y_base

    # ---- notch geometry ------------------------------------------------
    if junctions:
        notch_depth = height * 0.35          # how far the V cuts into the body
        notch_hw = max(abs(w) * 0.018, 120)  # half-width of the V at the edge
        # keep only junctions that fit between the body edges (not in the
        # arrowhead or too close to the blunt end)
        if strand == "+":
            valid = sorted(j for j in junctions
                           if x0 + notch_hw * 1.5 < j < x1 - aw - notch_hw)
        else:
            valid = sorted(j for j in junctions
                           if x0 + aw + notch_hw < j < x1 - notch_hw * 1.5)
    else:
        valid = []

    if not valid:
        # simple pentagon (no notches)
        if strand == "+":
            xs = [x0, x1 - aw, x1, x1 - aw, x0, x0]
            ys = [yb, yb, ym, yt, yt, yb]
        else:
            xs = [x0, x0 + aw, x1, x1, x0 + aw, x0]
            ys = [ym, yb, yb, yt, yt, ym]
        return xs, ys

    # ---- build polygon with V-notches ----------------------------------
    xs, ys = [], []
    if strand == "+":
        # bottom edge  left → right  (notches point *up* into body)
        xs.append(x0); ys.append(yb)
        for jx in valid:
            xs.extend([jx - notch_hw, jx, jx + notch_hw])
            ys.extend([yb, yb + notch_depth, yb])
        xs.append(x1 - aw); ys.append(yb)
        # arrowhead
        xs.append(x1); ys.append(ym)
        # top edge  right → left  (notches point *down* into body)
        xs.append(x1 - aw); ys.append(yt)
        for jx in reversed(valid):
            xs.extend([jx + notch_hw, jx, jx - notch_hw])
            ys.extend([yt, yt - notch_depth, yt])
        xs.append(x0); ys.append(yt)
        xs.append(x0); ys.append(yb)          # close
    else:  # strand "-"
        # arrowhead tip
        xs.append(x0); ys.append(ym)
        # bottom edge  arrowhead → right
        xs.append(x0 + aw); ys.append(yb)
        for jx in valid:
            xs.extend([jx - notch_hw, jx, jx + notch_hw])
            ys.extend([yb, yb + notch_depth, yb])
        xs.append(x1); ys.append(yb)
        # right side
        xs.append(x1); ys.append(yt)
        # top edge  right → arrowhead
        for jx in reversed(valid):
            xs.extend([jx + notch_hw, jx, jx - notch_hw])
            ys.extend([yt, yt - notch_depth, yt])
        xs.append(x0 + aw); ys.append(yt)
        xs.append(x0); ys.append(ym)          # close
    return xs, ys


def _hex_to_rgba(hexc, alpha):
    hexc = hexc.lstrip("#")
    r, g, b = int(hexc[0:2], 16), int(hexc[2:4], 16), int(hexc[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _get_coords(gene):
    return gene.get("start_plot", gene["start"]), gene.get("end_plot", gene["end"])


def add_gene(fig, gene, x_off, y_base, h, colour, border_clr, border_w,
             hover, show_legend, legend_group, line_dash="solid"):
    g_start, g_end = _get_coords(gene)
    x0 = g_start - x_off
    x1 = g_end - x_off

    # Compute junction x-positions for V-notch polygon
    junction_xs = None
    exon_coords = gene.get("exon_coords", [])
    n_exons_attr = gene.get("n_exons", 0)
    if len(exon_coords) >= 2:
        # Real exon boundaries available
        gene_start_raw = gene["start"]
        gene_end_raw = gene["end"]
        gene_span_raw = max(1, gene_end_raw - gene_start_raw)
        junction_xs = []
        for i in range(len(exon_coords) - 1):
            junction_genomic = (exon_coords[i][1] + exon_coords[i + 1][0]) / 2.0
            frac = (junction_genomic - gene_start_raw) / gene_span_raw
            frac = max(0.02, min(0.98, frac))
            junction_xs.append(x0 + frac * (x1 - x0))
    elif n_exons_attr and n_exons_attr >= 2:
        # No real CDS boundaries; synthesize evenly-spaced notches from
        # the Exons=N attribute (common for flanking_hit_span genes).
        gene_w = x1 - x0
        junction_xs = []
        for k in range(1, n_exons_attr):
            frac = k / n_exons_attr
            junction_xs.append(x0 + frac * gene_w)

    xs, ys = _arrow_xy(x0, x1, y_base, h, gene["strand"], junctions=junction_xs)
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        fill="toself",
        fillcolor=colour,
        line=dict(color=border_clr, width=border_w, dash=line_dash),
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
    if is_goi_flag:
        text = "* " + text
        fcolour = GOI_BORDER
        fsize = max(fsize, 10)
    fig.add_annotation(
        x=xc, y=y_base + h + h * 0.25,
        text=text, showarrow=False,
        font=dict(size=fsize, color=fcolour),
        textangle=-35,
        xanchor="center", yanchor="bottom",
    )



def draw_gap_break(fig, x_pos, y_pos, height, text, is_chrom_break=False):
    """Draw a visual break mark and text label for a compressed gap or chromosome boundary."""
    if is_chrom_break:
        # Draw a prominent vertical dashed line for chromosome boundaries
        fig.add_shape(
            type="line",
            x0=x_pos, x1=x_pos,
            y0=y_pos - 0.08, y1=y_pos + height + 0.08,
            line=dict(color="rgba(80,80,80,0.55)", width=2, dash="dot"),
        )
    else:
        fig.add_annotation(
            x=x_pos, y=y_pos + height / 2,
            text=f"<b>//</b><br>{text}",
            showarrow=False,
            font=dict(size=9, color="black"),
            xanchor="center", yanchor="middle",
        )


def _draw_chrom_labels(fig, track, ti, track_y, x_off, gene_h):
    """Draw chromosome labels below each chromosome segment within a track."""
    genes = track["genes"]
    if not genes:
        return
    yb = track_y[ti]

    # Group genes by chromosome (preserving plot order)
    from collections import OrderedDict
    chrom_segments = OrderedDict()
    for g in genes:
        ch = g["chrom"]
        if ch not in chrom_segments:
            chrom_segments[ch] = []
        chrom_segments[ch].append(g)

    if len(chrom_segments) <= 1:
        return  # Only one chromosome — no need for labels

    for ch, ch_genes in chrom_segments.items():
        # Find the center of this chromosome segment in plot coords
        xs = [g["start_plot"] - x_off for g in ch_genes] + \
             [g["end_plot"] - x_off for g in ch_genes]
        cx = (min(xs) + max(xs)) / 2
        # Shorten chromosome name for display
        short = ch
        if len(ch) > 12:
            short = ch[-10:]  # e.g. NC_045757.1 → _045757.1
        fig.add_annotation(
            x=cx, y=yb - 0.18,
            text=f"<b>{short}</b>",
            showarrow=False,
            font=dict(size=11, color="#333333"),
            xanchor="center", yanchor="top",
        )


def compress_track_coordinates(genes, threshold=50000, visual_gap=2000):
    """
    Compress large gaps between genes and add visual breaks between chromosomes.

    Genes are grouped by chromosome, with each chromosome's genes sorted by
    start position.  Chromosomes are ordered so that the one containing a GOI
    gene comes first, then remaining chromosomes ordered by descending gene
    count (most genes → most synteny evidence → shown first).

    Returns
    -------
    compressed_genes : list of dicts (with added 'start_plot', 'end_plot')
    breaks           : list of dicts (x, gap_size, text)
    """
    if not genes:
        return [], []

    # ---- Group by chromosome & order -----------------------------------
    from collections import defaultdict as _dd
    chrom_groups = _dd(list)
    for g in genes:
        chrom_groups[g["chrom"]].append(g)

    # Sort each chromosome group by start
    for chrom in chrom_groups:
        chrom_groups[chrom].sort(key=lambda g: g["start"])

    # Determine chromosome ordering: GOI chromosome first, then by gene count
    def _chrom_sort_key(chrom):
        has_goi = any(_is_goi_target_gene(g) or is_goi(g.get("name", "")) for g in chrom_groups[chrom])
        return (0 if has_goi else 1, -len(chrom_groups[chrom]))

    ordered_chroms = sorted(chrom_groups.keys(), key=_chrom_sort_key)

    # ---- Build linear sequence with compression ------------------------
    sorted_genes = []
    for chrom in ordered_chroms:
        sorted_genes.extend(chrom_groups[chrom])

    compressed = []
    breaks = []
    current_shift = 0
    CHROM_VISUAL_GAP = max(visual_gap * 5, 30000)  # wide gap between chromosomes

    for i, g in enumerate(sorted_genes):
        new_g = g.copy()

        if i > 0:
            prev = sorted_genes[i - 1]
            same_chrom = g["chrom"] == prev["chrom"]

            if same_chrom:
                gap = g["start"] - prev["end"]
                if gap > threshold:
                    remove = gap - visual_gap
                    current_shift += remove
                    prev_end_plot = prev["end"] - (current_shift - remove)
                    break_x = prev_end_plot + visual_gap / 2
                    breaks.append({
                        "x": break_x,
                        "gap_size": gap,
                        "text": (f"{gap / 1e6:.2f} Mb" if gap >= 1e6
                                 else f"{gap / 1e3:.0f} kb"),
                    })
            else:
                # Chromosome boundary — insert a visual chromosome-break gap
                prev_end_plot = prev["end"] - current_shift
                # Shift so the new chromosome starts CHROM_VISUAL_GAP after prev
                new_origin = prev_end_plot + CHROM_VISUAL_GAP
                actual_start = g["start"]
                current_shift = actual_start - new_origin

                break_x = prev_end_plot + CHROM_VISUAL_GAP / 2
                breaks.append({
                    "x": break_x,
                    "gap_size": 0,
                    "text": f"⧫ {g['chrom'][:16]}",
                    "is_chrom_break": True,
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
        if _is_goi_target_gene(g) or is_goi(g.get("name")) or is_goi(g.get("home_gene_id")):
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
    if not tree_file or not os.path.exists(tree_file):
        return

    try:
        if ETE3_AVAILABLE:
            t = Tree(tree_file)
        else:
            with open(tree_file) as fh:
                newick_str = fh.read().strip()
            t = _parse_newick(newick_str)
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
                        genome_pretty = f"{sp_name} ({genome_pretty})"
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
            text="<b>SynVoy GOI Phylogenetic Tree</b>",
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
    ap = argparse.ArgumentParser(description="SynVoy synteny plot")
    ap.add_argument("--home_bed",       required=True)
    ap.add_argument("--home_gff",       default=None)
    ap.add_argument("--query_bed",      default=None)
    ap.add_argument("--target_gffs",    nargs="*", default=[])
    ap.add_argument("--target_names",   nargs="*", default=[])
    ap.add_argument("--candidate_beds", nargs="*", default=[])
    ap.add_argument("--homology_tsvs",  nargs="*", default=[])
    ap.add_argument("--tree",           default=None)
    ap.add_argument("--sorted_genomes", default=None)
    ap.add_argument("--species_map",    default=None,
                    help="TSV mapping accession → species name")
    ap.add_argument("--gap_threshold",  type=int, default=50000, help="Min gap size to compress (bp)")
    ap.add_argument("--gap_visual_size",type=int, default=20000, help="Visual size of compressed gaps (bp)")
    ap.add_argument("--flank_fallback_bp",type=int, default=1000000, help="Fallback window if candidate genes miss GOI")
    ap.add_argument("--scale_bar_len",  type=int, default=10000, help="Length of the scale bar (bp)")
    ap.add_argument("--plot_width",     type=int, default=0, help="Total width of the output HTML plot (0=auto)")
    ap.add_argument("--plot_height",    type=int, default=0, help="Total height of the output HTML plot (0=auto)")
    ap.add_argument("--max_legend_entries", type=int, default=25, help="Maximum number of flanking genes to show in legend")
    ap.add_argument("--ribbon_alpha_dense", type=float, default=0.20, help="Alpha for flanking ribbons")
    ap.add_argument("--hide_goi_absent", action="store_true",
                    help="Hide target tracks with no GOI-like annotation when informative tracks exist")
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
            title="SynVoy Synteny Plot (Failed: empty home BED)",
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

    # Parse exon boundaries for home genes from the home GFF
    home_gene_names = {g["name"] for g in home_genes}
    home_exons = parse_home_gff_exons(args.home_gff, home_gene_names) if args.home_gff else {}
    for g in home_genes:
        if g["name"] in home_exons:
            g["exon_coords"] = home_exons[g["name"]]

    homology_map  = parse_homology_tsvs(args.homology_tsvs)

    goi_genome_colours, tree_target_order = parse_tree_clade_colours(args.tree)

    if args.sorted_genomes and os.path.exists(args.sorted_genomes):
        with open(args.sorted_genomes) as fh:
            tree_target_order = [line.strip().split("\t")[0] for line in fh if line.strip()]
            print(f"[plot] Overriding target order with {len(tree_target_order)} genomes from {args.sorted_genomes}")

    # -- 2. Build target tracks (matched by filename, not positional index)
    candidate_regions_by_genome = parse_candidate_regions(args.candidate_beds)

    target_tracks = []
    for gff_file in args.target_gffs:
        genome_id = clean_genome_name(
            os.path.basename(gff_file).replace(".gff", ""))
        genes_all = parse_target_gff(gff_file)
        candidate_regions = _match_regions_for_genome(candidate_regions_by_genome, genome_id)
        genes = filter_genes_to_candidate_regions(genes_all, candidate_regions)

        # Try restricting to GOI-containing candidate regions, but ONLY if
        # that doesn't discard most flanking context.  The GOI hit is often in
        # a tiny region separate from the main synteny block — blindly
        # restricting to GOI regions would throw away all flanking evidence.
        goi_candidate_regions = _candidate_regions_with_goi(candidate_regions, genes_all)
        if goi_candidate_regions:
            goi_only_genes = filter_genes_to_candidate_regions(genes_all, goi_candidate_regions)
            if len(goi_only_genes) >= max(3, len(genes) * 0.5):
                # GOI regions contain enough flanking context — use them.
                genes = goi_only_genes
                candidate_regions = goi_candidate_regions
            else:
                # GOI is isolated from flanking genes; keep ALL candidate
                # regions so flanking context is preserved.
                print(
                    f"[plot] {genome_id}: GOI-only regions have {len(goi_only_genes)} genes "
                    f"vs {len(genes)} from all regions — keeping all regions."
                )

        # If candidate regions exist but miss GOI, recover a GOI-centered context.
        if candidate_regions and not any(_is_goi_target_gene(g) for g in genes):
            fallback_genes = _select_goi_context_genes(genes_all, flank_bp=args.flank_fallback_bp)
            if fallback_genes:
                genes = fallback_genes
                print(
                    f"[plot] {genome_id}: candidate regions missed GOI; "
                    f"using GOI-centered fallback ({len(genes)} genes)."
                )

        # Richest-block fallback: if candidate regions yielded very few genes
        # (≤2), they likely cover low-quality scoring windows rather than the
        # block with the best synteny evidence.  Find the block containing
        # a GOI gene with the most flanking genes, or fall back to the globally
        # richest block.
        if len(genes) <= 2 and len(genes_all) > len(genes):
            import re as _re
            from collections import Counter
            block_counter = Counter()
            block_genes = defaultdict(list)
            for g in genes_all:
                gname = g.get("name", "")
                m = _re.search(r'_b(\d+)_', gname)
                if m:
                    bid = m.group(1)
                    block_counter[bid] += 1
                    block_genes[bid].append(g)

            if block_counter:
                # Blocks that contain at least one GOI gene
                goi_bids = set()
                for g in genes_all:
                    if _is_goi_target_gene(g):
                        m = _re.search(r'_b(\d+)_', g.get("name", ""))
                        if m:
                            goi_bids.add(m.group(1))

                if goi_bids:
                    # Pick richest block that contains a GOI
                    best_block = max(goi_bids, key=lambda bid: block_counter.get(bid, 0))
                else:
                    best_block = block_counter.most_common(1)[0][0]

                best_genes = block_genes[best_block]
                if len(best_genes) > len(genes):
                    genes = best_genes
                    print(
                        f"[plot] {genome_id}: candidate regions had ≤2 genes; "
                        f"using richest block b{best_block} ({len(genes)} genes)."
                    )

        # Focus on the most informative chromosomes for chromosome-level assemblies.
        # For scaffold/contig-level assemblies (many distinct contigs), flanking genes
        # may legitimately reside on different contigs — do NOT apply the restriction.
        all_chroms_in_gff = {g["chrom"] for g in genes_all}
        is_scaffold_assembly = len(all_chroms_in_gff) > 20

        if not is_scaffold_assembly and any(_is_goi_target_gene(g) for g in genes):
            goi_chroms = {g["chrom"] for g in genes if _is_goi_target_gene(g)}
            # Count flanking genes per chromosome (non-GOI genes)
            from collections import Counter as _Counter
            chrom_flank_counts = _Counter(
                g["chrom"] for g in genes if not _is_goi_target_gene(g)
            )

            # ── Synteny-aware GOI filter ──────────────────────────────
            # GOI hits on chromosomes without ANY flanking gene support
            # are almost certainly low-complexity/spurious matches (e.g.
            # proline-repeat hits matching the pro-peptide region).
            # Without synteny context, orthology cannot be established.
            # Drop them unless they are the ONLY GOI hits in this genome.
            unsupported_goi_chroms = {
                ch for ch in goi_chroms
                if chrom_flank_counts.get(ch, 0) == 0
            }
            supported_goi_chroms = goi_chroms - unsupported_goi_chroms

            if unsupported_goi_chroms and supported_goi_chroms:
                # There are GOI hits with synteny support elsewhere —
                # safe to drop the unsupported ones.
                n_dropped = sum(
                    1 for g in genes
                    if g["chrom"] in unsupported_goi_chroms
                )
                genes = [
                    g for g in genes
                    if g["chrom"] not in unsupported_goi_chroms
                ]
                print(
                    f"[plot] {genome_id}: dropped {n_dropped} unsupported GOI "
                    f"hits on {unsupported_goi_chroms} (no flanking genes = "
                    f"no synteny evidence)"
                )
                # Refresh GOI chroms after drop
                goi_chroms = supported_goi_chroms

            # Keep GOI chromosome(s) + any chromosome with >=3 flanking genes
            # (i.e. chromosomes with real synteny evidence)
            important_chroms = set(goi_chroms)
            for ch, cnt in chrom_flank_counts.items():
                if cnt >= 3:
                    important_chroms.add(ch)

            # Only filter if we'd still have enough genes
            filtered = [g for g in genes if g["chrom"] in important_chroms]
            if len(filtered) >= len(genes) * 0.3 or len(important_chroms) <= 3:
                if len(filtered) < len(genes):
                    dropped_chroms = {g["chrom"] for g in genes} - important_chroms
                    print(
                        f"[plot] {genome_id}: keeping {len(important_chroms)} informative "
                        f"chromosomes ({len(filtered)} genes), dropped {dropped_chroms}"
                    )
                genes = filtered

        # Don't skip target track if there are no genes found; we want to show it's empty
        genes.sort(key=lambda g: g["start"])
        # Use species name from mapping if available, format: "Species name (accession)"
        display = genome_id
        for acc, sp_name in species_map.items():
            if acc in genome_id:
                display = f"<i>{sp_name}</i> ({genome_id})"
                break
        gene_chroms = sorted({g["chrom"] for g in genes}) if genes else []
        if len(gene_chroms) == 1:
            target_chrom = gene_chroms[0]
        elif len(gene_chroms) > 1:
            target_chrom = f"{len(gene_chroms)} chr"
        else:
            target_chrom = candidate_regions[0][0] if candidate_regions else "unknown"
        
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

    for tt in target_tracks:
        tt["goi_status"] = _track_goi_status({"genes": tt.get("genes", [])})

    hidden_absent_tracks = 0
    informative_tracks = [t for t in target_tracks if t.get("goi_status") != "absent"]
    if args.hide_goi_absent and informative_tracks and len(informative_tracks) < len(target_tracks):
        hidden_absent_tracks = len(target_tracks) - len(informative_tracks)
        target_tracks = informative_tracks
        print(f"[plot] Hid {hidden_absent_tracks} GOI-absent target tracks from overview plot")

    ambiguous_track_count = sum(1 for t in target_tracks if t.get("goi_status") == "ambiguous")
    resolved_track_count = sum(1 for t in target_tracks if t.get("goi_status") == "resolved")

    # -- 3. Colour map ---------------------------------------------------

    gene_colours = assign_gene_colours(home_genes, query_intervals)

    # -- 4. Assemble track list & Compress -------------------------------
    # Each genome gets one horizontal track bar.  Genes from different
    # chromosomes are placed side-by-side with clear visual separators
    # (handled by compress_track_coordinates).

    home_chrom = home_genes[0]["chrom"]

    raw_tracks = [{
        "label":        f"Home genome ({home_chrom})",
        "genes":        home_genes,
        "is_home":      True,
        "genome_id":    "home",
        "goi_status":   "resolved",
    }]
    for tt in target_tracks:
        genes = tt["genes"]
        raw_tracks.append({
            "label":        tt['display_name'],
            "genes":        genes,
            "is_home":      False,
            "genome_id":    tt["genome_id"],
            "goi_status":   tt.get("goi_status", _track_goi_status({"genes": genes})),
        })

    all_tracks = []
    for track in raw_tracks:
        # 1. Compress
        c_genes, breaks = compress_track_coordinates(track["genes"], threshold=args.gap_threshold, visual_gap=args.gap_visual_size)
        track["genes"]  = c_genes
        track["breaks"] = breaks
        
        # 2. Find Anchor (GOI center) to align at x=0
        anchor = get_anchor_center(c_genes)
        track["offset"] = anchor  # This effectively centers the plot on the GOI
        
        all_tracks.append(track)

    n_tracks = len(all_tracks)

    # -- 5. Layout geometry -----------------------------------------------

    GENE_H      = 0.22          # gene arrow height  (flatter = gggenomes-like)
    TRACK_SPACE = 1.1           # vertical pitch between tracks
    RIBBON_GAP  = 0.08          # gap between gene arrow and ribbon edge

    fig = go.Figure()

    # Simple Y positions: uniform spacing, top to bottom.
    track_y = [(n_tracks - 1 - ti) * TRACK_SPACE for ti in range(n_tracks)]

    # -- 5a. Track background bands --------------------------------------
    # Split the background bar at chromosome breaks so each chromosome
    # gets its own distinct background rectangle.

    for ti, track in enumerate(all_tracks):
        yb = track_y[ti]
        x_off = track["offset"]
        if not track["genes"]:
            continue

        # Collect chromosome break x positions (in plot-offset coords)
        chrom_break_xs = sorted(
            brk["x"] - x_off for brk in track.get("breaks", [])
            if brk.get("is_chrom_break")
        )

        if not chrom_break_xs:
            # No chromosome breaks — single background rectangle
            x_min = min(g["start_plot"] for g in track["genes"]) - x_off - 1000
            x_max = max(g["end_plot"]   for g in track["genes"]) - x_off + 1000
            fig.add_shape(
                type="rect",
                x0=x_min, x1=x_max, y0=yb - 0.02, y1=yb + GENE_H + 0.02,
                fillcolor=TRACK_BG_CLR, line=dict(width=0), layer="below",
            )
        else:
            # Draw a separate background rectangle for each chromosome segment
            all_gene_xs = [(g["start_plot"] - x_off, g["end_plot"] - x_off)
                           for g in track["genes"]]
            # Create segment boundaries: [start, brk1, brk2, ..., end]
            boundaries = [min(s for s, e in all_gene_xs) - 1000]
            boundaries.extend(chrom_break_xs)
            boundaries.append(max(e for s, e in all_gene_xs) + 1000)

            for si in range(len(boundaries) - 1):
                seg_left  = boundaries[si]
                seg_right = boundaries[si + 1]
                # Inset from break positions to create visible gap
                inset = 2500
                seg_x0 = seg_left  + (inset if si > 0 else 0)
                seg_x1 = seg_right - (inset if si < len(boundaries) - 2 else 0)
                if seg_x0 < seg_x1:
                    fig.add_shape(
                        type="rect",
                        x0=seg_x0, x1=seg_x1,
                        y0=yb - 0.02, y1=yb + GENE_H + 0.02,
                        fillcolor=TRACK_BG_CLR, line=dict(width=0),
                        layer="below",
                    )

    # -- 5b. Ribbons (draw first so they sit behind genes) ---------------
    # Ribbons connect consecutive tracks (one genome per track).
    for ti in range(len(all_tracks) - 1):
        upper = all_tracks[ti]
        lower = all_tracks[ti + 1]
        y_u = track_y[ti]
        y_l = track_y[ti + 1]
        y_ribbon_top = y_u - RIBBON_GAP
        y_ribbon_bot = y_l + GENE_H + RIBBON_GAP

        for lg in lower["genes"]:
            home_id = lg.get("home_gene_id", "")
            if not home_id:
                continue
            ribbon_alpha = args.ribbon_alpha_dense
            for ug in upper["genes"]:
                u_name = ug["name"]
                u_home = ug.get("home_gene_id", u_name)
                match = (u_home == home_id or u_name == home_id
                         or (is_goi(u_home) and is_goi(home_id))
                         or (is_goi(u_name) and is_goi(home_id)))
                if match:
                    colour = gene_colours.get(home_id,
                             gene_colours.get(u_name, UNMATCHED_CLR))
                    if is_goi(home_id):
                        colour = _goi_colour_for_genome(
                            lower["genome_id"], goi_genome_colours)
                        if _is_goi_target_gene(lg) and not _is_resolved_goi_target_gene(lg):
                            ribbon_alpha = 0.10
                    add_ribbon(fig, ug, lg,
                               upper["offset"], lower["offset"],
                               y_ribbon_top, y_ribbon_bot,
                               colour, alpha=ribbon_alpha)

    # -- 5c. Gene arrows -------------------------------------------------
    legend_shown = set()

    for ti, track in enumerate(all_tracks):
        yb    = track_y[ti]
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
            goi_like = _is_goi_target_gene(gene)
            resolved_goi = _is_resolved_goi_target_gene(gene)
            ambiguous_goi = goi_like and not resolved_goi
            confidence = (gene.get("confidence") or "").upper()
            goi_f = _is_goi_target_gene(gene) if not track["is_home"] else (is_goi(name) or is_goi(home_id))

            # --- colour ---
            if resolved_goi:
                colour = _goi_colour_for_genome(track["genome_id"], goi_genome_colours)
                bclr = GOI_BORDER
                bw = 2.8 if confidence == "HIGH" else 2.2
                dash = "solid"
            elif ambiguous_goi:
                base_goi = _goi_colour_for_genome(track["genome_id"], goi_genome_colours)
                colour = _hex_to_rgba(base_goi, 0.32)
                bclr, bw, dash = GOI_BORDER, 1.8, "dash"
            elif home_id in gene_colours:
                colour = gene_colours[home_id]
                bclr, bw = "rgba(0,0,0,0.35)", 1
                dash = "dot" if confidence == "LOW" else "solid"
            elif name in gene_colours:
                colour = gene_colours[name]
                bclr, bw = "rgba(0,0,0,0.35)", 1
                dash = "dot" if confidence == "LOW" else "solid"
            else:
                colour = UNMATCHED_CLR
                bclr, bw = "rgba(0,0,0,0.15)", 0.5
                dash = "dot" if confidence == "LOW" else "solid"

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
                n_ex_home = len(gene.get("exon_coords", []))
                if n_ex_home <= 1:
                    n_ex_home = gene.get("n_exons", 0)
                if n_ex_home and n_ex_home > 1:
                    hover += f"<br>Exons: {n_ex_home}"
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
                n_ex = len(gene.get("exon_coords", []))
                if n_ex <= 1:
                    n_ex = gene.get("n_exons", 0)
                if n_ex and n_ex > 1:
                    hover += f"<br>Exons: {n_ex}"
                if confidence:
                    hover += f"<br>Confidence: {confidence}"
                if gene.get("goi_class"):
                    hover += f"<br>GOI class: {gene['goi_class'].replace('_', ' ')}"
                if gene.get("evidence_type"):
                    hover += f"<br>Evidence: {gene['evidence_type'].replace('_', ' ')}"
                if gene.get("synteny_context"):
                    hover += f"<br>Synteny: {gene['synteny_context'].replace('_', ' ')}"
                if gene.get("query_coverage") is not None:
                    hover += f"<br>Query coverage: {gene['query_coverage'] * 100:.1f}%"
                hover += (f"<br>{gene['chrom']}:{gene['start']:,}-{gene['end']:,}"
                          f"<br>Strand: {gene['strand']}")
                if resolved_goi:
                    hover += "<br><b>GENE OF INTEREST</b>"
                elif ambiguous_goi:
                    hover += "<br><b>GOI-LIKE / AMBIGUOUS</b>"

            # --- legend (one entry per home-gene name) ---
            lg_key = home_id if home_id else name
            show_leg = False
            if lg_key not in legend_shown:
                if goi_like or len(legend_shown) < args.max_legend_entries:
                    show_leg = True
                    legend_shown.add(lg_key)

            add_gene(fig, gene, x_off, yb, GENE_H, colour, bclr, bw,
                     hover, show_leg, lg_key, line_dash=dash)

    # -- 5c2. Absent-GOI placeholder for targets without GOI gene --------
    for ti, track in enumerate(all_tracks):
        if track["is_home"]:
            continue
        if track.get("goi_status") == "absent" and track["genes"]:
            yb = track_y[ti]
            # Draw a dashed-outline "?" box at x=0 (GOI center)
            dash_w = 2000
            fig.add_shape(
                type="rect",
                x0=-dash_w / 2, x1=dash_w / 2,
                y0=yb, y1=yb + GENE_H,
                fillcolor="rgba(227,26,28,0.08)",
                line=dict(color=GOI_COLOUR, width=1.5, dash="dash"),
            )
            fig.add_annotation(
                x=0, y=yb + GENE_H / 2,
                text="<b>?</b>",
                showarrow=False,
                font=dict(size=12, color=GOI_COLOUR),
                xanchor="center", yanchor="middle",
                hovertext=f"<b>GOI not found</b><br>{track['label']}",
            )

    # -- 5d. Gene labels -------------------------------------------------
    for ti, track in enumerate(all_tracks):
        yb    = track_y[ti]
        x_off = track["offset"]
        genes_in_track = track["genes"]
        n_genes = len(genes_in_track)

        # Build candidate labels with their x positions, then filter overlaps
        label_candidates = []  # [(x_center, label_text, gene, is_goi, priority)]
        for gene in genes_in_track:
            name    = gene["name"]
            home_id = gene.get("home_gene_id", name)
            goi_f = _is_goi_target_gene(gene) if not track["is_home"] else (is_goi(name) or is_goi(home_id))
            resolved_goi = _is_resolved_goi_target_gene(gene) if not track["is_home"] else goi_f

            # In dense tracks, only label GOI and matched flanking genes
            has_colour = (home_id in gene_colours or name in gene_colours)
            if n_genes > 12 and not goi_f and not has_colour:
                continue  # skip labels for unmatched genes in dense tracks

            if not track["is_home"]:
                if goi_f:
                    # GOI on target: show 'GOI #N' instead of 'copy_N'
                    label = clean_gene_label(name, keep_goi_prefix=True)
                    if not label or label == 'GOI':
                        label = clean_gene_label(home_id, keep_goi_prefix=True)
                    if not resolved_goi:
                        label = "~ " + label
                else:
                    label = clean_gene_label(_preferred_target_label(gene))
                    if not label and home_id:
                        label = clean_gene_label(home_id)
            else:
                label = _preferred_home_label(gene, home_products)

            g_start, g_end = _get_coords(gene)
            xc = (g_start + g_end) / 2 - x_off
            # Priority: GOI=0 (always show), coloured flanking=1, other=2
            priority = 0 if goi_f else (1 if has_colour else 2)
            label_candidates.append((xc, label, gene, goi_f, priority))

        # Sort by priority (highest first), then by x position
        label_candidates.sort(key=lambda c: (c[4], c[0]))

        # Collision detection: skip labels that would overlap with already-placed ones
        # Estimate label width in x-data units.
        # For rotated labels (-35°), the horizontal footprint is
        #   w_proj = w * cos(35°) ≈ 0.82 * w, but the diagonal sweep
        # still causes visual collisions — use a 70% factor as effective width.
        fsize = max(6, 11 - (n_genes // 6))
        char_width = 300 + (fsize * 25)
        import math
        rotation_factor = 0.70  # accounts for diagonal sweep of -35° text
        placed_ranges = []  # list of (x_left, x_right) of placed labels

        for xc, label, gene, goi_f, priority in label_candidates:
            g_start, g_end = _get_coords(gene)
            gw = g_end - g_start
            is_rotated = True  # labels always rotated
            est_width = int(len(label) * char_width * rotation_factor)
            lbl_left  = xc - est_width / 2
            lbl_right = xc + est_width / 2

            # Minimum clearance between labels
            margin = 350 if is_rotated else 450

            # Check if this label overlaps with any already placed
            overlaps = any(
                lbl_left < pr + margin and lbl_right > pl - margin
                for pl, pr in placed_ranges
            )
            if overlaps and not goi_f:
                continue  # skip non-GOI labels that overlap

            placed_ranges.append((lbl_left, lbl_right))
            add_label(fig, gene, x_off, yb, GENE_H, label,
                      fsize=fsize, is_goi_flag=goi_f)

    # -- 5e. Track labels (left margin) ----------------------------------
    for ti, track in enumerate(all_tracks):
        yb = track_y[ti]
        track_label = f"<b>{track['label']}</b>"
        if not track["is_home"] and track.get("goi_status") == "absent":
            track_label += "<br><span style='color:#d32f2f;font-size:9px'>✗ GOI absent</span>"
        elif not track["is_home"] and track.get("goi_status") == "ambiguous":
            track_label += "<br><span style='color:#d97706;font-size:9px'>~ GOI ambiguous</span>"
        fig.add_annotation(
            x=-0.01, y=yb + GENE_H / 2,
            text=track_label,
            showarrow=False,
            font=dict(size=11, color="black"),
            xref="paper", yref="y",
            xanchor="right", yanchor="middle",
        )

    # -- 5f. Gap breaks & chromosome labels ------------------------------
    for ti, track in enumerate(all_tracks):
        yb    = track_y[ti]
        x_off = track["offset"]
        for brk in track.get("breaks", []):
            draw_gap_break(fig, brk["x"] - x_off, yb, GENE_H, brk["text"],
                           is_chrom_break=brk.get("is_chrom_break", False))
        # Draw chromosome labels below each segment
        _draw_chrom_labels(fig, track, ti, track_y, x_off, GENE_H)

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

    fig_height = max(400, n_tracks * 140 + 60)

    # Compute Y range from actual track positions
    y_min_pos = min(track_y) if track_y else 0
    y_max_pos = max(track_y) if track_y else 0
    subtitle_bits = [
        "Genes coloured by homology group",
        "* = resolved GOI",
        "dashed GOI = ambiguous/tandem family member",
        "V-notches = exon boundaries",
        "ribbons connect orthologs",
        "// = compressed gaps",
    ]
    if hidden_absent_tracks:
        subtitle_bits.append(f"{hidden_absent_tracks} GOI-absent track(s) hidden")
    if ambiguous_track_count:
        subtitle_bits.append(f"{ambiguous_track_count} ambiguous track(s)")

    if args.plot_height > 0:
        fig_height = args.plot_height
    else:
        fig_height = max(450, n_tracks * 140 + 80)
        
    if args.plot_width > 0:
        fig_width = args.plot_width
    else:
        # Scale dynamically between 1800 and 5000 based on window width
        estimated_needed = max(2000, int((x_max - x_min) / 450))
        fig_width = min(5000, estimated_needed)

    fig.update_layout(
        title=dict(
            text=("<b>SynVoy Synteny Plot</b>"
                  f"<br><sup>{' | '.join(subtitle_bits)}</sup>"),
            x=0.5, font=dict(size=15),
        ),
        height=fig_height,
        width=fig_width,
        xaxis=dict(
            title="", # No title since numbers are relative/discontinuous
            showgrid=False, 
            zeroline=False,
            showticklabels=False, # Hide ticks as they are discontinuous
            range=x_range,
        ),
        yaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[y_min_pos - 0.6,
                   y_max_pos + GENE_H + 1.0],
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
    
    # Add Scale Bar
    # Place it in bottom right? Or top left?
    # Let's put it at bottom right
    scale_len = args.scale_bar_len
    sb_x1 = x_max
    sb_x0 = x_max - scale_len
    sb_y  = y_min_pos - 0.4
    
    fig.add_shape(
        type="line",
        x0=sb_x0, x1=sb_x1, y0=sb_y, y1=sb_y,
        line=dict(color="black", width=3),
    )
    fig.add_annotation(
        x=(sb_x0 + sb_x1)/2, y=sb_y - 0.1,
        text=f"<b>{_format_bp_label(scale_len)}</b>",
        showarrow=False,
        font=dict(size=10, color="black"),
        yanchor="top"
    )

    fig.write_html(args.output)
    print(f"Synteny plot saved to {args.output}")
    print(f"  Tracks: {n_tracks} ({n_tracks - 1} target genomes)")
    print(f"  GOI tracks: {resolved_track_count} resolved, {ambiguous_track_count} ambiguous")
    if hidden_absent_tracks:
        print(f"  Hidden absent tracks: {hidden_absent_tracks}")
    print(f"  Gap compression: active (>{args.gap_threshold} bp -> {args.gap_visual_size} bp visual)")

    # -- 7. Tree plot (separate HTML) ------------------------------------
    tree_output = args.output.replace("_synteny_plot.html", "_tree.html")
    if tree_output == args.output:
        tree_output = args.output.replace(".html", "_tree.html")
    _render_tree_html(args.tree, goi_genome_colours, tree_output,
                      species_map=species_map)


if __name__ == "__main__":
    main()
