#!/usr/bin/env python3
"""
plot_synteny.py  –  SVG-based synteny visualization for SynVoy

Layout
──────
  •  Home genome at top, target genomes below (ordered by phylogenetic distance)
  •  Gene models: exon blocks connected by intron lines, with directional arrows
  •  Smooth bezier-curve ribbons between homologous genes in adjacent tracks
  •  GOI highlighted with warm/red clade colours from the phylogenetic tree
  •  Flanking genes share a consistent colour derived from the home-genome name
  •  Interactive: hover tooltips, click-to-highlight orthologs, zoom controls

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
  --output          Self-contained interactive HTML file (SVG)
"""

import argparse
import collections
import colorsys
import json
import math
import os
import re
import sys
from collections import defaultdict
from html import escape as _html_escape
from urllib.parse import unquote

# Shared helpers (single source of truth for BED/GFF parsing).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import parse_bed, parse_gff_attributes as _parse_gff_attrs  # noqa: E402,F401

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

# Okabe-Ito / Wong 2011 colorblind-safe palette for publication SVG
PUB_PALETTE = [
    "#0072B2",  # blue
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#E69F00",  # orange
    "#F0E442",  # yellow
    "#666666",  # neutral gray
]
GOI_PUB_COLOUR = "#E64B35"  # Nature-red for publication GOI
GOI_PUB_BORDER = "#B71C1C"  # dark vermillion

GOI_COLOUR    = "#e31a1c"   # bright red (default for GOI)
GOI_BORDER    = "#8b0000"   # dark red
UNMATCHED_CLR = "#d9d9d9"   # light gray
TRACK_BG_CLR  = "#f3f5f8"   # very light blue-gray track background

# Anchor-grid / matrix cells show the identity number alone (it always fits a
# cell and reads cleanly). Query-coverage is a second number only worth the
# clutter when it's LOW — a short high-identity hit would otherwise masquerade
# as a full-length ortholog. Below this fraction we append "(cov%)".
COVERAGE_FLAG_THRESHOLD = 0.80


# ======================================================================
# Parsing helpers
# ======================================================================

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
                "model_status": attrs.get("ModelStatus", ""),
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
    MAX_GOI_PER_GENOME = _MAX_GOI_PER_GENOME
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

# Max GOI/toxin entries drawn per genome (set from --max_goi_per_genome in main()).
# Raise for tandem-array loci (e.g. a 31-copy GR1 toxin cluster) so the array isn't truncated.
_MAX_GOI_PER_GENOME = 10


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


def _synthesize_home_goi_gene(home_genes, query_intervals, home_gff_path):
    """Inject a synthetic GOI gene into `home_genes` if none is present.

    The home synteny-block BED only contains *flanking* genes — by design.
    For short queries (e.g. melittin, 26 aa) the GOI itself often sits
    inside a much larger container gene (e.g. LOC726866, 17 kb), so the
    `identify_goi_names()` size filter rejects everything and the home
    track ends up with a silent gap where the GOI should be.

    This helper scans `home_gff` for the gene that best overlaps the
    query span, then appends a fresh entry to `home_genes` so the
    downstream track-building code (and `identify_goi_names`) treats it
    like any other home gene. Mutates `home_genes` in place; returns the
    synthetic gene's `name` (or None if nothing could be synthesized).
    """
    if not query_intervals:
        return None
    if not home_gff_path or home_gff_path == "NO_GFF" or not os.path.exists(home_gff_path):
        return None

    chrom = query_intervals[0]["chrom"]
    qstart = min(q["start"] for q in query_intervals)
    qend = max(q["end"] for q in query_intervals)
    qstrand = query_intervals[0].get("strand", "+")
    # Use the true CDS coverage (sum of exon lengths, NOT min..max which
    # includes introns). For melittin: sum ≈ 210 bp; without this fix the
    # max_size grew to ~40 kb and the 17 kb container LOC726866 was
    # incorrectly accepted as "small enough" to be the GOI itself.
    q_cds_total = max(1, sum(q["end"] - q["start"] for q in query_intervals))
    max_size = min(max(q_cds_total * 10, 2000), 5000)

    # Skip if a small home gene already covers the query — `identify_goi_names`
    # will mark that one. We only inject when there's a true gap.
    for g in home_genes:
        if g["chrom"] != chrom:
            continue
        if (g["end"] - g["start"]) > max_size:
            continue
        if g["start"] < qend and g["end"] > qstart:
            return g["name"]

    # Find the best-overlapping gene in home_gff that's *small* enough to be
    # the GOI itself, not a container.
    best = None
    best_ov = 0
    try:
        with open(home_gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "gene" or p[0] != chrom:
                    continue
                try:
                    gs = int(p[3]) - 1
                    ge = int(p[4])
                except ValueError:
                    continue
                if (ge - gs) > max_size:
                    continue
                ov = max(0, min(ge, qend) - max(gs, qstart))
                if ov <= 0:
                    continue
                attrs = {}
                for kv in p[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k] = unquote(v)
                gid = attrs.get("ID", "") or attrs.get("Name", "") or "GOI"
                gname = attrs.get("Name") or attrs.get("gene") or gid
                if ov > best_ov:
                    best = {
                        "chrom": chrom, "start": gs, "end": ge,
                        "name": gid,
                        "strand": p[6] if p[6] in {"+", "-"} else qstrand,
                        "display_name": gname,
                    }
                    best_ov = ov
    except OSError as exc:
        print(f"[plot] WARN: could not scan home_gff for synthetic GOI gene: {exc}",
              file=sys.stderr)
        return None

    if best is None:
        # Fallback: synthesize purely from query_bed coords (no symbol).
        best = {
            "chrom": chrom, "start": qstart, "end": qend,
            "name": "GOI", "strand": qstrand, "display_name": "GOI",
        }

    home_genes.append(best)
    home_genes.sort(key=lambda g: g["start"])
    print(
        f"[plot] Synthesized home GOI gene '{best['display_name']}' at "
        f"{chrom}:{best['start']:,}-{best['end']:,} (strand {best['strand']}) "
        f"— home_bed had no row for the GOI itself.",
        file=sys.stderr,
    )
    return best["name"]


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
            m = re.match(r'copy_(\d+)', suffix)
            if m:
                return f"GOI #{m.group(1)}"
            m = re.match(r'(.*?)_copy_(\d+)', suffix)
            if m:
                return f"GOI #{m.group(2)}"
            return f"GOI {suffix}" if suffix else "GOI"
        return name[4:]
    return name


# ======================================================================
# Drawing helpers (shared)
# ======================================================================

def _hex_to_rgba(hexc, alpha):
    hexc = hexc.lstrip("#")
    r, g, b = int(hexc[0:2], 16), int(hexc[2:4], 16), int(hexc[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _darken_hex(hexc, factor=0.7):
    """Darken a hex colour by the given factor."""
    hexc = hexc.lstrip("#")
    r = int(int(hexc[0:2], 16) * factor)
    g = int(int(hexc[2:4], 16) * factor)
    b = int(int(hexc[4:6], 16) * factor)
    return f"#{min(r,255):02x}{min(g,255):02x}{min(b,255):02x}"


def _lerp_hex(a, b, t):
    """Linear-interpolate between two hex colours (t in [0,1])."""
    a = a.lstrip("#"); b = b.lstrip("#")
    if len(a) != 6 or len(b) != 6:
        return "#cccccc"
    ar, ag, ab = int(a[0:2], 16), int(a[2:4], 16), int(a[4:6], 16)
    br, bg, bb = int(b[0:2], 16), int(b[2:4], 16), int(b[4:6], 16)
    t = max(0.0, min(1.0, t))
    r = round(ar + (br - ar) * t)
    g = round(ag + (bg - ag) * t)
    bv = round(ab + (bb - ab) * t)
    return f"#{r:02x}{g:02x}{bv:02x}"


def _shade_by_identity(base, identity):
    """Lighten `base` toward white as identity drops (matrix-plot idiom).

    A 100 %-identity ortholog keeps the full homology colour; a divergent
    hit fades toward white so the reader can gauge similarity at a glance.
    """
    ident = max(0.0, min(100.0, float(identity or 0)))
    t = (100.0 - ident) / 100.0
    return _lerp_hex(base, "#ffffff", min(0.6, t * 0.7))


def _is_dark_hex(hexc):
    """Return True when a hex colour is dark enough for white text."""
    if not hexc:
        return False
    hexc = hexc.lstrip("#")
    if len(hexc) != 6:
        return False
    r = int(hexc[0:2], 16)
    g = int(hexc[2:4], 16)
    b = int(hexc[4:6], 16)
    luminance = (0.2126 * r) + (0.7152 * g) + (0.0722 * b)
    return luminance < 140


def _get_coords(gene):
    return gene.get("start_plot", gene["start"]), gene.get("end_plot", gene["end"])


def _dominant_chrom_span(items, goi_chrom=None):
    """Pick a single representative scaffold for a row's gene-position gutter.

    ``items`` is a list of ``(chrom, mid)`` pairs (one per recovered gene). On a
    fragmented assembly a neighbourhood can straddle >1 scaffold; rather than
    blanking the accession (the old behaviour, which left rows like *Bombus
    impatiens* showing only "1.51-1.68 Mb"), pick the scaffold carrying the GOI
    (if present) else the one carrying the most genes, and restrict lo/hi to that
    scaffold so the range stays meaningful. Returns
    ``(chrom, lo, hi, n_other_scaffolds)`` — ``n_other`` is how many *additional*
    scaffolds the row's genes fall on (0 = contiguous)."""
    items = [(c, m) for c, m in items if m is not None]
    if not items:
        return "", 0.0, 0.0, 0
    chrom_counts = collections.Counter(c for c, _ in items if c)
    if goi_chrom and goi_chrom in chrom_counts:
        primary = goi_chrom
    elif chrom_counts:
        primary = max(chrom_counts, key=lambda c: (chrom_counts[c], c))  # most genes; name tie-break
    else:
        primary = ""
    on = [m for c, m in items if c == primary] if primary else [m for _, m in items]
    lo, hi = (min(on), max(on)) if on else (0.0, 0.0)
    n_other = sum(1 for c in chrom_counts if c != primary)
    return primary, lo, hi, n_other


def _span_gutter_label(chrom, lo, hi, n_other=0):
    """Format a right-gutter location label, e.g. ``NC_037641.1: 12.37-12.44 Mb``.
    A single gene collapses to one coordinate; extra scaffolds add ``(+N)``."""
    if hi <= lo and not chrom:
        return ""
    rng = f"{lo/1e6:.2f}-{hi/1e6:.2f} Mb" if hi > lo else f"{lo/1e6:.2f} Mb"
    extra = f" (+{n_other})" if n_other else ""
    return f"{chrom + ': ' if chrom else ''}{rng}{extra}"


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
# Coordinate compression
# ======================================================================

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
    chrom_groups = defaultdict(list)
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
    # `max_end` is the furthest original end-coordinate reached so far on the
    # current chromosome — NOT just the previous gene's end. Genes overlap and
    # nest (a long readthrough/lncRNA can span several short neighbours), so
    # measuring the gap from the previous gene would compress the gap *after a
    # short nested gene* while the enclosing long gene still extends across it —
    # drawing that gene right over the compressed gap (visually implying it is
    # ~gap-bp longer than it is). Measuring from `max_end` only compresses
    # regions no gene spans, so no gene is ever drawn over a break.
    max_end = None
    CHROM_VISUAL_GAP = max(visual_gap * 5, 30000)  # wide gap between chromosomes

    for i, g in enumerate(sorted_genes):
        new_g = g.copy()

        if i > 0:
            prev = sorted_genes[i - 1]
            same_chrom = g["chrom"] == prev["chrom"]

            if same_chrom:
                gap = g["start"] - max_end
                if gap > threshold:
                    remove = gap - visual_gap
                    current_shift += remove
                    prev_end_plot = max_end - (current_shift - remove)
                    break_x = prev_end_plot + visual_gap / 2
                    breaks.append({
                        "x": break_x,
                        "gap_size": gap,
                        "text": (f"{gap / 1e6:.2f} Mb" if gap >= 1e6
                                 else f"{gap / 1e3:.0f} kb"),
                    })
            else:
                # Chromosome boundary — insert a visual chromosome-break gap
                prev_end_plot = max_end - current_shift
                # Shift so the new chromosome starts CHROM_VISUAL_GAP after prev
                new_origin = prev_end_plot + CHROM_VISUAL_GAP
                actual_start = g["start"]
                current_shift = actual_start - new_origin
                max_end = None  # reset furthest extent for the new chromosome

                break_x = prev_end_plot + CHROM_VISUAL_GAP / 2
                breaks.append({
                    "x": break_x,
                    "gap_size": 0,
                    "text": f"◆ {g['chrom'][:16]}",
                    "is_chrom_break": True,
                })

        new_g["start_plot"] = g["start"] - current_shift
        new_g["end_plot"]   = g["end"]   - current_shift
        compressed.append(new_g)
        if max_end is None or g["end"] > max_end:
            max_end = g["end"]

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


def _reverse_complement_track(track):
    """Reflect a track around its anchor (offset) so it reads in the opposite
    direction: mirror every gene's plotted position, flip its strand/arrow, and
    mirror its exon model within the gene. Used to align a natively
    reverse-oriented scaffold to the home reference (see _orient_tracks_to_home).
    """
    off = track.get("offset", 0.0)
    for g in track["genes"]:
        s_plot, e_plot = g["start_plot"], g["end_plot"]
        g["start_plot"] = 2 * off - e_plot
        g["end_plot"]   = 2 * off - s_plot
        g["strand"]     = "-" if g.get("strand", "+") == "+" else "+"
        # Mirror the exon/intron model within the gene's own native span so the
        # fine structure (and terminal-exon arrow tip) reads in the new direction.
        ex = g.get("exon_coords") or []
        if ex:
            gs, ge = g["start"], g["end"]
            g["exon_coords"] = [(gs + (ge - ee), gs + (ge - es)) for (es, ee) in reversed(ex)]
    for brk in track.get("breaks", []):
        brk["x"] = 2 * off - brk["x"]
    track["flipped"] = True


def _orient_tracks_to_home(all_tracks, min_anchors=4, min_concordance=0.15):
    """Reverse-complement (in plot space) any target track whose flanking-anchor
    order runs opposite to the home reference, so no row appears mirror-imaged.

    Orientation is decided from the *flanking* genes only (toxin/GOI copies have
    mixed strands and are unreliable): for genes shared with home, we measure the
    rank concordance between each track's plotted order and the home's. Negative
    concordance ⇒ the scaffold is assembled antiparallel to home ⇒ flip it.
    Tracks with too few shared anchors or an ambiguous (|concordance| small)
    signal are left untouched. Returns the number of tracks flipped.
    """
    if not all_tracks:
        return 0
    home = all_tracks[0]
    home_center = {}
    for g in home["genes"]:
        nm = g.get("name")
        if nm and nm not in home_center:
            home_center[nm] = (g["start_plot"] + g["end_plot"]) / 2.0

    def _concordance(pairs):
        conc = disc = 0
        for i in range(len(pairs)):
            ti, hi = pairs[i]
            for j in range(i + 1, len(pairs)):
                tj, hj = pairs[j]
                dt, dh = ti - tj, hi - hj
                if dt == 0 or dh == 0:
                    continue
                if (dt > 0) == (dh > 0):
                    conc += 1
                else:
                    disc += 1
        tot = conc + disc
        return (conc - disc) / tot if tot else 0.0

    n_flipped = 0
    for track in all_tracks[1:]:
        if track.get("is_home"):
            continue
        pairs = []
        for g in track["genes"]:
            if _is_goi_target_gene(g) or is_goi(g.get("name")):
                continue  # flanking only
            hid = g.get("home_gene_id", "")
            if hid in home_center:
                pairs.append(((g["start_plot"] + g["end_plot"]) / 2.0, home_center[hid]))
        if len(pairs) < min_anchors:
            continue
        if _concordance(pairs) < -min_concordance:
            _reverse_complement_track(track)
            n_flipped += 1
    return n_flipped


def _widen_sparse_plot(all_tracks, target_coverage=0.25, max_factor=4.0):
    """Inflate every gene's visual width around its center when the plot has
    pixel headroom — i.e. when the median per-track gene coverage of the plot
    range is small.

    Centers don't move, so ribbon endpoints stay aligned. Strand-split /
    bumping in `_assign_sub_tracks` picks up any overlaps the widening
    introduces. Dense plots (many genes per track) are untouched.

    Returns the applied factor (1.0 if no widening was done).
    """
    all_x = []
    for t in all_tracks:
        x_off = t["offset"]
        for g in t["genes"]:
            all_x.append(g["start_plot"] - x_off)
            all_x.append(g["end_plot"] - x_off)
    if not all_x:
        return 1.0

    plot_range = max(1, max(all_x) - min(all_x))
    coverages = []
    for t in all_tracks:
        if not t["genes"]:
            continue
        total_w = sum(g["end_plot"] - g["start_plot"] for g in t["genes"])
        coverages.append(total_w / plot_range)
    if not coverages:
        return 1.0

    coverages.sort()
    median_cov = coverages[len(coverages) // 2]
    if median_cov >= target_coverage:
        return 1.0
    factor = min(max_factor, target_coverage / max(median_cov, 1e-9))
    if factor <= 1.05:
        return 1.0

    for t in all_tracks:
        # Compressed-gap breaks bound each gene's segment: a gene must never be
        # widened ACROSS a break, or it would be drawn over the gap — implying
        # it is ~gap-bp longer than it really is. Clamp each widened gene to the
        # nearest break on either side of its centre (with a small margin so the
        # gap stays visible).
        bxs = sorted(b["x"] for b in t.get("breaks", []))
        for g in t["genes"]:
            os_, oe_ = g["start_plot"], g["end_plot"]   # original (pre-widen)
            c = (os_ + oe_) / 2.0
            half = (oe_ - os_) / 2.0 * factor
            ns, ne = c - half, c + half
            # Keep ≥ half of each adjacent gap, so a widened gene never reaches
            # (let alone crosses) the gap break beside it.
            left = max((bx for bx in bxs if bx <= os_), default=None)
            right = min((bx for bx in bxs if bx >= oe_), default=None)
            if left is not None:
                ns = max(ns, left + 0.5 * (os_ - left))
            if right is not None:
                ne = min(ne, right - 0.5 * (right - oe_))
            g["start_plot"], g["end_plot"] = ns, ne
    return factor


def _assign_sub_tracks(genes, x_off, min_gap=800):
    """Greedy interval scheduling: writes gene['_sub_track'] in-place.

    Sub-tracks are spent only when genes would visually overlap. No strand
    awareness — strand stays a property of each gene model, not of a row.
    """
    sorted_genes = sorted(genes, key=lambda g: g["start_plot"] - x_off)
    sub_ends = []  # rightmost x used by each sub-track so far
    for gene in sorted_genes:
        x0 = gene["start_plot"] - x_off
        x1 = gene["end_plot"]   - x_off
        placed = False
        for i, end_x in enumerate(sub_ends):
            if end_x + min_gap <= x0:
                gene["_sub_track"] = i
                sub_ends[i] = x1
                placed = True
                break
        if not placed:
            gene["_sub_track"] = len(sub_ends)
            sub_ends.append(x1)


# ======================================================================
# SVG Rendering Engine
# ======================================================================

def _svg_esc(text):
    """Escape text for safe embedding in SVG/HTML."""
    return _html_escape(str(text), quote=True)


def _resolve_css_vars(css):
    """Replace ``var(--name)`` / ``var(--name, fallback)`` with the literal
    values declared in the stylesheet's ``:root {}`` block.

    Many standalone-SVG renderers — cairosvg, librsvg, Illustrator, Inkscape's
    exporter — do not implement CSS custom properties, so an exported figure
    that keeps ``fill: var(--track-bg)`` renders those fills as *black*. (That
    is exactly why the static-view / publication SVGs looked like they had
    black target tracks outside a browser.) Resolving the variables to plain
    hex up front makes the static SVG render identically everywhere; browsers
    still receive valid CSS.
    """
    root = {}
    m = re.search(r":root\s*\{([^}]*)\}", css, re.DOTALL)
    if m:
        for decl in m.group(1).split(";"):
            if ":" in decl:
                k, v = decl.split(":", 1)
                k = k.strip()
                if k.startswith("--"):
                    root[k] = v.strip()
    if not root:
        return css

    var_re = re.compile(r"var\(\s*(--[\w-]+)\s*(?:,\s*([^()]*?))?\)")

    def _sub(expr):
        def repl(mm):
            val = root.get(mm.group(1).strip())
            if val is not None:
                return val
            return (mm.group(2) or "").strip() or "transparent"
        prev = None
        out = expr
        for _ in range(5):  # bounded passes resolve vars nested in vars
            out = var_re.sub(repl, out)
            if out == prev:
                break
            prev = out
        return out

    for k in list(root):
        root[k] = _sub(root[k])
    return _sub(css)


def _export_html_inline_svg(html_path, svg_path):
    """Write a standalone SVG that mirrors the interactive HTML's render.

    The interactive HTML embeds one inline `<svg>` plus a `<style>` block
    using CSS variables. Browsers resolve those when loading the HTML, but a
    raw SVG file has no such page context — without help, the standalone
    SVG renders mostly black because every `var(--…)` falls back to the
    default. This helper inlines the page CSS into the SVG inside CDATA so
    the resulting `.svg` is visually identical to the HTML view, drops
    cleanly into READMEs, and renders correctly in any SVG viewer.
    """
    with open(html_path) as fh:
        src = fh.read()
    style_blocks = re.findall(r'<style[^>]*>(.*?)</style>', src, re.DOTALL)
    # Resolve CSS custom properties to literal colours so the standalone SVG
    # renders correctly in viewers that don't support var() (cairosvg, etc.).
    combined_style = _resolve_css_vars("\n".join(style_blocks))
    m = re.search(r'(<svg[^>]*>)(.*?)(</svg>)', src, re.DOTALL)
    if not m:
        raise RuntimeError("no <svg>...</svg> element found in HTML output")
    svg_open, svg_inner, svg_close = m.group(1), m.group(2), m.group(3)
    if 'xmlns="http://www.w3.org/2000/svg"' not in svg_open:
        svg_open = svg_open.replace(
            '<svg', '<svg xmlns="http://www.w3.org/2000/svg"', 1
        )
    body = (svg_open
            + '\n<style type="text/css"><![CDATA[\n'
            + combined_style
            + '\n]]></style>\n'
            + svg_inner + svg_close)
    with open(svg_path, "w") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + body)


def _format_species_label(name):
    """Apply --common_names mode to a raw species string.

    Reads two module globals set by main() at startup:
      _synvoy_taxa       — sibling helper module (or None if init failed)
      _common_name_mode  — one of 'both', 'common', 'scientific', 'off'
    Falls back to the raw name unchanged when lookup is disabled or the
    helper is unavailable.
    """
    helper = globals().get("_synvoy_taxa")
    mode = globals().get("_common_name_mode", "off")
    if helper is None or mode == "off":
        return name.replace("_", " ") if name else ""
    try:
        return helper.label_for_species(name.replace("_", " "), mode=mode)
    except Exception:
        return name.replace("_", " ") if name else ""


def _svg_arrow_path(x0, x1, yb, h, strand, rx=3):
    """Generate SVG path 'd' attribute for a pentagon gene arrow with rounded back."""
    w = x1 - x0
    aw = min(w * 0.35, h * 0.8)
    if aw < 1:
        aw = min(w * 0.5, 1)
    ym = yb + h / 2
    yt = yb + h
    rx = min(rx, w * 0.2, h * 0.2)

    if strand == "+":
        # Rounded left edge, pointed right
        return (
            f"M{x0 + rx:.1f},{yb:.1f} "
            f"L{x1 - aw:.1f},{yb:.1f} "
            f"L{x1:.1f},{ym:.1f} "
            f"L{x1 - aw:.1f},{yt:.1f} "
            f"L{x0 + rx:.1f},{yt:.1f} "
            f"Q{x0:.1f},{yt:.1f} {x0:.1f},{yt - rx:.1f} "
            f"L{x0:.1f},{yb + rx:.1f} "
            f"Q{x0:.1f},{yb:.1f} {x0 + rx:.1f},{yb:.1f} Z"
        )
    else:
        # Pointed left, rounded right edge
        return (
            f"M{x0:.1f},{ym:.1f} "
            f"L{x0 + aw:.1f},{yb:.1f} "
            f"L{x1 - rx:.1f},{yb:.1f} "
            f"Q{x1:.1f},{yb:.1f} {x1:.1f},{yb + rx:.1f} "
            f"L{x1:.1f},{yt - rx:.1f} "
            f"Q{x1:.1f},{yt:.1f} {x1 - rx:.1f},{yt:.1f} "
            f"L{x0 + aw:.1f},{yt:.1f} Z"
        )


def _direction_chevrons(x0, x1, y, strand, colour, gene_h):
    """IGV-style direction chevrons along a multi-exon gene backbone.

    Multi-exon genes are drawn as exon blocks on a thin intron line; when the
    terminal exon is narrow its arrow-tip disappears, so the strand becomes
    invisible. A few small chevrons on the backbone make the coding direction
    obvious at any exon count. Returns a list of SVG path strings.
    """
    w = x1 - x0
    if w < 22:
        return []
    ch = min(4.0, gene_h * 0.18)
    n = max(1, int(w // 26.0))
    parts = []
    for k in range(n):
        cx = x0 + (k + 0.5) * (w / n)
        if strand == "-":
            d = f"M{cx + ch:.1f},{y - ch:.1f} L{cx - ch:.1f},{y:.1f} L{cx + ch:.1f},{y + ch:.1f}"
        else:
            d = f"M{cx - ch:.1f},{y - ch:.1f} L{cx + ch:.1f},{y:.1f} L{cx - ch:.1f},{y + ch:.1f}"
        parts.append(
            f'<path d="{d}" fill="none" stroke="{colour}" stroke-width="1.3" '
            f'stroke-linecap="round" stroke-linejoin="round" opacity="0.85" '
            f'pointer-events="none"/>'
        )
    return parts


def _svg_ribbon_path(ux0, ux1, uy_bot, lx0, lx1, ly_top):
    """Generate SVG path 'd' for a bezier-curve synteny ribbon."""
    cy = (uy_bot + ly_top) / 2
    return (
        f"M{ux0:.1f},{uy_bot:.1f} "
        f"C{ux0:.1f},{cy:.1f} {lx0:.1f},{cy:.1f} {lx0:.1f},{ly_top:.1f} "
        f"L{lx1:.1f},{ly_top:.1f} "
        f"C{lx1:.1f},{cy:.1f} {ux1:.1f},{cy:.1f} {ux1:.1f},{uy_bot:.1f} Z"
    )


def _build_tooltip_json(gene, track, home_products):
    """Build tooltip data dict for a gene, returned as escaped JSON string."""
    is_home = track.get("is_home", False)
    name = gene.get("name", "")
    home_id = gene.get("home_gene_id", name)
    goi_f = _is_goi_target_gene(gene) if not is_home else (is_goi(name) or is_goi(home_id))
    resolved = _is_resolved_goi_target_gene(gene) if not is_home else goi_f
    ambiguous = goi_f and not resolved if not is_home else False

    if is_home:
        cn = _preferred_home_label(gene, home_products)
        product = _lookup_product(name, home_products)
    else:
        cn = clean_gene_label(_preferred_target_label(gene))
        product = gene.get("target_product", "")

    n_ex = len(gene.get("exon_coords", []))
    if n_ex <= 1:
        n_ex = gene.get("n_exons", 0)

    data = {
        "name": cn,
        "product": product or "",
        "coords": f"{gene['chrom']}:{gene['start']:,}-{gene['end']:,}",
        "strand": gene.get("strand", "+"),
        "exons": n_ex if n_ex and n_ex > 1 else 0,
        "isHome": is_home,
    }

    if not is_home:
        if home_id:
            data["homolog"] = clean_gene_label(home_id)
        if "identity" in gene:
            data["identity"] = round(gene["identity"], 1)
        conf = (gene.get("confidence") or "").upper()
        if conf:
            data["confidence"] = conf
        gc = gene.get("goi_class", "")
        if gc:
            data["goiClass"] = gc.replace("_", " ")
        et = gene.get("evidence_type", "")
        if et:
            data["evidence"] = et.replace("_", " ")
        ms = gene.get("model_status", "")
        if ms:
            data["model"] = ms
        sc = gene.get("synteny_context", "")
        if sc:
            data["synteny"] = sc.replace("_", " ")
        qc = gene.get("query_coverage")
        if qc is not None:
            data["queryCov"] = round(qc * 100, 1)

    if resolved:
        data["goiTag"] = "GOI"
    elif ambiguous:
        data["goiTag"] = "GOI-like / ambiguous"
    elif goi_f and is_home:
        data["goiTag"] = "GOI"

    return _svg_esc(json.dumps(data, ensure_ascii=True))


def _gene_display_label(gene, track, home_products, goi_f, resolved_goi_f):
    """Return the human-readable label string for a gene, regardless of
    whether the label is currently emitted as on-canvas text. Used both
    when setting the gene-group `data-label` (for click-to-pin in JS)
    and when emitting the default label layer."""
    name = gene["name"]
    home_id = gene.get("home_gene_id", name)
    if not track["is_home"]:
        if goi_f:
            label = clean_gene_label(name, keep_goi_prefix=True)
            if not label or label == 'GOI':
                label = clean_gene_label(home_id, keep_goi_prefix=True)
            if not resolved_goi_f:
                label = "~ " + label
        else:
            label = clean_gene_label(_preferred_target_label(gene))
            if not label and home_id:
                label = clean_gene_label(home_id)
    else:
        label = _preferred_home_label(gene, home_products)
    return label or name


def render_synteny_html(all_tracks, gene_colours, goi_genome_colours,
                        home_products, args,
                        subtitle_bits, hidden_absent_tracks,
                        ambiguous_track_count, resolved_track_count,
                        force_home_labels=False):
    """Render synteny visualization as self-contained HTML with embedded SVG.

    `force_home_labels=True` adds an on-canvas label for every gene in the
    home track (used by the publication-SVG export). The interactive HTML
    leaves it False so the canvas stays clean — flanking gene names live in
    tooltips and click-to-pin labels."""

    n_tracks = len(all_tracks)

    # ---- Layout constants ----
    GENE_H        = 30
    SUB_TRACK_GAP = 10
    TRACK_MARGIN  = 85
    # LEFT_MARGIN will be computed dynamically below based on longest track label
    LEFT_MARGIN   = 220
    RIGHT_MARGIN  = 50
    # Header band must clear the title + the (long) descriptive subtitle so
    # neither the subtitle nor the centred "GOI" guide tag lands on the home
    # panel. The first track starts at TOP_MARGIN.
    TOP_MARGIN    = 104
    BOTTOM_MARGIN = 90  # More room for legend
    TRACK_PAD     = 10
    MIN_GENE_PX   = 4
    EXON_RX       = 3
    RIBBON_GAP    = min(8, TRACK_MARGIN * 0.1)

    # ---- Compute x range ----
    all_x_bp = []
    goi_x_bp = []      # centered-bp extents of GOI/toxin genes (for local zoom)
    for track in all_tracks:
        x_off = track["offset"]
        for g in track["genes"]:
            s_c = g["start_plot"] - x_off
            e_c = g["end_plot"] - x_off
            all_x_bp.append(s_c)
            all_x_bp.append(e_c)
            is_goi_g = (_is_goi_target_gene(g) if not track["is_home"]
                        else (is_goi(g.get("name")) or is_goi(g.get("home_gene_id", ""))))
            if is_goi_g:
                goi_x_bp.extend((s_c, e_c))

    # ---- Local GOI-window magnification (piecewise-linear 'broken' axis) ----
    # The genomic scale (px/bp) is multiplied by `goi_zoom` inside the band
    # spanned by GOI/toxin genes and left at 1x outside, so tandem toxin arrays
    # spread apart (stop overlapping) without stretching the flanking regions.
    # Every coordinate is routed through `_stretch` *before* the linear px map,
    # so genes, ribbons, breaks and the guide line all stay consistent.
    _zoom = max(1.0, float(getattr(args, "goi_zoom", 1.0) or 1.0))
    if _zoom > 1.0 and goi_x_bp:
        _gpad = 800.0  # bp of context kept at zoom on each side of the toxin band
        _w_lo = min(goi_x_bp) - _gpad
        _w_hi = max(goi_x_bp) + _gpad
        _w_span = max(1.0, _w_hi - _w_lo)

        def _stretch(b):
            if b <= _w_lo:
                return b
            if b >= _w_hi:
                return _w_lo + _w_span * _zoom + (b - _w_hi)
            return _w_lo + (b - _w_lo) * _zoom
    else:
        def _stretch(b):
            return b

    if not all_x_bp:
        x_min_bp, x_max_bp = -1000, 1000
        raw_range = 2000.0
    else:
        _sx = [_stretch(b) for b in all_x_bp]
        x_min_bp, x_max_bp = min(_sx), max(_sx)
        raw_range = max(1.0, max(all_x_bp) - min(all_x_bp))

    pad_bp = (x_max_bp - x_min_bp) * 0.05 + 5000
    x_min_bp -= pad_bp
    x_max_bp += pad_bp
    raw_range += 2 * pad_bp  # keep raw/stretched on the same padded footing

    # ---- Plot dimensions ----
    # Dynamic left margin: measure longest track label and allow more room
    try:
        max_lbl = 0
        for t in all_tracks:
            lbl = re.sub(r"<[^>]+>", "", t.get("label", "") or "")
            max_lbl = max(max_lbl, len(lbl))
        # approx char * px + padding
        LEFT_MARGIN = max(260, int(max_lbl * 8 + 60))
    except Exception:
        LEFT_MARGIN = 260
    if args.plot_width > 0:
        plot_w = max(800, args.plot_width)
    else:
        est = max(1200, int((x_max_bp - x_min_bp) / 350))
        plot_w = min(6000, est)

    bp_range = max(1, x_max_bp - x_min_bp)
    if _zoom > 1.0 and goi_x_bp:
        # Local magnification: anchor px/bp to the *flanking* (raw) scale so the
        # flanking regions keep their normal density, and grow the canvas to fit
        # the widened toxin band — the plot gets a bit wider rather than squashing
        # the flanking genes.
        available_w = plot_w - LEFT_MARGIN - RIGHT_MARGIN
        scale = available_w / raw_range            # px per bp at flanking scale
        plot_w_needed = LEFT_MARGIN + RIGHT_MARGIN + bp_range * scale
        if plot_w_needed > 9000:                    # hard cap; shrink scale to fit
            plot_w = 9000
            scale = (plot_w - LEFT_MARGIN - RIGHT_MARGIN) / bp_range
        else:
            plot_w = int(round(plot_w_needed))
        available_w = plot_w - LEFT_MARGIN - RIGHT_MARGIN
    else:
        available_w = plot_w - LEFT_MARGIN - RIGHT_MARGIN
        scale = available_w / bp_range  # px per bp

    # ---- Track heights & y positions ----
    track_heights = []
    for track in all_tracks:
        genes = track["genes"]
        n_sub = (max(g.get("_sub_track", 0) for g in genes) + 1) if genes else 1
        th = n_sub * GENE_H + max(0, n_sub - 1) * SUB_TRACK_GAP
        track_heights.append(th)

    track_y = []
    y_cursor = TOP_MARGIN
    for i in range(n_tracks):
        track_y.append(y_cursor)
        y_cursor += track_heights[i] + TRACK_MARGIN

    total_h = y_cursor - TRACK_MARGIN + BOTTOM_MARGIN
    if args.plot_height > 0:
        total_h = max(total_h, args.plot_height)

    # Optional caption block pinned below the plot content (e.g. conserved-locus
    # / no-toxin / off-contig notes). Reserve extra canvas height for it.
    _caption_lines = []
    if getattr(args, "caption_file", None) and os.path.exists(args.caption_file):
        _caption_lines = [ln.rstrip("\n") for ln in open(args.caption_file, encoding="utf-8")
                          if ln.strip()]
    _content_bottom = total_h
    total_h += (len(_caption_lines) * 15 + 26) if _caption_lines else 0

    # ---- Coordinate helpers (closures) ----
    def bp2px(bp_val):
        # bp_val is centered bp (gene start_plot/end_plot minus track offset);
        # apply the local GOI magnification, then the linear px map.
        return LEFT_MARGIN + (_stretch(bp_val) - x_min_bp) * scale

    def gene_px(gene, track):
        x_off = track["offset"]
        x0 = bp2px(gene["start_plot"] - x_off)
        x1 = bp2px(gene["end_plot"] - x_off)
        # GOI/toxin genes are tiny (~200 bp) at locus scale; give them a larger floor so
        # the exon/intron model (drawn relative to the box, needs w_px > 25) renders.
        floor = MIN_GENE_PX
        if getattr(args, "goi_min_px", 0) and _is_goi_target_gene(gene):
            floor = max(MIN_GENE_PX, args.goi_min_px)
        if x1 - x0 < floor:
            mid = (x0 + x1) / 2
            x0, x1 = mid - floor / 2, mid + floor / 2
        return x0, x1

    def gene_yb(ti, gene):
        sub = gene.get("_sub_track", 0)
        return track_y[ti] + sub * (GENE_H + SUB_TRACK_GAP)

    # ---- Build SVG elements ----
    svg_parts = []

    # ---- SVG defs (filters) ----
    svg_parts.append('<defs>')
    svg_parts.append(f"""
    <filter id="geneShadow" x="-4%" y="-15%" width="108%" height="140%">
      <feDropShadow dx="0" dy="1" stdDeviation="1.5" flood-opacity="0.10" flood-color="#000"/>
    </filter>
    <linearGradient id="geneGloss" x1="0%" y1="0%" x2="0%" y2="100%">
      <stop offset="0%" stop-color="#ffffff" stop-opacity="0.30"/>
      <stop offset="45%" stop-color="#ffffff" stop-opacity="0.05"/>
      <stop offset="55%" stop-color="#000000" stop-opacity="0.05"/>
      <stop offset="100%" stop-color="#000000" stop-opacity="0.25"/>
    </linearGradient>
    <!-- Diagonal-stripe pattern for ambiguous GOI / low-confidence rescue
         hits. Same idiom as plot_synteny_matrix.py so the two plots read alike. -->
    <pattern id="ambiguousGoi" patternUnits="userSpaceOnUse"
             width="7" height="7" patternTransform="rotate(45)">
      <rect width="7" height="7" fill="#fff0ec"/>
      <rect width="3" height="7" fill="{GOI_PUB_COLOUR}" opacity="0.75"/>
    </pattern>
    """)
    svg_parts.append('</defs>')

    # ---- Canvas background ----
    # Explicit white rect so the exported static SVG/PNG isn't transparent
    # (renders black) in viewers that ignore the CSS `background` property.
    svg_parts.append(
        f'<rect x="0" y="0" width="{plot_w}" height="{total_h}" fill="#ffffff"/>'
    )

    # ---- Track backgrounds ----
    for ti, track in enumerate(all_tracks):
        if not track["genes"]:
            continue
        yb = track_y[ti]
        th = track_heights[ti]
        x_off = track["offset"]

        gxs = []
        for g in track["genes"]:
            x0, x1 = gene_px(g, track)
            gxs.extend([x0, x1])
        if not gxs:
            continue

        x_left = min(gxs) - TRACK_PAD
        x_right = max(gxs) + TRACK_PAD

        chrom_break_xs = sorted(
            bp2px(brk["x"] - x_off)
            for brk in track.get("breaks", [])
            if brk.get("is_chrom_break")
        )

        if not chrom_break_xs:
            svg_parts.append(
                f'<rect x="{x_left:.1f}" y="{yb - TRACK_PAD:.1f}" '
                f'width="{x_right - x_left:.1f}" height="{th + 2*TRACK_PAD:.1f}" '
                f'class="track-bg track-item" data-track-idx="{ti}" rx="6"/>'
            )
        else:
            boundaries = [x_left] + chrom_break_xs + [x_right]
            for si in range(len(boundaries) - 1):
                inset = 10
                seg_x0 = boundaries[si] + (inset if si > 0 else 0)
                seg_x1 = boundaries[si + 1] - (inset if si < len(boundaries) - 2 else 0)
                if seg_x0 < seg_x1:
                    svg_parts.append(
                        f'<rect x="{seg_x0:.1f}" y="{yb - TRACK_PAD:.1f}" '
                        f'width="{seg_x1 - seg_x0:.1f}" height="{th + 2*TRACK_PAD:.1f}" '
                        f'class="track-bg track-item" data-track-idx="{ti}" rx="6"/>'
                    )

    # ---- GOI guide line ----
    # Every track is centred on its GOI (offset = anchor centre), so the GOI
    # sits at bp≈0 in every track and a single vertical guide aligns them all.
    # A faint dashed red line + a small tag makes the GOI instantly findable
    # down the stacked tracks.
    if all_tracks and all_tracks[0]["genes"]:
        x_guide = bp2px(0)
        y_top = track_y[0] - TRACK_PAD - 14
        y_bot = track_y[-1] + track_heights[-1] + TRACK_PAD
        svg_parts.append(
            f'<line class="goi-guide-line" x1="{x_guide:.1f}" y1="{y_top:.1f}" '
            f'x2="{x_guide:.1f}" y2="{y_bot:.1f}" stroke="{GOI_COLOUR}" '
            f'stroke-width="1" stroke-dasharray="4,4" opacity="0.45"/>'
        )
        svg_parts.append(
            f'<text class="goi-guide-tag" x="{x_guide:.1f}" y="{y_top - 2:.1f}" '
            f'text-anchor="middle" font-size="10" font-weight="700" '
            f'fill="{GOI_BORDER}" opacity="0.8">GOI</text>'
        )

    # ---- Ribbons (drawn first, behind genes) ----
    # ONE ribbon per orthology — not a quadratic fan-out. The GOI matches
    # GOI-to-GOI, and each genome can carry up to ~10 GOI copies, so connecting
    # every copy to every copy piled ~100 near-opaque ribbons per track pair
    # into a dark, muddy mass (the "black ribbons"). Instead draw at most ONE
    # GOI ribbon per track pair (best-resolved ↔ best-resolved) plus one
    # flanking ribbon per lower gene to its best upper match.
    def _gene_is_goi_in(track, g):
        if track["is_home"]:
            return is_goi(g.get("name")) or is_goi(g.get("home_gene_id", ""))
        return _is_goi_target_gene(g)

    def _best_goi_in(track):
        gois = [g for g in track["genes"] if _gene_is_goi_in(track, g)]
        if not gois:
            return None
        def _key(g):
            resolved = True if track["is_home"] else _is_resolved_goi_target_gene(g)
            return (1 if resolved else 0, g.get("identity", 0.0))
        return max(gois, key=_key)

    # Why two groups + group-opacity: every track is GOI-centred, so many
    # flanking ribbons funnel toward the same x. With per-ribbon alpha, those
    # overlaps COMPOUND — a dozen translucent ribbons of different hues stack
    # into a muddy near-black mass (the "black ribbons"). Drawing each flanking
    # ribbon at FULL opacity inside ONE group that carries a single `opacity`
    # makes the layer composite once: 10 overlapping ribbons look identical to
    # 1, so overlaps never darken. The GOI ribbon sits in its own group on top
    # and is painted with a vertical GRADIENT that flows from the upper gene's
    # colour into the lower gene's colour (so mismatched clade colours blend
    # *smoothly* instead of muddying).
    flank_parts = []
    goi_parts = []
    grad_defs = []
    for ti in range(len(all_tracks) - 1):
        upper = all_tracks[ti]
        lower = all_tracks[ti + 1]

        def _ribbon_geom(ug, lg, _ti=ti, _upper=upper, _lower=lower):
            y_ug_bot = gene_yb(_ti, ug) + GENE_H + RIBBON_GAP
            y_lg_top = gene_yb(_ti + 1, lg) - RIBBON_GAP
            ux0, ux1 = gene_px(ug, _upper)
            lx0, lx1 = gene_px(lg, _lower)
            return _svg_ribbon_path(ux0, ux1, y_ug_bot, lx0, lx1, y_lg_top), y_ug_bot, y_lg_top

        # Best upper FLANKING gene per home key (by identity).
        upper_by_key = {}
        for ug in upper["genes"]:
            if _gene_is_goi_in(upper, ug):
                continue
            key = ug.get("home_gene_id") or ug.get("name")
            if not key:
                continue
            cur = upper_by_key.get(key)
            if cur is None or ug.get("identity", 0) > cur.get("identity", 0):
                upper_by_key[key] = ug

        # Flanking ribbons: one per lower flanking gene -> best upper match.
        # Solid, lightened fill at FULL opacity; the group opacity below makes
        # the whole flanking layer translucent without per-overlap darkening.
        for lg in lower["genes"]:
            if _gene_is_goi_in(lower, lg):
                continue
            home_id = lg.get("home_gene_id", "")
            if not home_id:
                continue
            ug = upper_by_key.get(home_id) or upper_by_key.get(lg.get("name"))
            if ug is None:
                continue
            colour = gene_colours.get(home_id,
                     gene_colours.get(ug.get("name"), UNMATCHED_CLR))
            path_d, _, _ = _ribbon_geom(ug, lg)
            flank_parts.append(
                f'<path d="{path_d}" fill="{_lerp_hex(colour, "#ffffff", 0.28)}" '
                f'stroke="none" class="ribbon" data-homology="{_svg_esc(home_id)}" '
                f'data-upper-track="{ti}" data-lower-track="{ti + 1}"/>'
            )

        # Single GOI ribbon: best-resolved upper <-> best-resolved lower, drawn
        # with a vertical gradient (upper colour -> lower colour).
        ug_goi = _best_goi_in(upper)
        lg_goi = _best_goi_in(lower)
        if ug_goi is not None and lg_goi is not None:
            up_col = _goi_colour_for_genome(upper["genome_id"], goi_genome_colours)
            lo_col = _goi_colour_for_genome(lower["genome_id"], goi_genome_colours)
            path_d, y_ug_bot, y_lg_top = _ribbon_geom(ug_goi, lg_goi)
            gid = f"goiRib{ti}"
            grad_defs.append(
                f'<linearGradient id="{gid}" gradientUnits="userSpaceOnUse" '
                f'x1="0" y1="{y_ug_bot:.1f}" x2="0" y2="{y_lg_top:.1f}">'
                f'<stop offset="0" stop-color="{up_col}"/>'
                f'<stop offset="1" stop-color="{lo_col}"/></linearGradient>'
            )
            op = 0.62
            if not lower["is_home"] and not _is_resolved_goi_target_gene(lg_goi):
                op = 0.32
            goi_parts.append(
                f'<path d="{path_d}" fill="url(#{gid})" opacity="{op:.2f}" '
                f'stroke="{_hex_to_rgba(lo_col, 0.55)}" stroke-width="0.6" '
                f'class="ribbon" data-homology="GOI" '
                f'data-upper-track="{ti}" data-lower-track="{ti + 1}"/>'
            )

    if grad_defs:
        svg_parts.append('<defs>' + "".join(grad_defs) + '</defs>')
    # Flanking layer: full-opacity ribbons, single group opacity (no compounding).
    svg_parts.append('<g class="ribbons ribbons-flank" opacity="0.5">')
    svg_parts.extend(flank_parts)
    svg_parts.append('</g>')
    # GOI layer on top.
    svg_parts.append('<g class="ribbons ribbons-goi">')
    svg_parts.extend(goi_parts)
    svg_parts.append('</g>')

    # ---- Gene models ----
    svg_parts.append('<g class="genes">')
    legend_shown = set()

    for ti, track in enumerate(all_tracks):
        x_off = track["offset"]
        # Draw large genes first so small genes render on top
        sorted_genes = sorted(track["genes"],
                              key=lambda g: g["end_plot"] - g["start_plot"],
                              reverse=True)

        for gene in sorted_genes:
            yb = gene_yb(ti, gene)
            name = gene["name"]
            home_id = gene.get("home_gene_id", name)
            goi_like = _is_goi_target_gene(gene) if not track["is_home"] else (is_goi(name) or is_goi(home_id))
            resolved_goi = _is_resolved_goi_target_gene(gene) if not track["is_home"] else goi_like
            ambiguous_goi = goi_like and not resolved_goi if not track["is_home"] else False
            confidence = (gene.get("confidence") or "").upper()

            # --- colour ---
            if resolved_goi or (track["is_home"] and goi_like):
                colour = _goi_colour_for_genome(track["genome_id"], goi_genome_colours)
                bclr = GOI_BORDER
                bw = 2.2 if confidence == "HIGH" or track["is_home"] else 1.8
                dash = ""
            elif ambiguous_goi:
                # Ambiguous GOI: diagonal-stripe pattern fill + dashed border
                # to make 'low-confidence rescue hit' visible at a glance,
                # rather than just a slightly-different border weight on a
                # solid red block. Matches the matrix plot's idiom.
                colour = "url(#ambiguousGoi)"
                bclr = GOI_BORDER
                bw = 1.5
                dash = ' stroke-dasharray="6,3"'
            elif home_id in gene_colours:
                colour = gene_colours[home_id]
                bclr = _darken_hex(colour, 0.6)
                bw = 1
                dash = ' stroke-dasharray="2,2"' if confidence == "LOW" else ""
            elif name in gene_colours:
                colour = gene_colours[name]
                bclr = _darken_hex(colour, 0.6)
                bw = 1
                dash = ' stroke-dasharray="2,2"' if confidence == "LOW" else ""
            else:
                colour = UNMATCHED_CLR
                bclr = "#b0b0b0"
                bw = 1.0
                dash = ' stroke-dasharray="2,2"' if confidence == "LOW" else ""

            x0, x1 = gene_px(gene, track)
            w_px = x1 - x0
            strand = gene.get("strand", "+")

            tooltip_json = _build_tooltip_json(gene, track, home_products)
            goi_attr = ' data-is-goi="true"' if goi_like else ''
            hom_id_attr = _svg_esc(home_id)
            x0_attr = f'{x0:.1f}'
            x1_attr = f'{x1:.1f}'
            yb_attr = f'{yb:.1f}'

            # Pre-compute the gene's display label so click-to-pin in JS can
            # look it up without redoing the (track-aware) label resolution.
            display_label = _gene_display_label(
                gene, track, home_products, goi_like, resolved_goi
            )
            label_attr = _svg_esc(display_label)

            svg_parts.append(
                f'<g class="gene-group" data-homology="{hom_id_attr}" '
                f'data-track="{ti}" data-track-idx="{ti}" '
                f'data-x0="{x0_attr}" data-x1="{x1_attr}" data-yb="{yb_attr}" '
                f'data-fill="{colour}" data-identity="{gene.get("identity", 0.0)}" '
                f'data-label="{label_attr}" '
                f'data-tooltip=\'{tooltip_json}\'{goi_attr}>'
            )

            # Render gene body
            exon_coords = gene.get("exon_coords", [])
            n_exons_attr = gene.get("n_exons", 0)
            has_real_exons = len(exon_coords) >= 2 and w_px > 25
            has_synth_exons = (not has_real_exons and n_exons_attr
                               and n_exons_attr >= 2 and w_px > 25)

            if has_real_exons or has_synth_exons:
                # --- Exon/intron model ---
                mid_y = yb + GENE_H / 2

                # Intron backbone line
                svg_parts.append(
                    f'<line x1="{x0:.1f}" y1="{mid_y:.1f}" '
                    f'x2="{x1:.1f}" y2="{mid_y:.1f}" '
                    f'class="intron-line"{dash}/>'
                )

                if has_real_exons:
                    gene_s = gene["start"]
                    gene_e = gene["end"]
                    gene_span = max(1, gene_e - gene_s)
                    for ei, (es, ee) in enumerate(exon_coords):
                        frac_s = max(0, min(1, (es - gene_s) / gene_span))
                        frac_e = max(0, min(1, (ee - gene_s) / gene_span))
                        ex0 = x0 + frac_s * w_px
                        ex1 = x0 + frac_e * w_px
                        ew = max(2, ex1 - ex0)

                        # Last/first exon gets arrow tip
                        is_terminal = ((strand == "+" and ei == len(exon_coords) - 1) or
                                       (strand == "-" and ei == 0))
                        if is_terminal and ew > 10:
                            aw = min(ew * 0.3, GENE_H * 0.5)
                            if strand == "+":
                                d = (f"M{ex0:.1f},{yb:.1f} L{ex0 + ew - aw:.1f},{yb:.1f} "
                                     f"L{ex0 + ew:.1f},{mid_y:.1f} L{ex0 + ew - aw:.1f},{yb + GENE_H:.1f} "
                                     f"L{ex0:.1f},{yb + GENE_H:.1f} Z")
                            else:
                                d = (f"M{ex0:.1f},{mid_y:.1f} L{ex0 + aw:.1f},{yb:.1f} "
                                     f"L{ex0 + ew:.1f},{yb:.1f} L{ex0 + ew:.1f},{yb + GENE_H:.1f} "
                                     f"L{ex0 + aw:.1f},{yb + GENE_H:.1f} Z")
                            svg_parts.append(
                                f'<path d="{d}" fill="{colour}" stroke="{bclr}" '
                                f'stroke-width="{bw}" class="exon"{dash}/>'
                            )
                        else:
                            svg_parts.append(
                                f'<rect x="{ex0:.1f}" y="{yb:.1f}" width="{ew:.1f}" '
                                f'height="{GENE_H}" rx="{EXON_RX}" fill="{colour}" '
                                f'stroke="{bclr}" stroke-width="{bw}" class="exon"{dash}/>'
                            )
                else:
                    # Synthesized evenly-spaced exons
                    for k in range(n_exons_attr):
                        frac_s = k / n_exons_attr
                        frac_e = (k + 0.65) / n_exons_attr
                        ex0 = x0 + frac_s * w_px
                        ex1 = x0 + frac_e * w_px
                        ew = max(2, ex1 - ex0)

                        is_terminal = ((strand == "+" and k == n_exons_attr - 1) or
                                       (strand == "-" and k == 0))
                        if is_terminal and ew > 10:
                            aw = min(ew * 0.3, GENE_H * 0.5)
                            if strand == "+":
                                d = (f"M{ex0:.1f},{yb:.1f} L{ex0 + ew - aw:.1f},{yb:.1f} "
                                     f"L{ex0 + ew:.1f},{mid_y:.1f} L{ex0 + ew - aw:.1f},{yb + GENE_H:.1f} "
                                     f"L{ex0:.1f},{yb + GENE_H:.1f} Z")
                            else:
                                d = (f"M{ex0:.1f},{mid_y:.1f} L{ex0 + aw:.1f},{yb:.1f} "
                                     f"L{ex0 + ew:.1f},{yb:.1f} L{ex0 + ew:.1f},{yb + GENE_H:.1f} "
                                     f"L{ex0 + aw:.1f},{yb + GENE_H:.1f} Z")
                            svg_parts.append(
                                f'<path d="{d}" fill="{colour}" stroke="{bclr}" '
                                f'stroke-width="{bw}" class="exon"{dash}/>'
                            )
                        else:
                            svg_parts.append(
                                f'<rect x="{ex0:.1f}" y="{yb:.1f}" width="{ew:.1f}" '
                                f'height="{GENE_H}" rx="{EXON_RX}" fill="{colour}" '
                                f'stroke="{bclr}" stroke-width="{bw}" class="exon"{dash}/>'
                            )

                # Direction chevrons on the backbone: the strand stays obvious
                # even when the terminal exon is too narrow to show an arrow-tip.
                _chev_clr = "#475569" if str(colour).startswith("url") else _darken_hex(colour, 0.5)
                svg_parts.extend(
                    _direction_chevrons(x0, x1, mid_y, strand, _chev_clr, GENE_H)
                )
            else:
                # --- Single-block arrow gene (pentagon points along strand) ---
                path_d = _svg_arrow_path(x0, x1, yb, GENE_H, strand, rx=EXON_RX)
                svg_parts.append(
                    f'<path d="{path_d}" fill="{colour}" stroke="{bclr}" '
                    f'stroke-width="{bw}" class="exon"{dash}/>'
                )
                svg_parts.append(
                    f'<path d="{path_d}" fill="url(#geneGloss)" pointer-events="none"/>'
                )

            svg_parts.append('</g>')
    svg_parts.append('</g>')

    # ---- Absent-GOI placeholders ----
    for ti, track in enumerate(all_tracks):
        if track["is_home"]:
            continue
        if track.get("goi_status") == "absent" and track["genes"]:
            yb = track_y[ti]
            cx = bp2px(0)
            dash_w = max(MIN_GENE_PX, 2000 * scale)
            svg_parts.append(
                f'<rect x="{cx - dash_w/2:.1f}" y="{yb:.1f}" '
                f'width="{dash_w:.1f}" height="{GENE_H}" rx="4" '
                f'fill="rgba(227,26,28,0.06)" stroke="{GOI_COLOUR}" '
                f'stroke-width="1.5" stroke-dasharray="6,3"/>'
            )
            svg_parts.append(
                f'<text x="{cx:.1f}" y="{yb + GENE_H/2 + 5:.1f}" '
                f'text-anchor="middle" fill="{GOI_COLOUR}" font-size="14" '
                f'font-weight="700">?</text>'
            )

    # ---- Gene labels ----
    # Interactive HTML keeps the canvas clean: only GOIs are labelled by
    # default. Flanking-gene names live in the hover tooltip and on the
    # click-to-pin layer. Publication SVG (`force_home_labels=True`) adds
    # one label per home-track gene.
    svg_parts.append('<g class="gene-labels">')
    for ti, track in enumerate(all_tracks):
        x_off = track["offset"]
        n_genes = len(track["genes"])
        fsize = max(9, 13 - (n_genes // 6))
        is_home = track["is_home"]

        # Separate GOI hits from forced home-flanking labels. A divergent or
        # tandem GOI family can put many copies in one track; stacking one
        # rotated label per copy produced an unreadable pile (e.g. five
        # "★ ~ Melt" on top of each other). Instead label only the single best
        # copy (resolved > ambiguous, then highest identity) and tag it with
        # the copy count "×N" — the other copies stay drawn, just unlabelled,
        # with their details in the hover tooltip / click-to-pin layer.
        goi_genes, flank_force = [], []
        for gene in track["genes"]:
            name = gene["name"]
            home_id = gene.get("home_gene_id", name)
            goi_f = (_is_goi_target_gene(gene) if not is_home
                     else (is_goi(name) or is_goi(home_id)))
            if goi_f:
                resolved = (_is_resolved_goi_target_gene(gene) if not is_home
                            else True)
                goi_genes.append((gene, resolved))
            elif force_home_labels and is_home:
                flank_force.append(gene)

        label_jobs = []  # (gene, label_text, is_goi)
        if goi_genes:
            goi_genes.sort(key=lambda gr: (gr[1], gr[0].get("identity", 0.0)),
                           reverse=True)
            best_gene, best_resolved = goi_genes[0]
            label = _gene_display_label(best_gene, track, home_products,
                                        True, best_resolved)
            if len(goi_genes) > 1:
                label = f"{label} ×{len(goi_genes)}"
            label_jobs.append((best_gene, "★ " + label, True))
        for gene in flank_force:
            label_jobs.append(
                (gene, _gene_display_label(gene, track, home_products, False, False),
                 False))

        for gene, label, is_goi_lbl in label_jobs:
            g_start, g_end = _get_coords(gene)
            xc_px = bp2px((g_start + g_end) / 2 - x_off)
            yb_px = gene_yb(ti, gene)
            lbl_y = yb_px - 4
            lbl_class = "gene-label goi" if is_goi_lbl else "gene-label"
            # White halo (paint-order:stroke) so the label reads over ribbons
            # and adjacent gene arrows instead of disappearing into them.
            svg_parts.append(
                f'<text x="{xc_px:.1f}" y="{lbl_y:.1f}" '
                f'transform="rotate(-45 {xc_px:.1f} {lbl_y:.1f})" '
                f'class="{lbl_class} track-item" data-track-idx="{ti}" '
                f'font-size="{fsize}" '
                f'style="paint-order:stroke;stroke:#ffffff;stroke-width:2.6;'
                f'stroke-linejoin:round">{_svg_esc(label)}</text>'
            )
    svg_parts.append('</g>')

    # ---- Pinned-labels layer (populated by JS on click-to-pin) ----
    svg_parts.append('<g class="pinned-labels"></g>')

    # ---- Gap breaks & chromosome labels ----
    for ti, track in enumerate(all_tracks):
        yb = track_y[ti]
        x_off = track["offset"]
        th = track_heights[ti]

        # Compressed-gap labels used to sit at the track's vertical centre,
        # landing on top of the gene arrows and reading as if they belonged to
        # a gene rather than to the gap between two. Draw them ABOVE the track
        # instead, mark the exact break with a faint vertical line through the
        # gap, and skip a label that would collide with the previous one.
        last_lbl_x = -1.0e9
        for brk in sorted(track.get("breaks", []), key=lambda b: b["x"]):
            brk_px = bp2px(brk["x"] - x_off)
            if brk.get("is_chrom_break"):
                svg_parts.append(
                    f'<line x1="{brk_px:.1f}" y1="{yb - TRACK_PAD:.1f}" '
                    f'x2="{brk_px:.1f}" y2="{yb + th + TRACK_PAD:.1f}" '
                    f'class="chrom-break-line track-item" data-track-idx="{ti}"/>'
                )
            else:
                # faint vertical tick marking WHERE the gap was compressed
                svg_parts.append(
                    f'<line x1="{brk_px:.1f}" y1="{yb - 2:.1f}" '
                    f'x2="{brk_px:.1f}" y2="{yb + th + 2:.1f}" '
                    f'class="track-item" data-track-idx="{ti}" '
                    f'stroke="#cbd5e1" stroke-width="1" stroke-dasharray="2,3"/>'
                )
                if abs(brk_px - last_lbl_x) >= 48:
                    svg_parts.append(
                        f'<text x="{brk_px:.1f}" y="{yb - 6:.1f}" '
                        f'text-anchor="middle" '
                        f'class="break-label track-item" data-track-idx="{ti}">'
                        f'// {_svg_esc(brk["text"])}</text>'
                    )
                    last_lbl_x = brk_px

        # Chromosome labels below segments
        if track["genes"]:
            from collections import OrderedDict
            chrom_segs = OrderedDict()
            for g in track["genes"]:
                ch = g["chrom"]
                if ch not in chrom_segs:
                    chrom_segs[ch] = []
                chrom_segs[ch].append(g)

            if len(chrom_segs) > 1:
                for ch, ch_genes in chrom_segs.items():
                    xs = []
                    for g in ch_genes:
                        gx0, gx1 = gene_px(g, track)
                        xs.extend([gx0, gx1])
                    cx = (min(xs) + max(xs)) / 2
                    short = ch if len(ch) <= 14 else ch[-12:]
                    svg_parts.append(
                        f'<text x="{cx:.1f}" y="{yb + th + TRACK_PAD + 10:.1f}" '
                        f'text-anchor="middle" class="chrom-label track-item" data-track-idx="{ti}">{_svg_esc(short)}</text>'
                    )

    # ---- Scale bar ----
    # Wrapped in a group so the JS reflow can keep it pinned to the (shrinking)
    # canvas bottom when genomes are removed.
    scale_len_bp = args.scale_bar_len
    sb_x1_px = bp2px(x_max_bp - pad_bp * 0.5)
    sb_x0_px = sb_x1_px - scale_len_bp * scale
    sb_y = _content_bottom - BOTTOM_MARGIN + 20
    svg_parts.append('<g class="scale-bar-group">')
    svg_parts.append(
        f'<line x1="{sb_x0_px:.1f}" y1="{sb_y:.1f}" '
        f'x2="{sb_x1_px:.1f}" y2="{sb_y:.1f}" '
        f'class="scale-bar-line"/>'
    )
    # Ticks at ends
    svg_parts.append(
        f'<line x1="{sb_x0_px:.1f}" y1="{sb_y - 4:.1f}" '
        f'x2="{sb_x0_px:.1f}" y2="{sb_y + 4:.1f}" class="scale-bar-line"/>'
    )
    svg_parts.append(
        f'<line x1="{sb_x1_px:.1f}" y1="{sb_y - 4:.1f}" '
        f'x2="{sb_x1_px:.1f}" y2="{sb_y + 4:.1f}" class="scale-bar-line"/>'
    )
    svg_parts.append(
        f'<text x="{(sb_x0_px + sb_x1_px) / 2:.1f}" y="{sb_y + 18:.1f}" '
        f'text-anchor="middle" class="scale-bar-text">'
        f'{_format_bp_label(scale_len_bp)}</text>'
    )
    svg_parts.append('</g>')

    # ---- Track labels (left margin) ----
    # Layout:
    #   x=  8..26 → collapse-toggle button (own column, never under text)
    #   x=  LEFT_MARGIN-14 right edge → species, accession, status (right-aligned)
    TOGGLE_X      = 8
    TOGGLE_W      = 18
    TEXT_RIGHT_X  = LEFT_MARGIN - 14

    svg_parts.append('<g class="track-labels">')
    for ti, track in enumerate(all_tracks):
        yb = track_y[ti] + GENE_H / 2
        label = re.sub(r"<[^>]+>", "", track["label"])

        # Parse out accession if possible (assume label is "Species_name (Accession)")
        species = label
        acc = ""
        m = re.search(r"^(.*?)\s*\(([^)]+)\)$", label)
        if m:
            species = m.group(1).strip()
            acc = m.group(2).strip()

        # Collapse toggle for non-home tracks, in its own left-edge column so
        # the species name never overlaps it.
        if not track.get("is_home"):
            svg_parts.append(
                f'<g class="track-toggle" data-track-idx="{ti}">'
                f'<rect x="{TOGGLE_X}" y="{yb - TOGGLE_W/2:.1f}" '
                f'width="{TOGGLE_W}" height="{TOGGLE_W}" rx="3"/>'
                f'<text x="{TOGGLE_X + TOGGLE_W/2:.1f}" y="{yb:.1f}" '
                f'text-anchor="middle" dominant-baseline="central">▼</text>'
                f'</g>'
            )

        # (Clade-colour dot removed — the green/red per-genome swatch was
        # redundant with the tree plot and visually noisy.)

        # Species/common name (NCBI 'datasets' lookup if enabled).
        # data-track-idx on EVERY per-track label so hide/reflow toggles them
        # all (the species name, accession and status used to be left behind
        # when a track was collapsed because they carried no track index).
        species_label = _format_species_label(species)
        svg_parts.append(
            f'<text x="{TEXT_RIGHT_X}" y="{yb:.1f}" data-track-idx="{ti}" '
            f'text-anchor="end" class="track-label track-item" font-style="italic">'
            f'{_svg_esc(species_label)}</text>'
        )

        # Accession and span
        if acc:
            span_str = acc
            if track["genes"]:
                chroms = sorted({g["chrom"] for g in track["genes"]})
                chr_str = chroms[0] if len(chroms) == 1 else f"{len(chroms)} chr"
                span_str = f"{acc} • {chr_str}"

            svg_parts.append(
                f'<text x="{TEXT_RIGHT_X}" y="{yb + 14:.1f}" data-track-idx="{ti}" '
                f'text-anchor="end" class="chrom-label track-item">{_svg_esc(span_str)}</text>'
            )

        # GOI status indicators
        if not track["is_home"]:
            status = track.get("goi_status", "")
            if status == "absent":
                svg_parts.append(
                    f'<text x="{TEXT_RIGHT_X}" y="{yb + 28:.1f}" data-track-idx="{ti}" '
                    f'text-anchor="end" class="goi-status absent track-item">✗ GOI not found</text>'
                )
            elif status == "ambiguous":
                svg_parts.append(
                    f'<text x="{TEXT_RIGHT_X}" y="{yb + 28:.1f}" data-track-idx="{ti}" '
                    f'text-anchor="end" class="goi-status ambiguous track-item">⚠ Ambiguous orthology</text>'
                )
    svg_parts.append('</g>')

    # ---- Legend (Bottom left) ----
    # Removed as requested.

    # ---- Title & subtitle ----
    # Sits in the dedicated header band above TOP_MARGIN. The subtitle is the
    # descriptive legend line; it can be long, so its font shrinks to stay
    # within the canvas instead of spilling over the plot (the old fixed size
    # ran the text — and the centred 'GOI' tag — onto the home panel).
    title_x = plot_w / 2
    svg_parts.append(
        f'<text x="{title_x:.1f}" y="34" text-anchor="middle" '
        f'class="plot-title">SynVoy Synteny Plot</text>'
    )
    if subtitle_bits:
        sub_text = " · ".join(subtitle_bits)
        sub_fs = max(8.5, min(11.0, (plot_w - 60) / (len(sub_text) * 0.52)))
        svg_parts.append(
            f'<text x="{title_x:.1f}" y="60" text-anchor="middle" '
            f'class="plot-subtitle" style="font-size:{sub_fs:.1f}px">'
            f'{_svg_esc(sub_text)}</text>'
        )

    # ---- Caption block (bottom) ----
    if _caption_lines:
        svg_parts.append('<g class="plot-caption">')
        svg_parts.append(
            f'<line x1="{LEFT_MARGIN}" y1="{_content_bottom + 4:.1f}" '
            f'x2="{plot_w - RIGHT_MARGIN}" y2="{_content_bottom + 4:.1f}" '
            f'stroke="#dddddd" stroke-width="1"/>'
        )
        for i, line in enumerate(_caption_lines):
            weight = "600" if i == 0 else "400"
            svg_parts.append(
                f'<text x="{LEFT_MARGIN}" y="{_content_bottom + 18 + i * 15:.1f}" '
                f'style="font-size:10px;fill:#444;font-weight:{weight}">'
                f'{_svg_esc(line)}</text>'
            )
        svg_parts.append('</g>')

    # ---- Per-track geometry for the JS reflow ----
    # When a genome is removed, the JS reclaims its vertical space by re-laying
    # out the remaining tracks (translate + canvas resize) rather than leaving
    # a blank gap. It needs each track's original top + height and the margins.
    layout = {
        "topMargin":    TOP_MARGIN,
        "trackMargin":  TRACK_MARGIN,
        "bottomMargin": BOTTOM_MARGIN,
        "trackPad":     TRACK_PAD,
        "width":        plot_w,
        "guideX":       (bp2px(0) if (all_tracks and all_tracks[0]["genes"]) else None),
        "tracks": [{"idx": ti, "top": round(track_y[ti], 2),
                    "h": round(track_heights[ti], 2)}
                   for ti in range(n_tracks)],
    }

    # ---- Assemble full HTML ----
    svg_content = "\n".join(svg_parts)
    html = _assemble_full_html(svg_content, plot_w, total_h, json.dumps(layout))
    return html


def _render_tree_svg(tree_file, goi_genome_colours, output_path, species_map=None,
                     clade_count=4):
    """Render a horizontal dendrogram of the GOI phylogenetic tree as SVG HTML.

    The tree is midpoint-rooted and split into ``clade_count`` clades by
    cutting the K−1 longest non-root branches (see
    ``synvoy_tree.partition_clades``). Each clade gets one distinct colour
    from the colour-blind-safe ``CLADE_PALETTE``; leaves within the same
    cut subtree share that colour and so form visually contiguous groups —
    matching the matrix plot's idiom.

    The legacy per-genome ``goi_genome_colours`` is kept as a fallback for
    leaves that cannot be assigned to a clade (e.g. trees with <2 leaves).
    """
    if not tree_file or not os.path.exists(tree_file):
        return

    try:
        with open(tree_file) as fh:
            newick_str = fh.read()
        # Use the shared synvoy_tree parser → midpoint root → partition.
        # Side-effect: this is a separate code path from ete3, but the
        # implementation is dependency-free and produces equivalent layouts.
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import synvoy_tree as _stree
        raw = _stree.parse_newick_tree(newick_str)
        t = _stree.midpoint_root(raw) if raw is not None else None
    except Exception as exc:
        print(f"Warning: could not parse tree for rendering: {exc}")
        return

    if t is None:
        return

    leaves = list(t.leaves())
    if len(leaves) < 2:
        return

    # Clade partition — leaf name → 0-based clade id. K-longest-branches mode
    # (default K=4): cuts the K−1 deepest splits near the root so the
    # rendered colours track the natural major lineages.
    leaf_clade = _stree.partition_clades(t, target_k=clade_count)

    # Assign coordinates via recursive DFS using the rooted tree.
    node_coords = {}
    leaf_counter = [0]

    def _layout(node, x_offset):
        if node.is_leaf():
            y = leaf_counter[0]
            leaf_counter[0] += 1
            node_coords[id(node)] = (x_offset + node.dist, y)
        else:
            child_ys = []
            for child in node.children:
                _layout(child, x_offset + node.dist)
                child_ys.append(node_coords[id(child)][1])
            node_coords[id(node)] = (x_offset + node.dist,
                                     sum(child_ys) / len(child_ys))

    _layout(t, 0)

    n_leaves = len(leaves)
    max_x = max(c[0] for c in node_coords.values())
    if max_x <= 0:
        max_x = 1

    left_m = 40
    right_m = 320
    top_m = 60
    leaf_h = 50
    tree_w = 500
    total_w = left_m + tree_w + right_m
    total_h = top_m + n_leaves * leaf_h + 40
    x_scale = tree_w / max_x

    def tx(val):
        return left_m + val * x_scale

    def ty(val):
        return top_m + val * leaf_h + leaf_h / 2

    # Each non-root node inherits the clade colour of *any* descendant leaf
    # (since clades are contiguous by construction, all descendants share
    # the same id). This lets us colour internal branches consistently.
    def _node_clade_color(node):
        for leaf in _iter_leaves(node):
            cid = leaf_clade.get(leaf.name)
            if cid is not None:
                return _stree.color_for_clade(cid)
        return "#777"

    def _iter_leaves(node):
        if node.is_leaf():
            yield node
            return
        for c in node.children:
            yield from _iter_leaves(c)

    svg_parts = []

    # Branch lines, coloured by clade.
    for node in t.all_nodes():
        if node.parent is None:
            continue
        parent = node.parent
        px, py = node_coords[id(parent)]
        cx, cy = node_coords[id(node)]
        bcol = _node_clade_color(node)
        svg_parts.append(
            f'<polyline points="{tx(px):.1f},{ty(py):.1f} '
            f'{tx(px):.1f},{ty(cy):.1f} {tx(cx):.1f},{ty(cy):.1f}" '
            f'fill="none" stroke="{bcol}" stroke-width="1.5"/>'
        )

    # Clean, publication-ready leaf labels: the species name (from the GFF /
    # species-map), never the internal SynVoy locus id. A leaf name like
    # 'GOI_DCN|zebrafish_fna_b2_l2_fallback' becomes 'Danio rerio'; when a
    # species has several recovered copies they are disambiguated 'Danio rerio
    # · copy N' rather than leaking '_b2_l2_fallback' / '_b0_l1_exon_ann'.
    def _clean_species_for_leaf(name):
        if "|" not in name:
            return None  # home / query leaf
        sp = _stree.species_from_leaf(name) or name.split("|", 1)[1]
        if species_map:
            for acc, sp_name in species_map.items():
                if acc and (acc in sp or sp in acc):
                    return sp_name
        return sp.replace("_", " ")
    from collections import Counter, defaultdict as _dd
    _sp_total = Counter(s for s in (_clean_species_for_leaf(l.name) for l in leaves)
                        if s)
    _sp_seen = _dd(int)

    # Leaf nodes coloured by clade. Falls back to the legacy per-genome
    # palette only if clade assignment failed for this leaf.
    for leaf in leaves:
        lx, ly = node_coords[id(leaf)]
        cid = leaf_clade.get(leaf.name)
        if cid is not None:
            colour = _stree.color_for_clade(cid)
        else:
            gid = _genome_id_from_leaf(leaf.name)
            key = gid if gid else "home"
            colour = GOI_COLOUR
            if goi_genome_colours:
                if key in goi_genome_colours:
                    colour = goi_genome_colours[key]
                else:
                    for k, c in goi_genome_colours.items():
                        if k in key or key in k:
                            colour = c
                            break

        # Clean label (see _clean_species_for_leaf above).
        sp_label = _clean_species_for_leaf(leaf.name)
        if sp_label is None:
            label = f"{clean_gene_label(leaf.name)} (home)"
        else:
            _sp_seen[sp_label] += 1
            label = (f"{sp_label} · copy {_sp_seen[sp_label]}"
                     if _sp_total[sp_label] > 1 else sp_label)

        cpx, cpy = tx(lx), ty(ly)
        svg_parts.append(
            f'<circle cx="{cpx:.1f}" cy="{cpy:.1f}" r="7" '
            f'fill="{colour}" stroke="#333" stroke-width="1"/>'
        )
        svg_parts.append(
            f'<text x="{cpx + 14:.1f}" y="{cpy + 4:.1f}" '
            f'font-size="11" fill="#333">{_svg_esc(label)}</text>'
        )

    # Title
    svg_parts.append(
        f'<text x="{total_w / 2:.1f}" y="30" text-anchor="middle" '
        f'font-size="16" font-weight="700" fill="#1a1d26">'
        f'SynVoy GOI Phylogenetic Tree</text>'
    )

    # X-axis label
    svg_parts.append(
        f'<text x="{left_m + tree_w / 2:.1f}" y="{total_h - 10:.1f}" '
        f'text-anchor="middle" font-size="12" fill="#6b7280">'
        f'Evolutionary distance</text>'
    )

    svg_content = "\n".join(svg_parts)

    tree_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SynVoy GOI Phylogenetic Tree</title>
<style>
  body {{
    margin: 0; padding: 20px;
    background: #f8f9fb;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  }}
  .tree-container {{
    background: #fff;
    border-radius: 12px;
    box-shadow: 0 4px 20px rgba(0,0,0,0.08);
    padding: 20px;
    display: inline-block;
  }}
</style>
</head>
<body>
<div class="tree-container">
<svg width="{total_w}" height="{total_h}" xmlns="http://www.w3.org/2000/svg"
     style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
{svg_content}
</svg>
</div>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(tree_html)
    print(f"Tree plot saved to {output_path}")


# ======================================================================
# Publication SVG Renderer
# ======================================================================
#
# The publication SVG is the same content as the interactive HTML, but
# rendered with every home-genome gene labelled on the canvas (so the
# figure is self-describing in print). The narrow vertical layout the
# previous renderer produced has been removed — the user wanted the
# publication output to be a verbatim mirror of the HTML view.


def render_publication_svg(all_tracks, gene_colours, goi_genome_colours,
                           home_products, args,
                           subtitle_bits, hidden_absent_tracks,
                           ambiguous_track_count, resolved_track_count):
    """Render publication SVG = HTML plot with home-genome gene names forced on."""
    pub_html = render_synteny_html(
        all_tracks, gene_colours, goi_genome_colours,
        home_products, args,
        subtitle_bits, hidden_absent_tracks,
        ambiguous_track_count, resolved_track_count,
        force_home_labels=True,
    )

    style_blocks = re.findall(r"<style[^>]*>(.*?)</style>", pub_html, re.DOTALL)
    combined_style = "\n".join(style_blocks)
    m = re.search(r"(<svg[^>]*>)(.*?)(</svg>)", pub_html, re.DOTALL)
    if not m:
        raise RuntimeError("no <svg>...</svg> in publication HTML render")
    svg_open, svg_inner, svg_close = m.group(1), m.group(2), m.group(3)
    if 'xmlns="http://www.w3.org/2000/svg"' not in svg_open:
        svg_open = svg_open.replace(
            "<svg", '<svg xmlns="http://www.w3.org/2000/svg"', 1
        )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        + svg_open
        + "\n<style type=\"text/css\"><![CDATA[\n"
        + combined_style
        + "\n]]></style>\n"
        + svg_inner
        + svg_close
    )


# ======================================================================
# HTML + CSS + JS Templates
# ======================================================================

_CSS_TEMPLATE = """
:root {
  --bg: #f8f9fb;
  --surface: #ffffff;
  --track-bg: #f3f5f8;
  --track-border: #e8eaef;
  --text-primary: #1a1d26;
  --text-secondary: #555d6e;
  --text-muted: #8c95a6;
  --goi-color: #e31a1c;
  --goi-dark: #8b0000;
}
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
  background: var(--bg);
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
               'Helvetica Neue', Arial, sans-serif;
  color: var(--text-primary);
}
.toolbar {
  position: fixed; top: 12px; right: 16px; z-index: 100;
  display: flex; gap: 6px;
  background: rgba(255,255,255,0.92);
  backdrop-filter: blur(8px);
  border: 1px solid var(--track-border);
  border-radius: 8px; padding: 4px 6px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.08);
}
.toolbar button {
  width: 32px; height: 32px; border: none;
  background: transparent; border-radius: 6px;
  font-size: 18px; color: var(--text-secondary);
  cursor: pointer; display: flex; align-items: center;
  justify-content: center; transition: all 0.15s ease;
}
.toolbar button:hover {
  background: var(--track-bg); color: var(--text-primary);
}
.toolbar button.active {
  background: var(--track-bg); color: var(--text-primary);
  box-shadow: inset 0 0 0 1px var(--track-border);
}
.toolbar button { position: relative; }
/* Count badge on the track-manager button: how many genomes are hidden, so
   the restore control is discoverable after a row is removed. */
.tm-badge {
  position: absolute; top: -3px; right: -3px;
  min-width: 15px; height: 15px; padding: 0 3px;
  background: #e3342f; color: #fff;
  border-radius: 8px; font-size: 10px; font-weight: 700;
  line-height: 15px; text-align: center;
  box-shadow: 0 0 0 1.5px #fff;
}
.tm-badge[hidden] { display: none; }
.ribbons.ribbons-hidden { display: none; }
.plot-wrapper {
  width: 100%; overflow: auto; padding: 16px;
}
.zoom-container {
  transform-origin: left top;
  transition: transform 0.15s ease;
  display: inline-block;
}
.synteny-svg {
  background: var(--surface);
  border-radius: 12px;
  box-shadow: 0 4px 24px rgba(0,0,0,0.07), 0 1px 4px rgba(0,0,0,0.04);
  display: block;
}

/* Track backgrounds */
.track-bg {
  fill: var(--track-bg);
  stroke: var(--track-border);
  stroke-width: 0.5;
}

/* Gene groups */
.gene-group {
  cursor: pointer;
  filter: drop-shadow(0 1px 2px rgba(0,0,0,0.08));
  transition: filter 0.15s ease, opacity 0.2s ease;
}
.gene-group:hover {
  filter: drop-shadow(0 2px 6px rgba(0,0,0,0.15)) brightness(1.06);
}
.gene-group .exon {
  transition: filter 0.15s ease;
}
/* Right-click pinned: subtle dark outline so it's visually marked as
   "I asked for this label" without changing the gene fill. */
.gene-group.pinned .exon {
  stroke: #1a1d26;
}

/* Intron lines */
.intron-line {
  stroke: #94a3b8;
  stroke-width: 1.5;
}

/* Ribbons */
.ribbon {
  transition: opacity 0.2s ease;
}
.ribbon:hover {
  opacity: 0.45 !important;
}

/* Labels */
.gene-label {
  fill: var(--text-secondary);
  font-weight: 500;
  pointer-events: none;
}
.gene-label.goi {
  fill: var(--goi-dark);
  font-weight: 700;
}

/* Track labels */
.track-label {
  fill: var(--text-primary);
  font-size: 13px;
  font-weight: 600;
    pointer-events: none;
}
.track-toggle { cursor: pointer; }
.track-toggle rect {
    fill: var(--track-bg);
    stroke: var(--track-border);
    stroke-width: 1;
}
.track-toggle:hover rect {
    fill: var(--surface);
    stroke: var(--text-muted);
}
.track-toggle text {
    fill: var(--text-secondary);
    font-size: 12px;
    font-weight: 700;
    user-select: none;
    pointer-events: none;
}
.track-toggle:hover text { fill: var(--text-primary); }
.track-summary {
    display: none;
    fill: var(--text-secondary);
    font-size: 10px;
    font-style: italic;
    pointer-events: none;
}
.track-manager {
    position: fixed;
    top: 52px;
    right: 16px;
    z-index: 101;
    width: 280px;
    max-height: 55vh;
    overflow: auto;
    background: rgba(255,255,255,0.96);
    border: 1px solid var(--track-border);
    border-radius: 12px;
    box-shadow: 0 10px 34px rgba(0,0,0,0.14);
    padding: 12px;
    display: none;
}
.track-manager.visible {
    display: block;
}
.track-manager h3 {
    font-size: 13px;
    margin-bottom: 4px;
    color: var(--text-primary);
}
.track-manager .tm-hint {
    font-size: 10.5px;
    line-height: 1.35;
    color: var(--text-secondary);
    margin-bottom: 8px;
}
.track-manager label {
    display: flex;
    align-items: flex-start;
    gap: 8px;
    padding: 6px 4px;
    border-radius: 8px;
    font-size: 12px;
    color: var(--text-secondary);
}
.track-manager label:hover {
    background: var(--track-bg);
}
.track-manager input {
    margin-top: 2px;
}
.track-item-hidden {
    display: none !important;
}

/* GOI status indicators */
.goi-status {
  font-size: 10px;
  font-weight: 500;
}
.goi-status.absent { fill: #dc2626; }
.goi-status.ambiguous { fill: #d97706; }

/* Chromosome labels & breaks */
.chrom-label {
  fill: var(--text-muted);
  font-size: 10px;
  font-weight: 500;
    pointer-events: none;
}
.chrom-break-line {
  stroke: #94a3b8;
  stroke-width: 1.5;
  stroke-dasharray: 5,4;
    pointer-events: none;
}
.break-label {
  fill: var(--text-muted);
  font-size: 9px;
  font-weight: 600;
    pointer-events: none;
}

/* Scale bar */
.scale-bar-line {
  stroke: #6b7280;
  stroke-width: 2;
}
.scale-bar-text {
  fill: #6b7280;
  font-size: 12px;
  font-weight: 600;
}

/* Title */
.plot-title {
  fill: var(--text-primary);
  font-size: 18px;
  font-weight: 700;
}
.plot-subtitle {
  fill: var(--text-muted);
  font-size: 11px;
}

/* Pinned gene labels (added by JS on right-click) */
.pinned-label {
  fill: var(--text-primary);
  font-size: 11px;
  font-weight: 600;
  pointer-events: none;
}
.pinned-label-bg {
  fill: rgba(255,255,255,0.85);
  stroke: var(--track-border);
  stroke-width: 0.5;
  pointer-events: none;
}
.pinned-label-leader {
  stroke: var(--text-muted);
  stroke-width: 0.6;
  fill: none;
  pointer-events: none;
}

/* Highlight states (left-click follows a gene's orthologs across tracks) */
.gene-group.highlighted {
  filter: drop-shadow(0 2px 8px rgba(0,0,0,0.2)) brightness(1.1) !important;
  opacity: 1 !important;
}
.gene-group.dimmed {
  opacity: 0.18 !important;
  filter: saturate(0.3) !important;
}
.ribbon.highlighted { opacity: 0.5 !important; }
.ribbon.dimmed      { opacity: 0.03 !important; }

/* GOI pulse animation (one shot on load, draws the eye to the GOI). */
@keyframes goi-pulse {
  0%, 100% { filter: drop-shadow(0 1px 2px rgba(0,0,0,0.08)); }
  50%      { filter: drop-shadow(0 0 8px rgba(227,26,28,0.4)) brightness(1.1); }
}
.gene-group[data-is-goi="true"] {
  animation: goi-pulse 1.8s ease-in-out 3;
}

.track-bg.track-item-hidden,
.track-label.track-item-hidden,
.track-summary.track-item-hidden,
.track-toggle.track-item-hidden,
.chrom-label.track-item-hidden,
.chrom-break-line.track-item-hidden,
.break-label.track-item-hidden,
.goi-status.track-item-hidden,
.gene-label.track-item-hidden,
.gene-group.track-item-hidden {
    display: none !important;
}

/* Tooltip */
.tooltip {
  position: fixed;
  padding: 10px 14px;
  background: rgba(26,29,38,0.95);
  border: 1px solid rgba(255,255,255,0.1);
  border-radius: 10px;
  box-shadow: 0 8px 30px rgba(0,0,0,0.25);
  color: #e8eaef;
  font-size: 12px;
  line-height: 1.6;
  pointer-events: none;
  z-index: 1000;
  max-width: 360px;
  display: none;
  backdrop-filter: blur(12px);
  -webkit-backdrop-filter: blur(12px);
}
.tooltip .tt-name {
  font-size: 14px;
  font-weight: 700;
  color: #fff;
  margin-bottom: 6px;
  padding-bottom: 5px;
  border-bottom: 1px solid rgba(255,255,255,0.12);
}
.tooltip .tt-product {
  font-style: italic;
  color: #a5b4c8;
  margin-bottom: 6px;
}
.tooltip .tt-row {
  display: flex;
  justify-content: space-between;
  gap: 16px;
}
.tooltip .tt-label {
  color: #8c95a6;
  white-space: nowrap;
}
.tooltip .tt-value {
  color: #e8eaef;
  text-align: right;
  font-weight: 500;
}
.tooltip .tt-goi {
  color: #ff6b6b;
  font-weight: 700;
  margin-top: 6px;
  padding-top: 5px;
  border-top: 1px solid rgba(255,255,255,0.12);
}
"""

_JS_TEMPLATE = r"""
document.addEventListener('DOMContentLoaded', () => {
  const tooltip = document.getElementById('tooltip');
  const svg = document.querySelector('.synteny-svg');
  if (!svg) return;

    const trackManager = document.getElementById('track-manager');
    const trackManagerBtn = document.getElementById('track-manager-btn');
    const tmBadge = document.getElementById('track-manager-badge');
    const ribbonGroups = Array.from(svg.querySelectorAll('.ribbons'));
    const ribbonsGroup = ribbonGroups[0] || null;  // JS-rebuilt ribbons land here
    const guideLine = svg.querySelector('.goi-guide-line');
    const guideTag = svg.querySelector('.goi-guide-tag');
    const scaleBarGroup = svg.querySelector('.scale-bar-group');
    const LAYOUT = window.__SYNVOY_LAYOUT__ || null;
    const baseH = LAYOUT ? Math.max.apply(null, LAYOUT.tracks.map(t => t.top + t.h)) + LAYOUT.bottomMargin : 0;
    const trackCount = new Set(Array.from(svg.querySelectorAll('[data-track-idx]')).map(el => el.dataset.trackIdx).filter(Boolean)).size;
    const collapsedTracks = new Set();
    let ribbonsRebuilt = false;  // first removal swaps the static gradient ribbons for JS-rebuilt ones

  // ---- Zoom controls ----
  let zoom = 1;
  const container = document.querySelector('.zoom-container');
  const zoomIn = document.getElementById('zoom-in');
  const zoomOut = document.getElementById('zoom-out');
  const zoomReset = document.getElementById('zoom-reset');

  function applyZoom() {
    container.style.transform = 'scale(' + zoom + ')';
  }
  if (zoomIn) zoomIn.addEventListener('click', () => { zoom = Math.min(5, zoom * 1.25); applyZoom(); });
  if (zoomOut) zoomOut.addEventListener('click', () => { zoom = Math.max(0.2, zoom / 1.25); applyZoom(); });
  if (zoomReset) zoomReset.addEventListener('click', () => { zoom = 1; applyZoom(); });

  // ---- Ribbon visibility toggle ----
  const toggleRibbonsBtn = document.getElementById('toggle-ribbons');
  if (toggleRibbonsBtn && ribbonGroups.length) {
    toggleRibbonsBtn.addEventListener('click', () => {
      const hidden = !ribbonGroups[0].classList.contains('ribbons-hidden');
      ribbonGroups.forEach(g => g.classList.toggle('ribbons-hidden', hidden));
      toggleRibbonsBtn.classList.toggle('active', hidden);
      toggleRibbonsBtn.setAttribute('title', hidden
        ? 'Show synteny ribbons (connections between tracks)'
        : 'Hide synteny ribbons (connections between tracks)');
    });
  }

    // ---- Track manager / collapse controls ----
    function uniqueTrackIndices() {
        const indices = new Set();
        svg.querySelectorAll('[data-track-idx]').forEach(el => {
            if (el.dataset.trackIdx !== undefined && el.dataset.trackIdx !== '') indices.add(el.dataset.trackIdx);
        });
        return Array.from(indices).sort((a, b) => Number(a) - Number(b));
    }

    function trackLabelFor(idx) {
        const labelEl = svg.querySelector('.track-label[data-track-idx="' + idx + '"]');
        if (labelEl) return labelEl.textContent.trim();
        const summaryEl = svg.querySelector('.track-summary[data-track-idx="' + idx + '"]');
        if (summaryEl) return summaryEl.textContent.replace(/^⋯\\s*/, '').replace(/\\s*\\(collapsed\\)$/, '').trim();
        return 'Track ' + idx;
    }

    function visibleIndices() {
        return uniqueTrackIndices().filter(idx => !collapsedTracks.has(idx));
    }

    function updateBadge() {
        if (!tmBadge) return;
        const n = collapsedTracks.size;
        tmBadge.textContent = String(n);
        tmBadge.hidden = (n === 0);
    }

    function updateGuide(dyOf) {
        if (!guideLine || !LAYOUT) return;
        const vis = LAYOUT.tracks.filter(t => !collapsedTracks.has(String(t.idx)));
        if (!vis.length) return;
        const first = vis[0], last = vis[vis.length - 1];
        const pad = LAYOUT.trackPad || 10;
        const yTop = first.top + (dyOf[String(first.idx)] || 0) - pad - 14;
        const yBot = last.top + (dyOf[String(last.idx)] || 0) + last.h + pad;
        guideLine.setAttribute('y1', yTop.toFixed(1));
        guideLine.setAttribute('y2', yBot.toFixed(1));
        if (guideTag) guideTag.setAttribute('y', (yTop - 2).toFixed(1));
    }

    // Re-lay-out the visible tracks so a removed genome's row is reclaimed
    // (no blank gap) instead of merely hidden in place. Each track's elements
    // are translated to their new vertical slot, the canvas is resized, and
    // the ribbons are rebuilt between the now-adjacent visible tracks.
    function applyLayout(rebuild) {
        const dyOf = {};
        if (LAYOUT && LAYOUT.tracks) {
            let cursor = LAYOUT.topMargin;
            let lastBottom = LAYOUT.topMargin;
            LAYOUT.tracks.forEach(t => {
                const idx = String(t.idx);
                if (collapsedTracks.has(idx)) return;
                dyOf[idx] = cursor - t.top;
                lastBottom = cursor + t.h;
                cursor = lastBottom + LAYOUT.trackMargin;
            });
            LAYOUT.tracks.forEach(t => {
                const idx = String(t.idx);
                const hidden = collapsedTracks.has(idx);
                const dy = dyOf[idx] || 0;
                svg.querySelectorAll('[data-track-idx="' + idx + '"]').forEach(el => {
                    el.classList.toggle('track-item-hidden', hidden);
                    if (!hidden && dy) el.setAttribute('transform', 'translate(0,' + dy.toFixed(2) + ')');
                    else el.removeAttribute('transform');
                });
            });
            const newH = lastBottom + LAYOUT.bottomMargin;
            svg.setAttribute('height', newH);
            if (scaleBarGroup) {
                const dsb = newH - baseH;
                if (dsb) scaleBarGroup.setAttribute('transform', 'translate(0,' + dsb.toFixed(2) + ')');
                else scaleBarGroup.removeAttribute('transform');
            }
            updateGuide(dyOf);
        } else {
            // No layout payload: fall back to hide-in-place (older renders).
            uniqueTrackIndices().forEach(idx => {
                const hidden = collapsedTracks.has(idx);
                svg.querySelectorAll('[data-track-idx="' + idx + '"]').forEach(el => {
                    el.classList.toggle('track-item-hidden', hidden);
                });
            });
        }
        if (rebuild || ribbonsRebuilt) rebuildRibbons(dyOf);
        updateBadge();
    }

    function setTrackCollapsed(idx, collapsed) {
        idx = String(idx);
        if (idx === '0') return;  // the home reference always stays
        collapsedTracks[collapsed ? 'add' : 'delete'](idx);
        svg.querySelectorAll('.track-toggle[data-track-idx="' + idx + '"]').forEach(el => {
            const tEl = el.querySelector('text') || el;
            tEl.textContent = collapsed ? '▶' : '▼';
        });
        applyLayout(true);
    }

    function buildTrackManager() {
        if (!trackManager) return;
        const nHidden = collapsedTracks.size;
        trackManager.innerHTML = '<h3>Show / hide genomes</h3>'
            + '<p class="tm-hint">Uncheck to remove a genome (its row collapses and '
            + 'ribbons rebuild); re-check to restore it.'
            + (nHidden ? ' <b>' + nHidden + ' hidden.</b>' : '') + '</p>';
        uniqueTrackIndices().forEach(idx => {
            const row = document.createElement('label');
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.checked = !collapsedTracks.has(idx);
            if (idx === '0') {
                cb.checked = true;
                cb.disabled = true;
            }
            cb.addEventListener('change', () => { setTrackCollapsed(idx, !cb.checked); buildTrackManager(); });
            const span = document.createElement('span');
            span.textContent = trackLabelFor(idx);
            row.appendChild(cb);
            row.appendChild(span);
            trackManager.appendChild(row);
        });
    }

    if (trackManagerBtn && trackManager) {
        trackManagerBtn.addEventListener('click', () => {
            trackManager.classList.toggle('visible');
            if (trackManager.classList.contains('visible')) buildTrackManager();
        });
        document.addEventListener('click', (e) => {
            if (!trackManager.contains(e.target) && e.target !== trackManagerBtn) {
                trackManager.classList.remove('visible');
            }
        });
    }

    svg.querySelectorAll('.track-toggle').forEach(el => {
        el.addEventListener('click', (e) => {
            e.stopPropagation();
            const idx = el.dataset.trackIdx;
            setTrackCollapsed(idx, !collapsedTracks.has(idx));
            buildTrackManager();
        });
    });

    function buildRibbonPath(ux0, ux1, uy, lx0, lx1, ly) {
        const cy = (uy + ly) / 2;
        return 'M' + ux0.toFixed(2) + ',' + uy.toFixed(2) +
                     ' C' + ux0.toFixed(2) + ',' + cy.toFixed(2) + ' ' + lx0.toFixed(2) + ',' + cy.toFixed(2) +
                     ' ' + lx0.toFixed(2) + ',' + ly.toFixed(2) + ' L' + lx1.toFixed(2) + ',' + ly.toFixed(2) +
                     ' C' + lx1.toFixed(2) + ',' + cy.toFixed(2) + ' ' + ux1.toFixed(2) + ',' + cy.toFixed(2) +
                     ' ' + ux1.toFixed(2) + ',' + uy.toFixed(2) + ' Z';
    }

    function rebuildRibbons(dyOf) {
        if (!ribbonsGroup) return;
        dyOf = dyOf || {};
        // Clear EVERY ribbon group (the static render has a flank group AND a
        // GOI gradient group; clearing only the first left orphan GOI ribbons
        // dangling when a genome was removed). All rebuilt ribbons go into one.
        ribbonGroups.forEach(g => { g.innerHTML = ''; });
        ribbonsRebuilt = true;
        function emitRibbon(ug, lg, du, dl, isGoi, hom) {
            const ux0 = parseFloat(ug.dataset.x0), ux1 = parseFloat(ug.dataset.x1);
            const uy = parseFloat(ug.dataset.yb) + 11 + du;
            const lx0 = parseFloat(lg.dataset.x0), lx1 = parseFloat(lg.dataset.x1);
            const ly = parseFloat(lg.dataset.yb) - 11 + dl;
            const identity = Math.max(0, Math.min(100, parseFloat(lg.dataset.identity || '50')));
            const clr = lg.dataset.fill || '#909090';
            const alpha = isGoi ? 0.45 : (0.08 + (identity / 100) * 0.35);
            const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
            path.setAttribute('class', 'ribbon');
            path.setAttribute('data-homology', isGoi ? 'GOI' : hom);
            path.setAttribute('data-upper-track', String(upperIdx));
            path.setAttribute('data-lower-track', String(lowerIdx));
            path.setAttribute('d', buildRibbonPath(ux0, ux1, uy, lx0, lx1, ly));
            path.setAttribute('fill', rgbaFromHex(clr, alpha));
            path.setAttribute('stroke', rgbaFromHex(clr, Math.min(1, alpha * 1.6)));
            path.setAttribute('stroke-width', isGoi ? '0.6' : '0.3');
            ribbonsGroup.appendChild(path);
        }
        let upperIdx, lowerIdx;
        const bestGoi = genes => {
            // The GOI is matched SEMANTICALLY (data-is-goi), not by the homology
            // string — the home GOI's id (e.g. 'gene-Melt') differs from the
            // targets' ('GOI_Melt'), which used to break the home↔target GOI
            // ribbon after a reflow. Pick the single best (highest identity)
            // copy per track so multi-copy GOIs don't fan out into a dark mass.
            let best = null, bi = -1;
            genes.forEach(g => {
                if (g.dataset.isGoi !== 'true') return;
                const id = parseFloat(g.dataset.identity || '0');
                if (best === null || id > bi) { best = g; bi = id; }
            });
            return best;
        };
        const visible = visibleIndices();
        for (let vi = 0; vi < visible.length - 1; vi++) {
            upperIdx = visible[vi];
            lowerIdx = visible[vi + 1];
            const du = dyOf[upperIdx] || 0;
            const dl = dyOf[lowerIdx] || 0;
            const upperGenes = Array.from(svg.querySelectorAll('.gene-group[data-track-idx="' + upperIdx + '"]'));
            const lowerGenes = Array.from(svg.querySelectorAll('.gene-group[data-track-idx="' + lowerIdx + '"]'));
            // ---- flanking ribbons: match by homology id, GOI genes excluded ----
            const upperByHomology = new Map();
            upperGenes.forEach(g => {
                if (g.dataset.isGoi === 'true') return;
                const hom = g.dataset.homology;
                if (!upperByHomology.has(hom)) upperByHomology.set(hom, []);
                upperByHomology.get(hom).push(g);
            });
            lowerGenes.forEach(lg => {
                if (lg.dataset.isGoi === 'true') return;
                const matches = upperByHomology.get(lg.dataset.homology) || [];
                matches.forEach(ug => emitRibbon(ug, lg, du, dl, false, lg.dataset.homology));
            });
            // ---- one GOI ribbon, best upper copy ↔ best lower copy ----
            const ug = bestGoi(upperGenes), lg = bestGoi(lowerGenes);
            if (ug && lg) emitRibbon(ug, lg, du, dl, true, 'GOI');
        }
    }

    function rgbaFromHex(hex, alpha) {
        if (!hex) hex = '#909090';
        hex = String(hex).trim();
        // If input is already rgb/rgba, parse numbers and set alpha
        const rgbMatch = hex.match(/rgba?\s*\(\s*([0-9]+)\s*,\s*([0-9]+)\s*,\s*([0-9]+)(?:\s*,\s*([0-9\.]+))?\s*\)/i);
        if (rgbMatch) {
            const r = parseInt(rgbMatch[1], 10);
            const g = parseInt(rgbMatch[2], 10);
            const b = parseInt(rgbMatch[3], 10);
            return 'rgba(' + r + ',' + g + ',' + b + ',' + alpha.toFixed(3) + ')';
        }
        // Otherwise expect a hex string like #rrggbb
        const clean = hex.replace('#', '');
        if (clean.length === 3) {
            const r = parseInt(clean[0] + clean[0], 16);
            const g = parseInt(clean[1] + clean[1], 16);
            const b = parseInt(clean[2] + clean[2], 16);
            return 'rgba(' + r + ',' + g + ',' + b + ',' + alpha.toFixed(3) + ')';
        }
        const r = parseInt(clean.slice(0, 2), 16) || 0;
        const g = parseInt(clean.slice(2, 4), 16) || 0;
        const b = parseInt(clean.slice(4, 6), 16) || 0;
        return 'rgba(' + r + ',' + g + ',' + b + ',' + alpha.toFixed(3) + ')';
    }

    buildTrackManager();
    updateBadge();
    // Keep the server-rendered (gradient GOI) ribbons on first paint; they're
    // rebuilt only when the layout changes (a genome is removed / restored).

  // ---- Hover tooltips ----
  svg.querySelectorAll('.gene-group').forEach(el => {
    el.addEventListener('mouseenter', (e) => {
      try {
        const data = JSON.parse(el.dataset.tooltip);
        let html = '<div class="tt-name">' + esc(data.name) + '</div>';
        if (data.product) html += '<div class="tt-product">' + esc(data.product) + '</div>';
        if (data.coords) html += row('Coords', data.coords);
        if (data.strand) html += row('Strand', data.strand);
        if (data.exons) html += row('Exons', data.exons);
        if (data.homolog) html += row('Homolog', data.homolog);
        if (data.identity !== undefined) html += row('Identity', data.identity + '%');
        if (data.confidence) html += row('Confidence', data.confidence);
        if (data.goiClass) html += row('GOI class', data.goiClass);
        if (data.evidence) html += row('Evidence', data.evidence);
        if (data.model) html += row('Model', data.model);
        if (data.synteny) html += row('Synteny', data.synteny);
        if (data.queryCov !== undefined) html += row('Query cov', data.queryCov + '%');
        if (data.goiTag) html += '<div class="tt-goi">' + esc(data.goiTag) + '</div>';
        tooltip.innerHTML = html;
        tooltip.style.display = 'block';
        positionTooltip(e);
      } catch(err) {}
    });
    el.addEventListener('mousemove', positionTooltip);
    el.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
  });

  function esc(s) { const d = document.createElement('div'); d.textContent = s; return d.innerHTML; }
  function row(label, value) {
    return '<div class="tt-row"><span class="tt-label">' + label +
           '</span><span class="tt-value">' + esc(String(value)) + '</span></div>';
  }
  function positionTooltip(e) {
    let x = e.clientX + 16, y = e.clientY + 16;
    const r = tooltip.getBoundingClientRect();
    if (x + r.width > window.innerWidth - 8) x = e.clientX - r.width - 16;
    if (y + r.height > window.innerHeight - 8) y = e.clientY - r.height - 16;
    tooltip.style.left = Math.max(4, x) + 'px';
    tooltip.style.top = Math.max(4, y) + 'px';
  }

  // ---- Pinned gene labels (right-click) ----
  // Hover         → tooltip (above)
  // Left-click    → follow the gene's orthologs across tracks (below)
  // Right-click   → toggle a persistent label above that gene
  const pinnedLayer = svg.querySelector('.pinned-labels');
  const pinnedById = new Map();
  let pinSeq = 0;

  function geneKey(el) {
    if (!el.dataset.pinKey) el.dataset.pinKey = 'g' + (pinSeq++);
    return el.dataset.pinKey;
  }

  function pinLabel(el) {
    if (!pinnedLayer) return;
    const key = geneKey(el);
    if (pinnedById.has(key)) return;
    const x0 = parseFloat(el.dataset.x0);
    const x1 = parseFloat(el.dataset.x1);
    const yb = parseFloat(el.dataset.yb);
    const text = el.dataset.label || '';
    if (!text) return;
    const cx = (x0 + x1) / 2;
    const ly = yb - 6;

    const SVG_NS = 'http://www.w3.org/2000/svg';
    const group = document.createElementNS(SVG_NS, 'g');
    group.setAttribute('class', 'pinned-label-group');
    group.dataset.pinKey = key;

    const leader = document.createElementNS(SVG_NS, 'line');
    leader.setAttribute('x1', cx.toFixed(1));
    leader.setAttribute('y1', yb.toFixed(1));
    leader.setAttribute('x2', cx.toFixed(1));
    leader.setAttribute('y2', (ly + 2).toFixed(1));
    leader.setAttribute('class', 'pinned-label-leader');
    group.appendChild(leader);

    const tx = document.createElementNS(SVG_NS, 'text');
    tx.setAttribute('x', cx.toFixed(1));
    tx.setAttribute('y', ly.toFixed(1));
    tx.setAttribute('text-anchor', 'middle');
    tx.setAttribute('class', 'pinned-label');
    tx.textContent = text;
    pinnedLayer.appendChild(group);
    group.appendChild(tx);

    try {
      const bb = tx.getBBox();
      const pad = 3;
      const rect = document.createElementNS(SVG_NS, 'rect');
      rect.setAttribute('x', (bb.x - pad).toFixed(1));
      rect.setAttribute('y', (bb.y - pad).toFixed(1));
      rect.setAttribute('width',  (bb.width  + 2 * pad).toFixed(1));
      rect.setAttribute('height', (bb.height + 2 * pad).toFixed(1));
      rect.setAttribute('rx', '2');
      rect.setAttribute('class', 'pinned-label-bg');
      group.insertBefore(rect, tx);
    } catch(err) { /* getBBox may fail if not rendered yet */ }

    el.classList.add('pinned');
    pinnedById.set(key, group);
  }

  function unpinLabel(el) {
    const key = geneKey(el);
    const node = pinnedById.get(key);
    if (node) { node.remove(); pinnedById.delete(key); }
    el.classList.remove('pinned');
  }

  function togglePin(el) {
    if (pinnedById.has(geneKey(el))) unpinLabel(el);
    else pinLabel(el);
  }

  function clearAllPins() {
    pinnedById.forEach((node) => node.remove());
    pinnedById.clear();
    svg.querySelectorAll('.gene-group.pinned').forEach(g => g.classList.remove('pinned'));
  }

  svg.querySelectorAll('.gene-group').forEach(el => {
    el.addEventListener('contextmenu', (e) => {
      e.preventDefault();
      e.stopPropagation();
      togglePin(el);
    });
  });

  const clearPinsBtn = document.getElementById('clear-pins');
  if (clearPinsBtn) clearPinsBtn.addEventListener('click', clearAllPins);

  // ---- Left-click: follow a gene's orthologs across tracks ----
  let selectedHom = null;
  svg.querySelectorAll('.gene-group').forEach(el => {
    el.addEventListener('click', (e) => {
      e.stopPropagation();
      const hom = el.dataset.homology;
      if (selectedHom === hom) { clearHighlight(); return; }
      selectedHom = hom;
      highlightHomology(hom);
    });
  });

  document.addEventListener('click', (e) => {
    if (!e.target.closest('.gene-group')) clearHighlight();
  });

  function highlightHomology(homId) {
    svg.querySelectorAll('.gene-group').forEach(g => {
      if (g.dataset.homology === homId) {
        g.classList.add('highlighted'); g.classList.remove('dimmed');
      } else {
        g.classList.add('dimmed'); g.classList.remove('highlighted');
      }
    });
    svg.querySelectorAll('.ribbon').forEach(r => {
      if (r.dataset.homology === homId) {
        r.classList.add('highlighted'); r.classList.remove('dimmed');
      } else {
        r.classList.add('dimmed'); r.classList.remove('highlighted');
      }
    });
  }

  function clearHighlight() {
    selectedHom = null;
    svg.querySelectorAll('.highlighted, .dimmed').forEach(el => {
      el.classList.remove('highlighted', 'dimmed');
    });
  }
});
"""


def _assemble_full_html(svg_content, width, height, layout_json="null"):
    """Wrap SVG content in a complete self-contained HTML document."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SynVoy Synteny Plot</title>
<style>{_CSS_TEMPLATE}</style>
</head>
<body>
<div class="toolbar">
  <button id="zoom-in" title="Zoom in">+</button>
  <button id="zoom-out" title="Zoom out">−</button>
  <button id="zoom-reset" title="Reset zoom">⟲</button>
  <button id="clear-pins" title="Clear pinned labels (right-click a gene to pin its name)">⌫</button>
  <button id="toggle-ribbons" title="Hide/show synteny ribbons (connections between tracks)">≋</button>
  <button id="track-manager-btn" title="Show / hide genomes (restore removed rows here)">📋<span id="track-manager-badge" class="tm-badge" hidden></span></button>
</div>
<div id="track-manager" class="track-manager"></div>
<script>window.__SYNVOY_LAYOUT__ = {layout_json};</script>
<div class="plot-wrapper">
<div class="zoom-container">
<svg class="synteny-svg" width="{width}" height="{height}"
     xmlns="http://www.w3.org/2000/svg"
     style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
{svg_content}
</svg>
</div>
</div>
<div id="tooltip" class="tooltip"></div>
<script>{_JS_TEMPLATE}</script>
</body>
</html>"""


# ======================================================================
# Anchor-grid view  (matrix × synteny hybrid)
# ======================================================================
#
# A third figure, complementary to the ribbon plot and the matrix:
#   • columns  = home genes in genomic order (shared "anchor" lanes)
#   • rows     = home + each target species (phylogenetic order)
#   • each ortholog is drawn as a directional gene arrow inside its anchor
#     column, so homologues line up *vertically* across all species — the
#     synteny reads from the alignment itself, not from a tangle of ribbons.
# Identity drives fill shade, confidence drives the visual tier, strand vs.
# the home gene flags inversions, and a missing ortholog leaves a visible
# gap.  The GOI gets a full-height guide band so it can't be missed.  A small
# home-coordinate axis on top restores "where the genes actually are" and
# annotates large genomic gaps (the §18-bridged neighbourhood case).

def _goi_display_label(all_tracks):
    """Best short symbol for the GOI column (e.g. 'DCN'), else 'GOI'.

    Target GOI genes carry a SynVoy_Parent like 'GOI_DCN'; the synthetic home
    placeholder may instead be 'GOI_<chrom>_<pos>'. Prefer a short alphabetic
    symbol and reject accession/coordinate-style suffixes.
    """
    from collections import Counter
    cand = Counter()
    for tr in all_tracks:
        for g in tr.get("genes", []):
            for s in (g.get("home_gene_id", "") or "", g.get("name", "") or ""):
                if s.startswith("GOI_"):
                    suf = s[4:].split("|")[0]
                    if re.match(r"^[A-Za-z][A-Za-z0-9-]{0,11}$", suf):
                        cand[suf] += 1
    return cand.most_common(1)[0][0] if cand else "GOI"


def _grid_target_map(track, anchor_keys, goi_key):
    """Map a target track's genes onto anchor columns.

    Returns {anchor_key: {"best": gene, "n": copies}}. The best gene per anchor
    is the highest-confidence, highest-identity hit; `n` counts tandem copies.
    """
    buckets = defaultdict(list)
    for g in track.get("genes", []):
        if _is_goi_target_gene(g):
            buckets[goi_key].append(g)
            continue
        hid = g.get("home_gene_id", "") or g.get("name", "")
        if hid in anchor_keys:
            buckets[hid].append(g)
    out = {}
    for key, genes in buckets.items():
        genes.sort(key=lambda g: (_confidence_rank(g.get("confidence")),
                                   g.get("identity", 0.0)), reverse=True)
        out[key] = {"best": genes[0], "n": len(genes)}
    return out


# --- Species-tree (cladogram) panel + row ordering, shared by both grids -----

def _grid_collapsed_species_tree(tree_file):
    """Parse the GOI newick, midpoint-root it, and collapse to one leaf per
    species. Returns (rooted_tree_or_None, ordered_species_keys). A leaf's
    species key is ``synvoy_tree.species_from_leaf`` (e.g. 'zebrafish'), which
    matches a target track's ``genome_id`` so the tree aligns to the rows."""
    if not tree_file or not os.path.exists(tree_file):
        return None, []
    try:
        import synvoy_tree as _stree
        with open(tree_file) as fh:
            raw = _stree.parse_newick_tree(fh.read())
        rooted = _stree.midpoint_root(raw) if raw is not None else None
        if rooted is not None:
            rooted = _stree.collapse_to_one_leaf_per_species(
                rooted, _stree.species_from_leaf)
        order = []
        if rooted is not None:
            for lf in rooted.leaves():
                sp = _stree.species_from_leaf(lf.name) or lf.name
                if sp not in order:
                    order.append(sp)
        return rooted, order
    except Exception as exc:
        print(f"[plot] anchor-grid tree unavailable: {exc}", file=sys.stderr)
        return None, []


def _grid_species_tree(args):
    """Tree + species order for the grid cladograms.

    PREFERS a real NCBI-taxonomy SPECIES tree (built in ``main`` from the run's
    species list via ``synvoy_tree.build_taxonomy_tree_newick`` and stashed on
    ``args._species_tree_newick``) — the GOI ``--tree`` is a *gene* tree whose
    midpoint-rooted topology scrambles the species order (Apis species split
    apart, ants interleaved with bees), which is what made the cladogram look
    wrong. Falls back to the collapsed GOI gene tree only when taxonomy is
    unavailable. Returns (rooted_tree_or_None, ordered_species_keys, is_species).
    """
    nwk = getattr(args, "_species_tree_newick", None)
    if nwk:
        try:
            import synvoy_tree as _stree
            rooted = _stree.parse_newick_tree(nwk)  # NCBI topology is already rooted
            order = []
            if rooted is not None:
                for lf in rooted.leaves():
                    sp = lf.name
                    if sp and sp not in order:
                        order.append(sp)
            if rooted is not None and order:
                return rooted, order, True
        except Exception as exc:
            print(f"[plot] taxonomy species tree unavailable: {exc}", file=sys.stderr)
    rooted, order = _grid_collapsed_species_tree(getattr(args, "tree", None))
    return rooted, order, False


def _grid_match_species(genome_id, species_keys):
    """Best species key for a track genome_id (exact, else substring)."""
    if not genome_id:
        return None
    if genome_id in species_keys:
        return genome_id
    for sp in species_keys:
        if sp and (sp in genome_id or genome_id in sp):
            return sp
    return None


def _grid_order_targets(targets, species_order):
    """Stable-reorder target tracks to follow the tree's species order so the
    cladogram branches don't cross. Tracks whose species isn't in the tree keep
    their original relative order and trail the tree-ordered ones."""
    if not species_order:
        return list(targets)
    rank = {sp: i for i, sp in enumerate(species_order)}
    def _key(t):
        sp = _grid_match_species((t.get("genome_id") or ""), species_order)
        return (0, rank[sp]) if sp in rank else (1, 0)
    return sorted(targets, key=_key)


def _render_grid_tree_panel(rooted_tree, row_species, x0, x1, y_of_row,
                            col_line="#9aa3b2"):
    """Midpoint-rooted cladogram aligned to grid rows.

    row_species[ri] is the species key for row ri (None = home / unmapped);
    y_of_row(ri) returns row ri's vertical centre. Branches are right-angle
    ladder lines with x proportional to cumulative branch length.

    Crucially, EVERY grid row gets a connector. The ``--tree`` is the GOI
    *gene* tree, so its collapsed leaves only cover species in which the GOI
    was actually recovered — the home reference and any GOI-absent target have
    no leaf and used to be left with a blank tree column. Those *unplaced* rows
    are now attached to the root as a dashed basal polytomy (placement unknown,
    but visually connected), so the cladogram reads as covering all genomes.
    Returns a list of SVG strings (empty only when there's no tree at all)."""
    if rooted_tree is None or not list(rooted_tree.leaves()):
        return []
    import synvoy_tree as _stree
    sp_to_y = {sp: y_of_row(ri) for ri, sp in enumerate(row_species)
               if sp is not None}
    leaf_y = {}
    for lf in rooted_tree.leaves():
        sp = _stree.species_from_leaf(lf.name) or lf.name
        if sp in sp_to_y:
            leaf_y[id(lf)] = sp_to_y[sp]
    depth = {id(rooted_tree): 0.0}
    def _depth_walk(n):
        for c in n.children:
            depth[id(c)] = depth[id(n)] + (c.dist or 0.0)
            _depth_walk(c)
    _depth_walk(rooted_tree)
    max_depth = max(depth.values()) or 1.0
    span_x = (x1 - x0) - 6
    node_y = {}
    def _y_walk(n):
        if n.is_leaf():
            return [leaf_y[id(n)]] if id(n) in leaf_y else []
        ys = []
        for c in n.children:
            ys.extend(_y_walk(c))
        if ys:
            node_y[id(n)] = sum(ys) / len(ys)
        return ys
    _y_walk(rooted_tree)
    for lid, y in leaf_y.items():
        node_y[lid] = y
    def _x_for(n):
        return x0 + (depth[id(n)] / max_depth) * span_x
    out = ['<g class="grid-tree" fill="none">']
    def _draw(n):
        if n.is_leaf() or id(n) not in node_y:
            return
        nx = _x_for(n)
        child_ys = [node_y[id(c)] for c in n.children if id(c) in node_y]
        if child_ys:
            out.append(f'<line x1="{nx:.1f}" y1="{min(child_ys):.1f}" '
                       f'x2="{nx:.1f}" y2="{max(child_ys):.1f}" '
                       f'stroke="{col_line}" stroke-width="1.2"/>')
        for c in n.children:
            if id(c) not in node_y:
                continue
            cy, cx = node_y[id(c)], _x_for(c)
            out.append(f'<line x1="{nx:.1f}" y1="{cy:.1f}" x2="{cx:.1f}" '
                       f'y2="{cy:.1f}" stroke="{col_line}" stroke-width="1.2"/>')
            if c.is_leaf():
                # faint dotted leader from the branch tip to the grid edge so
                # the eye follows the row across the label gutter
                out.append(f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{x1:.1f}" '
                           f'y2="{cy:.1f}" stroke="{col_line}" stroke-width="0.7" '
                           f'stroke-dasharray="1,2" opacity="0.65"/>')
            _draw(c)
    _draw(rooted_tree)

    # ---- Attach UNPLACED rows (home + GOI-absent species) to the root ------
    # A row is "placed" when its species owns a tree leaf; everything else is
    # grafted onto the root as a dashed basal polytomy so no row is left blank.
    placed_species = {_stree.species_from_leaf(lf.name) or lf.name
                      for lf in rooted_tree.leaves() if id(lf) in leaf_y}
    unplaced_y = [y_of_row(ri) for ri, sp in enumerate(row_species)
                  if sp is None or sp not in placed_species]
    if unplaced_y:
        root_y = node_y.get(id(rooted_tree))
        if root_y is None:  # no row matched at all → synthesise a star root
            root_y = sum(unplaced_y) / len(unplaced_y)
        rx = _x_for(rooted_tree)
        vy0, vy1 = min([root_y] + unplaced_y), max([root_y] + unplaced_y)
        # extend the root's spine to span the unplaced rows (dashed)
        out.append(f'<line x1="{rx:.1f}" y1="{vy0:.1f}" x2="{rx:.1f}" '
                   f'y2="{vy1:.1f}" stroke="{col_line}" stroke-width="1.0" '
                   f'stroke-dasharray="3,2.5" opacity="0.7"/>')
        for uy in unplaced_y:
            out.append(f'<line x1="{rx:.1f}" y1="{uy:.1f}" x2="{x1:.1f}" '
                       f'y2="{uy:.1f}" stroke="{col_line}" stroke-width="1.0" '
                       f'stroke-dasharray="3,2.5" opacity="0.7"/>')

    out.append('</g>')
    return out


def render_anchor_grid(all_tracks, gene_colours, goi_genome_colours,
                       home_products, args, max_cols=50):
    """Render the anchor-grid (aligned-column) figure as self-contained HTML."""
    if not all_tracks:
        return _assemble_grid_html("", 400, 200, "SynVoy anchor grid")
    home_track = all_tracks[0]
    targets = all_tracks[1:]
    goi_key = "__GOI__"
    goi_label = _goi_display_label(all_tracks)

    # Species tree (collapsed to one leaf per genome) drives both the left
    # cladogram panel and the row order, so the branches read cleanly without
    # crossing. Falls back to no panel / original order when no tree is given.
    rooted_tree, species_order, _is_species_tree = _grid_species_tree(args)
    _tree_label = "Species tree" if _is_species_tree else "GOI phylogeny"
    targets = _grid_order_targets(targets, species_order)
    row_tracks = [home_track] + targets

    # ---- 1. Build ordered, de-duplicated anchors from the home track -------
    anchors = []          # list of dicts in home genomic order
    seen = set()
    for hg in sorted(home_track.get("genes", []), key=lambda g: g.get("start", 0)):
        nm = hg.get("name", "")
        is_g = is_goi(nm) or is_goi(hg.get("home_gene_id", "") or "")
        key = goi_key if is_g else nm
        if key in seen:
            continue
        seen.add(key)
        anchors.append({
            "key":     key,
            "is_goi":  is_g,
            "label":   goi_label if is_g else clean_gene_label(nm),
            "colour":  GOI_COLOUR if is_g else gene_colours.get(nm, UNMATCHED_CLR),
            "strand":  hg.get("strand", "+"),
            "start":   hg.get("start", 0),
            "end":     hg.get("end", 0),
            "product": _lookup_product(nm, home_products) if home_products else "",
        })
    anchor_keys = {a["key"] for a in anchors}

    # ---- 2. Map every target onto the anchor columns -----------------------
    target_maps = [_grid_target_map(t, anchor_keys, goi_key) for t in targets]

    # ---- 3. Focus + cap columns (GOI always kept) --------------------------
    # This is an ortholog-*alignment* view: a home gene recovered in no target
    # at all sits outside the aligned neighbourhood and only adds empty columns,
    # so drop it. Then cap the remainder by ortholog coverage.
    cov = {a["key"]: sum(1 for m in target_maps if a["key"] in m) for a in anchors}
    anchors = [a for a in anchors if a["is_goi"] or cov[a["key"]] > 0]
    if len(anchors) > max_cols:
        ordered = sorted(anchors, key=lambda a: (a["is_goi"], cov[a["key"]]), reverse=True)
        keep = {a["key"] for a in ordered[:max_cols]} | {goi_key}
        anchors = [a for a in anchors if a["key"] in keep]

    n_cols = len(anchors)
    n_rows = 1 + len(targets)

    # ---- 4. Layout geometry ------------------------------------------------
    GRID_BG = "#ffffff"
    GRID_HEADER = "#f5f7fb"
    GRID_FRAME = "#d8dee8"
    GRID_ROW_ALT = "#f8f9fc"
    GRID_ROW_HOME = "#eef2f7"
    GRID_AXIS = "#b6bfcd"
    GRID_AXIS_TEXT = "#6b7280"
    GRID_GAP = "#b45309"
    GRID_STRIPE = "#0f172a"

    # AXIS_H is deliberately tall: the home-coordinate axis carries slanted
    # leaders AND the "gap N Mb" rearrangement glyphs BELOW the axis line, and
    # those have to clear the top of the grid or the home row's opaque
    # background paints over them (the hidden-gap-label bug). Keeping the band
    # tall shifts the whole grid down so the gap annotation stays visible.
    MARGIN_TOP, HEADER_H, AXIS_H = 68, 130, 86
    ROW_H, ARROW_H = 42, 20
    COL_W, GOI_COL_W, RIGHT, LEGEND_H = 46, 68, 56, 132

    # labels follow the (reordered) row order: home first, then tree-sorted
    # targets — so labels[ri] always matches row_tracks[ri].
    labels = [_svg_esc(re.sub(r"<[^>]+>", "", (t.get("label") or "")).strip())
              for t in row_tracks]
    TREE_W = 144 if (rooted_tree is not None) else 0
    LABEL_W = max(190, min(360, int(max((len(l) for l in labels), default=12) * 7.0 + 40)))
    LEFT = TREE_W + LABEL_W

    col_x, col_w, cx = [], [], LEFT
    for a in anchors:
        w = GOI_COL_W if a["is_goi"] else COL_W
        col_x.append(cx); col_w.append(w); cx += w
    grid_x1 = cx
    grid_y0 = MARGIN_TOP + HEADER_H + AXIS_H
    grid_h = n_rows * ROW_H
    grid_w = grid_x1 - LEFT

    # Per-row genomic location ("chrom: lo–hi Mb"), shown in a right gutter like
    # the gene-position map — so the aligned grid also tells you WHERE each
    # genome's neighbourhood actually sits.
    home_chrom = next((g.get("chrom", "") for g in home_track.get("genes", [])), "")
    def _row_span(ri):
        if ri == 0:
            mids = [(a["start"] + a["end"]) / 2.0 for a in anchors if a["start"]]
            return home_chrom, (min(mids) if mids else 0), (max(mids) if mids else 0), 0
        tmap = target_maps[ri - 1]
        items, goi_chrom = [], ""
        for a in anchors:
            entry = tmap.get(a["key"])
            if not entry:
                continue
            gs, ge = _get_coords(entry["best"])
            ch = entry["best"].get("chrom", "")
            items.append((ch, (gs + ge) / 2.0))
            if a["is_goi"] and ch:
                goi_chrom = ch
        return _dominant_chrom_span(items, goi_chrom)
    def _span_text(ri):
        c, lo, hi, n_other = _row_span(ri)
        return _span_gutter_label(c, lo, hi, n_other)
    span_texts = [_span_text(ri) for ri in range(n_rows)]
    _max_span = max((len(s) for s in span_texts), default=0)
    SPAN_W = (max(150, int(_max_span * 6.0 + 26)) if _max_span else 10)

    # The caption + legend overflowed the canvas on narrow, few-column plots
    # (e.g. melittin). Size the canvas to the widest of the grid+span gutter, the
    # caption line, and the legend row so nothing spills outside.
    _caption = ("Columns = home gene order   ·   Rows = species "
                "(phylogenetic order)   ·   arrow points in coding strand")
    _legend_row_w = 560  # widest legend row (confidence tiers + GOI + no-ortholog)
    content_w = max(grid_w + SPAN_W, len(_caption) * 5.7, _legend_row_w)
    width = LEFT + content_w + RIGHT
    height = grid_y0 + grid_h + LEGEND_H
    header_y0 = MARGIN_TOP - 12
    header_h = HEADER_H + AXIS_H + 12

    def col_center(i):
        return col_x[i] + col_w[i] / 2

    P = []  # svg parts

    # ---- 5. Header backplate (behind axis + column labels) -----------------
    P.append(f'<rect x="{LEFT - 8:.1f}" y="{header_y0:.1f}" '
             f'width="{grid_w + 16:.1f}" height="{header_h:.1f}" '
             f'fill="{GRID_HEADER}" stroke="{GRID_FRAME}" stroke-width="1" rx="8"/>')

    # ---- 6. GOI guide band + column striping (drawn first, behind) ---------
    for i, a in enumerate(anchors):
        if a["is_goi"]:
            P.append(f'<rect x="{col_x[i]:.1f}" y="{grid_y0 - 6:.1f}" '
                     f'width="{col_w[i]:.1f}" height="{n_rows*ROW_H + 12:.1f}" '
                     f'fill="{GOI_COLOUR}" opacity="0.08" '
                     f'stroke="{GOI_BORDER}" stroke-width="0.6" stroke-opacity="0.3"/>')
        elif i % 2 == 0:
            P.append(f'<rect x="{col_x[i]:.1f}" y="{grid_y0:.1f}" '
                     f'width="{col_w[i]:.1f}" height="{n_rows*ROW_H:.1f}" '
                     f'fill="{GRID_STRIPE}" opacity="0.015"/>')

    # ---- 7. Home-coordinate axis (real positions + gap glyphs) -------------
    mids = [(a["start"] + a["end"]) / 2.0 for a in anchors if a["start"]]
    if len(mids) >= 2 and max(mids) > min(mids):
        pmin, pmax = min(mids), max(mids)
        span = pmax - pmin
        axis_y = MARGIN_TOP + HEADER_H + AXIS_H * 0.62
        gx0, gx1 = LEFT, grid_x1
        def real_x(pos):
            return gx0 + (pos - pmin) / span * (gx1 - gx0)
        P.append(f'<line x1="{gx0:.1f}" y1="{axis_y:.1f}" x2="{gx1:.1f}" '
                 f'y2="{axis_y:.1f}" stroke="{GRID_AXIS}" stroke-width="1.2"/>')
        P.append(f'<text x="{gx0:.1f}" y="{axis_y - 7:.1f}" font-size="9" '
                 f'fill="{GRID_AXIS_TEXT}">{pmin/1e6:.2f} Mb</text>')
        P.append(f'<text x="{gx1:.1f}" y="{axis_y - 7:.1f}" font-size="9" '
                 f'fill="{GRID_AXIS_TEXT}" text-anchor="end">{pmax/1e6:.2f} Mb</text>')
        gap_thresh = max(1.0e6, span * 0.18)
        for i, a in enumerate(anchors):
            if not a["start"]:
                continue
            rx = real_x((a["start"] + a["end"]) / 2.0)
            cc = col_center(i)
            col = GOI_COLOUR if a["is_goi"] else GRID_AXIS
            # slanted leader from true position down to its aligned column
            P.append(f'<path d="M{rx:.1f},{axis_y:.1f} L{cc:.1f},{grid_y0 - 4:.1f}" '
                     f'stroke="{col}" stroke-width="{1.4 if a["is_goi"] else 0.7:.1f}" '
                     f'fill="none" opacity="{0.85 if a["is_goi"] else 0.5}"/>')
            P.append(f'<circle cx="{rx:.1f}" cy="{axis_y:.1f}" '
                     f'r="{3.2 if a["is_goi"] else 2.1:.1f}" fill="{col}"/>')
        # annotate large genomic gaps between consecutive anchors
        for j in range(len(anchors) - 1):
            if not anchors[j]["start"] or not anchors[j + 1]["start"]:
                continue
            m0 = (anchors[j]["start"] + anchors[j]["end"]) / 2.0
            m1 = (anchors[j + 1]["start"] + anchors[j + 1]["end"]) / 2.0
            if m1 - m0 > gap_thresh:
                gxm = (real_x(m0) + real_x(m1)) / 2.0
                gap_label = f"gap {(m1 - m0)/1e6:.1f} Mb"
                slash_y0 = axis_y + 4
                slash_y1 = axis_y + 12
                P.append(f'<line x1="{gxm - 6:.1f}" y1="{slash_y0:.1f}" '
                         f'x2="{gxm - 1:.1f}" y2="{slash_y1:.1f}" '
                         f'stroke="{GRID_GAP}" stroke-width="1.2" stroke-linecap="round"/>')
                P.append(f'<line x1="{gxm + 1:.1f}" y1="{slash_y0:.1f}" '
                         f'x2="{gxm + 6:.1f}" y2="{slash_y1:.1f}" '
                         f'stroke="{GRID_GAP}" stroke-width="1.2" stroke-linecap="round"/>')
                gap_lbl_y = min(axis_y + 22, grid_y0 - 10)
                P.append(f'<text x="{gxm:.1f}" y="{gap_lbl_y:.1f}" font-size="8.5" '
                         f'fill="{GRID_GAP}" text-anchor="middle" font-weight="600">'
                         f'{gap_label}</text>')

    # ---- 8. Column header labels (rotated) ---------------------------------
    for i, a in enumerate(anchors):
        cc = col_center(i)
        ly = MARGIN_TOP + HEADER_H - 6
        cls = "goi" if a["is_goi"] else ""
        fill = GOI_BORDER if a["is_goi"] else "#42495a"
        weight = "700" if a["is_goi"] else "500"
        P.append(f'<text class="acol-lbl {cls}" x="{cc:.1f}" y="{ly:.1f}" '
                 f'font-size="10.5" fill="{fill}" font-weight="{weight}" '
                 f'transform="rotate(-55 {cc:.1f} {ly:.1f})">{_svg_esc(a["label"])}</text>')

    # ---- 9. Rows: home first, then targets ---------------------------------
    def draw_arrow(i, ri, base, strand, identity, conf, n_copies, inverted,
                   is_home, title, coverage=None):
        x0 = col_x[i] + 4
        x1 = col_x[i] + col_w[i] - 4
        yb = grid_y0 + ri * ROW_H + (ROW_H - ARROW_H) / 2
        conf = (conf or "").upper()
        if is_home:
            fill, stroke, dash, op = base, _darken_hex(base, 0.7), "", 0.95
        elif conf == "LOW":
            fill, stroke, dash, op = _lerp_hex(base, "#ffffff", 0.6), base, ' stroke-dasharray="2,2"', 0.85
        elif conf == "MEDIUM":
            fill, stroke, dash, op = _lerp_hex(base, "#ffffff", 0.42), base, ' stroke-dasharray="3.5,2"', 0.95
        else:  # HIGH / unknown
            fill, stroke, dash, op = _shade_by_identity(base, identity), _darken_hex(base, 0.65), "", 1.0
        d = _svg_arrow_path(x0, x1, yb, ARROW_H, strand)
        inner = [f'<title>{_svg_esc(title)}</title>',
                 f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="1"{dash} '
                 f'opacity="{op:.2f}"/>']
        if (not is_home) and identity >= 25:
            # Default to the identity number ALONE — one clean value that always
            # fits the cell. Query-coverage is appended (in parens) only when it
            # is low, since a short high-identity hit would otherwise read as a
            # full-length ortholog. Dark text with a white halo (paint-order:
            # stroke) reads on any fill shade; the font auto-shrinks so the
            # label can never spill past the arrow body (the old fixed size
            # clipped "100/100" to "00/100" in the narrow GOI column).
            if coverage is not None and coverage < COVERAGE_FLAG_THRESHOLD:
                num = f"{identity:.0f} ({coverage*100:.0f})"
            else:
                num = f"{identity:.0f}"
            avail = (x1 - x0) - 5
            fsz = max(6.5, min(10.5, avail / (len(num) * 0.60)))
            inner.append(
                f'<text x="{(x0+x1)/2:.1f}" y="{yb + ARROW_H/2 + 3.2:.1f}" '
                f'text-anchor="middle" font-size="{fsz:.1f}" fill="#15181f" '
                f'font-weight="{"700" if conf=="HIGH" else "600"}" '
                f'pointer-events="none" '
                f'style="paint-order:stroke;stroke:#ffffff;stroke-width:1.4">'
                f'{num}</text>'
            )
        if n_copies > 1:
            # Copy-count badge sits ABOVE the arrow (in the row's top padding),
            # clear of the identity number AND of the home row above it, so the
            # first target row's "×N" stays readable instead of tucking under
            # the home genome. White halo keeps it legible over the frame line.
            inner.append(
                f'<text x="{x1:.1f}" y="{yb - 3:.1f}" text-anchor="end" '
                f'font-size="8.5" fill="{_darken_hex(base, 0.55)}" font-weight="700" '
                f'pointer-events="none" '
                f'style="paint-order:stroke;stroke:#ffffff;stroke-width:1.8">'
                f'×{n_copies}</text>')
        P.append('<g class="acell">' + "".join(inner) + '</g>')

    def draw_absent(i, ri):
        cc = col_center(i)
        ymid = grid_y0 + ri * ROW_H + ROW_H / 2
        P.append(f'<circle cx="{cc:.1f}" cy="{ymid:.1f}" r="4.4" '
             f'fill="none" stroke="#d7dbe3" stroke-width="1" '
             f'stroke-dasharray="2,2"/>')

    for ri in range(n_rows):
        is_home = ri == 0
        track = home_track if is_home else targets[ri - 1]
        tmap = None if is_home else target_maps[ri - 1]
        row_y = grid_y0 + ri * ROW_H
        # row background
        if is_home:
            row_fill = GRID_ROW_HOME
        elif ri % 2 == 1:
            row_fill = GRID_ROW_ALT
        else:
            row_fill = GRID_BG
        P.append(f'<rect class="arow-bg" x="{LEFT:.1f}" y="{row_y:.1f}" '
                 f'width="{grid_w + SPAN_W:.1f}" height="{ROW_H:.1f}" '
                 f'fill="{row_fill}" opacity="1.0"/>')
        # genomic-location label in the right gutter ("chrom: lo–hi Mb")
        if span_texts[ri]:
            P.append(f'<text x="{grid_x1 + 10:.1f}" y="{row_y + ROW_H / 2 + 3:.1f}" '
                     f'font-size="9" fill="{GRID_AXIS_TEXT}">'
                     f'{_svg_esc(span_texts[ri])}</text>')
        # collinear thread between present columns
        present_cc = [col_center(i) for i, a in enumerate(anchors)
                      if is_home or (tmap and a["key"] in tmap)]
        if len(present_cc) >= 2:
            ymid = row_y + ROW_H / 2
            P.append(f'<line x1="{min(present_cc):.1f}" y1="{ymid:.1f}" '
                     f'x2="{max(present_cc):.1f}" y2="{ymid:.1f}" '
                     f'stroke="{GRID_FRAME}" stroke-width="1.2"/>')
        # species / row label (clade swatch removed — clade is conveyed by the
        # tree panel; the coloured left bar was redundant and distracting)
        lbl = labels[ri]
        if is_home:
            lbl = lbl or "Home"
        ly = row_y + ROW_H / 2 + 4
        P.append(f'<text class="arow-lbl" x="{TREE_W + 8:.1f}" y="{ly:.1f}" '
                 f'font-size="12" fill="#1a1d26" '
                 f'font-weight="{"700" if is_home else "600"}">'
                 f'{lbl}</text>')
        # cells
        for i, a in enumerate(anchors):
            if is_home:
                title = f'{a["label"]} (home reference) — {a["product"] or a["key"]}'
                draw_arrow(i, ri, a["colour"], a["strand"], 0.0, "", 1, False, True, title)
                continue
            entry = tmap.get(a["key"]) if tmap else None
            if not entry:
                draw_absent(i, ri)
                continue
            g = entry["best"]
            ident = g.get("identity", 0.0)
            conf = g.get("confidence", "")
            cov = g.get("query_coverage")
            tstrand = g.get("strand", "+")
            inverted = bool(tstrand) and bool(a["strand"]) and tstrand != a["strand"]
            tgt_label = clean_gene_label(_preferred_target_label(g))
            cov_txt = f", coverage {cov*100:.0f}%" if cov is not None else ""
            title = (f'{a["label"]} in {lbl} — identity {ident:.1f}%{cov_txt}, '
                     f'confidence {conf or "—"}'
                     + (", inverted" if inverted else "")
                     + (f', {entry["n"]}× copies' if entry["n"] > 1 else "")
                     + (f' — {tgt_label}' if tgt_label else ""))
            draw_arrow(i, ri, a["colour"], tstrand, ident, conf, entry["n"],
                       inverted, False, title, coverage=cov)

    # ---- 9b. Species-tree (cladogram) panel on the left --------------------
    if TREE_W:
        row_species = [None] + [
            _grid_match_species((t.get("genome_id") or ""), species_order)
            for t in targets]
        P.extend(_render_grid_tree_panel(
            rooted_tree, row_species, 12, TREE_W - 6,
            lambda ri: grid_y0 + ri * ROW_H + ROW_H / 2))
        P.append(f'<text x="{TREE_W / 2:.1f}" y="{grid_y0 - 9:.1f}" '
                 f'text-anchor="middle" font-size="9.5" fill="{GRID_AXIS_TEXT}" '
                 f'font-style="italic">{_tree_label}</text>')

    # ---- 10. Legend (two rows: keys on row 1, the long metric note on row 2,
    # so a wide note never runs off the right edge; wider symbol→text gaps) ---
    ly0 = grid_y0 + n_rows * ROW_H + 34
    lx = LEFT
    ARROW_W, SYM_GAP = 26, 9   # SYM_GAP: clear space between a symbol and its label
    P.append(f'<text x="{lx:.1f}" y="{ly0 - 16:.1f}" font-size="10.5" fill="#8c95a6">'
             f'{_svg_esc(_caption)}</text>')
    # Row 1 — confidence tiers (HIGH/MEDIUM/LOW arrow styles) + GOI + no-ortholog.
    cur = lx
    P.append(f'<text x="{cur:.1f}" y="{ly0 + 11:.1f}" font-size="10" '
             f'font-weight="700" fill="#42495a">Confidence:</text>')
    cur += 72
    for name, dash, fill in [
        ("HIGH", "", _shade_by_identity(GENE_PALETTE[0], 95)),
        ("MEDIUM", ' stroke-dasharray="3.5,2"', _lerp_hex(GENE_PALETTE[0], "#ffffff", 0.42)),
        ("LOW", ' stroke-dasharray="2,2"', _lerp_hex(GENE_PALETTE[0], "#ffffff", 0.6)),
    ]:
        d = _svg_arrow_path(cur, cur + ARROW_W, ly0, 14, "+")
        P.append(f'<path d="{d}" fill="{fill}" stroke="{_darken_hex(GENE_PALETTE[0],0.65)}" '
                 f'stroke-width="1"{dash}/>')
        P.append(f'<text x="{cur + ARROW_W + SYM_GAP:.1f}" y="{ly0 + 11:.1f}" font-size="10" '
                 f'fill="#42495a">{name}</text>')
        cur += ARROW_W + SYM_GAP + len(name) * 7 + 24
    # GOI swatch (a category, not a confidence tier).
    cur += 12
    d = _svg_arrow_path(cur, cur + ARROW_W, ly0, 14, "+")
    P.append(f'<path d="{d}" fill="{_shade_by_identity(GOI_COLOUR, 95)}" '
             f'stroke="{GOI_BORDER}" stroke-width="1"/>')
    P.append(f'<text x="{cur + ARROW_W + SYM_GAP:.1f}" y="{ly0 + 11:.1f}" font-size="10" '
             f'fill="#42495a">GOI</text>')
    cur += ARROW_W + SYM_GAP + 3 * 7 + 24
    # absent marker
    P.append(f'<circle cx="{cur + 5:.1f}" cy="{ly0 + 5:.1f}" r="4.5" '
             f'fill="none" stroke="#d7dbe3" stroke-dasharray="2,2"/>')
    P.append(f'<text x="{cur + 9.5 + SYM_GAP:.1f}" y="{ly0 + 9:.1f}" font-size="10" '
             f'fill="#42495a">no ortholog</text>')
    # Row 2 — metric note: each cell shows "%identity / %query-coverage"; shade = identity.
    P.append(f'<text x="{lx:.1f}" y="{ly0 + 33:.1f}" font-size="10" '
             f'fill="#8c95a6">number = % identity &#183; (n) = % coverage when '
             f'&lt; 80% &#183; shade = % identity</text>')

    # ---- 11. Title ---------------------------------------------------------
    title_txt = f"Anchor-grid synteny - GOI: {_svg_esc(goi_label)}"
    P.insert(0, f'<text class="grid-title" x="{LEFT:.1f}" y="30" font-size="17" '
                f'font-weight="700" fill="#1a1d26">{title_txt}</text>')
    P.insert(1, f'<text class="grid-subtitle" x="{LEFT:.1f}" y="48" font-size="11" '
                f'fill="{GRID_AXIS_TEXT}">'
                f'{n_cols} anchor genes x {n_rows} genomes | '
                f'orthologues aligned into shared columns</text>')
    # Explicit white canvas so the static-SVG/PNG export isn't transparent
    # (renders black) in viewers that ignore the CSS `background` property.
    P.insert(0, f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')

    # Frame around the grid area (drawn last so it sits above row fills).
    P.append(f'<rect x="{LEFT:.1f}" y="{grid_y0:.1f}" '
             f'width="{grid_w:.1f}" height="{grid_h:.1f}" '
             f'fill="none" stroke="{GRID_FRAME}" stroke-width="1"/>')

    return _assemble_grid_html("\n".join(P), width, height, "SynVoy anchor grid")


def render_anchor_grid_positional(all_tracks, gene_colours, goi_genome_colours,
                                  home_products, args):
    """Real-position variant of the anchor grid.

    Companion to :func:`render_anchor_grid`. Instead of evenly-spaced aligned
    columns, every genome is a horizontal *backbone* line carrying one dot at
    each gene's true (per-row normalised) position, coloured by its home
    ortholog. The aligned grid answers "is the ortholog present?"; this view
    answers "where does it actually sit, and is the neighbourhood
    rearranged?" — spacing, gaps and inversions the column grid hides. The
    species-tree panel and palette are shared with the aligned grid."""
    if not all_tracks:
        return _assemble_grid_html("", 400, 200, "SynVoy gene-position map")
    home_track = all_tracks[0]
    targets = all_tracks[1:]
    goi_key = "__GOI__"
    goi_label = _goi_display_label(all_tracks)

    rooted_tree, species_order, _is_species_tree = _grid_species_tree(args)
    _tree_label = "Species tree" if _is_species_tree else "GOI phylogeny"
    targets = _grid_order_targets(targets, species_order)
    row_tracks = [home_track] + targets

    # Home anchors → colour/label per key (same homology palette as the grid).
    anchors, seen = [], set()
    for hg in sorted(home_track.get("genes", []), key=lambda g: g.get("start", 0)):
        nm = hg.get("name", "")
        is_g = is_goi(nm) or is_goi(hg.get("home_gene_id", "") or "")
        key = goi_key if is_g else nm
        if key in seen:
            continue
        seen.add(key)
        anchors.append({"key": key, "is_goi": is_g,
                        "colour": GOI_COLOUR if is_g else gene_colours.get(nm, UNMATCHED_CLR),
                        "start": hg.get("start", 0), "end": hg.get("end", 0),
                        "strand": hg.get("strand", "+"),
                        "label": goi_label if is_g else clean_gene_label(nm)})
    anchor_keys = {a["key"] for a in anchors}
    anchor_colour = {a["key"]: a["colour"] for a in anchors}
    anchor_label = {a["key"]: a["label"] for a in anchors}
    anchor_strand = {a["key"]: a["strand"] for a in anchors}
    target_maps = [_grid_target_map(t, anchor_keys, goi_key) for t in targets]
    home_chrom = next((g.get("chrom", "") for g in home_track.get("genes", [])), "")

    # Per-row points: each is (x_real_bp, key, colour, is_goi, ident, conf,
    # strand, n_copies, label, chrom).
    def _row_points(ri):
        pts = []
        if ri == 0:
            for a in anchors:
                if not a["start"]:
                    continue
                pts.append({"x": (a["start"] + a["end"]) / 2.0, "key": a["key"],
                            "colour": a["colour"], "is_goi": a["is_goi"],
                            "ident": 100.0, "conf": "HIGH", "strand": a["strand"],
                            "n": 1, "label": a["label"], "chrom": home_chrom})
            return pts
        tmap = target_maps[ri - 1]
        for key, entry in tmap.items():
            g = entry["best"]
            gs, ge = _get_coords(g)
            is_g = (key == goi_key)
            pts.append({"x": (gs + ge) / 2.0, "key": key,
                        "colour": GOI_COLOUR if is_g else anchor_colour.get(key, UNMATCHED_CLR),
                        "is_goi": is_g, "ident": g.get("identity", 0.0),
                        "conf": (g.get("confidence") or "").upper(),
                        "strand": g.get("strand", "+"), "n": entry["n"],
                        "label": anchor_label.get(key, key) if not is_g else goi_label,
                        "chrom": g.get("chrom", "")})
        return pts

    row_pts = [_row_points(ri) for ri in range(len(row_tracks))]

    # ---- Layout -----------------------------------------------------------
    GRID_HEADER = "#f5f7fb"
    GRID_FRAME = "#d8dee8"
    GRID_ROW_ALT = "#f8f9fc"
    GRID_ROW_HOME = "#eef2f7"
    GRID_AXIS_TEXT = "#6b7280"
    GRID_BACKBONE = "#c2cad6"

    MARGIN_TOP, HEADER_H = 64, 30
    ROW_H, DOT_R, GOI_R, RIGHT, LEGEND_H = 40, 4.6, 6.4, 64, 96
    PLOT_W = 740

    # Right gutter for the per-row "chrom: lo–hi Mb" span label. Sized to the
    # LONGEST actual label (accession + range) so long scaffold/contig names
    # (e.g. 'WUUM01000001.1: 9.44–9.58 Mb') aren't clipped at the canvas edge.
    def _span_label_len(ri):
        pts = row_pts[ri]
        if not pts:
            return 0
        goi_ch = next((p["chrom"] for p in pts if p["is_goi"] and p["chrom"]), "")
        c, lo, hi, n_other = _dominant_chrom_span(
            [(p["chrom"], p["x"]) for p in pts], goi_ch)
        return len(_span_gutter_label(c, lo, hi, n_other))
    max_span_len = max((_span_label_len(ri) for ri in range(len(row_tracks))),
                       default=14)
    SPAN_W = max(150, int(max_span_len * 6.0 + 30))

    labels = [_svg_esc(re.sub(r"<[^>]+>", "", (t.get("label") or "")).strip())
              for t in row_tracks]
    TREE_W = 144 if (rooted_tree is not None) else 0
    LABEL_W = max(190, min(360, int(max((len(l) for l in labels), default=12) * 7.0 + 40)))
    LEFT = TREE_W + LABEL_W
    n_rows = len(row_tracks)

    grid_x0 = LEFT
    grid_x1 = LEFT + PLOT_W
    grid_y0 = MARGIN_TOP + HEADER_H
    grid_h = n_rows * ROW_H
    inset = 16  # keep edge dots off the frame
    bx0, bx1 = grid_x0 + inset, grid_x1 - inset

    _caption = ("Each genome normalised to its own neighbourhood   ·   dot = gene "
                "(colour = home ortholog)   ·   ◆ = GOI   ·   ring = inverted vs home")
    content_w = max(PLOT_W + SPAN_W, len(_caption) * 5.6, 640)
    width = LEFT + content_w + RIGHT
    height = grid_y0 + grid_h + LEGEND_H

    P = []
    P.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')

    # Header backplate over the plot band.
    P.append(f'<rect x="{grid_x0 - 8:.1f}" y="{MARGIN_TOP - 10:.1f}" '
             f'width="{PLOT_W + SPAN_W + 16:.1f}" height="{HEADER_H + grid_h + 10:.1f}" '
             f'fill="{GRID_HEADER}" stroke="{GRID_FRAME}" stroke-width="1" rx="8"/>')

    # Rows.
    for ri in range(n_rows):
        is_home = ri == 0
        track = row_tracks[ri]
        row_y = grid_y0 + ri * ROW_H
        cy = row_y + ROW_H / 2
        row_fill = GRID_ROW_HOME if is_home else (GRID_ROW_ALT if ri % 2 == 1 else "#ffffff")
        P.append(f'<rect x="{grid_x0:.1f}" y="{row_y:.1f}" width="{PLOT_W + SPAN_W:.1f}" '
                 f'height="{ROW_H:.1f}" fill="{row_fill}"/>')

        # Row label (clade swatch removed — clade is read from the tree panel).
        P.append(f'<text x="{TREE_W + 8:.1f}" y="{cy + 4:.1f}" font-size="12" '
                 f'fill="#1a1d26" font-weight="{"700" if is_home else "600"}">'
                 f'{labels[ri]}</text>')

        pts = sorted(row_pts[ri], key=lambda p: p["x"])
        if not pts:
            P.append(f'<text x="{(bx0 + bx1) / 2:.1f}" y="{cy + 4:.1f}" '
                     f'text-anchor="middle" font-size="10" fill="#aab2c0" '
                     f'font-style="italic">no orthologs mapped</text>')
            continue
        # On a fragmented assembly a row's genes can straddle >1 scaffold; mixing
        # their coordinates on one normalised line is meaningless, so keep the
        # dominant (or GOI-bearing) scaffold and note any dropped ones as "(+N)".
        goi_ch = next((p["chrom"] for p in pts if p["is_goi"] and p["chrom"]), "")
        prim, _plo, _phi, n_other = _dominant_chrom_span(
            [(p["chrom"], p["x"]) for p in pts], goi_ch)
        if prim:
            pts = [p for p in pts if p["chrom"] == prim] or pts
        xs = [p["x"] for p in pts]
        lo, hi = min(xs), max(xs)
        def _px(xr):
            if hi <= lo:
                return (bx0 + bx1) / 2
            return bx0 + (xr - lo) / (hi - lo) * (bx1 - bx0)

        # Backbone line through the gene dots.
        P.append(f'<line x1="{_px(lo):.1f}" y1="{cy:.1f}" x2="{_px(hi):.1f}" '
                 f'y2="{cy:.1f}" stroke="{GRID_BACKBONE}" stroke-width="1.6"/>')

        # Coordinate-span label in the right gutter.
        span_txt = _span_gutter_label(prim, lo, hi, n_other)
        if span_txt:
            P.append(f'<text x="{grid_x1 + 10:.1f}" y="{cy + 3.5:.1f}" font-size="9" '
                     f'fill="{GRID_AXIS_TEXT}">{_svg_esc(span_txt)}</text>')

        # Dots (GOI last so it sits on top).
        for p in sorted(pts, key=lambda q: q["is_goi"]):
            x = _px(p["x"])
            conf = p["conf"]
            inverted = (not p["is_goi"] and p["strand"] and
                        anchor_strand.get(p["key"]) and
                        p["strand"] != anchor_strand.get(p["key"]))
            if is_home:
                fill = p["colour"]
            elif conf == "LOW":
                fill = _lerp_hex(p["colour"], "#ffffff", 0.55)
            elif conf == "MEDIUM":
                fill = _lerp_hex(p["colour"], "#ffffff", 0.30)
            else:
                fill = _shade_by_identity(p["colour"], p["ident"])
            stroke = _darken_hex(p["colour"], 0.6)
            tip = (f'{p["label"]}'
                   + ('' if is_home else f' — identity {p["ident"]:.0f}%'
                      + (f', {conf}' if conf else '')
                      + (', inverted' if inverted else '')
                      + (f', {p["n"]}× copies' if p["n"] > 1 else '')))
            g = [f'<title>{_svg_esc(tip)}</title>']
            if p["is_goi"]:
                r = GOI_R
                g.append(f'<path d="M{x:.1f},{cy - r:.1f} L{x + r:.1f},{cy:.1f} '
                         f'L{x:.1f},{cy + r:.1f} L{x - r:.1f},{cy:.1f} Z" '
                         f'fill="{fill}" stroke="{GOI_BORDER}" stroke-width="1.3"/>')
            else:
                g.append(f'<circle cx="{x:.1f}" cy="{cy:.1f}" r="{DOT_R}" '
                         f'fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
                if inverted:
                    g.append(f'<circle cx="{x:.1f}" cy="{cy:.1f}" r="{DOT_R + 2.4}" '
                             f'fill="none" stroke="{stroke}" stroke-width="0.9" '
                             f'stroke-dasharray="1.5,1.5"/>')
            if p["n"] > 1:
                g.append(f'<text x="{x + DOT_R + 1:.1f}" y="{cy - DOT_R - 1:.1f}" '
                         f'font-size="7.5" fill="{stroke}" font-weight="700" '
                         f'style="paint-order:stroke;stroke:#fff;stroke-width:1.6">'
                         f'×{p["n"]}</text>')
            P.append('<g class="acell">' + "".join(g) + '</g>')

    # Tree panel.
    if TREE_W:
        row_species = [None] + [
            _grid_match_species((t.get("genome_id") or ""), species_order)
            for t in targets]
        P.extend(_render_grid_tree_panel(
            rooted_tree, row_species, 12, TREE_W - 6,
            lambda ri: grid_y0 + ri * ROW_H + ROW_H / 2))
        P.append(f'<text x="{TREE_W / 2:.1f}" y="{grid_y0 - 6:.1f}" '
                 f'text-anchor="middle" font-size="9.5" fill="{GRID_AXIS_TEXT}" '
                 f'font-style="italic">{_tree_label}</text>')

    # Frame.
    P.append(f'<rect x="{grid_x0:.1f}" y="{grid_y0:.1f}" width="{PLOT_W + SPAN_W:.1f}" '
             f'height="{grid_h:.1f}" fill="none" stroke="{GRID_FRAME}" stroke-width="1"/>')

    # Caption + legend.
    ly0 = grid_y0 + grid_h + 30
    P.append(f'<text x="{LEFT:.1f}" y="{ly0 - 12:.1f}" font-size="10.5" '
             f'fill="#8c95a6">{_svg_esc(_caption)}</text>')
    cur = LEFT
    P.append(f'<circle cx="{cur + 6:.1f}" cy="{ly0 + 6:.1f}" r="{DOT_R}" '
             f'fill="{_shade_by_identity(GENE_PALETTE[0], 95)}" '
             f'stroke="{_darken_hex(GENE_PALETTE[0], 0.6)}" stroke-width="1"/>')
    P.append(f'<text x="{cur + 20:.1f}" y="{ly0 + 9:.1f}" font-size="10" '
             f'fill="#42495a">flanking ortholog</text>')
    cur += 20 + len("flanking ortholog") * 6 + 26
    P.append(f'<path d="M{cur + 6:.1f},{ly0 + 6 - GOI_R:.1f} L{cur + 6 + GOI_R:.1f},'
             f'{ly0 + 6:.1f} L{cur + 6:.1f},{ly0 + 6 + GOI_R:.1f} '
             f'L{cur + 6 - GOI_R:.1f},{ly0 + 6:.1f} Z" '
             f'fill="{_shade_by_identity(GOI_COLOUR, 95)}" stroke="{GOI_BORDER}" '
             f'stroke-width="1.2"/>')
    P.append(f'<text x="{cur + 22:.1f}" y="{ly0 + 9:.1f}" font-size="10" '
             f'fill="#42495a">GOI</text>')
    cur += 22 + 3 * 7 + 26
    P.append(f'<circle cx="{cur + 6:.1f}" cy="{ly0 + 6:.1f}" r="{DOT_R}" '
             f'fill="{_lerp_hex(GENE_PALETTE[0], "#ffffff", 0.5)}" '
             f'stroke="{_darken_hex(GENE_PALETTE[0], 0.6)}" stroke-width="1"/>')
    P.append(f'<circle cx="{cur + 6:.1f}" cy="{ly0 + 6:.1f}" r="{DOT_R + 2.4}" '
             f'fill="none" stroke="{_darken_hex(GENE_PALETTE[0], 0.6)}" '
             f'stroke-width="0.9" stroke-dasharray="1.5,1.5"/>')
    P.append(f'<text x="{cur + 22:.1f}" y="{ly0 + 9:.1f}" font-size="10" '
             f'fill="#42495a">inverted vs home   ·   shade = % identity</text>')

    # Title.
    title_txt = f"Gene-position map - GOI: {_svg_esc(goi_label)}"
    P.insert(1, f'<text class="grid-title" x="{LEFT:.1f}" y="30" font-size="17" '
                f'font-weight="700" fill="#1a1d26">{title_txt}</text>')
    P.insert(2, f'<text class="grid-subtitle" x="{LEFT:.1f}" y="48" font-size="11" '
                f'fill="{GRID_AXIS_TEXT}">{len(anchors)} home genes x {n_rows} '
                f'genomes | true genomic positions per row</text>')

    return _assemble_grid_html("\n".join(P), width, height,
                               "SynVoy gene-position map")


def render_anchor_grid_threaded(all_tracks, gene_colours, goi_genome_colours,
                                home_products, args, max_cols=50):
    """Anchor grid where EVERY genome carries its own light real-position dot-line.

    The aligned-column anchor grid (:func:`render_anchor_grid`) answers "is the
    ortholog present, and in what order"; here each row ALSO gets the home-axis
    idiom — a faint line of dots at the genes' TRUE positions with slanted
    leaders down to their aligned column — so the real spacing and any
    rearrangement of every genome is shown unobtrusively right above its arrow
    row. Static SVG; shares the column layout, palette and tree panel with the
    aligned grid."""
    if not all_tracks:
        return _assemble_grid_html("", 400, 200, "SynVoy anchor positions")
    home_track = all_tracks[0]
    targets = all_tracks[1:]
    goi_key = "__GOI__"
    goi_label = _goi_display_label(all_tracks)

    rooted_tree, species_order, _is_species_tree = _grid_species_tree(args)
    _tree_label = "Species tree" if _is_species_tree else "GOI phylogeny"
    targets = _grid_order_targets(targets, species_order)

    # ---- anchors (home order, dedup, GOI keyed) ----
    anchors, seen = [], set()
    for hg in sorted(home_track.get("genes", []), key=lambda g: g.get("start", 0)):
        nm = hg.get("name", "")
        is_g = is_goi(nm) or is_goi(hg.get("home_gene_id", "") or "")
        key = goi_key if is_g else nm
        if key in seen:
            continue
        seen.add(key)
        anchors.append({"key": key, "is_goi": is_g,
                        "label": goi_label if is_g else clean_gene_label(nm),
                        "colour": GOI_COLOUR if is_g else gene_colours.get(nm, UNMATCHED_CLR),
                        "strand": hg.get("strand", "+"),
                        "start": hg.get("start", 0), "end": hg.get("end", 0)})
    anchor_keys = {a["key"] for a in anchors}
    target_maps = [_grid_target_map(t, anchor_keys, goi_key) for t in targets]

    # focus + cap columns (GOI always kept) — same policy as the aligned grid
    cov = {a["key"]: sum(1 for m in target_maps if a["key"] in m) for a in anchors}
    anchors = [a for a in anchors if a["is_goi"] or cov[a["key"]] > 0]
    if len(anchors) > max_cols:
        ordered = sorted(anchors, key=lambda a: (a["is_goi"], cov[a["key"]]), reverse=True)
        keep = {a["key"] for a in ordered[:max_cols]} | {goi_key}
        anchors = [a for a in anchors if a["key"] in keep]
    n_cols = len(anchors)
    row_tracks = [home_track] + targets
    n_rows = len(row_tracks)
    home_chrom = next((g.get("chrom", "") for g in home_track.get("genes", [])), "")

    # ---- per-row present points (col index + true midpoint + the gene) ----
    def _row_points(ri):
        pts = []
        if ri == 0:
            for i, a in enumerate(anchors):
                if not a["start"]:
                    continue
                pts.append({"i": i, "key": a["key"], "mid": (a["start"] + a["end"]) / 2.0,
                            "ident": 0.0, "conf": "", "strand": a["strand"], "n": 1,
                            "cov": None, "is_goi": a["is_goi"], "label": a["label"],
                            "chrom": home_chrom})
            return pts, home_chrom, 0
        tmap = target_maps[ri - 1]
        for i, a in enumerate(anchors):
            entry = tmap.get(a["key"])
            if not entry:
                continue
            g = entry["best"]
            gs, ge = _get_coords(g)
            pts.append({"i": i, "key": a["key"], "mid": (gs + ge) / 2.0,
                        "ident": g.get("identity", 0.0),
                        "conf": (g.get("confidence") or "").upper(),
                        "strand": g.get("strand", "+"), "n": entry["n"],
                        "cov": g.get("query_coverage"), "is_goi": a["is_goi"],
                        "label": a["label"], "chrom": g.get("chrom", "")})
        # On a fragmented assembly a row can straddle >1 scaffold; keep the
        # dominant (or GOI-bearing) one so the position line + label stay coherent.
        goi_ch = next((p["chrom"] for p in pts if p["is_goi"] and p["chrom"]), "")
        prim, _lo, _hi, n_other = _dominant_chrom_span(
            [(p["chrom"], p["mid"]) for p in pts], goi_ch)
        if prim:
            pts = [p for p in pts if p["chrom"] == prim] or pts
        return pts, prim, n_other
    row_pts, row_chrom, row_other = [], [], []
    for ri in range(n_rows):
        p, c, no = _row_points(ri)
        row_pts.append(p); row_chrom.append(c); row_other.append(no)

    # ---- geometry ----
    GRID_HEADER = "#f5f7fb"; GRID_FRAME = "#d8dee8"
    GRID_ROW_ALT = "#f8f9fc"; GRID_ROW_HOME = "#eef2f7"
    GRID_AXIS_TEXT = "#6b7280"; POS_LINE = "#c2cad6"; LEADER = "#c7cedb"
    MARGIN_TOP, HEADER_H = 60, 120
    POS_H, ARROW_LANE, ARROW_H = 26, 34, 18
    ROW_H = POS_H + ARROW_LANE
    COL_W, GOI_COL_W, RIGHT, LEGEND_H = 46, 66, 60, 120

    labels = [_svg_esc(re.sub(r"<[^>]+>", "", (t.get("label") or "")).strip())
              for t in row_tracks]
    TREE_W = 144 if (rooted_tree is not None) else 0
    LABEL_W = max(190, min(360, int(max((len(l) for l in labels), default=12) * 7.0 + 40)))
    LEFT = TREE_W + LABEL_W

    col_x, col_w, cx = [], [], LEFT
    for a in anchors:
        w = GOI_COL_W if a["is_goi"] else COL_W
        col_x.append(cx); col_w.append(w); cx += w
    grid_x1 = cx
    grid_y0 = MARGIN_TOP + HEADER_H
    grid_w = grid_x1 - LEFT
    grid_h = n_rows * ROW_H

    def col_center(i):
        return col_x[i] + col_w[i] / 2.0

    # per-row span label (chrom: lo-hi Mb) → adaptive right gutter
    def _span_text(ri):
        mids = [p["mid"] for p in row_pts[ri]]
        if not mids:
            return ""
        return _span_gutter_label(row_chrom[ri], min(mids), max(mids), row_other[ri])
    span_texts = [_span_text(ri) for ri in range(n_rows)]
    max_span = max((len(s) for s in span_texts), default=0)
    SPAN_W = (max(150, int(max_span * 6.0 + 26)) if max_span else 10)

    _caption = ("Columns = home gene order   ·   dots above each row = TRUE gene "
                "positions, leaders drop to the aligned column   ·   arrow points "
                "in coding strand")
    _legend_row_w = 560
    content_w = max(grid_w + SPAN_W, len(_caption) * 5.6, _legend_row_w)
    width = LEFT + content_w + RIGHT
    height = grid_y0 + grid_h + LEGEND_H

    P = []
    P.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')
    # header backplate behind column labels
    P.append(f'<rect x="{LEFT - 8:.1f}" y="{MARGIN_TOP - 12:.1f}" '
             f'width="{grid_w + SPAN_W + 16:.1f}" height="{HEADER_H + 12:.1f}" '
             f'fill="{GRID_HEADER}" stroke="{GRID_FRAME}" stroke-width="1" rx="8"/>')

    # GOI guide band (full height)
    for i, a in enumerate(anchors):
        if a["is_goi"]:
            P.append(f'<rect x="{col_x[i]:.1f}" y="{grid_y0 - 6:.1f}" '
                     f'width="{col_w[i]:.1f}" height="{grid_h + 12:.1f}" '
                     f'fill="{GOI_COLOUR}" opacity="0.07" stroke="{GOI_BORDER}" '
                     f'stroke-width="0.6" stroke-opacity="0.3"/>')
        elif i % 2 == 0:
            P.append(f'<rect x="{col_x[i]:.1f}" y="{grid_y0:.1f}" '
                     f'width="{col_w[i]:.1f}" height="{grid_h:.1f}" '
                     f'fill="#0f172a" opacity="0.015"/>')

    # column header labels (rotated)
    for i, a in enumerate(anchors):
        cc = col_center(i)
        ly = MARGIN_TOP + HEADER_H - 8
        fill = GOI_BORDER if a["is_goi"] else "#42495a"
        P.append(f'<text x="{cc:.1f}" y="{ly:.1f}" font-size="10.5" fill="{fill}" '
                 f'font-weight="{"700" if a["is_goi"] else "500"}" '
                 f'transform="rotate(-55 {cc:.1f} {ly:.1f})">{_svg_esc(a["label"])}</text>')

    def draw_arrow(i, base, strand, identity, conf, n_copies, is_home, yb, title, cov=None):
        x0 = col_x[i] + 4; x1 = col_x[i] + col_w[i] - 4
        conf = (conf or "").upper()
        if is_home:
            fill, stroke, dash, op = base, _darken_hex(base, 0.7), "", 0.95
        elif conf == "LOW":
            fill, stroke, dash, op = _lerp_hex(base, "#ffffff", 0.6), base, ' stroke-dasharray="2,2"', 0.85
        elif conf == "MEDIUM":
            fill, stroke, dash, op = _lerp_hex(base, "#ffffff", 0.42), base, ' stroke-dasharray="3.5,2"', 0.95
        else:
            fill, stroke, dash, op = _shade_by_identity(base, identity), _darken_hex(base, 0.65), "", 1.0
        d = _svg_arrow_path(x0, x1, yb, ARROW_H, strand)
        inner = [f'<title>{_svg_esc(title)}</title>',
                 f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="1"{dash} opacity="{op:.2f}"/>']
        if (not is_home) and identity >= 25:
            num = (f"{identity:.0f} ({cov*100:.0f})" if (cov is not None and cov < COVERAGE_FLAG_THRESHOLD)
                   else f"{identity:.0f}")
            avail = (x1 - x0) - 5
            fsz = max(6.5, min(10.0, avail / (len(num) * 0.60)))
            inner.append(f'<text x="{(x0+x1)/2:.1f}" y="{yb + ARROW_H/2 + 3.0:.1f}" '
                         f'text-anchor="middle" font-size="{fsz:.1f}" fill="#15181f" '
                         f'font-weight="{"700" if conf=="HIGH" else "600"}" pointer-events="none" '
                         f'style="paint-order:stroke;stroke:#ffffff;stroke-width:1.4">{num}</text>')
        if n_copies > 1:
            inner.append(f'<text x="{x1:.1f}" y="{yb - 2:.1f}" text-anchor="end" font-size="8.5" '
                         f'fill="{_darken_hex(base, 0.55)}" font-weight="700" pointer-events="none" '
                         f'style="paint-order:stroke;stroke:#ffffff;stroke-width:1.8">×{n_copies}</text>')
        P.append('<g class="acell">' + "".join(inner) + '</g>')

    # ---- rows ----
    for ri in range(n_rows):
        is_home = ri == 0
        track = row_tracks[ri]
        row_y = grid_y0 + ri * ROW_H
        pos_line_y = row_y + POS_H * 0.62
        arrow_top = row_y + POS_H
        arrow_yb = arrow_top + (ARROW_LANE - ARROW_H) / 2.0
        row_fill = GRID_ROW_HOME if is_home else (GRID_ROW_ALT if ri % 2 == 1 else "#ffffff")
        P.append(f'<rect x="{LEFT:.1f}" y="{row_y:.1f}" width="{grid_w + SPAN_W:.1f}" '
                 f'height="{ROW_H:.1f}" fill="{row_fill}"/>')

        pts = row_pts[ri]
        # --- position lane: line + dots at TRUE positions + slanted leaders ---
        if len(pts) >= 1:
            present_i = [p["i"] for p in pts]
            x_lo, x_hi = col_center(min(present_i)), col_center(max(present_i))
            mids = [p["mid"] for p in pts]
            lo, hi = min(mids), max(mids)
            def _posx(mid):
                if hi <= lo or x_hi <= x_lo:
                    return (x_lo + x_hi) / 2.0
                return x_lo + (mid - lo) / (hi - lo) * (x_hi - x_lo)
            if len(pts) >= 2:
                P.append(f'<line x1="{x_lo:.1f}" y1="{pos_line_y:.1f}" x2="{x_hi:.1f}" '
                         f'y2="{pos_line_y:.1f}" stroke="{POS_LINE}" stroke-width="1"/>')
            for p in pts:
                px = _posx(p["mid"]); cc = col_center(p["i"])
                lcol = GOI_COLOUR if p["is_goi"] else LEADER
                P.append(f'<path d="M{px:.1f},{pos_line_y + 2:.1f} L{cc:.1f},{arrow_top - 1:.1f}" '
                         f'stroke="{lcol}" stroke-width="{1.1 if p["is_goi"] else 0.6:.1f}" '
                         f'fill="none" opacity="{0.85 if p["is_goi"] else 0.55}"/>')
            for p in sorted(pts, key=lambda q: q["is_goi"]):
                px = _posx(p["mid"])
                base = GOI_COLOUR if p["is_goi"] else anchors[p["i"]]["colour"]
                if p["is_goi"]:
                    P.append(f'<circle cx="{px:.1f}" cy="{pos_line_y:.1f}" r="3.4" '
                             f'fill="{base}" stroke="{GOI_BORDER}" stroke-width="0.9"/>')
                else:
                    P.append(f'<circle cx="{px:.1f}" cy="{pos_line_y:.1f}" r="2.4" '
                             f'fill="{base}" stroke="{_darken_hex(base, 0.6)}" stroke-width="0.6"/>')

        # --- arrow lane ---
        for p in pts:
            a = anchors[p["i"]]
            if is_home:
                title = f'{a["label"]} (home reference)'
                draw_arrow(p["i"], a["colour"], a["strand"], 0.0, "", 1, True, arrow_yb, title)
            else:
                inv = p["strand"] and a["strand"] and p["strand"] != a["strand"]
                title = (f'{a["label"]} in {labels[ri]} — identity {p["ident"]:.1f}%'
                         + (", inverted" if inv else "")
                         + (f', {p["n"]}× copies' if p["n"] > 1 else ""))
                draw_arrow(p["i"], a["colour"], p["strand"], p["ident"], p["conf"],
                           p["n"], False, arrow_yb, title, cov=p["cov"])

        # row label
        ly = row_y + ROW_H / 2 + 4
        P.append(f'<text x="{TREE_W + 8:.1f}" y="{ly:.1f}" font-size="12" '
                 f'fill="#1a1d26" font-weight="{"700" if is_home else "600"}">'
                 f'{labels[ri] or "Home"}</text>')
        # span label
        if span_texts[ri]:
            P.append(f'<text x="{grid_x1 + 10:.1f}" y="{pos_line_y + 3:.1f}" font-size="9" '
                     f'fill="{GRID_AXIS_TEXT}">{_svg_esc(span_texts[ri])}</text>')

    # ---- tree panel ----
    if TREE_W:
        row_species = [None] + [
            _grid_match_species((t.get("genome_id") or ""), species_order) for t in targets]
        P.extend(_render_grid_tree_panel(
            rooted_tree, row_species, 12, TREE_W - 6,
            lambda ri: grid_y0 + ri * ROW_H + ROW_H / 2))
        P.append(f'<text x="{TREE_W / 2:.1f}" y="{grid_y0 - 9:.1f}" text-anchor="middle" '
                 f'font-size="9.5" fill="{GRID_AXIS_TEXT}" font-style="italic">{_tree_label}</text>')

    # ---- frame ----
    P.append(f'<rect x="{LEFT:.1f}" y="{grid_y0:.1f}" width="{grid_w:.1f}" height="{grid_h:.1f}" '
             f'fill="none" stroke="{GRID_FRAME}" stroke-width="1"/>')

    # ---- legend (two rows) ----
    ly0 = grid_y0 + grid_h + 34
    lx = LEFT
    ARROW_W, SYM_GAP = 26, 9
    P.append(f'<text x="{lx:.1f}" y="{ly0 - 16:.1f}" font-size="10.5" fill="#8c95a6">'
             f'{_svg_esc(_caption)}</text>')
    cur = lx
    P.append(f'<text x="{cur:.1f}" y="{ly0 + 11:.1f}" font-size="10" font-weight="700" '
             f'fill="#42495a">Confidence:</text>')
    cur += 72
    for name, dash, fill in [
        ("HIGH", "", _shade_by_identity(GENE_PALETTE[0], 95)),
        ("MEDIUM", ' stroke-dasharray="3.5,2"', _lerp_hex(GENE_PALETTE[0], "#ffffff", 0.42)),
        ("LOW", ' stroke-dasharray="2,2"', _lerp_hex(GENE_PALETTE[0], "#ffffff", 0.6)),
    ]:
        d = _svg_arrow_path(cur, cur + ARROW_W, ly0, 14, "+")
        P.append(f'<path d="{d}" fill="{fill}" stroke="{_darken_hex(GENE_PALETTE[0],0.65)}" '
                 f'stroke-width="1"{dash}/>')
        P.append(f'<text x="{cur + ARROW_W + SYM_GAP:.1f}" y="{ly0 + 11:.1f}" font-size="10" '
                 f'fill="#42495a">{name}</text>')
        cur += ARROW_W + SYM_GAP + len(name) * 7 + 24
    cur += 12
    d = _svg_arrow_path(cur, cur + ARROW_W, ly0, 14, "+")
    P.append(f'<path d="{d}" fill="{_shade_by_identity(GOI_COLOUR, 95)}" stroke="{GOI_BORDER}" '
             f'stroke-width="1"/>')
    P.append(f'<text x="{cur + ARROW_W + SYM_GAP:.1f}" y="{ly0 + 11:.1f}" font-size="10" '
             f'fill="#42495a">GOI</text>')
    cur += ARROW_W + SYM_GAP + 3 * 7 + 24
    P.append(f'<circle cx="{cur + 4:.1f}" cy="{ly0 + 5:.1f}" r="2.6" fill="{GENE_PALETTE[0]}"/>')
    P.append(f'<text x="{cur + 4 + 2.6 + SYM_GAP:.1f}" y="{ly0 + 9:.1f}" font-size="10" '
             f'fill="#42495a">true position (dot) → aligned column</text>')
    P.append(f'<text x="{lx:.1f}" y="{ly0 + 33:.1f}" font-size="10" fill="#8c95a6">'
             f'number = % identity &#183; (n) = % coverage when &lt; 80% &#183; '
             f'shade = % identity &#183; per-row dots normalised to that genome</text>')

    # ---- title ----
    title_txt = f"Anchor positions - GOI: {_svg_esc(goi_label)}"
    P.insert(1, f'<text class="grid-title" x="{LEFT:.1f}" y="30" font-size="17" '
                f'font-weight="700" fill="#1a1d26">{title_txt}</text>')
    P.insert(2, f'<text class="grid-subtitle" x="{LEFT:.1f}" y="48" font-size="11" '
                f'fill="{GRID_AXIS_TEXT}">{n_cols} anchor genes x {n_rows} genomes | '
                f'aligned columns + true positions per row</text>')

    return _assemble_grid_html("\n".join(P), width, height, "SynVoy anchor positions")


def _assemble_grid_html(svg_content, width, height, title):
    """Self-contained HTML for the anchor-grid figure (literal colours only,
    so the static-SVG export is correct in any viewer)."""
    css = """
:root {
    --bg-1: #f3f5fa;
    --bg-2: #eef1f6;
    --panel: #ffffff;
    --panel-border: #e1e6ef;
    --text: #1a1d26;
    --muted: #6b7280;
}
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
    background: radial-gradient(1200px 600px at 15% 10%, #f9fafc 0%, var(--bg-1) 45%, var(--bg-2) 100%);
    font-family: "IBM Plex Sans", "DejaVu Sans", "Liberation Sans", sans-serif;
    color: var(--text);
}
.plot-wrapper { width: 100%; overflow: auto; padding: 20px 24px 30px; }
.grid-svg {
    background: var(--panel);
    border-radius: 14px;
    border: 1px solid var(--panel-border);
    box-shadow: 0 12px 30px rgba(15, 23, 42, 0.12), 0 2px 8px rgba(15, 23, 42, 0.06);
    display: block;
}
.grid-svg text {
    font-family: "IBM Plex Sans", "DejaVu Sans", "Liberation Sans", sans-serif;
    letter-spacing: 0.1px;
}
.grid-title {
    font-family: "Source Serif 4", "Georgia", "Times New Roman", serif;
    letter-spacing: 0.2px;
}
.grid-subtitle {
    font-family: "IBM Plex Sans", "DejaVu Sans", "Liberation Sans", sans-serif;
    letter-spacing: 0.1px;
}
.acell { transition: opacity 0.15s ease; }
.acell:hover { opacity: 0.65; }
.acol-lbl { font-family: inherit; letter-spacing: 0.15px; }
@media (max-width: 1100px) {
    .plot-wrapper { padding: 12px; }
    .grid-svg { border-radius: 10px; }
}
"""
    return (
        '<!DOCTYPE html>\n<html lang="en"><head><meta charset="UTF-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">'
        f'<title>{_svg_esc(title)}</title><style>{css}</style></head><body>'
        '<div class="plot-wrapper">'
        f'<svg class="grid-svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">'
        f'{svg_content}</svg></div></body></html>'
    )


# ======================================================================
# Main
# ======================================================================

def main():
    ap = argparse.ArgumentParser(description="SynVoy synteny plot")
    ap.add_argument("--home_bed",       required=True)
    ap.add_argument("--home_gff",       default=None)
    ap.add_argument("--home_species",   default="",
                    help="Home-genome species name, used to label the home "
                         "track (e.g. 'Homo sapiens') instead of 'Home genome'.")
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
    ap.add_argument("--goi_min_px",     type=int, default=0,
                    help="Minimum on-screen width (px) for GOI/toxin genes so their exon "
                         "structure is visible at locus scale (0=off; ~30 shows notches). "
                         "Exons scale to the enlarged box; lengths become schematic for GOIs.")
    ap.add_argument("--caption_file",   default=None,
                    help="Path to a UTF-8 text file whose lines are rendered as a caption "
                         "block at the bottom of the figure (e.g. conserved-locus/no-toxin notes).")
    ap.add_argument("--max_goi_per_genome", type=int, default=10,
                    help="Max GOI/toxin genes drawn per genome (default 10). Raise for tandem-array "
                         "loci so large copy-number clusters (e.g. 31-copy GR1) aren't truncated.")
    ap.add_argument("--goi_zoom",       type=float, default=1.0,
                    help="Local magnification factor for the GOI/toxin window only: the genomic "
                         "scale (px/bp) is multiplied by this factor inside the band spanned by GOI "
                         "genes and left at 1x for flanking regions (piecewise-linear 'broken' axis). "
                         "Spreads tandem toxin arrays apart so they stop overlapping while keeping "
                         "flanking spacing compact (1.0=off; ~8-12 for dense arrays).")
    ap.add_argument("--orient_to_home", action="store_true",
                    help="Reverse-complement (in plot space) any target track whose flanking-anchor "
                         "order runs opposite to the home reference, so every track reads in the same "
                         "left-to-right orientation (no mirror-imaged 'flipped' rows). Reflects gene "
                         "positions, flips strands/arrows, and mirrors the exon model per gene.")
    ap.add_argument("--max_legend_entries", type=int, default=25, help="Maximum number of flanking genes to show in legend")
    ap.add_argument("--ribbon_alpha_dense", type=float, default=0.20, help="Alpha for flanking ribbons")
    ap.add_argument("--hide_goi_absent", action="store_true",
                    help="Hide target tracks with no GOI-like annotation when informative tracks exist")
    ap.add_argument("--pub_svg", action="store_true",
                    help="Also write a publication SVG: same as the interactive "
                         "HTML view, with every home-genome gene labelled.")
    ap.add_argument("--no_anchor_grid", dest="anchor_grid", action="store_false",
                    help="Skip the anchor-grid (aligned-column) figure. By "
                         "default an additional '*_anchor_grid.html'/'.svg' is "
                         "written alongside the ribbon plot and matrix.")
    ap.set_defaults(anchor_grid=True)
    # Legacy flags from the previous narrow-layout publication renderer.
    # Accepted but ignored so existing Nextflow invocations keep working.
    ap.add_argument("--pub_width", type=int, default=183, help=argparse.SUPPRESS)
    ap.add_argument("--pub_palette", default="okabe_ito",
                    choices=["okabe_ito", "tableau"], help=argparse.SUPPRESS)
    ap.add_argument("--output",         required=True)
    ap.add_argument("--common_names", choices=("both", "common", "scientific", "off"),
                    default="both",
                    help="Species label style. 'both' shows 'Scientific (common)'.")
    ap.add_argument("--common_names_tsv", default="",
                    help="Optional 2-column TSV (scientific<TAB>common) "
                         "overriding NCBI lookups.")
    ap.add_argument("--no_network", action="store_true",
                    help="Skip the NCBI 'datasets' CLI lookup for common names.")
    ap.add_argument("--clade_count", type=int, default=4,
                    help="Number of clades for tree-leaf colouring (default 4). "
                         "Iteratively splits the largest clade of the "
                         "midpoint-rooted tree (topology-driven). Each "
                         "clade is rendered with a distinct colour from "
                         "the colour-blind-safe CLADE_PALETTE.")
    args = ap.parse_args()

    global _MAX_GOI_PER_GENOME
    _MAX_GOI_PER_GENOME = max(1, args.max_goi_per_genome)

    # Initialize common-name resolver up front (cheap; just reads cache).
    if args.common_names != "off":
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import synvoy_taxa  # noqa: WPS433 — sibling helper module
            synvoy_taxa.init_lookup(
                common_names_tsv=args.common_names_tsv or None,
                allow_network=not args.no_network,
            )
            globals()["_synvoy_taxa"] = synvoy_taxa
        except Exception as exc:
            print(f"[plot] Common-name lookup unavailable: {exc}", file=sys.stderr)
            globals()["_synvoy_taxa"] = None
    globals()["_common_name_mode"] = args.common_names

    # -- 0. Load species mapping -----------------------------------------
    species_map = {}  # accession -> species name
    if args.species_map and os.path.exists(args.species_map):
        with open(args.species_map) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    species_map[parts[0]] = parts[1]
        print(f"[plot] Loaded species mapping for {len(species_map)} genomes")

    # -- 1. Parse inputs -------------------------------------------------

    home_genes = parse_bed(args.home_bed)
    if not home_genes:
        msg = f"ERROR: empty home BED: {args.home_bed}"
        print(msg, file=sys.stderr)
        error_html = f"""<!DOCTYPE html>
<html><head><title>SynVoy Error</title></head>
<body style="font-family:sans-serif;padding:40px;background:#f8f9fb;">
<h1 style="color:#dc2626;">SynVoy Synteny Plot Failed</h1>
<p style="color:#555;">{_html_escape(msg)}</p>
</body></html>"""
        with open(args.output, "w") as f:
            f.write(error_html)
        sys.exit(2)
    home_genes.sort(key=lambda g: g["start"])

    query_intervals = []
    if args.query_bed and os.path.exists(args.query_bed):
        for g in parse_bed(args.query_bed):
            query_intervals.append({"chrom": g["chrom"],
                                    "start": g["start"], "end": g["end"],
                                    "strand": g.get("strand", "+")})

    # The synteny-block BED that feeds `home_genes` carries flanking genes
    # only — by design, since `extract_flanking_genes.py` filters the GOI
    # itself out of the flanking set. That means the home track has every
    # neighbour but a *gap* where the GOI should be, which prevents the
    # legacy plot from drawing the GOI on the home row at all (severe bug).
    # Synthesize one from query_bed + home_gff before downstream track
    # construction sees `home_genes`.
    _synthesize_home_goi_gene(home_genes, query_intervals, args.home_gff)

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

    # -- NCBI-taxonomy species tree (correct species topology) --------------
    # The GOI --tree is a *gene* tree: its order scrambles the species (Apis
    # spp. split, ants among bees). Build a proper species tree once, from the
    # run's species list, and let every cladogram/row-order use it; fall back to
    # the gene tree only when taxonomy is unavailable.
    args._species_tree_newick = None
    try:
        import synvoy_tree as _stree
        sp_list = []
        _hs = (getattr(args, "home_species", "") or "").strip()
        if _hs and _hs.lower() not in ("home", "unknown", ""):
            sp_list.append(_hs)
        for gff in args.target_gffs:
            gid = clean_genome_name(os.path.basename(gff).replace(".gff", ""))
            name = next((spn for acc, spn in species_map.items() if acc in gid), None)
            if not name:
                name = re.sub(r"\.(fa|fna|fasta)$", "", gid).replace("_", " ")
            sp_list.append(name)
        args._species_tree_newick = _stree.build_taxonomy_tree_newick(
            sp_list, allow_network=not args.no_network)
        if args._species_tree_newick:
            order = [lf.name for lf in _stree.parse_newick_tree(args._species_tree_newick).leaves()]
            print(f"[plot] NCBI-taxonomy species tree built ({len(order)} species)")
            # Order ribbon-plot rows by the species tree too (consistent across plots).
            if order:
                tree_target_order = order
        else:
            print("[plot] taxonomy species tree unavailable — falling back to the GOI gene tree")
    except Exception as exc:
        print(f"[plot] species-tree build failed ({exc}); using GOI gene tree", file=sys.stderr)

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
            from collections import Counter
            block_counter = Counter()
            block_genes = defaultdict(list)
            for g in genes_all:
                gname = g.get("name", "")
                m = re.search(r'_b(\d+)_', gname)
                if m:
                    bid = m.group(1)
                    block_counter[bid] += 1
                    block_genes[bid].append(g)

            if block_counter:
                # Blocks that contain at least one GOI gene
                goi_bids = set()
                for g in genes_all:
                    if _is_goi_target_gene(g):
                        m = re.search(r'_b(\d+)_', g.get("name", ""))
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
            # are almost certainly low-complexity/spurious matches.
            # Without synteny context, orthology cannot be established.
            # Drop them unless they are the ONLY GOI hits in this genome.
            unsupported_goi_chroms = {
                ch for ch in goi_chroms
                if chrom_flank_counts.get(ch, 0) == 0
            }
            supported_goi_chroms = goi_chroms - unsupported_goi_chroms

            if unsupported_goi_chroms and supported_goi_chroms:
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
                display = f"{sp_name} ({genome_id})"
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

    home_chrom = home_genes[0]["chrom"]
    # Label the home track with the real species name when known, falling back
    # to the generic "Home genome" only when it's unset/unknown.
    _hs = (getattr(args, "home_species", "") or "").strip().replace("_", " ")
    if _hs and _hs.lower() not in ("home", "unknown", ""):
        home_label = f"{_format_species_label(_hs)} ({home_chrom})"
    else:
        home_label = f"Home genome ({home_chrom})"

    raw_tracks = [{
        "label":        home_label,
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

    # -- 4a. Orient every track to the home reference --------------------
    # Natively reverse-assembled scaffolds otherwise render mirror-imaged
    # (genes + arrows on the opposite side, ribbons crossing). Flip those so
    # all rows read left-to-right the same way. Anchor centres are preserved,
    # so the GOI guide line and ribbon endpoints stay aligned.
    if getattr(args, "orient_to_home", False):
        n_flipped = _orient_tracks_to_home(all_tracks)
        if n_flipped:
            print(f"[plot] Oriented {n_flipped} track(s) to home (reverse-complemented)")

    # -- 4b. Adaptive widening for sparse plots --------------------------
    # When most tracks have few genes spread across a wide bp range, the
    # default px-per-bp scale leaves gene models too narrow to read. Scale
    # each gene's visual width around its center; centers don't move so
    # ortholog ribbons stay aligned. No-op on dense plots.
    widen_factor = _widen_sparse_plot(all_tracks)
    if widen_factor > 1.0:
        print(f"[plot] Adaptive gene widening factor: {widen_factor:.2f}×")

    # -- 5. Sub-track assignment -----------------------------------------
    # Pure overlap-based bumping: a sub-track is spent only when genes
    # would collide on a single row.
    for track in all_tracks:
        _assign_sub_tracks(track["genes"], track["offset"])

    # -- 6. Build subtitle -----------------------------------------------

    subtitle_bits = [
        "Genes coloured by homology group",
        "★ = resolved GOI",
        "dashed = ambiguous",
        "exon blocks + intron lines",
        "ribbons connect orthologs",
        "// = compressed gaps",
    ]
    if hidden_absent_tracks:
        subtitle_bits.append(f"{hidden_absent_tracks} GOI-absent track(s) hidden")
    if ambiguous_track_count:
        subtitle_bits.append(f"{ambiguous_track_count} ambiguous track(s)")

    # -- 7. Render SVG ---------------------------------------------------

    html = render_synteny_html(
        all_tracks, gene_colours, goi_genome_colours,
        home_products, args,
        subtitle_bits, hidden_absent_tracks,
        ambiguous_track_count, resolved_track_count,
    )

    with open(args.output, "w") as f:
        f.write(html)
    print(f"Synteny plot (HTML) saved to {args.output}")

    # -- 7a-bis. Always export a static-SVG sibling of the interactive HTML.
    # This is a verbatim extraction of the inline <svg> with the page's CSS
    # CDATA-embedded — the result is visually identical to the HTML but
    # standalone. It's what the user actually wants to drop into a paper
    # or a README and not the narrow publication-format `--pub_svg` view.
    static_svg_path = args.output.replace(".html", "_view.svg")
    if static_svg_path == args.output:
        static_svg_path = args.output + ".view.svg"
    try:
        _export_html_inline_svg(args.output, static_svg_path)
        print(f"Static-view SVG saved to {static_svg_path}")
    except Exception as exc:
        print(f"  (could not export static-view SVG: {exc})", file=sys.stderr)

    # -- 7b. Render Publication SVG --------------------------------------
    # Same content as the interactive HTML but with every home-genome gene
    # labelled on the canvas, so the figure is self-describing in print.
    if args.pub_svg:
        pub_svg_content = render_publication_svg(
            all_tracks, gene_colours, goi_genome_colours,
            home_products, args,
            subtitle_bits, hidden_absent_tracks,
            ambiguous_track_count, resolved_track_count,
        )
        pub_output = args.output.replace(".html", ".svg")
        if pub_output == args.output:
            pub_output += ".svg"
        with open(pub_output, "w") as f:
            f.write(pub_svg_content)
        print(f"Publication SVG saved to {pub_output}")
    print(f"  Tracks: {n_tracks} ({n_tracks - 1} target genomes)")
    print(f"  GOI tracks: {resolved_track_count} resolved, {ambiguous_track_count} ambiguous")
    if hidden_absent_tracks:
        print(f"  Hidden absent tracks: {hidden_absent_tracks}")
    print(f"  Gap compression: active (>{args.gap_threshold} bp -> {args.gap_visual_size} bp visual)")

    # -- 7c. Anchor-grid view (matrix × synteny hybrid) ------------------
    # An additional figure: orthologues aligned into shared home-gene
    # columns, one row per species. Complements the ribbon plot (real
    # positions) and the matrix (presence/absence). Self-contained HTML +
    # a static SVG sibling.
    if getattr(args, "anchor_grid", True):
        try:
            grid_html = render_anchor_grid(
                all_tracks, gene_colours, goi_genome_colours,
                home_products, args,
            )
            grid_output = args.output.replace("_synteny_plot.html", "_anchor_grid.html")
            if grid_output == args.output:
                grid_output = args.output.replace(".html", "_anchor_grid.html")
            with open(grid_output, "w") as f:
                f.write(grid_html)
            print(f"Anchor-grid plot (HTML) saved to {grid_output}")
            grid_svg = grid_output.replace(".html", ".svg")
            try:
                _export_html_inline_svg(grid_output, grid_svg)
                print(f"Anchor-grid SVG saved to {grid_svg}")
            except Exception as exc:
                print(f"  (could not export anchor-grid SVG: {exc})", file=sys.stderr)
        except Exception as exc:
            print(f"  (anchor-grid render failed: {exc})", file=sys.stderr)

        # Companion real-position variant: same orthologs, true per-row
        # genomic positions (shows spacing / gaps / rearrangements the aligned
        # column grid hides). Written as '*_gene_positions.html'/'.svg'.
        try:
            pos_html = render_anchor_grid_positional(
                all_tracks, gene_colours, goi_genome_colours,
                home_products, args,
            )
            pos_output = args.output.replace("_synteny_plot.html", "_gene_positions.html")
            if pos_output == args.output:
                pos_output = args.output.replace(".html", "_gene_positions.html")
            with open(pos_output, "w") as f:
                f.write(pos_html)
            print(f"Gene-position map (HTML) saved to {pos_output}")
            try:
                _export_html_inline_svg(pos_output, pos_output.replace(".html", ".svg"))
                print(f"Gene-position map SVG saved to {pos_output.replace('.html', '.svg')}")
            except Exception as exc:
                print(f"  (could not export gene-position SVG: {exc})", file=sys.stderr)
        except Exception as exc:
            print(f"  (gene-position map render failed: {exc})", file=sys.stderr)

        # Threaded anchor grid: the aligned-column grid PLUS a light real-
        # position dot-line (with leaders down to each column) above every
        # genome's arrow row. Written as '*_anchor_positions.html'/'.svg'.
        try:
            thr_html = render_anchor_grid_threaded(
                all_tracks, gene_colours, goi_genome_colours,
                home_products, args,
            )
            thr_output = args.output.replace("_synteny_plot.html", "_anchor_positions.html")
            if thr_output == args.output:
                thr_output = args.output.replace(".html", "_anchor_positions.html")
            with open(thr_output, "w") as f:
                f.write(thr_html)
            print(f"Anchor-positions plot (HTML) saved to {thr_output}")
            try:
                _export_html_inline_svg(thr_output, thr_output.replace(".html", ".svg"))
                print(f"Anchor-positions SVG saved to {thr_output.replace('.html', '.svg')}")
            except Exception as exc:
                print(f"  (could not export anchor-positions SVG: {exc})", file=sys.stderr)
        except Exception as exc:
            print(f"  (anchor-positions render failed: {exc})", file=sys.stderr)

    # -- 8. Tree plot (separate HTML) ------------------------------------
    tree_output = args.output.replace("_synteny_plot.html", "_tree.html")
    if tree_output == args.output:
        tree_output = args.output.replace(".html", "_tree.html")
    _render_tree_svg(args.tree, goi_genome_colours, tree_output,
                     species_map=species_map, clade_count=args.clade_count)


if __name__ == "__main__":
    main()
