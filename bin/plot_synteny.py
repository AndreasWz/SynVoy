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
import colorsys
import json
import math
import os
import re
import sys
from collections import defaultdict
from html import escape as _html_escape
from urllib.parse import unquote

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
                    "text": f"◆ {g['chrom'][:16]}",
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
        for g in t["genes"]:
            c = (g["start_plot"] + g["end_plot"]) / 2.0
            half = (g["end_plot"] - g["start_plot"]) / 2.0 * factor
            g["start_plot"] = c - half
            g["end_plot"]   = c + half
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
    combined_style = "\n".join(style_blocks)
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
    TOP_MARGIN    = 72
    BOTTOM_MARGIN = 90  # More room for legend
    TRACK_PAD     = 10
    MIN_GENE_PX   = 4
    EXON_RX       = 3
    RIBBON_GAP    = min(8, TRACK_MARGIN * 0.1)

    # ---- Compute x range ----
    all_x_bp = []
    for track in all_tracks:
        x_off = track["offset"]
        for g in track["genes"]:
            all_x_bp.append(g["start_plot"] - x_off)
            all_x_bp.append(g["end_plot"] - x_off)

    if not all_x_bp:
        x_min_bp, x_max_bp = -1000, 1000
    else:
        x_min_bp, x_max_bp = min(all_x_bp), max(all_x_bp)

    pad_bp = (x_max_bp - x_min_bp) * 0.05 + 5000
    x_min_bp -= pad_bp
    x_max_bp += pad_bp

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

    available_w = plot_w - LEFT_MARGIN - RIGHT_MARGIN
    bp_range = max(1, x_max_bp - x_min_bp)
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

    # ---- Coordinate helpers (closures) ----
    def bp2px(bp_val):
        return LEFT_MARGIN + (bp_val - x_min_bp) * scale

    def gene_px(gene, track):
        x_off = track["offset"]
        x0 = bp2px(gene["start_plot"] - x_off)
        x1 = bp2px(gene["end_plot"] - x_off)
        if x1 - x0 < MIN_GENE_PX:
            mid = (x0 + x1) / 2
            x0, x1 = mid - MIN_GENE_PX / 2, mid + MIN_GENE_PX / 2
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

    # ---- Ribbons (drawn first, behind genes) ----
    svg_parts.append('<g class="ribbons">')
    for ti in range(len(all_tracks) - 1):
        upper = all_tracks[ti]
        lower = all_tracks[ti + 1]
        for lg in lower["genes"]:
            home_id = lg.get("home_gene_id", "")
            if not home_id:
                continue
            
            y_lg_top = gene_yb(ti + 1, lg) - RIBBON_GAP
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
                        ribbon_alpha = 0.40
                        if _is_goi_target_gene(lg) and not _is_resolved_goi_target_gene(lg):
                            ribbon_alpha = 0.15
                    else:
                        identity = lg.get("identity", 50.0)
                        ribbon_alpha = 0.08 + (min(identity, 100) / 100) * 0.35

                    y_ug_bot = gene_yb(ti, ug) + GENE_H + RIBBON_GAP

                    ux0, ux1 = gene_px(ug, upper)
                    lx0, lx1 = gene_px(lg, lower)
                    fill = _hex_to_rgba(colour, ribbon_alpha)
                    edge = _hex_to_rgba(colour, min(1.0, ribbon_alpha * 1.8))
                    path_d = _svg_ribbon_path(ux0, ux1, y_ug_bot, lx0, lx1, y_lg_top)
                    svg_parts.append(
                        f'<path d="{path_d}" fill="{fill}" stroke="{edge}" '
                        f'stroke-width="0.5" class="ribbon" '
                        f'data-homology="{_svg_esc(home_id)}" '
                        f'data-upper-track="{ti}" data-lower-track="{ti + 1}"/>'
                    )
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

        for gene in track["genes"]:
            name = gene["name"]
            home_id = gene.get("home_gene_id", name)
            goi_f = _is_goi_target_gene(gene) if not track["is_home"] else (is_goi(name) or is_goi(home_id))
            resolved_goi_f = _is_resolved_goi_target_gene(gene) if not track["is_home"] else goi_f

            home_force = force_home_labels and track["is_home"]
            if not (goi_f or home_force):
                continue

            label = _gene_display_label(gene, track, home_products, goi_f, resolved_goi_f)

            g_start, g_end = _get_coords(gene)
            xc_px = bp2px((g_start + g_end) / 2 - x_off)
            yb_px = gene_yb(ti, gene)
            lbl_y = yb_px - 4

            lbl_class = "gene-label"
            if goi_f:
                lbl_class += " goi"
                label = "★ " + label

            svg_parts.append(
                f'<text x="{xc_px:.1f}" y="{lbl_y:.1f}" '
                f'transform="rotate(-45 {xc_px:.1f} {lbl_y:.1f})" '
                f'class="{lbl_class} track-item" data-track-idx="{ti}" '
                f'font-size="{fsize}">{_svg_esc(label)}</text>'
            )
    svg_parts.append('</g>')

    # ---- Pinned-labels layer (populated by JS on click-to-pin) ----
    svg_parts.append('<g class="pinned-labels"></g>')

    # ---- Gap breaks & chromosome labels ----
    for ti, track in enumerate(all_tracks):
        yb = track_y[ti]
        x_off = track["offset"]
        th = track_heights[ti]

        for brk in track.get("breaks", []):
            brk_px = bp2px(brk["x"] - x_off)
            if brk.get("is_chrom_break"):
                svg_parts.append(
                    f'<line x1="{brk_px:.1f}" y1="{yb - TRACK_PAD:.1f}" '
                    f'x2="{brk_px:.1f}" y2="{yb + th + TRACK_PAD:.1f}" '
                    f'class="chrom-break-line track-item" data-track-idx="{ti}"/>'
                )
            else:
                svg_parts.append(
                    f'<text x="{brk_px:.1f}" y="{yb + th/2:.1f}" '
                    f'text-anchor="middle" dominant-baseline="central" '
                    f'class="break-label track-item" data-track-idx="{ti}">// {_svg_esc(brk["text"])}</text>'
                )

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
    scale_len_bp = args.scale_bar_len
    sb_x1_px = bp2px(x_max_bp - pad_bp * 0.5)
    sb_x0_px = sb_x1_px - scale_len_bp * scale
    sb_y = total_h - BOTTOM_MARGIN + 20
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

    # ---- Track labels (left margin) ----
    # Layout:
    #   x=  8..26 → collapse-toggle button (own column, never under text)
    #   x=  ...   → GOI clade dot, just left of the species label
    #   x=  LEFT_MARGIN-14 right edge → species, accession, status (right-aligned)
    TOGGLE_X      = 8
    TOGGLE_W      = 18
    TEXT_RIGHT_X  = LEFT_MARGIN - 14
    DOT_X         = LEFT_MARGIN - 4

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

        # GOI clade color indicator (right of all labels, before the plot)
        g_clr = goi_genome_colours.get(track.get("genome_id", ""), GOI_COLOUR)
        svg_parts.append(
            f'<circle cx="{DOT_X}" cy="{yb:.1f}" r="4" fill="{g_clr}" />'
        )

        # Species/common name (NCBI 'datasets' lookup if enabled).
        species_label = _format_species_label(species)
        svg_parts.append(
            f'<text x="{TEXT_RIGHT_X}" y="{yb:.1f}" '
            f'text-anchor="end" class="track-label" font-style="italic">{_svg_esc(species_label)}</text>'
        )

        # Accession and span
        if acc:
            span_str = acc
            if track["genes"]:
                chroms = sorted({g["chrom"] for g in track["genes"]})
                chr_str = chroms[0] if len(chroms) == 1 else f"{len(chroms)} chr"
                span_str = f"{acc} • {chr_str}"

            svg_parts.append(
                f'<text x="{TEXT_RIGHT_X}" y="{yb + 14:.1f}" '
                f'text-anchor="end" class="chrom-label">{_svg_esc(span_str)}</text>'
            )

        # GOI status indicators
        if not track["is_home"]:
            status = track.get("goi_status", "")
            if status == "absent":
                svg_parts.append(
                    f'<text x="{TEXT_RIGHT_X}" y="{yb + 28:.1f}" '
                    f'text-anchor="end" class="goi-status absent">✗ GOI not found</text>'
                )
            elif status == "ambiguous":
                svg_parts.append(
                    f'<text x="{TEXT_RIGHT_X}" y="{yb + 28:.1f}" '
                    f'text-anchor="end" class="goi-status ambiguous">⚠ Ambiguous orthology</text>'
                )
    svg_parts.append('</g>')

    # ---- Legend (Bottom left) ----
    # Removed as requested.

    # ---- Title & subtitle ----
    title_x = plot_w / 2
    svg_parts.append(
        f'<text x="{title_x:.1f}" y="28" text-anchor="middle" '
        f'class="plot-title">SynVoy Synteny Plot</text>'
    )
    if subtitle_bits:
        sub_text = " · ".join(subtitle_bits)
        svg_parts.append(
            f'<text x="{title_x:.1f}" y="48" text-anchor="middle" '
            f'class="plot-subtitle">{_svg_esc(sub_text)}</text>'
        )

    # ---- Assemble full HTML ----
    svg_content = "\n".join(svg_parts)
    html = _assemble_full_html(svg_content, plot_w, total_h)
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

        # Clean label
        label = leaf.name
        if "|" in label:
            parts = label.split("|")
            goi_part = parts[0]
            genome_part = parts[1] if len(parts) > 1 else ""
            genome_pretty = genome_part.replace("_fna_exon_ann", "").replace("_fna", "")
            if species_map:
                for acc, sp_name in species_map.items():
                    if acc in genome_pretty:
                        genome_pretty = f"{sp_name} ({genome_pretty})"
                        break
            label = f"{goi_part} | {genome_pretty}"
        else:
            label = f"{label} (home)"

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
    margin-bottom: 8px;
    color: var(--text-primary);
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
    const ribbonsGroup = svg.querySelector('.ribbons');
    const trackCount = new Set(Array.from(svg.querySelectorAll('[data-track-idx]')).map(el => el.dataset.trackIdx).filter(Boolean)).size;
    const collapsedTracks = new Set();

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

    function setTrackCollapsed(idx, collapsed) {
        const state = collapsed ? 'add' : 'delete';
        collapsedTracks[state](idx);
        svg.querySelectorAll('[data-track-idx="' + idx + '"]').forEach(el => {
            if (el.classList.contains('track-summary')) return;
            el.classList.toggle('track-item-hidden', collapsed);
        });
        svg.querySelectorAll('.track-summary[data-track-idx="' + idx + '"]').forEach(el => {
            el.classList.toggle('track-item-hidden', !collapsed);
        });
        svg.querySelectorAll('.track-toggle[data-track-idx="' + idx + '"]').forEach(el => {
            // The toggle is a <g> containing a <rect> background and a <text>
            // glyph. Update only the text child so the rect stays intact.
            const tEl = el.querySelector('text') || el;
            tEl.textContent = collapsed ? '▶' : '▼';
        });
        rebuildRibbons();
    }

    function buildTrackManager() {
        if (!trackManager) return;
        trackManager.innerHTML = '<h3>Manage tracks</h3>';
        uniqueTrackIndices().forEach(idx => {
            const row = document.createElement('label');
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.checked = !collapsedTracks.has(idx);
            if (idx === '0') {
                cb.checked = true;
                cb.disabled = true;
            }
            cb.addEventListener('change', () => setTrackCollapsed(idx, !cb.checked));
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

    function rebuildRibbons() {
        if (!ribbonsGroup) return;
        ribbonsGroup.innerHTML = '';
        const visible = uniqueTrackIndices().filter(idx => !collapsedTracks.has(idx));
        for (let vi = 0; vi < visible.length - 1; vi++) {
            const upperIdx = visible[vi];
            const lowerIdx = visible[vi + 1];
            const upperGenes = Array.from(svg.querySelectorAll('.gene-group[data-track-idx="' + upperIdx + '"]')).filter(el => !el.classList.contains('track-item-hidden'));
            const lowerGenes = Array.from(svg.querySelectorAll('.gene-group[data-track-idx="' + lowerIdx + '"]')).filter(el => !el.classList.contains('track-item-hidden'));
            const upperByHomology = new Map();
            upperGenes.forEach(g => {
                const hom = g.dataset.homology;
                if (!upperByHomology.has(hom)) upperByHomology.set(hom, []);
                upperByHomology.get(hom).push(g);
            });
            lowerGenes.forEach(lg => {
                const hom = lg.dataset.homology;
                const matches = upperByHomology.get(hom) || [];
                if (!matches.length) return;
                const lx0 = parseFloat(lg.dataset.x0);
                const lx1 = parseFloat(lg.dataset.x1);
                const ly = parseFloat(lg.dataset.yb) - 11;
                const identity = Math.max(0, Math.min(100, parseFloat(lg.dataset.identity || '50')));
                const isGoi = lg.dataset.isGoi === 'true';
                const clr = lg.dataset.fill || '#909090';
                const alpha = isGoi ? 0.40 : (0.08 + (identity / 100) * 0.35);
                matches.forEach(ug => {
                    const ux0 = parseFloat(ug.dataset.x0);
                    const ux1 = parseFloat(ug.dataset.x1);
                    const uy = parseFloat(ug.dataset.yb) + 11;
                    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                    path.setAttribute('class', 'ribbon');
                    path.setAttribute('data-homology', hom);
                    path.setAttribute('data-upper-track', String(upperIdx));
                    path.setAttribute('data-lower-track', String(lowerIdx));
                    path.setAttribute('d', buildRibbonPath(ux0, ux1, uy, lx0, lx1, ly));
                    path.setAttribute('fill', rgbaFromHex(clr, alpha));
                    path.setAttribute('stroke', rgbaFromHex(clr, Math.min(1, alpha * 1.6)));
                    path.setAttribute('stroke-width', isGoi ? '0.5' : '0.3');
                    ribbonsGroup.appendChild(path);
                });
            });
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
    rebuildRibbons();

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


def _assemble_full_html(svg_content, width, height):
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
  <button id="track-manager-btn" title="Manage tracks">📋</button>
</div>
<div id="track-manager" class="track-manager"></div>
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
    ap.add_argument("--pub_svg", action="store_true",
                    help="Also write a publication SVG: same as the interactive "
                         "HTML view, with every home-genome gene labelled.")
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

    # -- 8. Tree plot (separate HTML) ------------------------------------
    tree_output = args.output.replace("_synteny_plot.html", "_tree.html")
    if tree_output == args.output:
        tree_output = args.output.replace(".html", "_tree.html")
    _render_tree_svg(args.tree, goi_genome_colours, tree_output,
                     species_map=species_map, clade_count=args.clade_count)


if __name__ == "__main__":
    main()
