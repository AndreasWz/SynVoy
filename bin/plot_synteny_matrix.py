#!/usr/bin/env python3
"""
plot_synteny_matrix.py — phylogeny-anchored matrix view of SynVoy synteny.

Layout
──────
  ┌─tree──┬─species name──┬─col_1─┬─col_2─┬─...─┬─GOI─┬─...─┐  ← home (top)
  │       │               │       │       │     │     │     │
  │       │ target 1      │       │       │     │     │     │
  │       │ target 2      │       │       │     │     │     │
  │       │ ...           │       │       │     │     │     │

  - Columns are the ordered home flanking genes, with the GOI inserted at
    its correct genomic position (synthesized from --query_bed + --home_gff
    if not present in --home_bed). This fixes the long-standing bug where
    the GOI was missing from the home-genome row in the legacy plot.
  - Each target gene is mapped to a column via its `SynVoy_Parent` GFF attr
    or via the homology TSV (target_id → home_id).
  - Cell glyphs: directional arrow tinted by % identity; absent → muted gray.

Inputs (subset of plot_synteny.py — uses what is available, ignores the rest)
──────────────────────────────────────────────────────────────────────────────
  --home_bed        Synteny-block BED for the home genome (flanking genes).
  --home_gff        NCBI/RefSeq GFF for the home genome (gene symbols + GOI).
  --query_bed       BED with the query-gene location on the home genome.
  --target_gffs     SynVoy iterative-search GFFs (one per target genome).
  --target_names    Optional display names; otherwise derived from GFF stems.
  --homology_tsvs   Optional per-target homology TSVs (fallback for SynVoy_Parent).
  --tree            Newick tree of GOI leaves (used to order rows).
  --output          Output path (.html for interactive, .svg for static).

Design goals
────────────
  * Self-contained output (no external JS).
  * Linear-time rendering, ~600-line script.
  * Clear PDF/SVG export from the HTML by saving the embedded <svg> directly.
"""

import argparse
import os
import re
import sys
from collections import defaultdict, OrderedDict
from html import escape as _esc

# Local helper modules (siblings in bin/)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import synvoy_taxa  # noqa: E402
from synvoy_tree import (  # noqa: E402
    TreeNode,
    parse_newick_tree,
    parse_newick_leaf_order,
    midpoint_root,
    collapse_to_one_leaf_per_species,
    partition_clades,
    build_taxonomy_tree_newick,
    CLADE_PALETTE,
    species_from_leaf,
)

# ─────────────────────────────── parsers ────────────────────────────────────

def parse_bed(path):
    """BED → list of {chrom,start,end,name,strand}."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                start = int(p[1])
                end = int(p[2])
            except ValueError:
                continue
            name = p[3] if len(p) > 3 else f"{p[0]}:{start}-{end}"
            strand = p[5] if len(p) > 5 and p[5] in {"+", "-"} else "+"
            rows.append({"chrom": p[0], "start": start, "end": end,
                         "name": name, "strand": strand})
    return rows


def parse_gff_attrs(field):
    out = {}
    for kv in field.split(";"):
        kv = kv.strip()
        if not kv or "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def parse_home_gff_genes(path):
    """home GFF → dict[gene_id] = {symbol, product, chrom, start, end, strand}.

    Indexed by both `ID` ("gene-Melt") and `Name` ("Melt") so the BED-name
    lookups work regardless of which form the synteny-block BED uses.
    """
    by_id = {}
    if not path or not os.path.exists(path) or path == "NO_GFF":
        return by_id
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attrs = parse_gff_attrs(p[8])
            symbol = attrs.get("Name") or attrs.get("gene") or attrs.get("ID", "")
            product = attrs.get("description") or attrs.get("product") or ""
            try:
                start = int(p[3]) - 1
                end = int(p[4])
            except ValueError:
                continue
            entry = {"symbol": symbol, "product": product,
                     "chrom": p[0], "start": start, "end": end, "strand": p[6]}
            for key in (attrs.get("ID"), attrs.get("Name"), symbol):
                if key:
                    by_id[key] = entry
                    by_id[f"gene-{key}"] = entry
    return by_id


def parse_target_gff(path, candidate_regions=None):
    """SynVoy iterative GFF → list of cells indexed by SynVoy_Parent.

    `candidate_regions` is a list of (chrom, start, end) intervals from
    `regions.bed`; if provided, GFF rows outside these intervals are dropped.
    This is critical: without filtering, the matrix shows GOI hits from
    anywhere in the target genome (e.g. weak melittin matches on a different
    chromosome), which inflates the apparent number of orthologs.

    Returns: list of dicts with keys
        parent (home gene ID this cell maps to),
        role (goi|flanking),
        confidence, identity, strand, start, end, name, chrom,
        evidence_type, model_status, query_coverage, exons.
    """
    cells = []
    if not path or not os.path.exists(path):
        return cells
    accepted = {"gene", "mRNA", "transcript", "tandem_copy"}
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in accepted:
                continue
            attrs = parse_gff_attrs(p[8])
            parent = attrs.get("SynVoy_Parent", "")
            if not parent:
                continue
            try:
                ident = float(attrs.get("Identity", "0") or 0)
            except ValueError:
                ident = 0.0
            try:
                start = int(p[3]) - 1
                end = int(p[4])
            except ValueError:
                continue
            if candidate_regions and not _row_in_regions(p[0], start, end,
                                                          candidate_regions):
                continue
            cells.append({
                "parent": parent,
                "role": (attrs.get("SynVoyRole") or "").lower(),
                "confidence": (attrs.get("Confidence") or "").upper(),
                "identity": ident,
                "strand": p[6] if p[6] in {"+", "-"} else "+",
                "start": start, "end": end,
                "name": attrs.get("Name", ""),
                "chrom": p[0],
                "evidence_type": attrs.get("EvidenceType", ""),
                "model_status": attrs.get("ModelStatus", ""),
                "exons": attrs.get("Exons", ""),
            })
    return cells


def _row_in_regions(chrom, start, end, regions):
    """True if the [start, end) interval on `chrom` overlaps any region."""
    for (rchrom, rstart, rend) in regions:
        if rchrom != chrom:
            continue
        if start < rend and end > rstart:
            return True
    return False


def parse_candidate_regions(path):
    """regions.bed → list of (chrom, start, end). Empty if path missing."""
    regions = []
    if not path or not os.path.exists(path):
        return regions
    for row in parse_bed(path):
        regions.append((row["chrom"], row["start"], row["end"]))
    return regions


def parse_homology_tsv(path):
    """homology TSV → list of {target_id, home_id, role, confidence, identity}."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                p += [""] * (len(header) - len(p))
            try:
                ident = float(p[idx["identity"]] or 0)
            except (ValueError, KeyError):
                ident = 0.0
            rows.append({
                "target_id": p[idx.get("target_id", 0)],
                "home_id":   p[idx.get("home_id", 1)],
                "role":      (p[idx.get("role", 2)] or "").lower(),
                "confidence": (p[idx.get("confidence", 3)] or "").upper(),
                "identity":  ident,
            })
    return rows


# Tree helpers moved to bin/synvoy_tree.py — see imports at top of file.

def species_from_leaf(leaf):
    """`GOI_Melt|Apis_florea_fna_b0_l1_exon_ann` → `Apis_florea`."""
    if "|" not in leaf:
        return None
    tail = leaf.split("|", 1)[1]
    # Drop common sentinel suffixes piece by piece.
    suffixes = ("_exon_ann", "_fallback", "_flank_ann", "_full")
    for suf in suffixes:
        if tail.endswith(suf):
            tail = tail[: -len(suf)]
    # Drop trailing  "_b<N>_l<M>" or "_b<N>_fl<M>" segments.
    tail = re.sub(r"_(b|fl|l)\d+(_(b|fl|l)\d+)*$", "", tail)
    tail = re.sub(r"_(fa|fna|fasta)$", "", tail)
    return tail


# ─────────────────────────────── slot building ──────────────────────────────

def build_home_slots(home_bed_rows, query_bed_rows, home_gff_index):
    """Order home flanking genes by start; insert a GOI slot from query_bed.

    Returns (slots, goi_slot_id) where each slot is a dict:
        {id, label, start, end, strand, is_goi, product}
    """
    if not home_bed_rows:
        return [], None
    chrom = home_bed_rows[0]["chrom"]
    # Strip "gene-" prefix from labels; carry the raw name as id for joins.
    slots = []
    for row in sorted(home_bed_rows, key=lambda r: r["start"]):
        label, product = _label_and_product(row["name"], home_gff_index)
        slots.append({
            "id": row["name"], "label": label, "product": product,
            "start": row["start"], "end": row["end"],
            "strand": row["strand"], "is_goi": False,
        })

    # Synthesize GOI slot: union of query_bed intervals, look up gene symbol
    # in home_gff that overlaps this range. If query_bed is empty, scan
    # home_gff for a gene whose symbol appears in any home_bed name (no-op).
    goi_slot = _make_goi_slot(query_bed_rows, home_gff_index, chrom)
    if goi_slot is None:
        return slots, None

    # Drop flanking slots that are container genes for the GOI (≥50% of the
    # GOI is inside them and they are ≥3× larger). Without this, a 17kb locus
    # encompassing the GOI shows up as both a flanking column and the GOI
    # column, which confuses the matrix.
    goi_len = max(1, goi_slot["end"] - goi_slot["start"])
    cleaned = []
    for s in slots:
        ov = max(0, min(s["end"], goi_slot["end"]) - max(s["start"], goi_slot["start"]))
        s_len = max(1, s["end"] - s["start"])
        is_container = (ov / goi_len >= 0.5) and (s_len >= 3 * goi_len)
        if not is_container:
            cleaned.append(s)
    slots = cleaned

    # Insert GOI slot at correct genomic position.
    midpoint = (goi_slot["start"] + goi_slot["end"]) // 2
    insert_at = len(slots)
    for i, s in enumerate(slots):
        if midpoint < s["start"]:
            insert_at = i
            break
    slots.insert(insert_at, goi_slot)
    return slots, goi_slot["id"]


def _label_and_product(raw_name, home_gff_index):
    label = raw_name
    if label.startswith("gene-"):
        label = label[len("gene-"):]
    product = ""
    entry = home_gff_index.get(raw_name) or home_gff_index.get(label)
    if entry:
        if entry.get("symbol"):
            label = entry["symbol"]
        product = entry.get("product", "") or ""
    return label, product


def _make_goi_slot(query_bed_rows, home_gff_index, chrom):
    if not query_bed_rows:
        return None
    qrows = [r for r in query_bed_rows if r["chrom"] == chrom]
    if not qrows:
        qrows = query_bed_rows
    qstart = min(r["start"] for r in qrows)
    qend = max(r["end"] for r in qrows)
    qstrand = qrows[0]["strand"]
    # Find the gene in home_gff that best overlaps the query span. Score by
    # overlap fraction (overlap / gene_length) so a tightly-fitting gene like
    # Melt beats a 17kb container locus that happens to span the same range.
    best = None
    best_score = 0.0
    seen_ids = set()
    for entry in home_gff_index.values():
        if id(entry) in seen_ids:  # the index points the same dict at multiple keys
            continue
        seen_ids.add(id(entry))
        if entry.get("chrom") != chrom:
            continue
        ov = max(0, min(entry["end"], qend) - max(entry["start"], qstart))
        if ov <= 0:
            continue
        gene_len = max(1, entry["end"] - entry["start"])
        score = ov / gene_len  # prefer smaller, more specific genes
        if score > best_score:
            best = entry
            best_score = score
    if best:
        label = best.get("symbol") or "GOI"
        product = best.get("product", "") or ""
        return {
            "id": f"GOI_{label}",
            "label": label, "product": product or "gene of interest",
            "start": best["start"], "end": best["end"],
            "strand": best["strand"], "is_goi": True,
        }
    # Fall back to query_bed coords with a generic "GOI" label.
    return {
        "id": "GOI", "label": "GOI", "product": "gene of interest",
        "start": qstart, "end": qend, "strand": qstrand, "is_goi": True,
    }


def map_target_cells_to_slots(target_cells, slots, homology_rows):
    """For each slot, gather target cells whose `parent` (or homology home_id)
    matches the slot id. GOI slot collects every cell with role='goi'.
    """
    by_slot = defaultdict(list)
    slot_ids = {s["id"] for s in slots}
    goi_slot_id = next((s["id"] for s in slots if s["is_goi"]), None)

    # Index homology by target_id for fallback parent lookup.
    homology_by_tid = {h["target_id"]: h for h in homology_rows}

    for cell in target_cells:
        parent = cell["parent"]
        if cell["role"] == "goi":
            if goi_slot_id:
                by_slot[goi_slot_id].append(cell)
            continue
        if parent in slot_ids:
            by_slot[parent].append(cell)
            continue
        # Homology fallback: target_id in TSV may map to a home_id we know.
        cand = None
        for tid, h in homology_by_tid.items():
            if cell["name"] and cell["name"] in tid:
                cand = h["home_id"]
                break
        if cand and cand in slot_ids:
            by_slot[cand].append(cell)

    # For each slot keep the best cell (highest identity) but remember count.
    summary = {}
    for sid, cells in by_slot.items():
        cells.sort(key=lambda c: (c["confidence"] != "HIGH",
                                  c["confidence"] != "MEDIUM",
                                  -c["identity"]))
        summary[sid] = {"best": cells[0], "all": cells, "n": len(cells)}
    return summary


# ─────────────────────────────── rendering ──────────────────────────────────

# Color palette tuned for accessibility (Okabe-Ito-ish, plus muted neutrals).
COL_GOI         = "#D55E00"      # vermillion
COL_GOI_LIGHT   = "#FFB199"
COL_FLANKING    = "#0072B2"      # blue
COL_FLANKING_LO = "#A6CEE3"
COL_MISSING     = "#EFEFEF"
COL_MISSING_BG  = "#FAFAFA"
COL_GRID        = "#D9D9D9"
COL_TEXT        = "#222"
COL_TREE        = "#777"
COL_HEADER_BG   = "#F4F4F4"
COL_GOI_HEADER  = "#FFE7DD"

CONF_OPACITY = {"HIGH": 1.0, "MEDIUM": 0.75, "LOW": 0.50, "": 0.85}


def render_svg(slots, species_rows, home_label, goi_slot_id, rooted_tree=None):
    """Build an SVG string for the matrix.

    species_rows: ordered list of
        (display_label, summary_dict, is_home, tree_species_key, clade_color)
    rooted_tree:  midpoint-rooted TreeNode (or None for placeholder ladder).
    """
    cell_w, cell_h = 38, 28
    name_w = 240
    tree_w = 170
    header_h = 110
    gap_above = 8
    margin = 24

    n_rows = len(species_rows)
    n_cols = len(slots)
    matrix_w = n_cols * cell_w
    matrix_h = n_rows * cell_h

    total_w = margin + tree_w + name_w + matrix_w + margin
    total_h = margin + header_h + gap_above + matrix_h + 90  # 90 = two-row legend

    out = []
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_w} {total_h}" '
        f'width="{total_w}" height="{total_h}" '
        f'style="font-family: Helvetica, Arial, sans-serif; font-size: 11px; color: {COL_TEXT};">'
    )
    # Patterns for LOW-confidence cells (diagonal stripes signal "uncertain").
    out.append(
        '<defs>'
        '<pattern id="lowConfFlank" patternUnits="userSpaceOnUse" '
        f'width="6" height="6" patternTransform="rotate(45)">'
        f'<rect width="6" height="6" fill="#FFFFFF"/>'
        f'<rect width="2" height="6" fill="{COL_FLANKING}" opacity="0.40"/></pattern>'
        '<pattern id="lowConfGoi" patternUnits="userSpaceOnUse" '
        f'width="6" height="6" patternTransform="rotate(45)">'
        f'<rect width="6" height="6" fill="#FFFFFF"/>'
        f'<rect width="2" height="6" fill="{COL_GOI}" opacity="0.45"/></pattern>'
        '</defs>'
    )
    out.append(f'<rect width="100%" height="100%" fill="white"/>')

    # Define arrow marker via per-cell <polygon> (avoids marker reuse complications).
    matrix_x0 = margin + tree_w + name_w
    matrix_y0 = margin + header_h + gap_above

    # ── Column headers ────────────────────────────────────────────────────
    out.append(f'<g class="headers">')
    for ci, slot in enumerate(slots):
        x = matrix_x0 + ci * cell_w + cell_w / 2
        # Header background highlight for GOI column
        if slot["is_goi"]:
            out.append(
                f'<rect x="{matrix_x0 + ci * cell_w}" y="{margin}" '
                f'width="{cell_w}" height="{header_h + matrix_h + gap_above}" '
                f'fill="{COL_GOI_HEADER}" opacity="0.55"/>'
            )
        # Rotated label
        text_y = margin + header_h - 6
        rotation = -55
        weight = "700" if slot["is_goi"] else "500"
        fill = COL_GOI if slot["is_goi"] else COL_TEXT
        label = _esc(slot["label"])
        tooltip = _esc(slot.get("product", ""))
        out.append(
            f'<g transform="translate({x:.1f},{text_y}) rotate({rotation})">'
            f'<text text-anchor="start" font-weight="{weight}" fill="{fill}">'
            f'<title>{label} — {tooltip} ({slot["start"]:,}–{slot["end"]:,})</title>'
            f'{label}</text></g>'
        )
    out.append('</g>')

    # ── Matrix grid ──────────────────────────────────────────────────────
    out.append(f'<g class="grid">')
    out.append(
        f'<rect x="{matrix_x0}" y="{matrix_y0}" width="{matrix_w}" height="{matrix_h}" '
        f'fill="{COL_MISSING_BG}" stroke="{COL_GRID}" stroke-width="1"/>'
    )
    # Vertical grid lines
    for ci in range(1, n_cols):
        gx = matrix_x0 + ci * cell_w
        out.append(
            f'<line x1="{gx}" y1="{matrix_y0}" x2="{gx}" y2="{matrix_y0 + matrix_h}" '
            f'stroke="{COL_GRID}" stroke-width="0.5"/>'
        )
    # Horizontal grid lines
    for ri in range(1, n_rows):
        gy = matrix_y0 + ri * cell_h
        out.append(
            f'<line x1="{matrix_x0}" y1="{gy}" x2="{matrix_x0 + matrix_w}" y2="{gy}" '
            f'stroke="{COL_GRID}" stroke-width="0.5"/>'
        )
    out.append('</g>')

    # ── Species names + cells (rows) ─────────────────────────────────────
    out.append(f'<g class="rows">')
    for ri, row in enumerate(species_rows):
        display_name, summary, is_home, _sp_key, clade_color = row
        cy = matrix_y0 + ri * cell_h + cell_h / 2

        if is_home:
            out.append(
                f'<rect x="{margin + tree_w}" y="{matrix_y0 + ri * cell_h}" '
                f'width="{name_w + matrix_w}" height="{cell_h}" '
                f'fill="#FFFBF2" opacity="1"/>'
            )

        # Clade color swatch on the inside edge of the name column
        if clade_color:
            out.append(
                f'<rect x="{margin + tree_w + 2}" y="{matrix_y0 + ri * cell_h + 4}" '
                f'width="6" height="{cell_h - 8}" fill="{clade_color}" opacity="0.9"/>'
            )

        weight = "700" if is_home else "500"
        style = "italic" if not is_home else "normal"
        suffix = "  (home)" if is_home else ""
        # Species labels always render in the default text colour. The
        # clade is communicated by the coloured stripe to the right of the
        # name (above); colouring the species name itself made the label
        # hard to read on light clade hues and gave a noisy look.
        out.append(
            f'<text x="{margin + tree_w + name_w - 8}" y="{cy + 4}" '
            f'text-anchor="end" font-style="{style}" font-weight="{weight}" '
            f'fill="{COL_TEXT}">'
            f'{_esc(display_name + suffix)}</text>'
        )

        # Cells for each slot
        for ci, slot in enumerate(slots):
            x = matrix_x0 + ci * cell_w
            cell = (summary or {}).get(slot["id"])
            if cell:
                _draw_cell(out, x, matrix_y0 + ri * cell_h, cell_w, cell_h,
                           slot, cell, is_home,
                           n_copies=cell.get("_n_copies", 1))
            elif is_home:
                _draw_home_cell(out, x, matrix_y0 + ri * cell_h, cell_w, cell_h, slot)
    out.append('</g>')

    # ── Tree (real midpoint-rooted cladogram) ────────────────────────────
    out.append(_render_cladogram(rooted_tree, species_rows,
                                 margin, matrix_y0, tree_w, cell_h))

    # ── Legend ──────────────────────────────────────────────────────────
    out.append(_render_legend(margin, matrix_y0 + matrix_h + 22, COL_GOI, COL_FLANKING))

    out.append('</svg>')
    return "\n".join(out)


def _identity_color(base, identity):
    """Lerp from base color to its lightened variant by (1 - identity/100)."""
    ident = max(0.0, min(100.0, float(identity or 0)))
    t = (100.0 - ident) / 100.0
    return _lerp_hex(base, "#FFFFFF", min(0.55, t * 0.7))


def _lerp_hex(a, b, t):
    a = a.lstrip("#"); b = b.lstrip("#")
    ar, ag, ab = int(a[0:2], 16), int(a[2:4], 16), int(a[4:6], 16)
    br, bg, bb = int(b[0:2], 16), int(b[2:4], 16), int(b[4:6], 16)
    r = round(ar + (br - ar) * t)
    g = round(ag + (bg - ag) * t)
    bv = round(ab + (bb - ab) * t)
    return f"#{r:02X}{g:02X}{bv:02X}"


def _arrow_path(x, y, w, h, strand):
    """Arrow polygon points string; pointed in strand direction."""
    pad_x = 4
    pad_y = 5
    head = 7  # arrow-head depth
    left = x + pad_x
    right = x + w - pad_x
    top = y + pad_y
    bottom = y + h - pad_y
    if strand == "-":
        return (f"{right},{top} {left + head},{top} {left},{(top+bottom)/2:.1f} "
                f"{left + head},{bottom} {right},{bottom}")
    return (f"{left},{top} {right - head},{top} {right},{(top+bottom)/2:.1f} "
            f"{right - head},{bottom} {left},{bottom}")


def _draw_cell(out, x, y, w, h, slot, cell, is_home, n_copies=1):
    """Render one matrix cell. Confidence drives the visual distinction:
       HIGH   → solid fill, dark border, full opacity
       MEDIUM → solid fill, normal border, 80% opacity, slight desaturation
       LOW    → diagonal-stripe pattern fill, dashed border, lighter text
    """
    base = COL_GOI if slot["is_goi"] else COL_FLANKING
    conf = (cell.get("confidence") or "").upper()
    pts = _arrow_path(x, y, w, h, cell["strand"])
    ident = cell["identity"]

    # Three visually distinct tiers so the reader can tell HIGH from MEDIUM
    # from LOW at a glance — without reading numbers or hovering. This is
    # critical for the melittin case where the pipeline overcalls MEDIUM
    # rescue hits (~13 of them) but only one is biologically real (HIGH).
    if conf == "LOW":
        pattern = "url(#lowConfGoi)" if slot["is_goi"] else "url(#lowConfFlank)"
        fill = pattern
        stroke = base
        stroke_dash = ' stroke-dasharray="3,2"'
        opacity = 0.85
        txt_fill = COL_TEXT
    elif conf == "MEDIUM":
        # Visibly weaker than HIGH: heavily lightened fill (~50% toward white)
        # plus a dotted (not solid) border. Reads as "tentative".
        fill = _lerp_hex(base, "#FFFFFF", 0.45)
        stroke = base
        stroke_dash = ' stroke-dasharray="1,1.5"'
        opacity = 0.95
        txt_fill = COL_TEXT
    else:  # HIGH or unknown — treat as solid
        fill = _identity_color(base, ident)
        stroke = base
        stroke_dash = ""
        opacity = 1.0
        txt_fill = "#fff" if ident >= 70 else COL_TEXT

    if is_home:
        title = (f'{slot["label"]} (home reference) — '
                 f'{slot.get("product", "") or slot["id"]}')
    else:
        title = (
            f'{slot["label"]} in target — identity {ident:.1f}%, '
            f'confidence {conf or "—"}, role {cell["role"] or "—"}'
            + (f', evidence {cell["evidence_type"]}' if cell.get("evidence_type") else "")
            + (f', {n_copies}× copies' if n_copies > 1 else "")
        )
    out.append(
        f'<polygon points="{pts}" fill="{fill}" opacity="{opacity:.2f}" '
        f'stroke="{stroke}" stroke-width="0.9"{stroke_dash}>'
        f'<title>{_esc(title)}</title></polygon>'
    )
    # Suppress the identity number on the home row — it's the reference,
    # not a measurement. Showing "100" everywhere is misleading.
    if not is_home and ident >= 30:
        out.append(
            f'<text x="{x + w/2:.1f}" y="{y + h/2 + 3:.1f}" '
            f'text-anchor="middle" font-size="9" fill="{txt_fill}" '
            f'font-weight="{("700" if conf == "HIGH" else "500")}" '
            f'pointer-events="none">{ident:.0f}</text>'
        )
    if n_copies > 1:
        # Small "×N" badge in the upper-right corner.
        out.append(
            f'<text x="{x + w - 3:.1f}" y="{y + 9:.1f}" '
            f'text-anchor="end" font-size="8" fill="{base}" '
            f'font-weight="700" pointer-events="none">×{n_copies}</text>'
        )


def _draw_home_cell(out, x, y, w, h, slot):
    """Home (reference) cell — no identity number, just the strand arrow.

    Distinguishing from target cells avoids the misleading "100" number
    that previously read like a measurement; the home gene is the reference,
    not a hit.
    """
    base = COL_GOI if slot["is_goi"] else COL_FLANKING
    pts = _arrow_path(x, y, w, h, slot["strand"])
    title = f'{slot["label"]} (home reference) — {slot.get("product","")}'
    out.append(
        f'<polygon points="{pts}" fill="{base}" opacity="0.95" '
        f'stroke="{base}" stroke-width="0.8">'
        f'<title>{_esc(title)}</title></polygon>'
    )


def _render_cladogram(rooted_tree, species_rows, margin, matrix_y0, tree_w, cell_h):
    """Real midpoint-rooted cladogram aligned to matrix rows.

    Each leaf is positioned at its species row's vertical center. Internal
    nodes are placed at the mean y of their subtree leaves and at an x
    proportional to cumulative branch length from the root. Branches are
    drawn as right-angle "ladder" lines (horizontal from parent's x to
    child's x, then vertical to span sibling y-range) — the standard
    cladogram convention.

    Species not present in the tree (e.g. those whose only GOI evidence is
    a tandem_copy hit when the upstream filter still drops them) get a
    dashed leader line from their row to the right edge of the tree
    panel, so the reader can tell they are *missing*, not at depth zero.
    """
    if not species_rows:
        return ""

    out = ['<g class="tree" fill="none">']
    x0 = margin + 6
    x1 = margin + tree_w - 8

    if rooted_tree is None or not list(rooted_tree.leaves()):
        # No usable tree — vertical trunk stub with dashed leaders.
        out.append(_dashed_leaders(species_rows, x0, x1, matrix_y0, cell_h))
        out.append("</g>")
        return "\n".join(out)

    # Map each tree leaf-name → species key (matrix row key).
    leaf_to_species = {}
    for leaf in rooted_tree.leaves():
        sp = species_from_leaf(leaf.name) or leaf.name
        leaf_to_species[id(leaf)] = sp

    # Map species_key → row index for matrix alignment.
    row_y = {}
    for ri, row in enumerate(species_rows):
        sp_key = row[3]
        row_y[sp_key] = matrix_y0 + ri * cell_h + cell_h / 2

    # Reorder tree leaves to match the matrix's row order: build a mapping
    # tree_leaf_obj → y, dropping leaves not in the matrix (rare).
    leaf_y = {}
    for leaf in rooted_tree.leaves():
        sp = leaf_to_species[id(leaf)]
        if sp in row_y:
            leaf_y[id(leaf)] = row_y[sp]

    if not leaf_y:
        out.append(_dashed_leaders(species_rows, x0, x1, matrix_y0, cell_h))
        out.append("</g>")
        return "\n".join(out)

    # Compute internal-node depth (cumulative branch length from root).
    depth = {id(rooted_tree): 0.0}
    def _depth_walk(n):
        for c in n.children:
            depth[id(c)] = depth[id(n)] + (c.dist or 0.0)
            _depth_walk(c)
    _depth_walk(rooted_tree)
    max_depth = max(depth.values()) or 1.0
    span_x = x1 - x0 - 6  # leave 6px gutter at the tip

    # Internal-node y = mean y of descendant leaves that have row_y.
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
    # Leaves: their own y (when present)
    for leaf, y in leaf_y.items():
        node_y[leaf] = y

    def _x_for(n):
        return x0 + (depth[id(n)] / max_depth) * span_x

    # Draw branches.
    def _draw(n):
        if n.is_leaf() or id(n) not in node_y:
            return
        ny = node_y[id(n)]
        nx = _x_for(n)
        # Vertical span across drawable children
        child_ys = [node_y[id(c)] for c in n.children if id(c) in node_y]
        if child_ys:
            out.append(
                f'<line x1="{nx:.1f}" y1="{min(child_ys):.1f}" '
                f'x2="{nx:.1f}" y2="{max(child_ys):.1f}" '
                f'stroke="{COL_TREE}" stroke-width="1.2"/>'
            )
        for c in n.children:
            if id(c) not in node_y:
                continue
            cy = node_y[id(c)]
            cx = _x_for(c)
            out.append(
                f'<line x1="{nx:.1f}" y1="{cy:.1f}" '
                f'x2="{cx:.1f}" y2="{cy:.1f}" '
                f'stroke="{COL_TREE}" stroke-width="1.2"/>'
            )
            _draw(c)

    _draw(rooted_tree)

    # Dashed leaders for matrix rows whose species isn't in the tree.
    species_in_tree = set(sp for sp, _ in zip(
        (leaf_to_species[lid] for lid in leaf_y.keys()), leaf_y.keys()
    ))
    missing_rows = [(ri, row[3]) for ri, row in enumerate(species_rows)
                    if row[3] not in species_in_tree and not row[2]]
    for ri, _ in missing_rows:
        y = matrix_y0 + ri * cell_h + cell_h / 2
        out.append(
            f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" '
            f'stroke="{COL_TREE}" stroke-width="0.7" stroke-dasharray="2,3" '
            f'opacity="0.6"/>'
        )

    out.append("</g>")
    return "\n".join(out)


def _dashed_leaders(species_rows, x0, x1, matrix_y0, cell_h):
    parts = []
    for ri, _ in enumerate(species_rows):
        y = matrix_y0 + ri * cell_h + cell_h / 2
        parts.append(
            f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" '
            f'stroke="{COL_TREE}" stroke-width="0.7" stroke-dasharray="2,3" '
            f'opacity="0.6"/>'
        )
    return "\n".join(parts)


def _render_legend(x, y, c_goi, c_flank):
    """Two-row legend: confidence tiers on top, identity scale + role on bottom."""
    out = ['<g class="legend" font-size="10">']

    # Row 1: confidence tiers shown as flanking-coloured arrows.
    cx = x
    out.append(f'<text x="{cx}" y="{y + 10}" font-weight="600">Confidence:</text>')
    cx += 76
    # Match cell rendering: HIGH solid, MEDIUM lightened+dotted, LOW hatched+dashed.
    tiers = [
        ("HIGH",   c_flank,                    "",                       1.0),
        ("MEDIUM", _lerp_hex(c_flank, "#FFFFFF", 0.45),
                                              ' stroke-dasharray="1,1.5"', 0.95),
        ("LOW",    "url(#lowConfFlank)",       ' stroke-dasharray="3,2"',  0.85),
    ]
    for label, fill, dash, opacity in tiers:
        pts = _arrow_path(cx, y, 36, 16, "+")
        out.append(
            f'<polygon points="{pts}" fill="{fill}" opacity="{opacity:.2f}" '
            f'stroke="{c_flank}" stroke-width="0.9"{dash}/>'
        )
        out.append(f'<text x="{cx + 42}" y="{y + 11}">{label}</text>')
        cx += 92

    # Row 2: GOI vs flanking + absent + identity scale.
    y2 = y + 28
    cx = x
    items = [
        (c_goi,   "GOI"),
        (c_flank, "Flanking ortholog"),
        (COL_MISSING_BG, "Absent"),
    ]
    for color, label in items:
        out.append(
            f'<rect x="{cx}" y="{y2}" width="18" height="12" fill="{color}" '
            f'stroke="#888" stroke-width="0.5"/>'
        )
        out.append(f'<text x="{cx + 24}" y="{y2 + 10}">{_esc(label)}</text>')
        cx += 24 + 8 + max(60, len(label) * 7)

    # Identity gradient bar.
    bar_x = cx + 12
    bar_y = y2 - 1
    for k in range(20):
        t = k / 19.0
        col = _lerp_hex(c_flank, "#FFFFFF", min(0.55, (1 - t) * 0.7))
        out.append(f'<rect x="{bar_x + k * 6}" y="{bar_y}" width="6" height="14" fill="{col}"/>')
    out.append(f'<text x="{bar_x - 4}" y="{y2 + 26}" font-size="9">30%</text>')
    out.append(f'<text x="{bar_x + 105}" y="{y2 + 26}" font-size="9">100% identity</text>')

    out.append('</g>')
    return "\n".join(out)


# ─────────────────────────────── HTML wrapper ────────────────────────────────

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SynVoy synteny matrix — {title}</title>
<style>
  body {{ margin: 24px; font-family: Helvetica, Arial, sans-serif; color: #222; }}
  h1   {{ font-size: 18px; font-weight: 600; margin: 0 0 4px 0; }}
  p.sub {{ color: #666; margin: 0 0 20px 0; }}
  .wrap {{ overflow-x: auto; padding-bottom: 20px; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="sub">{subtitle}</p>
<div class="wrap">{svg}</div>
</body>
</html>
"""


def write_output(svg, title, subtitle, path):
    if path.endswith(".svg"):
        with open(path, "w") as fh:
            fh.write(svg)
        return
    with open(path, "w") as fh:
        fh.write(HTML_TEMPLATE.format(svg=svg, title=_esc(title),
                                      subtitle=_esc(subtitle)))


# ─────────────────────────────── orchestrator ───────────────────────────────

def _guess_candidate_bed(gff_path):
    """Given `Apis_cerana.fna.gff`, look for `Apis_cerana.fna.regions.bed` next
    to it. Return the path if it exists, else None.
    """
    if not gff_path:
        return None
    if gff_path.endswith(".gff"):
        candidate = gff_path[:-4] + ".regions.bed"
        if os.path.exists(candidate):
            return candidate
    return None


def derive_target_name(path):
    base = os.path.basename(path)
    # Strip the longest matching suffix so "Tetragonula_carbonaria.fa.gff"
    # becomes "Tetragonula_carbonaria", not "Tetragonula_carbonaria.fa".
    for suf in (".fasta.gff", ".fna.gff", ".fa.gff", ".gff"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return base.replace("_", " ")


def _parse_species_map(path):
    """Parse a SynVoy species-mapping TSV into {accession: species_name}.

    Format (one row per genome): accession<TAB>species<TAB>role[<TAB>...].
    Empty / missing path is tolerated and yields an empty mapping.
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            acc, species = parts[0].strip(), parts[1].strip()
            if acc and species:
                out[acc] = species
    return out


def _species_for_stem(stem, species_by_accession):
    """Look up a genome stem (e.g. `GCF_049350105.2.fna`) in the species map,
    progressively stripping FASTA suffixes (`.fna`/`.fa`/`.fasta`) so a stem
    matches the bare accession key the species map uses.
    """
    if not stem or not species_by_accession:
        return None
    candidates = [stem]
    for suf in (".fna", ".fa", ".fasta"):
        if stem.endswith(suf):
            candidates.append(stem[: -len(suf)])
            break
    for key in candidates:
        if key in species_by_accession:
            return species_by_accession[key]
    return None


def _matching_stem(path_or_name):
    """Canonical genome stem used to join GFFs, candidate BEDs, homology
    TSVs and display names across separate Nextflow channels.

    Nextflow does not guarantee positional alignment of independently
    collected channels, so these four lists can arrive in different orders
    for the same locus. Matching by stem fixes that.

    All of these resolve to the same stem:
        GCF_000001635.27.fna.gff           → GCF_000001635.27.fna
        GCF_000001635.27.fna.regions.bed   → GCF_000001635.27.fna
        GCF_000001635.27.fna.homology.tsv  → GCF_000001635.27.fna
        GCF_000001635.27.fna               → GCF_000001635.27.fna
    """
    base = os.path.basename(path_or_name or "")
    for suf in (".regions.bed", ".homology.tsv", ".gff"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base


def order_targets(target_specs, tree_leaves):
    """Order targets by phylogenetic tree leaf order; unknowns alphabetically
    at the end.

    Handles both leaf-name styles SynVoy uses:
      - GOI tree: 'GOI_Melt|Apis_florea_fna_b0_l1_exon_ann' (paralog-tagged).
        ``species_from_leaf`` strips the GOI prefix and locus tags down to
        'Apis_florea'.
      - Taxonomy tree: bare 'Apis_florea' — no '|' separator. Fall back to
        the leaf name itself in that case (``species_from_leaf`` returns
        None on these).
    """
    leaf_species_order = []
    seen = set()
    for leaf in tree_leaves:
        sp = species_from_leaf(leaf) or leaf
        if sp and sp not in seen:
            leaf_species_order.append(sp)
            seen.add(sp)
    def key_for(spec):
        underscore = spec["display"].replace(" ", "_")
        try:
            # Known leaf: stays in tree order. Empty-string secondary keeps
            # all keys 2-tuples so int-vs-tuple comparisons don't blow up.
            return (leaf_species_order.index(underscore), "")
        except ValueError:
            # Unknown: alphabetical at the bottom — deterministic, no hashing.
            return (10_000, spec["display"].lower())
    return sorted(target_specs, key=key_for)


def species_clade_colors(rooted_tree, clade_count):
    """Map species → clade color via partition_clades on the rooted tree.

    Tree leaves carry GOI-prefixed names; we extract the species component
    and dedupe (a species with multiple paralogs gets the color of its
    earliest-occurring leaf in left-to-right order). Partitioning uses
    the topology-driven iterative-split rule (see
    ``synvoy_tree.partition_clades``) with K = ``clade_count``.
    """
    if rooted_tree is None:
        return {}
    leaf_to_clade = partition_clades(rooted_tree, target_k=clade_count)
    species_to_clade = {}
    for leaf_name, clade_id in leaf_to_clade.items():
        sp = species_from_leaf(leaf_name)
        # Some leaves are bare GOI labels (e.g. "GOI_Melt") for the home query.
        if not sp:
            sp = leaf_name
        species_to_clade.setdefault(sp, clade_id)
    return {sp: CLADE_PALETTE[cid % len(CLADE_PALETTE)]
            for sp, cid in species_to_clade.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--home_bed", required=True)
    ap.add_argument("--query_bed", required=True)
    ap.add_argument("--home_gff", default="")
    ap.add_argument("--target_gffs", nargs="*", default=[])
    ap.add_argument("--target_names", nargs="*", default=[])
    ap.add_argument("--homology_tsvs", nargs="*", default=[])
    ap.add_argument("--candidate_beds", nargs="*", default=[],
                    help="Per-target regions.bed files. Cells outside these "
                         "intervals are dropped from the matrix — matches "
                         "the legacy plot's behaviour and prevents off-region "
                         "hits from inflating the apparent ortholog count.")
    ap.add_argument("--tree", default="")
    ap.add_argument("--species_map", default="",
                    help="TSV: accession<TAB>species<TAB>role. Used to render "
                         "target row labels with their real species name "
                         "instead of the raw genome accession.")
    ap.add_argument("--home_species", default="Home")
    ap.add_argument("--title", default="")
    ap.add_argument("--output", required=True)
    ap.add_argument("--clade_count", type=int, default=4,
                    help="Number of clades for row colouring (default 4). "
                         "Iteratively splits the largest clade of the "
                         "midpoint-rooted tree (topology-driven; "
                         "singleton-producing splits deprioritised). "
                         "Replaces the legacy --clade_depth_frac.")
    ap.add_argument("--tree_source", choices=("taxonomy", "goi"),
                    default="taxonomy",
                    help="Which tree to render alongside the matrix. "
                         "'taxonomy' (default) builds an NCBI-classification "
                         "tree from the species list so every searched "
                         "species gets a leaf — even ones with only weak "
                         "GOI hits. 'goi' uses the IQ-TREE locus tree at "
                         "--tree (legacy).")
    # Legacy flag kept so existing Nextflow invocations don't crash.
    ap.add_argument("--clade_depth_frac", type=float, default=None,
                    help=argparse.SUPPRESS)
    ap.add_argument("--common_names", choices=("both", "common", "scientific", "off"),
                    default="both",
                    help="Species label style. 'both' shows 'Scientific (common)'.")
    ap.add_argument("--common_names_tsv", default="",
                    help="Optional 2-column TSV (scientific<TAB>common) "
                         "overriding NCBI lookups.")
    ap.add_argument("--no_network", action="store_true",
                    help="Skip the NCBI 'datasets' CLI lookup for common names.")
    args = ap.parse_args()

    # Initialize common-name resolver up front (cheap; just reads cache).
    if args.common_names != "off":
        synvoy_taxa.init_lookup(common_names_tsv=args.common_names_tsv or None,
                                allow_network=not args.no_network)

    home_bed = parse_bed(args.home_bed)
    query_bed = parse_bed(args.query_bed)
    home_gff_index = parse_home_gff_genes(args.home_gff)

    if not home_bed:
        sys.exit(f"ERROR: no rows in --home_bed {args.home_bed}")

    slots, goi_slot_id = build_home_slots(home_bed, query_bed, home_gff_index)
    if not goi_slot_id:
        print(f"WARN: no GOI slot synthesized — query_bed empty or unmatched. "
              f"Home row will lack a GOI cell.", file=sys.stderr)

    # Pair target GFFs with their matching name, homology TSV and candidate
    # BED *by filename stem*. Naive positional zip is unsafe: main.nf collects
    # target_gffs / target_names / candidate_beds / homology_tsvs from four
    # independent Nextflow channels that are not guaranteed to be in the same
    # order, so a positional pairing silently mis-labels rows and renders the
    # wrong region/homology data per species.
    names_by_stem = {_matching_stem(n): n for n in args.target_names if n}
    beds_by_stem  = {_matching_stem(b): b for b in args.candidate_beds if b}
    homs_by_stem  = {_matching_stem(h): h for h in args.homology_tsvs if h}
    species_by_accession = _parse_species_map(args.species_map)

    target_specs = []
    for gff in args.target_gffs:
        stem = _matching_stem(gff)
        name = names_by_stem.get(stem)
        bed  = beds_by_stem.get(stem) or _guess_candidate_bed(gff)
        hom  = homs_by_stem.get(stem)
        # Prefer the species name from the SynVoy species map (real binomial),
        # otherwise fall back to the user-supplied target_name, otherwise to
        # the GFF stem with underscores → spaces.
        species_label = _species_for_stem(stem, species_by_accession)
        display = species_label or name or derive_target_name(gff).replace("_", " ")
        target_specs.append({"gff": gff, "display": display,
                             "homology": hom, "candidate_bed": bed})

    # ---- Build the tree the matrix renders alongside its rows ----
    # By default we use an NCBI-taxonomy species tree (ete3.NCBITaxa) so
    # every searched species — including those with only weak GOI hits and
    # therefore absent from the IQ-TREE locus tree — gets a leaf and a
    # clade colour. That makes the matrix's left-side cladogram show
    # "who is related to who" rather than per-locus sequence similarity.
    # `--tree_source goi` falls back to the legacy IQ-TREE locus tree.
    tree_source = getattr(args, "tree_source", "taxonomy")

    rooted_tree = None
    if tree_source == "taxonomy":
        target_sp = [spec["display"] for spec in target_specs]
        sp_list = [args.home_species] + target_sp
        tax_newick = build_taxonomy_tree_newick(
            sp_list, allow_network=not args.no_network)
        if tax_newick:
            raw_tax = parse_newick_tree(tax_newick)
            # NCBI topology has no real branch lengths; midpoint-root is a
            # no-op on unit branches but harmless.
            rooted_tree = midpoint_root(raw_tax) if raw_tax else None
        else:
            print("[matrix] Taxonomy tree unavailable (ete3 missing or no "
                  "taxids resolved); falling back to GOI tree.", file=sys.stderr)
            tree_source = "goi"

    if rooted_tree is None:
        # GOI-tree path (legacy): parse IQ-TREE newick, midpoint-root, and
        # collapse multi-paralog leaves down to one per species so the
        # cladogram renders one branch per matrix row.
        tree_text = _read_text(args.tree) if args.tree else ""
        raw_tree = parse_newick_tree(tree_text)
        rooted_tree = midpoint_root(raw_tree) if raw_tree else None
        if rooted_tree is not None:
            rooted_tree = collapse_to_one_leaf_per_species(
                rooted_tree, species_from_leaf)
        tree_leaves = ([n.name for n in rooted_tree.leaves() if n.name]
                       if rooted_tree else parse_newick_leaf_order(tree_text))
    else:
        tree_leaves = [n.name for n in rooted_tree.leaves() if n.name]

    target_specs = order_targets(target_specs, tree_leaves)
    sp_clade_color = species_clade_colors(rooted_tree, args.clade_count)

    # Build species_rows. Each entry: (display_label, cell_summary, is_home,
    # tree_species_key, clade_color).
    species_rows = []

    home_label = _label_species(args.home_species, args.common_names)
    home_clade = sp_clade_color.get(args.home_species.replace(" ", "_"))
    species_rows.append((home_label, _home_summary(slots, goi_slot_id), True,
                         args.home_species.replace(" ", "_"), home_clade))

    for spec in target_specs:
        regions = parse_candidate_regions(spec["candidate_bed"])
        cells = parse_target_gff(spec["gff"], candidate_regions=regions)
        homology = parse_homology_tsv(spec["homology"]) if spec["homology"] else []
        summary = map_target_cells_to_slots(cells, slots, homology)
        # Pass through both the best cell and the count of tandem hits per slot
        # so the renderer can show "×N" badges for multi-copy GOI cells.
        condensed = {}
        for sid, v in summary.items():
            best = v["best"].copy()
            best["_n_copies"] = v["n"]
            condensed[sid] = best
        sp_key = spec["display"].replace(" ", "_")
        label = _label_species(spec["display"], args.common_names)
        clade_color = sp_clade_color.get(sp_key)
        species_rows.append((label, condensed, False, sp_key, clade_color))

    title = args.title or f"{args.home_species} — synteny matrix"
    subtitle = (
        f'{len(slots)} home slots ({"GOI: " + next((s["label"] for s in slots if s["is_goi"]), "—")})'
        f' · {len(species_rows) - 1} target genomes'
    )
    svg = render_svg(slots, species_rows, args.home_species, goi_slot_id,
                     rooted_tree=rooted_tree)
    write_output(svg, title, subtitle, args.output)
    print(f"Wrote {args.output}", file=sys.stderr)


def _label_species(scientific, mode):
    """Wrap synvoy_taxa.label_for_species and tolerate the 'off' mode."""
    if mode == "off":
        return scientific.replace("_", " ")
    return synvoy_taxa.label_for_species(scientific.replace("_", " "), mode=mode)


def _home_summary(slots, goi_slot_id):
    """Home row: every slot present at 100% identity, HIGH confidence, with strand."""
    summary = {}
    for s in slots:
        summary[s["id"]] = {
            "parent": s["id"], "role": "goi" if s["is_goi"] else "flanking",
            "confidence": "HIGH", "identity": 100.0,
            "strand": s["strand"], "start": s["start"], "end": s["end"],
            "name": s["label"], "chrom": "",
            "evidence_type": "home", "model_status": "complete", "exons": "",
        }
    return summary


def _read_text(path):
    if not path or not os.path.exists(path):
        return ""
    with open(path) as fh:
        return fh.read()


if __name__ == "__main__":
    main()
