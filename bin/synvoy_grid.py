#!/usr/bin/env python3
"""The anchor grid: one row per genome, one column per home gene (import-only helper).

Everything is laid out in points at the size the figure is printed. The print
version states its width in millimetres (183 mm, text 7 / 6 / 5 pt); when a grid
is too large the spacing gives way first, then the columns with the fewest
orthologs are dropped and reported. Text is never made smaller than 5 pt. The
full version is the same drawing without a width limit, with the identity in
every arrow, the home-coordinate axis, the genomic location of every row and a
tooltip on every glyph; ``plot_synteny.py`` shows it enlarged on screen and adds
the gene-names table.

GOI column. SynVoy reports GOI *calls*. A call is a *copy* when it is the best
call of its genome or has at least medium confidence; every other call is a
*weak call*. The grid shows a small arrow per copy, in genomic order and
pointing in its coding direction, numbered along the chromosome when a genome
has several; a small pale mark per weak call; and ``//`` before what is not at
the locus. The exon structure of every copy is drawn on a line of its own,
5'->3', in a panel right of the grid, under the same number. When no genome
(the home genome included) has two copies there is no panel: the model is drawn
in the GOI column itself, in its genomic orientation.

    result = render(plot_synteny_module, tracks, gene_colours, home_products, args, style_for())

``plot_synteny`` is passed in as a module: it holds the gene-name registries,
and importing it a second time from here would give an empty copy of them.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, replace
from typing import Optional

import synvoy_legend as _legend
import synvoy_tree as _stree

PT_PER_MM = 72.0 / 25.4
FONT = "Arial, Helvetica, sans-serif"
GOI_KEY = "__GOI__"

INK = "#1a1a1a"
MUTED = "#555b66"
TREE_LINE = "#7b8494"
THREAD = "#c3c9d4"
STRIPE = "#f3f4f7"
HOME_BAND = "#e4e8ef"
AXIS = "#9aa3b2"
GAP_INK = "#b45309"
DASH = {"MEDIUM": "1.5,0.8", "LOW": "0.7,0.7", "AMBIGUOUS": "0.45,1.0"}

MIN_EXON = 0.9       # an exon is never drawn narrower
PAD = 3.0            # canvas margin
CELL_GAP = 2.2       # white space between two arrows of a row
MINI_W, MINI_H = 8.6, 5.4     # small arrow standing for one copy (room for its number)
MARK_W = 3.2         # small mark standing for one weak call
ITEM_GAP = 2.2       # between two glyphs of a GOI cell
LINE_UNITS = 6.0     # a line of a GOI cell holds six small arrows (a mark counts half)
SCALE_STEPS = (30, 60, 90, 150, 300, 600, 900, 1500, 3000, 6000, 9000, 15000)
NUMBER_MIN_IDENTITY = 25.0   # no number is written into a paler arrow than this
COV_LINE_H = 5.2     # room under an arrow for the coverage line "(44)" in 5 pt

_PS = None           # the plot_synteny module of the current render() call


@dataclass
class GridStyle:
    numbers: bool = False        # the identity printed in every flanking arrow
    full: bool = False           # title, home axis, location column, tooltips; no size limit
    width_mm: Optional[float] = 183.0
    max_height_mm: Optional[float] = 170.0
    fs_title: float = 8.0
    fs_species: float = 7.0
    fs_species_min: float = 6.0
    fs_col: float = 6.0
    fs_num: float = 6.0
    fs_small: float = 5.0
    col_w: float = 12.0
    col_w_min: float = 9.5
    col_w_max: float = 15.0
    row_h: float = 10.0
    row_h_min: float = 8.6
    arrow_h: float = 5.6
    lane_h: float = 6.4          # line pitch inside a GOI cell with several lines
    lane_h_min: float = 5.8
    goi_model_w: float = 60.0    # exons of the longest model
    goi_model_w_min: float = 36.0
    goi_model_w_max: float = 84.0
    tree_w: float = 28.0
    angle: float = 55.0
    max_models: int = 10         # most exon models drawn for one genome
    max_marks: int = 8           # most weak-call marks in one GOI cell; the rest is "+N"
    max_cols: int = 50


def style_for(numbers=False, full=False, **kw):
    """The style of the print version, or with ``full`` of the supplementary one."""
    st = GridStyle(numbers=bool(numbers or full), full=full)
    if st.numbers:      # a 5 pt number needs a taller, wider arrow
        st = replace(st, col_w=14.0, col_w_min=13.0, col_w_max=15.5,
                     row_h=11.0, row_h_min=10.0, arrow_h=7.2)
    if full:
        st = replace(st, width_mm=None, max_height_mm=None)
    return replace(st, **kw)


def tw(text, size, bold=False):
    return _PS.text_width(text, size, bold=bold)


def esc(text):
    return _PS._svg_esc(text)


def text_svg(x, y, size, label, fill=INK, bold=False, anchor="", extra=""):
    return '<text x="%.2f" y="%.2f" font-size="%g" fill="%s"%s%s%s>%s</text>' % (
        x, y, size, fill, ' font-weight="700"' if bold else "",
        f' text-anchor="{anchor}"' if anchor else "", extra, esc(label))


# ------------------------------------------------------------------ data ----

def build_anchors(home_track, gene_colours, goi_label, home_products):
    """Home genes in genomic order, the GOI as one column."""
    ps = _PS
    anchors, seen = [], set()
    for hg in sorted(home_track.get("genes", []), key=lambda g: g.get("start", 0)):
        nm = hg.get("name", "")
        is_g = ps.is_goi(nm) or ps.is_goi(hg.get("home_gene_id", "") or "")
        key = GOI_KEY if is_g else nm
        if key in seen:
            continue
        seen.add(key)
        anchors.append({
            "key": key, "is_goi": is_g,
            "label": goi_label if is_g else ps._home_gene_label(nm),
            "colour": ps.GOI_COLOUR if is_g else gene_colours.get(nm, ps.UNMATCHED_CLR),
            "strand": hg.get("strand", "+"), "start": hg.get("start", 0),
            "end": hg.get("end", 0), "chrom": hg.get("chrom", ""),
            "gene_id": "" if is_g else ps.clean_gene_label(nm),
            "product": ((ps._named(nm) or {}).get("full_name")
                        or (ps._lookup_product(nm, home_products) if home_products else "")),
        })
    return anchors


def row_label(track):
    return (track.get("species") or re.sub(r"<[^>]+>", "", track.get("label") or "")).strip()


class _Node:
    __slots__ = ("kids", "track", "h", "x", "y")

    def __init__(self, kids=None, track=None):
        self.kids, self.track, self.h, self.x, self.y = kids or [], track, 0, 0.0, 0.0


def species_tree_rows(args, tracks):
    """Row order and a pruned tree whose leaves are the rows.

    Returns (root or None, row order as track indices, is_species_tree). The
    home genome is row 0. A row is matched to a leaf by its genome id or by its
    species name, so a run whose genome ids are not species names ('chicken')
    still gets its tree. When the home genome has a leaf, the tree is rotated so
    that leaf comes first: the home row is then part of the cladogram instead of
    hanging from the root on a dashed line.
    """
    ps = _PS
    rooted, _order, is_species = ps._grid_species_tree(args)
    n = len(tracks)
    if rooted is None:
        return None, list(range(n)), False
    keys = {}
    for lf in rooted.leaves():
        for k in {lf.name, _stree.species_from_leaf(lf.name) or lf.name}:
            if k:
                keys.setdefault(k, lf)
    by_leaf = {}
    for ti, tr in enumerate(tracks):
        cands = [row_label(tr).replace(" ", "_")]
        gid = tr.get("genome_id") or ""
        if gid and gid != "home":
            cands.insert(0, gid)
        leaf = next((keys[c] for c in cands if c in keys), None)
        if leaf is None and gid and gid != "home":
            hit = ps._grid_match_species(gid, list(keys))
            leaf = keys.get(hit) if hit else None
        if leaf is not None:
            by_leaf.setdefault(id(leaf), []).append(ti)

    def prune(node):
        if node.is_leaf():
            tis = by_leaf.get(id(node), [])
            if not tis:
                return None
            if len(tis) == 1:
                return _Node(track=tis[0])
            return _Node(kids=[_Node(track=t) for t in tis])
        kids = [k for k in (prune(c) for c in node.children) if k is not None]
        if not kids:
            return None
        return kids[0] if len(kids) == 1 else _Node(kids=kids)

    root = prune(rooted)
    if root is None:
        return None, list(range(n)), is_species

    def rotate_home_first(node):
        if node.track is not None:
            return node.track == 0
        for i, k in enumerate(node.kids):
            if rotate_home_first(k):
                node.kids.insert(0, node.kids.pop(i))
                return True
        return False

    rotate_home_first(root)
    order = []

    def walk(node):
        if node.track is not None:
            order.append(node.track)
        for k in node.kids:
            walk(k)

    walk(root)
    placed = set(order)
    rest = [ti for ti in range(n) if ti not in placed]
    if 0 in rest:                      # home has no leaf: still the first row
        rest.remove(0)
        order = [0] + order
    elif order and order[0] != 0:      # a tree that cannot be rotated to the home leaf
        order.remove(0)
        order = [0] + order
    return root, order + rest, is_species


def model_blocks(gene):
    return sorted(gene.get("exon_coords") or []) or [(gene["start"], gene["end"])]


def cds_bp(gene):
    return sum(e - s + 1 for s, e in model_blocks(gene))


def caret_for(model_w):
    """Width of an intron caret: narrower when the models are drawn small."""
    return max(1.4, min(2.4, model_w / 25.0))


def model_geometry(gene, bp_pt, cap_bp=None, caret=2.4, flip=False):
    """(width, exon boxes, intron spans) of a GOI model: exons to scale, carets between.

    A hit-based call has no exon structure and is drawn as the span of its hits;
    one tandem-copy call spanned 5.6 kb next to a 227 bp gene. Such a bar is cut
    to `cap_bp`, the longest gene model of the figure, so it cannot set the scale.
    """
    lens = [e - s + 1 for s, e in model_blocks(gene)]
    if cap_bp and not _PS._is_spliced_model(gene) and sum(lens) > cap_bp:
        lens = [n * cap_bp / sum(lens) for n in lens]
    if flip:
        lens = lens[::-1]
    x, exons, introns = 0.0, [], []
    for i, n_bp in enumerate(lens):
        if i:
            introns.append((x, x + caret))
            x += caret
        w = max(MIN_EXON, n_bp * bp_pt)
        exons.append((x, x + w))
        x += w
    return x, exons, introns


def is_copy(gene, best):
    """A call counts as a gene copy: the best call of its genome, or at least medium."""
    return _PS._goi_is_copy(gene, best)


def split_calls(tmap, span_bp):
    """The GOI calls of one target row: (best, [(gene, elsewhere, is_copy)] in genomic order).

    Calls at the locus come first, in the order of the oriented chromosome; calls
    on another scaffold, or farther away than the neighbourhood is long, follow.
    """
    ps = _PS
    genes = (tmap.get(GOI_KEY) or {}).get("genes", [])
    if not genes:
        return None, []
    best = ps._goi_best([g for g in genes if not ps._is_fragment_goi(g)] or genes)
    _, lo, hi, _ = ps._dominant_chrom_span(
        [(e["best"].get("chrom", ""), (e["best"]["start"] + e["best"]["end"]) / 2.0)
         for e in tmap.values()], best.get("chrom", ""))
    reach = max(float(span_bp or 0), hi - lo)
    local = sorted((g for g in genes if ps._is_local_copy(g, best, reach)),
                   key=lambda g: g.get("start_plot", g["start"]))
    away = sorted((g for g in genes if not ps._is_local_copy(g, best, reach)),
                  key=lambda g: (g.get("chrom", ""), g["start"]))
    return best, ([(g, False, is_copy(g, best)) for g in local]
                  + [(g, True, is_copy(g, best)) for g in away])


def wrap_items(items):
    """Lines of a GOI cell: at most LINE_UNITS small arrows on one (a mark counts half)."""
    weight = {"mini": 1.0, "mark": 0.5, "break": 0.5, "full": 0.0, "count": 1.0}
    lines, used = [[]], 0.0
    for it in items:
        w = weight[it[0]]
        if lines[-1] and used + w > LINE_UNITS:
            carry = [lines[-1].pop()] if lines[-1][-1][0] == "break" else []
            lines.append(carry)
            used = 0.5 * len(carry)
        lines[-1].append(it)
        used += w
    return [ln for ln in lines if ln]


# ---------------------------------------------------------------- colours ----

def cell_style(base, identity, conf, is_home, numbers):
    """(fill, stroke, dash) of one flanking arrow."""
    ps = _PS
    conf = (conf or "").upper()
    dash = DASH.get(conf, "")
    if is_home:
        return base, ps._darken_hex(base, 0.7), ""
    if numbers:                 # pale fills: the black number must stay readable
        t = {"MEDIUM": 0.66, "LOW": 0.8, "AMBIGUOUS": 0.86}.get(conf, 0.5)
        return ps._lerp_hex(base, "#ffffff", t), ps._darken_hex(base, 0.7), dash
    if conf == "MEDIUM":
        return ps._lerp_hex(base, "#ffffff", 0.42), base, dash
    if conf in ("LOW", "AMBIGUOUS"):
        return ps._lerp_hex(base, "#ffffff", 0.62), base, dash
    return ps._shade_by_identity(base, identity), ps._darken_hex(base, 0.65), ""


def short_coverage(gene, numbers):
    """Query coverage in whole % where the figure has to state it, else None.

    Stated where the identity is printed and the ortholog covers less of the home
    protein than COVERAGE_FLAG_THRESHOLD: 100 % identity over 7 % of a protein is
    another finding than 100 % over all of it, and a bare flag hides which it is.
    """
    if not numbers or not gene:
        return None
    cov = gene.get("query_coverage")
    if cov is None or cov >= _PS.COVERAGE_FLAG_THRESHOLD:
        return None
    if (gene.get("identity") or 0.0) < NUMBER_MIN_IDENTITY:
        return None
    return _PS._coverage_pct(cov)


def goi_fill(identity, conf, is_home, colour):
    """(fill, dash) of one GOI glyph."""
    ps = _PS
    conf = (conf or "").upper()
    if is_home:
        return colour, ""
    if conf == "AMBIGUOUS":
        return ps._lerp_hex(colour, "#ffffff", 0.72), DASH["AMBIGUOUS"]
    if conf == "LOW":
        return ps._lerp_hex(colour, "#ffffff", 0.55), DASH["LOW"]
    if conf == "MEDIUM":
        return ps._lerp_hex(colour, "#ffffff", 0.36), DASH["MEDIUM"]
    return ps._shade_by_identity(colour, max(identity or 0, 55)), ""


# ---------------------------------------------------------------- tooltips ----

def goi_tooltip(gene, label, is_home, note=""):
    blocks = model_blocks(gene)
    ident = gene.get("identity", 0.0) or 0.0
    cov = gene.get("query_coverage")
    conf = (gene.get("confidence") or "").upper()
    return (f'{gene.get("name", "GOI")} in {label}: '
            f'{len(blocks)} exon{"s" if len(blocks) != 1 else ""}, CDS {cds_bp(gene):,} bp, '
            f'{gene.get("chrom", "")}:{blocks[0][0]:,}-{blocks[-1][1]:,} '
            f'({gene.get("genomic_strand", gene.get("strand", "+"))}'
            f'{", scaffold drawn reversed" if "genomic_strand" in gene else ""})'
            + (f', identity {ident:.0f}%' if (ident and not is_home) else '')
            + (f', coverage {_PS._coverage_pct(cov)}%' if (cov is not None and not is_home) else '')
            + (f', {conf}' if conf else '')
            + (f', {gene.get("evidence_type")}' if gene.get("evidence_type") else '')
            + (f' ({note})' if note else ''))


def titled(tooltip, body, on):
    """Wrap a glyph with its tooltip (full version only)."""
    if not on:
        return body
    return f'<g class="acell"><title>{esc(tooltip)}</title>{body}</g>'


# ----------------------------------------------------------------- render ----

def render(ps, all_tracks, gene_colours, home_products, args, st):
    """Draw one grid.

    Returns a dict: ``parts`` (SVG elements in points, no <svg> wrapper), ``width``
    and ``height`` in points, ``info``, ``source`` (one row per cell and call),
    ``legend_text`` (a draft figure legend) and ``legend_groups``.
    """
    global _PS
    _PS = ps
    goi_label = ps._goi_display_label(all_tracks)
    goi_colour, goi_border = ps.GOI_COLOUR, ps.GOI_BORDER
    break_w = tw("//", st.fs_small) + 1.0
    tips = st.full

    tree_root, order, is_species_tree = species_tree_rows(args, all_tracks)
    rows = [all_tracks[ti] for ti in order]            # rows[0] is the home genome
    home_track = rows[0]
    n_rows = len(rows)

    anchors = build_anchors(home_track, gene_colours, goi_label, home_products)
    keys = {a["key"] for a in anchors}
    span_bp = ps._home_locus_span(anchors)
    tmaps = [None] + [ps._grid_target_map(t, keys, GOI_KEY, span_bp) for t in rows[1:]]
    cov = {a["key"]: sum(1 for m in tmaps[1:] if a["key"] in m) for a in anchors}
    n_home_genes = sum(1 for a in anchors if not a["is_goi"])
    anchors = [a for a in anchors if a["is_goi"] or cov[a["key"]] > 0]
    n_unplaced_cols = n_home_genes - sum(1 for a in anchors if not a["is_goi"])
    capped = []
    if len([a for a in anchors if not a["is_goi"]]) > st.max_cols:
        ranked = sorted((a for a in anchors if not a["is_goi"]), key=lambda a: -cov[a["key"]])
        capped = ranked[st.max_cols:]
        gone = {a["key"] for a in capped}
        anchors = [a for a in anchors if a["key"] not in gone]

    # ---- GOI calls of every row ---------------------------------------------------
    home_goi = sorted([g for g in home_track.get("genes", [])
                       if ps.is_goi(g.get("name", "") or "")
                       or ps.is_goi(g.get("home_gene_id", "") or "")],
                      key=lambda g: g.get("start", 0))
    calls = [[(g, False, True) for g in home_goi]]      # per row: (gene, elsewhere, is_copy)
    bests = [None]
    for ri in range(1, n_rows):
        best, row_calls = split_calls(tmaps[ri], span_bp)
        bests.append(best)
        calls.append(row_calls)
    n_copies = [sum(1 for _, _, c in rc if c) for rc in calls]
    in_column = max(n_copies + [0]) <= 1     # no genome has two: the model goes in the column
    numbered = [(not in_column) and n > 1 for n in n_copies]

    cells = []           # per row: grid lines, panel lines, copies beyond max_models
    for ri in range(n_rows):
        items, seen_away, index = [], False, {}
        shown = [g for g, _, c in calls[ri] if c][:st.max_models]
        if bests[ri] is not None and not any(g is bests[ri] for g in shown) and shown:
            shown[-1] = bests[ri]
        shown_ids = {id(g) for g in shown}
        n = n_marks = 0
        for g, away, copy in calls[ri]:
            if not copy:
                n_marks += 1
                if n_marks > st.max_marks:
                    continue
            if away and not seen_away:
                items.append(("break", None, True))
                seen_away = True
            if copy:
                n += 1
                index[id(g)] = n
                items.append(("full" if in_column else "mini", g, away))
            else:
                items.append(("mark", g, away))
        if n_marks > st.max_marks:      # a family with many chance hits: the rest as a number
            items.append(("count", n_marks - st.max_marks, False))
        panel = [] if in_column else [(g, away) for g, away, c in calls[ri]
                                      if c and id(g) in shown_ids]
        cells.append({"grid": wrap_items(items), "panel": panel, "index": index,
                      "more": 0 if in_column else n_copies[ri] - len(panel)})
    has_goi_col = any(a["is_goi"] for a in anchors)
    full_drawn = [g for ri in range(n_rows) for g, _, c in calls[ri] if c
                  and (in_column or any(g is p for p, _ in cells[ri]["panel"]))]
    has_panel = (not in_column) and bool(full_drawn) and has_goi_col
    # The scale comes from the gene models (home gene, miniprot models), not from hit spans.
    max_cds = max([cds_bp(g) for g in ([g for g in full_drawn if ps._is_spliced_model(g)]
                                       or full_drawn)] + [1])
    mini_w = MINI_W + (2.2 if max(n_copies + [0]) >= 10 else 0.0)    # two-digit copy numbers

    def number_of(gene, ri):
        return "" if ri == 0 else ps._goi_number_label(gene)

    def item_w(it, bp_pt, caret, ri):
        if it[0] == "full":
            w = model_geometry(it[1], bp_pt, max_cds, caret)[0]
            num = number_of(it[1], ri)
            return w + ((2.0 + tw(num, st.fs_num, bold=True)) if num else 0.0)
        if it[0] == "count":
            return tw(f"+{it[1]}", st.fs_small)
        return {"mini": mini_w, "mark": MARK_W, "break": break_w}[it[0]]

    def line_w(line, bp_pt, caret, ri):
        return sum(item_w(it, bp_pt, caret, ri) for it in line) + ITEM_GAP * (len(line) - 1)

    def panel_line_w(gene, away, bp_pt, caret, ri):
        w = (break_w + ITEM_GAP if away else 0.0) + model_geometry(gene, bp_pt, max_cds, caret)[0]
        num = number_of(gene, ri)
        return w + ((2.0 + tw(num, st.fs_num, bold=True)) if num else 0.0)

    # ---- genomic location of every row (full version) ---------------------------
    def span_text(ri):
        if ri == 0:
            mids = [(a["start"] + a["end"]) / 2.0 for a in anchors if a["start"]]
            chrom = next((a["chrom"] for a in anchors if a.get("chrom")), "")
            return ps._span_gutter_label(chrom, min(mids) if mids else 0, max(mids) if mids else 0, 0)
        items, goi_chrom = [], ""
        for a in anchors:
            entry = tmaps[ri].get(a["key"])
            if not entry:
                continue
            ch = entry["best"].get("chrom", "")
            items.append((ch, (entry["best"]["start"] + entry["best"]["end"]) / 2.0))
            if a["is_goi"] and ch:
                goi_chrom = ch
        return ps._span_gutter_label(*ps._dominant_chrom_span(items, goi_chrom))
    spans = [span_text(ri) for ri in range(n_rows)] if st.full else []

    labels = [row_label(t) or ("Home" if ri == 0 else "") for ri, t in enumerate(rows)]
    tree_w = st.tree_w if tree_root is not None else 0.0
    cos_a, sin_a = math.cos(math.radians(st.angle)), math.sin(math.radians(st.angle))
    W = st.width_mm * PT_PER_MM if st.width_mm else float("inf")
    more_text = "+%d more"

    def measure(col_w, model_w, fs_species, cols):
        """Column positions and the width the canvas needs."""
        caret = caret_for(model_w)
        bp_pt = model_w / max_cds
        gutter = 2.4 if any(len(c["grid"]) > 1 for c in cells) else 0.0
        label_w = max(tw(l, fs_species, bold=(ri == 0)) for ri, l in enumerate(labels))
        left = PAD + tree_w + (3.0 if tree_w else 0.0) + label_w + 4.0
        goi_w = gutter + max([line_w(ln, bp_pt, caret, ri) for ri in range(n_rows)
                              for ln in cells[ri]["grid"]] + [MINI_W + 4.0]) + 6.0
        goi_horizontal = tw(goi_label, st.fs_species, bold=True) <= goi_w - 2
        xs, ws, x = [], [], left
        for a in cols:
            w = goi_w if a["is_goi"] else col_w
            xs.append(x)
            ws.append(w)
            x += w
        right = x
        for i, a in enumerate(cols):
            if a["is_goi"] and goi_horizontal:
                continue
            fs = st.fs_species if a["is_goi"] else st.fs_col
            right = max(right, xs[i] + ws[i] / 2 + tw(a["label"], fs, bold=a["is_goi"]) * cos_a + 1)
        panel_x0 = panel_w = 0.0
        head = []
        end = x
        if has_panel:
            panel_x0 = x + 8.0 + (5.5 if any(numbered) else 0.0)
            models_w = max([panel_line_w(g, away, bp_pt, caret, ri) for ri in range(n_rows)
                            for g, away in cells[ri]["panel"]]
                           + [tw(more_text % 99, st.fs_small) if any(c["more"] for c in cells) else 0.0])
            # The header must not be what sets the panel width: on two lines when it is wider.
            head = [f"{goi_label} gene models"]
            if tw(head[0], st.fs_col, bold=True) > models_w:
                head = [goi_label, "gene models"]
            panel_w = max([models_w] + [tw(h, st.fs_col, bold=True) for h in head])
            end = panel_x0 + panel_w
            right = max(right, end)
        loc_x0 = 0.0
        if spans:
            loc_x0 = end + 8.0
            right = max(right, loc_x0 + max(tw(s, st.fs_small) for s in spans))
        return {"left": left, "xs": xs, "ws": ws, "grid_x1": x, "need": right + PAD,
                "bp_pt": bp_pt, "caret": caret, "goi_horizontal": goi_horizontal,
                "label_w": label_w, "gutter": gutter, "panel_x0": panel_x0, "panel_w": panel_w,
                "panel_head": head, "loc_x0": loc_x0}

    # ---- fit the width: spacing gives way first, text last, then columns ----------
    col_w, model_w, fs_species = st.col_w, st.goi_model_w, st.fs_species
    cols, dropped = list(anchors), []
    n_flank = lambda: max(1, sum(1 for a in cols if not a["is_goi"]))  # noqa: E731
    m = measure(col_w, model_w, fs_species, cols)
    if m["need"] > W:
        col_w = max(st.col_w_min, col_w - (m["need"] - W) / n_flank())
        m = measure(col_w, model_w, fs_species, cols)
    for _ in range(8):
        if m["need"] <= W or model_w <= st.goi_model_w_min:
            break
        model_w = max(st.goi_model_w_min, model_w - max(2.0, m["need"] - W))
        m = measure(col_w, model_w, fs_species, cols)
    if m["need"] > W and fs_species > st.fs_species_min:
        fs_species = st.fs_species_min
        m = measure(col_w, model_w, fs_species, cols)
    while m["need"] > W and n_flank() > 2:
        gi = next((i for i, a in enumerate(cols) if a["is_goi"]), len(cols) // 2)
        worst = min((a for a in cols if not a["is_goi"]),
                    key=lambda a: (cov[a["key"]], -abs(cols.index(a) - gi)))
        cols.remove(worst)
        dropped.append(worst)
        m = measure(col_w, model_w, fs_species, cols)
    if m["need"] < W:
        col_w = min(st.col_w_max, col_w + (W - m["need"]) / n_flank())
        m = measure(col_w, model_w, fs_species, cols)
    if m["need"] < W:
        model_w = min(st.goi_model_w_max, model_w + (W - m["need"]))
        m = measure(col_w, model_w, fs_species, cols)
        if m["need"] > W:
            model_w = max(st.goi_model_w_min, model_w - (m["need"] - W))
            m = measure(col_w, model_w, fs_species, cols)
    left, xs, ws, grid_x1, bp_pt, caret = (m["left"], m["xs"], m["ws"], m["grid_x1"],
                                           m["bp_pt"], m["caret"])
    gutter = m["gutter"]
    goi_i = next((i for i, a in enumerate(cols) if a["is_goi"]), None)

    def has_axis_for(columns):
        pos = [(a["start"] + a["end"]) / 2.0 for a in columns if a["start"]]
        return st.full and len(pos) >= 2 and max(pos) > min(pos)

    # ---- legend (laid out now: a narrow figure is widened to hold it) --------------
    any_copies = any((not a["is_goi"]) and tmaps[ri].get(a["key"], {}).get("n", 1) > 1
                     for ri in range(1, n_rows) for a in cols)
    # Rows with a flanking ortholog whose coverage has to be stated (see short_coverage).
    short_rows = [ri > 0 and any((not a["is_goi"]) and a["key"] in tmaps[ri]
                                 and short_coverage(tmaps[ri][a["key"]]["best"], st.numbers) is not None
                                 for a in cols) for ri in range(n_rows)]
    kinds = {it[0] for c in cells for ln in c["grid"] for it in ln}
    tiers = {(g.get("confidence") or "").upper() for rc in calls[1:] for g, _, _ in rc}
    tiers |= {(tmaps[ri][a["key"]]["best"].get("confidence") or "").upper()
              for ri in range(1, n_rows) for a in cols if not a["is_goi"] and a["key"] in tmaps[ri]}

    def lay_out_legend(cov_mode):
        groups = build_legend(st, goi_label, kinds, full_drawn, in_column, any(numbered),
                              any(c["more"] for c in cells), any_copies, cov_mode,
                              tree_root is not None, is_species_tree, bool(spans),
                              has_axis_for(cols), tiers)
        legend_w = _legend.natural_width(groups, tw, st.fs_small)
        # A narrow grid is widened for its legend, but not beyond one journal column
        # (89 mm); past that the legend groups continue on a second band.
        width = min(W, max(m["need"], min(legend_w + 2 * PAD,
                                          89.0 * PT_PER_MM if st.width_mm else W)))
        parts, h = _legend.legend_svg(groups, PAD, 0.0, width - 2 * PAD, tw,
                                      fs=st.fs_small, ink=INK, muted=MUTED)
        return groups, width, parts, h

    # ---- heights: a row is as tall as the lines of its GOI cell need ---------------
    mids = [(a["start"] + a["end"]) / 2.0 for a in cols if a["start"]]
    has_axis = has_axis_for(cols)
    title_h = (st.fs_title + st.fs_small + 9.0) if st.full else 0.0
    axis_h = 24.0 if has_axis else 0.0
    header_h = st.fs_species + 5.0
    for i, a in enumerate(cols):
        if a["is_goi"] and m["goi_horizontal"]:
            continue
        fs = st.fs_species if a["is_goi"] else st.fs_col
        header_h = max(header_h, tw(a["label"], fs, bold=a["is_goi"]) * sin_a + fs * 0.75 * cos_a + 3.0)
    header_h = max(header_h, len(m["panel_head"]) * (st.fs_col + 1.0) + 3.0)
    label_y = PAD + title_h + header_h          # baseline band of the column labels
    grid_y0 = label_y + axis_h
    H_max = st.max_height_mm * PT_PER_MM if st.max_height_mm else float("inf")
    n_lines = [max(1, len(c["grid"]), len(c["panel"]) + (1 if c["more"] else 0)) for c in cells]

    def heights(row_h, lane_h, cov_line):
        """(height, room under the arrow line kept for the coverage) of every row."""
        arrow = min(st.arrow_h, row_h - 3.0)
        out = []
        for ri, n in enumerate(n_lines):
            base = max(row_h, n * lane_h + 2.4)
            room = max(0.0, COV_LINE_H - (base - arrow) / 2) if (cov_line and short_rows[ri]) else 0.0
            out.append((base + room, room))
        return out

    def fit_rows(cov_line, legend_h):
        below = 6.0 + (7.0 if (has_goi_col and full_drawn) else 0.0) + legend_h + PAD
        row_h, lane_h = st.row_h, st.lane_h

        def total():
            return grid_y0 + sum(h for h, _ in heights(row_h, lane_h, cov_line)) + below

        while total() > H_max and row_h > st.row_h_min:
            row_h = max(st.row_h_min, row_h - 0.1)
        while total() > H_max and lane_h > st.lane_h_min:
            lane_h = max(st.lane_h_min, lane_h - 0.1)
        return row_h, lane_h, below, total() <= H_max + 0.05

    # Low coverage is written under the arrow, "(44)". That line makes a row taller, so
    # a figure with a height limit that it would break falls back to a "*" on the number.
    cov_mode = "line" if any(short_rows) else ""
    groups, width, legend_parts, legend_h = lay_out_legend(cov_mode)
    row_h, lane_h, below, fits = fit_rows(cov_mode == "line", legend_h)
    if cov_mode == "line" and not fits:
        cov_mode = "star"
        groups, width, legend_parts, legend_h = lay_out_legend(cov_mode)
        row_h, lane_h, below, fits = fit_rows(False, legend_h)
    sized = heights(row_h, lane_h, cov_mode == "line")
    rh = [h for h, _ in sized]
    cov_room = [r for _, r in sized]
    tops = [grid_y0 + sum(rh[:ri]) for ri in range(n_rows)]
    arrow_h = min(st.arrow_h, row_h - 3.0)
    model_h = min(5.2, row_h - 4.4) if max(n_lines) == 1 else min(4.8, lane_h - 1.5)
    grid_y1 = grid_y0 + sum(rh)
    height = grid_y1 + below
    band_x1 = (m["panel_x0"] + m["panel_w"] + 2.0) if has_panel else grid_x1
    if spans:
        band_x1 = width - PAD

    def y_mid(ri):          # the coverage line, where a row has one, lies below this
        return tops[ri] + (rh[ri] - cov_room[ri]) / 2

    P = [f'<rect x="0" y="0" width="{width:.2f}" height="{height:.2f}" fill="#ffffff"/>']

    # ---- title (full version) --------------------------------------------------------
    if st.full:
        n_flank_cols = sum(1 for a in cols if not a["is_goi"])
        P.append(text_svg(left, PAD + st.fs_title, st.fs_title, f"Anchor-grid synteny - GOI: {goi_label}",
                          INK, bold=True, extra=' class="grid-title"'))
        P.append(text_svg(left, PAD + st.fs_title + st.fs_small + 2.5, st.fs_small,
                          f"{n_flank_cols} flanking genes x {n_rows} genomes | orthologs aligned "
                          f"into the columns of the home genes", MUTED, extra=' class="grid-subtitle"'))

    # ---- row bands, GOI column band -------------------------------------------------
    band_x0 = PAD + tree_w + (1.5 if tree_w else 0.0)
    for ri in range(n_rows):
        fill = HOME_BAND if ri == 0 else (STRIPE if ri % 2 == 0 else None)
        if fill:
            P.append(f'<rect x="{band_x0:.2f}" y="{tops[ri]:.2f}" '
                     f'width="{band_x1 - band_x0:.2f}" height="{rh[ri]:.2f}" fill="{fill}"/>')
    if goi_i is not None:
        P.append(f'<rect x="{xs[goi_i]:.2f}" y="{grid_y0:.2f}" width="{ws[goi_i]:.2f}" '
                 f'height="{grid_y1 - grid_y0:.2f}" fill="{goi_colour}" opacity="0.07"/>')

    # ---- home-coordinate axis: where the columns really lie (full version) ----------
    if has_axis:
        pmin, pmax = min(mids), max(mids)
        axis_y = label_y + 8.0
        ax0, ax1 = xs[0] + ws[0] / 2, xs[-1] + ws[-1] / 2

        def real_x(pos):
            return ax0 + (pos - pmin) / (pmax - pmin) * (ax1 - ax0)

        P.append(f'<line x1="{ax0:.2f}" y1="{axis_y:.2f}" x2="{ax1:.2f}" y2="{axis_y:.2f}" '
                 f'stroke="{AXIS}" stroke-width="0.5"/>')
        P.append(text_svg(ax0 - 2.0, axis_y + st.fs_small * 0.35, st.fs_small,
                          f"{pmin / 1e6:.2f} Mb", MUTED, anchor="end"))
        P.append(text_svg(ax1 + 2.0, axis_y + st.fs_small * 0.35, st.fs_small,
                          f"{pmax / 1e6:.2f} Mb", MUTED))
        for i, a in enumerate(cols):
            if not a["start"]:
                continue
            rx = real_x((a["start"] + a["end"]) / 2.0)
            col = goi_colour if a["is_goi"] else AXIS
            P.append(f'<path d="M{rx:.2f},{axis_y:.2f} L{xs[i] + ws[i] / 2:.2f},{grid_y0 - 1.5:.2f}" '
                     f'stroke="{col}" stroke-width="{0.6 if a["is_goi"] else 0.3}" fill="none" '
                     f'opacity="{0.9 if a["is_goi"] else 0.6}"/>')
            P.append(f'<circle cx="{rx:.2f}" cy="{axis_y:.2f}" r="{1.3 if a["is_goi"] else 0.9}" '
                     f'fill="{col}"/>')
        gap_min = max(1.0e6, (pmax - pmin) * 0.18)
        placed = [a for a in cols if a["start"]]
        for a, b in zip(placed, placed[1:]):
            m0, m1 = (a["start"] + a["end"]) / 2.0, (b["start"] + b["end"]) / 2.0
            if m1 - m0 > gap_min:
                gx = (real_x(m0) + real_x(m1)) / 2.0
                for dx in (-1.6, 0.6):
                    P.append(f'<line x1="{gx + dx:.2f}" y1="{axis_y + 1.8:.2f}" x2="{gx + dx + 1.4:.2f}" '
                             f'y2="{axis_y - 1.8:.2f}" stroke="{GAP_INK}" stroke-width="0.5"/>')
                P.append(text_svg(gx, axis_y - 3.0, st.fs_small, f"gap {(m1 - m0) / 1e6:.1f} Mb",
                                  GAP_INK, anchor="middle"))

    # ---- tree --------------------------------------------------------------------
    if tree_root is not None:
        P.append('<g class="grid-tree">' + "".join(
            draw_tree(tree_root, order, PAD, PAD + tree_w, y_mid)) + '</g>')
        if st.full:
            P.append(text_svg(PAD + tree_w / 2, grid_y0 - 3.0, st.fs_small,
                              "Species tree" if is_species_tree else "GOI phylogeny", MUTED,
                              anchor="middle", extra=' font-style="italic"'))

    def draw_cell(ri, x0, ym):
        """The GOI cell of one row: small arrows or the model, weak marks, '//'."""
        is_home = ri == 0
        lines, index = cells[ri]["grid"], cells[ri]["index"]
        n = len(lines)
        ys = [ym + (j - (n - 1) / 2) * lane_h for j in range(n)]
        at_locus = [ys[j] for j, ln in enumerate(lines)
                    if any(not it[2] and it[0] != "count" for it in ln)]
        if len(at_locus) > 1:          # the lines continue one stretch of chromosome
            P.append(f'<line x1="{x0 - 1.6:.2f}" y1="{at_locus[0]:.2f}" x2="{x0 - 1.6:.2f}" '
                     f'y2="{at_locus[-1]:.2f}" stroke="{THREAD}" stroke-width="0.5"/>')
        for j, line in enumerate(lines):
            y, x, pos = ys[j], x0, []
            for it in line:
                w = item_w(it, bp_pt, caret, ri)
                pos.append((x, w))
                x += w + ITEM_GAP
            local = [(px, pw) for (px, pw), it in zip(pos, line)
                     if not it[2] and it[0] != "count"]
            if local:                  # the chromosome line under what is at the locus
                P.append(f'<line x1="{x0 - (1.6 if len(at_locus) > 1 else 0.0):.2f}" y1="{y:.2f}" '
                         f'x2="{local[-1][0] + local[-1][1]:.2f}" y2="{y:.2f}" stroke="{THREAD}" '
                         f'stroke-width="0.4"/>')
            for (px, pw), (kind, g, away) in zip(pos, line):
                if kind == "break":
                    P.append(text_svg(px + pw / 2, y + st.fs_small * 0.35, st.fs_small, "//",
                                      MUTED, anchor="middle"))
                    continue
                if kind == "count":
                    P.append(titled(f"{g} further weak GOI calls in {labels[ri]} (source data)",
                                    text_svg(px, y + st.fs_small * 0.35, st.fs_small, f"+{g}", MUTED),
                                    tips))
                    continue
                fill, dash = goi_fill(g.get("identity"), g.get("confidence"), is_home, goi_colour)
                dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
                where = "elsewhere in the genome" if away else ""
                if kind == "full":
                    body, mw = draw_model(g, px, y, bp_pt, fill, goi_border, dash, model_h,
                                          max_cds, caret)
                    num = number_of(g, ri)
                    if num:
                        body += text_svg(px + mw + 2.0, y + st.fs_num * 0.35, st.fs_num, num, INK,
                                         bold=(g.get("confidence") or "").upper() == "HIGH")
                    P.append(titled(goi_tooltip(g, labels[ri], is_home, where), body, tips))
                elif kind == "mini":
                    strand = g.get("strand", "+")
                    d = ps._svg_arrow_path(px, px + pw, y - MINI_H / 2, MINI_H, strand, rx=0.5)
                    body = (f'<path d="{d}" fill="{fill}" stroke="{goi_border}" '
                            f'stroke-width="0.4"{dash_attr}/>')
                    if numbered[ri]:
                        body += text_svg(px + pw / 2 + (-0.7 if strand != "-" else 0.7),
                                         y + st.fs_small * 0.36, st.fs_small, str(index[id(g)]),
                                         "#ffffff" if ps._is_dark_hex(fill) else "#000000",
                                         anchor="middle")
                    P.append(titled(goi_tooltip(g, labels[ri], is_home, where), body, tips))
                else:                  # weak call: a box for a gene model, a bar for hits
                    mh = 3.2 if ps._is_spliced_model(g) else 1.7
                    body = (f'<rect x="{px:.2f}" y="{y - mh / 2:.2f}" width="{pw:.2f}" '
                            f'height="{mh:.2f}" rx="0.5" '
                            f'fill="{ps._lerp_hex(goi_colour, "#ffffff", 0.7)}" stroke="{goi_border}" '
                            f'stroke-width="0.35" stroke-dasharray="0.6,0.6"/>')
                    note = ("fragment" if ps._is_fragment_goi(g) else "weak call") + (
                        ", " + where if where else "")
                    P.append(titled(goi_tooltip(g, labels[ri], is_home, note), body, tips))

    def draw_panel(ri, x0, ym):
        """The gene models of one row: one line per copy, drawn 5'->3'."""
        is_home = ri == 0
        cell = cells[ri]
        n = len(cell["panel"]) + (1 if cell["more"] else 0)
        ys = [ym + (j - (n - 1) / 2) * lane_h for j in range(n)]
        for j, (g, away) in enumerate(cell["panel"]):
            y, x = ys[j], x0
            if numbered[ri]:
                P.append(text_svg(x0 - 2.4, y + st.fs_small * 0.35, st.fs_small,
                                  str(cell["index"][id(g)]), MUTED, anchor="end"))
            if away:
                P.append(text_svg(x + break_w / 2, y + st.fs_small * 0.35, st.fs_small, "//",
                                  MUTED, anchor="middle"))
                x += break_w + ITEM_GAP
            fill, dash = goi_fill(g.get("identity"), g.get("confidence"), is_home, goi_colour)
            body, mw = draw_model(g, x, y, bp_pt, fill, goi_border, dash, model_h, max_cds, caret,
                                  aligned=True)
            num = number_of(g, ri)
            if num:
                body += text_svg(x + mw + 2.0, y + st.fs_num * 0.35, st.fs_num, num, INK,
                                 bold=(g.get("confidence") or "").upper() == "HIGH")
            P.append(titled(goi_tooltip(g, labels[ri], is_home,
                                        "elsewhere in the genome" if away else ""), body, tips))
        if cell["more"]:
            P.append(text_svg(x0, ys[-1] + st.fs_small * 0.35, st.fs_small,
                              more_text % cell["more"], MUTED, extra=' font-style="italic"'))

    # ---- rows ----------------------------------------------------------------------
    source = []
    for ri, track in enumerate(rows):
        is_home = ri == 0
        tmap = tmaps[ri]
        ym = y_mid(ri)
        present = [xs[i] + ws[i] / 2 for i, a in enumerate(cols)
                   if not a["is_goi"] and (is_home or a["key"] in tmap)]
        thread = []
        if goi_i is not None and cells[ri]["grid"]:
            g0 = xs[goi_i] + 3.0 + gutter
            g1 = g0 + max(line_w(ln, bp_pt, caret, ri) for ln in cells[ri]["grid"])
            before = [p for p in present if p < g0]
            after = [p for p in present if p > g1]
            if before:
                thread.append((min(before), g0 - (1.6 if len(cells[ri]["grid"]) > 1 else 0.0)))
            if after:
                thread.append((g1 + 1.5, max(after)))
        elif len(present) >= 2:
            thread.append((min(present), max(present)))
        for sx0, sx1 in thread:
            P.append(f'<line x1="{sx0:.2f}" y1="{ym:.2f}" x2="{sx1:.2f}" y2="{ym:.2f}" '
                     f'stroke="{THREAD}" stroke-width="0.4"/>')
        italic = ' font-style="italic"' if labels[ri] not in ("Home genome", "Home") else ""
        P.append(text_svg(left - 4.0 - m["label_w"], ym + fs_species * 0.35, fs_species,
                          labels[ri], INK, bold=is_home, extra=italic + ' class="arow-lbl"'))
        if spans and spans[ri]:
            P.append(text_svg(m["loc_x0"], ym + st.fs_small * 0.35, st.fs_small, spans[ri], MUTED))
        for i, a in enumerate(cols):
            if a["is_goi"]:
                draw_cell(ri, xs[i] + 3.0 + gutter, ym)
                if has_panel:
                    draw_panel(ri, m["panel_x0"], ym)
                in_panel = {id(g) for g, _ in cells[ri]["panel"]}
                marked = {id(it[1]) for ln in cells[ri]["grid"] for it in ln if it[0] == "mark"}
                for g, away, copy in calls[ri]:
                    drawn = (("model" if (in_column or id(g) in in_panel) else "arrow") if copy
                             else "mark" if id(g) in marked else "counted")
                    source.append(source_row(ri, track, labels[ri], a, g, n_copies[ri], is_home, True,
                                             drawn, "elsewhere" if away else "locus",
                                             cells[ri]["index"].get(id(g), "")))
                if not calls[ri]:
                    source.append(source_row(ri, track, labels[ri], a, None, 0, False, True, "", "", ""))
                continue
            entry = None if is_home else tmap.get(a["key"])
            if not is_home and not entry:
                source.append(source_row(ri, track, labels[ri], a, None, 0, False, False, "", "", ""))
                continue
            g = None if is_home else entry["best"]
            ident = 0.0 if is_home else (g.get("identity") or 0.0)
            conf = "" if is_home else g.get("confidence", "")
            strand = a["strand"] if is_home else g.get("strand", "+")
            copies = 1 if is_home else entry["n"]
            cov_q = None if is_home else g.get("query_coverage")
            fill, stroke, dash = cell_style(a["colour"], ident, conf, is_home, st.numbers)
            x0, x1 = xs[i] + CELL_GAP / 2, xs[i] + ws[i] - CELL_GAP / 2
            yb = ym - arrow_h / 2
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            body = ""
            if copies > 1:
                d2 = ps._svg_arrow_path(x0 + 1.1, x1 + 1.1, yb - 1.1, arrow_h, strand, rx=0.6)
                body += (f'<path d="{d2}" fill="{ps._lerp_hex(fill, "#ffffff", 0.45)}" '
                         f'stroke="{stroke}" stroke-width="0.35"/>')
            d = ps._svg_arrow_path(x0, x1, yb, arrow_h, strand, rx=0.6)
            body += f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="0.5"{dash_attr}/>'
            if st.numbers and not is_home and ident >= NUMBER_MIN_IDENTITY:
                short = short_coverage(g, True)
                shift = -0.7 if strand != "-" else 0.7      # away from the arrow tip
                body += text_svg((x0 + x1) / 2 + shift, ym + st.fs_small * 0.36, st.fs_small,
                                 f"{ident:.0f}" + ("*" if (short is not None and cov_mode == "star")
                                                   else ""), "#000000", anchor="middle")
                if short is not None and cov_mode == "line":
                    body += text_svg((x0 + x1) / 2, ym + arrow_h / 2 + 0.8 + st.fs_small * 0.72,
                                     st.fs_small, f"({short})", INK, anchor="middle",
                                     extra=' class="cov-low"')
            if is_home:
                tip = f'{a["label"]} (home reference): ' + " · ".join(dict.fromkeys(
                    t for t in (a["product"], a["gene_id"]) if t))
            else:
                own = ps.clean_gene_label(ps._preferred_target_label(g))
                tip = (f'{a["label"]} in {labels[ri]}: identity {ident:.1f}%'
                       + (f', coverage {ps._coverage_pct(cov_q)}%' if cov_q is not None else '')
                       + f', confidence {conf or "-"}'
                       + (', inverted (opposite strand to the home gene)'
                          if strand != a["strand"] else '')
                       + (f', {copies} copies at the locus' if copies > 1 else '')
                       + (f' ({own})' if own else ''))
            P.append(titled(tip, body, tips))
            source.append(source_row(ri, track, labels[ri], a, g, copies, is_home, False,
                                     "arrow", "locus", ""))

    # ---- column labels -------------------------------------------------------------
    for i, a in enumerate(cols):
        cx = xs[i] + ws[i] / 2
        if a["is_goi"]:
            if m["goi_horizontal"]:
                P.append(text_svg(cx, label_y - 3.0, st.fs_species, a["label"], goi_border,
                                  bold=True, anchor="middle", extra=' class="acol-lbl goi"'))
                continue
            fs, fill, bold = st.fs_species, goi_border, True
        else:
            fs, fill, bold = st.fs_col, INK, False
        lx, ly = cx - 1.0, label_y - 2.5
        P.append(text_svg(lx, ly, fs, a["label"], fill, bold=bold,
                          extra=' class="acol-lbl%s" transform="rotate(-%g %.2f %.2f)"' % (
                              " goi" if a["is_goi"] else "", st.angle, lx, ly)))
    for n, line in enumerate(reversed(m["panel_head"])):
        P.append(text_svg(m["panel_x0"], label_y - 3.0 - n * (st.fs_col + 1.0), st.fs_col, line,
                          goi_border, bold=True))

    # ---- scale bar of the GOI models ---------------------------------------------
    y = grid_y1 + 6.0
    if has_goi_col and full_drawn:
        step = next((s for s in SCALE_STEPS if s * bp_pt >= 9.0), SCALE_STEPS[-1])
        bx = m["panel_x0"] if has_panel else xs[goi_i] + 3.0 + gutter
        label = f"{step:,} bp coding sequence"
        if bx + step * bp_pt + 2.0 + tw(label, st.fs_small) > width - PAD:
            bx = max(PAD, width - PAD - step * bp_pt - 2.0 - tw(label, st.fs_small))
        P.append(f'<line x1="{bx:.2f}" y1="{y - 1.5:.2f}" x2="{bx + step * bp_pt:.2f}" '
                 f'y2="{y - 1.5:.2f}" stroke="{INK}" stroke-width="0.6"/>')
        P.append(text_svg(bx + step * bp_pt + 2.0, y + 0.3, st.fs_small, label, MUTED))
        y += 7.0

    # ---- legend ----------------------------------------------------------------------
    P.append(f'<g class="figure-legend" transform="translate(0 {y:.2f})">' + "".join(legend_parts)
             + '</g>')

    sizes = sorted({float(s) for s in re.findall(r'font-size="([0-9.]+)"', "".join(P))})
    info = {
        "numbers": st.numbers, "full": st.full, "goi_form": "column" if in_column else "panel",
        "requested_width_mm": st.width_mm,
        "width_mm": round(width / PT_PER_MM, 1), "height_mm": round(height / PT_PER_MM, 1),
        "rows": n_rows, "flanking_columns": sum(1 for a in cols if not a["is_goi"]),
        "home_genes_without_ortholog": n_unplaced_cols,
        "dropped_columns": [a["label"] for a in dropped + capped],
        "text_pt": sizes, "species_pt": fs_species, "col_pitch_pt": round(col_w, 1),
        "row_pitch_pt": round(row_h, 1), "lane_pitch_pt": round(lane_h, 1),
        "model_width_pt": round(model_w, 1), "bp_per_pt": round(1.0 / bp_pt, 1),
        "max_lines": max(n_lines),
        "fits_height": height <= H_max + 0.05,
        # how low query coverage of a flanking ortholog is shown: "line" = the value in
        # brackets under the arrow, "star" = a * on the identity (no room for the line)
        "low_coverage_shown_as": cov_mode,
        "goi_label": goi_label, "home": labels[0],
        "goi_calls": sum(len(rc) for rc in calls[1:]), "goi_copies": sum(n_copies[1:]),
        "home_copies": n_copies[0], "most_copies_in_a_target": max(n_copies[1:] + [0]),
        "genomes_with_several_copies": sum(1 for n in n_copies if n > 1),
        "copies_not_drawn": sum(c["more"] for c in cells),
    }
    legend_text = build_legend_text(st, info, cols, anchors, is_species_tree, tree_root is not None,
                                    any_copies, full_drawn, kinds, in_column, cov_mode)
    return {"parts": P, "width": width, "height": height, "info": info, "source": source,
            "legend_text": legend_text, "legend_groups": groups,
            "column_keys": [a["key"] for a in cols if not a["is_goi"]]}


def svg_document(result, extra_attrs=""):
    """The stand-alone SVG of a render() result, sized in millimetres."""
    w, h = result["width"], result["height"]
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w / PT_PER_MM:.2f}mm" '
            f'height="{h / PT_PER_MM:.2f}mm" viewBox="0 0 {w:.2f} {h:.2f}" '
            f'font-family="{FONT}"{extra_attrs}>\n' + "\n".join(result["parts"]) + "\n</svg>\n")


def draw_model(gene, x, ym, bp_pt, fill, stroke, dash, h_full, cap_bp=None, caret=2.4,
               aligned=False):
    """One GOI model with its left edge at x; returns (svg, width).

    `aligned` draws the model 5'->3' whatever its strand, so that the first exons
    of all models of a panel line up.
    """
    ps = _PS
    strand = "+" if aligned else gene.get("strand", "+")
    width, exons, introns = model_geometry(gene, bp_pt, cap_bp, caret,
                                           flip=aligned and gene.get("strand", "+") == "-")
    spliced = ps._is_spliced_model(gene)
    h = h_full if spliced else h_full * 0.58
    yb = ym - h / 2
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    parts = []
    for x0, x1 in introns:
        if spliced:
            parts.append(f'<polyline points="{x + x0:.2f},{ym:.2f} {x + (x0 + x1) / 2:.2f},'
                         f'{yb - 1.3:.2f} {x + x1:.2f},{ym:.2f}" fill="none" stroke="{stroke}" '
                         f'stroke-width="0.4" stroke-linejoin="round"/>')
        else:
            parts.append(f'<line x1="{x + x0:.2f}" y1="{ym:.2f}" x2="{x + x1:.2f}" y2="{ym:.2f}" '
                         f'stroke="{stroke}" stroke-width="0.4"/>')
    tip = 0 if strand == "-" else len(exons) - 1
    for k, (x0, x1) in enumerate(exons):
        if k == tip and (x1 - x0) >= 1.6:
            d = ps._svg_arrow_path(x + x0, x + x1, yb, h, strand, rx=0.4)
            parts.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="0.4"{dash_attr}/>')
        else:
            parts.append(f'<rect x="{x + x0:.2f}" y="{yb:.2f}" width="{x1 - x0:.2f}" height="{h:.2f}" '
                         f'rx="0.3" fill="{fill}" stroke="{stroke}" stroke-width="0.4"{dash_attr}/>')
    return "".join(parts), width


def draw_tree(root, order, x0, x1, y_mid):
    """Cladogram with aligned tips. Rows without a leaf hang from the root, dashed."""
    row_of = {ti: ri for ri, ti in enumerate(order)}

    def height(node):
        node.h = 0 if node.track is not None else 1 + max(height(k) for k in node.kids)
        return node.h

    unit = (x1 - x0 - 2.0) / max(1, height(root))

    def place(node):
        node.x = x1 - node.h * unit
        if node.track is not None:
            node.y = y_mid(row_of[node.track])
        else:
            for k in node.kids:
                place(k)
            node.y = (node.kids[0].y + node.kids[-1].y) / 2

    place(root)
    out = [f'<g fill="none" stroke="{TREE_LINE}" stroke-width="0.5" stroke-linecap="square">']

    def draw(node):
        if node.track is not None:
            return
        ys = [k.y for k in node.kids]
        out.append(f'<line x1="{node.x:.2f}" y1="{min(ys):.2f}" x2="{node.x:.2f}" y2="{max(ys):.2f}"/>')
        for k in node.kids:
            out.append(f'<line x1="{node.x:.2f}" y1="{k.y:.2f}" x2="{k.x:.2f}" y2="{k.y:.2f}"/>')
            draw(k)

    draw(root)
    placed = set()

    def collect(node):
        if node.track is not None:
            placed.add(node.track)
        for k in node.kids:
            collect(k)

    collect(root)
    loose = [y_mid(ri) for ri, ti in enumerate(order) if ti not in placed]
    if loose:
        ys = [root.y] + loose
        out.append(f'<line class="tree-unplaced" x1="{root.x:.2f}" y1="{min(ys):.2f}" x2="{root.x:.2f}" '
                   f'y2="{max(ys):.2f}" stroke-dasharray="1.2,1.2"/>')
        for uy in loose:
            out.append(f'<line class="tree-unplaced" x1="{root.x:.2f}" y1="{uy:.2f}" x2="{x1:.2f}" '
                       f'y2="{uy:.2f}" stroke-dasharray="1.2,1.2"/>')
    out.append('</g>')
    return out


# ----------------------------------------------------------------- legend ----

def build_legend(st, goi_label, kinds, full_drawn, in_column, any_numbered, any_more, any_copies,
                 low_cov, has_tree, is_species_tree, has_location, has_axis, tiers):
    """The legend as titled groups (see synvoy_legend)."""
    ps = _PS
    fs = st.fs_small
    base = ps.GENE_PALETTE[0]
    goi_colour, goi_border = ps.GOI_COLOUR, ps.GOI_BORDER

    def arrow(fill, stroke, dash="", stacked=False, aw=9.0, h=4.6, label=""):
        def draw(x, y):
            out = ""
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            if stacked:
                d2 = ps._svg_arrow_path(x + 1.1, x + aw + 1.1, y - h / 2 - 1.1, h, "+", rx=0.6)
                out += (f'<path d="{d2}" fill="{ps._lerp_hex(fill, "#ffffff", 0.45)}" '
                        f'stroke="{stroke}" stroke-width="0.35"/>')
            d = ps._svg_arrow_path(x, x + aw, y - h / 2, h, "+", rx=0.6)
            out += f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="0.5"{dash_attr}/>'
            if label:
                out += text_svg(x + aw / 2 - 0.7, y + fs * 0.36, fs, label,
                                "#ffffff" if ps._is_dark_hex(fill) else "#000000", anchor="middle")
            return out
        return (aw + (1.1 if stacked else 0.0), draw)

    def model(gene, fill):
        bp_pt = 12.0 / cds_bp(gene)
        mw = model_geometry(gene, bp_pt)[0]
        return (mw, lambda x, y: draw_model(gene, x, y, bp_pt, fill, goi_border, "", 4.4)[0])

    def mark():
        return (MARK_W, lambda x, y: (
            f'<rect x="{x:.2f}" y="{y - 1.6:.2f}" width="{MARK_W:.2f}" height="3.2" rx="0.5" '
            f'fill="{ps._lerp_hex(goi_colour, "#ffffff", 0.7)}" stroke="{goi_border}" '
            f'stroke-width="0.35" stroke-dasharray="0.6,0.6"/>'))

    def ramp():
        def draw(x, y):
            return "".join(
                f'<rect x="{x + k * 3.2:.2f}" y="{y - 2.0:.2f}" width="3.0" height="4.0" '
                f'fill="{ps._shade_by_identity(base, ident)}" stroke="{ps._darken_hex(base, 0.65)}" '
                f'stroke-width="0.3"/>' for k, ident in enumerate((95, 65, 35)))
        return (9.4, draw)

    if st.numbers:
        fill_hi = ps._lerp_hex(base, "#ffffff", 0.5)
        fill_md = ps._lerp_hex(base, "#ffffff", 0.66)
        fill_lo = ps._lerp_hex(base, "#ffffff", 0.8)
        stroke = stroke_weak = ps._darken_hex(base, 0.7)
    else:
        fill_hi = ps._shade_by_identity(base, 95)
        fill_md = ps._lerp_hex(base, "#ffffff", 0.42)
        fill_lo = ps._lerp_hex(base, "#ffffff", 0.62)
        stroke, stroke_weak = ps._darken_hex(base, 0.65), base

    groups = []
    if st.full:
        layout = [("Columns", "home genes in genomic order; colour = home gene"),
                  ("Rows", "genomes" + (", ordered by the species tree (NCBI taxonomy)"
                                        if (has_tree and is_species_tree) else
                                        ", ordered by the GOI gene tree" if has_tree else "")),
                  ("Arrow", "ortholog placed in this neighbourhood; points in its coding "
                            "direction, scaffolds oriented like the home genome"),
                  ("Line", "from the first to the last gene placed in a row"),
                  ("Empty cell", "not placed here: no ortholog was placed in this neighbourhood, "
                                 "which is not evidence that the gene is absent; strong hits "
                                 "refused on synteny are listed in synvoy_report.json → "
                                 "rejected_candidates")]
        if has_axis:
            layout.append(("Top axis", "where the home genes lie on the home chromosome"))
        if has_location:
            layout.append(("Right", "scaffold and position of each row's neighbourhood"))
        groups.append(("Layout", layout))

    flank = [(arrow(fill_hi, stroke), "high confidence"),
             (arrow(fill_md, stroke_weak, DASH["MEDIUM"]), "medium confidence"),
             (arrow(fill_lo, stroke_weak, DASH["LOW"]), "low confidence")]
    if "AMBIGUOUS" in tiers:        # drawn, but must not read as a confident call
        flank.append((arrow(ps._lerp_hex(base, "#ffffff", 0.86 if st.numbers else 0.72), stroke_weak,
                            DASH["AMBIGUOUS"]), "ambiguous"))
    numbers = []
    if st.numbers:
        numbers.append(("87", "in an arrow: % identity to the home gene"))
        if low_cov == "line":
            numbers.append(("(44)", "under an arrow: query coverage in %, given when below 80 %"))
        elif low_cov == "star":
            numbers.append(("87*", "the same, query coverage below 80 %"))
    else:
        flank.append((ramp(), "paler = lower identity to the home gene"))
    if any_copies:
        flank.append((arrow(fill_hi, stroke, stacked=True), "more than one copy at the locus"))
    groups.append(("Flanking genes", flank))

    if full_drawn or "mini" in kinds or "mark" in kinds:
        key = {"start": 1, "end": 1400, "strand": "+", "evidence_type": "exon_annotation",
               "exon_coords": [(1, 180), (420, 520), (900, 1400)]}
        goi = []
        if "mini" in kinds:
            goi.append((arrow(goi_colour, goi_border, aw=MINI_W, h=MINI_H,
                              label="1" if any_numbered else ""),
                        "one copy" + ("; numbered along the chromosome" if any_numbered else "")))
        if full_drawn:
            goi.append((model(key, goi_colour),
                        "gene model: coding exons to scale, ∧ = intron"
                        + ("" if in_column else "; drawn 5′→3′")))
            if any(not ps._is_spliced_model(g) for g in full_drawn):
                hit = dict(key, evidence_type="fallback_hit_span", exon_coords=[(1, 300), (900, 1400)])
                goi.append((model(hit, ps._lerp_hex(goi_colour, "#ffffff", 0.36)),
                            "aligned hits without a gene model"))
        if "mark" in kinds:
            goi.append((mark(), "weak call (low or ambiguous confidence)"))
        if "break" in kinds:
            goi.append(("//", "what follows is elsewhere in the genome"))
        if "count" in kinds:
            goi.append(("+N", "further weak calls, not drawn"))
        if any_more:
            goi.append(("+N more", f"copies beyond the {st.max_models} drawn"))
        if tiers - {"HIGH", ""}:
            goi.append((None, "fill and outline: confidence, as for the flanking genes"))
        groups.append((f"GOI: {goi_label}", goi))
        if full_drawn:
            numbers.append(("78", "after a GOI model: % identity to the query; bold = high "
                                  "confidence"))
            numbers.append(("56 (44)", "identity (query coverage, when below 80 %)"))
    if numbers:
        groups.append(("Numbers", numbers))
    return groups


def build_legend_text(st, info, cols, anchors, is_species_tree, has_tree, any_copies, full_drawn,
                      kinds, in_column, low_cov):
    """The sentences the print figure does not carry, as a draft figure legend."""
    flank = [a for a in cols if not a["is_goi"]]
    placed = [a for a in anchors if a["start"]]
    where = ""
    if placed:
        lo, hi = min(a["start"] for a in placed), max(a["end"] for a in placed)
        chrom = next((a["chrom"] for a in placed if a.get("chrom")), "")
        where = f" ({chrom + ': ' if chrom else ''}{lo / 1e6:.2f}–{hi / 1e6:.2f} Mb)"
    goi = info["goi_label"]
    s = [f"{goi} and its neighbouring genes in {info['rows']} genomes."]
    s.append(f"Columns are the {len(flank)} genes flanking {goi} in {info['home']}{where}, "
             f"in genomic order; rows are genomes"
             + (", ordered by the NCBI-taxonomy species tree shown on the left (branch lengths "
                "not to scale)." if (has_tree and is_species_tree) else
                ", ordered by the gene tree of the GOI shown on the left." if has_tree else "."))
    notes = []
    if info["home_genes_without_ortholog"]:
        notes.append(f"{info['home_genes_without_ortholog']} home gene(s) placed in no other "
                     f"genome are not shown")
    if info["dropped_columns"]:
        notes.append(f"{len(info['dropped_columns'])} column(s) with the fewest orthologs were "
                     f"left out to fit the figure width ({', '.join(info['dropped_columns'])})")
    if notes:
        text = "; ".join(notes)
        s.append(text[0].upper() + text[1:] + ".")
    s.append("An arrow is the ortholog placed in that genome's neighbourhood and points in its "
             "coding direction, with each scaffold oriented like the home genome.")
    if st.numbers:
        s.append("Colour: the home gene the ortholog belongs to. Number: protein identity to the "
                 "home gene in %" + (
                     "; a number in brackets under an arrow is the query coverage in %, given "
                     "when it is below 80 %." if low_cov == "line" else
                     ", with * when the query coverage is below 80 %." if low_cov == "star"
                     else "."))
    else:
        s.append("Colour: the home gene the ortholog belongs to; paler fill = lower protein "
                 "identity to it.")
    s.append("Outline: solid = high, dashed = medium, dotted = low confidence.")
    if any_copies:
        s.append("A second arrow behind the first marks more than one copy at the locus.")
    s.append("The grey line joins the first and the last gene placed in a row; a gap on that "
             "line is a gene that was not placed in this neighbourhood, which is not evidence "
             "that the genome lacks it.")
    if full_drawn or "mini" in kinds:
        copies = ("A copy is the best GOI call of a genome and every further call of at least "
                  "medium confidence")
        if in_column:
            s.append(f"{goi} column: the gene model found in each genome, in its coding "
                     "direction.")
        else:
            s.append(f"{goi} column: one small arrow per copy, in genomic order and pointing in "
                     "its coding direction, numbered along the chromosome where a genome has "
                     "several. The panel on the right draws the exon structure of each copy on a "
                     "line of its own, 5′→3′, under the same number. " + copies + ".")
            if info["copies_not_drawn"]:
                s.append(f"At most {st.max_models} gene models are drawn for one genome; "
                         f"{info['copies_not_drawn']} further copies are in the source data.")
        if "mark" in kinds:
            s.append("A small pale mark is a weak call (low or ambiguous confidence): a box for a "
                     "gene model, a bar for aligned hits."
                     + (f" At most {st.max_marks} are drawn in a cell; +N counts the rest."
                        if "count" in kinds else ""))
        if "break" in kinds:
            s.append("Calls after // lie on another scaffold or farther from the locus than the "
                     "neighbourhood is long.")
        if full_drawn:
            s.append("Coding exons are to scale (bar), ∧ = intron (not to scale); thin boxes are "
                     "aligned hit segments without a gene model (drawn no longer than the longest "
                     "model). The numbers are the protein identities to the query in %, with the "
                     "query coverage in brackets when it is below 80 %.")
    if any(a["label"].endswith("*") for a in cols):
        s.append("Labels ending in * are abbreviations of the annotated product name, not "
                 "official gene symbols (gene-names table).")
    s.append("Values of every cell and the genomic location of every gene: source data.")
    return " ".join(s) + "\n"


SOURCE_COLUMNS = ("row", "species", "genome_id", "column", "home_gene_id", "state", "drawn_as",
                  "copy_number", "position", "identity_pct", "query_coverage_pct", "confidence",
                  "copies", "strand_drawn", "scaffold", "start", "end", "exons", "evidence")


def source_row(ri, track, label, anchor, gene, n, is_home, is_goi_col, drawn_as, position, index):
    row = dict.fromkeys(SOURCE_COLUMNS, "")
    row.update(row=ri + 1, species=label, genome_id=track.get("genome_id", ""),
               column=anchor["label"], home_gene_id=anchor.get("gene_id", ""),
               drawn_as=drawn_as, position=position, copy_number=index)
    if is_home and not is_goi_col:
        row.update(state="home", strand_drawn=anchor["strand"], scaffold=anchor.get("chrom", ""),
                   start=anchor["start"], end=anchor["end"], copies=1)
        return row
    if gene is None:
        row["state"] = "not placed"
        return row
    cov = gene.get("query_coverage")
    row.update(state="home" if is_home else "placed",
               identity_pct="" if is_home else f"{gene.get('identity') or 0:.1f}",
               query_coverage_pct="" if (is_home or cov is None) else str(_PS._coverage_pct(cov)),
               confidence=gene.get("confidence", "") or "", copies=n,
               strand_drawn=gene.get("strand", ""), scaffold=gene.get("chrom", ""),
               start=gene.get("start", ""), end=gene.get("end", ""),
               exons=len(gene.get("exon_coords") or []) if is_goi_col else "",
               evidence=gene.get("evidence_type", "") or gene.get("goi_model_source", "") or "")
    return row
