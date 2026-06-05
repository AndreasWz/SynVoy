#!/usr/bin/env python3
"""Tests for the anchor-grid figure (matrix × synteny hybrid) and the
CSS-variable resolver that fixes black tracks in exported SVGs.

The anchor grid aligns every recovered orthologue into shared home-gene
columns (one row per species), so the synteny reads from the vertical
alignment rather than a tangle of ribbons. These tests pin the contract:
the GOI gets its own emphasised column, identity numbers are drawn, strand
inversions are flagged, and home genes recovered in no target are dropped.
"""
import os
import sys
from types import SimpleNamespace
from xml.etree import ElementTree as ET

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402


def _gene(name, start, end, home_gene_id=None, identity=75.0, strand="+",
          confidence="HIGH", role=None):
    g = {
        "chrom": "chr1", "start": start, "end": end,
        "start_plot": start, "end_plot": end,
        "name": name, "home_gene_id": home_gene_id or name,
        "strand": strand, "identity": identity, "exon_coords": [],
        "confidence": confidence,
    }
    if role:
        g["role"] = role
    return g


def _args():
    return SimpleNamespace()


def _tracks():
    """Home (5 anchors incl. GOI) + two targets with partial recovery."""
    home = {
        "label": "Home genome (chr1)", "is_home": True, "genome_id": "home",
        "offset": 0, "breaks": [],
        "genes": [
            _gene("gene-AAA", 1000, 2000, strand="+"),
            _gene("gene-BBB", 3000, 4000, strand="-"),
            _gene("GOI_DCN", 5000, 6000, strand="-"),
            _gene("gene-CCC", 7000, 8000, strand="+"),
            _gene("gene-LONELY", 9000, 9500, strand="+"),  # recovered nowhere
        ],
    }
    cow = {
        "label": "Bos taurus (cow)", "is_home": False, "genome_id": "cow",
        "offset": 0, "breaks": [],
        "genes": [
            _gene("orthA", 100, 200, home_gene_id="gene-AAA", identity=99.0),
            # decorin recovered HIGH, same (minus) strand as home
            _gene("GOI_DCN|cow_b0_l1_exon_ann", 500, 600,
                  home_gene_id="GOI_DCN", identity=90.0, role="goi"),
            _gene("orthC", 700, 800, home_gene_id="gene-CCC", identity=88.0,
                  strand="-"),  # inverted vs home (+)
        ],
    }
    mouse = {
        "label": "Mus musculus (mouse)", "is_home": False, "genome_id": "mouse",
        "offset": 0, "breaks": [],
        "genes": [
            _gene("orthA2", 110, 210, home_gene_id="gene-AAA", identity=97.0),
            _gene("GOI_DCN|mouse_b0_l1_exon_ann", 510, 610,
                  home_gene_id="GOI_DCN", identity=80.0, role="goi"),
            # tandem second copy of the GOI
            _gene("GOI_DCN|mouse_b0_l2_exon_ann", 520, 620,
                  home_gene_id="GOI_DCN", identity=78.0, role="goi"),
        ],
    }
    return [home, cow, mouse]


def _render():
    gene_colours = {"gene-AAA": "#0072B2", "gene-BBB": "#009E73",
                    "gene-CCC": "#56B4E9", "gene-LONELY": "#666666"}
    goi_colours = {"home": "#E64B35", "cow": "#0072B2", "mouse": "#CC79A7"}
    return ps.render_anchor_grid(_tracks(), gene_colours, goi_colours, {}, _args())


# ---- structural / validity ------------------------------------------------

def test_anchor_grid_is_valid_html_with_one_svg():
    html = _render()
    low = html.lower()
    assert "<!doctype html>" in low
    assert html.count("<svg") == 1 and html.count("</svg>") == 1
    # The embedded SVG must parse as XML on its own.
    import re
    m = re.search(r"(<svg.*?</svg>)", html, re.DOTALL)
    assert m, "no <svg> in anchor-grid HTML"
    root = ET.fromstring(m.group(1))
    assert root.tag.endswith("svg")


def test_goi_column_label_and_guide_present():
    html = _render()
    # GOI symbol derived from 'GOI_DCN' -> 'DCN', used as the column label/title.
    assert ">DCN<" in html
    assert "Anchor-grid synteny" in html
    # Red GOI guide band drawn behind the columns.
    assert ps.GOI_COLOUR in html


# ---- alignment semantics --------------------------------------------------

def test_orthologs_carry_identity_numbers():
    html = _render()
    # cow decorin 90 %, mouse decorin 80 %, flanking 99/97/88 are all drawn.
    for ident in (">90<", ">80<", ">99<", ">88<"):
        assert ident in html, f"missing identity label {ident}"


def test_zero_coverage_home_gene_is_dropped():
    html = _render()
    # gene-LONELY is a home anchor recovered in no target -> dropped from the
    # alignment view; the others (incl. GOI) remain.
    assert ">LONELY<" not in html
    assert ">AAA<" in html and ">CCC<" in html


def test_inversion_is_flagged():
    # cow orthC is on '-' while home gene-CCC is '+': the inversion is conveyed
    # by drawing the arrow in the TARGET strand (so it points opposite the home
    # gene) and is annotated in the cell tooltip.
    html = _render()
    assert ", inverted" in html


def test_tandem_copy_badge():
    html = _render()
    # mouse has two GOI copies -> a '×2' badge.
    assert "×2" in html


def test_goi_display_label_prefers_alpha_symbol():
    # 'GOI_DCN' -> 'DCN'; coordinate-style 'GOI_NC_000012.12_91140483' rejected.
    tracks = [{"genes": [
        {"home_gene_id": "GOI_DCN", "name": "x"},
        {"home_gene_id": "GOI_NC_000012.12_91140483", "name": "y"},
    ]}]
    assert ps._goi_display_label(tracks) == "DCN"
    assert ps._goi_display_label([{"genes": []}]) == "GOI"


# ---- CSS-variable resolver (fixes black tracks in exported SVG) ------------

def test_resolve_css_vars_substitutes_root_literals():
    css = ":root { --track-bg: #f3f5f8; --goi-color: #e31a1c; }\n" \
          ".track-bg { fill: var(--track-bg); }\n" \
          ".goi { fill: var(--goi-color); }"
    out = ps._resolve_css_vars(css)
    assert "var(--" not in out
    assert "fill: #f3f5f8" in out
    assert "fill: #e31a1c" in out


def test_resolve_css_vars_uses_fallback_for_unknown_var():
    css = ":root { --a: #111111; }\n.x { fill: var(--missing, #abcdef); }"
    out = ps._resolve_css_vars(css)
    assert "#abcdef" in out and "var(--" not in out


def test_resolve_css_vars_noop_without_root():
    css = ".x { fill: var(--whatever); }"
    # No :root block -> nothing to resolve; returned unchanged.
    assert ps._resolve_css_vars(css) == css


def test_export_inline_svg_has_no_unresolved_vars(tmp_path):
    # End-to-end: the ribbon HTML's exported static SVG must not keep any
    # var(--…) (which renders black in cairosvg/Illustrator).
    html = tmp_path / "p_synteny_plot.html"
    html.write_text(
        "<html><head><style>:root{--track-bg:#f3f5f8;}"
        ".track-bg{fill:var(--track-bg);}</style></head><body>"
        '<svg width="10" height="10"><rect class="track-bg" '
        'width="10" height="10"/></svg></body></html>'
    )
    svg = tmp_path / "p_synteny_plot_view.svg"
    ps._export_html_inline_svg(str(html), str(svg))
    text = svg.read_text()
    assert "var(--" not in text
    assert "#f3f5f8" in text


def test_shade_by_identity_fades_to_white_for_divergent():
    base = "#0072B2"
    full = ps._shade_by_identity(base, 100.0)
    faded = ps._shade_by_identity(base, 20.0)
    assert full.lower() == base.lower()       # 100 % keeps the colour
    assert faded.lower() != base.lower()       # divergent fades
    # faded is lighter (closer to white) than full on every channel
    fr = tuple(int(faded.lstrip("#")[i:i+2], 16) for i in (0, 2, 4))
    br = tuple(int(base.lstrip("#")[i:i+2], 16) for i in (0, 2, 4))
    assert all(f >= b for f, b in zip(fr, br))
