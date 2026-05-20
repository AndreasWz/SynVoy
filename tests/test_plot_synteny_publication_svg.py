"""Contract test for the publication-SVG export.

The pub SVG used to be a narrow Nature-column-style layout produced by a
separate renderer; that renderer was removed when we unified the look with
the interactive HTML. The pub SVG is now exactly the HTML view's inline
SVG, but with every home-genome gene labelled on the canvas.

This test exercises that contract: the pub SVG mirrors the HTML render,
contains home gene labels, and stays a valid SVG (no <html> wrapper, no
embedded <script>).
"""

import os
import sys
from types import SimpleNamespace
from xml.etree import ElementTree as ET

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny


def _gene(name, start, end, home_gene_id=None, identity=75.0, strand="+"):
    return {
        "chrom": "chr1",
        "start": start,
        "end": end,
        "start_plot": start,
        "end_plot": end,
        "name": name,
        "home_gene_id": home_gene_id or name,
        "strand": strand,
        "identity": identity,
        "exon_coords": [],
    }


def _make_args():
    """Minimal arg namespace; the renderer reads these directly."""
    return SimpleNamespace(
        plot_width=0,
        plot_height=0,
        scale_bar_len=10000,
        # Legacy flags kept for Nextflow compatibility but ignored.
        pub_width=89,
        pub_palette="okabe_ito",
    )


def _make_tracks():
    return [
        {
            "label": "Apis mellifera (GCF_000001)",
            "genes": [
                _gene("geneA", 100, 180, identity=92.0),
                _gene("GOI_test", 240, 320, identity=99.0),
                _gene("hidden1", 360, 420, identity=45.0),
                _gene("hidden2", 460, 520, identity=45.0),
                _gene("geneB", 580, 660, identity=88.0, strand="-"),
            ],
            "is_home": True,
            "genome_id": "home",
            "goi_status": "resolved",
            "offset": 0,
            "breaks": [],
            "minus_strand_row": None,
        },
        {
            "label": "Bombus terrestris (GCF_000002)",
            "genes": [
                _gene("orthA", 120, 200, home_gene_id="geneA", identity=93.0),
                _gene("GOI_match", 250, 330, home_gene_id="GOI_test", identity=84.0),
            ],
            "is_home": False,
            "genome_id": "target1",
            "goi_status": "resolved",
            "offset": 0,
            "breaks": [],
            "minus_strand_row": None,
        },
    ]


def test_render_publication_svg_is_valid_xml_without_html_wrapper():
    svg = plot_synteny.render_publication_svg(
        _make_tracks(),
        {"geneA": "#0072B2", "geneB": "#009E73", "GOI_test": "#D55E00"},
        {"home": "#E64B35", "target1": "#0072B2"},
        {},
        _make_args(),
        ["test"], 0, 0, 1,
    )

    # Pub SVG is *not* an HTML document: no <html>, <body>, or <script>.
    lower = svg.lower()
    assert "<html" not in lower
    assert "<body" not in lower
    assert "<script" not in lower

    # Still parses as XML, with svg as the root element.
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")


def test_pub_svg_labels_every_home_gene_but_html_does_not():
    tracks = _make_tracks()
    gene_colours = {"geneA": "#0072B2", "geneB": "#009E73", "GOI_test": "#D55E00"}
    goi_colours  = {"home": "#E64B35", "target1": "#0072B2"}
    args = _make_args()

    pub_svg = plot_synteny.render_publication_svg(
        tracks, gene_colours, goi_colours, {}, args, ["test"], 0, 0, 1,
    )
    html = plot_synteny.render_synteny_html(
        tracks, gene_colours, goi_colours, {}, args, ["test"], 0, 0, 1,
        force_home_labels=False,
    )

    # The home track has 5 genes (one of which is the GOI). The pub SVG
    # should label all 5 — GOI as ".gene-label goi", flanking as plain
    # ".gene-label". The interactive HTML labels only the GOI.
    pub_goi    = pub_svg.count('class="gene-label goi')
    pub_flank  = pub_svg.count('class="gene-label track-item"')
    html_goi   = html.count('class="gene-label goi')
    html_flank = html.count('class="gene-label track-item"')

    # GOI label present in both. (Two GOI gene-groups in the fixture — one
    # home, one target.)
    assert pub_goi >= 2
    assert html_goi >= 2
    # Home flanking labels only in pub SVG.
    assert pub_flank >= 4   # geneA, hidden1, hidden2, geneB
    assert html_flank == 0


def test_pub_svg_carries_shared_style_defs():
    """The pub SVG mirrors the interactive HTML's look, including the
    shadow + gloss defs and the CSS-vars body. These exist in the inline
    SVG (geneGloss/geneShadow gradients) and in the embedded CSS
    (drop-shadow on .gene-group). The test guards against silent
    regressions where the shared styling drifts apart."""
    svg = plot_synteny.render_publication_svg(
        _make_tracks(),
        {"geneA": "#0072B2", "geneB": "#009E73", "GOI_test": "#D55E00"},
        {"home": "#E64B35", "target1": "#0072B2"},
        {},
        _make_args(),
        ["test"], 0, 0, 1,
    )
    assert 'id="geneGloss"' in svg
    assert 'id="geneShadow"' in svg
    assert "filter: drop-shadow" in svg
