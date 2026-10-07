#!/usr/bin/env python3
"""Regression tests for the 2026-09 figure pass (docs/VIZ_REVIEW_2026-09.md).

- the ribbon plot uses at most 3 lanes per genome, placing the GOI and
  ribbon-bearing genes first, and widening never adds a lane;
- the home GOI is the annotated gene (real strand, one transcript's CDS), not
  one '+' placeholder row per query hit;
- ModelStatus=fragment GOI models are hidden by default and cannot erase the
  flanking gene they overlap;
- GOI models in the anchor grid are drawn from their CDS blocks on one shared
  scale, and the location column prints genomic (not plot) coordinates.
"""
import os
import re
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402


def _g(name, start, end, **kw):
    g = {"chrom": "chr1", "start": start, "end": end, "start_plot": start,
         "end_plot": end, "name": name, "strand": "+", "exon_coords": []}
    g.update(kw)
    return g


# ---------------------------------------------------------------- lanes ----

def test_lanes_never_exceed_three_and_nothing_is_dropped():
    # Six genes stacked on the same interval would need six lanes.
    genes = [_g(f"g{i}", 1000, 5000) for i in range(6)]
    ps._assign_sub_tracks(genes, 0)
    assert {g["_sub_track"] for g in genes} <= {0, 1, 2}
    assert sum(1 for g in genes if g.get("_lane_overflow")) == 3
    assert len(genes) == 6


def test_priority_genes_get_lanes_before_others():
    goi = _g("GOI", 3000, 3500)
    others = [_g(f"o{i}", 1000, 6000) for i in range(3)]
    ps._assign_sub_tracks(others + [goi], 0, is_priority=lambda g: g["name"] == "GOI")
    assert not goi.get("_lane_overflow")
    assert goi["_sub_track"] == 0


def test_non_overlapping_genes_share_one_lane():
    genes = [_g("a", 0, 1000), _g("b", 2000, 3000), _g("c", 4000, 5000)]
    ps._assign_sub_tracks(genes, 0)
    assert {g["_sub_track"] for g in genes} == {0}


def test_widening_after_lane_assignment_keeps_lane_neighbours_apart():
    # Sparse track: widening would make neighbours overlap without the lane clamp.
    genes = [_g("a", 0, 1000), _g("b", 2000, 3000), _g("c", 4000, 5000)]
    track = {"genes": genes, "offset": 0, "breaks": []}
    ps._assign_sub_tracks(genes, 0)
    factor = ps._widen_sparse_plot([track], target_coverage=2.0, max_factor=4.0)
    assert factor > 1.0
    spans = sorted((g["start_plot"], g["end_plot"]) for g in genes)
    for (s0, e0), (s1, e1) in zip(spans, spans[1:]):
        assert e0 < s1


# ------------------------------------------------------ home GOI models ----

HOME_GFF = """\
##gff-version 3
chr1\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene-CONT;Name=CONT
chr1\t.\tmRNA\t1000\t9000\t.\t+\t.\tID=rna-CONT;Parent=gene-CONT
chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds-CONT;Parent=rna-CONT
chr1\t.\tCDS\t8800\t9000\t.\t+\t0\tID=cds-CONT;Parent=rna-CONT
chr1\t.\tgene\t3000\t4600\t.\t-\t.\tID=gene-TOX;Name=Tox
chr1\t.\tmRNA\t3000\t4600\t.\t-\t.\tID=rna-TOX1;Parent=gene-TOX
chr1\t.\tCDS\t3000\t3100\t.\t-\t0\tID=cds-TOX1;Parent=rna-TOX1
chr1\t.\tCDS\t4500\t4600\t.\t-\t0\tID=cds-TOX1;Parent=rna-TOX1
chr1\t.\tmRNA\t3000\t4600\t.\t-\t.\tID=rna-TOX2;Parent=gene-TOX
chr1\t.\tCDS\t3000\t3100\t.\t-\t0\tID=cds-TOX2;Parent=rna-TOX2
chr1\t.\tgene\t12000\t13000\t.\t+\t.\tID=gene-PARA;Name=Para
chr1\t.\tmRNA\t12000\t13000\t.\t+\t.\tID=rna-PARA;Parent=gene-PARA
chr1\t.\tCDS\t12000\t12300\t.\t+\t0\tID=cds-PARA;Parent=rna-PARA
"""


def _home_setup(tmp_path):
    gff = tmp_path / "home.gff"
    gff.write_text(HOME_GFF)
    hits = [  # BED half-open, as in locus_*.bed
        {"chrom": "chr1", "start": 3010, "end": 3100, "strand": "+"},
        {"chrom": "chr1", "start": 4510, "end": 4600, "strand": "+"},
        {"chrom": "chr1", "start": 12050, "end": 12250, "strand": "+"},
    ]
    home = [_g("gene-FLANK", 20000, 21000)] + [
        _g(f"GOI_chr1_{q['start']}", q["start"], q["end"]) for q in hits]
    return str(gff), hits, home


def test_home_goi_resolves_to_annotated_genes(tmp_path):
    gff, hits, home = _home_setup(tmp_path)
    names = ps._resolve_home_goi_models(home, hits, gff)
    assert names == ["gene-PARA", "gene-TOX"]
    by = {g["name"]: g for g in home}
    # Placeholders are gone; the flanking gene stays; the container gene whose
    # intron holds the hits (no CDS overlap) is not chosen.
    assert not any(n.startswith("GOI_") for n in by)
    assert "gene-FLANK" in by and "gene-CONT" not in by
    tox = by["gene-TOX"]
    assert tox["strand"] == "-"
    # Transcript covering both hits wins over the one-exon isoform.
    assert tox["exon_coords"] == [(3000, 3100), (4500, 4600)]
    assert tox["goi_model_source"] == "rna-TOX1"


def test_home_goi_without_gff_leaves_placeholders(tmp_path):
    _, hits, home = _home_setup(tmp_path)
    before = [g["name"] for g in home]
    assert ps._resolve_home_goi_models(home, hits, "NO_GFF") == []
    assert [g["name"] for g in home] == before


# ------------------------------------------------------------ fragments ----

TARGET_GFF = """\
##gff-version 3
t1\tSynVoy\tmRNA\t1000\t3000\t.\t+\t.\tID=gene-FL|t_b0_fl1_flank_ann;Name=gene-FL;SynVoy_Parent=gene-FL;SynVoyRole=flanking;Identity=90;Confidence=HIGH;ModelStatus=complete;Exons=1
t1\tSynVoy\tmRNA\t1500\t1600\t.\t+\t.\tID=GOI_X|t_b0_l1_extra1;Name=GOI_X;SynVoy_Parent=GOI_X;SynVoyRole=goi;Identity=60;EvidenceType=rescued_exon;Confidence=LOW;ModelStatus=fragment;QueryCoverage=0.04;Exons=1
t1\tSynVoy\tmRNA\t5000\t6000\t.\t-\t.\tID=GOI_X|t_b0_l1_exon_ann;Name=GOI_X;SynVoy_Parent=GOI_X;SynVoyRole=goi;Identity=70;EvidenceType=exon_annotation;Confidence=HIGH;ModelStatus=complete;QueryCoverage=0.95;Exons=2
t1\tSynVoy\tCDS\t5000\t5200\t.\t-\t0\tParent=GOI_X|t_b0_l1_exon_ann
t1\tSynVoy\tCDS\t5800\t6000\t.\t-\t0\tParent=GOI_X|t_b0_l1_exon_ann
"""


@pytest.fixture
def target_gff(tmp_path):
    p = tmp_path / "t.fna.gff"
    p.write_text(TARGET_GFF)
    yield str(p)
    ps._SHOW_GOI_FRAGMENTS = False


def test_fragments_hidden_by_default_and_flanking_kept(target_gff):
    ps._SHOW_GOI_FRAGMENTS = False
    names = sorted(g["name"] for g in ps.parse_target_gff(target_gff))
    # Without the filter, dedup (GOI always wins) would erase gene-FL.
    assert names == ["GOI_X", "gene-FL"]
    goi = [g for g in ps.parse_target_gff(target_gff) if g["name"] == "GOI_X"]
    assert goi[0]["exon_coords"] == [(5000, 5200), (5800, 6000)]


def test_fragments_shown_in_variant(target_gff):
    ps._SHOW_GOI_FRAGMENTS = True
    genes = ps.parse_target_gff(target_gff)
    assert any(ps._is_fragment_goi(g) for g in genes)


def test_variant_path_suffix():
    assert ps._variant_path("x_anchor_grid.svg", "_with_fragments") == \
        "x_anchor_grid_with_fragments.svg"
    assert ps._variant_path("x_anchor_grid.svg", "") == "x_anchor_grid.svg"


# ------------------------------------------------------ GOI cell models ----

def _model(strand="-"):
    return {"start": 1000, "end": 3000, "strand": strand,
            "evidence_type": "exon_annotation",
            "exon_coords": [(1000, 1099), (2001, 2300), (2901, 3000)]}


def test_cds_geometry_is_to_scale_with_fixed_introns():
    w, exons, introns, pointing = ps._goi_model_geometry(_model(), "cds", 0.1)
    assert [round(x1 - x0, 6) for x0, x1 in exons] == [10.0, 30.0, 10.0]
    assert all(round(x1 - x0, 6) == ps.GOI_INTRON_PX for x0, x1 in introns)
    assert pointing == "-"
    assert w == pytest.approx(50 + 2 * ps.GOI_INTRON_PX)


def test_cds_aligned_draws_minus_strand_five_prime_first():
    _, exons, _, pointing = ps._goi_model_geometry(_model("-"), "cds_aligned", 0.1)
    # 5' exon of a '-' gene is the genomically LAST block (100 bp) ... all
    # blocks are symmetric here, so check the middle one stays in the middle
    # and the tip points right.
    assert [round(x1 - x0, 6) for x0, x1 in exons] == [10.0, 30.0, 10.0]
    assert pointing == "+"
    asym = dict(_model("-"), exon_coords=[(1000, 1099), (2001, 2300)])
    _, ex2, _, _ = ps._goi_model_geometry(asym, "cds_aligned", 0.1)
    assert round(ex2[0][1] - ex2[0][0], 6) == 30.0


def test_genomic_geometry_keeps_intron_lengths():
    w, exons, introns, _ = ps._goi_model_geometry(_model(), "genomic", 0.1)
    assert w == pytest.approx(200.1)
    assert introns[0][1] - introns[0][0] == pytest.approx(90.1)


def test_column_scale_is_shared_and_longest_model_gets_budget():
    short = dict(_model(), exon_coords=[(1000, 1099)])
    long_ = _model()
    bp_px, col_w = ps._goi_column_scale([[short], [long_]], "cds")
    assert bp_px == pytest.approx(ps.GOI_MODEL_MAX_W / 500)
    assert col_w <= ps.GOI_COL_MAX_W


def test_hit_chains_are_not_drawn_as_gene_models():
    hit = dict(_model(), evidence_type="fallback_hit_span")
    svg_model, _ = ps._goi_model_svg(_model(), 0, 20, "cds", 0.1, "#f00", "#800", "")
    svg_hit, _ = ps._goi_model_svg(hit, 0, 20, "cds", 0.1, "#f00", "#800", "")
    assert "<polyline" in svg_model          # intron carets
    assert "<polyline" not in svg_hit        # plain connector lines
    h_model = float(re.search(r'height="([\d.]+)"', svg_model).group(1))
    h_hit = float(re.search(r'height="([\d.]+)"', svg_hit).group(1))
    assert h_hit < h_model


def test_ambiguous_goi_is_not_styled_as_high():
    high, _ = ps._goi_copy_fill(80, "HIGH", False)
    amb, dash = ps._goi_copy_fill(80, "AMBIGUOUS", False)
    assert amb != high and "dasharray" in dash


def test_grid_location_column_uses_genomic_coordinates():
    home = {"label": "Home", "species": "Apis mellifera", "is_home": True,
            "genome_id": "home", "offset": 0, "breaks": [],
            "genes": [_g("gene-A", 1000, 2000), _g("GOI_T", 3000, 3500)]}
    # A reflected target: plot coordinates negative, genomic coordinates at 5 Mb.
    tgt_genes = [
        _g("a", 5_000_000, 5_001_000, home_gene_id="gene-A", identity=90.0,
           confidence="HIGH", start_plot=-5_001_000, end_plot=-5_000_000),
        _g("GOI_T|x", 5_010_000, 5_010_500, home_gene_id="GOI_T", role="goi",
           identity=80.0, confidence="HIGH", start_plot=-5_010_500, end_plot=-5_010_000),
    ]
    tgt = {"label": "T", "species": "Bombus terrestris", "is_home": False,
           "genome_id": "bter", "offset": 0, "breaks": [], "genes": tgt_genes}
    html = ps.render_anchor_grid([home, tgt], {"gene-A": "#4e79a7"}, {}, {},
                                 SimpleNamespace())
    assert "chr1: 5.00-5.01 Mb" in html
    assert ": -" not in html
    # Species label: binomial only, italic.
    assert re.search(r'font-style="italic"[^>]*>Bombus terrestris<', html)


def test_text_width_matches_arial_metrics():
    assert ps.text_width("Apis", 10) == pytest.approx((667 + 556 + 222 + 500) / 100)
    assert ps.text_width("Apis", 10, bold=True) > ps.text_width("Apis", 10)


# ------------------------------------------------------ orientation ----

def _orient_fixture():
    """Home: 4 flanking genes (+,-,+,-) around a GOI. Target: the same genes on
    TWO scaffolds, each read the wrong way round (reversed order, flipped strands),
    laid out one after the other as compress_track_coordinates does."""
    home_genes = [_g("gA", 0, 100, strand="+"), _g("gB", 200, 300, strand="-"),
                  _g("GOI_T", 400, 500, strand="+"),
                  _g("gC", 600, 700, strand="+"), _g("gD", 800, 900, strand="-"),
                  _g("gE", 1000, 1100, strand="+"), _g("gF", 1200, 1300, strand="-")]
    home = {"is_home": True, "genes": home_genes, "offset": 450, "breaks": []}

    def t(name, hid, s, e, strand, chrom, **kw):
        return _g(name, s, e, chrom=chrom, home_gene_id=hid, strand=strand, **kw)
    # scaffold s1 (GOI side): home order A B GOI C read backwards -> C GOI B A
    s1 = [t("c", "gC", 0, 100, "-", "s1"),
          t("goi", "GOI_T", 200, 300, "-", "s1", role="goi",
            exon_coords=[(200, 230), (260, 300)]),
          t("b", "gB", 400, 500, "+", "s1"), t("a", "gA", 600, 700, "-", "s1")]
    # scaffold s2: home order D E F read backwards -> F E D
    s2 = [t("f", "gF", 30000, 30100, "+", "s2"), t("e", "gE", 30200, 30300, "-", "s2"),
          t("d", "gD", 30400, 30500, "+", "s2")]
    tgt = {"is_home": False, "genes": s1 + s2, "offset": 250,
           "breaks": [{"x": 15000, "gap_size": 0, "is_chrom_break": True}]}
    return [home, tgt]


def test_each_reversed_scaffold_is_flipped_independently():
    tracks = _orient_fixture()
    n = ps._orient_tracks_to_home(tracks)
    tgt = tracks[1]
    assert n == 2 and sorted(tgt["flipped_chroms"]) == ["s1", "s2"]
    home_strand = {g["name"]: g["strand"] for g in tracks[0]["genes"]}
    for g in tgt["genes"]:
        if g.get("role") != "goi":
            assert g["strand"] == home_strand[g["home_gene_id"]]
            assert g["genomic_strand"] != g["strand"]
    # Order now follows home within each scaffold.
    s2 = sorted((g for g in tgt["genes"] if g["chrom"] == "s2"), key=lambda g: g["start_plot"])
    assert [g["home_gene_id"] for g in s2] == ["gD", "gE", "gF"]
    # Scaffold stays in its own plot extent; the chromosome break does not move.
    assert min(g["start_plot"] for g in s2) == 30000
    assert tgt["breaks"][0]["x"] == 15000
    # The GOI exon model is mirrored within the gene and the track re-centred on it.
    goi = next(g for g in tgt["genes"] if g.get("role") == "goi")
    assert goi["exon_coords"] == [(200, 240), (270, 300)]
    assert tgt["offset"] == pytest.approx((goi["start_plot"] + goi["end_plot"]) / 2)


def test_partial_inversion_is_not_flipped():
    tracks = _orient_fixture()
    tgt = tracks[1]
    # Put s2 back in home orientation except one gene (1 of 3 opposite).
    for g, (hid, strand, pos) in zip(
            [g for g in tgt["genes"] if g["chrom"] == "s2"],
            [("gD", "-", 30000), ("gE", "+", 30200), ("gF", "+", 30400)]):
        g.update(home_gene_id=hid, strand=strand, start_plot=pos, end_plot=pos + 100)
    ps._orient_tracks_to_home(tracks)
    assert "s2" not in tgt.get("flipped_chroms", [])


# ------------------------------------------- ribbon: no invented exons ----

def _ribbon_args():
    return SimpleNamespace(plot_width=0, plot_height=0, scale_bar_len=10000,
                           pub_width=89, pub_palette="okabe_ito")


def _ribbon_tracks(target_gene):
    home = {"label": "Home", "is_home": True, "genome_id": "home", "goi_status": "resolved",
            "offset": 0, "breaks": [],
            "genes": [_g("gene-FL", 0, 6000, home_gene_id="gene-FL", identity=100.0),
                      _g("GOI_T", 8000, 8300, home_gene_id="GOI_T", identity=100.0)]}
    tgt = {"label": "T", "is_home": False, "genome_id": "t", "goi_status": "resolved",
           "offset": 0, "breaks": [],
           "genes": [target_gene,
                     _g("GOI_T|t", 8000, 8300, home_gene_id="GOI_T", role="goi",
                        identity=80.0, confidence="HIGH")]}
    return [home, tgt]


def _gene_groups(html, track):
    """Inner SVG of each gene group drawn on one track, in document order."""
    return re.findall(r'<g class="gene-group"[^>]*data-track="%d"[^>]*>(.*?)</g>' % track,
                      html, re.S)


def test_ribbon_does_not_invent_exons_from_exon_count():
    fl = _g("gene-FL|t_flank_ann", 0, 6000, home_gene_id="gene-FL", identity=90.0,
            n_exons=6, evidence_type="flanking_miniprot")
    html = ps.render_synteny_html(_ribbon_tracks(fl), {"gene-FL": "#4e79a7"}, {}, {},
                                  _ribbon_args(), [], 0, 0, 1)
    flank = _gene_groups(html, 1)[0]
    # Exons=6 but no CDS rows: one arrow, no intron backbone.
    assert flank.count('class="exon"') == 1
    assert 'class="intron-line"' not in flank


def test_ribbon_draws_real_exons_and_thin_hit_chains():
    model = _g("gene-FL|t_flank_ann", 0, 6000, home_gene_id="gene-FL", identity=90.0,
               evidence_type="flanking_miniprot",
               exon_coords=[(0, 1000), (2500, 3500), (5000, 6000)])
    html = ps.render_synteny_html(_ribbon_tracks(model), {"gene-FL": "#4e79a7"}, {}, {},
                                  _ribbon_args(), [], 0, 0, 1)
    assert _gene_groups(html, 1)[0].count('class="exon"') == 3

    hits = dict(model, evidence_type="flanking_hit_span")
    html2 = ps.render_synteny_html(_ribbon_tracks(hits), {"gene-FL": "#4e79a7"}, {}, {},
                                   _ribbon_args(), [], 0, 0, 1)
    chain = _gene_groups(html2, 1)[0]
    assert chain.count('class="exon"') == 3
    heights = [float(h) for h in re.findall(r'<rect[^>]*height="([\d.]+)"[^>]*class="exon"', chain)]
    assert heights and max(heights) < 30   # thinner than a spliced model's GENE_H (30)


# ------------------------------------------------- many GOI models in a cell ----

def _goi(start, conf="LOW", ident=40.0, chrom="chr1", status="complete"):
    return _g(f"GOI_T|m{start}", start, start + 300, chrom=chrom, role="goi",
              home_gene_id="GOI_T", confidence=conf, identity=ident,
              model_status=status, evidence_type="exon_annotation",
              exon_coords=[(start, start + 100), (start + 200, start + 300)])


def test_cell_with_few_models_draws_all_of_them():
    genes = [_goi(1000), _goi(2000), _goi(3000)]
    models, n, frags, nf = ps._goi_cell_items(genes)
    assert len(models) == 3 and n == 3 and not frags and nf == 0


def test_cell_with_many_models_draws_best_and_counts_all():
    genes = [_goi(s) for s in range(1000, 46000, 1000)] + [_goi(99000, conf="HIGH", ident=70.0)]
    genes += [_goi(s, status="fragment") for s in (50000, 51000)]
    models, n, frags, nf = ps._goi_cell_items(genes)
    assert n == 46 and [g["confidence"] for g in models] == ["HIGH"]
    assert nf == 2 and len(frags) == 2


def test_grid_prints_count_for_collapsed_cell():
    home = {"label": "Home", "species": "Apis mellifera", "is_home": True,
            "genome_id": "home", "offset": 0, "breaks": [],
            "genes": [_g("gene-A", 1000, 2000), _g("GOI_T", 3000, 3500)]}
    tgt_genes = [_g("a", 1000, 2000, home_gene_id="gene-A", identity=90.0, confidence="HIGH")]
    tgt_genes += [_goi(s) for s in range(5000, 12000, 1000)]
    tgt = {"label": "T", "species": "Bombus terrestris", "is_home": False,
           "genome_id": "bter", "offset": 0, "breaks": [], "genes": tgt_genes}
    html = ps.render_anchor_grid([home, tgt], {"gene-A": "#4e79a7"}, {}, {}, SimpleNamespace())
    assert ">×7</text>" in html
    assert "best of 7 GOI models in this neighbourhood (7 LOW)" in html


def test_capped_goi_models_inside_the_view_are_restored():
    view = [_g("fl1", 1000, 2000), _g("fl2", 9000, 10000)]
    inside, outside, other_chrom = _goi(5000), _goi(20000), _goi(5000, chrom="chr9")
    assert ps._restore_goi_in_view(view, [inside, outside, other_chrom]) == [inside]


def test_cap_keeps_the_best_goi_models():
    genes = [_goi(1000, conf="HIGH", ident=80.0), _goi(2000), _goi(3000), _g("fl", 0, 100)]
    kept, dropped = ps._cap_goi_genes(genes, 1)
    assert [g["confidence"] for g in kept if g.get("role") == "goi"] == ["HIGH"]
    assert len(dropped) == 2 and any(g["name"] == "fl" for g in kept)
