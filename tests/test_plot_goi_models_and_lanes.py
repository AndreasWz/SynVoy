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
import synvoy_grid as sg  # noqa: E402

sg._PS = ps      # the grid module is handed plot_synteny by render(); set here for its helpers


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
    w, exons, introns = sg.model_geometry(_model(), 0.1, caret=2.4)
    assert [round(x1 - x0, 6) for x0, x1 in exons] == [10.0, 30.0, 10.0]
    assert all(round(x1 - x0, 6) == 2.4 for x0, x1 in introns)
    assert w == pytest.approx(50 + 2 * 2.4)


def test_panel_model_is_drawn_five_prime_first():
    # The 5' exon of a '-' gene is its genomically LAST block.
    asym = dict(_model("-"), exon_coords=[(1000, 1099), (2001, 2300)])
    _, exons, _ = sg.model_geometry(asym, 0.1, flip=True)
    assert round(exons[0][1] - exons[0][0], 6) == 30.0
    svg, _ = sg.draw_model(asym, 0, 20, 0.1, "#f00", "#800", "", 5.0, aligned=True)
    assert svg.count("<path") == 1          # one arrow tip, on the last exon drawn


def test_a_tiny_exon_keeps_a_visible_width():
    micro = dict(_model(), exon_coords=[(1000, 1002), (2001, 2300)])
    _, exons, _ = sg.model_geometry(micro, 0.01)
    assert exons[0][1] - exons[0][0] == sg.MIN_EXON


def test_hit_span_cannot_set_the_scale():
    # A hit-based call has no exon structure; its span is cut to the longest gene model.
    hit = dict(_model(), evidence_type="tandem_copy", exon_coords=[(1000, 6600)])
    assert sg.model_geometry(hit, 0.1, cap_bp=500)[0] == pytest.approx(50.0)
    assert sg.model_geometry(_model(), 0.1, cap_bp=100)[0] == pytest.approx(50 + 2 * 2.4)


def test_hit_chains_are_not_drawn_as_gene_models():
    hit = dict(_model(), evidence_type="fallback_hit_span")
    svg_model, _ = sg.draw_model(_model(), 0, 20, 0.1, "#f00", "#800", "", 5.0)
    svg_hit, _ = sg.draw_model(hit, 0, 20, 0.1, "#f00", "#800", "", 5.0)
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


def _grid_pair(target_goi, home_goi=None):
    home = {"label": "Home", "species": "Apis mellifera", "is_home": True,
            "genome_id": "home", "offset": 0, "breaks": [],
            "genes": [_g("gene-A", 1000, 2000)] + (home_goi or [_g("GOI_T", 3000, 3500)])}
    tgt_genes = [_g("a", 1000, 2000, home_gene_id="gene-A", identity=90.0, confidence="HIGH")]
    tgt = {"label": "T", "species": "Bombus terrestris", "is_home": False,
           "genome_id": "bter", "offset": 0, "breaks": [], "genes": tgt_genes + target_goi}
    return [home, tgt]


def _grid(tracks, **style):
    return sg.render(ps, tracks, {"gene-A": "#4e79a7"}, {}, SimpleNamespace(), sg.style_for(**style))


def _drawn(result, row=2):
    return [r["drawn_as"] for r in result["source"] if r["row"] == row and r["column"] == "T"]


def test_copy_is_the_best_call_or_at_least_medium():
    best, medium, low = _goi(1000, "LOW", 60.0), _goi(2000, "MEDIUM"), _goi(3000, "LOW")
    frag = _goi(4000, "HIGH", status="fragment")
    assert ps._goi_is_copy(best, best)                  # the best call, whatever its tier
    assert ps._goi_is_copy(medium, best)
    assert not ps._goi_is_copy(low, best)
    assert not ps._goi_is_copy(frag, best)              # a fragment is never a copy


def test_one_copy_per_genome_puts_the_model_in_the_goi_column():
    res = _grid(_grid_pair([_goi(5000, "HIGH", 70.0)]))
    assert res["info"]["goi_form"] == "column"
    assert _drawn(res) == ["model"]
    assert "gene models" not in "".join(res["parts"])      # no panel header


def test_two_copies_get_numbered_arrows_and_a_panel():
    res = _grid(_grid_pair([_goi(5000, "HIGH", 70.0), _goi(7000, "MEDIUM", 55.0)]))
    assert res["info"]["goi_form"] == "panel"
    assert res["info"]["most_copies_in_a_target"] == 2
    rows = [r for r in res["source"] if r["row"] == 2 and r["column"] == "T"]
    assert [(r["drawn_as"], r["copy_number"]) for r in rows] == [("model", 1), ("model", 2)]
    assert rows[0]["start"] < rows[1]["start"]             # numbered along the chromosome
    svg = "".join(res["parts"])
    assert "T gene models" in svg or ">gene models<" in svg


def test_two_goi_genes_in_the_home_genome_also_need_the_panel():
    home_goi = [_g("GOI_T", 3000, 3500), _g("GOI_T2", 4000, 4500, home_gene_id="GOI_T")]
    tracks = _grid_pair([_goi(5000, "HIGH", 70.0)], home_goi=home_goi)
    ps._GOI_NAMES.update({"GOI_T", "GOI_T2"})
    try:
        res = _grid(tracks)
    finally:
        ps._GOI_NAMES.difference_update({"GOI_T", "GOI_T2"})
    assert res["info"]["home_copies"] == 2 and res["info"]["goi_form"] == "panel"


def test_weak_calls_are_marks_and_never_counted_as_copies():
    # Seven LOW calls: the best one is the copy, the other six are weak calls.
    res = _grid(_grid_pair([_goi(s) for s in range(5000, 12000, 1000)]))
    assert res["info"]["goi_copies"] == 1 and res["info"]["goi_calls"] == 7
    assert sorted(_drawn(res)) == ["mark"] * 6 + ["model"]
    assert "×7" not in "".join(res["parts"])


def test_many_weak_calls_are_capped_and_the_rest_is_counted():
    genes = [_goi(s) for s in range(1000, 46000, 1000)] + [_goi(99000, conf="HIGH", ident=70.0)]
    res = _grid(_grid_pair(genes))
    drawn = _drawn(res)
    assert drawn.count("model") == 1 and drawn.count("mark") == sg.GridStyle().max_marks
    assert drawn.count("counted") == 45 - sg.GridStyle().max_marks
    assert f">+{45 - sg.GridStyle().max_marks}<" in "".join(res["parts"])


def test_a_long_tandem_array_draws_the_first_models_and_says_how_many_more():
    genes = [_goi(s, "HIGH", 70.0) for s in range(5000, 17000, 1000)]      # 12 copies
    res = _grid(_grid_pair(genes), max_models=10)
    assert res["info"]["goi_copies"] == 12 and res["info"]["copies_not_drawn"] == 2
    drawn = _drawn(res)
    assert drawn.count("model") == 10 and drawn.count("arrow") == 2    # every copy keeps its arrow
    assert ">+2 more<" in "".join(res["parts"])


def test_full_grid_tooltip_names_weak_calls():
    html = ps.render_anchor_grid(_grid_pair([_goi(s) for s in range(5000, 12000, 1000)]),
                                 {"gene-A": "#4e79a7"}, {}, {}, SimpleNamespace())
    assert html.count("(weak call") == 6       # three at the locus, three elsewhere
    assert html.count("(weak call, elsewhere in the genome)") == 3
    assert "×7" not in html


def test_capped_goi_models_inside_the_view_are_restored():
    view = [_g("fl1", 1000, 2000), _g("fl2", 9000, 10000)]
    inside, outside, other_chrom = _goi(5000), _goi(20000), _goi(5000, chrom="chr9")
    assert ps._restore_goi_in_view(view, [inside, outside, other_chrom]) == [inside]


def test_cap_keeps_the_best_goi_models():
    genes = [_goi(1000, conf="HIGH", ident=80.0), _goi(2000), _goi(3000), _g("fl", 0, 100)]
    kept, dropped = ps._cap_goi_genes(genes, 1)
    assert [g["confidence"] for g in kept if g.get("role") == "goi"] == ["HIGH"]
    assert len(dropped) == 2 and any(g["name"] == "fl" for g in kept)
