"""A searched genome in which no GOI call was made is a result, and the figures say so.

Until 2026-10-09 the pipeline left such genomes out of every figure
(``hide_goi_absent_tracks = true``): on 26 finished loci, 16 of 90 searched genomes
were not in the anchor grid, among them one with three blocks of flanking genes. Now
every searched genome is drawn; these tests pin what a reader gets for such a genome.
"""
import os
import re
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402
import synvoy_grid as sg  # noqa: E402

REPO = os.path.join(os.path.dirname(__file__), "..")


def _g(name, start, end, **kw):
    g = {"chrom": "chr1", "start": start, "end": end, "start_plot": start,
         "end_plot": end, "name": name, "strand": "+", "exon_coords": []}
    g.update(kw)
    return g


def _goi(start):
    return _g(f"GOI_T|m{start}", start, start + 300, role="goi", home_gene_id="GOI_T",
              confidence="HIGH", identity=70.0, evidence_type="exon_annotation",
              exon_coords=[(start, start + 100), (start + 200, start + 300)])


def _tracks(n_flank=6, n_targets=3):
    names = [f"gene-F{k}" for k in range(n_flank)]
    home_genes = [_g(n, 1000 + 2000 * k, 2000 + 2000 * k) for k, n in enumerate(names)]
    home_genes.insert(n_flank // 2, _g("GOI_T", 6100, 6400))
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True,
               "genome_id": "home", "offset": 0, "breaks": [], "genes": home_genes}]
    for t in range(n_targets):
        genes = [_g(f"t{t}_{k}", 1000 + 2000 * k, 2000 + 2000 * k, home_gene_id=n,
                    identity=80.0, confidence="HIGH", query_coverage=0.95)
                 for k, n in enumerate(names)]
        genes.append(_goi(7000))
        tracks.append({"label": f"T{t}", "species": f"Bombus species{t}", "is_home": False,
                       "genome_id": f"g{t}", "offset": 0, "breaks": [], "genes": genes})
    return tracks, {n: ps.GENE_PALETTE[k % len(ps.GENE_PALETTE)] for k, n in enumerate(names)}


def _without_goi():
    """Target 2 keeps two flanking genes and no GOI; in target 3 nothing was placed."""
    tracks, colours = _tracks()
    tracks[2]["genes"] = [g for g in tracks[2]["genes"] if g.get("role") != "goi"][:2]
    tracks[3]["genes"] = []
    return tracks, colours


def _dashes(svg):
    return len(re.findall(r">%s</text>" % re.escape(sg.NO_CALL_MARK), svg))


def test_every_searched_genome_keeps_its_row():
    svg, legend, source, info = ps.render_anchor_grid_print(*_without_goi(), {}, SimpleNamespace())
    assert info["rows"] == 4
    assert {r["species"] for r in source} == {"Apis mellifera", "Bombus species0",
                                              "Bombus species1", "Bombus species2"}
    # the flanking genes that were placed stay in the row; the rest reads "not placed"
    placed = [r for r in source if r["species"] == "Bombus species1" and r["state"] == "placed"]
    assert len(placed) == 2
    assert all(r["state"] == "not placed" for r in source if r["species"] == "Bombus species2")


def test_goi_cell_says_that_no_model_was_placed():
    with_calls, _, _, _ = ps.render_anchor_grid_print(*_tracks(), {}, SimpleNamespace())
    assert _dashes(with_calls) == 0                      # nothing to say: no mark, no legend entry
    svg, _, _, _ = ps.render_anchor_grid_print(*_without_goi(), {}, SimpleNamespace())
    assert _dashes(svg) == 3                             # two cells and the legend entry
    assert "no model placed (not evidence of absence)" in svg


def test_legend_text_names_the_genomes_and_what_was_placed():
    _, legend, _, info = ps.render_anchor_grid_print(*_without_goi(), {}, SimpleNamespace())
    assert info["home_flanking_genes"] == 6
    assert [(m["species"], m["genome_id"], m["flanking_placed"])
            for m in info["rows_without_goi"]] == [("Bombus species1", "g1", 2),
                                                   ("Bombus species2", "g2", 0)]
    assert "No T model was placed in 2 of the 3 searched genomes" in legend
    assert "Bombus species1 (2 of 6 flanking genes placed)" in legend
    assert "Bombus species2 (0 of 6 flanking genes placed)" in legend
    assert "not evidence that the gene is absent" in legend
    _, legend, _, info = ps.render_anchor_grid_print(*_tracks(), {}, SimpleNamespace())
    assert info["rows_without_goi"] == [] and "No T model" not in legend


def test_full_grid_says_it_in_the_subtitle_and_where_the_location_would_stand():
    html = ps.render_anchor_grid(*_without_goi(), {}, {}, SimpleNamespace())
    assert "no T model in 2 genomes" in html
    assert html.count(sg.NOTHING_PLACED) == 1            # the row in which nothing was placed
    assert "No T model was placed in Bombus species2" in html      # tooltip of the dash
    html = ps.render_anchor_grid(*_tracks(), {}, {}, SimpleNamespace())
    assert "no T model" not in html and sg.NOTHING_PLACED not in html


def test_all_genomes_without_a_call_still_gives_a_figure():
    tracks, colours = _tracks()
    for track in tracks[1:]:
        track["genes"] = [g for g in track["genes"] if g.get("role") != "goi"]
    svg, legend, _, info = ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace())
    assert len(info["rows_without_goi"]) == 3 and _dashes(svg) == 4
    assert "No T model was placed in 3 of the 3 searched genomes" in legend


def test_genomes_without_a_call_go_below_the_others_in_the_same_order():
    """The ribbon plot links each track to the next one: a nearly empty track between
    two others would cut the ribbons that join them."""
    tracks = [{"genome_id": name, "goi_status": status} for name, status in (
        ("a", "resolved"), ("b", "absent"), ("c", "ambiguous"), ("d", "absent"), ("e", "resolved"))]
    assert [t["genome_id"] for t in ps._absent_tracks_last(tracks)] == ["a", "c", "e", "b", "d"]
    assert ps._absent_tracks_last([]) == []


def test_pipeline_and_redraw_script_draw_every_genome_by_default():
    config = open(os.path.join(REPO, "nextflow.config")).read()
    assert re.search(r"^\s*hide_goi_absent_tracks\s*=\s*false\b", config, re.M)
    replot = open(os.path.join(REPO, "scripts", "replot.sh")).read()
    assert "--hide_goi_absent" not in replot             # same default as the pipeline
