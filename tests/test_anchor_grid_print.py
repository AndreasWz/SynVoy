"""The print version of the anchor grid, its legend and the shared legend builder.

The screen figures used to be laid out in pixels and scaled to the page, which
put their text at 2-5 pt in a 183 mm figure. ``bin/synvoy_grid.py`` lays the grid
out in points at final size; these tests pin what a journal checks (size in mm,
text of 5-7 pt) and what the figure must not hide (dropped columns, copies).
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402
import synvoy_grid as sg  # noqa: E402
import synvoy_legend as sl  # noqa: E402


def _g(name, start, end, **kw):
    g = {"chrom": "chr1", "start": start, "end": end, "start_plot": start,
         "end_plot": end, "name": name, "strand": "+", "exon_coords": []}
    g.update(kw)
    return g


def _goi(start, conf="HIGH", ident=70.0, chrom="chr1", **kw):
    return _g(f"GOI_T|m{start}", start, start + 300, chrom=chrom, role="goi",
              home_gene_id="GOI_T", confidence=conf, identity=ident,
              evidence_type="exon_annotation",
              exon_coords=[(start, start + 100), (start + 200, start + 300)], **kw)


def _tracks(n_flank=6, target_goi=None, n_targets=2):
    names = [f"gene-F{k}" for k in range(n_flank)]
    home_genes = [_g(n, 1000 + 2000 * k, 2000 + 2000 * k) for k, n in enumerate(names)]
    home_genes.insert(n_flank // 2, _g("GOI_T", 1500 + 2000 * (n_flank // 2 - 1) + 600,
                                       1500 + 2000 * (n_flank // 2 - 1) + 900))
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True,
               "genome_id": "home", "offset": 0, "breaks": [], "genes": home_genes}]
    for t in range(n_targets):
        genes = [_g(f"t{t}_{k}", 1000 + 2000 * k, 2000 + 2000 * k, home_gene_id=n,
                    identity=80.0 - 5 * t, confidence="HIGH", query_coverage=0.95)
                 for k, n in enumerate(names)]
        genes += [dict(g) for g in (target_goi if target_goi is not None else [_goi(7000)])]
        tracks.append({"label": f"T{t}", "species": f"Bombus species{t}", "is_home": False,
                       "genome_id": f"g{t}", "offset": 0, "breaks": [], "genes": genes})
    return tracks, {n: ps.GENE_PALETTE[k % len(ps.GENE_PALETTE)] for k, n in enumerate(names)}


def _print(tracks, colours, **args):
    return ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace(**args))


def _mm(svg, attr):
    return float(re.search(rf'\b{attr}="([0-9.]+)mm"', svg).group(1))


# ------------------------------------------------------------ what a journal checks ----

def test_print_svg_states_its_size_in_mm_and_parses():
    svg, legend, source, info = _print(*_tracks())
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert _mm(svg, "width") <= 183.0 and _mm(svg, "height") <= 170.0
    assert info["width_mm"] == pytest.approx(_mm(svg, "width"), abs=0.06)
    assert "Arial" in svg


def test_print_text_is_five_to_seven_points():
    for kw in ({}, {"print_grid_numbers": True}):
        svg, _, _, info = _print(*_tracks(n_flank=30, n_targets=8), **kw)
        sizes = {float(x) for x in re.findall(r'font-size="([0-9.]+)"', svg)}
        assert sizes <= {5.0, 6.0, 7.0}, sizes
        assert info["text_pt"] == sorted(sizes)


def test_print_figure_carries_no_title_no_location_no_tooltip():
    svg, legend, source, _ = _print(*_tracks())
    assert "Anchor-grid synteny" not in svg and "<title>" not in svg
    assert "chr1:" not in svg
    # ... they are in the legend text and the source table instead.
    assert "Columns are the 6 genes flanking T in Apis mellifera" in legend
    assert {r["scaffold"] for r in source if r["state"] == "placed"} == {"chr1"}


def test_narrow_figure_drops_columns_and_says_which():
    tracks, colours = _tracks(n_flank=30, n_targets=3)
    # one flanking gene is found in a single genome: the first to go
    tracks[2]["genes"] = [g for g in tracks[2]["genes"] if g.get("home_gene_id") != "gene-F29"]
    tracks[3]["genes"] = [g for g in tracks[3]["genes"] if g.get("home_gene_id") != "gene-F29"]
    svg, legend, _, info = _print(tracks, colours, print_width_mm=89.0)
    assert _mm(svg, "width") <= 89.0 + 0.01
    assert info["dropped_columns"] and info["dropped_columns"][0] == "F29"
    assert info["flanking_columns"] == 30 - len(info["dropped_columns"])
    assert "left out to fit the figure width" in legend and "F29" in legend
    assert {float(x) for x in re.findall(r'font-size="([0-9.]+)"', svg)} <= {5.0, 6.0, 7.0}


def test_many_rows_are_packed_before_the_height_limit_is_passed():
    svg, _, _, info = _print(*_tracks(n_targets=40))
    assert info["rows"] == 41 and info["fits_height"]
    assert info["row_pitch_pt"] < sg.GridStyle().row_h


# ------------------------------------------------------------------- the GOI column ----

def test_goi_form_follows_the_copy_number():
    assert _print(*_tracks())[3]["goi_form"] == "column"
    two = [_goi(7000), _goi(9000, "MEDIUM", 55.0)]
    info = _print(*_tracks(target_goi=two))[3]
    assert info["goi_form"] == "panel" and info["genomes_with_several_copies"] == 2


def test_copy_elsewhere_follows_a_break_mark():
    calls = [_goi(7000), _goi(500_000, "MEDIUM", 50.0, chrom="chr9")]
    svg, _, source, _ = _print(*_tracks(target_goi=calls))
    rows = [r for r in source if r["row"] == 2 and r["column"] == "T"]
    assert [r["position"] for r in rows] == ["locus", "elsewhere"]
    assert svg.count(">//</text>") >= 2            # in the grid and in the panel, per target


def test_numbers_are_an_option_of_the_print_grid():
    tracks, colours = _tracks()
    tracks[1]["genes"][0]["query_coverage"] = 0.4
    plain = _print(tracks, colours)[0]
    numbered = _print(tracks, colours, print_grid_numbers=True)[0]
    assert ">80<" not in plain and "(40)" not in plain
    assert ">80<" in numbered and ">(40)<" in numbered      # (40) = query coverage, below 80 %


def _low_coverage_tracks(n_targets):
    tracks, colours = _tracks(n_targets=n_targets)
    for track in tracks[1:]:
        track["genes"][0].update(identity=60.0, query_coverage=0.4)
    return tracks, colours


def test_low_coverage_is_written_under_the_arrow_and_makes_the_row_taller():
    plain = _print(*_tracks(n_targets=3), print_grid_numbers=True)
    low = _print(*_low_coverage_tracks(3), print_grid_numbers=True)
    assert low[0].count('class="cov-low"') == 3 and low[0].count(">(40)<") == 3
    assert ">(44)<" in low[0] and ">(44)<" not in plain[0]         # the legend key
    assert "*<" not in low[0]
    assert low[3]["low_coverage_shown_as"] == "line" and plain[3]["low_coverage_shown_as"] == ""
    extra_pt = (low[3]["height_mm"] - plain[3]["height_mm"]) * sg.PT_PER_MM
    # three rows, each at most COV_LINE_H taller, plus one more legend line (5 pt x 1.5)
    assert 3 * 2.0 < extra_pt <= 3 * sg.COV_LINE_H + 7.5 + 0.3
    assert "in brackets under an arrow is the query coverage" in low[1]
    assert {r["query_coverage_pct"] for r in low[2]
            if r["column"] == "F0" and r["state"] == "placed"} == {"40"}


def test_coverage_line_gives_way_to_a_flag_when_the_height_limit_is_hit():
    tracks, colours = _low_coverage_tracks(10)
    roomy = _print(tracks, colours, print_grid_numbers=True)
    assert roomy[3]["low_coverage_shown_as"] == "line" and roomy[3]["fits_height"]
    tight = _print(tracks, colours, print_grid_numbers=True,
                   print_height_mm=roomy[3]["height_mm"] - 8.0)
    assert tight[3]["low_coverage_shown_as"] == "star" and tight[3]["fits_height"]
    assert 'class="cov-low"' not in tight[0] and tight[0].count(">60*<") == 10
    assert "with * when the query coverage is below 80 %" in tight[1]
    assert "query coverage below 80 %" in tight[0]             # the legend key of the flag


def test_no_coverage_without_an_identity_number():
    tracks, colours = _tracks(n_targets=2)
    tracks[1]["genes"][0].update(identity=20.0, query_coverage=0.4)     # too pale for a number
    svg, _, _, info = _print(tracks, colours, print_grid_numbers=True)
    assert "(40)" not in svg and info["low_coverage_shown_as"] == ""


@pytest.mark.parametrize("cov, pct", [(0.07, 7), (0.404, 40), (0.796, 79), (0.7999, 79),
                                       (0.80, 80), (0.996, 100)])
def test_a_flagged_coverage_is_never_printed_as_the_threshold(cov, pct):
    assert ps._coverage_pct(cov) == pct


def test_source_table_has_one_row_per_cell_and_call():
    tracks, colours = _tracks(target_goi=[_goi(7000), _goi(9000, "LOW", 30.0)])
    _, _, source, _ = _print(tracks, colours)
    assert set(source[0]) == set(sg.SOURCE_COLUMNS)
    goi_rows = [r for r in source if r["column"] == "T"]
    assert len(goi_rows) == 1 + 2 * 2                      # home gene + two calls per target
    assert [r["drawn_as"] for r in goi_rows if r["row"] == 2] == ["model", "mark"]
    assert len([r for r in source if r["column"] != "T"]) == 3 * 6


# --------------------------------------------------------------------------- legend ----

def test_grid_legend_is_grouped():
    tracks, colours = _tracks(target_goi=[_goi(7000), _goi(9000, "LOW", 30.0)])
    res = sg.render(ps, tracks, colours, {}, SimpleNamespace(), sg.style_for())
    titles = [t for t, _ in res["legend_groups"]]
    assert titles == ["Flanking genes", "GOI: T", "Numbers"]
    full = sg.render(ps, tracks, colours, {}, SimpleNamespace(), sg.style_for(full=True))
    assert [t for t, _ in full["legend_groups"]][0] == "Layout"
    goi_texts = [text for _, text in dict(res["legend_groups"])["GOI: T"]]
    assert any("weak call" in t for t in goi_texts)


def test_legend_groups_sit_side_by_side_and_wrap_to_a_second_band():
    groups = [("One", [(None, "alpha beta gamma")]), ("Two", [("×N", "copies")]),
              ("Three", [((10.0, lambda x, y: f'<rect x="{x}" y="{y}"/>'), "a box")])]
    wide, h_wide = sl.legend_svg(groups, 0, 0, 1000, ps.text_width, fs=10)
    narrow, h_narrow = sl.legend_svg(groups, 0, 0, 120, ps.text_width, fs=10)
    assert h_narrow > h_wide
    heads = [float(x) for x in re.findall(r'<text x="([0-9.]+)"[^>]*font-weight="700"[^>]*>(?:One|Two|Three)<',
                                         "".join(wide))]
    assert heads == sorted(heads) and len(set(heads)) == 3     # three columns
    assert sl.natural_width(groups, ps.text_width, fs=10) <= 1000
    assert "<rect" in "".join(wide)


def test_legend_wraps_long_text_and_skips_empty_groups():
    long = "word " * 40
    parts, height = sl.legend_svg([("Empty", []), ("Text", [("Key", long)])], 0, 0, 300,
                                  ps.text_width, fs=10)
    lines = [t for t in re.findall(r">([^<]+)</text>", "".join(parts)) if t.startswith("word")]
    assert len(lines) > 1 and all(ps.text_width(t, 10) <= 300 for t in lines)
    assert "Empty" not in "".join(parts)
    assert sl.legend_svg([], 0, 0, 300, ps.text_width) == ([], 0.0)
    assert sl.legend_text([("Text", [("Key", "meaning")])]) == "Text: Key = meaning."


def test_every_screen_figure_has_a_grouped_legend():
    tracks, colours = _tracks()
    args = SimpleNamespace()
    for html in (ps.render_anchor_grid(tracks, colours, {}, {}, args),
                 ps.render_anchor_grid_positional(tracks, colours, {}, {}, args),
                 ps.render_anchor_grid_threaded(tracks, colours, {}, {}, args)):
        texts = re.findall(r">([^<]+)</text>", html)
        assert "Layout" in texts and "Flanking genes" in texts and "GOI: T" in texts
        assert not any(" · " in t and len(t) > 60 for t in texts), "a dotted list is back"


# ------------------------------------------------------------------------- the tree ----

def test_rows_are_matched_to_the_species_tree_by_species_name():
    track = {"genome_id": "chicken", "species": "Gallus gallus"}
    keys = ["Bos_taurus", "Gallus_gallus"]
    assert ps._track_species_key(track, keys) == "Gallus_gallus"
    assert ps._track_species_key({"genome_id": "Bos_taurus", "species": ""}, keys) == "Bos_taurus"
    assert ps._track_species_key({"genome_id": "zebrafish", "species": "Danio rerio"}, keys) is None


def test_species_tree_puts_the_home_genome_into_the_cladogram():
    tracks, colours = _tracks(n_targets=3)
    for t, name in zip(tracks[1:], ("Bombus terrestris", "Vespa crabro", "Bombus impatiens")):
        t["species"] = name
    args = SimpleNamespace(_species_tree_newick=(
        "((Apis_mellifera,(Bombus_terrestris,Bombus_impatiens)),Vespa_crabro);"))
    res = sg.render(ps, tracks, colours, {}, args, sg.style_for())
    svg = "".join(res["parts"])
    assert 'class="grid-tree"' in svg and 'class="tree-unplaced"' not in svg
    order = [r["species"] for r in res["source"] if r["column"] == "F0"]
    assert order == ["Apis mellifera", "Bombus terrestris", "Bombus impatiens", "Vespa crabro"]
