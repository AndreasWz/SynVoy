"""The clade version of the print anchor grid (synvoy_grid.clade_style / goi_clades):
GOI copies, their arrows and their exon models in the panel, coloured by their clade in the
GOI gene tree, written only when it shows more than the species tree beside the rows."""
import os
import re
import sys
import types

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import synvoy_grid as sg  # noqa: E402
import synvoy_tree as stree  # noqa: E402

# Two well separated groups: the home gene with A and B, and C with D.
TREE = ("((GOI_X:0.1,(GOI_X|gA_fna_b0_l1:0.1,GOI_X|gB_fna_b0_l1:0.1):0.1):1.0,"
        "(GOI_X|gC_fna_b0_l1:0.1,GOI_X|gD_fna_b0_l1:0.1):1.0);")


def _args(tmp_path, tree=TREE):
    path = tmp_path / "goi_tree.nwk"
    path.write_text(tree)
    return types.SimpleNamespace(tree=str(path), clade_count=2)


def _call(gff_id):
    return ({"gff_id": gff_id}, False, True)


def _rows(*genomes):
    """Home row, then one row per genome with one copy each."""
    return [[({"gff_id": ""}, False, True)]] + [[_call(f"GOI_X|{g}_fna_b0_l1")] for g in genomes]


def _clade_of(tmp_path, leaf):
    t = stree.midpoint_root(stree.parse_newick_tree(open(_args(tmp_path).tree).read()))
    return stree.partition_clades(t, target_k=2)[leaf]


def test_leaf_key_keeps_the_part_before_the_bar():
    """Tandem copies of one model differ only before the '|'; keying on the part after it
    gave every copy of an array the clade of one of them."""
    a = sg._leaf_key("GOI_copy_1|Pogonomyrmex_barbatus_fna_b0_r2")
    b = sg._leaf_key("GOI_copy_2|Pogonomyrmex_barbatus_fna_b0_r2")
    assert a != b
    assert sg._leaf_key("GOI_X|Sp_fna_b5_l1_exon_ann") == sg._leaf_key("GOI_X|Sp_fna_b5_l1")
    assert sg._leaf_key("GOI_X") == ""


def test_each_copy_gets_the_clade_of_its_own_leaf(tmp_path):
    tree = ("((GOI_X:0.1,(GOI_copy_1|gA_fna_b0_l1:0.1,GOI_X|gB_fna_b0_l1:0.1):0.1):1.0,"
            "(GOI_copy_2|gA_fna_b0_l1:0.1,GOI_X|gD_fna_b0_l1:0.1):1.0);")
    args = _args(tmp_path, tree)
    calls = [[({"gff_id": ""}, False, True)],
             [_call("GOI_copy_1|gA_fna_b0_l1"), _call("GOI_copy_2|gA_fna_b0_l1")]]
    clade_of, info = sg.goi_clades(args, calls)
    first, second = (clade_of[id(g)] for g, _a, _c in calls[1])
    assert first != second
    assert info["rows_with_two_placed"] == 1 and info["informative"]


def test_single_copies_in_species_order_say_nothing_new(tmp_path):
    clade_of, info = sg.goi_clades(_args(tmp_path), _rows("gA", "gB", "gC", "gD"))
    assert info["copies_matched"] == info["copies"] == 5
    assert info["rows_with_two_placed"] == 0 and info["rows_out_of_pattern"] == 0
    assert not info["informative"]


def test_a_copy_that_breaks_the_pattern_of_its_neighbours_is_informative(tmp_path):
    """Rows A, C, B, D: C sits between two rows of the other clade."""
    assert _clade_of(tmp_path, "GOI_X|gA_fna_b0_l1") == _clade_of(tmp_path, "GOI_X|gB_fna_b0_l1")
    _clade_of_c = _clade_of(tmp_path, "GOI_X|gC_fna_b0_l1")
    assert _clade_of_c != _clade_of(tmp_path, "GOI_X|gA_fna_b0_l1")
    _map, info = sg.goi_clades(_args(tmp_path), _rows("gA", "gC", "gB", "gD"))
    assert info["rows_out_of_pattern"] >= 1 and info["informative"]


def test_no_tree_means_no_clades(tmp_path):
    clade_of, info = sg.goi_clades(_args(tmp_path, "(GOI_X:0.1,GOI_X|gA_fna_b0_l1:0.1);"),
                                   _rows("gA"))
    assert clade_of == {} and info["n_clades"] == 0 and not info["informative"]


def _grid_tracks():
    """Home + two targets; target g0 carries two copies, g1 one, all leaves of TREE2."""
    def gene(name, start, end, **kw):
        return dict({"chrom": "chr1", "start": start, "end": end, "start_plot": start,
                     "end_plot": end, "name": name, "strand": "+", "exon_coords": []}, **kw)

    def goi(start, gff_id):
        return gene(f"GOI_X|m{start}", start, start + 300, role="goi", home_gene_id="GOI_X",
                    confidence="HIGH", identity=70.0, evidence_type="exon_annotation",
                    gff_id=gff_id, exon_coords=[(start, start + 100), (start + 200, start + 300)])

    names = [f"gene-F{k}" for k in range(4)]
    home = [gene(n, 1000 + 2000 * k, 2000 + 2000 * k) for k, n in enumerate(names)]
    home.insert(2, gene("GOI_X", 4600, 4900))
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home",
               "offset": 0, "breaks": [], "genes": home}]
    for t, gois in ((0, [goi(4600, "GOI_X|g0_fna_b0_l1"), goi(5200, "GOI_X|g0_fna_b0_l2")]),
                    (1, [goi(4600, "GOI_X|g1_fna_b0_l1")])):
        genes = [gene(f"t{t}_{k}", 1000 + 2000 * k, 2000 + 2000 * k, home_gene_id=n,
                      identity=80.0, confidence="HIGH", query_coverage=0.95)
                 for k, n in enumerate(names)] + gois
        tracks.append({"label": f"T{t}", "species": f"Bombus species{t}", "is_home": False,
                       "genome_id": f"g{t}", "offset": 0, "breaks": [], "genes": genes})
    return tracks, {n: "#4477aa" for n in names}


TREE2 = ("((GOI_X:0.1,(GOI_X|g0_fna_b0_l1:0.1,GOI_X|g1_fna_b0_l1:0.1):0.1):1.0,"
         "GOI_X|g0_fna_b0_l2:1.2);")


def test_clade_grid_draws_the_exon_panel_in_clade_colours(tmp_path):
    """The numbered copies of the GOI column have their exon models on the right, as in
    the print grid, each filled with the colour of its clade (user, 2026-10-14)."""
    import plot_synteny as ps
    tracks, colours = _grid_tracks()
    res = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), sg.clade_style())
    svg = sg.svg_document(res)
    assert res["info"]["goi_form"] == "panel"
    assert "gene models" in svg
    goi_rows = [r for r in res["source"]
                if r["column"] == res["info"]["goi_label"] and r["state"] == "placed"]
    assert len(goi_rows) == 3
    assert {r["drawn_as"] for r in goi_rows} == {"model"}
    assert sorted(r["clade"] for r in goi_rows if r["species"] == "Bombus species0") == [1, 2]
    for cid in (0, 1):      # an exon of a panel model is a rect in that clade's colour
        assert re.search(rf'<rect [^>]*fill="{stree.color_for_clade(cid)}"', svg)
    assert "exon structure of each copy" in res["legend_text"]


def test_compact_grid_draws_one_stack_per_side_with_the_count(tmp_path):
    """Compact version: the flanking genes of each side are one grey stack with how many
    were placed (home row: all); the source data still lists every gene."""
    import plot_synteny as ps
    tracks, colours = _grid_tracks()
    tracks[2]["genes"] = [g for g in tracks[2]["genes"] if g.get("home_gene_id") != "gene-F3"]
    res = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), sg.compact_style(inline=True))
    svg = sg.svg_document(res)
    # home 2|2, target 0 2|2, target 1 lost F3 on the right: 2|1
    assert sorted(re.findall(r">×(\d+)</text>", svg)) == ["1", "2", "2", "2", "2", "2"]
    assert res["info"]["compact"] and res["info"]["flanking_columns"] == 2
    stack_rows = [r for r in res["source"] if r["drawn_as"] == "stack"]
    assert len(stack_rows) == 4 + 4 + 3        # one row per placed gene, every genome
    assert any(r["state"] == "not placed" and r["column"] == "F3" for r in res["source"])
    assert "one grey stack per side" in res["legend_text"]


def test_compact_inline_puts_every_copy_in_the_column(tmp_path):
    """Inline: no panel; every copy's exon model stands in the GOI column, so a genome with
    two copies needs one line instead of two (or two with one model to a line)."""
    import plot_synteny as ps
    tracks, colours = _grid_tracks()
    res = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), sg.compact_style(inline=True))
    assert res["info"]["goi_form"] == "column" and res["info"]["max_lines"] == 1
    assert "gene models" not in sg.svg_document(res)
    goi = [r for r in res["source"] if r["column"] == res["info"]["goi_label"] and r["state"] == "placed"]
    assert {r["drawn_as"] for r in goi} == {"model"}
    one = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2),
                    sg.compact_style(inline=True, inline_per_line=1))
    assert one["info"]["max_lines"] == 2


def test_compact_panel_keeps_the_exon_panel(tmp_path):
    import plot_synteny as ps
    tracks, colours = _grid_tracks()
    res = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), sg.compact_style(inline=False))
    assert res["info"]["goi_form"] == "panel" and "gene models" in sg.svg_document(res)


def test_clade_tree_collapses_each_clade_to_one_tip():
    """((a,b),(c,(d,e))) cut into clades {a,b}, {c}, {d,e} -> (0, (1, 2)); a polytomy's
    leftover siblings that share a clade become one tip; pruning keeps the drawn clades."""
    t = stree.parse_newick_tree("((a:1,b:1):1,(c:1,(d:1,e:1):1):1);")
    leaf_clade = {"a": 0, "b": 0, "c": 1, "d": 2, "e": 2}
    assert sg._clade_tree(t, leaf_clade) == (0, (1, 2))
    poly = stree.parse_newick_tree("(a:1,b:1,c:1,(d:1,e:1):1);")
    assert sg._clade_tree(poly, {"a": 0, "b": 0, "c": 0, "d": 1, "e": 1}) == (0, 1)
    assert sg._prune_clade_tree((0, (1, 2)), {0, 2}) == (0, 2)
    assert sg.clade_tips((0, (1, 2))) == [0, 1, 2]


def test_legend_of_a_clade_figure_carries_the_clade_tree(tmp_path):
    import plot_synteny as ps
    tracks, colours = _grid_tracks()
    for style in (sg.clade_style(), sg.compact_style(inline=True), sg.clade_muted_style()):
        res = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), style)
        assert 'class="clade-tree"' in sg.svg_document(res)
        assert "small tree in the legend" in res["legend_text"]
    plain = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), sg.style_for())
    assert 'class="clade-tree"' not in sg.svg_document(plain)


def test_clades_muted_keeps_the_flanking_colours_but_paler(tmp_path):
    """Version 'clades_muted': the print grid's home-gene colours, muted; GOI copies by clade."""
    import plot_synteny as ps
    tracks, colours = _grid_tracks()          # every flanking gene in #4477aa
    res = sg.render(ps, tracks, colours, {}, _args(tmp_path, TREE2), sg.clade_muted_style())
    svg = sg.svg_document(res)
    muted = ps._lerp_hex("#4477aa", sg.MUTE_TOWARDS, sg.MUTE_SHARE)
    assert f'fill="{muted}"' in svg and 'fill="#4477aa"' not in svg
    assert not re.search(rf'<path [^>]*fill="{sg.FLANK_GREY}"', svg)     # no grey arrow
    assert res["info"]["goi_form"] == "panel"
    assert re.search(rf'<rect [^>]*fill="{stree.color_for_clade(1)}"', svg)


def test_gene_palette_is_large_distinct_and_free_of_the_goi_reds():
    """30 pale flanking-gene colours (Tol's light scheme, extended), no two closer than
    CIEDE2000 9, none near the GOI red, neighbours in the order genes take them clearly
    different, none dark (user, 2026-10-14: distinct, and pale so the GOI stands out)."""
    import plot_synteny as ps
    sys.path.insert(0, os.path.join(REPO, "scripts"))
    import make_gene_palette as mp
    pal = ps.GENE_PALETTE
    assert len(pal) >= 30 and len(set(pal)) == len(pal)
    labs = [mp.hex_to_lab(h) for h in pal]
    assert min(mp.de2000(a, b) for i, a in enumerate(labs) for b in labs[i + 1:]) >= 9
    assert min(mp.de2000(a, b) for a, b in zip(labs, labs[1:])) >= 11
    assert min(lab[0] for lab in labs) >= 64            # pale: CIELAB lightness
    assert pal[0] == mp.TOL_LIGHT[0] and set(pal[:7]) == set(mp.TOL_LIGHT)    # Tol's colours first
    assert pal == mp.make_palette(len(pal))             # the script reproduces the constant
    red = mp.hex_to_lab(ps.GOI_COLOUR)
    assert min(mp.de2000(a, red) for a in labs) > 30


def test_a_gene_holding_the_goi_in_an_intron_is_not_a_goi(monkeypatch):
    """Melittin lies in an intron of MESR3: the span overlap tied the two and the figures
    drew MESR3 as a second home melittin copy. The CDS owner of the hit decides."""
    import plot_synteny as ps
    monkeypatch.setattr(ps, "_GOI_NAMES", set())
    home = [{"chrom": "c", "start": 20000, "end": 37000, "name": "gene-HOST"},    # 17 kb host
            {"chrom": "c", "start": 22900, "end": 24400, "name": "gene-Melt"},    # in its intron
            {"chrom": "c", "start": 17500, "end": 19300, "name": "gene-NEXT"}]
    hit = [{"chrom": "c", "start": 23000, "end": 24350, "strand": "-"}]
    ps.identify_goi_names(home, hit)               # the span rule alone: both
    assert ps._GOI_NAMES == {"gene-Melt", "gene-HOST"}
    ps.set_goi_names(home, hit, ["gene-Melt"])     # with the resolved CDS owner: only that
    assert ps._GOI_NAMES == {"gene-Melt"}
    ps.set_goi_names(home + [{"chrom": "c", "start": 1, "end": 9, "name": "GOI_c_1"}], hit, ["gene-Melt"])
    assert ps._GOI_NAMES == {"gene-Melt", "GOI_c_1"}
    ps.set_goi_names(home, hit, [])                # nothing resolved: the span rule
    assert "gene-Melt" in ps._GOI_NAMES


def test_grid_variants_option():
    import argparse
    import plot_synteny as ps
    assert ps._grid_variant_list("all") == list(ps.GRID_VARIANTS)
    assert ps._grid_variant_list("compact, clades_muted") == ["compact", "clades_muted"]
    with pytest.raises(argparse.ArgumentTypeError):
        ps._grid_variant_list("compact,inline")


def test_clade_style_is_narrower_with_smaller_column_text():
    st, base = sg.clade_style(), sg.style_for()
    assert st.clade and not base.clade
    assert st.col_w < base.col_w and st.fs_col < base.fs_col
    assert "clade" in sg.SOURCE_COLUMNS_CLADES and "clade" not in sg.SOURCE_COLUMNS
