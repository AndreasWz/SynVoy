"""The clade version of the print anchor grid (synvoy_grid.clade_style / goi_clades):
GOI copies coloured by their clade in the GOI gene tree, written only when it shows more than
the species tree beside the rows."""
import os
import sys
import types

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


def test_clade_style_is_narrower_with_smaller_column_text():
    st, base = sg.clade_style(), sg.style_for()
    assert st.clade and not base.clade
    assert st.col_w < base.col_w and st.fs_col < base.fs_col
    assert "clade" in sg.SOURCE_COLUMNS_CLADES and "clade" not in sg.SOURCE_COLUMNS
