"""Tests for collapse_to_one_leaf_per_species — the matrix-cladogram fix.

The matrix view aligns one tree branch per species row. Multi-copy GOI trees
(tandem paralogs, gene-family expansions) used to be drawn with every leaf
mapped to its species' row y, producing a tangled mess of degenerate
branches. The collapse helper prunes to one leaf per species and contracts
single-child internals so the cladogram has the same shape but one tip per
matrix row.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "bin"))

from synvoy_tree import (  # noqa: E402
    collapse_to_one_leaf_per_species,
    parse_newick_tree,
    species_from_leaf,
)


def _leaf_names(root):
    return [l.name for l in root.leaves()]


def _max_depth(root):
    best = 0.0
    def _walk(n, d):
        nonlocal best
        if n.is_leaf():
            best = max(best, d)
            return
        for c in n.children:
            _walk(c, d + (c.dist or 0.0))
    _walk(root, 0.0)
    return best


def test_collapse_drops_duplicate_species_keeping_first():
    # Apis_florea appears 3× (extra1/extra2/extra3); only the first should survive.
    text = ("(((GOI_Melt|Apis_florea_fna_b0_l1_extra1:0.1,"
            "GOI_Melt|Apis_florea_fna_b0_l1_extra2:0.2):0.05,"
            "GOI_Melt|Bombus_terrestris_fna_b0_l1_fallback:0.3):0.1,"
            "GOI_Melt|Apis_florea_fna_b0_l1_extra3:0.4);")
    root = parse_newick_tree(text)
    collapsed = collapse_to_one_leaf_per_species(root, species_from_leaf)
    names = _leaf_names(collapsed)
    # Two distinct species → two leaves.
    assert len(names) == 2
    # First Apis_florea leaf in left-to-right order is _extra1.
    assert any("extra1" in n for n in names)
    assert all("extra2" not in n and "extra3" not in n for n in names)
    assert any("Bombus_terrestris" in n for n in names)


def test_collapse_drops_unmappable_leaves():
    # GOI_Melt (no |species suffix) cannot be mapped to a row.
    text = ("(GOI_Melt:0.0,"
            "(GOI_Melt|Apis_florea_fna_b0_l1_exon_ann:0.2,"
            "GOI_Melt|Bombus_terrestris_fna_b0_l1_fallback:0.3):0.1);")
    root = parse_newick_tree(text)
    collapsed = collapse_to_one_leaf_per_species(root, species_from_leaf)
    names = _leaf_names(collapsed)
    assert len(names) == 2
    assert all(name and "|" in name for name in names)


def test_collapse_contracts_single_child_internals_preserving_depth():
    # After dropping all but one Apis_florea, the immediate parent has only
    # one surviving child — that node must be collapsed and its branch length
    # absorbed so the leaf's depth from the root stays the same.
    text = ("((GOI_Melt|Apis_florea_fna_b0_l1_extra1:0.1,"
            "GOI_Melt|Apis_florea_fna_b0_l1_extra2:0.1):0.5,"
            "GOI_Melt|Bombus_terrestris_fna_b0_l1_fallback:0.7);")
    root = parse_newick_tree(text)
    apis_depth_before = 0.5 + 0.1  # parent edge + leaf edge
    collapsed = collapse_to_one_leaf_per_species(root, species_from_leaf)
    # Only two leaves left, attached directly to the new root.
    assert len(collapsed.children) == 2
    apis_leaf = next(l for l in collapsed.leaves() if "Apis_florea" in l.name)
    # Single-child contraction: leaf.dist now includes the absorbed parent edge.
    assert apis_leaf.dist == pytest.approx(apis_depth_before)
    # Total tree depth unchanged for this leaf.
    assert _max_depth(collapsed) == pytest.approx(0.7)


def test_collapse_no_op_on_single_copy_tree():
    # If every species already has exactly one leaf, topology is unchanged.
    text = ("(GOI_Melt|Apis_florea_fna_b0_l1_exon_ann:0.1,"
            "(GOI_Melt|Bombus_terrestris_fna_b0_l1_fallback:0.2,"
            "GOI_Melt|Solenopsis_invicta_fna_b1_l1_fallback:0.3):0.05);")
    root = parse_newick_tree(text)
    collapsed = collapse_to_one_leaf_per_species(root, species_from_leaf)
    assert sorted(_leaf_names(collapsed)) == sorted(_leaf_names(root))


def test_collapse_returns_none_when_everything_pruned():
    # Tree of only home-reference leaves with no species mapping.
    text = "(GOI_Melt:0.1,GOI_Melt:0.2);"
    root = parse_newick_tree(text)
    assert collapse_to_one_leaf_per_species(root, species_from_leaf) is None


def test_collapse_realistic_melittin_tree_one_leaf_per_species():
    # Tandem-rich tree similar to what melittin produces: many copies
    # per species. After collapsing, every surviving leaf maps to a unique
    # species.
    text = (
        "((((GOI_copy_2|Colletes_gigas_fa_b0_l1:0.5,"
        "GOI_copy_2|Megachile_rotundata_fna_b0_l1:1.3):1.1,"
        "GOI_copy_2|Cardiocondyla_obscurior_fna_b0_l1:0.2):0.1,"
        "(GOI_Melt|Nomia_melanderi_fna_b0_l2_fallback:0.4,"
        "GOI_copy_3|Tetramorium_bicarinatum_fna_b0_l1:0.3):0.3):0.2,"
        "(GOI_copy_1|Cardiocondyla_obscurior_fna_b0_l1:0.5,"
        "GOI_Melt|Apis_florea_fna_b0_l1_exon_ann:0.0):0.1);"
    )
    root = parse_newick_tree(text)
    collapsed = collapse_to_one_leaf_per_species(root, species_from_leaf)
    species = [species_from_leaf(l.name) for l in collapsed.leaves()]
    assert len(species) == len(set(species))
    assert set(species) == {
        "Colletes_gigas", "Megachile_rotundata", "Cardiocondyla_obscurior",
        "Nomia_melanderi", "Tetramorium_bicarinatum", "Apis_florea",
    }
