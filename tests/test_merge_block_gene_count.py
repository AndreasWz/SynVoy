"""A merged block counts each flanking gene once.

merge_synteny_blocks() joins blocks whose padded spans overlap. It used to ADD their
gene counts and keep the first part's gene list. Two loci of one multi-copy flanking
gene less than 2 x padding apart therefore read as a 2-gene block and passed
``min_block_genes = 2``, the rule that a neighbourhood needs two different anchors.
Replaying block building on the saved genome-wide hits of two local runs (firefly
luciferase against two flies, yeast STE2; 2026-10-06): 89 of the 544 blocks that passed
the filter held a single distinct flanking gene, and each of them got a full window
search. The gene list of the later parts was also lost, so the "block seeded by the
GOI's own home gene" exception and the per-genome cap never saw those genes.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402


def _block(start, end, genes, chrom="chr1", **kw):
    b = {"chrom": chrom, "start": start, "end": end, "genes": list(genes),
         "genes_count": len(genes), "loci_count": len(genes),
         "collinear_chain_len": 1, "collinear_direction": "+", "bridged": False}
    b.update(kw)
    return b


def _hit(query, start, end, chrom="chr1"):
    return {"query": query, "chrom": chrom, "start": start, "end": end, "strand": "+",
            "pident": 60.0, "alnlen": 100, "evalue": 1e-20, "bits": 120.0,
            "qstart": 1, "qend": 100}


def test_two_copies_of_one_gene_are_one_gene():
    merged = isr.merge_synteny_blocks(
        [_block(100_000, 101_000, ["gene-A"]), _block(300_000, 301_000, ["gene-A"])],
        padding=150_000)
    assert len(merged) == 1
    assert merged[0]["genes"] == ["gene-A"]
    assert merged[0]["genes_count"] == 1


def test_different_genes_are_united():
    merged = isr.merge_synteny_blocks(
        [_block(100_000, 101_000, ["gene-A", "gene-B"]),
         _block(300_000, 301_000, ["gene-B", "gene-C"])],
        padding=150_000)
    assert len(merged) == 1
    assert merged[0]["genes"] == ["gene-A", "gene-B", "gene-C"]
    assert merged[0]["genes_count"] == 3


def test_gene_list_of_every_part_survives_a_chain_of_merges():
    merged = isr.merge_synteny_blocks(
        [_block(100_000, 101_000, ["gene-A"]), _block(250_000, 251_000, ["gene-A"]),
         _block(400_000, 401_000, ["gene-LY6E"])],
        padding=150_000)
    assert len(merged) == 1
    # The later part's gene is what the "seeded by the GOI's own home gene" check reads.
    assert "gene-LY6E" in merged[0]["genes"]
    assert merged[0]["genes_count"] == 2


def test_blocks_without_a_gene_list_keep_the_summed_count():
    merged = isr.merge_synteny_blocks(
        [{"chrom": "chr1", "start": 1000, "end": 2000, "genes_count": 2},
         {"chrom": "chr1", "start": 2100, "end": 3000, "genes_count": 3}], padding=200)
    assert merged[0]["genes_count"] == 5


def test_unmerged_blocks_are_untouched():
    blocks = [_block(100_000, 101_000, ["gene-A"]),
              _block(900_000, 901_000, ["gene-A"]),
              _block(100_000, 101_000, ["gene-A"], chrom="chr2")]
    merged = isr.merge_synteny_blocks(blocks, padding=150_000)
    assert len(merged) == 3
    assert all(b["genes_count"] == 1 for b in merged)


def test_a_duplicated_flanking_gene_does_not_pass_the_two_gene_filter():
    # Two loci of gene-A 200 kb apart: farther than cluster_distance (two blocks),
    # closer than 2 x region_padding (merged again). Before the fix the merged block
    # reported genes_count 2.
    hits = [_hit("gene-A", 100_000, 100_900), _hit("gene-A", 300_000, 300_900),
            _hit("gene-B", 5_000_000, 5_000_900), _hit("gene-C", 5_020_000, 5_020_900)]
    blocks = isr.identify_synteny_blocks(hits, max_intron=20_000, cluster_distance=150_000)
    assert len(blocks) == 3
    merged = isr.merge_synteny_blocks(blocks, 150_000)
    by_start = {b["start"]: b for b in merged}
    assert len(merged) == 2
    assert by_start[100_000]["genes_count"] == 1        # one gene, two copies
    assert by_start[5_000_000]["genes_count"] == 2      # a real two-gene neighbourhood
    kept = [b for b in merged if b["genes_count"] >= 2]
    assert [b["genes"] for b in kept] == [["gene-B", "gene-C"]]
