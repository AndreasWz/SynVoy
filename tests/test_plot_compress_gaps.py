#!/usr/bin/env python3
"""Regression: a gene must never be drawn over a compressed gap.

`compress_track_coordinates` shrinks large *intergenic* gaps. The gap used to
be measured from the immediately-previous gene's end, but genes overlap and
nest (a long readthrough/lncRNA can span several short neighbours). Measuring
from the previous gene let the gap *after a short nested gene* be compressed
while the enclosing long gene still extended across it — so the long gene was
drawn right over the break, implying it was ~gap-bp longer than it really is.
The fix measures gaps from the running furthest extent (`max_end`).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402


def _g(start, end, chrom="c1", name=None, strand="+"):
    return {"chrom": chrom, "start": start, "end": end,
            "name": name or f"g{start}", "strand": strand,
            "home_gene_id": name or f"g{start}", "identity": 90.0}


def _no_gene_spans_any_break(compressed, breaks):
    for brk in breaks:
        bx = brk["x"]
        for g in compressed:
            # A break must sit in empty space, not inside a gene's drawn body.
            if g["start_plot"] < bx < g["end_plot"]:
                return False, (brk, g)
    return True, None


def test_long_gene_nesting_short_neighbour_then_gap():
    # A short gene, a long readthrough that spans it, a gene nested inside the
    # long one, then a far gene that opens a real gap *after* the long gene.
    genes = [
        _g(1_000, 2_000, name="short_A"),
        _g(1_500, 50_000, name="long_readthrough"),
        _g(3_000, 3_500, name="nested_C"),
        _g(80_000, 81_000, name="far_D"),
    ]
    compressed, breaks = ps.compress_track_coordinates(
        genes, threshold=10_000, visual_gap=2_000)
    ok, offender = _no_gene_spans_any_break(compressed, breaks)
    assert ok, f"a gene is drawn over a compressed gap: {offender}"
    # The real gap (50 kb -> 80 kb) is still compressed: exactly one break.
    assert len([b for b in breaks if not b.get("is_chrom_break")]) == 1
    # The long readthrough keeps its full width (it was never split).
    lr = next(g for g in compressed if g["name"] == "long_readthrough")
    assert lr["end_plot"] - lr["start_plot"] == 50_000 - 1_500


def test_simple_intergenic_gap_still_compressed():
    genes = [_g(1_000, 2_000, name="A"), _g(500_000, 501_000, name="B")]
    compressed, breaks = ps.compress_track_coordinates(
        genes, threshold=10_000, visual_gap=2_000)
    gaps = [b for b in breaks if not b.get("is_chrom_break")]
    assert len(gaps) == 1 and gaps[0]["text"].endswith("kb")
    # B is pulled left to ~visual_gap past A, not left at its raw 500 kb.
    b = next(g for g in compressed if g["name"] == "B")
    assert b["start_plot"] < 10_000
    ok, _ = _no_gene_spans_any_break(compressed, breaks)
    assert ok


def test_no_break_when_overlapping_genes_have_no_real_gap():
    # Dense, overlapping genes with no gap above threshold -> no gap breaks.
    genes = [_g(1_000, 9_000, name="A"), _g(2_000, 3_000, name="B"),
             _g(4_000, 12_000, name="C"), _g(11_000, 13_000, name="D")]
    compressed, breaks = ps.compress_track_coordinates(
        genes, threshold=10_000, visual_gap=2_000)
    assert [b for b in breaks if not b.get("is_chrom_break")] == []
    ok, offender = _no_gene_spans_any_break(compressed, breaks)
    assert ok, offender


def test_chromosome_break_not_inside_a_gene():
    genes = [_g(1_000, 2_000, chrom="c1", name="A"),
             _g(1_500, 40_000, chrom="c1", name="long"),
             _g(5_000, 6_000, chrom="c2", name="B")]
    compressed, breaks = ps.compress_track_coordinates(
        genes, threshold=10_000, visual_gap=2_000)
    assert any(b.get("is_chrom_break") for b in breaks)
    ok, offender = _no_gene_spans_any_break(compressed, breaks)
    assert ok, offender
