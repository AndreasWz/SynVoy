"""§1m collinearity-aware seed placement (docs/TODO.md §1m — the proper fix for
the flanking-gap "fumble", replacing the §17 hull rescue).

Covers the parts that were actually load-bearing in the design:
  - _longest_collinear_run: longest monotonic run over home-ranks (either direction).
  - _can_bridge: only bridge a large gap when the flanking collinearly continue the
    block's home-order run, with enough anchors and a bounded rank jump.
  - build_home_rank: rank flanking by first appearance, skip GOI proxies.
  - identify_synteny_blocks:
      * gap-bridging (#1): the decorin-like case — two clusters bracketing a gap
        (where the GOI sits) become ONE block spanning the gap.
      * collinear scoring (#3): an ordered block outranks a larger scrambled cluster.
      * STRICT no-op when home_rank is absent (legacy behaviour for existing callers).
"""
import os
import sys

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402


def _hit(query, chrom, start, end=None):
    return {"query": query, "chrom": chrom, "start": start,
            "end": end if end is not None else start + 10000}


# --------------------------------------------------------------------------- #
# _longest_collinear_run
# --------------------------------------------------------------------------- #
def test_collinear_run_increasing():
    assert isr._longest_collinear_run([0, 1, 2, 3]) == (4, "+")


def test_collinear_run_decreasing():
    length, direction = isr._longest_collinear_run([7, 6, 5, 4])
    assert length == 4 and direction == "-"


def test_collinear_run_scrambled():
    # [10, 3, 50, 1, 30]: longest non-increasing (reading as decreasing) = [10,3,1] = 3
    length, _ = isr._longest_collinear_run([10, 3, 50, 1, 30])
    assert length == 3


def test_collinear_run_ignores_none():
    # GOI proxies (None) are skipped, leaving [0,1,2].
    assert isr._longest_collinear_run([0, None, 1, None, 2]) == (3, "+")


def test_collinear_run_trivial():
    assert isr._longest_collinear_run([]) == (0, "+")
    assert isr._longest_collinear_run([5]) == (1, "+")


# --------------------------------------------------------------------------- #
# _can_bridge
# --------------------------------------------------------------------------- #
def _block(ranks_with_pos, home_rank):
    """Build a minimal block (list of loci) from (gene_id, start) pairs."""
    return [{"query": gid, "chrom": "chr1", "start": pos, "end": pos + 1000}
            for gid, pos in ranks_with_pos]


def test_can_bridge_forward_continuation():
    home_rank = {"g0": 0, "g1": 1, "g2": 2, "g3": 3, "g4": 4}
    block = _block([("g0", 0), ("g1", 100), ("g2", 200), ("g3", 300)], home_rank)
    cand = {"query": "g4", "chrom": "chr1", "start": 9999, "end": 10999}
    assert isr._can_bridge(block, cand, home_rank, max_rank_gap=5, min_anchors=3)


def test_can_bridge_rejects_too_few_anchors():
    home_rank = {"g0": 0, "g1": 1, "g4": 4}
    block = _block([("g0", 0), ("g1", 100)], home_rank)  # only 2 anchors
    cand = {"query": "g4", "chrom": "chr1", "start": 9999}
    assert not isr._can_bridge(block, cand, home_rank, max_rank_gap=5, min_anchors=3)


def test_can_bridge_rejects_unknown_candidate():
    home_rank = {"g0": 0, "g1": 1, "g2": 2}
    block = _block([("g0", 0), ("g1", 100), ("g2", 200)], home_rank)
    cand = {"query": "GOI_X", "chrom": "chr1", "start": 9999}  # no rank
    assert not isr._can_bridge(block, cand, home_rank, max_rank_gap=5, min_anchors=3)


def test_can_bridge_rejects_rank_jump_too_big():
    home_rank = {"g0": 0, "g1": 1, "g2": 2, "g3": 3, "gX": 20}
    block = _block([("g0", 0), ("g1", 100), ("g2", 200), ("g3", 300)], home_rank)
    cand = {"query": "gX", "chrom": "chr1", "start": 9999}  # rank 20, frontier 3
    assert not isr._can_bridge(block, cand, home_rank, max_rank_gap=5, min_anchors=3)


def test_can_bridge_rejects_direction_break():
    # Block is increasing; a candidate BELOW the frontier (would invert) is not bridged.
    home_rank = {"g2": 2, "g3": 3, "g4": 4, "g5": 5, "g1": 1}
    block = _block([("g2", 0), ("g3", 100), ("g4", 200), ("g5", 300)], home_rank)
    cand = {"query": "g1", "chrom": "chr1", "start": 9999}
    assert not isr._can_bridge(block, cand, home_rank, max_rank_gap=5, min_anchors=3)


def test_can_bridge_accepts_uniformly_inverted_block():
    """A wholly inverted neighbourhood bridges — only a direction REVERSAL is refused.

    The docstring used to claim "inversions are deliberately not bridged", which is
    wrong: _can_bridge takes the direction from the BLOCK's own anchors, so a block
    running g5,g4,g3 continues correctly into g2. Pinning it so the wrong reading
    cannot come back (docs corrected 2026-07-21).
    """
    home_rank = {"g1": 1, "g2": 2, "g3": 3, "g4": 4, "g5": 5}
    block = _block([("g5", 0), ("g4", 100), ("g3", 200)], home_rank)
    cand = {"query": "g2", "chrom": "chr1", "start": 9999, "end": 10999}
    assert isr._can_bridge(block, cand, home_rank, max_rank_gap=5, min_anchors=3)


def test_bridge_spans_gap_for_inverted_neighbourhood():
    """End-to-end counterpart: inverted anchors either side of a 2.8 Mb gap = ONE block."""
    home_rank = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}
    hits = [{"query": q, "chrom": "chr1", "start": s, "end": s + 5000,
             "strand": "+", "bits": 100}
            for q, s in [("F", 100_000), ("E", 150_000), ("D", 200_000),
                         ("C", 3_000_000), ("B", 3_050_000), ("A", 3_100_000)]]
    blocks = isr.identify_synteny_blocks(
        hits, max_intron=20_000, cluster_distance=150_000, home_rank=home_rank,
        bridge_max_gap=6_000_000, bridge_max_rank_gap=5, bridge_min_anchors=3)
    assert len(blocks) == 1
    assert blocks[0]["bridged"] is True
    assert blocks[0]["collinear_chain_len"] == 6
    assert blocks[0]["collinear_direction"] == "-"


# --------------------------------------------------------------------------- #
# build_home_rank
# --------------------------------------------------------------------------- #
def test_build_home_rank(tmp_path):
    fa = tmp_path / "initial_db.faa"
    fa.write_text(
        ">gene-A\nMAAAA\n"
        ">gene-B|exon_1\nMBBBB\n"   # same base as next, must not double-rank
        ">gene-B\nMBBBB\n"
        ">gene-C\nMCCCC\n"
        ">GOI_X\nMGGGG\n"           # GOI proxy: skipped
        ">GOI_X|exon_1\nMG\n"       # GOI exon: skipped
    )
    rank = isr.build_home_rank(str(fa))
    assert rank == {"gene-A": 0, "gene-B": 1, "gene-C": 2}
    assert "GOI_X" not in rank


# --------------------------------------------------------------------------- #
# identify_synteny_blocks — #1 gap-bridging
# --------------------------------------------------------------------------- #
def _decorin_like_hits():
    """Cluster A (ranks 0-3) at ~1.0-1.3 Mb, then a ~2.7 Mb gap (where the GOI sits),
    then cluster B (ranks 4-7) at ~4.0-4.3 Mb. All on chr1, collinear (0..7)."""
    return [
        _hit("g0", "chr1", 1_000_000), _hit("g1", "chr1", 1_100_000),
        _hit("g2", "chr1", 1_200_000), _hit("g3", "chr1", 1_300_000),
        _hit("g4", "chr1", 4_000_000), _hit("g5", "chr1", 4_100_000),
        _hit("g6", "chr1", 4_200_000), _hit("g7", "chr1", 4_300_000),
    ]


HOME_RANK_8 = {f"g{i}": i for i in range(8)}


def test_bridge_spans_gap_when_collinear():
    blocks = isr.identify_synteny_blocks(
        _decorin_like_hits(), cluster_distance=300_000,
        home_rank=HOME_RANK_8, bridge_max_gap=6_000_000,
        bridge_max_rank_gap=5, bridge_min_anchors=3,
    )
    assert len(blocks) == 1                      # A and B bridged into one block
    b = blocks[0]
    assert b["start"] == 1_000_000               # spans...
    assert b["end"] == 4_310_000                 # ...the whole neighbourhood
    assert b["start"] <= 2_650_000 <= b["end"]   # the GOI gap is now inside the window
    assert b["bridged"] is True
    assert b["collinear_chain_len"] == 8


def test_no_bridge_when_gap_exceeds_max():
    blocks = isr.identify_synteny_blocks(
        _decorin_like_hits(), cluster_distance=300_000,
        home_rank=HOME_RANK_8, bridge_max_gap=2_000_000,  # < 2.7 Mb gap
        bridge_max_rank_gap=5, bridge_min_anchors=3,
    )
    assert len(blocks) == 2                       # gap too large to bridge


def test_no_bridge_when_not_collinear():
    # Cluster B is a paralog cluster: ranks 20-23, not a continuation of A's 0-3.
    home_rank = {"g0": 0, "g1": 1, "g2": 2, "g3": 3,
                 "g4": 20, "g5": 21, "g6": 22, "g7": 23}
    blocks = isr.identify_synteny_blocks(
        _decorin_like_hits(), cluster_distance=300_000,
        home_rank=home_rank, bridge_max_gap=6_000_000,
        bridge_max_rank_gap=5, bridge_min_anchors=3,
    )
    assert len(blocks) == 2                        # not bridged (rank discontinuity)


# --------------------------------------------------------------------------- #
# identify_synteny_blocks — #3 collinear scoring
# --------------------------------------------------------------------------- #
def test_collinear_block_outranks_larger_scrambled_block():
    # Block X (chr1): 4 collinear genes (ranks 0-3) -> chain 4.
    # Block Y (chr2): 5 genes but scrambled ranks -> chain 3, genes_count 5.
    hits = [
        _hit("g0", "chr1", 1_000_000), _hit("g1", "chr1", 1_100_000),
        _hit("g2", "chr1", 1_200_000), _hit("g3", "chr1", 1_300_000),
        _hit("y0", "chr2", 1_000_000), _hit("y1", "chr2", 1_100_000),
        _hit("y2", "chr2", 1_200_000), _hit("y3", "chr2", 1_300_000),
        _hit("y4", "chr2", 1_400_000),
    ]
    home_rank = {"g0": 0, "g1": 1, "g2": 2, "g3": 3,
                 "y0": 10, "y1": 3, "y2": 50, "y3": 1, "y4": 30}
    blocks = isr.identify_synteny_blocks(
        hits, cluster_distance=300_000, home_rank=home_rank,
        bridge_max_gap=0,  # bridging off; this test is purely about scoring
    )
    assert blocks[0]["chrom"] == "chr1"           # collinear block ranked first...
    assert blocks[0]["collinear_chain_len"] == 4
    # ...despite the scrambled block having MORE genes.
    scrambled = next(b for b in blocks if b["chrom"] == "chr2")
    assert scrambled["genes_count"] == 5
    assert scrambled["collinear_chain_len"] == 3


# --------------------------------------------------------------------------- #
# Legacy no-op: without home_rank, behaviour is unchanged
# --------------------------------------------------------------------------- #
def test_legacy_no_home_rank_does_not_bridge():
    blocks = isr.identify_synteny_blocks(_decorin_like_hits(), cluster_distance=300_000)
    assert len(blocks) == 2                        # legacy: gap splits into two blocks
    for b in blocks:
        assert b["collinear_chain_len"] == 0       # scoring inert
        assert b["bridged"] is False


def test_legacy_sorts_by_gene_count():
    # Larger cluster must come first under legacy gene-count sort.
    hits = [
        _hit("a", "chr1", 1000), _hit("b", "chr1", 1100),          # 2 genes
        _hit("c", "chr1", 50000), _hit("d", "chr1", 50100), _hit("e", "chr1", 50200),  # 3
    ]
    blocks = isr.identify_synteny_blocks(hits, cluster_distance=500)
    assert blocks[0]["genes_count"] == 3
