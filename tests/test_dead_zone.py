"""The dead zone (profile handoff, constraint C2): the dispersed-GOI
rescue leaves to a block only the GOI hits inside its GOI zone (block +- the GOI margin, the core a
block keeps GOI hits in), and it decides that against the blocks themselves, built without it
(iterative_search_runner.select_dispersed_goi_seeds, goi_zones). Until 2026-11-16 it skipped every
hit within the GOI margin of an anchor, assuming that anchor formed a block; a lone anchor forms
none, and a strong hit beside it was searched by nobody.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import iterative_search_runner as isr  # noqa: E402


def _goi(chrom, start):
    return {"query": "GOI_q", "chrom": chrom, "start": start, "end": start + 300,
            "pident": 60.0, "evalue": 1e-20, "alnlen": 80}


def _anchor(gene, chrom, start):
    return {"query": gene, "chrom": chrom, "start": start, "end": start + 1000,
            "pident": 90.0, "evalue": 1e-50, "alnlen": 300, "bits": 500.0}


ANCHORS = [_anchor(f"gene-a{i}", "c1", 1_000_000 + i * 200_000) for i in range(4)]
# window margin 60 kb, GOI margin 20 kb; anchor gene-a1 spans 1,200,000-1,201,000 and forms a
# block of its own here: its GOI zone is 1,180,000-1,221,000
ZONES = isr.goi_zones([{"chrom": "c1", "start": 1_200_000, "end": 1_201_000}], 20_000)
NEAR = _goi("c1", 1_210_000)          # 9 kb from it: inside the GOI zone
BETWEEN = _goi("c1", 1_241_000)       # 40 kb from it: beyond the GOI zone, inside the window margin
FAR = _goi("c1", 1_500_000)           # 99 kb from the nearest anchor


def _seeds(hits, **kw):
    drops, in_zone = {}, []
    seeds = isr.select_dispersed_goi_seeds(hits, ANCHORS, region_padding=60_000, envelope_margin=500_000,
                                           drops=drops, in_zone_out=in_zone, **kw)
    return [h["start"] for h in seeds], drops, [h["start"] for h in in_zone]


def test_a_hit_inside_a_goi_zone_is_left_to_its_block():
    starts, drops, in_zone = _seeds([NEAR, BETWEEN, FAR], goi_zones=ZONES)
    assert starts == [1_241_000, 1_500_000]
    assert drops == {"dispersed_in_goi_zone": 1} and in_zone == [1_210_000]


def test_a_hit_beside_an_anchor_without_a_block_is_a_seed():
    """No block holds gene-a1 (no zone): the hit 9 kb from it is a seed (rule envelope: c1 carries
    four anchors). Before, it was skipped as 'next to an anchor' and lay in no searched zone."""
    drops, in_zone = {}, []
    seeds = isr.select_dispersed_goi_seeds([NEAR, BETWEEN], ANCHORS, region_padding=60_000,
                                           envelope_margin=500_000, drops=drops, in_zone_out=in_zone)
    assert [(s["start"], s["rescue_rule"], s["rescue_anchors"]) for s in seeds] == [(1_210_000, "envelope", 4)]
    assert seeds[0]["rescue_anchor_bp"] == 9_000
    assert drops == {} and in_zone == []          # the two hits share one bin (60 kb): one seed


def test_the_search_uses_one_goi_margin_for_the_core_the_zones_and_the_near_radius():
    """The GOI core of process_region_block, the zones the rescue leaves to the blocks and the
    near-neighbour radius are the same number (goi_core_margin), so no hit lies between them."""
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    assert re.search(r"def goi_core_margin\(args\) -> int:\n(?:    .*\n)*?"
                     r"    return int\(getattr\(args, 'goi_margin', None\) or args\.gap_search_window\)", src)
    region = src.split("def process_region_block(", 1)[1].split("\ndef ", 1)[0]
    assert "core_margin = goi_core_margin(args)" in region
    selector = src.split("def rescue_selector(", 1)[1].split("\ndef ", 1)[0]
    assert "margin = goi_core_margin(args)" in selector
    assert "zones = goi_zones(blocks, margin)" in selector and "near_radius=margin" in selector
    axis = src.split("def genome_axis_blocks(", 1)[1].split("\ndef ", 1)[0]
    assert "goi_core_margin(args)" in axis
    assert "skip_radius" not in src
