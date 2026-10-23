"""The dead-zone fix (profile handoff, constraint C2): the dispersed-GOI rescue leaves to the
blocks only the GOI hits within the GOI margin of an anchor, the margin a block keeps GOI hits
in (iterative_search_runner.select_dispersed_goi_seeds, skip_radius). Before, the radius was
the window margin, and a strong GOI hit between the two margins was used by neither path.
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
# window margin 60 kb, GOI margin 20 kb; anchor gene-a1 spans 1,200,000-1,201,000
NEAR = _goi("c1", 1_210_000)          # 9 kb from it: inside the GOI margin
BETWEEN = _goi("c1", 1_241_000)       # 40 kb from it: beyond the GOI margin, inside the window margin
FAR = _goi("c1", 1_500_000)           # 99 kb from the nearest anchor


def _seeds(hits, **kw):
    drops, near = {}, []
    seeds = isr.select_dispersed_goi_seeds(hits, ANCHORS, region_padding=60_000, envelope_margin=500_000,
                                           drops=drops, near_out=near, **kw)
    return [h["start"] for h in seeds], drops, [h["start"] for h in near]


def test_a_hit_between_the_goi_margin_and_the_window_margin_is_rescued():
    starts, drops, near = _seeds([NEAR, BETWEEN, FAR], skip_radius=20_000)
    assert starts == [1_241_000, 1_500_000]
    assert drops == {"dispersed_near_anchor": 1} and near == [1_210_000]


def test_without_a_skip_radius_the_window_margin_decides_as_before():
    starts, drops, near = _seeds([NEAR, BETWEEN, FAR])
    assert starts == [1_500_000]
    assert drops == {"dispersed_near_anchor": 2} and near == [1_210_000, 1_241_000]


def test_the_search_passes_the_margin_of_its_goi_core():
    """The skip radius and the GOI core (block +- core_margin in process_region_block) are the
    same expression, so no hit lies between them."""
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    core = re.search(r"core_margin = (.+)\n", src).group(1).strip()
    assert core == "getattr(args, 'goi_margin', None) or args.gap_search_window"
    assert f"skip_radius={core})" in src
