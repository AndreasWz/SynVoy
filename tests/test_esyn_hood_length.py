"""Review fix RF-1b (batch of 2026-10-15): E_syn gets the length its run was counted in.

The collinear run k of a block is counted over the neighbourhood (block +- padding); the
DNA searched can be smaller (derived distances: block +- search margin). E_syn is computed
for the neighbourhood: equal to before where the two lengths are equal (yeast, honeybee),
larger where the searched part is smaller (large-pitch loci).
"""
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402


def _args():
    ranks = {f"g{i}": i for i in range(1, 11)}
    return SimpleNamespace(genome_ctx={"genome_len": 100_000_000, "places": {g: 2 for g in ranks}},
                           home_rank=ranks, goi_ranks={}, home_geometry={"pitch": 30_000}, expansion=1.0)


def _esyn(k, w, hood):
    a = _args()
    terms = isr.block_place_terms(a, 0, 1_000_000, 1_000_000 + w, 1_000_000, 1_000_000 + w, k, [],
                                  hood_len=hood)
    return terms["esyn"], terms["window_len"]


def test_equal_lengths_give_the_old_value():
    for k in (2, 3, 4):
        new, wl = _esyn(k, 1_600_000, 1_600_000)
        old = synvoy_evalue.e_syn(k, [2] * 10, 1_600_000, 100_000_000)
        assert new == old and wl == 1_600_000


def test_a_smaller_searched_window_no_longer_shrinks_esyn():
    for k in (2, 4):
        searched, _ = _esyn(k, 1_600_000, None)            # the old computation
        hood, wl = _esyn(k, 1_600_000, 4_800_000)
        assert hood > searched and wl == 1_600_000         # tblastn rescaling keeps the window
        assert hood == synvoy_evalue.e_syn(k, [2] * 10, 4_800_000, 100_000_000)


def test_the_runner_passes_the_neighbourhood():
    src = open(os.path.join(os.path.dirname(__file__), "..", "bin", "iterative_search_runner.py")).read()
    assert "hood_len=hood_end - hood_start," in src
