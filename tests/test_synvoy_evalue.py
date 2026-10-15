"""bin/synvoy_evalue.py: the block E-value E_syn, places and the per-candidate E-value of
the label rework (tmp/formula_check_2026-10-14/FINAL_for_vscode.md, section 3.2)."""
import itertools
import math
import os
import random
import sys

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import synvoy_evalue as ev  # noqa: E402


def _monotone_run(perm):
    def lis(seq):
        best = []
        for x in seq:
            i = next((j for j, b in enumerate(best) if b >= x), len(best))
            best[i:i + 1] = [x]
        return len(best)
    return max(lis(perm), lis([-x for x in perm]))


@pytest.mark.parametrize("n", range(1, 8))
def test_monotone_run_probability_is_exact(n):
    """Against brute force over every permutation of n <= 7."""
    perms = list(itertools.permutations(range(n)))
    for k in range(1, n + 2):
        brute = sum(_monotone_run(p) >= k for p in perms) / len(perms)
        assert ev.p_monotone_run(k, n) == pytest.approx(brute, abs=1e-12), (k, n)


def test_erdos_szekeres():
    """Any sequence of more than (k-1)^2 distinct values has a monotone run of k."""
    assert ev.p_monotone_run(4, 10) == 1.0
    assert ev.p_monotone_run(4, 9) < 1.0


def test_poisson_binomial():
    d = ev.poisson_binomial([0.5, 0.5])
    assert d == pytest.approx([0.25, 0.5, 0.25])
    assert sum(ev.poisson_binomial([0.1, 0.7, 0.3, 0.9])) == pytest.approx(1.0)


def test_e_syn_falls_with_the_run_and_rises_with_places():
    base = dict(window_len=50_000, genome_len=14_000_000)
    e = [ev.e_syn(k, [1] * 20, **base) for k in (2, 3, 4, 6, 10)]
    assert all(a > b for a, b in zip(e, e[1:]))
    assert ev.e_syn(4, [30] * 20, **base) > ev.e_syn(4, [1] * 20, **base)
    assert ev.e_syn(0, [1] * 20, **base) == pytest.approx(280.0)


def test_e_syn_matches_a_simulation():
    """Random windows of a random genome: the expected count of windows with a run >= k."""
    rng = random.Random(5)
    G, W, places, k = 2_000_000, 100_000, [2] * 12, 3
    hits = []
    for gene, m in enumerate(places):
        for _ in range(m):
            hits.append((rng.randrange(G), gene))
    hits.sort()
    windows = 0
    for start in range(0, G - W + 1, W):        # non-overlapping windows
        genes = []
        for pos, gene in hits:
            if start <= pos < start + W and gene not in genes:
                genes.append(gene)
        windows += _monotone_run(genes) >= k if genes else 0
    # one draw of a 20-window genome: compare with the expectation over draws instead
    expect = ev.e_syn(k, places, W, G)
    sims = []
    for rep in range(300):
        hits = sorted((rng.randrange(G), g) for g, m in enumerate(places) for _ in range(m))
        n = 0
        for start in range(0, G - W + 1, W):
            genes = []
            for pos, gene in hits:
                if start <= pos < start + W and gene not in genes:
                    genes.append(gene)
            n += _monotone_run(genes) >= k if genes else 0
        sims.append(n)
    assert sum(sims) / len(sims) == pytest.approx(expect, rel=0.25)


def test_places_per_gene_merges_nearby_hits():
    hits = [{"query": "gene-a|exon_1", "chrom": "c1", "start": 100, "end": 400},
            {"query": "gene-a|exon_2", "chrom": "c1", "start": 3000, "end": 3300},
            {"query": "gene-a", "chrom": "c1", "start": 90000, "end": 90300},
            {"query": "gene-a", "chrom": "c2", "start": 100, "end": 300},
            {"query": "gene-b", "chrom": "c1", "start": 5000, "end": 5300}]
    assert ev.places_per_gene(hits) == {"gene-a": 3, "gene-b": 1}


def test_slot_between_the_nearest_neighbours_plus_a_pitch():
    placed = {-2: (1000, 2000), -1: (5000, 6000), 1: (9000, 9500), 3: (20000, 21000)}
    assert ev.slot_span(placed, 1000) == (5000, 10000)
    inverted = {-1: (9000, 9500), 1: (5000, 6000)}
    assert ev.slot_span(inverted, 0) == (6000, 9000)
    assert ev.slot_span({-1: (5000, 6000)}, 1000) is None      # one side only


def test_place_is_the_smallest_that_holds_the_candidate():
    slot, core, window = (5000, 10000), (3000, 30000), (0, 60000)
    assert ev.place_of((6000, 7000), True, slot, core, window, 10**7) == ("slot", 5000)
    assert ev.place_of((9000, 12000), True, slot, core, window, 10**7) == ("core", 27000)
    assert ev.place_of((40000, 41000), True, slot, core, window, 10**7) == ("window", 60000)
    assert ev.place_of((6000, 7000), False, slot, core, window, 10**7) == ("genome", 10**7)
    assert ev.place_of((6000, 7000), True, None, core, window, 10**7) == ("core", 27000)


def test_search_space_is_one_frame_of_the_place():
    assert ev.search_space(100, 30_000) == 100 * 10_000


def test_split_at_stops_keeps_pieces_apart():
    assert ev.split_at_stops("MKT*AAAAAAAA*BBBBB**CCCCCCC") == ["AAAAAAAA", "BBBBB", "CCCCCCC"]


def test_blastp_evalue_is_linear_in_the_search_space():
    if not ev._blastp():
        pytest.skip("blastp not installed")
    rng = random.Random(3)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    q = "".join(rng.choice(aa) for _ in range(250))
    x = list(q)
    for i in rng.sample(range(250), 150):
        x[i] = rng.choice(aa)
    x = "".join(x)
    e6, b6 = ev.blastp_evalue(q, [x], 10 ** 6)
    e8, b8 = ev.blastp_evalue(q, [x], 10 ** 8)
    assert b6 == b8 and e8 == pytest.approx(100 * e6, rel=1e-3)
    junk = "".join(rng.choice(aa) for _ in range(250))
    e_min, _b = ev.blastp_evalue(q, [junk, x], 10 ** 6)
    assert e_min == pytest.approx(e6, rel=1e-6)            # the best subject counts
    assert ev.blastp_evalue(q, [], 10 ** 6) == (math.inf, 0.0)
