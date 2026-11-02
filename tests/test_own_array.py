"""Fix-phase item 10: own-array specificity and ownership without a margin.

HomeProteome.best_other excluded only q's own gene, so in a tandem array a sibling copy could be the best other home
protein: with a copy identical to q every significant target copy was AMBIGUOUS not_specific, at a margin of at most 0
bits. A home protein in q's own home locus that
passes the own-gene rule against q (>= 95 % identity over >= 90 % of q) is now no competitor for a model x, unless x favours
it over q by the seed's bar of 10 bits: at 450 aa a 96 % tandem neighbour is told apart by tens of bits. Real paralogs in the
neighbourhood, below that identity, stay competitors. The rule and its values exist already: no new value.
"""
import os
import random
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
pytest.importorskip("parasail")
import synvoy_specificity as ss  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(n, seed):
    rng = random.Random(seed)
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _mutate(seq, k, seed):
    """``seq`` with ``k`` substitutions at distinct positions."""
    rng = random.Random(seed)
    out = list(seq)
    for i in rng.sample(range(1, len(seq)), k):
        out[i] = rng.choice([a for a in AA if a != out[i]])
    return "".join(out)


Q = _protein(120, 1)


def _home(genes):
    """HomeProteome of (id, seq, chrom, start) with 1-kb genes."""
    return ss.HomeProteome([(g, s) for g, s, _c, _p in genes], {g: (c, p, p + 999) for g, _s, c, p in genes})


def _array(*copies, other=()):
    genes = [("gene-q", Q, "c1", 10_000)]
    for i, seq in enumerate(copies):
        genes.append((f"gene-a{i}", seq, "c1", 12_000 + 2_000 * i))
    genes.extend(other)
    # A far protein, so the proteome is never only the array.
    genes.append(("gene-far", _protein(200, 77), "c5", 1_000))
    return _home(genes)


def _spec(home, x, n=5, own="gene-q"):
    bits = ss.raw_to_bits(home.scores(x)[0][0])
    return ss.specificity(home, "GOI_q", Q, x, bits, None, n, own_gene_id=own)


@pytest.mark.parametrize("diffs,passes", [(0, True), (1, True), (5, True), (6, True), (7, False), (8, False), (20, False)])
def test_the_own_gene_rule_at_120_aa(diffs, passes):
    """95 % of 120 aa: up to 6 differences. Copies with 8-20 differences are real competitors."""
    assert ss.passes_own_gene_rule(Q, _mutate(Q, diffs, 3)) is passes


def test_a_short_copy_does_not_pass():
    assert not ss.passes_own_gene_rule(Q, Q[:100])                   # 83 % of q
    assert ss.passes_own_gene_rule(Q, Q[:115])                       # 96 % of q


def test_an_identical_copy_in_the_array_is_no_competitor():
    home = _array(Q)
    spec = _spec(home, _mutate(Q, 10, 4))
    assert spec["best_other"] != "gene-a0"
    assert not spec["real_competitor"] or spec["margin"] > 50


def test_a_near_identical_copy_is_no_competitor():
    home = _array(_mutate(Q, 3, 5))
    assert _spec(home, _mutate(Q, 10, 4))["best_other"] != "gene-a0"


def _q_bits(home, q, x):
    import parasail
    return ss.raw_to_bits(float(parasail.sw_striped_sat(x, q, 12, 1, parasail.blosum62).score))


def test_a_target_copy_closer_to_a_near_identical_neighbour_is_told_apart():
    """The real shape: tubulin beta, 444 aa against a 96.4 % tandem neighbour. x, the target's copy of the neighbour, is
    tens of bits closer to it than to q: the neighbour stays p*, and x stays not_specific."""
    long_q = _protein(450, 21)
    neighbour = _mutate(long_q, 16, 22)              # 96.4 %
    assert ss.passes_own_gene_rule(long_q, neighbour)
    home = _home([("gene-q", long_q, "c1", 10_000), ("gene-n", neighbour, "c1", 12_000),
                  ("gene-far", _protein(200, 77), "c5", 1_000)])
    x = neighbour                                     # the neighbour's ortholog, identical to it
    bits = _q_bits(home, long_q, x)
    spec = ss.specificity(home, "GOI_q", long_q, x, bits, None, 5, own_gene_id="gene-q")
    assert spec["best_other"] == "gene-n" and spec["margin"] <= -ss.SEED_MIN_MARGIN_BITS and spec["real_competitor"]
    # A target copy equally far from both is not told apart: the neighbour is left out.
    between = long_q[:225] + neighbour[225:]
    bits = _q_bits(home, long_q, between)
    assert ss.specificity(home, "GOI_q", long_q, between, bits, None, 5, own_gene_id="gene-q")["best_other"] != "gene-n"


def test_a_copy_outside_n_neighbours_but_inside_the_hull_is_no_competitor():
    """The flanking walk skips relatives, so this locus's hull spans an array that n genes on each side do not."""
    copies = [_protein(120, 100 + i) for i in range(4)] + [Q]
    home = _array(*copies)
    x = _mutate(Q, 10, 4)
    bits = ss.raw_to_bits(home.scores(x)[0][0])
    assert ss.specificity(home, "GOI_q", Q, x, bits, None, 2, own_gene_id="gene-q")["best_other"] == "gene-a4"
    hull = {"gene-q", "gene-a0", "gene-a1", "gene-a2", "gene-a3", "gene-a4"}
    assert ss.specificity(home, "GOI_q", Q, x, bits, hull, 2, own_gene_id="gene-q")["best_other"] != "gene-a4"


def test_a_hull_without_qs_own_gene_adds_nothing():
    """A run from a secondary home locus: a copy in that hull is not in q's own home locus."""
    home = _array(other=[("gene-far-copy", Q, "c7", 5_000)])
    x = _mutate(Q, 10, 4)
    bits = ss.raw_to_bits(home.scores(x)[0][0])
    assert ss.specificity(home, "GOI_q", Q, x, bits, {"gene-far-copy"}, 2, own_gene_id="gene-q")["best_other"] == "gene-far-copy"


def test_the_rule_also_reads_qs_own_gene():
    """The hull rescue's q is the query file: 92 % of the home protein. A copy of the home gene is an own copy all the same."""
    variant = _mutate(Q, 10, 31)
    home = _array(Q)
    assert not ss.passes_own_gene_rule(variant, Q)
    assert home.own_copies("GOI_v", variant, 5, "gene-q") == {"gene-q", "gene-a0"}


def test_a_model_closer_to_the_home_gene_than_to_q_is_not_told_apart_from_the_gene():
    """q is another species' protein, 75 % of the home gene: x, close to the home gene and its copy, is tens of bits above
    q for both. That says nothing about the copy against the gene: x does not favour the copy, so the copy stays out."""
    variant = _mutate(Q, 30, 41)
    copy = _mutate(Q, 3, 42)
    assert ss.passes_own_gene_rule(Q, copy)
    home = _home([("gene-q", Q, "c1", 10_000), ("gene-c", copy, "c1", 12_000), ("gene-far", _protein(200, 77), "c5", 1_000)])
    x = _mutate(Q, 2, 43)
    bits = _q_bits(home, variant, x)
    assert ss.raw_to_bits(dict((r, sc) for sc, r in home.scores(x))["gene-c"]) - bits >= ss.SEED_MIN_MARGIN_BITS
    spec = ss.specificity(home, "GOI_v", variant, x, bits, None, 5, own_gene_id="gene-q")
    assert spec["best_other"] == "gene-far"
    # x's own copy of the neighbour, favoured over the gene as well, stays its competitor.
    far_copy = _mutate(_protein(450, 21), 16, 22)
    long_q = _protein(450, 21)
    home = _home([("gene-q", long_q, "c1", 10_000), ("gene-c", far_copy, "c1", 12_000), ("gene-far", _protein(200, 77), "c5", 1_000)])
    variant = _mutate(long_q, 100, 44)
    bits = _q_bits(home, variant, far_copy)
    assert ss.specificity(home, "GOI_v", variant, far_copy, bits, None, 5, own_gene_id="gene-q")["best_other"] == "gene-c"


def test_a_copy_with_8_to_20_differences_stays_a_competitor():
    for k in (8, 13, 20):
        home = _array(_mutate(Q, k, 6))
        spec = _spec(home, _mutate(Q, 10, 4))
        assert spec["best_other"] == "gene-a0" and spec["real_competitor"], k


def test_a_paralog_in_the_neighbourhood_stays_a_competitor():
    para = _mutate(Q, 40, 8)
    home = _array(Q, para)
    spec = _spec(home, _mutate(para, 3, 9))
    assert spec["best_other"] == "gene-a1" and spec["real_competitor"]


def test_a_copy_outside_qs_own_locus_stays_a_competitor():
    """The rule reads q's own home locus: an identical copy on another chromosome is a place position can tell from q."""
    home = _array(other=[("gene-far-copy", Q, "c7", 5_000)])
    spec = _spec(home, _mutate(Q, 10, 4))
    assert spec["best_other"] == "gene-far-copy" and abs(spec["margin"]) < 1


def test_a_copy_beyond_n_neighbours_stays_a_competitor():
    copies = [_protein(120, 100 + i) for i in range(4)] + [Q]        # 4 unrelated genes, then the copy
    home = _array(*copies)
    assert _spec(home, _mutate(Q, 10, 4), n=4)["best_other"] == "gene-a4"
    assert _spec(home, _mutate(Q, 10, 4), n=5)["best_other"] != "gene-a4"


def test_without_a_flank_size_only_the_own_gene_is_left_out():
    home = _array(Q)
    assert _spec(home, _mutate(Q, 10, 4), n=0)["best_other"] == "gene-a0"


def test_without_an_own_gene_nothing_is_left_out():
    """A query from another species that no home protein matches (a decoy too): every protein competes."""
    query = _protein(120, 2)
    home = _home([("gene-a", _mutate(query, 40, 5), "c1", 1000), ("gene-b", _protein(150, 9), "c1", 3000)])
    assert home.own_copies("GOI_q", query, 5, None) == set()
    x = _mutate(query, 40, 5)
    bits = ss.raw_to_bits(home.scores(x)[0][0])
    assert ss.specificity(home, "GOI_q", query, x, bits, None, 5)["best_other"] == "gene-a"


def test_a_decoy_keeps_every_competitor():
    """A decoy's q is its own reversed sequence (annotate_goi_specificity): no home protein passes the rule for it."""
    decoy = Q[::-1]
    home = _array(Q)
    assert home.own_copies("GOI_q~rev", decoy, 5, None) == set()


def test_the_copies_are_found_once_per_query():
    home = _array(Q)
    calls = []
    real = ss.passes_own_gene_rule
    try:
        ss.passes_own_gene_rule = lambda a, b: (calls.append(1), real(a, b))[1]
        first = home.own_copies("GOI_q", Q, 5, "gene-q")
        n = len(calls)
        assert home.own_copies("GOI_q", Q, 5, "gene-q") == first and len(calls) == n
    finally:
        ss.passes_own_gene_rule = real
    assert first == {"gene-q", "gene-a0"}
