"""Pairwise helpers in bin/sequence_utils.py: is-it-related E-value, aligned query
span, and the score of ONE copy of a gene among the hits of a locus."""
import random
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bin"))

pytest.importorskip("parasail")

from sequence_utils import one_copy_score, sw_align_span, sw_pair_evalue  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice(AA) for _ in range(n))


def _mutate(seq, share, seed):
    rnd = random.Random(seed)
    out = list(seq)
    for i in rnd.sample(range(len(seq)), int(len(seq) * share)):
        out[i] = rnd.choice(AA.replace(out[i], ""))
    return "".join(out)


# ------------------------------------------------------------------- sw_pair_evalue

def test_related_proteins_have_a_tiny_evalue():
    a = _protein(300, 1)
    e, score, ident = sw_pair_evalue(a, _mutate(a, 0.5, 2))
    assert e < 1e-20 and score > 300 and 40 < ident < 65


def test_unrelated_proteins_are_at_chance_whatever_their_lengths():
    """The measure this replaced (common subsequence / longer length) flagged about half
    of all unrelated pairs of similar length. The E-value does not depend on length."""
    flagged = 0
    for i in range(60):
        n = 80 + 37 * (i % 9)
        e, _score, _ident = sw_pair_evalue(_protein(n, 100 + i), _protein(n + i % 5, 500 + i))
        flagged += e <= 1e-5
    assert flagged == 0


def test_a_shared_domain_in_otherwise_unrelated_proteins_is_found():
    domain = _protein(90, 7)
    a = _protein(200, 8) + domain + _protein(150, 9)
    b = _protein(60, 10) + _mutate(domain, 0.3, 11) + _protein(400, 12)
    e, _score, _ident = sw_pair_evalue(a, b)
    assert e < 1e-10


def test_evalue_is_symmetric_and_empty_input_is_infinite():
    a, b = _protein(150, 3), _mutate(_protein(150, 3), 0.4, 4)
    assert sw_pair_evalue(a, b)[0] == pytest.approx(sw_pair_evalue(b, a)[0])
    assert sw_pair_evalue("", a) == (float("inf"), 0.0, 0.0)
    assert sw_pair_evalue(a, "") == (float("inf"), 0.0, 0.0)


def test_evalue_grows_with_the_search_space():
    """Same alignment, longer unrelated tail: the same score is less surprising."""
    a = _protein(60, 21)
    b = _mutate(a, 0.45, 22)
    short = sw_pair_evalue(a, b)
    long_ = sw_pair_evalue(a, b + _protein(3000, 23))
    assert long_[1] >= short[1]
    assert long_[0] > short[0] * 10


# --------------------------------------------------------------------- sw_align_span

def test_span_of_a_full_length_match():
    a = _protein(200, 31)
    r = sw_align_span(a, a)
    assert (r["q_start"], r["q_end"]) == (1, 200)
    assert r["identity"] == 100.0 and r["query_cov"] == 1.0


def test_span_reports_the_part_of_the_query_that_aligned():
    """Subject = residues 51-150 of the query between unrelated flanks."""
    q = _protein(200, 32)
    subject = _protein(40, 33) + q[50:150] + _protein(70, 34)
    r = sw_align_span(q, subject)
    assert abs(r["q_start"] - 51) <= 3 and abs(r["q_end"] - 150) <= 3
    assert r["query_cov"] == pytest.approx(0.5, abs=0.03)
    assert r["identity"] > 95


def test_span_is_measured_on_the_query_not_on_the_subject():
    q = _protein(100, 35)
    long_subject = _protein(300, 36) + q + _protein(300, 37)
    r = sw_align_span(q, long_subject)
    assert (r["q_start"], r["q_end"]) == (1, 100) and r["query_cov"] == 1.0
    r2 = sw_align_span(long_subject, q)
    assert r2["query_cov"] == pytest.approx(100 / 700, abs=0.01)
    assert abs(r2["q_start"] - 301) <= 2


def test_span_with_a_gap_in_the_subject_counts_the_skipped_query_residues():
    """A deletion in the subject does not shorten the stretch of the query that is covered."""
    q = _protein(160, 38)
    subject = q[:70] + q[90:]
    r = sw_align_span(q, subject)
    assert (r["q_start"], r["q_end"]) == (1, 160)
    assert r["query_cov"] == 1.0 and r["identity"] == 100.0


def test_span_of_empty_input():
    assert sw_align_span("", "MKT")["query_cov"] == 0.0
    assert sw_align_span("MKT", "")["q_end"] == 0


def test_span_values_are_plain_python_numbers():
    r = sw_align_span(_protein(50, 39), _protein(50, 39))
    assert type(r["q_start"]) is int and type(r["q_end"]) is int
    assert type(r["identity"]) is float and type(r["query_cov"]) is float


# -------------------------------------------------------------------- one_copy_score

def _hit(qstart, qend, bits):
    return {"qstart": qstart, "qend": qend, "bits": bits}


def test_one_copy_of_a_three_exon_gene_adds_its_exons():
    assert one_copy_score([_hit(1, 100, 80), _hit(101, 200, 70), _hit(201, 300, 60)]) == 210


def test_a_tandem_array_scores_as_one_copy_not_as_the_sum():
    copies = [_hit(1, 300, 400), _hit(1, 300, 380), _hit(5, 298, 350)]
    assert one_copy_score(copies) == 400


def test_weak_relatives_in_a_row_do_not_add_up_to_a_strong_locus():
    """Five genes that each match the same domain of the query: one domain, not five."""
    hits = [_hit(120, 210, 45), _hit(118, 205, 42), _hit(125, 212, 40),
            _hit(119, 208, 39), _hit(121, 209, 38)]
    assert one_copy_score(hits) == 45


def test_best_combination_is_found_not_the_first():
    hits = [_hit(1, 150, 100), _hit(1, 60, 70), _hit(61, 150, 70), _hit(151, 300, 90)]
    assert one_copy_score(hits) == 70 + 70 + 90


def test_small_overlap_between_neighbouring_exons_is_allowed():
    assert one_copy_score([_hit(1, 105, 50), _hit(100, 200, 50)]) == 100
    assert one_copy_score([_hit(1, 105, 50), _hit(100, 200, 50)], max_query_overlap=0) == 50
    assert one_copy_score([_hit(1, 130, 50), _hit(100, 200, 50)]) == 50


def test_no_hits_score_zero_and_order_does_not_matter():
    assert one_copy_score([]) == 0.0
    hits = [_hit(201, 300, 60), _hit(1, 100, 80), _hit(70, 230, 75), _hit(101, 200, 70)]
    assert one_copy_score(hits) == one_copy_score(list(reversed(hits))) == 210
