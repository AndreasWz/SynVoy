"""Shared home-hit helpers (batch of 2026-10-15, S0; no output change): load_home_hits (moved
from split_loci.load_raw_hits, strand-blind as before, now also with strand and E-value) and
one query-continuation predicate in coding order (continues_query, in_coding_order), for the
home hit chaining, the adaptive flank walk and the GOI hit chain."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import split_loci as sl  # noqa: E402
from sequence_utils import (QUERY_CONTINUATION_TOL, continues_query, in_coding_order,  # noqa: E402
                            load_home_hits)


def _m8(path, rows):
    with open(path, "w") as fh:
        for q, c, qs, qe, ss, se, e, b in rows:
            fh.write(f"{q}\t{c}\t50.0\t40\t10\t0\t{qs}\t{qe}\t{ss}\t{se}\t{e}\t{b}\n")
    return str(path)


def test_strand_and_evalue_are_kept_and_the_dedup_is_strand_blind(tmp_path):
    path = _m8(tmp_path / "h.m8", [
        ("q", "c1", 1, 40, 1000, 1119, 1e-20, 80.0),     # plus
        ("q", "c1", 41, 88, 5000, 4857, 1e-10, 60.0),    # minus (subject coordinates backwards)
        ("q", "c1", 1, 40, 1110, 1001, 1e-5, 30.0),      # minus, on the first hit's DNA: dropped
        ("q", "c1", 90, 150, 9000, 9182, 0.5, 70.0),     # E above the bar: dropped
    ])
    hits = load_home_hits([path], 0.01)
    assert [(h["start"], h["end"], h["strand"], h["evalue"], h["bits"]) for h in hits] == [
        (999, 1119, "+", 1e-20, 80.0), (4856, 5000, "-", 1e-10, 60.0)]
    assert [(h["qstart"], h["qend"]) for h in hits] == [(1, 40), (41, 88)]
    # split_loci reads the same hits as before (the keys it used are unchanged)
    old = sl.load_raw_hits([path], 0.01)
    assert [{k: h[k] for k in ("chrom", "start", "end", "bits", "qstart", "qend")} for h in old] == \
           [{k: h[k] for k in ("chrom", "start", "end", "bits", "qstart", "qend")} for h in hits]


def _hit(qs, qe, start=0, end=0):
    return {"qstart": qs, "qend": qe, "start": start, "end": end}


def test_the_predicate_at_its_boundaries():
    a = _hit(1, 40)
    tol = QUERY_CONTINUATION_TOL
    assert tol == 10
    assert continues_query(a, _hit(41, 88))                      # the next exon
    assert continues_query(a, _hit(40 - tol, 60))                # overlap of exactly tol: allowed
    assert not continues_query(a, _hit(40 - tol - 1, 60))        # one residue more: refused
    assert not continues_query(a, _hit(35, 40))                  # inside the end of a: refused
    assert not continues_query(a, _hit(30, 40))                  # (the end test does this)
    assert not continues_query(a, _hit(1, 40))                   # the same query part again
    assert not continues_query(_hit(41, 88), a)                  # backwards in the query


def test_both_strands_in_coding_order():
    exons = [(1, 40), (41, 88), (90, 150)]
    plus = [_hit(qs, qe, start, start + 3 * (qe - qs + 1))
            for (qs, qe), start in zip(exons, (1000, 3000, 7000))]
    # On the minus strand the first exon in coding order lies downstream.
    minus = [_hit(qs, qe, start, start + 3 * (qe - qs + 1))
             for (qs, qe), start in zip(exons, (7000, 3000, 1000))]
    for hits, strand in ((plus, "+"), (minus, "-")):
        ordered = in_coding_order(list(reversed(hits)), strand)
        assert [(h["qstart"], h["qend"]) for h in ordered] == exons
        assert all(continues_query(x, y) for x, y in zip(ordered, ordered[1:]))
    # In genomic order the minus strand reads backwards in the query: nothing continues.
    genomic = sorted(minus, key=lambda h: h["start"])
    assert not any(continues_query(x, y) for x, y in zip(genomic, genomic[1:]))
