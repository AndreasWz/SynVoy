"""The default tandem detector (annotate_goi_exons.detect_tandem_duplications, without
--tandem_query_restart): what it calls copies, pinned on synthetic hits.

Two hits that each cover >= 35 % of the query make an array; hits less than 5 kb apart are one
copy. So two exons of one gene more than 5 kb apart are two copies, and two copies less than
5 kb apart are no array (miniprot's models find such a pair; the search uses these copies only
where miniprot has none). docs/ALGORITHM.md § 7 step 3 says so.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import annotate_goi_exons as age  # noqa: E402

DNA = "GCT" * 15000          # 45 kb without a stop codon in any forward frame
QUERY = "M" + "A" * 199      # 200 aa


def _hit(qs, qe, gs, ge):
    return {"qstart": qs, "qend": qe, "gstart": gs, "gend": ge, "strand": "+", "pident": 70.0,
            "evalue": 1e-20, "query": "GOI_X"}


def _detect(hits):
    return age.detect_tandem_duplications(hits, QUERY, DNA, "chr1")


def test_two_exons_of_one_gene_far_apart_are_two_copies():
    is_tandem, copies = _detect([_hit(1, 80, 1000, 1240), _hit(121, 200, 9300, 9540)])
    assert is_tandem is True
    assert [(c["gstart"], c["gend"]) for c in copies] == [(1000, 1240), (9300, 9540)]
    assert [c["id"] for c in copies] == ["GOI_copy_1", "GOI_copy_2"]


def test_two_copies_less_than_5_kb_apart_are_no_array():
    assert _detect([_hit(1, 180, 1000, 1540), _hit(5, 190, 4600, 5155)]) == (False, [])


def test_two_copies_more_than_5_kb_apart_are_two_copies():
    is_tandem, copies = _detect([_hit(1, 180, 1000, 1540), _hit(5, 190, 9600, 10155)])
    assert is_tandem is True
    assert [(c["gstart"], c["gend"]) for c in copies] == [(1000, 1540), (9600, 10155)]


def test_the_detector_takes_no_distance_argument():
    """max_intergenic_distance was accepted and never read."""
    import inspect
    params = inspect.signature(age.detect_tandem_duplications).parameters
    assert "max_intergenic_distance" not in params
    assert list(params) == ["hits", "query_seq", "chrom_seq", "chrom_name", "query_restart", "max_gap"]
