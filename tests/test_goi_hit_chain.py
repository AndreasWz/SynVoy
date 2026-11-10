"""The GOI's hit chain when miniprot gives no model (fallback_hit_span), and batch G1's link rule.

G1's chain that followed the query (--goi_chain_query_aware, 2026-10-15) was removed on
2026-11-10: switched on, it refused a same-frame continuation of the query and undid a gain on
the curated truth (main.nf removedParams; tests/test_removed_params_preview.py). The chain is the
one whose query centres rise with the most hits, as it was with the switch off. G1's link rule
(_query_chain_link) stays: the leftover join of batch G3B uses it (tests/test_join_leftover_hits.py).
These cases were in tests/test_goi_chain_query_aware.py; query 300 aa, all synthetic."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402


def _h(qs, qe, gs, ge, bits, strand="+"):
    return {"qstart": qs, "qend": qe, "gstart": gs, "gend": ge, "bits": bits,
            "alnlen": qe - qs + 1, "strand": strand}


def _q(chain):
    return [(h["qstart"], h["qend"]) for h in chain]


def _chain(hits, strand="+", max_gap=50000):
    ordered = sorted(hits, key=lambda h: h["gstart"])
    return isr._longest_monotonic_query_chain(ordered, strand, max_gap=max_gap)


def test_the_chain_is_the_one_whose_query_centres_rise():
    hits = [_h(50, 150, 1000, 1303, 200), _h(90, 110, 2303, 2366, 40)]
    assert _q(_chain(hits)) == [(50, 150), (90, 110)]
    hits = [_h(100, 114, 1000, 1045, 20), _h(150, 164, 2045, 2090, 20), _h(200, 214, 3090, 3135, 20),
            _h(20, 270, 9000, 9753, 400)]
    assert _q(_chain(hits)) == [(100, 114), (150, 164), (200, 214)]


def test_the_chain_has_no_query_aware_switch_any_more():
    import inspect
    assert "query_aware" not in inspect.signature(isr._longest_monotonic_query_chain).parameters
    assert not hasattr(isr, "_query_aware_chain") and not hasattr(isr, "merge_one_exon_pieces")


def test_the_chain_query_attribute():
    assert isr.chain_query_attr([_h(41, 88, 2120, 2264, 80), _h(1, 40, 1000, 1120, 70)]) == \
        "ChainQuery=1-40,41-88"


def _link(a, b, strand="+", max_gap=50000):
    return isr._query_chain_link(a, b, strand, max_gap)


def test_g1_link_rule():
    assert _link(_h(10, 60, 1000, 1153, 80), _h(12, 62, 4153, 4306, 70)) == (False, False)  # same part twice
    assert _link(_h(50, 150, 1000, 1303, 200), _h(90, 110, 2303, 2366, 40)) == (False, False)  # inside
    assert _link(_h(1, 40, 1000, 1120, 70), _h(250, 290, 1420, 1543, 60)) == (False, False)  # 209 aa in 300 bp
    assert _link(_h(100, 130, 1000, 1093, 50), _h(131, 160, 1069, 1159, 45)) == (False, False)  # 24 bp overlap
    assert _link(_h(1, 40, 1000, 1120, 70), _h(41, 88, 2120, 2264, 80)) == (True, False)  # a next exon
    assert _link(_h(41, 88, 4183, 4327, 80, "-"), _h(1, 40, 5327, 5447, 70, "-"), "-") == (True, False)
    a = _h(10, 50, 1000, 1123, 60)
    assert _link(a, _h(52, 90, 1126, 1243, 55)) == (True, True)    # one codon skipped, same frame
    assert _link(a, _h(52, 90, 1127, 1244, 55)) == (False, False)  # another frame


def test_a_hit_candidate_is_judged_on_its_hits():
    gff = ["c\tfallback_hits\tmRNA\t1\t300\t90\t+\t.\tID=x;EvidenceType=fallback_hit_span",
           "c\tfallback_hits\tCDS\t1\t300\t.\t+\t0\tID=x_CDS1;Parent=x"]
    plain = {"gff": gff, "record": {"seq": "M"}}
    isr.prepare_label_proteins([plain], "ATG" + "GCT" * 120, 0)
    assert "*" not in plain["label_protein"]                          # one row: one piece
