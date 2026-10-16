"""The GOI's hit chain follows the query (batch of 2026-10-15, G1; --goi_chain_query_aware, off
by default). With the switch a link must continue the query (strand-aware, 10 aa tolerance),
the residues it skips must fit into the gap, and a gap under 40 bp joins only two pieces of one
exon (same frame), which become one CDS row; the chain with the most bits wins. The judge still
reads the original hits. Cases of the gene-model spec, query 300 aa; all synthetic."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402


def _h(qs, qe, gs, ge, bits, strand="+"):
    return {"qstart": qs, "qend": qe, "gstart": gs, "gend": ge, "bits": bits,
            "alnlen": qe - qs + 1, "strand": strand}


def _chain(hits, strand="+", query_aware=True, max_gap=50000):
    ordered = sorted(hits, key=lambda h: h["gstart"])
    return isr._longest_monotonic_query_chain(ordered, strand, max_gap=max_gap, query_aware=query_aware)


def _q(chain):
    return [(h["qstart"], h["qend"]) for h in chain]


def test_a_the_same_query_part_twice_is_one_hit():
    hits = [_h(10, 60, 1000, 1153, 80), _h(12, 62, 4153, 4306, 70)]
    assert _q(_chain(hits)) == [(10, 60)]


def test_b_one_strong_hit_beats_three_weak_pieces():
    # Three 20-bit 15-aa hits that chain among themselves, inside the query part the 400-bit
    # hit covers: the 400-bit hit alone (today the chain with the most hits wins).
    hits = [_h(100, 114, 1000, 1045, 20), _h(150, 164, 2045, 2090, 20), _h(200, 214, 3090, 3135, 20),
            _h(20, 270, 9000, 9753, 400)]
    assert _q(_chain(hits)) == [(20, 270)]
    assert _q(_chain(hits, query_aware=False)) == [(100, 114), (150, 164), (200, 214)]


def test_c_a_hit_inside_another_is_not_a_next_exon():
    hits = [_h(50, 150, 1000, 1303, 200), _h(90, 110, 2303, 2366, 40)]
    assert _q(_chain(hits)) == [(50, 150)]


def test_d_skipped_residues_must_fit_into_the_gap():
    hits = [_h(1, 40, 1000, 1120, 70), _h(250, 290, 1420, 1543, 60)]   # 209 aa skipped, 300 bp
    assert _q(_chain(hits)) == [(1, 40)]


def test_e_an_overlap_that_is_no_codon_shift_is_refused():
    hits = [_h(100, 130, 1000, 1093, 50), _h(131, 160, 1069, 1159, 45)]  # 24 bp genomic overlap
    assert _q(_chain(hits)) == [(100, 130)]


def test_a_real_three_exon_gene_on_both_strands():
    plus = [_h(1, 40, 1000, 1120, 70), _h(41, 88, 2120, 2264, 80), _h(90, 150, 5264, 5447, 90)]
    assert _q(_chain(plus)) == [(1, 40), (41, 88), (90, 150)]
    minus = [_h(90, 150, 1000, 1183, 90, "-"), _h(41, 88, 4183, 4327, 80, "-"),
             _h(1, 40, 5327, 5447, 70, "-")]
    assert _q(_chain(minus, "-")) == [(90, 150), (41, 88), (1, 40)]   # genomic order


def test_one_exon_broken_for_one_codon_becomes_one_cds():
    a, b = _h(10, 50, 1000, 1123, 60), _h(52, 90, 1126, 1243, 55)    # one codon skipped, same frame
    chain = _chain([a, b])
    assert _q(chain) == [(10, 50), (52, 90)]
    assert isr.merge_one_exon_pieces(chain, "+") == [(1000, 1243)]
    shifted = _h(52, 90, 1127, 1244, 55)                             # another frame: refused
    assert _q(_chain([a, shifted])) == [(10, 50)]


def test_no_link_gives_the_best_single_hit_and_the_attribute():
    assert _q(_chain([_h(1, 40, 1000, 1120, 30), _h(1, 40, 90000, 90120, 50)])) == [(1, 40)]
    assert isr.chain_query_attr([_h(41, 88, 2120, 2264, 80), _h(1, 40, 1000, 1120, 70)]) == \
        "ChainQuery=1-40,41-88"


def test_off_by_default_today_s_chain():
    hits = [_h(50, 150, 1000, 1303, 200), _h(90, 110, 2303, 2366, 40)]
    assert _q(_chain(hits, query_aware=False)) == [(50, 150), (90, 110)]   # centres rise: chained


def test_the_judge_reads_the_hits_also_when_cds_rows_were_merged():
    subseq = "ATG" + "GCT" * 100 + "TAA" + "GCT" * 20
    gff = ["c\tfallback_hits\tmRNA\t1\t300\t90\t+\t.\tID=x;EvidenceType=fallback_hit_span",
           "c\tfallback_hits\tCDS\t1\t300\t.\t+\t0\tID=x_CDS1;Parent=x"]
    label_rows = ["c\tfallback_hits\tCDS\t1\t150\t.\t+\t0\tParent=x",
                  "c\tfallback_hits\tCDS\t154\t300\t.\t+\t0\tParent=x"]
    cand = {"gff": gff, "record": {"seq": "M"}, "label_rows": label_rows}
    isr.prepare_label_proteins([cand], subseq, 0)
    assert cand["label_protein"].count("*") == 1                     # two pieces, joined by '*'
    plain = {"gff": gff, "record": {"seq": "M"}, "label_rows": None}
    isr.prepare_label_proteins([plain], subseq, 0)
    assert "*" not in plain["label_protein"]                          # one row: one piece
