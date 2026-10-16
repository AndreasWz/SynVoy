"""Tandem copies of the GOI are query restarts (batch of 2026-10-15, G2; --tandem_query_restart,
off by default). Copies are built per strand from hits that continue the query, a copy counts at
>= 0.35 union query coverage with a best hit of >= 20 aa or >= 30 bits (both under 150 aa), copies
that share genome are one, and it is a tandem when two counted copies overlap in the query by
more than half. The spec's cases (scripts/replay_tandem.py, with bit scores added) and the
duplicate-hit case; all synthetic."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import annotate_goi_exons as ag  # noqa: E402
import iterative_search_runner as isr  # noqa: E402

DNA = "GCT" * 20000


def H(q1, q2, g1, g2, pid=60, ev=1e-10, strand="+", bits=80.0):
    return dict(qstart=q1, qend=q2, gstart=g1, gend=g2, pident=pid, evalue=ev, strand=strand,
                chrom="c", bits=bits)


def run(hits, qlen, restart=True):
    return ag.detect_tandem_duplications(hits, "M" * qlen, DNA, "c", query_restart=restart)


def test_1_and_2_two_exons_of_one_gene_are_one_copy():
    for intron in (6000, 3000):
        hits = [H(1, 40, 1000, 1120), H(41, 80, 1120 + intron, 1240 + intron)]
        assert run(hits, 80)[0] is False


def test_3_two_exons_on_opposite_strands_are_no_tandem():
    hits = [H(1, 40, 1000, 1120, strand="+"), H(41, 80, 9000, 9120, strand="-")]
    assert run(hits, 80)[0] is False
    assert run(hits, 80, restart=False)[0] is True                  # today: two hits >= 35 %


def test_4_three_short_hits_on_one_motif_are_too_small():
    hits = [H(30, 41, 1000, 1036, pid=40, strand="+", bits=25), H(31, 42, 20000, 20036, pid=26, strand="-", bits=20),
            H(30, 41, 27000, 27036, pid=45, strand="+", bits=26)]
    assert run(hits, 736)[0] is False


def test_5_two_full_length_copies_are_a_tandem():
    hits = [H(1, 80, 1000, 1240), H(1, 80, 4240, 4480)]
    tandem, copies = run(hits, 80)
    assert tandem is True and [(c["gstart"], c["gend"]) for c in copies] == [(1000, 1240), (4240, 4480)]
    assert all(len(c["hits"]) == 1 and c["qcov"] == 1.0 and c["id"].startswith("GOI_copy_") for c in copies)
    assert all((h["gend"] - h["gstart"]) % 3 == 0 for c in copies for h in c["hits"])


def test_a_duplicate_hit_does_not_fake_a_tandem():
    # Two programs' hits of one copy: the same place (dedup), or 60 % shared genome (one copy).
    assert run([H(1, 80, 1000, 1240, bits=90), H(1, 78, 1006, 1240, bits=85)], 80)[0] is False
    assert run([H(1, 80, 1000, 1240, bits=90), H(1, 80, 1096, 1336, bits=85)], 80)[0] is False


def test_copies_of_two_exons_each_and_an_inverted_copy():
    hits = [H(1, 40, 1000, 1120), H(41, 80, 3120, 3240), H(1, 40, 9000, 9120), H(41, 80, 11120, 11240)]
    tandem, copies = run(hits, 80)
    assert tandem is True and [len(c["hits"]) for c in copies] == [2, 2]
    inv = [H(1, 80, 1000, 1240, strand="+"), H(1, 80, 6000, 6240, strand="-")]
    assert run(inv, 80)[0] is True


def test_a_short_query_copy_needs_both_length_and_bits():
    hits = [H(1, 25, 1000, 1075, bits=20), H(1, 25, 4000, 4075, bits=20)]     # 25 aa, 20 bits
    assert run(hits, 70)[0] is False
    strong = [H(1, 25, 1000, 1075, bits=35), H(1, 25, 4000, 4075, bits=35)]
    assert run(strong, 70)[0] is True


def test_the_judge_reads_a_copy_s_hits():
    subseq = "ATG" + "GCT" * 100 + "TAA" + "GCT" * 20
    gff = ["c\ttandem_copy\tmRNA\t1\t300\t90\t+\t.\tID=x;EvidenceType=tandem_copy",
           "c\ttandem_copy\tCDS\t1\t150\t.\t+\t0\tID=x_CDS1;Parent=x",
           "c\ttandem_copy\tCDS\t154\t300\t.\t+\t0\tID=x_CDS2;Parent=x"]
    cand = {"gff": gff, "record": {"seq": "MAAA"}, "hit_based": True}
    isr.prepare_label_proteins([cand], subseq, 0)
    assert cand["label_protein"].count("*") == 1
    single = {"gff": ["c\ttandem_copy\tgene\t1\t300\t90\t+\t.\tID=y;EvidenceType=tandem_copy"],
              "record": {"seq": "MAAA"}}
    isr.prepare_label_proteins([single], subseq, 0)
    assert single["label_protein"] == "MAAA"                           # the default: one span
    rows_only = {"gff": gff, "record": {"seq": "MAAA"}}                 # not built from hits
    isr.prepare_label_proteins([rows_only], subseq, 0)
    assert rows_only["label_protein"] == "MAAA"
