"""Fix-phase item 1: two hits in one stop-free frame are one piece of the label protein.

hit_label_protein joined a hit-built candidate's hits with '*', so split_at_stops never joined two hits (review
fix RF-1a), even when they lay in one open reading frame. orf19.6528 / C. parapsilosis: two hits on the minus
strand, 24 bp apart and 879 bp end to end, in one ORF, the truth gene (897 bp with its stop). Never joined across
a stop: that is what RF-1a is for.
"""
import os
import shutil
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402
import synvoy_label  # noqa: E402
from sequence_utils import reverse_complement, translate  # noqa: E402

CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
         "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT", "*": "TAA"}
ORF = "MKLVTPSGHRDEAFCYWNQIMKLVTPSGHRDEAFCYWNQI"


def _dna(protein):
    return "".join(CODON[a] for a in protein)


def _rows(spans, strand="+"):
    return [f"c1\tfallback_hits\tCDS\t{a}\t{b}\t.\t{strand}\t0\tID=x" for a, b in spans]


def _label(spans, seq, strand="+", w_start=0):
    return isr.hit_label_protein(_rows(spans, strand), seq, w_start)


def _pieces(spans, seq, strand="+"):
    return _label(spans, seq, strand).split("*")


# Plus strand: codon i (0-based) of ORF is bases 3i+1 .. 3i+3.
def _codons(i, j):
    return 3 * i + 1, 3 * j


def test_two_hits_in_one_open_reading_frame_are_one_piece():
    seq = _dna(ORF)
    got = _label([_codons(0, 10), _codons(14, 30)], seq)
    assert got == ORF[:30] and "*" not in got


def test_the_gap_is_translated_not_deleted():
    seq = _dna(ORF)
    got = _label([_codons(0, 10), _codons(14, 30)], seq)
    assert got[10:14] == ORF[10:14]


def test_a_stop_in_the_gap_keeps_them_apart():
    prot = ORF[:12] + "*" + ORF[13:]
    seq = _dna(prot)
    assert _pieces([_codons(0, 10), _codons(14, 30)], seq) == [prot[:10], prot[14:30]]


def test_a_stop_inside_the_second_hit_keeps_them_apart():
    prot = ORF[:20] + "*" + ORF[21:]
    seq = _dna(prot)
    assert _pieces([_codons(0, 10), _codons(14, 30)], seq) == [prot[:10], prot[14:20], prot[21:30]]


def test_another_frame_keeps_them_apart():
    """The second hit starts one base later: no stop anywhere, but not one frame."""
    seq = _dna(ORF)
    first, (s, e) = _codons(0, 10), _codons(14, 30)
    off = seq[s:e + 1]
    assert _label([first, (s + 1, e + 1)], seq) == ORF[:10] + "*" + translate(off[:len(off) - len(off) % 3])


def test_overlapping_hits_in_one_frame_are_one_piece():
    seq = _dna(ORF)
    assert _label([_codons(0, 12), _codons(9, 30)], seq) == ORF[:30]


def test_a_hit_inside_another_is_one_piece():
    seq = _dna(ORF)
    assert _label([_codons(0, 30), _codons(9, 15)], seq) == ORF[:30]


def test_overlapping_hits_in_another_frame_stay_apart():
    seq = _dna(ORF)
    s, e = _codons(9, 20)
    off = seq[s + 1:e + 2]
    assert _label([_codons(0, 12), (s + 2, e + 2)], seq) == ORF[:12] + "*" + translate(off[:len(off) - len(off) % 3])


def test_three_hits_join_as_far_as_the_frame_is_open():
    prot = ORF[:26] + "*" + ORF[27:]
    seq = _dna(prot)
    assert _pieces([_codons(0, 8), _codons(10, 20), _codons(28, 38)], seq) == [prot[:20], prot[28:38]]


def test_only_consecutive_hits_join():
    """A hit in another frame between two in one frame keeps them apart: A, B and C stay three pieces."""
    seq = _dna(ORF)
    a, c = _codons(0, 8), _codons(20, 30)
    s, e = _codons(10, 16)
    got = _label([a, (s + 1, e + 1), c], seq)
    assert got.count("*") == 2 and got.startswith(ORF[:8] + "*") and got.endswith("*" + ORF[20:30])


def test_a_gap_of_ns_keeps_them_apart():
    """An assembly gap: whether the frame is open there is not known, and the length of the N run should not decide."""
    seq = _dna(ORF[:10]) + "N" * 99 + _dna(ORF[14:30])
    assert _pieces([(1, 30), (130, 177)], seq) == [ORF[:10], ORF[14:30]]


def test_the_two_strands_are_mirrors_with_a_row_inside_a_row():
    """O holds I (another frame), N follows O in O's frame. In transcript order O, I, N: I lies between, so nothing joins,
    on both strands. The minus strand's order used to be (start, end) descending: I, O, N, and O joined N."""
    n = 200
    prot = (ORF * 6)[:n]
    plus = _dna(prot)
    o, i, nn = (1, 300), (273, 299), (304, 360)
    mirror = lambda r: (3 * n - r[1] + 1, 3 * n - r[0] + 1)
    got_plus = _label([o, i, nn], plus)
    got_minus = _label([mirror(o), mirror(i), mirror(nn)], reverse_complement(plus), "-")
    assert got_plus == got_minus
    assert got_plus.split("*")[0] == prot[:100] and len(got_plus.split("*")) >= 3


def test_a_hit_with_its_own_stop_joins_nothing():
    prot = ORF[:5] + "*" + ORF[6:]
    seq = _dna(prot)
    assert _pieces([_codons(0, 10), _codons(14, 30)], seq) == [prot[:5], prot[6:10], prot[14:30]]


# Minus strand: the genome holds the reverse complement. Rows are genomic.
def _minus(i, j, total):
    return 3 * (total - j) + 1, 3 * (total - i)


def test_the_minus_strand_joins_in_transcript_order():
    seq = reverse_complement(_dna(ORF))
    n = len(ORF)
    assert _label([_minus(14, 30, n), _minus(0, 10, n)], seq, "-") == ORF[:30]
    assert _label([_minus(0, 10, n), _minus(14, 30, n)], seq, "-") == ORF[:30]


def test_the_minus_strand_keeps_a_stop_and_another_frame_apart():
    prot = ORF[:12] + "*" + ORF[13:]
    seq = reverse_complement(_dna(prot))
    n = len(prot)
    assert _pieces([_minus(0, 10, n), _minus(14, 30, n)], seq, "-") == [prot[:10], prot[14:30]]
    seq = reverse_complement(_dna(ORF))
    s, e = _minus(14, 30, n)
    off = reverse_complement(seq[s - 2:e - 1])
    assert _label([_minus(0, 10, n), (s - 1, e - 1)], seq, "-") == ORF[:10] + "*" + translate(off[:len(off) - len(off) % 3])


def test_the_minus_strands_frame_reads_from_the_row_ends():
    """A row cut one base short at its 3' end (a low genomic coordinate) keeps its frame: the frame reads from the transcript's
    start, the row's end coordinate."""
    seq = reverse_complement(_dna(ORF))
    n = len(ORF)
    a, (s, e) = _minus(0, 10, n), _minus(14, 30, n)
    assert _label([a, (s + 1, e)], seq, "-") == ORF[:29]


def test_the_case_of_the_item():
    """The shape of orf19.6528 / C. parapsilosis: minus strand, 24 bp apart, 56 and 229 codons, one ORF."""
    prot = (ORF * 8)[:56 + 8 + 229]
    seq = reverse_complement(_dna(prot))
    n = len(prot)
    a, b = _minus(0, 56, n), _minus(64, 293, n)
    assert a[0] - b[1] == 25 and b[1] - b[0] + 1 == 687 and a[1] - a[0] + 1 == 168
    got = _label([a, b], seq, "-")
    assert got == prot and len(synvoy_evalue.split_at_stops(got)) == 1


# The two hits themselves: Candida parapsilosis cpar_Chr_7:753241-754119, as the assembly reads it, and the home
# protein of orf19.6528. On the minus strand the hits are 753952-754119 and 753241-753927.
CPAR_START = 753241
CPAR_WINDOW = (
    "ATGATACCCTTTTAATCCTTGGCCTAATGCATTTACAACCGAATCAGTCTTGATCTTGGCAATAGCGTCAACTTCAAAATCAAAAGAAGAGAGATCTCTTTCTTCC"
    "CTTGCTGTTGCTGCTGTTGCTGCTGTTGCTGTTGTCGTACTTAAATTTCCTGATTGGTAGAATAGTACGACAAAATACCCATTAAACAAATCTACATCATAATCTT"
    "CATAAGCGTTATCAAATTCAGCATCTTTATCATCATTAGTTGTAGATTTCGTGTTATTTTCCCCACTACGCTTGATATCTTGTCCATTATCTTTCAATTGCTTAGC"
    "GCCTCTATTATCTTTTTGTCGTGATGGGTTTTGATCAAGGTTAACTTTTGAGATTATACATTTAACAGTACCGAAATCAACAATTGACCAATCATTAGCATCATCG"
    "GTGTCACTACTTGACTCGTCAGAGTCGGTGGCTTCGGTAATTGTATTGTATGCAATTAATGAATATATCTTTAATTGATCCAAGGTGATGAACGACGTTCCATCTT"
    "GCGAGTCTGCGTTACCCCCTGTATAAACTGATGGTGTTGCTCTAGCGTGGTCGTATTCTTCACCTGTTACTGTGGGCATTGTTGTTGTTGTTGGAGTACCAGCTTT"
    "GGAACTATGCAAGGAGTGCAATTGCAGTTGTGGCTGCAGTTGCAGCTGCAGTTGATGATGATCAACACGCTCCGTGGATTGACTCTTTGTAGTGAGTGAGATTAAT"
    "GGATTGCCATGGTTGGATAATAATGAAATTGAAATTGGATGGAAATTGGGTACTGAGGGAGTTGATGATGTTGTTGTAGATCCCGGTGACAGTAGTGCTGATTTGA"
    "GCAAGTTTGATATTTGATTCTGCTGCAACAT"
)
ORF19_6528 = (
    "MLQQKQISKILQQGLTPITLNSKSNSQSATTIPTSTHSIKSISLLSKQGIPLITVSSTTDSGSGIIDQDYKIYSIIAYNSLNKPNENDNDNDNDNSKNTNKLKSKENE"
    "TWTIIQFEPGLKCMIERCLDMYLVLYYDGRVDNKDRDAGTDNKNNKLEGAANGDAGNVKTKETRQLEVEVDAFVKLKMDGLVLSLHDGLLGYNNE"
)
CPAR_HITS = [(753952, 754119), (753241, 753927)]


def _cpar_label(spans):
    return isr.hit_label_protein(_rows(spans, "-"), CPAR_WINDOW, CPAR_START - 1)


def test_the_two_real_hits_are_one_piece():
    """Translated through the 24-bp gap: 56 + 8 + 229 codons, no stop. The gap's residues are in the piece."""
    assert len(CPAR_WINDOW) == 879
    got = _cpar_label(CPAR_HITS)
    assert len(got) == 293 and "*" not in got and got.startswith("MLQQNQISNL")
    assert got[56:64] == "TERVDHHQ"
    assert synvoy_evalue.split_at_stops(got) == [got]


def test_each_real_hit_alone_reads_its_own_residues():
    got = _cpar_label(CPAR_HITS)
    assert _cpar_label(CPAR_HITS[:1]) == got[:56]
    assert _cpar_label(CPAR_HITS[1:]) == got[64:]


@pytest.mark.skipif(shutil.which("blastp") is None, reason="needs blastp")
def test_what_the_join_does_to_E_on_the_real_hits():
    """E in the run's slot (7,330 bp). The main hit alone reads 3.35e-4, 1.005e-3 with the 3 significant blocks, just
    past the MEDIUM bar of 1e-3 (AMBIGUOUS on LRZ), and the piece 9.34e-5, 2.8e-4 with them (MEDIUM). E is that of one
    alignment and none crosses the gap: the 9.34e-5 is the first hit's own alignment, 2.05e-3 alone. The gain is
    composition statistics on the longer subject. Allow a blastp version some room, but not across the piece's bar."""
    got = _cpar_label(CPAR_HITS)
    space = synvoy_evalue.search_space(len(ORF19_6528), 7330)
    e_first = synvoy_evalue.blastp_evalue(ORF19_6528, [got[:56]], space)[0]
    e_main = synvoy_evalue.blastp_evalue(ORF19_6528, [got[64:]], space)[0]
    e_piece = synvoy_evalue.blastp_evalue(ORF19_6528, [got], space)[0]
    assert e_main == pytest.approx(3.35e-4, rel=0.25)
    assert e_first == pytest.approx(2.05e-3, rel=0.25)
    assert e_piece == pytest.approx(9.34e-5, rel=0.25)
    assert e_piece < e_main < e_first
    assert 3 * e_piece < synvoy_label.MEDIUM_EVALUE


def test_a_window_offset_is_honoured():
    seq = "N" * 100 + _dna(ORF)
    s1, e1 = _codons(0, 10)
    s2, e2 = _codons(14, 30)
    assert isr.hit_label_protein(_rows([(s1 + 1100, e1 + 1100), (s2 + 1100, e2 + 1100)]), seq, 1000) == ORF[:30]


def test_a_single_hit_reads_as_before():
    seq = _dna(ORF)
    assert _label([_codons(0, 10)], seq) == ORF[:10]


def test_no_cds_row_gives_no_protein():
    assert isr.hit_label_protein(["c1\tx\tmRNA\t1\t30\t.\t+\t.\tID=x"], _dna(ORF), 0) == ""


def test_the_candidate_is_still_judged_on_the_hits_themselves(monkeypatch):
    """label_subjects: the joined piece is added to the window hits, which stay subjects."""
    seq = _dna(ORF)
    rows = _rows([_codons(0, 10), _codons(14, 30)])
    cand = {"start": 1, "end": 90, "label_protein": isr.hit_label_protein(rows, seq, 0)}
    hits = {"GOI_x": [{"gstart": 0, "gend": 30, "strand": "+"}, {"gstart": 42, "gend": 90, "strand": "+"}]}
    subjects = isr.label_subjects(cand, "GOI_x", hits, seq, 0)
    assert subjects == [ORF[:30], ORF[:10], ORF[14:30]]


def test_translate_is_what_the_stop_test_reads():
    assert translate("TAA") == "*"
