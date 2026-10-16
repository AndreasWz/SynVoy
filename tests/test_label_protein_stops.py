"""Review fix RF-1a (batch of 2026-10-15): the label protein keeps its stops.

Candidates built from search hits (hit chains, rescued exons, raw hits) carry a record
protein with the stops deleted; the place E-value was computed on that protein, so pieces
on either side of a stop were judged as one sequence. The label protein is now translated
again from the CDS rows with the stops kept and the rows joined with '*', and
split_at_stops judges the pieces. The hull rescue's judge gets its stops back too.
"""
import os
import random
import sys

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
import synvoy_evalue  # noqa: E402
from sequence_utils import reverse_complement, translate  # noqa: E402

CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'}
W_START = 1000


def _dna(protein):
    return "".join(CODON[a] for a in protein)


def _window(strand):
    """A window with two hits of one chain; the first holds an internal stop."""
    rnd = random.Random(5)
    filler = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    hit1 = _dna("MKTAYIAKQR*QISFVKSHF")          # 20 codons, stop at codon 11
    hit2 = _dna("SRQLEERLGLIEVQAPILSR")
    if strand == "-":
        hit1, hit2 = reverse_complement(hit2), reverse_complement(hit1)
    left, mid = filler(300), filler(150)
    seq = left + hit1 + mid + hit2 + filler(300)
    a = (W_START + len(left) + 1, W_START + len(left) + len(hit1))
    b = (a[1] + len(mid) + 1, a[1] + len(mid) + len(hit2))
    return seq, a, b


def _candidate(strand, evidence="fallback_hit_span"):
    seq, a, b = _window(strand)
    rows = [f"c1\tfallback_hits\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=x_CDS{i};Parent=x"
            for i, (s, e) in enumerate((a, b), 1)]
    mrna = (f"c1\tfallback_hits\tmRNA\t{a[0]}\t{b[1]}\t60.0\t{strand}\t.\tID=GOI_q|g_b0_l1_fallback;"
            f"SynVoy_Parent=GOI_q;SynVoyRole=goi;EvidenceType={evidence}")
    record = translate(seq[a[0] - W_START - 1:a[1] - W_START] if strand == "+" else
                       reverse_complement(seq[b[0] - W_START - 1:b[1] - W_START])).replace("*", "")
    return seq, {"start": a[0], "end": b[1], "gff": [mrna] + rows,
                 "record": {"id": "GOI_q|g_b0_l1_fallback", "seq": record}}


def test_a_hit_chain_is_judged_on_the_pieces_between_stops():
    for strand in ("+", "-"):
        seq, cand = _candidate(strand)
        before = cand["record"]["seq"]
        isr.prepare_label_proteins([cand], seq, W_START)
        assert cand["label_protein"] == "MKTAYIAKQR*QISFVKSHF*SRQLEERLGLIEVQAPILSR", strand
        assert cand["record"]["seq"] == before and "*" not in before   # the record is left alone
        assert synvoy_evalue.split_at_stops(cand["label_protein"]) == [
            "MKTAYIAKQR", "QISFVKSHF", "SRQLEERLGLIEVQAPILSR"]


def test_rescued_exons_and_raw_hits_keep_their_stops_too():
    for evidence in ("rescued_exon", "raw_hit"):
        seq, cand = _candidate("+", evidence)
        isr.prepare_label_proteins([cand], seq, W_START)
        assert "*" in cand["label_protein"], evidence


def test_the_place_evalue_sees_only_pieces(monkeypatch):
    seq, cand = _candidate("+")
    isr.prepare_label_proteins([cand], seq, W_START)
    seen = []
    monkeypatch.setattr(synvoy_evalue, "blastp_evalue",
                        lambda q, subjects, sp: (seen.extend(subjects) or (1.0, 0.0)))
    terms = {"significant": True, "slot": None, "core": (W_START, W_START + len(seq)),
             "window": (W_START, W_START + len(seq)), "genome_len": 10 ** 7, "esyn": 1e-3,
             "window_len": len(seq)}
    isr.annotate_place_evalues([cand], {"GOI_q": {"seq": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSR"}},
                               {}, terms, seq, W_START)
    assert seen and all("*" not in s for s in seen)
    assert "MKTAYIAKQRQISFVKSHF" not in seen            # never joined across the stop


def test_exon_models_and_tandem_copies_are_unchanged():
    seq, cand = _candidate("+", "tandem_copy")
    isr.prepare_label_proteins([cand], seq, W_START)
    assert cand["label_protein"] == cand["record"]["seq"]


def test_the_hull_judge_gets_the_stops():
    window = "AAA" + _dna("MKTAYIAKQR*QISFVKSHF") + "TAA" + "GGG"
    rows = [["w", "miniprot", "CDS", "4", str(3 + 63), ".", "+", "0", "Parent=m"]]
    assert hull._translate_model(window, rows) == "MKTAYIAKQRXQISFVKSHF"
    assert hull._translate_model(window, rows, keep_stops=True) == "MKTAYIAKQR*QISFVKSHF"
    src = open(os.path.join(BIN, "rescue_goi_hull.py")).read()
    assert "label_protein=_translate_model(window, cds_rows, keep_stops=True)" in src
