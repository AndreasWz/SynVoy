"""Review fix RF-1e (batch of 2026-10-15): one protein for Bits, specificity and the E-value.

The place E-value is the best over the candidate's label protein (in pieces between stops)
and its hit regions; Bits= and the specificity used the record protein, which for a hit
chain is the stop-deleted join of its hits. Now Bits= is the best of the subjects the E is
computed on, and the specificity runs on that same subject (cand['judged_protein']).
"""
import os
import random
import sys
from types import SimpleNamespace

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
from sequence_utils import sw_bits  # noqa: E402

CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'}
HOME = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV"


def _chain():
    """A hit chain over the home protein's first 40 residues with a stop after residue 20."""
    rnd = random.Random(9)
    filler = "".join(rnd.choice("ACGT") for _ in range(300))
    hit = "".join(CODON[a] for a in HOME[:20] + "*" + HOME[21:40])
    seq = filler + hit + filler
    s, e = len(filler) + 1, len(filler) + len(hit)
    mrna = (f"c1\tfallback_hits\tmRNA\t{s}\t{e}\t60.0\t+\t.\tID=GOI_q|g_b0_l1_fallback;"
            f"SynVoy_Parent=GOI_q;SynVoyRole=goi;EvidenceType=fallback_hit_span")
    cds = f"c1\tfallback_hits\tCDS\t{s}\t{e}\t.\t+\t0\tID=x_CDS1;Parent=x"
    record = (HOME[:20] + HOME[21:40])                    # the stop-deleted join
    return seq, {"start": s, "end": e, "gff": [mrna, cds], "record": {"id": "GOI_q|g_b0_l1_fallback", "seq": record}}


def test_bits_and_specificity_use_the_best_piece(monkeypatch):
    seq, cand = _chain()
    by_parent = {"GOI_q": {"seq": HOME}}
    isr.prepare_label_proteins([cand], seq, 0)
    isr.prepare_label_subjects([cand], {}, seq, 0)
    assert cand["label_subjects"] == [HOME[:20], HOME[21:40]]
    isr.score_goi_candidates([cand], by_parent)
    assert cand["judged_protein"] == max(cand["label_subjects"], key=lambda s: sw_bits(HOME, s))
    assert cand["bits"] == sw_bits(HOME, cand["judged_protein"]) < sw_bits(HOME, cand["record"]["seq"])
    seen = {}

    def fake_spec(proteome, parent, q, x, bits, hood, n_flanking=0, own_gene_id=None):
        seen["x"], seen["bits"] = x, bits
        return {"margin": 1.0, "real_competitor": False}
    monkeypatch.setattr(isr, "_home_proteome", lambda args: (object(), None))
    monkeypatch.setattr(isr.synvoy_specificity, "specificity", fake_spec)
    monkeypatch.setattr(isr.synvoy_specificity, "gff_attrs", lambda spec: "SpecMargin=1.0")
    isr.annotate_goi_specificity([cand], by_parent, SimpleNamespace())
    assert seen == {"x": cand["judged_protein"], "bits": cand["bits"]}


def test_without_prepared_subjects_the_record_protein_is_scored():
    _seq, cand = _chain()
    isr.score_goi_candidates([cand], {"GOI_q": {"seq": HOME}})
    assert cand["bits"] == sw_bits(HOME, cand["record"]["seq"])


def test_the_hull_scores_the_best_piece(monkeypatch):
    seen = {}
    monkeypatch.setattr(hull.synvoy_evalue, "blastp_evalue", lambda q, s, sp: (1e-30, 200.0))

    def fake_spec(proteome, qid, q, x, bits, hood, n_flanking=0, own_gene_id=None):
        seen["x"] = x
        return {"margin": 50.0, "real_competitor": False}
    monkeypatch.setattr(hull.synvoy_specificity, "specificity", fake_spec)
    monkeypatch.setattr(hull.synvoy_specificity, "gff_attrs", lambda spec: "SpecMargin=50")
    label = HOME[:20] + "*" + HOME[21:40]
    rows = [["w", "miniprot", "CDS", "1", "123", ".", "+", "0", "Parent=m"]]
    _conf, _key, _attrs, bits = hull.judge_hull_model("q", HOME, HOME[:20] + "X" + HOME[21:40], rows,
                                                      10 ** 7, object(), label_protein=label, coverage=0.7)
    assert seen["x"] in (HOME[:20], HOME[21:40]) and bits == sw_bits(HOME, seen["x"])
