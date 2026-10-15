"""A bit score for every GOI candidate, and the choice among models by it (label rework,
steps 1b and 1c, 2026-10-14).

Every GOI candidate of the search gets the Smith-Waterman bit score of its protein against
the home protein (``Bits=``; BLOSUM62, BLAST gaps 11/1, lambda 0.267, K 0.041). Overlapping
candidates in a window are chosen by it (then the shorter span), the report keeps the
member of a duplicate group with the best label and then the best bit score, and the
figures pick their representative GOI model the same way. Until then only hit chains had
a bit score and every choice went by identity: on the Candida pilot that picked the model
on the curated gene in 16 of 31 genome pairs with several models, bits in 29.
"""
import math
import os
import sys

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
parasail = pytest.importorskip("parasail")
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import plot_synteny as ps  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
from sequence_utils import SW_K_11_1, SW_LAMBDA_11_1, sw_bits  # noqa: E402

HOME = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWEAQ"


def _mutate(seq, every):
    """Change every `every`-th residue (A <-> W), so the bit score drops in steps."""
    return "".join(("W" if c != "W" else "A") if i % every == 0 else c for i, c in enumerate(seq))


# ---------------------------------------------------------------- sequence_utils.sw_bits

def test_bits_follow_karlin_altschul_at_blast_11_1():
    score = parasail.sw_striped_32(HOME, HOME, 12, 1, parasail.blosum62).score
    expect = (SW_LAMBDA_11_1 * score - math.log(SW_K_11_1)) / math.log(2)
    assert sw_bits(HOME, HOME) == pytest.approx(expect)


def test_bits_fall_with_divergence_and_are_zero_for_nothing():
    full = sw_bits(HOME, HOME)
    assert full > sw_bits(HOME, _mutate(HOME, 7)) > sw_bits(HOME, _mutate(HOME, 3)) > 0
    assert sw_bits(HOME, "") == 0.0 and sw_bits("", HOME) == 0.0


# ---------------------------------------------------------------- the search (runner)

def _cand(mid, start, end, seq, parent="GOI_q", ev="exon_annotation", ident=50.0, flank=0):
    attrs = (f"ID={parent}|g_b0_{mid};Name={parent};SynVoy_Parent={parent};Identity={ident};"
             f"SynVoyRole=goi;EvidenceType={ev};Confidence=LOW")
    return {"start": start, "end": end, "score": ident, "flanking_support": flank,
            "record": {"id": f"{parent}|g_b0_{mid}", "seq": seq},
            "gff": [f"chr1\t{ev}\tmRNA\t{start}\t{end}\t{ident}\t+\t.\t{attrs}"]}


def test_every_candidate_gets_bits_against_the_home_protein():
    cands = [_cand("l1", 100, 900, HOME), _cand("l1_extra1", 1000, 1300, HOME[40:120], ev="rescued_exon"),
             _cand("l2", 2000, 2600, HOME, parent="GOI_other")]
    n = isr.score_goi_candidates(cands, {"GOI_q": {"id": "GOI_q", "seq": HOME}})
    assert n == 2
    assert cands[0]["bits"] == pytest.approx(sw_bits(HOME, HOME))
    assert cands[1]["bits"] == pytest.approx(sw_bits(HOME, HOME[40:120]))
    assert cands[0]["gff"][0].endswith(f";Bits={cands[0]['bits']:.1f}")
    assert "Bits=" not in cands[2]["gff"][0] and "bits" not in cands[2]   # no home protein known


def test_choice_goes_by_bits_not_identity():
    strong = _cand("l1", 100, 900, HOME, ident=35.0)                  # lower identity, more bits
    weak = _cand("l2", 120, 400, _mutate(HOME[:90], 4), ident=62.0)
    cands = [weak, strong]
    isr.score_goi_candidates(cands, {"GOI_q": {"seq": HOME}})
    assert strong["bits"] > weak["bits"]
    cands.sort(key=isr.goi_candidate_order)
    assert cands[0] is strong


def test_equal_bits_the_shorter_span_wins():
    """A chain over two copies scores what its better copy scores: the copy must win."""
    copy = _cand("l1", 1000, 1600, HOME)
    chain = _cand("l1_fallback", 1000, 9000, HOME, ev="fallback_hit_span")
    cands = [chain, copy]
    isr.score_goi_candidates(cands, {"GOI_q": {"seq": HOME}})
    assert chain["bits"] == copy["bits"]
    cands.sort(key=isr.goi_candidate_order)
    assert cands[0] is copy


def test_order_then_flanking_support_start_and_id():
    a = dict(_cand("a", 500, 900, HOME), bits=80.0, flanking_support=3)
    b = dict(_cand("b", 100, 500, HOME), bits=80.0, flanking_support=7)
    c = dict(_cand("c", 50, 450, HOME), bits=80.0, flanking_support=7)
    assert sorted([a, b, c], key=isr.goi_candidate_order) == [c, b, a]


# ---------------------------------------------------------------- the report

def _rec(start, end, identity, bits="", confidence="HIGH", source="locus_1"):
    return {"genome": "g1", "source": source, "chrom": "chr1", "start": start, "end": end,
            "confidence": confidence, "goi_class": "confident_goi", "identity": str(identity),
            "query_cov": "0.9", "bits": bits, "target_gene": "", "model_status": "complete",
            "mrna_id": f"m{start}", "context_promoted": ""}


def test_report_keeps_the_label_then_the_bits():
    hi_ident = _rec(1000, 2000, 71.0, bits="95.0")
    hi_bits = _rec(1010, 1990, 48.0, bits="160.0")
    out = gr.dedupe_goi_annotations([hi_ident, hi_bits])
    assert len(out["records"]) == 1
    assert out["records"][0]["identity"] == "48.0" and out["records"][0]["bits"] == "160.0"


def test_report_label_still_comes_first():
    medium = _rec(1000, 2000, 30.0, bits="300.0", confidence="MEDIUM")
    high = _rec(1010, 1990, 60.0, bits="90.0")
    rep = gr.dedupe_goi_annotations([medium, high])["records"][0]
    assert rep["confidence"] == "HIGH" and rep["bits"] == "90.0"


def test_report_without_bits_falls_back_to_identity():
    """Runs written before 2026-10-14 carry no Bits: the old rule applies unchanged."""
    a = _rec(1000, 2000, 71.0)
    b = _rec(1010, 1990, 48.0)
    rep = gr.dedupe_goi_annotations([b, a])["records"][0]
    assert rep["identity"] == "71.0" and rep["bits"] == ""


def test_report_reads_bits_from_the_gff(tmp_path):
    gff = tmp_path / "locus_1__g1.fna.gff"
    gff.write_text("##gff-version 3\n" + "\t".join([
        "chr1", "exon_annotation", "mRNA", "100", "900", "50.0", "+", ".",
        "ID=GOI_q|g_b0_l1;Name=GOI_q;SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=HIGH;"
        "GOIClass=confident_goi;EvidenceType=exon_annotation;Identity=50.0;Bits=123.4"]) + "\n")
    recs = gr.collect_goi_annotations([str(gff)])
    assert recs and recs[0]["bits"] == "123.4"


# ---------------------------------------------------------------- the figures

def _gene(identity, bits=None, confidence="HIGH", start=100):
    return {"identity": identity, "bits": bits, "confidence": confidence, "start": start,
            "query_coverage": 0.9}


def test_figures_pick_the_model_with_more_bits():
    a, b = _gene(71.0, 95.0), _gene(48.0, 160.0, start=200)
    assert ps._goi_best([a, b]) is b
    assert ps._goi_best([_gene(71.0), _gene(48.0, start=200)])["identity"] == 71.0
    assert ps._goi_best([_gene(30.0, 300.0, "MEDIUM"), _gene(60.0, 90.0)])["confidence"] == "HIGH"


# ---------------------------------------------------------------- the hull rescue

def test_hull_rescue_writes_bits():
    mrna = ["w", "miniprot", "mRNA", "11", "910", ".", "+", ".", "ID=MP1;Identity=0.81"]
    cds = [["w", "miniprot", "CDS", "11", "910", ".", "+", "0", "Parent=MP1"]]
    rows = hull._build_hull_gff_rows(mrna, cds, "chr1", 1000, "q", "GOI_q|g_chr1_hull_rescue",
                                     "HIGH", query_coverage=0.95, bits=210.3)
    assert all("Bits=210.3" in r for r in rows[:2])
    rows = hull._build_hull_gff_rows(mrna, cds, "chr1", 1000, "q", "GOI_q|g_chr1_hull_rescue", "HIGH")
    assert not any("Bits=" in r for r in rows)
