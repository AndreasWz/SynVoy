"""Review fix RF-1f (batch of 2026-10-15): the judge and the seed bar use full precision.

The runner writes EPlace / Eeff with 3 significant digits and SpecMargin at 0.1 bit, and
used to read those strings back: an E of 1.0004e-6 read as 1e-06 (HIGH at the bar of 1e-6),
a margin of 9.96 bits as 10.0 (a seed at the bar of 10). The numbers now travel with the
candidate (judge_numbers) and per model id to the genome end; the strings are for the GFF.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402
import synvoy_label  # noqa: E402

ROW = ("c1\texon_annotation\tmRNA\t100\t900\t60.0\t+\t.\tID=GOI_q|g_b0_l1_exon_ann;SynVoy_Parent=GOI_q;"
       "SynVoyRole=goi;EvidenceType=exon_annotation;Confidence=LOW;GOIClass=x;InferenceReason=awaiting_label;"
       "QueryCoverage=0.900;Place=slot;BlockEsyn=1e-05;EPlace=1e-06;EPlaceTblastn=inf;"
       "BestOtherHome=gene-p;SpecMargin=10.0;RealCompetitor=1;BestOtherInNeighbourhood=0")
NUMBERS = {"EPlace": 1.0004e-6, "EPlaceTblastn": float("inf"), "BlockEsyn": 1e-5,
           "SpecMargin": 9.96, "RealCompetitor": True, "BestOtherInNeighbourhood": False}


def _attrs(line):
    return isr._parse_gff_attr_string(line.split("\t")[8])


def test_the_string_and_the_number_disagree_at_the_bars():
    assert synvoy_label.HIGH_EVALUE == 1e-6 and synvoy_label.HIGH_MARGIN_BITS == 10.0
    assert isr._judge_attrs(_attrs(ROW), "EPlace")[0] == "HIGH"                 # the rounded strings
    assert isr._judge_attrs(_attrs(ROW), "EPlace", dict(NUMBERS))[0] == "MEDIUM"
    assert isr._judge_attrs(_attrs(ROW), "EPlace", dict(NUMBERS, EPlace=0.9e-6))[0] == "MEDIUM"  # margin 9.96 < 10
    assert isr._judge_attrs(_attrs(ROW), "EPlace", dict(NUMBERS, EPlace=0.9e-6, SpecMargin=10.0))[0] == "HIGH"


def test_the_window_judge_uses_the_candidate_numbers():
    cand = {"start": 100, "end": 900, "gff": [ROW], "record": {"id": "GOI_q|g_b0_l1_exon_ann"},
            "judge_numbers": dict(NUMBERS)}
    isr.judge_candidates([cand])
    assert cand["rank_key"][0] == -synvoy_label.LABEL_RANK["MEDIUM"]


def test_the_genome_end_computes_eeff_from_the_numbers_and_judges_on_them():
    ctx = {"block_esyn": {0: 1e-5}, "judge_numbers": {"GOI_q|g_b0_l1_exon_ann": dict(NUMBERS)}}
    lines = isr.finalize_place_evalues([ROW], ctx)
    numbers = ctx["judge_numbers"]["GOI_q|g_b0_l1_exon_ann"]
    assert numbers["Eeff"] == 1.0004e-6 and "Eeff=1e-06" in lines[0]          # the string is rounded
    out, counts = isr.relabel_goi_lines(lines, ctx["judge_numbers"])
    assert "Confidence=MEDIUM" in out[0] and counts == {"MEDIUM/significant_specific": 1}
    out_strings, _ = isr.relabel_goi_lines(lines)
    assert "Confidence=HIGH" in out_strings[0]                                   # what it was before


def test_the_seed_bar_reads_the_full_margin():
    from types import SimpleNamespace
    args = SimpleNamespace(genome_ctx={"judge_numbers": {"m1": {"SpecMargin": 9.96, "RealCompetitor": True}}})
    n = isr._numbers_of(args, "m1")
    assert not isr._seed_is_specific({"spec_margin": n["SpecMargin"], "real_competitor": n["RealCompetitor"]})
    assert isr._seed_is_specific({"spec_margin": "10.0", "real_competitor": "1"})      # the old string
    assert isr._numbers_of(args, "m2") == {}
