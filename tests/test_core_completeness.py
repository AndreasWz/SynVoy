"""B3 signal_anchor, core completeness for HIGH (preset plan of 2026-10-15, PR-10).

With core_completeness a GOI gene model carries CoreCoverage: the share of the home
protein's core (the protein without its signal peptide and GPI tail) that the model
covers; the judge then reads the larger of QueryCoverage and CoreCoverage. The judge's
constant (0.7) does not change; MEDIUM needs no completeness and is unaffected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import iterative_search_runner as isr  # noqa: E402
import resolve_effective_params as rep  # noqa: E402

# UniProt Q16553 (LY6E): signal peptide 1-20, GPI tail of 32 aa (query_features).
LY6E = ("MKIFLPVLLAALLGVERASSLMCFSCLNQKSNLYCLKPTICSDQDNYCVTVSASAGIGNLVTFGHSLSKTCSPACPIPEGVNVG"
        "VASMGISCCQSFLCNFSAADGGLRASVTLLGAGLLLSLLPALLRFGP")


def test_the_core_of_ly6e_leaves_out_its_signal_peptide_and_gpi_tail():
    assert isr.query_core(LY6E) == (21, len(LY6E) - 32)
    assert isr.query_core("ACDEFGHIKLMNPQRSTVWY" * 5) == (1, 100)       # neither: the whole protein


def test_core_coverage_is_a_union_clipped_to_the_core():
    first, last = isr.query_core(LY6E)
    core_len = last - first + 1
    assert isr.core_coverage([(first, last)], LY6E) == 1.0
    assert isr.core_coverage([(1, 60)], LY6E) == pytest.approx((60 - first + 1) / core_len)
    assert isr.core_coverage([(30, 50), (40, 70), (1, 5)], LY6E) == pytest.approx(41 / core_len)
    assert isr.core_coverage([(1, 10)], LY6E) == 0.0
    assert isr.core_coverage([(1, 10)], "") is None


ATTRS = {"Eeff": "1e-20", "Place": "slot", "BlockEsyn": "1e-4", "SpecMargin": "50",
         "RealCompetitor": "0", "BestOtherInNeighbourhood": "0", "QueryCoverage": "0.600",
         "CoreCoverage": "1.000"}


def test_the_judge_reads_the_core_only_with_the_flag(monkeypatch):
    monkeypatch.setattr(isr, "CORE_COMPLETENESS", False)
    assert isr._judge_attrs(dict(ATTRS), "Eeff") == ("MEDIUM", "significant_specific")
    monkeypatch.setattr(isr, "CORE_COMPLETENESS", True)
    assert isr._judge_attrs(dict(ATTRS), "Eeff") == ("HIGH", "significant_specific")
    # the larger of the two: a model complete on the query stays complete
    assert isr._judge_attrs(dict(ATTRS, QueryCoverage="0.900", CoreCoverage="0.500"), "Eeff")[0] == "HIGH"
    # without the attribute nothing changes
    plain = {k: v for k, v in ATTRS.items() if k != "CoreCoverage"}
    assert isr._judge_attrs(plain, "Eeff")[0] == "MEDIUM"


def test_core_coverage_is_written_only_when_given():
    with_core = isr._goi_feature_attrs({"ID": "x"}, "exon_annotation", query_cov=0.6, core_cov=0.95)
    assert with_core["CoreCoverage"] == "0.950"
    assert "CoreCoverage" not in isr._goi_feature_attrs({"ID": "x"}, "exon_annotation", query_cov=0.6)


def test_a_model_found_with_the_relative_gets_core_coverage_from_the_home_alignment():
    pytest.importorskip("parasail")
    first, last = isr.query_core(LY6E)
    protein = LY6E[first - 1:last]                      # exactly the core
    gff = ("c1\texon_annotation\tmRNA\t1\t300\t60.0\t+\t.\tID=GOI_L~rel|g_b0_l1_exon_ann;"
           "Name=GOI_L~rel;SynVoy_Parent=GOI_L~rel;Identity=60.0;SynVoyRole=goi;"
           "EvidenceType=exon_annotation;Exons=3;BlockFlankingSupport=5")
    cand = {"start": 1, "end": 300, "gff": [gff], "record": {"id": "GOI_L~rel|g_b0_l1_exon_ann", "seq": protein}}
    by_parent = {"GOI_L": {"id": "GOI_L", "seq": LY6E}, "GOI_L~rel": {"id": "GOI_L~rel", "seq": protein}}
    out = isr.rebase_relative_candidates([cand], ["GOI_L~rel"], by_parent, core=True)
    attrs = isr._parse_gff_attr_string(out[0]["gff"][0].split("\t")[8])
    assert float(attrs["QueryCoverage"]) == pytest.approx((last - first + 1) / len(LY6E), abs=0.01)
    assert float(attrs["CoreCoverage"]) == pytest.approx(1.0, abs=0.01)
    plain = isr.rebase_relative_candidates([cand], ["GOI_L~rel"], by_parent)
    assert "CoreCoverage" not in plain[0]["gff"][0]


def test_b3_sets_the_key_and_legacy_does_not():
    defaults = {"core_completeness": False}
    f = {"sp": True, "sp_len": 20}
    assert rep.resolve(defaults, auto_preset=None, modules_spec="auto", features=f)["settings"] == \
        {"core_completeness": True}
    assert rep.resolve(defaults, auto_preset={"preset": "preset_single_copy"}, modules_spec="legacy",
                       features=f)["settings"] == defaults


def test_the_flag_is_set_for_each_genome_task():
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    body = src[src.index("def process_single_genome"):]
    assert "CORE_COMPLETENESS = bool(getattr(args, 'core_completeness', False))" in body[:6000]
    assert "core=bool(getattr(args, 'core_completeness', False))" in src
