"""Review fixes RF-3 (batch of 2026-10-15), decoys and the tblastn E source.

* F3: each decoy case (home GOI x kind) is chosen among its own candidates; pooled, an
  overlapping decoy call of another kind suppressed it (the rate was a lower bound).
* F6: a decoy row carries PlaceBp and its E-value in every nested place that holds it.
* F18: no home locus span (goi=) means no slot: the runner says so.
* F5: a tblastn row counts for EPlaceTblastn when it overlaps the candidate in its frame,
  not only when it lies wholly inside.
"""
import os
import sys

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402


def _cand(parent, start, end):
    row = (f"c1\texon_annotation\tmRNA\t{start}\t{end}\t.\t+\t.\tID={parent}|g_b0_x{start};"
           f"SynVoy_Parent={parent};SynVoyRole=goi;EvidenceType=exon_annotation")
    return {"start": start, "end": end, "gff": [row], "record": {"id": f"{parent}|g_b0_x{start}", "seq": "MK"}}


def test_each_decoy_case_is_its_own_choice():
    cands = [_cand("GOI_q~shuf1", 100, 900), _cand("GOI_q~rev", 120, 880), _cand("GOI_q~shuf1", 2000, 2500)]
    cases = isr.decoy_cases(cands)
    assert [[c["record"]["id"] for c in case] for case in cases] == [
        ["GOI_q~rev|g_b0_x120"], ["GOI_q~shuf1|g_b0_x100", "GOI_q~shuf1|g_b0_x2000"]]
    src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    assert "for case in decoy_cases(decoy_candidates):" in src


def test_nested_places_of_a_decoy_row():
    terms = {"significant": True, "slot": (1000, 2000), "core": (500, 5000), "window": (0, 10000),
             "genome_len": 1_000_000}
    q = 100
    attrs = dict(kv.split("=") for kv in isr.decoy_place_attrs((1200, 1800), 1000, 1e-6, q, terms).strip(";").split(";"))
    assert attrs["PlaceBp"] == "1000" and attrs["EPlaceSlot"] == "1e-06"
    ratio = lambda bp: synvoy_evalue.search_space(q, bp) / synvoy_evalue.search_space(q, 1000)
    for key, bp in (("EPlaceCore", 4500), ("EPlaceWindow", 10000), ("EPlaceGenome", 1_000_000)):
        assert float(attrs[key]) == float(f"{1e-6 * ratio(bp):.3g}"), key
    outside = isr.decoy_place_attrs((3000, 3500), 4500, 1e-6, q, terms)
    assert "EPlaceSlot" not in outside and "EPlaceCore" in outside
    insignificant = isr.decoy_place_attrs((1200, 1800), 1_000_000, 1e-6, q, dict(terms, significant=False))
    assert "EPlaceCore" not in insignificant and "EPlaceGenome=1e-06" in insignificant


def test_a_missing_home_locus_span_is_announced():
    src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    assert "if not args.goi_ranks:" in src and "the slot of every block is disabled" in src


def test_a_tblastn_row_reaching_past_the_candidate_counts(monkeypatch):
    monkeypatch.setattr(synvoy_evalue, "blastp_evalue", lambda q, s, sp: (float("inf"), 0.0))
    cand = _cand("GOI_q", 1001, 1600)
    cand["label_subjects"] = []
    terms = {"significant": True, "slot": None, "core": (0, 10000), "window": (0, 10000),
             "genome_len": 10 ** 6, "esyn": 1e-3, "window_len": 10000}
    rows = [{"query": "GOI_q", "start": 900, "end": 1300, "strand": "+", "evalue": 1e-8}]   # starts before
    isr.annotate_place_evalues([cand], {"GOI_q": {"seq": "M" * 100}}, {}, terms, "A" * 10000, 0, rows)
    assert "EPlaceTblastn=inf" not in cand["gff"][0]
    away = _cand("GOI_q", 5001, 5600)
    away["label_subjects"] = []
    isr.annotate_place_evalues([away], {"GOI_q": {"seq": "M" * 100}}, {}, terms, "A" * 10000, 0, rows)
    assert "EPlaceTblastn=inf" in away["gff"][0]
