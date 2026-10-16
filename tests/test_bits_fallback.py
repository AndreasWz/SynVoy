"""Review fix RF-3 (S01-11): without parasail every GOI candidate type is ordered by the same
fallback, the best bit score of its own search hits. Before, a hit chain kept the bits of its
best hit while every other type had 0, so a chain won every choice."""
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402


def _cand(evidence, start, end, bits=None):
    rid = f"GOI_q|g_b0_{evidence}"
    row = (f"c1\tsynvoy\tmRNA\t{start}\t{end}\t50.0\t+\t.\tID={rid};Name=GOI_q;"
           f"SynVoy_Parent=GOI_q;SynVoyRole=goi;EvidenceType={evidence}")
    c = {"start": start, "end": end, "gff": [row], "record": {"id": rid, "seq": "MKV" * 20},
         "label_protein": "MKV" * 20}
    if bits is not None:
        c["bits"] = bits
    return c


def test_every_type_falls_back_to_its_own_hits(monkeypatch, caplog):
    def no_parasail(_a, _b):
        raise RuntimeError("parasail-python is required for Smith-Waterman alignment")
    monkeypatch.setattr(isr, "sw_bits", no_parasail)
    monkeypatch.setattr(isr, "_SW_BITS_FAILED", [False])
    chain = _cand("fallback_hit_span", 100, 400, bits=45.0)    # created with its best hit's bits
    model = _cand("exon_annotation", 1000, 1600)
    copy = _cand("tandem_copy", 3000, 3300)                     # no hit of its own in the span
    hits = {"GOI_q": [{"gstart": 99, "gend": 399, "strand": "+", "bits": 45.0},
                      {"gstart": 999, "gend": 1599, "strand": "+", "bits": 80.0},
                      {"gstart": 1200, "gend": 1500, "strand": "+", "bits": 30.0}]}
    cands = [chain, model, copy]
    isr.prepare_label_subjects(cands, hits, "ATG" * 2000, 0)
    with caplog.at_level(logging.WARNING):
        isr.score_goi_candidates(cands, {"GOI_q": {"seq": "MKV" * 20}})
    assert (chain["bits"], model["bits"], copy["bits"]) == (45.0, 80.0, 0.0)
    assert all(c["bits_from_hits"] for c in cands)
    assert not any("Bits=" in c["gff"][0] for c in cands)
    assert sorted(cands, key=isr.goi_candidate_order)[0] is model
    assert sum("ordered by the best bit score of their own search hits" in r.getMessage()
               for r in caplog.records) == 1

    # No specificity terms from bits of another unit.
    def must_not_run(*_a, **_k):
        raise AssertionError("specificity computed from fallback bits")
    monkeypatch.setattr(isr, "_home_proteome", lambda _args: (object(), None))
    monkeypatch.setattr(isr.synvoy_specificity, "specificity", must_not_run)
    assert isr.annotate_goi_specificity(cands, {"GOI_q": {"seq": "MKV" * 20}}, None) == 0


def test_with_parasail_the_bits_are_smith_waterman_bits():
    chain = _cand("fallback_hit_span", 100, 400, bits=45.0)
    hits = {"GOI_q": [{"gstart": 99, "gend": 399, "strand": "+", "bits": 45.0}]}
    isr.prepare_label_subjects([chain], hits, "ATG" * 2000, 0)
    isr.score_goi_candidates([chain], {"GOI_q": {"seq": "MKV" * 20}})
    assert chain["bits"] == isr.sw_bits("MKV" * 20, "MKV" * 20) and "bits_from_hits" not in chain
    assert f"Bits={chain['bits']:.1f}" in chain["gff"][0]
