"""Decoys of the relative's query (review fix RF-3, decided 2026-10-15): in decoy mode the
second GOI query of a window (``<parent>~rel``, the model of a relative searched earlier) gets
a reversed and two shuffled copies too, ids ``<parent>~rev~rel`` etc. Like a real ``~rel``
candidate, a model they find is measured against the home decoy of the same kind and filed
under it, so it competes in that decoy case and is written to decoys/ with SearchQuery=."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402

HOME = "MKVLAAGIVLLWSTRPQEDCHNFY" * 3
REL = "MKVLSAGLVLLWSTKPQEDCHNWY" * 3


def _queries():
    return [{"id": "GOI_q", "seq": HOME},
            {"id": "GOI_q~rel", "seq": REL, "model_id": "GOI_q|G_fna_b0_l1"},
            {"id": "GOI_copy_1|G_fna_b0_l1", "seq": "MKKKKKKKKKKKKKKKK"}]


def test_relative_query_gets_three_decoys_home_decoys_unchanged():
    decoys = isr.decoy_queries(_queries(), {"GOI_q"})
    ids = [d["id"] for d in decoys]
    assert ids == ["GOI_q~rev", "GOI_q~shuf1", "GOI_q~shuf2",
                   "GOI_q~rev~rel", "GOI_q~shuf1~rel", "GOI_q~shuf2~rel"]
    by_id = {d["id"]: d["seq"] for d in decoys}
    assert by_id["GOI_q~rev~rel"] == REL[::-1] and sorted(by_id["GOI_q~shuf1~rel"]) == sorted(REL)
    # The home decoys are those of the code before (same ids, same seeds).
    home_only = isr.decoy_queries([_queries()[0]], {"GOI_q"})
    assert [(d["id"], d["seq"]) for d in home_only] == [(d["id"], d["seq"]) for d in decoys[:3]]
    # A seeded copy's own query still gets none (its hits would count as real).
    assert not any(d["id"].startswith("GOI_copy_1") for d in decoys)
    assert all(isr.is_decoy_parent(isr.extract_base_gene_id(i)) for i in ids)
    assert not isr.is_decoy_parent("GOI_q~rel") and not isr.is_decoy_parent("GOI_q")


def test_a_relative_decoy_model_is_filed_under_the_home_decoy():
    decoys = isr.decoy_queries(_queries(), {"GOI_q"})
    by_parent = {d["id"]: d for d in decoys}
    protein = REL[::-1][5:60]
    row = ("c1\tminiprot\tmRNA\t1000\t1165\t80.0\t+\t.\tID=GOI_q~rev~rel|G2_fna_b0_l1;"
           "Name=GOI_q~rev~rel;SynVoy_Parent=GOI_q~rev~rel;Type=exon_annotation;Identity=80.0;"
           "SynVoyRole=goi;EvidenceType=exon_annotation;Confidence=LOW;Exons=1")
    cand = {"start": 1000, "end": 1165, "gff": [row],
            "record": {"id": "GOI_q~rev~rel|G2_fna_b0_l1", "seq": protein, "description": ""}}
    out = isr.rebase_relative_candidates([cand], ["GOI_q~rev~rel"], by_parent)
    assert len(out) == 1
    attrs = isr._parse_gff_attr_string(out[0]["gff"][0].split("\t")[8])
    assert attrs["SynVoy_Parent"] == "GOI_q~rev" and attrs["ID"] == "GOI_q~rev|G2_fna_b0_r1"
    assert attrs["SearchQuery"] == "GOI_q~rev~rel"
    assert isr.is_decoy_parent(isr._candidate_parent(out[0]))
    # Same decoy case as a model of the home decoy of that kind.
    home_row = row.replace("GOI_q~rev~rel", "GOI_q~rev").replace("1000\t1165", "5000\t5165")
    home_cand = {"start": 5000, "end": 5165, "gff": [home_row], "record": {"id": "x", "seq": protein}}
    cases = isr.decoy_cases([out[0], home_cand])
    assert len(cases) == 1 and len(cases[0]) == 2
    # Where the home decoy has a model of its own on the same DNA, the relative's is dropped.
    over = dict(home_cand, start=1010, end=1160,
                gff=[home_row.replace("5000\t5165", "1010\t1160")])
    assert isr.rebase_relative_candidates([over, cand], ["GOI_q~rev", "GOI_q~rev~rel"],
                                          by_parent) == [over]
