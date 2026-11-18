"""Decoy shuffles are seeded by the decoy kind and the protein sequence, not by the gene's name:
a renamed gene with the same protein got other ~shuf1/~shuf2 decoys, and the
relative's decoys were one fixed permutation per name applied to whatever relative a window had.
Real calls never depended on it (decoys are written apart). All synthetic."""
import os
import sys

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
OTHER = "MSDNELQKLLAEAKRLGLSQEEFVKKAWELLGEPAPAAPAAGTSSAAEDKKAAA"


def _by_kind(decoys):
    return {d["id"].split("~", 1)[1]: d for d in decoys}


def test_the_same_protein_under_two_names_gets_the_same_decoys():
    a = _by_kind(isr.decoy_queries([{"id": "GOI_alpha", "seq": SEQ}], {"GOI_alpha"}))
    b = _by_kind(isr.decoy_queries([{"id": "GOI_beta", "seq": SEQ}], {"GOI_beta"}))
    assert [a[k]["seq"] for k in ("rev", "shuf1", "shuf2")] == [b[k]["seq"] for k in ("rev", "shuf1", "shuf2")]
    assert a["shuf1"]["id"] == "GOI_alpha~shuf1" and b["shuf1"]["id"] == "GOI_beta~shuf1"


def test_the_kinds_differ_and_keep_the_composition():
    d = _by_kind(isr.decoy_queries([{"id": "GOI_q", "seq": SEQ}], {"GOI_q"}))
    assert d["rev"]["seq"] == SEQ[::-1]
    assert d["shuf1"]["seq"] != d["shuf2"]["seq"]
    assert sorted(d["shuf1"]["seq"]) == sorted(SEQ) == sorted(d["shuf2"]["seq"])
    assert d["shuf1"]["seq"] != SEQ


def test_a_rerun_gives_the_same_decoys():
    first = isr.decoy_queries([{"id": "GOI_q", "seq": SEQ}], {"GOI_q"})
    again = isr.decoy_queries([{"id": "GOI_q", "seq": SEQ}], {"GOI_q"})
    assert first == again


def test_the_relatives_shuffles_follow_its_own_sequence():
    def rel(name, seq):
        return _by_kind(isr.decoy_queries([{"id": f"{name}~rel", "seq": seq}], {name}))
    one, two = rel("GOI_q", SEQ), rel("GOI_q", OTHER)
    assert one["shuf1~rel"]["seq"] != two["shuf1~rel"]["seq"]
    renamed = rel("GOI_r", SEQ)
    assert renamed["shuf1~rel"]["seq"] == one["shuf1~rel"]["seq"]
    assert renamed["shuf1~rel"]["id"] == "GOI_r~shuf1~rel"


def test_the_relatives_tag_is_in_the_seed():
    home = _by_kind(isr.decoy_queries([{"id": "GOI_q", "seq": SEQ}], {"GOI_q"}))
    relative = _by_kind(isr.decoy_queries([{"id": "GOI_q~rel", "seq": SEQ}], {"GOI_q"}))
    assert home["shuf1"]["seq"] != relative["shuf1~rel"]["seq"]
    assert home["rev"]["seq"] == relative["rev~rel"]["seq"]
