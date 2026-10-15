"""bin/synvoy_label.py: the one judging function of the label rework (plan section 3.4),
its ranking and the frame guard (section 3.1)."""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import synvoy_label as L  # noqa: E402


def j(e, place="slot", sig=True, margin=50.0, real=True, near=False, qcov=0.9):
    return L.judge(e, place, sig, margin, real, near, qcov)


def test_the_five_rows_of_the_table():
    assert j(1e-20) == ("HIGH", "significant_specific")
    assert j(1e-4) == ("MEDIUM", "significant_specific")
    assert j(1e-20, margin=-20.0) == ("AMBIGUOUS", "not_specific")
    assert j(0.05) == ("AMBIGUOUS", "gene_order_only")
    assert j(5.0) == ("LOW", "not_significant")


def test_high_needs_the_10_bit_margin_and_a_complete_model():
    assert j(1e-20, margin=9.9) == ("MEDIUM", "significant_specific")
    assert j(1e-20, margin=-50.0, real=False) == ("HIGH", "significant_specific")   # no real competitor
    assert j(1e-20, qcov=0.69) == ("MEDIUM", "significant_specific")
    assert j(1e-20, qcov=None) == ("MEDIUM", "significant_specific")


def test_thresholds_are_inclusive():
    assert j(1e-6)[0] == "HIGH" and j(1e-3)[0] == "MEDIUM" and j(1.0)[0] == "AMBIGUOUS"


def test_gene_order_only_needs_the_slot_of_a_significant_block():
    assert j(0.05, place="core")[0] == "LOW"
    assert j(0.05, sig=False)[0] == "LOW"
    assert j(0.05, margin=-10.0) == ("LOW", "not_significant")       # not specific


def test_positional_tie_break_and_its_exception():
    # within +-3.3 bits of a real competitor: position decides in slot/core of a significant block
    assert j(1e-5, margin=1.0, place="slot") == ("MEDIUM", "significant_specific")
    assert j(1e-5, margin=-3.0, place="core") == ("MEDIUM", "significant_specific")
    assert j(1e-5, margin=1.0, place="window") == ("AMBIGUOUS", "not_specific")
    assert j(1e-5, margin=1.0, sig=False, place="genome") == ("AMBIGUOUS", "not_specific")
    # ... unless the competitor lies in the gene's own home neighbourhood (or that is unknown)
    assert j(1e-5, margin=1.0, near=True) == ("AMBIGUOUS", "not_specific")
    assert j(1e-5, margin=1.0, near=None) == ("AMBIGUOUS", "not_specific")
    # a margin of >= 3.3 bits is specific wherever the competitor lies
    assert j(1e-5, margin=3.3, near=True) == ("MEDIUM", "significant_specific")
    assert j(1e-5, margin=-3.3, place="slot") == ("AMBIGUOUS", "not_specific")


def test_without_specificity_terms_nothing_is_held_back():
    assert L.judge(1e-20, "slot", True, None, None, None, 0.9) == ("HIGH", "significant_specific")


def test_missing_e_is_not_significant():
    assert j(None) == ("LOW", "not_significant") and j(float("nan"))[0] == "LOW"


def test_rank_key():
    k = L.rank_key
    assert k("MEDIUM", 1e-4, 1e-3, 30, 900, 5) < k("AMBIGUOUS", 1e-30, 1e-30, 300, 900, 5)   # label first
    assert k("HIGH", 1e-20, 1e-3, 30, 900, 5) < k("HIGH", 1e-8, 1e-30, 300, 900, 5)          # then E
    # below 1e-10 the E-values are equal: the stronger block, then bits, then the shorter span
    assert k("HIGH", 1e-50, 1e-5, 30, 900, 5) > k("HIGH", 1e-12, 1e-9, 30, 900, 5)
    assert k("HIGH", 1e-50, 1e-9, 30, 900, 5) > k("HIGH", 1e-12, 1e-9, 99, 900, 5)
    assert k("HIGH", 1e-50, 1e-9, 99, 2000, 5) > k("HIGH", 1e-12, 1e-9, 99, 900, 5)


def test_frame_guard():
    assert L.frame_guard("ATGGCTGCTTAA") == ("MAA", "ok")
    assert L.frame_guard("CATGGCTGCTGCT") == ("MAAA", "fixed")        # one extra base in front
    assert L.frame_guard("TAAATAAATAAAGGGCCC")[1] == "broken"         # a stop inside every frame
    assert L.frame_guard("atggctgct") == ("MAA", "ok")
