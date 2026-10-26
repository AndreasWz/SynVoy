"""Whose model a seeded tandem copy is (label rework REPORT finding 6; fixed behind
``--tandem_copy_parent``, default ``seed`` until the measured switch of the preset plan).

``seed``: a single-segment tandem copy is recorded as
``GOI_copy_<n>|<genome>_b..._l..._by_<parent>``; once it seeds the next wave, ``GOI_copy_<n>``
is a GOI of its own and its models are judged against the seed. ``home``: it is recorded under
the home GOI and is one of its models.
"""
import os
import random
import sys

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice(AA) for _ in range(n))


def test_record_ids_of_the_two_modes():
    assert isr.tandem_copy_record_id("GOI_copy_2", "GOI_X", "GenA_fna", 3, 1) == \
        "GOI_copy_2|GenA_fna_b3_l1_by_GOI_X"
    assert isr.tandem_copy_record_id("GOI_copy_2", "GOI_X", "GenA_fna", 3, 1, "home") == \
        "GOI_X|GenA_fna_b3_l1_tc2"
    # a copy found with the relative query is renamed to the home GOI by the rebase step
    assert isr.tandem_copy_record_id("GOI_copy_1", "GOI_X~rel", "GenA_fna", 0, 2, "home") == \
        "GOI_X~rel|GenA_fna_b0_l2_tc1"


def _db(copy_id, copy_seq, home_seq):
    return {"GOI_X": {"id": "GOI_X", "seq": home_seq},
            copy_id: {"id": copy_id, "seq": copy_seq}}


def test_a_seeded_copy_is_a_goi_of_its_own_in_seed_mode():
    home = _protein(120, 1)
    copy = home[:40] + _protein(20, 2) + home[60:100]
    copy_id = isr.tandem_copy_record_id("GOI_copy_1", "GOI_X", "GenA_fna", 0, 1, "seed")
    queries, by_parent, alias_of = isr.select_goi_queries(_db(copy_id, copy, home), {"GOI_X"})
    assert "GOI_copy_1" in by_parent                      # REPORT finding 6
    assert by_parent["GOI_copy_1"]["seq"] == copy


def test_a_seeded_copy_is_a_model_of_the_home_goi_in_home_mode():
    home = _protein(120, 1)
    copy = home[:40] + _protein(20, 2) + home[60:100]
    copy_id = isr.tandem_copy_record_id("GOI_copy_1", "GOI_X", "GenA_fna", 0, 1, "home")
    queries, by_parent, alias_of = isr.select_goi_queries(_db(copy_id, copy, home), {"GOI_X"})
    assert not any(p.startswith("GOI_copy_") for p in by_parent)
    assert by_parent["GOI_X"]["seq"] == home               # the home protein stays the query
    assert alias_of == {copy_id: "GOI_X~rel"}               # the copy is the relative query
    assert isr._model_genome_tag(copy_id) == "GenA_fna"


def test_the_runner_takes_the_mode_from_its_option():
    src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    assert "getattr(args, 'tandem_copy_parent', 'seed')" in src
    assert '--tandem_copy_parent' in src
