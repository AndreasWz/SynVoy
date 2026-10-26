"""One id per tandem copy (2026-10-23).

detect_tandem_duplications counts its copies ``GOI_copy_<n>`` from 1 for each parent, block and
locus, and the seed form of the record id was ``GOI_copy_<n>|<genome>_b<k>_l<j>``: copies of
different parents (the home GOI, a seeded copy, a home tandem copy, a decoy) at the same block,
locus and n shared one id. The judge's full-precision terms are kept per id until the genome's
end (judge_numbers), so Eeff, the label and the seed rule of such a copy came from whichever was
judged last: the real copies are judged after the decoys, so a decoy copy could be labelled
with a real copy's E-value. Now ``GOI_copy_<n>|<genome>_b<k>_l<j>_by_<parent>``; the base id
stays ``GOI_copy_<n>`` (one pool of next-wave GOIs, as before).
"""
import logging
import math
import os
import sys
import types

import pytest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import synvoy_tree  # noqa: E402
from sequence_utils import block_index_of  # noqa: E402

PARENTS = ["GOI_X", "GOI_orf19.3940.1", "GOI_copy_2", "GOI_X~rev", "GOI_X~shuf1"]
TAGS = [("GenA_fna", "GenA"), ("Candida_b2_fna", "Candida_b2")]
HOME = "MKVLAAGIVLLWSTRPQEDCHNFY" * 3
REL = "MKVLSAGLVLLWSTKPQEDCHNWY" * 3


def _old_seed_id(copy_name, tag, block, locus):
    """The seed form before 2026-10-23 (for the mutation checks)."""
    return f"{copy_name}|{tag}_b{block}_l{locus}"


def _copy_row(model_id, parent, start=1000, end=1299, terms=""):
    attrs = isr._goi_feature_attrs(
        {"ID": model_id, "Name": model_id.split("|")[0], "SynVoy_Parent": parent, "Type": "tandem_copy"},
        evidence_type="tandem_copy", identity=80.0, exon_count=1, query_cov=0.9, flanking_support=6)
    return f"c1\ttandem_copy\tgene\t{start}\t{end}\t80.0\t+\t.\t{isr._compose_gff_attrs(attrs)}{terms}"


@pytest.mark.parametrize("tag, species", TAGS)
def test_copies_of_five_parents_at_one_place_have_five_ids(tag, species):
    ids = [isr.tandem_copy_record_id("GOI_copy_1", p, tag, 3, 1) for p in PARENTS]
    assert len(set(ids)) == 5
    assert len({_old_seed_id("GOI_copy_1", tag, 3, 1) for _p in PARENTS}) == 1     # the bug
    for parent, model_id in zip(PARENTS, ids):
        assert model_id == f"GOI_copy_1|{tag}_b3_l1_by_{parent}"
        assert isr.extract_base_gene_id(model_id) == "GOI_copy_1"
        assert block_index_of(model_id, tag) == 3
        assert isr._model_genome_tag(model_id) == tag
        assert synvoy_tree.species_from_leaf(model_id) == species


@pytest.mark.parametrize("tag, species", TAGS)
def test_species_of_the_home_copy_form_and_of_a_numbered_fallback(tag, species):
    home = isr.tandem_copy_record_id("GOI_copy_2", "GOI_orf19.3940.1", tag, 0, 4, "home")
    assert home == f"GOI_orf19.3940.1|{tag}_b0_l4_tc2"
    for model_id in (home, f"GOI_X|{tag}_b0_l1_fallback_3", f"GOI_X|{tag}_b5_fl2_flank_hits",
                     _old_seed_id("GOI_copy_1", tag, 0, 1), f"GOI_X|{tag}_b0_l1m2_exon_ann"):
        assert synvoy_tree.species_from_leaf(model_id) == species
    assert block_index_of(home, tag) == 0 and isr._model_genome_tag(home) == tag


def test_a_parent_that_would_mislead_the_readers_is_logged_once(caplog):
    isr._WARNED_COPY_PARENTS.discard("GOI_b1_l2")
    with caplog.at_level(logging.WARNING, logger=isr.logger.name):
        model_id = isr.tandem_copy_record_id("GOI_copy_1", "GOI_b1_l2", "GenA_fna", 0, 1)
        isr.tandem_copy_record_id("GOI_copy_2", "GOI_b1_l2", "GenA_fna", 0, 1)
        isr.tandem_copy_record_id("GOI_copy_1", "GOI_X", "GenA_fna", 0, 1)
    assert [r.getMessage() for r in caplog.records].count(
        "tandem copy parent 'GOI_b1_l2' holds whitespace, '|' or a block-like token: the genome, "
        "block or species read from its copies' ids can be wrong") == 1
    assert len(caplog.records) == 1
    assert isr._model_genome_tag(model_id) != "GenA_fna"        # why it is logged


def test_copies_of_a_relatives_query_are_renamed_to_the_home_parent():
    by_parent = {"GOI_X": {"id": "GOI_X", "seq": HOME},
                 "GOI_X~rel": {"id": "GOI_X~rel", "seq": REL, "model_id": "GOI_X|GenB_fna_b0_l1_exon_ann"},
                 "GOI_X~rev": {"id": "GOI_X~rev", "seq": HOME[::-1]},
                 "GOI_X~rev~rel": {"id": "GOI_X~rev~rel", "seq": REL[::-1]}}
    out = {}
    for parent, home_parent, protein in (("GOI_X~rel", "GOI_X", REL[5:60]),
                                         ("GOI_X~rev~rel", "GOI_X~rev", REL[::-1][5:60])):
        copy_id = isr.tandem_copy_record_id("GOI_copy_1", parent, "GenA_fna", 2, 3)
        assert copy_id == f"GOI_copy_1|GenA_fna_b2_l3_by_{parent}"
        cand = {"start": 1000, "end": 1165, "gff": [_copy_row(copy_id, parent, 1000, 1165)],
                "record": {"id": copy_id, "seq": protein, "description": ""}}
        (rebased,) = isr.rebase_relative_candidates([cand], [parent], by_parent)
        a = isr._parse_gff_attr_string(rebased["gff"][0].split("\t")[8])
        assert a["ID"] == rebased["record"]["id"] == f"GOI_copy_1|GenA_fna_b2_r3_by_{home_parent}"
        assert a["SynVoy_Parent"] == home_parent and a["Name"] == "GOI_copy_1"
        assert isr.extract_base_gene_id(a["ID"]) == "GOI_copy_1"
        assert block_index_of(a["ID"], "GenA_fna") == 2 and isr._model_genome_tag(a["ID"]) == "GenA_fna"
        assert synvoy_tree.species_from_leaf(a["ID"]) == "GenA"
        out[home_parent] = a["ID"]
    assert out["GOI_X"] != out["GOI_X~rev"]     # a real and a decoy relative's copy


def test_copies_of_two_parents_are_one_goi_of_the_next_wave():
    a = isr.tandem_copy_record_id("GOI_copy_1", "GOI_X", "GenA_fna", 0, 1)
    b = isr.tandem_copy_record_id("GOI_copy_1", "GOI_copy_2", "GenA_fna", 0, 1)
    db = {"GOI_X": {"id": "GOI_X", "seq": HOME},
          a: {"id": a, "seq": HOME[:50]}, b: {"id": b, "seq": HOME[:60]}}
    queries, by_parent, alias_of = isr.select_goi_queries(db, {"GOI_X"})
    assert sorted(by_parent) == ["GOI_X", "GOI_copy_1"]
    assert by_parent["GOI_copy_1"]["id"] == b                     # the longest of the pool
    assert sorted(q["id"] for q in queries) == ["GOI_X", b] and alias_of == {}
    collapsed = isr.collapse_goi_queries_by_parent([db[a], dict(db[b], seq=HOME[:50])])
    assert list(collapsed) == ["GOI_copy_1"] and collapsed["GOI_copy_1"]["id"] == a   # tie: smallest id


# ------------------------------------------------------------- the judge-number table
def _judged(model_id, parent, eplace):
    row = _copy_row(model_id, parent, terms=f";Place=slot;BlockEsyn=0.0001;EPlace={eplace:.3g};"
                                            "EPlaceTblastn=inf;SpecMargin=40.0;RealCompetitor=0")
    numbers = {"EPlace": eplace, "EPlaceTblastn": math.inf, "BlockEsyn": 1e-4,
               "SpecMargin": 40.0, "RealCompetitor": False, "BestOtherInNeighbourhood": None}
    return {"start": 1000, "end": 1299, "gff": [row], "record": {"id": model_id, "seq": "MK"},
            "judge_numbers": numbers}


def _genome_end(cands_in_order, ctx):
    """record_judge_numbers per window in the runner's order (decoy cases first, the real
    candidates last), then each row through finalize_place_evalues and relabel_goi_lines."""
    numbers, writers = ctx.setdefault("judge_numbers", {}), {}
    repeats = sum(isr.record_judge_numbers(numbers, [c], writers) for c in cands_in_order)
    rows = []
    for c in cands_in_order:
        lines, _ = isr.relabel_goi_lines(isr.finalize_place_evalues(c["gff"], ctx), numbers)
        rows.append(isr._parse_gff_attr_string(lines[0].split("\t")[8]))
    return repeats, writers, rows


def test_a_real_copy_and_a_decoy_copy_keep_their_own_numbers_and_labels():
    real_id = isr.tandem_copy_record_id("GOI_copy_1", "GOI_X", "GenA_fna", 0, 1)
    decoy_id = isr.tandem_copy_record_id("GOI_copy_1", "GOI_X~rev", "GenA_fna", 0, 1)
    decoy, real = _judged(decoy_id, "GOI_X~rev", 5.0), _judged(real_id, "GOI_X", 1e-8)
    ctx = {"block_esyn": {0: 1e-4}}
    repeats, writers, (d, r) = _genome_end([decoy, real], ctx)
    assert repeats == 0 and writers == {decoy_id: ["GOI_X~rev"], real_id: ["GOI_X"]}
    assert ctx["judge_numbers"][decoy_id] is decoy["judge_numbers"]
    assert ctx["judge_numbers"][real_id] is real["judge_numbers"]
    assert (d["Eeff"], d["Confidence"], d["Basis"]) == ("5", "LOW", "not_significant")
    assert (r["Eeff"], r["Confidence"], r["Basis"]) == ("1e-08", "HIGH", "significant_specific")
    assert decoy["judge_numbers"]["Eeff"] == 5.0 and real["judge_numbers"]["Eeff"] == 1e-8
    # writing a candidate's own entry again is no repeat
    assert isr.record_judge_numbers(ctx["judge_numbers"], [real]) == 0


def test_with_the_old_seed_ids_the_decoy_copy_read_the_real_copys_numbers():
    """Mutation check: the old seed return gives both copies one id, the helper counts the
    repeat, and the decoy row shows its own EPlace next to the real copy's Eeff and label."""
    old = _old_seed_id("GOI_copy_1", "GenA_fna", 0, 1)
    decoy, real = _judged(old, "GOI_X~rev", 5.0), _judged(old, "GOI_X", 1e-8)
    repeats, writers, (d, r) = _genome_end([decoy, real], {"block_esyn": {0: 1e-4}})
    assert repeats == 1 and writers == {old: ["GOI_X~rev", "GOI_X"]}
    assert (d["EPlace"], d["Eeff"], d["Confidence"]) == ("5", "1e-08", "HIGH")


def test_the_genome_line_and_one_line_per_repeated_id(caplog):
    args = types.SimpleNamespace(genome_ctx={"judge_repeats": 1, "judge_writers": {
        "GOI_copy_1|G_fna_b0_l1": ["GOI_X~rev", "GOI_X"], "GOI_X|G_fna_b0_l1_exon_ann": ["GOI_X"]}})
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.log_judge_number_repeats(args, "G.fna")
        isr.log_judge_number_repeats(types.SimpleNamespace(), "H.fna")    # ended before its blocks
    assert [(r.levelname, r.getMessage()) for r in caplog.records] == [
        ("WARNING", "[G.fna] judge_numbers: 1 repeated model ids"),
        ("WARNING", "[G.fna] judge_numbers: GOI_copy_1|G_fna_b0_l1 written by 2 candidates, "
                    "SynVoy_Parent GOI_X~rev, GOI_X"),
        ("INFO", "[H.fna] judge_numbers: 0 repeated model ids")]


def test_the_search_writes_the_table_through_the_helper_and_logs_every_genome():
    src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    window = src[src.index("    def _judge_window(cands):"):]
    window = window[:window.index("        return greedy_select(cands)")]
    assert "record_judge_numbers(" in window and "] = c['judge_numbers']" not in window
    assert "    finally:\n        log_judge_number_repeats(args, genome_name)\n" in src
