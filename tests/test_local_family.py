"""A1 local_family (preset plan of 2026-10-15, PR-6): GOI loci per query in a window, and
the copies of the family in the target's slot.

* goi_candidate_regions_per_query replaces the hard-coded 2 GOI loci kept per query in a
  window; A1 sets it to the home copies + 1 (2 to 10), legacy and the defaults keep 2.
* report_copies_in_slot (A1 turns it on): CopiesInSlot / SlotRank on the GOI models and,
  per genome, home copies : copies in the slot in the report. Decides nothing.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import generate_report  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
import select_modules as sm  # noqa: E402

DEFAULTS = {"goi_relative_evalue": 1e-5, "goi_candidate_regions_per_query": 2,
            "report_copies_in_slot": False}


def _rel(n):
    return [{"gene_id": f"gene-{i}", "rank": i + 1, "evalue_seg_no": 1e-5} for i in range(n)]


# ============================================================ how many
def test_home_copies_from_the_coverage_or_the_family_members_nearby():
    assert sm.home_copies({"primary_copies": 2.34, "relatives_nearby": []}) == 2
    assert sm.home_copies({"primary_copies": 2.5, "relatives_nearby": None}) == 3   # half up
    assert sm.home_copies({"primary_copies": 1.0, "relatives_nearby": _rel(3)}) == 4
    assert sm.home_copies({"primary_copies": None, "relatives_nearby": _rel(1)}) == 2
    assert sm.home_copies({"primary_copies": None, "relatives_nearby": None}) is None


def test_a1_keeps_home_copies_plus_one_loci_between_2_and_10():
    by_id = {m.id: m for m in sm.MODULES}
    for copies, rel, want in ((1.6, 0, 3), (2.34, 0, 3), (1.0, 5, 7), (40.0, 0, 10), (None, None, 2)):
        f = {"primary_copies": copies, "relatives_nearby": _rel(rel) if rel is not None else None}
        settings, source = sm.merge([by_id["A1"]], f)
        assert settings["goi_candidate_regions_per_query"] == want, (copies, rel)
        assert settings["report_copies_in_slot"] is True
        assert source["goi_candidate_regions_per_query"] == "module:A1"


def test_the_report_switch_is_one_rule_of_a1():
    """One rule per key (profile handoff, stage 1): no other module can switch it off."""
    rule = sm.RULES_BY_KEY["report_copies_in_slot"]
    assert rule.owners == ("A1",) and rule.evidence == "report only"
    assert sm.merge([sm.module("B1"), sm.module("A1")], {})[0]["report_copies_in_slot"] is True
    assert "report_copies_in_slot" not in sm.merge([sm.module("B1")], {})[0]


def test_legacy_and_the_defaults_keep_two_loci_and_no_report():
    f = {"primary_copies": 3.0, "relatives_nearby": _rel(2)}
    legacy = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_tandem_family"},
                         modules_spec="legacy", features=f)
    assert legacy["settings"]["goi_candidate_regions_per_query"] == 2
    assert legacy["settings"]["report_copies_in_slot"] is False
    assert legacy["settings"]["goi_relative_evalue"] == 1e-3            # the old preset
    auto = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="auto", features=f)
    assert auto["settings"]["goi_candidate_regions_per_query"] == 4
    assert auto["settings"]["report_copies_in_slot"] is True
    assert auto["home_copies"] == legacy["home_copies"] == 3


def test_a_value_the_user_set_wins_over_a1():
    f = {"primary_copies": 3.0, "relatives_nearby": None}
    res = rep.resolve(dict(DEFAULTS, goi_candidate_regions_per_query=2), auto_preset=None,
                      modules_spec="auto", features=f, user_set={"goi_candidate_regions_per_query": 2})
    assert res["settings"]["goi_candidate_regions_per_query"] == 2
    assert res["user_set_kept"] == ["goi_candidate_regions_per_query"]


def test_the_runner_reads_the_key_where_the_fixed_two_was():
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    # The cap and its count live in helpers since stage 2 (cap_goi_parent_loci,
    # cap_flanking_loci); process_region_block calls them.
    i, j = src.index("def cap_goi_parent_loci"), src.index("def cap_flanking_loci")
    assert "sorted(parent_loci, key=rank, reverse=True)" in src[i:j]
    assert "getattr(args, 'goi_candidate_regions_per_query', 2)" in src[i:j]
    assert "[:2]" not in src[i:j]
    # the cap on the loci of a flanking gene stays
    assert re.search(r"\)\[:2\]", src[j:src.index("def count_fallback_refusal")])
    block = src[src.index("def process_region_block"):]
    assert "cap_goi_parent_loci(args, parent_id, parent_loci, _locus_rank)" in block
    assert "cap_flanking_loci(args, parent_loci)" in block


# ============================================================ the slot
def test_orientation_of_the_slot():
    assert isr.slot_inverted({-1: (100, 200), 1: (500, 600)}) is False
    assert isr.slot_inverted({-2: (900, 1000), 1: (100, 200)}) is True
    assert isr.slot_inverted({-1: (100, 600), 1: (500, 700)}) is None
    assert isr.slot_inverted({1: (100, 200)}) is None


def _row(model, start, end, conf, place, kind="mRNA", genome="GenA_fna", block=0):
    return (f"c1\texon_annotation\t{kind}\t{start}\t{end}\t.\t+\t.\t"
            f"ID=GOI_X|{genome}_b{block}_{model};SynVoy_Parent=GOI_X;SynVoyRole=goi;"
            f"Confidence={conf};Place={place}")


def _attrs(line):
    return dict(kv.split("=", 1) for kv in line.split("\t")[8].split(";"))


def test_copies_in_the_slot_are_counted_and_ranked_from_the_left_neighbour():
    lines = [
        _row("l1_exon_ann", 5000, 6000, "HIGH", "slot"),
        "c1\texon_annotation\tCDS\t5000\t6000\t.\t+\t0\tID=x_CDS1;Parent=GOI_X|GenA_fna_b0_l1_exon_ann",
        _row("l2_exon_ann", 3000, 4000, "AMBIGUOUS", "slot"),
        _row("l3_exon_ann", 7000, 8000, "LOW", "slot"),           # not a copy
        _row("l4_exon_ann", 9000, 9900, "MEDIUM", "core"),        # outside the slot
        _row("l1_exon_ann", 100, 900, "HIGH", "slot", block=1),   # a block without a slot
        _row("l1_exon_ann", 100, 900, "HIGH", "slot", genome="GenB_fna"),
    ]
    slots = {0: ((2000, 8500), False), 1: (None, None)}
    out = isr.annotate_slot_copies(lines, "GenA_fna", slots)
    a = [_attrs(x) if x.split("\t")[2] == "mRNA" else None for x in out]
    assert (a[0]["CopiesInSlot"], a[0]["SlotRank"]) == ("2", "2")
    assert (a[2]["CopiesInSlot"], a[2]["SlotRank"]) == ("2", "1")
    assert a[3]["CopiesInSlot"] == "2" and "SlotRank" not in a[3]
    assert a[4]["CopiesInSlot"] == "2" and "SlotRank" not in a[4]
    assert "CopiesInSlot" not in a[5] and "CopiesInSlot" not in a[6]
    assert out[1] == lines[1]
    inverted = isr.annotate_slot_copies(lines, "GenA_fna", {0: ((2000, 8500), True)})
    assert _attrs(inverted[0])["SlotRank"] == "1" and _attrs(inverted[2])["SlotRank"] == "2"
    assert isr.annotate_slot_copies(lines, "GenA_fna", {}) == lines


def test_the_report_pairs_home_copies_with_the_slot_of_the_best_call():
    recs = [
        {"genome": "GenA", "chrom": "c1", "start": 5000, "end": 6000, "confidence": "HIGH",
         "eeff": "1e-30", "block_esyn": "1e-9", "bits": "300", "copies_in_slot": "3", "slot_rank": "2"},
        {"genome": "GenA", "chrom": "c9", "start": 10, "end": 900, "confidence": "MEDIUM",
         "eeff": "1e-6", "block_esyn": "1e-3", "bits": "60", "copies_in_slot": "1", "slot_rank": "1"},
        {"genome": "GenB", "chrom": "c1", "start": 1, "end": 2, "confidence": "HIGH", "eeff": "1e-9"},
    ]
    out = generate_report.copies_in_slot_by_genome(recs, 2)
    assert set(out) == {"GenA"}
    assert out["GenA"]["unit"] == "2:3" and out["GenA"]["slot_of"]["chrom"] == "c1"
    assert generate_report.copies_in_slot_by_genome(recs, None)["GenA"]["unit"] == "?:3"
    assert generate_report.copies_in_slot_by_genome(recs[2:], 2) == {}
