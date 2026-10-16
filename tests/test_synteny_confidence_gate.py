"""Where synteny ORDER enters a GOI call's label since 2026-10-15.

Until then a separate gate in the identity-based classifier demoted a HIGH call to MEDIUM
when its block's flanking genes were not in home order (classify_high_min_collinear, the
biglycan-as-decorin signature). The label is synvoy_label.judge's now, and order enters it
through the block's E_syn: the longest run of anchors in home order decides whether a block
is significant (E_syn <= 0.1), and only a significant block gives a call its slot or core
as the place of its E-value. A scrambled paralog neighbourhood is therefore not a place
that predicts the gene. Rewritten from the collinear-gate tests on 2026-10-15.
"""
import importlib.util
import os
import sys

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")
sys.path.insert(0, BIN)
spec = importlib.util.spec_from_file_location(
    "iterative_search_runner", os.path.join(BIN, "iterative_search_runner.py"))
isr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(isr)

import synvoy_evalue  # noqa: E402


def test_longest_collinear_run_helper():
    assert isr._longest_collinear_run([1, 2, 3, 4])[0] == 4        # ordered
    assert isr._longest_collinear_run([4, 3, 2, 1])[0] == 4        # reverse-ordered (still collinear)
    assert isr._longest_collinear_run([1, 9, 2, 8, 3])[0] == 3     # scrambled → run of 3 (1,2,3)
    assert isr._longest_collinear_run([5, 1, 4, 2, 3])[0] <= 3     # scrambled


def test_order_decides_whether_a_block_is_significant():
    """Twenty single-place anchors, a 300-kb window in a 12-Mb genome: six of them in home
    order are a significant block, the same six scrambled (longest run 2) are not."""
    places = [1] * 20
    ordered = synvoy_evalue.e_syn(isr._longest_collinear_run([1, 2, 3, 4, 5, 6])[0],
                                  places, 300_000, 12_000_000)
    scrambled = synvoy_evalue.e_syn(isr._longest_collinear_run([3, 1, 6, 2, 5, 4])[0],
                                    places, 300_000, 12_000_000)
    assert ordered <= synvoy_evalue.SIGNIFICANT_BLOCK_ESYN < scrambled


def test_every_goi_row_waits_for_the_judge():
    """The label written when a model is built is provisional; relabel_goi_lines replaces
    it with the judge's. A row the judge cannot label stays LOW (never counted)."""
    for evidence in ("exon_annotation", "fallback_hit_span", "rescued_exon", "raw_hit"):
        assert isr._provisional_goi_label(evidence) == (
            "LOW", "ambiguous_goi_family_member", "awaiting_label")
    assert isr._provisional_goi_label("tandem_copy") == ("LOW", "tandem_goi_copy", "awaiting_label")
    attrs = isr._goi_feature_attrs({"ID": "m1"}, evidence_type="exon_annotation", identity=95.0,
                                   exon_count=4, query_cov=0.98, flanking_support=9,
                                   collinear_support=9)
    assert (attrs["Confidence"], attrs["InferenceReason"]) == ("LOW", "awaiting_label")
    assert (attrs["ModelStatus"], attrs["BlockCollinearSupport"]) == ("complete", "9")


def test_a_row_without_eeff_is_counted_as_unjudged():
    row = ("chr1\tSynVoy\tmRNA\t1\t900\t.\t+\t.\t"
           "ID=GOI_x|G_b0_l1;SynVoyRole=goi;Confidence=LOW;GOIClass=ambiguous_goi_family_member;"
           "InferenceReason=awaiting_label")
    out, counts = isr.relabel_goi_lines([row])
    assert out == [row] and counts == {"LOW/unjudged": 1}
