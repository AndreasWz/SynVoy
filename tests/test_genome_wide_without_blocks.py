"""PR-9 (H1, --goi_genome_wide_without_blocks, off by default): the genome-wide path of FINAL
step 5 also in genomes without a block at E_syn <= 0.1, and RetrocopyCandidate=1 on a call that
runs across a home exon junction without an intron. All synthetic."""
import logging
import os
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402


def test_home_exon_junctions_are_the_last_residues_of_all_exons_but_the_last():
    assert isr.home_exon_junctions([17, 39, 73]) == [17, 56]
    assert isr.home_exon_junctions([120]) == [] and isr.home_exon_junctions([]) == []


def test_a_piece_spans_a_junction_with_ten_residues_on_each_side():
    junctions = [17, 56]
    assert isr.spans_a_junction([(5, 30)], junctions)                   # 8..17 and 18..27 aligned
    assert not isr.spans_a_junction([(9, 30)], junctions)               # only 9 residues before
    assert not isr.spans_a_junction([(5, 26)], junctions)               # only 9 residues after
    assert not isr.spans_a_junction([(1, 17), (18, 56), (57, 129)], junctions)   # one per exon
    assert isr.spans_a_junction([(1, 17), (40, 80)], junctions)         # exons 2-3 contiguous
    assert not isr.spans_a_junction([(5, 30)], [])                      # single-exon home gene


def _args(**kw):
    a = types.SimpleNamespace(genome_ctx={"block_esyn": {0: 0.5}}, home_goi_ids={"GOI_q"})
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def test_without_a_significant_block_the_pass_runs_only_with_the_switch(monkeypatch, caplog):
    monkeypatch.setattr(isr, "_home_proteome", lambda a: (None, None))
    with caplog.at_level(logging.INFO):
        assert isr.genome_wide_goi_pass(_args(), [], {}, {}, "g", "u", 1, None, "", 5000) == ([], [])
    assert "not run (--goi_genome_wide_without_blocks runs it)" in caplog.text
    caplog.clear()
    with caplog.at_level(logging.INFO):
        isr.genome_wide_goi_pass(_args(goi_genome_wide_without_blocks=True), [], {}, {}, "g", "u", 1,
                                 None, "", 5000)
    assert "no home proteome" in caplog.text                             # past the block gate
