"""Every search block of a genome has a number of its own (finding 2 of the plot-selection work,
2026-10-19). Extra blocks were numbered from the number of models found so far (proxy: blocks
+ 1000 + models, sweep: blocks + 2000 + models), so an extra block that yielded no model gave
the next extra block of the same kind its number: that block's E_syn, context and table row
replaced its own, and its E_syn left SigBlocks. next_block_number keeps the old number when it
is free and moves a taken one up to the next free number."""
import csv
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402

TAG = "G_fna"


def _args(tmp_path):
    ranks = {f"g{i}": i for i in range(1, 11)}
    return SimpleNamespace(output_dir=str(tmp_path),
                           genome_ctx={"genome_len": 20_000_000, "places": {g: 2 for g in ranks},
                                       "block_esyn": {}},
                           home_rank=ranks, goi_ranks={}, home_geometry={"pitch": 20_000},
                           expansion=1.0)


def test_a_free_number_stays_a_taken_one_moves_up():
    used = {0, 1, 2}
    assert isr.next_block_number(used, 1005) == 1005
    assert isr.next_block_number(used, 1005) == 1006
    assert isr.next_block_number(used, 1005) == 1007
    assert isr.next_block_number(used, 2) == 3
    assert used == {0, 1, 2, 3, 1005, 1006, 1007}


def _two_extra_blocks(tmp_path, offset, numbers):
    """Two extra blocks of one kind, each significant, searched while the model count stays 0
    (the first yields no model), numbered by ``numbers``; one GOI row in the second."""
    args = _args(tmp_path)
    n_blocks = 2
    wanted = n_blocks + offset + 0
    got = []
    for k, chrom in enumerate(("c7", "c8")):
        b = numbers(wanted)
        got.append(b)
        start = 1_000_000 * (k + 1)
        isr.record_block_row(args, b, {"chrom": chrom, "start": start, "end": start + 100_000,
                                       "genes": ["gene-GOI_PROXY"]}, start - 50_000, start + 150_000, 1, 7)
        isr.block_place_terms(args, b, start - 50_000, start + 150_000, start - 50_000,
                              start + 150_000, 7, [], hood_len=200_000)
    line = (f"c8\tx\tmRNA\t2000001\t2000900\t.\t+\t.\tID=GOI_q|{TAG}_b{got[1]}_l1_exon_ann;"
            f"SynVoy_Parent=GOI_q;SynVoyRole=goi;Place=core;EPlace=1e-08;EPlaceTblastn=inf")
    out = isr.finalize_place_evalues([line], args.genome_ctx)
    isr._write_block_table(args, "", "G.fna", out, args.genome_ctx)
    with open(tmp_path / "blocks" / "G.fna.blocks.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    sig = isr._parse_gff_attr_string(out[0].split("\t")[8])["SigBlocks"]
    return got, rows, sig


def test_two_proxies_keep_two_numbers_two_rows_and_both_esyn(tmp_path):
    used = set(range(2))
    got, rows, sig = _two_extra_blocks(tmp_path / "new", 1000, lambda w: isr.next_block_number(used, w))
    assert got == [1002, 1003]
    assert [r["block"] for r in rows] == ["1002", "1003"] and [r["chrom"] for r in rows] == ["c7", "c8"]
    assert sig == "2"
    old, rows, sig = _two_extra_blocks(tmp_path / "old", 1000, lambda w: w)   # numbering before the fix
    assert old == [1002, 1002] and len(rows) == 1 and rows[0]["chrom"] == "c8" and sig == "1"


def test_two_sweeps_keep_two_numbers(tmp_path):
    used = set(range(2))
    got, rows, sig = _two_extra_blocks(tmp_path, 2000, lambda w: isr.next_block_number(used, w))
    assert got == [2002, 2003] and len(rows) == 2 and sig == "2"


def test_genome_wide_windows_keep_their_gaps_and_move_up(monkeypatch):
    places = [{"parent": "GOI_q", "chrom": "c1", "start": s, "end": s + 900, "hits": 1,
               "hit_list": [{"pident": 50.0}]} for s in (1000, 5000, 9000)]
    windows = {1000: (500, 2500), 5000: (5, 5), 9000: (8500, 10500)}     # the 2nd has no window
    called = []
    monkeypatch.setattr(isr, "_home_proteome", lambda args: ("proteome", "hood"))
    monkeypatch.setattr(isr, "home_goi_parents", lambda args: ["GOI_q"])
    monkeypatch.setattr(isr, "select_goi_queries", lambda *a, **k: (None, {}, {}))
    monkeypatch.setattr(isr, "genome_wide_places", lambda *a, **k: [dict(p) for p in places])
    monkeypatch.setattr(isr, "_hit_pieces", lambda genome_seqs, h: [])
    monkeypatch.setattr(isr, "choose_genome_wide_places", lambda ps, cap_of: list(ps))
    monkeypatch.setattr(isr, "genome_wide_window", lambda p, margin, n, w: windows[p["start"]])
    monkeypatch.setattr(isr, "_write_genome_wide_table", lambda *a, **k: None)

    def fake_block(number, block, *a, **k):
        called.append(number)
        return [], []
    monkeypatch.setattr(isr, "process_region_block", fake_block)
    args = SimpleNamespace(genome_ctx={"block_esyn": {0: 1e-5}, "genome_len": 10_000_000, "windows": []},
                           gap_search_window=1000, goi_home_span=0, expansion=1.0, goi_fallback_max_gap=0,
                           max_intron=5000, home_goi_ids=None, home_goi_gene_id=None, n_flanking_genes=10)

    def run(**kw):
        called.clear()
        isr.genome_wide_goi_pass(args, [], {"c1": "A" * 20_000}, {}, "G.fna", "u", 1, None, "", 5003, **kw)
        return list(called)
    assert run() == [5003, 5005]                                     # numbering before the fix
    assert run(block_numbers={0, 1, 2}) == [5003, 5005]              # no collision: unchanged
    assert run(block_numbers={0, 1, 2, 5005}) == [5003, 5006]        # the taken 5005 moves up


def test_the_three_call_sites_use_the_helper():
    src = open(os.path.join(os.path.dirname(__file__), "..", "bin", "iterative_search_runner.py")).read()
    assert "_block_numbers = set(range(len(synteny_blocks)))" in src
    assert ("_proxy_idx = next_block_number(\n"
            "                            _block_numbers, len(synteny_blocks) + 1000 + len(all_genes))") in src
    assert ("_sw_idx = next_block_number(\n"
            "                        _block_numbers, len(synteny_blocks) + 2000 + len(all_genes))") in src
    assert "block_numbers=_block_numbers)" in src
    assert "next_block_number(block_numbers, first_idx + k)" in src
    assert "len(synteny_blocks) + 1000 + len(all_genes)\n" not in src.replace(
        "_block_numbers, len(synteny_blocks) + 1000 + len(all_genes))", "")
