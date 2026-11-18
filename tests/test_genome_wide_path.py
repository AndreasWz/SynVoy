"""The genome-wide path (FINAL step 5; --goi_genome_wide_path): places of GOI
hits outside every searched GOI zone and every written row, tried when the judge's HIGH terms
hold in the genome place, best first per home GOI; windows clipped at the zones and judged in the
genome place. The answer rules are in test_genome_wide_answer.py. All synthetic."""
import os
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402


def _hit(q, chrom, s, e, strand="+", pident=30.0):
    return {"query": q, "chrom": chrom, "start": s, "end": e, "strand": strand, "pident": pident,
            "qstart": 1, "qend": (e - s) // 3, "evalue": 1e-5, "bits": 40.0}


def test_places_are_the_goi_hits_outside_every_goi_zone_and_every_row():
    hits = [_hit("GOI_q", "c1", 1000, 1300), _hit("GOI_q|exon_2", "c1", 4000, 4300),   # one place
            _hit("GOI_q", "c1", 90000, 90300),                                         # a second
            _hit("GOI_q|Sp_fna_b0_l1_exon_ann", "c2", 500, 800),                        # a model hit
            _hit("GOI_q", "c1", 200000, 200300),                                       # in a GOI zone
            _hit("gene-FLANK", "c1", 5000, 5300),                                      # flanking
            _hit("GOI_q~rev", "c3", 100, 400),                                         # a decoy
            _hit("GOI_q", "c4", 7000, 7300)]                                           # under a row
    places = isr.genome_wide_places(hits, [("c1", 150000, 250000)], {"GOI_q"}, max_gap=10000,
                                    records=[("c4", 6900, 8000)])
    assert [(p["chrom"], p["start"], p["end"], p["hits"]) for p in places] == \
        [("c1", 1000, 4300, 2), ("c1", 90000, 90300, 1), ("c2", 500, 800, 1)]
    assert {p["parent"] for p in places} == {"GOI_q"}


def test_a_hit_between_a_goi_zone_and_the_window_edge_is_a_place():
    """Bug list item 4: the window is searched over block +- 7 gene steps, the GOI zone is block +-
    2; a hit between the two was skipped as "inside a window" and is a place now."""
    zone, window = ("c1", 20000, 30000), ("c1", 5000, 45000)
    hits = [_hit("GOI_q", "c1", 38000, 38600)]
    assert isr.genome_wide_places(hits, [window], {"GOI_q"}, max_gap=10000) == []
    places = isr.genome_wide_places(hits, [zone], {"GOI_q"}, max_gap=10000)
    assert [(p["start"], p["end"]) for p in places] == [(38000, 38600)]


def test_no_place_spans_a_goi_zone():
    """Two hits on either side of a narrow GOI zone, closer than max_gap, are two places: a place
    across the zone had a window that held the block's zone (searched again in the genome place)."""
    hits = [_hit("GOI_q", "c1", 99000, 99300), _hit("GOI_q", "c1", 146000, 146300)]
    places = isr.genome_wide_places(hits, [("c1", 100000, 145000)], {"GOI_q"}, max_gap=50000)
    assert [(p["start"], p["end"]) for p in places] == [(99000, 99300), (146000, 146300)]
    for p in places:
        lo, hi = isr.genome_wide_window(p, 48000, 1_000_000, [("c1", 100000, 145000)])
        assert not isr.overlaps_any("c1", lo, hi, [("c1", 100000, 145000)])
    joined = isr.genome_wide_places(hits, [], {"GOI_q"}, max_gap=50000)
    assert [(p["start"], p["end"], p["hits"]) for p in joined] == [(99000, 146300, 2)]


def test_hits_of_another_gene_need_their_own_home_parent():
    hits = [_hit("GOI_copy_1|Sp_fna_b0_l1", "c1", 1000, 1300), _hit("GOI_a", "c1", 5000, 5300)]
    assert [p["parent"] for p in isr.genome_wide_places(hits, [], {"GOI_q"}, 100)] == ["GOI_q", "GOI_q"]
    two = isr.genome_wide_places(hits, [], {"GOI_a", "GOI_b"}, 100)
    assert [(p["parent"], p["start"]) for p in two] == [("GOI_a", 5000)]       # copy: no parent


def test_the_window_stays_outside_every_goi_zone():
    place = {"chrom": "c1", "start": 10000, "end": 11000}
    assert isr.genome_wide_window(place, 5000, 100000, []) == (5000, 16000)
    assert isr.genome_wide_window(place, 5000, 12000, []) == (5000, 12000)               # chromosome end
    zones = [("c1", 0, 8000), ("c1", 13000, 20000), ("c2", 0, 100000)]
    assert isr.genome_wide_window(place, 5000, 100000, zones) == (8000, 13000)
    lo, hi = isr.genome_wide_window(place, 5000, 100000, [("c1", 10500, 30000)])   # a later zone on it
    assert hi <= lo


def test_genome_wide_windows_are_no_blocks():
    ctx = {"genome_len": 10_000_000, "places": {}, "block_esyn": {0: 1e-5}}
    args = types.SimpleNamespace(genome_ctx=ctx, home_rank={}, home_geometry={}, expansion=1.0)
    terms = isr.block_place_terms(args, 7, 100, 9000, 100, 9000, 3, [], genome_wide=True)
    assert terms["significant"] is False and terms["slot"] is None and ctx["block_esyn"] == {0: 1e-5}
    assert synvoy_evalue.place_of((200, 400), terms["significant"], terms["slot"], terms["core"],
                                  terms["window"], terms["genome_len"])[0] == "genome"


def test_hit_pieces_are_split_at_stops_on_both_strands():
    seq = "ATGGCT" + "TAA" + "GCTGCTGCTGCTGCT"
    assert isr._hit_pieces({"c": seq}, {"chrom": "c", "start": 0, "end": 24, "strand": "+"}) == ["AAAAA"]
    rc = isr.reverse_complement(seq)
    assert isr._hit_pieces({"c": rc}, {"chrom": "c", "start": 0, "end": 24, "strand": "-"}) == ["AAAAA"]


def test_the_pass_needs_the_home_proteome(monkeypatch, caplog):
    """A significant block is not needed any more (the genome place does not depend on the
    blocks); the home proteome is: without it there is no specificity and no answer."""
    import logging
    args = types.SimpleNamespace(genome_ctx={"block_esyn": {0: 0.5}}, home_goi_ids={"GOI_q"})
    monkeypatch.setattr(isr, "_home_proteome", lambda a: (None, None))
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        assert isr.genome_wide_goi_pass(args, [], {}, {}, "g", "u", 1, None, "", 5000) == ([], [])
    assert "no home proteome" in caplog.text and "not run" in caplog.text


def test_no_place_outside_the_goi_zones_is_said_and_tabled(monkeypatch, tmp_path, caplog):
    """Fix item 6: with every GOI hit inside a searched GOI zone the pass says so and writes the
    header-only table (a pass without a significant block runs as well)."""
    import logging
    args = types.SimpleNamespace(genome_ctx={"block_esyn": {0: 0.5}, "goi_cores": [("c1", 0, 1000)]},
                                 home_goi_ids={"GOI_q"}, gap_search_window=500, max_intron=1000,
                                 goi_home_span=0, expansion=1.0, output_dir=str(tmp_path))
    monkeypatch.setattr(isr, "_home_proteome", lambda a: (object(), None))
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        out = isr.genome_wide_goi_pass(args, [_hit("GOI_q", "c1", 100, 400)], {}, {}, "g", "u", 1,
                                       None, "", 5000)
    assert out == ([], [])
    assert "[g] genome-wide path: 0 place(s) outside every searched GOI zone" in caplog.text
    table = tmp_path / "hits" / "g.genome_wide.tsv"
    assert table.read_text().splitlines() == ["\t".join(isr.GENOME_WIDE_COLUMNS)]


def test_the_report_lists_the_places_and_the_decoy_places_apart(tmp_path):
    path = tmp_path / "locus_1_G.fna.genome_wide.tsv"
    rows = [dict(genome="G.fna", parent="GOI_q", chrom="c1", start=100, end=400, hits=2,
                 best_hit_identity=87.5, evalue_genome="2e-80", bits="260.1", coverage="0.910",
                 retrocopy_hits="0", best_other="gene-p", margin="47.0", real_competitor=1, specific=1,
                 rank=1, decision="modelled", block=5003, model_decision="answer", near_block="",
                 significant_blocks=0, median_flank_identity="71.2"),
            dict(genome="G.fna", parent="GOI_q", chrom="c2", start=100, end=400, hits=1,
                 best_hit_identity=40.0, evalue_genome="1e-05", coverage="0.400", decision="evalue",
                 significant_blocks=0, median_flank_identity="71.2"),
            dict(genome="G.fna", parent="GOI_q", decoy="rev", chrom="c3", start=10, end=90, hits=1,
                 best_hit_identity=30.0, evalue_genome="0.3", coverage="0.100", decision="evalue",
                 significant_blocks=0, median_flank_identity="71.2")]
    with open(path, "w") as fh:
        fh.write("\t".join(isr.GENOME_WIDE_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in isr.GENOME_WIDE_COLUMNS) + "\n")
    records = gr.collect_genome_wide_places([str(path)])
    block = gr.summarize_genome_wide_places(records)
    assert block["total"] == 2 and block["by_decision"] == {"modelled": 1, "evalue": 1}
    assert block["by_model_decision"] == {"answer": 1}
    assert [r["evalue_genome"] for r in block["places"]] == ["2e-80", "1e-05"]
    assert block["places"][0]["coverage"] == "0.910" and block["places"][0]["block"] == "5003"
    assert block["decoy_places"]["total"] == 1 and block["decoy_places"]["places"][0]["decoy"] == "rev"
    assert "decoy_mode" in block["note"]
