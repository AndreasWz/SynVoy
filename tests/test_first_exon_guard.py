"""Do not invent a start where a short home first exon is missing (batch of 2026-10-15, G3A;
--goi_first_exon_guard, off by default). When a model's first CDS begins at query position p > 3
and p - 1 is within 5 aa of the home first exon's length (an exon of <= 30 aa), the start walk
would find an ATG in the intron; the model keeps miniprot's start and says FirstExon=missing,
MissingQueryStart=1-<p-1>; the same at the stop end. Home GOI and decoy parents only, never the
relative's query; single-exon homes unchanged. All synthetic."""
import json
import os
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import annotate_goi_exons as ag  # noqa: E402
import iterative_search_runner as isr  # noqa: E402

LY6E_LIKE = {"first_aa": 17, "last_aa": 73, "num_exons": 3}


def _cds(qs, qe, gs, ge, phase=0):
    return {"qstart": qs, "qend": qe, "gstart": gs, "gend": ge, "phase": phase, "strand": "+"}


def test_a_model_without_the_short_first_exon_is_flagged():
    cds = [_cds(19, 57, 3000, 3117), _cds(58, 130, 4000, 4219)]
    assert ag.missing_terminal_exons(cds, 130, LY6E_LIKE) == {"start": (1, 18)}
    assert ag.missing_terminal_exons(cds, 130, dict(LY6E_LIKE, first_aa=40)) == {}   # walks as today
    assert ag.missing_terminal_exons(cds, 130, dict(LY6E_LIKE, num_exons=1)) == {}   # one exon
    assert ag.missing_terminal_exons(cds, 130, None) == {}
    assert ag.missing_terminal_exons([_cds(30, 57, 3000, 3084)], 130, LY6E_LIKE) == {}   # 29 vs 17


def test_the_window_for_tiny_exons_and_the_stop_end():
    tiny = {"first_aa": 6, "last_aa": 12, "num_exons": 4}
    assert ag.missing_terminal_exons([_cds(5, 90, 1000, 1258)], 100, tiny)["start"] == (1, 4)
    assert "start" not in ag.missing_terminal_exons([_cds(4, 90, 1000, 1261)], 100, tiny)   # 3 < 4
    end = ag.missing_terminal_exons([_cds(1, 88, 1000, 1264)], 100, tiny)
    assert end == {"end": (89, 100)}                                    # 12 unaligned vs 12


def test_the_start_walk_is_skipped_only_when_asked():
    seq = "C" * 30 + "ATG" + "GCT" * 10 + "GCT" * 40 + "TAA" + "C" * 30
    edge = 30 + 3 + 30
    for walk, expect_start in ((True, 30), (False, edge)):
        cds = [{"gstart": edge, "gend": edge + 120, "phase": 0}]
        ag.refine_model_termini(cds, seq, "+", walk_start=walk, walk_stop=False)
        assert cds[0]["gstart"] == expect_start


def test_build_model_flags_and_keeps_miniprot_s_start():
    seq = "C" * 30 + "ATG" + "GCT" * 10 + "GCT" * 40 + "TAA" + "C" * 30
    edge = 63
    best = {"strand": "+", "identity": 0.8, "rank": 1, "score": 100, "protein": "A" * 40,
            "cds_list": [{"gstart": edge, "gend": edge + 120, "strand": "+", "phase": 0,
                          "qstart": 19, "qend": 58}]}
    exons, prot = ag._build_model(json.loads(json.dumps(best)), "M" * 58, seq, "c",
                                  home_structure=LY6E_LIKE)
    assert exons[0]["gstart"] == edge and exons[0]["missing_query_start"] == (1, 18)
    assert isr.missing_exon_attrs(exons) == ";MissingQueryStart=1-18;FirstExon=missing"
    exons2, _p = ag._build_model(json.loads(json.dumps(best)), "M" * 58, seq, "c")
    assert exons2[0]["gstart"] == 30 and "missing_query_start" not in exons2[0]   # walked, as today
    assert isr.missing_exon_attrs(exons2) == ""


def test_which_parents_get_the_guard(tmp_path):
    info = tmp_path / "goi_info.json"
    info.write_text(json.dumps({"method": "gff_annotation", "exons": [
        {"id": "GOI_q|exon_2", "exon_num": 2, "length_aa": 39},
        {"id": "GOI_q|exon_1", "exon_num": 1, "length_aa": 17},
        {"id": "GOI_q|exon_3", "exon_num": 3, "length_aa": 73},
        {"id": "GOI_copy_1", "exon_num": 4, "length_aa": 90}]}))
    assert isr.read_home_exons(str(info)) == ("gff_annotation", [17, 39, 73])
    assert isr.read_home_exons(str(tmp_path / "missing.json")) == (None, [])
    assert isr.home_exon_terms("gff_annotation", [17, 39, 73]) == (3, LY6E_LIKE)
    assert isr.home_exon_terms("hit_annotation", [17, 39, 73]) == (3, None)       # not from the GFF
    assert isr.home_exon_terms("gff_annotation", [40]) == (1, None)               # one exon
    args = types.SimpleNamespace(goi_first_exon_guard=True, home_goi_ids={"GOI_q", "GOI_q|exon_1"},
                                 home_goi_structure=LY6E_LIKE)
    assert isr.first_exon_guard_structure(args, "GOI_q") == LY6E_LIKE
    assert isr.first_exon_guard_structure(args, "GOI_q~rev") == LY6E_LIKE         # a decoy
    assert isr.first_exon_guard_structure(args, "GOI_q~rel") is None               # the relative
    assert isr.first_exon_guard_structure(args, "GOI_copy_1|g_b0_l1") is None      # a seeded copy
    args.goi_first_exon_guard = False
    assert isr.first_exon_guard_structure(args, "GOI_q") is None
