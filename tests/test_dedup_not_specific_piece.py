"""Fix-phase item 2: dedup lets a counted piece override its not-specific whole model (the luciferase runs).

In dedupe_goi_annotations a counted piece contained in a longer model of the same cluster won on confidence, so a
paralog's exon at MEDIUM or HIGH beat its whole gene, which is AMBIGUOUS with Basis=not_specific: Round 2's luciferase
runs counted 6 calls, 3 of them such pieces. The whole gene's verdict now stands, unless the fused-model rule prefers
the piece: a wide model that spans at least twice the piece while adding under 15 points of query coverage. A piece is a
model the pipeline calls a fragment, on the wide model's strand, shorter, covering no more of the query: span containment
alone also took complete genes inside an intron and a model's rescue twin at the same coordinates.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import generate_report as gr  # noqa: E402


def _rec(start, end, conf, basis="significant_specific", cov="0.700", mid=None, goi_class="probable_goi",
         status="partial", strand="+", evidence="exon_annotation", cds=None):
    return {"genome": "G", "chrom": "c1", "start": start, "end": end, "strand": strand, "evidence": evidence,
            "cds": cds or [], "confidence": conf,
            "basis": basis, "goi_class": goi_class, "identity": "60.0", "query_cov": cov, "bits": "100",
            "target_gene": "", "model_status": status, "mrna_id": mid or f"GOI_x|m{start}", "source": "locus_1",
            "eeff": "1e-10", "place": "core", "block_esyn": "1e-5"}


def _piece(start, end, conf="MEDIUM", cov="0.300", **kw):
    return _rec(start, end, conf, cov=cov, status="fragment", **kw)


def _one(*recs):
    out = gr.dedupe_goi_annotations(list(recs))["records"]
    assert len(out) == 1
    return out[0]


# The real shapes of the luciferase runs, spans and coverages from Round 2's records; all three pieces are fragments.
LUCIFERASE = [
    (_rec(51107366, 51125346, "AMBIGUOUS", "not_specific", "0.978"), _piece(51116785, 51118110, "HIGH", cov="0.709")),
    (_rec(25214948, 25215970, "AMBIGUOUS", "not_specific", "0.627"), _piece(25215740, 25215970, "MEDIUM", cov="0.131")),
    (_rec(120728192, 120729890, "AMBIGUOUS", "not_specific", "0.882"), _piece(120729141, 120729875, "MEDIUM", cov="0.382")),
]


@pytest.mark.parametrize("wide,piece", LUCIFERASE)
def test_the_luciferase_pieces_leave_the_count(wide, piece):
    r = _one(dict(wide), dict(piece))
    assert (r["start"], r["end"], r["confidence"], r["basis"]) == (wide["start"], wide["end"], "AMBIGUOUS", "not_specific")
    assert r["overruled_counted_piece"] is True


@pytest.mark.parametrize("wide,piece", LUCIFERASE)
def test_the_order_of_the_records_does_not_matter(wide, piece):
    """Same start: the piece may come first (the cluster's order is start, then end)."""
    w, p = dict(wide), dict(piece)
    p["start"] = w["start"]
    p["end"] = w["start"] + (piece["end"] - piece["start"])
    a, b = _one(dict(w), dict(p)), _one(dict(p), dict(w))
    assert (a["start"], a["end"], a["confidence"]) == (b["start"], b["end"], b["confidence"]) == (w["start"], w["end"], "AMBIGUOUS")


def test_two_pieces_in_one_gene_leave_together():
    r = _one(_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.9"), _piece(2000, 2500, "HIGH", cov="0.4"),
             _piece(6000, 6500, "MEDIUM", cov="0.3"))
    assert (r["start"], r["confidence"]) == (1000, "AMBIGUOUS")


def _mirrored(recs, size=1_000_000):
    return [dict(r, start=size - r["end"], end=size - r["start"]) for r in recs]


def _counted(recs):
    return sorted(r["confidence"] for r in gr.dedupe_goi_annotations([dict(x) for x in recs])["records"]
                  if r["confidence"] in ("HIGH", "MEDIUM"))


def test_a_complete_gene_in_an_intron_stays_a_call():
    """Real from a search run: a GOI ortholog on the other strand inside the first intron of a paralog's ortholog."""
    wide = _rec(1000, 19000, "AMBIGUOUS", "not_specific", "0.95")
    inner = _rec(9000, 10500, "HIGH", cov="0.75", status="complete", strand="-")
    assert _counted([wide, inner]) == ["HIGH"]
    assert _counted([wide, dict(inner, strand="+")]) == ["HIGH"]


def test_a_model_and_its_rescue_twin_stay_a_call():
    """The plant run: Brassica FT, a main-path MEDIUM that holds by the positional tie-break, and the hull rescue's row
    at the same coordinates, which reads AMBIGUOUS not_specific because the rescue judges in the genome place."""
    main = _rec(22975663, 22977990, "MEDIUM", cov="1.000", status="complete")
    twin = _rec(22975663, 22977990, "AMBIGUOUS", "not_specific", "1.000", status="rescue")
    assert _counted([main, twin]) == ["MEDIUM"] and _counted([twin, main]) == ["MEDIUM"]
    # The same for a fragment-status twin at the same coordinates: no longer, so no piece.
    frag = dict(main, model_status="fragment")
    assert _counted([frag, twin]) == ["MEDIUM"] and _counted([twin, frag]) == ["MEDIUM"]


def test_a_near_identical_counted_model_stays_a_call():
    assert _counted([_rec(1000, 5000, "HIGH", cov="0.9", status="complete"),
                     _rec(1000, 5003, "AMBIGUOUS", "not_specific", "0.9")]) == ["HIGH"]


def test_a_wide_model_that_covers_less_than_its_piece_is_no_whole_gene():
    assert _counted([_rec(1000, 2000, "AMBIGUOUS", "not_specific", "0.50"), _piece(1100, 1800, "MEDIUM", cov="0.80")]) == ["MEDIUM"]
    assert _counted([_rec(1000, 2000, "AMBIGUOUS", "not_specific", "0.80"), _piece(1100, 1800, "MEDIUM", cov="0.50")]) == []


def test_three_copies_under_a_chain_stay_calls():
    assert _counted([_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.95"),
                     _rec(1000, 2000, "HIGH", cov="0.75"), _rec(4500, 5500, "HIGH", cov="0.75"),
                     _rec(8000, 9000, "HIGH", cov="0.75")]) == ["HIGH", "HIGH", "HIGH"]


@pytest.mark.parametrize("evidence", ["fallback_hit_span", "tandem_copy", "synteny_hull_rescue"])
def test_a_chain_a_tandem_row_or_a_hull_row_overrules_nothing(evidence):
    """The real shape: a fallback chain over four curated copies, a rescued exon on one of them. Not_specific says
    something about the chain, not about the copy. A hull row's verdict can come from its place alone."""
    wide = _rec(1000, 8400, "AMBIGUOUS", "not_specific", "0.737", evidence=evidence)
    assert _counted([wide, _piece(4900, 5000, "MEDIUM", cov="0.463", evidence="rescued_exon")]) == ["MEDIUM"]


def test_a_rescued_hit_can_overrule():
    """The luciferase record: a wide rescued exon (cov 0.627) over a MEDIUM one (cov 0.131)."""
    wide = _rec(25214948, 25215970, "AMBIGUOUS", "not_specific", "0.627", status="fragment", evidence="rescued_exon")
    piece = _piece(25215740, 25215970, "MEDIUM", cov="0.131", evidence="rescued_exon")
    assert _counted([wide, piece]) == []


def test_a_piece_in_the_wide_models_intron_stays_a_call():
    wide = _rec(1000, 19000, "AMBIGUOUS", "not_specific", "0.95", cds=[(1000, 1300), (18000, 19000)])
    in_intron = _piece(9000, 9300, "HIGH", cov="0.3", cds=[(9000, 9300)])
    on_exon = _piece(18200, 18500, "HIGH", cov="0.3", cds=[(18200, 18500)])
    assert _counted([wide, in_intron]) == ["HIGH"]
    assert _counted([wide, on_exon]) == []


def test_a_fragment_on_the_other_strand_is_no_piece():
    wide = _rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.95")
    assert _counted([wide, _piece(2000, 2500, "HIGH", cov="0.3", strand="-")]) == ["HIGH"]
    assert _counted([wide, _piece(2000, 2500, "HIGH", cov="0.3", strand="+")]) == []
    assert _counted([wide, _piece(2000, 2500, "HIGH", cov="0.3", strand=".")]) == []     # unknown: no test


def test_a_near_identical_fragment_is_no_piece():
    """A few base pairs longer does not make a model the whole gene of its own twin: the judge's order decides."""
    assert _counted([_piece(1000, 2000, "MEDIUM", cov="0.3"),
                     _rec(1000, 2003, "AMBIGUOUS", "not_specific", "0.3", status="fragment")]) == ["MEDIUM"]


def test_a_counted_model_merged_under_the_piece_is_judged_on_its_own():
    """Same start, so the short ones come first: a fragment, a complete model near-identical to it, then the wide one.
    The complete model must not vanish into the wide model's record."""
    frag = _piece(1000, 1500, "MEDIUM", cov="0.3")
    whole = _rec(1000, 1560, "MEDIUM", cov="0.35", status="complete")
    wide = _rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.95")
    assert _counted([frag, whole, wide]) == ["MEDIUM"]
    assert _counted(_mirrored([frag, whole, wide])) == ["MEDIUM"]


def test_a_piece_goes_back_when_the_wide_model_loses_the_record():
    """Wide + a fragment + a later compact model the fused rule prefers: the wide model is then no gene, and the first
    piece must not vanish into a record 1.5 kb away. The same on the mirrored scaffold."""
    recs = [_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.80"), _piece(2000, 2500, "MEDIUM", cov="0.30"),
            _rec(4000, 4500, "MEDIUM", cov="0.85")]
    assert _counted(recs) == ["MEDIUM", "MEDIUM"]
    assert _counted(_mirrored(recs)) == ["MEDIUM", "MEDIUM"]
    out = gr.dedupe_goi_annotations([dict(x) for x in recs])["records"]
    assert sorted((r["start"], r["n_merged_hits"]) for r in out) == [(2000, 1), (4000, 2)]
    assert not any(r.get("overruled_counted_piece") for r in out)


def test_a_better_wide_model_that_does_not_hold_the_piece_gives_it_back():
    w1 = dict(_rec(1000, 5000, "AMBIGUOUS", "not_specific", "0.8"), eeff="1e-5")
    w2 = dict(_rec(1400, 6000, "AMBIGUOUS", "not_specific", "0.8"), eeff="1e-30")
    piece = _piece(1100, 1500, "HIGH", cov="0.3")
    assert _counted([w1, piece, w2]) == ["HIGH"]
    assert _counted(_mirrored([w1, piece, w2])) == ["HIGH"]


def test_the_flag_sits_on_the_record_that_overruled():
    out = gr.dedupe_goi_annotations([_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.9"), _piece(2000, 2500)])["records"]
    assert [r.get("overruled_counted_piece") for r in out] == [True]
    for field in ("strand", "evidence", "cds"):
        assert field not in out[0]     # the record's fields do not change


def test_the_fused_model_rule_prefers_the_piece():
    """The wide model spans 4 times the piece and adds 5 points: it chained in foreign sequence."""
    r = _one(_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.90"), _piece(2000, 4000, "HIGH", cov="0.85"))
    assert (r["start"], r["confidence"], r["superseded_wide_model"]) == (2000, "HIGH", True)
    assert "overruled_counted_piece" not in r


def test_the_rule_reads_only_not_specific():
    """A gene-order-only wide model says nothing about another home gene: the piece wins, as before."""
    r = _one(_rec(1000, 9000, "AMBIGUOUS", "gene_order_only", "0.90"), _piece(2000, 2500, "HIGH", cov="0.40"))
    assert (r["start"], r["confidence"]) == (2000, "HIGH") and "overruled_counted_piece" not in r


def test_a_counted_wide_model_keeps_its_old_order():
    r = _one(_rec(1000, 9000, "MEDIUM", cov="0.90"), _piece(2000, 2500, "HIGH", cov="0.40"))
    assert (r["start"], r["confidence"]) == (2000, "HIGH")


def test_an_ambiguous_piece_is_no_counted_piece():
    """Two AMBIGUOUS models: the judge's order decides, as before, and nothing is flagged overruled."""
    r = _one(_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.90"), _rec(2000, 2500, "AMBIGUOUS", "gene_order_only", "0.40"))
    assert "overruled_counted_piece" not in r


def test_a_piece_outside_the_gene_stays_a_call():
    out = gr.dedupe_goi_annotations([_rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.90"), _piece(20000, 21000, "HIGH", cov="0.3")])
    assert sorted((r["start"], r["confidence"]) for r in out["records"]) == [(1000, "AMBIGUOUS"), (20000, "HIGH")]


def test_a_partial_overlap_is_no_containment():
    out = gr.dedupe_goi_annotations([_rec(1000, 3000, "AMBIGUOUS", "not_specific", "0.90"), _rec(2800, 5000, "HIGH", cov="0.8")])
    assert len(out["records"]) == 2


def test_other_genomes_and_chromosomes_stay_apart():
    a = _rec(1000, 9000, "AMBIGUOUS", "not_specific", "0.90")
    b = _piece(2000, 2500, "HIGH", cov="0.4")
    b["chrom"] = "c2"
    assert len(gr.dedupe_goi_annotations([a, b])["records"]) == 2


def test_a_hull_twin_that_outranks_its_exon_twin_still_overrules_the_piece():
    """Fix item 44 (9 Nov 2026): c7 tested the piece against the cluster's representative only, and a hull rescue
    (synteny_hull_rescue) is no whole gene, so when the hull twin of a gene model outranked it, the counted piece
    under both stayed counted. The piece is now tested against the representative's whole-model twin."""
    w = _rec(1000, 5000, "AMBIGUOUS", basis="not_specific", cov="0.900", mid="GOI_x|w",
             goi_class="syntenic_candidate_unconfirmed", cds=[(1000, 5000)])
    w["eeff"] = "1e-5"
    hull = dict(_rec(1000, 5000, "AMBIGUOUS", basis="not_specific", cov="0.900", mid="GOI_x|h",
                     goi_class="synteny_hull_rescue", evidence="synteny_hull_rescue", cds=[(1000, 5000)]),
                eeff="1e-8")
    piece = _piece(2000, 2600, cds=[(2000, 2600)])
    for order in ((w, hull, piece), (piece, hull, w), (hull, w, piece)):
        recs = gr.dedupe_goi_annotations([dict(r) for r in order])["records"]
        assert len(recs) == 1, recs
        assert recs[0]["confidence"] == "AMBIGUOUS" and recs[0].get("overruled_counted_piece"), recs[0]
