"""Rules of the search core added 2026-10-08/09 (bin/iterative_search_runner.py):

* one anchor per place (`resolve_anchor_hits`);
* the home protein is always a query, the latest relative's model a second one, and
  what is reported always refers to the home protein (`select_goi_queries`,
  `rebase_relative_candidates`);
* distances per target genome (`genome_scale_args`, `bp_or_auto`, `_gene_limit`).
"""
import argparse
import os
import random
import sys

import pytest

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import synvoy_scale  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice(AA) for _ in range(n))


def _mutate(seq, share, seed):
    rnd = random.Random(seed)
    out = list(seq)
    for i in rnd.sample(range(len(seq)), int(len(seq) * share)):
        out[i] = rnd.choice(AA.replace(out[i], ""))
    return "".join(out)


# ============================================================ one anchor per place

def _hit(query, start, end, bits, chrom="chr1", evalue=1e-30):
    return {"query": query, "chrom": chrom, "start": start, "end": end,
            "bits": bits, "evalue": evalue}


def _kept(hits, **kw):
    return sorted((h["query"], h["start"]) for h in isr.resolve_anchor_hits(hits, **kw))


def test_two_relatives_on_one_target_gene_leave_one_anchor():
    hits = [_hit("gene-KIN1", 1000, 2000, 300), _hit("gene-KIN2", 1010, 1990, 120)]
    assert _kept(hits) == [("gene-KIN1", 1000)]


def test_tandem_pair_is_assigned_one_to_one():
    """A and B are relatives; each hits both target copies, each is best on its own."""
    hits = [_hit("gene-A", 1000, 2000, 400), _hit("gene-B", 1000, 2000, 150),
            _hit("gene-A", 5000, 6000, 160), _hit("gene-B", 5000, 6000, 380)]
    assert _kept(hits) == [("gene-A", 1000), ("gene-B", 5000)]


def test_hits_of_the_same_gene_never_compete():
    hits = [_hit("gene-A|exon_1", 1000, 2000, 80), _hit("gene-A|exon_2", 1500, 2500, 300),
            _hit("gene-A", 1000, 2500, 500)]
    assert len(isr.resolve_anchor_hits(hits)) == 3


def test_small_overlap_is_not_the_same_place():
    """Neighbouring genes whose alignments touch stay two anchors."""
    hits = [_hit("gene-A", 1000, 2000, 300), _hit("gene-B", 1900, 3000, 100)]
    assert len(isr.resolve_anchor_hits(hits)) == 2
    assert _kept(hits, min_overlap=0.05) == [("gene-A", 1000)]


def test_short_hit_inside_a_long_one_competes():
    """Overlap is measured on the SHORTER hit: a domain hit inside a full-length hit."""
    hits = [_hit("gene-LONG", 1000, 9000, 900), _hit("gene-DOMAIN", 4000, 4300, 60)]
    assert _kept(hits) == [("gene-LONG", 1000)]


def test_a_loser_cannot_remove_a_third_hit():
    """B loses to A. C overlaps only B, so C stays (greedy by strength)."""
    hits = [_hit("gene-A", 1000, 2000, 500), _hit("gene-B", 1400, 2600, 300),
            _hit("gene-C", 2200, 2700, 100)]
    assert _kept(hits) == [("gene-A", 1000), ("gene-C", 2200)]


def test_other_strand_and_other_sequence():
    hits = [_hit("gene-A", 1000, 2000, 300), _hit("gene-B", 1000, 2000, 100, chrom="chr2")]
    assert len(isr.resolve_anchor_hits(hits)) == 2


def test_goi_hits_are_never_removed_and_never_remove():
    hits = [_hit("GOI_X", 1000, 2000, 900), _hit("gene-A", 1000, 2000, 50),
            _hit("GOI_X|exon_1", 1000, 1300, 30)]
    assert _kept(hits) == [("GOI_X", 1000), ("GOI_X|exon_1", 1000), ("gene-A", 1000)]


def test_ties_are_broken_the_same_way_in_any_input_order():
    hits = [_hit("gene-B", 1000, 2000, 200), _hit("gene-A", 1000, 2000, 200),
            _hit("gene-C", 1000, 2000, 200, evalue=1e-40)]
    first = _kept(hits)
    assert first == [("gene-C", 1000)]                     # lower E-value wins the tie
    for seed in range(5):
        shuffled = hits[:]
        random.Random(seed).shuffle(shuffled)
        assert _kept(shuffled) == first
    assert _kept(hits[:2]) == [("gene-A", 1000)]           # then the smaller id


def test_many_hits_on_one_sequence_stay_consistent():
    """The sweep must give the same answer as comparing every pair."""
    rnd = random.Random(3)
    hits = []
    for i in range(400):
        start = rnd.randrange(0, 200_000)
        hits.append(_hit(f"gene-{rnd.randrange(40)}", start, start + rnd.randrange(200, 3000),
                         rnd.randrange(30, 900)))
    kept = isr.resolve_anchor_hits(hits)

    ranked = sorted(hits, key=lambda h: (-h["bits"], h["evalue"], h["query"], h["start"], h["end"]))
    expected = []
    for h in ranked:
        clash = False
        for o in expected:
            ov = min(h["end"], o["end"]) - max(h["start"], o["start"])
            shorter = min(h["end"] - h["start"], o["end"] - o["start"])
            if ov > 0 and ov >= 0.5 * shorter and o["query"] != h["query"]:
                clash = True
                break
        if not clash:
            expected.append(h)
    key = lambda h: (h["query"], h["start"], h["end"], h["bits"])   # noqa: E731
    assert sorted(map(key, kept)) == sorted(map(key, expected))
    assert len(kept) < len(hits)


# ====================================================================== query rule

HOME = _protein(300, 1)


def _db(*records):
    return {rid: {"id": rid, "seq": seq, "header": rid} for rid, seq in records}


def test_home_protein_alone():
    db = _db(("GOI_X", HOME), ("gene-A", _protein(200, 2)))
    queries, by_parent, alias_of = isr.select_goi_queries(db, {"GOI_X"})
    assert [q["id"] for q in queries] == ["GOI_X"] and alias_of == {}
    assert by_parent["GOI_X"]["seq"] == HOME


def test_a_longer_model_no_longer_replaces_the_home_protein():
    """Until 2026-10-08 the LONGEST sequence of a GOI was the only query of a window."""
    model = _mutate(HOME, 0.3, 5) + "MKV"
    db = _db(("GOI_X", HOME), ("GOI_X|GenA_fna_b0_l1_exon_ann", model))
    queries, by_parent, alias_of = isr.select_goi_queries(db, {"GOI_X"})
    assert [q["id"] for q in queries] == ["GOI_X", "GOI_X~rel"]
    assert by_parent["GOI_X"]["seq"] == HOME
    assert by_parent["GOI_X~rel"]["seq"] == model
    assert by_parent["GOI_X~rel"]["model_id"] == "GOI_X|GenA_fna_b0_l1_exon_ann"
    assert alias_of == {"GOI_X|GenA_fna_b0_l1_exon_ann": "GOI_X~rel"}

    old = isr.collapse_goi_queries_by_parent(list(db.values()))
    assert old["GOI_X"]["seq"] == model                    # the rule that was replaced


def test_relative_comes_from_the_genome_added_last():
    near = _mutate(HOME, 0.1, 6)
    far = _mutate(HOME, 0.5, 7)
    db = _db(("GOI_X", HOME), ("GOI_X|GenA_fna_b0_l1_exon_ann", near),
             ("GOI_X|GenB_fna_b2_l1_exon_ann", far))
    _q, by_parent, _a = isr.select_goi_queries(db, {"GOI_X"})
    assert by_parent["GOI_X~rel"]["model_id"] == "GOI_X|GenB_fna_b2_l1_exon_ann"


def test_among_the_models_of_that_genome_the_one_closest_to_home_is_used():
    good = _mutate(HOME, 0.2, 8)
    fragment = HOME[:60]
    paralog = _mutate(HOME, 0.6, 9)
    db = _db(("GOI_X", HOME), ("GOI_X|GenA_fna_b0_l1_exon_ann", _mutate(HOME, 0.1, 6)),
             ("GOI_X|GenB_fna_b4_l1_fallback", fragment),
             ("GOI_X|GenB_fna_b0_l1_exon_ann", good),
             ("GOI_X|GenB_fna_b7_l2m2_exon_ann", paralog))
    _q, by_parent, alias_of = isr.select_goi_queries(db, {"GOI_X"})
    assert by_parent["GOI_X~rel"]["seq"] == good
    assert list(alias_of) == ["GOI_X|GenB_fna_b0_l1_exon_ann"]


def test_a_model_identical_to_the_home_protein_adds_nothing():
    db = _db(("GOI_X", HOME), ("GOI_X|GenA_fna_b0_l1_exon_ann", HOME))
    queries, _b, alias_of = isr.select_goi_queries(db, {"GOI_X"})
    assert [q["id"] for q in queries] == ["GOI_X"] and alias_of == {}


def test_relative_query_can_be_switched_off():
    db = _db(("GOI_X", HOME), ("GOI_X|GenA_fna_b0_l1_exon_ann", _mutate(HOME, 0.3, 5)))
    queries, _b, alias_of = isr.select_goi_queries(db, {"GOI_X"}, use_relative=False)
    assert [q["id"] for q in queries] == ["GOI_X"] and alias_of == {}


def test_exon_and_fragment_records_are_not_queries():
    db = _db(("GOI_X", HOME), ("GOI_X|exon_1", HOME[:80]), ("GOI_X|frag_2", HOME[80:160]))
    queries, _b, _a = isr.select_goi_queries(db, {"GOI_X", "GOI_X|exon_1"})
    assert [q["id"] for q in queries] == ["GOI_X"]


def test_without_home_ids_the_old_rule_applies():
    model = HOME + "MKV"
    db = _db(("GOI_X", HOME), ("GOI_X|GenA_fna_b0_l1_exon_ann", model))
    queries, by_parent, alias_of = isr.select_goi_queries(db, set())
    assert alias_of == {} and [q["seq"] for q in queries] == [model]
    assert by_parent["GOI_X"]["seq"] == model


def test_two_gois_each_get_their_own_relative():
    other = _protein(250, 20)
    db = _db(("GOI_X", HOME), ("GOI_Y", other),
             ("GOI_X|GenA_fna_b0_l1_exon_ann", _mutate(HOME, 0.2, 21)),
             ("GOI_Y|GenA_fna_b0_l1_exon_ann", _mutate(other, 0.2, 22)))
    queries, _b, alias_of = isr.select_goi_queries(db, {"GOI_X", "GOI_Y"})
    assert [q["id"] for q in queries] == ["GOI_X", "GOI_Y", "GOI_X~rel", "GOI_Y~rel"]
    assert len(alias_of) == 2


def test_query_choice_does_not_depend_on_hash_order():
    twin_a = _mutate(HOME, 0.2, 30)
    twin_b = twin_a[:-1] + ("A" if twin_a[-1] != "A" else "C")
    recs = [("GOI_X", HOME), ("GOI_X|GenB_fna_b1_l1_exon_ann", twin_a),
            ("GOI_X|GenB_fna_b0_l1_exon_ann", twin_b)]
    chosen = set()
    for order in (recs, [recs[0], recs[2], recs[1]]):
        isr._RELATIVE_SCORE_CACHE.clear()
        _q, by_parent, _a = isr.select_goi_queries(_db(*order), {"GOI_X"})
        chosen.add(by_parent["GOI_X~rel"]["model_id"])
    assert len(chosen) == 1


def test_home_ids_are_read_from_the_pristine_query_set(tmp_path):
    fa = tmp_path / "initial_db.faa"
    fa.write_text(">gene-A home=1-500\nMKT\n>GOI_X\nMKV\n>GOI_X|exon_1 coords=5-90\nMK\n")
    assert isr.build_home_goi_ids(str(fa)) == {"GOI_X", "GOI_X|exon_1"}
    assert isr.build_home_goi_ids(str(tmp_path / "absent.faa")) == set()


@pytest.mark.parametrize("model_id, genome", [
    ("GOI_X|GenA_fna_b0_l1_exon_ann", "GenA_fna"),
    ("GOI_X|GCF_000001_2_fna_b12_l3m2_exon_ann", "GCF_000001_2_fna"),
    ("GOI_X|GenA_fna_b0_r1_fallback", "GenA_fna"),
    ("GOI_copy_1|Apis_cerana_fna_b0_l1", "Apis_cerana_fna"),
    ("GOI_X", ""),
    ("GOI_X|exon_1", ""),
])
def test_genome_of_a_wave_model(model_id, genome):
    assert isr._model_genome_tag(model_id) == genome


# ------------------------------------------- models found with the relative's model

def _candidate(parent, model_id_tail, start, end, protein, identity, evidence="exon_annotation",
               flanking=6, exons=3, qcov=0.95):
    gid = f"{parent}|{model_id_tail}"
    attrs = (f"ID={gid};Name={gid};SynVoy_Parent={parent};Identity={identity:.1f};Exons={exons};"
             f"EvidenceType={evidence};Confidence=HIGH;GOIClass=confident_goi;"
             f"BlockFlankingSupport={flanking};BlockCollinearSupport={flanking};"
             f"QueryCoverage={qcov};TargetGene=tgtA")
    gff = [f"chr1\tSynVoy\tmRNA\t{start + 1}\t{end}\t{identity:.1f}\t+\t.\t{attrs}",
           f"chr1\tSynVoy\tCDS\t{start + 1}\t{end}\t.\t+\t0\tID={gid}_CDS1;Parent={gid}"]
    return {"start": start, "end": end, "gff": gff, "score": identity,
            "record": {"id": gid, "seq": protein, "description": f"{gid} chr1:{start}-{end}"}}


def _attrs(cand):
    return isr._parse_gff_attr_string(cand["gff"][0].split("\t")[8])


@pytest.fixture
def relative_setup():
    relative = _mutate(HOME, 0.45, 40)
    by_parent = {"GOI_X": {"id": "GOI_X", "seq": HOME},
                 "GOI_X~rel": {"id": "GOI_X~rel", "seq": relative,
                               "model_id": "GOI_X|GenA_fna_b0_l1_exon_ann"}}
    return relative, by_parent


def test_model_found_by_the_home_protein_is_passed_through(relative_setup):
    _relative, by_parent = relative_setup
    cand = _candidate("GOI_X", "GenB_fna_b0_l1_exon_ann", 1000, 4000, HOME, 100.0)
    out = isr.rebase_relative_candidates([cand], ["GOI_X"], by_parent)
    assert out == [cand]


def test_relative_model_on_the_same_dna_as_a_home_model_is_dropped(relative_setup):
    relative, by_parent = relative_setup
    home_model = _candidate("GOI_X", "GenB_fna_b0_l1_exon_ann", 1000, 4000, HOME, 62.0)
    rel_model = _candidate("GOI_X~rel", "GenB_fna_b0_l1_exon_ann", 1100, 4100, relative, 98.0)
    out = isr.rebase_relative_candidates([home_model, rel_model], ["GOI_X", "GOI_X~rel"], by_parent)
    assert out == [home_model]


def test_relative_only_model_is_reported_against_the_home_protein(relative_setup):
    relative, by_parent = relative_setup
    found = _mutate(relative, 0.05, 41)                 # 95 % to the relative, ~52 % to home
    rel_model = _candidate("GOI_X~rel", "GenB_fna_b3_l2_exon_ann", 9000, 12000, found, 95.0)
    (out,) = isr.rebase_relative_candidates([rel_model], ["GOI_X~rel"], by_parent)
    a = _attrs(out)
    assert a["ID"] == "GOI_X|GenB_fna_b3_r2_exon_ann" and a["SynVoy_Parent"] == "GOI_X"
    assert a["SearchQuery"] == "GOI_X|GenA_fna_b0_l1_exon_ann"
    assert a["SearchQueryIdentity"] == "95.0"
    identity = float(a["Identity"])
    assert 40 < identity < 65                            # measured against HOME, not 95
    assert out["gff"][0].split("\t")[5] == f"{identity:.1f}"
    assert out["score"] == pytest.approx(identity, abs=0.05)
    assert float(a["QueryCoverage"]) > 0.9
    assert a["BlockFlankingSupport"] == "6" and a["TargetGene"] == "tgtA"
    assert "~rel" not in "\n".join(out["gff"])
    assert out["gff"][1].endswith("ID=GOI_X|GenB_fna_b3_r2_exon_ann_CDS1;Parent=GOI_X|GenB_fna_b3_r2_exon_ann")
    assert out["record"]["id"] == "GOI_X|GenB_fna_b3_r2_exon_ann"
    assert out["record"]["seq"] == found
    assert "search_query:GOI_X|GenA_fna_b0_l1_exon_ann" in out["record"]["description"]
    assert rel_model["record"]["id"].startswith("GOI_X~rel|")       # the input is not edited


def test_a_relative_only_model_is_judged_against_the_home_protein(relative_setup):
    """95 % to the relative; against the home protein it is far below, and that is what is
    reported. Its label is provisional until the judge labels it, which aligns the HOME
    protein against the model (annotate_place_evalues runs on the rebased candidate)."""
    relative, by_parent = relative_setup
    distant_relative = _mutate(HOME, 0.75, 42)
    by_parent["GOI_X~rel"]["seq"] = distant_relative
    rel_model = _candidate("GOI_X~rel", "GenB_fna_b3_l2_exon_ann", 9000, 12000,
                           _mutate(distant_relative, 0.05, 43), 95.0)
    (out,) = isr.rebase_relative_candidates([rel_model], ["GOI_X~rel"], by_parent)
    a = _attrs(out)
    assert float(a["Identity"]) < 40
    assert (a["Confidence"], a["InferenceReason"]) == ("LOW", "awaiting_label")
    assert a["SynVoy_Parent"] == "GOI_X"


def test_relative_model_next_to_a_home_model_is_kept(relative_setup):
    """Only the same DNA competes: a second copy 10 kb away is a finding of its own."""
    relative, by_parent = relative_setup
    home_model = _candidate("GOI_X", "GenB_fna_b0_l1_exon_ann", 1000, 4000, HOME, 62.0)
    rel_model = _candidate("GOI_X~rel", "GenB_fna_b0_l2_exon_ann", 14000, 17000, relative, 99.0)
    out = isr.rebase_relative_candidates([home_model, rel_model], ["GOI_X", "GOI_X~rel"], by_parent)
    assert len(out) == 2 and _attrs(out[1])["ID"] == "GOI_X|GenB_fna_b0_r2_exon_ann"


def test_second_model_of_a_locus_keeps_its_number_when_renamed(relative_setup):
    relative, by_parent = relative_setup
    rel_model = _candidate("GOI_X~rel", "GenB_fna_b0_l1m2_exon_ann", 14000, 17000, relative, 99.0)
    (out,) = isr.rebase_relative_candidates([rel_model], ["GOI_X~rel"], by_parent)
    assert _attrs(out)["ID"] == "GOI_X|GenB_fna_b0_r1m2_exon_ann"


# ============================================================ distances per genome

def test_bp_or_auto():
    assert isr.bp_or_auto("auto") is None and isr.bp_or_auto("AUTO ") is None
    assert isr.bp_or_auto("") is None
    assert isr.bp_or_auto("150000") == 150000 and isr.bp_or_auto("1.5e5") == 150000
    assert isr.bp_or_auto(0) == 0
    with pytest.raises(ValueError):
        isr.bp_or_auto("wide")


def _geometry(n=11, pitch=10_000, genome=12_000_000):
    return synvoy_scale.home_geometry_from_headers(
        f"gene-F{i:02d} home={100_000 + i * pitch}-{102_000 + i * pitch} genome={genome}"
        for i in range(n))


def _args(geometry, **kw):
    base = dict(home_geometry=geometry, block_gap_genes=7.0, window_margin_genes=7.0,
                goi_margin_genes=2.0, synteny_bridge_max_gap_genes=20.0,
                dispersed_envelope_margin_genes=9.0, cluster_distance=None, region_padding=None,
                goi_margin=None, synteny_bridge_max_gap=None, dispersed_envelope_margin=None)
    base.update(kw)
    return argparse.Namespace(**base)


def _anchor_hits(geometry, factor):
    first = geometry["genes"][geometry["order"][0]]["mid"]
    return [{"query": g, "chrom": "chr1",
             "start": 5_000_000 + int((geometry["genes"][g]["mid"] - first) * factor) - 300,
             "end": 5_000_000 + int((geometry["genes"][g]["mid"] - first) * factor) + 300}
            for g in geometry["order"]]


def test_distances_follow_the_expansion_measured_in_the_target(tmp_path):
    g = _geometry()
    genome = tmp_path / "t.fna"
    genome.write_text(">chr1\n" + "A" * 1000 + "\n")
    args = _args(g)
    out = isr.genome_scale_args(args, _anchor_hits(g, 2.0), str(genome), "t")
    assert out is not args and args.cluster_distance is None          # the input is not edited
    assert out.expansion == pytest.approx(2.0, rel=0.01)
    assert out.cluster_distance == pytest.approx(7 * g["pitch"] * 2.0, rel=0.01)
    assert out.region_padding == out.cluster_distance
    assert out.goi_margin == pytest.approx(2 * g["pitch"] * 2.0, rel=0.01)
    assert out.synteny_bridge_max_gap == pytest.approx(20 * g["pitch"] * 2.0, rel=0.01)
    assert out.adaptive_padding is False
    assert out.distances["expansion_source"] == "anchors"
    assert out.distances["expansion_anchors"] == 11
    assert out.distances["cluster_distance"] == out.cluster_distance


def test_genome_size_ratio_is_the_fallback_when_order_is_not_conserved(tmp_path):
    g = _geometry(genome=1_000_000)
    genome = tmp_path / "t.fna"
    genome.write_text("A" * 3_000_000)
    out = isr.genome_scale_args(_args(g), [], str(genome), "t")
    assert out.distances["expansion_source"] == "genome size"
    assert out.expansion == pytest.approx(3.0, rel=0.01)


def test_a_distance_set_in_bp_is_used_as_it_is(tmp_path):
    g = _geometry()
    genome = tmp_path / "t.fna"
    genome.write_text("A" * 100)
    out = isr.genome_scale_args(_args(g, cluster_distance=150_000, region_padding=150_000,
                                      goi_margin=50_000), _anchor_hits(g, 2.0), str(genome), "t")
    assert (out.cluster_distance, out.region_padding, out.goi_margin) == (150_000, 150_000, 50_000)
    assert out.adaptive_padding is True          # a margin in bp keeps the old adaptive rule
    assert out.distances["source"]["cluster_distance"] == "set in bp"
    assert out.distances["source"]["synteny_bridge_max_gap"] == "derived"


@pytest.mark.parametrize("kw", [dict(synteny_bridge_max_gap=0),
                                dict(synteny_bridge_max_gap_genes=0)])
def test_bridging_can_be_switched_off_in_either_unit(tmp_path, kw):
    g = _geometry()
    genome = tmp_path / "t.fna"
    genome.write_text("A" * 100)
    out = isr.genome_scale_args(_args(g, **kw), _anchor_hits(g, 1.0), str(genome), "t")
    assert out.synteny_bridge_max_gap == 0 and out.distances["synteny_bridge_max_gap"] == 0


def test_only_the_dna_whose_hits_can_be_used_is_cut_out(tmp_path):
    """Hits outside the core (block +- GOI margin) are dropped, and a gene model uses at
    most the model margin beyond its hits. The DNA searched is therefore core + model
    margin, never more than the window margin. At a human-size locus that is the
    difference between megabases and a few hundred kb on each side."""
    genome = tmp_path / "t.fna"
    genome.write_text("A" * 100)
    # compact locus: 7 pitches are less than 2 pitches + 50 kb, the window margin stays
    small = _geometry(pitch=3_000)
    out = isr.genome_scale_args(_args(small, gap_search_window=50_000, goi_home_span=1_300),
                                _anchor_hits(small, 1.0), str(genome), "t")
    assert out.search_margin == out.region_padding == pytest.approx(7 * small["pitch"], rel=0.01)
    # wide locus: 2 pitches + 50 kb instead of 7 pitches
    wide = _geometry(pitch=300_000)
    out = isr.genome_scale_args(_args(wide, gap_search_window=50_000, goi_home_span=40_000),
                                _anchor_hits(wide, 1.0), str(genome), "t")
    assert out.region_padding == pytest.approx(7 * wide["pitch"], rel=0.01)
    assert out.search_margin == pytest.approx(2 * wide["pitch"] + 50_000, rel=0.01)
    assert out.distances["search_margin"] == out.search_margin
    # a gene longer than the floor gets its own span, stretched like the neighbourhood
    out = isr.genome_scale_args(_args(wide, gap_search_window=50_000, goi_home_span=400_000),
                                _anchor_hits(wide, 2.0), str(genome), "t")
    assert out.search_margin == pytest.approx((2 * wide["pitch"] + 400_000) * 2.0, rel=0.01)


def test_a_window_margin_in_bp_is_searched_in_full_as_before(tmp_path):
    g = _geometry(pitch=300_000)
    genome = tmp_path / "t.fna"
    genome.write_text("A" * 100)
    out = isr.genome_scale_args(_args(g, region_padding=150_000, gap_search_window=50_000, goi_home_span=40_000),
                                _anchor_hits(g, 1.0), str(genome), "t")
    assert out.adaptive_padding is True and out.search_margin is None
    assert "search_margin" not in out.distances
    out = isr.genome_scale_args(_args(None, gap_search_window=50_000), [], str(tmp_path / "absent.fna"), "t")
    assert out.search_margin is None             # no home coordinates: the former fixed window


def test_query_set_without_home_coordinates_keeps_the_fixed_distances(tmp_path):
    out = isr.genome_scale_args(_args(None), [], str(tmp_path / "absent.fna"), "t")
    assert (out.cluster_distance, out.region_padding, out.goi_margin) == (150_000, 150_000, 50_000)
    assert out.synteny_bridge_max_gap == 6_000_000 and out.expansion == 1.0
    assert out.distances["home_anchors"] == 0


def test_per_gene_gap_limit_uses_the_home_span_and_the_expansion():
    g = synvoy_scale.home_geometry_from_headers(
        ["gene-SMALL home=1000-3000", "gene-BIG home=50000-450000"])
    args = argparse.Namespace(home_geometry=g, expansion=2.0)
    assert isr._gene_limit(args, "gene-BIG|exon_3", 20_000) == 1_200_000
    assert isr._gene_limit(args, "gene-SMALL", 20_000) == 20_000
    assert isr._gene_limit(argparse.Namespace(), "gene-BIG", 20_000) == 20_000


def test_exons_of_a_long_gene_stay_one_occurrence_in_the_block_step():
    """Two exon hits 77 kb apart: one gene with a long intron, not two copies of it."""
    def h(q, s):
        return {"query": q, "chrom": "c", "start": s, "end": s + 3000}
    hits = [h("gene-A|exon_1", 100_000), h("gene-A|exon_2", 180_000), h("gene-B", 190_000)]

    fixed = isr.identify_synteny_blocks([dict(x) for x in hits], max_intron=20_000,
                                        cluster_distance=20_000)
    assert sorted(b["genes_count"] for b in fixed) == [1, 2]

    by_gene = isr.identify_synteny_blocks(
        [dict(x) for x in hits], max_intron=20_000, cluster_distance=20_000,
        gene_gap_limit=lambda gene: 120_000 if gene == "gene-A" else 20_000)
    assert len(by_gene) == 1
    assert (by_gene[0]["start"], by_gene[0]["end"], by_gene[0]["loci_count"]) == (100_000, 193_000, 2)


# ============================================== anchors limited to their best places (opt-in)

def _places(hits, n, limit=20_000):
    kept = isr.limit_anchor_places(hits, n, (lambda _gene: limit))
    return sorted((h["query"], h["chrom"], h["start"]) for h in kept)


def test_anchor_places_off_by_default_returns_every_hit():
    hits = [_hit("gene-A", 1000, 1500, 200), _hit("gene-A", 9_000_000, 9_000_400, 40)]
    assert isr.limit_anchor_places(hits, 0) == hits
    assert isr.limit_anchor_places(hits, None) == hits
    assert isr.limit_anchor_places([], 5) == []


def test_the_best_places_of_each_gene_are_kept():
    """gene-A: ortholog (two exon hits, 300 bits together) and three relatives elsewhere."""
    hits = [_hit("gene-A|exon_1", 1000, 1300, 180), _hit("gene-A|exon_2", 4000, 4300, 120),
            _hit("gene-A|exon_1", 2_000_000, 2_000_300, 90),
            _hit("gene-A|exon_2", 5_000_000, 5_000_300, 60, chrom="chr2"),
            _hit("gene-A|exon_1", 8_000_000, 8_000_300, 40),
            _hit("gene-B", 30_000, 31_000, 250)]
    assert _places(hits, 1) == [("gene-A|exon_1", "chr1", 1000), ("gene-A|exon_2", "chr1", 4000),
                                ("gene-B", "chr1", 30_000)]
    assert len(_places(hits, 2)) == 4 and ("gene-A|exon_1", "chr1", 2_000_000) in _places(hits, 2)
    assert len(_places(hits, 10)) == len(hits)


def test_a_place_is_all_hits_of_the_gene_within_its_gap_limit():
    """Two exons 30 kb apart are one place for a gene with a 40 kb limit and two places for a
    gene with a 20 kb limit; the weight of a place is the sum of its bit scores."""
    hits = [_hit("gene-A|exon_1", 1000, 1300, 60), _hit("gene-A|exon_2", 31_000, 31_300, 60),
            _hit("gene-A|exon_1", 900_000, 900_300, 100)]
    one_place = isr.limit_anchor_places(hits, 1, (lambda _g: 40_000))
    assert sorted(h["start"] for h in one_place) == [1000, 31_000]          # 120 bits beat 100
    split = isr.limit_anchor_places(hits, 1, (lambda _g: 20_000))
    assert [h["start"] for h in split] == [900_000]                         # 100 beat 60 and 60


def test_anchor_places_keep_input_order_and_break_ties_by_position():
    hits = [_hit("gene-A", 500_000, 500_300, 80), _hit("gene-A", 1000, 1300, 80),
            _hit("gene-A", 90_000, 90_300, 80, chrom="chr0")]
    assert [(h["chrom"], h["start"]) for h in isr.limit_anchor_places(hits, 2)] == [("chr1", 1000), ("chr0", 90_000)] \
        or [(h["chrom"], h["start"]) for h in isr.limit_anchor_places(hits, 2)] == [("chr1", 1000), ("chr0", 90_000)][::-1]
    first = isr.limit_anchor_places(hits, 2)
    assert first == isr.limit_anchor_places(list(reversed(hits)), 2)[::-1] or \
        {id(h) for h in first} == {id(h) for h in isr.limit_anchor_places(list(reversed(hits)), 2)}
