"""Distances in gene pitches (bin/synvoy_scale.py).

The search states "how far apart" in steps between neighbouring anchors of the home
locus and converts them to base pairs per target genome. These tests pin the three
parts: the geometry read from the query set, the expansion measured in a target, and
the conversion with its fall-backs.
"""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bin"))

import synvoy_scale as sc  # noqa: E402


def _headers(n=11, pitch=10_000, gene_len=2_000, genome=12_000_000, start=100_000):
    """n flanking genes, one every `pitch` bp, as PREPARE_INITIAL_DB writes them."""
    out = []
    for i in range(n):
        s = start + i * pitch
        out.append(f"gene-F{i:02d} home={s}-{s + gene_len} genome={genome}")
    return out


def _geometry(**kw):
    return sc.home_geometry_from_headers(_headers(**kw))


def _hits(geometry, factor=1.0, chrom="chr1", offset=5_000_000, skip=()):
    """One hit per anchor, the neighbourhood stretched by `factor`."""
    first = geometry["genes"][geometry["order"][0]]["mid"]
    hits = []
    for gene in geometry["order"]:
        if gene in skip:
            continue
        mid = offset + int((geometry["genes"][gene]["mid"] - first) * factor)
        hits.append({"query": gene, "chrom": chrom, "start": mid - 300, "end": mid + 300})
    return hits


# --------------------------------------------------------------------------- geometry

def test_geometry_pitch_is_span_over_steps():
    g = _geometry(n=11, pitch=10_000, gene_len=2_000)
    assert g["n"] == 11
    assert g["span"] == 10 * 10_000 + 2_000
    assert g["pitch"] == pytest.approx(g["span"] / 10)
    assert g["max_gap"] == 8_000
    assert g["home_genome_bp"] == 12_000_000
    assert g["order"][0] == "gene-F00" and g["order"][-1] == "gene-F10"


def test_geometry_ignores_the_goi_and_its_exons():
    headers = _headers(n=4) + ["GOI_X home=1-9000000 genome=5", "exon_1 coords=1-90"]
    g = sc.home_geometry_from_headers(headers)
    assert g["n"] == 4 and "GOI_X" not in g["genes"]


def test_geometry_from_exon_records_gives_span_and_longest_intron():
    headers = [
        "gene-A|exon_1 coords=1000-1200",
        "gene-A|exon_2 coords=5200-5400",
        "gene-A|exon_3 coords=5600-5900",
        "gene-B|exon_1 coords=20000-20500",
    ]
    g = sc.home_geometry_from_headers(headers)
    assert g["genes"]["gene-A"]["start"] == 1000 and g["genes"]["gene-A"]["end"] == 5900
    assert g["genes"]["gene-A"]["max_intron"] == 4000
    assert g["genes"]["gene-B"]["max_intron"] == 0
    assert g["home_genome_bp"] is None


def test_gene_header_wins_over_exon_records_for_the_span():
    headers = [
        "gene-A home=900-6000 genome=100",
        "gene-A|exon_1 coords=1000-1200",
        "gene-A|exon_2 coords=5200-5400",
        "gene-B home=9000-9500",
    ]
    g = sc.home_geometry_from_headers(headers)
    assert (g["genes"]["gene-A"]["start"], g["genes"]["gene-A"]["end"]) == (900, 6000)
    assert g["genes"]["gene-A"]["max_intron"] == 4000


@pytest.mark.parametrize("headers", [
    [],
    ["gene-A home=10-500"],                       # one anchor: no step to measure
    ["gene-A", "gene-B some description"],        # query set of an older run
    ["gene-A home=500-10", "gene-B home=7-3"],    # unusable coordinates
])
def test_no_geometry_without_two_located_anchors(headers):
    assert sc.home_geometry_from_headers(headers) is None


def test_geometry_from_a_fasta_file_and_from_a_missing_file(tmp_path):
    fa = tmp_path / "db.faa"
    fa.write_text("".join(f">{h}\nMKT\n" for h in _headers(n=3)))
    assert sc.home_geometry(str(fa))["n"] == 3
    assert sc.home_geometry(str(tmp_path / "absent.faa")) is None


# -------------------------------------------------------------------------- expansion

@pytest.mark.parametrize("factor", [0.5, 1.0, 2.0, 4.0])
def test_expansion_is_the_stretch_of_the_neighbourhood(factor):
    g = _geometry()
    e, anchors, source = sc.target_expansion(g, _hits(g, factor))
    assert source == "anchors"
    assert e == pytest.approx(factor, rel=0.01)
    assert anchors == 11        # all eleven lie in one run


def test_expansion_survives_one_anchor_moved_far_away():
    g = _geometry()
    hits = _hits(g, 2.0)
    hits[5]["start"] += 40_000_000
    hits[5]["end"] += 40_000_000
    e, anchors, source = sc.target_expansion(g, hits)
    assert source == "anchors" and e == pytest.approx(2.0, rel=0.01)
    assert anchors == 10        # the run steps over the missing one


def test_scrambled_genome_falls_back_to_the_genome_size_ratio():
    """Anchors megabases apart say nothing about the size of the neighbourhood."""
    g = _geometry()
    hits = [{"query": gene, "chrom": "chr1", "start": 1_000_000 + i * 3_000_000,
             "end": 1_000_600 + i * 3_000_000} for i, gene in enumerate(g["order"])]
    e, anchors, source = sc.target_expansion(g, hits, size_ratio=1.7)
    assert (e, anchors, source) == (1.7, 1, "genome size")      # no two anchors form a run


def test_few_conserved_anchors_among_many_are_not_trusted():
    """Three anchors in a row plus eight scattered: a run of three, below the minimum."""
    g = _geometry()
    near = _hits(g, 1.0)[:3]
    far = [{"query": gene, "chrom": "chr1", "start": 50_000_000 + i * 5_000_000,
            "end": 50_000_600 + i * 5_000_000} for i, gene in enumerate(g["order"][3:])]
    e, anchors, source = sc.target_expansion(g, near + far, size_ratio=0.8)
    assert source == "genome size" and e == 0.8 and anchors == 3


def test_anchors_on_different_sequences_give_no_run():
    g = _geometry()
    hits = [dict(h, chrom=f"scaf{i}") for i, h in enumerate(_hits(g))]
    assert sc.target_expansion(g, hits) == (1.0, 1, "none")
    assert sc.target_expansion(g, []) == (1.0, 0, "none")


def test_genome_size_ratio_is_kept_inside_the_bounds():
    assert sc.target_expansion(None, [], size_ratio=40.0)[0] == sc.EXPANSION_BOUNDS[1]
    assert sc.target_expansion(None, [], size_ratio=0.01)[0] == sc.EXPANSION_BOUNDS[0]
    assert sc.target_expansion(None, []) == (1.0, 0, "none")


def test_closest_copy_of_a_duplicated_anchor_is_used():
    """A second, distant hit of an anchor (a relative elsewhere) must not widen the estimate."""
    g = _geometry()
    hits = _hits(g, 1.0)
    extra = [dict(h, start=h["start"] + 30_000_000, end=h["end"] + 30_000_000) for h in hits[::2]]
    e, _pairs, source = sc.target_expansion(g, hits + extra)
    assert source == "anchors" and e == pytest.approx(1.0, rel=0.01)


def test_hits_of_exon_queries_count_for_their_gene():
    g = _geometry()
    hits = [dict(h, query=h["query"] + "|exon_2") for h in _hits(g, 1.5)]
    e, _pairs, source = sc.target_expansion(g, hits)
    assert source == "anchors" and e == pytest.approx(1.5, rel=0.01)


# ------------------------------------------------------------------------- conversion

def test_distances_scale_with_pitch_and_expansion():
    g = _geometry()                      # pitch 10.2 kb, largest gap 8 kb
    d1 = sc.derive_distances(g, 1.0)
    d2 = sc.derive_distances(g, 2.0)
    assert d1["cluster_distance"] == round(7 * g["pitch"])
    assert d1["region_padding"] == round(7 * g["pitch"])
    assert d1["goi_margin"] == round(2 * g["pitch"])
    assert d1["synteny_bridge_max_gap"] == round(20 * g["pitch"])
    assert d1["dispersed_envelope_margin"] == round(9 * g["pitch"])
    for key in ("cluster_distance", "region_padding", "goi_margin",
                "synteny_bridge_max_gap", "dispersed_envelope_margin"):
        assert d2[key] == pytest.approx(2 * d1[key], abs=1)
        assert d1["source"][key] == "derived"
    assert d1["home_anchors"] == 11 and d1["home_pitch_bp"] == round(g["pitch"])


def test_block_gap_never_below_the_largest_home_gap():
    """One gene desert at home: an unchanged neighbourhood must stay one block."""
    headers = _headers(n=21, pitch=10_000) + ["gene-FAR home=600000-602000"]
    g = sc.home_geometry_from_headers(headers)
    d = sc.derive_distances(g, 1.0)
    assert g["max_gap"] == 600_000 - (100_000 + 20 * 10_000 + 2_000)
    assert d["cluster_distance"] == round(1.25 * g["max_gap"])
    assert d["cluster_distance"] > 7 * g["pitch"]
    assert d["region_padding"] == round(7 * g["pitch"])        # only the block gap is raised


def test_a_value_in_bp_overrides_the_pitch_for_that_distance_only():
    g = _geometry()
    d = sc.derive_distances(g, 3.0, fixed_bp={"cluster_distance": 150_000, "region_padding": None,
                                              "gap_search_window": 50_000})
    assert d["cluster_distance"] == 150_000 and d["source"]["cluster_distance"] == "set in bp"
    assert d["goi_margin"] == 50_000 and d["source"]["goi_margin"] == "set in bp"
    assert d["region_padding"] == round(7 * g["pitch"] * 3.0)
    assert d["source"]["region_padding"] == "derived"


def test_gene_units_can_be_changed():
    g = _geometry()
    d = sc.derive_distances(g, 1.0, genes={"block_gap_genes": 3, "bridge_max_gap_genes": 40,
                                           "window_margin_genes": None})
    assert d["cluster_distance"] == round(max(3 * g["pitch"], 1.25 * g["max_gap"]))
    assert d["synteny_bridge_max_gap"] == round(40 * g["pitch"])
    assert d["region_padding"] == round(7 * g["pitch"])        # None = the default


def test_without_geometry_every_distance_keeps_its_former_value():
    d = sc.derive_distances(None, 2.5)
    assert d["cluster_distance"] == 150_000
    assert d["region_padding"] == 150_000
    assert d["goi_margin"] == 50_000
    assert d["synteny_bridge_max_gap"] == 6_000_000
    assert d["dispersed_envelope_margin"] == 200_000
    assert set(d["source"].values()) == {"fixed (no home coordinates)"}
    assert d["home_pitch_bp"] is None and d["home_anchors"] == 0


def test_gene_gap_limit_follows_the_home_span_of_the_gene():
    headers = ["gene-SMALL home=1000-3000", "gene-BIG home=50000-450000"]
    g = sc.home_geometry_from_headers(headers)
    assert sc.gene_gap_limit(g, "gene-SMALL", 1.0, 20_000) == 20_000          # floor
    assert sc.gene_gap_limit(g, "gene-BIG", 1.0, 20_000) == 600_000           # 1.5 x 400 kb
    assert sc.gene_gap_limit(g, "gene-BIG", 0.5, 20_000) == 300_000
    assert sc.gene_gap_limit(g, "gene-UNKNOWN", 1.0, 20_000) == 20_000
    assert sc.gene_gap_limit(None, "gene-BIG", 1.0, 20_000) == 20_000


# ------------------------------------------------------------------- result GFF header

def test_header_round_trip(tmp_path):
    d = sc.derive_distances(_geometry(), 1.26)
    d["expansion_source"] = "anchors"
    d["expansion_anchors"] = 11
    gff = tmp_path / "g.gff"
    gff.write_text("##gff-version 3\n" + sc.format_gff_header(d) + "\nchr1\tSynVoy\tmRNA\t1\t9\t.\t+\t.\tID=x\n")
    back = sc.parse_gff_header(str(gff))
    assert back == json.loads(json.dumps(d))
    assert back["cluster_distance"] == d["cluster_distance"]


def test_header_line_is_one_line_and_stable():
    d = sc.derive_distances(_geometry(), 1.0)
    line = sc.format_gff_header(d)
    assert line.startswith("##synvoy-distances {") and "\n" not in line
    assert line == sc.format_gff_header(dict(reversed(list(d.items()))))     # key order free


@pytest.mark.parametrize("text", [
    "##gff-version 3\nchr1\tx\tmRNA\t1\t9\t.\t+\t.\tID=x\n",             # GFF of an older run
    "chr1\tx\tmRNA\t1\t9\t.\t+\t.\tID=x\n##synvoy-distances {\"cluster_distance\":5}\n",  # not in the head
    "##synvoy-distances {broken\n",
    "",
])
def test_no_header_means_none(tmp_path, text):
    gff = tmp_path / "g.gff"
    gff.write_text(text)
    assert sc.parse_gff_header(str(gff)) is None


def test_missing_gff_means_none(tmp_path):
    assert sc.parse_gff_header(str(tmp_path / "absent.gff")) is None


def test_genes_between_counts_only_genes_entirely_inside():
    genes = [(100, 200), (300, 400), (450, 700), (800, 900)]
    assert sc.genes_between(genes, 200, 800) == 2
    assert sc.genes_between(genes, 150, 800) == 2      # the first gene starts before
    assert sc.genes_between(genes, 200, 650) == 1      # the third gene ends after
    assert sc.genes_between(genes, 900, 2000) == 0


# ------------------------------------------------ expansion: like for like, and only where
# ------------------------------------------------ the gene order is conserved (2026-10-09)

def _exon_headers(n=8, pitch=6_000, exons=((0, 300), (1_500, 1_800), (3_600, 4_000))):
    """n flanking genes of three exons each, one gene every `pitch` bp: compact genes with
    introns, as at the fly Defensin locus where the estimate came out a quarter too small."""
    out = []
    for i in range(n):
        s = 200_000 + i * pitch
        out.append(f"gene-M{i:02d} home={s}-{s + exons[-1][1]} genome=140000000")
        for k, (a, b) in enumerate(exons, start=1):
            out.append(f"gene-M{i:02d}|exon_{k} coords={s + a}-{s + b}")
    return out


def _exon_hits(headers, factor=1.0, offset=9_000_000, chrom="chr2"):
    """One hit per exon record, the whole neighbourhood stretched by `factor`."""
    hits = []
    for h in headers:
        rid, coords = h.split()[0], h.split()[1]
        if "|exon_" not in rid:
            continue
        a, b = (int(x) for x in coords.split("=")[1].split("-"))
        mid = offset + int(((a + b) // 2 - 200_000) * factor)
        hits.append({"query": rid, "chrom": chrom, "start": mid - 100, "end": mid + 100})
    return hits


def test_geometry_remembers_where_each_query_record_lies():
    g = sc.home_geometry_from_headers(_exon_headers(n=3))
    assert g["records"]["gene-M00|exon_1"] == 200_150
    assert g["records"]["gene-M00|exon_3"] == 203_800
    assert g["records"]["gene-M00"] == g["genes"]["gene-M00"]["mid"]
    assert not any(r.startswith("GOI") for r in g["records"])


@pytest.mark.parametrize("factor", [0.5, 1.0, 2.5])
def test_exon_level_hits_give_the_true_stretch(factor):
    """Both distances between the same two exons. Dividing the distance of the nearest
    exons by the distance of the two gene midpoints gave 0.74 for a neighbourhood that is
    exactly as long as at home (D. simulans, Defensin locus)."""
    headers = _exon_headers()
    g = sc.home_geometry_from_headers(headers)
    e, _pairs, source = sc.target_expansion(g, _exon_hits(headers, factor))
    assert source == "anchors"
    assert e == pytest.approx(factor, rel=0.02)


def test_gene_midpoints_instead_of_exon_positions_shorten_the_run_at_home_only_at_its_ends():
    """Without the positions of the exon records every hit counts for its gene midpoint. The
    measure of the first version (nearest exons of two neighbours against their gene
    midpoints, pair by pair) came out a quarter too small that way. Along a run the offsets
    cancel except at its two ends, so the quotient is off by a few percent; with the record
    positions it is exact (the test above)."""
    headers = _exon_headers()
    g = sc.home_geometry_from_headers(headers)
    by_gene = dict(g, records={})
    e_by_gene, _n, _s = sc.target_expansion(by_gene, _exon_hits(headers, 1.0))
    assert 0.85 < e_by_gene < 1.0
    e_exact, _n, _s = sc.target_expansion(g, _exon_hits(headers, 1.0))
    assert abs(e_exact - 1.0) < abs(e_by_gene - 1.0)


def _points(positions):
    """[(target position, home rank)] on one sequence -> points for longest_anchor_chain,
    home positions 10 kb apart."""
    return {"chr1": [(t, 100_000 + r * 10_000, r) for t, r in positions]}


def test_chain_counts_anchors_in_a_row_in_either_direction():
    forward = _points([(1_000_000 + i * 10_000, i) for i in range(8)])
    assert sc.longest_anchor_chain(forward) == 8
    inverted = _points([(1_000_000 + i * 10_000, 7 - i) for i in range(8)])
    assert sc.longest_anchor_chain(inverted) == 8


def test_chain_skips_up_to_three_ranks_and_no_more():
    with_gaps = _points([(1_000_000, 0), (1_030_000, 3), (1_050_000, 5), (1_060_000, 6)])
    assert sc.longest_anchor_chain(with_gaps) == 4
    too_far = _points([(1_000_000, 0), (1_040_000, 4), (1_050_000, 5)])
    assert sc.longest_anchor_chain(too_far) == 2


def test_chain_needs_comparable_distances_and_one_direction():
    far_apart = _points([(1_000_000 + i * 2_000_000, i) for i in range(8)])      # 200 x home
    assert sc.longest_anchor_chain(far_apart) == 1
    zigzag = _points([(1_000_000, 0), (1_010_000, 2), (1_020_000, 1), (1_030_000, 3), (1_040_000, 2)])
    assert sc.longest_anchor_chain(zigzag) <= 3


def test_chain_is_per_sequence_and_the_search_limit_changes_nothing():
    pts = {"a": [(1_000_000 + i * 10_000, 100_000 + i * 10_000, i) for i in range(4)],
           "b": [(5_000_000 + i * 10_000, 100_000 + (i + 4) * 10_000, i + 4) for i in range(4)]}
    assert sc.longest_anchor_chain(pts) == 4
    long_run = _points([(1_000_000 + i * 12_000, i) for i in range(12)])
    assert sc.longest_anchor_chain(long_run) == sc.longest_anchor_chain(long_run, max_home_step=40_000) == 12
    assert sc.longest_anchor_chain({}) == 0


def test_chance_pairs_without_a_run_of_five_are_not_trusted():
    """Firefly against fruit fly, 2026-10-09: 11 pairs of anchors lay within eight times
    their home distance of each other by chance, scattered over the chromosome arms. That
    passed the former rule and gave a factor of 3.7, a block gap of 956 kb and a search
    window that needed 10 GB. No five of them were in a row."""
    g = _geometry(n=20)
    order = g["order"]
    hits = []
    # six pairs of neighbours, each pair at 4-7 times its home distance, pairs megabases apart
    for k, i in enumerate(range(0, 12, 2)):
        base = 2_000_000 + k * 6_000_000
        gap = 10_000 * (4 + k % 4)
        hits.append({"query": order[i], "chrom": "arm", "start": base, "end": base + 600})
        hits.append({"query": order[i + 1], "chrom": "arm", "start": base + gap, "end": base + gap + 600})
    e, anchors, source = sc.target_expansion(g, hits, size_ratio=0.31)
    assert anchors == 2 and source == "genome size" and e == 0.31


def test_four_anchors_in_a_row_are_not_enough_five_are():
    g = _geometry(n=12)
    assert sc.target_expansion(g, _hits(g, 3.0)[:4], size_ratio=0.5)[2] == "genome size"
    e, _pairs, source = sc.target_expansion(g, _hits(g, 3.0)[:5], size_ratio=0.5)
    assert source == "anchors" and e == pytest.approx(3.0, rel=0.01)



# ------------------------------------------------- the expansion is a quotient of two lengths

def _patchy_hits(geometry, steps, chrom="chr1", offset=5_000_000):
    """One hit per anchor; `steps[i]` is the factor by which the distance from anchor i to
    anchor i + 1 has grown."""
    order = geometry["order"]
    mids = [geometry["genes"][g]["mid"] for g in order]
    pos, hits = offset, []
    for i, gene in enumerate(order):
        if i:
            pos += int((mids[i] - mids[i - 1]) * steps[i - 1])
        hits.append({"query": gene, "chrom": chrom, "start": pos - 300, "end": pos + 300})
    return hits


def test_a_few_large_insertions_count_in_full():
    """Most neighbours keep their distance and three insertions make the neighbourhood 2.5
    times as long. The median of the pair ratios, the measure until 2026-10-11, was about 1
    here (A. lyrata at the FT locus: 1.5 for a neighbourhood 4.9 times as long)."""
    g = _geometry()                                   # 11 anchors, 10 steps
    steps = [1, 1, 6, 1, 1, 6, 1, 1, 6, 1]
    e, anchors, source = sc.target_expansion(g, _patchy_hits(g, steps))
    assert source == "anchors" and anchors == 11
    assert e == pytest.approx(sum(steps) / len(steps), rel=0.01)      # 2.5


def test_an_insertion_larger_than_the_bound_is_stepped_over():
    """A step of 20 times its home length is outside the bounds, and so is a step of two
    ranks across it (10.5 times). A step of three ranks crosses it (7.3 times), which
    leaves out two anchors next to the insertion."""
    g = _geometry()
    steps = [1, 1, 1, 1, 20, 1, 1, 1, 1, 1]
    e, anchors, source = sc.target_expansion(g, _patchy_hits(g, steps))
    assert source == "anchors"
    assert anchors == 9
    assert sc.EXPANSION_BOUNDS[0] <= e <= sc.EXPANSION_BOUNDS[1]
    assert e == pytest.approx(29 / 10, rel=0.02)      # the whole stretch, the insertion included


def test_an_insertion_that_no_step_can_cross_ends_the_run():
    g = _geometry()
    steps = [1, 1, 1, 1, 1, 60, 1, 1, 1, 1]           # 600 kb between two halves
    e, anchors, source = sc.target_expansion(g, _patchy_hits(g, steps), size_ratio=0.9)
    assert anchors == 6 and source == "anchors"       # the longer half: six anchors, unstretched
    assert e == pytest.approx(1.0, rel=0.01)


def test_the_estimate_does_not_depend_on_the_order_of_the_hits():
    import random
    g = _geometry(n=14)
    steps = [1, 3, 1, 5, 1, 1, 7, 2, 1, 4, 1, 6, 1]
    hits = _patchy_hits(g, steps) + _patchy_hits(g, list(reversed(steps)), chrom="chr9")
    expected = sc.target_expansion(g, hits)
    for seed in range(8):
        shuffled = hits[:]
        random.Random(seed).shuffle(shuffled)
        assert sc.target_expansion(g, shuffled) == expected


def test_anchor_run_returns_both_lengths():
    pts = _points([(1_000_000, 0), (1_020_000, 1), (1_050_000, 2), (1_060_000, 3), (1_100_000, 4)])
    anchors, in_target, at_home = sc.anchor_run(pts)
    assert (anchors, in_target, at_home) == (5, 100_000, 40_000)
    assert sc.anchor_run({}) == (0, 0, 0)
    assert sc.anchor_run({"x": [(5, 100, 0)]}) == (1, 0, 0)


def test_of_two_runs_with_as_many_anchors_the_longer_one_at_home_is_used():
    # five neighbours (40 kb at home) unstretched, and five every second gene (80 kb at home) doubled
    short = [(1_000_000 + i * 10_000, 100_000 + i * 10_000, i) for i in range(5)]
    long_ = [(9_000_000 + i * 40_000, 100_000 + (2 * i + 6) * 10_000, 2 * i + 6) for i in range(5)]
    assert sc.anchor_run({"a": short, "b": long_}) == (5, 160_000, 80_000)
    assert sc.anchor_run({"b": short, "a": long_}) == (5, 160_000, 80_000)
