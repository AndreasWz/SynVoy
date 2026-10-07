"""Figure fixes from the 2026-10-06 review of the latest runs.

- The grid's "xN" badge above a flanking arrow counted every model named after that
  home gene anywhere in the genome. A yeast flanking gene with nine paralog hits on
  six chromosomes was drawn as "x9" above its single-copy ortholog.
- The home GOI label of the synteny plot climbs at 45 degrees into a header of fixed
  height: a long product name ran through the subtitle and off the top edge.
- The home row printed its scaffold twice ("NC_037641.1 . NC_037641.1") and carried
  UniProt's strain qualifier in the species label.
"""
from __future__ import annotations

import os
import re
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402
import sequence_utils as su  # noqa: E402


def _g(name, start, end, chrom="chr1", **kw):
    g = {"chrom": chrom, "start": start, "end": end, "start_plot": start,
         "end_plot": end, "name": name, "strand": "+", "exon_coords": []}
    g.update(kw)
    return g


def _flank(name, start, chrom="chr1", identity=90.0, home="gene-A"):
    return _g(name, start, start + 900, chrom=chrom, home_gene_id=home,
              identity=identity, confidence="HIGH")


def _goi(name, start, chrom="chr1", identity=80.0):
    return _g(name, start, start + 500, chrom=chrom, home_gene_id="GOI_T", role="goi",
              identity=identity, confidence="HIGH")


# ------------------------------------------------------------ grid copies ----

def _paralog_rich_track():
    """gene-A: the ortholog, one tandem copy beside it, and paralog hits elsewhere."""
    return {"genes": [
        _flank("a_ortholog", 100_000, identity=100.0),
        _flank("a_tandem", 104_000, identity=88.0),
        _flank("a_same_chrom_far", 5_000_000, identity=48.0),
        _flank("a_chr2", 600_000, chrom="chr2", identity=39.0),
        _flank("a_chr3", 800_000, chrom="chr3", identity=33.0),
        _flank("a_chr4", 40_000, chrom="chr4", identity=30.0),
        _goi("GOI_T|x", 120_000),
        _flank("b", 140_000, home="gene-B"),
    ]}


def test_flanking_copy_count_ignores_hits_elsewhere_in_the_genome():
    tmap = ps._grid_target_map(_paralog_rich_track(), {"gene-A", "gene-B", "GOI"}, "GOI",
                               locus_span_bp=50_000)
    assert tmap["gene-A"]["best"]["name"] == "a_ortholog"
    assert tmap["gene-A"]["n"] == 2          # the ortholog and its tandem copy
    assert len(tmap["gene-A"]["genes"]) == 6  # nothing is dropped from the track
    assert tmap["gene-B"]["n"] == 1


def test_copy_count_without_a_home_span_falls_back_to_the_row_stretch():
    # The row's best genes span 100-140 kb, so a copy 4 kb away still counts.
    tmap = ps._grid_target_map(_paralog_rich_track(), {"gene-A", "gene-B", "GOI"}, "GOI")
    assert tmap["gene-A"]["n"] == 2


def test_goi_column_still_counts_every_model_of_the_view():
    track = {"genes": [_goi("GOI_T|1", 120_000), _goi("GOI_T|2", 126_000, identity=70.0),
                       _goi("GOI_T|3", 900_000, chrom="chr2", identity=60.0),
                       _flank("a", 100_000)]}
    tmap = ps._grid_target_map(track, {"gene-A", "GOI"}, "GOI", locus_span_bp=50_000)
    assert tmap["GOI"]["n"] == 3


def test_local_copy_needs_the_same_scaffold_and_the_reach():
    best = _flank("best", 100_000)
    assert ps._is_local_copy(best, best, 0)
    assert ps._is_local_copy(_flank("near", 130_000), best, 50_000)
    assert not ps._is_local_copy(_flank("far", 400_000), best, 50_000)
    assert not ps._is_local_copy(_flank("other", 100_500, chrom="chr9"), best, 50_000)


def test_home_locus_span_ignores_unplaced_anchors():
    anchors = [{"start": 1000, "end": 2000}, {"start": 0, "end": 0}, {"start": 61_000, "end": 62_000}]
    assert ps._home_locus_span(anchors) == 61_000
    assert ps._home_locus_span([]) == 0


def test_rendered_grid_prints_the_local_copy_count():
    home = {"label": "Home", "species": "Saccharomyces cerevisiae", "is_home": True,
            "genome_id": "home", "offset": 0, "breaks": [],
            "genes": [_g("gene-A", 1000, 2000), _g("GOI_T", 3000, 3500),
                      _g("gene-B", 50_000, 51_000)]}
    tgt = dict(_paralog_rich_track(), label="T", species="Saccharomyces paradoxus",
               is_home=False, genome_id="spar", offset=0, breaks=[])
    html = ps.render_anchor_grid([home, tgt], {"gene-A": "#4e79a7", "gene-B": "#59a14f"},
                                 {}, {}, SimpleNamespace())
    # Two of the nine models are at the locus: a second arrow behind the first, the
    # count in the tooltip, and no other cell claims copies.
    assert re.findall(r"(\d+) copies at the locus", html) == ["2"]
    assert "more than one copy at the locus" in html


# ------------------------------------------------------- synteny plot labels ----

def _synteny_tracks(home_goi_name):
    def gene(name, start, end, home_gene_id=None, identity=75.0, chrom="chr1"):
        return {"chrom": chrom, "start": start, "end": end, "start_plot": start,
                "end_plot": end, "name": name, "home_gene_id": home_gene_id or name,
                "strand": "+", "identity": identity, "exon_coords": []}
    scaf = "NW_022170249.1"   # the pipeline labels the home track with its scaffold
    return [
        {"label": f"Photinus pyralis ({scaf})", "is_home": True, "genome_id": "home",
         "goi_status": "resolved", "offset": 0, "breaks": [], "minus_strand_row": None,
         "genes": [gene("geneA", 100, 180, chrom=scaf),
                   gene(home_goi_name, 240, 320, identity=99.0, chrom=scaf),
                   gene("geneB", 580, 660, chrom=scaf)]},
        {"label": "Drosophila melanogaster (Drosophila_melanogaster)", "is_home": False,
         "genome_id": "target1", "goi_status": "resolved", "offset": 0, "breaks": [],
         "minus_strand_row": None,
         "genes": [gene("orthA", 120, 200, home_gene_id="geneA"),
                   gene("GOI_match", 250, 330, home_gene_id=home_goi_name)]},
    ]


def _render(tracks, products=None):
    args = SimpleNamespace(plot_width=0, plot_height=0, scale_bar_len=10000,
                           pub_width=89, pub_palette="okabe_ito")
    return ps.render_synteny_html(tracks, {"geneA": "#0072B2", "geneB": "#009E73"},
                                  {"home": "#E64B35", "target1": "#0072B2"},
                                  products or {}, args, ["test"], 0, 0, 1)


def _home_goi_label(html):
    m = re.search(r'<text x="[\d.]+" y="([\d.]+)" transform="rotate\(-45 [^"]*"\s+'
                  r'class="gene-label goi track-item" data-track-idx="0"\s+'
                  r'font-size="(\d+)"[^>]*>([^<]*)</text>', html)
    assert m, "home GOI label not found"
    return float(m.group(1)), int(m.group(2)), m.group(3)


def test_long_home_goi_label_stays_below_the_subtitle():
    name = "GOI_4-coumarate--CoA_ligase_1-like"
    y, size, text = _home_goi_label(_render(_synteny_tracks(name)))
    assert y > 200                          # the header grew: the default puts it at 100
    top = y - ps._rotated_label_rise(text, size, True)
    assert top >= 66, f"label top at {top:.0f} px runs into the subtitle (baseline 60)"


def test_header_grows_only_as_far_as_the_label_needs():
    short = ps._top_margin_for_labels(_synteny_tracks("GOI_Melt")[0], {}, False, 40)
    long_ = ps._top_margin_for_labels(
        _synteny_tracks("GOI_4-coumarate--CoA_ligase_1-like")[0], {}, False, 40)
    assert short < 120 < long_
    # A track with no label at all asks for nothing: the default header stays.
    no_goi = {"is_home": True, "genes": [g for g in _synteny_tracks("GOI_Melt")[0]["genes"]
                                         if not g["name"].startswith("GOI_")]}
    assert ps._top_margin_for_labels(no_goi, {}, False, 40) == 0


def test_label_in_a_lower_lane_needs_less_header():
    track = _synteny_tracks("GOI_4-coumarate--CoA_ligase_1-like")[0]
    lane0 = ps._top_margin_for_labels(track, {}, False, 40)
    track["genes"][1]["_sub_track"] = 2
    assert ps._top_margin_for_labels(track, {}, False, 40) == lane0 - 80


def test_overlong_label_is_clipped_with_an_ellipsis():
    assert ps._clip_label("apyrase") == "apyrase"
    long_name = "x" * 80
    clipped = ps._clip_label(long_name)
    assert len(clipped) == ps.GOI_LABEL_MAX_CHARS and clipped.endswith("…")


def test_label_jobs_tag_the_best_goi_with_the_copy_count():
    track = _synteny_tracks("GOI_Melt")[1]
    track["genes"].append(dict(track["genes"][1], name="GOI_second", identity=40.0,
                               start=400, end=480, start_plot=400, end_plot=480))
    jobs = ps._track_label_jobs(track, {}, False)
    assert len(jobs) == 1 and jobs[0][2] is True
    # The second call has no confidence of at least medium: it is a weak call, not a copy.
    assert jobs[0][0]["name"] == "GOI_match" and jobs[0][1].endswith("+1 weak")
    assert "×" not in jobs[0][1]
    track["genes"][1]["confidence"], track["genes"][-1]["confidence"] = "HIGH", "MEDIUM"
    assert ps._track_label_jobs(track, {}, False)[0][1].endswith("×2")


def test_home_row_prints_its_scaffold_once():
    html = _render(_synteny_tracks("GOI_Melt"))
    assert ">NW_022170249.1</text>" in html
    assert "NW_022170249.1 • " not in html
    # A target row still shows genome and scaffold.
    assert "Drosophila_melanogaster • chr1" in html


# ------------------------------------------------------------ species names ----

def test_strain_qualifier_is_stripped_from_a_species_name():
    assert su.strip_species_qualifier(
        "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)") == "Saccharomyces cerevisiae"
    assert su.strip_species_qualifier("Escherichia coli (strain K12)") == "Escherichia coli"
    for kept in ("Apis mellifera", "Foo bar (nickname)", "Drosophila (Sophophora) melanogaster", ""):
        assert su.strip_species_qualifier(kept) == kept
    assert su.strip_species_qualifier(None) is None


def test_guide_tag_is_dropped_only_when_a_label_runs_through_it():
    # Label anchored 30 px left of the tag and 22 px below it: at 45 degrees it crosses the tag.
    box = (289, 66, 311, 76)
    assert ps._rotated_label_crosses(270, 100, "★ 4-coumarate--CoA ligase 1-like ×5", 9, True, box)
    # Anchored right of the tag, the label climbs away from it.
    assert not ps._rotated_label_crosses(320, 100, "★ Melt", 13, True, box)
    # Too short to reach the tag.
    assert not ps._rotated_label_crosses(200, 100, "★ M", 9, True, box)
