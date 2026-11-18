"""An answer of the genome-wide path in the figures (review round 3, 2026-11-17).

With --goi_genome_wide_path a genome without a counted call can keep one answer: a strong hit
without synteny, written MEDIUM with GOIClass=genome_wide_strong_hit, which the report never
counts. The figures read only the label and drew it as a counted call: the row's status was
'resolved', the anchor grid listed it among the counted calls outside the stretch drawn, the
ribbon plot counted it '+1 elsewhere', the copy rule took it as the genome's best call, and the
block choice drew its window as a call block. Now it is no counted call, no copy and never the
best call; it is drawn outlined with its own legend entry, its row status is strong_hit, the
grid's source table gives it its own state, and the block choice draws no block for it.
Synthetic ids and coordinates only."""
import contextlib
import io
import os
import re
import sys
from types import SimpleNamespace

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import plot_synteny as ps  # noqa: E402
import synvoy_label  # noqa: E402

NAMES = [f"gene-F{k}" for k in range(3)]
COLOURS = {n: ps.GENE_PALETTE[k] for k, n in enumerate(NAMES)}
OUTLINE = f'fill="#ffffff" stroke="{ps.GOI_BORDER}"'


def _flank(chrom, s, e, gid):
    return (f"{chrom}\tflanking_miniprot\tmRNA\t{s}\t{e}\t80\t+\t.\tID={gid}|G_fna_b0_{s}_flank_ann;Name={gid};"
            f"SynVoy_Parent={gid};Identity=80;Confidence=HIGH;SynVoyRole=flanking;"
            f"EvidenceType=flanking_miniprot")


def _goi(chrom, s, e, ident, conf, block=0, extra=""):
    return (f"{chrom}\texon_annotation\tmRNA\t{s}\t{e}\t{ident}\t+\t.\tID=GOI_T|G_fna_b{block}_l1_{s}_exon_ann;"
            f"Name=GOI_T;SynVoy_Parent=GOI_T;Identity={ident};Confidence={conf};SynVoyRole=goi;"
            f"EvidenceType=exon_annotation;ModelStatus=complete;QueryCoverage=0.95{extra}")


ANSWER = ";GOIClass=genome_wide_strong_hit;GenomeWide=1;Place=genome"


def _gff(tmp_path, answer_at=2_500_000, answer_block=5003):
    """Three flanking genes and a LOW model at the locus on chr1; the answer of the path at
    ``answer_at`` on the same scaffold, in its own window (block ``answer_block``)."""
    rows = [_flank("chr1", 10_000 * k, 10_000 * k + 900, n) for k, n in enumerate(NAMES)]
    rows.append(_goi("chr1", 15_000, 16_000, 40.0, "LOW"))
    rows.append(_goi("chr1", answer_at, answer_at + 1_000, 85.0, "MEDIUM", block=answer_block, extra=ANSWER))
    path = tmp_path / "G.fna.gff"
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return str(path)


def _home():
    home = [{"chrom": "h", "start": 10_000 * k, "end": 10_000 * k + 900, "start_plot": 10_000 * k,
             "end_plot": 10_000 * k + 900, "name": n, "strand": "+", "exon_coords": []}
            for k, n in enumerate(NAMES)]
    home.append({"chrom": "h", "start": 15_000, "end": 16_000, "start_plot": 15_000, "end_plot": 16_000,
                 "name": "GOI_T", "strand": "+", "exon_coords": []})
    return {"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home", "offset": 0,
            "breaks": [], "genes": home, "goi_status": "resolved"}


def _target(genes, outside=(), strong=()):
    for g in genes:
        g.setdefault("start_plot", g["start"])
        g.setdefault("end_plot", g["end"])
    t = {"label": "T0", "species": "Bombus x", "is_home": False, "genome_id": "G", "offset": 0, "breaks": [],
         "genes": genes, "goi_outside_view": list(outside), "goi_strong_outside_view": list(strong)}
    t["goi_status"] = ps._track_goi_status({"genes": genes + list(outside) + list(strong)})
    return t


def _view(gff, region_end=40_000):
    with contextlib.redirect_stdout(io.StringIO()):
        genes, _ = ps.select_target_genes(gff, "G", {"G": [("chr1", 0, region_end, "Reg1_G3_CHIGH_S0.9")]},
                                          1_000_000)
    return genes


def _model(conf="MEDIUM", cls="", start=1_000):
    return {"chrom": "c1", "start": start, "end": start + 900, "name": "GOI_T", "home_gene_id": "GOI_T",
            "role": "goi", "confidence": conf, "goi_class": cls, "identity": 80.0, "strand": "+",
            "exon_coords": [], "gff_id": f"GOI_T|m{start}"}


# ------------------------------------------------------------------------- the rule
def test_an_answer_is_no_counted_call_no_copy_and_never_the_best_call():
    assert ps.GENOME_WIDE_CLASS == synvoy_label.GENOME_WIDE_CLASS
    answer, medium, low = _model(cls=ps.GENOME_WIDE_CLASS), _model(start=5_000), _model("LOW", start=9_000)
    assert ps._is_genome_wide_answer(answer) and not ps._is_genome_wide_answer(medium)
    assert not ps._is_resolved_goi_target_gene(answer) and ps._is_resolved_goi_target_gene(medium)
    assert ps._goi_best_call([answer, low]) is low
    assert not ps._goi_is_copy(answer, answer)                 # not even as the only call
    assert ps._goi_copies_and_weak([answer, low]) == [(answer, False), (low, True)]
    assert ps._track_goi_status({"genes": [answer, low]}) == "strong_hit"
    assert ps._track_goi_status({"genes": [answer, low, medium]}) == "resolved"
    assert ps._track_goi_status({"genes": [low]}) == "ambiguous"


# ------------------------------------------------------------------------- the anchor grid
def test_the_grid_lists_an_answer_outside_the_view_outlined_and_counts_it_nowhere(tmp_path):
    """The reviewer's case: a LOW model at the locus, the answer 2 Mb away on the same scaffold.
    Before: 'counted calls outside the view' held the answer, the row was resolved, the source
    table gave it as an 'elsewhere ... model', the legend text called it a counted call and the
    ribbon label read '+1 weak +1 elsewhere'."""
    gff = _gff(tmp_path)
    genes = _view(gff)
    assert [(g["start"], g["confidence"]) for g in genes if ps._is_goi_target_gene(g)] == [(15_000, "LOW")]
    assert ps._counted_calls_outside_view(gff, genes) == []
    strong = ps._strong_hits_outside_view(gff, genes)
    assert [(g["start"], g["confidence"], g["outside_view"]) for g in strong] == [(2_500_000, "MEDIUM", True)]
    tracks = [_home(), _target(genes, strong=strong)]
    assert tracks[1]["goi_status"] == "strong_hit"
    svg, legend, source, info = ps.render_anchor_grid_print(tracks, COLOURS, {}, SimpleNamespace())
    rows = {r["start"]: r for r in source if r["row"] == 2 and r["column"] == "T"}
    assert (rows[2_500_000]["position"], rows[2_500_000]["drawn_as"], rows[2_500_000]["state"]) == (
        "elsewhere", "outline", ps.STRONG_HIT_NOTE)
    assert (rows[15_000]["drawn_as"], rows[15_000]["state"]) == ("model", "placed")   # the best call
    assert info["goi_strong_hits"] == 1 and info["goi_calls_outside_view"] == 0
    assert "counted call (high or medium confidence) is listed" not in legend
    assert "An outlined arrow in the T column is a strong hit without synteny" in legend
    assert OUTLINE in svg
    # the full grid carries the same glyph with its tooltip and legend entry
    html = ps.render_anchor_grid(tracks, COLOURS, {}, {}, SimpleNamespace())
    assert OUTLINE in html and ps.STRONG_HIT_NOTE in html
    # the ribbon plot's label of the row counts neither a weak call nor a call elsewhere
    assert ps._track_label_jobs(tracks[1], {}) == []


def test_an_answer_in_the_view_is_drawn_outlined_in_the_grid(tmp_path):
    """An answer next to the neighbourhood (inside a block's window, outside its GOI zone) is in
    the view: the grid draws it outlined at the locus; the LOW model stays the best call."""
    gff = _gff(tmp_path, answer_at=30_000)
    genes = _view(gff)
    assert ps._strong_hits_outside_view(gff, genes) == []
    tracks = [_home(), _target(genes)]
    svg, legend, source, info = ps.render_anchor_grid_print(tracks, COLOURS, {}, SimpleNamespace())
    rows = {r["start"]: r for r in source if r["row"] == 2 and r["column"] == "T"}
    assert (rows[30_000]["position"], rows[30_000]["drawn_as"]) == ("locus", "outline")
    assert rows[15_000]["drawn_as"] == "model"
    assert info["goi_strong_hits"] == 1 and OUTLINE in svg


# ------------------------------------------------------------------------- the ribbon plot
def test_the_ribbon_plot_draws_an_answer_outlined_with_its_status_and_no_goi_ribbon(tmp_path):
    gff = _gff(tmp_path, answer_at=30_000)
    genes = [g for g in _view(gff) if not (ps._is_goi_target_gene(g) and g["confidence"] == "LOW")]
    tracks = [_home(), _target(genes)]
    assert tracks[1]["goi_status"] == "strong_hit"
    args = SimpleNamespace(plot_width=0, plot_height=0, scale_bar_len=10000)
    html = ps.render_synteny_html(tracks, COLOURS, {}, args, [], 0, 0, 0)
    group = re.search(r'<g class="gene-group" data-homology="GOI_T" data-track="1"[^>]*>', html).group(0)
    assert 'data-fill="#ffffff"' in group and 'data-goi-resolved' not in group
    assert ps.STRONG_HIT_NOTE in html                    # tooltip tag and legend entry
    assert "Strong hit without synteny, not counted" in html
    assert 'data-homology="GOI" data-upper-track="0"' not in html     # no GOI ribbon to it
    # a counted call in the same row is the row's GOI, with its ribbon
    answer = next(g for g in genes if ps._is_genome_wide_answer(g))
    medium = dict(answer, start=20_000, end=20_900, start_plot=20_000, end_plot=20_900,
                  goi_class="", gff_id="GOI_T|G_fna_b0_l1_20000_exon_ann")
    tracks = [_home(), _target(genes + [medium])]
    assert tracks[1]["goi_status"] == "resolved"
    html = ps.render_synteny_html(tracks, COLOURS, {}, args, [], 0, 0, 1)
    assert len(re.findall(r'data-homology="GOI" data-upper-track="0"', html)) == 1


# ------------------------------------------------------------------------- the gene-position map
def test_the_gene_position_map_draws_an_answer_alone_outlined(tmp_path):
    gff = _gff(tmp_path, answer_at=30_000)
    genes = [g for g in _view(gff) if not (ps._is_goi_target_gene(g) and g["confidence"] == "LOW")]
    html = ps.render_anchor_grid_positional([_home(), _target(genes)], COLOURS, {}, {}, SimpleNamespace())
    assert f'fill="#ffffff" stroke="{ps.GOI_BORDER}" stroke-width="1.3"' in html
    assert ps.STRONG_HIT_NOTE in html


# ------------------------------------------------------------------------- the block choice
def test_the_block_choice_draws_no_block_for_an_answer(tmp_path):
    """The answer's window (its block number in the model id) was drawn as a call block and took
    one of the row's block slots; its genome then had a 'counted' block."""
    gff = _gff(tmp_path)
    table = {0: {"esyn": 0.001, "significant": "1", "chrom": "chr1", "hood_start": 0, "hood_end": 40_000}}
    choice = ps.choose_blocks(ps.choice_models(gff), ps.genome_tag(gff), table)
    assert choice.drawn == {0: "predicted_no_call"} and choice.predicted == {0}
