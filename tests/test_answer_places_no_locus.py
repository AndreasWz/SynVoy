"""An answer of the genome-wide path places no locus (review round 4 of the bug-fix series, 2026-11-18).

An answer (GOIClass=genome_wide_strong_hit, MEDIUM, never counted) lies outside every searched GOI
zone: it has no synteny. Three steps still took it as the genome's GOI:

- the region step (cluster_grs.py) took it as a GOI interval: it got the GOI-overlap bonus, went
  first, could be injected as a goi_anchor region, and its score raised the floor of the genome's
  other regions, so two-gene neighbourhoods the search had placed were not written;
- the figures' view choice (plot_synteny.select_target_genes) took it as the row's GOI: a
  GOI-centred view of the answer alone when the regions held no other GOI model, its scaffold kept
  as a GOI scaffold;
- an answer that was a row's only GOI call in view set the row's place in the anchor grid
  (_grid_target_map): the location column read the answer's scaffold, the reach in which a
  flanking gene's other models count as copies widened to it, and the gene-position map kept only
  the points on its scaffold, so the flanking genes were not drawn.

Now none of them reads it: the regions and the view are those of a run without the path, the
answer is drawn where it lies inside that view, and the anchor grid lists it elsewhere. The grid's
legend sentence also names the report count of answers inside a block's searched window.
Synthetic ids and coordinates only."""
import contextlib
import csv
import html
import io
import os
import re
import sys
from types import SimpleNamespace

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import cluster_grs as cg  # noqa: E402
import plot_synteny as ps  # noqa: E402
import synvoy_label  # noqa: E402

FLANKS = [f"F{i}" for i in range(1, 6)]
COLOURS = {n: ps.GENE_PALETTE[i] for i, n in enumerate(FLANKS)}
OUTLINED_DIAMOND = f'fill="#ffffff" stroke="{ps.GOI_BORDER}" stroke-width="1.3"'
ANSWER_ATTRS = ";GOIClass=genome_wide_strong_hit;GenomeWide=1;Place=genome"


def _gene(name, chrom, start, conf="HIGH", ident=80.0, role="flanking", cls="", gid=None, length=600):
    return {"chrom": chrom, "start": start, "end": start + length - 1, "start_plot": start,
            "end_plot": start + length - 1, "name": name, "home_gene_id": name, "role": role,
            "confidence": conf, "identity": ident, "goi_class": cls, "strand": "+", "exon_coords": [],
            "gff_id": gid or f"{name}|G_fna_b0_{start}_ann", "query_coverage": 0.95,
            "evidence_type": "exon_annotation" if role == "goi" else "flanking_miniprot"}


def _answer(chrom, start, **kw):
    return _gene("GOI_X", chrom, start, conf="MEDIUM", ident=85.0, role="goi", cls=ps.GENOME_WIDE_CLASS,
                 gid=f"GOI_X|G_fna_b5003_l1_{start}_exon_ann", **kw)


def _home():
    genes = [{"chrom": "chrH", "start": i * 10_000, "end": i * 10_000 + 2_000, "start_plot": i * 10_000,
              "end_plot": i * 10_000 + 2_000, "name": n, "strand": "+", "exon_coords": []}
             for i, n in enumerate(FLANKS, 1)]
    genes.append({"chrom": "chrH", "start": 35_000, "end": 36_000, "start_plot": 35_000, "end_plot": 36_000,
                  "name": "GOI_X", "strand": "+", "exon_coords": []})
    return {"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home", "offset": 0,
            "breaks": [], "genes": genes, "goi_status": "resolved"}


def _flanking(chrom="chrA"):
    """Five flanking genes at chrA:0.11-0.15 Mb, in home order."""
    return [_gene(n, chrom, 110_001 + 10_000 * k) for k, n in enumerate(FLANKS)]


def _target(genes, outside=(), strong=()):
    t = {"label": "T0", "species": "Bombus x", "is_home": False, "genome_id": "G", "offset": 0, "breaks": [],
         "genes": genes, "goi_outside_view": list(outside), "goi_strong_outside_view": list(strong)}
    t["goi_status"] = ps._track_goi_status({"genes": genes + list(outside) + list(strong)})
    return t


def _figures(tracks):
    """(full grid HTML, gene-position map HTML, print-grid source rows of the target, info)."""
    full = html.unescape(ps.render_anchor_grid(tracks, COLOURS, {}, {}, SimpleNamespace()))
    pos = ps.render_anchor_grid_positional(tracks, COLOURS, {}, {}, SimpleNamespace())
    _svg, legend, source, info = ps.render_anchor_grid_print(tracks, COLOURS, {}, SimpleNamespace())
    return full, pos, [r for r in source if r["row"] == 2], legend, info


def _labels(text):
    return sorted(set(re.findall(r">([^<>]*Mb[^<>]*)<", text)))


def _text(svg):
    """The text of an SVG or HTML page, its elements joined by single spaces (legend lines wrap)."""
    return re.sub(r"\s+", " ", " ".join(re.findall(r">([^<>]+)<", html.unescape(svg))))


def _target_dot_titles(pos):
    """The tooltips of the target row's points in the gene-position map (identity 80 or 85 %)."""
    return [t for t in re.findall(r"<title>([^<]*)</title>", html.unescape(pos))
            if " — identity " in t and ("identity 80%" in t or "identity 85%" in t)]


# ------------------------------------------------------------- the row's place in the figures ----

def test_an_answer_alone_on_another_scaffold_does_not_place_the_row():
    """Five flanking genes on chrA and the row's only GOI call in view an answer on chrB. Before: the
    location column read 'chrB: 2.50 Mb (+1)', and the gene-position map drew the answer alone
    under that label, without the five flanking genes."""
    tracks = [_home(), _target(_flanking() + [_answer("chrB", 2_500_001)])]
    full, pos, source, _legend, _info = _figures(tracks)
    assert "chrA: 0.11-0.15 Mb" in _labels(full) and not any("chrB" in x for x in _labels(full))
    # the gene-position map: the row on its flanking genes' scaffold, the answer on another one
    # counted in (+1) and named in the label's tooltip, not drawn
    assert "chrA: 0.11-0.15 Mb (+1)" in _labels(pos)
    titles = _target_dot_titles(pos)
    assert sorted(t.split(" — ")[0] for t in titles) == FLANKS
    assert OUTLINED_DIAMOND not in pos
    tip = next(t for t in re.findall(r"<title>([^<]*)</title>", html.unescape(pos)) if "chrB:2,500,001" in t)
    assert ps.STRONG_HIT_NOTE in tip and "not drawn" in tip
    # the print grid lists the answer after //, the flanking genes at the locus
    goi = [r for r in source if r["column"] == "X"]
    assert [(r["position"], r["drawn_as"], r["state"]) for r in goi] == [("elsewhere", "outline", ps.STRONG_HIT_NOTE)]
    assert sorted(r["column"] for r in source if r["state"] == "placed" and r["column"] != "X") == FLANKS


def test_an_answer_far_on_the_same_scaffold_does_not_stretch_the_row():
    """The same answer on chrA, 2.4 Mb from the flanking genes. Before: both labels read
    'chrA: 0.11-2.50 Mb' and the gene-position map's row was scaled to that stretch."""
    tracks = [_home(), _target(_flanking() + [_answer("chrA", 2_500_001)])]
    full, pos, source, _legend, _info = _figures(tracks)
    assert "chrA: 0.11-0.15 Mb" in _labels(full)
    assert "chrA: 0.11-0.15 Mb" in _labels(pos)
    assert not any("2.50" in x for x in _labels(full) + _labels(pos))
    assert OUTLINED_DIAMOND not in pos and len(_target_dot_titles(pos)) == 5
    assert any(ps.STRONG_HIT_NOTE in t and "chrA:2,500,001" in t
               for t in re.findall(r"<title>([^<]*)</title>", html.unescape(pos)))
    goi = [r for r in source if r["column"] == "X"]
    assert [(r["position"], r["drawn_as"]) for r in goi] == [("elsewhere", "outline")]


def test_an_answer_at_the_locus_is_drawn_where_it_lies():
    """An answer 10 kb beyond the last flanking gene is at the locus (within the neighbourhood's
    length): the gene-position map draws it outlined, and the row's label stays the stretch of
    its flanking genes."""
    tracks = [_home(), _target(_flanking() + [_answer("chrA", 160_001)])]
    _full, pos, source, _legend, _info = _figures(tracks)
    assert OUTLINED_DIAMOND in pos and len(_target_dot_titles(pos)) == 6
    assert "chrA: 0.11-0.15 Mb" in _labels(pos)
    assert "the row's only GOI call, drawn where it lies in the neighbourhood" in _text(pos)
    goi = [r for r in source if r["column"] == "X"]
    assert [(r["position"], r["drawn_as"]) for r in goi] == [("locus", "outline")]


def test_an_answer_does_not_widen_the_reach_of_a_flanking_copy():
    """A weaker model of F1 at 2,495 kb, next to the answer at 2,500 kb: with the answer in view it
    counted as a second copy at the locus (x2), as beside a counted call; now as with no call."""
    def row(goi):
        genes = [_gene(n, "chr1", 10_000 * k) for k, n in enumerate(FLANKS[:3])]
        genes.append(_gene("F1", "chr1", 2_495_000, conf="MEDIUM", ident=45.0, gid="F1|G_fna_b9_x_ann"))
        return _target(genes + goi)
    keys = set(FLANKS)
    with_answer = ps._grid_target_map(row([_answer("chr1", 2_500_000)]), keys, "__GOI__")
    medium = _gene("GOI_X", "chr1", 2_500_000, conf="MEDIUM", ident=85.0, role="goi", cls="probable_goi")
    with_call = ps._grid_target_map(row([medium]), keys, "__GOI__")
    without = ps._grid_target_map(row([]), keys, "__GOI__")
    assert (with_answer["F1"]["n"], with_call["F1"]["n"], without["F1"]["n"]) == (1, 2, 1)
    assert ps._row_place(with_answer, "__GOI__")[:3] == ps._row_place(without, "__GOI__")[:3]


def test_the_gene_position_map_names_an_answer_outside_the_view():
    """An answer outside the view chosen without it is no point of the gene-position map; the
    tooltip of the row's location names it, as the anchor grid lists it."""
    outside = dict(_answer("chrB", 2_500_001), outside_view=True)
    tracks = [_home(), _target(_flanking(), strong=[outside])]
    _full, pos, source, _legend, _info = _figures(tracks)
    assert "chrA: 0.11-0.15 Mb" in _labels(pos) and len(_target_dot_titles(pos)) == 5
    (tip,) = [t for t in re.findall(r"<title>([^<]*)</title>", html.unescape(pos)) if "chrB:2,500,001" in t]
    assert tip == (f"X chrB:2,500,001-2,500,600, {ps.STRONG_HIT_NOTE}: not drawn, it lies outside this "
                   "neighbourhood (anchor grid)")
    assert [(r["position"], r["drawn_as"]) for r in source if r["column"] == "X"] == [("elsewhere", "outline")]


def test_a_row_with_only_an_answer_has_no_place():
    """Nothing of the neighbourhood placed, only an answer: the location column says so and the
    gene-position map names the answer in a tooltip, instead of placing the row on it."""
    tracks = [_home(), _target([_answer("chrB", 2_500_001)])]
    full, pos, source, _legend, _info = _figures(tracks)
    assert "no flanking gene placed" in full and not any("chrB" in x for x in _labels(full))
    assert "no orthologs mapped" in pos and OUTLINED_DIAMOND not in pos
    assert any(ps.STRONG_HIT_NOTE in t for t in re.findall(r"<title>([^<]*)</title>", html.unescape(pos)))
    goi = [r for r in source if r["column"] == "X"]
    assert [(r["position"], r["drawn_as"]) for r in goi] == [("elsewhere", "outline")]


# ------------------------------------------------------------------------ the legend sentence ----

def test_the_legend_sentence_names_the_report_count_of_answers_near_a_block():
    """The report counts an answer inside a block's searched window (NearBlock) apart, in
    summary.near_block_strong_hits; the sentence pointed to genome_wide_strong_hits only."""
    near = dict(_answer("chrA", 160_001), near_block="1", outside_view=True)
    far = dict(_answer("chrC", 900_001), outside_view=True)

    def legend(*answers):
        tracks = [_home()] + [dict(_target(_flanking(), strong=[a]), genome_id=f"G{i}", label=f"T{i}")
                              for i, a in enumerate(answers)]
        return ps.render_anchor_grid_print(tracks, COLOURS, {}, SimpleNamespace())[1]
    assert ("(1 in this figure; synvoy_report.json, summary.genome_wide_strong_hits)." in legend(far))
    assert ("(1 in this figure, inside a block's searched window, outside its GOI zone; synvoy_report.json, "
            "summary.near_block_strong_hits)." in legend(near))
    assert ("(2 in this figure, 1 of them inside a block's searched window, outside its GOI zone; "
            "synvoy_report.json, summary.genome_wide_strong_hits and summary.near_block_strong_hits)."
            in legend(far, near))


def test_the_figures_read_near_block_from_the_gff(tmp_path):
    gff = tmp_path / "G.fna.gff"
    gff.write_text("##gff-version 3\nchrA\tminiprot\tmRNA\t1000\t2000\t85\t+\t.\tID=GOI_X|G_fna_b5003_l1_exon_ann;"
                   "Name=GOI_X;SynVoy_Parent=GOI_X;Identity=85;Confidence=MEDIUM;SynVoyRole=goi;"
                   "EvidenceType=exon_annotation;ModelStatus=complete" + ANSWER_ATTRS + ";NearBlock=7\n")
    with contextlib.redirect_stdout(io.StringIO()):
        (g,) = ps.parse_target_gff(str(gff), cap=False)
    assert g["near_block"] == "7" and ps._is_genome_wide_answer(g)


# ------------------------------------------------------------------- the figures' view choice ----

def _flank_row(chrom, s, gid):
    return (f"{chrom}\tflanking_miniprot\tmRNA\t{s}\t{s + 599}\t80\t+\t.\tID={gid}|G_fna_b0_{s}_flank_ann;"
            f"Name={gid};SynVoy_Parent={gid};Identity=80;Confidence=HIGH;SynVoyRole=flanking;"
            f"EvidenceType=flanking_miniprot")


def _goi_row(chrom, s, conf, extra="", block=0, ident=85):
    return (f"{chrom}\texon_annotation\tmRNA\t{s}\t{s + 599}\t{ident}\t+\t.\tID=GOI_X|G_fna_b{block}_l1_{s}_exon_ann;"
            f"Name=GOI_X;SynVoy_Parent=GOI_X;Identity={ident};Confidence={conf};SynVoyRole=goi;"
            f"EvidenceType=exon_annotation;ModelStatus=complete;QueryCoverage=0.95{extra}")


def _gff(path, rows):
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return str(path)


def _view(gff, regions):
    with contextlib.redirect_stdout(io.StringIO()) as out:
        genes, _ = ps.select_target_genes(gff, "G", {"G": regions}, 1_000_000)
    return genes, out.getvalue()


def test_the_view_is_chosen_without_an_answer(tmp_path):
    """The regions hold the five flanking genes and no GOI model; the answer lies on chrB. Before:
    'candidate regions missed GOI; using GOI-centered fallback', and the row showed the answer
    alone. Now the row shows its neighbourhood, as without the path, and the grid lists the
    answer."""
    flanks = [_flank_row("chrA", 110_001 + 10_000 * k, n) for k, n in enumerate(FLANKS)]
    answer = _goi_row("chrB", 2_500_001, "MEDIUM", ANSWER_ATTRS, block=5003)
    regions = [("chrA", 110_000, 150_600, "Reg1_FLANKonly5_CMEDIUM_S0.40")]
    on = _gff(tmp_path / "on.fna.gff", flanks + [answer])
    off = _gff(tmp_path / "off.fna.gff", flanks)
    view, log = _view(on, regions)
    assert "GOI-centered fallback" not in log
    assert [(g["chrom"], g["start"]) for g in view] == [(g["chrom"], g["start"]) for g in _view(off, regions)[0]]
    assert [g["start"] for g in ps._strong_hits_outside_view(on, view)] == [2_500_001]


def test_an_answer_takes_no_part_in_the_choice_and_is_drawn_inside_the_view(tmp_path):
    """A LOW model at the locus, a flanking gene's paralog copy and the answer on chrB, and a
    region around each: the GOI regions, the GOI scaffolds kept and the GOI cap are those of a run
    without the path; the answer is in the view only where it lies between the view's models."""
    flanks = [_flank_row("chrA", 110_001 + 10_000 * k, n) for k, n in enumerate(FLANKS)]
    low = _goi_row("chrA", 135_001, "LOW", ident=40)
    copy_b = _flank_row("chrB", 2_400_001, "F1")
    regions = [("chrA", 110_000, 150_600, "Reg1_G5_CHIGH_S0.90"),
               ("chrB", 2_380_000, 2_600_000, "Reg2_G1_CLOW_S0.10")]
    for at, inside in ((2_500_001, False), (125_001, True)):
        answer = _goi_row("chrB" if at > 1_000_000 else "chrA", at, "MEDIUM", ANSWER_ATTRS, block=5003)
        on = _gff(tmp_path / f"on{at}.fna.gff", flanks + [low, copy_b, answer])
        off = _gff(tmp_path / f"off{at}.fna.gff", flanks + [low, copy_b])
        view, _ = _view(on, regions)
        rest = [(g["chrom"], g["start"]) for g in view if g.get("goi_class") != synvoy_label.GENOME_WIDE_CLASS]
        assert rest == [(g["chrom"], g["start"]) for g in _view(off, regions)[0]]
        assert any(ps._is_genome_wide_answer(g) for g in view) is inside
        assert len(ps._strong_hits_outside_view(on, view)) == (0 if inside else 1)


def test_the_track_is_not_centred_on_an_answer():
    flanking = [dict(g, start_plot=g["start"], end_plot=g["end"]) for g in _flanking()]
    centre = ps.get_anchor_center(flanking)
    assert ps.get_anchor_center(flanking + [_answer("chrA", 160_001)]) == centre
    call = _gene("GOI_X", "chrA", 131_001, conf="MEDIUM", role="goi", cls="probable_goi")
    assert ps.get_anchor_center(flanking + [call, _answer("chrA", 160_001)]) == (131_001 + 131_600) / 2


# ------------------------------------------------------------------------------ the region step ----

def _home_bed(path, n=30):
    path.write_text("".join(f"chrH\t{i * 10_000}\t{i * 10_000 + 2_000}\tF{i}\t0\t+\n" for i in range(1, n + 1)))
    return path


def _hit(q, chrom, s, ident=80.0, ev=1e-30):
    return f"{q}\t{chrom}\t{ident}\t300\t0\t0\t1\t100\t{s}\t{s + 600}\t{ev}\t200\n"


def _regions(tmp_path, tag, gff_rows, hits):
    """cluster_grs.main() on one genome: (regions.bed lines, scores rows, stderr)."""
    d = tmp_path / tag
    d.mkdir()
    (d / "g.m8").write_text("".join(hits))
    gff = d / "G.fna.gff"
    gff.write_text("##gff-version 3\n" + "\n".join(gff_rows) + "\n")
    argv = ["cluster_grs.py", "--hits", str(d / "g.m8"), "--synteny_bed", str(_home_bed(d / "home.bed")),
            "--target_gff", str(gff), "--genome_name", "G.fna", "--output", str(d / "g.regions.bed"),
            "--scores_output", str(d / "g.scores.tsv"), "--cluster_distance", "50000",
            "--goi_overlap_bonus", "0.15"]
    old, sys.argv = sys.argv, argv
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            cg.main()
    finally:
        sys.argv = old
    with open(d / "g.scores.tsv", newline="") as fh:
        scores = list(csv.DictReader(fh, delimiter="\t"))
    return (d / "g.regions.bed").read_text().splitlines(), scores, err.getvalue()


def _neighbourhood_hits():
    """A four-gene neighbourhood (F1-F4) and three two-gene ones on chrA, of 30 home genes."""
    hits = [_hit(f"F{i}", "chrA", 100_000 + 10_000 * i) for i in range(1, 5)]
    for k, (a, b) in enumerate(((10, 11), (12, 13), (14, 15))):
        base = 1_000_000 * (k + 1)
        hits += [_hit(f"F{a}", "chrA", base), _hit(f"F{b}", "chrA", base + 10_000)]
    return hits


ANSWER_ROW = ("chrB\texon_annotation\tmRNA\t5000001\t5000600\t85\t+\t.\tID=GOI_X|G_fna_b5003_l1_exon_ann;"
              "Name=GOI_X;SynVoy_Parent=GOI_X;Identity=85;Confidence=MEDIUM;SynVoyRole=goi;"
              "EvidenceType=exon_annotation;ModelStatus=complete;QueryCoverage=0.95" + ANSWER_ATTRS)
LOW_ROW = ("chrA\tfallback_hits\tmRNA\t3500001\t3500200\t30\t+\t.\tID=GOI_X|G_fna_b0_r1_fallback;Name=GOI_X;"
           "SynVoy_Parent=GOI_X;Identity=30;Confidence=LOW;SynVoyRole=goi;EvidenceType=fallback_hit_span")


def test_the_region_step_leaves_an_answer_out(tmp_path):
    """The fly_def Anopheles case in small: a LOW GOI model, the answer on chrB with the GOI hits
    there, a four-gene and three two-gene neighbourhoods. As a GOI interval the answer's hits
    formed a region of their own (score 0.15, the bonus alone) that raised the floor to
    0.3 x 0.15 = 0.045, so the two-gene neighbourhoods (0.04) were not written. Now the regions
    are those without the answer."""
    hits = _neighbourhood_hits() + [_hit("GOI_X", "chrB", 5_000_001, ident=85.0, ev=1e-80),
                                    _hit("GOI_X", "chrB", 5_000_301, ident=85.0, ev=1e-60)]
    on, scores_on, err = _regions(tmp_path, "on", [LOW_ROW, ANSWER_ROW], hits)
    off, scores_off, _ = _regions(tmp_path, "off", [LOW_ROW], hits)
    assert on == off
    assert not any(line.startswith("chrB") for line in on) and len(on) == 4
    assert [r["goi_overlap"] for r in scores_on] == ["False"] * 4
    assert "answer(s) of the genome-wide path (genome_wide_strong_hit) out of the GOI intervals" in err
    gff = str(tmp_path / "on" / "G.fna.gff")
    assert cg.load_goi_intervals_from_gff(gff) == []
    assert cg.load_goi_block_anchors_from_gff(gff, genome_tag="G_fna") == []
    # what a region holds still names it (goi_model_class), and it is no counted model
    (model,) = [m for m in cg.load_goi_model_labels_from_gff(gff) if m["goi_class"]]
    assert model["goi_class"] == synvoy_label.GENOME_WIDE_CLASS and not model["counted"]


def test_raw_goi_hits_stand_in_as_without_the_path(tmp_path):
    """No GOI model but the answer: raw GOI hits stand in for the models, as in a run without the
    path, and the answer's search window does not widen the region they give."""
    hits = _neighbourhood_hits() + [_hit("GOI_X", "chrB", 5_000_001, ident=85.0, ev=1e-80)]
    on, _, err = _regions(tmp_path, "on", [ANSWER_ROW], hits)
    off, _, _ = _regions(tmp_path, "off", [], hits)
    assert on == off
    assert "GOI interval(s) from hit file" in err
