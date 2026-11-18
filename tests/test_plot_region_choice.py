"""Which stretch of a target genome a row of the figures shows, and that no counted call is lost.

Without block tables a target row shows the candidate regions of the region scorer
(plot_synteny.select_target_genes). Three rules of that choice left out what the figures
should show:

- the regions with a GOI model replaced all others whenever they held half the genes, so a
  strong neighbourhood without a GOI model (cluster_grs names it Reg<k>_FLANKonly<N>_..., a
  homeolog that lost the gene, or the neighbourhood of a GOI too divergent to model) vanished,
  and no line said so;
- with the models on 20 sequences or fewer, only the GOI sequences and those with three or more
  flanking genes were kept, so a fragmented neighbourhood lost every contig with one or two
  neighbours;
- the GOI cap that steers the view (max_goi_per_genome) and the drop of a GOI on a sequence
  without flanking genes took counted calls (HIGH or MEDIUM) out of every figure.

Now the strong neighbourhoods are drawn, every region left out is named on the console, a
sequence that holds a neighbour placed nowhere else is kept, and the counted calls outside the
view are listed in the anchor grid and counted in the ribbon plot ('+N elsewhere'). In the grid
they are placed as the calls in view are: at the locus when they lie within reach of it (a copy of
a tandem array beyond the cap), after '//' on another scaffold or farther away,
and a scaffold drawn reversed mirrors them with it.
"""
import contextlib
import csv
import io
import os
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))

import plot_synteny as ps  # noqa: E402
import synvoy_grid as sg  # noqa: E402


def _flank(chrom, s, e, gid, ident=80.0, conf="HIGH"):
    return (f"{chrom}\tflanking_miniprot\tmRNA\t{s}\t{e}\t{ident}\t+\t.\t"
            f"ID={gid}|G_fna_{chrom}_{s}_flank_ann;Name={gid};SynVoy_Parent={gid};Identity={ident};"
            f"Confidence={conf};SynVoyRole=flanking;EvidenceType=flanking_miniprot")


def _goi(chrom, s, e, ident=70.0, conf="HIGH"):
    return (f"{chrom}\texon_annotation\tmRNA\t{s}\t{e}\t{ident}\t+\t.\t"
            f"ID=GOI_X|G_fna_{chrom}_{s}_exon_ann;Name=GOI_X;SynVoy_Parent=GOI_X;"
            f"Identity={ident};Confidence={conf};SynVoyRole=goi;EvidenceType=exon_annotation;"
            f"ModelStatus=complete;QueryCoverage=0.95")


def _gff(tmp_path, rows, name="G.fna.gff"):
    path = tmp_path / name
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return str(path)


def _select(tmp_path, rows, regions, cap=10, monkeypatch=None):
    """select_target_genes on a synthetic row: ({chrom: models drawn}, console)."""
    if monkeypatch is not None:
        monkeypatch.setattr(ps, "_MAX_GOI_PER_GENOME", cap)
    gff = _gff(tmp_path, rows)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        genes, _regions = ps.select_target_genes(gff, "G", {"G": regions}, 1_000_000)
    drawn = {}
    for g in genes:
        drawn[g["chrom"]] = drawn.get(g["chrom"], 0) + 1
    return drawn, out.getvalue(), gff, genes


def _homeolog(region2_name):
    rows = [_goi("chr1", 50_000, 52_000)]
    rows += [_flank("chr1", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(1, 7)]
    rows += [_flank("chr2", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(1, 6)]
    regions = [("chr1", 0, 100_000, "Apis_x|Reg1_G7_CHIGH_S0.80"), ("chr2", 0, 100_000, region2_name)]
    return rows, regions


# ---------------------------------------------------------- regions without a GOI model ----

def test_a_strong_neighbourhood_without_a_goi_model_is_drawn(tmp_path, monkeypatch):
    rows, regions = _homeolog("Apis_x|Reg2_FLANKonly5_CMEDIUM_S0.40")
    drawn, log, _, _ = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 7, "chr2": 5}
    assert "left out" not in log


def test_a_region_without_a_goi_model_is_named_when_it_is_left_out(tmp_path, monkeypatch):
    rows, regions = _homeolog("Apis_x|Reg2_G5_CLOW_S0.12")
    regions.append(("chr4", 0, 5000, "Apis_x|Reg3_DISP_CMEDIUM_S0.45"))     # holds no model
    drawn, log, _, _ = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 7}
    assert ("[plot] G: kept the 1 region(s) with a GOI model; left out "
            "Reg2_G5_CLOW_S0.12 chr2:1-100000 (5 models)\n") in log


def test_candidate_regions_keep_their_name(tmp_path):
    named = tmp_path / "G.fna.regions.bed"
    named.write_text("chr1\t0\t100000\tApis_x|Reg1_G7_CHIGH_S0.80\t0.80\t+\n"
                     "chr2\t500\t900\tReg2_FLANKonly6_CMEDIUM_S0.30\t0.30\t-\n")
    bare = tmp_path / "H.fna.regions.bed"
    bare.write_text("chr9\t10\t20\n")
    regions = ps.parse_candidate_regions([str(named), str(bare)])
    assert regions["G"] == [("chr1", 0, 100000, "Apis_x|Reg1_G7_CHIGH_S0.80"),
                            ("chr2", 500, 900, "Reg2_FLANKonly6_CMEDIUM_S0.30")]
    assert regions["H"] == [("chr9", 10, 20, "")]
    gene = {"chrom": "chr2", "start": 600, "end": 700}
    assert ps.filter_genes_to_candidate_regions([gene], regions["G"]) == [gene]
    assert ps._region_overlaps_gene(regions["G"][1], gene)
    assert [ps._is_strong_no_goi_region(r) for r in regions["G"] + regions["H"]] == [
        False, True, False]


@pytest.mark.parametrize("name,strong", [
    ("Reg3_FLANKonly7_CMEDIUM_S0.42", True), ("Apis_mellifera|Reg1_FLANKonly12_CHIGH_S1.20", True),
    ("Reg2_G5_CLOW_S0.12", False), ("Reg1_GOI_anchor_CHIGH_S0.90", False),
    ("Reg4_DISP_CMEDIUM_S0.30", False), ("FLANKonly", False), ("", False)])
def test_the_strong_neighbourhood_name_is_the_region_scorers(name, strong):
    """The figures recognise the class by the name cluster_grs.py writes into the BED."""
    assert ps._is_strong_no_goi_region(("c", 0, 1, name)) is strong
    src = open(os.path.join(REPO, "bin", "cluster_grs.py")).read()
    assert 'f"Reg{i+1}_FLANKonly{best[\'high_flanking_count\']}_"' in src
    assert 'region_class = "goi_missing_but_strong_synteny"' in src
    i, count, conf, score = 2, 7, "MEDIUM", 0.4231
    built = f"Reg{i+1}_FLANKonly{count}_" f"C{conf}_S{score:.2f}"
    assert ps._is_strong_no_goi_region(("c", 0, 1, built))


# ----------------------------------------------- neighbours on sequences of their own ----

def test_a_fragmented_neighbourhood_keeps_every_contig_with_a_neighbour(tmp_path, monkeypatch):
    rows = [_goi("ctg1", 5000, 6000), _flank("ctg1", 1000, 1900, "gene-F3"),
            _flank("ctg1", 8000, 8900, "gene-F4"),
            _flank("ctg2", 1000, 1900, "gene-F1"), _flank("ctg2", 4000, 4900, "gene-F2"),
            _flank("ctg3", 1000, 1900, "gene-F5"),
            _flank("ctg4", 1000, 1900, "gene-F6"), _flank("ctg4", 4000, 4900, "gene-F7")]
    regions = [(c, 0, 10_000, f"Reg{k}_G2_CLOW_S0.10") for k, c in
               enumerate(("ctg1", "ctg2", "ctg3", "ctg4"), 1)]
    drawn, log, _, _ = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert drawn == {"ctg1": 3, "ctg2": 2, "ctg3": 1, "ctg4": 2}
    assert "dropped" not in log


def test_a_sequence_that_only_repeats_placed_neighbours_is_still_dropped(tmp_path, monkeypatch):
    """chr9 holds a second hit of F2, which is placed beside the GOI; chr7 holds F6, which is
    placed nowhere else. No candidate regions (a genome without any): the sequence rule alone."""
    rows = [_goi("chr1", 50_000, 52_000)]
    rows += [_flank("chr1", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(1, 6)]
    rows += [_flank("chr9", 30_000, 30_900, "gene-F2"), _flank("chr7", 70_000, 70_900, "gene-F6")]
    drawn, log, _, _ = _select(tmp_path, rows, [], monkeypatch=monkeypatch)
    assert drawn == {"chr1": 6, "chr7": 1}
    assert "[plot] G: keeping 2 informative chromosomes (7 genes), dropped chr9" in log


def test_more_than_twenty_sequences_are_drawn_unfiltered_as_before(tmp_path, monkeypatch):
    rows = [_goi("ctg0", 5000, 6000), _flank("ctg0", 1000, 1900, "gene-F1")]
    rows += [_flank(f"ctg{k}", 1000, 1900, "gene-F1") for k in range(1, 21)]   # 21 sequences
    regions = [(f"ctg{k}", 0, 10_000, f"Reg{k + 1}_G1_CLOW_S0.10") for k in range(21)]
    drawn, log, _, _ = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert len(drawn) == 21 and "dropped" not in log


# ------------------------------------------------- counted calls outside the view ----

def _capped_rows():
    """Ten HIGH copies at the locus on chr1 and two MEDIUM calls in a second neighbourhood."""
    rows = [_goi("chr1", 40_000 + 3000 * k, 41_000 + 3000 * k, ident=80.0 - k) for k in range(1, 11)]
    rows += [_flank("chr1", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(1, 4)]
    rows += [_goi("chr3", 50_000, 51_000, ident=60.0, conf="MEDIUM"),
             _goi("chr3", 55_000, 56_000, ident=58.0, conf="MEDIUM")]
    rows += [_flank("chr3", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(1, 5)]
    return rows, [("chr1", 0, 100_000, "Reg1_G13_CHIGH_S0.90"), ("chr3", 0, 100_000, "Reg2_G6_CHIGH_S0.70")]


def test_counted_calls_below_the_goi_cap_are_listed_outside_the_view(tmp_path, monkeypatch):
    rows, regions = _capped_rows()
    drawn, log, gff, genes = _select(tmp_path, rows, regions, cap=10, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 13}
    assert "left out Reg2_G6_CHIGH_S0.70 chr3:1-100000 (6 models, incl. 2 GOI models below the GOI cap)" in log
    outside = ps._counted_calls_outside_view(gff, genes)
    assert [(g["chrom"], g["start"], g["confidence"]) for g in outside] == [
        ("chr3", 50_000, "MEDIUM"), ("chr3", 55_000, "MEDIUM")]
    assert all(g["outside_view"] for g in outside)
    # with a cap that holds every call, both neighbourhoods are drawn and nothing is outside
    drawn, _, gff, genes = _select(tmp_path, rows, regions, cap=60, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 13, "chr3": 6}
    assert ps._counted_calls_outside_view(gff, genes) == []


def test_a_counted_call_on_a_sequence_without_neighbours_is_listed(tmp_path, monkeypatch):
    rows = [_goi("chr1", 50_000, 52_000)]
    rows += [_flank("chr1", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(1, 6)]
    rows += [_goi("chr5", 900_000, 902_000, ident=88.0), _goi("chr6", 10_000, 11_000, ident=30.0, conf="LOW")]
    regions = [("chr1", 0, 100_000, "Reg1_G6_CHIGH_S0.80"), ("chr5", 890_000, 910_000, "Reg2_DISP_CMEDIUM_S0.30"),
               ("chr6", 0, 20_000, "Reg3_G1_CLOW_S0.05")]
    drawn, log, gff, genes = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 6}
    assert "dropped 2 unsupported GOI hits on chr5, chr6" in log
    outside = ps._counted_calls_outside_view(gff, genes)      # the LOW call is no counted call
    assert [(g["chrom"], g["confidence"]) for g in outside] == [("chr5", "HIGH")]


def _grid_tracks(outside):
    names = [f"gene-F{k}" for k in range(6)]
    home = [{"chrom": "chr1", "start": 1000 + 2000 * k, "end": 2000 + 2000 * k, "start_plot": 1000 + 2000 * k,
             "end_plot": 2000 + 2000 * k, "name": n, "strand": "+", "exon_coords": []} for k, n in enumerate(names)]
    home.insert(3, {"chrom": "chr1", "start": 6600, "end": 6900, "start_plot": 6600, "end_plot": 6900,
                    "name": "GOI_T", "strand": "+", "exon_coords": []})

    def goi(start, conf="HIGH", ident=70.0, chrom="chr1", **kw):
        g = {"chrom": chrom, "start": start, "end": start + 300, "start_plot": start, "end_plot": start + 300,
             "name": f"GOI_T|m{chrom}_{start}", "gff_id": f"GOI_T|m{chrom}_{start}", "strand": "+",
             "role": "goi", "home_gene_id": "GOI_T", "confidence": conf, "identity": ident,
             "evidence_type": "exon_annotation", "query_coverage": 0.95,
             "exon_coords": [(start, start + 100), (start + 200, start + 300)]}
        g.update(kw)
        return g
    genes = [{"chrom": "chr1", "start": 1000 + 2000 * k, "end": 2000 + 2000 * k, "start_plot": 1000 + 2000 * k,
              "end_plot": 2000 + 2000 * k, "name": f"t_{k}", "strand": "+", "exon_coords": [],
              "home_gene_id": n, "identity": 80.0, "confidence": "HIGH", "query_coverage": 0.95}
             for k, n in enumerate(names)] + [goi(7000)]
    target = {"label": "T0", "species": "Bombus species0", "is_home": False, "genome_id": "g0", "offset": 0,
              "breaks": [], "genes": genes,
              "goi_outside_view": [goi(s, "MEDIUM", 60.0, chrom="chr3", outside_view=True) for s in outside]}
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home",
               "offset": 0, "breaks": [], "genes": home}, target]
    return tracks, {n: ps.GENE_PALETTE[k] for k, n in enumerate(names)}


def test_the_anchor_grid_lists_counted_calls_outside_the_view_after_the_break():
    svg, legend, source, info = ps.render_anchor_grid_print(*_grid_tracks([50_000, 55_000]), {},
                                                            SimpleNamespace())
    rows = [r for r in source if r["row"] == 2 and r["column"] == "T"]
    assert [(r["position"], r["scaffold"], r["drawn_as"]) for r in rows] == [
        ("locus", "chr1", "model"), ("elsewhere", "chr3", "model"), ("elsewhere", "chr3", "model")]
    assert info["goi_calls_outside_view"] == 2 and info["goi_copies"] == 3
    assert ">//</text>" in svg
    assert "Calls after // lie on another scaffold or farther from the locus than the neighbourhood is long." in legend
    assert ("2 counted calls (high or medium confidence) are listed that the ribbon plot and the "
            "gene-position map do not draw: they lie outside the stretch of their genome those "
            "figures show.") in legend
    # the full grid's tooltip says where the call lies
    full = ps.render_anchor_grid(*_grid_tracks([50_000]), {}, {}, SimpleNamespace())
    assert "outside the stretch the ribbon plot draws for this genome" in full
    # without such calls the grid is drawn as before
    plain_svg, plain_legend, plain_source, plain_info = ps.render_anchor_grid_print(
        *_grid_tracks([]), {}, SimpleNamespace())
    assert plain_info["goi_calls_outside_view"] == 0 and ">//</text>" not in plain_svg
    assert "counted call" not in plain_legend


def test_the_ribbon_plot_counts_them_and_the_genome_has_the_gene():
    tracks, _ = _grid_tracks([50_000, 55_000])
    assert ps._track_label_jobs(tracks[1], {})[0][1] == "+2 elsewhere"
    tracks[1]["genes"].append(dict(tracks[1]["genes"][-1], start=9000, end=9300, start_plot=9000,
                                   end_plot=9300, confidence="LOW", identity=30.0))
    assert ps._track_label_jobs(tracks[1], {})[0][1] == "+1 weak +2 elsewhere"
    texts = [t for _, rows in ps._ribbon_legend_groups("T", elsewhere=True) for _, t in rows]
    assert "counted calls outside the stretch drawn; listed in the anchor grid" in texts
    assert not any("outside the stretch" in t for _, rows in ps._ribbon_legend_groups("T")
                   for _, t in rows)
    # a genome whose only counted calls lie outside the view is no "GOI not placed" genome
    weak = {"genes": [dict(tracks[1]["genes"][-1])]}
    assert ps._track_goi_status(weak) == "ambiguous"
    assert ps._track_goi_status({"genes": weak["genes"] + tracks[1]["goi_outside_view"]}) == "resolved"


def test_a_weak_call_in_view_is_no_copy_when_the_genome_has_a_counted_call_elsewhere():
    """The copy rule names the best call of the genome: a LOW call beside the flanking genes is
    a weak call when a counted call of the genome lies outside the view, and the row's copy
    without one, as before."""
    tracks, colours = _grid_tracks([50_000])
    tracks[1]["genes"][-1].update(confidence="LOW", identity=30.0)
    _svg, _legend, source, _info = ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace())
    rows = [r for r in source if r["row"] == 2 and r["column"] == "T"]
    assert [(r["position"], r["confidence"], r["drawn_as"]) for r in rows] == [
        ("locus", "LOW", "mark"), ("elsewhere", "MEDIUM", "model")]
    assert ps._track_label_jobs(tracks[1], {})[0][1] == "+1 weak +1 elsewhere"
    html = ps.render_anchor_grid_positional(tracks, colours, {}, {}, SimpleNamespace())
    assert "1 weak call(s), 1 counted call(s) outside this neighbourhood (anchor grid)" in html
    tracks[1]["goi_outside_view"] = []
    _svg, _legend, source, _info = ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace())
    assert [r["drawn_as"] for r in source if r["row"] == 2 and r["column"] == "T"] == ["model"]
    assert ps._track_label_jobs(tracks[1], {}) == []


# ------------------------------------------------------------ the figures, end to end ----

HOME = [("FLA", 1000, 2000), ("FLB", 4000, 5000), ("QRY1", 7000, 8200), ("FLC", 10000, 11000),
        ("FLD", 13000, 14000)]


def test_a_counted_call_outside_the_view_reaches_the_figures(tmp_path):
    """plot_synteny.py on a run whose target holds a HIGH call on a scaffold without neighbours."""
    (tmp_path / "synteny_block_locus_1.bed").write_text("".join(
        f"hchr1\t{s - 1}\t{e}\tgene-{g}\t.\t+\t{g}\n" for g, s, e in HOME if g != "QRY1"))
    (tmp_path / "locus_1.bed").write_text("hchr1\t7100\t8100\tgene_loc\t1e-40\t+\n")
    home = ["##gff-version 3"]
    for g, s, e in HOME:
        home += [f"hchr1\tRefSeq\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{g};Name={g}",
                 f"hchr1\tRefSeq\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{g};Parent=gene-{g}",
                 f"hchr1\tRefSeq\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{g};Parent=rna-{g}"]
    (tmp_path / "home.gff").write_text("\n".join(home) + "\n")
    rows = []
    for g, s, e in HOME:
        if g == "QRY1":
            rows.append(_goi("tchr5", s + 500_000, e + 500_000).replace("GOI_X", "GOI_QRY1"))
        else:
            rows.append(_flank("tchr5", s + 500_000, e + 500_000, f"gene-{g}"))
    rows.append(_goi("tchr9", 900_000, 901_200, ident=91.0).replace("GOI_X", "GOI_QRY1"))
    gff = _gff(tmp_path, rows, "GCA_000000009.1.fna.gff")
    (tmp_path / "GCA_000000009.1.fna.regions.bed").write_text(
        "tchr5\t500000\t515000\tReg1_G5_CHIGH_S0.90\t0.90\t+\n"
        "tchr9\t899000\t902000\tReg2_DISP_CMEDIUM_S0.30\t0.30\t+\n")
    out = tmp_path / "synteny_block_locus_1_synteny_plot.html"
    proc = subprocess.run(
        [sys.executable, os.path.join(REPO, "bin", "plot_synteny.py"),
         "--home_bed", str(tmp_path / "synteny_block_locus_1.bed"),
         "--query_bed", str(tmp_path / "locus_1.bed"), "--home_gff", str(tmp_path / "home.gff"),
         "--target_gffs", gff, "--candidate_beds", str(tmp_path / "GCA_000000009.1.fna.regions.bed"),
         "--output", str(out), "--no_network", "--common_names", "off", "--no_fragment_variant"],
        capture_output=True, text=True, cwd=tmp_path,
        env=dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1"))
    assert proc.returncode == 0, proc.stderr
    assert ("[plot] GCA_000000009.1: 1 counted GOI call(s) outside the neighbourhood drawn: "
            "listed in the anchor grid\n") in proc.stdout
    with open(tmp_path / "synteny_block_locus_1_anchor_grid_print.source.tsv", newline="") as fh:
        goi = [r for r in csv.DictReader(fh, delimiter="\t")
               if r["genome_id"] == "GCA_000000009.1" and r["column"] == "QRY1"]
    assert [(r["position"], r["scaffold"], r["confidence"]) for r in goi] == [
        ("locus", "tchr5", "HIGH"), ("elsewhere", "tchr9", "HIGH")]
    assert "+1 elsewhere" in out.read_text()
    assert re.search(r"1 counted call \(high or medium confidence\) is listed that the ribbon plot", (
        tmp_path / "synteny_block_locus_1_anchor_grid_print.legend.txt").read_text())


# ------------------------------ what is named left out, and where off-view calls go ----

def test_a_region_drawn_through_a_kept_one_is_not_named_left_out(tmp_path, monkeypatch):
    """cluster_grs widens a GOI region to its search block, so a weaker region can lie inside
    it: its models are drawn, and the log named it left out with them."""
    rows = [_goi("chr1", 50_000, 52_000)] + [_flank("chr1", 10_000 * k, 10_000 * k + 900, f"gene-F{k}")
                                             for k in range(1, 9)]
    regions = [("chr1", 0, 100_000, "Reg1_G9_CHIGH_S0.9"), ("chr1", 65_000, 95_000, "Reg2_G3_CLOW_S0.2")]
    drawn, log, _, _ = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 9} and "left out" not in log
    # a region reaching beyond the kept one is named with the models it alone holds
    rows.append(_flank("chr1", 130_000, 130_900, "gene-F9"))
    regions[1] = ("chr1", 65_000, 140_000, "Reg2_G4_CLOW_S0.2")
    drawn, log, _, _ = _select(tmp_path, rows, regions, monkeypatch=monkeypatch)
    assert drawn == {"chr1": 9}
    assert "left out Reg2_G4_CLOW_S0.2 chr1:65001-140000 (1 model)" in log


def test_a_capped_copy_beside_the_drawn_array_is_listed_at_the_locus(tmp_path, monkeypatch):
    """Twelve copies at the end of a neighbourhood with the GOI cap at 10: the two weakest lie
    3 kb past the last copy drawn. They were drawn after '//', listed as 'elsewhere', and the
    legend said they lie on another scaffold or farther away."""
    rows = [_flank("chr1", 10_000 * k, 10_000 * k + 900, f"gene-F{k}") for k in range(0, 3)]
    rows += [_goi("chr1", 40_000 + 3000 * k, 41_000 + 3000 * k, ident=90.0 - k,
                  conf="HIGH" if k < 10 else "MEDIUM") for k in range(12)]
    regions = [("chr1", 0, 100_000, "Reg1_G3_CHIGH_S0.9")]
    _drawn, _log, gff, genes = _select(tmp_path, rows, regions, cap=10, monkeypatch=monkeypatch)
    outside = ps._counted_calls_outside_view(gff, genes)
    assert [(g["start"], g["confidence"]) for g in outside] == [(70_000, "MEDIUM"), (73_000, "MEDIUM")]
    names = [f"gene-F{k}" for k in range(3)]
    home = [{"chrom": "h", "start": 10_000 * k, "end": 10_000 * k + 900, "start_plot": 10_000 * k,
             "end_plot": 10_000 * k + 900, "name": n, "strand": "+", "exon_coords": []} for k, n in enumerate(names)]
    home.append({"chrom": "h", "start": 40_000, "end": 41_000, "start_plot": 40_000, "end_plot": 41_000,
                 "name": "GOI_X", "strand": "+", "exon_coords": []})
    for g in genes:
        g.setdefault("start_plot", g["start"])
        g.setdefault("end_plot", g["end"])
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home",
               "offset": 0, "breaks": [], "genes": home},
              {"label": "T0", "species": "Bombus x", "is_home": False, "genome_id": "G", "offset": 0,
               "breaks": [], "genes": genes, "goi_outside_view": outside}]
    svg, legend, source, info = ps.render_anchor_grid_print(
        tracks, {n: ps.GENE_PALETTE[k] for k, n in enumerate(names)}, {}, SimpleNamespace())
    goi = [(r["position"], int(r["start"]), r["confidence"]) for r in source
           if r["row"] == 2 and r["column"] == "X"]
    assert goi == [("locus", 40_000 + 3000 * k, "HIGH" if k < 10 else "MEDIUM") for k in range(12)]
    assert ">//</text>" not in svg and "Calls after //" not in legend
    assert info["goi_calls_outside_view"] == 2
    assert "2 counted calls (high or medium confidence) are listed that the ribbon plot" in legend
    # the ribbon plot, which does not draw them, still counts them
    assert ps._track_label_jobs(tracks[1], {})[0][1] == "×10 +2 elsewhere"


def _flipped_calls():
    """Three copies in view on a scaffold drawn reversed (plot positions mirrored), and an
    off-view copy beyond the highest genomic position (no plot position)."""
    def g(start, **kw):
        d = {"chrom": "c", "start": start, "end": start + 300, "strand": "-", "genomic_strand": "+",
             "name": f"GOI_X|m{start}", "gff_id": f"GOI_X|m{start}", "home_gene_id": "GOI_X",
             "role": "goi", "confidence": "HIGH", "identity": 80.0, "evidence_type": "exon_annotation"}
        d.update(kw)
        return d
    drawn = [g(s, start_plot=100_000 - s - 300, end_plot=100_000 - s) for s in (10_000, 12_000, 14_000)]
    off = g(16_000, confidence="MEDIUM", outside_view=True)
    return drawn, off


def test_off_view_calls_on_a_reversed_scaffold_take_its_order(monkeypatch):
    monkeypatch.setattr(sg, "_PS", ps)
    drawn, off = _flipped_calls()
    tmap = {sg.GOI_KEY: {"best": drawn[0], "genes": drawn, "n": 3}}
    _best, calls = sg.split_calls(tmap, 50_000, [off])
    # drawn right to left in genomic terms: 14 kb, 12 kb, 10 kb; the off-view copy at 16 kb comes first
    assert [(g["start"], away) for g, away, _c in calls] == [
        (16_000, False), (14_000, False), (12_000, False), (10_000, False)]
    # on a scaffold drawn forwards it follows the copies
    for d in drawn + [off]:
        d.pop("genomic_strand")
        d["strand"] = "+"
    for d in drawn:
        d["start_plot"], d["end_plot"] = d["start"], d["end"]
    _best, calls = sg.split_calls(tmap, 50_000, [off])
    assert [g["start"] for g, _a, _c in calls] == [10_000, 12_000, 14_000, 16_000]
    # far away on the same scaffold: after //
    far = dict(off, start=900_000, end=900_300)
    _best, calls = sg.split_calls(tmap, 5_000, [far])
    assert [(g["start"], away) for g, away, _c in calls][-1] == (900_000, True)


def test_off_view_calls_with_no_call_in_view_are_placed_by_the_flanking_genes(monkeypatch):
    monkeypatch.setattr(sg, "_PS", ps)
    near = {"chrom": "c", "start": 25_000, "end": 25_300, "strand": "+", "gff_id": "GOI_X|a",
            "home_gene_id": "GOI_X", "role": "goi", "confidence": "HIGH", "identity": 80.0,
            "evidence_type": "exon_annotation", "outside_view": True}
    other = dict(near, chrom="d", start=5_000, end=5_300, gff_id="GOI_X|b", confidence="MEDIUM")
    flank = [{"chrom": "c", "start": s, "end": s + 900, "start_plot": s, "end_plot": s + 900}
             for s in (10_000, 20_000)]
    tmap = {f"gene-F{k}": {"best": f, "genes": [f], "n": 1} for k, f in enumerate(flank)}
    best, calls = sg.split_calls(tmap, 20_000, [other, near])
    assert best is near
    assert [(g["gff_id"], away) for g, away, _c in calls] == [("GOI_X|a", False), ("GOI_X|b", True)]
    assert sg.split_calls({}, 20_000, [near])[1] == [(near, True, True)]   # nothing drawn: no locus


def test_off_view_calls_are_mirrored_with_their_scaffold():
    """_orient_tracks_to_home flipped only the genes in view, so an off-view copy on a scaffold
    drawn reversed kept its genomic strand and pointed against the copies beside it."""
    home = [{"chrom": "h", "start": 1000 * k, "end": 1000 * k + 500, "start_plot": 1000 * k,
             "end_plot": 1000 * k + 500, "name": f"gene-F{k}", "strand": "+"} for k in range(4)]
    target = [{"chrom": "t", "start": 10_000 - 1000 * k, "end": 10_500 - 1000 * k,
               "start_plot": 10_000 - 1000 * k, "end_plot": 10_500 - 1000 * k, "name": f"t{k}",
               "home_gene_id": f"gene-F{k}", "strand": "-"} for k in range(4)]
    off = {"chrom": "t", "start": 50_000, "end": 51_000, "strand": "-", "name": "GOI_X|m",
           "home_gene_id": "GOI_X", "role": "goi", "exon_coords": [(50_000, 50_200), (50_700, 51_000)],
           "outside_view": True}
    elsewhere = dict(off, chrom="u", exon_coords=[(50_000, 50_200)])
    tracks = [{"is_home": True, "genes": home},
              {"is_home": False, "genes": target, "breaks": [], "goi_outside_view": [off, elsewhere]}]
    assert ps._orient_tracks_to_home(tracks) == 1
    assert tracks[1]["flipped_chroms"] == ["t"]
    assert (off["strand"], off["genomic_strand"]) == ("+", "-")
    assert off["exon_coords"] == [(50_000, 50_300), (50_800, 51_000)]
    assert elsewhere["strand"] == "-" and "genomic_strand" not in elsewhere   # another scaffold


def test_an_off_view_call_on_a_reversed_scaffold_points_as_it_is_drawn(tmp_path):
    """plot_synteny.py end to end: the neighbourhood reversed on tchr5 and a second HIGH call
    2 Mb away on the same scaffold, both on '-'."""
    loci = [("FLA", 1000, 2000), ("FLB", 4000, 5000), ("QRY1", 7000, 8200), ("FLC", 10000, 11000),
            ("FLD", 13000, 14000)]
    (tmp_path / "synteny_block_locus_1.bed").write_text("".join(
        f"hchr1\t{s - 1}\t{e}\tgene-{g}\t.\t+\t{g}\n" for g, s, e in loci if g != "QRY1"))
    (tmp_path / "locus_1.bed").write_text("hchr1\t7100\t8100\tgene_loc\t1e-40\t+\n")
    home = ["##gff-version 3"]
    for g, s, e in loci:
        home += [f"hchr1\tRefSeq\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{g};Name={g}",
                 f"hchr1\tRefSeq\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{g};Parent=gene-{g}",
                 f"hchr1\tRefSeq\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{g};Parent=rna-{g}"]
    (tmp_path / "home.gff").write_text("\n".join(home) + "\n")
    rows = []
    for g, s, e in loci:
        ts, te = 520_000 - e, 520_000 - s
        row = (_goi("tchr5", ts, te).replace("GOI_X", "GOI_QRY1") if g == "QRY1"
               else _flank("tchr5", ts, te, f"gene-{g}"))
        rows.append(row.replace("\t+\t.\t", "\t-\t.\t"))
    rows.append(_goi("tchr5", 2_500_000, 2_501_200, ident=91.0).replace("GOI_X", "GOI_QRY1")
                .replace("\t+\t.\t", "\t-\t.\t"))
    gff = _gff(tmp_path, rows, "GCA_000000009.1.fna.gff")
    (tmp_path / "GCA_000000009.1.fna.regions.bed").write_text(
        "tchr5\t500000\t520000\tReg1_G5_CHIGH_S0.90\t0.90\t-\n")
    out = tmp_path / "synteny_block_locus_1_synteny_plot.html"
    proc = subprocess.run(
        [sys.executable, os.path.join(REPO, "bin", "plot_synteny.py"),
         "--home_bed", str(tmp_path / "synteny_block_locus_1.bed"),
         "--query_bed", str(tmp_path / "locus_1.bed"), "--home_gff", str(tmp_path / "home.gff"),
         "--target_gffs", gff, "--candidate_beds", str(tmp_path / "GCA_000000009.1.fna.regions.bed"),
         "--output", str(out), "--no_network", "--common_names", "off", "--no_fragment_variant"],
        capture_output=True, text=True, cwd=tmp_path,
        env=dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1"))
    assert proc.returncode == 0, proc.stderr
    assert "Oriented 1 scaffold(s) to home" in proc.stdout
    with open(tmp_path / "synteny_block_locus_1_anchor_grid_print.source.tsv", newline="") as fh:
        goi = [r for r in csv.DictReader(fh, delimiter="\t")
               if r["genome_id"] == "GCA_000000009.1" and r["column"] == "QRY1"]
    assert [(r["position"], r["start"], r["strand_drawn"]) for r in goi] == [
        ("locus", "511800", "+"), ("elsewhere", "2500000", "+")]
    html = (tmp_path / "synteny_block_locus_1_anchor_grid.html").read_text()
    tips = [t for t in re.findall(r"<title>([^<]*)</title>", html) if "2,500,0" in t or "2500000" in t]
    assert tips and all("scaffold drawn reversed" in t for t in tips)
