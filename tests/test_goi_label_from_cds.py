"""The figures name the GOI after a home gene only when the query hits fall on that gene's CDS.

Without an annotated gene under the hits' CDS (a curated seed, an ORF, a gene in an intron) the
home GOI is the gene the span overlap finds (identify_goi_names), and that can be an lncRNA over
the hit. The home row's placeholder and the GOI's own id were already guarded (freeze fix e3);
the GOI column, the figure titles and the legends (_goi_display_label) and the ribbon plot's
home label still took that gene's name. Now they take the GOI's own name, as its target models
carry it, and the annotated gene names the GOI only on the CDS path.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402

ROWS = {"gene-NB1": {"label": "NB1", "source": "gff_symbol", "full_name": "", "ncbi_gene_id": ""},
        "gene-lncX": {"label": "lncX", "source": "gff_symbol", "full_name": "long non-coding RNA X",
                      "ncbi_gene_id": ""},
        "gene-NB2": {"label": "NB2", "source": "gff_symbol", "full_name": "", "ncbi_gene_id": ""}}


@pytest.fixture
def state(monkeypatch):
    monkeypatch.setattr(ps, "_GENE_NAMES", dict(ROWS))
    monkeypatch.setattr(ps, "_GOI_NAMES", set())
    monkeypatch.setattr(ps, "_GOI_FROM_CDS", [False])
    monkeypatch.setattr(ps, "_GOI_QUERY_NAME", [None])
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [None])


def _home():
    return [{"chrom": "c", "start": 1000, "end": 2000, "name": "gene-NB1", "display_name": "NB1"},
            # an lncRNA over the hit span: no CDS, in the home BED as a neighbour
            {"chrom": "c", "start": 4800, "end": 6200, "name": "gene-lncX", "display_name": "lncX"},
            {"chrom": "c", "start": 9000, "end": 9900, "name": "gene-NB2", "display_name": "NB2"}]


HITS = [{"chrom": "c", "start": 5000, "end": 5600, "strand": "+"}]


def _tracks(home, target_names=("GOI_SEED7",)):
    genes = [{"name": f"{n}|Apis_florea_fna_b0_l1_exon_ann", "home_gene_id": n} for n in target_names]
    return [{"is_home": True, "genes": home}, {"is_home": False, "genes": genes}]


def test_a_gene_the_span_overlap_found_does_not_name_the_goi(state):
    home = _home()
    ps.set_goi_names(home, HITS, [])              # nothing under the hits' CDS: the span rule
    assert ps._GOI_NAMES == {"gene-lncX"} and ps._GOI_FROM_CDS[0] is False
    assert ps._goi_display_label(_tracks(home)) == "SEED7"


def test_the_gene_under_the_hits_cds_names_the_goi(state):
    home = _home()
    ps.set_goi_names(home, HITS, ["gene-NB1"])
    assert ps._GOI_FROM_CDS[0] is True
    assert ps._goi_display_label(_tracks(home)) == "NB1"


def test_without_a_model_name_the_goi_is_goi(state):
    home = _home()
    ps.set_goi_names(home, HITS, [])
    assert ps._goi_display_label(_tracks(home, target_names=())) == "GOI"


def test_the_ribbon_plot_names_the_goi_as_the_other_figures(state, monkeypatch):
    home = [dict(g) for g in _home()]
    ps.set_goi_names(home, HITS, [])
    tracks = _tracks(home)
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [ps._goi_display_label(tracks)])
    (gene, text), = ps._track_label_jobs(tracks[0], {})
    assert gene["name"] == "gene-lncX" and text == "★ SEED7"
    # on the CDS path the gene's own label, as before
    monkeypatch.setattr(ps, "_GOI_NAMES", {"gene-NB1"})
    monkeypatch.setattr(ps, "_GOI_FROM_CDS", [True])
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [ps._goi_display_label(tracks)])
    (gene, text), = ps._track_label_jobs(tracks[0], {})
    assert gene["name"] == "gene-NB1" and text == "★ NB1"


# ---- review of item 42: the Gene names table and <locus>_gene_names.tsv -------------------

def _own_name_state(monkeypatch, home):
    ps.set_goi_names(home, HITS, [])
    tracks = _tracks(home)
    monkeypatch.setattr(ps, "_GOI_QUERY_NAME", [ps._goi_query_name(tracks)])
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [ps._goi_display_label(tracks)])
    return tracks


def test_the_names_table_names_the_goi_as_the_goi_column(state, monkeypatch):
    """The GOI column read 'SEED7' while the names table under it listed the lncRNA the span
    overlap found as the GOI row, in GOI red: the figure named the GOI two ways."""
    home = _home()
    _own_name_state(monkeypatch, home)
    rows = ps._gene_name_rows(home, {"gene-NB1": "#4e79a7", "gene-NB2": "#59a14f"})
    assert [(r["label"], r["is_goi"]) for r in rows] == [("NB1", False), ("SEED7", True), ("NB2", False)]
    goi = rows[1]
    assert goi["full"] == "query hits on no annotated CDS (home row: lncX)"
    assert goi["source"] == "query name" and goi["colour"] == ps.GOI_COLOUR
    assert not any("long non-coding" in r["full"] for r in rows)
    # on the CDS path the annotated gene is the GOI row, and the lncRNA an ordinary gene
    ps.set_goi_names(home, HITS, ["gene-NB1"])
    rows = ps._gene_name_rows(home, {"gene-NB1": "#4e79a7", "gene-NB2": "#59a14f"})
    assert [(r["label"], r["is_goi"]) for r in rows] == [("NB1", True), ("lncX", False), ("NB2", False)]


def test_a_bare_query_hit_without_a_name_keeps_its_row(state, monkeypatch):
    home = [g for g in _home() if g["name"] != "gene-lncX"]
    home.insert(1, {"chrom": "c", "start": 5000, "end": 5600, "name": "GOI_c_5000", "display_name": ""})
    ps.set_goi_names(home, HITS, [])
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [ps._goi_display_label(_tracks(home, ()))])
    rows = ps._gene_name_rows(home, {})
    assert [(r["label"], r["full"]) for r in rows if r["is_goi"]] == [
        ("GOI", "query hit; no annotated home gene")]


def test_the_names_tsv_has_the_goi_row_and_reads_it_back(tmp_path, state, monkeypatch):
    from types import SimpleNamespace
    home = _home()
    monkeypatch.setattr(ps, "_GENE_NAMES_ORDER", [g["name"] for g in home])
    monkeypatch.setattr(ps, "_GENE_NAMES_NOTE", ["skipped (--no_network)"])
    monkeypatch.setattr(ps, "_GOI_NAME_ROWS", {})
    tracks = _own_name_state(monkeypatch, home)
    out = tmp_path / "locus_1_gene_names.tsv"
    args = SimpleNamespace(gene_names_out=str(out), output=str(tmp_path / "locus_1_synteny_plot.html"),
                           home_gff=None, gene_names_tsv="", no_network=True, gene_name_lookup=False)
    ps._write_gene_names_tsv(args)
    lines = out.read_text().splitlines()
    assert "GOI_SEED7\tSEED7\t\tquery_name\t" in lines
    assert "gene-lncX\tlncX\tlong non-coding RNA X\tgff_symbol\t" in lines   # the gene keeps its own row
    # passed back unchanged, the table gives the same name; an edited label renames the GOI
    args.gene_names_tsv = str(out)
    ps._load_gene_names(args, home)
    assert ps._goi_display_label(tracks) == "SEED7"
    edited = tmp_path / "edited.tsv"
    edited.write_text("gene_id\tlabel\tfull_name\tsource\nGOI_SEED7\tSeed-7\tseed toxin 7\tuser\n")
    args.gene_names_tsv = str(edited)
    ps._load_gene_names(args, home)
    assert ps._goi_display_label(tracks) == "Seed-7"
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", ["Seed-7"])
    goi, = [r for r in ps._gene_name_rows(home, {}) if r["is_goi"]]
    assert (goi["label"], goi["full"], goi["source"]) == (
        "Seed-7", "seed toxin 7 (home row: lncX)", "names table")
    # on the CDS path no GOI row: the annotated gene's own row names the GOI
    ps.set_goi_names(home, HITS, ["gene-NB1"])
    ps._write_gene_names_tsv(args)
    assert not any(ln.startswith("GOI_") for ln in out.read_text().splitlines())


def test_the_figures_and_the_names_tsv_name_the_goi_one_way(tmp_path):
    """plot_synteny.py end to end: a curated seed GOI_SEED7 whose home hit lies inside an lncRNA."""
    import re
    import subprocess

    repo = os.path.join(os.path.dirname(__file__), "..")
    flanks = [("FLA", 1000, 2000), ("FLB", 4000, 5000), ("FLC", 10000, 11000), ("FLD", 13000, 14000)]

    def inputs(d):
        d.mkdir()
        bed = [f"hchr1\t{s - 1}\t{e}\tgene-{g}\t.\t+\t{g}\n" for g, s, e in flanks]
        bed.insert(2, "hchr1\t6499\t8800\tgene-LNC9\t.\t+\tLNC9\n")
        (d / "synteny_block_locus_1.bed").write_text("".join(bed))
        (d / "locus_1.bed").write_text("hchr1\t7100\t8100\tgene_loc\t1e-40\t+\n")
        home = ["##gff-version 3"]
        for g, s, e in flanks:
            home += [f"hchr1\tRefSeq\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{g};Name={g}",
                     f"hchr1\tRefSeq\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{g};Parent=gene-{g};product={g} protein",
                     f"hchr1\tRefSeq\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{g};Parent=rna-{g}"]
        home += ["hchr1\tRefSeq\tgene\t6500\t8800\t.\t+\t.\tID=gene-LNC9;Name=LNC9;gene_biotype=lncRNA;"
                 "description=long non-coding RNA 9",
                 "hchr1\tRefSeq\tlnc_RNA\t6500\t8800\t.\t+\t.\tID=rna-LNC9;Parent=gene-LNC9",
                 "hchr1\tRefSeq\texon\t6500\t8800\t.\t+\t.\tID=exon-LNC9;Parent=rna-LNC9"]
        (d / "home.gff").write_text("\n".join(home) + "\n")
        rows = [f"tchr5\tflanking_miniprot\tmRNA\t{s + 500_000}\t{e + 500_000}\t80\t+\t.\t"
                f"ID=gene-{g}|G_fna_tchr5_{s}_flank_ann;Name=gene-{g};SynVoy_Parent=gene-{g};Identity=80;"
                f"Confidence=HIGH;SynVoyRole=flanking;EvidenceType=flanking_miniprot" for g, s, e in flanks]
        rows.append("tchr5\texon_annotation\tmRNA\t507100\t508100\t70\t+\t.\tID=GOI_SEED7|G_fna_tchr5_507100_"
                    "exon_ann;Name=GOI_SEED7;SynVoy_Parent=GOI_SEED7;Identity=70;Confidence=HIGH;SynVoyRole=goi;"
                    "EvidenceType=exon_annotation;ModelStatus=complete;QueryCoverage=0.95")
        (d / "G.fna.gff").write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
        (d / "G.fna.regions.bed").write_text("tchr5\t500000\t515000\tReg1_G5_CHIGH_S0.90\t0.90\t+\n")

    def draw(d, *extra):
        proc = subprocess.run(
            [sys.executable, os.path.join(repo, "bin", "plot_synteny.py"),
             "--home_bed", str(d / "synteny_block_locus_1.bed"), "--query_bed", str(d / "locus_1.bed"),
             "--home_gff", str(d / "home.gff"), "--target_gffs", str(d / "G.fna.gff"),
             "--candidate_beds", str(d / "G.fna.regions.bed"),
             "--output", str(d / "synteny_block_locus_1_synteny_plot.html"), "--no_network",
             "--common_names", "off", "--no_fragment_variant", *extra],
            capture_output=True, text=True, cwd=d,
            env=dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1"))
        assert proc.returncode == 0, proc.stderr

    first = tmp_path / "first"
    inputs(first)
    draw(first)
    tsv = (first / "synteny_block_locus_1_gene_names.tsv").read_text()
    assert "GOI_SEED7\tSEED7\t\tquery_name\t" in tsv.splitlines()
    for page in ("anchor_grid", "gene_positions"):
        html = (first / f"synteny_block_locus_1_{page}.html").read_text()
        table = html[html.index('class="gene-names-table"'):]
        titles = re.findall(r"<title>([^<]*)</title>", table)
        assert "SEED7 · query hits on no annotated CDS (home row: LNC9) · query name" in titles
        assert not any(t.startswith("LNC9") for t in titles)
    # a re-draw from the run's own table names everything as the first drawing did
    again = tmp_path / "again"
    inputs(again)
    draw(again, "--gene_names_tsv", str(first / "synteny_block_locus_1_gene_names.tsv"))
    assert (again / "synteny_block_locus_1_gene_names.tsv").read_text() == tsv
    assert ((again / "synteny_block_locus_1_anchor_grid.svg").read_text().replace(str(again), "")
            == (first / "synteny_block_locus_1_anchor_grid.svg").read_text().replace(str(first), ""))


def test_the_tree_page_names_the_home_leaf_as_the_other_figures(tmp_path, state, monkeypatch):
    """An edited GOI_<name> row of the names table renames the GOI in the tree page's home
    leaf too; without an edit the leaf reads as before."""
    import re
    nwk = tmp_path / "locus_1_tree.nwk"
    nwk.write_text("(GOI_SEED7:0.1,(GOI_SEED7|A_fna_b0_l1:0.1,GOI_SEED7|B_fna_b0_l1:0.1):0.2,"
                   "GOI_SEED7|C_fna_b0_l1:0.3);\n")
    out = tmp_path / "synteny_block_locus_1_tree.html"
    monkeypatch.setattr(ps, "_GOI_QUERY_NAME", ["SEED7"])

    def home_leaf(display):
        monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [display])
        ps._render_tree_svg(str(nwk), {}, str(out), species_map={}, gene_legend=False)
        return [t for t in re.findall(r"<text[^>]*>([^<]*)</text>", out.read_text()) if "(home)" in t]
    assert home_leaf("SEED7") == ["SEED7 (home)"]
    assert home_leaf("Seed-7") == ["Seed-7 (home)"]
    monkeypatch.setattr(ps, "_GOI_FROM_CDS", [True])      # the CDS path: the leaf's own gene
    assert home_leaf("Seed-7") == ["SEED7 (home)"]
