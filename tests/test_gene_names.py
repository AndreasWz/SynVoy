"""Gene names in the figures.

The figures used to print whatever identifier a gene happened to carry: the raw
gene ID in the anchor grid (`Dmel_CG4491` although the GFF says `noc`), a
transcript accession in the home track (`XM_006567126.3`), `LOC...` nearly
everywhere for genomes annotated without symbols. Every figure now labels a home
gene from one names table: the symbol in the home GFF, else the current symbol at
NCBI Gene, else an abbreviation of the product name (marked, because it is not an
official symbol), else the ID. The table is also drawn under each figure.
"""
from __future__ import annotations

import io
import json
import os
import re
import sys
import urllib.error
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402
import plot_synteny_matrix as pm  # noqa: E402
import synvoy_gene_names as gn  # noqa: E402


# ---- the abbreviation rule -------------------------------------------------

@pytest.mark.parametrize("product, label", [
    # every word gives its first letter; numbers stay
    ("pancreatic triacylglycerol lipase", "PTL"),
    ("coiled-coil domain-containing protein 170", "CCDC170"),
    ("ceramide-1-phosphate transfer protein", "C1PT"),
    ("peptidylglycine alpha-hydroxylating monooxygenase", "PAHM"),
    ("cyclin-dependent kinase 8", "CDK8"),
    ("structural maintenance of chromosomes protein 6", "SMC6"),
    # "like" gives L where it stands
    ("zinc finger protein 629-like", "ZF629L"),
    ("transcription factor CP2-like protein 1", "TFCP2L1"),
    # an existing abbreviation or a designator is kept whole
    ("5'-AMP-activated protein kinase subunit beta-2", "5AMPAKSB2"),
    ("ADP-ribosylation factor-like protein 6", "ADPRFL6"),
    ("vegetative cell wall protein gp1", "VCWGP1"),
    ("28S ribosomal protein S14, mitochondrial", "28SRS14"),
    # a one- or two-word name would be too short: the first word gives three letters
    ("arginine kinase", "ARGK"),
    ("cytochrome b5", "CYTB5"),
    ("apolipoprotein D", "APOD"),
    ("luciferin 4-monooxygenase", "LUC4M"),
    # a gene symbol at the end of the name is the label
    ("protein LTV1 homolog", "LTV1"),
    ("E3 ubiquitin-protein ligase MYLIP", "MYLIP"),
    ("probable C-mannosyltransferase DPY19L1", "DPY19L1"),
    ("segmentation protein Runt", "Runt"),
    ("neurogenic locus Notch protein", "Notch"),
    # a short one-word name is kept as it is
    ("protein lozenge", "lozenge"),
    ("apyrase", "apyrase"),
    ("formin-2", "formin-2"),
    ("titin-like", "titin-like"),
    ("melittin", "melittin"),
])
def test_abbreviation_rule(product, label):
    assert gn.abbreviate_product(product) == label


@pytest.mark.parametrize("product", [
    "uncharacterized LOC726866",
    "uncharacterized protein LOC116158585 isoform X1",
    "hypothetical protein",
    "",
])
def test_a_product_that_names_nothing_gives_no_label(product):
    assert gn.abbreviate_product(product) == ""


def test_transcript_and_compartment_qualifiers_do_not_reach_the_label():
    assert gn.abbreviate_product(
        "pancreatic triacylglycerol lipase, transcript variant X2") == "PTL"
    assert gn.clean_product("alcohol dehydrogenase, isoform C") == "alcohol dehydrogenase"
    assert gn.clean_product("malate dehydrogenase, mitochondrial") == "malate dehydrogenase"
    assert gn.clean_product("decorin [Source:HGNC Symbol;Acc:HGNC:2705]") == "decorin"


def test_a_designator_is_not_taken_for_a_gene_symbol():
    # "S14", "B9", "IIH" end a name without being its symbol.
    assert gn.abbreviate_product("annexin B9") == "annexin B9"
    assert gn.abbreviate_product("general transcription factor IIH subunit 3") == "GTFIIHS3"
    assert gn.abbreviate_product("casein kinase I-like") == "CKIL"
    # a chemical locant does not make a word a symbol
    assert gn.abbreviate_product("protein-L-isoaspartate O-methyltransferase") == "LIOM"


def test_a_one_word_name_is_not_cut_to_three_letters():
    assert gn.abbreviate_product("junctophilin-1") == "junctophilin-1"
    assert gn.abbreviate_product("tetraspanin-12") == "tetraspanin-12"
    assert gn.abbreviate_product("4-nitrophenylphosphatase") == "4NIT"      # too long to keep


def test_labels_never_exceed_the_column_header_budget():
    long_name = ("NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 2 "
                 "regulatory subunit gamma delta epsilon")
    assert len(gn.abbreviate_product(long_name)) <= gn.MAX_LABEL_CHARS
    assert len(gn.abbreviate_product("gastrula zinc finger protein XlCGF46.1-like")) \
        <= gn.MAX_LABEL_CHARS


@pytest.mark.parametrize("label, generic", [
    ("LOC726817", True), ("gene-LOC726817", True), ("XM_006567126.3", True),
    ("Dmel_CG4491", True), ("XY_0012345", True), ("ENSG00000011465", True),
    ("g123.t1", True), ("", True),
    ("Melt", False), ("Adh", False), ("ACT1", False), ("CG33090", False),
    ("G6PD", False), ("alc", False),
])
def test_identifiers_are_told_from_gene_symbols(label, generic):
    assert gn.is_generic_id(label) is generic


# ---- the home GFF ----------------------------------------------------------

GFF = "\n".join([
    "##gff-version 3",
    # NCBI gene with a symbol
    "chr1\tGnomon\tgene\t100\t900\t.\t+\t.\tID=gene-Melt;Dbxref=GeneID:406130;Name=Melt;"
    "description=melittin;gene=Melt;gene_biotype=protein_coding",
    # NCBI gene without one: the product is on the transcript
    "chr1\tGnomon\tgene\t1000\t1900\t.\t+\t.\tID=gene-LOC726817;Dbxref=BEEBASE:GB44098,"
    "GeneID:726817;Name=LOC726817;gene=LOC726817;gene_biotype=protein_coding",
    "chr1\tGnomon\tmRNA\t1000\t1900\t.\t+\t.\tID=rna-XM_006567126.3;Parent=gene-LOC726817;"
    "Name=XM_006567126.3;gene=LOC726817;product=pancreatic triacylglycerol lipase%2C "
    "transcript variant X2",
    # nothing names this one
    "chr1\tGnomon\tgene\t2000\t2900\t.\t-\t.\tID=gene-LOC100578368;Dbxref=GeneID:100578368;"
    "Name=LOC100578368;gene=LOC100578368",
    "chr1\tGnomon\tmRNA\t2000\t2900\t.\t-\t.\tID=rna-XM_1.1;Parent=gene-LOC100578368;"
    "product=uncharacterized LOC100578368%2C transcript variant X1",
    # FlyBase style: ID and locus_tag are the annotation ID, gene= is the symbol
    "chr1\tRefSeq\tgene\t3000\t3900\t.\t+\t.\tID=gene-Dmel_CG4491;Dbxref=FLYBASE:FBgn0005771,"
    "GeneID:34889;Name=noc;description=no ocelli;gene=noc;locus_tag=Dmel_CG4491",
    # yeast ORF without a gene name: Name repeats the locus tag
    "chr1\tRefSeq\tgene\t4000\t4900\t.\t+\t.\tID=gene-YFL019C;Dbxref=GeneID:850525;"
    "Name=YFL019C;locus_tag=YFL019C",
]) + "\n"


@pytest.fixture
def gff(tmp_path):
    path = tmp_path / "home.gff"
    path.write_text(GFF)
    return str(path)


def test_gff_index_separates_symbols_from_identifiers(gff):
    info = gn.load_gff_gene_info(gff)
    assert info["gene-Melt"] == {"symbol": "Melt", "product": "melittin",
                                 "ncbi_gene_id": "406130"}
    assert info["gene-LOC726817"]["symbol"] == ""
    assert info["gene-LOC726817"]["product"].startswith("pancreatic triacylglycerol lipase")
    assert info["gene-LOC726817"]["ncbi_gene_id"] == "726817"
    assert info["gene-Dmel_CG4491"]["symbol"] == "noc"
    assert info["gene-YFL019C"]["symbol"] == ""          # Name == locus_tag
    assert gn.load_gff_gene_info(gff, wanted=["gene-Melt"]).keys() == {"gene-Melt"}
    assert gn.load_gff_gene_info("NO_GFF") == {}


# ---- resolution ------------------------------------------------------------

IDS = ["gene-Melt", "gene-LOC726817", "gene-LOC100578368", "gene-Dmel_CG4491", "gene-YFL019C"]


def _no_fetch(ids):
    raise AssertionError(f"NCBI must not be asked ({ids})")


def test_sources_are_used_in_order(gff):
    asked = []

    def fetch(ids):
        asked.extend(ids)
        return {"726817": ("LOC726817", "pancreatic triacylglycerol lipase"),
                "100578368": ("MESR3", "misexpression suppressor of ras 3"),
                "850525": ("YFL019C", "uncharacterized protein")}

    table, note = gn.resolve_gene_names(IDS, gn.load_gff_gene_info(gff), fetch=fetch)
    # only the genes the GFF leaves unnamed go to NCBI
    assert asked == ["726817", "850525", "100578368"]
    assert note.startswith("ok")
    assert (table["gene-Melt"]["label"], table["gene-Melt"]["source"]) == ("Melt", gn.SOURCE_GFF)
    assert table["gene-Dmel_CG4491"]["label"] == "noc"
    assert (table["gene-LOC100578368"]["label"],
            table["gene-LOC100578368"]["source"]) == ("MESR3", gn.SOURCE_NCBI)
    assert table["gene-LOC100578368"]["full_name"] == "misexpression suppressor of ras 3"
    # NCBI has no symbol either: abbreviation of the product name, marked in the figures
    assert (table["gene-LOC726817"]["label"],
            table["gene-LOC726817"]["source"]) == ("PTL", gn.SOURCE_PRODUCT)
    assert gn.display_label(table["gene-LOC726817"]) == "PTL" + gn.DERIVED_MARK
    assert gn.display_label(table["gene-LOC100578368"]) == "MESR3"
    # a "symbol" that repeats the ID names nothing
    assert (table["gene-YFL019C"]["label"],
            table["gene-YFL019C"]["source"]) == ("YFL019C", gn.SOURCE_ID)
    assert table["gene-YFL019C"]["full_name"] == ""


def test_ncbi_is_not_asked_when_every_gene_has_a_symbol(gff):
    table, note = gn.resolve_gene_names(["gene-Melt", "gene-Dmel_CG4491"],
                                        gn.load_gff_gene_info(gff), fetch=_no_fetch)
    assert note.startswith("not needed")
    assert {r["source"] for r in table.values()} == {gn.SOURCE_GFF}


def test_offline_run_falls_back_to_the_gff_and_says_so(gff):
    table, note = gn.resolve_gene_names(IDS, gn.load_gff_gene_info(gff),
                                        allow_network=False, fetch=_no_fetch)
    assert note.startswith("skipped")
    assert table["gene-LOC726817"]["label"] == "PTL"
    assert table["gene-LOC100578368"]["source"] == gn.SOURCE_ID
    assert table["gene-LOC100578368"]["label"] == "LOC100578368"


def test_a_failed_lookup_is_reported_not_fatal(gff):
    def fetch(ids):
        raise OSError("NCBI Gene lookup failed: timed out")

    table, note = gn.resolve_gene_names(IDS, gn.load_gff_gene_info(gff), fetch=fetch)
    assert note.startswith("FAILED") and "timed out" in note
    assert table["gene-LOC726817"]["label"] == "PTL"


def test_names_table_overrides_every_other_source(gff):
    overrides = {"gene-LOC726817": {"label": "Ptl", "full_name": "my lipase",
                                    "source": gn.SOURCE_USER, "ncbi_gene_id": ""},
                 "gene-LOC100578368": {"label": "Unk1", "full_name": "",
                                       "source": gn.SOURCE_USER, "ncbi_gene_id": ""}}
    table, note = gn.resolve_gene_names(
        ["gene-Melt", "gene-LOC726817", "gene-LOC100578368"],
        gn.load_gff_gene_info(gff), overrides=overrides, fetch=_no_fetch)
    assert table["gene-LOC726817"]["label"] == "Ptl"
    assert gn.display_label(table["gene-LOC726817"]) == "Ptl"       # a user label is not marked
    assert note.startswith("not needed")


def test_two_genes_never_share_a_derived_label():
    info = {f"gene-LOC{i}": {"symbol": "", "product": "protein lozenge", "ncbi_gene_id": ""}
            for i in (1, 2)}
    info["gene-LOC3"] = {"symbol": "", "product": "pancreatic triacylglycerol lipase",
                         "ncbi_gene_id": ""}
    table, _ = gn.resolve_gene_names(["gene-LOC1", "gene-LOC2", "gene-LOC3"], info,
                                     allow_network=False)
    assert [table[g]["label"] for g in ("gene-LOC1", "gene-LOC2", "gene-LOC3")] == \
        ["lozenge-1", "lozenge-2", "PTL"]


# ---- NCBI Gene -------------------------------------------------------------

def test_ncbi_reports_are_parsed(monkeypatch):
    seen = {}

    def fake_urlopen(req, timeout=None):
        seen["url"], seen["timeout"] = req.full_url, timeout
        return io.BytesIO(json.dumps({"reports": [
            {"gene": {"gene_id": "409662", "symbol": "alc",
                      "description": "5'-AMP-activated protein kinase subunit beta-1"}},
            {"gene": {"gene_id": "726817", "symbol": "LOC726817",
                      "description": "pancreatic triacylglycerol lipase"}},
        ]}).encode())

    monkeypatch.setattr(gn.urllib.request, "urlopen", fake_urlopen)
    out = gn.fetch_ncbi_gene_names(["726817", "409662", "not-an-id"])
    assert out["409662"] == ("alc", "5'-AMP-activated protein kinase subunit beta-1")
    assert out["726817"][0] == "LOC726817"
    assert seen["url"].startswith(gn.NCBI_GENE_URL + "409662,726817")
    assert seen["timeout"] == gn.NCBI_TIMEOUT_S


def test_an_unreachable_ncbi_raises_oserror(monkeypatch):
    def fake_urlopen(req, timeout=None):
        raise urllib.error.URLError("no route to host")

    monkeypatch.setattr(gn.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(OSError, match="NCBI Gene lookup failed"):
        gn.fetch_ncbi_gene_names(["726817"])


# ---- the names table file --------------------------------------------------

def test_names_table_round_trip(tmp_path, gff):
    table, note = gn.resolve_gene_names(IDS, gn.load_gff_gene_info(gff), allow_network=False)
    path = tmp_path / "gene_names.tsv"
    gn.write_names_tsv(str(path), table, IDS, lookup_note=note)
    text = path.read_text()
    assert text.splitlines()[-len(IDS) - 1].split("\t") == list(gn.TSV_COLUMNS)
    assert "# NCBI Gene lookup: skipped" in text
    assert gn.read_names_tsv(str(path)) == table
    # the saved table reproduces the labels without any lookup
    again, _ = gn.resolve_gene_names(IDS, {}, overrides=gn.read_names_tsv(str(path)),
                                     fetch=_no_fetch)
    assert again == table


def test_a_hand_written_two_column_table_is_enough(tmp_path):
    path = tmp_path / "names.tsv"
    path.write_text("# my names\ngene_id\tlabel\ngene-LOC123456\tApy\n\nbroken line\n")
    assert gn.read_names_tsv(str(path)) == {
        "gene-LOC123456": {"label": "Apy", "full_name": "", "source": gn.SOURCE_USER,
                           "ncbi_gene_id": ""}}
    assert gn.read_names_tsv(str(tmp_path / "missing.tsv")) == {}


# ---- the table under the figures -------------------------------------------

def _width(text, size, bold=False):
    return len(text or "") * size * 0.55


def _rows():
    named = {"label": "PTL", "full_name": "pancreatic triacylglycerol lipase",
             "source": gn.SOURCE_PRODUCT}
    official = {"label": "Melt", "full_name": "melittin", "source": gn.SOURCE_GFF}
    unnamed = {"label": "LOC100578368", "full_name": "", "source": gn.SOURCE_ID}
    return [gn.table_row(named, "LOC726817", "#4e79a7"),
            gn.table_row(official, "Melt", "#e31a1c", is_goi=True),
            gn.table_row(unnamed, "LOC100578368", "#59a14f")]


def test_table_lists_label_name_id_and_source():
    svg, height = gn.names_table_svg(_rows(), 24, 500, 1200, _width)
    text = "\n".join(svg)
    assert height > 0 and 'class="gene-names-table"' in text
    for expected in ("Gene names", ">PTL*<", "pancreatic triacylglycerol lipase", ">LOC726817<",
                     "product name", ">Melt<", "home GFF", "no name annotated"):
        assert expected in text, expected
    assert gn.DERIVED_NOTE.replace("'", "&#x27;") in text or gn.DERIVED_NOTE in text


def test_the_mark_is_explained_only_when_it_is_used():
    svg, _ = gn.names_table_svg(_rows()[1:], 24, 500, 1200, _width)
    assert "not an official gene symbol" not in "\n".join(svg)


def test_a_table_of_bare_ids_is_not_drawn():
    assert gn.names_table_svg(_rows()[2:], 24, 500, 1200, _width) == ([], 0)
    assert gn.names_table_svg([], 24, 500, 1200, _width) == ([], 0)


def test_a_skipped_or_failed_lookup_is_written_under_the_table():
    failed = "\n".join(gn.names_table_svg(_rows(), 24, 500, 1200, _width,
                                          lookup_note="FAILED (timed out)")[0])
    assert "NCBI Gene was not queried (not reachable)" in failed
    ok = "\n".join(gn.names_table_svg(_rows(), 24, 500, 1200, _width,
                                      lookup_note="ok (3 of 3 records returned)")[0])
    assert "NCBI Gene was not queried" not in ok


def test_one_entry_always_fits_the_canvas():
    rows = _rows()
    rows[0]["full"] = "a very long product name " * 12
    svg, _ = gn.names_table_svg(rows, 24, 500, 420, _width)
    xs = [float(x) for x in re.findall(r'<text x="([0-9.]+)"', "\n".join(svg))]
    assert max(xs) < 24 + 420


# ---- plot_synteny.py -------------------------------------------------------

NAMES = {
    "gene-AAA": {"label": "PTL", "full_name": "pancreatic triacylglycerol lipase",
                 "source": gn.SOURCE_PRODUCT, "ncbi_gene_id": "1"},
    "gene-BBB": {"label": "alc", "full_name": "alicorn", "source": gn.SOURCE_NCBI,
                 "ncbi_gene_id": "2"},
    "gene-CCC": {"label": "LOC3", "full_name": "", "source": gn.SOURCE_ID, "ncbi_gene_id": "3"},
    "gene-Melt": {"label": "Melt", "full_name": "melittin", "source": gn.SOURCE_GFF,
                  "ncbi_gene_id": "4"},
}


@pytest.fixture
def names(monkeypatch):
    monkeypatch.setattr(ps, "_GENE_NAMES", dict(NAMES))
    monkeypatch.setattr(ps, "_GENE_NAMES_NOTE", ["ok"])
    monkeypatch.setattr(ps, "_GOI_NAMES", {"gene-Melt"})
    return NAMES


def _g(name, start, end, **kw):
    g = {"chrom": "chr1", "start": start, "end": end, "start_plot": start, "end_plot": end,
         "name": name, "home_gene_id": name, "strand": "+", "identity": 80.0,
         "exon_coords": [], "confidence": "HIGH"}
    g.update(kw)
    return g


def _tracks():
    home = {"label": "Apis mellifera", "is_home": True, "genome_id": "home", "offset": 0,
            "breaks": [], "genes": [_g("gene-AAA", 1000, 2000), _g("gene-BBB", 3000, 4000),
                                    _g("gene-Melt", 5000, 6000), _g("gene-CCC", 7000, 8000)]}
    target = {"label": "Bombus terrestris", "is_home": False, "genome_id": "bter", "offset": 0,
              "breaks": [], "genes": [
                  _g("gene-AAA|bter_fl1", 100, 200, home_gene_id="gene-AAA"),
                  _g("gene-BBB|bter_fl2", 300, 400, home_gene_id="gene-BBB"),
                  _g("GOI_Melt|bter_l1", 500, 600, home_gene_id="GOI_Melt", role="goi"),
                  _g("gene-CCC|bter_fl3", 700, 800, home_gene_id="gene-CCC")]}
    return [home, target]


COLOURS = {"gene-AAA": "#4e79a7", "gene-BBB": "#b07aa1", "gene-CCC": "#59a14f"}


def _texts(html):
    return re.findall(r"<text[^>]*>([^<]*)</text>", html)


def test_home_gene_label_accepts_every_spelling_of_the_id(names):
    assert ps._home_gene_label("gene-AAA") == "PTL*"
    assert ps._home_gene_label("AAA") == "PTL*"
    assert ps._home_gene_label("gene-AAA|bter_fl1") == "PTL*"
    assert ps._home_gene_label("GOI_Melt") == "Melt"
    assert ps._home_gene_label("gene-ZZZ") == "ZZZ"            # not in the table: bare ID
    assert ps.clean_gene_label("GOI_Melt", keep_goi_prefix=True) == "GOI Melt"
    assert ps.clean_gene_label("GOI_copy_2", keep_goi_prefix=True) == "GOI #2"


@pytest.mark.parametrize("render", [ps.render_anchor_grid, ps.render_anchor_grid_threaded])
def test_grid_columns_carry_names_not_ids(names, render):
    html = render(_tracks(), COLOURS, {}, {}, SimpleNamespace())
    # the figure itself, without the names table (which lists the IDs on purpose)
    texts = _texts(html.split('class="gene-names-table"')[0])
    for label in ("PTL*", "alc", "Melt", "LOC3"):
        assert label in texts, label
    assert "AAA" not in texts and "BBB" not in texts
    assert any(t.endswith("GOI: Melt") for t in texts)


@pytest.mark.parametrize("render", [ps.render_anchor_grid, ps.render_anchor_grid_positional,
                                    ps.render_anchor_grid_threaded])
def test_every_grid_figure_has_the_names_table(names, render):
    html = render(_tracks(), COLOURS, {}, {}, SimpleNamespace())
    assert html.count('class="gene-names-table"') == 1
    texts = _texts(html)
    for expected in ("Gene names", "pancreatic triacylglycerol lipase", "alicorn", "melittin",
                     "AAA", "NCBI Gene", "product name", "home GFF", "no name annotated"):
        assert expected in texts, expected
    assert any("not an official gene symbol" in t for t in texts)
    # the canvas grew by the table: nothing is drawn below its lower edge
    height = float(re.search(r'<svg class="grid-svg" width="[0-9.]+" height="([0-9.]+)"',
                             html).group(1))
    table = html[html.index('class="gene-names-table"'):]
    assert max(float(y) for y in re.findall(r' y="([0-9.]+)"', table)) < height


@pytest.mark.parametrize("render", [ps.render_anchor_grid, ps.render_anchor_grid_positional,
                                    ps.render_anchor_grid_threaded])
def test_no_gene_legend_switch_removes_the_table(names, render):
    html = render(_tracks(), COLOURS, {}, {}, SimpleNamespace(gene_legend=False))
    assert "gene-names-table" not in html
    assert "PTL*" in _texts(html) or render is ps.render_anchor_grid_positional


def test_without_a_names_table_the_figure_is_drawn_as_before(monkeypatch):
    monkeypatch.setattr(ps, "_GENE_NAMES", {})
    html = ps.render_anchor_grid(_tracks(), COLOURS, {}, {}, SimpleNamespace())
    assert "gene-names-table" not in html
    assert "AAA" in _texts(html)


def test_target_gene_keeps_its_own_annotation_name(names):
    own_symbol = _g("x", 1, 2, home_gene_id="gene-AAA", target_gene="Lip3",
                    target_product="lipase 3")
    assert ps._preferred_target_label(own_symbol) == "Lip3"
    own_product = _g("x", 1, 2, home_gene_id="gene-AAA", target_gene="LOC6732354",
                     target_product="heterogeneous nuclear ribonucleoprotein A1")
    assert ps._preferred_target_label(own_product) == "HNRA1*"
    # no name of its own: the label of the home gene it is the ortholog of
    unnamed = _g("gene-BBB|bter_fl2", 1, 2, home_gene_id="gene-BBB",
                 target_gene="LOC6732354", target_product="uncharacterized LOC6732354")
    assert ps._preferred_target_label(unnamed) == "alc"


def test_tooltip_gives_label_full_name_and_id(names):
    home, target = _tracks()
    tip = json.loads(re.sub(r"&quot;", '"', ps._build_tooltip_json(home["genes"][0], home, {})))
    assert (tip["name"], tip["product"], tip["geneId"]) == \
        ("PTL*", "pancreatic triacylglycerol lipase", "AAA")
    tip = json.loads(re.sub(r"&quot;", '"', ps._build_tooltip_json(target["genes"][1], target, {})))
    assert tip["homolog"] == "alc (BBB)"


def _ribbon_args(**kw):
    base = dict(plot_width=0, plot_height=0, scale_bar_len=10000, goi_zoom=1.0, goi_min_px=0,
                caption_file=None, ribbon_alpha_dense=0.2, max_legend_entries=25)
    base.update(kw)
    return SimpleNamespace(**base)


def _ribbon_tracks():
    tracks = _tracks()
    for t in tracks:
        t.update({"chrom": "chr1", "species": t["label"], "gap_breaks": []})
        for g in t["genes"]:
            g.setdefault("_sub_track", 0)
    return tracks


def test_interactive_ribbon_plot_has_no_table(names):
    # Label, full name and gene ID are in the tooltip; the table is for print only.
    html = ps.render_synteny_html(_ribbon_tracks(), COLOURS, {}, {}, _ribbon_args(), [], 0, 0, 1)
    assert 'class="gene-names-table"' not in html
    # The footer holds the legend only; it follows the canvas bottom on reflow.
    assert html.count('class="figure-legend"') == 1
    layout = json.loads(re.search(r"__SYNVOY_LAYOUT__ = (\{.*?\});</script>", html).group(1))
    assert 0 < layout["footerH"] < 200
    tooltips = " ".join(re.findall(r"data-tooltip='([^']*)'", html))
    assert "pancreatic triacylglycerol lipase" in tooltips


def test_publication_ribbon_plot_has_the_table_and_keeps_it_on_reflow(names):
    html = ps.render_synteny_html(_ribbon_tracks(), COLOURS, {}, {}, _ribbon_args(), [], 0, 0, 1,
                                  force_home_labels=True)
    assert html.count('class="gene-names-table"') == 1
    assert "pancreatic triacylglycerol lipase" in _texts(html)
    layout = json.loads(re.search(r"__SYNVOY_LAYOUT__ = (\{.*?\});</script>", html).group(1))
    assert layout["footerH"] > 0
    assert ".plot-caption, .gene-names-table" in html      # the reflow script moves it
    svg_h = float(re.search(r'<svg class="synteny-svg" width="\d+" height="([0-9.]+)"',
                            html).group(1))
    table = html[html.index('class="gene-names-table"'):html.index("window.__SYNVOY", 0)
                 if html.index("window.__SYNVOY") > html.index('class="gene-names-table"')
                 else len(html)]
    assert max(float(y) for y in re.findall(r' y="([0-9.]+)"', table)) < svg_h

    off = ps.render_synteny_html(_ribbon_tracks(), COLOURS, {}, {},
                                 _ribbon_args(gene_legend=False), [], 0, 0, 1)
    assert "gene-names-table" not in off.split("<script>")[0].split("</style>")[-1]


def test_load_gene_names_end_to_end(tmp_path, gff, monkeypatch, capsys):
    monkeypatch.setattr(ps, "_GENE_NAMES", {})
    monkeypatch.setattr(ps, "_GENE_NAMES_NOTE", [""])
    home_genes = [
        {"name": "gene-Melt", "start": 100, "end": 900, "display_name": "Melt"},
        # the BED's display column holds a transcript accession for this gene
        {"name": "gene-LOC726817", "start": 1000, "end": 1900, "display_name": "XM_006567126.3"},
        {"name": "gene-LOC100578368", "start": 2000, "end": 2900, "display_name": ""},
        # not in the GFF (a predicted gene): the BED's own name is all there is
        {"name": "g42", "start": 5000, "end": 5900, "display_name": "arginine kinase"},
        {"name": "GOI_chr1_100", "start": 100, "end": 200, "display_name": ""},
    ]
    out = tmp_path / "locus_1_synteny_plot.html"
    args = SimpleNamespace(home_gff=gff, gene_names_tsv="", gene_names_out="",
                           no_network=True, gene_name_lookup=False, output=str(out))
    ps._load_gene_names(args, home_genes)
    assert ps._home_gene_label("gene-LOC726817") == "PTL*"      # not XM_006567126.3
    assert ps._home_gene_label("gene-LOC100578368") == "LOC100578368"
    assert ps._home_gene_label("g42") == "ARGK*"
    assert "GOI_chr1_100" not in ps._GENE_NAMES
    saved = tmp_path / "locus_1_gene_names.tsv"
    assert saved.exists() and "gene-LOC726817\tPTL\t" in saved.read_text()
    assert "NCBI Gene lookup: skipped" in capsys.readouterr().out

    # an edited table decides the label, whichever way it spells the gene ID
    edited = tmp_path / "edited.tsv"
    edited.write_text("gene_id\tlabel\tfull_name\tsource\nLOC726817\tPtl\tlipase\tuser\n")
    args.gene_names_tsv = str(edited)
    args.gene_names_out = str(tmp_path / "second.tsv")
    ps._load_gene_names(args, home_genes)
    assert ps._home_gene_label("gene-LOC726817") == "Ptl"
    assert (tmp_path / "second.tsv").exists()


def test_goi_title_label_comes_from_the_annotated_home_gene(names):
    assert ps._goi_display_label(_tracks()) == "Melt"


def test_long_column_labels_get_a_taller_header():
    short = [{"label": "alc", "is_goi": False}]
    long_ = [{"label": "lncRNA:CR44733-extra-long", "is_goi": False}]
    assert ps._rotated_header_height(long_, 17, 55) > 150 > ps._rotated_header_height(short, 17, 55)


# ---- plot_synteny_matrix.py ------------------------------------------------

def test_matrix_slots_and_table_use_the_names(gff):
    home_bed = [{"chrom": "chr1", "start": 999, "end": 1900, "name": "gene-LOC726817",
                 "strand": "+", "display_name": "XM_006567126.3"},
                {"chrom": "chr1", "start": 2999, "end": 3900, "name": "gene-Dmel_CG4491",
                 "strand": "+", "display_name": "noc"}]
    query_bed = [{"chrom": "chr1", "start": 150, "end": 400, "name": "q", "strand": "+",
                  "display_name": ""}]
    args = SimpleNamespace(home_gff=gff, gene_names_tsv="", no_network=True)
    index = pm.parse_home_gff_genes(gff)
    table, note = pm.load_gene_names(args, home_bed, query_bed, index)
    slots, goi_id = pm.build_home_slots(home_bed, query_bed, index, table)
    assert [s["label"] for s in slots] == ["Melt", "PTL*", "noc"]
    assert goi_id == "GOI_Melt"                        # the join key keeps the GFF symbol
    rows = [("Apis mellifera", pm._home_summary(slots, goi_id), True, "Apis_mellifera", None)]
    svg = pm.render_svg(slots, rows, "Apis mellifera", goi_id, names=table, lookup_note=note)
    assert svg.count('class="gene-names-table"') == 1
    assert "pancreatic triacylglycerol lipase" in svg and "NCBI Gene was not queried" in svg
    assert "gene-names-table" not in pm.render_svg(slots, rows, "Apis mellifera", goi_id)
