"""Fix-phase item 19: a GOI whose name starts with a digit keeps its records.

extract_id's locus_tag pattern ^([A-Z]{2,6}_\\d+) has no boundary after the digits and GOI matches
[A-Z]{2,6}, so GOI_14-3-3epsilon and GOI_14-3-3epsilon|exon_1 both became GOI_14. PREPARE_INITIAL_DB then saw
no exon, counted every record as a full GOI and built fallback fragments from the last one, an exon. Read back,
those fragments carried that id too, so the search kept the last fragment as the GOI's home protein (50 of 180 aa
on a real gene). The new first pattern takes the whole first word of a GOI_ header and leaves every other header to
the old patterns. The figures' name of the GOI took a letter first only, so such a GOI read "GOI" (or a relative's
name) where GOI_Ly6e read Ly6e.
"""
import os
import random
import re
import string
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import build_home_paralog_panel as panel  # noqa: E402
import iterative_search_runner as runner  # noqa: E402
import plot_synteny as ps  # noqa: E402
import sequence_utils as su  # noqa: E402

DIGIT_FIRST = ["14", "14-3-3epsilon", "1700001C19Rik", "1_venom", "14.1", "5-HT2A", "4E-BP1", "18S", "2_1"]
FULL = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRL"


def old_extract_id(header):
    """extract_id as it was before the fix: the same function without the new pattern."""
    keep = su.ID_PATTERNS
    su.ID_PATTERNS = [p for p in keep if p[1] != "goi"]
    try:
        return su.extract_id(header)
    finally:
        su.ID_PATTERNS = keep


@pytest.mark.parametrize("header,expected", [
    ("GOI_14-3-3epsilon|exon_1", "GOI_14-3-3epsilon|exon_1"),
    ("GOI_14-3-3epsilon", "GOI_14-3-3epsilon"),
    ("GOI_1700001C19Rik|exon_2", "GOI_1700001C19Rik|exon_2"),
    ("GOI_1_venom|exon_1", "GOI_1_venom|exon_1"),
    ("GOI_14.1", "GOI_14.1"),
    ("GOI_Ly6e|exon_1", "GOI_Ly6e|exon_1"),
    ("GOI_14-3-3epsilon|exon_1 a description", "GOI_14-3-3epsilon|exon_1"),
])
def test_a_goi_record_keeps_its_whole_first_word(header, expected):
    assert su.extract_id(header) == expected


@pytest.mark.parametrize("name", DIGIT_FIRST)
def test_the_old_patterns_cut_these_names(name):
    """What the fix is for: the same headers without the new pattern lose their exon number."""
    cut = old_extract_id(f"GOI_{name}|exon_1")
    assert cut != f"GOI_{name}|exon_1" and re.fullmatch(r"GOI_\d+", cut)
    assert su.extract_id(f"GOI_{name}|exon_1") == f"GOI_{name}|exon_1"


@pytest.mark.parametrize("header", [
    "GOI_Ly6e", "GOI_Ly6e|exon_1", "GOI_copy_3", "GOI_copy_1|Apis_florea_fna_b0_l1_by_GOI_Melt",
    "GOI_rescue_Apis_florea_locus_block1_2", "GOI_Melt|Apis_florea_fna_b0_l1_exon_ann",
    "GOI_Melt|Colletes_gigas_WUUM01000001.1_hull_rescue", "GOI_NEIGHBOR_gene-ABC1", "GOI_lcl|ORF1",
    "GOI_lcl|NC_123.1_cds_XP_456.1_7", "GOI_12345", "GOI_12345 a description", "GOI_x~shuf1~rel",
])
def test_every_other_goi_header_gives_what_it_gave(header):
    assert su.extract_id(header) == old_extract_id(header)


def test_no_other_goi_name_changes():
    """Any GOI_ header whose name does not start with a digit, or is only digits, as before."""
    rng = random.Random(19)
    letters = string.ascii_letters + "._-|~"
    rest = string.ascii_letters + string.digits + "._-|~"
    for _ in range(5000):
        name = rng.choice(letters) + "".join(rng.choice(rest) for _ in range(rng.randint(0, 20)))
        header = f"GOI_{name}" + rng.choice(["", " a description"])
        assert su.extract_id(header) == old_extract_id(header), header
    for n in (1, 14, 12345):
        assert su.extract_id(f"GOI_{n}") == old_extract_id(f"GOI_{n}") == f"GOI_{n}"


def _fixture_headers():
    """Headers of the repository's own material: the tracked example query, the FASTA text in the tests'
    sources, one header per kind of id the patterns know."""
    out = [line[1:].rstrip("\n") for line in open(ROOT / "assets" / "example_queries" / "P60615.faa")
           if line.startswith(">")]
    for path in (ROOT / "tests").glob("*.py"):
        if path.name == Path(__file__).name:
            continue
        out += re.findall(r"""(?:^|["'n])>([A-Za-z0-9_|.\-][^"'\\\n]*)""", path.read_text(errors="replace"), flags=re.M)
    out += ["XP_123456.1 a protein", "NP_123456.1", "AAA12345.1 hypothetical protein", "gi|12345|ref|XP_1.1|",
            "NC_123456.1 chromosome 1", "NZ_CP123456.1", "lcl|NC_123456.1_cds_XP_123456.1_1 [gene=x]", "P12345", "Q9NXB0",
            "A0A1B2C3D4", "sp|P12345|NAME_HUMAN desc", "tr|A0A123|NAME", "ENSP00000123456.2", "ENSMUSG00000000001",
            "ENST00000123456", "FBgn0000001", "FBpp0000001", "WBGene00000001", "AT1G01010.1", "1_234 # 1000 # 2000 # 1 # ID=1_234",
            "ABC_12345", "ABC_12345.1", "ABC_12345 locus tag style", "SUPER_10", "SUPER_10_unloc_1", "gene_12345", "gene_12-like",
            "protein_12", "prot_7-x", "gene-LOC123456789|exon_1", "LOC123456789 desc", "random_header_with no known pattern",
            "14-3-3epsilon", "1_venom", "1700001C19Rik", "goi_14-3-3", "GOI", "GOI-14", "xGOI_14-3-3", "", "   "]
    return out


@pytest.mark.parametrize("header,expected", [
    ("XP_123456.1 a protein", "XP_123456.1"), ("AAA12345.1 hypothetical protein", "AAA12345.1"),
    ("gi|12345|ref|XP_1.1|", "12345"), ("NC_123456.1 chromosome 1", "NC_123456.1"), ("sp|P12345|NAME_HUMAN desc", "P12345"),
    ("ENSP00000123456.2", "ENSP00000123456.2"), ("FBgn0000001", "FBgn0000001"), ("WBGene00000001", "WBGene00000001"),
    ("AT1G01010.1", "AT1G01010.1"), ("1_234 # 1000 # 2000 # 1 # ID=1_234", "1_234"), ("ABC_12345 locus tag style", "ABC_12345"),
    ("gene_12345", "gene_12345"), ("protein_12", "protein_12"), ("gene-LOC123456789|exon_1", "gene-LOC123456789|exon_1"),
    ("14-3-3epsilon", "14-3-3epsilon"), ("1_venom", "1_venom"),
])
def test_other_headers_keep_their_ids(header, expected):
    """The ids as they stood before the fix, written out: the patterns behind the new one are the old ones."""
    assert su.extract_id(header) == expected


def test_the_new_pattern_fires_on_no_other_header():
    headers = _fixture_headers()
    assert len(headers) > 100
    for h in headers:
        if not h.strip().startswith("GOI_") or not h.strip()[4:5].isdigit():
            assert su.extract_id(h) == old_extract_id(h), h
    rng = random.Random(1919)
    alphabet = string.ascii_letters + string.digits + "._-|~ #"
    n = 0
    while n < 20000:
        h = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 24)))
        if h.strip().startswith("GOI_"):
            continue
        n += 1
        assert su.extract_id(h) == old_extract_id(h), h


def _prepare_initial_db(tmp_path, goi_exons, flanking_faa):
    """The module's own script, with Nextflow's substitutions done by hand."""
    nf = (ROOT / "modules" / "prepare_initial_db.nf").read_text()
    body = nf.split('"""')[1]
    for key, value in (("${projectDir}", str(ROOT)), ("${flanking_faa}", str(flanking_faa)),
                       ("${goi_exons}", str(goi_exons)), ("${locus_id}", "locus_1")):
        body = body.replace(key, value)
    script = tmp_path / "prepare.py"
    script.write_text(textwrap.dedent(body).replace("\\\\", "\\"))
    r = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    return r.stdout, tmp_path / "initial_db_locus_1.faa"


def _goi_records(name):
    return [(f"GOI_{name}", FULL), (f"GOI_{name}|exon_1", FULL[:40]), (f"GOI_{name}|exon_2", FULL[40:70]),
            (f"GOI_{name}|exon_3", FULL[70:])]


@pytest.mark.parametrize("name", DIGIT_FIRST + ["Ly6e"])
def test_prepare_initial_db_sees_the_exons(tmp_path, name):
    goi = tmp_path / "goi_exons.faa"
    su.write_fasta(_goi_records(name), str(goi))
    flank = tmp_path / "flank.faa"
    flank.write_text(">gene-A|exon_1 home=1-100 coords=1-100\nMKVLAAGIVGLLLAQWERTY\n")
    out, db = _prepare_initial_db(tmp_path, goi, flank)
    assert "- GOI (full): 1" in out and "- GOI exons: 3" in out and "- Fallback fragments: 0" in out
    assert f"Added full GOI: GOI_{name} ({len(FULL)} aa)" in out
    assert "generating fallback fragments" not in out
    ids = [cid for _h, cid, _s in su.parse_fasta(str(db))]
    assert len(ids) == len(set(ids)), ids
    assert {f"GOI_{name}", f"GOI_{name}|exon_1", f"GOI_{name}|exon_2", f"GOI_{name}|exon_3"} <= set(ids)


@pytest.mark.parametrize("name", DIGIT_FIRST)
def test_the_search_and_the_panel_take_the_full_protein(tmp_path, name):
    goi = tmp_path / "goi_exons.faa"
    su.write_fasta(_goi_records(name), str(goi))
    seqs = {cid: seq for _h, cid, seq in su.parse_fasta(str(goi))}      # as the runner's base_db_sequences
    assert seqs[f"GOI_{name}"] == FULL
    assert len(seqs) == 4
    assert runner.is_full_length_goi_query_id(f"GOI_{name}")
    assert not runner.is_full_length_goi_query_id(f"GOI_{name}|exon_1")
    assert runner.extract_base_gene_id(f"GOI_{name}|exon_1") == f"GOI_{name}"
    assert runner.build_home_goi_ids(str(goi)) == set(seqs)
    _queries, by_parent, _aliases = runner.select_goi_queries({k: {"id": k, "seq": v} for k, v in seqs.items()},
                                                              runner.build_home_goi_ids(str(goi)))
    assert by_parent[f"GOI_{name}"]["seq"] == FULL
    assert panel.home_goi_record(str(goi)) == (f"GOI_{name}", FULL)


def _target_model(name):
    return {"home_gene_id": f"GOI_{name}", "name": f"GOI_{name}|Apis_florea_fna_b0_l1_exon_ann"}


@pytest.mark.parametrize("name", ["5-HT2A", "4E-BP1", "18S", "14", "Ly6e"])
def test_the_figures_take_a_digit_first_name_from_the_models(name):
    """Without a named home gene (a run without a home GFF) the title and the GOI column read the models' name."""
    assert ps._goi_display_label([{"genes": [_target_model(name)]}]) == name


@pytest.mark.parametrize("placeholder", ["GOI_1_4100", "GOI_2L_4100", "GOI_NC_000012.12_91140483", "GOI_copy_3"])
def test_a_home_placeholder_or_a_copy_is_no_figure_name(placeholder):
    assert ps._goi_display_label([{"genes": [{"home_gene_id": placeholder, "name": placeholder}]}]) == "GOI"


@pytest.mark.parametrize("name", ["14.1", "1_venom", "14-3-3epsilon", "1700001C19Rik", "Ly6e.1", "Ly6e_1", "Serpina3n-long"])
def test_a_name_with_a_dot_an_underscore_or_over_12_characters_is_no_figure_name(name):
    """The rule's limit, first character aside: such a name reads "GOI" without a named home gene."""
    assert ps._goi_display_label([{"genes": [_target_model(name)]}]) == "GOI"


@pytest.mark.parametrize("goi,relative", [("5-HT2A", "5-HT2B"), ("Ly6e", "Ly6d")])
def test_of_two_home_genes_the_figures_name_the_one_the_models_carry(goi, relative, monkeypatch):
    names = {f"gene-{g}": {"label": g, "full": "", "id": f"gene-{g}", "source": "gff_symbol"} for g in (goi, relative)}
    monkeypatch.setattr(ps, "_GENE_NAMES", names)
    monkeypatch.setattr(ps, "_GOI_NAMES", set(names))
    home = {"genes": [{"name": f"gene-{relative}"}, {"name": f"gene-{goi}"}]}
    assert ps._goi_display_label([home, {"genes": [_target_model(goi)]}]) == goi


CODON = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG',
         'Y': 'TAT', 'V': 'GTT'}
INTRON = 'GT' + 'ACGTTGCA' * 30 + 'AG'
SPACER = 'TAACTAGTTGA' * 40


@pytest.mark.parametrize("gene,query_id", [("14-3-3epsilon", "14-3-3epsilon"), ("1700001C19Rik", "some_query")])
def test_a_digit_first_gene_of_the_home_gff_reaches_the_query_set_whole(tmp_path, gene, query_id):
    """With a home GFF the name is the gene's in the GFF, whatever the query's id: ANNOTATE_GOI end to end, then
    PREPARE_INITIAL_DB on its output."""
    pytest.importorskip("parasail")
    e1, e2 = FULL[:50], FULL[50:]
    x1 = (len(SPACER), len(SPACER) + 3 * len(e1))
    x2 = (x1[1] + len(INTRON), x1[1] + len(INTRON) + 3 * len(e2) + 3)
    dna = SPACER + ''.join(CODON[a] for a in e1) + INTRON + ''.join(CODON[a] for a in e2) + 'TAA' + SPACER
    (tmp_path / "home.fna").write_text(f">chrA\n{dna}\n")
    (tmp_path / "q.faa").write_text(f">{query_id}\n{FULL}\n")

    def row(kind, span, attrs, phase='.'):
        return f"chrA\tsrc\t{kind}\t{span[0] + 1}\t{span[1]}\t.\t+\t{phase}\t{attrs}\n"
    (tmp_path / "home.gff").write_text(
        "##gff-version 3\n" + row('gene', (x1[0], x2[1]), f"ID=gene-{gene};Name={gene};gene={gene}")
        + row('mRNA', (x1[0], x2[1]), f"ID=rna-1;Parent=gene-{gene}")
        + row('CDS', x1, "ID=cds-1;Parent=rna-1", '0') + row('CDS', x2, "ID=cds-1;Parent=rna-1", '0'))
    (tmp_path / "hits.m8").write_text(f"{query_id}\tchrA\t100\t50\t0\t0\t1\t50\t{x1[0] + 1}\t{x1[1]}\t1e-20\t100\n")
    out = {k: tmp_path / f"out.{k}" for k in ("faa", "bed", "json")}
    r = subprocess.run([sys.executable, str(ROOT / "bin" / "annotate_goi_exons.py"), "--query", str(tmp_path / "q.faa"),
                        "--genome", str(tmp_path / "home.fna"), "--gff", str(tmp_path / "home.gff"),
                        "--blast_hits", str(tmp_path / "hits.m8"), "--output_exons", str(out["faa"]),
                        "--output_bed", str(out["bed"]), "--output_info", str(out["json"])],
                       capture_output=True, text=True, env=dict(os.environ))
    assert r.returncode == 0, r.stderr
    ids = [cid for _h, cid, _s in su.parse_fasta(str(out["faa"]))]
    assert ids == [f"GOI_{gene}", f"GOI_{gene}|exon_1", f"GOI_{gene}|exon_2"]
    flank = tmp_path / "flank.faa"
    flank.write_text(">gene-A|exon_1 home=1-100 coords=1-100\nMKVLAAGIVGLLLAQWERTY\n")
    stdout, db = _prepare_initial_db(tmp_path, out["faa"], flank)
    assert "- GOI (full): 1" in stdout and "- GOI exons: 2" in stdout and "- Fallback fragments: 0" in stdout
    assert {cid: seq for _h, cid, seq in su.parse_fasta(str(db))}[f"GOI_{gene}"] == FULL
