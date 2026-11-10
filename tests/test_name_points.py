"""Fix items 28 and 34 (9 Nov 2026): smaller points about names. None of them moves a gate set: each needs a name or an
input no gate set has.

- 28a: a GFF3-reserved character (, ; = & %) in the GOI's name: SynVoy writes it into the ID, Parent and SynVoy_Parent
  of every model, where ',' makes a Parent a list and ';' ends the attribute. normalize_query.py and ANNOTATE_GOI stop.
- 28b: a GOI named copy_<n> read as a tandem copy (ANNOTATE_GOI stops), a tandem copy's bare number made flanking genes
  '1' or 'gene-1' GOI proxies, and the ribbon plot read any '<x>_copy_<n>' as 'GOI #n'.
- 28c: ANNOTATE_GOI's --query_id was unquoted and --gene_name in plain single quotes (lib/ShellArgs.groovy now), and
  UniProt's symbol search was unquoted and took its first hit as it came (Su(H) gave Su(dx)).
- 28d: the structural discovery dropped every GOI template whose name holds 'exon' (Hexon).
- 28e: the ncbi_cds pattern could never match, so an NCBI CDS download's header (lcl|..._cds_XP_..._1) kept its first word;
  a query id with a database tag (ref|X|, gb|X|, ENA|X|X.1, lcl|X) kept its '|', and the GOI's name read as the tag.
  The common forms read to their accession; any '|' left stops normalize_query.py and ANNOTATE_GOI.
- 28f: LOC123456 (NCBI Gene's provisional symbol) read as a protein accession; efetch's error text was written as FASTA.
- 34: figures stripped only NCBI's 'gene-' prefix (not Ensembl's 'gene:'), and the ribbon plot cut a scaffold label to
  its last 12 characters.
"""
import json
import os
import shutil
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import sequence_utils as su  # noqa: E402

NEXTFLOW = shutil.which("nextflow")
PROT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV"


# --------------------------------------------------------------------------- 28a

@pytest.mark.parametrize("name,chars", [("my,gene", "','"), ("my;gene", "';'"), ("a=b", "'='"), ("a&b", "'&'"),
                                        ("50%", "'%'")])
def test_a_gff_reserved_character_is_named(name, chars):
    assert chars in su.id_problem(name, whitespace=False, gff=True)
    assert su.id_problem(name, whitespace=False) is None          # only where it goes into GFF column 9


@pytest.mark.parametrize("qid", ["my,gene", "my;gene", "a=b"])
def test_normalize_query_stops_on_a_reserved_character(tmp_path, qid):
    (tmp_path / "q.faa").write_text(f">{qid} a protein\n{PROT}\n")
    p = subprocess.run([sys.executable, os.path.join(BIN, "normalize_query.py"), "--input", str(tmp_path / "q.faa"),
                        "--output", str(tmp_path / "out.faa")], capture_output=True, text=True)
    assert p.returncode == 1 and qid in p.stderr and "reserved in GFF3" in p.stderr, p.stderr


def _annotate(tmp_path, query_id):
    (tmp_path / "q.faa").write_text(f">q\n{PROT}\n")
    (tmp_path / "h.fna").write_text(">c1\n" + "ACGT" * 50 + "\n")
    return subprocess.run([sys.executable, os.path.join(BIN, "annotate_goi_exons.py"), "--query", str(tmp_path / "q.faa"),
                           "--genome", str(tmp_path / "h.fna"), "--gff", "NO_GFF", "--query_id", query_id,
                           "--output_exons", str(tmp_path / "e.faa"), "--output_bed", str(tmp_path / "e.bed"),
                           "--output_info", str(tmp_path / "i.json")], capture_output=True, text=True, cwd=tmp_path)


@pytest.mark.parametrize("query_id,shown", [("a,b", "reserved in GFF3"), ("a b", "whitespace"),
                                            ("copy_1", "tandem copy"), ("ref|XP_1.1|", "'|'")])
def test_annotate_goi_stops_on_a_goi_name_the_search_cannot_carry(tmp_path, query_id, shown):
    """--query_id (easy mode) reaches ANNOTATE_GOI without normalize_query's check of the header."""
    if shutil.which("miniprot") is None:
        pytest.skip("miniprot not on PATH")
    p = _annotate(tmp_path, query_id)
    assert p.returncode == 3 and shown in p.stderr and repr(query_id) in p.stderr, p.stderr
    assert "Traceback" not in p.stderr


# --------------------------------------------------------------------------- 28b

def test_a_tandem_copy_gives_no_goi_proxy():
    import iterative_search_runner as isr
    ids = {"GOI_copy_1|G_fna_b0_l1_by_GOI_x", "GOI_copy_1", "GOI_x", "gene-1", "1", "gene-x", "gene-y"}
    assert isr.goi_proxy_parents(ids) == {"gene-x"}


def test_only_a_tandem_copy_reads_goi_number():
    import plot_synteny as ps
    assert ps.clean_gene_label("GOI_copy_2", keep_goi_prefix=True) == "GOI #2"
    assert ps.clean_gene_label("GOI_Mel_copy_2", keep_goi_prefix=True) != "GOI #2"


# --------------------------------------------------------------------------- 28c

def test_annotate_goi_quotes_its_name_arguments():
    nf = open(os.path.join(ROOT, "modules", "annotate_goi.nf")).read()
    assert '"--query_id ${ShellArgs.q(query_id)}"' in nf and '"--gene_name ${ShellArgs.q(gene_name)}"' in nf
    res = open(os.path.join(ROOT, "modules", "resolve_effective_params.nf")).read()
    assert "'--run_resource ' + ShellArgs.q(" in res
    rq = open(os.path.join(ROOT, "modules", "resolve_query.nf")).read()      # the same easy-mode input, earlier
    assert "resolved_input_path=${ShellArgs.q(resolved_input)}" in rq and '"--species ${ShellArgs.q(species_override)}"' in rq
    assert 'resolved_input_path="${resolved_input}"' not in rq


@pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
def test_shell_args_keep_every_byte(tmp_path):
    values = ["14-3-3(epsilon)", "5'-nucleotidase", "a b", "$(echo x)", "`x`", "a;b", "plain"]
    (tmp_path / "lib").mkdir()
    shutil.copy(os.path.join(ROOT, "lib", "ShellArgs.groovy"), tmp_path / "lib" / "ShellArgs.groovy")
    body = "\n".join(f"    println 'RESULT ' + ShellArgs.q({json.dumps(v)})" for v in values)
    (tmp_path / "t.nf").write_text("workflow {\n" + body.replace("$", "\\$") + "\n}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path, capture_output=True,
                          text=True, timeout=300)
    quoted = [line[len("RESULT "):] for line in proc.stdout.splitlines() if line.startswith("RESULT ")]
    assert len(quoted) == len(values), proc.stdout + proc.stderr
    for v, q in zip(values, quoted):
        assert subprocess.run(["bash", "-c", "printf %s " + q], capture_output=True, text=True).stdout == v
    assert quoted[-1] == "'plain'"                  # a value without ' renders as before


def test_uniprot_symbols_are_quoted_and_the_hit_checked(monkeypatch, tmp_path):
    import resolve_gene_input as rgi
    asked = []

    class Resp:
        def __init__(self, data):
            self.data = data

        def read(self):
            return json.dumps(self.data).encode()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def urlopen(req, timeout=30):
        asked.append(req.full_url)
        hit = {"primaryAccession": "B7YZY9", "genes": [{"geneName": {"value": "Su(dx)"}}]}
        return Resp({"results": [hit]})
    monkeypatch.setattr(rgi.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(rgi, "fetch_uniprot", lambda acc, out: {"accession": acc})
    assert rgi.search_uniprot("Su(H)", "Drosophila melanogaster", str(tmp_path)) is None
    assert "gene_exact%3A%22Su%28H%29%22" in asked[0]


def test_a_uniprot_hit_carrying_the_symbol_is_taken(monkeypatch, tmp_path):
    import resolve_gene_input as rgi
    hit = {"primaryAccession": "P28159", "genes": [{"geneName": {"value": "Su(H)"},
                                                    "synonyms": [{"value": "Suppressor of Hairless"}]}]}
    monkeypatch.setattr(rgi, "fetch_uniprot", lambda acc, out: {"accession": acc})

    class Resp:
        def read(self):
            return json.dumps({"results": [hit]}).encode()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    monkeypatch.setattr(rgi.urllib.request, "urlopen", lambda req, timeout=30: Resp())
    assert rgi.search_uniprot("su(h)", "Drosophila melanogaster", str(tmp_path)) == {"accession": "P28159"}


# --------------------------------------------------------------------------- 28d

def test_structural_templates_are_full_length_goi_proteins():
    import iterative_search_runner as isr
    goi = [("GOI_Hexon", "M" * 300), ("GOI_Hexon|exon_1", "M" * 40), ("GOI_Hexon|frag_1_50_300", "M" * 50)]
    assert isr.full_length_goi_templates(goi) == [("GOI_Hexon", "M" * 300)]
    assert isr.full_length_goi_templates(goi[1:]) == goi[1:]      # none full length: all, as before


# --------------------------------------------------------------------------- 28e

@pytest.mark.parametrize("header,expected", [
    ("lcl|NC_000001.1_cds_XP_100001.1_1 [gene=LOC1] [protein=x]", "XP_100001.1"),
    ("lcl|NC_1.1_cds_NP_000537.3_12", "NP_000537.3"),
    ("lcl|NC_000001.1_prot_XP_100001.1_1 [gene=LOC1]", "XP_100001.1"),   # the "FASTA Protein" download
    ("lcl|CM012345.1_cds_KAF1234567.1_1 [gene=x]", "KAF1234567.1"),            # a GenBank protein
    ("lcl|CM1.1_prot_AAA12345.1_3", "AAA12345.1"),
    ("lcl|X_cds_XP_1.1_1|y", "lcl|X_cds_XP_1.1_1|y"),       # not the download's form: its first word
    ("cds_XP_1.1_x", "cds_XP_1.1_x"),
])
def test_an_ncbi_cds_download_names_its_protein(header, expected):
    assert su.extract_id(header) == expected


@pytest.mark.parametrize("header,expected", [
    ("ref|XP_100002.1| x", "XP_100002.1"), ("gb|AAA12345.1| x", "AAA12345.1"),
    ("gb|AAA12345.1|LOCUS1 x", "AAA12345.1"), ("emb|CAB12345.1|", "CAB12345.1"), ("dbj|BAA12345.1|", "BAA12345.1"),
    ("ENA|CAA12345|CAA12345.1 Homo sapiens", "CAA12345.1"), ("lcl|ORF1 d", "ORF1"),
    # the forms that were read before stay as they were
    ("sp|P12345|LY6E_HUMAN d", "P12345"), ("gi|12345|ref|XP_1.1| x", "12345"), ("P12345|LY6E_HUMAN", "P12345"),
    # no known form: the first word, '|' and all (normalize_query.py then stops)
    ("jgi|Aspni7|1234|fgenesh1_pg.1_#_1", "jgi|Aspni7|1234|fgenesh1_pg.1_#_1"), ("lcl|a|b", "lcl|a|b"),
])
def test_a_tagged_id_reads_as_its_accession(header, expected):
    assert su.extract_id(header) == expected


@pytest.mark.parametrize("header,written", [("ref|XP_100002.1| x", "XP_100002.1"),
                                            ("ENA|CAA12345|CAA12345.1 d", "CAA12345.1")])
def test_normalize_query_writes_the_accession_of_a_tagged_id(tmp_path, header, written):
    (tmp_path / "q.faa").write_text(f">{header}\n{PROT}\n")
    p = subprocess.run([sys.executable, os.path.join(BIN, "normalize_query.py"), "--input", str(tmp_path / "q.faa"),
                        "--output", str(tmp_path / "out.faa")], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "out.faa").read_text().splitlines()[0] == f">{written}"


@pytest.mark.parametrize("qid", ["jgi|Aspni7|1234|fgenesh1_pg.1_#_1", "lcl|a|b", "000000F|arrow"])
def test_normalize_query_stops_on_a_pipe_left_in_the_id(tmp_path, qid):
    """A GOI's name is read up to its first '|': 'jgi', 'lcl' or '000000F' would name the GOI."""
    (tmp_path / "q.faa").write_text(f">{qid} a protein\n{PROT}\n")
    p = subprocess.run([sys.executable, os.path.join(BIN, "normalize_query.py"), "--input", str(tmp_path / "q.faa"),
                        "--output", str(tmp_path / "out.faa")], capture_output=True, text=True)
    assert p.returncode == 1 and qid in p.stderr and "'|'" in p.stderr, p.stderr
    assert qid.replace("|", "_") in p.stderr                    # the rename hint


# --------------------------------------------------------------------------- 28f

def test_a_loc_id_is_a_symbol_and_accessions_keep_their_form():
    import resolve_gene_input as rgi
    for sym in ("LOC123456", "LOC1234567"):          # LOC1234567: three letters and seven digits
        assert rgi.detect_input_type(sym)[0] == "symbol", sym
    for acc in ("AAB12345", "AAB12345.1", "KAF1234567.1", "XP_031329057"):
        assert rgi.detect_input_type(acc)[0] == "ncbi", acc


def test_efetchs_error_text_is_no_fasta(monkeypatch, tmp_path):
    import resolve_gene_input as rgi
    monkeypatch.setattr(rgi.subprocess, "run", lambda *a, **k: types.SimpleNamespace(
        stdout="Error: F a i l e d  t o  u n d e r s t a n d  i d :  L O C 1 2 3 4 5 6", stderr="", returncode=0))
    assert rgi.fetch_ncbi("LOC123456", str(tmp_path)) is None
    assert not list(tmp_path.iterdir())


# --------------------------------------------------------------------------- 34

def test_one_gene_prefix_rule_in_the_figures():
    import plot_synteny as ps
    assert ps.clean_gene_label("gene:ENSG00000141510") == "ENSG00000141510"
    assert ps.clean_gene_label("gene-LOC412898") == "LOC412898"
    assert ps.clean_gene_label("FBgn0000055") == "FBgn0000055"
    ps._GENE_NAMES.clear()
    ps._GENE_NAMES["gene:ENSG00000141510"] = {"label": "TP53"}
    try:
        assert ps._named("ENSG00000141510") == {"label": "TP53"}
    finally:
        ps._GENE_NAMES.clear()


def test_a_long_scaffold_label_keeps_both_ends():
    src = open(os.path.join(BIN, "plot_synteny.py")).read()
    assert 'short = ch if len(ch) <= 20 else ch[:8] + "\\u2026" + ch[-8:]' in src
    assert "ch[-12:]" not in src


# --------------------------------------------------------------------------- the tree page's species map (C1's 28b)

@pytest.mark.parametrize("leaf_species,expected", [
    ("GCA_963969585_1", "Drosophila helvetica"),        # an easy-mode accession key
    ("GCF_011750605_1_ASM", "Drosophila kikkawai"),     # the genome file named beyond the accession
    ("Candida_tropicalis", "Candida tropicalis"),       # a species-stem key with its extension
    ("Apis_cerana_japonica", "Apis cerana japonica"),   # the longer key first
    ("Apis_cerana", "Apis cerana"),
    ("Bombus_terrestris", None),
])
def test_tree_leaves_meet_their_species(leaf_species, expected):
    import plot_synteny as ps
    species_map = {"GCA_963969585.1": "Drosophila helvetica", "GCF_011750605.1": "Drosophila kikkawai",
                   "Candida_tropicalis.fna": "Candida tropicalis", "Apis_cerana": "Apis cerana",
                   "Apis_cerana_japonica": "Apis cerana japonica"}
    assert ps.leaf_species_name(leaf_species, species_map) == expected
