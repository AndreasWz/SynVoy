"""Fix-phase item 25b: a file SynVoy wrote is read by the first word of each header.

SynVoy builds ids as <id>|<part>: flanking exons (gene|exon_2), fragments (gene|frag_1_10_20), models
(parent|<genome>_b0_l1_exon_ann), panel entries (locus_1|gene). extract_id's accession patterns stop at '|', so
with a home GFF whose gene ids are FlyBase, TAIR or bare RefSeq or Ensembl accessions, FBgn0000001|exon_2 read as
FBgn0000001: the exon records of a gene merged in PREPARE_INITIAL_DB, and every model of one parent merged when a
model FASTA was read back. parse_fasta now takes the first word; the ladder of known id forms stays where a user's
header enters the pipeline, normalize_query.py, which writes the query's id once.
"""
import ast
import random
import string
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import build_home_paralog_panel as panel  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import sequence_utils as su  # noqa: E402
import synvoy_scale  # noqa: E402
import synvoy_specificity  # noqa: E402

# (first word of the header, what the ladder reads at the query entry)
FORMS = [
    # The forms the ladder cut at '|'
    ("FBgn0000001|exon_2", "FBgn0000001"),
    ("FBpp0000001|Target_fna_b0_l1_exon_ann", "FBpp0000001"),
    ("AT1G01010|exon_1", "AT1G01010"),
    ("AT1G01010.1|frag_1_10_20", "AT1G01010.1"),
    ("XP_123456.1|exon_3", "XP_123456.1"),
    ("NP_123456|Target_fna_b0_l1_flanking_miniprot", "NP_123456"),
    ("ENSG00000123456|exon_1", "ENSG00000123456"),
    ("ENSP00000123456.2|slide_1_50_200", "ENSP00000123456.2"),
    ("ENST00000123456|exon_1", "ENST00000123456"),
    ("WBGene00000001|exon_1", "WBGene00000001"),
    ("AAA12345.1|exon_2", "AAA12345.1"),
    ("NC_000001.11|exon_1", "NC_000001.11"),
    ("P12345|LY6E_HUMAN", "P12345"),
    ("Q9NXB0|exon_1", "Q9NXB0"),
    ("sp|P12345|LY6E_HUMAN", "P12345"),
    ("gi|12345|ref|XP_1.1|", "12345"),
    # The forms the generic patterns cut (25a): whole in both readings now
    ("FUN_000001|exon_2", "FUN_000001|exon_2"),
    ("NT_033779.5__9376937__9477759_g46", "NT_033779.5__9376937__9477759_g46"),
    ("NT_033779.5__9376937__9477759_g46|exon_1", "NT_033779.5__9376937__9477759_g46|exon_1"),
    ("NODE_1_length_5000_cov_1.0", "NODE_1_length_5000_cov_1.0"),
    ("SUPER_10_unloc_1", "SUPER_10_unloc_1"),
    ("ABC_12345.1", "ABC_12345.1"),
    ("gene_12345_b|exon_1", "gene_12345_b|exon_1"),
    ("GCA_018908115.1.fna_structdisc_0", "GCA_018908115.1.fna_structdisc_0"),
    # The same in both readings on the parent too
    ("GOI_14|exon_1", "GOI_14|exon_1"),
    ("GOI_14-3-3epsilon|exon_1", "GOI_14-3-3epsilon|exon_1"),
    ("GOI_Ly6e", "GOI_Ly6e"),
    ("GOI_copy_1|Apis_b0_l1_by_GOI_Melt", "GOI_copy_1|Apis_b0_l1_by_GOI_Melt"),
    ("GOI_lcl|ORF1", "GOI_lcl|ORF1"),
    ("gene-LOC123456|exon_1", "gene-LOC123456|exon_1"),
    ("locus_1|gene-X", "locus_1|gene-X"),
    ("home_gene-BGN|gene-BGN", "home_gene-BGN|gene-BGN"),
    ("XP_123456.1", "XP_123456.1"),
    ("FBgn0000001", "FBgn0000001"),
    ("1_234", "1_234"),
]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")
    return str(path)


def test_a_synvoy_file_is_read_by_the_whole_first_word(tmp_path):
    path = write_fasta(tmp_path / "x.faa", [(f"{word} home=1-2 d", "MKV" + "A" * i) for i, (word, _l) in enumerate(FORMS)])
    got = list(su.parse_fasta(path))
    assert [g[1] for g in got] == [word for word, _l in FORMS]
    assert [g[0] for g in got] == [f"{word} home=1-2 d" for word, _l in FORMS]      # the raw header stays


def test_a_users_header_goes_through_the_ladder(tmp_path):
    path = write_fasta(tmp_path / "q.faa", [(word, "MKV") for word, _l in FORMS])
    got = list(su.parse_fasta(path, user_header=True))
    assert [g[1] for g in got] == [ladder for _w, ladder in FORMS]


def test_the_exon_records_of_one_gene_stay_apart(tmp_path):
    for gene in ("FBgn0000001", "AT1G01010", "XP_123456.1", "ENSG00000123456", "FUN_000001", "gene-a"):
        path = write_fasta(tmp_path / "x.faa", [(gene, "M" * 9), (f"{gene}|exon_1", "M" * 3), (f"{gene}|exon_2", "M" * 6)])
        assert su.parse_fasta_to_dict(path) == {gene: "M" * 9, f"{gene}|exon_1": "M" * 3, f"{gene}|exon_2": "M" * 6}


def test_two_models_of_one_parent_stay_apart(tmp_path):
    path = write_fasta(tmp_path / "models.faa", [("XP_1.1|T_fna_b0_l1_flanking_miniprot", "MK"),
                                                 ("XP_1.1|T_fna_b1_l1_flanking_miniprot", "MKV")])
    assert len(su.parse_fasta_to_dict(path)) == 2


def test_random_ids_read_back_as_written(tmp_path):
    """Any id without whitespace, written with write_fasta, reads back as it is."""
    rng = random.Random(25)
    alphabet = string.ascii_letters + string.digits + "._-|~:;=,#()[]"
    ids = sorted({"".join(rng.choice(alphabet) for _ in range(rng.randint(1, 30))) for _ in range(5000)})
    path = tmp_path / "x.faa"
    su.write_fasta([(i, "MKV") for i in ids], str(path))
    assert [g[1] for g in su.parse_fasta(str(path))] == ids


# The query entry

def _normalize(tmp_path, header):
    q = write_fasta(tmp_path / "q.faa", [(header, "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKA")])
    out = tmp_path / "n.faa"
    subprocess.run([sys.executable, str(ROOT / "bin" / "normalize_query.py"), "--input", q, "--output", str(out)],
                   check=True, capture_output=True, text=True)
    return out


@pytest.mark.parametrize("header,expected", [
    ("P12345|LY6E_HUMAN", "P12345"),
    ("P12345|LY6E_HUMAN a protein", "P12345"),
    ("sp|P12345|LY6E_HUMAN Lymphocyte antigen OS=Homo sapiens", "P12345"),
    ("tr|A0A1B2C3D4|A0A1B2C3D4_APIME x", "A0A1B2C3D4"),
    ("gi|12345|ref|XP_1.1| x", "12345"),
    ("XP_123456.1 a protein", "XP_123456.1"),
    ("XY_0000001 locus=x start=y", "XY_0000001"),
    ("FAM_Genus_species_XY_0000002", "FAM_Genus_species_XY_0000002"),
    ("orf19.1562", "orf19.1562"),
    ("LOC123456789", "LOC123456789"),
    ("lcl|ORF1 d", "lcl|ORF1"),
    # What item 25 changes at this entry: the word stays whole
    ("ABC_12345.1", "ABC_12345.1"),
    ("XY_0000001_exon2", "XY_0000001_exon2"),
    ("FUN_000001|exon_2", "FUN_000001|exon_2"),
])
def test_normalize_query_writes_the_id_of_the_ladder(tmp_path, header, expected):
    assert _normalize(tmp_path, header).read_text().splitlines()[0] == f">{expected}"


@pytest.mark.parametrize("header,expected", [
    ("P12345|LY6E_HUMAN", "P12345"), ("sp|P12345|LY6E_HUMAN d", "P12345"), ("gi|12345|ref|XP_1.1| d", "12345"),
    ("XP_123456.1 d", "XP_123456.1"), ("lcl|ORF1", "lcl|ORF1"), ("ABC_12345.1", "ABC_12345.1"),
])
def test_every_later_step_reads_the_id_normalize_query_wrote(tmp_path, header, expected):
    """The query's id is settled once: a later reader gets the same id, ladder or not."""
    out = str(_normalize(tmp_path, header))
    assert [g[1] for g in su.parse_fasta(out)] == [expected]
    assert [g[1] for g in su.parse_fasta(out, user_header=True)] == [expected]


# Only the query entry uses the ladder.

def _calls(tree, name):
    return [n for n in ast.walk(tree) if isinstance(n, ast.Call)
            and getattr(n.func, "id", getattr(n.func, "attr", None)) == name]


def test_nothing_but_the_entries_of_a_users_header_uses_the_ladder():
    entries, extract = [], []
    for path in sorted((ROOT / "bin").glob("*.py")):
        tree = ast.parse(path.read_text())
        for call in _calls(tree, "parse_fasta"):
            if any(k.arg == "user_header" for k in call.keywords):
                entries.append(path.name)
        if path.name != "sequence_utils.py" and _calls(tree, "extract_id"):
            extract.append(path.name)
    assert entries == ["fragment_query.py", "normalize_query.py"]
    assert extract == []
    for path in sorted((ROOT / "modules").glob("*.nf")) + sorted((ROOT / "subworkflows").glob("*.nf")) + [ROOT / "main.nf"]:
        text = path.read_text()
        assert "user_header" not in text and "extract_id" not in text, path.name


def test_the_query_reaches_the_later_steps_only_through_normalize_query():
    """main.nf hands the raw query to NORMALIZE_QUERY alone: every other process takes its output."""
    text = (ROOT / "main.nf").read_text()
    uses = [line.strip() for line in text.splitlines() if "raw_gene_ch" in line and not line.strip().startswith("//")]
    assert len(uses) == 3 and sum("NORMALIZE_QUERY(raw_gene_ch)" in u for u in uses) == 1
    assert sum("raw_gene_ch =" in u for u in uses) == 2
    assert text.count("fromPath(params.query)") == 1 and "path(params.query)" not in text


# The readers of SynVoy's own files

def _prepare_initial_db(tmp_path, flanking_faa, goi_exons):
    """The module's own script, with Nextflow's substitutions done by hand."""
    nf = (ROOT / "modules" / "prepare_initial_db.nf").read_text()
    body = nf.split('"""')[1]
    for key, value in (("${projectDir}", str(ROOT)), ("${flanking_faa}", str(flanking_faa)),
                       ("${goi_exons}", str(goi_exons)), ("${locus_id}", "locus_1")):
        body = body.replace(key, value)
    script = tmp_path / "prepare.py"
    script.write_text(textwrap.dedent(body).replace("\\\\", "\\"))
    done = subprocess.run([sys.executable, str(script)], cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    return list(su.parse_fasta(str(tmp_path / "initial_db_locus_1.faa")))


@pytest.mark.parametrize("gene", ["FBgn0000001", "AT1G01010", "XP_123456.1", "ENSG00000123456", "FUN_000001",
                                  "NT_033779.5__9376937__9477759_g46"])
def test_prepare_initial_db_keeps_each_exon_record(tmp_path, gene):
    flank = write_fasta(tmp_path / "flank.faa", [
        (f"{gene}|exon_1 parent={gene} exon=1/2 coords=100-130 strand=+", "MKVLA" + "A" * 5),
        (f"{gene}|exon_2 parent={gene} exon=2/2 coords=200-260 strand=+", "GGGGG" + "G" * 15),
        (f"{gene} full_length_protein exons=2 strand=+ home=100-260 genome=5000 goi=500-600", "MKVLA" + "A" * 5 + "G" * 20)])
    goi = write_fasta(tmp_path / "goi.faa", [("GOI_x", "MKT" * 20), ("GOI_x|exon_1", "MKT" * 20)])
    got = _prepare_initial_db(tmp_path, flank, goi)
    assert [g[1] for g in got] == [f"{gene}|exon_1", f"{gene}|exon_2", gene, "GOI_x", "GOI_x|exon_1"]
    # the home fields stay behind the id, and the search reads its geometry from them
    assert got[2][0].endswith("home=100-260 genome=5000 goi=500-600") and got[1][0].endswith("coords=200-260")


def test_the_searchs_distances_see_every_gene_of_a_query_set(tmp_path):
    """Written by PREPARE_INITIAL_DB and read by synvoy_scale: 3 anchors, each exon its own record."""
    genes = ["NT_033779.5__9376937__9477759_g45", "NT_033779.5__9376937__9477759_g46", "NT_033779.5__9376937__9477759_g47"]
    flank = []
    for k, g in enumerate(genes):
        lo = 1000 * (k + 1)
        flank += [(f"{g}|exon_1 parent={g} exon=1/1 coords={lo}-{lo + 300} strand=+", "MKVLAGGGGG" + "A" * 10),
                  (f"{g} full_length_protein exons=1 strand=+ home={lo}-{lo + 300} genome=50000 goi=9000-9100", "MKVLAGGGGG" + "A" * 10)]
    _prepare_initial_db(tmp_path, write_fasta(tmp_path / "flank.faa", flank),
                        write_fasta(tmp_path / "goi.faa", [("GOI_x", "MKT" * 20), ("GOI_x|exon_1", "MKT" * 20)]))
    geo = synvoy_scale.home_geometry(str(tmp_path / "initial_db_locus_1.faa"))
    assert geo["order"] == genes and geo["n"] == 3 and f"{genes[1]}|exon_1" in geo["records"]


def test_the_search_reads_its_query_set_by_the_whole_id(tmp_path):
    db = write_fasta(tmp_path / "db.faa", [("FBgn1|exon_1", "M" * 10), ("FBgn1|exon_2", "M" * 20), ("FBgn1", "M" * 30),
                                           ("XP_2.1|exon_1", "K" * 10), ("XP_2.1", "K" * 10), ("GOI_x", "A" * 30), ("GOI_x|exon_1", "A" * 30)])
    assert isr.build_home_goi_ids(db) == {"GOI_x", "GOI_x|exon_1"}
    assert isr.build_home_rank(db) == {"FBgn1": 0, "XP_2.1": 1}
    # the anchors of a gene's parts: extract_base_gene_id reads the part before the first '|'
    assert {isr.extract_base_gene_id(g[1]) for g in su.parse_fasta(db)} == {"FBgn1", "XP_2.1", "GOI_x"}


def test_the_candidate_models_of_a_region_are_found_by_their_whole_id(tmp_path):
    """GPU layers: the model ids of the GFF are looked up in the model FASTA."""
    faa = write_fasta(tmp_path / "r.faa", [("XP_1.1|T_fna_b0_l1_exon_ann", "M" * 12), ("XP_1.1|T_fna_b1_l1_exon_ann", "K" * 12)])
    gff = tmp_path / "r.gff"
    gff.write_text("c1\tSynVoy\tmRNA\t1\t36\t.\t+\t.\tID=XP_1.1|T_fna_b0_l1_exon_ann;SynVoyRole=goi\n"
                   "c1\tSynVoy\tmRNA\t50\t86\t.\t+\t.\tID=XP_1.1|T_fna_b1_l1_exon_ann;SynVoyRole=goi\n")
    assert isr._select_goi_candidates_from_region(str(gff), faa) == [("XP_1.1|T_fna_b0_l1_exon_ann", "M" * 12),
                                                                     ("XP_1.1|T_fna_b1_l1_exon_ann", "K" * 12)]


def test_the_home_proteome_keeps_every_predicted_gene(tmp_path):
    """Predicted genes of a home genome without a GFF, on a RefSeq-style scaffold (the fly gate case),
    and bare accessions, as the search's specificity check and the panel read them."""
    ids = ["NT_033779.5__9376937__9477759_g45", "NT_033779.5__9376937__9477759_g46", "XP_1.1", "XP_1.2", "gene-a"]
    (tmp_path / "home_proteome.faa").write_text("".join(f">{i} # 1 # 30 # 1 # ID=1\nMKV{'A' * k}\n" for k, i in enumerate(ids)))
    (tmp_path / "home_genes.tsv").write_text("gene_id\tchrom\tstart\tend\tstrand\n" + "".join(f"{i}\tc1\t1\t30\t+\n" for i in ids))
    prot = synvoy_specificity.HomeProteome.load(str(tmp_path))
    assert [r[0] for r in prot.records] == ids
    assert su.parse_fasta_to_dict(str(tmp_path / "home_proteome.faa"))[ids[1]] == "MKVA"


def test_the_panel_reads_the_home_goi_by_its_whole_id(tmp_path):
    path = write_fasta(tmp_path / "goi_exons.faa", [("GOI_14-3-3e", "M" * 40), ("GOI_14-3-3e|exon_1", "M" * 10), ("GOI_copy_1", "K" * 50)])
    assert panel.home_goi_record(path) == ("GOI_14-3-3e", "M" * 40)


def test_the_smith_waterman_child_reports_each_query_by_its_whole_id(tmp_path):
    pytest.importorskip("parasail")
    code = {"M": "ATG", "K": "AAA", "T": "ACT", "A": "GCT", "Y": "TAT", "I": "ATT", "Q": "CAA", "R": "CGT", "S": "TCT",
            "F": "TTT", "V": "GTT", "H": "CAT", "L": "CTG", "E": "GAA", "G": "GGT", "P": "CCT", "D": "GAT", "N": "AAT",
            "W": "TGG", "C": "TGT"}
    rng = random.Random(25)
    prot1 = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80))
    prot2 = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80))
    dna = "".join(rng.choice("ACGT") for _ in range(300)) + "".join(code[a] for a in prot1) + \
        "".join(rng.choice("ACGT") for _ in range(300)) + "".join(code[a] for a in prot2)
    query = write_fasta(tmp_path / "q.faa", [("FBgn0000001|exon_1", prot1), ("FBgn0000001|exon_2", prot2)])
    target = write_fasta(tmp_path / "t.fna", [("region_seq", dna)])
    out = tmp_path / "hits.m8"
    subprocess.run([sys.executable, str(ROOT / "bin" / "smith_waterman_search.py"), "--query", query, "--target", target,
                    "--output", str(out), "--method", "parasail"], check=True, capture_output=True, text=True)
    rows = [line.split("\t") for line in out.read_text().splitlines() if line.strip()]
    assert {r[0] for r in rows} == {"FBgn0000001|exon_1", "FBgn0000001|exon_2"}
    assert {r[1] for r in rows} == {"region_seq"}


def test_a_home_genome_without_a_gff_has_one_query_per_predicted_gene(tmp_path):
    """The fly gate case: 14 anchors read as 9 because 6 predicted genes read as one id."""
    genes = [f"NT_033779.5__9376937__9477759_g{i}" for i in range(40, 46)]
    db = write_fasta(tmp_path / "db.faa", [(f"{g}|exon_1", "M" * 10) for g in genes] + [(g, "M" * 10) for g in genes])
    assert set(isr.build_home_rank(db)) == set(genes)
