"""Fix-phase item 25a: a genome's sequence name is its header's first word.

The aligners, samtools and the GFF name a sequence by that word; the genome readers keyed it by
extract_id's cut. So a target assembly with SPAdes-style names gave nothing, exit 0 and no warning: the hits
were on NODE_1_length_5000_cov_1.0, the genome held NODE_1, and process_region_block returned nothing at
`chrom not in genome_seqs`. SUPER_10_unloc_1 overwrote SUPER_10. Predicted genes of a home genome without a GFF
(<scaffold>__<start>__<end>_gN) all read as the scaffold's accession, and FUN_000001|exon_2 as FUN_000001.
The patterns locus_tag, generic_gene and generic_protein are gone, two records with one name stop the task, and a
hit or block on a sequence the genome does not hold stops the search.
"""
import random
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import iterative_search_runner as isr  # noqa: E402
import rescue_goi_hull  # noqa: E402
import sequence_utils as su  # noqa: E402
from sequence_utils import StepFailed  # noqa: E402

# (header, the name). The word goes whole, whatever it looks like.
NAMES = [
    ("NODE_1_length_5000_cov_1.0", "NODE_1_length_5000_cov_1.0"),
    ("NODE_12_length_100_cov_5.2 a description", "NODE_12_length_100_cov_5.2"),
    ("SUPER_10", "SUPER_10"),
    ("SUPER_10_unloc_1", "SUPER_10_unloc_1"),
    ("NC_000001.11 Homo sapiens chromosome 1", "NC_000001.11"),
    ("NT_033779.5", "NT_033779.5"),
    ("NC_000007.14_chr7", "NC_000007.14_chr7"),
    ("NC_000006.12:1-5000", "NC_000006.12:1-5000"),
    ("NC_000003.12|chr3", "NC_000003.12|chr3"),
    ("gi|12345|ref|NC_000001.11|", "gi|12345|ref|NC_000001.11|"),
    ("ref|NC_000002.12|", "ref|NC_000002.12|"),
    ("sp|P12345|CHR5", "sp|P12345|CHR5"),
    ("lcl|contig1", "lcl|contig1"),
    ("000000F|arrow", "000000F|arrow"),
    ("trailing|", "trailing|"),
    ("A0A1B2C3D4", "A0A1B2C3D4"),
    ("chr9\ttab description", "chr9"),
    ("  chr10 leading space", "chr10"),
    ("chr11 no-break", "chr11 no-break"),     # part of the word for samtools and miniprot too
    ("1 dna:chromosome chromosome:GRCh38:1:1:248956422:1", "1"),
]


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")
    return str(path)


def _distinct(names):
    """One record per name, with a sequence of its own (lengths 3, 4, ...)."""
    return [(h, "ACGT"[i % 4] * (i + 3)) for i, (h, _n) in enumerate(names)]


@pytest.mark.parametrize("header,name", NAMES)
def test_first_word_is_the_whole_word(header, name):
    assert su.first_word(header) == name


def test_first_word_of_an_empty_header():
    assert su.first_word("") == su.first_word("   ") == "unknown"


READERS = {
    "load_genome": su.load_genome,
    "load_genome_subset": lambda path: su.load_genome_subset(path, [n for _h, n in NAMES]),
    "get_sequence_lengths": su.get_sequence_lengths,
    "stream_fasta_lengths": su.stream_fasta_lengths,
}


@pytest.mark.parametrize("reader", sorted(READERS))
def test_every_genome_reader_keys_by_the_first_word(reader, tmp_path):
    records = _distinct(NAMES)
    got = READERS[reader](write_fasta(tmp_path / "g.fna", records))
    assert list(got) == [n for _h, n in NAMES]
    for (_h, seq), (_h2, name) in zip(records, NAMES):
        assert got[name] in (seq, len(seq))


def test_the_readers_read_a_gzipped_genome(tmp_path):
    import gzip
    path = tmp_path / "g.fna.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">NODE_1_length_4_cov_1.0 d\nACGT\n>SUPER_10\nAC\n")
    assert su.load_genome(str(path)) == {"NODE_1_length_4_cov_1.0": "ACGT", "SUPER_10": "AC"}
    assert su.stream_fasta_lengths(str(path)) == {"NODE_1_length_4_cov_1.0": 4, "SUPER_10": 2}


def test_an_unlocalized_scaffold_no_longer_overwrites_its_chromosome(tmp_path):
    path = write_fasta(tmp_path / "g.fna", [("SUPER_10", "A" * 10), ("SUPER_10_unloc_1", "C" * 4)])
    assert su.load_genome(path) == {"SUPER_10": "A" * 10, "SUPER_10_unloc_1": "C" * 4}
    assert su.stream_fasta_lengths(path) == {"SUPER_10": 10, "SUPER_10_unloc_1": 4}
    assert su.load_genome_subset(path, ["SUPER_10"]) == {"SUPER_10": "A" * 10}


# Two records with one name stop the task, whichever reader sees them first.

DUPLICATE = [("chr1 first copy", "ACGT"), ("chr2", "GG"), ("chr1 second copy", "TTTT")]


@pytest.mark.parametrize("reader", sorted(READERS))
def test_a_name_used_twice_stops_the_task_and_names_both_headers_and_the_file(reader, tmp_path):
    path = write_fasta(tmp_path / "g.fna", DUPLICATE)
    call = (lambda p: su.load_genome_subset(p, ["chr1"])) if reader == "load_genome_subset" else READERS[reader]
    with pytest.raises(StepFailed) as err:
        call(path)
    text = str(err.value)
    assert ">chr1 first copy" in text and ">chr1 second copy" in text and path in text and "chr1" in text


def test_a_name_used_twice_further_down_stops_a_subset_read_too(tmp_path):
    path = write_fasta(tmp_path / "g.fna", [("chr1", "ACGT"), ("chr2", "GG"), ("chr3", "T"), ("chr2 again", "CC")])
    with pytest.raises(StepFailed):
        su.load_genome_subset(path, ["chr1"])


def test_records_that_differ_after_the_first_word_are_the_same_name(tmp_path):
    path = write_fasta(tmp_path / "g.fna", [("chr1 a", "AC"), ("chr1\tb", "GT")])
    with pytest.raises(StepFailed):
        su.load_genome(path)


def test_the_stop_passes_the_fallbacks_of_the_search():
    """StepFailed is a BaseException, so `except Exception` around a block does not hold it."""
    assert not issubclass(StepFailed, Exception)


def test_the_quality_check_stops_at_a_name_used_twice(tmp_path):
    """The first task that reads a plain-file target: the run ends before any search."""
    bad = write_fasta(tmp_path / "bad.fna", DUPLICATE)
    good = write_fasta(tmp_path / "good.fna", [("NODE_1_length_4_cov_1.0", "ACGT"), ("NODE_2_length_2_cov_1.0", "GG")])
    script = str(ROOT / "bin" / "assess_genome_quality.py")
    out = tmp_path / "qc.json"
    done = subprocess.run([sys.executable, script, "--genome", bad, "--output", str(out)],
                          capture_output=True, text=True)
    assert done.returncode == 1
    assert ">chr1 first copy" in done.stdout and ">chr1 second copy" in done.stdout and bad in done.stdout
    assert not out.exists()
    done = subprocess.run([sys.executable, script, "--genome", good, "--output", str(out)],
                          capture_output=True, text=True)
    assert done.returncode == 0 and out.exists()


# The hull rescue's two ways to a window agree (the fallback read extract_id's cut).

def test_the_hull_rescues_fallback_reads_the_first_word(tmp_path):
    path = write_fasta(tmp_path / "g.fna", [("NODE_1_length_8_cov_1.0", "ACGTACGT"), ("NODE_1_length_4_cov_1.0", "TTTT")])
    get = rescue_goi_hull._extract_window
    assert get(path, "NODE_1_length_8_cov_1.0", 2, 4, samtools_bin="no-such-samtools") == "CGT"
    assert get(path, "NODE_1_length_4_cov_1.0", 1, 4, samtools_bin="no-such-samtools") == "TTTT"
    assert get(path, "NODE_1", 1, 4, samtools_bin="no-such-samtools") == ""


@pytest.mark.skipif(not shutil.which("samtools"), reason="samtools not on PATH")
@pytest.mark.parametrize("header,name", [n for n in NAMES if not n[0].startswith(" ")])
def test_the_readers_name_a_sequence_as_samtools_does(header, name, tmp_path):
    records = _distinct(NAMES)
    path = write_fasta(tmp_path / "g.fna", [r for r in records if not r[0].startswith(" ")])
    subprocess.run(["samtools", "faidx", path], check=True)
    with open(path + ".fai", encoding="utf-8") as fh:
        fai = [line.split("\t")[0] for line in fh]
    assert name in fai
    assert fai == list(su.stream_fasta_lengths(path))
    seq = dict(su.iter_genome(path))[name]
    assert rescue_goi_hull._extract_window(path, name, 1, len(seq)) == seq
    assert rescue_goi_hull._extract_window(path, name, 1, len(seq), samtools_bin="no-such-samtools") == seq


@pytest.mark.skipif(not shutil.which("mmseqs"), reason="mmseqs not on PATH")
def test_mmseqs2_keeps_these_names_and_rewrites_those(tmp_path):
    """The premise of the stop below, on the installed MMseqs2: the names of the patterns that are gone
    come back whole, and a database tag does not."""
    rng = random.Random(25)
    code = {"M": "ATG", "K": "AAA", "T": "ACT", "A": "GCT", "Y": "TAT", "I": "ATT", "Q": "CAA", "R": "CGT", "S": "TCT",
            "F": "TTT", "V": "GTT", "H": "CAT", "L": "CTG", "E": "GAA", "G": "GGT", "P": "CCT", "D": "GAT", "N": "AAT",
            "W": "TGG", "C": "TGT"}
    prot = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(150))
    cds = "".join(code[a] for a in prot)
    names = {"NODE_1_length_5000_cov_1.0": "NODE_1_length_5000_cov_1.0", "SUPER_10_unloc_1": "SUPER_10_unloc_1",
             "000000F|arrow": "000000F|arrow", "lcl|contig1": "contig1"}
    genome = write_fasta(tmp_path / "g.fna", [(n, "".join(rng.choice("ACGT") for _ in range(400)) + cds
                                               + "".join(rng.choice("ACGT") for _ in range(400))) for n in names])
    query = write_fasta(tmp_path / "q.faa", [("q1", prot)])
    out = tmp_path / "hits.m8"
    subprocess.run(["mmseqs", "easy-search", query, genome, str(out), str(tmp_path / "tmp"), "--search-type", "3",
                    "--format-output", "target", "--threads", "1", "-v", "0"], check=True, capture_output=True)
    assert {line.strip() for line in out.read_text().splitlines()} == set(names.values())


# A hit or block on a sequence the genome does not hold.

def test_a_block_on_a_sequence_the_genome_does_not_hold_stops_the_search():
    with pytest.raises(StepFailed) as err:
        isr.process_region_block(3, {"chrom": "NODE_1_length_5000_cov_1.0"}, [], {"NODE_1": "ACGT"}, {},
                                 "Target.fna", SimpleNamespace(), "u", 1)
    text = str(err.value)
    assert "NODE_1_length_5000_cov_1.0" in text and "Target.fna" in text and "NODE_1" in text


def test_the_message_names_the_sequence_and_the_genomes_names():
    """A backstop since fix item 29 (every MMseqs2 table is read by its FASTA headers): it names the
    sequence, the genome and three of its names, and no longer blames MMseqs2's rewriting."""
    text = str(su.sequence_not_in_genome("contig1", "Target.fna", ["lcl|contig1", "lcl|contig2", "lcl|contig3", "lcl|contig4"]))
    assert "contig1" in text and "Target.fna" in text and "lcl|contig1" in text and "first word" in text
    assert "4 sequence names" in text and "lcl|contig4" not in text


def test_the_off_chromosome_path_stops_too():
    """By its source: the check of the table's hits fires first, so this line is a backstop."""
    src = open(ROOT / "bin" / "iterative_search_runner.py").read()
    at = src.index("if off_chrom not in genome_seqs:")
    assert src[at:at + 130].split("\n")[1].strip() == "raise sequence_not_in_genome(off_chrom, genome_name, genome_seqs)"


def _search_args(tmp_path):
    (tmp_path / "hits").mkdir(exist_ok=True)
    return SimpleNamespace(output_dir=str(tmp_path), min_identity=20.0, min_length=10, evalue=1.0)


def test_a_hit_on_a_sequence_the_genome_does_not_hold_stops_the_search(tmp_path, monkeypatch):
    """Before any block forms: a genome with too few anchors for a block is told too. The real call,
    with the MMseqs2 search and its table replaced; lcl|contig1 is what MMseqs2's own target column
    reports as contig1 (since fix item 29 the table is read by its headers, so this is a backstop)."""
    genome = write_fasta(tmp_path / "Target.fna", [("lcl|contig1", "ACGT" * 50)])
    db = write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV"), ("GOI_x", "MKT")])
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", lambda **k: (True, "", False))
    monkeypatch.setattr(isr, "parse_hits", lambda *a, **k: [{"query": "gene-a|exon_1", "chrom": "contig1", "start": 0, "end": 30}])
    name, genes, extras, error = isr.process_single_genome(genome, db, _search_args(tmp_path), None, "", 1)
    assert name == "Target.fna" and genes == [] and extras == []
    assert "contig1" in error and "lcl|contig1" in error and "first word" in error


def test_a_hit_on_the_genomes_own_name_passes_that_check(tmp_path, monkeypatch):
    genome = write_fasta(tmp_path / "Target.fna", [("NODE_1_length_200_cov_1.0", "ACGT" * 50)])
    db = write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV")])
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", lambda **k: (True, "", False))
    monkeypatch.setattr(isr, "parse_hits", lambda *a, **k: [{"query": "gene-a|exon_1", "chrom": "NODE_1_length_200_cov_1.0",
                                                             "start": 0, "end": 30}])
    monkeypatch.setattr(isr.synvoy_evalue, "places_per_gene", lambda hits: (_ for _ in ()).throw(RuntimeError("past the check")))
    _name, _g, _e, error = isr.process_single_genome(genome, db, _search_args(tmp_path), None, "", 1)
    assert error == "past the check"


def test_a_name_used_twice_stops_the_search_before_mmseqs2_runs(tmp_path, monkeypatch):
    genome = write_fasta(tmp_path / "Target.fna", DUPLICATE)
    db = write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV")])

    def must_not_run(**_kw):
        raise AssertionError("MMseqs2 ran on a genome with a name used twice")
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", must_not_run)
    _name, genes, _e, error = isr.process_single_genome(genome, db, _search_args(tmp_path), None, "", 1)
    assert genes == [] and ">chr1 first copy" in error and ">chr1 second copy" in error and genome in error


# The ladder: what is left for a user's header.

@pytest.mark.parametrize("kind", ["locus_tag", "generic_gene", "generic_protein"])
def test_the_generic_patterns_are_gone(kind):
    assert kind not in [k for _p, k in su.ID_PATTERNS]


@pytest.mark.parametrize("header,expected", [
    # Cut on the parent, whole now
    ("NODE_1_length_5000_cov_1.0", "NODE_1_length_5000_cov_1.0"),
    ("SUPER_10_unloc_1", "SUPER_10_unloc_1"),
    ("FUN_000001|exon_2", "FUN_000001|exon_2"),
    ("NT_033779.5__9376937__9477759_g46", "NT_033779.5__9376937__9477759_g46"),
    ("NT_033779.5__9376937__9477759_g46|exon_1 parent=x", "NT_033779.5__9376937__9477759_g46|exon_1"),
    ("ABC_12345.1", "ABC_12345.1"),
    ("XY_0000001.1 a description", "XY_0000001.1"),
    ("gene_12345_b", "gene_12345_b"),
    ("protein_12345.2", "protein_12345.2"),
    ("prot_12|x", "prot_12|x"),
    ("GCA_018908115.1.fna_structdisc_0", "GCA_018908115.1.fna_structdisc_0"),
    # The same on the parent
    ("ABC_12345 locus tag style", "ABC_12345"),
    ("XY_0000001 locus=x start=y", "XY_0000001"),
    ("gene_12345", "gene_12345"),
    ("protein_12345 d", "protein_12345"),
    ("XP_123456.1 a protein", "XP_123456.1"),
    ("NC_000001.11 Homo sapiens", "NC_000001.11"),
    ("P12345 d", "P12345"),
    ("sp|P12345|LY6E_HUMAN d", "P12345"),
    ("P12345|LY6E_HUMAN", "P12345"),
    ("gi|12345|ref|XP_1.1|", "12345"),
    ("1_234 # 1000 # 2000 # 1 # ID=1_234;partial=00", "1_234"),
    ("LOC143834063 d", "LOC143834063"),
    ("AAA12345.1 d", "AAA12345.1"),
])
def test_the_ladder_keeps_its_known_forms(header, expected):
    assert su.extract_id(header) == expected


def test_no_generic_pattern_comes_back_in_another_form():
    """A word of the old patterns' shape, whatever follows the digits, reads whole."""
    for word in ("ABCDEF_1x", "AB_1-2", "XY_12_34", "gene_1_2", "gene_1.2", "protein_1x", "prot_1_", "GOI_1-2", "FUN_1|exon_1"):
        assert su.extract_id(word) == word
