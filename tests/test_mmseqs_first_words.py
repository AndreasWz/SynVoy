"""Fix item 29 (9 Nov 2026): every MMseqs2 table the pipeline parses names its rows by the first word of the
FASTA headers (the qheader and theader columns), not by MMseqs2's own query and target columns.

MMseqs2 18-8cc5c rewrites names in those columns (measured on 31 Oct, fix-phase item 25): a name that
starts with a database tag by the accession inside (lcl|X, ref|X|, sp|X|Y, gnl|db|X, gi|1|ref|X| as X;
gnl|a|1 and gnl|b|1 both as 1), a leading uc or consensus_ dropped (ucon_5 as on_5, consensus_24 as 24,
which can be another record's name). tblastn, miniprot and samtools keep the whole first word. Since
2026-11-02 such names stopped the run (fix-phase d4); now they work. The four live tables: LOCATE_GENE,
the search's genome-wide table, its window table, and the annotation transfer of a run without a home GFF.
"""
import os
import random
import shutil
import sys
from types import SimpleNamespace

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import sequence_utils as su  # noqa: E402

needs_mmseqs = pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
         'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA', 'R': 'CGT', 'S': 'TCT',
         'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT'}


def _inputs(tmp_path):
    """Four contigs whose names MMseqs2 rewrites (consensus_24 and the real contig 24 read alike), and
    three queries, two with names it rewrites."""
    rng = random.Random(29)
    prots = {q: "M" + "".join(rng.choice(AA) for _ in range(150)) for q in ("ucp1|exon_1", "gene-b|exon_1", "GOI_x")}
    rnd = lambda k: "".join(rng.choice("ACGT") for _ in range(k))  # noqa: E731
    contigs = [("lcl|contig1", "ucp1|exon_1"), ("consensus_24", "gene-b|exon_1"), ("24", None), ("ucon_5", "GOI_x")]
    genome, queries = tmp_path / "Target.fna", tmp_path / "db.faa"
    genome.write_text("".join(f">{name} some description\n" + rnd(400)
                              + ("".join(CODON[a] for a in prots[prot]) if prot else rnd(453)) + rnd(400) + "\n"
                              for name, prot in contigs))
    queries.write_text("".join(f">{k}\n{v}\n" for k, v in prots.items()))
    return str(genome), str(queries)


EXPECTED = {("ucp1|exon_1", "lcl|contig1"), ("gene-b|exon_1", "consensus_24"), ("GOI_x", "ucon_5")}


# --------------------------------------------------------------------------- the row rule

def test_a_row_is_named_by_the_first_words_of_its_headers():
    row = "\t".join(["p1", "contig1", "90.0", "100", "0", "0", "1", "100", "401", "700", "1e-30", "200",
                     "ucp1|exon_1", "lcl|contig1 some description"]) + "\n"
    assert su.mmseqs_first_words(row) == ["ucp1|exon_1", "lcl|contig1", "90.0", "100", "0", "0", "1", "100",
                                          "401", "700", "1e-30", "200"]


def test_a_tab_inside_the_target_header_stays_in_it():
    row = "\t".join(["GOI_x", "1", "90.0", "100", "0", "0", "1", "100", "1", "300", "1e-30", "200",
                     "GOI_x", "gnl|a|1 a\tdescription\twith tabs"]) + "\n"
    assert su.mmseqs_first_words(row)[:2] == ["GOI_x", "gnl|a|1"]


def test_a_tab_inside_the_query_header_stops():
    """The target column would be the first word of the query's description: MMseqs2's own target name
    (contig1) is no part of it, so the row stops instead of naming the wrong target."""
    row = "\t".join(["GOI_x", "contig1", "90.0", "100", "0", "0", "1", "100", "1", "300", "1e-30", "200",
                     "GOI_x\tdescription", "lcl|contig1 d"]) + "\n"
    with pytest.raises(su.StepFailed, match="misaligned"):
        su.mmseqs_first_words(row)


@pytest.mark.parametrize("name,word", [("X", "lcl|X"), ("1", "gnl|a|1"), ("on_5", "ucon_5"), ("24", "consensus_24"),
                                       ("P1", "sp|P1|NAME"), ("X", "gi|1|ref|X|"), ("", "pir||X"), ("q", "q")])
def test_every_measured_rewrite_is_a_part_of_the_word(name, word):
    row = "\t".join(["q", name, "90.0", "100", "0", "0", "1", "100", "1", "300", "1e-30", "200", "q d", f"{word} d"])
    assert su.mmseqs_first_words(row)[:2] == ["q", word]


def test_a_row_without_the_headers_stops():
    with pytest.raises(su.StepFailed, match="12 BLAST-tabular ones and the two headers"):
        su.mmseqs_first_words("\t".join(["q", "t", "90.0", "100", "0", "0", "1", "100", "1", "300", "1e-30", "200"]))


def test_the_table_keeps_its_rows_and_columns(tmp_path):
    raw, out = tmp_path / "h.raw", tmp_path / "h.m8"
    rows = [["q1", "c1", "90.0", "100", "0", "0", "1", "100", "1", "300", "1e-30", "200", "q1 x", "lcl|c1 y"],
            ["q2", "1", "80.0", "50", "1", "0", "5", "55", "301", "151", "1e-5", "60", "q2", "gnl|b|1"]]
    raw.write_text("".join("\t".join(r) + "\n" for r in rows))
    assert su.mmseqs_table_first_words(str(raw), str(out)) == 2
    assert out.read_text().splitlines() == ["\t".join(["q1", "lcl|c1"] + rows[0][2:12]),
                                            "\t".join(["q2", "gnl|b|1"] + rows[1][2:12])]


# --------------------------------------------------------------------------- the four tables, with MMseqs2

@needs_mmseqs
def test_the_genome_table_names_what_the_fasta_names(tmp_path):
    import iterative_search_runner as isr
    genome, queries = _inputs(tmp_path)
    args = SimpleNamespace(mmseqs_sens=7.5, mmseqs_split_memory_limit="2G", evalue=1e-3, mmseqs_verbosity=1)
    out = str(tmp_path / "hits.m8")
    ok, _det, _rf = isr.run_mmseqs_easy_search_with_retries(db_path=queries, genome_path=genome, hits_file=out,
                                                            tmp_dir=str(tmp_path / "tmp"), args=args,
                                                            threads_per_job=1, genome_name="Target.fna")
    assert ok
    rows = [line.split("\t") for line in open(out)]
    assert {(r[0], r[1]) for r in rows} == EXPECTED
    assert all(len(r) == 12 for r in rows) and not os.path.exists(out + ".raw")


@needs_mmseqs
def test_the_window_table_names_what_the_fasta_names(tmp_path):
    """The window search's queries and target are SynVoy's own (GOI_*, region_seq): changing it is
    defensive; it must read the same table as before."""
    import iterative_search_runner as isr
    genome, queries = _inputs(tmp_path)
    args = SimpleNamespace(mmseqs_sens=7.5, mmseqs_split_memory_limit="2G", evalue=1e-3, mmseqs_verbosity=1,
                           aug_relaxed_evalue_cap=1.0, aug_relaxed_evalue_mult=1.0, window_mmseqs_min_length=30)
    out = str(tmp_path / "aug.m8")
    ok, _det = isr.run_augmented_mmseqs_with_retries(variants_fasta=queries, region_fasta=genome, aug_hits_file=out,
                                                     aug_tmp_dir=str(tmp_path / "tmpw"), args=args, threads=1,
                                                     genome_name="Target.fna")
    assert ok
    assert {(r.split("\t")[0], r.split("\t")[1]) for r in open(out)} == EXPECTED


@needs_mmseqs
def test_the_annotation_transfer_names_what_the_fasta_names(tmp_path):
    import borrow_annotations as ba
    genome, queries = _inputs(tmp_path)
    hits = ba.mmseqs_hits(queries, genome)
    assert {(h["query"], h["chrom"]) for h in hits} == EXPECTED


def test_locate_gene_reads_the_headers_before_its_first_reader():
    nf = open(os.path.join(ROOT, "modules", "locate_gene.nf")).read()
    assert "evalue,bits,qheader,theader" in nf
    rewrite = nf.index("mmseqs_first_words.py --input hits_mmseqs.raw --output hits_mmseqs.m8")
    assert nf.index("mmseqs easy-search $gene $home_genome hits_mmseqs.raw") < rewrite
    for reader in ("cat hits_blast.txt hits_mmseqs.m8", "--hits hits_mmseqs.m8", "--mmseqs hits_mmseqs.m8"):
        assert rewrite < nf.index(reader), reader


def test_every_table_the_pipeline_parses_asks_for_the_headers():
    """No MMseqs2 call that the pipeline parses writes MMseqs2's own names alone."""
    for f in ("iterative_search_runner.py", "borrow_annotations.py"):
        src = open(os.path.join(BIN, f)).read()
        assert "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits\"" not in src, f
        assert "MMSEQS_FORMAT_OUTPUT" in src, f
