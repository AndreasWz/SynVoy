"""The home gene taken from the annotation must be the query's own gene.

Found 2026-10-08 on two human queries:

* Q8N1Y9 has no gene symbol at UniProt, only the protein name "Putative uncharacterized
  protein FLJ37218". The name lookup split that name into words; "uncharacterized" matched
  12,484 human genes and the first one in the file, a non-coding gene on another chromosome,
  became the home gene. The run searched ten genomes for an 11-residue translation of that
  gene's span.
* TP53 was read from the first transcript in the file, the 354-residue isoform, although the
  query is the 393-residue protein. The protein of the home gene model replaces the query in
  the query set, so the run searched with the isoform.

A gene found by name (or by position) is now compared with the query, and the transcript
that agrees best with the query is the one used.
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import annotate_goi_exons as ag  # noqa: E402
from sequence_utils import load_genome_subset, reverse_complement  # noqa: E402

CODON = {  # one fixed codon per residue, none of them a stop
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA',
    'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
    'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
}
STOP = 'TAA'
INTRON = 'GT' + 'ACGTTGCA' * 30 + 'AG'      # stops in every frame, so no read-through
SPACER = 'TAACTAGTTGA' * 40                 # intergenic, stops in every frame


def _protein(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return 'M' + ''.join(rng.choice('ARNDCQEGHILKFPSTWYV') for _ in range(length - 1))


def _dna(protein: str) -> str:
    return ''.join(CODON[a] for a in protein)


def _layout(tmp_path: Path, minus: bool = False):
    """A genome with one coding gene of two transcripts, and two decoys.

    chrA   ALPHA1: three coding exons. Transcript t_short comes FIRST in the file and
           holds the last exon only; t_full holds all three.
    chrB   a non-coding gene described as "uncharacterized LOC900", first in the file;
           and an unrelated coding gene whose description contains "uncharacterized".
    """
    full = _protein(120, seed=1)
    e1, e2, e3 = full[:40], full[40:80], full[80:]
    other = _protein(110, seed=2)

    gene_dna = _dna(e1) + INTRON + _dna(e2) + INTRON + _dna(e3) + STOP
    chr_a = SPACER + gene_dna + SPACER
    g_start = len(SPACER)                          # 0-based
    x1 = (g_start, g_start + 120)
    x2 = (x1[1] + len(INTRON), x1[1] + len(INTRON) + 120)
    x3 = (x2[1] + len(INTRON), x2[1] + len(INTRON) + 120 + 3)
    strand = '+'
    if minus:
        total = len(chr_a)
        chr_a = reverse_complement(chr_a)
        x1, x2, x3 = [(total - e, total - s) for s, e in (x1, x2, x3)]
        strand = '-'
    gene_span = (min(x1[0], x3[0]), max(x1[1], x3[1]))

    noncoding = 'ACGGTCA' * 60
    other_dna = _dna(other) + STOP
    chr_b = SPACER + noncoding + SPACER + other_dna + SPACER
    n_span = (len(SPACER), len(SPACER) + len(noncoding))
    o_start = n_span[1] + len(SPACER)
    o_span = (o_start, o_start + len(other_dna))

    fasta = tmp_path / "home.fna"
    fasta.write_text(f">chrB\n{chr_b}\n>chrA\n{chr_a}\n")

    def row(chrom, kind, span, strand_, attrs, phase='.'):
        return f"{chrom}\tsrc\t{kind}\t{span[0] + 1}\t{span[1]}\t.\t{strand_}\t{phase}\t{attrs}\n"

    lines = ["##gff-version 3\n"]
    # chrB first: the decoys precede the real gene in the file
    lines.append(row('chrB', 'gene', n_span, '+',
                     'ID=gene-LOC900;Name=LOC900;description=uncharacterized LOC900;gene_biotype=lncRNA'))
    lines.append(row('chrB', 'gene', o_span, '+',
                     'ID=gene-LOC901;Name=LOC901;description=uncharacterized protein LOC901'))
    lines.append(row('chrB', 'mRNA', o_span, '+', 'ID=rna-o1;Parent=gene-LOC901'))
    lines.append(row('chrB', 'CDS', o_span, '+', 'ID=cds-o1;Parent=rna-o1', phase='0'))
    lines.append(row('chrA', 'gene', gene_span, strand,
                     'ID=gene-ALPHA1;Name=ALPHA1;description=alpha one protein'))
    # the short isoform FIRST
    lines.append(row('chrA', 'mRNA', x3, strand, 'ID=rna-a_short;Parent=gene-ALPHA1'))
    lines.append(row('chrA', 'CDS', x3, strand, 'ID=cds-a_short;Parent=rna-a_short', phase='0'))
    lines.append(row('chrA', 'mRNA', gene_span, strand, 'ID=rna-b_full;Parent=gene-ALPHA1'))
    for x in sorted((x1, x2, x3)):
        lines.append(row('chrA', 'CDS', x, strand, 'ID=cds-b_full;Parent=rna-b_full', phase='0'))
    gff = tmp_path / "home.gff"
    gff.write_text(''.join(lines))

    # one hit row inside the gene, in the 12-column format of LOCATE_GENE
    hits = tmp_path / "hits.tsv"
    hs, he = sorted(x2)
    tcoords = (hs + 1, he) if strand == '+' else (he, hs + 1)
    hits.write_text(f"q\tchrA\t100\t40\t0\t0\t41\t80\t{tcoords[0]}\t{tcoords[1]}\t1e-20\t80\n")
    return dict(fasta=str(fasta), gff=str(gff), hits=str(hits), full=full, short=e3,
                other=other, exons=sorted((x1, x2, x3)), strand=strand)


def _model_protein(match, lay):
    seq = load_genome_subset(lay['fasta'], {match['chrom']})[match['chrom']]
    return ag._cds_protein(match['cds_parts'], seq, match['strand'])


@pytest.mark.parametrize("minus", [False, True])
def test_exact_name_takes_the_transcript_that_is_the_query(tmp_path, minus):
    lay = _layout(tmp_path, minus=minus)
    match = ag._search_gff_by_name(lay['gff'], ['ALPHA1'], query_seq=lay['full'],
                                   genome_file=lay['fasta'])
    assert match is not None and match['gene_id'] == 'gene-ALPHA1'
    assert len(match['cds_parts']) == 3, "the first transcript in the file is the short isoform"
    assert _model_protein(match, lay) == lay['full']


def test_short_isoform_query_takes_the_short_transcript(tmp_path):
    """The choice follows the query, not the transcript length."""
    lay = _layout(tmp_path)
    match = ag._search_gff_by_name(lay['gff'], ['ALPHA1'], query_seq=lay['short'],
                                   genome_file=lay['fasta'])
    assert match is not None and len(match['cds_parts']) == 1
    assert _model_protein(match, lay) == lay['short']


def test_a_word_of_the_protein_name_does_not_select_a_gene(tmp_path):
    """The Q8N1Y9 case: only a descriptive protein name, whose words occur in the
    description of unrelated genes. No gene may be taken."""
    lay = _layout(tmp_path)
    names = ['Putative uncharacterized protein FLJ37218', 'Q8N1Y9']
    assert ag._search_gff_by_name(lay['gff'], names, query_seq=lay['full'],
                                  genome_file=lay['fasta']) is None


def test_word_match_is_still_taken_without_a_query(tmp_path):
    """Without a query there is nothing to compare with; the old behaviour is kept for
    callers that pass names only. The pipeline always passes the query."""
    lay = _layout(tmp_path)
    match = ag._search_gff_by_name(lay['gff'], ['Putative uncharacterized protein FLJ37218'])
    assert match is not None and match['gene_id'] == 'gene-LOC900'


def test_wrong_gene_with_the_right_name_is_rejected(tmp_path):
    """A name can be right and the gene still wrong (other species, reused symbol)."""
    lay = _layout(tmp_path)
    assert ag._search_gff_by_name(lay['gff'], ['ALPHA1'], query_seq=lay['other'],
                                  genome_file=lay['fasta']) is None


def test_match_falls_back_to_the_gene_at_the_hits(tmp_path):
    """Name lookup fails on the decoys, the positional lookup then finds the gene and
    again takes the transcript that is the query."""
    lay = _layout(tmp_path)
    match = ag.match_goi_in_gff('FLJ37218 uncharacterized', lay['full'], lay['gff'],
                                lay['fasta'], lay['hits'], None, 1000)
    assert match is not None and match['gene_id'] == 'gene-ALPHA1'
    assert len(match['cds_parts']) == 3
    assert _model_protein(match, lay) == lay['full']


def test_gene_row_without_cds_passes_only_if_its_span_is_the_query(tmp_path):
    protein = _protein(90, seed=5)
    chrom = SPACER + _dna(protein) + STOP + SPACER
    span = (len(SPACER), len(SPACER) + 3 * len(protein) + 3)
    fasta = tmp_path / "g.fna"
    fasta.write_text(f">c1\n{chrom}\n")
    gff = tmp_path / "g.gff"
    gff.write_text("##gff-version 3\n"
                   f"c1\tsrc\tgene\t{span[0] + 1}\t{span[1]}\t.\t+\t.\tID=gene-ORF1;Name=ORF1\n")
    match = ag._search_gff_by_name(str(gff), ['ORF1'], query_seq=protein, genome_file=str(fasta))
    assert match is not None and match['gene_id'] == 'gene-ORF1'
    assert ag._search_gff_by_name(str(gff), ['ORF1'], query_seq=_protein(90, seed=6),
                                  genome_file=str(fasta)) is None


def test_a_few_residues_of_the_query_are_not_the_gene():
    """The score is taken as a fraction of the QUERY's self-score. Taken of the shorter
    sequence, a six-residue junk translation scores 1.0 against any query containing it."""
    pytest.importorskip("parasail")   # without it the step falls back to a word-sharing measure
    query = _protein(200, seed=9)
    frac, ident = ag._query_match(query, query[50:56])
    assert ident == pytest.approx(100.0)
    assert frac < 0.1
    assert not ag._is_home_gene(frac, ident)
    assert ag._is_home_gene(*ag._query_match(query, query))
    # an isoform lacking a tenth of the query is still the gene
    assert ag._is_home_gene(*ag._query_match(query, query[20:]))


def test_load_genome_subset_keeps_only_what_is_asked(tmp_path):
    fasta = tmp_path / "x.fna"
    fasta.write_text(">a desc\nACGT\nAC\n>b\nGG\n>c\nTT\n")
    assert load_genome_subset(str(fasta), {'b'}) == {'b': 'GG'}
    assert load_genome_subset(str(fasta), {'a', 'c'}) == {'a': 'ACGTAC', 'c': 'TT'}
    assert load_genome_subset(str(fasta), set()) == {}
    assert load_genome_subset(str(fasta), {'missing'}) == {}
