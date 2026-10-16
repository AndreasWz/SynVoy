"""The home gene by position uses the name path's check (batch of 2026-10-15, HG-3).

Until then the gene by position was the best shared-3-mer score above 30 among the genes within
100 kb of ANY home hit (tblastn reports to E 10): about 75-80 % identity at equal length, so a
short-isoform query or a query from another species fell to miniprot, and a close paralog could
be taken in silence where the query's own gene is not annotated. Now: windows from hits with
E <= search_evalue, the Smith-Waterman check of the name path (_query_match + _is_home_gene),
genes that overlap a hit first, a 0.05 margin over the runner-up gene, and the evidence in
goi_info.json. Synthetic genomes only."""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import annotate_goi_exons as ag  # noqa: E402

CODON = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA',
         'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
         'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT'}
INTRON = 'GT' + 'ACGTTGCA' * 30 + 'AG'
SPACER = 'TAACTAGTTGA' * 40
AA = 'ARNDCQEGHILKFPSTWYV'


def _protein(length, seed):
    rng = random.Random(seed)
    return 'M' + ''.join(rng.choice(AA) for _ in range(length - 1))


def _mutate(protein, identity, seed):
    rng = random.Random(seed)
    pos = rng.sample(range(1, len(protein)), int(round((1 - identity) * (len(protein) - 1))))
    out = list(protein)
    for i in pos:
        out[i] = rng.choice([a for a in AA if a != out[i]])
    return ''.join(out)


def _layout(tmp_path, genes, hits, gap_kb=20):
    """chrA holds ``genes`` [(gene_id, protein)] in this order, each as two exons, ``gap_kb``
    apart; chrB holds an unrelated gene. ``hits``: [(gene_index or 'B', exon, evalue)], a
    LOCATE-style hit row on that exon. Returns paths and exon spans."""
    gap = SPACER * max(1, (gap_kb * 1000) // len(SPACER))
    chrom, spans, rows = SPACER, {}, []
    for gid, prot in genes:
        half = len(prot) // 2
        e1 = (len(chrom), len(chrom) + 3 * half)
        chrom += ''.join(CODON[a] for a in prot[:half]) + INTRON
        e2 = (len(chrom), len(chrom) + 3 * (len(prot) - half) + 3)
        chrom += ''.join(CODON[a] for a in prot[half:]) + 'TAA' + gap
        spans[gid] = (e1, e2)
    other = _protein(150, seed=99)
    chr_b = SPACER + ''.join(CODON[a] for a in other) + 'TAA' + SPACER
    o_span = (len(SPACER), len(SPACER) + 3 * len(other) + 3)
    fasta = tmp_path / "home.fna"
    fasta.write_text(f">chrA\n{chrom}\n>chrB\n{chr_b}\n")

    def row(c, kind, span, attrs, phase='.'):
        return f"{c}\tsrc\t{kind}\t{span[0] + 1}\t{span[1]}\t.\t+\t{phase}\t{attrs}\n"

    lines = ["##gff-version 3\n"]
    for gid, (e1, e2) in spans.items():
        lines.append(row('chrA', 'gene', (e1[0], e2[1]), f'ID={gid};Name={gid}'))
        lines.append(row('chrA', 'mRNA', (e1[0], e2[1]), f'ID=rna-{gid};Parent={gid}'))
        for e in (e1, e2):
            lines.append(row('chrA', 'CDS', e, f'ID=cds-{gid};Parent=rna-{gid}', phase='0'))
    lines.append(row('chrB', 'gene', o_span, 'ID=gene-OTHER;Name=OTHER'))
    lines.append(row('chrB', 'mRNA', o_span, 'ID=rna-OTHER;Parent=gene-OTHER'))
    lines.append(row('chrB', 'CDS', o_span, 'ID=cds-OTHER;Parent=rna-OTHER', phase='0'))
    gff = tmp_path / "home.gff"
    gff.write_text(''.join(lines))
    with open(tmp_path / "hits.tsv", "w") as fh:
        for where, exon, ev in hits:
            if where == 'B':
                c, (a, b) = 'chrB', o_span
            else:
                c, (a, b) = 'chrA', spans[genes[where][0]][exon]
            fh.write(f"q\t{c}\t90.0\t40\t2\t0\t1\t40\t{a + 1}\t{b}\t{ev}\t80\n")
    return dict(fasta=str(fasta), gff=str(gff), hits=str(tmp_path / "hits.tsv"), spans=spans)


def _find(lay, query, max_evalue=0.01, window=100000):
    regions = ag._parse_hit_regions(lay['hits'], None)
    return ag._search_gff_by_region_and_sequence(lay['gff'], lay['fasta'], query, regions,
                                                 window, max_evalue=max_evalue)


@pytest.fixture(autouse=True)
def _need_parasail():
    pytest.importorskip("parasail")


def test_a_half_length_isoform_query_is_found(tmp_path):
    """(a) The query is the first half of an annotated gene: the 3-mer rule scored it about 25
    (below 30) and fell to miniprot; the name path's check accepts it."""
    gene = _protein(240, seed=1)
    lay = _layout(tmp_path, [("gene-G", gene)], [(0, 0, 1e-30)])
    query = gene[:120]
    assert ag._quick_protein_identity(query, gene) < 30        # the old rule fails it
    m = _find(lay, query)
    assert m and m['gene_id'] == 'gene-G' and m['found_by'] == 'position'
    assert m['overlaps_hit'] is True and m['transcript_id'] == 'rna-gene-G'
    assert m['identity'] == pytest.approx(100.0) and m['home_gene_ambiguous'] is False


def test_a_query_from_another_species_is_found(tmp_path, capsys):
    """(c) 70 % identity: accepted by the check, and named in a warning (not the query itself)."""
    gene = _protein(200, seed=2)
    query = _mutate(gene, 0.70, seed=3)
    lay = _layout(tmp_path, [("gene-G", gene)], [(0, 0, 1e-15), (0, 1, 1e-12)])
    m = _find(lay, query)
    assert m and m['gene_id'] == 'gene-G'
    assert 60 <= m['identity'] < 80 and m['score_fraction'] >= 0.5
    assert "is not the query itself" in capsys.readouterr().out


def test_a_paralog_where_the_own_gene_is_missing_is_not_taken_in_silence(tmp_path, capsys):
    """(b) The query's own gene is not annotated; an 85 % paralog 20 kb away is the only gene
    that agrees: taken, but the run says it is not the query itself."""
    query = _protein(200, seed=4)
    lay = _layout(tmp_path, [("gene-UNREL", _protein(200, seed=5)),
                             ("gene-PARA", _mutate(query, 0.85, seed=6))],
                  [(0, 0, 1e-20), (1, 0, 1e-10)])
    m = _find(lay, query)
    out = capsys.readouterr().out
    assert m is None or m['home_gene_ambiguous'] or "is not the query itself" in out
    assert m['gene_id'] == 'gene-PARA' and "is not the query itself" in out


def test_two_close_relatives_flag_the_choice(tmp_path, capsys):
    """The knock-out case: two relatives agree with the query within 0.05 of score fraction:
    the better one is taken, home_gene_ambiguous is set, and the runner-up is recorded."""
    query = _protein(200, seed=7)
    a, b = _mutate(query, 0.82, seed=8), _mutate(query, 0.80, seed=9)
    lay = _layout(tmp_path, [("gene-A", a), ("gene-B", b)], [(0, 0, 1e-12), (1, 0, 1e-11)])
    m = _find(lay, query)
    assert m['home_gene_ambiguous'] is True
    assert {m['gene_id'], m['runner_up_gene_id']} == {"gene-A", "gene-B"}
    assert m['score_fraction'] - m['runner_up_score_fraction'] < ag.HOME_GENE_POSITION_MARGIN
    assert "the home gene is ambiguous" in capsys.readouterr().out
    # The query's own gene annotated next to them: taken, with a clear margin.
    lay2 = _layout(tmp_path, [("gene-A", a), ("gene-Q", query), ("gene-B", b)],
                   [(0, 0, 1e-12), (1, 0, 1e-30), (2, 0, 1e-11)])
    m2 = _find(lay2, query)
    assert m2['gene_id'] == 'gene-Q' and m2['home_gene_ambiguous'] is False


def test_only_significant_hits_make_windows(tmp_path):
    """A hit at E 5 next to an unrelated gene made it a candidate; with search_evalue it is not."""
    gene = _protein(200, seed=10)
    lay = _layout(tmp_path, [("gene-G", gene)], [('B', 0, 5.0)])
    assert _find(lay, gene, max_evalue=0.01) is None
    assert _find(lay, gene, max_evalue=10) is None            # unrelated: rejected by the check


def test_a_gene_that_overlaps_a_hit_comes_first(tmp_path):
    """A nearby gene that agrees slightly better loses to the gene the hits fall on."""
    query = _protein(200, seed=11)
    on_hit, near = _mutate(query, 0.90, seed=12), _mutate(query, 0.97, seed=13)
    lay = _layout(tmp_path, [("gene-ON", on_hit), ("gene-NEAR", near)], [(0, 0, 1e-20), (0, 1, 1e-18)])
    m = _find(lay, query)
    assert m['gene_id'] == 'gene-ON' and m['overlaps_hit'] is True
    assert m['runner_up_gene_id'] == 'gene-NEAR'


def test_main_writes_the_evidence_and_keeps_the_method(tmp_path, monkeypatch):
    gene = _protein(240, seed=14)
    lay = _layout(tmp_path, [("gene-G", gene)], [(0, 0, 1e-30)])
    (tmp_path / "q.faa").write_text(f">q\n{gene[:120]}\n")
    monkeypatch.setattr(ag, "fetch_uniprot_names", lambda _q: [])
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(sys, "argv", ["annotate_goi_exons.py", "--query", str(tmp_path / "q.faa"),
                                      "--genome", lay['fasta'], "--gff", lay['gff'],
                                      "--blast_hits", lay['hits'], "--query_id", "NOSUCHNAME",
                                      "--max_evalue", "0.01",
                                      "--output_exons", str(out / "e.faa"), "--output_bed", str(out / "e.bed"),
                                      "--output_info", str(out / "info.json")])
    ag.main()
    info = json.loads((out / "info.json").read_text())
    assert info['method'] == 'gff_annotation'             # split_loci keys the primary locus on it
    assert info['found_by'] == 'position' and info['gene_id'] == 'gene-G'
    assert info['transcript_id'] == 'rna-gene-G' and info['home_gene_ambiguous'] is False
    for key in ('identity', 'score_fraction', 'overlaps_hit', 'runner_up_gene_id', 'runner_up_score_fraction'):
        assert key in info


def test_the_name_path_records_its_evidence(tmp_path):
    gene = _protein(200, seed=15)
    lay = _layout(tmp_path, [("gene-G", gene), ("gene-P", _mutate(gene, 0.75, seed=16))], [(0, 0, 1e-30)])
    m = ag._search_gff_by_name(lay['gff'], ['gene-G', 'gene-P'], query_seq=gene, genome_file=lay['fasta'])
    assert m['gene_id'] == 'gene-G' and m['found_by'] == 'name'
    assert m['transcript_id'] == 'rna-gene-G' and m['identity'] == pytest.approx(100.0)
    assert m['runner_up_gene_id'] == 'gene-P' and m['runner_up_score_fraction'] < m['score_fraction']


def test_without_parasail_the_old_rule_is_called_explicitly(tmp_path, monkeypatch, capsys):
    import builtins
    real_import = builtins.__import__

    def no_parasail(name, *a, **k):
        if name == "parasail":
            raise ImportError("no parasail in this test")
        return real_import(name, *a, **k)
    gene = _protein(200, seed=17)
    lay = _layout(tmp_path, [("gene-G", gene)], [(0, 0, 1e-30)])
    monkeypatch.setattr(builtins, "__import__", no_parasail)
    m = _find(lay, gene)
    assert m['gene_id'] == 'gene-G' and m['found_by'] == 'position_3mer'
    assert "the rule before 2026-10-15" in capsys.readouterr().out
