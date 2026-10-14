"""Anchor selection at the home locus (bin/extract_flanking_genes.py, 2026-10-08/09).

* A neighbour is skipped as an anchor only when it is a RELATIVE of the gene of interest
  (pairwise alignment E-value), not when it merely has a similar length.
* Every anchor protein carries its home coordinates and the home genome length in the
  FASTA header; the search derives its distances from them.
* The marker row of the gene of interest in the home BED is not a query.
"""
import os
import random
import subprocess
import sys

import pytest

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))
sys.path.insert(0, BIN)

pytest.importorskip("parasail")

import extract_flanking_genes as efg  # noqa: E402
import synvoy_scale  # noqa: E402

SCRIPT = os.path.join(BIN, "extract_flanking_genes.py")
AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
         "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
         "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


def _protein(n, seed):
    rnd = random.Random(seed)
    return "M" + "".join(rnd.choice(AA) for _ in range(n - 1))


def _mutate(seq, share, seed):
    rnd = random.Random(seed)
    out = list(seq)
    for i in rnd.sample(range(1, len(seq)), int(len(seq) * share)):
        out[i] = rnd.choice(AA.replace(out[i], ""))
    return "".join(out)


# ------------------------------------------------------------------ the relatedness test

GOI = _protein(240, 1)


def test_relative_of_the_goi_is_recognised():
    assert efg.is_goi_relative(_mutate(GOI, 0.55, 2), [GOI], 1e-5)


def test_similar_length_alone_does_not_make_a_relative():
    """The measure that was replaced flagged about half of such proteins (it scores two
    unrelated proteins of equal length near its own threshold of 35)."""
    unrelated = [_protein(240 + (i % 7) - 3, 50 + i) for i in range(40)]
    assert not any(efg.is_goi_relative(p, [GOI], 1e-5) for p in unrelated)
    old_measure = sum(efg._lcs_similarity(p, GOI) >= 35.0 for p in unrelated)
    assert old_measure >= 10


def test_a_looser_evalue_reaches_more_distant_relatives():
    distant = _mutate(GOI, 0.80, 15)
    e = efg.sw_pair_evalue(distant, GOI)[0]
    assert 1e-5 < e < 1e-3, e                    # between the default and the tandem preset
    assert not efg.is_goi_relative(distant, [GOI], 1e-5)
    assert efg.is_goi_relative(distant, [GOI], 1e-3)


def test_zero_switches_the_test_off_and_empty_input_is_no_relative():
    assert not efg.is_goi_relative(GOI, [GOI], 0)
    assert not efg.is_goi_relative("", [GOI], 1e-5)
    assert not efg.is_goi_relative(GOI, [], 1e-5)


def test_any_of_several_goi_sequences_counts():
    other = _protein(200, 9)
    assert efg.is_goi_relative(_mutate(other, 0.3, 10), [GOI, other], 1e-5)


def test_without_an_aligner_the_old_measure_is_used_with_a_warning(monkeypatch, capsys):
    def no_aligner(_a, _b):
        raise RuntimeError("parasail-python is required")
    monkeypatch.setattr(efg, "sw_pair_evalue", no_aligner)
    monkeypatch.setattr(efg.is_goi_relative, "_warned", False, raising=False)
    assert efg.is_goi_relative(GOI, [GOI], 1e-5)
    assert not efg.is_goi_relative("MKT" * 30, [GOI], 1e-5)
    assert capsys.readouterr().err.count("falls back to the common-subsequence measure") == 1


# ------------------------------------------------------------------------- end to end

def _locus(tmp_path, strand_of_goi="+"):
    """chr1 with nine single-exon genes, the gene of interest in the middle.

    Gene 3 is a relative of the gene of interest; gene 5 has its length and nothing else."""
    rnd = random.Random(7)
    proteins = [_protein(180 + 15 * i, 100 + i) for i in range(9)]
    proteins[4] = GOI
    proteins[3] = _mutate(GOI, 0.5, 11)
    proteins[5] = _protein(len(GOI), 12)
    comp = str.maketrans("ACGT", "TGCA")
    seq, genes = [], []
    pos = 0
    for i, prot in enumerate(proteins):
        spacer = "".join(rnd.choice("ACGT") for _ in range(3000))
        seq.append(spacer)
        pos += len(spacer)
        cds = "".join(CODON[a] for a in prot) + "TAA"
        strand = strand_of_goi if i == 4 else "+"
        seq.append(cds if strand == "+" else cds.translate(comp)[::-1])
        genes.append((f"G{i}", pos, pos + len(cds), strand))
        pos += len(cds)
    seq.append("".join(rnd.choice("ACGT") for _ in range(3000)))
    genome = "".join(seq)
    (tmp_path / "home.fna").write_text(">chr1\n" + genome + "\n>chrUn\n" + "ACGT" * 250 + "\n")
    rows = ["##gff-version 3"]
    for name, start, end, strand in genes:
        rows.append(f"chr1\tt\tgene\t{start + 1}\t{end}\t.\t{strand}\t.\tID=gene-{name};Name={name};gene_biotype=protein_coding")
        rows.append(f"chr1\tt\tmRNA\t{start + 1}\t{end}\t.\t{strand}\t.\tID=rna-{name};Parent=gene-{name}")
        rows.append(f"chr1\tt\tCDS\t{start + 1}\t{end}\t.\t{strand}\t0\tID=cds-{name};Parent=rna-{name}")
    (tmp_path / "home.gff").write_text("\n".join(rows) + "\n")
    goi = genes[4]
    (tmp_path / "locus.bed").write_text(f"chr1\t{goi[1]}\t{goi[2]}\tgene_loc\t0\t{goi[3]}\t500\n")
    (tmp_path / "goi.faa").write_text(f">GOI_X\n{GOI}\n")
    return genes, len(genome) + 1000


def _run(tmp_path, *extra):
    out_bed, out_faa = tmp_path / "out.bed", tmp_path / "out.faa"
    r = subprocess.run([sys.executable, SCRIPT, "--bed", str(tmp_path / "locus.bed"),
                        "--gff", str(tmp_path / "home.gff"), "--genome", str(tmp_path / "home.fna"),
                        "--n_flank", "3", "--min_size", "300", "--exon_mode", "false",
                        "--goi_faa", str(tmp_path / "goi.faa"),
                        "--out_bed", str(out_bed), "--out_faa", str(out_faa), *extra],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    headers = [line[1:].strip() for line in out_faa.read_text().splitlines() if line.startswith(">")]
    bed = [line.split("\t") for line in out_bed.read_text().splitlines() if line.strip()]
    return headers, bed, r.stdout


def _ids(headers):
    return [h.split()[0] for h in headers]


def test_relative_is_skipped_and_the_walk_goes_one_gene_further(tmp_path):
    _locus(tmp_path)
    headers, _bed, _out = _run(tmp_path)
    ids = _ids(headers)
    assert "gene-G3" not in ids                         # the relative
    assert "gene-G5" in ids                             # same length, unrelated: an anchor
    assert ids == ["gene-G0", "gene-G1", "gene-G2", "gene-G5", "gene-G6", "gene-G7"]


def test_with_the_test_switched_off_the_relative_is_an_anchor(tmp_path):
    _locus(tmp_path)
    headers, _bed, _out = _run(tmp_path, "--goi_relative_evalue", "0")
    assert _ids(headers) == ["gene-G1", "gene-G2", "gene-G3", "gene-G5", "gene-G6", "gene-G7"]


def test_the_gene_of_interest_is_never_its_own_neighbour(tmp_path):
    """Switching the neighbour test off must not bring the gene at the GOI's own position
    back as an anchor: it would count as one flanking gene next to every true call."""
    genes, _bp = _locus(tmp_path)
    for extra in ((), ("--goi_relative_evalue", "0"), ("--goi_relative_evalue", "1e-300")):
        headers, bed, _out = _run(tmp_path, *extra)
        assert "gene-G4" not in _ids(headers), extra
        assert [row[3] for row in bed if row[3].startswith("GOI_")] == [f"GOI_chr1_{genes[4][1]}"]


def test_nearest_gene_is_an_anchor_when_the_locus_is_between_genes(tmp_path):
    """A locus that overlaps no gene (an unannotated ORF seed) sits between G6 and G7; G6 is
    the nearest gene. It used to be skipped like the gene of interest (bug hunt A7)."""
    genes, _bp = _locus(tmp_path)
    g6_end = genes[6][2]
    (tmp_path / "locus.bed").write_text(f"chr1\t{g6_end + 1000}\t{g6_end + 1300}\torf_loc\t0\t+\t500\n")
    headers, _bed, _out = _run(tmp_path)
    ids = _ids(headers)
    assert "gene-G6" in ids
    assert ids == ["gene-G2", "gene-G5", "gene-G6", "gene-G7", "gene-G8"]


def test_anchor_headers_carry_home_coordinates_and_genome_length(tmp_path):
    genes, genome_bp = _locus(tmp_path)
    headers, _bed, _out = _run(tmp_path)
    by_name = {f"gene-{name}": (start, end) for name, start, end, _s in genes}
    for h in headers:
        gene = h.split()[0]
        assert f"home={by_name[gene][0]}-{by_name[gene][1]}" in h.split()
        assert f"genome={genome_bp}" in h.split()

    geometry = synvoy_scale.home_geometry(str(tmp_path / "out.faa"))
    assert geometry["n"] == 6 and geometry["home_genome_bp"] == genome_bp
    assert geometry["order"] == ["gene-G0", "gene-G1", "gene-G2", "gene-G5", "gene-G6", "gene-G7"]
    assert geometry["span"] == by_name["gene-G7"][1] - by_name["gene-G0"][0]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_goi_marker_stays_in_the_bed_and_is_not_a_query(tmp_path, strand):
    """On the plus strand the translation of the marker interval is the gene once more:
    as a FASTA record it was a second query for the gene of interest."""
    genes, _bp = _locus(tmp_path, strand_of_goi=strand)
    headers, bed, out = _run(tmp_path)
    marker = f"GOI_chr1_{genes[4][1]}"
    assert "Injected GOI pseudo-gene" in out
    assert [row[3] for row in bed if row[3].startswith("GOI_")] == [marker]
    assert not any(h.startswith("GOI_") for h in headers)
    assert len(bed) == len(headers) + 1
