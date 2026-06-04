"""§1n — name-lookup home-locus resolver tests."""
import json
import os
import subprocess
import sys

import pytest

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import resolve_home_locus as rhl  # noqa: E402
from sequence_utils import translate  # noqa: E402

# Deterministic, stop-free codon per amino acid (back-translation for fixtures).
CODON = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA',
    'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTT', 'K': 'AAA',
    'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG',
    'Y': 'TAT', 'V': 'GTT',
}


def bt(aa: str) -> str:
    return "".join(CODON[a] for a in aa)


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


DCN_AA = "MAWKPGFDHKEFQYRSTVGI"        # 20 aa
LUM_AA = "VVVVVKKKKKDDDDDGGGGG"        # paralog-ish, dissimilar
MINUS_AA = "ITKPWAGFDHMEFQYRSTVL"      # minus-strand gene


@pytest.fixture
def home(tmp_path):
    """Build a synthetic home genome + GFF.

    chr1 layout (1-based):
      31-60   DCN exon1 (+)          151-210  LUM CDS (+)
      91-120  DCN exon2 (+)          251-310  MINUS CDS (gene on '-' strand)
    """
    dcn_e1 = bt(DCN_AA[:10])   # 30 bp
    dcn_e2 = bt(DCN_AA[10:])   # 30 bp
    lum = bt(LUM_AA)           # 60 bp
    minus_plus_strand = revcomp(bt(MINUS_AA))  # genome carries RC; gene is on '-'

    seq = list("A" * 360)
    def place(s, start1):  # start1 is 1-based inclusive
        seq[start1 - 1:start1 - 1 + len(s)] = list(s)
    place(dcn_e1, 31)
    place(dcn_e2, 91)
    place(lum, 151)
    place(minus_plus_strand, 251)
    genome_fa = tmp_path / "home.fna"
    genome_fa.write_text(">chr1\n" + "".join(seq) + "\n")

    gff = tmp_path / "home.gff"
    gff.write_text(
        "\n".join([
            "##gff-version 3",
            "chr1\tsrc\tgene\t31\t120\t.\t+\t.\tID=gene-DCN;Name=DCN;Alias=PG40,SLRR1B",
            "chr1\tsrc\tmRNA\t31\t120\t.\t+\t.\tID=rna-DCN;Parent=gene-DCN",
            "chr1\tsrc\tCDS\t31\t60\t.\t+\t0\tID=cds-DCN1;Parent=rna-DCN",
            "chr1\tsrc\tCDS\t91\t120\t.\t+\t0\tID=cds-DCN2;Parent=rna-DCN",
            "chr1\tsrc\tgene\t151\t210\t.\t+\t.\tID=gene-LUM;Name=LUM",
            "chr1\tsrc\tmRNA\t151\t210\t.\t+\t.\tID=rna-LUM;Parent=gene-LUM",
            "chr1\tsrc\tCDS\t151\t210\t.\t+\t0\tID=cds-LUM;Parent=rna-LUM",
            "chr1\tsrc\tgene\t251\t310\t.\t-\t.\tID=gene-MIN;Name=MINUSG",
            "chr1\tsrc\tmRNA\t251\t310\t.\t-\t.\tID=rna-MIN;Parent=gene-MIN",
            "chr1\tsrc\tCDS\t251\t310\t.\t-\t0\tID=cds-MIN;Parent=rna-MIN",
            "",
        ])
    )
    qdcn = tmp_path / "query_dcn.faa"
    qdcn.write_text(">DCN\n" + DCN_AA + "\n")
    return {"genome": str(genome_fa), "gff": str(gff), "query_dcn": str(qdcn), "dir": tmp_path}


def run(home, symbol, query=None, ident=60.0, pad=0):
    out_bed = home["dir"] / f"locus_{symbol or 'none'}.bed"
    out_status = home["dir"] / f"status_{symbol or 'none'}.json"
    rc = subprocess.run(
        [sys.executable, os.path.join(BIN, "resolve_home_locus.py"),
         "--query", query or home["query_dcn"],
         "--gene_symbol", symbol,
         "--home_gff", home["gff"],
         "--home_genome", home["genome"],
         "--out_bed", str(out_bed),
         "--out_status", str(out_status),
         "--min_consistency_identity", str(ident),
         "--pad", str(pad)],
        capture_output=True, text=True,
    )
    assert rc.returncode == 0, rc.stderr
    status = json.loads(out_status.read_text())
    bed = out_bed.read_text().strip()
    return status, bed


# ---- unit-level ----

def test_find_gene_features_by_name(home):
    feats = rhl.find_gene_features(home["gff"], "DCN")
    assert len(feats) == 1 and feats[0]["attributes"]["ID"] == "gene-DCN"


def test_find_gene_features_case_insensitive(home):
    assert len(rhl.find_gene_features(home["gff"], "dcn")) == 1


def test_find_gene_features_by_alias(home):
    assert len(rhl.find_gene_features(home["gff"], "SLRR1B")) == 1


def test_translate_cds_multi_exon(home):
    from sequence_utils import load_genome
    genome = load_genome(home["genome"])
    feats = rhl.find_gene_features(home["gff"], "DCN")
    cds_by_gene = rhl.collect_cds_by_gene(home["gff"], {"gene-DCN"})
    protein = rhl.translate_cds(cds_by_gene["gene-DCN"], genome, "+")
    assert protein == DCN_AA


def test_translate_cds_minus_strand(home):
    from sequence_utils import load_genome
    genome = load_genome(home["genome"])
    cds_by_gene = rhl.collect_cds_by_gene(home["gff"], {"gene-MIN"})
    protein = rhl.translate_cds(cds_by_gene["gene-MIN"], genome, "-")
    assert protein == MINUS_AA


# ---- end-to-end CLI ----

def test_matched_emits_single_locus_bed(home):
    status, bed = run(home, "DCN")
    assert status["status"] == "matched" and status["matched"] is True
    assert status["consistency_identity"] >= 99.0
    cols = bed.split("\t")
    assert cols[0] == "chr1"
    assert int(cols[1]) == 30 and int(cols[2]) == 120   # 0-based start, end
    assert cols[3] == "gene-DCN" and cols[5] == "+"


def test_matched_with_pad(home):
    status, bed = run(home, "DCN", pad=10)
    cols = bed.split("\t")
    assert int(cols[1]) == 20 and int(cols[2]) == 130


def test_not_in_gff_falls_back(home):
    status, bed = run(home, "ZZZZ")
    assert status["status"] == "not_in_gff" and status["matched"] is False
    assert bed == ""


def test_no_symbol_falls_back(home):
    status, bed = run(home, "")
    assert status["status"] == "no_symbol" and bed == ""


def test_failed_consistency_falls_back(home):
    # Look up LUM but pass the DCN query → paralog mismatch → must NOT anchor.
    status, bed = run(home, "LUM", query=home["query_dcn"], ident=60.0)
    assert status["status"] == "failed_consistency" and status["matched"] is False
    assert bed == ""
    assert status["consistency_identity"] < 60.0
    # but it DID locate the gene (for the diagnostic)
    assert status["gene_id"] == "gene-LUM" and status["chrom"] == "chr1"


def test_minus_strand_gene_matched(home):
    qmin = home["dir"] / "qmin.faa"
    qmin.write_text(">M\n" + MINUS_AA + "\n")
    status, bed = run(home, "MINUSG", query=str(qmin))
    assert status["status"] == "matched"
    cols = bed.split("\t")
    assert cols[5] == "-" and int(cols[1]) == 250 and int(cols[2]) == 310
