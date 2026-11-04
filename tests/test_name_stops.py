"""Fix-phase d4 (2 Nov; the loud halves of items 29, 30, 31, 32 and 36): names that lost results in silence stop the run.

- 29: a flanking gene id, or a home sequence name, that MMseqs2 reports as another name (lcl|X as X, gnl|db|X as X, ucon_5 as
  on_5 ...): the flanking gene's hits were filed under a name nothing else used and its model was lost; a home sequence gave two
  loci at one place. Since item 25a a target named so stops the search; now these stop too. Making them work waits.
- 30: a '|' in a home sequence name when the home genes are predicted (no GFF): their ids hold the name and are read up to it.
- 31: a home FASTA and GFF that share no seqid stop; a partial overlap is a report warning (the per-row skipping stays: the
  easy-mode chromosome filter removes records from the FASTA alone, and those are not counted).
- 32: a gene id with whitespace (ID=ORF%201 reads "ORF 1") at the three writers of gene records (gff_to_faa.py,
  extract_flanking_genes.py, borrow_annotations.py): "ORF 1" and "ORF 2" read as one id.
- 36: a query id holding |exon_, |frag_ or |slide_ stops in normalize_query.py: the search had no full-length query.
"""
import json
import os
import random
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import annotate_goi_exons as ag  # noqa: E402
import sequence_utils as su  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT'}


# --------------------------------------------------------------------------- the rule

@pytest.mark.parametrize("word,reported", [
    ("lcl|X", "X"), ("ref|NC_1.1|", "NC_1.1"), ("sp|P12345|Y_HUMAN", "P12345"), ("tr|A0A1|B", "A0A1"),
    ("gb|X1.1|", "X1.1"), ("pdb|1ABC|A", "1ABC"), ("bbs|123", "123"), ("cl|X", "X"), ("gnl|db|X", "X"),
    ("pat|US|X", "X"), ("pir||X", "X"), ("prf||X", "X"), ("gi|1|ref|X|", "X"), ("gi|2", ""),
    ("ucon_5", "on_5"), ("uc001aaa.3", "001aaa.3"), ("consensus_24", "24"),
    ("emb|X|", "emb|X|"), ("dbj|X|", "dbj|X|"), ("NC_000012.12", "NC_000012.12"), ("ctg1|arrow", "ctg1|arrow")])
def test_the_names_mmseqs2_reports(word, reported):
    """Item 25's measurement on MMseqs2 18-8cc5c."""
    assert su.mmseqs_reported_name(word) == reported


@pytest.mark.parametrize("name", ["gene-LOC408474", "NC_000012.12", "NT_033777.3", "FBgn0000508", "AT1G01010",
                                  "orf19.1562", "gene-orf19.1562", "NODE_1_length_5000_cov_1.0", "Ca22chr1A",
                                  "SUPER_10_unloc_1", "scaffold_12"])
def test_ordinary_names_pass(name):
    assert su.id_problem(name, pipe=True, mmseqs=True) is None


def test_each_class_is_named():
    assert "whitespace" in su.id_problem("ORF 1")
    assert "'|'" in su.id_problem("ctg|arrow", pipe=True, mmseqs=False)
    assert su.id_problem("ctg|arrow") is None               # '|' alone is no MMseqs2 rewrite
    assert "MMseqs2 as 'X'" in su.id_problem("lcl|X")
    assert su.id_problem("lcl|X", mmseqs=False) is None


def test_check_gene_ids_names_the_id_and_both_files():
    with pytest.raises(su.StepFailed) as stop:
        su.check_gene_ids(["gene-a", "ORF 1"], "home.gff", "home_proteome.faa")
    msg = str(stop.value)
    assert "'ORF 1'" in msg and "home.gff" in msg and "home_proteome.faa" in msg and "Rename" in msg
    su.check_gene_ids(["lcl|X"], "home.gff", "out.faa")       # MMseqs2 only when asked
    with pytest.raises(su.StepFailed):
        su.check_gene_ids(["lcl|X"], "home.gff", "out.faa", mmseqs=True)


# --------------------------------------------------------------------------- the home genome (ANNOTATE_GOI)

def _home(tmp_path, names, gff_seqids=None):
    fna = tmp_path / "home.fna"
    fna.write_text("".join(f">{n} some text\nACGTACGTAC\n" for n in names))
    gff = None
    if gff_seqids is not None:
        gff = tmp_path / "home.gff"
        gff.write_text("##gff-version 3\n" + "".join(
            f"{s}\tsrc\tregion\t1\t10\t.\t+\t.\tID={s}\n{s}\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-{i}\n"
            for i, s in enumerate(gff_seqids)))
    return str(fna), (str(gff) if gff else None)


def test_a_home_name_mmseqs2_rewrites_stops(tmp_path):
    fna, gff = _home(tmp_path, ["chr1", "lcl|chr2"], ["chr1", "lcl|chr2"])
    with pytest.raises(su.StepFailed) as stop:
        ag.check_home_names(fna, gff)
    assert "'lcl|chr2'" in str(stop.value) and fna in str(stop.value) and "Rename it" in str(stop.value)


def test_a_pipe_in_a_home_name_stops_only_without_a_gff(tmp_path):
    fna, gff = _home(tmp_path, ["000000F|arrow", "000001F|arrow"], ["000000F|arrow", "000001F|arrow"])
    assert ag.check_home_names(fna, gff) is None
    with pytest.raises(su.StepFailed) as stop:
        ag.check_home_names(fna, None)
    assert "'000000F|arrow'" in str(stop.value) and "predicted" in str(stop.value)


def test_a_gff_that_shares_no_seqid_stops(tmp_path):
    fna, gff = _home(tmp_path, ["NC_1.1|foo", "NC_2.1|foo"], ["NC_1.1", "NC_2.1"])
    with pytest.raises(su.StepFailed) as stop:
        ag.check_home_names(fna, gff)
    msg = str(stop.value)
    assert "share no sequence name" in msg and "NC_1.1" in msg and "NC_1.1|foo" in msg


def test_a_partial_overlap_is_counted_not_stopped(tmp_path, capsys):
    fna, gff = _home(tmp_path, ["LG_1", "LG_2_pilon"], ["LG_1", "LG_2"])
    assert ag.check_home_names(fna, gff) == {"seqids": 1, "of": 2, "examples": ["LG_2"]}
    assert "1 of the home GFF's 2 sequences with genes are not in the home FASTA" in capsys.readouterr().err


def test_records_the_chromosome_filter_removed_are_not_counted(tmp_path):
    fna, gff = _home(tmp_path, ["NC_1.1"], ["NC_1.1", "NW_7.1"])
    assert ag.check_home_names(fna, gff, removed_names=["NW_7.1"]) is None
    assert ag.check_home_names(fna, gff)["examples"] == ["NW_7.1"]


def test_a_seqid_without_genes_does_not_count(tmp_path):
    fna, gff = _home(tmp_path, ["chr1"], ["chr1"])
    with open(gff, "a") as fh:
        fh.write("chrM\tsrc\tregion\t1\t10\t.\t+\t.\tID=chrM\n")
    assert ag.check_home_names(fna, gff) is None


# --------------------------------------------------------------------------- the writers

def _flank_inputs(tmp_path, flank_id):
    rng = random.Random(1)
    goi = "M" + "".join(rng.choice(AA) for _ in range(99))
    seq = list("".join(rng.choice("ACGT") for _ in range(9000)))
    gff = ["##gff-version 3"]
    for i, (gid, prot) in enumerate([(flank_id, None), ("gene-GOI", goi), ("gene-b", None)]):
        prot = prot or "M" + "".join(rng.choice(AA) for _ in range(99))
        s = 1500 * (i + 1)
        cds = "".join(CODON[a] for a in prot) + "TAA"
        seq[s - 1:s - 1 + len(cds)] = cds
        e = s + len(cds) - 1
        gff += [f"c\tt\tgene\t{s}\t{e}\t.\t+\t.\tID={gid};gene_biotype=protein_coding",
                f"c\tt\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{i};Parent={gid}",
                f"c\tt\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{i};Parent=rna-{i}"]
    (tmp_path / "g.fna").write_text(">c\n" + "".join(seq) + "\n")
    (tmp_path / "g.gff").write_text("\n".join(gff) + "\n")
    (tmp_path / "goi.faa").write_text(f">GOI_q\n{goi}\n")
    (tmp_path / "hits.tsv").write_text(f"q\tc\t100.0\t100\t0\t0\t1\t100\t3000\t3299\t1e-30\t200\n")
    (tmp_path / "locus_1.bed").write_text("c\t2999\t3302\tgene_loc\t1e-30\t+\t150\n")


def _flank(tmp_path):
    return subprocess.run([sys.executable, os.path.join(BIN, "extract_flanking_genes.py"),
                           "--bed", str(tmp_path / "locus_1.bed"), "--gff", str(tmp_path / "g.gff"),
                           "--genome", str(tmp_path / "g.fna"), "--n_flank", "2", "--min_size", "0",
                           "--exon_mode", "true", "--goi_faa", str(tmp_path / "goi.faa"),
                           "--goi_relative_evalue", "0.00001", "--home_hits", str(tmp_path / "hits.tsv"),
                           "--max_evalue", "0.01", "--out_bed", str(tmp_path / "out.bed"),
                           "--out_faa", str(tmp_path / "out.faa")], capture_output=True, text=True, cwd=tmp_path)


def test_flanking_genes_with_ordinary_ids_run(tmp_path):
    _flank_inputs(tmp_path, "gene-a")
    p = _flank(tmp_path)
    assert p.returncode == 0, p.stderr
    assert "gene-a" in (tmp_path / "out.bed").read_text()


@pytest.mark.parametrize("flank_id,shown", [("lcl|gene-a", "'lcl|gene-a'"), ("ucp1", "'ucp1'"),
                                            ("ORF 1", "'ORF 1'")])
def test_a_flanking_gene_id_that_loses_its_hits_stops(tmp_path, flank_id, shown):
    _flank_inputs(tmp_path, flank_id)
    p = _flank(tmp_path)
    assert p.returncode != 0
    assert shown in p.stderr and "g.gff" in p.stderr and "out.faa" in p.stderr


def test_gff_to_faa_stops_on_a_gene_id_with_whitespace(tmp_path):
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(300))
    prot = "".join(CODON[a] for a in "MKTAYIAKQR") + "TAA"
    seq = seq[:30] + prot + seq[30 + len(prot):]
    (tmp_path / "g.fna").write_text(f">c\n{seq}\n")
    (tmp_path / "g.gff").write_text(
        f"c\tt\tgene\t31\t63\t.\t+\t.\tID=ORF 1\nc\tt\tmRNA\t31\t63\t.\t+\t.\tID=r1;Parent=ORF 1\n"
        f"c\tt\tCDS\t31\t63\t.\t+\t0\tID=c1;Parent=r1\n")
    p = subprocess.run([sys.executable, os.path.join(BIN, "gff_to_faa.py"), "--gff", str(tmp_path / "g.gff"),
                        "--genome", str(tmp_path / "g.fna"), "--output", str(tmp_path / "p.faa")],
                       capture_output=True, text=True)
    assert p.returncode != 0 and "'ORF 1'" in p.stderr and "whitespace" in p.stderr


def test_borrow_annotations_checks_both_protein_files():
    """The third writer: the home proteome of a run without a home GFF, and the window proteins."""
    src = open(os.path.join(BIN, "borrow_annotations.py")).read()
    assert src.count("check_gene_ids(") == 2


# --------------------------------------------------------------------------- the query

def _normalize(tmp_path, qid):
    (tmp_path / "q.faa").write_text(f">{qid} a protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n")
    return subprocess.run([sys.executable, os.path.join(BIN, "normalize_query.py"), "--input", str(tmp_path / "q.faa"),
                           "--output", str(tmp_path / "n.faa")], capture_output=True, text=True)


@pytest.mark.parametrize("qid", ["my_gene|exon_2", "FUN_1|exon_2", "x|frag_1_30_120", "x|slide_1_50_200", "x|frag_a"])
def test_a_query_id_of_a_piece_stops(tmp_path, qid):
    p = _normalize(tmp_path, qid)
    assert p.returncode == 1 and qid in p.stderr and "Rename the query" in p.stderr
    assert not (tmp_path / "n.faa").exists()


@pytest.mark.parametrize("qid", ["my_gene", "XP_006563422.1", "sp|P01501|MEL_APIME", "exon_2", "my|gene"])
def test_other_query_ids_run(tmp_path, qid):
    assert _normalize(tmp_path, qid).returncode == 0


# --------------------------------------------------------------------------- the report

def test_the_report_warns_on_a_partial_overlap(tmp_path):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n")
    gi = tmp_path / "goi_info.json"
    gi.write_text(json.dumps({"home_gff_unmatched": {"seqids": 2, "of": 9, "examples": ["LG_2", "LG_3"]}}))
    out = tmp_path / "report.json"
    subprocess.run([sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
                    "--goi_info", str(gi), "--output", str(out), "--allow-empty"], check=True, capture_output=True)
    rep = json.loads(out.read_text())
    assert rep["home_gff_unmatched"]["seqids"] == 2
    assert any("2 of the home GFF's 9 sequences" in w and "LG_2" in w for w in rep["warnings"])


def test_a_gff_whose_genes_all_lie_on_removed_records_does_not_stop(tmp_path):
    """Review of the session: leave the removed records out before asking whether the GFF shares a seqid."""
    fna, gff = _home(tmp_path, ["NC_1.1"], ["NW_7.1"])
    assert ag.check_home_names(fna, gff, removed_names=["NW_7.1"]) is None
