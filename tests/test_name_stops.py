"""Fix-phase d4 (2 Nov; the loud halves of items 29, 30, 31, 32 and 36): names that lost results in silence stop the run.

- 29: a flanking gene id, or a home sequence name, that MMseqs2 reports as another name (lcl|X as X, gnl|db|X as X, ucon_5 as
  on_5 ...): the flanking gene's hits were filed under a name nothing else used and its model was lost; a home sequence gave two
  loci at one place. d4 stopped them; since 9 Nov (fix item 29) every MMseqs2 table is read by its FASTA headers and these
  names work (tests/test_mmseqs_first_words.py). A '|' in a flanking gene id still stops: SynVoy reads a gene id up to it.
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

@pytest.mark.parametrize("name", ["gene-LOC100004", "NC_000012.12", "NT_033777.3", "FBgn0000508", "AT1G01010",
                                  "orf19.1562", "gene-orf19.1562", "NODE_1_length_5000_cov_1.0", "Ca22chr1A",
                                  "SUPER_10_unloc_1", "scaffold_12"])
def test_ordinary_names_pass(name):
    assert su.id_problem(name, pipe=True) is None


def test_each_class_is_named():
    assert "whitespace" in su.id_problem("ORF 1")
    assert "'|'" in su.id_problem("ctg|arrow", pipe=True)
    assert su.id_problem("ctg|arrow") is None               # '|' only where SynVoy cuts at it
    for name in ("lcl|X", "ucon_5", "consensus_24", "gi|1|ref|X|"):   # MMseqs2 rewrites: no problem since item 29
        assert su.id_problem(name) is None


def test_check_gene_ids_names_the_id_and_both_files():
    with pytest.raises(su.StepFailed) as stop:
        su.check_gene_ids(["gene-a", "ORF 1"], "home.gff", "home_proteome.faa")
    msg = str(stop.value)
    assert "'ORF 1'" in msg and "home.gff" in msg and "home_proteome.faa" in msg and "Rename" in msg
    su.check_gene_ids(["lcl|X"], "home.gff", "out.faa")       # '|' only when asked
    with pytest.raises(su.StepFailed, match="first '\\|'"):
        su.check_gene_ids(["lcl|X"], "home.gff", "out.faa", pipe=True)


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


def test_a_home_name_mmseqs2_rewrites_runs_with_a_gff(tmp_path):
    """Fix item 29: LOCATE_GENE reads its table by the FASTA headers, so lcl|chr2 is one locus, not two.
    Without a GFF the '|' still stops (item 30: the predicted genes' ids hold the name)."""
    fna, gff = _home(tmp_path, ["chr1", "lcl|chr2"], ["chr1", "lcl|chr2"])
    assert ag.check_home_names(fna, gff) is None
    with pytest.raises(su.StepFailed) as stop:
        ag.check_home_names(fna, None)
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


def test_a_flanking_gene_id_mmseqs2_rewrites_runs(tmp_path):
    """Fix item 29: ucp1 (MMseqs2 reports it as p1) stopped since 2026-11-02; its hits are now read by the
    FASTA header, so it runs."""
    _flank_inputs(tmp_path, "ucp1")
    p = _flank(tmp_path)
    assert p.returncode == 0, p.stderr
    assert "ucp1" in (tmp_path / "out.bed").read_text()


@pytest.mark.parametrize("flank_id,shown", [("lcl|gene-a", "'lcl|gene-a'"), ("ctg|1", "'ctg|1'"),
                                            ("ORF 1", "'ORF 1'")])
def test_a_flanking_gene_id_that_loses_its_hits_stops(tmp_path, flank_id, shown):
    _flank_inputs(tmp_path, flank_id)
    p = _flank(tmp_path)
    # A stop ends with exit 3 and its message, no traceback (fix item d4: an uncaught StepFailed, exit 1).
    assert p.returncode == 3 and "Traceback" not in p.stderr, p.stderr
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
    assert p.returncode == 3 and "'ORF 1'" in p.stderr and "whitespace" in p.stderr
    assert "Traceback" not in p.stderr


def test_annotate_goi_exons_stops_with_exit_3_and_no_traceback(tmp_path):
    """Fix item d4: ANNOTATE_GOI's name stops ended as an uncaught StepFailed (exit 1 with a traceback;
    the loudness session's two constructed runs)."""
    (tmp_path / "q.faa").write_text(">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n")
    (tmp_path / "h.fna").write_text(">c1|arrow\n" + "ACGT" * 50 + "\n")
    p = subprocess.run([sys.executable, os.path.join(BIN, "annotate_goi_exons.py"), "--query", str(tmp_path / "q.faa"),
                        "--genome", str(tmp_path / "h.fna"), "--gff", "NO_GFF",
                        "--output_exons", str(tmp_path / "e.faa"), "--output_bed", str(tmp_path / "e.bed"),
                        "--output_info", str(tmp_path / "i.json")], capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 3, p.stderr
    assert "Traceback" not in p.stderr and "'c1|arrow'" in p.stderr and p.stderr.startswith("ERROR: ")


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


@pytest.mark.parametrize("qid", ["my_gene", "XP_100003.1", "sp|P01501|MEL_APIME", "exon_2"])
def test_other_query_ids_run(tmp_path, qid):
    assert _normalize(tmp_path, qid).returncode == 0


def test_a_pipe_left_in_a_query_id_stops(tmp_path):
    """my|gene ran until freeze fix e2, and the GOI was named 'my' (tests/test_name_points.py)."""
    p = _normalize(tmp_path, "my|gene")
    assert p.returncode == 1 and "'|'" in p.stderr and "my_gene" in p.stderr


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
