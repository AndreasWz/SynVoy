"""Every length floor of the search counts aligned query residues (bug list item 6, 2026-11-16).

MMseqs2's translated search, run without a backtrace as the search runs it, writes the target
span in nucleotides in its alignment-length column; tblastn and Smith-Waterman write alignment
columns. The floors are set in residues, so an MMseqs2 hit was held to about a third of each
(the rescue's 50 acted like 17). They now read qend - qstart + 1
(iterative_search_runner.aligned_query_residues); a hit keeps the engine's own alnlen.
"""
import os
import random
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import iterative_search_runner as isr  # noqa: E402
from sequence_utils import MMSEQS_FORMAT_OUTPUT, mmseqs_table_first_words  # noqa: E402


def _m8(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


# query target pident alnlen mismatch gapopen qstart qend tstart tend evalue bits
MMSEQS_ROW = ["GOI_q", "c1", 60.0, 60, 0, 0, 1, 20, 1001, 1060, 1e-5, 40.0]   # alnlen: target nt
TBLASTN_ROW = ["GOI_q", "c1", 60.0, 21, 0, 1, 1, 20, 2001, 2063, 1e-5, 40.0]  # alnlen: columns


def test_the_floor_reads_query_residues_for_every_engine(tmp_path):
    p = tmp_path / "h.m8"
    _m8(p, [MMSEQS_ROW, TBLASTN_ROW])
    drops = {}
    # 20 residues each: below a floor of 25 both are dropped (the MMseqs2 row's 60 was kept before)
    assert isr.parse_hits(str(p), 10.0, 25, 0.01, drops=drops) == []
    assert drops == {"len_floor": 2}
    hits = isr.parse_hits(str(p), 10.0, 20, 0.01)
    assert len(hits) == 2
    assert [isr.aligned_query_residues(h) for h in hits] == [20, 20]
    assert sorted(h["alnlen"] for h in hits) == [21, 60]           # the engine's own column kept


def test_the_short_query_floor_and_the_long_query_count_read_residues(tmp_path):
    p = tmp_path / "h.m8"
    # 12 residues, 36 nt: on a 200-aa query below the floor of 15 and counted as a short drop
    _m8(p, [["GOI_q", "c1", 60.0, 36, 0, 0, 5, 16, 1001, 1036, 1e-5, 30.0]])
    drops = {}
    assert isr.parse_hits(str(p), 10.0, 15, 0.01, query_lengths={"GOI_q": 200}, drops=drops) == []
    assert drops == {"len_floor": 1, "len_floor_lt15_long_query": 1}
    # a 16-aa query: its floor is max(12, round(0.6 x 16)) = 12 residues, which the hit reaches
    assert len(isr.parse_hits(str(p), 10.0, 15, 0.01, query_lengths={"GOI_q": 16})) == 1


def test_the_rescue_length_test_reads_query_residues():
    strong = {"query": "GOI_q", "chrom": "c1", "start": 5_000_000, "end": 5_000_150, "strand": "+",
              "qstart": 1, "qend": 50, "pident": 80.0, "evalue": 1e-30, "alnlen": 150, "bits": 120.0}
    drops = {}
    assert isr.select_dispersed_goi_seeds([strong], [], min_alnlen=60, drops=drops) == []
    assert drops == {"dispersed_alnlen": 1}
    drops, rejected = {}, []
    isr.select_dispersed_goi_seeds([strong], [], min_alnlen=50, drops=drops, rejected_out=rejected)
    assert "dispersed_alnlen" not in drops and len(rejected) == 1   # past the length test
    assert rejected[0]["alnlen"] == 150                              # the record keeps alnlen


def test_pooled_identity_weights_by_query_residues():
    mmseqs = {"pident": 80.0, "alnlen": 300, "qstart": 1, "qend": 100}       # 300 nt
    tblastn = {"pident": 40.0, "alnlen": 100, "qstart": 101, "qend": 200}    # 100 columns
    assert isr._pooled_identity([mmseqs, tblastn]) == 60.0                    # was 70.0


def test_the_helper_falls_back_to_alnlen_without_query_coordinates():
    assert isr.aligned_query_residues({"alnlen": 42}) == 42
    assert isr.aligned_query_residues({"alnlen": 42, "qstart": "", "qend": ""}) == 42
    assert isr.aligned_query_residues({}) == 0
    assert isr.aligned_query_residues({"qstart": 30, "qend": 11, "alnlen": 90}) == 20


def test_the_hit_chain_tie_break_reads_query_residues():
    # two single hits 60 kb apart, equal bits: the one with more query residues is returned
    a = {"gstart": 0, "gend": 300, "qstart": 1, "qend": 40, "bits": 50.0, "alnlen": 300}
    b = {"gstart": 60_000, "gend": 60_150, "qstart": 41, "qend": 90, "bits": 50.0, "alnlen": 50}
    assert isr._longest_monotonic_query_chain([a, b], "+", max_gap=20_000) == [b]


def test_no_length_floor_of_the_runner_reads_alnlen():
    """By the source: a comparison or sum on an alnlen read is gone from the runner (the copies of
    hit dicts and the helper's own fallback stay)."""
    with open(os.path.join(ROOT, "bin", "iterative_search_runner.py")) as fh:
        src = fh.read()
    head, rest = src.split("def aligned_query_residues(", 1)
    src = head + rest.split("\ndef ", 1)[1]                  # the helper's own fallback stays
    for bad in ("get('alnlen', 0) <", "get('alnlen', 0)) <", "int(h.get('alnlen'",
                "int(ordered_hits[0].get('alnlen'", 'int(h.get("alnlen"', "alnlen < effective_min_length",
                "alnlen >= effective_min_length"):
        assert bad not in src, bad


# ------------------------------------------------------------------ the engines themselves
AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
         "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
         "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


def _synthetic(tmp_path):
    """A 200-aa protein and a contig that encodes it with 10 extra codons in the middle."""
    rng = random.Random(7)
    q = "M" + "".join(rng.choice(AA) for _ in range(199))
    target = q[:100] + "GSGSGSGSGS" + q[100:]
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    contig = flank(2000) + "".join(CODON[a] for a in target) + flank(2000)
    (tmp_path / "q.faa").write_text(f">GOI_q\n{q}\n")
    (tmp_path / "t.fna").write_text(f">ctg1\n{contig}\n")
    return q


def _run(cmd, cwd):
    done = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-500:]


@pytest.mark.skipif(not shutil.which("mmseqs"), reason="mmseqs not on PATH")
def test_mmseqs_as_the_search_runs_it_writes_the_target_span(tmp_path):
    _synthetic(tmp_path)
    for cmd in (["mmseqs", "createdb", "q.faa", "qdb", "-v", "1"],
                ["mmseqs", "createdb", "t.fna", "tdb", "-v", "1"],
                ["mmseqs", "search", "qdb", "tdb", "res", "tmp", "--search-type", "2", "-s", "7.5",
                 "-e", "10", "--threads", "1", "-v", "1"],
                ["mmseqs", "convertalis", "qdb", "tdb", "res", "mm.raw", "--format-output",
                 MMSEQS_FORMAT_OUTPUT, "-v", "1"]):
        _run(cmd, tmp_path)
    mmseqs_table_first_words(str(tmp_path / "mm.raw"), str(tmp_path / "mm.m8"))
    rows = [ln.split("\t") for ln in (tmp_path / "mm.m8").read_text().splitlines() if ln]
    assert rows
    best = max(rows, key=lambda r: float(r[11]))
    alnlen, qs, qe, ts, te = (int(best[i]) for i in (3, 6, 7, 8, 9))
    assert alnlen == abs(te - ts) + 1                    # nucleotides of the target
    assert alnlen > 2 * (qe - qs + 1)
    hits = isr.parse_hits(str(tmp_path / "mm.m8"), 0.0, 1, 10.0)
    assert max(isr.aligned_query_residues(h) for h in hits) == qe - qs + 1
    # a floor one residue above the hit's query span drops it, whatever alnlen says
    assert not [h for h in isr.parse_hits(str(tmp_path / "mm.m8"), 0.0, qe - qs + 2, 10.0)
                if (h["qstart"], h["qend"]) == (qs, qe)]


@pytest.mark.skipif(not (shutil.which("tblastn") and shutil.which("makeblastdb")),
                    reason="BLAST+ not on PATH")
def test_tblastn_writes_alignment_columns_and_the_floor_reads_residues(tmp_path):
    _synthetic(tmp_path)
    _run(["makeblastdb", "-in", "t.fna", "-dbtype", "nucl", "-out", "tdb"], tmp_path)
    _run(["tblastn", "-query", "q.faa", "-db", "tdb", "-seg", "no", "-evalue", "10",
          "-outfmt", isr.TBLASTN_TABULAR, "-out", "tb.m8"], tmp_path)
    rows = [ln.split("\t") for ln in (tmp_path / "tb.m8").read_text().splitlines() if ln]
    assert rows
    best = max(rows, key=lambda r: float(r[11]))
    alnlen, qs, qe = int(best[3]), int(best[6]), int(best[7])
    assert alnlen >= qe - qs + 1                          # query residues plus gap columns
    hits = isr.parse_hits(str(tmp_path / "tb.m8"), 0.0, 1, 10.0)
    assert max(isr.aligned_query_residues(h) for h in hits) == qe - qs + 1
