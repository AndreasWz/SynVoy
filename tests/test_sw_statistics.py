"""Smith-Waterman hits carry real Karlin-Altschul statistics (2026-09-11).

Before: every SW hit was written with E=0.001 and its RAW score as the bit score, so the
~6 chance hits SW returns per window (best alignment per reading frame, homolog or not)
cleared every E-value filter and outranked real tblastn/MMseqs2 bit scores in the
spatial dedup. Now: BLAST gap costs 11/1 (parasail 12/1), bit score and E-value for the
window's six-frame search space, so the runner's existing E <= 1 filter works.
"""
from __future__ import annotations

import math
import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import smith_waterman_search as sws  # noqa: E402

parasail = pytest.importorskip("parasail")

CODON = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA',
         'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
         'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT'}
QUERY = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"


def test_blast_gap_costs_converted_for_parasail():
    assert sws.to_parasail_gaps(11, 1) == (12, 1)
    # one deleted residue between exact matches costs BLAST's 11 + 1*1 = 12
    q, t = "ACDEFGHIK", "ACDEGHIK"          # F deleted; matches score 47 in BLOSUM62
    p_open, p_ext = sws.to_parasail_gaps(11, 1)
    assert 47 - parasail.sw_striped_32(q, t, p_open, p_ext, parasail.blosum62).score == 12


def test_bit_score_matches_blast_for_11_1():
    # BLAST reports "43.1 bits (100)" for raw score 100 with BLOSUM62 11/1.
    bits, _ = sws.karlin_altschul(100, 70, 1_000_000, 11, 1)
    assert round(bits, 1) == 43.1


def test_evalue_scales_with_search_space():
    _, e1 = sws.karlin_altschul(50, 70, 100_000, 11, 1)
    _, e2 = sws.karlin_altschul(50, 70, 200_000, 11, 1)
    assert e2 == pytest.approx(2 * e1)
    assert e1 == pytest.approx(0.041 * 70 * 100_000 * math.exp(-0.267 * 50))


def test_unsupported_gap_costs_fail_loud():
    with pytest.raises(ValueError):
        sws.karlin_altschul(50, 70, 1000, 7, 3)


def _run(tmp_path, query, window, name):
    q = tmp_path / f"{name}.faa"
    t = tmp_path / f"{name}.fna"
    o = tmp_path / f"{name}.m8"
    q.write_text(f">GOI_q\n{query}\n")
    t.write_text(f">region_seq\n{window}\n")
    subprocess.run([sys.executable, str(ROOT / "bin" / "smith_waterman_search.py"),
                    "--query", str(q), "--target", str(t), "--output", str(o),
                    "--min_score", "20", "--min_identity", "10.0", "--method", "parasail"],
                   check=True, capture_output=True, text=True)
    rows = [line.rstrip("\n").split("\t") for line in o.read_text().splitlines() if line.strip()]
    return [(float(r[10]), float(r[11]), r) for r in rows]


def test_planted_homolog_significant_chance_hits_not(tmp_path):
    rng = random.Random(17)
    dna = ''.join(rng.choice("ACGT") for _ in range(60_000))
    gene = ''.join(CODON[a] for a in QUERY)
    window = dna[:30_000] + gene + dna[30_000:]

    hits = _run(tmp_path, QUERY, window, "planted")
    best_e, best_bits, best = min(hits, key=lambda h: h[0])
    assert best_e < 1e-10
    assert int(best[8]) == 30_001                 # tstart, 1-based, forward frame
    assert best_bits < float(best[3]) * 5         # a bit score, not the raw score

    shuffled = list(QUERY)
    rng.shuffle(shuffled)
    null = _run(tmp_path, ''.join(shuffled), dna, "null")
    # SW finds a best alignment in every frame, homolog or not; none of them survives the
    # runner's E <= 1 filter, and since 2026-10-29 (fix-phase item 18) the child no longer
    # aligns them with a traceback or writes them (min_reportable_score).
    assert all(e > 1.0 for e, _b, _r in null)
    assert null == []


# ----------------------------------------------------------- long frames in pieces (2026-10-09)
# The traceback needs about 4 bytes per cell (query length x frame length). One search
# window of 13.7 Mb took 10 GB for a 550-aa query and the kernel ended the run.

def test_a_frame_within_the_budget_is_one_piece():
    assert sws.protein_chunks("A" * 1000, 70) == [(0, "A" * 1000)]
    assert sws.protein_chunks("", 70) == [(0, "")]


def test_pieces_cover_the_frame_and_overlap_by_twice_the_query():
    protein = ''.join(random.Random(3).choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(10_000))
    qlen = 50
    pieces = sws.protein_chunks(protein, qlen, max_cells=50 * 1_000)      # pieces of 1,000 residues
    assert len(pieces) > 5
    assert all(len(p) * qlen <= 50 * 1_000 for _s, p in pieces)
    assert pieces[0][0] == 0 and pieces[-1][0] + len(pieces[-1][1]) == len(protein)
    for (s1, p1), (s2, _p2) in zip(pieces, pieces[1:]):
        assert s1 + len(p1) - s2 == 2 * qlen                              # the overlap
    # every stretch of up to 2 x query length lies entirely inside one piece
    for start in range(0, len(protein) - 2 * qlen, 37):
        assert any(s <= start and start + 2 * qlen <= s + len(p) for s, p in pieces)


def _run_main(tmp_path, monkeypatch, query, window, name, max_cells=None, extra=()):
    q, t, o = tmp_path / f"{name}.faa", tmp_path / f"{name}.fna", tmp_path / f"{name}.m8"
    q.write_text(f">GOI_q\n{query}\n")
    t.write_text(f">region_seq\n{window}\n")
    if max_cells is not None:
        monkeypatch.setattr(sws, "SW_MAX_CELLS", max_cells)
    monkeypatch.setattr(sys, "argv", ["smith_waterman_search.py", "--query", str(q), "--target", str(t),
                                      "--output", str(o), "--min_score", "20", "--min_identity", "10.0",
                                      "--method", "parasail", *extra])
    sws.main()
    return [line.rstrip("\n").split("\t") for line in o.read_text().splitlines() if line.strip()]


@pytest.mark.parametrize("position", [30_000, 4_190, 8_399, 12_600])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_pieces_find_the_same_gene_at_the_same_place(tmp_path, monkeypatch, position, strand):
    """A gene anywhere in the window, also across the border of two pieces, on either
    strand, comes back with the coordinates, identity and E-value of the undivided run."""
    rng = random.Random(23)
    dna = ''.join(rng.choice("ACGT") for _ in range(60_000))
    gene = ''.join(CODON[a] for a in QUERY)
    if strand == "-":
        gene = gene[::-1].translate(str.maketrans("ACGT", "TGCA"))
    window = dna[:position] + gene + dna[position:]

    whole = _run_main(tmp_path, monkeypatch, QUERY, window, "whole")
    # 70 aa x 1,400 residues per piece = pieces of 4.2 kb of DNA, overlapping by 420 bp
    pieces = _run_main(tmp_path, monkeypatch, QUERY, window, "pieces", max_cells=len(QUERY) * 1_400)

    best_whole = min(whole, key=lambda r: float(r[10]))
    best_pieces = min(pieces, key=lambda r: float(r[10]))
    assert float(best_whole[10]) < 1e-10
    assert best_pieces == best_whole
    lo, hi = sorted((int(best_whole[8]), int(best_whole[9])))
    assert (lo, hi) == (position + 1, position + len(gene))
    # Only pieces that can give a kept hit are aligned with a traceback (fix-phase item 18), so
    # the pieces add no line with E <= 1 and the gene is written once.
    assert [r for r in pieces if float(r[10]) <= 1] == [r for r in whole if float(r[10]) <= 1]
    assert len({(r[8], r[9]) for r in pieces}) == len(pieces)   # nothing twice


def test_a_query_too_long_for_the_traceback_is_skipped_with_a_message(tmp_path, monkeypatch, capsys):
    """Pieces cannot be shorter than four query lengths, so memory grows with the square of
    the query length; beyond the limit the query is left to the other search programs."""
    monkeypatch.setattr(sws, "SW_MAX_QUERY_LEN", 50)
    rng = random.Random(5)
    dna = ''.join(rng.choice("ACGT") for _ in range(3_000))
    gene = ''.join(CODON[a] for a in QUERY)
    rows = _run_main(tmp_path, monkeypatch, QUERY, dna[:1_500] + gene + dna[1_500:], "long")   # 70 aa > 50
    assert rows == []
    assert "Smith-Waterman skipped for GOI_q: 70 aa" in capsys.readouterr().err
    monkeypatch.setattr(sws, "SW_MAX_QUERY_LEN", 10_000)
    assert _run_main(tmp_path, monkeypatch, QUERY, dna[:1_500] + gene + dna[1_500:], "ok")


# ---------------------------------------------------------------- several alignments per piece
def _gene(protein):
    return ''.join(CODON[a] for a in protein)


def test_two_copies_in_one_frame_are_both_reported_with_three_alignments(tmp_path, monkeypatch):
    """The best alignment of a piece hides a second copy in the same frame: at the default
    one of the two is reported, at 3 both, each at its own coordinates."""
    rng = random.Random(41)
    dna = ''.join(rng.choice("ACGT") for _ in range(9_000))
    gene = _gene(QUERY)
    first, second = 2_400, 2_400 + len(gene) + 900          # same frame: 900 is a multiple of 3
    window = dna[:first] + gene + dna[first:first + 900] + gene + dna[first + 900:]
    copies = {(first + 1, first + len(gene)), (second + 1, second + len(gene))}

    def strong(rows):
        return {(int(r[8]), int(r[9])) for r in rows if float(r[10]) < 1e-10}

    assert len(strong(_run_main(tmp_path, monkeypatch, QUERY, window, "one"))) == 1
    three = _run_main(tmp_path, monkeypatch, QUERY, window, "three",
                      extra=("--max_alignments_per_piece", "3"))
    assert strong(three) == copies


def test_three_alignments_add_nothing_without_a_real_hit(tmp_path, monkeypatch):
    """Without an alignment at E <= 1 a piece gives nothing at 3, as at 1 (since 2026-10-29 its
    best alignment, E > 1, is not aligned with a traceback or written at all)."""
    window = ''.join(random.Random(43).choice("ACGT") for _ in range(30_000))
    one = _run_main(tmp_path, monkeypatch, QUERY, window, "one")
    three = _run_main(tmp_path, monkeypatch, QUERY, window, "three",
                      extra=("--max_alignments_per_piece", "3"))
    assert one == [] and three == []


# Pieces of 1,400 residues overlapping by 140: at these positions half of the copy lies in
# the last 35 residues of a piece (frame residues 1,365-1,434 or 2,625-2,694).
@pytest.mark.parametrize("strand,position", [("+", 4_095), ("-", 55_905), ("+", 7_875)])
def test_a_copy_across_a_piece_border_is_reported_once_and_whole(tmp_path, monkeypatch, position, strand):
    """Overlapping pieces see the same copy, one of them cut at its edge. At 3 the copy is
    reported once, with the coordinates of the undivided run; without the edge rule it came
    back as two halves (a part of it scores E <= 1; chance alignments of the other frames
    over the same DNA do not)."""
    rng = random.Random(23)
    dna = ''.join(rng.choice("ACGT") for _ in range(60_000))
    gene = _gene(QUERY)
    if strand == "-":
        gene = gene[::-1].translate(str.maketrans("ACGT", "TGCA"))
    window = dna[:position] + gene + dna[position:]
    rows = _run_main(tmp_path, monkeypatch, QUERY, window, "pieces", max_cells=len(QUERY) * 1_400,
                     extra=("--max_alignments_per_piece", "3"))
    at_gene = [r for r in rows if float(r[10]) <= 1
               and min(int(r[8]), int(r[9])) <= position + len(gene) and max(int(r[8]), int(r[9])) > position]
    assert len(at_gene) == 1
    assert sorted((int(at_gene[0][8]), int(at_gene[0][9]))) == [position + 1, position + len(gene)]
    assert float(at_gene[0][10]) < 1e-10


def test_the_option_must_be_positive(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        _run_main(tmp_path, monkeypatch, QUERY, "ACGT" * 100, "zero", extra=("--max_alignments_per_piece", "0"))
