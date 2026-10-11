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
    assert null, "SW always returns a best hit per frame"
    assert all(e > 1.0 for e, _b, _r in null)     # none survive the runner's E <= 1 filter


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


def _run_main(tmp_path, monkeypatch, query, window, name, max_cells=None):
    q, t, o = tmp_path / f"{name}.faa", tmp_path / f"{name}.fna", tmp_path / f"{name}.m8"
    q.write_text(f">GOI_q\n{query}\n")
    t.write_text(f">region_seq\n{window}\n")
    if max_cells is not None:
        monkeypatch.setattr(sws, "SW_MAX_CELLS", max_cells)
    monkeypatch.setattr(sys, "argv", ["smith_waterman_search.py", "--query", str(q), "--target", str(t),
                                      "--output", str(o), "--min_score", "20", "--min_identity", "10.0",
                                      "--method", "parasail"])
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
    assert len(pieces) > len(whole)                      # the best hit of every piece is written
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
