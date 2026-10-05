"""normalize_query.py is the first script that touches user input.

It used to pass a protein query through verbatim, so anything that was formatting
rather than sequence reached the search: lowercase letters, the gaps of an aligned
FASTA, the digits and spaces of numbered GenBank-style lines, a trailing ``*``.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import normalize_query as nq  # noqa: E402

MELITTIN = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"
MELITTIN_CDS = (
    "ATGAAATTCTTAGTCAACGTTGCCCTTGTTTTTATGGTCGTATACATTTCTTACATCTATGCGGCCCCTGAACCGGAACCG"
    "GCACCAGAGCCAGAGGCGGAGGCAGACGCGGAGGCAGATCCGGAAGCGGGAATTGGAGCAGTTCTGAAGGTATTAACCACA"
    "GGATTGCCCGCCCTCATAAGTTGGATTAAACGTAAGAGGCAACAGGGTTAA"
)


def run(tmp_path: Path, text: str, *extra: str):
    src = tmp_path / "in.fa"
    src.write_bytes(text.encode("utf-8"))
    out = tmp_path / "out.faa"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "normalize_query.py"),
         "--input", str(src), "--output", str(out), *extra],
        capture_output=True, text=True)
    seq = None
    if out.exists():
        lines = out.read_text().splitlines()
        seq = "".join(l for l in lines if not l.startswith(">"))
    return proc, seq


# ───────────────────────── cleaning ─────────────────────────

def test_plain_protein_is_unchanged(tmp_path):
    proc, seq = run(tmp_path, f">q1 melittin\n{MELITTIN}\n")
    assert proc.returncode == 0 and seq == MELITTIN
    assert "Removed" not in proc.stderr


def test_lowercase_is_folded_to_upper(tmp_path):
    proc, seq = run(tmp_path, f">q1\n{MELITTIN.lower()}\n")
    assert proc.returncode == 0 and seq == MELITTIN


def test_numbered_and_spaced_lines_are_cleaned(tmp_path):
    body = f"        1 {MELITTIN[:10].lower()} {MELITTIN[10:20].lower()}\n       21 {MELITTIN[20:]}\n"
    proc, seq = run(tmp_path, ">q1\n" + body)
    assert proc.returncode == 0 and seq == MELITTIN
    assert "formatting character" in proc.stderr


def test_alignment_gaps_are_removed(tmp_path):
    gapped = MELITTIN[:30] + "---" + MELITTIN[30:50] + "..." + MELITTIN[50:]
    proc, seq = run(tmp_path, f">q1\n{gapped}\n")
    assert proc.returncode == 0 and seq == MELITTIN


def test_trailing_stop_symbol_is_dropped(tmp_path):
    proc, seq = run(tmp_path, f">q1\n{MELITTIN}*\n")
    assert proc.returncode == 0 and seq == MELITTIN


def test_internal_stop_is_removed_with_a_warning(tmp_path):
    proc, seq = run(tmp_path, f">q1\n{MELITTIN[:40]}*{MELITTIN[40:]}\n")
    assert proc.returncode == 0 and seq == MELITTIN
    assert "internal stop" in proc.stderr


def test_windows_line_endings_and_wrapping(tmp_path):
    proc, seq = run(tmp_path, f">q1\r\n{MELITTIN[:35]}\r\n{MELITTIN[35:]}\r\n")
    assert proc.returncode == 0 and seq == MELITTIN


def test_ambiguity_codes_and_non_ascii_header_are_accepted(tmp_path):
    odd = MELITTIN[:30] + "XBZUO" + MELITTIN[30:]
    proc, seq = run(tmp_path, f">q1 café α-toxin\n{odd}\n")
    assert proc.returncode == 0 and seq == odd


# ───────────────────────── rejections ─────────────────────────

@pytest.mark.parametrize("text", ["", "   \n\n", MELITTIN + "\n"],
                         ids=["empty", "whitespace", "no_header"])
def test_no_record_is_a_clean_error(tmp_path, text):
    proc, seq = run(tmp_path, text)
    assert proc.returncode == 1 and seq is None
    assert "Traceback" not in proc.stderr and "No sequences found" in proc.stderr


def test_header_without_sequence_is_a_clean_error(tmp_path):
    proc, seq = run(tmp_path, ">q1\n\n")
    assert proc.returncode == 1 and seq is None
    assert "no sequence" in proc.stderr


def test_non_sequence_characters_are_rejected(tmp_path):
    proc, seq = run(tmp_path, f">q1\n{MELITTIN[:30]}<br/>{MELITTIN[30:]}\n")
    assert proc.returncode == 1 and seq is None
    assert "not amino-acid" in proc.stderr and "'<'" in proc.stderr


def test_short_query_is_rejected_unless_overridden(tmp_path):
    proc, seq = run(tmp_path, ">q1\nMKFLVNVALV\n")
    assert proc.returncode == 2 and seq is None
    assert "below the minimum" in proc.stderr
    proc, seq = run(tmp_path, ">q1\nMKFLVNVALV\n", "--min_length", "0")
    assert proc.returncode == 0 and seq == "MKFLVNVALV"


def test_length_is_checked_after_cleaning(tmp_path):
    # 25 residues padded with gaps to 45 characters must still be "too short".
    proc, _ = run(tmp_path, ">q1\n" + MELITTIN[:25] + "-" * 20 + "\n")
    assert proc.returncode == 2


# ───────────────────────── nucleotide input ─────────────────────────

def test_dna_is_translated(tmp_path):
    proc, seq = run(tmp_path, f">q1\n{MELITTIN_CDS}\n")
    assert proc.returncode == 0 and seq == MELITTIN
    assert "nucleotide" in proc.stderr


def test_rna_and_lowercase_dna_are_translated(tmp_path):
    proc, seq = run(tmp_path, f">q1\n{MELITTIN_CDS.replace('T', 'U')}\n")
    assert proc.returncode == 0 and seq == MELITTIN
    proc, seq = run(tmp_path, f">q1\n{MELITTIN_CDS.lower()}\n")
    assert proc.returncode == 0 and seq == MELITTIN


def test_gly_ser_rich_peptide_is_not_mistaken_for_dna():
    # Every letter is also an IUPAC nucleotide code (G, S, A, K, D, ...).
    peptide = "GSGSGSGSGSKDGSGSGSGSGSGSKDGSGSGSGSGS"
    assert all(ch in nq.NUC_ALPHABET for ch in peptide)
    assert not nq.is_nucleotide(peptide)
    assert nq.is_nucleotide("ACGTACGTNNACGTACGTRYACGT")
    assert not nq.is_nucleotide("")
    assert not nq.is_nucleotide(MELITTIN)


# ───────────────────────── multi-record ─────────────────────────

def test_multi_record_uses_the_first_and_says_so(tmp_path):
    proc, seq = run(tmp_path, f">first\n{MELITTIN}\n>second\n{MELITTIN[::-1]}\n")
    assert proc.returncode == 0 and seq == MELITTIN
    assert "2 sequences found" in proc.stderr and "'first'" in proc.stderr


def test_clean_sequence_reports_what_it_removed():
    assert nq.clean_sequence("mk f-l.1v\t") == ("MKFLV", 5)
    assert nq.clean_sequence("MKFLV") == ("MKFLV", 0)
