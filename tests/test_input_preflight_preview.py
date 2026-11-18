"""What the start of a run checks before any task (bug list of 16 Nov 2026).

Runs ``nextflow run main.nf -preview`` (builds the DAG, runs no task), about 8 s per run; skipped
without Nextflow on PATH (CI has none). The checks:

* item 31 / F4a: a query FASTA with fewer sequence characters than --min_query_length stops the
  run at its start (NORMALIZE_QUERY would stop it for certain: it reads the first record and a DNA
  query translates to at most a third of its length). It warned with a fixed "30 aa minimum".
* item 31 / F4b: the warning for --max_genomes 1 or 2 says what fewer targets change; it described
  a classifier that no longer exists.
* F3: easy mode with --home_genome, --home_gff or --target_genomes stops with the hint to add
  --mode pro; it downloaded other genomes and ignored them without a word.
* item 39b: a listed target genome or annotation that matches nothing stops the run; it was left
  out without a word.
* item 22: --max_flanking_distance is checked as a number (from the command line it is text).
* items 21 and 27: the removed settings stop or warn (main.nf removedParams()).
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")

pytestmark = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")

PROTEIN_60 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV"


def _run(tmp_path: Path, *args: str, query: str = PROTEIN_60) -> tuple[int, str]:
    (tmp_path / "q.faa").write_text(f">q\n{query}\n")
    (tmp_path / "h.fna").write_text(">c1\nACGTACGTACGTACGTACGTACGTACGT\n")
    (tmp_path / "t1.fna").write_text(">c1\nACGTACGTACGTACGTACGTACGTACGT\n")
    (tmp_path / "t1.gff").write_text("##gff-version 3\nc1\t.\tgene\t1\t20\t.\t+\t.\tID=g1\n")
    noconda = tmp_path / "noconda.config"
    noconda.write_text("conda.enabled = false\nprocess.conda = null\nprocess.executor = 'local'\n")
    proc = subprocess.run(
        [NEXTFLOW, "run", str(ROOT / "main.nf"), "-preview", "-ansi-log", "false",
         "-c", str(noconda), "--outdir", "out", *args],
        cwd=tmp_path, capture_output=True, text=True, timeout=300)
    return proc.returncode, re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout + proc.stderr)


def _pro(tmp_path: Path, *args: str, query: str = PROTEIN_60) -> tuple[int, str]:
    return _run(tmp_path, "--mode", "pro", "--query", "q.faa", "--home_genome", "h.fna",
                "--home_species", "Testus_testus", *args, query=query)


# ----------------------------------------------------------------------- the query length
def test_a_query_shorter_than_the_minimum_stops_at_the_start(tmp_path):
    code, text = _pro(tmp_path, query="MKTAYIAKQRQISFV")
    assert code != 0
    assert ("Query FASTA 'q.faa' holds 15 sequence characters, fewer than --min_query_length "
            "(20 aa), so the run would stop at NORMALIZE_QUERY.") in text
    assert "30 aa minimum" not in text


def test_a_lower_minimum_lets_the_same_query_start(tmp_path):
    code, text = _pro(tmp_path, "--min_query_length", "10", query="MKTAYIAKQRQISFV")
    assert code == 0, text[-2000:]
    assert "fewer than --min_query_length" not in text


def test_a_query_above_the_minimum_is_not_flagged(tmp_path):
    code, text = _pro(tmp_path, query=PROTEIN_60[:25])
    assert code == 0, text[-2000:]
    assert "fewer than --min_query_length" not in text and "sequence characters" not in text


# ----------------------------------------------------------------------- easy mode
def test_few_genomes_warn_with_what_they_change(tmp_path):
    code, text = _run(tmp_path, "--mode", "easy", "--query_id", "P01501", "--max_genomes", "2")
    assert code == 0, text[-2000:]
    assert ("--max_genomes 2: each call is labelled from the evidence in its own genome, so the labels "
            "mean the same with one target as with many") in text
    assert "consensus across species" not in text


def test_easy_mode_with_pro_mode_inputs_stops(tmp_path):
    code, text = _run(tmp_path, "--mode", "easy", "--query_id", "P04637", "--home_genome", "h.fna",
                      "--target_genomes", "t1.fna")
    assert code != 0
    assert ("--home_genome, --target_genomes are pro-mode inputs: add --mode pro (easy mode, the default, "
            "downloads the home and the target genomes itself and would ignore them).") in text
    code, text = _run(tmp_path, "--query_id", "P04637", "--home_gff", "h.gff")
    assert code != 0 and "--home_gff is a pro-mode input: add --mode pro" in text


# ----------------------------------------------------------------------- the target lists
def test_a_listed_genome_that_does_not_exist_stops(tmp_path):
    code, text = _pro(tmp_path, "--target_genomes", "t1.fna,missing.fna")
    assert code != 0
    assert "--target_genomes names 1 file(s), folder(s) or pattern(s) that match nothing: missing.fna." in text


def test_a_pattern_that_matches_nothing_stops(tmp_path):
    code, text = _pro(tmp_path, "--target_genomes", "t1.fna,nomatch/*.fna")
    assert code != 0 and "that match nothing: nomatch/*.fna." in text
    code, text = _pro(tmp_path, "--target_genomes", "missing.fna")             # nothing at all: as before
    assert code != 0 and "No target genomes found for --target_genomes 'missing.fna'" in text


def test_existing_targets_and_empty_list_elements_pass(tmp_path):
    code, text = _pro(tmp_path, "--target_genomes", "t1.fna,,")
    assert code == 0, text[-2000:]
    assert "match nothing" not in text


def test_a_listed_annotation_that_does_not_exist_stops(tmp_path):
    code, text = _pro(tmp_path, "--target_genomes", "t1.fna", "--target_gffs", "t1.gff,missing.gff")
    assert code != 0
    assert "--target_gffs names 1 file(s), folder(s) or pattern(s) that match nothing: missing.gff." in text
    code, text = _pro(tmp_path, "--target_genomes", "t1.fna", "--target_gffs", "missing.gff")
    assert code != 0 and "No target annotations found for --target_gffs 'missing.gff'" in text


# ----------------------------------------------------------------------- parameters
def test_the_flank_walk_limit_is_checked_as_a_number(tmp_path):
    code, text = _pro(tmp_path, "--max_flanking_distance", "-5")
    assert code != 0
    assert "Invalid --max_flanking_distance (-5). Must be >= 0 (bp; 0 = no limit)." in text
    code, text = _pro(tmp_path, "--max_flanking_distance", "500000")
    assert code == 0, text[-2000:]


def test_the_removed_settings_warn_and_the_parasail_switch_stops(tmp_path):
    code, text = _pro(tmp_path, "--disable_paralog_check", "true", "--paralog_confusion_min_gap", "3",
                      "--structural_max_length", "500", "--discovery_paralog_margin", "2",
                      "--allow_missing_smith_waterman", "false")
    assert code == 0, text[-2000:]
    for name in ("disable_paralog_check", "paralog_confusion_min_gap", "structural_max_length",
                 "discovery_paralog_margin"):
        assert f"--{name} no longer exists. It had no effect (removed on 2026-11-16" in text, name
    assert ("--allow_missing_smith_waterman no longer exists. Set to false, its old default, it changes "
            "nothing") in text
    code, text = _pro(tmp_path, "--allow_missing_smith_waterman", "true")
    assert code != 0
    assert ("--allow_missing_smith_waterman no longer exists: removed on 2026-11-16: switched on, it let "
            "the search start without parasail") in text
