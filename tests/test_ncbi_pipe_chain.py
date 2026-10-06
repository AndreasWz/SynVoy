"""The NCBI query chain (esearch | efetch | xtract) must not hang on its own stderr.

run_piped_command() reads only the LAST stage. The stderr of the earlier stages was
a pipe that nobody drained, so a stage that wrote more than the pipe buffer (~64 kB)
blocked for good. entrez-direct echoes the whole request for every call it has to
retry: 23.5 kB for one 1,235-assembly query on a connection where each request
failed once (2026-10-05). Each stage now writes to a temporary file, and the tail
is reported when the chain fails or times out.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import fetch_related_genomes as fr  # noqa: E402

PASS_THROUGH = [sys.executable, "-c", "import sys; sys.stdout.write(sys.stdin.read())"]


def _stage(code: str):
    return [sys.executable, "-c", code]


def test_chain_survives_an_upstream_stage_that_floods_stderr():
    noisy = _stage("import sys; sys.stderr.write('x' * 300000); sys.stderr.flush(); print('payload')")
    assert fr.run_piped_command([noisy, PASS_THROUGH], timeout=30) == "payload"


def test_chain_returns_the_last_stage_output():
    first = _stage("print('a\\tb')")
    upper = _stage("import sys; sys.stdout.write(sys.stdin.read().upper())")
    assert fr.run_piped_command([first, upper], timeout=30) == "A\tB"


def test_timeout_reports_what_the_stalled_stage_said(capsys):
    stuck = _stage("import sys, time; sys.stderr.write('curl: (56) connection dropped\\n');"
                   " sys.stderr.flush(); time.sleep(60)")
    assert fr.run_piped_command([stuck, PASS_THROUGH], timeout=1) is None
    err = capsys.readouterr().err
    assert "timed out after 1s" in err
    assert "curl: (56) connection dropped" in err


def test_failing_last_stage_returns_none_and_reports_its_stderr(capsys):
    fail = _stage("import sys; sys.stderr.write('xtract: bad pattern\\n'); sys.exit(3)")
    assert fr.run_piped_command([_stage("print('x')"), fail], timeout=30) is None
    assert "xtract: bad pattern" in capsys.readouterr().err
