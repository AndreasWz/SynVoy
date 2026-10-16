"""Review fix RF-1c (batch of 2026-10-15): a failing blastp is never silent.

synvoy_evalue.blastp_evalue returned E = inf on a nonzero exit, an empty or unparsable
output and a missing blastp, without a word, and the candidate was judged LOW. Each
failure is now a warning (and counted); the runner stops at the start without blastp.
"""
import json
import logging
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import synvoy_evalue as ev  # noqa: E402


def _fake_run(returncode=0, stdout="", stderr=""):
    return lambda *a, **k: SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


def test_each_kind_of_failure_warns(monkeypatch, caplog):
    monkeypatch.setattr(ev, "_blastp", lambda: "/usr/bin/blastp")
    cases = [(_fake_run(1, "", "BLAST engine error: bad"), "exit 1"),
             (_fake_run(0, ""), "no output"),
             (_fake_run(0, "{not json"), "unparsable")]
    for run, words in cases:
        monkeypatch.setattr(ev.subprocess, "run", run)
        before = ev.BLASTP_FAILURES["count"]
        with caplog.at_level(logging.WARNING, logger="synvoy_evalue"):
            assert ev.blastp_evalue("MKTAYIAKQR", ["MKTAYIAKQR"], 1000) == (float("inf"), 0.0)
        assert ev.BLASTP_FAILURES["count"] == before + 1
        assert words in caplog.text
    monkeypatch.setattr(ev, "_blastp", lambda: None)
    with caplog.at_level(logging.WARNING, logger="synvoy_evalue"):
        ev.blastp_evalue("MKTAYIAKQR", ["MKTAYIAKQR"], 1000)
    assert "blastp not found" in caplog.text


def test_no_alignment_is_not_a_failure(monkeypatch, caplog):
    monkeypatch.setattr(ev, "_blastp", lambda: "/usr/bin/blastp")
    empty = json.dumps({"BlastOutput2": [{"report": {"results": {"bl2seq": [{"hits": []}]}}}]})
    monkeypatch.setattr(ev.subprocess, "run", _fake_run(0, empty))
    before = ev.BLASTP_FAILURES["count"]
    with caplog.at_level(logging.WARNING, logger="synvoy_evalue"):
        assert ev.blastp_evalue("MKTAYIAKQR", ["WWWWWWWWWW"], 1000) == (float("inf"), 0.0)
    assert ev.BLASTP_FAILURES["count"] == before and "failed" not in caplog.text


def test_the_runner_stops_without_blastp():
    src = open(os.path.join(os.path.dirname(__file__), "..", "bin", "iterative_search_runner.py")).read()
    assert "if synvoy_evalue._blastp() is None:" in src
