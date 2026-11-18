"""Regression tests for the fail-loud parasail / Smith-Waterman guards.

SW is the load-bearing tier for divergent GOI recovery (e.g. re-finding a 26-40 %
myrmicitoxin from a melittin query). A missing-parasail env — classically a `.venv`
shadowing the conda env on PATH — used to SILENTLY disable SW and crash the §1m panel,
producing wrong results with no loud signal (see the 2026-07-01 melittin validation).
These tests pin the guard behaviour so it can't regress back to silent-degrade.
"""
import importlib.util
import os
import shutil
import subprocess
import sys

import pytest

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")


def _load(mod_name, filename):
    spec = importlib.util.spec_from_file_location(mod_name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


isr = _load("iterative_search_runner", "iterative_search_runner.py")


# ── parasail_requirement_error: every search needs parasail ─────────────────
# --allow_missing_smith_waterman let the search start
# without parasail, but the labels' bit scores and specificity terms, the hull rescue and the
# ownership panel need it too, so such a run lost them and stopped at the panel after the whole
# search. The option is gone (main.nf removedParams, offOnly), and the search stops at its start.

def test_with_parasail_nothing_stops():
    assert isr.parasail_requirement_error(True) is None
    assert isr.parasail_requirement_error(True, "/home/u/.venv") is None


def test_without_parasail_the_search_stops_with_the_reason_and_the_fix():
    msg = isr.parasail_requirement_error(False)
    assert "parasail is not importable" in msg and "Every search needs it" in msg
    assert "conda install -c bioconda parasail-python" in msg
    assert "VIRTUAL_ENV" not in msg and "allow_missing" not in msg


def test_an_active_virtualenv_is_named():
    assert "VIRTUAL_ENV=/home/u/.venv" in isr.parasail_requirement_error(False, "/home/u/.venv")


def test_the_degraded_switch_is_gone():
    import argparse
    assert not hasattr(isr, "resolve_smith_waterman_enabled")
    with open(os.path.join(BIN, "iterative_search_runner.py")) as fh:
        text = fh.read()
    assert "--allow_missing_smith_waterman" not in text
    main = text[text.index("def main():"):]
    # the check runs whatever --enable_smith_waterman and --sw_method say
    assert "parasail_requirement_error(has_parasail_available()" in main


@pytest.mark.skipif(not (shutil.which("blastp") and shutil.which("miniprot")),
                    reason="the runner checks blastp and miniprot before parasail")
def test_the_search_stops_before_its_first_genome_without_parasail(tmp_path):
    """The runner under a python that cannot import parasail: exit 1 with the reason, also with
    --enable_smith_waterman false (the bit scores and the specificity still need it)."""
    shim = tmp_path / "shim"
    (shim / "parasail").mkdir(parents=True)
    (shim / "parasail" / "__init__.py").write_text("raise ImportError('shadowed for test')\n")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(shim) + os.pathsep + BIN + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("VIRTUAL_ENV", None)
    db = tmp_path / "initial.faa"
    db.write_text(">GOI_q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n")
    (tmp_path / "genomes").mkdir()
    (tmp_path / "sorted.txt").write_text("")
    proc = subprocess.run(
        [sys.executable, os.path.join(BIN, "iterative_search_runner.py"),
         "--initial_db", str(db), "--sorted_genomes", str(tmp_path / "sorted.txt"),
         "--genomes_dir", str(tmp_path / "genomes"), "--output_dir", str(tmp_path / "out"),
         "--enable_smith_waterman", "false"],
        capture_output=True, text=True, env=env, timeout=300)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 1, out[-2000:]
    assert "parasail is not importable in this environment. Every search needs it" in out


# ── build_home_paralog_panel: §1m panel fails loud without parasail ──────────

def test_paralog_panel_exits_loud_when_parasail_missing(tmp_path):
    """Run the panel script under a python env that cannot import parasail and
    confirm it exits non-zero with an actionable message (not a silent empty panel)."""
    # A shim package dir whose `parasail/__init__.py` raises ImportError, placed first
    # on sys.path so `import parasail` fails exactly like the .venv-shadow scenario.
    shim = tmp_path / "shim"
    (shim / "parasail").mkdir(parents=True)
    (shim / "parasail" / "__init__.py").write_text("raise ImportError('shadowed for test')\n")

    gff = tmp_path / "home.gff"
    gff.write_text("##gff-version 3\nchr1\t.\tgene\t1\t99\t.\t+\t.\tID=gene-A\n")
    prot = tmp_path / "home.faa"
    prot.write_text(">gene-A\nACDEFGHIKLMNPQRSTVWY\n")
    bed = tmp_path / "locus_1.bed"
    bed.write_text("chr1\t1\t99\tlocus_1\n")
    query = tmp_path / "q.faa"
    query.write_text(">q\nACDEFGHIKLMNPQRSTVWY\n")

    env = dict(os.environ)
    env["PYTHONPATH"] = str(shim) + os.pathsep + BIN + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, os.path.join(BIN, "build_home_paralog_panel.py"),
         "--home_gff", str(gff), "--home_proteome", str(prot),
         "--locus_beds", str(bed), "--query", str(query),
         "--output", str(tmp_path / "panel.faa")],
        capture_output=True, text=True, env=env,
    )
    assert proc.returncode != 0, "panel must fail loud, not exit 0 with an empty panel"
    assert "parasail" in (proc.stderr + proc.stdout).lower()
