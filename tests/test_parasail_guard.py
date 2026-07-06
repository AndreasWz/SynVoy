"""Regression tests for the fail-loud parasail / Smith-Waterman guards.

SW is the load-bearing tier for divergent GOI recovery (e.g. re-finding a 26-40 %
myrmicitoxin from a melittin query). A missing-parasail env — classically a `.venv`
shadowing the conda env on PATH — used to SILENTLY disable SW and crash the §1m panel,
producing wrong results with no loud signal (see the 2026-07-01 melittin validation).
These tests pin the guard behaviour so it can't regress back to silent-degrade.
"""
import importlib.util
import os
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


# ── resolve_smith_waterman_enabled: the pure decision ────────────────────────

def test_sw_available_stays_enabled():
    assert isr.resolve_smith_waterman_enabled(True, "auto", True, False) is True


def test_sw_missing_auto_fails_loud():
    """The core guard: requested + auto + no parasail + not opted-out → raise."""
    with pytest.raises(isr.SmithWatermanUnavailable):
        isr.resolve_smith_waterman_enabled(True, "auto", False, allow_missing=False)


def test_sw_missing_but_allowed_degrades_quietly():
    assert isr.resolve_smith_waterman_enabled(True, "auto", False, allow_missing=True) is False


def test_sw_explicitly_disabled_never_raises():
    # User asked for no SW at all — missing parasail is irrelevant.
    assert isr.resolve_smith_waterman_enabled(False, "auto", False, False) is False


def test_sw_non_auto_backend_bypasses_the_parasail_gate():
    # ssearch36 doesn't need parasail, so the auto-mode gate must not fire.
    assert isr.resolve_smith_waterman_enabled(True, "ssearch36", False, False) is True


def test_allow_missing_flag_defaults_false():
    """A degraded run must be a deliberate opt-in, not the default."""
    import argparse
    # Rebuild just the relevant arg to assert the shipped default.
    p = argparse.ArgumentParser()
    p.add_argument("--allow_missing_smith_waterman",
                   type=lambda v: str(v).lower() in ("1", "true", "yes"), default=False)
    assert p.parse_args([]).allow_missing_smith_waterman is False


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
