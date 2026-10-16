"""A removed parameter is not accepted and ignored (main.nf removedParams()).

Runs ``nextflow run main.nf -preview`` (builds the DAG, runs no task) on a one-sequence
pro-mode input, about 8 s per run; skipped without Nextflow on PATH (CI has none).

* ``replaced``: the old value asks for something the pipeline can no longer do, so the
  run stops and the message names what does it now;
* ``noEffect``: the run goes on with a warning.

The parameters of the identity-based GOI label (removed 2026-10-15) are the cases.
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


def _preview(tmp_path: Path, *extra: str) -> tuple[int, str]:
    (tmp_path / "q.faa").write_text(">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n")
    (tmp_path / "h.fna").write_text(">c1\nACGTACGTACGTACGTACGTACGTACGT\n")
    noconda = tmp_path / "noconda.config"
    noconda.write_text("conda.enabled = false\nprocess.conda = null\nprocess.executor = 'local'\n")
    proc = subprocess.run(
        [NEXTFLOW, "run", str(ROOT / "main.nf"), "-preview", "-ansi-log", "false",
         "-c", str(noconda), "--mode", "pro", "--query", "q.faa", "--home_genome", "h.fna",
         "--home_species", "Testus_testus", "--outdir", "out", *extra],
        cwd=tmp_path, capture_output=True, text=True, timeout=300)
    text = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout + proc.stderr)
    return proc.returncode, text


def test_a_replaced_parameter_stops_the_run(tmp_path):
    code, text = _preview(tmp_path, "--classify_high_min_identity", "50")
    assert code != 0
    assert ("--classify_high_min_identity no longer exists: a GOI call's label comes from "
            "one judge") in text
    assert "docs/ALGORITHM.md 7.8" in text


def test_a_parameter_without_effect_only_warns(tmp_path):
    code, text = _preview(tmp_path, "--distance_autotune_far_pct", "40",
                          "--strict_goi_family", "true")
    assert code == 0, text[-2000:]
    assert ("--distance_autotune_far_pct no longer exists. It had no effect; remove it "
            "from the command line or config.") in text
    assert ("--strict_goi_family no longer exists. It had no effect (the family-name check "
            "is annotation only now") in text
    assert "All parameters validated" in text
