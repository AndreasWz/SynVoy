"""Counts tied to n_flanking_genes (profile handoff, stage 6, constraint C9): the bridge limit in
pitches (2n), the HIGH flanking genes of a strong-synteny region (n / 2, rounded up) and those in
the hull of the hull rescue (round(0.4 n)) are 'auto' by default, resolved in Nextflow
(lib/FlankingCounts.groovy) because the scripts take numbers. At the default n = 10 they are 20, 5
and 4, the fixed defaults before, so such a run passes the same command lines. The Nextflow runs
are skipped without Nextflow on PATH (CI has none)."""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")

# param -> (module, flag of its script)
COUNTS = {
    "synteny_bridge_max_gap_genes": ("iterative_search.nf", "--synteny_bridge_max_gap_genes"),
    "strong_synteny_min_flanking": ("cluster_regions.nf", "--strong_synteny_min_flanking"),
    "goi_hull_min_flanking": ("rescue_goi_hull.nf", "--min_flanking"),
}


def test_defaults_are_auto_and_the_modules_resolve_them():
    config = (ROOT / "nextflow.config").read_text()
    groovy = (ROOT / "lib" / "FlankingCounts.groovy").read_text()
    for name, (module, flag) in COUNTS.items():
        assert re.search(rf"^\s*{name}\s*=\s*'auto'", config, re.M), name
        assert f"'{name}'" in groovy, name
        text = (ROOT / "modules" / module).read_text()
        assert f"{flag} ${{FlankingCounts.value(params, '{name}')}}" in text, module
    nf_files = [ROOT / "main.nf", *(ROOT / "modules").glob("*.nf"), *(ROOT / "subworkflows").glob("*.nf")]
    for nf in nf_files:
        text = nf.read_text()
        for name in COUNTS:
            assert f"params.{name}" not in text, f"{nf.name} reads the raw {name}"   # never the string


@needs_nextflow
def test_the_resolved_values(tmp_path):
    (tmp_path / "lib").mkdir()
    shutil.copy(ROOT / "lib" / "FlankingCounts.groovy", tmp_path / "lib" / "FlankingCounts.groovy")
    (tmp_path / "t.nf").write_text(
        "workflow {\n"
        "    [[10, 'auto'], ['10', ' auto '], [20, 'auto'], ['5', 'auto'], [3, 'auto'], [15, 'auto'],\n"
        "     [1, 'auto'], [20, '7'], [20, 0]].each { c ->\n"
        "        def p = [n_flanking_genes: c[0]]\n"
        "        FlankingCounts.NAMES.each { name -> p[name] = c[1] }\n"
        "        println \"RESULT ${c[0]}|${FlankingCounts.NAMES.collect { FlankingCounts.value(p, it) }.join('|')}\"\n"
        "    }\n"
        "}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=300)
    rows = [line[len("RESULT "):] for line in proc.stdout.splitlines() if line.startswith("RESULT ")]
    # bridge 2n | strong ceil(n / 2) | hull round(0.4 n); a number passes as given
    assert rows == ["10|20|5|4", "10|20|5|4", "20|40|10|8", "5|10|3|2", "3|6|2|1", "15|30|8|6",
                    "1|2|1|0", "20|7|7|7", "20|0|0|0"], proc.stdout + proc.stderr


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
    return proc.returncode, re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout + proc.stderr)


@needs_nextflow
@pytest.mark.parametrize("name,value,ok", [
    ("synteny_bridge_max_gap_genes", "auto", True),
    ("strong_synteny_min_flanking", "3", True),
    ("goi_hull_min_flanking", "-1", False),
    ("synteny_bridge_max_gap_genes", "x", False),
])
def test_validation_takes_auto_and_numbers(tmp_path, name, value, ok):
    code, text = _preview(tmp_path, f"--{name}", value)
    assert (f"Invalid --{name}" not in text) is ok, text[-2000:]
    if not ok:
        assert code != 0, text[-2000:]
