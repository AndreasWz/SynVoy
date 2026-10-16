"""locus_gap_genes follows n_flanking_genes (batch of 2026-10-15, HG-1C): 'auto' (the default)
is floor(n_flanking_genes / 2), 5 at the default of 10, resolved in Nextflow (lib/LocusGap.groovy)
because split_loci.py and profile_hits.py take a number. The Nextflow runs are skipped without
Nextflow on PATH (CI has none)."""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")


def test_default_is_auto_and_both_modules_resolve_it():
    config = (ROOT / "nextflow.config").read_text()
    assert re.search(r"^\s*locus_gap_genes\s*=\s*'auto'", config, re.M)
    for module in ("split_loci.nf", "locate_gene.nf"):
        text = (ROOT / "modules" / module).read_text()
        assert "--locus_gap_genes ${LocusGap.genes(params)}" in text
        assert "params.locus_gap_genes" not in text          # never the raw string


@needs_nextflow
def test_the_resolved_values(tmp_path):
    (tmp_path / "lib").mkdir()
    shutil.copy(ROOT / "lib" / "LocusGap.groovy", tmp_path / "lib" / "LocusGap.groovy")
    (tmp_path / "t.nf").write_text(
        "workflow {\n"
        "    [[10, 'auto'], ['7', 'auto'], [1, 'auto'], [10, 4], [10, '6'], [10, ' auto ']].each { c ->\n"
        "        def p = [n_flanking_genes: c[0], locus_gap_genes: c[1]]\n"
        "        println \"RESULT ${c[0]}|${c[1]}|${LocusGap.genes(p)}|${LocusGap.isAuto(p)}\"\n"
        "    }\n"
        "}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=300)
    rows = [line[len("RESULT "):].split("|") for line in proc.stdout.splitlines()
            if line.startswith("RESULT ")]
    assert [(n, g.strip(), int(v), a) for n, g, v, a in rows] == [
        ("10", "auto", 5, "true"), ("7", "auto", 3, "true"), ("1", "auto", 0, "true"),
        ("10", "4", 4, "false"), ("10", "6", 6, "false"), ("10", "auto", 5, "true")], proc.stdout + proc.stderr


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
@pytest.mark.parametrize("value,ok", [("auto", True), ("3", True), ("-1", False)])
def test_validation_takes_auto_and_numbers(tmp_path, value, ok):
    code, text = _preview(tmp_path, "--locus_gap_genes", value)
    assert ("Invalid --locus_gap_genes" not in text) is ok, text[-2000:]
    if not ok:
        assert code != 0 and "Must be 'auto' or >= 0" in text
