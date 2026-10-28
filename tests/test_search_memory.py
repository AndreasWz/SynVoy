"""The memory request of ITERATIVE_SEARCH (profile stage 5b; it changes no result).

One task searches every target genome; the genomes of a wave are searched at once, each by a
worker with its own MMseqs2 run (at most mmseqs_split_memory_limit) and the whole genome in
memory. The request was a fixed 10 GB, and a retry after an out-of-memory kill asked for the same
10 GB again. Now (lib/SearchMemory.groovy) iterative_search_memory 'auto' is max(10 GB, the need of
the largest genomes searched at once), x the attempt, and under the local executor at most the
executor's memory, which Nextflow refuses to exceed. A size the user sets is the base.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")


def _groovy(tmp_path, body):
    (tmp_path / "lib").mkdir(exist_ok=True)
    shutil.copy(ROOT / "lib" / "SearchMemory.groovy", tmp_path / "lib" / "SearchMemory.groovy")
    (tmp_path / "t.nf").write_text("workflow {\n" + body + "\n}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=300)
    rows = [line[len("RESULT "):] for line in proc.stdout.splitlines() if line.startswith("RESULT ")]
    return rows, proc.stdout + proc.stderr


def _sparse(folder: Path, sizes):
    folder.mkdir(parents=True, exist_ok=True)
    for i, size in enumerate(sizes):
        with open(folder / f"g{i}.fna", "wb") as fh:
            fh.truncate(int(size))
    return folder


@needs_nextflow
def test_the_largest_wave_is_the_runners(tmp_path):
    rows, log = _groovy(tmp_path, "    println \"RESULT ${(1..40).collect { n -> SearchMemory.largestWave(n, true, false) }.join(',')}\"\n"
                                  "    println \"RESULT ${SearchMemory.largestWave(12, true, true)},${SearchMemory.largestWave(12, false, false)},${SearchMemory.largestWave(3, false, false)}\"")
    runner = [max(len(w) for w in isr.assign_waves_by_rank(list(range(n)))) for n in range(1, 41)]
    assert rows[0] == ",".join(map(str, runner)), log
    assert rows[1] == "12,5,3", log                     # one wave; the old binning's waves of 5


@needs_nextflow
def test_the_request_for_attempts_one_to_three(tmp_path):
    """A 15 GB genome next to a small one: base = 8 GB split + 15 GB + 1 + 1 GB = 25 GB, x the
    attempt where nothing caps it (a scheduler), and at most the local executor's memory."""
    g = _sparse(tmp_path / "genomes", [15 * 2**30, 400 * 2**20])
    body = f"""
    def p = [iterative_search_memory: 'auto', mmseqs_split_memory_limit: '8G', rank_wave_binning: true,
             disable_wavefront: false, iterative_search_cpus: 2]
    def g = new File('{g}')
    (1..3).each {{ a -> println "RESULT slurm ${{SearchMemory.request(p, g, [attempt: a, executor: 'slurm', cpus: 2]).memory.toBytes()}}" }}
    def local = SearchMemory.localMemory().toBytes()
    (1..3).each {{ a -> println "RESULT local ${{SearchMemory.request(p, g, [attempt: a, executor: 'local', cpus: 2]).memory.toBytes() == Math.min(local, a * 25L * (1L << 30))}}" }}
    println "RESULT line ${{SearchMemory.request(p, g, [attempt: 2, executor: 'slurm', cpus: 2]).line}}"
    """
    rows, log = _groovy(tmp_path, body)
    gib = 2**30
    assert rows[:3] == [f"slurm {a * 25 * gib}" for a in (1, 2, 3)], log
    assert rows[3:6] == ["local true"] * 3, log
    assert rows[6] == ("line ITERATIVE_SEARCH memory: 50 GB (attempt 2 x 25 GB; max(10 GB, 25 GB for 1 "
                       "genome(s) at once, the largest 15 GB))"), log


@needs_nextflow
def test_small_genomes_keep_the_old_ten_gigabytes(tmp_path):
    """Candida-sized targets at the low_mem split: the need is below the floor, so attempt 1 asks
    for 10 GB, as before."""
    g = _sparse(tmp_path / "genomes", [14.6e6] * 5)
    body = f"""
    def p = [iterative_search_memory: 'auto', mmseqs_split_memory_limit: '2G', rank_wave_binning: true,
             disable_wavefront: false, iterative_search_cpus: 4]
    println "RESULT ${{SearchMemory.request(p, new File('{g}'), [attempt: 1, executor: 'slurm', cpus: 4]).memory}}"
    println "RESULT ${{SearchMemory.need(p, new File('{g}'), 4).at_once}}"
    """
    rows, log = _groovy(tmp_path, body)
    assert rows == ["10 GB", "2"], log


@needs_nextflow
def test_genomes_searched_at_once_and_a_value_the_user_set(tmp_path):
    """Three 2 GB genomes of one wave at split 16G need 3 x (16 + 2 + 1) + 1 = 58 GB (LRZ measured
    37-38 GiB with these sizes); gzip counts 4 x; a size the user set is the base."""
    g = _sparse(tmp_path / "genomes", [2 * 2**30] * 9)
    gz = tmp_path / "gz"
    gz.mkdir()
    with open(gz / "a.fna.gz", "wb") as fh:
        fh.truncate(2**30)
    (gz / "notes.txt").write_text("x" * 10**6)
    body = f"""
    def p = [iterative_search_memory: 'auto', mmseqs_split_memory_limit: '16G', rank_wave_binning: true,
             disable_wavefront: false, iterative_search_cpus: 16]
    println "RESULT ${{SearchMemory.request(p, new File('{g}'), [attempt: 1, executor: 'slurm', cpus: 16]).memory}}"
    println "RESULT ${{SearchMemory.targetBytes(new File('{gz}'))}}"
    def u = p + [iterative_search_memory: '6 GB']
    println "RESULT ${{(1..3).collect {{ a -> SearchMemory.request(u, new File('{g}'), [attempt: a, executor: 'slurm', cpus: 16]).memory }}.join(',')}}"
    """
    rows, log = _groovy(tmp_path, body)
    assert rows == ["58 GB", f"[{4 * 2**30}]", "6 GB,12 GB,18 GB"], log


def test_the_module_sets_the_memory_and_the_config_does_not():
    nf = (ROOT / "modules" / "iterative_search.nf").read_text()
    assert "memory { def m = SearchMemory.request(params, genomes, task); log.info m.line; m.memory }" in nf
    config = (ROOT / "nextflow.config").read_text()
    block = config[config.index("withName: 'ITERATIVE_SEARCH' {", config.index("// Process-specific")):]
    block = block[:block.index("}")]
    assert not re.search(r"^\s*memory\s*=", block, re.M), block
    assert re.search(r"^\s*iterative_search_memory = 'auto'$", config, re.M)
    # the memory tiers keep their sizes: a size the user (or a profile) sets is the base
    assert "params.iterative_search_memory = '8 GB'" in config and "params.iterative_search_memory = '4 GB'" in config


@needs_nextflow
@pytest.mark.parametrize("value,ok", [("auto", True), ("24 GB", True), ("6.GB", True), ("lots", False)])
def test_validation(tmp_path, value, ok):
    (tmp_path / "q.faa").write_text(">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n")
    (tmp_path / "h.fna").write_text(">c1\nACGTACGTACGTACGTACGTACGTACGT\n")
    (tmp_path / "noconda.config").write_text("conda.enabled = false\nprocess.conda = null\nprocess.executor = 'local'\n")
    proc = subprocess.run(
        [NEXTFLOW, "run", str(ROOT / "main.nf"), "-preview", "-ansi-log", "false", "-c", "noconda.config",
         "--mode", "pro", "--query", "q.faa", "--home_genome", "h.fna", "--home_species", "Testus_testus",
         "--outdir", "out", "--iterative_search_memory", value],
        cwd=tmp_path, capture_output=True, text=True, timeout=300)
    text = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout + proc.stderr)
    assert ("Invalid --iterative_search_memory" not in text) is ok, text[-2000:]
