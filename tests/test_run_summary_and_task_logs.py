"""The end summary and the task logs of a run (items 32 and 38 of the bug list of 16 Nov 2026).

Item 32: the summary listed regions/*.regions.bed at the top level of regions/, but CLUSTER_REGIONS
publishes regions/<locus>/, so the line never printed; every *_tree.nwk was "(GOI phylogeny)", also
the placeholder of --skip_tree (which the launcher's default low_mem sets); and the main figure,
the anchor grid, was not named. lib/RunSummary.groovy lists the region tables per home locus and
tells a tree from a placeholder (a tree has its <locus>_tree.info.json, not placeholder_tree).
Review of the series (fix round 1): the info must be this run's. A placeholder tree published into a
reused output folder sat next to an earlier run's info and was named "(GOI phylogeny)"; an info that
the caller's fromThisRun refuses, or that is older than its tree, now counts as absent.

Item 38: logs/ was filled by a recursive walk of the whole work folder, so it held the tasks of
earlier runs (and walked into work/conda/), with __2, __3 counters in the file system's listing
order. lib/TaskLogs.groovy takes only the task folders (<2 hex>/<30 hex>/) whose .command.run this
session wrote, names each folder by process, first output and task hash, and copies non-empty logs.

The Groovy checks run the classes in a small Nextflow script; skipped without Nextflow on PATH.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
HEX30 = "0123456789abcdef0123456789abcd"


def _groovy(tmp_path: Path, body: str) -> tuple[list[str], str]:
    (tmp_path / "lib").mkdir(exist_ok=True)
    for name in ("RunSummary.groovy", "TaskLogs.groovy"):
        shutil.copy(ROOT / "lib" / name, tmp_path / "lib" / name)
    (tmp_path / "t.nf").write_text("workflow {\n" + body + "\n}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=300)
    rows = [line[len("RESULT "):] for line in proc.stdout.splitlines() if line.startswith("RESULT ")]
    return rows, proc.stdout + proc.stderr


@needs_nextflow
def test_the_region_tables_are_listed_per_home_locus(tmp_path):
    out = tmp_path / "out"
    for rel in ("regions/locus_2/b.regions.bed", "regions/locus_1/a.regions.bed", "regions/locus_1/c.txt",
                "regions/x.regions.bed"):
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        (out / rel).write_text("c1\t1\t2\n")
    rows, log = _groovy(tmp_path, f"""
    def beds = RunSummary.regionBeds(new File('{out}'), {{ f -> true }})
    println "RESULT ${{beds.collect {{ b -> b.parentFile.name + '/' + b.name }}.join(',')}}"
    println "RESULT ${{RunSummary.regionBeds(new File('{out}'), {{ f -> f.name.startsWith('a') }}).size()}}"
    println "RESULT ${{RunSummary.regionBeds(new File('{tmp_path}/none'), {{ f -> true }}).size()}}"
""")
    assert rows == ["locus_1/a.regions.bed,locus_2/b.regions.bed", "1", "0"], log


@needs_nextflow
def test_a_placeholder_tree_is_told_from_a_tree(tmp_path):
    for locus, info in (("locus_1", None), ("locus_2", '{"bootstrap": {"status": "placeholder_tree"}}'),
                        ("locus_3", '{"bootstrap": {"status": "ultrafast"}}'), ("locus_4", "{not json")):
        (tmp_path / f"{locus}_tree.nwk").write_text("(a,b,c);\n")
        if info is not None:
            (tmp_path / f"{locus}_tree.info.json").write_text(info)
    rows, log = _groovy(tmp_path, "    println \"RESULT ${(1..4).collect { i -> RunSummary.treeKind(new File('"
                        + str(tmp_path) + "/locus_' + i + '_tree.nwk')) }.join(',')}\"")
    assert rows == ["placeholder,placeholder,tree,placeholder"], log


@needs_nextflow
def test_an_earlier_runs_info_does_not_make_a_placeholder_a_tree(tmp_path):
    """A placeholder written by this run next to the info of an earlier run's real tree (the output folder
    reused, with or without -resume) is a placeholder; the info of this run's tree, written together with it,
    keeps it a tree."""
    now = time.time()
    ultrafast = '{"bootstrap": {"status": "ultrafast"}}'
    cases = (("locus_1", now, now - 2 * 86400),     # fresh placeholder, info two days old
             ("locus_2", now, now - 60),            # info one minute older than its tree
             ("locus_3", now, now - 1),             # both from one task: within the slack
             ("locus_4", now - 5, now))             # info written after its tree
    for locus, nwk_time, info_time in cases:
        nwk, info = tmp_path / f"{locus}_tree.nwk", tmp_path / f"{locus}_tree.info.json"
        nwk.write_text("(GOI_placeholder:0.0);\n")
        info.write_text(ultrafast)
        os.utime(nwk, (nwk_time, nwk_time))
        os.utime(info, (info_time, info_time))
    since = int((now - 3600) * 1000)                  # this "run" started an hour ago
    rows, log = _groovy(tmp_path, f"""
    println "RESULT ${{RunSummary.treeKind(new File('{tmp_path}/locus_1_tree.nwk'))}}"
    def kinds = {{ Closure keep -> (1..4).collect {{ i -> RunSummary.treeKind(new File('{tmp_path}/locus_' + i + '_tree.nwk'), keep) }}.join(',') }}
    println "RESULT ${{kinds.call({{ File f -> true }})}}"
    println "RESULT ${{kinds.call({{ File f -> f.lastModified() >= {since}L }})}}"
""")
    assert rows[:1] == ["placeholder"], log          # was "tree": the earlier run's info was read
    # keep everything (as fromThisRun does on -resume), or fromThisRun on a fresh run (every info but
    # locus_1's is of this run): the age against the tree decides the same
    assert rows[1:] == ["placeholder,placeholder,tree,tree"] * 2, log


@needs_nextflow
def test_the_callers_keep_refuses_an_info(tmp_path):
    (tmp_path / "locus_1_tree.nwk").write_text("(a,b,c);\n")
    (tmp_path / "locus_1_tree.info.json").write_text('{"bootstrap": {"status": "ultrafast"}}')
    rows, log = _groovy(tmp_path, f"""
    def nwk = new File('{tmp_path}/locus_1_tree.nwk')
    println "RESULT ${{RunSummary.treeKind(nwk, {{ File f -> true }})}},${{RunSummary.treeKind(nwk, {{ File f -> !f.name.endsWith('.info.json') }})}}"
""")
    assert rows == ["tree,placeholder"], log


def _task(work: Path, top: str, rest: str, name: str, output: str | None, log: str = "out\n",
          err: str = "", age_s: float = 0.0) -> Path:
    task = work / top / rest
    task.mkdir(parents=True)
    header = ["#!/bin/bash", "### ---", f"### name: '{name}'"]
    if output:
        header += ["### outputs:", f"### - '{output}'"]
    (task / ".command.run").write_text("\n".join(header + ["### ...", "set -e"]) + "\n")
    (task / ".command.log").write_text(log)
    (task / ".command.err").write_text(err)
    (task / ".command.sh").write_text("#!/bin/bash -ue\necho hi\n")
    if age_s:
        old = time.time() - age_s
        os.utime(task / ".command.run", (old, old))
    return task


@needs_nextflow
def test_only_this_sessions_tasks_reach_logs(tmp_path):
    work = tmp_path / "work"
    _task(work, "ab", HEX30, "PLOT_SYNTENY (1)", "x_synteny_plot.html", log="plot\n")
    _task(work, "cd", HEX30, "PHYLO_SORT (sorting)", "sorted_genomes.txt", log="", err="warn\n")
    _task(work, "ef", HEX30, "GENERATE_REPORT", None, log="", err="")            # empty logs: no folder
    _task(work, "01", HEX30, "OLD_RUN (1)", "old.txt", age_s=86400)            # an earlier session's task
    _task(work / "conda" / "env-x", "ab", HEX30, "CONDA_INSIDE (1)", "c.txt")  # never entered
    since = int((time.time() - 3600) * 1000)
    rows, log = _groovy(tmp_path, f"""
    def n1 = TaskLogs.collect(new File('{work}'), new File('{tmp_path}/logs1'), {since}L)
    def n2 = TaskLogs.collect(new File('{work}'), new File('{tmp_path}/logs2'), {since}L)
    println "RESULT ${{n1}},${{n2}}"
""")
    assert rows == ["2,2"], log
    first = sorted(p.name for p in (tmp_path / "logs1").iterdir())
    assert first == sorted(p.name for p in (tmp_path / "logs2").iterdir())
    assert first == ["PHYLO_SORT__sorting___sorted_genomes.txt__cd012345",
                     "PLOT_SYNTENY__1___x_synteny_plot.html__ab012345"], first
    plot = tmp_path / "logs1" / "PLOT_SYNTENY__1___x_synteny_plot.html__ab012345"
    assert sorted(p.name for p in plot.iterdir()) == ["script.sh", "stdout.log"]       # empty stderr left out
    assert (plot / "stdout.log").read_text() == "plot\n"
    sort = tmp_path / "logs1" / "PHYLO_SORT__sorting___sorted_genomes.txt__cd012345"
    assert sorted(p.name for p in sort.iterdir()) == ["script.sh", "stderr.log"]


def test_main_nf_uses_the_classes():
    main = (ROOT / "main.nf").read_text()
    assert "eachDirRecurse" not in main and "def collectTaskLogs" not in main
    assert "TaskLogs.collect(new File(workflow.workDir.toString())" in main
    assert "RunSummary.regionBeds(outdir" in main and "RunSummary.treeKind(t, { f -> fromThisRun(f) })" in main
    assert "(main figure)" in main and "(placeholder, " in main
    assert not re.search(r"regions_dir\.listFiles\(\)", main)


def test_replot_run_still_finds_the_log_folders():
    """scripts/replot_run.sh reads the home species from logs/PLOT_SYNTENY*/ and logs/PHYLO_SORT*/."""
    from fnmatch import fnmatch
    text = (ROOT / "scripts" / "replot_run.sh").read_text()
    assert '"$RUN"/logs/PLOT_SYNTENY*/script.sh' in text and '"$RUN"/logs/PHYLO_SORT*/script.sh' in text
    assert fnmatch("PLOT_SYNTENY__1___x_synteny_plot.html__ab012345", "PLOT_SYNTENY*")
    assert fnmatch("PHYLO_SORT__sorting___sorted_genomes.txt__cd012345", "PHYLO_SORT*")
