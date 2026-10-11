"""How far a running search is, read from what the pipeline writes while it runs.

Three sources, all produced by the pipeline version in this checkout:

* its console output: ``[OK  ] LABEL  message`` when a stage is complete (written by
  ``main.nf``) and ``[PROCESS ab/123456] NAME (tag)`` when a task is submitted (written
  by Nextflow; the stage lines alone are enough to tell the steps apart);
* the log of the running genome search (``Wave 2/3 [genomes 2-4/10]`` ...);
* the Nextflow trace file, for the number of finished and failed tasks.

Reading is tolerant by design: text that does not match is ignored, and the worst
outcome of a changed log line is a step that shows less detail. It can never fail a
search. ``tests/test_progress.py`` checks that ``main.nf`` still writes the marker
lines used here.
"""

import re
from pathlib import Path

STEP_KEYS = ("locate", "neighbours", "search", "check", "report")

_ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")
# "[PROCESS ab/123456] NAME (tag)" is what current Nextflow prints without its live
# table; "[ab/123456] Submitted process > NAME (tag)" is the older form.
_PROCESS_RE = re.compile(r"\[(?:PROCESS )?([0-9a-f]{2})/([0-9a-f]{6})\] "
                         r"(?:(?:Submitted|Cached) process > |process > )?([A-Za-z0-9_:]+)")
_OK_RE = re.compile(r"\[OK  \] ([A-Z_]+)[ \t]+([^\[\r\n]*)")
_LOCI_RE = re.compile(r"identified (\d+) locus/loci")
_WAVE_START_RE = re.compile(r"Wave (\d+)/(\d+) \[genomes (\d+)-(\d+)/(\d+)\]")
_WAVE_DONE_RE = re.compile(r"Wave (\d+) complete")

# Marker strings that must exist in main.nf (checked by a test).
MAIN_NF_MARKERS = ("[OK  ]", "'SPLIT_LOCI'.padRight", "'QC_FILTER'.padRight",
                   "'ITERATIVE_SEARCH'.padRight", "'GENERATE_REPORT'.padRight",
                   "identified ${count} locus/loci")
# ... and in bin/iterative_search_runner.py, for the genome counts.
SEARCH_LOG_MARKERS = ("Wave {i+1}/{len(waves)} ", "Wave {i+1} complete")


def read_tail(path, max_bytes):
    """Last ``max_bytes`` of a text file ('' when it cannot be read)."""
    try:
        with open(str(path), "rb") as fh:
            fh.seek(0, 2)
            size = fh.tell()
            fh.seek(max(0, size - max_bytes))
            return fh.read().decode("utf-8", "replace")
    except OSError:
        return ""


def parse_console(text):
    """Submitted tasks and completed stages found in the console output."""
    text = _ANSI_RE.sub("", text or "")
    submitted = [(m.group(1), m.group(2), m.group(3).split(":")[-1])
                 for m in _PROCESS_RE.finditer(text)]
    completed = [(m.group(1), m.group(2).strip()) for m in _OK_RE.finditer(text)]
    loci = None
    for label, message in completed:
        if label == "SPLIT_LOCI":
            match = _LOCI_RE.search(message)
            if match:
                loci = int(match.group(1))
    return {"submitted": submitted, "completed": completed, "loci": loci,
            "success": "[SUCCESS]" in text}


def parse_search_log(text):
    """(genomes searched, genomes in total) from the log of one genome-search task, or
    None when the log does not say."""
    ends, total = {}, None
    for match in _WAVE_START_RE.finditer(text or ""):
        ends[int(match.group(1))] = int(match.group(4))
        total = int(match.group(5))
    if total is None:
        return None
    done = 0
    for match in _WAVE_DONE_RE.finditer(text):
        done = max(done, ends.get(int(match.group(1)), 0))
    return min(done, total), total


def _search_detail(console, work_dir):
    """Genome counts of the search task that is running right now."""
    if not work_dir:
        return None
    for first, second, name in console["submitted"]:
        if name != "ITERATIVE_SEARCH":
            continue
        folder = Path(work_dir) / first
        try:
            matches = [p for p in folder.iterdir() if p.name.startswith(second)]
        except OSError:
            continue
        for task_dir in matches:
            if (task_dir / ".exitcode").exists():
                continue                                   # this one has finished
            counts = parse_search_log(read_tail(task_dir / ".command.log", 400000))
            if counts:
                return counts
    return None


def parse_trace(paths):
    """(tasks finished, tasks failed) over the trace files of all attempts."""
    done = failed = 0
    for path in paths:
        try:
            lines = Path(path).read_text(errors="replace").splitlines()
        except OSError:
            continue
        if not lines:
            continue
        header = lines[0].split("\t")
        if "status" not in header:
            continue
        column = header.index("status")
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) <= column:
                continue
            status = fields[column].strip().upper()
            if status in ("COMPLETED", "CACHED"):
                done += 1
            elif status in ("FAILED", "ABORTED"):
                failed += 1
    return done, failed


def build_progress(console_text, work_dir=None, trace_paths=()):
    """The ``progress`` object a visitor sees for a running search."""
    console = parse_console(console_text)
    names = {name for _, _, name in console["submitted"]}
    labels = [label for label, _ in console["completed"]]
    loci = console["loci"]
    searches_done = labels.count("ITERATIVE_SEARCH")

    done = {
        "locate": "SPLIT_LOCI" in labels or bool(names & {"EXTRACT_FLANKING", "ITERATIVE_SEARCH"}),
        "neighbours": "ITERATIVE_SEARCH" in names or "QC_FILTER" in labels,
        "search": "GENERATE_REPORT" in names or (searches_done >= max(1, loci or 1)),
        "check": "GENERATE_REPORT" in names or "GENERATE_REPORT" in labels,
        "report": console["success"],
    }
    # A later step that is complete implies the earlier ones are.
    reached = False
    for key in reversed(STEP_KEYS):
        reached = reached or done[key]
        done[key] = reached

    steps, active_given = [], False
    for key in STEP_KEYS:
        if done[key]:
            state = "done"
        elif not active_given:
            state, active_given = "active", True
        else:
            state = "waiting"
        steps.append({"key": key, "state": state, "detail": {}})

    by_key = {s["key"]: s for s in steps}
    if loci is not None:
        by_key["locate"]["detail"]["loci"] = loci
    search = by_key["search"]
    if loci and loci > 1:
        search["detail"]["loci_done"] = min(searches_done, loci)
        search["detail"]["loci_total"] = loci
    if search["state"] == "active":
        counts = _search_detail(console, work_dir)
        if counts:
            search["detail"]["genomes_done"], search["detail"]["genomes_total"] = counts

    tasks_done, tasks_failed = parse_trace(trace_paths)
    return {"steps": steps, "tasks_done": tasks_done, "tasks_failed": tasks_failed}


def initial_progress():
    return build_progress("")


def finished_progress(progress):
    """Every step ticked, for a search that has ended successfully."""
    progress = dict(progress or initial_progress())
    progress["steps"] = [dict(step, state="done") for step in progress.get("steps", [])]
    return progress
