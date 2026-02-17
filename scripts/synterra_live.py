#!/usr/bin/env python3
"""
Run SynTerra with a controlled, modern terminal UI.

This wrapper executes Nextflow with ANSI task table disabled and renders a
single live status line with spinner + compact state.
"""

from __future__ import annotations

import argparse
import os
import queue
import re
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from shutil import get_terminal_size, which


ANSI_RE = re.compile(r"\x1B\[[0-?]*[ -/]*[@-~]")
STATUS_RE = re.compile(r"\[(RUN |OK  |WARN|FAIL)\]\s+([A-Z0-9_]+)\s*(.*)")
PHASE_RE = re.compile(r"^Phase\s+(\d+)\s*\|\s*(.+)$")
PROC_RE = re.compile(r"process\s+>\s+([A-Z0-9_]+)\b")
TABLE_PROC_RE = re.compile(r"^\[[^]]+\]\s+(?:process\s+>\s+)?([A-Z0-9_]+)\b")
SUBMIT_RE = re.compile(r"Submitted process >\s+([A-Z0-9_]+)\b")
PROC_COUNTS_RE = re.compile(r"([0-9]+)\s+of\s+([0-9]+)")
DURATION_RE = re.compile(r"Duration:\s*(.+)$")
TASKS_RE = re.compile(r"Tasks Completed:\s*([0-9]+)")
ERROR_RE = re.compile(r"^ERROR\s*~\s*(.+)$")

SPINNER = ["|", "/", "-", "\\"]

PHASE_BY_TASK = {
    "RESOLVE_GENE_INPUT": "P0 Setup",
    "FETCH_QUERY_FROM_ID": "P0 Setup",
    "FETCH_HOME_GENOME": "P0 Setup",
    "FETCH_RELATED_GENOMES": "P0 Setup",
    "NORMALIZE_QUERY": "P0 Setup",
    "STAGE_GENOMES": "P0 Setup",
    "LOCATE_GENE": "P1 Gene Localization in Home Genome",
    "ANNOTATE_GOI": "P1 Gene Localization in Home Genome",
    "SPLIT_LOCI": "P1 Gene Localization in Home Genome",
    "PREPARE_HOME": "P1 Gene Localization in Home Genome",
    "PREPARE_HOME_PROTEOME": "P1 Gene Localization in Home Genome",
    "BORROW_ANNOT": "P1 Gene Localization in Home Genome",
    "BORROW_ANNOTATIONS": "P1 Gene Localization in Home Genome",
    "EXTRACT_FLANKING": "P1 Gene Localization in Home Genome",
    "PREPARE_DB": "P1 Gene Localization in Home Genome",
    "PREPARE_INITIAL_DB": "P1 Gene Localization in Home Genome",
    "PHYLO_SORT": "P2 Phylogenetic Ordering and Iterative Search",
    "GENOME_QC": "P2 Phylogenetic Ordering and Iterative Search",
    "ASSESS_GENOME_QUALITY": "P2 Phylogenetic Ordering and Iterative Search",
    "ITERATIVE_SEARCH": "P2 Phylogenetic Ordering and Iterative Search",
    "CLUSTER_REGIONS": "P3 Region Clustering",
    "COMPUTE_TREE": "P4 Phylogenetics and Visualization",
    "PLOT_SYNTENY": "P4 Phylogenetics and Visualization",
    "GENERATE_REPORT": "P5 Report Generation",
}


def strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def set_phase_for_task(state: "LiveState", task: str) -> None:
    phase = PHASE_BY_TASK.get(task)
    if phase:
        state.phase = phase


class LiveState:
    def __init__(self) -> None:
        self.phase = "Init"
        self.current_task = "Starting"
        self.task_detail = ""
        self.progress_label = ""
        self.progress_counts = ""
        self.last_info = ""
        self.last_error = ""
        self.pipeline_duration = ""
        self.pipeline_tasks = ""
        self.success = False
        self.events: list[str] = []

    def add_event(self, text: str) -> None:
        if not text:
            return
        self.events.append(text)
        if len(self.events) > 6:
            self.events = self.events[-6:]


def should_skip_line(line: str) -> bool:
    if not line:
        return True
    if line.startswith("executor >"):
        return True
    if line.startswith("WARN: The operator `first` is useless"):
        return True
    if line.startswith("curl: (6) Could not resolve host: www.nextflow.io"):
        return True
    if line.startswith("N E X T F L O W"):
        return True
    if line.startswith("Launching `main.nf`"):
        return True
    # old-style task table lines
    if line.startswith("[") and "process >" in line:
        return True
    if line.startswith("[-        ]"):
        return True
    return False


def update_state(state: LiveState, line: str) -> None:
    if not line:
        return

    m = PHASE_RE.search(line)
    if m:
        state.phase = f"P{m.group(1)} {m.group(2)}"
        state.last_info = line
        state.add_event(line)
        return

    m = STATUS_RE.search(line)
    if m:
        level, task, detail = m.group(1).strip(), m.group(2), m.group(3).strip()
        if level == "RUN":
            state.current_task = task
            state.task_detail = detail
            set_phase_for_task(state, task)
        elif level == "OK":
            msg = f"{task} done"
            if detail:
                msg = f"{msg}: {detail}"
            state.last_info = msg
            state.add_event(msg)
        elif level in {"WARN", "FAIL"}:
            msg = f"{level} {task}"
            if detail:
                msg = f"{msg}: {detail}"
            state.last_info = msg
            state.add_event(msg)
        return

    m = PROC_RE.search(line)
    if m:
        proc = m.group(1)
        state.progress_label = proc
        state.current_task = proc
        set_phase_for_task(state, proc)
        c = PROC_COUNTS_RE.search(line)
        if c:
            state.progress_counts = f"{c.group(1)}/{c.group(2)}"
        elif "[100%]" in line:
            state.progress_counts = "100%"
        return

    m = TABLE_PROC_RE.search(line)
    if m:
        proc = m.group(1)
        state.progress_label = proc
        state.current_task = proc
        set_phase_for_task(state, proc)
        c = PROC_COUNTS_RE.search(line)
        if c:
            state.progress_counts = f"{c.group(1)}/{c.group(2)}"
        elif "[100%]" in line:
            state.progress_counts = "100%"
        return

    m = SUBMIT_RE.search(line)
    if m:
        proc = m.group(1)
        state.current_task = proc
        state.progress_label = proc
        set_phase_for_task(state, proc)
        return

    m = DURATION_RE.search(line)
    if m:
        state.pipeline_duration = m.group(1).strip()
        return

    m = TASKS_RE.search(line)
    if m:
        state.pipeline_tasks = m.group(1)
        return

    m = ERROR_RE.search(line)
    if m:
        state.last_error = m.group(1).strip()
        return

    if "Pipeline completed successfully" in line:
        state.success = True
        state.last_info = "Pipeline completed successfully"
        state.add_event(state.last_info)
        return

    if line.startswith("[CLUSTER]") or line.startswith("[PLOT]") or line.startswith("[SEARCH]") or line.startswith("[TREE]"):
        state.last_info = line
        state.add_event(line)
        return


def compact(text: str, width: int) -> str:
    if len(text) <= width:
        return text
    if width <= 3:
        return text[:width]
    return text[: width - 3] + "..."


def render_line(state: LiveState, start_time: float, spinner_idx: int) -> str:
    elapsed = int(time.time() - start_time)
    mm, ss = divmod(elapsed, 60)
    hh, mm = divmod(mm, 60)
    tstr = f"{hh:02d}:{mm:02d}:{ss:02d}"
    spin = SPINNER[spinner_idx % len(SPINNER)]

    core = f"{spin} {tstr}  {state.phase}  |  {state.current_task}"
    if state.progress_label:
        core += f" [{state.progress_label}"
        if state.progress_counts:
            core += f" {state.progress_counts}"
        core += "]"
    if state.task_detail:
        core += f"  -  {state.task_detail}"
    elif state.last_info:
        core += f"  -  {state.last_info}"
    return core


def print_final(state: LiveState, exit_code: int, log_path: Path) -> None:
    if exit_code == 0 and state.success:
        print("Pipeline finished successfully.")
    elif exit_code == 0:
        print("Pipeline finished (no explicit success banner detected).")
    else:
        print("Pipeline failed.")

    if state.pipeline_duration:
        print(f"Duration: {state.pipeline_duration}")
    if state.pipeline_tasks:
        print(f"Tasks Completed: {state.pipeline_tasks}")
    if state.last_error:
        print(f"Last Error: {state.last_error}")
    print(f"Raw log: {log_path}")

    if state.events:
        print("Recent events:")
        for evt in state.events[-5:]:
            print(f"- {evt}")


def resolve_nextflow_bin(project_root: Path) -> str:
    """Resolve nextflow executable from PATH or project-local wrapper."""
    nf_in_path = which("nextflow")
    if nf_in_path:
        return nf_in_path

    local_nf = project_root / "nextflow"
    if local_nf.exists() and local_nf.is_file():
        return str(local_nf)

    raise FileNotFoundError(
        "Could not find `nextflow` on PATH or `./nextflow` in project root."
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run SynTerra with controlled live logging",
        epilog=(
            "All unknown arguments are forwarded to `nextflow run main.nf`.\n"
            "Example: ./synterra --gene P01501 --mode easy --outdir results -resume"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Bypass live UI and print raw nextflow output",
    )
    ns, forwarded = parser.parse_known_args()

    project_root = Path(__file__).resolve().parents[1]
    main_nf = project_root / "main.nf"
    if not main_nf.exists():
        print(f"Cannot find pipeline entrypoint: {main_nf}", file=sys.stderr)
        return 2

    try:
        nextflow_bin = resolve_nextflow_bin(project_root)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    # Ensure project script is used and disable Nextflow ANSI task table
    cmd = [nextflow_bin, "run", str(main_nf), "-ansi-log", "false", *forwarded]

    logs_dir = Path(".synterra_logs")
    logs_dir.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    raw_log = logs_dir / f"run_{stamp}.log"

    env = os.environ.copy()
    env.setdefault("NXF_DISABLE_CHECK_LATEST", "true")

    proc = subprocess.Popen(
        cmd,
        cwd=str(project_root),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )

    if proc.stdout is None:
        print("Failed to capture Nextflow output stream.", file=sys.stderr)
        return 2

    q: queue.Queue[str | None] = queue.Queue()

    def reader() -> None:
        for raw in proc.stdout:
            q.put(raw.rstrip("\n"))
        q.put(None)

    t = threading.Thread(target=reader, daemon=True)
    t.start()

    state = LiveState()
    start_time = time.time()
    spinner_idx = 0
    done = False
    interactive_ui = (not ns.raw) and sys.stdout.isatty()
    last_non_tty_line = ""
    last_non_tty_emit = 0.0

    with raw_log.open("w", encoding="utf-8") as raw_out:
        while True:
            drained = False
            while True:
                try:
                    item = q.get_nowait()
                except queue.Empty:
                    break

                drained = True
                if item is None:
                    done = True
                    break

                raw_out.write(item + "\n")
                raw_out.flush()
                clean = strip_ansi(item).strip()

                if ns.raw:
                    print(clean)
                    continue

                update_state(state, clean)
                if should_skip_line(clean):
                    continue

            if interactive_ui:
                width = get_terminal_size((120, 30)).columns
                line = compact(render_line(state, start_time, spinner_idx), max(20, width - 1))
                sys.stdout.write("\r" + line.ljust(width - 1))
                sys.stdout.flush()
                spinner_idx += 1
            elif not ns.raw:
                # Non-interactive output: print occasional snapshots, not a spammy spinner stream.
                line = render_line(state, start_time, spinner_idx)
                now = time.time()
                if line != last_non_tty_line and (now - last_non_tty_emit) >= 1.0:
                    print(line)
                    last_non_tty_line = line
                    last_non_tty_emit = now
                spinner_idx += 1

            if done and proc.poll() is not None:
                break

            # Keep spinner moving even when no new lines arrive.
            time.sleep(0.08 if not drained else 0.02)

    exit_code = proc.wait()
    if interactive_ui:
        width = get_terminal_size((120, 30)).columns
        sys.stdout.write("\r" + " " * (width - 1) + "\r")
        sys.stdout.flush()
    print_final(state, exit_code, raw_log)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
