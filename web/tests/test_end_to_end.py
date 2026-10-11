"""The worker with the REAL pipeline. Not part of the normal test run: it needs the
pipeline's software, real genomes and several minutes.

    SYNVOY_WEB_E2E_CONFIG=var/settings.local.json \\
    SYNVOY_WEB_E2E_QUERY=../local_data/demo/query/PGK1.faa \\
    SYNVOY_WEB_E2E_HOME=saccharomyces_cerevisiae \\
        python -m pytest tests/test_end_to_end.py -s

``CONFIG`` is a settings file of a working installation (its genome library and conda
set-up are used; searches and database go to a temporary folder). ``QUERY`` is a
protein FASTA file, ``HOME`` the library id of the genome it comes from. The three
closest genomes of the library are searched.

What it shows, on real processes (Nextflow, its tasks, the tools they start):
  1. cancelling a running search stops every process and leaves no working files;
  2. a search whose processes are all killed (a crash, a reboot) is resumed by the
     worker and comes to a result, reusing the steps that had finished.
"""

import json
import os
import signal
import time
from pathlib import Path

import pytest

from synvoy_web import db, jobs, procs
from synvoy_web.config import load_settings
from synvoy_web.db import Store
from synvoy_web.jobs import Context
from synvoy_web.library import LibraryWatcher
from synvoy_web.runner import JobDir, Worker

CONFIG = os.environ.get("SYNVOY_WEB_E2E_CONFIG")
QUERY = os.environ.get("SYNVOY_WEB_E2E_QUERY")
HOME = os.environ.get("SYNVOY_WEB_E2E_HOME")

pytestmark = pytest.mark.skipif(
    not (CONFIG and QUERY and HOME),
    reason="set SYNVOY_WEB_E2E_CONFIG, _QUERY and _HOME to run the real pipeline")


@pytest.fixture(scope="module")
def installation(tmp_path_factory):
    base = load_settings(CONFIG, env={})
    values = base.as_dict()
    values.update({
        "data_dir": str(tmp_path_factory.mktemp("e2e")),
        "library_file": str(base.library_path),
        "conda_cache_dir": str(base.conda_cache_path),      # the built environment is reused
        "poll_seconds": 1.0, "max_running": 1, "max_attempts": 2,
        "max_active_per_client": 5, "max_submissions_per_day": 50,
    })
    settings = load_settings(env={}, overrides=values)
    store = Store(settings.db_path)
    ctx = Context(settings, store, LibraryWatcher(settings.library_path))
    logs = []
    worker = Worker(ctx, log=lambda message: (logs.append(message), print("   worker:", message)))
    worker.lock()
    ready = worker.preflight()
    assert ready["ok"], ready["problems"]
    library = ctx.library()
    home = library.get(HOME)
    assert home is not None and home.can_be_home, "%s is not a starting genome of the library" % HOME
    targets = [g.id for g in library.targets_for(home)
               if g.size_mb <= settings.max_genome_mb][:3]
    assert targets, "the library has no genome to search"
    yield {"ctx": ctx, "worker": worker, "settings": settings, "store": store, "logs": logs,
           "targets": targets, "sequence": Path(QUERY).read_text()}
    for job in store.jobs_in_state(db.RUNNING, db.QUEUED):
        procs.terminate(job["id"], 20)
    store.close()


def submit(inst, name):
    payload = {"query": {"sequence": inst["sequence"], "name": name, "home": HOME},
               "targets": inst["targets"]}
    inst["store"].set_meta("worker_heartbeat", int(time.time()))
    return jobs.submit(inst["ctx"], payload, "127.0.0.1")


def tick_until(inst, predicate, timeout, what):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        inst["worker"].tick()
        if predicate():
            return
        time.sleep(1.0)
    raise AssertionError("timed out after %d s waiting for %s" % (timeout, what))


def wait_until_searching(inst, job_id, seconds):
    """Tick until the genome search of this job has been at work for ``seconds``: the
    step is active and the search has processes of its own (launcher, Nextflow, the
    task and the tools it starts). Looking earlier would catch a task that has only
    just been handed over."""
    tick_until(inst, lambda: in_search_step(inst, job_id), 900, "the genome search to start")
    since = time.monotonic()
    tick_until(inst, lambda: in_search_step(inst, job_id) and time.monotonic() - since > seconds
               and len(procs.job_pids(job_id)) >= 3, 300, "the genome search to be at work")


def in_search_step(inst, job_id):
    job = inst["store"].get_job(job_id)
    if job["state"] in db.FINAL_STATES:            # do not wait for what cannot come
        tail = JobDir(inst["settings"].jobs_dir, job_id).console
        raise AssertionError("the search ended early: %s %s %s\n%s" % (
            job["state"], job["error_code"], job["error_text"],
            tail.read_text(errors="replace")[-1500:] if tail.is_file() else ""))
    steps = {s["key"]: s["state"] for s in (job["progress"] or {}).get("steps", [])}
    return job["state"] == db.RUNNING and steps.get("search") == "active"


def test_cancel_stops_the_real_pipeline(installation):
    inst = installation
    outcome = submit(inst, "e2e_cancel")
    wait_until_searching(inst, outcome["id"], 15)
    running = procs.job_pids(outcome["id"])
    print("   %d processes belong to the running search" % len(running))
    assert len(running) >= 3                       # launcher, Nextflow, at least one task
    jd = JobDir(inst["settings"].jobs_dir, outcome["id"])
    assert (jd.run / "work").is_dir()

    assert jobs.cancel(inst["ctx"], outcome["id"], outcome["token"])["cancelled"] is True
    started = time.monotonic()
    tick_until(inst, lambda: inst["store"].get_job(outcome["id"])["state"] == db.CANCELLED, 120,
               "the search to be cancelled")
    print("   cancelled %.0f s after the request" % (time.monotonic() - started))
    time.sleep(2)
    assert procs.job_pids(outcome["id"]) == []     # nothing survived
    for pid in running:
        assert not os.path.exists("/proc/%d" % pid) or \
            b"SYNVOY_WEB_JOB" not in open("/proc/%d/environ" % pid, "rb").read()
    assert not jd.run.exists() and not jd.out.exists() and not jd.input.exists()
    assert jd.console.is_file()


def test_killed_search_is_resumed_to_a_result(installation):
    inst = installation
    outcome = submit(inst, "e2e_resume")
    job_id = outcome["id"]
    wait_until_searching(inst, job_id, 25)
    before = inst["store"].get_job(job_id)["progress"]["tasks_done"]
    victims = procs.job_pids(job_id)
    print("   killing %d processes (%d tasks had finished)" % (len(victims), before))
    assert len(victims) >= 3                       # launcher, Nextflow, the search and its tools
    for pid in victims:                            # no warning, as a crash or a reboot would
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    time.sleep(2)
    for pid in procs.job_pids(job_id):             # children started in that instant
        os.kill(pid, signal.SIGKILL)

    tick_until(inst, lambda: inst["store"].get_job(job_id)["attempt"] == 2, 120,
               "the worker to start the search again")
    assert any("interrupted" in line for line in inst["logs"])
    tick_until(inst, lambda: inst["store"].get_job(job_id)["state"] in db.FINAL_STATES, 3600,
               "the resumed search to end")
    job = inst["store"].get_job(job_id)
    assert job["state"] == db.DONE, (job["error_code"], job["error_text"])
    assert job["attempt"] == 2

    jd = JobDir(inst["settings"].jobs_dir, job_id)
    console = jd.console.read_text(errors="replace")
    if any("from the beginning" in line for line in inst["logs"]):
        # Someone changed the pipeline's code between the two attempts. The worker then
        # reuses nothing, which is right, but it is not what this test is about.
        assert "-resume" not in console
        pytest.skip("the checkout changed between the attempts, so the search was started "
                    "from the beginning; run this test while nobody edits the pipeline")
    assert console.count("===== attempt ") == 2 and "-resume" in console
    cached = sum(1 for line in jd.trace(2).read_text().splitlines() if "\tCACHED\t" in line)
    print("   %d finished steps were reused by the second attempt" % cached)
    # The work of the first attempt was not lost. (Two steps of slack for what ended in
    # the instant of the kill. With a query file that differed between the attempts
    # this was 5 of 13: every step that reads the query ran again.)
    assert cached >= max(1, before - 2), "%d steps had finished, %d were reused" % (before, cached)

    result = json.loads((jd.result / "result.json").read_text())
    print("   result:", {k: result["summary"][k] for k in ("genomes", "found", "candidate", "none")})
    assert result["summary"]["genomes"] == len(inst["targets"])
    assert result["self_check"]["agrees"] is True
    assert result["loci"] and result["loci"][0]["files"].get("grid")
    assert procs.job_pids(job_id) == []
    assert not jd.run.exists()
