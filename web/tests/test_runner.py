"""The worker under the conditions that matter for a public service: success, every
kind of failure, cancelling, the time limit, a restart in the middle of a search, a
full disk, and cleaning up afterwards."""

import json
import os
import re
import time
from pathlib import Path

import pytest

from conftest import PROTEIN, run_until
from synvoy_web import db, jobs, procs
from synvoy_web.query import QueryError
from synvoy_web.runner import JobDir, Worker, WorkerError, remove_tree


def submit(env, name="ok", address="10.0.0.1", **kwargs):
    env.beat()
    return jobs.submit(env.ctx, env.payload(name=name, **kwargs), address)


def state(env, job_id):
    return env.store.get_job(job_id)["state"]


def wait_for(env, worker, job_id, wanted, timeout=30.0):
    run_until(worker, lambda: state(env, job_id) in wanted, timeout)
    return env.store.get_job(job_id)


def folder(env, job_id):
    return JobDir(env.settings.jobs_dir, job_id)


# ── the normal course ───────────────────────────────────────────────────────

def test_search_runs_to_a_result(env):
    worker = env.worker()
    job_id = submit(env)["id"]
    job = wait_for(env, worker, job_id, (db.DONE, db.FAILED))
    assert job["state"] == db.DONE, job
    assert job["attempt"] == 1 and job["error_code"] is None
    assert job["expires_at"] > job["finished_at"]
    assert all(step["state"] == "done" for step in job["progress"]["steps"])

    jd = folder(env, job_id)
    result = json.loads((jd.result / "result.json").read_text())
    assert result["summary"]["genomes"] == 3
    assert result["summary"]["found"] == 2 and result["summary"]["candidate"] == 1
    assert result["summary"]["failed_tasks"] == 1          # the ignored task of the trace
    assert result["self_check"]["agrees"] is True
    assert [g["species"] for g in result["genomes"]] == ["Alpha beta", "Epsilon zeta", "Eta theta"]
    assert (jd.result / "files" / "synvoy_report.json").is_file()

    # Working files are gone, the logs are kept, the library is untouched.
    assert not jd.run.exists() and not jd.out.exists() and not jd.input.exists()
    assert jd.console.is_file() and jd.trace(1).is_file()
    for name in ("Homus_primus.fna", "Homus_primus.gff", "Alpha_beta.fna"):
        assert (env.genomes_dir / name).is_file()
    assert procs.job_pids(job_id) == []


def test_command_line_and_inputs(env):
    worker = env.worker()
    job = env.store.claim(submit(env)["id"], int(time.time()))
    jd = folder(env, job["id"])
    inputs = worker.prepare_inputs(jd, job)
    argv = worker.command(jd, job, inputs, resume=False)
    assert argv[:2] == env.settings.pipeline_command
    for flag, value in (("--mode", "pro"), ("--home_species", "Homus primus"),
                        ("--target_genomes", str(jd.input / "targets")),
                        ("--outdir", str(jd.out)), ("-work-dir", str(jd.work)),
                        ("-profile", "auto,low_mem")):
        assert argv[argv.index(flag) + 1] == value
    assert "-resume" not in argv and "--target_gffs" not in argv
    assert "-resume" in worker.command(jd, job, inputs, resume=True)
    # File names are the species: the pipeline uses them as genome identifiers.
    assert sorted(p.name for p in (jd.input / "targets").iterdir()) == [
        "Alpha_beta.fna", "Epsilon_zeta.fna", "Eta_theta.fna"]
    assert (jd.input / "home" / "Homus_primus.gff").is_symlink()
    assert (jd.input / "query.faa").read_text().startswith(">ok\n")


def test_target_annotations_only_when_switched_on(make_env):
    env = make_env(use_target_annotations=True)
    worker = env.worker()
    job = env.store.claim(submit(env, targets=["gamma_delta", "alpha_beta"])["id"], int(time.time()))
    jd = folder(env, job["id"])
    argv = worker.command(jd, job, worker.prepare_inputs(jd, job), resume=False)
    assert argv[argv.index("--target_gffs") + 1] == str(jd.input / "target_gffs")
    assert [p.name for p in (jd.input / "target_gffs").iterdir()] == ["Gamma_delta.gff"]


def test_service_settings_are_not_passed_to_the_pipeline(env, monkeypatch):
    monkeypatch.setenv("SYNVOY_WEB_SMTP_PASSWORD", "secret")
    monkeypatch.setenv("SYNVOY_WEB_CONFIG", "/somewhere/settings.json")
    monkeypatch.setenv("PYTHONPATH", "/somewhere")
    child = env.worker().child_env("abc")
    assert "SYNVOY_WEB_SMTP_PASSWORD" not in child and "SYNVOY_WEB_CONFIG" not in child
    assert "PYTHONPATH" not in child
    assert child[procs.MARKER] == "abc"
    assert child["NXF_CONDA_CACHEDIR"] == str(env.settings.conda_cache_path)


# ── failures ────────────────────────────────────────────────────────────────

def test_failed_task_is_reported_with_its_step(env):
    worker = env.worker()
    job = wait_for(env, worker, submit(env, "fail")["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.FAILED and job["error_code"] == "pipeline_failed"
    assert "searching the genomes" in job["error_text"]
    assert not folder(env, job["id"]).run.exists()
    assert folder(env, job["id"]).console.is_file()         # kept for the operator


def test_gene_not_in_home_genome_is_remembered(env):
    worker = env.worker()
    first = submit(env, "nogene")
    job = wait_for(env, worker, first["id"], (db.DONE, db.FAILED))
    assert job["error_code"] == "gene_not_in_home" and "Homus primus" in job["error_text"]
    again = submit(env, "nogene", address="10.9.9.9")
    assert again["id"] == first["id"] and again["reused"] is True   # no second run


def test_other_failures_are_run_again(env):
    worker = env.worker()
    first = submit(env, "fail")
    wait_for(env, worker, first["id"], (db.FAILED,))
    again = submit(env, "fail")
    assert again["id"] != first["id"] and again["reused"] is False


def test_success_without_report_is_a_failure(env):
    worker = env.worker()
    job = wait_for(env, worker, submit(env, "noreport")["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.FAILED and job["error_code"] == "no_result"


def test_environment_that_cannot_be_set_up(make_env, tmp_path):
    env = make_env(conda_env="does_not_exist", conda_base=str(tmp_path / "no_conda"))
    worker = env.worker()
    env.beat()
    job_id = jobs.submit(env.ctx, env.payload(), "10.0.0.1")["id"]
    worker._preflight_due = time.monotonic() + 3600        # let the search start anyway
    job = wait_for(env, worker, job_id, (db.DONE, db.FAILED))
    assert job["error_code"] == "server_setup"
    assert worker._preflight_due == 0                       # the environment is checked again


def test_genome_replaced_after_submission(env):
    worker = env.worker()
    job_id = submit(env)["id"]
    data = json.loads(env.library_file.read_text())
    for genome in data["genomes"]:
        if genome["id"] == "alpha_beta":
            genome["checksum"] = "sha256:another-assembly"
    env.library_file.write_text(json.dumps(data))
    job = wait_for(env, worker, job_id, (db.DONE, db.FAILED))
    assert job["error_code"] == "genome_changed" and "Alpha beta" in job["error_text"]


def test_code_updated_under_a_running_search_is_flagged(make_env, tmp_path, monkeypatch):
    repo = tmp_path / "checkout"
    (repo / "bin").mkdir(parents=True)
    (repo / "VERSION").write_text("9.9.9\n")
    (repo / "bin" / "tool.py").write_text("print(1)\n")
    env = make_env(repo_dir=str(repo))
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "2")
    worker = env.worker()
    clean = wait_for(env, worker, submit(env, "ok_clean")["id"], (db.DONE, db.FAILED))
    assert json.loads((folder(env, clean["id"]).result / "result.json").read_text())["notes"] == []

    outcome = running_search(env, worker, "slow_mixed")
    (repo / "bin" / "tool.py").write_text("print(2)\n")            # an update in mid-search
    job = wait_for(env, worker, outcome["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.DONE and job["cache_key"] == "mixed-code:" + job["id"]
    notes = json.loads((folder(env, job["id"]).result / "result.json").read_text())["notes"]
    assert [n["kind"] for n in notes] == ["code_changed"] and "updated while" in notes[0]["text"]
    assert any("code changed while" in line for line in env.logs)
    # The same question asked again is computed afresh, with the new code.
    env.ctx._version = None
    again = submit(env, "slow_mixed", address="10.0.0.7")
    assert again["id"] != job["id"] and again["reused"] is False
    jobs.cancel(env.ctx, again["id"], again["token"])
    wait_for(env, worker, again["id"], (db.CANCELLED, db.DONE))


def test_one_bad_search_does_not_stop_the_others(env, monkeypatch):
    worker = env.worker()
    bad = submit(env, "slow_bad")["id"]
    good = submit(env, "ok_good", address="10.0.0.2")["id"]
    real = Worker.read_progress

    def broken(self, jd):
        if jd.path.name == bad:
            raise RuntimeError("boom")
        return real(self, jd)

    monkeypatch.setattr(Worker, "read_progress", broken)
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "60")
    job = wait_for(env, worker, bad, (db.FAILED,))
    assert job["error_code"] == "internal_error"
    assert procs.job_pids(bad) == []                        # it was stopped, not abandoned
    assert wait_for(env, worker, good, (db.DONE, db.FAILED))["state"] == db.DONE


# ── stopping ────────────────────────────────────────────────────────────────

def running_search(env, worker, name="slow", **kwargs):
    outcome = submit(env, name, **kwargs)
    run_until(worker, lambda: state(env, outcome["id"]) == db.RUNNING
              and bool(procs.job_pids(outcome["id"])))
    return outcome


def test_cancel_stops_every_process(env, monkeypatch):
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "120")
    worker = env.worker()
    outcome = running_search(env, worker)
    assert len(procs.job_pids(outcome["id"])) >= 2          # launcher + pipeline
    assert jobs.cancel(env.ctx, outcome["id"], outcome["token"]) == {
        "cancelled": True, "others_waiting": 0}
    job = wait_for(env, worker, outcome["id"], (db.CANCELLED, db.DONE, db.FAILED))
    assert job["state"] == db.CANCELLED
    assert procs.job_pids(outcome["id"]) == []
    assert not folder(env, outcome["id"]).run.exists()


def test_cancel_needs_the_token_of_the_submitter(env, monkeypatch):
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "120")
    worker = env.worker()
    outcome = running_search(env, worker)
    for token in ("", None, "not-the-token"):
        with pytest.raises(QueryError) as caught:
            jobs.cancel(env.ctx, outcome["id"], token)
        assert caught.value.status == 403
    assert state(env, outcome["id"]) == db.RUNNING
    jobs.cancel(env.ctx, outcome["id"], outcome["token"])
    wait_for(env, worker, outcome["id"], (db.CANCELLED,))


def test_search_shared_by_two_visitors_stops_only_when_both_leave(env, monkeypatch):
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "120")
    worker = env.worker()
    first = running_search(env, worker)
    second = submit(env, "slow", address="10.0.0.2")
    assert second["id"] == first["id"] and second["reused"] is True
    assert jobs.cancel(env.ctx, first["id"], first["token"]) == {
        "cancelled": False, "others_waiting": 1}
    worker.tick()
    assert state(env, first["id"]) == db.RUNNING
    assert jobs.cancel(env.ctx, first["id"], second["token"])["cancelled"] is True
    wait_for(env, worker, first["id"], (db.CANCELLED,))


def test_cancel_while_waiting_never_starts_it(make_env, monkeypatch):
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "120")
    env = make_env(max_running=1)
    worker = env.worker()
    first = running_search(env, worker, "slow_1")
    second = submit(env, "slow_2")
    worker.tick()
    assert state(env, second["id"]) == db.QUEUED
    assert jobs.job_view(env.ctx, env.store.get_job(second["id"]))["queue"] == {
        "ahead": 0, "running": 1, "genomes": []}
    jobs.cancel(env.ctx, second["id"], second["token"])
    assert state(env, second["id"]) == db.CANCELLED
    jobs.cancel(env.ctx, first["id"], first["token"])
    wait_for(env, worker, first["id"], (db.CANCELLED,))
    worker.tick()
    assert state(env, second["id"]) == db.CANCELLED
    assert not folder(env, second["id"]).path.exists()      # it never got a folder


def test_time_limit_ends_a_search_that_ignores_requests_to_stop(make_env):
    env = make_env(max_runtime_hours=0.0003)                # about one second
    worker = env.worker()
    outcome = submit(env, "hang")
    job = wait_for(env, worker, outcome["id"], (db.FAILED, db.DONE), timeout=40)
    assert job["error_code"] == "timeout" and "time limit" in job["error_text"]
    assert procs.job_pids(outcome["id"]) == []              # TERM was ignored; KILL was not


def test_process_left_behind_by_a_finished_search_is_stopped(env):
    worker = env.worker()
    job = wait_for(env, worker, submit(env, "leak")["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.DONE
    assert procs.job_pids(job["id"]) == []


# ── restarts ────────────────────────────────────────────────────────────────

def test_restarted_worker_takes_over_a_running_search(env, monkeypatch):
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "3")
    first = env.worker()
    outcome = running_search(env, first)
    del first                                               # the worker process is gone
    second = env.worker()                                   # a new one knows nothing in memory
    assert second._handles == {}
    second.tick()
    assert state(env, outcome["id"]) == db.RUNNING
    job = wait_for(env, second, outcome["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.DONE and job["attempt"] == 1  # not started a second time


def test_search_that_died_with_the_machine_is_resumed(env):
    worker = env.worker()
    outcome = submit(env, "crashonce")
    job = wait_for(env, worker, outcome["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.DONE and job["attempt"] == 2
    console = folder(env, outcome["id"]).console.read_text()
    assert console.count("===== attempt ") == 2 and "-resume" in console
    assert any("interrupted" in line for line in env.logs)
    # Nextflow knows an input file by path, size and modification time. Both attempts
    # must have seen the very same query file, or the finished steps are run again.
    seen = re.findall(r"query file: (\d+ bytes, modified \d+)", console)
    assert len(seen) == 2 and seen[0] == seen[1], seen
    assert re.findall(r"working folder has seen: (\d+)", console) == ["1", "2"]   # it was kept


def test_search_started_again_after_a_code_update_begins_from_scratch(make_env, tmp_path):
    """Killed under the old code, started again under the new: nothing the first attempt
    computed may be reused, and the result is one of the new code only."""
    repo = tmp_path / "checkout"
    (repo / "bin").mkdir(parents=True)
    (repo / "VERSION").write_text("9.9.9\n")
    (repo / "bin" / "tool.py").write_text("print(1)\n")
    env = make_env(repo_dir=str(repo))
    worker = env.worker()
    outcome = submit(env, "crashfirst")
    first_key = env.store.get_job(outcome["id"])["cache_key"]
    run_until(worker, lambda: env.store.get_job(outcome["id"])["attempt"] == 1)    # launched once
    (repo / "bin" / "tool.py").write_text("print(2)\n")            # the update, before attempt 2
    job = wait_for(env, worker, outcome["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.DONE and job["attempt"] == 2
    console = folder(env, job["id"]).console.read_text()
    assert console.count("===== attempt ") == 2 and "-resume" not in console
    assert re.findall(r"working folder has seen: (\d+)", console) == ["1", "1"]   # it was emptied
    assert any("from the beginning" in line for line in env.logs)
    assert json.loads(job["version"])["fingerprint"] == jobs.code_fingerprint(repo)
    assert job["cache_key"] != first_key and not job["cache_key"].startswith("mixed-code")
    assert json.loads((folder(env, job["id"]).result / "result.json").read_text())["notes"] == []


def test_search_that_keeps_dying_fails_with_a_reason(env):
    worker = env.worker()
    job = wait_for(env, worker, submit(env, "crash")["id"], (db.DONE, db.FAILED))
    assert job["state"] == db.FAILED and job["error_code"] == "interrupted"
    assert job["attempt"] == env.settings.max_attempts


def test_only_one_worker_per_data_folder(env):
    first, second = env.worker(), env.worker()
    first.lock()
    with pytest.raises(WorkerError):
        second.lock()


def test_worker_that_is_not_ready_starts_nothing(make_env):
    env = make_env()
    env.library_file.write_text(json.dumps({"groups": [], "genomes": []}))
    worker = env.worker()
    assert worker.preflight()["ok"] is False
    accepting, code, _ = jobs.service_state(env.ctx)
    assert accepting is False and code in ("library_empty", "worker_not_ready")


# ── disk ────────────────────────────────────────────────────────────────────

def test_search_waits_when_the_disk_is_short(make_env):
    env = make_env(min_free_gb=5.0, disk_factor=4.0)
    free = {"bytes": 1 * 10 ** 9}
    worker = env.worker(disk_free=lambda path: free["bytes"])
    job_id = submit(env)["id"]
    for _ in range(3):
        worker.tick()
    assert state(env, job_id) == db.QUEUED
    assert env.store.get_meta("paused")["code"] == "low_disk"
    accepting, code, message = jobs.service_state(env.ctx)
    assert accepting is False and code == "paused" and "disk" in message
    free["bytes"] = 500 * 10 ** 9
    job = wait_for(env, worker, job_id, (db.DONE, db.FAILED))
    assert job["state"] == db.DONE
    assert not env.store.get_meta("paused")


def test_disk_need_grows_with_the_genomes(make_env):
    env = make_env(disk_factor=4.0)
    small = env.store.get_job(submit(env, "ok_a", targets=["epsilon_zeta"])["id"])
    large = env.store.get_job(submit(env, "ok_b", targets=["gamma_delta", "alpha_beta"])["id"])
    assert small["need_mb"] == 4 * (12 + 20)
    assert large["need_mb"] == 4 * (12 + 900 + 300)


# ── cleaning up ─────────────────────────────────────────────────────────────

def test_results_expire(env):
    worker = env.worker()
    job = wait_for(env, worker, submit(env)["id"], (db.DONE,))
    env.store.update_job(job["id"], expires_at=int(time.time()) - 1)
    worker.collect_garbage()
    assert state(env, job["id"]) == db.EXPIRED
    assert not folder(env, job["id"]).path.exists()
    # The sequence goes with the result; the name of the search stays.
    kept = env.store.get_job(job["id"])["request"]
    assert kept["query"]["fasta"] == "" and kept["query"]["label"] == "ok"
    assert PROTEIN not in "\n".join(env.store._conn().iterdump())
    # The same search can be started again and runs afresh.
    again = submit(env)
    assert again["id"] != job["id"] and again["reused"] is False


def test_nothing_links_a_finished_search_to_a_visitor(make_env):
    env = make_env(smtp_host="mail.invalid", smtp_from="synvoy@example.org",
                   public_base_url="https://synvoy.example.org")
    worker = env.worker()
    sent = []
    import synvoy_web.notify as notify
    original = notify.send
    notify.send = lambda settings, message, log=None: sent.append(message) or True
    try:
        env.beat()
        outcome = jobs.submit(env.ctx, env.payload(email="someone@example.org"), "203.0.113.9")
        rows = env.store._conn().execute("SELECT client, email FROM submitters").fetchall()
        assert rows[0]["email"] == "someone@example.org" and len(rows[0]["client"]) == 32
        wait_for(env, worker, outcome["id"], (db.DONE,))
    finally:
        notify.send = original
    assert len(sent) == 1 and sent[0]["To"] == "someone@example.org"
    assert "https://synvoy.example.org/result.html?id=" + outcome["id"] in sent[0].get_content()
    rows = env.store._conn().execute("SELECT client, email FROM submitters").fetchall()
    assert [(r["client"], r["email"]) for r in rows] == [("", None)]
    assert "someone@example.org" not in "\n".join(env.store._conn().iterdump())


def test_failed_and_cancelled_searches_lose_their_sequence_later(make_env):
    env = make_env(failed_retention_days=1)
    worker = env.worker()
    job = wait_for(env, worker, submit(env, "fail")["id"], (db.FAILED,))
    worker.collect_garbage()
    assert env.store.get_job(job["id"])["request"]["query"]["fasta"]      # still within a day
    env.store.update_job(job["id"], finished_at=int(time.time()) - 2 * 86400)
    worker.collect_garbage()
    assert env.store.get_job(job["id"])["request"]["query"]["fasta"] == ""
    assert not folder(env, job["id"]).path.exists()


def test_leftovers_of_a_crash_are_removed(env):
    worker = env.worker()
    job = wait_for(env, worker, submit(env)["id"], (db.DONE,))
    jd = folder(env, job["id"])
    (jd.run / "work").mkdir(parents=True)                   # as if clean-up had been cut short
    (jd.run / "work" / "big.tmp").write_text("x")
    orphan = env.settings.jobs_dir / "Zz0000000000000000000a"
    orphan.mkdir()
    old = time.time() - 7200
    os.utime(str(orphan), (old, old))
    fresh = env.settings.jobs_dir / "Zz0000000000000000000b"
    fresh.mkdir()
    foreign = env.settings.jobs_dir / "notes"
    foreign.mkdir()
    worker.collect_garbage()
    assert not jd.run.exists() and (jd.result / "result.json").is_file()
    assert not orphan.exists()
    assert fresh.exists()                                   # may belong to a search being made
    assert foreign.exists()                                 # not ours to delete


def test_remove_tree_never_leaves_its_root(tmp_path):
    root = tmp_path / "jobs"
    inside = root / "a" / "b"
    inside.mkdir(parents=True)
    outside = tmp_path / "precious"
    outside.mkdir()
    (outside / "genome.fna").write_text(">x\nACGT\n")
    with pytest.raises(RuntimeError):
        remove_tree(outside, root)
    with pytest.raises(RuntimeError):
        remove_tree(root, root)
    with pytest.raises(RuntimeError):
        remove_tree(root / "a" / ".." / ".." / "precious", root)
    link = root / "a" / "link"
    link.symlink_to(outside)
    remove_tree(link, root)                                 # removes the link only
    assert not link.exists() and (outside / "genome.fna").is_file()
    (inside / "linked").symlink_to(outside)
    locked = inside / "locked"
    locked.mkdir()
    (locked / "file").write_text("x")
    os.chmod(str(locked), 0o500)
    remove_tree(root / "a", root)
    assert not (root / "a").exists() and (outside / "genome.fna").is_file()


def test_search_folder_rejects_names_that_are_not_ids(env):
    for name in ("../x", "short", "a" * 22 + "/..", ""):
        with pytest.raises(ValueError):
            JobDir(env.settings.jobs_dir, name)


def test_progress_is_visible_while_running(env, monkeypatch):
    monkeypatch.setenv("FAKE_SLOW_SECONDS", "120")
    worker = env.worker()
    outcome = running_search(env, worker)

    def search_active():
        progress = env.store.get_job(outcome["id"])["progress"] or {}
        steps = {s["key"]: s for s in progress.get("steps", [])}
        return steps.get("search", {}).get("state") == "active" \
            and steps["search"]["detail"].get("genomes_total") == 3

    run_until(worker, search_active)
    steps = {s["key"]: s for s in env.store.get_job(outcome["id"])["progress"]["steps"]}
    assert steps["locate"]["state"] == "done" and steps["locate"]["detail"] == {"loci": 1}
    assert steps["search"]["detail"] == {"genomes_done": 1, "genomes_total": 3}
    assert steps["check"]["state"] == "waiting"
    jobs.cancel(env.ctx, outcome["id"], outcome["token"])
    wait_for(env, worker, outcome["id"], (db.CANCELLED,))
