"""The environment a search runs in. A search that picks up the wrong Python fails in
a way that looks like a problem with the gene, so the launcher makes sure it cannot,
the start-up check looks at what a task would really get, and a failure of this kind
stops the queue instead of failing one waiting search after the other.

The tests use a make-believe conda (a few shell scripts), so they need no real one.
"""

import os
import stat
import subprocess
import time
from pathlib import Path

import pytest

from conftest import run_until
from synvoy_web import cli, db, jobs
from synvoy_web.runner import TASK_PROBE, WRAPPER


def script(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\n" + text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def fake(tmp_path):
    """A conda base with one environment, a virtualenv, and a PATH in which the
    virtualenv comes first. ``conda activate`` behaves like the real one when another
    environment is active: it puts the new environment where the old one stood, that
    is BEHIND whatever was in front."""
    base = tmp_path / "conda"
    env = base / "envs" / "synvoy_env"
    script(env / "bin" / "python3", 'echo "conda-python $*"\n')
    script(base / "bin" / "python3", 'echo "base-python"\n')
    (base / "etc" / "profile.d").mkdir(parents=True)
    (base / "etc" / "profile.d" / "conda.sh").write_text(
        'conda() {\n'
        '  if [ "$1" = "activate" ]; then\n'
        '    export CONDA_PREFIX="%s/envs/$2"\n'
        '    [ -d "$CONDA_PREFIX" ] || return 1\n'
        '    export PATH="${PATH//%s\\/bin/$CONDA_PREFIX/bin}"\n'      # in place of the base
        '  fi\n'
        '}\n' % (base, str(base).replace("/", "\\/")))
    venv = tmp_path / "project" / ".venv"
    (venv / "pyvenv.cfg").parent.mkdir(parents=True)
    (venv / "pyvenv.cfg").write_text("home = /usr/bin\n")
    script(venv / "bin" / "python3", 'echo "venv-python"\n')
    other = tmp_path / "plain"
    script(other / "bin" / "tool", 'echo tool\n')
    state = tmp_path / "state"
    state.mkdir()
    (base / "condabin").mkdir()                              # where the conda command lives
    path = os.pathsep.join([str(venv / "bin"), str(other / "bin"), str(base / "bin"),
                            str(base / "condabin"), "/usr/bin", "/bin"])
    return {"base": base, "env": env, "venv": venv, "state": state, "path": path, "tmp": tmp_path}


def launch(fake, command, virtual_env=None, conda_env="synvoy_env"):
    env = {"PATH": fake["path"], "SYNVOY_WEB_CONDA_ENV": conda_env,
           "SYNVOY_WEB_CONDA_BASE": str(fake["base"]), "HOME": str(fake["tmp"])}
    if virtual_env:
        env["VIRTUAL_ENV"] = virtual_env
    done = subprocess.run(["bash", str(WRAPPER), str(fake["state"])] + command,
                          capture_output=True, text=True, env=env, timeout=30)
    code_file = fake["state"] / "exit_code"
    return done, (code_file.read_text().strip() if code_file.exists() else None)


@pytest.mark.parametrize("named", [False, True])
def test_search_gets_the_conda_python_whatever_is_in_front(fake, named):
    """The run that failed for real: a virtualenv first in PATH, VIRTUAL_ENV not set."""
    done, code = launch(fake, ["bash", "-c", 'python3; echo "PATH=$PATH"; echo "VENV=[${VIRTUAL_ENV:-}]"'],
                        virtual_env=str(fake["venv"]) if named else None)
    assert done.stdout.splitlines()[0].strip() == "conda-python", done.stdout + done.stderr
    path = [l for l in done.stdout.splitlines() if l.startswith("PATH=")][0][5:].split(os.pathsep)
    assert path[0] == str(fake["env"] / "bin")              # the environment comes first
    assert str(fake["venv"] / "bin") not in path            # the virtualenv is gone
    assert str(fake["tmp"] / "plain" / "bin") in path       # ordinary folders stay
    assert "VENV=[]" in done.stdout and code == "0"


def test_virtualenv_named_only_by_the_variable_is_dropped_too(fake):
    bare = fake["tmp"] / "bare_venv"                         # no pyvenv.cfg in it
    script(bare / "bin" / "python3", 'echo "venv-python"\n')
    fake["path"] = str(bare / "bin") + os.pathsep + fake["path"]
    done, _ = launch(fake, ["bash", "-c", 'python3; echo "PATH=$PATH"'], virtual_env=str(bare))
    lines = done.stdout.splitlines()
    assert lines[0].strip() == "conda-python"
    assert str(bare / "bin") not in lines[1][5:].split(os.pathsep)


def test_exit_code_of_the_search_is_recorded(fake):
    done, code = launch(fake, ["bash", "-c", "exit 7"])
    assert done.returncode == 7 and code == "7"
    assert not (fake["state"] / "exit_code.tmp").exists()
    done, code = launch(fake, ["bash", "-c", "echo fine"], conda_env="")     # no conda wanted
    assert (done.returncode, code, done.stdout.strip()) == (0, "0", "fine")


def test_environment_that_cannot_be_activated_is_code_96(fake):
    done, code = launch(fake, ["python3"], conda_env="no_such_env")
    assert done.returncode == 96 and code == "96" and "cannot activate" in done.stdout
    fake["base"] = fake["tmp"] / "nowhere"
    done, code = launch(fake, ["python3"])
    assert done.returncode == 96 and code == "96" and "conda not found" in done.stdout


def make_task_env(fake, python_text):
    """What Nextflow builds for the tasks, plus the two conda entry points a task uses."""
    task_env = fake["tmp"] / "cache" / "env-abc123"
    script(task_env / "bin" / "python3", python_text)
    script(fake["base"] / "condabin" / "conda",
           'echo "{"; echo "  \\"conda_prefix\\": \\"%s\\","; echo "}"\n' % fake["base"])
    # `source <base>/bin/activate <env>`: the new environment replaces the active one in place.
    (fake["base"] / "bin" / "activate").write_text(
        'export PATH="${PATH//${CONDA_PREFIX//\\//\\\\/}\\/bin/$1/bin}"\nexport CONDA_PREFIX="$1"\n')
    return task_env


def test_probe_sees_what_a_task_would_run(fake):
    good = 'if [ "$1" = "-c" ]; then echo "%s"; else echo task-python; fi\n'
    task_env = make_task_env(fake, "")
    script(task_env / "bin" / "python3", good % (task_env / "bin" / "python3"))
    done, _ = launch(fake, ["bash", str(TASK_PROBE), str(task_env)])
    assert done.returncode == 0 and done.stdout.strip() == str(task_env / "bin" / "python3")

    # The same probe WITHOUT the launcher's clean-up: the virtualenv stays in front,
    # and the probe says so. This is the state that failed two real searches.
    env = {"PATH": fake["path"].replace(str(fake["base"] / "bin"), str(fake["env"] / "bin")
                                        + os.pathsep + str(fake["base"] / "bin")),
           "CONDA_PREFIX": str(fake["env"]), "HOME": str(fake["tmp"])}
    raw = subprocess.run(["bash", str(TASK_PROBE), str(task_env)], capture_output=True,
                         text=True, env=env, timeout=30)
    assert raw.stdout.strip() == "venv-python"               # not the task's Python

    script(task_env / "bin" / "python3", 'echo "No module named parasail" >&2; exit 1\n')
    done, _ = launch(fake, ["bash", str(TASK_PROBE), str(task_env)])
    assert done.returncode == 1


# ── what the worker does with it ────────────────────────────────────────────

def test_start_up_check_of_the_worker_looks_at_the_tasks_python(make_env, fake, monkeypatch):
    """The check the worker runs before it takes searches, on its real path: Nextflow,
    conda, and the Python a task would get."""
    script(fake["env"] / "bin" / "nextflow", 'echo "  version 26.04.0 build 1"\n')
    monkeypatch.setenv("PATH", fake["path"])
    cache = fake["tmp"] / "cache"
    checkout = fake["tmp"] / "checkout"                      # all the check wants to see of it
    checkout.mkdir()
    (checkout / "main.nf").write_text("// stand-in\n")
    env = make_env(pipeline_command=[], conda_env="synvoy_env", conda_base=str(fake["base"]),
                   conda_cache_dir=str(cache), repo_dir=str(checkout))
    worker = env.worker()
    ready = worker.preflight()                               # nothing built yet: no objection
    assert ready["ok"] and ready["nextflow"] == "26.04.0", ready

    task_env = make_task_env(fake, "")
    python = task_env / "bin" / "python3"
    script(python, 'if [ "$1" = "-c" ]; then echo "%s"; fi\n' % python)
    assert worker.preflight()["ok"]

    script(python, 'echo "No module named parasail" >&2; exit 1\n')
    ready = worker.preflight()
    assert not ready["ok"] and "parasail" in ready["problems"][0]
    assert jobs.service_state(env.ctx)[0] is False           # and no search is taken
    assert any(line.startswith("NOT READY") for line in env.logs)

    script(python, 'echo "/somewhere/else/bin/python3"\n')   # answers, but is another Python
    ready = worker.preflight()
    assert not ready["ok"] and "instead of the Python" in ready["problems"][0]

    script(python, 'if [ "$1" = "-c" ]; then echo "%s"; fi\n' % python)
    assert worker.preflight()["ok"]                          # repaired: ready again


def submit(env, name, address="10.0.0.1"):
    env.beat()
    return jobs.submit(env.ctx, env.payload(name=name), address)


def test_broken_software_pauses_the_queue(env):
    worker = env.worker()
    first = submit(env, "noparasail")
    run_until(worker, lambda: env.store.get_job(first["id"])["state"] == db.FAILED)
    job = env.store.get_job(first["id"])
    assert job["error_code"] == "server_setup" and "operator" in job["error_text"]
    paused = env.store.get_meta("paused")
    assert paused["code"] == "server_setup" and paused["until"] > time.time()
    assert jobs.service_state(env.ctx)[:2] == (False, "paused")
    assert any("PAUSED until" in line and "resume" in line for line in env.logs)

    env.store.set_meta("paused", None)                       # as if repaired, to queue one
    second = submit(env, "ok_after")
    worker.pause_after_setup_failure(first["id"])
    for _ in range(5):
        worker.tick()
    assert env.store.get_job(second["id"])["state"] == db.QUEUED       # it waits, it does not fail
    assert env.store.get_job(second["id"])["attempt"] == 0

    paused = env.store.get_meta("paused")                    # the pause ends by itself
    env.store.set_meta("paused", dict(paused, until=int(time.time()) - 1))
    run_until(worker, lambda: env.store.get_job(second["id"])["state"] == db.DONE)
    assert env.store.get_meta("paused") is None
    assert jobs.service_state(env.ctx)[0] is True


def test_operator_can_pause_and_resume(env, tmp_path, capsys):
    import json
    config = tmp_path / "settings.json"
    config.write_text(json.dumps({k: v for k, v in env.settings.as_dict().items()
                                  if k in ("data_dir", "site_dir", "repo_dir", "conda_env",
                                           "pipeline_command", "min_free_gb", "disk_factor")}))
    worker = env.worker()
    assert cli.main(["--config", str(config), "pause", "Back", "at", "noon."]) == 0
    env.beat()
    accepting, code, message = jobs.service_state(env.ctx)
    assert (accepting, code, message) == (False, "paused", "Back at noon.")
    with pytest.raises(Exception) as caught:
        jobs.submit(env.ctx, env.payload(), "10.0.0.1")
    assert caught.value.status == 503 and caught.value.message == "Back at noon."

    env.store.set_meta("paused", None)
    waiting = submit(env, "ok_waiting")
    cli.main(["--config", str(config), "pause"])
    for _ in range(5):
        worker.tick()                                        # a full disk check must not reopen it
    assert env.store.get_job(waiting["id"])["state"] == db.QUEUED
    assert env.store.get_meta("paused")["code"] == "operator"
    assert cli.main(["--config", str(config), "resume"]) == 0
    run_until(worker, lambda: env.store.get_job(waiting["id"])["state"] == db.DONE)
    assert cli.main(["--config", str(config), "resume"]) == 0
    assert "was not paused" in capsys.readouterr().out
