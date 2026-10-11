"""Finding and stopping the processes of one search.

Every process started for a search inherits ``SYNVOY_WEB_JOB=<id>`` in its
environment. Looking for that marker under /proc finds all of them: the launcher,
Nextflow, each task and whatever a task started, also after they changed process
group or lost their parent because the worker was restarted. A signal to the process
group alone would miss those.
"""

import os
import signal
import time

MARKER = "SYNVOY_WEB_JOB"


def job_pids(job_id, proc_root="/proc"):
    """Process ids that belong to the search. Empty when none is left (a zombie has no
    environment any more and is not counted). None when /proc cannot be read at all."""
    needle = ("%s=%s" % (MARKER, job_id)).encode("utf-8")
    try:
        names = os.listdir(proc_root)
    except OSError:
        return None
    own = os.getpid()
    found = []
    for name in names:
        if not name.isdigit() or int(name) == own:
            continue
        try:
            with open(os.path.join(proc_root, name, "environ"), "rb") as fh:
                data = fh.read()
        except OSError:
            continue        # gone in the meantime, or another user's process
        if needle in data.split(b"\0"):
            found.append(int(name))
    return found


def _signal(pids, sig):
    for pid in pids:
        try:
            os.kill(pid, sig)
        except (ProcessLookupError, PermissionError):
            pass


def is_alive(job_id, handle=None):
    """True while any process of the search exists. ``handle`` is the worker's own
    ``Popen`` object when it started the search itself (cheaper than reading /proc)."""
    if handle is not None and handle.poll() is None:
        return True
    pids = job_pids(job_id)
    if pids is None:                      # no /proc: all we have is the handle
        return False
    return bool(pids)


def terminate(job_id, grace, handle=None, sleep=time.sleep, kill_wait=10.0):
    """Stop every process of the search: ask (TERM), wait up to ``grace`` seconds, then
    force (KILL). Returns the number of processes still there afterwards (0 = clean)."""
    pids = job_pids(job_id)
    if pids is None:
        # No /proc: fall back to the process group of the launcher we started.
        if handle is None or handle.poll() is not None:
            return 0
        for sig, wait in ((signal.SIGTERM, grace), (signal.SIGKILL, kill_wait)):
            try:
                os.killpg(handle.pid, sig)
            except (ProcessLookupError, PermissionError):
                return 0
            deadline = time.monotonic() + wait
            while time.monotonic() < deadline:
                if handle.poll() is not None:
                    return 0
                sleep(0.2)
        return 1
    if not pids:
        return 0
    _signal(pids, signal.SIGTERM)
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        if handle is not None:
            handle.poll()                 # reap our child so it does not linger as a zombie
        pids = job_pids(job_id) or []
        if not pids:
            return 0
        sleep(0.2)
    # Whatever is left did not react; new processes may also have appeared meanwhile.
    deadline = time.monotonic() + kill_wait
    while True:
        pids = job_pids(job_id) or []
        if not pids:
            return 0
        _signal(pids, signal.SIGKILL)
        if handle is not None:
            handle.poll()
        if time.monotonic() >= deadline:
            return len(pids)
        sleep(0.2)
