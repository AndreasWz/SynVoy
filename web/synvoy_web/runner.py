"""The worker: takes queued searches, runs the pipeline for each, watches it and turns
its output into a result. One worker per data folder.

How it stays reliable
---------------------
* The worker keeps nothing in memory that it cannot rebuild. The state of a search
  is in the database and in its folder; whether it is still running is read from
  /proc (``procs.py``), its outcome from ``state/exit_code`` written by
  ``run_job.sh``. A restarted worker therefore carries on with searches that are
  still running, and resumes (``-resume``) those that died with the machine.
* A search is only started when the disk has room for it; otherwise it waits.
* Every search has a time limit. Stopping one (limit, cancel) ends every process it
  started, checked afterwards.
* One search that makes the worker raise does not stop the others: after three
  errors in a row it is marked failed and the queue moves on.
* Working folders are deleted when a search ends. Results expire.

Folder of a search: ``<data_dir>/jobs/<id>/``

    input/    the query and links to the library genomes (file name = species)
    run/      where Nextflow is started: its cache, its log, ``work/``
    out/      the pipeline's --outdir
    logs/     console output and the trace file of each attempt
    state/    pid, start time, exit code
    result/   the bundle shown on the result page (``results.py``)
"""

import fcntl
import json
import os
import re
import shlex
import shutil
import stat
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

from . import db, notify, procs, results
from . import progress as progress_mod
from .fsutil import remove_tree, write_atomic
from .jobs import JOB_ID_RE, cache_key, iso, job_version, pipeline_version

WRAPPER = Path(__file__).with_name("run_job.sh")
TASK_PROBE = Path(__file__).with_name("probe_task_env.sh")
ATTEMPT_MARK = "===== attempt "
SETUP_EXIT_CODE = 96
MAX_CONSOLE_BYTES = 5 * 1024 * 1024

# (pattern in the logs, error code, message for the visitor). First match wins.
_FAILURES = (
    (re.compile(r"Gene not found in home genome"), "gene_not_in_home",
     "The protein could not be located in the genome of {home}. If it comes from "
     "another species, the two may be too different for this search."),
    (re.compile(r"aa, below the|No sequences found in query|not amino-acid"), "query_rejected",
     "The pipeline did not accept the query sequence."),
    (re.compile(r"No space left on device"), "no_disk_space",
     "The server ran out of disk space during the search. Please start it again later."),
    (re.compile(r"OutOfMemoryError|Cannot fit databases into|exit status 137|"
                r"terminated with an error exit status \(137\)"), "out_of_memory",
     "A step of the search was stopped by the system, most likely because it needed more "
     "memory than this server has for one search."),
    (re.compile(r"Failed to create Conda environment|conda: command not found|"
                r"parasail is not importable|ModuleNotFoundError|No module named"), "server_setup",
     "The server could not prepare its software for the search. The operator has to fix this."),
)
_FAILED_PROCESS_RE = re.compile(r"Error executing process > '([A-Za-z0-9_:]+)")
_STEP_NAMES = {
    "locate": "locating the gene in its home genome",
    "neighbours": "reading the neighbouring genes",
    "search": "searching the genomes",
    "check": "checking the calls",
    "report": "writing figures and summary",
}
_PROCESS_STEP = (
    (("ITERATIVE_SEARCH",), "search"),
    (("EXTRACT_FLANKING", "PREPARE_INITIAL_DB", "PHYLO_SORT", "FILTER_SORTED_GENOMES"), "neighbours"),
    (("GENERATE_REPORT", "PLOT_SYNTENY", "PLOT_SYNTENY_MATRIX", "COMPUTE_TREE"), "report"),
    (("CLUSTER_REGIONS", "RESCUE", "OWNERSHIP", "PARALOG", "PHYLO_PLACEMENT", "STAGE_REGION"), "check"),
)


class WorkerError(RuntimeError):
    """The worker cannot run at all (second worker, unusable data folder)."""


class JobProblem(Exception):
    """A search cannot be started for a reason that is its own."""

    def __init__(self, code, message):
        Exception.__init__(self, message)
        self.code = code
        self.message = message


class JobDir:
    def __init__(self, jobs_dir, job_id):
        if not JOB_ID_RE.match(job_id):
            raise ValueError("not a search id: %r" % job_id)
        self.root = Path(jobs_dir)
        self.path = self.root / job_id
        self.input = self.path / "input"
        self.run = self.path / "run"
        self.work = self.run / "work"
        self.out = self.path / "out"
        self.logs = self.path / "logs"
        self.state = self.path / "state"
        self.result = self.path / "result"
        self.console = self.logs / "console.log"

    def trace(self, attempt):
        return self.logs / ("trace_%d.tsv" % attempt)

    def traces(self):
        return sorted(self.logs.glob("trace_*.tsv"))

    def _read_int(self, name):
        try:
            return int((self.state / name).read_text().strip())
        except (OSError, ValueError):
            return None

    def exit_code(self):
        return self._read_int("exit_code")

    def launched_at(self):
        return self._read_int("launched_at")


def last_attempt(text):
    """The part of the console log written by the latest attempt."""
    index = text.rfind(ATTEMPT_MARK)
    return text if index < 0 else text[index:]


class Worker:
    GC_EVERY = 600
    PREFLIGHT_EVERY = 3600
    PREFLIGHT_RETRY = 60
    MAX_ERRORS = 3
    SETUP_PAUSE_SECONDS = 1800      # after a search failed on broken software
    AUTO_PAUSES = ("low_disk", "genome_storage")

    def __init__(self, ctx, log=None, disk_free=None, sleep=time.sleep):
        self.ctx = ctx
        self.settings = ctx.settings
        self.store = ctx.store
        self.log = log or self._default_log
        self.disk_free = disk_free or (lambda path: shutil.disk_usage(str(path)).free)
        self.sleep = sleep
        self.genomes = ctx.genomes  # genomes fetched from NCBI (genomes.py)
        self.genomes.log = self.log
        self._fetch_thread = None   # one genome is fetched at a time
        self._fetch_stop = threading.Event()
        self._taxonomy_due = 0.0    # when to look for genomes without a lineage again
        self._handles = {}          # id -> Popen, for searches this process started
        self._errors = {}           # id -> errors in a row while watching it
        self._last_gc = 0
        self._last_forced_gc = 0
        self._preflight_due = 0
        self._lock_file = None
        self.settings.jobs_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _default_log(message):
        sys.stderr.write("%s worker: %s\n" % (time.strftime("%Y-%m-%dT%H:%M:%S"), message))
        sys.stderr.flush()

    def jobdir(self, job_id):
        return JobDir(self.settings.jobs_dir, job_id)

    # ── start-up ────────────────────────────────────────────────────────────
    def lock(self):
        """Two workers on one data folder would start the same search twice."""
        path = Path(self.settings.data_dir) / "worker.lock"
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = open(str(path), "a+")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            handle.close()
            raise WorkerError("another worker is already running on %s" % self.settings.data_dir)
        handle.seek(0)
        handle.truncate()
        handle.write("%d\n" % os.getpid())
        handle.flush()
        self._lock_file = handle

    def child_env(self, job_id=None):
        """Environment of the pipeline. The service's own settings (which may hold the
        mail password) are not passed on."""
        env = {k: v for k, v in os.environ.items() if not k.startswith("SYNVOY_WEB_")}
        for key in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP"):
            env.pop(key, None)
        if job_id:
            env[procs.MARKER] = job_id
        env["SYNVOY_WEB_CONDA_ENV"] = self.settings.conda_env
        env["SYNVOY_WEB_CONDA_BASE"] = self.settings.conda_base
        env["NXF_CONDA_CACHEDIR"] = str(self.settings.conda_cache_path)
        env["NXF_ANSI_LOG"] = "false"
        env["NXF_DISABLE_CHECK_LATEST"] = "true"
        env["NO_COLOR"] = "1"
        return env

    def _run_wrapped(self, argv, timeout):
        """Run a short command the way a search is run (same environment set-up)."""
        scratch = Path(self.settings.data_dir) / "preflight"
        scratch.mkdir(parents=True, exist_ok=True)
        try:
            done = subprocess.run(["bash", str(WRAPPER), str(scratch)] + list(argv),
                                  capture_output=True, text=True, timeout=timeout,
                                  env=self.child_env(), stdin=subprocess.DEVNULL,
                                  cwd=str(scratch))
            return done.returncode, (done.stdout + done.stderr).strip()
        except subprocess.TimeoutExpired:
            return 124, "no answer within %d s" % timeout
        except OSError as exc:
            return 127, str(exc)

    def preflight(self):
        """Can a search be run at all? The answer is stored for the API, which stops
        taking searches while it is no."""
        settings = self.settings
        problems = []
        data_dir = Path(settings.data_dir)
        try:
            probe = data_dir / ".write_test"
            probe.write_text("ok")
            probe.unlink()
        except OSError as exc:
            problems.append("the data folder %s is not writable: %s" % (data_dir, exc))
        library = self.ctx.library()
        if self.ctx.libraries.error:
            problems.append(self.ctx.libraries.error)
        elif not library.homes() and not self.ctx.genomes.enabled:
            # With the genome store on, a visitor can name a starting genome at NCBI.
            problems.append("no genome with an annotation on this server: register one with "
                            "`library add --gff`, or set genome_store_gb so that genomes "
                            "are fetched from NCBI")
        nextflow = ""
        if not settings.pipeline_command:
            if not (Path(settings.repo_dir) / "main.nf").is_file():
                problems.append("no main.nf in %s (repo_dir)" % settings.repo_dir)
            code, output = self._run_wrapped(["nextflow", "-version"], 180)
            if code != 0:
                problems.append("Nextflow cannot be started: %s" % (output[-400:] or "no output"))
            else:
                match = re.search(r"version\s+([0-9][^\s]*)", output)
                nextflow = match.group(1) if match else ""
            if settings.conda_env:
                # The pipeline builds the environment of its tasks with conda.
                code, output = self._run_wrapped(["conda", "--version"], 60)
                if code != 0:
                    problems.append("conda is not on the PATH of a search: %s" % output[-200:])
                problem = self.probe_task_environment()
                if problem:
                    problems.append(problem)
        result = {"ok": not problems, "problems": problems, "nextflow": nextflow,
                  "checked_at": self.ctx.now()}
        self.store.set_meta("worker_preflight", result)
        wait = self.PREFLIGHT_EVERY if result["ok"] else self.PREFLIGHT_RETRY
        self._preflight_due = time.monotonic() + wait
        for problem in problems:
            self.log("NOT READY: " + problem)
        for problem in library.problems:
            self.log("genome library: " + problem)
        return result

    def probe_task_environment(self):
        """Would a task of the pipeline get the right Python? '' when it would, or when
        the pipeline has not built its environment yet (the first search does that)."""
        try:
            built = [p for p in self.settings.conda_cache_path.glob("env-*")
                     if (p / "bin" / "python3").exists()]
        except OSError:
            built = []
        if not built:
            return ""
        newest = max(built, key=lambda p: p.stat().st_mtime)
        code, output = self._run_wrapped(["bash", str(TASK_PROBE), str(newest)], 120)
        lines = output.strip().splitlines()
        used = lines[-1] if lines else ""
        if code != 0:
            return ("the tasks of a search cannot load the pipeline's Python packages "
                    "(parasail): %s" % (output[-300:] or "no output"))
        if not used.startswith(str(newest)):
            return ("the tasks of a search would run %s instead of the Python of the "
                    "pipeline's environment; something in PATH is in front of it" % used)
        return ""

    # ── the loop ────────────────────────────────────────────────────────────
    def run_forever(self, stop):
        """``stop`` is a threading.Event. Searches that are running when the worker
        stops are left alone; the next worker takes them over."""
        self.lock()
        # This is the only worker, so no download of a worker is alive. One started by
        # hand (`library fetch`) may be.
        self.genomes.recover(alive=self.fetch_owner_alive)
        self.preflight()
        self.log("ready (data folder %s)" % self.settings.data_dir)
        while not stop.is_set():
            try:
                self.tick()
            except Exception:                              # the loop must not die
                self.log("error in the main loop:\n" + traceback.format_exc())
            stop.wait(self.settings.poll_seconds)
        self.stop_fetching()
        self.log("stopping; %d search(es) keep running" % self.store.count_state(db.RUNNING))

    def tick(self):
        now = self.ctx.now()
        self.store.set_meta("worker_heartbeat", now)
        if time.monotonic() >= self._preflight_due:
            self.preflight()
        for job in self.store.jobs_in_state(db.RUNNING):
            self._guarded(self.supervise, job)
        self.poll_fetch()
        self.start_waiting()
        if time.monotonic() - self._last_gc >= self.GC_EVERY:
            self._last_gc = time.monotonic()
            self.collect_garbage()

    def _guarded(self, function, job):
        job_id = job["id"]
        try:
            function(job)
            self._errors.pop(job_id, None)
        except Exception:
            count = self._errors.get(job_id, 0) + 1
            self._errors[job_id] = count
            self.log("error while handling %s (%d of %d):\n%s"
                     % (job_id, count, self.MAX_ERRORS, traceback.format_exc()))
            if count >= self.MAX_ERRORS:
                try:
                    procs.terminate(job_id, self.settings.kill_grace_seconds,
                                    self._handles.get(job_id), self.sleep)
                    self.finalise(job, db.FAILED, "internal_error",
                                  "The server could not complete this search because of an "
                                  "internal error.")
                except Exception:
                    self.log("could not mark %s as failed:\n%s" % (job_id, traceback.format_exc()))

    # ── genomes that are not on the server yet ──────────────────────────────
    def poll_fetch(self):
        """Fetch the next genome somebody waits for. One at a time, in a thread of this
        process: a download is input and output only, and what has arrived of it is
        continued by the next worker after a restart."""
        thread = self._fetch_thread
        if thread is not None and thread.is_alive():
            return
        self._fetch_thread = None
        if not self.genomes.enabled or self._fetch_stop.is_set():
            return
        row = self.genomes.next_wanted()
        if row is not None:
            thread = threading.Thread(
                target=self.genomes.fetch, name="genome-fetch", daemon=True,
                args=(row["accession"], self._fetch_stop.is_set, "worker:%d" % os.getpid()))
        elif time.monotonic() >= self._taxonomy_due:
            # Nothing to fetch: a genome that was stored while NCBI's taxonomy did not
            # answer gets its lineage now. Looked at every ten minutes.
            self._taxonomy_due = time.monotonic() + 600
            if not self.genomes.without_taxonomy():
                return
            thread = threading.Thread(target=self.genomes.repair_taxonomy,
                                      name="genome-taxonomy", daemon=True)
        else:
            return
        self._fetch_thread = thread
        thread.start()

    @staticmethod
    def fetch_owner_alive(owner):
        kind, _, pid = str(owner).partition(":")
        if kind != "cli" or not pid.isdigit():
            return False
        try:
            os.kill(int(pid), 0)
            return True
        except OSError:
            return False

    def stop_fetching(self, wait=15.0):
        self._fetch_stop.set()
        thread = self._fetch_thread
        if thread is not None:
            thread.join(wait)

    # ── starting ────────────────────────────────────────────────────────────
    def blocked_reason(self, job, running_need_mb):
        """Why this search cannot be started right now (None: it can)."""
        settings = self.settings
        free = self.disk_free(settings.data_dir)
        need = (job["need_mb"] + running_need_mb) * 1e6 + settings.min_free_gb * 1e9
        if free < need:
            return "low_disk"
        library = self.ctx.library()
        request = job["request"] or {}
        for ref in [request.get("home") or {}] + list(request.get("targets") or []):
            genome = library.get(ref.get("id", ""))
            if genome is not None and not genome.available:
                return "genome_storage"
        return None

    def start_waiting(self):
        settings = self.settings
        preflight = self.store.get_meta("worker_preflight") or {}
        if not preflight.get("ok"):
            return
        paused = self.store.get_meta("paused") or {}
        if paused.get("code") == "operator":
            return                                   # closed by hand; `resume` opens it
        if paused.get("code") == "server_setup":
            if self.ctx.now() < int(paused.get("until") or 0):
                return
            self.log("trying the next search after the set-up failure")
            self.store.set_meta("paused", None)
            paused = {}
        running = self.store.jobs_in_state(db.RUNNING)
        running_need = sum(j["need_mb"] for j in running)
        count = len(running)
        blocked = None
        for job in self.store.jobs_in_state(db.QUEUED):
            if count >= settings.max_running:
                break
            state, detail = self.genomes.readiness(job["request"])
            if state == "failed":
                self.finalise(job, db.FAILED, "genome_download_failed", detail)
                continue
            if state == "waiting":
                continue                 # its genomes are on their way; others may start
            reason = self.blocked_reason(job, running_need)
            if reason:
                blocked = blocked or reason
                continue
            claimed = self.store.claim(job["id"], self.ctx.now())
            if claimed is None:
                continue
            self._guarded(self.launch, claimed)
            count += 1
            running_need += claimed["need_mb"]
        if blocked and count == 0:
            # Nothing runs and the next search cannot start: waiting alone will not help.
            if blocked == "low_disk" and time.monotonic() - self._last_forced_gc > 60:
                self._last_forced_gc = time.monotonic()
                self.collect_garbage()
            messages = {
                "low_disk": "This server is short of disk space. Searches are paused until "
                            "space has been freed.",
                "genome_storage": "The genome storage of this server cannot be read at the "
                                  "moment. Searches are paused.",
            }
            if paused.get("code") != blocked:
                self.log("PAUSED: " + blocked)
            self.store.set_meta("paused", {"code": blocked, "message": messages[blocked]})
        elif (self.store.get_meta("paused") or {}).get("code") in self.AUTO_PAUSES:
            self.log("no longer paused")
            self.store.set_meta("paused", None)

    def pause_after_setup_failure(self, job_id):
        """A search failed because the pipeline's software could not be used. The next
        one would fail the same way, so the queue waits: for the operator, or for
        ``SETUP_PAUSE_SECONDS``, whichever comes first."""
        until = self.ctx.now() + self.SETUP_PAUSE_SECONDS
        self.store.set_meta("paused", {
            "code": "server_setup", "until": until,
            "message": "The search service is being repaired. Searches that are waiting "
                       "keep their place."})
        self._preflight_due = 0
        self.log("PAUSED until %s: the search %s failed because the pipeline's software could "
                 "not be used (see its logs). After the repair: python -m synvoy_web resume"
                 % (iso(until), job_id))

    def prepare_inputs(self, jd, job):
        """Write the query and link the library genomes into the search folder."""
        library = self.ctx.library()
        request = job["request"]
        for folder in (jd.input / "home", jd.input / "targets", jd.run, jd.logs, jd.state):
            folder.mkdir(parents=True, exist_ok=True)
        # Nextflow knows an input file by path, size and modification time. A second
        # attempt must be handed the very same query file, or `-resume` runs every
        # step that reads it again, which is nearly all of them. So the file keeps the
        # time the search was created, however often it is written.
        query_file = jd.input / "query.faa"
        try:
            unchanged = query_file.read_text() == request["query"]["fasta"]
        except (OSError, ValueError):
            unchanged = False
        if not unchanged:
            write_atomic(query_file, request["query"]["fasta"])
        os.utime(str(query_file), (job["created_at"], job["created_at"]))

        def genome_of(ref):
            genome = library.get(ref["id"])
            if genome is None or genome.identity() != ref["identity"]:
                raise JobProblem(
                    "genome_changed",
                    "The genome of %s was replaced or removed on the server after this "
                    "search was started. Please start the search again." % ref["species"])
            if not genome.available:
                raise JobProblem("genome_unavailable",
                                 "The genome of %s cannot be read at the moment." % ref["species"])
            return genome

        def link(source, target):
            if target.is_symlink() or target.exists():
                target.unlink()
            os.symlink(source, str(target))

        home = genome_of(request["home"])
        if not home.gff:
            raise JobProblem("genome_changed", "The annotation of %s is no longer on the "
                             "server. Please start the search again." % home.species)
        home_fasta = jd.input / "home" / (home.stem + Path(home.fasta).suffix)
        home_gff = jd.input / "home" / (home.stem + ".gff")
        link(home.fasta, home_fasta)
        link(home.gff, home_gff)
        with_gffs = False
        for ref in request["targets"]:
            genome = genome_of(ref)
            link(genome.fasta, jd.input / "targets" / (genome.stem + Path(genome.fasta).suffix))
            if self.settings.use_target_annotations and genome.gff:
                (jd.input / "target_gffs").mkdir(exist_ok=True)
                link(genome.gff, jd.input / "target_gffs" / (genome.stem + ".gff"))
                with_gffs = True
        return home_fasta, home_gff, with_gffs

    def command(self, jd, job, inputs, resume):
        settings = self.settings
        home_fasta, home_gff, with_gffs = inputs
        argv = list(settings.pipeline_command) or [
            "nextflow", "run", str(Path(settings.repo_dir) / "main.nf")]
        argv += ["-profile", settings.nextflow_profile]
        if settings.nextflow_config:
            argv += ["-c", settings.nextflow_config]
        argv += ["-ansi-log", "false",
                 "-work-dir", str(jd.work),
                 "-with-trace", str(jd.trace(job["attempt"])),
                 "--mode", "pro",
                 "--query", str(jd.input / "query.faa"),
                 "--home_genome", str(home_fasta),
                 "--home_gff", str(home_gff),
                 "--home_species", job["request"]["home"]["species"],
                 "--target_genomes", str(jd.input / "targets"),
                 "--outdir", str(jd.out)]
        if with_gffs:
            argv += ["--target_gffs", str(jd.input / "target_gffs")]
        argv += list(settings.pipeline_args)
        if resume:
            argv.append("-resume")
        return argv

    def launch(self, job):
        """Start the pipeline for a search that has just been claimed."""
        job_id = job["id"]
        jd = self.jobdir(job_id)
        now = self.ctx.now()
        updates = {}
        # A genome that was fetched for this search is only now known for what it is.
        request, settled = self.genomes.settle(job["request"], self.ctx.library())
        if settled:
            job = dict(job, request=request)
            updates["request"] = request
        try:
            inputs = self.prepare_inputs(jd, job)
        except JobProblem as problem:
            self.finalise(job, db.FAILED, problem.code, problem.message)
            return
        version = pipeline_version(self.settings.repo_dir, self.settings.version_label)
        try:
            queued_with = json.loads(job["version"]).get("fingerprint")
        except (TypeError, ValueError, AttributeError):
            queued_with = None
        resume = job["attempt"] > 1 and (jd.run / ".nextflow").is_dir()
        if settled:
            updates["cache_key"] = cache_key(version, request, self.settings)
        if version["fingerprint"] != queued_with:
            # The checkout changed while the search waited: record what really runs.
            updates["version"] = json.dumps(version, sort_keys=True)
            updates["cache_key"] = cache_key(version, job["request"], self.settings)
            if resume:
                # An earlier attempt ran with the older code. Nothing it computed may
                # go into a result that carries the new version: begin again.
                for stale in (jd.run, jd.out):
                    remove_tree(stale, self.settings.jobs_dir)
                jd.run.mkdir(parents=True, exist_ok=True)
                resume = False
                self.log("%s is started from the beginning: the pipeline's code changed "
                         "since its first attempt" % job_id)
        argv = self.command(jd, job, inputs, resume)
        for name in ("exit_code", "exit_code.tmp"):
            try:
                (jd.state / name).unlink()
            except OSError:
                pass
        trace = jd.trace(job["attempt"])
        if trace.exists():
            trace.unlink()
        try:
            with open(str(jd.console), "ab") as console:
                console.write(("\n%s%d, %s =====\n%s\n\n" % (
                    ATTEMPT_MARK, job["attempt"], iso(now),
                    " ".join(shlex.quote(a) for a in argv))).encode("utf-8"))
                console.flush()
                process = subprocess.Popen(
                    ["bash", str(WRAPPER), str(jd.state)] + argv,
                    cwd=str(jd.run), env=self.child_env(job_id), stdin=subprocess.DEVNULL,
                    stdout=console, stderr=subprocess.STDOUT, start_new_session=True,
                    close_fds=True)
        except OSError as exc:
            self.log("could not start %s: %s" % (job_id, exc))
            self.finalise(job, db.FAILED, "internal_error",
                          "The server could not start this search because of an internal error.")
            return
        self._handles[job_id] = process
        if job["attempt"] == 1:              # one use of each genome per search
            self.genomes.record_use(
                [ref["id"] for ref in [request.get("home") or {}] + list(request.get("targets") or [])
                 if ref.get("id")])
        write_atomic(jd.state / "pid", "%d\n" % process.pid)
        write_atomic(jd.state / "launched_at", "%d\n" % now)
        self.store.update_job(job_id, progress=progress_mod.initial_progress(), **updates)
        self.log("started %s (attempt %d%s, pid %d)"
                 % (job_id, job["attempt"], ", resumed" if resume else "", process.pid))

    # ── watching ────────────────────────────────────────────────────────────
    def supervise(self, job):
        job_id = job["id"]
        jd = self.jobdir(job_id)
        handle = self._handles.get(job_id)
        code = jd.exit_code()
        if code is None:
            alive = procs.is_alive(job_id, handle)
            if not alive:
                code = jd.exit_code()            # written between the two looks?
            if code is None and alive:
                self._watch_running(job, jd, handle)
                return
            if code is None:
                self.interrupted(job, jd)
                return
        self.finish(job, jd, code)

    def _watch_running(self, job, jd, handle):
        settings = self.settings
        now = self.ctx.now()
        if job["cancel_requested"]:
            self.stop(job, handle)
            self.finalise(job, db.CANCELLED)
            return
        launched = jd.launched_at() or job["started_at"] or now
        if now - launched > settings.max_runtime_hours * 3600:
            self.stop(job, handle)
            hours = settings.max_runtime_hours
            self.finalise(job, db.FAILED, "timeout",
                          "The search was stopped after %s hours, the time limit of this "
                          "server. Fewer or smaller genomes take less time."
                          % ("%g" % hours))
            return
        current = self.read_progress(jd)
        if current != job["progress"]:
            self.store.update_job(job["id"], progress=current)

    def read_progress(self, jd):
        text = last_attempt(progress_mod.read_tail(jd.console, 2000000))
        return progress_mod.build_progress(text, jd.work, jd.traces())

    def stop(self, job, handle):
        left = procs.terminate(job["id"], self.settings.kill_grace_seconds, handle, self.sleep)
        if left:
            self.log("WARNING: %d process(es) of %s could not be stopped" % (left, job["id"]))
        return left

    def interrupted(self, job, jd):
        """No exit code and no process left: the machine or the search was killed."""
        self.stop(job, self._handles.pop(job["id"], None))
        if job["cancel_requested"]:
            self.finalise(job, db.CANCELLED)
        elif job["attempt"] < self.settings.max_attempts:
            self.store.update_job(job["id"], state=db.QUEUED, progress=None)
            self.log("%s was interrupted; queued again (attempt %d of %d was used)"
                     % (job["id"], job["attempt"], self.settings.max_attempts))
        else:
            self.finalise(job, db.FAILED, "interrupted",
                          "The search was interrupted %d times by a restart of the server. "
                          "Please start it again." % job["attempt"])

    def finish(self, job, jd, code):
        """The pipeline has ended with ``code``."""
        handle = self._handles.get(job["id"])
        if handle is not None:
            try:
                handle.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
        self.stop(job, handle)                    # nothing should be left; make sure
        if code == 0:
            current = self.read_progress(jd)
            meta = self.result_meta(job, current)
            if self.code_changed_since_launch(job):
                # The checkout was updated under a running search: later steps ran
                # newer scripts than earlier ones. Say so, and never hand this result
                # to someone else as the answer to the same question.
                meta["notes"] = [{"kind": "code_changed", "text": (
                    "The software of this server was updated while this search was running, "
                    "so parts of the result may come from the older version. Start the "
                    "search again for a clean result.")}]
                self.store.update_job(job["id"], cache_key="mixed-code:" + job["id"])
                self.log("WARNING: the pipeline's code changed while %s was running" % job["id"])
            try:
                results.build_bundle(jd.out, jd.result, meta)
            except (results.ResultError, OSError, ValueError) as exc:
                self.log("%s finished without a usable result: %s" % (job["id"], exc))
                self.finalise(job, db.FAILED, "no_result",
                              "The search ran to the end but left no readable result.")
                return
            self.finalise(job, db.DONE, progress=progress_mod.finished_progress(current))
            return
        if job["cancel_requested"]:
            self.finalise(job, db.CANCELLED)
            return
        error_code, message = self.classify_failure(jd, job, code)
        self.finalise(job, db.FAILED, error_code, message, progress=self.read_progress(jd))
        if error_code == "server_setup":
            self.pause_after_setup_failure(job["id"])

    def code_changed_since_launch(self, job):
        try:
            launched_with = json.loads(job["version"]).get("fingerprint")
        except (TypeError, ValueError, AttributeError):
            return False
        return bool(launched_with) and \
            pipeline_version(self.settings.repo_dir, self.settings.version_label)["fingerprint"] != launched_with

    def result_meta(self, job, current):
        request = job["request"]
        now = self.ctx.now()
        started = job["started_at"] or now
        return {
            "id": job["id"],
            "scope": "",
            "query": {k: request["query"].get(k) for k in
                      ("source", "accession", "name", "gene", "organism", "length", "label",
                       "fasta")},
            "home": {k: request["home"].get(k) for k in
                     ("id", "species", "common_name", "assembly")},
            "targets": [{k: t.get(k) for k in ("stem", "species", "common_name")}
                        for t in request["targets"]],
            "pipeline": dict(job_version(job), profile=self.settings.nextflow_profile),
            "timing": {"started_at": iso(started), "finished_at": iso(now),
                       "runtime_s": max(0, now - started)},
            "failed_tasks": (current or {}).get("tasks_failed", 0),
        }

    def classify_failure(self, jd, job, code):
        """(error code, message for the visitor) for a pipeline that ended with an error."""
        if code == SETUP_EXIT_CODE:
            return "server_setup", ("The server could not prepare its software for the search. "
                                    "The operator has to fix this.")
        text = last_attempt(progress_mod.read_tail(jd.console, 400000))
        text += "\n" + progress_mod.read_tail(jd.run / ".nextflow.log", 400000)
        home = (job["request"].get("home") or {}).get("species") or "the starting species"
        for pattern, error_code, message in _FAILURES:
            if pattern.search(text):
                return error_code, message.format(home=home)
        match = _FAILED_PROCESS_RE.search(text)
        if match:
            process = match.group(1).split(":")[-1]
            for names, step in _PROCESS_STEP:
                if any(name in process for name in names):
                    return "pipeline_failed", ("The search stopped with an error while %s."
                                               % _STEP_NAMES[step])
            return "pipeline_failed", "The search stopped with an error in the step %s." % process
        return "pipeline_failed", "The search stopped with an error."

    # ── ending ──────────────────────────────────────────────────────────────
    def finalise(self, job, state, error_code=None, error_text=None, progress=None):
        job_id = job["id"]
        now = self.ctx.now()
        started = job["started_at"] or now
        fields = {"state": state, "finished_at": now, "runtime_s": max(0, now - started),
                  "cancel_requested": False, "error_code": error_code, "error_text": error_text}
        if progress is not None:
            fields["progress"] = progress
        if state == db.DONE:
            fields["expires_at"] = now + self.settings.retention_days * 86400
        self.store.update_job(job_id, **fields)
        self._handles.pop(job_id, None)
        self._errors.pop(job_id, None)
        self.log("%s -> %s%s" % (job_id, state, " (%s)" % error_code if error_code else ""))
        try:
            self.clean_up(job_id)
        except Exception:
            self.log("could not clean up %s (will be retried):\n%s"
                     % (job_id, traceback.format_exc()))
        self.send_mail(job, state, error_text)

    def clean_up(self, job_id):
        """Delete what a finished search no longer needs. The result and the logs stay."""
        jd = self.jobdir(job_id)
        nextflow_log = jd.run / ".nextflow.log"
        if nextflow_log.is_file():
            jd.logs.mkdir(parents=True, exist_ok=True)
            (jd.logs / "nextflow.log").write_text(
                progress_mod.read_tail(nextflow_log, MAX_CONSOLE_BYTES))
        if jd.console.is_file() and jd.console.stat().st_size > MAX_CONSOLE_BYTES:
            jd.console.write_text(progress_mod.read_tail(jd.console, MAX_CONSOLE_BYTES))
        for folder in (jd.run, jd.out, jd.input, jd.state):
            remove_tree(folder, self.settings.jobs_dir)

    def send_mail(self, job, state, error_text):
        try:
            addresses = self.store.emails_for(job["id"])
            if addresses and self.settings.email_enabled and state in (db.DONE, db.FAILED):
                request = job["request"]
                title = "%s in %d genome%s" % (
                    request["query"].get("name") or request["query"].get("label") or "search",
                    len(request["targets"]), "" if len(request["targets"]) == 1 else "s")
                for address in addresses:
                    message = notify.compose(self.settings, address, job["id"], title,
                                             "done" if state == db.DONE else "failed",
                                             error_text or "")
                    notify.send(self.settings, message, self.log)
        except Exception:
            self.log("e-mail step failed for %s:\n%s" % (job["id"], traceback.format_exc()))
        finally:
            # The search is over: nothing links it to a visitor from here on.
            self.store.forget_visitors(job["id"])

    # ── housekeeping ────────────────────────────────────────────────────────
    @staticmethod
    def without_sequence(job):
        """The stored request with the query sequence taken out. What stays (name,
        genomes) is enough to say which search a link once pointed to."""
        request = json.loads(json.dumps(job["request"] or {}))
        if isinstance(request.get("query"), dict):
            request["query"]["fasta"] = ""
        return request

    def collect_garbage(self):
        """Expire old results, remove what crashes left behind, trim small tables."""
        settings, store = self.settings, self.store
        now = self.ctx.now()
        for job in store.jobs_in_state(db.DONE):
            if job["expires_at"] is not None and job["expires_at"] <= now:
                remove_tree(self.jobdir(job["id"]).path, settings.jobs_dir)
                store.update_job(job["id"], state=db.EXPIRED, request=self.without_sequence(job))
                self.log("%s expired" % job["id"])
        old = now - settings.failed_retention_days * 86400
        for job in store.jobs_in_state(db.FAILED, db.CANCELLED):
            jd = self.jobdir(job["id"])
            if (job["finished_at"] or 0) > old:
                continue
            if jd.path.exists():
                remove_tree(jd.path, settings.jobs_dir)
            if ((job["request"] or {}).get("query") or {}).get("fasta"):
                store.update_job(job["id"], request=self.without_sequence(job))
        # Working folders of searches that are over (a crash between "done" and clean-up).
        for job in store.jobs_in_state(db.DONE, db.FAILED, db.CANCELLED):
            jd = self.jobdir(job["id"])
            if jd.run.exists() or jd.out.exists() or jd.input.exists():
                try:
                    self.clean_up(job["id"])
                except Exception:
                    self.log("clean-up of %s failed:\n%s" % (job["id"], traceback.format_exc()))
        # Rows nobody can reach any more.
        ancient = now - 400 * 86400
        for job in store.jobs_in_state(db.EXPIRED, db.FAILED, db.CANCELLED):
            if (job["finished_at"] or job["created_at"]) <= ancient:
                store.delete_job(job["id"])
        # Folders without a database row (older than an hour, so a search that is being
        # submitted right now is not touched). Anything not named like a search is left alone.
        try:
            entries = list(settings.jobs_dir.iterdir())
        except OSError:
            entries = []
        for entry in entries:
            if not JOB_ID_RE.match(entry.name) or store.get_job(entry.name) is not None:
                continue
            try:
                age = time.time() - entry.stat().st_mtime
            except OSError:
                continue
            if age > 3600:
                self.log("removing folder without a search: %s" % entry.name)
                remove_tree(entry, settings.jobs_dir)
        store.prune_events(now - 2 * 86400)
        store.prune_lookups(now - 60 * 86400)
        try:
            self.genomes.housekeeping()
        except Exception:
            self.log("clean-up of the genome store failed:\n%s" % traceback.format_exc())
