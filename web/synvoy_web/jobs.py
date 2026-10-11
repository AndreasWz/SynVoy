"""Accepting a search: checks, limits, the result cache, and what a visitor may see.

Nothing here starts a process. A search that passes is written to the database as
``queued``; the worker (``runner.py``) picks it up.
"""

import hashlib
import json
import re
import secrets
import subprocess
import threading
import time
from pathlib import Path

from . import db
from .genomes import GenomeStore, StoreProblem
from .library import Library
from .ncbi import normalise_accession
from .query import QueryError, Resolver, classify_accession, clean_protein, safe_name

JOB_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{21}$")
EMAIL_RE = re.compile(r"^[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]{1,64}@[A-Za-z0-9-]{1,63}(\.[A-Za-z0-9-]{1,63})+$")

# A search that failed for one of these reasons fails the same way every time, so an
# identical request gets the stored answer instead of the same hours of computing.
PERMANENT_ERRORS = ("gene_not_in_home", "query_rejected")

# Files whose content decides what a search returns. Their hash is part of the cache
# key, so a result is reused only for the code that produced it, committed or not.
_CODE_PATHS = ("VERSION", "main.nf", "nextflow.config", "environment.yml",
               "bin", "modules", "subworkflows", "conf")
_CODE_SUFFIXES = (".py", ".nf", ".config", ".yml", ".yaml", ".sh", "")


def iso(ts):
    if ts is None:
        return None
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(int(ts)))


def _git(repo, *args):
    try:
        out = subprocess.run(["git", "-C", str(repo)] + list(args), capture_output=True,
                             text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def code_fingerprint(repo):
    """Hash over the pipeline's code and settings files (names and contents)."""
    repo = Path(repo)
    digest = hashlib.sha256()
    files = []
    for entry in _CODE_PATHS:
        path = repo / entry
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            for sub in sorted(path.rglob("*")):
                if sub.is_file() and "__pycache__" not in sub.parts \
                        and sub.suffix in _CODE_SUFFIXES:
                    files.append(sub)
    for path in files:
        digest.update(str(path.relative_to(repo)).encode("utf-8") + b"\0")
        try:
            digest.update(path.read_bytes())
        except OSError:
            digest.update(b"<unreadable>")
        digest.update(b"\0")
    return digest.hexdigest()[:16]


PRE_RELEASE = "pre-release"


def pipeline_version(repo, label=""):
    """What is shown on every result page, plus the fingerprint used for the cache.

    A version number is shown only for a release: the checkout stands exactly on the
    git tag of its VERSION file (``v1.0.0`` or ``1.0.0``) and has no local changes.
    Everything else is "pre-release" and is told apart by its commit; the number in
    VERSION alone says nothing about what is running. ``label`` (the setting
    ``version_label``) is shown instead, for an installation without git history."""
    repo = Path(repo)
    try:
        version = (repo / "VERSION").read_text().strip() or "unknown"
    except OSError:
        version = "unknown"
    commit = _git(repo, "rev-parse", "--short=10", "HEAD") or ""
    status = _git(repo, "status", "--porcelain", "--untracked-files=no")
    tag = _git(repo, "describe", "--tags", "--exact-match", "HEAD") or ""
    release = version != "unknown" and tag in (version, "v" + version) and status == ""
    return {
        "version": version,
        "commit": commit,
        "modified": bool(status),
        "release": release,
        "label": str(label or "").strip() or (version if release else PRE_RELEASE),
        "fingerprint": code_fingerprint(repo),
    }


def cache_key(version, request, settings):
    """Same key = same question to the same code with the same settings."""
    config_hash = ""
    if settings.nextflow_config:
        try:
            config_hash = hashlib.sha256(Path(settings.nextflow_config).read_bytes()).hexdigest()
        except OSError:
            config_hash = "unreadable"
    material = {
        "v": 1,
        "code": version["fingerprint"],
        "run": [settings.nextflow_profile, config_hash, list(settings.pipeline_args),
                bool(settings.use_target_annotations)],
        "query": request["query"]["fasta"],
        "home": request["home"]["identity"],
        "targets": sorted(t["identity"] for t in request["targets"]),
    }
    text = json.dumps(material, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def new_job_id():
    while True:
        job_id = secrets.token_urlsafe(16)
        if JOB_ID_RE.match(job_id):
            return job_id


def token_hash(token):
    return hashlib.sha256(("synvoy-web-token:" + token).encode("utf-8")).hexdigest()


def wrap_fasta(header, sequence, width=60):
    lines = [">" + header]
    lines += [sequence[i:i + width] for i in range(0, len(sequence), width)]
    return "\n".join(lines) + "\n"


class Context:
    """Everything a request or the worker needs, in one place."""

    def __init__(self, settings, store, libraries, resolver=None, clock=time.time, genomes=None):
        self.settings = settings
        self.store = store
        self.libraries = libraries
        self.clock = clock
        self.resolver = resolver or Resolver(store, clock=clock)
        self.genomes = genomes or GenomeStore(settings, store, clock=clock)
        self.genomes.own_files = self._library_files
        self._secret = None
        self._version = None
        self._version_at = 0.0
        self._lock = threading.Lock()

    def library(self):
        """The genomes a search may use: those the operator registered, and those
        fetched from NCBI that are on the server right now."""
        return self.genomes.merged(self.libraries.current())

    def _library_files(self):
        """The files of the genomes the operator registered, and whether that list is
        complete (it is not while the library file cannot be read). The genome store
        leaves them alone, also when they lie in its folder."""
        library = self.libraries.current()
        files = [path for genome in library.genomes for path in (genome.fasta, genome.gff) if path]
        return files, not self.libraries.error

    def now(self):
        return int(self.clock())

    def version(self, max_age=60.0):
        """Pipeline version, re-read once a minute so an updated checkout is noticed."""
        with self._lock:
            if self._version is None or time.monotonic() - self._version_at > max_age:
                self._version = pipeline_version(self.settings.repo_dir,
                                                 self.settings.version_label)
                self._version_at = time.monotonic()
            return self._version

    def client_key(self, address):
        """Stable, non-reversible stand-in for a visitor's address."""
        if self._secret is None:
            self._secret = self.settings.secret()
        text = self._secret + "|" + (address or "unknown")
        return hashlib.sha256(text.encode("utf-8")).hexdigest()[:32]


# ── is the service taking searches? ─────────────────────────────────────────

def service_state(ctx):
    """``(accepting, code, message)``. New searches are refused, with a reason a
    visitor can read, whenever one could not be run right now."""
    settings, store = ctx.settings, ctx.store
    if not settings.accepting:
        return False, "closed", (settings.closed_message
                                 or "New searches are switched off at the moment.")
    if ctx.libraries.error:
        return False, "library_error", "The genome library of this server cannot be read."
    if not ctx.library().genomes and not ctx.genomes.enabled:
        return False, "library_empty", "This server has no genomes to search yet."
    preflight = store.get_meta("worker_preflight")
    if preflight is not None and not preflight.get("ok"):
        return False, "worker_not_ready", "The search service is being repaired."
    beat = store.get_meta("worker_heartbeat")
    if beat is None or ctx.now() - int(beat) > settings.worker_stale_seconds:
        return False, "worker_down", "The search service is not running at the moment."
    paused = store.get_meta("paused")
    if paused:
        return False, "paused", str(paused.get("message") or "Searches are paused.")
    return True, "", ""


def status_view(ctx):
    settings = ctx.settings
    accepting, code, message = service_state(ctx)
    version = ctx.version()
    return {
        "service": "synvoy-web",
        "api": "v1",
        "accepting": accepting,
        "reason": code,
        "message": message,
        "pipeline": shown_version(version),
        "limits": {
            "max_targets": settings.max_targets,
            "max_genome_mb": settings.max_genome_mb,
            "max_total_mb": settings.max_total_mb,
            "min_query_length": settings.min_query_length,
            "max_query_length": settings.max_query_length,
            "max_active_per_client": settings.max_active_per_client,
            "max_submissions_per_day": settings.max_submissions_per_day,
            "max_runtime_hours": settings.max_runtime_hours,
            "retention_days": settings.retention_days,
            # May a visitor name a genome that is not on the server?
            "downloads": ctx.genomes.enabled,
            "max_fetches_per_day": settings.max_fetches_per_day,
        },
        "queue": {"waiting": ctx.store.count_state(db.QUEUED),
                  "running": ctx.store.count_state(db.RUNNING)},
        "features": {"email": settings.email_enabled},
        "operator": {"name": settings.operator_name, "contact": settings.contact_email},
    }


# ── look-up of an accession ─────────────────────────────────────────────────

def _check_length(length, settings, field):
    if length < settings.min_query_length:
        raise QueryError("sequence_too_short",
                         "This protein has %d residues. The shortest accepted is %d."
                         % (length, settings.min_query_length), field)
    if length > settings.max_query_length:
        raise QueryError("sequence_too_long",
                         "This protein has %d residues. The longest accepted is %d."
                         % (length, settings.max_query_length), field)


def _resolve(ctx, accession, client):
    """Resolver call with the per-visitor limit on look-ups that leave the server.
    An answer that is already stored costs nothing and is not counted."""
    _, normalised = classify_accession(accession)
    record = ctx.resolver.cached(normalised)
    if record is not None:
        return record
    now = ctx.now()
    if ctx.store.count_events("lookup", client, now - 3600) >= ctx.settings.max_lookups_per_hour:
        raise QueryError("too_many_lookups",
                         "Too many look-ups from your address in the last hour. "
                         "Try again later, or paste the sequence.", "accession", status=429)
    ctx.store.add_event("lookup", client, now)
    record, _ = ctx.resolver.resolve(normalised)
    return record


def lookup_view(ctx, payload, address):
    """Answer for the search form: what the accession is, and which genome on this
    server it belongs to. For pasted text: whether the sequence is acceptable."""
    if not isinstance(payload, dict) or set(payload) - {"accession", "sequence"}:
        raise QueryError("bad_request", "The request is not in the expected form.")
    settings = ctx.settings
    if "sequence" in payload and "accession" not in payload:
        header, sequence = clean_protein(payload["sequence"], settings.min_query_length,
                                         settings.max_query_length)
        return {"source": "sequence", "length": len(sequence),
                "name": safe_name(header.split()[0] if header else "", "query"),
                "home": None, "warnings": _length_warnings(len(sequence))}
    if "accession" not in payload or "sequence" in payload:
        raise QueryError("bad_request", "Send either an accession or a sequence.")
    record = _resolve(ctx, payload["accession"], ctx.client_key(address))
    home = ctx.library().home_for(record.get("taxon_id"), record.get("organism"))
    view = {
        "source": record["source"],
        "accession": record["accession"],
        "name": record["name"],
        "gene": record["gene"],
        "organism": record["organism"],
        "organism_common": record["organism_common"],
        "taxon_id": record.get("taxon_id"),
        "length": record["length"],
        "home": home.public() if home else None,
        "warnings": _length_warnings(record["length"]),
        "acceptable": True,
        "problem": "",
    }
    try:
        _check_length(record["length"], settings, "accession")
    except QueryError as exc:
        view["acceptable"], view["problem"] = False, exc.message
    return view


def _length_warnings(length):
    # Measured on the validation runs: proteins under 100 residues gave few confident calls.
    return ["short_protein"] if length < 100 else []


# ── building and checking a request ─────────────────────────────────────────

def _genome_ref(genome, with_gff=False):
    """A genome as a search remembers it. ``store`` marks a genome of the genome store:
    the worker waits for it to be on the server and records what it finally is."""
    ref = {"id": genome.id, "species": genome.species, "common_name": genome.common_name,
           "assembly": genome.assembly, "stem": genome.stem, "size_bp": genome.size_bp,
           "identity": genome.identity()}
    if genome.origin == "ncbi":
        ref.update({"store": True, "accession": genome.accession, "with_gff": bool(with_gff)})
    return ref


def _is_item(value):
    """A genome in a request: its id on this server, or ``{"accession": "GCF_..."}``."""
    return isinstance(value, str) or (
        isinstance(value, dict) and set(value) == {"accession"}
        and isinstance(value["accession"], str))


def _pick(ctx, library, item, client, field, with_gff=False):
    """The genome a request names, or None when the server does not know the id. A
    genome named by accession that is not here is put on the list of genomes to
    fetch; it comes back marked as not available yet."""
    if isinstance(item, str):
        return library.get(item)
    accession = normalise_accession(item["accession"])
    if not accession:
        raise QueryError("genome_unknown", "%r is not an NCBI assembly accession (it looks "
                         "like GCF_000146045.2)." % item["accession"][:40], field)
    for genome in library.genomes:              # registered by the operator, or in the store
        if genome.assembly == accession and genome.available and (genome.gff or not with_gff):
            return genome
    try:
        row = ctx.genomes.want(accession, with_gff=with_gff, client=client)
    except StoreProblem as exc:
        raise QueryError(exc.code, exc.message, field, status=exc.status)
    return ctx.genomes.as_genome(row)


def targets_by_taxon(ctx, library, taxon_id, address):
    """``(ids of the available genomes, closest relatives of an NCBI taxon first; ids
    of those that are ticked by default)``: the species list when the starting genome
    is not on the server (yet). The lineage is asked of NCBI once and kept."""
    key = "taxon:%d" % taxon_id
    lineage = ctx.store.get_lookup(key, ctx.now() - 180 * 86400)
    if lineage is None:
        client, now = ctx.client_key(address), ctx.now()
        if ctx.store.count_events("lookup", client, now - 3600) >= ctx.settings.max_lookups_per_hour:
            lineage = []
        else:
            ctx.store.add_event("lookup", client, now)
            try:
                _, lineage, _ = ctx.genomes.taxonomy("", taxon_id=taxon_id)
                ctx.store.put_lookup(key, lineage, now)
            except Exception:                    # NCBI does not answer: the list stays unordered
                lineage = []

    class _Taxon:
        pass
    anchor = _Taxon()
    anchor.lineage = lineage
    genomes = sorted((g for g in library.genomes if g.available),
                     key=lambda g: (-Library.shared_depth(anchor, g), g.species.lower()))
    return [g.id for g in genomes], [g.id for g in Library.related(anchor, genomes)]


def _waits_for_fetch(genome):
    return genome.origin == "ncbi" and not genome.available


def find_view(ctx, payload, address):
    """Assemblies at NCBI for a species name, a taxon id or an accession, for the
    search form: which are on this server already, which could be fetched."""
    if not isinstance(payload, dict) or set(payload) - {"term", "purpose"} \
            or not isinstance(payload.get("term"), str) \
            or payload.get("purpose", "target") not in ("target", "home"):
        raise QueryError("bad_request", "The request is not in the expected form.")
    term = re.sub(r"\s+", " ", payload["term"]).strip()
    if not re.match(r"^[A-Za-z0-9][A-Za-z0-9 ._\-]{1,79}$", term):
        raise QueryError("term_invalid", "Type a species name (for example Tetramorium "
                         "bicarinatum) or an assembly accession (GCF_...).", "term")
    client, now = ctx.client_key(address), ctx.now()
    if ctx.store.count_events("lookup", client, now - 3600) >= ctx.settings.max_lookups_per_hour:
        raise QueryError("too_many_lookups", "Too many look-ups from your address in the last "
                         "hour. Try again later.", "term", status=429)
    ctx.store.add_event("lookup", client, now)
    try:
        candidates = ctx.genomes.find(term, payload.get("purpose", "target"), ctx.library())
    except StoreProblem as exc:
        raise QueryError(exc.code, exc.message, "term", status=exc.status)
    return {"term": term, "candidates": candidates, "downloads": ctx.genomes.enabled}


def build_request(ctx, payload, client):
    """Validate the body of a submission and return the stored form of the request."""
    settings, library = ctx.settings, ctx.library()
    if not isinstance(payload, dict):
        raise QueryError("bad_request", "The request is not in the expected form.")
    extra = set(payload) - {"query", "targets", "email"}
    if extra:
        raise QueryError("bad_request", "Unknown field: %s." % ", ".join(sorted(extra)))
    query_in = payload.get("query")
    if not isinstance(query_in, dict) or set(query_in) - {"accession", "sequence", "name", "home"}:
        raise QueryError("bad_request", "The request is not in the expected form.", "query")
    has_acc, has_seq = "accession" in query_in, "sequence" in query_in
    if has_acc == has_seq:
        raise QueryError("bad_request", "Give either an accession or a sequence.", "query")
    home_in = query_in.get("home")
    if home_in is not None and not _is_item(home_in):
        raise QueryError("bad_request", "The request is not in the expected form.", "home")
    # The cheap checks on the genomes come before anything is asked of NCBI.
    target_items = payload.get("targets")
    if not isinstance(target_items, list) or not all(_is_item(t) for t in target_items):
        raise QueryError("bad_request", "Choose the genomes to search.", "targets")
    if not target_items:
        raise QueryError("no_targets", "Choose at least one genome to search.", "targets")
    if len(target_items) > settings.max_targets:
        raise QueryError("too_many_targets",
                         "Choose at most %d genomes per search." % settings.max_targets, "targets")

    if has_acc:
        record = _resolve(ctx, query_in["accession"], client)
        _check_length(record["length"], settings, "accession")
        label = safe_name(record["gene"] or record["accession"], record["accession"])
        query = {"source": record["source"], "accession": record["accession"],
                 "name": record["name"] or record["accession"], "gene": record["gene"],
                 "organism": record["organism"], "length": record["length"], "label": label,
                 "fasta": wrap_fasta(record["header"] or record["accession"], record["sequence"])}
        home = _pick(ctx, library, home_in, client, "home", with_gff=True) if home_in \
            else library.home_for(record.get("taxon_id"), record.get("organism"))
        if home is None and not home_in:
            raise QueryError(
                "home_not_in_library",
                "The genome of %s is not on this server. Choose a related species as the "
                "starting genome%s." % (record["organism"] or "this protein's species",
                                        ", or have its genome fetched from NCBI"
                                        if ctx.genomes.enabled else ""), "home")
    else:
        header, sequence = clean_protein(query_in["sequence"], settings.min_query_length,
                                         settings.max_query_length)
        name_in = query_in.get("name")
        if name_in is not None and not isinstance(name_in, str):
            raise QueryError("bad_request", "The request is not in the expected form.", "name")
        label = safe_name(name_in or (header.split()[0] if header else ""), "query")
        query = {"source": "sequence", "accession": "", "name": label, "gene": "",
                 "organism": "", "length": len(sequence), "label": label,
                 "fasta": wrap_fasta(label, sequence)}
        if not home_in:
            raise QueryError("home_missing",
                             "Choose the species the protein comes from.", "home")
        home = _pick(ctx, library, home_in, client, "home", with_gff=True)

    if home is None:
        raise QueryError("home_unknown", "The starting genome is not on this server.", "home")
    if not home.available and not _waits_for_fetch(home):
        raise QueryError("genome_unavailable",
                         "The genome of %s cannot be used at the moment." % home.species, "home",
                         status=503)
    if not home.gff:
        raise QueryError("home_without_annotation",
                         "%s has no gene annotation on this server, so it cannot be the "
                         "starting genome." % home.species, "home")

    targets = []
    for item in target_items:
        genome = _pick(ctx, library, item, client, "targets")
        if genome is None:
            raise QueryError("target_unknown", "One of the chosen genomes is not on this server.",
                             "targets")
        if genome.id == home.id:
            raise QueryError("target_is_home",
                             "%s is the starting genome and cannot also be searched."
                             % home.species, "targets")
        if genome.id in {g.id for g in targets}:
            raise QueryError("bad_request", "A genome was chosen twice.", "targets")
        if not genome.available and not _waits_for_fetch(genome):
            raise QueryError("genome_unavailable",
                             "The genome of %s cannot be used at the moment." % genome.species,
                             "targets", status=503)
        if genome.size_mb > settings.max_genome_mb:
            raise QueryError("genome_too_large",
                             "The genome of %s (%d Mb) is above the limit of %d Mb per genome."
                             % (genome.species, round(genome.size_mb), settings.max_genome_mb),
                             "targets")
        targets.append(genome)
    # The pipeline knows a genome by its file name, which is made from the species.
    stems = {}
    for genome in [home] + targets:
        if genome.stem in stems:
            raise QueryError("same_species_twice",
                             "Two of the chosen genomes are %s. One search takes one genome "
                             "per species." % genome.species, "targets")
        stems[genome.stem] = genome.id
    total_mb = sum(g.size_mb for g in targets)
    if total_mb > settings.max_total_mb:
        raise QueryError("selection_too_large",
                         "The chosen genomes add up to %d Mb. The limit per search is %d Mb."
                         % (round(total_mb), settings.max_total_mb), "targets")
    # Closest relatives of the starting genome first. A genome that is still to be
    # fetched has no lineage yet and comes last.
    ordered = sorted(targets, key=lambda g: (-Library.shared_depth(home, g), g.species.lower()))
    largest = sorted((g.size_mb for g in targets), reverse=True)[:3]
    need_mb = int(settings.disk_factor * (home.size_mb + sum(largest)))
    return {"query": query, "home": _genome_ref(home, with_gff=True),
            "targets": [_genome_ref(g) for g in ordered],
            "total_mb": round(total_mb, 1)}, need_mb


def check_email(ctx, value):
    if value in (None, ""):
        return None
    if not ctx.settings.email_enabled:
        raise QueryError("email_off", "This server does not send e-mails.", "email")
    if not isinstance(value, str) or len(value) > 254 or not EMAIL_RE.match(value):
        raise QueryError("email_invalid", "This does not look like an e-mail address.", "email")
    return value


def _reusable(prior, now, settings, result_exists):
    if prior is None:
        return False
    if prior["state"] in db.ACTIVE_STATES:
        return True
    if prior["state"] == db.DONE:
        alive = prior["expires_at"] is None or prior["expires_at"] > now
        return alive and result_exists(prior["id"])
    if prior["state"] == db.FAILED and prior["error_code"] in PERMANENT_ERRORS:
        return (prior["finished_at"] or 0) > now - settings.failed_retention_days * 86400
    return False


def submit(ctx, payload, address):
    """Accept a search. Returns ``{"id", "token", "reused", ...}`` or raises QueryError."""
    settings, store = ctx.settings, ctx.store
    client = ctx.client_key(address)
    request, need_mb = build_request(ctx, payload, client)
    email = check_email(ctx, payload.get("email"))
    version = ctx.version()
    key = cache_key(version, request, settings)
    token = secrets.token_urlsafe(24)
    now = ctx.now()

    def result_exists(job_id):
        return (settings.jobs_dir / job_id / "result" / "result.json").is_file()

    with store.transaction():
        prior = store.newest_with_key(key)
        if _reusable(prior, now, settings, result_exists):
            store.add_submitter(prior["id"], token_hash(token), client, email, now)
            if prior["state"] == db.DONE:
                keep_until = now + settings.retention_days * 86400
                if (prior["expires_at"] or 0) < keep_until:
                    store.update_job(prior["id"], expires_at=keep_until)
                store.forget_emails(prior["id"])     # nothing left to announce
            return {"id": prior["id"], "token": token, "reused": True, "state": prior["state"]}

        accepting, code, message = service_state(ctx)
        if not accepting:
            raise QueryError(code, message, status=503)
        if store.count_events("submit", client, now - 86400) >= settings.max_submissions_per_day:
            raise QueryError("too_many_submissions",
                             "You have started %d searches in the last 24 hours, which is the "
                             "limit. Try again tomorrow." % settings.max_submissions_per_day,
                             status=429)
        active = store.active_for_client(client)
        if active >= settings.max_active_per_client:
            raise QueryError("too_many_active",
                             "You already have %d search%s waiting or running. Wait until one "
                             "has finished." % (active, "" if active == 1 else "es"), status=429)
        if store.count_state(db.QUEUED) >= settings.max_queue:
            raise QueryError("queue_full",
                             "The waiting list is full at the moment. Try again later.",
                             status=503)
        job_id = new_job_id()
        store.insert_job({
            "id": job_id, "cache_key": key, "state": db.QUEUED, "created_at": now,
            "queued_at": now, "request": request, "need_mb": need_mb,
            "version": json.dumps(version, sort_keys=True),
        })
        store.add_submitter(job_id, token_hash(token), client, email, now)
        store.add_event("submit", client, now)
    return {"id": job_id, "token": token, "reused": False, "state": db.QUEUED}


def cancel(ctx, job_id, token):
    """Withdraw one visitor's interest. The search itself stops only when nobody else
    is waiting for it."""
    store = ctx.store
    if not isinstance(token, str) or not token:
        raise QueryError("token_missing", "Only the browser that started this search can "
                         "cancel it.", status=403)
    with store.transaction():
        job = store.get_job(job_id)
        if job is None:
            raise QueryError("not_found", "There is no search with this link.", status=404)
        if job["state"] not in db.ACTIVE_STATES:
            raise QueryError("not_active", "This search is no longer waiting or running.",
                             status=409)
        remaining = store.withdraw(job_id, token_hash(token))
        if remaining is None:
            raise QueryError("not_yours", "Only the browser that started this search can "
                             "cancel it.", status=403)
        if remaining > 0:
            return {"cancelled": False, "others_waiting": remaining}
        if job["state"] == db.QUEUED:
            store.update_job(job_id, state=db.CANCELLED, finished_at=ctx.now())
            store.forget_visitors(job_id)
        else:
            store.update_job(job_id, cancel_requested=True)
        return {"cancelled": True, "others_waiting": 0}


# ── what a visitor sees of a search ─────────────────────────────────────────

def shown_version(version):
    """The part of ``pipeline_version`` that pages show. A record made before releases
    were told apart has no label: it was not a release."""
    version = version if isinstance(version, dict) else {}
    return {"label": version.get("label") or PRE_RELEASE, "release": bool(version.get("release")),
            "commit": version.get("commit") or "", "modified": bool(version.get("modified"))}


def job_version(job):
    try:
        version = json.loads(job["version"])
    except (TypeError, ValueError):
        version = {}
    return shown_version(version)


def job_view(ctx, job):
    request = job["request"] or {}
    query = request.get("query") or {}
    view = {
        "id": job["id"],
        "state": job["state"],
        "created_at": iso(job["created_at"]),
        "started_at": iso(job["started_at"]),
        "finished_at": iso(job["finished_at"]),
        "expires_at": iso(job["expires_at"]),
        "attempt": job["attempt"],
        "runtime_s": job["runtime_s"],
        "query": {k: query.get(k) for k in
                  ("source", "accession", "name", "gene", "organism", "length", "label")},
        "home": {k: (request.get("home") or {}).get(k) for k in
                 ("id", "species", "common_name", "assembly")},
        "targets": [{k: t.get(k) for k in ("id", "species", "common_name")}
                    for t in request.get("targets") or []],
        "total_mb": request.get("total_mb"),
        "pipeline": job_version(job),
        "queue": None,
        "progress": None,
        "error": None,
        "cancel_requested": job["cancel_requested"],
    }
    if job["state"] == db.QUEUED:
        view["queue"] = {"ahead": ctx.store.queue_ahead(job),
                         "running": ctx.store.count_state(db.RUNNING),
                         "genomes": ctx.genomes.waiting(request)}
    if job["state"] == db.RUNNING:
        view["progress"] = job["progress"]
    if job["state"] == db.FAILED:
        view["error"] = {"code": job["error_code"] or "failed",
                         "message": job["error_text"] or "The search stopped with an error."}
    return view
