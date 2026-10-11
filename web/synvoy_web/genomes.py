"""Genomes fetched from NCBI and kept in a folder of limited size: the genome store.

The operator gives the store a size (``genome_store_gb``) and usually fills part of it
(``library fetch``). The rest fills as visitors ask for genomes that are not on the
server. Every search that starts counts as one use of each of its genomes. When a new
genome needs room, the genomes used least in the last ``store_score_days`` days go
first; one that a waiting or running search needs, or that the operator has pinned, is
never removed.

A genome of the store is known by its assembly accession with version, which names one
fixed set of sequences. It can be removed and fetched again without changing a result,
and a result that used it stays valid.

States of a genome (table ``store_genomes``):

    wanted     somebody asked for it; the worker will fetch it
    fetching   being downloaded and checked
    ready      on the server, usable in searches
    failed     could not be fetched; the reason is kept for some hours
    evicted    was on the server and has been removed to make room

The genomes of the library file (``library add``) are not part of the store: they stay
where the operator put them and are never removed.
"""

import os
import shutil
import threading
import time
import traceback
from pathlib import Path
from urllib.parse import quote

from . import db as dbm
from .fsutil import remove_tree, tree_bytes
from .library import Genome, Library, fasta_stats, file_stem
from .ncbi import Ncbi, NcbiError, normalise_accession, unpack
from .net import FetchError
from .taxonomy import lookup_taxonomy

WANTED, FETCHING, READY, FAILED, EVICTED = "wanted", "fetching", "ready", "failed", "evicted"
GROUP = {"id": "ncbi", "name": "Added from NCBI",
         "note": "Fetched on request. The least used are removed when room is needed."}
INCOMING = ".incoming"
MAX_ATTEMPTS = 3                 # tries for one genome before it counts as failed
FAILED_FOR = 6 * 3600            # a failure is remembered this long before a new try
ROOM_RETRY = 300                 # no room right now: look again after this many seconds
ROOM_TRIES = 24
UNWANTED_AFTER = 3600            # a request nobody waits for any more is dropped
STALE_FETCH = 900                # a download that has not reported for this long is dead
UNKNOWN_FRESH = 60               # files in no record are measured again after this long


class StoreProblem(Exception):
    """A genome cannot be had. ``message`` is written for the visitor."""

    def __init__(self, code, message, status=400):
        Exception.__init__(self, message)
        self.code, self.message, self.status = code, message, status


def store_id(accession):
    """The library id of a store genome: 'GCF_000146045.2' -> 'gcf_000146045_2'."""
    return accession.lower().replace(".", "_")


def needed_bytes(size_bp, with_gff):
    """Disk a genome takes once unpacked. The annotation of a gene-dense genome can be
    as large as its sequence (yeast: 12.5 MB for 12.3 MB)."""
    return int(size_bp * (1.02 + (1.1 if with_gff else 0.0)))


def store_refs(request):
    """The store genomes a stored search request names."""
    request = request or {}
    refs = [request.get("home") or {}] + list(request.get("targets") or [])
    return [ref for ref in refs if ref.get("store") and ref.get("accession")]


def view(row):
    """What a page may know about a genome that is on its way."""
    return {"accession": row["accession"], "species": row["species"], "state": row["state"],
            "done_bytes": row["done_bytes"], "total_bytes": row["total_bytes"]}


class GenomeStore:
    def __init__(self, settings, store, ncbi=None, clock=time.time, disk_free=None,
                 taxonomy=lookup_taxonomy, log=None):
        self.settings = settings
        self.db = store
        self.ncbi = ncbi or Ncbi()
        self.clock = clock
        self.disk_free = disk_free or (lambda path: shutil.disk_usage(str(path)).free)
        self.taxonomy = taxonomy
        self.log = log or (lambda message: None)
        self._lock = threading.Lock()
        self._merged = (None, None, None)       # (library of the file, generation, both)
        self._unknown_told = ()                 # what housekeeping last reported as unknown
        self._unknown_sizes = ((), 0.0, [])     # (names, when measured, [(name, bytes)])
        # The files of the genomes the operator registered, and whether that list is
        # complete: () -> (paths, complete). Set by the service (jobs.Context).
        self.own_files = lambda: ([], True)

    # ── what the store is ───────────────────────────────────────────────────
    @property
    def enabled(self):
        return self.settings.genome_store_gb > 0

    @property
    def root(self):
        return self.settings.genome_store_path

    @property
    def budget(self):
        return int(self.settings.genome_store_gb * 1e9)

    def _now(self):
        return int(self.clock())

    def _changed(self):
        self.db.set_meta("store_generation", int(self.db.get_meta("store_generation", 0)) + 1)

    def _paths(self, row):
        folder = self.root / row["accession"]
        return folder / (row["stem"] + ".fna"), folder / (row["stem"] + ".gff")

    def _present(self, row, need_gff=False):
        fasta, gff = self._paths(row)
        return fasta.is_file() and (gff.is_file() or not (need_gff or row["with_gff"]))

    def usable(self, row, need_gff=False):
        return (row is not None and row["state"] == READY and (row["with_gff"] or not need_gff)
                and self._present(row, need_gff))

    def as_genome(self, row):
        fasta, gff = self._paths(row)
        genome = Genome()
        genome.id, genome.species, genome.stem = row["id"], row["species"], row["stem"]
        genome.common_name, genome.group = row["common_name"], GROUP["id"]
        genome.assembly, genome.accession, genome.origin = row["accession"], row["accession"], "ncbi"
        genome.taxon_id, genome.lineage = row["taxon_id"], list(row["lineage"])
        genome.size_bp, genome.n_sequences, genome.n50 = row["size_bp"], row["n_sequences"], row["n50"]
        genome.fasta = str(fasta)
        genome.gff = str(gff) if row["with_gff"] else ""
        genome.checksum = ("gff:" + (row.get("gff_md5") or "")) if row["with_gff"] else ""
        genome.available = row["state"] == READY and self._present(row)
        genome.problem = "" if genome.available else "not on the server at the moment"
        return genome

    def merged(self, base):
        """The library of the file plus the genomes of the store that are ready. A
        genome the operator registered wins over the same assembly in the store."""
        generation = self.db.get_meta("store_generation", 0)
        with self._lock:
            if self._merged[0] is base and self._merged[1] == generation:
                return self._merged[2]
        taken = {g.assembly for g in base.genomes if g.assembly}
        rows = [r for r in self.db.store_genomes(READY) if r["accession"] not in taken]
        if not rows and not self.enabled:
            merged = base
        else:
            groups = list(base.groups)
            if rows and not any(g["id"] == GROUP["id"] for g in groups):
                groups.append(dict(GROUP))
            merged = Library(groups, list(base.genomes) + [self.as_genome(r) for r in rows],
                             path=base.path, problems=base.problems)
        with self._lock:
            self._merged = (base, generation, merged)
        return merged

    def unknown(self):
        """What lies in the store's folder without a record: [(name, bytes)], largest
        first. It is left from a database that was lost or replaced, or was put there by
        hand. The store never removes it by itself (``sweep`` does, when the operator
        asks), and counts it against its size, so that the folder as a whole stays
        within the space that was set."""
        try:
            entries = sorted(path.name for path in self.root.iterdir())
        except OSError:
            return []
        skip = {r["accession"] for r in self.db.store_genomes(READY, WANTED, FETCHING)}
        skip.add(INCOMING)
        skip.update(self._operators_own())
        names = tuple(name for name in entries if name not in skip)
        measured, at, found = self._unknown_sizes
        if names != measured or time.monotonic() - at > UNKNOWN_FRESH:
            found = []
            for name in names:
                path = self.root / name
                if path.is_symlink() or path.exists():      # not removed meanwhile
                    found.append((name, tree_bytes(path)))
            found.sort(key=lambda item: (-item[1], item[0]))
            self._unknown_sizes = (names, time.monotonic(), found)
        return list(found)

    def _operators_own(self):
        """Names in the store's folder that hold a file of a genome the operator
        registered (``library add``). They are the operator's, whatever else is true."""
        names = set()
        roots = {os.path.abspath(str(self.root)), os.path.realpath(str(self.root))}
        for path in self.own_files()[0]:
            for form in (os.path.abspath(str(path)), os.path.realpath(str(path))):
                for root in roots:
                    if form.startswith(root + os.sep):
                        names.add(form[len(root) + 1:].split(os.sep)[0])
        return names

    def unknown_bytes(self):
        return sum(size for _, size in self.unknown())

    def sweep(self):
        """Delete what ``unknown`` lists. Only ever called for the operator."""
        if not self.own_files()[1]:
            raise StoreProblem("library_unreadable", "The library file cannot be read, so it "
                               "is not known which files are the operator's own. Nothing was "
                               "deleted.", status=500)
        removed = []
        with self.db.transaction():              # no download is put in place meanwhile
            for name, size in self.unknown():
                path = self.root / name
                if path.is_symlink() or not path.is_dir():
                    path.unlink()
                else:
                    remove_tree(path, self.root)
                removed.append((name, size))
                self.log("genome store: deleted %s (%d MB), which was in no record"
                         % (name, size // 1000000))
        return removed

    def usage(self):
        rows = self.db.store_genomes()
        ready = [r for r in rows if r["state"] == READY]
        coming = [r for r in rows if r["state"] in (WANTED, FETCHING)]
        used = sum(r["disk_bytes"] for r in ready)
        reserved = sum(needed_bytes(r["size_bp"], r["with_gff"]) for r in coming)
        unknown = self.unknown()
        stray = sum(size for _, size in unknown)
        return {"enabled": self.enabled, "budget": self.budget, "used": used, "reserved": reserved,
                "unknown": stray, "unknown_names": [name for name, _ in unknown],
                "free": max(0, self.budget - used - reserved - stray), "ready": len(ready),
                "pinned": sum(1 for r in ready if r["pinned"]), "coming": len(coming)}

    # ── who needs what ──────────────────────────────────────────────────────
    def needed_now(self):
        """Accessions that waiting or running searches name."""
        needed = set()
        for job in self.db.jobs_in_state(dbm.QUEUED, dbm.RUNNING):
            needed.update(ref["accession"] for ref in store_refs(job["request"]))
        return needed

    def record_use(self, genome_ids):
        """One use of every genome of a search that starts."""
        self.db.add_genome_uses(sorted(set(genome_ids)), self._now())

    def use_counts(self):
        since = self._now() - self.settings.store_score_days * 86400
        return self.db.genome_use_counts(since)

    # ── asking for a genome ─────────────────────────────────────────────────
    def describe(self, report, purpose="target", library=None):
        """A candidate from NCBI as the search form shows it: is it here already, and
        if not, could it be fetched."""
        row = self.db.get_store_genome(report["accession"])
        here = None
        if library is not None:
            for genome in library.genomes:
                if genome.assembly == report["accession"] and genome.available \
                        and (genome.gff or purpose != "home"):
                    here = genome.id
                    break
        problem = ""
        size_mb = report["size_bp"] / 1e6
        if here is None:
            if not self.enabled:
                problem = "This server does not fetch genomes."
            elif size_mb > self.settings.max_genome_mb:
                problem = "%d Mb, above the limit of %d Mb per genome." % (
                    round(size_mb), self.settings.max_genome_mb)
            elif purpose == "home" and not report["annotated"]:
                problem = "NCBI has no gene annotation for this assembly."
            elif needed_bytes(report["size_bp"], purpose == "home") > self.budget:
                problem = "Larger than the space this server keeps for genomes."
        return {
            "accession": report["accession"], "species": report["species"],
            "taxon_id": report["taxon_id"],
            "common_name": report["common_name"], "assembly": report["assembly_name"],
            "level": report["level"], "reference": bool(report["category"]),
            "size_mb": round(size_mb, 1), "annotated": report["annotated"],
            "on_server": here, "state": row["state"] if row else None,
            "usable": not problem, "problem": problem,
        }

    def find(self, term, purpose="target", library=None, limit=6):
        """Assemblies at NCBI for a species name, a taxon id or an accession."""
        accession = normalise_accession(term)
        try:
            reports = [self.ncbi.assembly(accession)] if accession \
                else self.ncbi.search(term, limit=limit)
        except NcbiError as exc:
            if exc.kind == "not_found":
                return []
            raise StoreProblem("ncbi_unavailable", "NCBI does not answer at the moment. "
                               "Try again in a few minutes.", status=502)
        return [self.describe(r, purpose, library) for r in reports]

    def _locked_bytes(self, rows, protect, skip):
        """Space that cannot be freed: pinned genomes, genomes a search needs, and
        genomes that are on their way."""
        total = 0
        for row in rows:
            if row["accession"] == skip:
                continue
            if row["state"] == READY and (row["pinned"] or row["accession"] in protect):
                total += row["disk_bytes"]
            elif row["state"] in (WANTED, FETCHING):
                total += needed_bytes(row["size_bp"], row["with_gff"])
        return total

    def want(self, accession, with_gff=False, client="", origin="visitor", pinned=False):
        """Make sure a genome is on the server or on its way, and return its row.
        ``origin``: 'visitor' (counted against the visitor's daily limit), 'search' (a
        waiting search needs it again) or 'operator'. Raises StoreProblem."""
        accession = normalise_accession(accession)
        if not accession:
            raise StoreProblem("genome_unknown", "This is not an NCBI assembly accession "
                               "(it looks like GCF_000146045.2).")
        row = self.db.get_store_genome(accession)
        if self.usable(row, with_gff):
            if pinned and not row["pinned"]:
                self.db.update_store_genome(accession, pinned=True)
                row = self.db.get_store_genome(accession)
            return row
        if not self.enabled:
            raise StoreProblem("downloads_off", "This server does not fetch genomes. Choose "
                               "from the genomes it keeps.")
        now = self._now()
        if row and row["state"] in (WANTED, FETCHING) and (row["with_gff"] or not with_gff):
            return row
        if row and row["state"] == FAILED and now - row["updated_at"] < FAILED_FOR \
                and origin != "operator":
            raise StoreProblem("genome_download_failed",
                               "The genome of %s could not be fetched from NCBI a short while "
                               "ago (%s). Try again later." % (row["species"], row["error_text"]),
                               status=502)
        if row is not None and origin == "search":
            # It was here before: what it is, is on record. No question to NCBI from
            # inside the worker's loop.
            report = {"accession": accession, "species": row["species"],
                      "common_name": row["common_name"], "assembly_name": row["assembly_name"],
                      "level": row["level"], "taxon_id": row["taxon_id"],
                      "size_bp": row["size_bp"], "annotated": row["annotated"]}
        else:
            try:
                report = self.ncbi.assembly(accession)
            except NcbiError as exc:
                if exc.kind == "not_found":
                    raise StoreProblem("genome_unknown", "NCBI has no assembly %s." % accession)
                raise StoreProblem("ncbi_unavailable", "NCBI does not answer at the moment. "
                                   "Try again in a few minutes.", status=502)
        size_mb = report["size_bp"] / 1e6
        if size_mb > self.settings.max_genome_mb and origin != "operator":
            raise StoreProblem("genome_too_large",
                               "The genome of %s (%d Mb) is above the limit of %d Mb per genome."
                               % (report["species"], round(size_mb), self.settings.max_genome_mb))
        if with_gff and not report["annotated"]:
            raise StoreProblem("home_without_annotation",
                               "NCBI has no gene annotation for %s of %s, so it cannot be the "
                               "starting genome." % (accession, report["species"]))
        need = needed_bytes(report["size_bp"], with_gff)
        stray = self.unknown_bytes()             # measured before the database is held
        with self.db.transaction():
            row = self.db.get_store_genome(accession)        # it may have changed meanwhile
            rows = self.db.store_genomes()
            protect = self.needed_now()
            locked = self._locked_bytes(rows, protect, accession) + stray
            if need > self.budget - locked:
                raise StoreProblem("no_room",
                                   "This server has no room for the genome of %s at the moment "
                                   "(%d Mb). Try again later, or choose from the genomes it "
                                   "keeps." % (report["species"], round(need / 1e6)), status=503)
            if row is not None and row["state"] == READY:
                self._changed()                  # it leaves the list until it is back
            new_request = row is None or row["state"] in (EVICTED, FAILED)
            if origin == "visitor" and new_request:
                since = now - 86400
                if self.db.count_events("fetch", client, since) >= self.settings.max_fetches_per_day:
                    raise StoreProblem(
                        "too_many_fetches",
                        "You have had %d genomes fetched in the last 24 hours, which is the "
                        "limit. Choose from the genomes on the server, or try again tomorrow."
                        % self.settings.max_fetches_per_day, status=429)
                self.db.add_event("fetch", client, now)
            fields = {
                "state": WANTED, "species": report["species"], "stem": file_stem(report["species"]),
                "common_name": report["common_name"], "assembly_name": report["assembly_name"],
                "level": report["level"], "taxon_id": report["taxon_id"],
                "size_bp": report["size_bp"], "annotated": report["annotated"],
                "with_gff": bool(with_gff or (row and row["with_gff"])),
                "pinned": bool(pinned or (row and row["pinned"])),
                "attempts": 0, "retry_at": 0, "error_code": None, "error_text": None,
                "done_bytes": 0, "total_bytes": 0, "owner": "", "updated_at": now,
            }
            if row is None:
                self.db.insert_store_genome(dict(
                    fields, accession=accession, id=store_id(accession), added_at=now,
                    origin="operator" if origin == "operator" else "visitor"))
            else:
                if origin == "operator":
                    fields["origin"] = "operator"
                self.db.update_store_genome(accession, **fields)
            return self.db.get_store_genome(accession)

    def readiness(self, request):
        """Can a stored search start? ``("ready", [])``, ``("waiting", [views])`` or
        ``("failed", message)``. A genome that was removed while the search waited is
        asked for again."""
        waiting = []
        for ref in store_refs(request):
            need_gff = bool(ref.get("with_gff"))
            row = self.db.get_store_genome(ref["accession"])
            if self.usable(row, need_gff):
                continue
            if row is None or row["state"] in (EVICTED, READY) \
                    or (row["state"] in (WANTED, FETCHING) and need_gff and not row["with_gff"]):
                try:
                    row = self.want(ref["accession"], need_gff, origin="search")
                except StoreProblem as exc:
                    return "failed", exc.message
            if row["state"] == FAILED:
                return "failed", ("The genome of %s could not be fetched from NCBI (%s)."
                                  % (ref.get("species") or row["species"],
                                     row["error_text"] or "no reason recorded"))
            waiting.append(view(row))
        return ("waiting", waiting) if waiting else ("ready", [])

    def waiting(self, request):
        """The genomes a stored search still waits for. Changes nothing."""
        out = []
        for ref in store_refs(request):
            row = self.db.get_store_genome(ref["accession"])
            if self.usable(row, bool(ref.get("with_gff"))):
                continue
            out.append(view(row) if row else {
                "accession": ref["accession"], "species": ref.get("species") or "",
                "state": WANTED, "done_bytes": 0, "total_bytes": 0})
        return out

    def settle(self, request, library):
        """The request with every store genome described as it now is on the server
        (what it is could not be known before it was fetched)."""
        request = dict(request or {})
        changed = False

        def current(ref):
            genome = library.get(ref.get("id", "")) if ref.get("store") else None
            if genome is None:
                return ref
            fresh = dict(ref, identity=genome.identity(), size_bp=genome.size_bp,
                         common_name=genome.common_name, stem=genome.stem)
            return fresh

        home = current(request.get("home") or {})
        targets = [current(ref) for ref in request.get("targets") or []]
        if home != request.get("home") or targets != list(request.get("targets") or []):
            changed = True
            request["home"], request["targets"] = home, targets
        return request, changed

    # ── making room ─────────────────────────────────────────────────────────
    def eviction_order(self, rows, protect):
        """Ready genomes that may be removed, the first to go first: fewest uses in
        the scoring window, then longest unused, then longest on the server."""
        counts = self.use_counts()
        candidates = [r for r in rows if r["state"] == READY and not r["pinned"]
                      and r["accession"] not in protect]
        return sorted(candidates, key=lambda r: (counts.get(r["id"], (0, 0))[0],
                                                 counts.get(r["id"], (0, 0))[1],
                                                 r["ready_at"] or 0, r["accession"]))

    def evict(self, row, reason="room is needed"):
        remove_tree(self.root / row["accession"], self.root)
        self.db.update_store_genome(row["accession"], state=EVICTED, disk_bytes=0,
                                    updated_at=self._now(), owner="")
        self._changed()
        self.log("genome store: removed %s (%s, %d MB): %s"
                 % (row["accession"], row["species"], row["disk_bytes"] // 1000000, reason))

    def make_room(self, need, keep):
        """Free ``need`` bytes within the store's size. Raises StoreProblem when that
        cannot be done without touching a genome that is pinned or needed."""
        stray = self.unknown_bytes()             # takes room, and is not the store's to delete
        with self.db.transaction():
            rows = self.db.store_genomes(READY)
            used = stray + sum(r["disk_bytes"] for r in rows if r["accession"] != keep)
            order = self.eviction_order(rows, self.needed_now() | {keep})
            while used + need > self.budget:
                if not order:
                    raise StoreProblem("no_room", "the genome store is full of genomes that "
                                       "are pinned or needed by searches"
                                       + (", and of files that are in no record" if stray else ""),
                                       status=503)
                victim = order.pop(0)
                self.evict(victim)
                used -= victim["disk_bytes"]
        free = self.disk_free(self.root) - self.settings.min_free_gb * 1e9
        if free < need * 1.4:                   # packed and unpacked lie side by side for a while
            raise StoreProblem("low_disk", "the disk of the genome store is nearly full",
                               status=503)

    # ── fetching ────────────────────────────────────────────────────────────
    def next_wanted(self):
        """The genome to fetch next: one that a search waits for, else one the
        operator asked for. None when there is nothing to do right now."""
        now = self._now()
        needed = self.needed_now()
        rows = [r for r in self.db.store_genomes(WANTED) if r["retry_at"] <= now]
        for row in rows:
            if row["accession"] in needed:
                return row
        for row in rows:
            if row["origin"] == "operator":
                return row
        return None

    def fetch(self, accession, stopped=None, owner=""):
        """Download, check and unpack one genome that is wanted. Every outcome is
        written to its row; nothing is raised. Returns the row."""
        with self.db.transaction():              # one fetcher per genome, across processes
            row = self.db.get_store_genome(accession)
            if row is None or row["state"] != WANTED:
                return row
            attempts = row["attempts"] + 1
            self.db.update_store_genome(accession, state=FETCHING, owner=owner, attempts=attempts,
                                        done_bytes=0, total_bytes=0, updated_at=self._now())
        deadline = time.monotonic() + self.settings.fetch_timeout_minutes * 60
        timed_out = [False]
        beat = [time.monotonic()]

        def halted():
            if time.monotonic() > deadline:
                timed_out[0] = True
                return True
            if time.monotonic() - beat[0] >= 5.0:    # still at work (also while unpacking)
                beat[0] = time.monotonic()
                self.db.update_store_genome(accession, updated_at=self._now())
            return bool(stopped is not None and stopped())

        self.log("genome store: fetching %s (%s, %d Mb)%s" % (
            accession, row["species"], row["size_bp"] // 1000000,
            " with its annotation" if row["with_gff"] else ""))
        try:
            self._fetch_files(row, halted)
            self.log("genome store: %s is ready" % accession)
        except NcbiError as exc:
            now = self._now()
            if exc.kind == "stopped" and not timed_out[0]:
                # The worker is stopping: what has arrived stays, the next worker goes on.
                self.db.update_store_genome(accession, state=WANTED, attempts=attempts - 1,
                                            owner="", updated_at=now)
            elif exc.kind == "disk" and attempts < ROOM_TRIES:
                self.db.update_store_genome(accession, state=WANTED, owner="", updated_at=now,
                                            retry_at=now + ROOM_RETRY, error_code="low_disk",
                                            error_text="the disk of the server is full")
                self.log("genome store: PAUSED fetching %s: %s" % (accession, exc))
            elif (exc.kind == "unavailable" or timed_out[0]) and attempts < MAX_ATTEMPTS:
                wait = 60 * 4 ** attempts
                self.db.update_store_genome(accession, state=WANTED, owner="", updated_at=now,
                                            retry_at=now + wait, error_code=exc.kind,
                                            error_text=self._reason(exc, timed_out[0]))
                self.log("genome store: %s not fetched (%s); next try in %d min"
                         % (accession, exc, wait // 60))
            else:
                self._failed(accession, exc.kind, self._reason(exc, timed_out[0]))
        except StoreProblem as exc:
            now = self._now()
            if attempts < ROOM_TRIES:
                self.db.update_store_genome(accession, state=WANTED, owner="", updated_at=now,
                                            retry_at=now + ROOM_RETRY, error_code=exc.code,
                                            error_text=exc.message)
            else:
                self._failed(accession, exc.code, "the server had no room for it")
        except Exception:
            self.log("genome store: error while fetching %s:\n%s" % (accession, traceback.format_exc()))
            self._failed(accession, "internal_error", "an error on the server")
        return self.db.get_store_genome(accession)

    @staticmethod
    def _reason(exc, timed_out):
        if timed_out:
            return "the download took too long"
        return {"not_found": "NCBI no longer has its files",
                "disk": "the server had no room for it",
                "checksum": "the files arrived damaged",
                "bad_answer": "NCBI sent something unexpected",
                "unavailable": "NCBI did not answer"}.get(exc.kind, exc.kind)

    def _failed(self, accession, code, text):
        self.db.update_store_genome(accession, state=FAILED, owner="", error_code=code,
                                    error_text=text, updated_at=self._now())
        remove_tree(self.root / INCOMING / accession, self.root)
        self.log("genome store: %s FAILED: %s" % (accession, text))

    def _fetch_files(self, row, halted):
        accession, stem = row["accession"], row["stem"]
        self.root.mkdir(parents=True, exist_ok=True)
        incoming = self.root / INCOMING / accession
        final = self.root / accession
        report = self.ncbi.assembly(accession)
        files = self.ncbi.files(accession)
        if row["with_gff"] and not files["gff"]:
            raise NcbiError("not_found", "no annotation file for %s" % accession)
        self.make_room(needed_bytes(report["size_bp"], row["with_gff"]), accession)
        incoming.mkdir(parents=True, exist_ok=True)

        last = [0.0]

        def progress(done, total):
            if time.monotonic() - last[0] >= 1.0:
                last[0] = time.monotonic()
                self.db.update_store_genome(accession, done_bytes=done, total_bytes=total or 0,
                                            updated_at=self._now())

        wanted = [("fna", files["fasta"], b">")]
        if row["with_gff"]:
            wanted.append(("gff", files["gff"], b"#"))
        for suffix, name, first in wanted:
            target = incoming / ("%s.%s" % (stem, suffix))
            kept = final / target.name
            if target.is_file():
                continue
            if kept.is_file() and suffix == "fna" and row["ready_at"]:
                # The sequence is here already: this genome was fetched and checked
                # before, and only its annotation is new. A file of the same name that
                # the store did not put there itself (no earlier record) is not trusted.
                os.replace(str(kept), str(target))
                continue
            packed = incoming / ("genomic.%s.gz" % suffix)
            self.ncbi.download(files["base"] + quote(name, safe=""), packed, files["md5"][name],
                               progress=progress, stopped=halted)
            unpack(packed, target, first, stopped=halted)
            packed.unlink()

        total, n_sequences, n50 = fasta_stats(str(incoming / (stem + ".fna")))
        if total < 1:
            raise NcbiError("bad_answer", "the genome file holds no sequence")
        taxon_id, lineage, common = report["taxon_id"], [], report["common_name"]
        try:
            taxon_id, lineage, found = self.taxonomy(report["species"], taxon_id=taxon_id)
            common = common or found
        except (FetchError, ValueError) as exc:
            self.log("genome store: no taxonomy for %s (%s); it is listed last" % (accession, exc))
        disk = sum(p.stat().st_size for p in incoming.iterdir() if p.is_file())
        with self.db.transaction():
            remove_tree(final, self.root)
            os.replace(str(incoming), str(final))
            self.db.update_store_genome(
                accession, state=READY, ready_at=self._now(), updated_at=self._now(), owner="",
                size_bp=total, n_sequences=n_sequences, n50=n50, disk_bytes=disk,
                taxon_id=taxon_id, lineage=lineage, common_name=common or "",
                gff_md5=files["md5"].get(files["gff"] or "", "") if row["with_gff"] else "",
                error_code=None, error_text=None, retry_at=0)
            self._changed()

    # ── keeping it tidy ─────────────────────────────────────────────────────
    def recover(self, alive=lambda owner: False):
        """After a restart: a download whose owner is gone is wanted again. What has
        arrived of it stays and is continued."""
        for row in self.db.store_genomes(FETCHING):
            if not alive(row["owner"]):
                self.db.update_store_genome(row["accession"], state=WANTED, owner="",
                                            attempts=max(0, row["attempts"] - 1),
                                            updated_at=self._now())
                self.log("genome store: the download of %s was interrupted; it will go on"
                         % row["accession"])

    def housekeeping(self):
        now = self._now()
        self.db.prune_genome_uses(now - 2 * self.settings.store_score_days * 86400)
        needed = self.needed_now()
        keep_incoming = set()
        for row in self.db.store_genomes():
            state, accession = row["state"], row["accession"]
            if state == FETCHING and now - row["updated_at"] > STALE_FETCH:
                # Whoever fetched it is gone (a command that was killed): it goes on later.
                self.db.update_store_genome(accession, state=WANTED, owner="", updated_at=now,
                                            attempts=max(0, row["attempts"] - 1))
                keep_incoming.add(accession)
            elif state in (WANTED, FETCHING):
                idle = now - row["updated_at"] > UNWANTED_AFTER
                if state == WANTED and row["origin"] != "operator" and accession not in needed and idle:
                    # Asked for by a search that has been cancelled since.
                    self.db.delete_store_genome(accession)
                    remove_tree(self.root / accession, self.root)
                else:
                    keep_incoming.add(accession)
            elif state == FAILED and now - row["updated_at"] > 4 * FAILED_FOR:
                self.db.delete_store_genome(accession)
            elif state == READY and not self._present(row):
                # The files are gone (a disk was swapped): say so instead of failing searches.
                self.db.update_store_genome(accession, state=EVICTED, disk_bytes=0, updated_at=now)
                self._changed()
                self.log("genome store: the files of %s are missing; it is no longer listed"
                         % accession)
        incoming = self.root / INCOMING
        if incoming.is_dir():
            for path in incoming.iterdir():
                if path.name not in keep_incoming:
                    remove_tree(path, self.root)
        unknown = self.unknown()
        names = tuple(name for name, _ in unknown)
        if names != self._unknown_told:          # once, not at every round
            self._unknown_told = names
            if names:
                self.log("genome store: NOT IN ANY RECORD: %d item(s), %d MB, in %s (%s). They "
                         "count against the store's size and are never deleted automatically; "
                         "`library sweep` lists them, `library sweep --yes` deletes them."
                         % (len(names), sum(size for _, size in unknown) // 1000000, self.root,
                            ", ".join(names[:6]) + (", ..." if len(names) > 6 else "")))
        if self.enabled:
            rows = self.db.store_genomes(READY)
            used = sum(r["disk_bytes"] for r in rows)
            for victim in self.eviction_order(rows, needed):
                if used <= self.budget:
                    break
                self.evict(victim, "the store is over its size")
                used -= victim["disk_bytes"]
