"""SQLite store: searches, who asked for them, rate-limit events, look-up cache, and
the catalogue of genomes fetched from NCBI.

One file in write-ahead mode. The API (several threads) and the worker (another
process) share it. Every change that reads before it writes runs inside
``BEGIN IMMEDIATE``, so two writers cannot interleave and a limit cannot be passed
by two requests arriving together.

Keep the file on a local disk: write-ahead logging does not work on network shares.
"""

import json
import sqlite3
import threading
from contextlib import contextmanager
from pathlib import Path

SCHEMA_VERSION = 2

# A search is in exactly one of these states.
QUEUED, RUNNING, DONE, FAILED, CANCELLED, EXPIRED = (
    "queued", "running", "done", "failed", "cancelled", "expired")
ACTIVE_STATES = (QUEUED, RUNNING)
FINAL_STATES = (DONE, FAILED, CANCELLED, EXPIRED)
ALL_STATES = ACTIVE_STATES + FINAL_STATES

_SCHEMA_V1 = """
CREATE TABLE jobs (
    id               TEXT PRIMARY KEY,
    cache_key        TEXT NOT NULL,
    state            TEXT NOT NULL,
    created_at       INTEGER NOT NULL,
    queued_at        INTEGER NOT NULL,
    started_at       INTEGER,
    finished_at      INTEGER,
    expires_at       INTEGER,
    attempt          INTEGER NOT NULL DEFAULT 0,
    cancel_requested INTEGER NOT NULL DEFAULT 0,
    request          TEXT NOT NULL,
    version          TEXT NOT NULL,
    need_mb          INTEGER NOT NULL DEFAULT 0,
    progress         TEXT,
    error_code       TEXT,
    error_text       TEXT,
    runtime_s        INTEGER
);
CREATE INDEX jobs_state ON jobs(state, queued_at);
CREATE INDEX jobs_cache ON jobs(cache_key);

CREATE TABLE submitters (
    token_hash TEXT PRIMARY KEY,
    job_id     TEXT NOT NULL REFERENCES jobs(id) ON DELETE CASCADE,
    client     TEXT NOT NULL,
    email      TEXT,
    created_at INTEGER NOT NULL,
    withdrawn  INTEGER NOT NULL DEFAULT 0
);
CREATE INDEX submitters_job ON submitters(job_id);
CREATE INDEX submitters_client ON submitters(client);

CREATE TABLE events (
    kind   TEXT NOT NULL,
    client TEXT NOT NULL,
    at     INTEGER NOT NULL
);
CREATE INDEX events_lookup ON events(kind, client, at);

CREATE TABLE lookups (
    accession  TEXT PRIMARY KEY,
    fetched_at INTEGER NOT NULL,
    payload    TEXT NOT NULL
);

CREATE TABLE meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""

# Genomes fetched from NCBI into the genome store (genomes.py), and every use of a
# genome by a search (what decides which genome goes first when space is needed).
_SCHEMA_V2 = """
CREATE TABLE store_genomes (
    accession     TEXT PRIMARY KEY,
    id            TEXT NOT NULL UNIQUE,
    state         TEXT NOT NULL,
    species       TEXT NOT NULL,
    stem          TEXT NOT NULL,
    common_name   TEXT NOT NULL DEFAULT '',
    assembly_name TEXT NOT NULL DEFAULT '',
    level         TEXT NOT NULL DEFAULT '',
    taxon_id      INTEGER,
    lineage       TEXT NOT NULL DEFAULT '[]',
    size_bp       INTEGER NOT NULL,
    n_sequences   INTEGER NOT NULL DEFAULT 0,
    n50           INTEGER NOT NULL DEFAULT 0,
    annotated     INTEGER NOT NULL DEFAULT 0,
    with_gff      INTEGER NOT NULL DEFAULT 0,
    disk_bytes    INTEGER NOT NULL DEFAULT 0,
    gff_md5       TEXT NOT NULL DEFAULT '',
    pinned        INTEGER NOT NULL DEFAULT 0,
    origin        TEXT NOT NULL DEFAULT 'visitor',
    added_at      INTEGER NOT NULL,
    ready_at      INTEGER,
    attempts      INTEGER NOT NULL DEFAULT 0,
    retry_at      INTEGER NOT NULL DEFAULT 0,
    error_code    TEXT,
    error_text    TEXT,
    done_bytes    INTEGER NOT NULL DEFAULT 0,
    total_bytes   INTEGER NOT NULL DEFAULT 0,
    owner         TEXT NOT NULL DEFAULT '',
    updated_at    INTEGER NOT NULL
);
CREATE INDEX store_genomes_state ON store_genomes(state, added_at);

CREATE TABLE genome_uses (
    genome_id TEXT NOT NULL,
    at        INTEGER NOT NULL
);
CREATE INDEX genome_uses_genome ON genome_uses(genome_id, at);
"""

_GENOME_COLUMNS = ("accession", "id", "state", "species", "stem", "common_name", "assembly_name",
                   "level", "taxon_id", "lineage", "size_bp", "n_sequences", "n50", "annotated",
                   "with_gff", "disk_bytes", "gff_md5", "pinned", "origin", "added_at", "ready_at",
                   "attempts", "retry_at", "error_code", "error_text", "done_bytes",
                   "total_bytes", "owner", "updated_at")
_GENOME_FLAGS = ("annotated", "with_gff", "pinned")

_JOB_COLUMNS = ("id", "cache_key", "state", "created_at", "queued_at", "started_at",
                "finished_at", "expires_at", "attempt", "cancel_requested", "request",
                "version", "need_mb", "progress", "error_code", "error_text", "runtime_s")
_JSON_COLUMNS = ("request", "progress")
_UPDATABLE = set(_JOB_COLUMNS) - {"id", "created_at"}


class StoreError(RuntimeError):
    """The database cannot be used as it is (wrong version, damaged file)."""


class Store:
    def __init__(self, path):
        self.path = str(path)
        self._local = threading.local()
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._migrate()

    # ── connections ─────────────────────────────────────────────────────────
    def _conn(self):
        conn = getattr(self._local, "conn", None)
        if conn is None:
            # isolation_level=None: no implicit transactions, we say BEGIN ourselves.
            conn = sqlite3.connect(self.path, timeout=30, isolation_level=None)
            conn.row_factory = sqlite3.Row
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA synchronous=FULL")
            conn.execute("PRAGMA foreign_keys=ON")
            self._local.conn = conn
            self._local.depth = 0
        return conn

    def close(self):
        conn = getattr(self._local, "conn", None)
        if conn is not None:
            conn.close()
            self._local.conn = None

    @contextmanager
    def transaction(self):
        """All statements inside succeed together or not at all. May be nested."""
        conn = self._conn()
        if self._local.depth:
            self._local.depth += 1
            try:
                yield conn
            finally:
                self._local.depth -= 1
            return
        conn.execute("BEGIN IMMEDIATE")
        self._local.depth = 1
        try:
            yield conn
        except BaseException:
            conn.execute("ROLLBACK")
            raise
        else:
            conn.execute("COMMIT")
        finally:
            self._local.depth = 0

    def _migrate(self):
        conn = self._conn()
        with self.transaction():
            version = conn.execute("PRAGMA user_version").fetchone()[0]
            if version > SCHEMA_VERSION:
                raise StoreError(
                    "the database %s was written by a newer version of the web service "
                    "(schema %d, this code knows %d)" % (self.path, version, SCHEMA_VERSION))
            # executescript would end the transaction; run the statements one by one.
            for number, schema in ((1, _SCHEMA_V1), (2, _SCHEMA_V2)):
                if version < number:
                    for statement in schema.split(";"):
                        if statement.strip():
                            conn.execute(statement)
            if version < SCHEMA_VERSION:
                conn.execute("PRAGMA user_version=%d" % SCHEMA_VERSION)

    def quick_check(self):
        """'ok' when SQLite finds the file intact."""
        return self._conn().execute("PRAGMA quick_check").fetchone()[0]

    # ── searches ────────────────────────────────────────────────────────────
    @staticmethod
    def _job(row):
        if row is None:
            return None
        job = dict(row)
        for key in _JSON_COLUMNS:
            if job.get(key):
                try:
                    job[key] = json.loads(job[key])
                except ValueError:
                    job[key] = None
        job["cancel_requested"] = bool(job["cancel_requested"])
        return job

    def insert_job(self, job):
        values = dict(job)
        for key in _JSON_COLUMNS:
            if values.get(key) is not None:
                values[key] = json.dumps(values[key], sort_keys=True)
        columns = [c for c in _JOB_COLUMNS if c in values]
        self._conn().execute(
            "INSERT INTO jobs (%s) VALUES (%s)" % (", ".join(columns), ", ".join("?" * len(columns))),
            [values[c] for c in columns])

    def get_job(self, job_id):
        row = self._conn().execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
        return self._job(row)

    def update_job(self, job_id, **fields):
        if not fields:
            return
        unknown = set(fields) - _UPDATABLE
        if unknown:
            raise ValueError("cannot update column(s) %s" % ", ".join(sorted(unknown)))
        values = dict(fields)
        for key in _JSON_COLUMNS:
            if values.get(key) is not None:
                values[key] = json.dumps(values[key], sort_keys=True)
        if "cancel_requested" in values:
            values["cancel_requested"] = int(bool(values["cancel_requested"]))
        columns = sorted(values)
        self._conn().execute(
            "UPDATE jobs SET %s WHERE id = ?" % ", ".join("%s = ?" % c for c in columns),
            [values[c] for c in columns] + [job_id])

    def delete_job(self, job_id):
        self._conn().execute("DELETE FROM jobs WHERE id = ?", (job_id,))

    def jobs_in_state(self, *states):
        marks = ", ".join("?" * len(states))
        rows = self._conn().execute(
            "SELECT * FROM jobs WHERE state IN (%s) ORDER BY queued_at, rowid" % marks,
            states).fetchall()
        return [self._job(r) for r in rows]

    def count_state(self, state):
        return self._conn().execute(
            "SELECT COUNT(*) FROM jobs WHERE state = ?", (state,)).fetchone()[0]

    def newest_with_key(self, cache_key):
        """The most recent search with this key, whatever its state."""
        row = self._conn().execute(
            "SELECT * FROM jobs WHERE cache_key = ? ORDER BY created_at DESC, rowid DESC LIMIT 1",
            (cache_key,)).fetchone()
        return self._job(row)

    def queue_ahead(self, job):
        """How many queued searches are started before this one."""
        return self._conn().execute(
            "SELECT COUNT(*) FROM jobs WHERE state = ? AND (queued_at < ? OR "
            "(queued_at = ? AND rowid < (SELECT rowid FROM jobs WHERE id = ?)))",
            (QUEUED, job["queued_at"], job["queued_at"], job["id"])).fetchone()[0]

    def claim(self, job_id, now):
        """Move one queued search to RUNNING and return it. None when it is not queued
        any more (cancelled a moment ago, or taken already)."""
        with self.transaction() as conn:
            cursor = conn.execute(
                "UPDATE jobs SET state = ?, attempt = attempt + 1, cancel_requested = 0, "
                "started_at = COALESCE(started_at, ?) WHERE id = ? AND state = ?",
                (RUNNING, now, job_id, QUEUED))
            if cursor.rowcount != 1:
                return None
            return self.get_job(job_id)

    # ── who asked ───────────────────────────────────────────────────────────
    def add_submitter(self, job_id, token_hash, client, email, now):
        self._conn().execute(
            "INSERT INTO submitters (token_hash, job_id, client, email, created_at) "
            "VALUES (?, ?, ?, ?, ?)", (token_hash, job_id, client, email or None, now))

    def withdraw(self, job_id, token_hash):
        """Mark one submitter as gone. Returns how many still want the search, or None
        when the token does not belong to this search."""
        with self.transaction() as conn:
            row = conn.execute(
                "SELECT withdrawn FROM submitters WHERE token_hash = ? AND job_id = ?",
                (token_hash, job_id)).fetchone()
            if row is None:
                return None
            conn.execute(
                "UPDATE submitters SET withdrawn = 1, email = NULL WHERE token_hash = ?",
                (token_hash,))
            return conn.execute(
                "SELECT COUNT(*) FROM submitters WHERE job_id = ? AND withdrawn = 0",
                (job_id,)).fetchone()[0]

    def active_for_client(self, client):
        return self._conn().execute(
            "SELECT COUNT(DISTINCT j.id) FROM jobs j JOIN submitters s ON s.job_id = j.id "
            "WHERE s.client = ? AND s.withdrawn = 0 AND j.state IN (?, ?)",
            (client, QUEUED, RUNNING)).fetchone()[0]

    def emails_for(self, job_id):
        rows = self._conn().execute(
            "SELECT DISTINCT email FROM submitters WHERE job_id = ? AND withdrawn = 0 "
            "AND email IS NOT NULL", (job_id,)).fetchall()
        return [r[0] for r in rows]

    def forget_emails(self, job_id):
        self._conn().execute("UPDATE submitters SET email = NULL WHERE job_id = ?", (job_id,))

    def forget_visitors(self, job_id):
        """Drop what links a finished search to who asked for it."""
        self._conn().execute(
            "UPDATE submitters SET client = '', email = NULL WHERE job_id = ?", (job_id,))

    # ── rate limits ─────────────────────────────────────────────────────────
    def add_event(self, kind, client, now):
        self._conn().execute(
            "INSERT INTO events (kind, client, at) VALUES (?, ?, ?)", (kind, client, now))

    def count_events(self, kind, client, since):
        return self._conn().execute(
            "SELECT COUNT(*) FROM events WHERE kind = ? AND client = ? AND at >= ?",
            (kind, client, since)).fetchone()[0]

    def prune_events(self, before):
        self._conn().execute("DELETE FROM events WHERE at < ?", (before,))

    # ── accession look-ups ──────────────────────────────────────────────────
    def get_lookup(self, accession, not_before):
        row = self._conn().execute(
            "SELECT payload FROM lookups WHERE accession = ? AND fetched_at >= ?",
            (accession, not_before)).fetchone()
        if row is None:
            return None
        try:
            return json.loads(row[0])
        except ValueError:
            return None

    def put_lookup(self, accession, payload, now):
        self._conn().execute(
            "INSERT OR REPLACE INTO lookups (accession, fetched_at, payload) VALUES (?, ?, ?)",
            (accession, now, json.dumps(payload, sort_keys=True)))

    def prune_lookups(self, before):
        self._conn().execute("DELETE FROM lookups WHERE fetched_at < ?", (before,))

    # ── genomes fetched from NCBI ───────────────────────────────────────────
    @staticmethod
    def _genome(row):
        if row is None:
            return None
        genome = dict(row)
        try:
            genome["lineage"] = [int(x) for x in json.loads(genome.get("lineage") or "[]")]
        except (ValueError, TypeError):
            genome["lineage"] = []
        for key in _GENOME_FLAGS:
            genome[key] = bool(genome[key])
        return genome

    @staticmethod
    def _genome_values(fields):
        unknown = set(fields) - set(_GENOME_COLUMNS)
        if unknown:
            raise ValueError("unknown genome column(s) %s" % ", ".join(sorted(unknown)))
        values = dict(fields)
        if "lineage" in values:
            values["lineage"] = json.dumps(list(values["lineage"] or []))
        for key in _GENOME_FLAGS:
            if key in values:
                values[key] = int(bool(values[key]))
        return values

    def get_store_genome(self, accession):
        row = self._conn().execute(
            "SELECT * FROM store_genomes WHERE accession = ?", (accession,)).fetchone()
        return self._genome(row)

    def store_genomes(self, *states):
        """Genomes of the store in the given states (all when none is named), oldest first."""
        sql, args = "SELECT * FROM store_genomes", ()
        if states:
            sql += " WHERE state IN (%s)" % ", ".join("?" * len(states))
            args = states
        rows = self._conn().execute(sql + " ORDER BY added_at, accession", args).fetchall()
        return [self._genome(r) for r in rows]

    def insert_store_genome(self, fields):
        values = self._genome_values(fields)
        columns = sorted(values)
        self._conn().execute(
            "INSERT INTO store_genomes (%s) VALUES (%s)"
            % (", ".join(columns), ", ".join("?" * len(columns))), [values[c] for c in columns])

    def update_store_genome(self, accession, **fields):
        values = self._genome_values(fields)
        values.pop("accession", None)
        if not values:
            return 0
        columns = sorted(values)
        cursor = self._conn().execute(
            "UPDATE store_genomes SET %s WHERE accession = ?" % ", ".join("%s = ?" % c for c in columns),
            [values[c] for c in columns] + [accession])
        return cursor.rowcount

    def delete_store_genome(self, accession):
        self._conn().execute("DELETE FROM store_genomes WHERE accession = ?", (accession,))

    def add_genome_uses(self, genome_ids, now):
        self._conn().executemany("INSERT INTO genome_uses (genome_id, at) VALUES (?, ?)",
                                 [(genome_id, now) for genome_id in genome_ids])

    def genome_use_counts(self, since=0):
        """{genome id: (number of uses since ``since``, time of the last use ever)}."""
        rows = self._conn().execute(
            "SELECT genome_id, SUM(CASE WHEN at >= ? THEN 1 ELSE 0 END), MAX(at) "
            "FROM genome_uses GROUP BY genome_id", (since,)).fetchall()
        return {r[0]: (int(r[1] or 0), int(r[2] or 0)) for r in rows}

    def prune_genome_uses(self, before):
        self._conn().execute("DELETE FROM genome_uses WHERE at < ?", (before,))

    # ── small facts shared by API and worker ────────────────────────────────
    def get_meta(self, key, default=None):
        row = self._conn().execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
        if row is None:
            return default
        try:
            return json.loads(row[0])
        except ValueError:
            return default

    def set_meta(self, key, value):
        self._conn().execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)", (key, json.dumps(value)))
