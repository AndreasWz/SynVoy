"""Settings of the web service.

Values come from three places, later ones winning:

1. the defaults below,
2. a JSON file (``--config FILE`` or ``$SYNVOY_WEB_CONFIG``),
3. environment variables ``SYNVOY_WEB_<KEY>`` (lists and booleans as JSON).

An unknown key is an error. A limit that was mistyped must not stay at its default
without anyone noticing.
"""

import json
import os
import secrets
from pathlib import Path

WEB_DIR = Path(__file__).resolve().parent.parent
REPO_DIR = WEB_DIR.parent

DEFAULTS = {
    # ── where things are ────────────────────────────────────────────────────
    "data_dir": str(WEB_DIR / "var"),      # database, job folders, caches (local disk)
    "site_dir": str(WEB_DIR / "site"),     # the static pages
    "library_file": "",                    # genome library; default <data_dir>/library.json
    "repo_dir": str(REPO_DIR),             # the SynVoy checkout that is run
    # ── genomes fetched from NCBI on request (the genome store) ─────────────
    "genome_store_dir": "",                # default <data_dir>/genomes
    "genome_store_gb": 0.0,                # space the store may take; 0 = nothing is fetched
    "max_fetches_per_day": 5,              # genomes one visitor may have fetched per day
    "fetch_timeout_minutes": 180,          # one genome download is given up after this
    "store_score_days": 90,                # uses in this time decide which genome goes first
    # ── how a search is started ─────────────────────────────────────────────
    "conda_base": "",                      # e.g. /opt/miniforge3; empty = use PATH as it is
    "conda_env": "synvoy_env",
    "conda_cache_dir": "",                 # default <data_dir>/conda_cache (shared by all searches)
    "nextflow_profile": "auto,low_mem",
    "nextflow_config": "",                 # extra `-c FILE`
    "pipeline_args": [],                   # appended to every run
    "pipeline_command": [],                # replaces `nextflow run main.nf` (tests)
    "use_target_annotations": False,       # pass library GFFs of the targets as --target_gffs
    # ── http ────────────────────────────────────────────────────────────────
    "host": "127.0.0.1",
    "port": 8080,
    "trust_proxy": False,                  # read the client address from X-Forwarded-For
    "cors_origins": [],                    # other origins allowed to call the API
    "public_base_url": "",                 # used in e-mails, e.g. https://synvoy.example.org
    # ── what a visitor may ask for ──────────────────────────────────────────
    "accepting": True,                     # operator switch: False closes new searches
    "closed_message": "",
    "max_targets": 10,
    "max_genome_mb": 1000,                 # largest single genome that may be searched
    "max_total_mb": 6000,                  # sum over the selected genomes
    "min_query_length": 30,
    "max_query_length": 3000,
    "max_active_per_client": 2,            # queued + running searches per visitor
    "max_submissions_per_day": 10,
    "max_lookups_per_hour": 60,            # accession look-ups that reach UniProt / NCBI
    "max_queue": 20,
    # ── how searches are run ────────────────────────────────────────────────
    "max_running": 1,
    "max_runtime_hours": 12.0,
    "max_attempts": 2,                     # a search interrupted by a restart is resumed once
    "retention_days": 30,
    "failed_retention_days": 7,
    "min_free_gb": 20.0,                   # never start a search below this much free disk
    "disk_factor": 4.0,                    # working space per searched base pair (see README)
    "poll_seconds": 3.0,
    "kill_grace_seconds": 30.0,
    "worker_stale_seconds": 120,
    # ── e-mail (off unless smtp_host is set) ────────────────────────────────
    "smtp_host": "",
    "smtp_port": 587,
    "smtp_user": "",
    "smtp_password": "",                   # prefer $SYNVOY_WEB_SMTP_PASSWORD over the file
    "smtp_from": "",
    "smtp_starttls": True,
    # ── shown on the pages ──────────────────────────────────────────────────
    "version_label": "",                   # shown as the SynVoy version; empty = read from git (a release tag, else "pre-release")
    "operator_name": "",
    "contact_email": "",
}

_PATH_KEYS = ("data_dir", "site_dir", "library_file", "repo_dir", "conda_base",
              "conda_cache_dir", "nextflow_config", "genome_store_dir")


class ConfigError(ValueError):
    """A setting is missing, unknown or out of range."""


def _coerce(key, raw, default):
    """Turn the text of an environment variable into the type of the default."""
    if isinstance(default, bool):
        low = raw.strip().lower()
        if low in ("1", "true", "yes", "on"):
            return True
        if low in ("0", "false", "no", "off", ""):
            return False
        raise ConfigError("%s: expected true or false, got %r" % (key, raw))
    if isinstance(default, int):
        try:
            return int(raw)
        except ValueError:
            raise ConfigError("%s: expected a whole number, got %r" % (key, raw))
    if isinstance(default, float):
        try:
            return float(raw)
        except ValueError:
            raise ConfigError("%s: expected a number, got %r" % (key, raw))
    if isinstance(default, list):
        try:
            value = json.loads(raw)
        except ValueError:
            raise ConfigError("%s: expected a JSON list, got %r" % (key, raw))
        if not isinstance(value, list):
            raise ConfigError("%s: expected a JSON list, got %r" % (key, raw))
        return value
    return raw


def _check_type(key, value, default):
    if isinstance(default, bool):
        ok = isinstance(value, bool)
    elif isinstance(default, int):
        ok = isinstance(value, int) and not isinstance(value, bool)
    elif isinstance(default, float):
        ok = isinstance(value, (int, float)) and not isinstance(value, bool)
    elif isinstance(default, list):
        ok = isinstance(value, list) and all(isinstance(v, str) for v in value)
    else:
        ok = isinstance(value, str)
    if not ok:
        raise ConfigError("%s: expected %s, got %r" % (key, type(default).__name__, value))
    return float(value) if isinstance(default, float) else value


class Settings:
    """Validated settings. Read attributes, e.g. ``settings.max_targets``."""

    def __init__(self, values):
        self._values = dict(values)
        for key, value in self._values.items():
            setattr(self, key, value)

    def as_dict(self):
        return dict(self._values)

    # Derived locations -------------------------------------------------------
    @property
    def db_path(self):
        return Path(self.data_dir) / "synvoy_web.sqlite3"

    @property
    def jobs_dir(self):
        return Path(self.data_dir) / "jobs"

    @property
    def library_path(self):
        return Path(self.library_file) if self.library_file else Path(self.data_dir) / "library.json"

    @property
    def conda_cache_path(self):
        if self.conda_cache_dir:
            return Path(self.conda_cache_dir)
        return Path(self.data_dir) / "conda_cache"

    @property
    def genome_store_path(self):
        if self.genome_store_dir:
            return Path(self.genome_store_dir)
        return Path(self.data_dir) / "genomes"

    @property
    def email_enabled(self):
        return bool(self.smtp_host and self.smtp_from and self.public_base_url)

    def secret(self):
        """Salt for the hashed client keys. Made once and kept in the data folder, so
        visitor addresses are never stored and the keys survive a restart."""
        path = Path(self.data_dir) / "secret"
        try:
            text = path.read_text().strip()
            if len(text) >= 32:
                return text
        except OSError:
            pass
        text = secrets.token_hex(32)
        path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(fd, "w") as fh:
            fh.write(text + "\n")
        return text


def _validate(values):
    def need(cond, message):
        if not cond:
            raise ConfigError(message)

    need(1 <= values["port"] <= 65535, "port: must be between 1 and 65535")
    need(values["max_targets"] >= 1, "max_targets: must be at least 1")
    need(values["max_genome_mb"] >= 1, "max_genome_mb: must be at least 1")
    need(values["max_total_mb"] >= values["max_genome_mb"],
         "max_total_mb: must not be smaller than max_genome_mb")
    need(values["min_query_length"] >= 1, "min_query_length: must be at least 1")
    need(values["max_query_length"] >= values["min_query_length"],
         "max_query_length: must not be smaller than min_query_length")
    need(values["max_active_per_client"] >= 1, "max_active_per_client: must be at least 1")
    need(values["max_submissions_per_day"] >= 1, "max_submissions_per_day: must be at least 1")
    need(values["max_lookups_per_hour"] >= 1, "max_lookups_per_hour: must be at least 1")
    need(values["max_queue"] >= 1, "max_queue: must be at least 1")
    need(values["max_running"] >= 1, "max_running: must be at least 1")
    need(values["max_runtime_hours"] > 0, "max_runtime_hours: must be above 0")
    need(values["max_attempts"] >= 1, "max_attempts: must be at least 1")
    need(values["retention_days"] >= 1, "retention_days: must be at least 1")
    need(values["failed_retention_days"] >= 1, "failed_retention_days: must be at least 1")
    need(values["min_free_gb"] >= 0, "min_free_gb: must not be negative")
    need(values["genome_store_gb"] >= 0, "genome_store_gb: must not be negative")
    need(values["max_fetches_per_day"] >= 1, "max_fetches_per_day: must be at least 1")
    need(values["fetch_timeout_minutes"] >= 1, "fetch_timeout_minutes: must be at least 1")
    need(values["store_score_days"] >= 1, "store_score_days: must be at least 1")
    need(values["disk_factor"] >= 0, "disk_factor: must not be negative")
    need(values["poll_seconds"] > 0, "poll_seconds: must be above 0")
    need(values["kill_grace_seconds"] >= 0, "kill_grace_seconds: must not be negative")
    need(values["worker_stale_seconds"] >= 10, "worker_stale_seconds: must be at least 10")
    for origin in values["cors_origins"]:
        need(origin.startswith(("http://", "https://")) and not origin.endswith("/"),
             "cors_origins: %r must look like https://host (no trailing slash)" % origin)
    base = values["public_base_url"]
    need(not base or base.startswith(("http://", "https://")),
         "public_base_url: must start with http:// or https://")
    if values["smtp_host"]:
        need(values["smtp_from"], "smtp_from: needed when smtp_host is set")
        need(base, "public_base_url: needed when smtp_host is set (the e-mail links to the result)")
    # Profile names go on a command line; keep them to what Nextflow profiles look like.
    profile = values["nextflow_profile"]
    need(profile and all(part and part.replace("_", "").isalnum() for part in profile.split(",")),
         "nextflow_profile: expected names separated by commas, got %r" % profile)


def load_settings(path=None, env=None, overrides=None):
    """Build the settings. ``overrides`` (a dict) wins over file and environment; the
    tests and the command line use it."""
    env = os.environ if env is None else env
    values = dict(DEFAULTS)
    path = path or env.get("SYNVOY_WEB_CONFIG") or ""
    if path:
        try:
            with open(path) as fh:
                data = json.load(fh)
        except OSError as exc:
            raise ConfigError("cannot read the settings file %s: %s" % (path, exc))
        except ValueError as exc:
            raise ConfigError("the settings file %s is not valid JSON: %s" % (path, exc))
        if not isinstance(data, dict):
            raise ConfigError("the settings file %s must hold one JSON object" % path)
        for key, value in data.items():
            if key.startswith("_"):          # "_comment" entries are allowed
                continue
            if key not in DEFAULTS:
                raise ConfigError("unknown setting %r in %s" % (key, path))
            values[key] = _check_type(key, value, DEFAULTS[key])
    for key, default in DEFAULTS.items():
        raw = env.get("SYNVOY_WEB_" + key.upper())
        if raw is not None:
            values[key] = _check_type(key, _coerce(key, raw, default), default)
    for key, value in (overrides or {}).items():
        if key not in DEFAULTS:
            raise ConfigError("unknown setting %r" % key)
        values[key] = _check_type(key, value, DEFAULTS[key])
    for key in _PATH_KEYS:
        if values[key]:
            values[key] = str(Path(values[key]).expanduser().resolve())
    _validate(values)
    return Settings(values)
