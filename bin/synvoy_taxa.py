"""Common-name lookup for SynVoy plots — no hardcoded species list.

Resolution order:
  1. user-supplied TSV (--common_names_tsv)            — highest priority.
  2. on-disk cache (default: ~/.cache/synvoy/common_names.tsv).
  3. NCBI `datasets summary taxonomy taxon "<name>"`   — populates the cache.
  4. fall back to the scientific name unchanged.

The `datasets` CLI is already in the SynVoy env (used elsewhere by
`fetch_related_genomes.py`); we shell out to it lazily on cache miss and
parse `curator_common_name`. If `datasets` is missing or offline, lookups
silently return None and the caller falls back to the scientific name.

The cache is intentionally human-readable so users can hand-edit when NCBI
returns nothing useful (e.g. niche species without curator names).
"""

import csv
import json
import os
import re
import shutil
import subprocess
import time
from typing import Optional

DEFAULT_CACHE = os.path.join(
    os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")),
    "synvoy", "common_names.tsv",
)
DATASETS_TIMEOUT_S = 8
# A lookup NCBI did not answer is asked once more after a pause, at most RETRY_BUDGET times per
# plot run (fix of 9 Nov: a run started two minutes after boot recorded 'unreachable' where a later
# run recorded 'none'). The budget bounds what a dead network adds: RETRY_BUDGET x (pause + timeout).
RETRY_BUDGET = 2
RETRY_PAUSE_S = 2

_OVERRIDES: dict = {}
_CACHE: dict = {}
_CACHE_PATH: str = ""
_DATASETS_BIN: Optional[str] = None
_NEGATIVE: set = set()  # names we already tried and got nothing for, this run
_RETRIES_LEFT = RETRY_BUDGET
# Where each name of this run came from (fix-phase d7): table, cache, ncbi, none (NCBI answered
# without a common name), unreachable (datasets failed or ran past DATASETS_TIMEOUT_S), offline (no
# lookup: --no_network or no datasets). A figure label then says what it rests on: the per-user
# cache makes a label differ between machines, and a lookup that did not answer drops it.
SOURCES: dict = {}


def _normalise(name: str) -> str:
    if not name:
        return ""
    s = name.replace("_", " ").strip()
    parts = s.split(" ", 1)
    if not parts:
        return ""
    head = parts[0].capitalize()
    rest = parts[1].lower() if len(parts) > 1 else ""
    return (head + " " + rest).strip()


def init_lookup(common_names_tsv: Optional[str] = None,
                cache_path: Optional[str] = None,
                allow_network: bool = True) -> None:
    """Configure the lookup. Call once at plot startup.

    common_names_tsv: optional 2-column TSV `<scientific>\\t<common>` —
                     overrides everything else.
    cache_path:      where to read/write resolved names. Defaults to
                     ~/.cache/synvoy/common_names.tsv.
    allow_network:   if False, skip the NCBI lookup (cache + overrides only).
                     Use this in CI/offline runs.
    """
    global _OVERRIDES, _CACHE, _CACHE_PATH, _DATASETS_BIN, _NEGATIVE, _RETRIES_LEFT
    _OVERRIDES = {}
    _CACHE = {}
    _NEGATIVE = set()
    _RETRIES_LEFT = RETRY_BUDGET
    _CACHE_PATH = cache_path or DEFAULT_CACHE

    if common_names_tsv and os.path.exists(common_names_tsv):
        for sci, common in _read_tsv(common_names_tsv):
            _OVERRIDES[_normalise(sci)] = common

    if os.path.exists(_CACHE_PATH):
        for sci, common in _read_tsv(_CACHE_PATH):
            _CACHE[_normalise(sci)] = common

    _DATASETS_BIN = shutil.which("datasets") if allow_network else None


def common_name(scientific: str) -> Optional[str]:
    """Return the common name for `scientific`, or None if no mapping found."""
    global _RETRIES_LEFT
    norm = _normalise(scientific)
    if not norm:
        return None
    if norm in _OVERRIDES:
        SOURCES.setdefault(norm, ("table", _OVERRIDES[norm]))
        return _OVERRIDES[norm] or None
    if norm in _CACHE:
        SOURCES.setdefault(norm, ("cache", _CACHE[norm]))
        return _CACHE[norm] or None
    if norm in _NEGATIVE:
        return None
    if _DATASETS_BIN:
        cn, answered = _query_datasets_cli(norm)
        if not answered and _RETRIES_LEFT > 0:
            _RETRIES_LEFT -= 1
            time.sleep(RETRY_PAUSE_S)
            cn, answered = _query_datasets_cli(norm)
        if cn:
            _CACHE[norm] = cn
            _persist_cache_entry(norm, cn)
            SOURCES[norm] = ("ncbi", cn)
            return cn
        SOURCES[norm] = ("none" if answered else "unreachable", "")
    else:
        SOURCES[norm] = ("offline", "")
    _NEGATIVE.add(norm)
    return None


def write_sources(path: str) -> None:
    """The run's species names and where each came from (SOURCES), as a TSV beside the figures."""
    with open(path, "w", newline="") as fh:
        fh.write("# SynVoy species names: the common name each figure label carries, and its source "
                 "(table, cache = the per-user cache of an earlier lookup, ncbi, none = NCBI has none, "
                 "unreachable = NCBI did not answer, offline = not looked up).\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")   # LF like the header line
        w.writerow(["scientific_name", "common_name", "source"])
        for sci in sorted(SOURCES):
            w.writerow([sci, SOURCES[sci][1], SOURCES[sci][0]])


def label_for_species(scientific: str, mode: str = "both") -> str:
    """Render a species label per `mode`:
       - "scientific": always the scientific name
       - "common":     common name if known, else scientific
       - "both":       "Scientific (common)" if known, else just scientific
    """
    sci = _normalise(scientific) or scientific or ""
    if mode == "scientific":
        return sci
    cn = common_name(sci)
    if mode == "common":
        return cn or sci
    if cn and cn.lower() != sci.lower():
        return f"{sci} ({cn})"
    return sci


# ─────────────────────────────── helpers ───────────────────────────────────

def _read_tsv(path: str):
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or not row[0].strip() or row[0].startswith("#"):
                continue
            sci = row[0].strip()
            common = row[1].strip() if len(row) > 1 else ""
            yield sci, common


# NCBI's datasets exits 1 when it rejects a name ("The taxonomy name ... is not recognized", "not
# exact"): that is an answer, not a failed lookup (first review of fix-phase d7).
DATASETS_REJECTED = re.compile(r"not recognized|not exact|no .*found|not found", re.IGNORECASE)


def _query_datasets_cli(name: str):
    """Shell out to `datasets summary taxonomy taxon ...` and parse the
    response. Returns (the curator common name or None, whether NCBI answered)."""
    if not _DATASETS_BIN:
        return None, False
    try:
        proc = subprocess.run(
            [_DATASETS_BIN, "summary", "taxonomy", "taxon", name,
             "--as-json-lines"],
            capture_output=True, text=True,
            timeout=DATASETS_TIMEOUT_S,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None, False
    if proc.returncode != 0:
        return None, bool(DATASETS_REJECTED.search(proc.stderr or ""))
    if not proc.stdout.strip():
        return None, True
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            doc = json.loads(line)
        except json.JSONDecodeError:
            continue
        tax = doc.get("taxonomy") or {}
        cn = (tax.get("curator_common_name")
              or _first_common_name(tax.get("other_names", {})))
        if cn:
            return cn.strip(), True
    return None, True


def _first_common_name(other_names) -> Optional[str]:
    if isinstance(other_names, dict):
        for key in ("genbank_common_name", "common_name", "common_names"):
            v = other_names.get(key)
            if isinstance(v, str) and v.strip():
                return v
            if isinstance(v, list) and v:
                return v[0]
    return None


def _persist_cache_entry(scientific: str, common: str) -> None:
    """Append a single row to the on-disk cache."""
    try:
        os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
        with open(_CACHE_PATH, "a", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerow([scientific, common])
    except OSError:
        pass  # cache is best-effort; not having it is not a fatal error
