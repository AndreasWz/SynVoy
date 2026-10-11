"""Shared fixtures for the tests of the web service.

No test touches the network or the real pipeline (``test_end_to_end.py`` does, and only
when asked). Searches are run with ``fake_pipeline.py``; accession look-ups answer from
canned UniProt / NCBI replies.
"""

import gzip
import hashlib
import io
import json
import os
import re
import sys
import time
from pathlib import Path

import pytest

WEB_DIR = Path(__file__).resolve().parents[1]
REPO_DIR = WEB_DIR.parent
if str(WEB_DIR) not in sys.path:
    sys.path.insert(0, str(WEB_DIR))

from synvoy_web import procs  # noqa: E402
from synvoy_web.config import load_settings  # noqa: E402
from synvoy_web.db import Store  # noqa: E402
from synvoy_web.genomes import GenomeStore  # noqa: E402
from synvoy_web.jobs import Context  # noqa: E402
from synvoy_web.library import LibraryWatcher  # noqa: E402
from synvoy_web.ncbi import Ncbi  # noqa: E402
from synvoy_web.net import FetchError  # noqa: E402
from synvoy_web.query import Resolver  # noqa: E402
from synvoy_web.runner import Worker  # noqa: E402

FAKE_PIPELINE = Path(__file__).with_name("fake_pipeline.py")
PROTEIN = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"

UNIPROT_JSON = {
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "primaryAccession": "P01501",
    "organism": {"scientificName": "Homus primus (strain T1)", "commonName": "Test animal",
                 "taxonId": 99010},
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Melittin"}}},
    "genes": [{"geneName": {"value": "MELT"}}],
    "sequence": {"value": PROTEIN, "length": len(PROTEIN)},
}
UNIPROT_FASTA = (">sp|P01501|MEL_TEST Melittin OS=Homus primus OX=99010 GN=MELT PE=1 SV=1\n"
                 + PROTEIN[:60] + "\n" + PROTEIN[60:] + "\n")


class FakeFetch:
    """Stands in for ``net.http_get``. ``replies`` maps a URL fragment to bytes or to an
    exception to raise. Every URL asked for is recorded in ``calls``."""

    def __init__(self, replies=None):
        self.replies = dict(replies or {})
        self.calls = []

    def __call__(self, url, timeout=15, max_bytes=2000000, accept="*/*"):
        self.calls.append(url)
        for fragment, reply in self.replies.items():
            if fragment in url:
                if isinstance(reply, Exception):
                    raise reply
                return reply
        raise FetchError("not_found", "no canned reply for %s" % url)


class _Stream:
    """A response of ``net.http_stream`` over bytes, which can break off part-way."""

    def __init__(self, data, partial, fail_after=None, chunk=None, delay=0.0):
        self.data, self.position, self.partial = data, 0, partial
        self.fail_after, self.chunk, self.delay = fail_after, chunk, delay
        self.status = 206 if partial else 200

    def read(self, size):
        if self.delay:
            time.sleep(self.delay)
        if self.fail_after is not None and self.position >= self.fail_after:
            raise ConnectionResetError("connection reset by the test")
        end = self.position + (min(size, self.chunk) if self.chunk else size)
        if self.fail_after is not None:
            end = min(end, self.fail_after)
        block = self.data[self.position:end]
        self.position += len(block)
        return block

    def close(self):
        pass


class FakeNcbi:
    """NCBI as the genome store sees it: genomes made up here, served through stand-ins
    for ``net.http_get`` and ``net.http_stream``. It can be down, drop a connection
    part-way, ignore the request to resume, or send a damaged file."""

    def __init__(self):
        self.assemblies = {}
        self.calls = []                 # every address asked for, with the offset of a stream
        self.down = False               # nothing answers
        self.drops = {}                 # file name -> [bytes after which the connection breaks]
        self.damaged = {}               # file name -> how many times it arrives damaged
        self.ignore_range = False       # the server sends the whole file whatever is asked
        self.chunk = None               # bytes per read (None: as many as asked for)
        self.delay = 0.0                # seconds per read: a slow connection
        self.lineages = {}              # taxon id -> lineage

    def add(self, accession, species, size=4000, annotated=False, taxon=None, name="ASM1",
            level="Chromosome", category="reference genome", lineage=None, reported=None):
        """``reported``: the size NCBI states, when it should differ from the file made here."""
        half = size // 2
        fasta = (">chr1 %s\n%s\n>chr2\n%s\n" % (species, "ACGT" * (half // 4),
                                                   "TTGA" * ((size - half) // 4))).encode()
        gff = b"##gff-version 3\nchr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=gene-A\n" if annotated else None
        folder = "%s_%s" % (accession, name)
        taxon = taxon or (9000 + len(self.assemblies))
        self.lineages[taxon] = list(lineage or [1, 2, taxon])
        self.assemblies[accession] = {
            "folder": folder, "fasta": gzip.compress(fasta, mtime=0), "plain": fasta,
            "gff": gzip.compress(gff, mtime=0) if gff else None, "plain_gff": gff,
            "report": {
                "accession": accession, "current_accession": accession,
                "organism": {"tax_id": taxon, "organism_name": species},
                "assembly_info": {"assembly_level": level, "assembly_status": "current",
                                  "assembly_name": name, "refseq_category": category},
                "assembly_stats": {"total_sequence_length": str(
                    reported or (4 * (half // 4) + 4 * ((size - half) // 4))),
                                   "contig_n50": half, "scaffold_n50": half,
                                   "number_of_scaffolds": 2},
                **({"annotation_info": {"name": "Annotation 1"}} if annotated else {}),
            },
        }
        return self.assemblies[accession]

    def _files(self, entry):
        files = {entry["folder"] + "_genomic.fna.gz": entry["fasta"]}
        if entry["gff"]:
            files[entry["folder"] + "_genomic.gff.gz"] = entry["gff"]
        return files

    def fetch(self, url, timeout=15, max_bytes=2000000, accept="*/*"):
        self.calls.append(url)
        if self.down:
            raise FetchError("unavailable", "HTTP 503")
        match = re.search(r"/genome/accession/([^/]+)/dataset_report", url)
        if match:
            entry = self.assemblies.get(match.group(1))
            return json.dumps({"reports": [entry["report"]] if entry else [],
                               "total_count": 1 if entry else 0}).encode()
        match = re.search(r"/genome/taxon/([^/]+)/dataset_report", url)
        if match:
            term = match.group(1).replace("%20", " ").lower()
            found = [e["report"] for e in self.assemblies.values()
                     if e["report"]["organism"]["organism_name"].lower().startswith(term)
                     or str(e["report"]["organism"]["tax_id"]) == term]
            if "reference_only=true" in url:
                found = [r for r in found if r["assembly_info"]["refseq_category"]]
            return json.dumps({"reports": found, "total_count": len(found)}).encode()
        for accession, entry in self.assemblies.items():
            if url.endswith("/" + entry["folder"] + "/md5checksums.txt"):
                return "".join("%s  ./%s\n" % (hashlib.md5(data).hexdigest(), name)
                               for name, data in sorted(self._files(entry).items())).encode()
            digits = accession.split("_")[1].split(".")[0]
            if url.endswith("/%s/%s/%s/%s/" % (accession[:3], digits[:3], digits[3:6], digits[6:9])):
                return ('<a href="/genomes/all/">Parent</a> <a href="%s/">%s/</a>'
                        % (entry["folder"], entry["folder"])).encode()
        raise FetchError("not_found", "HTTP 404")

    def stream(self, url, offset=0, timeout=60):
        self.calls.append("%s @%d" % (url, offset))
        if self.down:
            raise FetchError("unavailable", "HTTP 503")
        name = url.rsplit("/", 1)[1]
        for entry in self.assemblies.values():
            data = self._files(entry).get(name)
            if data is None:
                continue
            if self.damaged.get(name, 0) > 0:
                self.damaged[name] -= 1
                data = data[:-8] + bytes(8)
            if offset >= len(data) > 0:
                raise FetchError("complete", "HTTP 416")
            partial = bool(offset) and not self.ignore_range
            start = offset if partial else 0
            drops = self.drops.get(name) or []
            fail_after = drops.pop(0) if drops else None
            return _Stream(data[start:], partial, fail_after, self.chunk, self.delay), len(data)
        raise FetchError("not_found", "HTTP 404")

    def taxonomy(self, species, taxon_id=None, fetch=None):
        if self.down or taxon_id not in self.lineages:
            raise FetchError("unavailable", "no taxonomy in this test")
        return taxon_id, list(self.lineages[taxon_id]), ""

    def client(self):
        return Ncbi(fetch=self.fetch, stream=self.stream, sleep=lambda seconds: None)


def uniprot_replies(accession="P01501", data=None, fasta=None):
    data = dict(UNIPROT_JSON if data is None else data)
    return {"%s.json" % accession: json.dumps(data).encode("utf-8"),
            "%s.fasta" % accession: (UNIPROT_FASTA if fasta is None else fasta).encode("utf-8")}


def library_data(root):
    """A small library: one starting genome and five targets, files in ``root``."""
    root = Path(root)
    genomes = [
        # id, species, lineage, size_bp, with annotation
        ("homus_primus", "Homus primus", [1, 2, 3, 99010], 12000000, True),
        ("alpha_beta", "Alpha beta", [1, 2, 3, 11], 300000000, False),
        ("gamma_delta", "Gamma delta", [1, 2, 4, 12], 900000000, True),
        ("epsilon_zeta", "Epsilon zeta", [1, 5, 6, 13], 20000000, False),
        ("eta_theta", "Eta theta", [], 15000000, False),
        ("big_one", "Magnus maximus", [1, 2, 3, 14], 1500000000, False),
    ]
    entries = []
    for genome_id, species, lineage, size, annotated in genomes:
        stem = species.replace(" ", "_")
        fasta = root / (stem + ".fna")
        fasta.write_text(">chr1\nACGTACGTACGT\n>chr2\nACGT\n")
        gff = ""
        if annotated:
            gff = str(root / (stem + ".gff"))
            Path(gff).write_text("##gff-version 3\n")
        entries.append({"id": genome_id, "species": species, "common_name": "", "group": "testers",
                        "assembly": "GCF_%s.1" % genome_id, "taxon_id": lineage[-1] if lineage else None,
                        "lineage": lineage, "fasta": str(fasta), "gff": gff, "size_bp": size,
                        "n_sequences": 2, "n50": 12, "checksum": "sha256:" + genome_id})
    return {"groups": [{"id": "testers", "name": "Test animals", "note": "made up"}],
            "genomes": entries}


class Env:
    """One isolated installation: data folder, library, pages, context."""

    def __init__(self, tmp_path, **overrides):
        self.root = Path(tmp_path)
        self.genomes_dir = self.root / "genomes"
        self.genomes_dir.mkdir(exist_ok=True)
        self.site_dir = self.root / "site"
        (self.site_dir / "assets").mkdir(parents=True, exist_ok=True)
        (self.site_dir / "index.html").write_text("<!doctype html><title>home</title>")
        (self.site_dir / "404.html").write_text("<!doctype html><title>not found</title>")
        (self.site_dir / "assets" / "site.css").write_text("body{margin:0}")
        self.data_dir = self.root / "var"
        self.data_dir.mkdir(exist_ok=True)
        self.library_file = self.data_dir / "library.json"
        self.library_file.write_text(json.dumps(library_data(self.genomes_dir)))
        values = {
            "data_dir": str(self.data_dir), "site_dir": str(self.site_dir),
            "repo_dir": str(REPO_DIR), "conda_env": "", "conda_base": "",
            "pipeline_command": [sys.executable, str(FAKE_PIPELINE)],
            "min_free_gb": 0.0, "disk_factor": 0.0, "poll_seconds": 0.05,
            "kill_grace_seconds": 1.0, "max_submissions_per_day": 100,
            "max_active_per_client": 50,
        }
        values.update(overrides)
        self.settings = load_settings(env={}, overrides=values)
        self.store = Store(self.settings.db_path)
        self.fetch = FakeFetch(uniprot_replies())
        self.now = None                       # set to a number to freeze the clock
        self.ncbi = FakeNcbi()                # the genome store asks this, never the network
        self.store_disk_free = 10 ** 13
        self.genome_logs = []                 # what the store reports when no worker runs
        self.genomes = GenomeStore(self.settings, self.store, ncbi=self.ncbi.client(),
                                   clock=self.clock, taxonomy=self.ncbi.taxonomy,
                                   disk_free=lambda path: self.store_disk_free,
                                   log=self.genome_logs.append)
        self.ctx = Context(self.settings, self.store,
                           LibraryWatcher(self.settings.library_path, check_every=0.0),
                           resolver=Resolver(self.store, fetch=self.fetch, clock=self.clock),
                           clock=self.clock, genomes=self.genomes)
        self.logs = []
        self.workers = []

    def clock(self):
        return time.time() if self.now is None else self.now

    def worker(self, **kwargs):
        kwargs.setdefault("log", self.logs.append)
        kwargs.setdefault("disk_free", lambda path: 10 ** 13)
        worker = Worker(self.ctx, **kwargs)
        worker.test_logs = self.logs
        self.workers.append(worker)
        return worker

    def beat(self):
        """Pretend a healthy worker is running (for tests without one)."""
        self.store.set_meta("worker_heartbeat", int(self.clock()))
        self.store.set_meta("worker_preflight", {"ok": True, "problems": []})

    def payload(self, name="ok", targets=("alpha_beta", "epsilon_zeta", "eta_theta"),
                home="homus_primus", **extra):
        body = {"query": {"sequence": PROTEIN, "name": name, "home": home},
                "targets": list(targets)}
        body.update(extra)
        return body

    def job_ids(self):
        return [row[0] for row in self.store._conn().execute("SELECT id FROM jobs").fetchall()]

    def close(self):
        """Stop anything a test left running, so one failing test cannot disturb the rest."""
        for job_id in self.job_ids():
            for worker in self.workers:
                handle = worker._handles.get(job_id)
                if handle is not None:
                    procs.terminate(job_id, 0.5, handle)
            procs.terminate(job_id, 0.5)
        self.store.close()


@pytest.fixture
def env(tmp_path):
    environment = Env(tmp_path)
    yield environment
    environment.close()


@pytest.fixture
def make_env(tmp_path):
    """Factory for installations with their own settings: ``make_env(max_queue=1)``."""
    made = []

    def factory(**overrides):
        folder = tmp_path / ("env%d" % len(made))
        folder.mkdir()
        environment = Env(folder, **overrides)
        made.append(environment)
        return environment

    yield factory
    for environment in made:
        environment.close()


def run_until(worker, predicate, timeout=30.0, pause=0.05):
    """Tick the worker until ``predicate()`` is true."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        worker.tick()
        if predicate():
            return
        time.sleep(pause)
    raise AssertionError("not reached within %.0f s; worker log:\n%s"
                         % (timeout, "\n".join(getattr(worker, "test_logs", [])[-30:])))


def call(app, method, path, body=None, headers=None, query="", remote="10.0.0.1"):
    """One request against the WSGI app. Returns ``(status, headers, body bytes)``."""
    environ = {
        "REQUEST_METHOD": method, "PATH_INFO": path, "QUERY_STRING": query,
        "SERVER_NAME": "test", "SERVER_PORT": "80", "SERVER_PROTOCOL": "HTTP/1.1",
        "REMOTE_ADDR": remote, "wsgi.version": (1, 0), "wsgi.url_scheme": "http",
        "wsgi.errors": sys.stderr, "wsgi.multithread": True, "wsgi.multiprocess": False,
        "wsgi.run_once": False,
    }
    raw = b""
    if body is not None:
        raw = body if isinstance(body, bytes) else json.dumps(body).encode("utf-8")
        environ["CONTENT_TYPE"] = "application/json"
        environ["CONTENT_LENGTH"] = str(len(raw))
    environ["wsgi.input"] = io.BytesIO(raw)
    for key, value in (headers or {}).items():
        if key.lower() == "content-type":
            environ["CONTENT_TYPE"] = value
        elif key.lower() == "content-length":
            if value is None:
                environ.pop("CONTENT_LENGTH", None)
            else:
                environ["CONTENT_LENGTH"] = value
        else:
            environ["HTTP_" + key.upper().replace("-", "_")] = value
    captured = {}

    def start_response(status, response_headers, exc_info=None):
        captured["status"] = int(status.split()[0])
        captured["headers"] = dict(response_headers)

    chunks = app(environ, start_response)
    data = b"".join(chunks)
    if hasattr(chunks, "close"):
        chunks.close()
    return captured["status"], captured["headers"], data


def call_json(app, method, path, body=None, **kwargs):
    status, headers, data = call(app, method, path, body, **kwargs)
    return status, json.loads(data.decode("utf-8")) if data else None
