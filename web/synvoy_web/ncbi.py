"""NCBI genome assemblies: what one is, where its files are, and getting them.

Used when a visitor or the operator asks for a genome that is not on the server.
Everything is read over https from two fixed hosts:

    api.ncbi.nlm.nih.gov/datasets/v2     what an assembly is: species, size, annotation
    ftp.ncbi.nlm.nih.gov/genomes/all     its files: <name>_genomic.fna.gz, .gff.gz and
                                         md5checksums.txt

A file is fetched in a way that can be resumed after a broken connection, is compared
with NCBI's MD5 sum, and is unpacked only after that.

An assembly accession with its version (``GCF_000146045.2``) names one fixed set of
sequences for ever, which is why it can stand for the genome in the result cache.
"""

import gzip
import hashlib
import http.client
import json
import os
import re
import time
import zlib
from urllib.parse import quote

from .net import FetchError, http_get, http_stream

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
FILES = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
ACCESSION_RE = re.compile(r"^(GC[AF])_(\d{9})\.(\d{1,3})$")
_LEVELS = {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}
_OPEN_NAMES = ("sp.", "sp", "cf.", "cf", "aff.", "aff", "x")
BLOCK = 1 << 20


class NcbiError(Exception):
    """``kind``: not_found, unavailable, bad_answer, checksum, stopped, disk."""

    def __init__(self, kind, detail=""):
        Exception.__init__(self, "%s: %s" % (kind, detail) if detail else kind)
        self.kind = kind
        self.detail = detail


def normalise_accession(text):
    """'gcf_000146045.2 ' -> 'GCF_000146045.2'; None when it is not an assembly
    accession with a version."""
    if not isinstance(text, str):
        return None
    value = text.strip().upper()
    return value if ACCESSION_RE.match(value) else None


def species_name(organism):
    """A name that can be shown and turned into a file name: 'Saccharomyces cerevisiae
    S288C' -> 'Saccharomyces cerevisiae'; an open name keeps its label
    ('Probolomyrmex sp. GAGA-0580')."""
    text = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", " ", str(organism or ""))
    words = re.sub(r"[^A-Za-z0-9 .\-]+", " ", text).split()
    if len(words) >= 2 and words[1].lower() not in _OPEN_NAMES:
        words = words[:2]
    name = " ".join(words[:4])[:80].strip(" .-")
    return name if re.match(r"^[A-Za-z][A-Za-z0-9 .\-]{2,79}$", name) else ""


def _int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def read_report(report):
    """One entry of a Datasets ``dataset_report`` as a plain dict, or None when it is
    not usable (no accession with a version, no species name, no size)."""
    if not isinstance(report, dict):
        return None
    accession = normalise_accession(report.get("current_accession") or report.get("accession"))
    organism = report.get("organism") or {}
    info = report.get("assembly_info") or {}
    stats = report.get("assembly_stats") or {}
    species = species_name(organism.get("organism_name"))
    size = _int(stats.get("total_sequence_length"))
    if not accession or not species or size < 1:
        return None
    return {
        "accession": accession,
        "species": species,
        "organism": str(organism.get("organism_name") or "")[:120],
        "common_name": str(organism.get("common_name") or "")[:80],
        "taxon_id": _int(organism.get("tax_id")) or None,
        "assembly_name": str(info.get("assembly_name") or "")[:80],
        "level": str(info.get("assembly_level") or "")[:40],
        "category": str(info.get("refseq_category") or "")[:40],
        "current": str(info.get("assembly_status") or "current") == "current",
        "size_bp": size,
        "contig_n50": _int(stats.get("contig_n50")),
        "scaffold_n50": _int(stats.get("scaffold_n50")),
        "sequences": _int(stats.get("number_of_scaffolds")) or _int(stats.get("number_of_contigs")),
        "annotated": bool(report.get("annotation_info")),
    }


def rank(assembly):
    """Order for a list of candidates: reference genomes first, then the more complete
    assembly, then one with an annotation, then the larger N50."""
    return (0 if assembly["category"] else 1,
            _LEVELS.get(assembly["level"].lower(), 4),
            0 if assembly["annotated"] else 1,
            -assembly["contig_n50"], assembly["accession"])


class Ncbi:
    """``fetch`` and ``stream`` are the two functions of ``net``; tests replace them."""

    def __init__(self, fetch=http_get, stream=http_stream, sleep=time.sleep):
        self.fetch = fetch
        self.stream = stream
        self.sleep = sleep

    # ── what an assembly is ─────────────────────────────────────────────────
    def _reports(self, url):
        try:
            data = json.loads(self.fetch(url, timeout=30, accept="application/json").decode("utf-8"))
        except FetchError as exc:
            raise NcbiError("not_found" if exc.kind == "not_found" else "unavailable", str(exc))
        except (ValueError, UnicodeDecodeError):
            raise NcbiError("bad_answer", "NCBI sent something that is not JSON")
        reports = (data or {}).get("reports") if isinstance(data, dict) else None
        return [r for r in (read_report(x) for x in reports or []) if r]

    def assembly(self, accession):
        accession = normalise_accession(accession)
        if not accession:
            raise NcbiError("not_found", "not an assembly accession")
        for found in self._reports("%s/accession/%s/dataset_report" % (API, accession)):
            if found["accession"] == accession:
                return found
        raise NcbiError("not_found", "NCBI has no assembly %s" % accession)

    def search(self, term, limit=6):
        """Assemblies for a species name or a taxon id, the best first. One per
        species when NCBI has named reference genomes for the taxon."""
        term = re.sub(r"\s+", " ", str(term or "")).strip()
        if not re.match(r"^[A-Za-z0-9][A-Za-z0-9 .\-]{1,79}$", term):
            raise NcbiError("not_found", "not a species name")
        base = "%s/taxon/%s/dataset_report?filters.assembly_version=current&page_size=%d" % (
            API, quote(term, safe=""), max(limit * 3, 20))
        found = self._reports(base + "&filters.reference_only=true")
        if not found:
            found = self._reports(base)
        return sorted(found, key=rank)[:limit]

    # ── where its files are ─────────────────────────────────────────────────
    def files(self, accession):
        """``{"base", "fasta", "gff", "md5": {file name: md5}}`` for one assembly.
        ``gff`` is None when NCBI has no annotation file for it."""
        accession = normalise_accession(accession)
        match = ACCESSION_RE.match(accession or "")
        if not match:
            raise NcbiError("not_found", "not an assembly accession")
        digits = match.group(2)
        parent = "%s/%s/%s/%s/%s/" % (FILES, match.group(1), digits[0:3], digits[3:6], digits[6:9])
        try:
            listing = self.fetch(parent, timeout=30).decode("utf-8", "replace")
        except FetchError as exc:
            raise NcbiError("not_found" if exc.kind == "not_found" else "unavailable", str(exc))
        folders = [name for name in re.findall(r'href="([A-Za-z0-9_.\-#+]+)/"', listing)
                   if name.startswith(accession + "_")]
        if len(folders) != 1:
            raise NcbiError("not_found", "NCBI's file server has no folder for %s" % accession)
        folder = folders[0]
        base = parent + quote(folder, safe="") + "/"
        try:
            sums = self.fetch(base + "md5checksums.txt", timeout=30).decode("utf-8", "replace")
        except FetchError as exc:
            raise NcbiError("not_found" if exc.kind == "not_found" else "unavailable", str(exc))
        md5 = {}
        for line in sums.splitlines():
            parts = line.split()
            if len(parts) == 2 and re.match(r"^[0-9a-f]{32}$", parts[0]):
                md5[parts[1].lstrip("./")] = parts[0]
        fasta, gff = folder + "_genomic.fna.gz", folder + "_genomic.gff.gz"
        if fasta not in md5:
            raise NcbiError("bad_answer", "NCBI lists no checksum for %s" % fasta)
        return {"base": base, "fasta": fasta, "gff": gff if gff in md5 else None, "md5": md5}

    # ── getting one file ────────────────────────────────────────────────────
    def download(self, url, dest, md5, progress=None, stopped=None, attempts=8):
        """Fetch ``url`` to ``dest``, resuming after a broken connection, and compare
        it with ``md5``. Returns the size in bytes. ``progress(done, total)`` is called
        as bytes arrive; ``stopped()`` returning True ends the download."""
        part = str(dest) + ".part"
        wrong = 0
        problem = "no attempt made"
        for attempt in range(attempts):
            if attempt:
                self.sleep(min(60, 2 ** attempt))
            try:
                have = os.path.getsize(part)
            except OSError:
                have = 0
            try:
                response, total = self.stream(url, offset=have)
            except FetchError as exc:
                if exc.kind == "not_found":
                    raise NcbiError("not_found", "%s: %s" % (url, exc))
                if exc.kind != "complete":
                    problem = str(exc)
                    continue
                response, total = None, have
            try:
                if response is not None:
                    if have and not getattr(response, "partial", False):
                        have = 0                 # the server sent the whole file again
                    with open(part, "ab" if have else "wb") as out:
                        while True:
                            if stopped is not None and stopped():
                                raise NcbiError("stopped", "the download was stopped")
                            block = response.read(BLOCK)
                            if not block:
                                break
                            try:
                                out.write(block)
                            except OSError as exc:
                                raise NcbiError("disk", "cannot write %s: %s" % (part, exc))
                            have += len(block)
                            if progress is not None:
                                progress(have, total)
            except NcbiError:
                raise
            except (OSError, ValueError, EOFError, http.client.HTTPException) as exc:
                # the connection was lost: carry on from what has arrived
                problem = "connection lost after %d bytes: %s" % (have, exc)
                continue
            finally:
                if response is not None:
                    try:
                        response.close()
                    except OSError:
                        pass
            if total is not None and have < total:
                problem = "only %d of %d bytes arrived" % (have, total)
                continue
            if file_md5(part) != md5:
                os.unlink(part)                  # a damaged file cannot be resumed
                wrong += 1
                problem = "the file does not match NCBI's checksum"
                if wrong >= 2:
                    raise NcbiError("checksum", "%s: %s, twice" % (url, problem))
                continue
            os.replace(part, str(dest))
            return have
        raise NcbiError("unavailable", "%s: %s" % (url, problem))


def file_md5(path):
    digest = hashlib.md5(usedforsecurity=False)  # the sums NCBI publishes are MD5
    with open(str(path), "rb") as fh:
        while True:
            block = fh.read(BLOCK)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def unpack(source, dest, first_byte, stopped=None):
    """Write the content of the gzip file ``source`` to ``dest`` (through a temporary
    file). The content must start with ``first_byte`` (b'>' for FASTA, b'#' for GFF).
    Returns the unpacked size."""
    temporary = str(dest) + ".tmp"
    name = os.path.basename(str(source))
    size = 0
    try:
        with gzip.open(str(source), "rb") as packed, open(temporary, "wb") as out:
            while True:
                if stopped is not None and stopped():
                    raise NcbiError("stopped", "the download was stopped")
                try:
                    block = packed.read(BLOCK)
                except (OSError, EOFError, zlib.error) as exc:
                    raise NcbiError("bad_answer", "%s cannot be unpacked: %s" % (name, exc))
                if not block:
                    break
                if size == 0 and block[:1] != first_byte:
                    raise NcbiError("bad_answer", "%s does not hold what its name says" % name)
                out.write(block)
                size += len(block)
            out.flush()
            os.fsync(out.fileno())
    except NcbiError:
        _remove(temporary)
        raise
    except OSError as exc:
        _remove(temporary)
        raise NcbiError("disk", "cannot write %s: %s" % (dest, exc))
    if size == 0:
        _remove(temporary)
        raise NcbiError("bad_answer", "%s is empty" % name)
    os.replace(temporary, str(dest))
    return size


def _remove(path):
    try:
        os.unlink(path)
    except OSError:
        pass
