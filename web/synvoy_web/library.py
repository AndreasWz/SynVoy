"""The genome library: the genomes kept on the server that a visitor may search.

The library is one JSON file (see ``deploy/library.example.json``). Paths in it are
absolute or relative to the file. A genome with an annotation (``gff``) can be the
starting genome of a search ("home"); every genome can be a target.

Names are checked strictly because they travel: the species name becomes a file name
in the search folder and a command-line value, and the pipeline uses the file name as
the genome's identifier in every output.
"""

import hashlib
import json
import os
import re
import threading
import time
from pathlib import Path

ID_RE = re.compile(r"^[a-z0-9][a-z0-9_]{1,63}$")
SPECIES_RE = re.compile(r"^[A-Za-z][A-Za-z0-9 .\-]{2,79}$")
TEXT_RE = re.compile(r"^[^\x00-\x1f<>]{0,120}$")
FASTA_SUFFIXES = (".fna", ".fa", ".fasta")
GFF_SUFFIXES = (".gff", ".gff3")


class LibraryError(ValueError):
    """The library file is wrong in a way that must be fixed by the operator."""


def file_stem(species):
    """File name (without suffix) used for a species in a search folder. Letters,
    digits and underscores only: the pipeline strips suffixes at dots and rejects
    file names with spaces or shell characters."""
    stem = re.sub(r"[^A-Za-z0-9]+", "_", species).strip("_")
    return stem


def fasta_stats(path):
    """(total bases, number of sequences, N50) of a FASTA file, read as a stream."""
    lengths = []
    current = 0
    seen_header = False
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if seen_header:
                    lengths.append(current)
                seen_header = True
                current = 0
            elif seen_header:
                current += len(line.strip())
    if seen_header:
        lengths.append(current)
    total = sum(lengths)
    n50 = 0
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= total:
            n50 = length
            break
    return total, len(lengths), n50


def file_checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(1 << 20)
            if not block:
                break
            digest.update(block)
    return "sha256:" + digest.hexdigest()


class Genome:
    __slots__ = ("id", "species", "common_name", "group", "assembly", "taxon_id", "lineage",
                 "fasta", "gff", "size_bp", "n_sequences", "n50", "checksum", "stem",
                 "available", "problem", "origin", "accession")

    def __init__(self):
        self.origin = "library"       # "ncbi": fetched into the genome store (genomes.py)
        self.accession = ""           # its NCBI assembly accession, for a store genome

    @property
    def size_mb(self):
        return self.size_bp / 1e6

    @property
    def can_be_home(self):
        return bool(self.gff) and self.available

    def identity(self):
        """What makes this genome 'the same genome' for the result cache. An assembly
        accession with its version names one fixed set of sequences, so a genome of
        the store keeps its identity when it is removed and fetched again; its
        annotation can be re-issued by NCBI and is told apart by its checksum."""
        if self.origin == "ncbi":
            return ["ncbi", self.accession, self.checksum]
        return [self.id, self.assembly, self.size_bp, self.checksum]

    def public(self):
        return {
            "id": self.id,
            "species": self.species,
            "common_name": self.common_name,
            "group": self.group,
            "assembly": self.assembly,
            "size_mb": round(self.size_mb, 1),
            "sequences": self.n_sequences,
            "n50": self.n50,
            "can_be_home": self.can_be_home,
            "available": self.available,
            "origin": self.origin,
        }


class Library:
    def __init__(self, groups, genomes, path=None, problems=()):
        self.groups = groups                       # [{id, name, note}]
        self.genomes = genomes                     # [Genome], in file order
        self.by_id = {g.id: g for g in genomes}
        self.by_stem = {g.stem: g for g in genomes}
        self.path = path
        self.problems = list(problems)             # files that are missing right now

    def get(self, genome_id):
        return self.by_id.get(genome_id)

    def homes(self):
        return [g for g in self.genomes if g.can_be_home]

    @staticmethod
    def shared_depth(a, b):
        """Number of taxonomy levels two genomes have in common, counted from the root."""
        depth = 0
        for x, y in zip(a.lineage, b.lineage):
            if x != y:
                break
            depth += 1
        return depth

    def targets_for(self, home):
        """Every other available genome, closest relatives of ``home`` first."""
        others = [g for g in self.genomes if g.id != home.id and g.available]
        return sorted(others, key=lambda g: (-self.shared_depth(home, g), g.species.lower()))

    @classmethod
    def related(cls, anchor, genomes):
        """The genomes near a starting genome: those that share at least half of its
        lineage. The search form uses it for genomes fetched from NCBI, which are
        whatever visitors asked for: a server that holds ant genomes does not tick
        them for a yeast gene. Without a lineage for the starting genome nothing can
        be told, and all are taken."""
        if not anchor.lineage:
            return list(genomes)
        need = max(1, (len(anchor.lineage) + 1) // 2)
        return [g for g in genomes if cls.shared_depth(anchor, g) >= need]

    def home_for(self, taxon_id, organism):
        """The library genome a protein of this organism belongs to, or None.

        Tried in order: the same taxon; the same species name; the same genus and
        species epithet (so "Apis mellifera ligustica" finds "Apis mellifera").
        """
        homes = self.homes()
        if taxon_id:
            for genome in homes:
                if genome.taxon_id == taxon_id:
                    return genome
        name = strip_qualifier(organism or "").lower()
        if not name:
            return None
        for genome in homes:
            if genome.species.lower() == name:
                return genome
        binomial = " ".join(name.split()[:2])
        if len(binomial.split()) == 2:
            for genome in homes:
                if " ".join(genome.species.lower().split()[:2]) == binomial:
                    return genome
        return None

    def public(self):
        return {"groups": self.groups, "genomes": [g.public() for g in self.genomes]}


def strip_qualifier(name):
    """Species name without a strain or isolate in brackets:
    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)" -> "Saccharomyces cerevisiae"."""
    previous = None
    while previous != name:
        previous = name
        name = re.sub(r"\s*\([^()]*\)", "", name)
    return re.sub(r"\s+", " ", name).strip()


def _text(entry, key, where, required=False):
    value = entry.get(key, "")
    if value is None:
        value = ""
    if not isinstance(value, str) or not TEXT_RE.match(value):
        raise LibraryError("%s: %r is not usable text for %s" % (where, value, key))
    if required and not value.strip():
        raise LibraryError("%s: %s is missing" % (where, key))
    return value.strip()


def _file_problem(path, suffixes, what):
    if not path.lower().endswith(suffixes):
        return "%s %s must end in %s" % (what, path, " / ".join(suffixes))
    if not os.path.isfile(path):
        return "%s %s does not exist" % (what, path)
    if not os.access(path, os.R_OK):
        return "%s %s cannot be read" % (what, path)
    return ""


def parse_library(data, base_dir, strict=True):
    """Build a Library from the decoded JSON.

    A mistake in the file itself (unknown key, duplicate id, bad name) always raises.
    A file that is missing on disk raises only with ``strict``; otherwise that genome
    is kept but marked unavailable, so one unmounted disk does not take every genome
    off the site.
    """
    if not isinstance(data, dict):
        raise LibraryError("the library must be one JSON object")
    unknown = set(data) - {"groups", "genomes"} - {k for k in data if k.startswith("_")}
    if unknown:
        raise LibraryError("unknown key(s) at the top of the library: %s" % ", ".join(sorted(unknown)))

    groups = []
    group_ids = set()
    for index, entry in enumerate(data.get("groups") or []):
        where = "group %d" % (index + 1)
        if not isinstance(entry, dict):
            raise LibraryError("%s: expected an object" % where)
        group_id = entry.get("id")
        if not isinstance(group_id, str) or not ID_RE.match(group_id):
            raise LibraryError("%s: id %r must be lower-case letters, digits, underscores" % (where, group_id))
        if group_id in group_ids:
            raise LibraryError("%s: id %r is used twice" % (where, group_id))
        group_ids.add(group_id)
        groups.append({"id": group_id,
                       "name": _text(entry, "name", where, required=True),
                       "note": _text(entry, "note", where)})

    genomes = []
    seen_ids, seen_stems, problems = set(), {}, []
    allowed = {"id", "species", "common_name", "group", "assembly", "taxon_id", "lineage",
               "fasta", "gff", "size_bp", "n_sequences", "n50", "checksum"}
    for index, entry in enumerate(data.get("genomes") or []):
        where = "genome %d" % (index + 1)
        if not isinstance(entry, dict):
            raise LibraryError("%s: expected an object" % where)
        extra = {k for k in entry if not k.startswith("_")} - allowed
        if extra:
            raise LibraryError("%s: unknown key(s) %s" % (where, ", ".join(sorted(extra))))
        genome = Genome()
        genome.id = entry.get("id")
        if not isinstance(genome.id, str) or not ID_RE.match(genome.id):
            raise LibraryError("%s: id %r must be lower-case letters, digits, underscores" % (where, genome.id))
        where = "genome %r" % genome.id
        if genome.id in seen_ids:
            raise LibraryError("%s: id is used twice" % where)
        seen_ids.add(genome.id)

        genome.species = entry.get("species")
        if not isinstance(genome.species, str) or not SPECIES_RE.match(genome.species):
            raise LibraryError(
                "%s: species %r must be 3-80 characters: letters, digits, spaces, dots, hyphens"
                % (where, genome.species))
        genome.species = re.sub(r"\s+", " ", genome.species).strip()
        genome.stem = file_stem(genome.species)
        if genome.stem in seen_stems:
            raise LibraryError("%s: its file name %r is the same as that of genome %r; "
                               "two genomes of one species need different species labels"
                               % (where, genome.stem, seen_stems[genome.stem]))
        seen_stems[genome.stem] = genome.id

        genome.common_name = _text(entry, "common_name", where)
        genome.assembly = _text(entry, "assembly", where)
        genome.checksum = _text(entry, "checksum", where)
        genome.group = entry.get("group")
        if genome.group not in group_ids:
            raise LibraryError("%s: group %r is not defined under \"groups\"" % (where, genome.group))

        taxon = entry.get("taxon_id")
        if taxon is not None and (not isinstance(taxon, int) or isinstance(taxon, bool) or taxon < 1):
            raise LibraryError("%s: taxon_id must be a positive whole number" % where)
        genome.taxon_id = taxon
        lineage = entry.get("lineage") or []
        if not isinstance(lineage, list) or not all(
                isinstance(x, int) and not isinstance(x, bool) for x in lineage):
            raise LibraryError("%s: lineage must be a list of taxon ids, root first" % where)
        genome.lineage = list(lineage)
        if taxon and genome.lineage and genome.lineage[-1] != taxon:
            genome.lineage.append(taxon)

        for key in ("size_bp", "n_sequences", "n50"):
            value = entry.get(key, 0)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise LibraryError("%s: %s must be a whole number" % (where, key))
            setattr(genome, key, value)
        if genome.size_bp < 1:
            raise LibraryError("%s: size_bp is missing (add the genome with `library add`)" % where)

        genome.available, genome.problem = True, ""
        for key, suffixes, what, required in (("fasta", FASTA_SUFFIXES, "genome file", True),
                                              ("gff", GFF_SUFFIXES, "annotation file", False)):
            raw = entry.get(key) or ""
            if not isinstance(raw, str):
                raise LibraryError("%s: %s must be a path" % (where, key))
            if not raw:
                if required:
                    raise LibraryError("%s: %s is missing" % (where, key))
                setattr(genome, key, "")
                continue
            path = raw if os.path.isabs(raw) else os.path.join(base_dir, raw)
            path = os.path.abspath(path)
            setattr(genome, key, path)
            problem = _file_problem(path, suffixes, what)
            if problem:
                if strict:
                    raise LibraryError("%s: %s" % (where, problem))
                genome.available, genome.problem = False, problem
                problems.append("%s: %s" % (where, problem))
        genomes.append(genome)
    return Library(groups, genomes, problems=problems)


def load_library(path, strict=True):
    path = Path(path)
    try:
        with open(str(path)) as fh:
            data = json.load(fh)
    except OSError as exc:
        raise LibraryError("cannot read the genome library %s: %s" % (path, exc))
    except ValueError as exc:
        raise LibraryError("the genome library %s is not valid JSON: %s" % (path, exc))
    library = parse_library(data, str(path.parent), strict=strict)
    library.path = str(path)
    return library


def save_library_data(path, data):
    """Write the library file so that a crash cannot leave half a file behind."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with open(str(temporary), "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=False)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(str(temporary), str(path))


class LibraryWatcher:
    """Hands out the current library and re-reads the file when it has changed, so the
    operator can add a genome without restarting anything. A file that has become
    invalid, or that was there and is gone, keeps the last good library in use and
    reports the error. No file at all is an empty library and no error: nothing has
    been registered yet, and genomes can come from the genome store alone
    (``missing`` is then true)."""

    def __init__(self, path, check_every=5.0, reload_every=300.0, clock=time.monotonic):
        self.path = str(path)
        self.check_every = check_every
        self.reload_every = reload_every      # files can vanish without the library changing
        self._clock = clock
        self._lock = threading.Lock()
        self._library = None
        self._stamp = None
        self._checked = None
        self._loaded = None
        self._had_file = False
        self.error = ""
        self.missing = False

    def _file_stamp(self):
        try:
            st = os.stat(self.path)
            return (st.st_mtime_ns, st.st_size)
        except OSError:
            return None

    def current(self):
        with self._lock:
            now = self._clock()
            if self._library is not None and self._checked is not None \
                    and now - self._checked < self.check_every:
                return self._library
            self._checked = now
            stamp = self._file_stamp()
            fresh = self._loaded is not None and now - self._loaded < self.reload_every
            if self._library is not None and stamp == self._stamp and fresh:
                return self._library
            try:
                library = load_library(self.path, strict=False)
            except LibraryError as exc:
                if stamp is None and not self._had_file:
                    self.error, self.missing = "", True      # never been there: empty
                    self._library = Library([], [], path=self.path)
                    return self._library
                self.error = str(exc)
                if self._library is None:
                    self._library = Library([], [], path=self.path)
                return self._library
            self.error, self.missing, self._had_file = "", False, True
            self._library = library
            self._stamp = stamp
            self._loaded = now
            return library
