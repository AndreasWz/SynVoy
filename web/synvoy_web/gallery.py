"""Finished runs as examples on the start page.

``python -m synvoy_web gallery build --spec FILE`` turns pipeline output folders into
result bundles under ``site/data/examples/<slug>/`` and writes ``site/data/gallery.json``.
These are plain files: the site shows them with or without a running server.

The examples are chosen by the operator. Nothing a visitor submits gets here.

Spec (paths may be relative to the spec file):

    {"examples": [{
        "slug": "melittin",                      required, lower-case
        "run_dir": "runs/melittin",              required, a pipeline --outdir
        "figures_dir": "redrawn/melittin",       optional, figures taken from here instead
        "console_log": "runs/melittin.log",      optional, the run time is read from it
        "name": "Melittin",                      required, shown as the title
        "about": "Honeybee venom peptide",       one line under the title
        "accession": "P01501", "gene": "MELT",
        "scope": "bees, wasps and ants",         "... in 19 <scope>"
        "home": {"species": "Apis mellifera", "common_name": "western honeybee",
                 "assembly": "GCF_003254395.2"},
        "runtime_s": 2640, "label": "1.0.0", "commit": "abc1234",
                                                 label: the release that made the run; left
                                                 out, the example reads "pre-release"
        "rename": {"internal-id-17": "MELT"},    optional: a run started under an internal
                                                 label is published under another name
                                                 (replaced in every text file)
        "must_not_contain": ["internal-id-17"]   optional: words that no published file of
                                                 this example may hold (the build stops
                                                 and names the file)
    }]}

A search that was run on this server can be shown as an example without its output
folder, as long as its result has not expired:

    {"slug": "ste2", "job": "<id of the search>", "about": "Yeast pheromone receptor",
     "scope": "yeasts"}
"""

import json
import os
import re
import shutil
import tempfile
import time
import zipfile
from pathlib import Path

from . import results

SLUG_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{0,40}$")
_DURATION_RE = re.compile(r"Duration:\s*(?:(\d+)d\s*)?(?:(\d+)h\s*)?(?:(\d+)m\s*)?(?:(\d+)s)?")
# Things that should not be in files published for everyone: paths of the machine the
# run was made on. Reported, not removed; the operator decides.
_PATH_RE = re.compile(rb"/(?:home|Users|root|dss|mnt|scratch|data)/[A-Za-z0-9._/-]{3,}")
_KEYS = {"slug", "run_dir", "figures_dir", "console_log", "name", "about", "accession", "gene",
         "scope", "home", "runtime_s", "label", "commit", "job", "must_not_contain", "rename"}
_JOB_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{21}$")


class GalleryError(ValueError):
    pass


def read_duration(path):
    """Seconds from the 'Duration: 3m 56s' line of a console log, or None."""
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return None
    best = None
    for match in _DURATION_RE.finditer(text):
        days, hours, minutes, seconds = (int(g) if g else 0 for g in match.groups())
        total = days * 86400 + hours * 3600 + minutes * 60 + seconds
        if total:
            best = total
    return best


def _query_sequence(run_dir):
    path = Path(run_dir) / "intermediate" / "query" / "normalized_query.faa"
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return ""
    return "".join(l.strip() for l in lines if l and not l.startswith(">"))


def _query_fasta(run_dir, name):
    """The query of a run under the example's own name. The header a run was started
    with is often an internal label, which is not for publishing."""
    sequence = _query_sequence(run_dir)
    if not sequence:
        return ""
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "query"
    return ">%s\n%s\n" % (label, "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60)))


# What results.build_bundle reads of a run folder (relative patterns).
_RUN_PARTS = ("synvoy_report.json", "locus_*_tree.nwk", "synteny_block_locus_*",
              "intermediate/split_loci/locus_*.bed", "intermediate/flanking/synteny_block_locus_*.bed",
              "intermediate/initial_db/initial_db_locus_*.faa", "intermediate/query/normalized_query.faa",
              "intermediate/locate_gene/effective_params.json",
              "plot_inputs_synteny_block_locus_*/*", "rescue/locus_*/*")


def relabelled_copy(source, target, mapping, patterns=_RUN_PARTS):
    """Copy the parts of a run folder that go into a bundle, with every word of
    ``mapping`` replaced in the text files. Returns the number of replacements."""
    source, target = Path(source), Path(target)
    pairs = [(old.encode("utf-8"), new.encode("utf-8")) for old, new in sorted(mapping.items())]
    replaced = 0
    for pattern in patterns:
        for path in sorted(source.glob(pattern)):
            if not path.is_file() or path.is_symlink():
                continue
            data = path.read_bytes()
            if b"\0" not in data[:4096]:                # text
                for old, new in pairs:
                    replaced += data.count(old)
                    data = data.replace(old, new)
            destination = target / path.relative_to(source)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
    return replaced


def scan_for_words(folder, words):
    """Published files that hold one of ``words`` (also inside the archives)."""
    hits = []
    needles = [(w, w.encode("utf-8")) for w in words if w]
    for path in sorted(Path(folder).rglob("*")):
        if not path.is_file():
            continue
        relative = str(path.relative_to(folder))
        try:
            if path.suffix == ".zip":
                with zipfile.ZipFile(str(path)) as archive:
                    parts = [(relative + ":" + n, archive.read(n)) for n in archive.namelist()]
                    parts += [(relative + ":" + n, n.encode("utf-8")) for n in archive.namelist()]
            else:
                parts = [(relative, path.read_bytes())]
        except (OSError, zipfile.BadZipFile):
            continue
        parts.append((relative, relative.encode("utf-8")))
        for name, data in parts:
            for word, needle in needles:
                if needle in data:
                    hits.append((name, word))
    return sorted(set(hits))


def scan_for_paths(folder):
    """Files in a bundle that contain what looks like a path of the run's machine."""
    hits = []
    for path in sorted(Path(folder).rglob("*")):
        if not path.is_file() or path.suffix == ".zip":
            continue
        try:
            match = _PATH_RE.search(path.read_bytes())
        except OSError:
            continue
        if match:
            hits.append((str(path.relative_to(folder)), match.group(0).decode("utf-8", "replace")))
    return hits


def _from_search(example, jobs_dir, bundle):
    """Copy the result of a search run on this server and relabel it as an example."""
    slug, job_id = example["slug"], example["job"]
    if not isinstance(job_id, str) or not _JOB_ID_RE.match(job_id) or jobs_dir is None:
        raise GalleryError("example %r: %r is not the id of a search" % (slug, job_id))
    source = Path(jobs_dir) / job_id / "result"
    try:
        result = json.loads((source / "result.json").read_text())
    except (OSError, ValueError):
        raise GalleryError("example %r: the search %s.. has no result (expired?)"
                           % (slug, job_id[:6]))
    if bundle.exists():
        shutil.rmtree(str(bundle))
    shutil.copytree(str(source), str(bundle))
    result["id"] = slug                       # the link of the search stays private
    result["scope"] = example.get("scope") or ""
    query = result.setdefault("query", {})
    if example.get("name"):
        query["name"] = example["name"]
    query["about"] = example.get("about") or ""
    query["source"] = "example"
    temporary = bundle / "result.json.tmp"
    temporary.write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    os.replace(str(temporary), str(bundle / "result.json"))
    return result


def _finish(result, slug, example, bundle, log):
    """Checks on a finished bundle, and its card for the start page."""
    if not result["self_check"]["agrees"]:
        log("WARNING %s: the counts differ from the report's own headline" % slug)
    for name, found in scan_for_paths(bundle):
        log("WARNING %s: %s contains a local path (%s)" % (slug, name, found))
    words = example.get("must_not_contain") or []
    if not isinstance(words, list) or not all(isinstance(w, str) for w in words):
        raise GalleryError("example %r: must_not_contain must be a list of words" % slug)
    hits = scan_for_words(bundle, words)
    if hits:
        shutil.rmtree(str(bundle))
        raise GalleryError("example %r is not published: %s" % (slug, "; ".join(
            "%s contains %r" % hit for hit in hits[:8])))
    entry = results.card(result)
    entry.update({"slug": slug, "about": example.get("about") or ""})
    summary = result["summary"]
    log("%-12s %d genomes: %d found, %d candidate, %d no call"
        % (slug, summary["genomes"], summary["found"], summary["candidate"], summary["none"]))
    return entry


def build_gallery(spec_path, site_dir, log=print, jobs_dir=None):
    """Build every example of the spec. Returns the list written to gallery.json."""
    spec_path = Path(spec_path)
    try:
        spec = json.loads(spec_path.read_text())
    except OSError as exc:
        raise GalleryError("cannot read %s: %s" % (spec_path, exc))
    except ValueError as exc:
        raise GalleryError("%s is not valid JSON: %s" % (spec_path, exc))
    examples = spec.get("examples") if isinstance(spec, dict) else None
    if not isinstance(examples, list) or not examples:
        raise GalleryError('%s needs a list "examples"' % spec_path)

    def locate(value):
        path = Path(value)
        return path if path.is_absolute() else (spec_path.parent / path)

    data_dir = Path(site_dir) / "data"
    target_root = data_dir / "examples"
    target_root.mkdir(parents=True, exist_ok=True)
    cards, slugs = [], set()
    for index, example in enumerate(examples):
        where = "example %d" % (index + 1)
        if not isinstance(example, dict):
            raise GalleryError("%s: expected an object" % where)
        unknown = set(example) - _KEYS
        if unknown:
            raise GalleryError("%s: unknown key(s) %s" % (where, ", ".join(sorted(unknown))))
        slug = example.get("slug")
        if not isinstance(slug, str) or not SLUG_RE.match(slug) or slug in slugs:
            raise GalleryError("%s: slug %r must be unique, lower-case letters, digits, - or _"
                               % (where, slug))
        slugs.add(slug)
        bundle = target_root / slug
        if example.get("job"):
            result = _from_search(example, jobs_dir, bundle)
            cards.append(_finish(result, slug, example, bundle, log))
            continue
        if not example.get("name") or not example.get("run_dir"):
            raise GalleryError("example %r: name and run_dir are required" % slug)
        run_dir = locate(example["run_dir"])
        if not (run_dir / results.REPORT_NAME).is_file():
            raise GalleryError("example %r: no %s in %s" % (slug, results.REPORT_NAME, run_dir))
        runtime = example.get("runtime_s")
        if runtime is None and example.get("console_log"):
            runtime = read_duration(locate(example["console_log"]))
        meta = {
            "id": slug,
            "scope": example.get("scope") or "",
            "query": {"source": "example", "name": example["name"],
                      "accession": example.get("accession") or "",
                      "gene": example.get("gene") or "", "label": example["name"],
                      "length": len(_query_sequence(run_dir)) or None,
                      "about": example.get("about") or "",
                      "fasta": _query_fasta(run_dir, example.get("gene") or example["name"])},
            "home": example.get("home") or {},
            "targets": [],
            "pipeline": {"label": str(example.get("label") or "") or "pre-release",
                         "release": bool(example.get("label")),
                         "commit": example.get("commit") or "", "modified": False},
            "timing": {"runtime_s": runtime},
        }
        figures = locate(example["figures_dir"]) if example.get("figures_dir") else None
        if figures is not None and not figures.is_dir():
            raise GalleryError("example %r: figures_dir %s does not exist" % (slug, figures))
        rename = example.get("rename") or {}
        if not isinstance(rename, dict) or not all(
                isinstance(k, str) and isinstance(v, str) and len(k) >= 4 and
                re.match(r"^[A-Za-z0-9._-]{1,40}$", v) for k, v in rename.items()):
            raise GalleryError("example %r: rename maps words of at least 4 characters to "
                               "plain names" % slug)
        scratch = None
        try:
            if rename:
                scratch = Path(tempfile.mkdtemp(prefix=".relabel-", dir=str(data_dir)))
                count = relabelled_copy(run_dir, scratch / "run", rename)
                if figures is not None:
                    count += relabelled_copy(figures, scratch / "figures", rename,
                                             ("synteny_block_locus_*",))
                    figures = scratch / "figures"
                log("%-12s relabelled in %d place(s): %s" % (slug, count, ", ".join(
                    "%s -> %s" % pair for pair in sorted(rename.items()))))
                run_dir = scratch / "run"
            result = results.build_bundle(run_dir, bundle, meta, figures_dir=figures)
        finally:
            if scratch is not None:
                shutil.rmtree(str(scratch), ignore_errors=True)
        cards.append(_finish(result, slug, example, bundle, log))

    # Examples that are no longer in the spec are removed, so the folder matches the list.
    for entry in target_root.iterdir():
        if entry.is_dir() and entry.name not in slugs and SLUG_RE.match(entry.name):
            shutil.rmtree(str(entry))
    listing = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "examples": cards}
    temporary = data_dir / "gallery.json.tmp"
    temporary.write_text(json.dumps(listing, indent=1, sort_keys=True) + "\n")
    os.replace(str(temporary), str(data_dir / "gallery.json"))
    return cards
