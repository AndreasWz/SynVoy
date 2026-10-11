"""Command line of the web service.

    python -m synvoy_web run                    pages + API + worker in one process (trying it out)
    python -m synvoy_web serve                  pages + API
    python -m synvoy_web worker                 runs the searches
    python -m synvoy_web check                  settings, library, database, pipeline start-up
    python -m synvoy_web library list|add|remove|check        genomes the operator registered
    python -m synvoy_web library status                       the genome store: space, uses
    python -m synvoy_web library fetch NAME_OR_ACCESSION...   get genomes from NCBI into it
    python -m synvoy_web library pin|unpin|evict ACCESSION
    python -m synvoy_web gallery build --spec FILE
    python -m synvoy_web jobs list|show ID|cancel ID
    python -m synvoy_web pause ["text for visitors"]    no new searches, none started
    python -m synvoy_web resume                         open again

Every command takes ``--config FILE`` (or ``$SYNVOY_WEB_CONFIG``).
"""

import argparse
import json
import os
import signal
import sys
import threading
import time
from pathlib import Path

from . import db
from .config import ConfigError, load_settings
from .gallery import GalleryError, build_gallery
from .genomes import EVICTED, FAILED, FETCHING, READY, WANTED, StoreProblem
from .jobs import iso
from .library import (LibraryError, SPECIES_RE, ID_RE, fasta_stats, file_checksum, file_stem,
                      load_library, save_library_data)
from .ncbi import NcbiError, normalise_accession
from .net import FetchError
from .runner import Worker, WorkerError
from .server import build_app, build_context, make_server


def _say(message=""):
    sys.stdout.write(message + "\n")
    sys.stdout.flush()


def _stop_on_signals():
    stop = threading.Event()

    def handler(signum, frame):
        stop.set()

    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGINT, handler)
    return stop


# ── services ────────────────────────────────────────────────────────────────

def cmd_serve(settings, args):
    ctx = build_context(settings)
    server = make_server(build_app(ctx=ctx), settings.host, settings.port)
    stop = _stop_on_signals()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    _say("SynVoy web: http://%s:%d/  (pages from %s)" % (settings.host, settings.port,
                                                        settings.site_dir))
    stop.wait()
    server.shutdown()
    return 0


def cmd_worker(settings, args):
    ctx = build_context(settings)
    worker = Worker(ctx)
    stop = _stop_on_signals()
    try:
        worker.run_forever(stop)
    except WorkerError as exc:
        _say("ERROR: %s" % exc)
        return 1
    return 0


def cmd_run(settings, args):
    ctx = build_context(settings)
    worker = Worker(ctx)
    try:
        worker.lock()
    except WorkerError as exc:
        _say("ERROR: %s" % exc)
        return 1
    server = make_server(build_app(ctx=ctx), settings.host, settings.port)
    stop = _stop_on_signals()
    threading.Thread(target=server.serve_forever, daemon=True).start()
    _say("SynVoy web: http://%s:%d/" % (settings.host, settings.port))
    worker.preflight()
    worker.log("ready (data folder %s)" % settings.data_dir)
    while not stop.is_set():
        try:
            worker.tick()
        except Exception as exc:                 # keep serving; the next tick tries again
            worker.log("error in the main loop: %r" % exc)
        stop.wait(settings.poll_seconds)
    server.shutdown()
    return 0


def cmd_check(settings, args):
    """Everything an operator wants to know before opening the site. Exit code 1 when
    something has to be fixed."""
    failed = False

    def line(ok, text):
        _say("%s %s" % ("[ok]  " if ok else "[FAIL]", text))

    try:
        library = load_library(settings.library_path, strict=True)
        homes = library.homes()
        line(bool(homes), "genome library: %d genome(s), %d usable as starting genome  (%s)"
             % (len(library.genomes), len(homes), settings.library_path))
        failed = failed or not homes
        for genome in library.genomes:
            if genome.size_mb > settings.max_genome_mb:
                _say("       note: %s (%d Mb) is above max_genome_mb and cannot be searched"
                     % (genome.species, round(genome.size_mb)))
            if not genome.lineage:
                _say("       note: %s has no taxonomy; it is listed last, not by relatedness"
                     % genome.species)
    except LibraryError as exc:
        line(False, "genome library: %s" % exc)
        failed = True
    try:
        ctx = build_context(settings)
        state = ctx.store.quick_check()
        line(state == "ok", "database %s: %s" % (settings.db_path, state))
        failed = failed or state != "ok"
    except Exception as exc:
        line(False, "database %s: %s" % (settings.db_path, exc))
        return 1
    site = Path(settings.site_dir)
    line((site / "index.html").is_file(), "pages in %s" % site)
    failed = failed or not (site / "index.html").is_file()
    worker = Worker(ctx)
    result = worker.preflight()
    line(result["ok"], "a search can be started%s"
         % (" (Nextflow %s)" % result["nextflow"] if result.get("nextflow") else ""))
    for problem in result["problems"]:
        _say("       " + problem)
    failed = failed or not result["ok"]
    free_gb = worker.disk_free(settings.data_dir) / 1e9
    line(free_gb >= settings.min_free_gb, "free disk in the data folder: %.0f GB (searches "
         "start above %.0f GB plus their own need)" % (free_gb, settings.min_free_gb))
    failed = failed or free_gb < settings.min_free_gb
    store = ctx.genomes
    if store.enabled:
        usage = store.usage()
        place = store.root if store.root.exists() else Path(settings.data_dir)
        disk_gb = worker.disk_free(place) / 1e9
        still_gb = (usage["budget"] - usage["used"]) / 1e9
        _say("       genome store: %.1f of %.1f GB used, %d genome(s)  (%s)"
             % (usage["used"] / 1e9, usage["budget"] / 1e9, usage["ready"], store.root))
        if disk_gb - settings.min_free_gb < still_gb:
            _say("       note: the store may grow by %.0f GB, the disk has %.0f GB free (%.0f GB "
                 "are kept clear for searches): lower genome_store_gb or free the disk"
                 % (still_gb, disk_gb, settings.min_free_gb))
        if usage["unknown_names"]:
            _say("       note: %d item(s), %.1f GB, in the store's folder are in no record and "
                 "count against its size: see `library sweep`"
                 % (len(usage["unknown_names"]), usage["unknown"] / 1e9))
        own = sorted(store._operators_own())
        if own:
            _say("       note: genomes of the library lie in the store's folder (%s). They are "
                 "left alone; a folder of their own is clearer" % ", ".join(own[:4]))
    else:
        _say("       genome store: off (genome_store_gb is 0); only the genomes of the library "
             "can be searched")
    if settings.smtp_host:
        line(settings.email_enabled, "e-mail: configured (%s)" % settings.smtp_host)
    else:
        _say("       e-mail: off (no smtp_host)")
    if not settings.operator_name or not settings.contact_email:
        _say("       note: operator_name / contact_email are empty; imprint and privacy page "
             "need them")
    return 1 if failed else 0


# ── genome library ──────────────────────────────────────────────────────────

def _library_data(path):
    try:
        with open(str(path)) as fh:
            data = json.load(fh)
    except FileNotFoundError:
        return {"groups": [], "genomes": []}
    except (OSError, ValueError) as exc:
        raise LibraryError("cannot read %s: %s" % (path, exc))
    if not isinstance(data, dict):
        raise LibraryError("%s must hold one JSON object" % path)
    data.setdefault("groups", [])
    data.setdefault("genomes", [])
    return data


def _store_status(ctx):
    store, settings = ctx.genomes, ctx.settings
    usage = store.usage()
    if not store.enabled:
        _say("genome store: off (genome_store_gb is 0). Genomes are not fetched from NCBI.")
    else:
        _say("genome store: %s" % store.root)
        _say("  %.2f of %.2f GB used (%d %%), %d genome(s), %d pinned; %.2f GB on its way"
             % (usage["used"] / 1e9, usage["budget"] / 1e9,
                round(100.0 * usage["used"] / usage["budget"]), usage["ready"], usage["pinned"],
                usage["reserved"] / 1e9))
    if usage["unknown_names"]:
        _say("  NOT IN ANY RECORD: %d item(s), %.2f GB, in the store's folder: %s"
             % (len(usage["unknown_names"]), usage["unknown"] / 1e9,
                ", ".join(usage["unknown_names"][:6])
                + (", ..." if len(usage["unknown_names"]) > 6 else "")))
        _say("    They count against the store's size and are never deleted automatically.")
        _say("    `library sweep` lists them, `library sweep --yes` deletes them.")
    rows = ctx.store.store_genomes()
    counts = store.use_counts()
    if rows:
        _say("  %-17s %-30s %9s  %-8s %5s  %-10s %s" % (
            "accession", "species", "size", "state", "uses", "last used", ""))
    for row in rows:
        uses, last = counts.get(row["id"], (0, 0))
        notes = [n for n in ("pinned" if row["pinned"] else "",
                             "with annotation" if row["with_gff"] else "",
                             row["error_text"] or "") if n]
        _say("  %-17s %-30s %6.0f Mb  %-8s %5d  %-10s %s" % (
            row["accession"], row["species"][:30], row["size_bp"] / 1e6, row["state"], uses,
            time.strftime("%Y-%m-%d", time.gmtime(last)) if last else "never", ", ".join(notes)))
    order = store.eviction_order(ctx.store.store_genomes(READY), store.needed_now())
    if order:
        _say("  first to go when room is needed: " + ", ".join(
            "%s (%s)" % (r["accession"], r["species"]) for r in order[:5]))
    if rows:
        _say("  uses are counted over the last %d days" % settings.store_score_days)
    return 0


def _fetch_one(ctx, term, args):
    """Get one genome into the store, in this process. Returns True when it is there."""
    store = ctx.genomes
    accession = normalise_accession(term)
    if not accession:
        try:
            found = store.ncbi.search(term, limit=1)
        except NcbiError as exc:
            _say("%s: NCBI did not answer (%s)" % (term, exc))
            return False
        if not found:
            _say("%s: NCBI has no assembly under this name" % term)
            return False
        accession = found[0]["accession"]
        _say("%s: %s (%s, %s, %.0f Mb%s)" % (
            term, accession, found[0]["assembly_name"], found[0]["level"] or "level unknown",
            found[0]["size_bp"] / 1e6, ", reference genome" if found[0]["category"] else ""))
    try:
        row = store.want(accession, with_gff=args.annotation, origin="operator", pinned=args.pin)
    except StoreProblem as exc:
        _say("%s: %s" % (accession, exc.message))
        return False
    owner = "cli:%d" % os.getpid()
    tries = 0
    while True:
        if row["state"] == READY:
            _say("%s: on the server (%s, %.0f Mb%s)" % (
                accession, row["species"], row["size_bp"] / 1e6,
                ", with annotation" if row["with_gff"] else ""))
            return True
        if row["state"] in (FAILED, EVICTED):
            _say("%s: NOT fetched: %s" % (accession, row["error_text"] or row["state"]))
            return False
        if row["state"] == FETCHING:          # the worker (or another command) has it
            time.sleep(2.0)
            row = ctx.store.get_store_genome(accession)
            continue
        if row["retry_at"] > ctx.now():       # it failed a moment ago and would wait
            tries += 1
            if tries > 3:
                _say("%s: NOT fetched: %s" % (accession, row["error_text"] or "no reason recorded"))
                return False
            _say("%s: %s; trying again" % (accession, row["error_text"] or "not fetched"))
            ctx.store.update_store_genome(accession, retry_at=0)
        row = store.fetch(accession, owner=owner)


def _library_store(settings, args):
    ctx = build_context(settings)
    store = ctx.genomes
    store.log = _say
    if args.action == "status":
        return _store_status(ctx)
    if args.action == "fetch":
        if not store.enabled:
            _say("The genome store is off: set genome_store_gb in the settings first.")
            return 1
        terms = list(args.terms)
        if args.from_file:
            try:
                terms += [line.split("#")[0].strip() for line in
                          Path(args.from_file).read_text().splitlines()]
            except OSError as exc:
                _say("ERROR: cannot read %s: %s" % (args.from_file, exc))
                return 1
        terms = [t for t in terms if t]
        if not terms:
            _say("Name at least one species or assembly accession.")
            return 2
        missing = 0
        for term in terms:
            usage = store.usage()
            share = 100.0 * (usage["used"] + usage["reserved"] + usage["unknown"]) / usage["budget"]
            if args.fill_to is not None and share >= args.fill_to:
                _say("The store is %d %% full, which is the --fill-to mark. Stopping before %s."
                     % (round(share), term))
                break
            if not _fetch_one(ctx, term, args):
                missing += 1
        _store_status(ctx)
        return 1 if missing else 0
    if args.action == "sweep":
        unknown = store.unknown()
        if not unknown:
            _say("Everything in %s is on record. Nothing to delete." % store.root)
            return 0
        for name, size in unknown:
            _say("  %-44s %10.1f MB" % (name, size / 1e6))
        total = sum(size for _, size in unknown) / 1e9
        if not args.yes:
            _say("%d item(s), %.2f GB, in %s are in no record of the genome store. They were "
                 "left by a database that was lost or replaced, or put there by hand."
                 % (len(unknown), total, store.root))
            _say("`library sweep --yes` deletes exactly these. Genomes registered with "
                 "`library add` are never among them.")
            return 0
        try:
            removed = store.sweep()
        except StoreProblem as exc:
            _say("ERROR: %s" % exc.message)
            return 1
        _say("Deleted %d item(s), %.2f GB." % (len(removed), sum(s for _, s in removed) / 1e9))
        return 0
    accession = normalise_accession(args.accession)
    row = ctx.store.get_store_genome(accession) if accession else None
    if row is None:
        _say("no genome %r in the store" % args.accession)
        return 1
    if args.action in ("pin", "unpin"):
        ctx.store.update_store_genome(accession, pinned=args.action == "pin")
        _say("%s (%s) %s" % (accession, row["species"],
                             "is pinned: it is never removed to make room"
                             if args.action == "pin" else "is no longer pinned"))
        return 0
    if args.action == "evict":
        if accession in store.needed_now():
            _say("%s is needed by a search that is waiting or running. Not removed." % accession)
            return 1
        if row["state"] != READY:
            _say("%s is not on the server (state: %s)" % (accession, row["state"]))
            return 1
        store.evict(row, "removed by the operator")
        return 0
    return 2


def cmd_library(settings, args):
    path = settings.library_path
    if args.action in ("status", "fetch", "pin", "unpin", "evict", "sweep"):
        return _library_store(settings, args)
    try:
        if args.action == "list":
            library = load_library(path, strict=False)
            for group in library.groups:
                _say("%s  (%s)" % (group["name"], group["id"]))
                for genome in library.genomes:
                    if genome.group != group["id"]:
                        continue
                    flags = []
                    flags.append("home+target" if genome.gff else "target only")
                    if not genome.available:
                        flags.append("UNAVAILABLE: " + genome.problem)
                    _say("  %-28s %-34s %8.1f Mb  %s  %s" % (
                        genome.id, genome.species, genome.size_mb, genome.assembly or "-",
                        ", ".join(flags)))
            return 0
        if args.action == "check":
            library = load_library(path, strict=True)
            _say("ok: %d genome(s) in %d group(s)" % (len(library.genomes), len(library.groups)))
            return 0
        if args.action == "remove":
            data = _library_data(path)
            before = len(data["genomes"])
            data["genomes"] = [g for g in data["genomes"] if g.get("id") != args.id]
            if len(data["genomes"]) == before:
                _say("no genome with the id %r" % args.id)
                return 1
            save_library_data(path, data)
            _say("removed %s from the library (the files were not touched)" % args.id)
            return 0
        if args.action == "add":
            return _library_add(settings, args, path)
    except LibraryError as exc:
        _say("ERROR: %s" % exc)
        return 1
    return 2


def _library_add(settings, args, path):
    if not SPECIES_RE.match(args.species):
        raise LibraryError("species %r: 3-80 characters, letters, digits, spaces, dots, hyphens"
                           % args.species)
    genome_id = args.id or file_stem(args.species).lower()
    if not ID_RE.match(genome_id):
        raise LibraryError("id %r: lower-case letters, digits, underscores" % genome_id)
    fasta = Path(args.fasta).expanduser().resolve()
    if not fasta.is_file():
        raise LibraryError("genome file %s does not exist" % fasta)
    gff = Path(args.gff).expanduser().resolve() if args.gff else None
    if gff is not None and not gff.is_file():
        raise LibraryError("annotation file %s does not exist" % gff)
    data = _library_data(path)
    if any(g.get("id") == genome_id for g in data["genomes"]) and not args.replace:
        raise LibraryError("the library already has a genome %r (use --replace)" % genome_id)
    if not any(g.get("id") == args.group for g in data["groups"]):
        if not args.group_name:
            raise LibraryError("group %r does not exist yet; give --group-name to create it"
                               % args.group)
        data["groups"].append({"id": args.group, "name": args.group_name,
                               "note": args.group_note or ""})

    _say("reading %s ..." % fasta.name)
    size_bp, n_sequences, n50 = fasta_stats(str(fasta))
    if size_bp == 0:
        raise LibraryError("%s holds no sequence" % fasta)
    checksum = file_checksum(str(fasta))
    taxon_id, lineage, common = args.taxon_id, [], args.common_name or ""
    if args.lineage:
        try:
            lineage = [int(x) for x in args.lineage.split(",") if x.strip()]
        except ValueError:
            raise LibraryError("--lineage: taxon ids separated by commas")
    elif not args.no_network:
        from .taxonomy import lookup_taxonomy
        try:
            taxon_id, lineage, found_common = lookup_taxonomy(args.species, taxon_id)
            common = common or found_common
        except FetchError as exc:
            _say("WARNING: no taxonomy from NCBI (%s). The genome is added without it and "
                 "will be listed last instead of by relatedness. Add it again with "
                 "--taxon-id to fix this." % exc)
    entry = {
        "id": genome_id, "species": args.species, "common_name": common, "group": args.group,
        "assembly": args.assembly or "", "taxon_id": taxon_id, "lineage": lineage,
        "fasta": str(fasta), "gff": str(gff) if gff else "",
        "size_bp": size_bp, "n_sequences": n_sequences, "n50": n50, "checksum": checksum,
    }
    previous = json.loads(json.dumps(data))
    data["genomes"] = [g for g in data["genomes"] if g.get("id") != genome_id] + [entry]
    save_library_data(path, data)
    try:
        load_library(path, strict=True)
    except LibraryError:
        save_library_data(path, previous)          # leave the library as it was
        raise
    _say("added %s: %s, %.1f Mb in %d sequence(s), N50 %.2f Mb, %s"
         % (genome_id, args.species, size_bp / 1e6, n_sequences, n50 / 1e6,
            "with annotation (can be a starting genome)" if gff else "no annotation (target only)"))
    if size_bp / 1e6 > settings.max_genome_mb:
        _say("NOTE: larger than max_genome_mb (%d Mb); visitors cannot search it until the "
             "limit is raised." % settings.max_genome_mb)
    return 0


# ── examples ────────────────────────────────────────────────────────────────

def cmd_gallery(settings, args):
    try:
        cards = build_gallery(args.spec, args.site or settings.site_dir, log=_say,
                              jobs_dir=settings.jobs_dir)
    except GalleryError as exc:
        _say("ERROR: %s" % exc)
        return 1
    _say("wrote %d example(s) to %s" % (len(cards), Path(args.site or settings.site_dir) / "data"))
    return 0


# ── searches (operator) ─────────────────────────────────────────────────────

def cmd_jobs(settings, args):
    ctx = build_context(settings)
    store = ctx.store
    if args.action == "list":
        states = [args.state] if args.state else list(db.ACTIVE_STATES) + [db.DONE, db.FAILED]
        rows = store.jobs_in_state(*states)
        for job in rows[-args.limit:]:
            request = job["request"] or {}
            _say("%s  %-9s  %s  %-22s %2d genome(s)  %s" % (
                job["id"], job["state"], iso(job["created_at"]),
                (request.get("query") or {}).get("label", "")[:22],
                len(request.get("targets") or []), job["error_code"] or ""))
        _say("%d search(es)" % len(rows))
        return 0
    job = store.get_job(args.id)
    if job is None:
        _say("no search with the id %r" % args.id)
        return 1
    if args.action == "show":
        job["request"]["query"].pop("fasta", None)
        _say(json.dumps(job, indent=2, sort_keys=True))
        _say("folder: %s" % (settings.jobs_dir / job["id"]))
        return 0
    if args.action == "cancel":
        with store.transaction():
            job = store.get_job(args.id)
            if job["state"] == db.QUEUED:
                store.update_job(args.id, state=db.CANCELLED, finished_at=int(time.time()))
                _say("cancelled (it had not started)")
            elif job["state"] == db.RUNNING:
                store.update_job(args.id, cancel_requested=True)
                _say("the worker will stop it within a few seconds")
            else:
                _say("nothing to do: the search is %s" % job["state"])
        return 0
    return 2


def cmd_pause(settings, args):
    """Close the service without stopping anything: running searches finish, waiting
    ones keep their place, new ones are refused with the given text."""
    ctx = build_context(settings)
    text = " ".join(args.message).strip() or "New searches are switched off at the moment."
    ctx.store.set_meta("paused", {"code": "operator", "message": text})
    _say("paused. Visitors are told: %s" % text)
    return 0


def cmd_resume(settings, args):
    ctx = build_context(settings)
    paused = ctx.store.get_meta("paused")
    ctx.store.set_meta("paused", None)
    _say("open again (was: %s)" % paused.get("code") if paused else "the service was not paused")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(prog="python -m synvoy_web", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", help="settings file (JSON); default $SYNVOY_WEB_CONFIG")
    commands = parser.add_subparsers(dest="command")
    for name in ("serve", "run"):
        sub = commands.add_parser(name)
        sub.add_argument("--host")
        sub.add_argument("--port", type=int)
    commands.add_parser("worker")
    commands.add_parser("check")
    commands.add_parser("pause").add_argument("message", nargs="*")
    commands.add_parser("resume")

    library = commands.add_parser("library").add_subparsers(dest="action")
    library.add_parser("list")
    library.add_parser("check")
    remove = library.add_parser("remove")
    remove.add_argument("id")
    add = library.add_parser("add")
    add.add_argument("--species", required=True, help='e.g. "Apis mellifera"')
    add.add_argument("--fasta", required=True, help="genome, .fna / .fa / .fasta (not compressed)")
    add.add_argument("--gff", help="annotation; with it the genome can be a starting genome")
    add.add_argument("--group", required=True, help="group id, e.g. bees")
    add.add_argument("--group-name", help="name of the group if it is new")
    add.add_argument("--group-note", help="one line under the group name")
    add.add_argument("--id", help="default: from the species name")
    add.add_argument("--common-name")
    add.add_argument("--assembly", help="accession, e.g. GCF_003254395.2")
    add.add_argument("--taxon-id", type=int)
    add.add_argument("--lineage", help="taxon ids root to species, comma separated (no look-up)")
    add.add_argument("--no-network", action="store_true", help="do not ask NCBI for the taxonomy")
    add.add_argument("--replace", action="store_true")
    library.add_parser("status", help="the genome store: space used, genomes, uses")
    fetch = library.add_parser("fetch", help="get genomes from NCBI into the genome store")
    fetch.add_argument("terms", nargs="*", metavar="NAME_OR_ACCESSION",
                       help='a species name ("Apis mellifera") or an assembly accession')
    fetch.add_argument("--from-file", help="one name or accession per line")
    fetch.add_argument("--annotation", action="store_true",
                       help="also get the annotation, so the genome can be a starting genome")
    fetch.add_argument("--pin", action="store_true", help="never remove it to make room")
    fetch.add_argument("--fill-to", type=float, metavar="PERCENT",
                       help="stop when the store is this full (e.g. 50)")
    for name in ("pin", "unpin", "evict"):
        library.add_parser(name).add_argument("accession")
    sweep = library.add_parser("sweep", help="list what lies in the genome store's folder "
                                             "without a record; --yes deletes it")
    sweep.add_argument("--yes", action="store_true", help="delete what is listed")

    gallery = commands.add_parser("gallery").add_subparsers(dest="action")
    build = gallery.add_parser("build")
    build.add_argument("--spec", required=True)
    build.add_argument("--site", help="default: site_dir of the settings")

    jobs = commands.add_parser("jobs").add_subparsers(dest="action")
    listing = jobs.add_parser("list")
    listing.add_argument("--state", choices=db.ALL_STATES)
    listing.add_argument("--limit", type=int, default=50)
    for name in ("show", "cancel"):
        jobs.add_parser(name).add_argument("id")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command or (args.command in ("library", "gallery", "jobs") and not args.action):
        parser.print_help()
        return 2
    overrides = {}
    if getattr(args, "host", None):
        overrides["host"] = args.host
    if getattr(args, "port", None):
        overrides["port"] = args.port
    try:
        settings = load_settings(args.config, overrides=overrides)
    except ConfigError as exc:
        _say("ERROR in the settings: %s" % exc)
        return 2
    handler = {"serve": cmd_serve, "worker": cmd_worker, "run": cmd_run, "check": cmd_check,
               "library": cmd_library, "gallery": cmd_gallery, "jobs": cmd_jobs,
               "pause": cmd_pause, "resume": cmd_resume}[args.command]
    return handler(settings, args)
