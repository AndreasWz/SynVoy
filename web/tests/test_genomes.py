"""The genome store: genomes fetched from NCBI on request, kept in a folder of limited
size, counted by use, and removed least-used first when room is needed.

NCBI is a stand-in here (``conftest.FakeNcbi``): it serves genomes made up in the test
through the same two functions the real client uses, and can be told to be down, to
drop a connection part-way, to ignore a request to resume, or to send a damaged file.
"""

import json
import os
import time

import pytest

from conftest import call_json, run_until
from synvoy_web import db, genomes, jobs, ncbi
from synvoy_web.genomes import EVICTED, FAILED, FETCHING, READY, WANTED, StoreProblem
from synvoy_web.library import LibraryWatcher
from synvoy_web.net import FetchError
from synvoy_web.query import QueryError

A, B, C, D = "GCF_000000001.1", "GCA_000000002.1", "GCA_000000003.2", "GCF_000000004.1"
ROOM_FOR_THREE = 13000 / 1e9          # the made-up genomes take about 4 kB each


def stocked(env, accessions=(A, B, C), **kwargs):
    """Genomes at the stand-in NCBI; returns the store."""
    names = {A: "Formica alpha", B: "Lasius beta", C: "Atta gamma", D: "Apis delta"}
    for accession in accessions:
        env.ncbi.add(accession, names[accession], **kwargs)
    return env.genomes


def fetched(env, accession, **want):
    env.genomes.want(accession, origin="operator", **want)
    return env.genomes.fetch(accession)


# ── what NCBI says ──────────────────────────────────────────────────────────

REAL_REPORT = {
    "accession": "GCF_000146045.2", "current_accession": "GCF_000146045.2",
    "paired_accession": "GCA_000146045.2", "source_database": "SOURCE_DATABASE_REFSEQ",
    "organism": {"tax_id": 559292, "organism_name": "Saccharomyces cerevisiae S288C",
                 "infraspecific_names": {"strain": "S288C"}},
    "assembly_info": {"assembly_level": "Complete Genome", "assembly_status": "current",
                      "assembly_name": "R64", "refseq_category": "reference genome"},
    "assembly_stats": {"total_sequence_length": "12071326", "contig_n50": 924431,
                       "scaffold_n50": 924431, "number_of_scaffolds": 16},
    "annotation_info": {"name": "Annotation submitted by SGD"},
}


def test_report_of_an_assembly_is_read():
    """The shape NCBI's service answered with on 2026-10-11."""
    found = ncbi.read_report(REAL_REPORT)
    assert found == {
        "accession": "GCF_000146045.2", "species": "Saccharomyces cerevisiae",
        "organism": "Saccharomyces cerevisiae S288C", "common_name": "", "taxon_id": 559292,
        "assembly_name": "R64", "level": "Complete Genome", "category": "reference genome",
        "current": True, "size_bp": 12071326, "contig_n50": 924431, "scaffold_n50": 924431,
        "sequences": 16, "annotated": True}
    for broken in (None, {}, {"accession": "GCF_1"}, dict(REAL_REPORT, organism={}),
                   dict(REAL_REPORT, assembly_stats={"total_sequence_length": "many"})):
        assert ncbi.read_report(broken) is None


@pytest.mark.parametrize("organism, expected", [
    ("Saccharomyces cerevisiae S288C", "Saccharomyces cerevisiae"),
    ("Apis mellifera ligustica", "Apis mellifera"),
    ("Probolomyrmex sp. GAGA-0580", "Probolomyrmex sp. GAGA-0580"),
    ("Tetramorium bicarinatum", "Tetramorium bicarinatum"),
    ("Candida albicans (strain SC5314)", "Candida albicans"),
    ("[Candida] auris", "auris"),
    ("x; rm -rf /", "x rm"), ("", ""), ("A", ""),
])
def test_species_name_is_made_safe_for_a_file_name(organism, expected):
    assert ncbi.species_name(organism) == expected


@pytest.mark.parametrize("text, expected", [
    ("GCF_000146045.2", "GCF_000146045.2"), (" gca_928718305.1 ", "GCA_928718305.1"),
    ("GCF_000146045", None), ("GCF_000146045.2/../../etc", None), ("NC_001133.9", None),
    ("GCF_00014604.2", None), (None, None), (5, None), ("", None),
])
def test_only_assembly_accessions_with_a_version_are_taken(text, expected):
    assert ncbi.normalise_accession(text) == expected


def test_candidates_are_ordered_reference_and_complete_first(env):
    env.ncbi.add(A, "Formica alpha", level="Scaffold", category="")
    env.ncbi.add(B, "Formica alpha", level="Chromosome", category="")
    env.ncbi.add(C, "Formica alpha", level="Contig", category="reference genome", annotated=True)
    client = env.ncbi.client()
    # NCBI names a reference genome: that one alone is offered.
    assert [a["accession"] for a in client.search("Formica alpha")] == [C]
    env.ncbi.assemblies[C]["report"]["assembly_info"]["refseq_category"] = ""
    assert [a["accession"] for a in client.search("formica alpha")] == [B, A, C]
    assert client.search("Nomen nescio") == []
    with pytest.raises(ncbi.NcbiError):
        client.search("x'; drop table")                       # not a name: never sent


# ── fetching ────────────────────────────────────────────────────────────────

def test_genome_is_fetched_checked_and_listed(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    row = store.want(A, client="visitor-1")
    assert row["state"] == WANTED and row["id"] == "gcf_000000001_1"
    assert env.ctx.library().get(row["id"]) is None            # not listed before it is here
    row = store.fetch(A)
    assert row["state"] == READY and row["error_text"] is None
    folder = store.root / A
    fasta = folder / "Formica_alpha.fna"
    assert fasta.read_bytes() == env.ncbi.assemblies[A]["plain"]
    assert sorted(p.name for p in folder.iterdir()) == ["Formica_alpha.fna"]   # nothing packed left
    assert not (store.root / genomes.INCOMING / A).exists()
    assert (row["size_bp"], row["n_sequences"], row["disk_bytes"]) == (4000, 2, fasta.stat().st_size)
    genome = env.ctx.library().get(row["id"])
    assert genome.available and genome.origin == "ncbi" and genome.group == "ncbi"
    assert genome.species == "Formica alpha" and genome.stem == "Formica_alpha"
    assert genome.lineage == env.ncbi.lineages[row["taxon_id"]]
    assert not genome.can_be_home                              # no annotation was fetched
    # The accession with its version is the genome's identity, wherever the file lies.
    assert genome.identity() == ["ncbi", A, ""]
    assert "ncbi" in [g["id"] for g in env.ctx.library().groups]
    assert env.ctx.library().get("homus_primus").origin == "library"


def test_annotation_comes_with_a_genome_that_is_to_be_a_starting_genome(make_env):
    env = make_env(genome_store_gb=0.001)
    env.ncbi.add(A, "Formica alpha", annotated=True)
    env.ncbi.add(B, "Lasius beta")
    with pytest.raises(StoreProblem) as refused:               # NCBI has no annotation for it
        env.genomes.want(B, with_gff=True)
    assert refused.value.code == "home_without_annotation"
    row = fetched(env, A, with_gff=True)
    assert row["state"] == READY and row["with_gff"]
    genome = env.ctx.library().get(row["id"])
    assert genome.can_be_home and open(genome.gff).read().startswith("##gff-version 3")
    assert genome.identity()[2].startswith("gff:") and len(genome.identity()[2]) == 36


def test_annotation_is_added_to_a_genome_that_is_already_here(make_env):
    env = make_env(genome_store_gb=0.001)
    env.ncbi.add(A, "Formica alpha", annotated=True)
    assert fetched(env, A)["with_gff"] is False
    before = len([c for c in env.ncbi.calls if c.endswith("_genomic.fna.gz @0")])
    row = env.genomes.want(A, with_gff=True, origin="operator")
    assert row["state"] == WANTED
    assert env.ctx.library().get(row["id"]) is None            # off the list while it changes
    row = env.genomes.fetch(A)
    assert row["state"] == READY and row["with_gff"] and env.ctx.library().get(row["id"]).can_be_home
    assert len([c for c in env.ncbi.calls if c.endswith("_genomic.fna.gz @0")]) == before   # not twice


def test_broken_connection_is_continued_where_it_stopped(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env, size=400000)
    name = env.ncbi.assemblies[A]["folder"] + "_genomic.fna.gz"
    size = len(env.ncbi.assemblies[A]["fasta"])
    env.ncbi.drops[name] = [100, 150]                          # two breaks
    assert fetched(env, A)["state"] == READY
    assert (store.root / A / "Formica_alpha.fna").read_bytes() == env.ncbi.assemblies[A]["plain"]
    offsets = [int(c.rsplit("@", 1)[1]) for c in env.ncbi.calls if name in c]
    assert offsets == [0, 100, 250] and size > 250             # each try starts where the last ended


def test_server_that_cannot_resume_sends_the_file_again(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env, size=400000)
    name = env.ncbi.assemblies[A]["folder"] + "_genomic.fna.gz"
    env.ncbi.drops[name] = [120]
    env.ncbi.ignore_range = True
    assert fetched(env, A)["state"] == READY
    assert (store.root / A / "Formica_alpha.fna").read_bytes() == env.ncbi.assemblies[A]["plain"]


def test_damaged_file_is_fetched_again_and_refused_the_second_time(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    name = env.ncbi.assemblies[A]["folder"] + "_genomic.fna.gz"
    env.ncbi.damaged[name] = 1
    assert fetched(env, A)["state"] == READY                   # the second copy was sound
    assert (store.root / A / "Formica_alpha.fna").read_bytes() == env.ncbi.assemblies[A]["plain"]

    name = env.ncbi.assemblies[B]["folder"] + "_genomic.fna.gz"
    env.ncbi.damaged[name] = 5
    row = fetched(env, B)
    assert row["state"] == FAILED and row["error_code"] == "checksum"
    assert row["error_text"] == "the files arrived damaged"
    assert not (store.root / B).exists() and not (store.root / genomes.INCOMING / B).exists()
    assert env.ctx.library().get(row["id"]) is None
    with pytest.raises(StoreProblem) as again:                 # a visitor is told, not made to wait
        store.want(B, client="v")
    assert again.value.code == "genome_download_failed" and "damaged" in again.value.message


def test_file_that_is_not_what_its_name_says_is_refused(make_env):
    env = make_env(genome_store_gb=0.001)
    stocked(env)
    import gzip
    entry = env.ncbi.assemblies[A]
    entry["fasta"] = gzip.compress(b"<html>Service unavailable</html>", mtime=0)
    row = fetched(env, A)
    assert row["state"] == FAILED and row["error_code"] == "bad_answer"
    # The same for an annotation: nothing else would notice a page of text in its place.
    env.ncbi.add(D, "Apis delta", annotated=True)
    env.ncbi.assemblies[D]["gff"] = gzip.compress(b"<html>Service unavailable</html>", mtime=0)
    row = fetched(env, D, with_gff=True)
    assert row["state"] == FAILED and row["error_code"] == "bad_answer"
    assert not (env.genomes.root / D).exists()


def test_ncbi_being_down_is_tried_again_later_and_then_given_up(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    store.want(A, origin="operator")
    env.ncbi.down = True
    env.now = 1000000
    row = store.fetch(A)
    assert row["state"] == WANTED and row["attempts"] == 1 and row["retry_at"] == env.now + 240
    assert store.next_wanted() is None                         # not before its time
    env.now = row["retry_at"]
    assert store.next_wanted()["accession"] == A
    row = store.fetch(A)
    assert row["state"] == WANTED and row["retry_at"] == env.now + 960
    env.now = row["retry_at"]
    row = store.fetch(A)
    assert row["state"] == FAILED and row["error_text"] == "NCBI did not answer"
    env.ncbi.down = False
    env.now += genomes.FAILED_FOR + 1                          # later, a new request tries afresh
    assert store.want(A, client="v")["state"] == WANTED
    assert store.fetch(A)["state"] == READY


def test_download_that_takes_too_long_is_given_up(make_env, monkeypatch):
    env = make_env(genome_store_gb=0.001, fetch_timeout_minutes=1)
    store = stocked(env)
    store.want(A, origin="operator")
    clock = [0.0]
    monkeypatch.setattr(genomes.time, "monotonic", lambda: clock.__setitem__(0, clock[0] + 40) or clock[0])
    row = store.fetch(A)
    assert row["state"] == WANTED and row["error_text"] == "the download took too long"


def test_stopping_the_worker_keeps_what_has_arrived(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env, size=400000)
    env.ncbi.chunk = 64                                        # the file arrives in small pieces
    name = env.ncbi.assemblies[A]["folder"] + "_genomic.fna.gz"
    store.want(A, origin="operator")
    seen = []
    row = store.fetch(A, stopped=lambda: bool(seen.append(1)) or len(seen) > 1)
    assert row["state"] == WANTED and row["attempts"] == 0     # not counted as a failed try
    part = store.root / genomes.INCOMING / A / "genomic.fna.gz.part"
    assert part.is_file() and 0 < part.stat().st_size < len(env.ncbi.assemblies[A]["fasta"])
    kept = part.stat().st_size
    assert store.fetch(A)["state"] == READY                    # the next worker goes on
    assert ("%s @%d" % (name, kept)) in " ".join(env.ncbi.calls)


def test_download_cut_off_by_a_restart_is_wanted_again(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    store.want(A, origin="operator")
    env.store.update_store_genome(A, state=FETCHING, owner="worker:1", attempts=1)
    assert store.fetch(A)["state"] == FETCHING                 # somebody has it: hands off
    assert not [c for c in env.ncbi.calls if "genomic" in c]
    store.recover(alive=lambda owner: owner == "worker:1")     # its owner lives: left alone
    assert env.store.get_store_genome(A)["state"] == FETCHING
    store.recover()
    row = env.store.get_store_genome(A)
    assert (row["state"], row["attempts"], row["owner"]) == (WANTED, 0, "")
    assert any("interrupted" in line for line in env.genome_logs)
    assert store.fetch(A)["state"] == READY


# ── room ────────────────────────────────────────────────────────────────────

def test_least_used_genome_goes_first_when_room_is_needed(make_env):
    env = make_env(genome_store_gb=ROOM_FOR_THREE)
    store = stocked(env, (A, B, C, D))
    env.now = 2000000
    for accession in (A, B, C):
        assert fetched(env, accession)["state"] == READY
    ids = {acc: env.store.get_store_genome(acc)["id"] for acc in (A, B, C)}
    store.record_use([ids[A], ids[A], ids[C]])                 # one search used A and C
    env.now += 10
    store.record_use([ids[A]])
    # B was never used: it goes. A (2 uses) and C (1 use) stay.
    assert [r["accession"] for r in store.eviction_order(env.store.store_genomes(READY), set())] == [B, C, A]
    assert fetched(env, D)["state"] == READY
    states = {acc: env.store.get_store_genome(acc)["state"] for acc in (A, B, C, D)}
    assert states == {A: READY, B: EVICTED, C: READY, D: READY}
    assert not (store.root / B).exists()
    assert env.ctx.library().get(ids[B]) is None and env.ctx.library().get(ids[A]) is not None
    assert any("removed %s" % B in line for line in env.genome_logs)
    usage = store.usage()
    assert usage["ready"] == 3 and usage["used"] <= usage["budget"]


def test_old_uses_no_longer_count(make_env):
    env = make_env(genome_store_gb=ROOM_FOR_THREE, store_score_days=30)
    store = stocked(env)
    env.now = 5000000
    for accession in (A, B):
        fetched(env, accession)
    ids = {acc: env.store.get_store_genome(acc)["id"] for acc in (A, B)}
    store.record_use([ids[A]] * 5)                             # much used, long ago
    env.now += 40 * 86400
    store.record_use([ids[B]])                                 # used once, recently
    order = [r["accession"] for r in store.eviction_order(env.store.store_genomes(READY), set())]
    assert order == [A, B]


def test_pinned_and_needed_genomes_are_never_removed(make_env):
    env = make_env(genome_store_gb=ROOM_FOR_THREE)
    store = stocked(env, (A, B, C, D))
    for accession in (A, B, C):
        fetched(env, accession)
    # C is the one in use; A and B were never used and would be the first to go.
    store.record_use([env.store.get_store_genome(C)["id"]])
    assert [r["accession"] for r in store.eviction_order(env.store.store_genomes(READY), set())][-1] == C
    env.store.update_store_genome(A, pinned=True)
    env.beat()
    waiting = jobs.submit(env.ctx, env.payload(targets=[{"accession": B}]), "10.0.0.1")
    assert genomes.store_refs(env.store.get_job(waiting["id"])["request"])[0]["accession"] == B
    fetched(env, D)                                            # C is the only one that may go
    states = {acc: env.store.get_store_genome(acc)["state"] for acc in (A, B, C, D)}
    assert states == {A: READY, B: READY, C: EVICTED, D: READY}
    # Nothing can go any more: a further genome is refused at once, with a reason.
    env.store.update_store_genome(D, pinned=True)
    env.ncbi.add("GCA_000000009.1", "Novum genus")
    with pytest.raises(StoreProblem) as full:
        store.want("GCA_000000009.1", client="v")
    assert full.value.code == "no_room" and full.value.status == 503
    assert env.store.get_store_genome("GCA_000000009.1") is None


def test_genome_larger_than_the_limits_is_refused_before_any_download(make_env):
    env = make_env(genome_store_gb=0.001, max_genome_mb=1)
    env.ncbi.add(A, "Formica alpha", reported=2500000)
    env.ncbi.add(B, "Lasius beta", reported=900000)            # under the limit, over the store
    with pytest.raises(StoreProblem) as too_large:
        env.genomes.want(A, client="v")
    assert too_large.value.code == "genome_too_large" and "2 Mb" in too_large.value.message
    env.settings.genome_store_gb = 0.0005
    with pytest.raises(StoreProblem) as no_room:
        env.genomes.want(B, client="v")
    assert no_room.value.code == "no_room"
    assert not [c for c in env.ncbi.calls if "genomic.fna.gz" in c]
    assert env.store.store_genomes() == []


def test_full_disk_waits_instead_of_failing(make_env):
    env = make_env(genome_store_gb=0.001, min_free_gb=1.0)
    store = stocked(env)
    store.want(A, origin="operator")
    env.store_disk_free = int(1.0e9) + 100                     # the disk, not the store, is full
    env.now = 3000000
    row = store.fetch(A)
    assert row["state"] == WANTED and row["error_code"] == "low_disk"
    assert row["retry_at"] == env.now + genomes.ROOM_RETRY
    env.store_disk_free = 10 ** 13
    env.now = row["retry_at"]
    assert store.fetch(A)["state"] == READY


def test_store_that_is_switched_off_fetches_nothing(env):
    stocked(env)
    assert env.genomes.enabled is False
    with pytest.raises(StoreProblem) as off:
        env.genomes.want(A, client="v")
    assert off.value.code == "downloads_off"
    env.beat()
    with pytest.raises(QueryError) as refused:
        jobs.submit(env.ctx, env.payload(targets=[{"accession": A}]), "10.0.0.1")
    assert refused.value.code == "downloads_off"
    assert env.genomes.merged(env.ctx.libraries.current()) is env.ctx.libraries.current()


def test_one_visitor_cannot_have_many_genomes_fetched(make_env):
    env = make_env(genome_store_gb=0.001, max_fetches_per_day=2)
    store = stocked(env, (A, B, C))
    store.want(A, client="v1")
    store.want(A, client="v1")                                 # asking again costs nothing
    store.want(B, client="v1")
    with pytest.raises(StoreProblem) as limit:
        store.want(C, client="v1")
    assert limit.value.code == "too_many_fetches" and limit.value.status == 429
    assert store.want(C, client="v2")["state"] == WANTED       # another visitor may


# ── searches that need a genome which is not here ───────────────────────────

def submit(env, targets, name="ok", address="10.0.0.1", **extra):
    env.beat()
    return jobs.submit(env.ctx, env.payload(name=name, targets=targets, **extra), address)


def test_search_waits_for_its_genome_and_then_runs(make_env):
    env = make_env(genome_store_gb=0.001)
    stocked(env)
    worker = env.worker()
    outcome = submit(env, ["alpha_beta", {"accession": A}])
    job = env.store.get_job(outcome["id"])
    first_key = job["cache_key"]
    refs = {t["species"]: t for t in job["request"]["targets"]}
    assert refs["Formica alpha"]["store"] and refs["Formica alpha"]["identity"] == ["ncbi", A, ""]
    assert "store" not in refs["Alpha beta"]
    waiting = jobs.job_view(env.ctx, job)["queue"]["genomes"]
    assert [(g["accession"], g["species"], g["state"]) for g in waiting] == [(A, "Formica alpha", WANTED)]

    # The same question asked meanwhile joins this search.
    again = submit(env, ["alpha_beta", {"accession": A}], address="10.0.0.2")
    assert again["id"] == outcome["id"] and again["reused"] is True

    run_until(worker, lambda: env.store.get_job(outcome["id"])["state"] in (db.DONE, db.FAILED))
    job = env.store.get_job(outcome["id"])
    assert job["state"] == db.DONE, job["error_text"]
    row = env.store.get_store_genome(A)
    assert row["state"] == READY
    # The pipeline was handed the fetched file under the species' name.
    console = (env.settings.jobs_dir / job["id"] / "logs" / "console.log").read_text()
    assert "--target_genomes" in console
    result = json.loads((env.settings.jobs_dir / job["id"] / "result" / "result.json").read_text())
    assert sorted(g["id"] for g in result["genomes"]) == ["Alpha_beta", "Formica_alpha"]
    # One use each of the starting genome and the two targets.
    counts = env.store.genome_use_counts()
    assert {k: v[0] for k, v in counts.items()} == {"homus_primus": 1, "alpha_beta": 1, row["id"]: 1}
    # Asked again, by accession or from the list: the same result, no new run.
    for targets in (["alpha_beta", {"accession": A}], ["alpha_beta", row["id"]]):
        third = submit(env, targets, address="10.0.0.3")
        assert third["id"] == outcome["id"] and third["reused"] is True
    assert env.store.get_job(outcome["id"])["cache_key"] == first_key


def test_other_searches_start_while_one_waits_for_its_genome(make_env):
    env = make_env(genome_store_gb=0.001)
    stocked(env)
    worker = env.worker()
    waiting = submit(env, [{"accession": A}], name="ok_first")
    env.ncbi.down = True                                       # the download cannot happen yet
    behind = submit(env, ["alpha_beta"], name="ok_second", address="10.0.0.2")
    run_until(worker, lambda: env.store.get_job(behind["id"])["state"] == db.DONE)
    assert env.store.get_job(waiting["id"])["state"] == db.QUEUED


def test_search_fails_with_the_reason_when_its_genome_cannot_be_fetched(make_env):
    env = make_env(genome_store_gb=0.001)
    stocked(env)
    name = env.ncbi.assemblies[A]["folder"] + "_genomic.fna.gz"
    env.ncbi.damaged[name] = 9
    worker = env.worker()
    outcome = submit(env, [{"accession": A}])
    run_until(worker, lambda: env.store.get_job(outcome["id"])["state"] in (db.DONE, db.FAILED))
    job = env.store.get_job(outcome["id"])
    assert (job["state"], job["error_code"]) == (db.FAILED, "genome_download_failed")
    assert "Formica alpha" in job["error_text"] and "damaged" in job["error_text"]
    assert job["attempt"] == 0                                 # the pipeline was never started


def test_genome_removed_while_a_search_waited_is_fetched_again(make_env):
    env = make_env(genome_store_gb=ROOM_FOR_THREE)
    store = stocked(env)
    row = fetched(env, A)
    outcome = submit(env, [row["id"]])                         # chosen from the list
    assert genomes.store_refs(env.store.get_job(outcome["id"])["request"])[0]["accession"] == A
    store.evict(row, "the test says so")                       # as if room had been needed
    assert env.ctx.library().get(row["id"]) is None
    asked = len(env.ncbi.calls)
    worker = env.worker()
    run_until(worker, lambda: env.store.get_job(outcome["id"])["state"] in (db.DONE, db.FAILED))
    assert env.store.get_job(outcome["id"])["state"] == db.DONE
    assert env.store.get_store_genome(A)["state"] == READY
    # What the genome is was on record: NCBI was asked for its files only.
    assert not [c for c in env.ncbi.calls[asked:] if "dataset_report" in c][1:]


def test_starting_genome_can_come_from_ncbi(make_env):
    env = make_env(genome_store_gb=0.001)
    env.ncbi.add(A, "Formica alpha", annotated=True)
    worker = env.worker()
    outcome = submit(env, ["alpha_beta"], home={"accession": A})
    request = env.store.get_job(outcome["id"])["request"]
    assert request["home"]["store"] and request["home"]["with_gff"] is True
    run_until(worker, lambda: env.store.get_job(outcome["id"])["state"] in (db.DONE, db.FAILED))
    job = env.store.get_job(outcome["id"])
    assert job["state"] == db.DONE, job["error_text"]
    console = (env.settings.jobs_dir / job["id"] / "logs" / "console.log").read_text()
    assert "Formica_alpha.gff" in console and "--home_species 'Formica alpha'" in console
    assert job["request"]["home"]["identity"][2].startswith("gff:")     # known once fetched


def without_library(env):
    """The same installation before any genome was registered with ``library add``."""
    env.library_file.unlink()
    env.ctx.libraries = LibraryWatcher(env.settings.library_path, check_every=0.0)
    return env


def test_server_can_run_on_the_genome_store_alone(make_env):
    """No genome registered by hand: the store is filled by the operator's
    ``library fetch`` and by what visitors name."""
    env = without_library(make_env(genome_store_gb=0.001))
    env.ncbi.add(A, "Formica alpha", annotated=True)
    env.ncbi.add(B, "Lasius beta")
    worker = env.worker()
    ready = worker.preflight()
    assert ready["ok"], ready["problems"]
    worker.tick()
    assert jobs.service_state(env.ctx)[0] is True         # searches are taken
    assert env.ctx.library().genomes == []
    outcome = jobs.submit(env.ctx, env.payload(targets=[{"accession": B}], home={"accession": A}),
                          "10.0.0.1")
    run_until(worker, lambda: env.store.get_job(outcome["id"])["state"] in (db.DONE, db.FAILED))
    job = env.store.get_job(outcome["id"])
    assert job["state"] == db.DONE, job["error_text"]
    assert sorted(g.species for g in env.ctx.library().genomes) == ["Formica alpha", "Lasius beta"]
    assert env.genomes.sweep() == []                     # and the store knows what is its own


def test_server_without_any_source_of_genomes_says_so(make_env):
    env = without_library(make_env())                    # no library, store off
    ready = env.worker().preflight()
    assert not ready["ok"] and "genome_store_gb" in " ".join(ready["problems"])
    env.beat()
    assert jobs.service_state(env.ctx)[:2] == (False, "library_empty")


def test_what_a_request_may_and_may_not_name(make_env):
    env = make_env(genome_store_gb=0.001)
    stocked(env)
    env.ncbi.add(D, "Alpha beta")                              # the species of a library genome
    env.beat()

    def refused(**payload):
        with pytest.raises(QueryError) as caught:
            jobs.submit(env.ctx, env.payload(**payload), "10.0.0.1")
        return caught.value.code

    assert refused(targets=[{"accession": "GCF_000000001"}]) == "genome_unknown"      # no version
    assert refused(targets=[{"accession": "GCF_000000001.1/../x"}]) == "genome_unknown"
    assert refused(targets=[{"accession": "GCF_999999999.1"}]) == "genome_unknown"    # not at NCBI
    assert refused(targets=[{"accession": A, "path": "/etc"}]) == "bad_request"
    assert refused(targets=[{"file": "/etc/passwd"}]) == "bad_request"
    assert refused(targets=[{"accession": A}, {"accession": A}]) == "bad_request"     # twice
    assert refused(targets=["alpha_beta", {"accession": D}]) == "same_species_twice"
    assert refused(targets=[{"accession": A}], home={"accession": B}) == "home_without_annotation"
    assert env.store.store_genomes() == [] or all(
        r["accession"] in (A, B, D) for r in env.store.store_genomes())


def test_request_nobody_waits_for_is_dropped(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    env.now = 4000000
    outcome = submit(env, [{"accession": A}])
    store.want(B, origin="operator")
    (store.root / genomes.INCOMING / "GCA_000000077.1").mkdir(parents=True)   # left by a crash
    jobs.cancel(env.ctx, outcome["id"], outcome["token"])
    assert store.next_wanted()["accession"] == B               # nobody waits for A any more
    worker = env.worker()
    run_until(worker, lambda: env.store.get_job(outcome["id"])["state"] == db.CANCELLED
              and env.store.get_store_genome(B)["state"] == READY)
    store.housekeeping()
    assert env.store.get_store_genome(A)["state"] == WANTED    # an hour of grace
    env.now += genomes.UNWANTED_AFTER + 1
    store.housekeeping()
    assert env.store.get_store_genome(A) is None               # the visitor's request is gone
    assert env.store.get_store_genome(B)["state"] == READY     # what the operator asked for is here
    assert not (store.root / genomes.INCOMING / "GCA_000000077.1").exists()
    assert store.next_wanted() is None


def test_download_whose_owner_died_is_taken_up_again(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    env.now = 6000000
    store.want(A, origin="operator")
    env.store.update_store_genome(A, state=FETCHING, owner="cli:999999", attempts=1,
                                  updated_at=env.now)
    store.housekeeping()
    assert env.store.get_store_genome(A)["state"] == FETCHING  # it may still be at work
    env.now += genomes.STALE_FETCH + 1
    store.housekeeping()
    assert env.store.get_store_genome(A)["state"] == WANTED
    from synvoy_web.runner import Worker
    assert Worker.fetch_owner_alive("cli:%d" % os.getpid()) is True
    assert Worker.fetch_owner_alive("cli:999999999") is False
    assert Worker.fetch_owner_alive("worker:%d" % os.getpid()) is False     # an older worker is gone


def test_store_made_smaller_is_brought_under_its_size(make_env):
    env = make_env(genome_store_gb=ROOM_FOR_THREE)
    store = stocked(env)
    for accession in (A, B, C):
        fetched(env, accession)
    store.record_use([env.store.get_store_genome(C)["id"]])
    env.settings.genome_store_gb = 5000 / 1e9                  # room for one
    store.housekeeping()
    states = {acc: env.store.get_store_genome(acc)["state"] for acc in (A, B, C)}
    assert states == {A: EVICTED, B: EVICTED, C: READY}


def test_files_that_vanished_take_the_genome_off_the_list(make_env):
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    row = fetched(env, A)
    os.unlink(str(store.root / A / "Formica_alpha.fna"))
    assert env.genomes.usable(env.store.get_store_genome(A)) is False
    store.housekeeping()
    assert env.store.get_store_genome(A)["state"] == EVICTED
    assert env.ctx.library().get(row["id"]) is None


def test_genome_stored_without_its_taxonomy_gets_it_later(make_env):
    """NCBI's taxonomy did not answer while the genome was fetched: the genome is kept,
    listed last, and placed among its relatives as soon as the lineage is known."""
    env = make_env(genome_store_gb=0.001)
    store = stocked(env, (A, B), lineage=[1, 2, 3, 77])
    real, down = store.taxonomy, [True]

    def taxonomy(species, taxon_id=None):
        if down[0]:
            raise FetchError("not_found", "HTTP 400")
        return real(species, taxon_id=taxon_id)

    store.taxonomy = taxonomy
    row = fetched(env, A)
    assert row["state"] == READY and row["lineage"] == []
    assert env.ctx.library().get(row["id"]).lineage == []
    assert [r["accession"] for r in store.without_taxonomy()] == [A]
    assert store.repair_taxonomy() == 0                         # still no answer: nothing lost
    down[0] = False
    worker = env.worker()                                       # the worker does it by itself
    run_until(worker, lambda: not store.without_taxonomy())
    assert env.store.get_store_genome(A)["lineage"] == [1, 2, 3, 77]
    assert env.ctx.library().get(row["id"]).lineage == [1, 2, 3, 77]
    assert any("taxonomy of %s" % A in line for line in env.logs)


# ── files in the store's folder that are in no record ───────────────────────

def stray_folder(store, name="GCA_999999999.1", kilobytes=6):
    """What a lost database leaves behind: a genome folder nobody has a record of."""
    folder = store.root / name
    folder.mkdir(parents=True)
    (folder / "Old_genome.fna").write_bytes(b">s\n" + b"ACGT" * (250 * kilobytes))
    return folder


def test_files_in_no_record_count_against_the_size_of_the_store(make_env):
    """The store cannot know what they are, so it never deletes them by itself. It
    must still stay within the space it was given."""
    env = make_env(genome_store_gb=ROOM_FOR_THREE)
    store = stocked(env, (A, B, C, D))
    assert fetched(env, A)["state"] == READY
    stray = stray_folder(store)                                 # 6 kB
    usage = store.usage()
    assert usage["unknown_names"] == ["GCA_999999999.1"] and usage["unknown"] == 6003
    assert usage["free"] == usage["budget"] - usage["used"] - 6003
    # 13 kB in all: without the 6 kB nobody knows, a second genome would fit next to A.
    assert fetched(env, B)["state"] == READY
    assert env.store.get_store_genome(A)["state"] == EVICTED
    on_disk = sum(p.stat().st_size for p in store.root.rglob("*") if p.is_file())
    assert on_disk <= store.budget
    # Nothing the store does by itself removes them, and it says so once, not every round.
    store.housekeeping()
    store.housekeeping()
    assert (stray / "Old_genome.fna").is_file()
    told = [line for line in env.genome_logs if "NOT IN ANY RECORD" in line]
    assert len(told) == 1 and "GCA_999999999.1" in told[0] and "library sweep" in told[0]
    # When they take the room that is left, a request is refused before anything is fetched.
    stray_folder(store, "left_by_hand", kilobytes=5)
    with pytest.raises(StoreProblem) as refused:
        store.want(C)
    assert refused.value.code == "no_room"
    assert env.store.get_store_genome(C) is None
    # The operator's command deletes exactly these.
    assert sorted(name for name, _ in store.sweep()) == ["GCA_999999999.1", "left_by_hand"]
    assert not stray.exists() and store.usage()["unknown"] == 0
    assert (store.root / B / "Lasius_beta.fna").is_file()       # what is on record stays
    assert store.want(C)["state"] == WANTED


def test_file_the_store_did_not_fetch_is_not_taken_for_the_genome(make_env):
    """A folder with the genome's name and a sequence file in it, but no record that
    the store ever fetched it: the genome is fetched and checked like any other."""
    env = make_env(genome_store_gb=0.001)
    store = stocked(env)
    folder = store.root / A
    folder.mkdir(parents=True)
    (folder / "Formica_alpha.fna").write_bytes(b">something_else\nAAAA\n")
    assert fetched(env, A)["state"] == READY
    assert (folder / "Formica_alpha.fna").read_bytes() == env.ncbi.assemblies[A]["plain"]


def test_genomes_the_operator_registered_are_never_taken_for_strays(make_env):
    """The operator may keep registered genomes in the folder the store uses."""
    env = make_env(genome_store_gb=ROOM_FOR_THREE)
    env.settings.genome_store_dir = str(env.genomes_dir)        # where the library's files lie
    store = stocked(env)
    own = sorted(p.name for p in env.genomes_dir.iterdir())
    assert "Homus_primus.fna" in own and "Homus_primus.gff" in own
    stray_folder(store)
    assert store.usage()["unknown_names"] == ["GCA_999999999.1"]
    assert [name for name, _ in store.sweep()] == ["GCA_999999999.1"]
    assert sorted(p.name for p in env.genomes_dir.iterdir()) == own
    # While the library file cannot be read, nobody knows what is the operator's:
    # nothing is deleted then.
    stray_folder(store)
    env.library_file.write_text("{ this is not a library")
    with pytest.raises(StoreProblem) as refused:
        store.sweep()
    assert refused.value.code == "library_unreadable"
    assert (store.root / "GCA_999999999.1").is_dir()
    assert sorted(p.name for p in env.genomes_dir.iterdir()) == sorted(own + ["GCA_999999999.1"])


# ── through the API ─────────────────────────────────────────────────────────

@pytest.fixture
def store_app(make_env):
    from synvoy_web.api import App
    env = make_env(genome_store_gb=0.001, max_genome_mb=1000)
    env.beat()
    return App(env.ctx), env


def test_form_can_ask_what_ncbi_has(store_app):
    app, env = store_app
    env.ncbi.add(A, "Formica alpha", annotated=True)
    env.ncbi.add(B, "Formica alpha", category="", level="Scaffold")
    env.ncbi.add(C, "Magnus giganteus", reported=2000000000)
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": "Formica alpha"})
    assert status == 200 and data["downloads"] is True
    assert [c["accession"] for c in data["candidates"]] == [A]             # the reference genome
    first = data["candidates"][0]
    assert (first["species"], first["usable"], first["on_server"], first["annotated"],
            first["reference"]) == ("Formica alpha", True, None, True, True)
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": B.lower()})
    assert [c["accession"] for c in data["candidates"]] == [B]
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": B, "purpose": "home"})
    assert data["candidates"][0]["usable"] is False and "annotation" in data["candidates"][0]["problem"]
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": "Magnus giganteus"})
    assert data["candidates"][0]["usable"] is False and "2000 Mb" in data["candidates"][0]["problem"]
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": "Nomen nescio"})
    assert status == 200 and data["candidates"] == []
    for bad in ({"term": "<script>"}, {"term": ""}, {"term": "x" * 200}, {"name": "x"},
                {"term": "Formica", "purpose": "everything"}, ["Formica"]):
        assert call_json(app, "POST", "/api/v1/genomes/find", bad)[0] == 400
    env.ncbi.down = True
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": "Formica alpha"})
    assert status == 502 and data["error"]["code"] == "ncbi_unavailable"


def test_genome_already_on_the_server_is_pointed_to(store_app):
    app, env = store_app
    env.ncbi.add(A, "Formica alpha")
    row = fetched(env, A)
    _, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": A})
    assert data["candidates"][0]["on_server"] == row["id"] and data["candidates"][0]["usable"]
    _, listing = call_json(app, "GET", "/api/v1/library")
    entry = [g for g in listing["genomes"] if g["id"] == row["id"]][0]
    assert entry["origin"] == "ncbi" and entry["group"] == "ncbi" and entry["available"]
    assert "fasta" not in entry and str(env.genomes.root) not in json.dumps(listing)


def test_finding_genomes_is_limited_per_visitor(make_env):
    from synvoy_web.api import App
    env = make_env(genome_store_gb=0.001, max_lookups_per_hour=2)
    env.beat()
    app = App(env.ctx)
    env.ncbi.add(A, "Formica alpha")
    for _ in range(2):
        assert call_json(app, "POST", "/api/v1/genomes/find", {"term": "Formica alpha"})[0] == 200
    status, data = call_json(app, "POST", "/api/v1/genomes/find", {"term": "Formica alpha"})
    assert status == 429 and data["error"]["code"] == "too_many_lookups"


def test_status_says_whether_genomes_are_fetched(store_app, make_env):
    from synvoy_web.api import App
    with_store, _ = store_app
    limits = call_json(with_store, "GET", "/api/v1/status")[1]["limits"]
    assert limits["downloads"] is True and limits["max_fetches_per_day"] == 5
    without = make_env()
    without.beat()
    assert call_json(App(without.ctx), "GET", "/api/v1/status")[1]["limits"]["downloads"] is False


def test_only_related_genomes_are_ticked_by_default(store_app):
    """A server can hold genomes of quite different groups. For a starting genome the
    form ticks the close ones; the others are listed and can be ticked."""
    app, env = store_app
    env.ncbi.add(A, "Formica alpha", lineage=[1, 2, 3, 77])    # close to Homus primus [1, 2, 3, 99010]
    env.ncbi.add(B, "Saccharomyces remotus", lineage=[1, 50, 51, 52])
    for accession in (A, B):
        fetched(env, accession)
    _, data = call_json(app, "GET", "/api/v1/library", query="home=homus_primus")
    assert set(data["targets"]) >= {"gcf_000000001_1", "gca_000000002_1", "alpha_beta", "eta_theta"}
    assert "gcf_000000001_1" in data["related"] and "alpha_beta" in data["related"]
    assert "gca_000000002_1" not in data["related"]            # shares one level of four
    assert "eta_theta" not in data["related"]                  # no lineage on record
    assert data["related"] == [g for g in data["targets"] if g in set(data["related"])]   # same order


def test_species_list_can_be_ordered_for_a_species_that_is_not_here(store_app):
    app, env = store_app
    env.ncbi.lineages[777] = [1, 5, 6, 777]                    # closest to Epsilon zeta
    status, data = call_json(app, "GET", "/api/v1/library", query="taxon=777")
    assert status == 200 and data["targets"][0] == "epsilon_zeta"
    assert data["related"] == ["epsilon_zeta"]                 # only this one shares half the lineage
    assert set(data["targets"]) == {g["id"] for g in data["genomes"] if g["available"]}
    asked = len(env.ncbi.calls)
    call_json(app, "GET", "/api/v1/library", query="taxon=777")
    assert call_json(app, "GET", "/api/v1/library", query="taxon=abc")[0] == 400
    # An unknown taxon still gives the list, in name order.
    status, data = call_json(app, "GET", "/api/v1/library", query="taxon=31337")
    assert status == 200 and len(data["targets"]) == len([g for g in data["genomes"] if g["available"]])
    assert asked == len(env.ncbi.calls)
