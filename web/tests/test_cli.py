"""The operator's commands: genome library, checks, searches."""

import json
import sys
import time

import pytest

from conftest import FAKE_PIPELINE, PROTEIN, REPO_DIR, WEB_DIR, library_data
from synvoy_web import cli, db, jobs
from synvoy_web.config import load_settings
from synvoy_web.server import build_context


@pytest.fixture
def config(tmp_path):
    """A settings file for an empty installation that runs the stand-in pipeline."""
    data = tmp_path / "var"
    data.mkdir()
    path = tmp_path / "settings.json"
    path.write_text(json.dumps({
        "data_dir": str(data), "site_dir": str(WEB_DIR / "site"), "repo_dir": str(REPO_DIR),
        "conda_env": "", "pipeline_command": [sys.executable, str(FAKE_PIPELINE)],
        "min_free_gb": 0.0, "disk_factor": 0.0}))
    return path


def run(config, *argv):
    return cli.main(["--config", str(config)] + list(argv))


def genome_file(tmp_path, name, text=">chr1\nACGTACGT\n>chr2\nACGT\n"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_library_add_list_check_remove(config, tmp_path, capsys):
    fasta = genome_file(tmp_path, "Homus_primus.fna")
    gff = genome_file(tmp_path, "Homus_primus.gff", "##gff-version 3\n")
    assert run(config, "library", "add", "--species", "Homus primus", "--fasta", fasta, "--gff", gff,
               "--group", "testers", "--group-name", "Test animals", "--assembly", "GCF_1.1",
               "--lineage", "1,2,3", "--taxon-id", "3", "--common-name", "test animal") == 0
    assert run(config, "library", "add", "--species", "Alpha beta", "--no-network", "--group",
               "testers", "--fasta", genome_file(tmp_path, "Alpha_beta.fa")) == 0
    settings = load_settings(str(config), env={})
    data = json.loads(settings.library_path.read_text())
    first = data["genomes"][0]
    assert (first["id"], first["size_bp"], first["n_sequences"], first["n50"]) == ("homus_primus", 12, 2, 8)
    assert first["checksum"].startswith("sha256:") and first["lineage"] == [1, 2, 3]
    assert data["genomes"][1]["gff"] == "" and data["genomes"][1]["lineage"] == []
    capsys.readouterr()
    assert run(config, "library", "list") == 0
    listing = capsys.readouterr().out
    assert "home+target" in listing and "target only" in listing and "Test animals" in listing
    assert run(config, "library", "check") == 0
    assert run(config, "library", "remove", "alpha_beta") == 0
    assert run(config, "library", "remove", "alpha_beta") == 1
    assert len(json.loads(settings.library_path.read_text())["genomes"]) == 1
    assert (tmp_path / "Alpha_beta.fa").is_file()            # files are never touched


def test_library_add_refuses_what_would_break_the_library(config, tmp_path, capsys):
    fasta = genome_file(tmp_path, "g.fna")
    base = ["library", "add", "--no-network", "--group", "testers", "--fasta", fasta]
    assert run(config, *base, "--species", "Homus primus") == 1                 # group not there yet
    assert run(config, *base, "--species", "Homus primus", "--group-name", "Testers") == 0
    before = load_settings(str(config), env={}).library_path.read_text()
    assert run(config, *base, "--species", "Homus primus") == 1                 # same id again
    assert run(config, *base, "--species", "Homus-primus", "--id", "other") == 1  # same file name
    assert run(config, *base, "--species", "x; rm -rf /") == 1
    assert run(config, *base, "--species", "Beta gamma", "--id", "Bad Id") == 1
    assert run(config, "library", "add", "--no-network", "--group", "testers", "--species",
               "Beta gamma", "--fasta", str(tmp_path / "missing.fna")) == 1
    assert run(config, *base[:-1], genome_file(tmp_path, "empty.fna", ""), "--species", "Beta gamma") == 1
    assert run(config, *base[:-1], genome_file(tmp_path, "packed.fna.gz", "x"), "--species",
               "Beta gamma") == 1
    assert load_settings(str(config), env={}).library_path.read_text() == before   # left as it was
    assert run(config, *base, "--species", "Homus primus", "--replace") == 0
    assert "ERROR" in capsys.readouterr().out


def test_check_reports_what_is_missing(config, tmp_path, capsys):
    assert run(config, "check") == 1                          # no genome library yet
    assert "[FAIL] genome library" in capsys.readouterr().out
    settings = load_settings(str(config), env={})
    settings.library_path.write_text(json.dumps(library_data(tmp_path)))
    assert run(config, "check") == 0
    out = capsys.readouterr().out
    assert "[ok]   a search can be started" in out and "[ok]   database" in out
    assert "Magnus maximus (1500 Mb) is above max_genome_mb" in out
    assert "Eta theta has no taxonomy" in out


def test_searches_can_be_listed_and_cancelled(config, tmp_path, capsys):
    settings = load_settings(str(config), env={})
    settings.library_path.write_text(json.dumps(library_data(tmp_path)))
    ctx = build_context(settings)
    ctx.store.set_meta("worker_heartbeat", int(time.time()))
    outcome = jobs.submit(ctx, {"query": {"sequence": PROTEIN, "name": "mine", "home": "homus_primus"},
                                "targets": ["alpha_beta"]}, "10.0.0.1")
    capsys.readouterr()
    assert run(config, "jobs", "list") == 0
    assert outcome["id"] in capsys.readouterr().out
    assert run(config, "jobs", "show", outcome["id"]) == 0
    shown = capsys.readouterr().out
    assert '"label": "mine"' in shown and PROTEIN[:30] not in shown
    assert run(config, "jobs", "cancel", outcome["id"]) == 0
    assert ctx.store.get_job(outcome["id"])["state"] == db.CANCELLED
    assert run(config, "jobs", "cancel", outcome["id"]) == 0  # nothing to do, no error
    assert run(config, "jobs", "show", "Aa0000000000000000000a") == 1


def test_wrong_use_is_reported(config, tmp_path, capsys):
    assert cli.main([]) == 2
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"max_targetz": 1}))
    assert cli.main(["--config", str(bad), "check"]) == 2
    assert "unknown setting" in capsys.readouterr().out
    assert run(config, "gallery", "build", "--spec", str(tmp_path / "none.json")) == 1


# ── the genome store ────────────────────────────────────────────────────────

@pytest.fixture
def store_cli(make_env, monkeypatch, tmp_path):
    """The command line on an installation whose NCBI is the stand-in."""
    env = make_env(genome_store_gb=13000 / 1e9)
    monkeypatch.setattr(cli, "build_context", lambda settings: env.ctx)
    monkeypatch.setattr(cli.time, "sleep", lambda seconds: None)
    config = tmp_path / "store.json"
    config.write_text(json.dumps({"data_dir": str(env.data_dir)}))
    env.ncbi.add("GCF_000000001.1", "Formica alpha", annotated=True)
    env.ncbi.add("GCA_000000002.1", "Lasius beta")
    env.ncbi.add("GCA_000000003.1", "Atta gamma")
    env.ncbi.add("GCA_000000004.1", "Apis delta")
    return env, lambda *argv: cli.main(["--config", str(config), "library"] + list(argv))


def test_operator_fills_the_store_by_name_or_accession(store_cli, capsys):
    env, library = store_cli
    assert library("fetch", "Formica alpha", "GCA_000000002.1", "--annotation") == 1
    out = capsys.readouterr().out
    # By name: the reference genome NCBI names. Lasius has no annotation at NCBI.
    assert "Formica alpha: GCF_000000001.1 (ASM1, Chromosome" in out and "reference genome" in out
    assert "GCF_000000001.1: on the server (Formica alpha, 0 Mb, with annotation)" in out
    assert "GCA_000000002.1: NCBI has no gene annotation" in out
    assert env.store.get_store_genome("GCF_000000001.1")["origin"] == "operator"
    assert env.ctx.library().get("gcf_000000001_1").can_be_home
    assert library("fetch", "GCA_000000002.1", "--pin") == 0
    assert env.store.get_store_genome("GCA_000000002.1")["pinned"] is True
    assert library("fetch", "Nomen nescio") == 1
    assert "NCBI has no assembly under this name" in capsys.readouterr().out
    assert library("fetch") == 2


def test_store_can_be_filled_to_a_share_of_its_size(store_cli, tmp_path, capsys):
    env, library = store_cli
    wanted = tmp_path / "genomes.txt"
    wanted.write_text("GCF_000000001.1   # ant\n\nGCA_000000002.1\nGCA_000000003.1\nGCA_000000004.1\n")
    assert library("fetch", "--from-file", str(wanted), "--fill-to", "50") == 0
    out = capsys.readouterr().out
    assert "which is the --fill-to mark. Stopping before GCA_000000003.1" in out
    states = {r["accession"]: r["state"] for r in env.store.store_genomes()}
    assert states == {"GCF_000000001.1": "ready", "GCA_000000002.1": "ready"}


def test_status_lists_space_uses_and_what_goes_first(store_cli, capsys):
    env, library = store_cli
    library("fetch", "GCF_000000001.1", "GCA_000000002.1", "GCA_000000003.1")
    env.genomes.record_use(["gca_000000002_1", "gcf_000000001_1"])     # one search
    env.genomes.record_use(["gca_000000002_1"])                        # another
    assert library("pin", "GCA_000000003.1") == 0
    capsys.readouterr()
    assert library("status") == 0
    out = capsys.readouterr().out
    assert "3 genome(s), 1 pinned" in out and "of 0.00 GB used" in out
    line = [l for l in out.splitlines() if "GCA_000000002.1" in l and "Lasius" in l][0]
    assert " ready " in line and "     2  " in line
    # The pinned one is not on the list; the unused... there is none: least used first.
    assert "first to go when room is needed: GCF_000000001.1 (Formica alpha), GCA_000000002.1" in out
    assert library("unpin", "GCA_000000003.1") == 0
    assert library("pin", "GCA_999999999.1") == 1


def test_operator_can_remove_a_genome_unless_a_search_needs_it(store_cli, capsys):
    env, library = store_cli
    library("fetch", "GCF_000000001.1", "GCA_000000002.1")
    env.beat()
    from synvoy_web import jobs
    jobs.submit(env.ctx, env.payload(targets=["gca_000000002_1"]), "10.0.0.1")
    assert library("evict", "GCA_000000002.1") == 1
    assert "needed by a search" in capsys.readouterr().out
    assert library("evict", "GCF_000000001.1") == 0
    assert env.store.get_store_genome("GCF_000000001.1")["state"] == "evicted"
    assert library("evict", "GCF_000000001.1") == 1            # no longer here


def test_sweep_lists_first_and_deletes_only_when_told(store_cli, capsys):
    env, library = store_cli
    library("fetch", "GCF_000000001.1")
    stray = env.genomes.root / "GCA_999999999.1"
    stray.mkdir()
    (stray / "Old.fna").write_bytes(b">s\n" + b"ACGT" * 1000)
    capsys.readouterr()
    assert library("status") == 0
    out = capsys.readouterr().out
    assert "NOT IN ANY RECORD: 1 item(s)" in out and "GCA_999999999.1" in out
    assert library("sweep") == 0
    out = capsys.readouterr().out
    assert "GCA_999999999.1" in out and "`library sweep --yes` deletes exactly these" in out
    assert stray.is_dir()                                      # listed, not deleted
    assert library("sweep", "--yes") == 0
    assert "Deleted 1 item(s)" in capsys.readouterr().out
    assert not stray.exists()
    assert (env.genomes.root / "GCF_000000001.1").is_dir()     # the genome on record stays
    assert library("sweep", "--yes") == 0
    assert "Nothing to delete" in capsys.readouterr().out


def test_check_says_how_the_genome_store_stands(config, tmp_path, capsys):
    run(config, "library", "add", "--species", "Homus primus", "--group", "testers",
        "--fasta", genome_file(tmp_path, "h.fna"), "--gff", genome_file(tmp_path, "h.gff", "##gff-version 3\n"),
        "--no-network")
    capsys.readouterr()
    run(config, "check")
    assert "genome store: off" in capsys.readouterr().out
    settings = json.loads(config.read_text())
    settings["genome_store_gb"] = 10.0 ** 7                    # more than any disk has
    config.write_text(json.dumps(settings))
    run(config, "check")
    out = capsys.readouterr().out
    assert "genome store: 0.0 of 10000000.0 GB used, 0 genome(s)" in out
    assert "lower genome_store_gb or free the disk" in out
    assert "in no record" not in out
    left = tmp_path / "var" / "genomes" / "GCA_999999999.1"
    left.mkdir(parents=True)
    (left / "Old.fna").write_text(">s\nACGT\n")
    run(config, "check")
    assert "1 item(s), 0.0 GB, in the store's folder are in no record" in capsys.readouterr().out


def test_check_accepts_an_installation_that_has_only_the_genome_store(config, capsys):
    assert run(config, "check") == 1                     # no library and no store: nothing to search
    out = capsys.readouterr().out
    assert "[FAIL] genome library: no file" in out and "genome_store_gb" in out
    settings = json.loads(config.read_text())
    settings["genome_store_gb"] = 0.001
    config.write_text(json.dumps(settings))
    assert run(config, "check") == 0
    out = capsys.readouterr().out
    assert "[FAIL]" not in out
    assert "[ok]   genome library: no file" in out
    assert "this server has the genomes its store fetches from NCBI, and no others" in out
    assert "no genome on this server can be a starting genome yet" in out


def test_store_commands_say_when_the_store_is_off(make_env, monkeypatch, tmp_path, capsys):
    env = make_env()
    monkeypatch.setattr(cli, "build_context", lambda settings: env.ctx)
    config = tmp_path / "off.json"
    config.write_text(json.dumps({"data_dir": str(env.data_dir)}))
    assert cli.main(["--config", str(config), "library", "fetch", "GCF_000000001.1"]) == 1
    assert cli.main(["--config", str(config), "library", "status"]) == 0
    assert "genome store: off" in capsys.readouterr().out
