"""Settings, the genome library, sequences and accessions: everything that comes in
from the operator or from a visitor before a search exists."""

import json
import os

import pytest

from conftest import PROTEIN, UNIPROT_JSON, FakeFetch, library_data, uniprot_replies
from synvoy_web import jobs
from synvoy_web.config import ConfigError, load_settings
from synvoy_web.db import Store, StoreError
from synvoy_web.library import (LibraryError, LibraryWatcher, fasta_stats, file_stem,
                                load_library, parse_library, save_library_data, strip_qualifier)
from synvoy_web.net import FetchError, http_get
from synvoy_web.query import (QueryError, Resolver, classify_accession, clean_protein,
                              safe_name)
from synvoy_web.taxonomy import lookup_taxonomy


# ── settings ────────────────────────────────────────────────────────────────

def test_defaults_are_valid():
    settings = load_settings(env={})
    assert settings.max_targets == 10 and settings.max_running == 1
    assert settings.email_enabled is False
    assert settings.library_path.name == "library.json"


def test_file_then_environment_then_overrides(tmp_path):
    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"_comment": "x", "max_targets": 5, "port": 9000}))
    env = {"SYNVOY_WEB_CONFIG": str(path), "SYNVOY_WEB_PORT": "9001",
           "SYNVOY_WEB_TRUST_PROXY": "true", "SYNVOY_WEB_CORS_ORIGINS": '["https://a.example"]',
           "SYNVOY_WEB_MIN_FREE_GB": "3"}
    settings = load_settings(env=env, overrides={"max_queue": 7})
    assert (settings.max_targets, settings.port, settings.max_queue) == (5, 9001, 7)
    assert settings.trust_proxy is True and settings.cors_origins == ["https://a.example"]
    assert settings.min_free_gb == 3.0


@pytest.mark.parametrize("data", [
    {"max_targetz": 5},                         # a typo must not pass silently
    {"max_targets": "5"}, {"max_targets": 0}, {"port": 70000}, {"accepting": "yes"},
    {"max_total_mb": 10, "max_genome_mb": 100}, {"cors_origins": ["example.org"]},
    {"cors_origins": ["https://example.org/"]}, {"nextflow_profile": "auto; rm -rf /"},
    {"smtp_host": "mail.example.org"},          # mail without sender and site address
    {"public_base_url": "example.org"}, {"pipeline_args": "x"}, {"max_runtime_hours": 0},
])
def test_wrong_settings_are_refused(tmp_path, data):
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ConfigError):
        load_settings(str(path), env={})


def test_broken_settings_file(tmp_path):
    path = tmp_path / "settings.json"
    path.write_text("{not json")
    with pytest.raises(ConfigError):
        load_settings(str(path), env={})
    with pytest.raises(ConfigError):
        load_settings(str(tmp_path / "missing.json"), env={})
    with pytest.raises(ConfigError):
        load_settings(env={"SYNVOY_WEB_PORT": "abc"})


def test_secret_is_made_once_and_kept_private(tmp_path):
    settings = load_settings(env={}, overrides={"data_dir": str(tmp_path)})
    first = settings.secret()
    assert len(first) >= 32 and settings.secret() == first
    assert (os.stat(str(tmp_path / "secret")).st_mode & 0o077) == 0


def test_database_from_a_newer_version_is_refused(tmp_path):
    store = Store(tmp_path / "x.sqlite3")
    store._conn().execute("PRAGMA user_version=99")
    store.close()
    with pytest.raises(StoreError):
        Store(tmp_path / "x.sqlite3")


def test_database_changes_roll_back_together(tmp_path):
    store = Store(tmp_path / "x.sqlite3")
    with pytest.raises(RuntimeError):
        with store.transaction():
            store.set_meta("a", 1)
            raise RuntimeError("stop")
    assert store.get_meta("a") is None
    with store.transaction():
        with store.transaction():               # nesting is allowed
            store.set_meta("a", 2)
    assert store.get_meta("a") == 2 and store.quick_check() == "ok"


# ── genome library ──────────────────────────────────────────────────────────

def test_library_loads_and_orders(tmp_path):
    library = parse_library(library_data(tmp_path), str(tmp_path))
    home = library.get("homus_primus")
    assert home.can_be_home and home.stem == "Homus_primus" and home.size_mb == 12.0
    assert [g.id for g in library.homes()] == ["homus_primus", "gamma_delta"]
    assert [g.id for g in library.targets_for(home)] == [
        "alpha_beta", "big_one", "gamma_delta", "epsilon_zeta", "eta_theta"]
    assert library.home_for(99010, "") is home
    assert library.home_for(None, "Homus primus (strain T1 / ATCC 1)") is home
    assert library.home_for(None, "Homus primus subsp. tertius") is home      # same binomial
    assert library.home_for(None, "Alpha beta") is None                        # no annotation
    assert library.home_for(None, "Homus") is None and library.home_for(None, "") is None
    assert "fasta" not in library.public()["genomes"][0]


def mutated(tmp_path, change):
    data = library_data(tmp_path)
    change(data)
    return data


@pytest.mark.parametrize("change", [
    lambda d: d["genomes"][1].update(id="homus_primus"),               # same id twice
    lambda d: d["genomes"][1].update(id="Alpha"),                      # upper case
    lambda d: d["genomes"][1].update(species="Alpha; rm -rf /"),
    lambda d: d["genomes"][1].update(species="A"),
    lambda d: d["genomes"][1].update(species="Homus-primus"),          # same file name as the home
    lambda d: d["genomes"][1].update(group="nobody"),
    lambda d: d["genomes"][1].update(size_bp=0),
    lambda d: d["genomes"][1].update(taxon_id="11"),
    lambda d: d["genomes"][1].update(lineage="1,2"),
    lambda d: d["genomes"][1].update(common_name="<script>"),
    lambda d: d["genomes"][1].update(surprise=1),
    lambda d: d["genomes"][1].pop("fasta"),
    lambda d: d["groups"].append({"id": "testers", "name": "again"}),
    lambda d: d.update(extra=1),
])
def test_mistakes_in_the_library_file_are_refused(tmp_path, change):
    with pytest.raises(LibraryError):
        parse_library(mutated(tmp_path, change), str(tmp_path), strict=False)


def test_missing_file_disables_one_genome_only(tmp_path):
    data = library_data(tmp_path)
    os.unlink(data["genomes"][1]["fasta"])
    with pytest.raises(LibraryError):
        parse_library(data, str(tmp_path), strict=True)
    library = parse_library(data, str(tmp_path), strict=False)
    assert library.get("alpha_beta").available is False
    assert library.get("homus_primus").available is True
    assert len(library.problems) == 1 and "does not exist" in library.problems[0]
    assert "alpha_beta" not in [g.id for g in library.targets_for(library.get("homus_primus"))]


def test_relative_paths_and_file_types(tmp_path):
    data = library_data(tmp_path)
    data["genomes"][0]["fasta"] = "Homus_primus.fna"
    assert parse_library(data, str(tmp_path)).get("homus_primus").fasta == str(
        tmp_path / "Homus_primus.fna")
    packed = tmp_path / "Alpha_beta.fna.gz"
    packed.write_text("x")
    data["genomes"][1]["fasta"] = str(packed)
    with pytest.raises(LibraryError):
        parse_library(data, str(tmp_path))


def test_file_names_made_from_species():
    assert file_stem("Apis mellifera") == "Apis_mellifera"
    assert file_stem("Candida sp. JCM-15000") == "Candida_sp_JCM_15000"
    assert strip_qualifier("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)") == \
        "Saccharomyces cerevisiae"
    assert strip_qualifier("Homo sapiens") == "Homo sapiens"


def test_fasta_statistics(tmp_path):
    path = tmp_path / "g.fna"
    path.write_text(">a\nACGT\nACGT\n>b desc\nAC\n\n>c\nACGTAC\n")
    assert fasta_stats(str(path)) == (16, 3, 8)
    empty = tmp_path / "e.fna"
    empty.write_text("")
    assert fasta_stats(str(empty)) == (0, 0, 0)


def test_library_is_reread_when_the_file_changes(tmp_path):
    path = tmp_path / "library.json"
    data = library_data(tmp_path)
    save_library_data(path, data)
    now = [0.0]
    watcher = LibraryWatcher(path, check_every=5.0, clock=lambda: now[0])
    assert len(watcher.current().genomes) == 6
    data["genomes"] = data["genomes"][:2]
    save_library_data(path, data)
    assert len(watcher.current().genomes) == 6          # not looked at again yet
    now[0] = 10.0
    assert len(watcher.current().genomes) == 2
    path.write_text("{broken")                          # an editing accident
    now[0] = 20.0
    assert len(watcher.current().genomes) == 2          # the last good library stays in use
    assert "not valid JSON" in watcher.error
    save_library_data(path, data)
    now[0] = 30.0
    assert watcher.current() and watcher.error == ""


def test_no_library_file_is_an_empty_library_and_a_vanished_one_is_a_fault(tmp_path):
    """Nothing registered yet is not a fault: genomes can come from the genome store
    alone. A file that was there and is gone is one."""
    path = tmp_path / "library.json"
    now = [0.0]
    watcher = LibraryWatcher(path, check_every=5.0, clock=lambda: now[0])
    assert watcher.current().genomes == [] and watcher.error == "" and watcher.missing
    with pytest.raises(LibraryError):
        load_library(path)                               # who asks for the file itself is told
    save_library_data(path, library_data(tmp_path))
    now[0] = 10.0
    assert len(watcher.current().genomes) == 6 and not watcher.missing
    path.unlink()
    now[0] = 20.0
    assert len(watcher.current().genomes) == 6           # the last good library stays in use
    assert "cannot read" in watcher.error


def test_taxonomy_lookup_reads_the_lineage():
    xml = (b"<TaxaSet><Taxon><TaxId>7460</TaxId><ScientificName>Apis mellifera</ScientificName>"
           b"<OtherNames><GenbankCommonName>honey bee</GenbankCommonName></OtherNames>"
           b"<LineageEx><Taxon><TaxId>131567</TaxId></Taxon><Taxon><TaxId>7459</TaxId></Taxon>"
           b"</LineageEx></Taxon></TaxaSet>")
    fetch = FakeFetch({"esearch": b'{"esearchresult":{"idlist":["7460"]}}', "efetch": xml})
    assert lookup_taxonomy("Apis mellifera", fetch=fetch) == (7460, [131567, 7459, 7460], "honey bee")
    with pytest.raises(FetchError):
        lookup_taxonomy("No such", fetch=FakeFetch({"esearch": b'{"esearchresult":{"idlist":[]}}'}))
    # NCBI's newer service is asked first. It counts from the root of all (1) and leaves
    # the taxon itself out; the lineage is the same list either way.
    newer = json.dumps({"taxonomy_nodes": [{"query": ["7460"], "taxonomy": {
        "tax_id": 7460, "organism_name": "Apis mellifera", "genbank_common_name": "honey bee",
        "lineage": [1, 131567, 7459]}}]}).encode()
    fetch = FakeFetch({"datasets/v2/taxonomy/taxon/7460": newer})
    assert lookup_taxonomy("", taxon_id=7460, fetch=fetch) == (7460, [131567, 7459, 7460], "honey bee")
    assert not any("efetch" in url for url in fetch.calls)
    unknown = json.dumps({"taxonomy_nodes": [{"errors": [{"reason": "not recognized"}]}]}).encode()
    fetch = FakeFetch({"datasets/v2/taxonomy/taxon/7460": unknown, "efetch": xml})
    assert lookup_taxonomy("", taxon_id=7460, fetch=fetch)[1] == [131567, 7459, 7460]


def test_taxonomy_lookup_asks_again_when_ncbi_stumbles():
    """The older service answers "HTTP 400" now and then for a request that works a
    second later (seen on 2026-10-11 for two of four yeasts)."""
    xml = (b"<TaxaSet><Taxon><TaxId>27288</TaxId><LineageEx><Taxon><TaxId>4893</TaxId></Taxon>"
           b"</LineageEx></Taxon></TaxaSet>")
    answers = [FetchError("not_found", "HTTP 400"), FetchError("unavailable", "HTTP 502"), xml]
    asked, pauses = [], []

    def fetch(url, timeout=15, max_bytes=2000000, accept="*/*"):
        asked.append(url)
        if "datasets" in url:
            raise FetchError("timeout", "no answer")
        answer = answers.pop(0)
        if isinstance(answer, Exception):
            raise answer
        return answer

    assert lookup_taxonomy("", taxon_id=27288, fetch=fetch, sleep=pauses.append) == (27288, [4893, 27288], "")
    assert len([u for u in asked if "efetch" in u]) == 3 and len(pauses) == 2

    def gone(url, timeout=15, max_bytes=2000000, accept="*/*"):
        asked.append(url)
        raise FetchError("not_found", "HTTP 404")

    del asked[:]
    with pytest.raises(FetchError):
        lookup_taxonomy("", taxon_id=27288, fetch=gone, sleep=pauses.append)
    assert len([u for u in asked if "efetch" in u]) == 1 and len(pauses) == 2     # not asked again


# ── sequences ───────────────────────────────────────────────────────────────

def test_sequence_is_cleaned_like_the_pipeline_does():
    text = ">sp|X|Y some protein\n  1 mkflv nvalv\n 11 FMVVY-ISYIY.\n" + PROTEIN[20:] + "*\n"
    header, sequence = clean_protein(text, 30, 3000)
    assert header == "sp|X|Y some protein" and sequence == PROTEIN
    assert clean_protein(PROTEIN.lower(), 30, 3000) == ("", PROTEIN)


@pytest.mark.parametrize("text, code", [
    ("", "sequence_missing"), (">only a header\n", "sequence_missing"), (None, "sequence_missing"),
    (PROTEIN[:29], "sequence_too_short"), (PROTEIN * 50, "sequence_too_long"),
    ("x" * 70000, "sequence_too_long"),
    (PROTEIN + "!", "sequence_invalid"), (PROTEIN[:30] + "*" + PROTEIN[30:], "sequence_invalid"),
    ("ACGTTGCAAC" * 10, "sequence_is_dna"), ("ACGUUGCAAC" * 10, "sequence_is_dna"),
    (">a\n" + PROTEIN + "\n>b\n" + PROTEIN, "several_sequences"),
])
def test_sequences_that_are_refused(text, code):
    with pytest.raises(QueryError) as caught:
        clean_protein(text, 30, 3000)
    assert caught.value.code == code and caught.value.message


def test_a_glycine_serine_peptide_is_not_dna():
    assert clean_protein("GSGSGSGSGS" * 4 + "AAAA", 30, 3000)[1].startswith("GSGS")


def test_alphabets_agree_with_the_pipeline():
    """The site must accept exactly the residues ``bin/normalize_query.py`` accepts."""
    import importlib.util
    from conftest import REPO_DIR
    from synvoy_web import query
    spec = importlib.util.spec_from_file_location(
        "normalize_query", str(REPO_DIR / "bin" / "normalize_query.py"))
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        pytest.skip("bin/normalize_query.py cannot be imported here: %s" % exc)
    assert set(query.PROTEIN_ALPHABET) == set(module.PROTEIN_ALPHABET)
    assert set(query.NUC_ALPHABET) == set(module.NUC_ALPHABET)
    assert set(query._FORMATTING) == set(module._FORMATTING)
    for text in (PROTEIN, "ACGTTGCAAC" * 10, "GSGSGSGS", "ACGTNNNNRY"):
        assert query.is_nucleotide(text) == module.is_nucleotide(text)


def test_names_are_made_safe():
    assert safe_name("my protein <script>") == "my_protein_script"
    assert safe_name("../../etc/passwd") == "etc_passwd"
    assert safe_name("") == "query" and safe_name("---") == "query"
    assert len(safe_name("a" * 100)) == 40


# ── accessions ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("text, expected", [
    ("P01501", ("uniprot", "P01501")), (" p01501.2 ", ("uniprot", "P01501")),
    ("A0A6M3Z554", ("uniprot", "A0A6M3Z554")), ("D6VTK4", ("uniprot", "D6VTK4")),
    ("XP_031329057.1", ("ncbi", "XP_031329057.1")), ("xp_031329057", ("ncbi", "XP_031329057")),
    ("NP_001011617.1", ("ncbi", "NP_001011617.1")), ("KAF1234567.1", ("ncbi", "KAF1234567.1")),
])
def test_accessions_are_recognised(text, expected):
    assert classify_accession(text) == expected


@pytest.mark.parametrize("text", [
    "", "   ", "melittin", "P0150", "https://evil.example/x", "P01501/../../x", "P01501?x=1",
    "XP_", "A" * 30, None, 7, "P01501 P01502",
])
def test_everything_else_is_not_an_accession(text):
    with pytest.raises(QueryError):
        classify_accession(text)


def resolver(tmp_path, replies):
    store = Store(tmp_path / "r.sqlite3")
    fetch = FakeFetch(replies)
    return Resolver(store, fetch=fetch), fetch


def test_uniprot_record_is_read_and_kept(tmp_path):
    res, fetch = resolver(tmp_path, uniprot_replies())
    record, cached = res.resolve("P01501")
    assert cached is False and record["sequence"] == PROTEIN and record["gene"] == "MELT"
    assert record["organism"] == "Homus primus (strain T1)" and record["taxon_id"] == 99010
    assert record["header"].startswith("sp|P01501|MEL_TEST") and record["reviewed"] is True
    assert all(url.startswith("https://rest.uniprot.org/uniprotkb/P01501.") for url in fetch.calls)
    assert res.resolve("p01501")[1] is True and len(fetch.calls) == 2


def test_uniprot_answers_that_cannot_be_used(tmp_path):
    inactive = dict(UNIPROT_JSON, entryType="Inactive")
    res, _ = resolver(tmp_path, uniprot_replies(data=inactive))
    with pytest.raises(QueryError) as caught:
        res.resolve("P01501")
    assert caught.value.code == "accession_inactive"

    for replies, code, status in (
            ({}, "accession_not_found", 404),
            ({"P01501.json": FetchError("timeout")}, "lookup_unavailable", 503),
            ({"P01501.json": FetchError("unavailable")}, "lookup_unavailable", 503),
            ({"P01501.json": b"<html>maintenance</html>"}, "lookup_unavailable", 503),
            (uniprot_replies(fasta=">x\n" + PROTEIN[:-1] + "A\n"), "lookup_unavailable", 503),
            (uniprot_replies(fasta="no fasta at all"), "lookup_unavailable", 503)):
        res, _ = resolver(_fresh(tmp_path), replies)
        with pytest.raises(QueryError) as caught:
            res.resolve("P01501")
        assert (caught.value.code, caught.value.status) == (code, status)
        assert res.cached("P01501") is None         # a failure is not remembered


_counter = [0]


def _fresh(tmp_path):
    _counter[0] += 1
    folder = tmp_path / ("r%d" % _counter[0])
    folder.mkdir()
    return folder


def test_ncbi_record_is_read(tmp_path):
    summary = {"result": {"uids": ["123"], "123": {
        "title": "luciferin 4-monooxygenase [Photinus pyralis]", "organism": "Photinus pyralis",
        "taxid": 7054, "accessionversion": "XP_031329057.1"}}}
    replies = {"esummary": json.dumps(summary).encode("utf-8"),
               "efetch": (">XP_031329057.1 luciferin 4-monooxygenase [Photinus pyralis]\n"
                          + PROTEIN + "\n").encode("utf-8")}
    res, fetch = resolver(tmp_path, replies)
    record, _ = res.resolve("XP_031329057")
    assert record["source"] == "ncbi" and record["accession"] == "XP_031329057.1"
    assert record["name"] == "luciferin 4-monooxygenase" and record["organism"] == "Photinus pyralis"
    assert record["taxon_id"] == 7054 and record["length"] == len(PROTEIN)
    assert all("eutils.ncbi.nlm.nih.gov" in url for url in fetch.calls)

    res, _ = resolver(_fresh(tmp_path), {"esummary": b'{"result":{"uids":[]}}'})
    with pytest.raises(QueryError) as caught:
        res.resolve("XP_000000001")
    assert caught.value.code == "accession_not_found"


def test_requests_go_to_the_two_databases_only():
    for url in ("http://rest.uniprot.org/x", "https://example.org/x", "https://169.254.169.254/",
                "file:///etc/passwd", "https://rest.uniprot.org.evil.example/x"):
        with pytest.raises(FetchError) as caught:
            http_get(url)
        assert caught.value.kind == "refused"


# ── the cache key ───────────────────────────────────────────────────────────

def test_cache_key_follows_what_decides_the_result(env, tmp_path):
    env.beat()
    request, _ = jobs.build_request(env.ctx, env.payload(), "c")
    version = {"fingerprint": "abc"}
    key = jobs.cache_key(version, request, env.settings)
    assert key == jobs.cache_key(version, json.loads(json.dumps(request)), env.settings)
    assert key != jobs.cache_key({"fingerprint": "abd"}, request, env.settings)
    other, _ = jobs.build_request(env.ctx, env.payload(name="other"), "c")
    assert key != jobs.cache_key(version, other, env.settings)
    fewer, _ = jobs.build_request(env.ctx, env.payload(targets=["alpha_beta"]), "c")
    assert key != jobs.cache_key(version, fewer, env.settings)
    profile = load_settings(env={}, overrides=dict(env.settings.as_dict(), nextflow_profile="standard"))
    assert key != jobs.cache_key(version, request, profile)


def test_code_fingerprint_changes_with_the_code(tmp_path):
    repo = tmp_path / "repo"
    (repo / "bin").mkdir(parents=True)
    (repo / "VERSION").write_text("1.0\n")
    (repo / "bin" / "tool.py").write_text("print(1)\n")
    (repo / "docs").mkdir()
    (repo / "docs" / "notes.md").write_text("x")
    first = jobs.code_fingerprint(repo)
    (repo / "docs" / "notes.md").write_text("changed")           # not code
    assert jobs.code_fingerprint(repo) == first
    (repo / "bin" / "tool.py").write_text("print(2)\n")
    assert jobs.code_fingerprint(repo) != first
    version = jobs.pipeline_version(repo)
    assert version["version"] == "1.0" and len(version["fingerprint"]) == 16


@pytest.mark.parametrize("tag, status, label, expected", [
    ("v1.0.0", "", "", ("1.0.0", True)),          # the checkout is the tagged release
    ("1.0.0", "", "", ("1.0.0", True)),
    ("", "", "", ("pre-release", False)),         # no tag here: the number says nothing
    ("v0.9.0", "", "", ("pre-release", False)),   # another release than VERSION names
    ("v1.0.0", " M bin/x.py", "", ("pre-release", False)),      # changed since the tag
    (None, None, "", ("pre-release", False)),     # no git history at all
    (None, None, "1.0.0 (archive)", ("1.0.0 (archive)", False)),   # the operator says so
])
def test_version_number_is_shown_for_a_tagged_release_only(tmp_path, monkeypatch, tag, status,
                                                         label, expected):
    """VERSION alone is a number in a file. What runs is a release only when the
    checkout stands on that release's tag and has not been changed since."""
    (tmp_path / "VERSION").write_text("1.0.0\n")
    answers = {"rev-parse": None if tag is None else "abc1234def", "status": status,
               "describe": tag or None}
    monkeypatch.setattr(jobs, "_git", lambda repo, command, *rest: answers[command])
    version = jobs.pipeline_version(tmp_path, label)
    assert (version["label"], version["release"]) == expected
    assert jobs.shown_version(version)["label"] == expected[0]
    assert "version" not in jobs.shown_version(version)       # the bare number is not shown
    # A search stored before releases were told apart was not one.
    assert jobs.shown_version({"version": "2.0.0", "commit": "abc"}) == {
        "label": "pre-release", "release": False, "commit": "abc", "modified": False}


def test_search_ids_cannot_be_guessed_or_misread():
    ids = {jobs.new_job_id() for _ in range(200)}
    assert len(ids) == 200 and all(jobs.JOB_ID_RE.match(i) for i in ids)
    assert all(not i.startswith(("-", "_")) for i in ids)
