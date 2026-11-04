"""Fix-phase d7 (2 Nov; item 23): a value that can change a result says where it came from when the network decides it.

- ANNOTATE_GOI's UniProt lookup (its names can decide which home gene is the GOI): goi_info.json records whether UniProt
  answered (uniprot_lookup), and the report warns when it did not. Absent when no lookup was tried (no UniProt accession).
- The figures' common names: each with its source (table, cache, ncbi, none, unreachable, offline) in
  <base>_species_names.tsv; a lookup NCBI did not answer is a warning.
- PHYLO_SORT's NCBI lookups: a failure is recorded in sort_provenance.json, and the report says how the genomes were
  ordered (genome_order) and warns when NCBI did not answer.
- Easy mode's download: download_record.tsv names a genome or a GFF that did not arrive; the report warns (downloads).
"""
import json
import os
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)


# --------------------------------------------------------------------------- UniProt

def test_a_failed_uniprot_lookup_is_recorded(monkeypatch):
    import annotate_goi_exons as ag
    monkeypatch.setattr(ag, "UNIPROT_LOOKUP", {})

    def offline(*a, **k):
        raise OSError("Network is unreachable")
    monkeypatch.setattr(ag.urllib.request, "urlopen", offline)
    assert ag.fetch_uniprot_names("P01501") == []
    assert ag.UNIPROT_LOOKUP == {"query_id": "P01501", "status": "failed", "error": "Network is unreachable"}


def test_no_lookup_for_an_id_that_is_no_accession(monkeypatch):
    import annotate_goi_exons as ag
    monkeypatch.setattr(ag, "UNIPROT_LOOKUP", {})
    assert ag.fetch_uniprot_names("orf19.1562") == [] and ag.UNIPROT_LOOKUP == {}


def test_an_answer_is_recorded_with_its_names(monkeypatch):
    import io
    import annotate_goi_exons as ag
    monkeypatch.setattr(ag, "UNIPROT_LOOKUP", {})
    body = json.dumps({"genes": [{"geneName": {"value": "MELT"}}],
                       "proteinDescription": {"recommendedName": {"fullName": {"value": "Melittin"}}}}).encode()

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    monkeypatch.setattr(ag.urllib.request, "urlopen", lambda *a, **k: Resp(body))
    assert ag.fetch_uniprot_names("P01501") == ["MELT", "Melittin"]
    assert ag.UNIPROT_LOOKUP["status"] == "ok" and ag.UNIPROT_LOOKUP["names"] == ["MELT", "Melittin"]


def _report(tmp_path, *extra):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True, exist_ok=True)
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n")
    out = tmp_path / "report.json"
    p = subprocess.run([sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
                        "--output", str(out), "--allow-empty", *extra], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return json.loads(out.read_text())


def test_the_report_warns_when_uniprot_did_not_answer(tmp_path):
    gi = tmp_path / "goi_info.json"
    gi.write_text(json.dumps({"uniprot_lookup": {"query_id": "P01501", "status": "failed", "error": "timed out"}}))
    rep = _report(tmp_path, "--goi_info", str(gi))
    assert any("UniProt did not answer for P01501 (timed out)" in w for w in rep["warnings"])
    gi.write_text(json.dumps({"uniprot_lookup": {"query_id": "P01501", "status": "ok", "names": []}}))
    assert not any("UniProt" in w for w in _report(tmp_path, "--goi_info", str(gi))["warnings"])


# --------------------------------------------------------------------------- common names

@pytest.fixture
def taxa(monkeypatch, tmp_path):
    import synvoy_taxa as st
    for name in ("_OVERRIDES", "_CACHE"):
        monkeypatch.setattr(st, name, {})
    monkeypatch.setattr(st, "_NEGATIVE", set())
    monkeypatch.setattr(st, "SOURCES", {})
    (tmp_path / "cache.tsv").write_text("Apis mellifera\thoney bee\n")
    (tmp_path / "table.tsv").write_text("Bombus terrestris\tbuff-tailed bumblebee\n")
    monkeypatch.setattr(st.shutil, "which", lambda b: "/usr/bin/datasets")
    st.init_lookup(common_names_tsv=str(tmp_path / "table.tsv"), cache_path=str(tmp_path / "cache.tsv"))
    return st


def test_each_name_has_its_source(taxa, monkeypatch, tmp_path):
    def run(cmd, *a, **k):
        name = cmd[4]
        if name == "Candida albicans":
            raise subprocess.TimeoutExpired(cmd, 8)
        if name == "Nomen nudum":
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")
        return types.SimpleNamespace(returncode=0, stderr="", stdout=json.dumps(
            {"taxonomy": {"curator_common_name": "fruit fly"}}) + "\n")
    monkeypatch.setattr(taxa.subprocess, "run", run)
    for sci in ("Apis mellifera", "Bombus terrestris", "Drosophila melanogaster", "Candida albicans", "Nomen nudum"):
        taxa.common_name(sci)
    assert taxa.SOURCES == {"Apis mellifera": ("cache", "honey bee"),
                            "Bombus terrestris": ("table", "buff-tailed bumblebee"),
                            "Drosophila melanogaster": ("ncbi", "fruit fly"),
                            "Candida albicans": ("unreachable", ""), "Nomen nudum": ("none", "")}
    taxa.write_sources(str(tmp_path / "x_species_names.tsv"))
    rows = (tmp_path / "x_species_names.tsv").read_text().splitlines()
    assert rows[1] == "scientific_name\tcommon_name\tsource" and "Candida albicans\t\tunreachable" in rows


def test_offline_names_say_so(taxa, monkeypatch):
    taxa.init_lookup(allow_network=False, cache_path="/nonexistent/cache.tsv")
    taxa.common_name("Candida albicans")
    assert taxa.SOURCES["Candida albicans"] == ("offline", "")


def test_the_plot_writes_the_sources_and_warns(taxa, monkeypatch, tmp_path, capsys):
    import plot_synteny as ps
    monkeypatch.setattr(taxa.subprocess, "run", lambda *a, **k: (_ for _ in ()).throw(OSError("down")))
    taxa.common_name("Candida albicans")
    monkeypatch.setitem(ps.__dict__, "_synvoy_taxa", taxa)
    ps._write_species_name_sources(str(tmp_path / "synteny_block_locus_1_synteny_plot.html"))
    assert (tmp_path / "synteny_block_locus_1_species_names.tsv").exists()
    assert "NCBI did not answer for the common name of Candida albicans" in capsys.readouterr().err


def test_the_module_publishes_the_sources():
    nf = open(os.path.join(ROOT, "modules", "plot_synteny.nf")).read()
    assert 'path "*_species_names.tsv", emit: species_names, optional: true' in nf


# --------------------------------------------------------------------------- PHYLO_SORT

def test_a_failed_datasets_lookup_is_in_the_provenance(monkeypatch):
    import phylo_sort as ps
    monkeypatch.setattr(ps, "PROVENANCE", {})
    monkeypatch.setattr(ps.subprocess, "run",
                        lambda *a, **k: types.SimpleNamespace(returncode=1, stdout="", stderr="Network unreachable"))
    assert ps._fetch_lineage_from_accessions(["GCA_000001.1"]) == {}
    assert ps.PROVENANCE["network_failures"] == ["datasets summary genome: exit 1: Network unreachable"]


def test_the_report_says_how_the_genomes_were_ordered(tmp_path):
    d = tmp_path / "sp" / "1"
    d.mkdir(parents=True)
    (d / "sort_provenance.json").write_text(json.dumps({"strategy": "minhash", "network_failures": ["datasets x: timed out"]}))
    (tmp_path / "sp" / "2").mkdir()
    (tmp_path / "sp" / "2" / "sort_provenance.json").write_text(json.dumps({"strategy": "ete3_taxonomy"}))
    rep = _report(tmp_path, "--sort_provenance_dir", str(tmp_path / "sp"))
    assert rep["genome_order"] == [{"strategy": "ete3_taxonomy"},
                                   {"strategy": "minhash", "network_failures": ["datasets x: timed out"]}]
    assert any("NCBI did not answer while the target genomes were ordered" in w and "ordered by minhash" in w
               for w in rep["warnings"])


# --------------------------------------------------------------------------- downloads

def test_the_download_record(monkeypatch, tmp_path, capsys):
    import fetch_related_genomes as fr
    monkeypatch.setattr(fr, "DOWNLOAD_RECORD", [])
    monkeypatch.setattr(fr, "_HAS_ROBUST_DOWNLOADER", True)
    gff = tmp_path / "GCA_2.1.gff"
    gff.write_text("##gff-version 3\nc\tx\tgene\t1\t9\t.\t+\t.\tID=g\n")

    def robust(acc, output_path, fna_name, gff_name, max_retries):
        return {"GCA_1.1": (None, None), "GCA_2.1": (tmp_path / fna_name, gff), "GCA_3.1": (tmp_path / fna_name, None)}[acc]
    monkeypatch.setattr(fr, "_download_robust", robust, raising=False)
    for acc, sp in (("GCA_1.1", "A a"), ("GCA_2.1", "B b"), ("GCA_3.1", "C c")):
        fr._download_one(acc, str(tmp_path), sp)
    fr.write_download_record(str(tmp_path))
    rows = (tmp_path / "download_record.tsv").read_text().splitlines()
    assert rows == ["accession\tspecies\tgenome\tgff", "GCA_1.1\tA a\tfailed\t", "GCA_2.1\tB b\tok\tok",
                    "GCA_3.1\tC c\tok\tnone"]
    err = capsys.readouterr().err
    assert "1 chosen genome(s) could not be downloaded" in err and "GCA_3.1" in err
    rep = _report(tmp_path, "--download_record", str(tmp_path / "download_record.tsv"))
    assert [d["genome"] for d in rep["downloads"]] == ["failed", "ok", "ok"]
    assert any("could not be downloaded and were not searched: GCA_1.1" in w for w in rep["warnings"])
    assert any("came without a GFF: GCA_3.1" in w for w in rep["warnings"])


def test_a_pro_mode_report_has_no_download_record(tmp_path):
    rep = _report(tmp_path, "--download_record", str(tmp_path / "staged_genomes" / "download_record.tsv"))
    assert "downloads" not in rep


# --------------------------------------------------------------------------- answers that say no (first review)

def test_a_name_ncbi_rejects_is_an_answer(taxa, monkeypatch):
    monkeypatch.setattr(taxa.subprocess, "run", lambda cmd, *a, **k: types.SimpleNamespace(
        returncode=1, stdout="", stderr="Error: The taxonomy name 'cow' is not recognized."))
    taxa.common_name("cow")
    assert taxa.SOURCES["Cow"] == ("none", "")


def test_a_rejection_is_no_network_failure_in_phylo_sort(monkeypatch):
    import phylo_sort as ps
    monkeypatch.setattr(ps, "PROVENANCE", {})
    ps._network_failure("datasets summary taxonomy cow", "exit 1: Error: The taxonomy name 'cow' is not recognized.")
    assert "network_failures" not in ps.PROVENANCE
    ps._network_failure("datasets summary taxonomy cow", "Command timed out after 30 seconds")
    assert ps.PROVENANCE["network_failures"] == ["datasets summary taxonomy cow: Command timed out after 30 seconds"]


def test_a_uniprot_404_is_not_found(monkeypatch):
    import annotate_goi_exons as ag
    monkeypatch.setattr(ag, "UNIPROT_LOOKUP", {})

    def gone(*a, **k):
        raise ag.urllib.error.HTTPError("u", 404, "Not Found", {}, None)
    monkeypatch.setattr(ag.urllib.request, "urlopen", gone)
    ag.fetch_uniprot_names("P99999")
    assert ag.UNIPROT_LOOKUP["status"] == "not_found"
