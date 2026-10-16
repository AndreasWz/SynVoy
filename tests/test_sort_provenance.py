"""Review fix RF-3 (provenance): phylo_sort.py --provenance records which strategy ordered the
target genomes (and so the search waves) and, for the taxonomy strategy, the database ETE3
opened, with its md5: the laptop and LRZ ordered the waves differently, and nothing said why."""
import hashlib
import json
import os
import sys

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)

import phylo_sort  # noqa: E402


def _genomes(tmp_path):
    tdir = tmp_path / "targets"
    tdir.mkdir()
    (tdir / "Bgenus_beta.fna").write_text(">b\n" + "ACGTTGCA" * 400 + "\n")
    (tdir / "Agenus_alpha.fna").write_text(">a\n" + "ACGTACGT" * 400 + "\n")
    home = tmp_path / "home.fna"
    home.write_text(">h\n" + "ACGTACGT" * 400 + "\n")
    return tdir, home


def _run(monkeypatch, tmp_path, tdir, home, taxdb="NO_DB"):
    out, prov = tmp_path / "sorted_genomes.txt", tmp_path / "sort_provenance.json"
    monkeypatch.setattr(phylo_sort, "PROVENANCE", {})
    monkeypatch.setattr(sys, "argv", ["phylo_sort.py", "--home", "Hgenus_home", "--home_fasta", str(home),
                                      "--targets_dir", str(tdir), "--taxdb", taxdb,
                                      "--output", str(out), "--provenance", str(prov)])
    phylo_sort.main()
    return out.read_text().splitlines(), json.loads(prov.read_text())


def test_minhash_order_is_recorded(monkeypatch, tmp_path):
    monkeypatch.setattr(phylo_sort, "NCBITaxa", None)
    monkeypatch.setattr(phylo_sort, "_datasets_cli_available", lambda: False)
    tdir, home = _genomes(tmp_path)
    rows, prov = _run(monkeypatch, tmp_path, tdir, home)
    assert prov["strategy"] == "minhash" and prov["minhash"]["k"] == 21
    assert prov["order"] == [r.split("\t") for r in rows] and rows[0].startswith("Agenus_alpha.fna")
    assert "taxonomy_db" not in prov and prov["taxdb"] == "NO_DB"


def test_taxonomy_database_path_and_md5_are_recorded(monkeypatch, tmp_path):
    db = tmp_path / "taxa.sqlite"
    db.write_bytes(b"not really sqlite, only fingerprinted")
    lineages = {1: [1, 2, 3, 10], 2: [1, 2, 3, 11], 3: [1, 2, 20, 21]}
    names = {"Hgenus home": 1, "Agenus alpha": 2, "Bgenus beta": 3}

    class FakeNCBITaxa:
        def __init__(self, dbfile=None):
            self.dbfile = dbfile or str(db)

        def get_name_translator(self, wanted):
            return {w: [names[w]] for w in wanted if w in names}

        def get_lineage(self, taxid):
            return lineages[taxid]

    monkeypatch.setattr(phylo_sort, "NCBITaxa", FakeNCBITaxa)
    tdir, home = _genomes(tmp_path)
    rows, prov = _run(monkeypatch, tmp_path, tdir, home, taxdb=str(db))
    assert prov["strategy"] == "ete3_taxonomy" and prov["home_taxid"] == 1
    assert prov["taxonomy_db"]["path"] == str(db)
    assert prov["taxonomy_db"]["md5"] == hashlib.md5(db.read_bytes()).hexdigest()
    assert prov["taxonomy_db"]["bytes"] == db.stat().st_size
    assert [r.split("\t")[0] for r in rows] == ["Agenus_alpha.fna", "Bgenus_beta.fna"]
    assert prov["targets_without_taxid"] == []


def test_a_database_that_cannot_place_the_home_is_not_reported_as_used(monkeypatch, tmp_path):
    db = tmp_path / "taxa.sqlite"
    db.write_bytes(b"x")

    class EmptyNCBITaxa:
        def __init__(self, dbfile=None):
            self.dbfile = str(db)

        def get_name_translator(self, wanted):
            return {}

    monkeypatch.setattr(phylo_sort, "NCBITaxa", EmptyNCBITaxa)
    monkeypatch.setattr(phylo_sort, "_datasets_cli_available", lambda: False)
    tdir, home = _genomes(tmp_path)
    _, prov = _run(monkeypatch, tmp_path, tdir, home, taxdb=str(db))
    assert prov["strategy"] == "minhash"
    assert "taxonomy_db" not in prov and prov["taxonomy_db_not_used"]["path"] == str(db)


def test_module_passes_the_provenance_file():
    nf = open(os.path.join(BIN, "..", "modules", "phylo_sort.nf")).read()
    assert nf.count("--provenance sort_provenance.json") == 3
    assert 'path "sort_provenance.json", emit: provenance, optional: true' in nf
