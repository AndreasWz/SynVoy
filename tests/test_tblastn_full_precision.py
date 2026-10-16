"""Review fix RF-3, F4: the window tblastn runs once into an archive; blast_formatter writes the
tabular output the search reads (byte-identical to a direct run) and a JSON output whose
E-values reach the label's tblastn E source at full precision (tabular output prints E with
3 decimals from 0.0009 to 0.1, 2 below 1, 1 below 10)."""
import os
import random
import shutil
import subprocess
import sys

import pytest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402

CODONS = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
          "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
          "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


@pytest.fixture
def blast_path(monkeypatch):
    monkeypatch.setenv("PATH", os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""))
    if not all(shutil.which(b) for b in ("tblastn", "makeblastdb", "blast_formatter")):
        pytest.skip("BLAST+ (tblastn, makeblastdb, blast_formatter) not on PATH")


def _window(tmp_path):
    """A 30-kb random window holding the query once on each strand (one copy mutated), a
    query and a decoy-like reversed query; the BLAST database of the window."""
    rng = random.Random(11)
    aa = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(150))
    mutated = "".join(c if rng.random() > 0.55 else rng.choice("ACDEFGHIKLMNPQRSTVWY") for c in aa)
    dna = [rng.choice("ACGT") for _ in range(30000)]
    plus = "".join(CODONS[c] for c in aa)
    minus = "".join(CODONS[c] for c in mutated)
    minus = minus[::-1].translate(str.maketrans("ACGT", "TGCA"))
    dna[4000:4000 + len(plus)] = plus
    dna[20000:20000 + len(minus)] = minus
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "".join(dna) + "\n")
    query = tmp_path / "q.faa"
    query.write_text(f">GOI_x\n{aa}\n>GOI_copy_2|Genome_fna_b0_l1~rev\n{aa[::-1]}\n")
    db = str(tmp_path / "db")
    subprocess.run(["makeblastdb", "-in", str(region), "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    cmd = ["tblastn", "-query", str(query), "-db", db, "-evalue", "1000", "-seg", "no",
           "-max_target_seqs", "5000"]
    return cmd, {"GOI_x": len(aa), "GOI_copy_2|Genome_fna_b0_l1~rev": len(aa)}


def test_tabular_identical_and_every_row_has_its_full_e(blast_path, tmp_path):
    cmd, qlens = _window(tmp_path)
    direct = tmp_path / "direct.m8"
    subprocess.run(cmd + ["-outfmt", isr.TBLASTN_TABULAR, "-out", str(direct)], check=True)
    tab = tmp_path / "hits.m8"
    json_path = isr.run_tblastn_tabular_and_json(cmd, str(tab))
    assert json_path and os.path.exists(json_path)
    assert tab.read_bytes() == direct.read_bytes()
    full = synvoy_evalue.tblastn_json_evalues(json_path)
    hits = isr.parse_hits(str(tab), 0.0, 1, float("inf"), query_lengths=qlens)
    assert len(hits) > 5 and {h["strand"] for h in hits} == {"+", "-"}
    assert {h["query"] for h in hits} == set(qlens)
    finer = 0
    for h in hits:
        key = (h["query"], h["chrom"], h["qstart"], h["qend"], h["start"], h["end"], h["strand"])
        assert key in full, key
        e = full[key]
        assert abs(e - h["evalue"]) <= max(0.051 * e, 0.0005 + 1e-12) or h["evalue"] == 0.0, (e, h["evalue"])
        finer += e != h["evalue"]
    assert finer > 0   # the tabular strings are rounded; the JSON values are not


def test_without_blast_formatter_one_direct_run(blast_path, tmp_path, monkeypatch):
    cmd, _ = _window(tmp_path)
    direct = tmp_path / "direct.m8"
    subprocess.run(cmd + ["-outfmt", isr.TBLASTN_TABULAR, "-out", str(direct)], check=True)
    real_which = shutil.which
    monkeypatch.setattr(isr.shutil, "which", lambda b: None if b == "blast_formatter" else real_which(b))
    tab = tmp_path / "hits.m8"
    assert isr.run_tblastn_tabular_and_json(cmd, str(tab)) is None
    assert tab.read_bytes() == direct.read_bytes()
    assert not os.path.exists(f"{tab}.json")


def test_json_reader_keys_and_bad_input(tmp_path):
    js = tmp_path / "t.json"
    js.write_text('{"BlastOutput2": [{"report": {"results": {"search": {"query_title": "GOI_q desc",'
                  ' "hits": [{"description": [{"id": "gnl|BL_ORD_ID|0", "title": "region_seq"}],'
                  ' "hsps": [{"evalue": 0.00123456, "query_from": 3, "query_to": 40,'
                  ' "hit_from": 900, "hit_to": 790, "hit_frame": -2},'
                  ' {"evalue": 10.4946, "query_from": 1, "query_to": 20,'
                  ' "hit_from": 100, "hit_to": 159, "hit_frame": 1}]}]}}}}]}')
    full = synvoy_evalue.tblastn_json_evalues(str(js))
    assert full == {("GOI_q", "region_seq", 3, 40, 789, 900, "-"): 0.00123456,
                    ("GOI_q", "region_seq", 1, 20, 99, 159, "+"): 10.4946}
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    assert synvoy_evalue.tblastn_json_evalues(str(bad)) == {}
    assert synvoy_evalue.tblastn_json_evalues(str(tmp_path / "missing.json")) == {}
