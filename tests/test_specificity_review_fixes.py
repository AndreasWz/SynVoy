"""Review fixes RF-2 in bin/synvoy_specificity.py (all synthetic):

* one-sided neighbourhood (S2-F3): with flanking genes on one side only, the hull of the
  flanking genes ended before q and its adjacent relatives; it now reaches the home locus;
* q's own home locus in every locus run: a run from a secondary home locus knew only its own
  neighbourhood, so a p* next to q's own gene allowed the positional tie-break;
* no silent fallback (S2-F5, SEG-4): a failed blastp call of the family test is logged and
  counted and gives E = inf (it fell through to plain Karlin-Altschul in silence);
* parsed like synvoy_evalue (SEG-5): ``bl2seq`` or ``search`` blocks;
* direction (SEG-3): q is the query, also in the walk's relative test (blastp_cbs).
"""
import logging
import os
import random
import stat
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import extract_flanking_genes as efg  # noqa: E402
import synvoy_specificity as spec  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(rng, n=150):
    return "".join(rng.choice(AA) for _ in range(n))


def _mutate(rng, seq, frac):
    return "".join(rng.choice(AA) if rng.random() < frac else c for c in seq)


# ------------------------------------------------------------------ neighbourhoods
def _proteome(rng):
    q = _protein(rng)
    paralog = _mutate(rng, q, 0.35)
    records = [("gene-Q", q), ("gene-P", paralog)] + [(f"gene-{i}", _protein(rng)) for i in range(8)]
    coords = {"gene-Q": ("chrA", 10_000, 11_000), "gene-P": ("chrA", 12_000, 13_000),
              "gene-0": ("chrA", 1_000, 2_000), "gene-1": ("chrA", 30_000, 31_000),
              "gene-2": ("chrA", 40_000, 41_000), "gene-3": ("chrA", 50_000, 51_000),
              "gene-4": ("chrB", 1_000, 2_000), "gene-5": ("chrB", 5_000, 6_000),
              "gene-6": ("chrB", 9_000, 10_000), "gene-7": ("chrB", 13_000, 14_000)}
    return spec.HomeProteome(records, coords), q, paralog


def test_the_hull_reaches_the_home_locus():
    proteome, _q, _p = _proteome(random.Random(3))
    # Flanking genes on the right of the locus only (q at a scaffold end, say).
    one_sided = proteome.neighbourhood(["gene-1", "gene-2", "gene-3"])
    assert "gene-P" not in one_sided and "gene-Q" not in one_sided
    hood = proteome.neighbourhood(["gene-1", "gene-2", "gene-3"], locus_span=(10_000, 11_000))
    assert {"gene-Q", "gene-P", "gene-1", "gene-2", "gene-3"} <= hood and "gene-0" not in hood


def test_q_s_own_locus_counts_in_a_secondary_locus_run():
    rng = random.Random(4)
    proteome, q, paralog = _proteome(rng)
    assert proteome.own_gene("GOI_q", q) == "gene-Q"
    assert proteome.own_neighbourhood("GOI_q", q, 1) == {"gene-0", "gene-Q", "gene-P"}
    x = _mutate(rng, paralog, 0.05)                 # a model of the paralog, found from locus 2
    secondary = proteome.neighbourhood(["gene-4", "gene-5", "gene-7"], locus_span=(9_000, 10_000))
    assert "gene-P" not in secondary
    before = spec.specificity(proteome, "GOI_q", q, x, 50.0, secondary)
    after = spec.specificity(proteome, "GOI_q", q, x, 50.0, secondary, n_flanking=1)
    assert before["best_other"] == after["best_other"] == "gene-P"
    assert before["in_neighbourhood"] is False and after["in_neighbourhood"] is True
    # Without coordinates for this locus, q's own locus still gives an answer.
    assert spec.specificity(proteome, "GOI_q", q, x, 50.0, None, n_flanking=1)["in_neighbourhood"] is True


def test_no_own_gene_no_own_locus():
    rng = random.Random(5)
    proteome, _q, _p = _proteome(rng)
    stranger = _protein(rng)
    assert proteome.own_gene("GOI_s", stranger) is None
    assert proteome.own_neighbourhood("GOI_s", stranger, 2) is None


# ------------------------------------------------------------------ the family test's blastp
def _fake_blastp(tmp_path, body):
    script = tmp_path / "blastp"
    script.write_text("#!/bin/sh\n" + body + "\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(script)


@pytest.mark.parametrize("body,reason", [
    ("echo boom >&2; exit 3", "exit 3: boom"),
    ("exit 0", "no output"),
    ("echo '{not json'", "unparsable JSON output"),
    ("echo '{\"BlastOutput2\": [{\"report\": {\"results\": {}}}]}'", "unparsable JSON output"),
    ("echo '{\"other\": 1}'", "unparsable JSON output"),
])
def test_a_failed_call_is_logged_counted_and_inf(tmp_path, monkeypatch, caplog, body, reason):
    monkeypatch.setattr(spec, "_blastp", lambda: _fake_blastp(tmp_path, body))
    monkeypatch.setattr(spec, "FAMILY_TEST_FAILURES", {"count": 0})
    with caplog.at_level(logging.WARNING):
        assert spec.pairwise_evalue("MKV" * 20, "MKV" * 20) == float("inf")
    assert spec.FAMILY_TEST_FAILURES["count"] == 1
    assert any(reason in r.getMessage() and "no family" in r.getMessage() for r in caplog.records)


@pytest.mark.parametrize("payload,expected", [
    ('{"BlastOutput2": [{"report": {"results": {"bl2seq": [{"hits": [{"hsps": [{"evalue": 2e-5}, '
     '{"evalue": 3e-9}]}]}]}}}]}', 3e-9),
    ('{"BlastOutput2": [{"report": {"results": {"search": {"hits": [{"hsps": [{"evalue": 1e-7}]}]}}}}]}', 1e-7),
    ('{"BlastOutput2": [{"report": {"results": {"bl2seq": [{"hits": []}]}}}]}', float("inf")),
])
def test_bl2seq_and_search_blocks_are_read(tmp_path, monkeypatch, payload, expected):
    monkeypatch.setattr(spec, "_blastp", lambda: _fake_blastp(tmp_path, f"echo '{payload}'"))
    monkeypatch.setattr(spec, "FAMILY_TEST_FAILURES", {"count": 0})
    assert spec.pairwise_evalue("MKV" * 20, "MKV" * 20) == expected
    assert spec.FAMILY_TEST_FAILURES["count"] == 0          # no hit is an answer, not a failure


def test_the_goi_is_the_query_in_the_walk_s_relative_test(monkeypatch):
    seen = []
    monkeypatch.setattr(spec, "pairwise_evalue", lambda a, b: seen.append((a, b)) or 1.0)
    assert not efg.is_goi_relative("NEIGHBOUR", ["GOISEQ"], 1e-5, "blastp_cbs")
    assert seen == [("GOISEQ", "NEIGHBOUR")]


def test_the_modules_pass_n_flanking_genes():
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "modules")
    for module in ("iterative_search.nf", "rescue_goi_hull.nf"):
        assert "--n_flanking_genes ${params.n_flanking_genes}" in open(os.path.join(root, module)).read()
