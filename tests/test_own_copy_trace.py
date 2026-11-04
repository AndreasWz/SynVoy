"""Fix-phase d2 (2 Nov): the trace of an own copy set aside.

Since item 10 (c9) a home protein in q's own home locus that passes the own-gene rule is left out of p* for a model x,
unless x favours it over q and q's own gene by 10 bits. Nothing in the result said so. When the best home protein
against x other than q's own gene is such a copy, specificity() now returns it (own_copy_set_aside) and by how many
bits x favours it (own_copy_bits, signed); gff_attrs writes OwnCopySetAside and OwnCopyBits; the report's records carry
them, and summary.own_copy_set_aside_goi counts the records (AMBIGUOUS or better) that have one. No label moves.
"""
import json
import os
import random
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
pytest.importorskip("parasail")
import synvoy_specificity as ss  # noqa: E402
from sequence_utils import sw_bits  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(n, seed):
    rng = random.Random(seed)
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _mutate(seq, k, seed, lo=1, hi=None):
    rng = random.Random(seed)
    out = list(seq)
    for i in rng.sample(range(lo, hi or len(seq)), k):
        out[i] = rng.choice([a for a in AA if a != out[i]])
    return "".join(out)


Q = _protein(300, 1)
COPY = _mutate(Q, 6, 2, 1, 150)               # 98 % of q, differences in the first half only
FAR = _protein(200, 77)


def _home(*extra):
    genes = [("gene-q", Q, "c1", 10_000), ("gene-a0", COPY, "c1", 12_000)] + list(extra)
    genes.append(("gene-far", FAR, "c5", 1_000))
    return ss.HomeProteome([(g, s) for g, s, _c, _p in genes], {g: (c, p, p + 999) for g, _s, c, p in genes})


def _spec(home, x, own="gene-q"):
    return ss.specificity(home, "GOI_q", Q, x, sw_bits(Q, x), None, 5, own_gene_id=own)


def _closer_to_copy(k):
    """A target model that takes ``k`` of the copy's six differences: it favours the copy by a few bits."""
    diffs = [i for i in range(len(Q)) if Q[i] != COPY[i]]
    out = list(_mutate(Q, 20, 9, 150, 300))
    for i in diffs[:k]:
        out[i] = COPY[i]
    return "".join(out)


def test_a_copy_set_aside_is_named_with_its_bits():
    home = _home()
    x = _closer_to_copy(4)
    spec = _spec(home, x)
    favour = sw_bits(COPY, x) - sw_bits(Q, x)
    assert 0 < favour < ss.SEED_MIN_MARGIN_BITS
    assert spec["best_other"] != "gene-a0"
    assert spec["own_copy_set_aside"] == "gene-a0"
    assert spec["own_copy_bits"] == pytest.approx(favour, abs=0.05)


def test_the_bits_are_signed():
    """A model closer to q than to the copy, with the copy still the best other protein: a negative value."""
    home = _home()
    x = _mutate(Q, 20, 9, 150, 300)
    spec = _spec(home, x)
    assert spec["own_copy_set_aside"] == "gene-a0"
    assert spec["own_copy_bits"] < 0


def test_a_copy_told_apart_is_p_star_and_no_trace():
    """A model that favours the copy by >= 10 bits is the copy's: the copy stays p*, nothing was set aside."""
    far_copy = _mutate(Q, 14, 3, 1, 150)
    assert ss.passes_own_gene_rule(Q, far_copy)
    home = ss.HomeProteome([("gene-q", Q), ("gene-a0", far_copy), ("gene-far", FAR)],
                           {"gene-q": ("c1", 10_000, 10_999), "gene-a0": ("c1", 12_000, 12_999),
                            "gene-far": ("c5", 1_000, 1_999)})
    x = far_copy
    assert sw_bits(far_copy, x) - sw_bits(Q, x) >= ss.SEED_MIN_MARGIN_BITS
    spec = _spec(home, x)
    assert spec["best_other"] == "gene-a0"
    assert spec["own_copy_set_aside"] == "" and spec["own_copy_bits"] is None


def test_a_paralog_that_beats_the_copy_leaves_no_trace():
    """The best other protein is a real paralog, not the copy: p* did not change, nothing to trace."""
    x = _closer_to_copy(4)
    paralog = x                                  # a home protein identical to the model
    home = _home(("gene-par", paralog, "c3", 1_000))
    spec = _spec(home, x)
    assert spec["best_other"] == "gene-par"
    assert spec["own_copy_set_aside"] == ""


def test_without_an_own_gene_nothing_is_set_aside():
    home = ss.HomeProteome([("gene-a0", COPY), ("gene-far", FAR)],
                           {"gene-a0": ("c1", 12_000, 12_999), "gene-far": ("c5", 1_000, 1_999)})
    spec = ss.specificity(home, "GOI_q", _protein(150, 5), _closer_to_copy(4), 100.0, None, 5)
    assert spec["own_copy_set_aside"] == ""


def test_the_label_terms_do_not_move():
    """The trace adds keys only: every other term is what the specificity gave before."""
    home = _home()
    x = _closer_to_copy(4)
    spec = _spec(home, x)
    p_id, raw = home.best_other(x, home.own_copies("GOI_q", Q, 5, "gene-q", None))
    assert spec["best_other"] == p_id
    assert spec["margin"] == pytest.approx(sw_bits(Q, x) - (ss.raw_to_bits(raw) if p_id else 0.0))


def test_gff_attrs_write_the_trace():
    base = {"best_other": "gene-far", "margin": 120.0, "real_competitor": False, "in_neighbourhood": None}
    assert ss.gff_attrs(dict(base, own_copy_set_aside="", own_copy_bits=None)) == \
        "BestOtherHome=gene-far;SpecMargin=120.0;RealCompetitor=0"
    assert ss.gff_attrs(base) == "BestOtherHome=gene-far;SpecMargin=120.0;RealCompetitor=0"
    assert ss.gff_attrs(dict(base, own_copy_set_aside="gene-a0", own_copy_bits=6.24)).endswith(
        ";OwnCopySetAside=gene-a0;OwnCopyBits=+6.2")
    assert ss.gff_attrs(dict(base, own_copy_set_aside="gene-a0", own_copy_bits=-3.0)).endswith(";OwnCopyBits=-3.0")


def _gff_row(start, end, mid, conf, extra=""):
    attrs = (f"ID={mid};Name=GOI_G;SynVoy_Parent=GOI_G;SynVoyRole=goi;Confidence={conf};GOIClass=confident_goi;"
             f"EvidenceType=exon_annotation;Identity=60.0;QueryCoverage=0.900;ModelStatus=complete;"
             f"Basis=significant_specific{extra}")
    return "\t".join(["c1", "exon_annotation", "mRNA", str(start), str(end), "60", "+", ".", attrs]) + "\n"


def _report(tmp_path, rows):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n" + "".join(rows))
    out = tmp_path / "report.json"
    proc = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "generate_report.py"), "--results_dir", str(stage),
                           "--output", str(out), "--allow-empty"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return json.loads(out.read_text())


def test_the_report_carries_and_counts_the_trace(tmp_path):
    rep = _report(tmp_path, [
        _gff_row(1000, 2000, "GOI_G|m1", "HIGH", ";OwnCopySetAside=gene-a0;OwnCopyBits=+6.2"),
        _gff_row(5000, 6000, "GOI_G|m2", "AMBIGUOUS", ";OwnCopySetAside=gene-a0;OwnCopyBits=-1.0"),
        _gff_row(9000, 10000, "GOI_G|m3", "MEDIUM"),
        _gff_row(20000, 21000, "GOI_G|m4", "LOW", ";OwnCopySetAside=gene-a0;OwnCopyBits=+2.0")])
    recs = {r["start"]: r for r in rep["goi_dedup"]["records"]}
    assert (recs[1000]["own_copy_set_aside"], recs[1000]["own_copy_bits"]) == ("gene-a0", "+6.2")
    assert recs[5000]["own_copy_bits"] == "-1.0"
    assert "own_copy_set_aside" not in recs[9000] and "own_copy_bits" not in recs[9000]
    assert 20000 not in recs                                   # LOW: not a record, not counted
    assert rep["summary"]["own_copy_set_aside_goi"] == 2


def test_a_report_without_the_trace_counts_zero(tmp_path):
    rep = _report(tmp_path, [_gff_row(1000, 2000, "GOI_G|m1", "HIGH")])
    assert rep["summary"]["own_copy_set_aside_goi"] == 0
    assert all("own_copy_set_aside" not in r for r in rep["goi_dedup"]["records"])


def test_a_copy_tied_with_a_paralog_leaves_no_trace():
    """Review of the session: a set-aside copy with the same score as a paralog changed nothing (p* is the paralog at
    the same bits), so no trace is written, whatever the ids' order."""
    x = _closer_to_copy(4)
    for par_id in ("gene-b", "gene-0"):
        home = _home((par_id, COPY, "c3", 1_000))
        spec = _spec(home, x)
        assert spec["best_other"] == par_id and spec["own_copy_set_aside"] == ""
