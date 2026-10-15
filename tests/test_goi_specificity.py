"""Specificity of GOI models against the home proteome, and what step 2 of the label
rework (2026-10-14) does with it.

bin/synvoy_specificity.py: p* (the best other home protein, exhaustive Smith-Waterman),
the margin bits(q, x) - bits(p*, x), whether p* is a real competitor (E <= 0.01 over the
proteome, or a home family member of q at a composition-corrected pairwise E <= 1e-3), and
whether p* lies in the home neighbourhood of the locus. iterative_search_runner writes them
on every kept GOI model; the window search deduplicates its hits per query (A5).
"""
import os
import random
import subprocess
import sys
import types

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
pytest.importorskip("parasail")
import iterative_search_runner as isr  # noqa: E402
import synvoy_specificity as spec  # noqa: E402
from sequence_utils import sw_bits  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _random(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice(AA) for _ in range(n))


def _mutate(seq, frac, seed):
    rng = random.Random(seed)
    out = list(seq)
    for i in rng.sample(range(len(seq)), int(len(seq) * frac)):
        out[i] = rng.choice([a for a in AA if a != out[i]])
    return "".join(out)


Q = _random(220, 1)
PARALOG = _mutate(Q, 0.45, 2)          # a home relative of q
ORTHOLOG = _mutate(Q, 0.25, 3)         # what an ortholog of q looks like
PARALOG_COPY = _mutate(PARALOG, 0.1, 4)  # the relative's ortholog in a target genome


def _write_proteome(folder, with_q=True, with_coords=True, paralog_chrom="chr1"):
    recs = [("gene-par", PARALOG)] + [(f"gene-u{i}", _random(300, 100 + i)) for i in range(20)]
    if with_q:
        recs.append(("gene-q", Q))
    recs += [("gene-f1", _random(250, 50)), ("gene-f2", _random(250, 51))]
    with open(os.path.join(folder, spec.PROTEOME_FILE), "w") as fh:
        for rid, seq in recs:
            fh.write(f">{rid}\n{seq}\n")
    if with_coords:
        rows = [("gene-f1", "chr1", 1000, 2000), ("gene-q", "chr1", 5000, 6000),
                ("gene-par", paralog_chrom, 8000, 9000), ("gene-f2", "chr1", 12000, 13000)]
        rows += [(f"gene-u{i}", "chr2", 100000 * (i + 1), 100000 * (i + 1) + 900) for i in range(20)]
        with open(os.path.join(folder, spec.GENES_FILE), "w") as fh:
            fh.write("gene_id\tchrom\tstart\tend\tstrand\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\t+\n")
    return spec.HomeProteome.load(folder)


# ---------------------------------------------------------------- the helper

def test_q_is_its_own_gene_and_never_p_star(tmp_path):
    hp = _write_proteome(str(tmp_path))
    assert hp.own_gene("GOI_q", Q) == "gene-q"
    s = spec.specificity(hp, "GOI_q", Q, ORTHOLOG, sw_bits(Q, ORTHOLOG), hp.neighbourhood(["gene-f1", "gene-f2"]))
    assert s["best_other"] == "gene-par"
    assert s["margin"] == pytest.approx(sw_bits(Q, ORTHOLOG) - sw_bits(PARALOG, ORTHOLOG), abs=0.01)
    assert s["margin"] > spec.SEED_MIN_MARGIN_BITS
    assert s["real_competitor"] is True and s["seed_specific"] is True


def test_a_relatives_copy_is_not_specific(tmp_path):
    hp = _write_proteome(str(tmp_path))
    s = spec.specificity(hp, "GOI_q", Q, PARALOG_COPY, sw_bits(Q, PARALOG_COPY), None)
    assert s["best_other"] == "gene-par" and s["margin"] < 0
    assert s["real_competitor"] is True and s["seed_specific"] is False


def test_a_chance_match_has_no_real_competitor(tmp_path):
    hp = _write_proteome(str(tmp_path))
    junk = _random(60, 999)
    s = spec.specificity(hp, "GOI_q", Q, junk, sw_bits(Q, junk), None)
    assert s["real_competitor"] is False and s["seed_specific"] is True


def test_without_q_in_the_proteome_nothing_is_excluded(tmp_path):
    hp = _write_proteome(str(tmp_path), with_q=False)
    assert hp.own_gene("GOI_q", Q) is None
    s = spec.specificity(hp, "GOI_q", Q, ORTHOLOG, sw_bits(Q, ORTHOLOG), None)
    assert s["best_other"] == "gene-par"


def test_neighbourhood_is_the_span_of_the_flanking_genes(tmp_path):
    hp = _write_proteome(str(tmp_path))
    hood = hp.neighbourhood(["gene-f1", "gene-f2"])
    assert hood == {"gene-f1", "gene-q", "gene-par", "gene-f2"}
    s = spec.specificity(hp, "GOI_q", Q, ORTHOLOG, sw_bits(Q, ORTHOLOG), hood)
    assert s["in_neighbourhood"] is True
    other = _write_proteome(str(tmp_path), paralog_chrom="chr3")
    hood = other.neighbourhood(["gene-f1", "gene-f2"])
    assert "gene-par" not in hood
    assert spec.specificity(other, "GOI_q", Q, ORTHOLOG, sw_bits(Q, ORTHOLOG), hood)["in_neighbourhood"] is False


def test_no_coordinates_no_neighbourhood(tmp_path):
    hp = _write_proteome(str(tmp_path), with_coords=False)
    assert hp.neighbourhood(["gene-f1"]) is None
    s = spec.specificity(hp, "GOI_q", Q, ORTHOLOG, sw_bits(Q, ORTHOLOG), None)
    assert s["in_neighbourhood"] is None
    assert "BestOtherInNeighbourhood" not in spec.gff_attrs(s)


def test_family_test_uses_blastp_composition_statistics():
    """q and its relative are one family; the E-value is blastp's own, composition-corrected,
    read at full precision from its JSON output."""
    blastp = spec._blastp()
    if not blastp:
        pytest.skip("blastp not installed")
    assert spec.pairwise_evalue(Q, PARALOG) < spec.FAMILY_MEMBER_EVALUE
    a, b = Q[:120], _mutate(Q[:120], 0.55, 11)
    import json, tempfile
    with tempfile.TemporaryDirectory() as tmp:
        fa, fb = os.path.join(tmp, "a.faa"), os.path.join(tmp, "b.faa")
        open(fa, "w").write(f">a\n{a}\n")
        open(fb, "w").write(f">b\n{b}\n")
        out = subprocess.run([blastp, "-query", fa, "-subject", fb, "-comp_based_stats", "2",
                              "-seg", spec.FAMILY_TEST_SEG, "-outfmt", "15"],
                             capture_output=True, text=True, check=True).stdout
    hits = json.loads(out)["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["hits"]
    expect = min(h["evalue"] for hit in hits for h in hit["hsps"]) if hits else float("inf")
    assert spec.pairwise_evalue(a, b) == expect


def test_family_test_without_blastp_falls_back_to_karlin_altschul(monkeypatch):
    import math
    import parasail
    monkeypatch.setattr(spec, "_blastp", lambda: None)
    score = parasail.sw_striped_sat(Q, PARALOG, 12, 1, parasail.blosum62).score
    expect = spec.SW_K_11_1 * len(Q) * len(PARALOG) * math.exp(-spec.SW_LAMBDA_11_1 * score)
    assert spec.pairwise_evalue(Q, PARALOG) == pytest.approx(expect)


def test_gff_attrs():
    s = {"best_other": "gene-par", "margin": 12.345, "real_competitor": True, "in_neighbourhood": False}
    assert spec.gff_attrs(s) == ("BestOtherHome=gene-par;SpecMargin=12.3;RealCompetitor=1;"
                                 "BestOtherInNeighbourhood=0")


def test_missing_proteome_is_none(tmp_path):
    assert spec.HomeProteome.load(str(tmp_path)) is None
    assert spec.HomeProteome.load(None) is None


# ---------------------------------------------------------------- the search

def _cand(mid, seq, parent="GOI_q"):
    attrs = f"ID={parent}|g_b0_{mid};Name={parent};SynVoy_Parent={parent};SynVoyRole=goi;Confidence=MEDIUM"
    return {"start": 100, "end": 900, "record": {"id": f"{parent}|g_b0_{mid}", "seq": seq},
            "gff": [f"chr1\texon_annotation\tmRNA\t100\t900\t50\t+\t.\t{attrs}"]}


def test_runner_writes_the_terms_on_kept_models(tmp_path):
    _write_proteome(str(tmp_path))
    isr._HOME_PROTEOME_CACHE.clear()
    args = types.SimpleNamespace(home_db_dir=str(tmp_path), home_flanking_ids=["gene-f1", "gene-f2"])
    by_parent = {"GOI_q": {"id": "GOI_q", "seq": Q}}
    cands = [_cand("l1", ORTHOLOG), _cand("l2", PARALOG_COPY)]
    isr.score_goi_candidates(cands, by_parent)
    assert isr.annotate_goi_specificity(cands, by_parent, args) == 2
    a1 = isr._parse_gff_attr_string(cands[0]["gff"][0].split("\t")[8])
    a2 = isr._parse_gff_attr_string(cands[1]["gff"][0].split("\t")[8])
    assert a1["BestOtherHome"] == "gene-par" and float(a1["SpecMargin"]) > 10 and a1["RealCompetitor"] == "1"
    assert float(a2["SpecMargin"]) < 0 and a2["BestOtherInNeighbourhood"] == "1"
    isr._HOME_PROTEOME_CACHE.clear()


def test_runner_without_proteome_writes_nothing(tmp_path):
    isr._HOME_PROTEOME_CACHE.clear()
    args = types.SimpleNamespace(home_db_dir=str(tmp_path), home_flanking_ids=[])
    cands = [_cand("l1", ORTHOLOG)]
    isr.score_goi_candidates(cands, {"GOI_q": {"seq": Q}})
    assert isr.annotate_goi_specificity(cands, {"GOI_q": {"seq": Q}}, args) == 0
    assert "SpecMargin" not in cands[0]["gff"][0]
    isr._HOME_PROTEOME_CACHE.clear()


def test_window_hits_are_deduplicated_per_query(tmp_path, monkeypatch):
    """A5: the relative's hit at a paralog no longer removes the home protein's hit there."""
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")

    def fake_run(cmd, *a, **k):
        if cmd[0] == "tblastn":
            with open(cmd[cmd.index("-out") + 1], "w") as fh:
                fh.write("GOI_q\tregion_seq\t40.0\t50\t30\t0\t1\t50\t100\t249\t1e-05\t60.0\n")
                fh.write("GOI_q~rel\tregion_seq\t90.0\t50\t5\t0\t1\t50\t100\t249\t1e-20\t110.0\n")
                fh.write("GOI_q\tregion_seq\t38.0\t50\t31\t0\t1\t50\t110\t259\t1e-04\t50.0\n")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(isr.subprocess, "run", fake_run)
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", lambda **k: (False, "skipped"))
    args = types.SimpleNamespace(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_evalue_mult=1000.0,
                                 aug_relaxed_evalue_cap=10.0, aug_relaxed_parse_evalue_mult=10.0,
                                 aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6,
                                 min_identity=10.0, aug_relaxed_length_min=15, aug_relaxed_length_div=2,
                                 min_length=10, enable_smith_waterman=False, deterministic_goi_search=True)
    hits = isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": Q}, {"id": "GOI_q~rel", "seq": PARALOG}],
                                    "g1", args, "t1", 1)
    kept = sorted((h["query"], h["start"]) for h in hits)
    # one hit per query at that place: the home protein's best and the relative's
    assert kept == [("GOI_q", 99), ("GOI_q~rel", 99)]


# ---------------------------------------------------------------- the home proteome's coordinates

def test_gff_to_faa_writes_gene_coordinates(tmp_path):
    genome = tmp_path / "g.fna"
    genome.write_text(">chr1\n" + "ATG" + "GCT" * 40 + "TAA" + "N" * 100 + "\n")
    gff = tmp_path / "g.gff"
    gff.write_text("\n".join([
        "chr1\tx\tgene\t1\t126\t.\t+\t.\tID=gene-a",
        "chr1\tx\tmRNA\t1\t126\t.\t+\t.\tID=rna-a;Parent=gene-a",
        "chr1\tx\tCDS\t1\t126\t.\t+\t0\tID=cds-a;Parent=rna-a"]) + "\n")
    out = tmp_path / "p.faa"
    tsv = tmp_path / "genes.tsv"
    subprocess.run([sys.executable, os.path.join(REPO, "bin", "gff_to_faa.py"), "--gff", str(gff),
                    "--genome", str(genome), "--output", str(out), "--coords_tsv", str(tsv)], check=True)
    assert tsv.read_text().splitlines() == ["gene_id\tchrom\tstart\tend\tstrand", "gene-a\tchr1\t1\t126\t+"]


def test_prepare_home_puts_proteome_and_coordinates_into_the_db_folder():
    text = open(os.path.join(REPO, "modules", "prepare_home.nf")).read()
    assert "--coords_tsv home_genes.tsv" in text
    assert f"home_proteome_db/{spec.PROTEOME_FILE}" in text and f"home_proteome_db/{spec.GENES_FILE}" in text
