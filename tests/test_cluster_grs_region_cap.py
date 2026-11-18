#!/usr/bin/env python3
"""The region step's selection and its regions (cluster_grs.py):

- the adaptive cap writes every cluster that passes a floor, in its order, up to the cap; it
  stopped at the first cluster below both floors, and a weak GOI cluster ordered first hid a
  4-gene neighbourhood;
- the score floor is relative to the best score of all clusters, not the first one's;
- a chain of hits longer than the largest neighbourhood the search considers
  (2 x home span x expansion) is split at its largest gaps;
- the chain's reach is measured from its furthest end so far;
- raw GOI hits are GOI intervals only where the target GFF holds no GOI model of any label.
"""
import contextlib
import csv
import io
import json
import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

import cluster_grs as cg  # noqa: E402
import synvoy_scale  # noqa: E402


def _hit_line(q, chrom, s, e, strand="+", ident=80.0, ev=1e-30):
    ts, te = (s, e) if strand == "+" else (e, s)
    return f"{q}\t{chrom}\t{ident}\t300\t0\t0\t1\t100\t{ts}\t{te}\t{ev}\t200\n"


def _home_bed(tmp):
    bed = tmp / "home.bed"
    with open(bed, "w") as fh:
        for i in range(1, 11):
            fh.write(f"chrH\t{i * 10000}\t{i * 10000 + 2000}\tF{i}\t0\t{'+' if i % 2 else '-'}\n")
    return bed


def _three_clusters(with_weak_goi=True):
    rows = []
    # chrB: F1..F8 in home order and the GOI hit.
    for i in range(1, 9):
        rows.append(_hit_line(f"F{i}", "chrB", 100000 + i * 10000, 100000 + i * 10000 + 600,
                              "+" if i % 2 else "-"))
    rows.append(_hit_line("GOI_X", "chrB", 145000, 145600))
    if with_weak_goi:
        # chrA: a weak GOI cluster, one GOI hit and one flanking hit.
        rows.append(_hit_line("GOI_X", "chrA", 500000, 500600, ident=45.0, ev=1e-5))
        rows.append(_hit_line("F9", "chrA", 510000, 510600))
    # chrC: a 4-gene neighbourhood without a GOI hit.
    for j, i in enumerate((2, 3, 4, 5)):
        rows.append(_hit_line(f"F{i}", "chrC", 1000000 + j * 10000, 1000000 + j * 10000 + 600,
                              "+" if i % 2 else "-"))
    return rows


def _run(tmp, rows, extra=(), gff=None):
    m8 = tmp / "g.m8"
    m8.write_text("".join(rows))
    regions, scores = tmp / "g.regions.bed", tmp / "g.scores.tsv"
    argv = ["cluster_grs.py", "--hits", str(m8), "--synteny_bed", str(_home_bed(tmp)),
            "--output", str(regions), "--scores_output", str(scores),
            "--cluster_distance", "50000", "--goi_overlap_bonus", "0.15", *extra]
    if gff:
        argv += ["--target_gff", str(gff)]
    old = sys.argv
    sys.argv = argv
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            cg.main()
    finally:
        sys.argv = old
    with open(scores, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t")), err.getvalue()


def test_a_weak_goi_cluster_no_longer_hides_a_neighbourhood(tmp_path):
    """The weak chrA GOI cluster passes neither floor: it is not written, and the 4-gene chrC
    neighbourhood behind it in the order is (the loop stopped at chrA and wrote chrB alone)."""
    rows, log = _run(tmp_path, _three_clusters())
    assert [r["chrom"] for r in rows] == ["chrB", "chrC"], log
    chr_c = next(r for r in rows if r["chrom"] == "chrC")
    assert chr_c["selection_reason"] in ("unique_gene_floor", "score_floor")
    assert int(chr_c["unique_genes"]) == 4


def test_without_the_weak_cluster_the_same_neighbourhood_is_written(tmp_path):
    rows, _ = _run(tmp_path, _three_clusters(with_weak_goi=False))
    assert [r["chrom"] for r in rows] == ["chrB", "chrC"]


def test_the_floor_comes_from_the_best_score_of_all_clusters(tmp_path):
    """A non-GOI cluster scoring above the GOI cluster ordered first sets the floor: a cluster
    under 0.3 x that score with fewer than 3 genes is not written."""
    rows = [_hit_line("GOI_X", "chrA", 500000, 500600, ident=45.0, ev=1e-5),
            _hit_line("F9", "chrA", 510000, 510600)]
    for i in range(1, 11):                          # chrB: all ten in order, no GOI
        rows.append(_hit_line(f"F{i}", "chrB", 100000 + i * 10000, 100000 + i * 10000 + 600,
                              "+" if i % 2 else "-"))
    rows.append(_hit_line("F1", "chrD", 900000, 900600))      # chrD: one gene
    out, log = _run(tmp_path, rows)
    best = max(float(r["score"]) for r in out)
    assert any(r["chrom"] == "chrB" for r in out)
    assert f"best_score={best:.2f}" in log
    assert all(r["chrom"] != "chrD" for r in out)


def test_the_adaptive_cap_still_caps(tmp_path):
    rows, _ = _run(tmp_path, _three_clusters(), extra=("--adaptive_max_regions", "1"))
    assert [r["chrom"] for r in rows] == ["chrB"]


def test_the_user_cap_is_unchanged(tmp_path):
    rows, _ = _run(tmp_path, _three_clusters(), extra=("--max_regions", "2"))
    assert [r["chrom"] for r in rows] == ["chrB", "chrA"]
    assert {r["selection_reason"] for r in rows} == {"user_cap"}


def test_the_first_cluster_is_written_when_none_passes(tmp_path):
    rows, _ = _run(tmp_path, [_hit_line("F1", "chrD", 900000, 900600)],
                   extra=("--adaptive_score_floor_abs", "0.99", "--adaptive_unique_gene_floor", "5"))
    assert len(rows) == 1
    assert rows[0]["selection_reason"] == "adaptive_backstop"


def _chain(n=300, step=40000):
    return [{"query": f"F{1 + i % 10}", "chrom": "chrZ", "start": i * step, "end": i * step + 500,
             "strand": "+", "identity": 30.0, "evalue": 1.0} for i in range(n)]


def test_a_long_chain_is_split_within_the_cap():
    clusters = cg.cluster_hits_proximity(_chain(), {}, 50000)
    assert len(clusters) == 1                       # one 12 Mb chain without a cap
    cap = cg.max_region_span({"home_span_bp": 100000, "expansion": 2.5})
    assert cap == 500000
    pieces = cg.split_long_clusters(clusters, cap)
    assert len(pieces) > 1
    assert sum(len(p) for p in pieces) == 300
    for p in pieces:
        assert max(h["end"] for h in p) - min(h["start"] for h in p) <= cap
    flat = [h for p in pieces for h in p]
    assert flat == clusters[0]                      # order kept, no hit lost or doubled


def test_the_split_cuts_the_largest_gap_first():
    hits = [{"query": "F1", "chrom": "c", "start": s, "end": s + 100, "strand": "+"}
            for s in (0, 1000, 2000, 50000, 51000, 60000)]
    pieces = cg.split_long_clusters([hits], 55000)
    assert [[h["start"] for h in p] for p in pieces] == [[0, 1000, 2000], [50000, 51000, 60000]]


def test_no_cap_without_the_search_distances():
    clusters = cg.cluster_hits_proximity(_chain(), {}, 50000)
    assert cg.max_region_span(None) is None
    assert cg.max_region_span({"expansion": 2.0}) is None
    assert cg.split_long_clusters(clusters, None) == clusters


def test_the_cap_is_read_from_the_target_gff(tmp_path):
    gff = tmp_path / "t.gff"
    record = {"home_span_bp": 100000, "expansion": 1.0, "cluster_distance": 50000}
    gff.write_text(f"{synvoy_scale.GFF_HEADER} {json.dumps(record)}\n")
    rows = []
    for i in range(30):                             # 30 hits 40 kb apart: 1.16 Mb chain
        rows.append(_hit_line(f"F{1 + i % 10}", "chrZ", i * 40000 + 1, i * 40000 + 500))
    _, log = _run(tmp_path, rows, gff=gff)
    assert "split" in log and "200000 bp" in log


def test_the_reach_is_the_furthest_end():
    long_hit = {"query": "F1", "chrom": "c", "start": 0, "end": 100000, "strand": "+"}
    inside = {"query": "F2", "chrom": "c", "start": 10000, "end": 10500, "strand": "+"}
    after = {"query": "F3", "chrom": "c", "start": 120000, "end": 120500, "strand": "+"}
    clusters = cg.cluster_hits_proximity([long_hit, inside, after], {}, 50000)
    assert len(clusters) == 1                       # 20 kb after the long hit's end


def _gff_row(chrom, start, end, conf, gid="GOI_X|g_m1"):
    return "\t".join([chrom, "SynVoy", "mRNA", str(start), str(end), ".", "+", ".",
                      f"ID={gid};SynVoy_Parent=GOI_X;SynVoyRole=goi;Confidence={conf};"
                      "GOIClass=ambiguous_goi_family_member;EvidenceType=fallback_hit_span"]) + "\n"


def test_raw_goi_hits_are_no_goi_interval_beside_a_low_model(tmp_path):
    """The search's only GOI model is LOW on chrA: the raw GOI hit there no longer makes the
    weak chrA cluster a GOI region ordered first."""
    gff = tmp_path / "t.gff"
    gff.write_text("##gff-version 3\n" + _gff_row("chrA", 500001, 500600, "LOW"))
    rows, log = _run(tmp_path, _three_clusters(), gff=gff)
    assert all(r["goi_overlap"] == "False" for r in rows)
    assert "GOI hits are not used as GOI intervals" in log
    chr_a = [r for r in rows if r["chrom"] == "chrA"]
    assert all(r["goi_model_confidence"] == "LOW" for r in chr_a)


def test_raw_goi_hits_stay_goi_intervals_without_any_goi_model(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("##gff-version 3\n")
    rows, log = _run(tmp_path, _three_clusters(), gff=gff)
    assert [(r["chrom"], r["goi_overlap"]) for r in rows] == [("chrB", "True"), ("chrC", "False")]
    assert "Derived 2 GOI interval(s) from hit file" in log


def test_goi_model_columns_in_the_scores_table(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("##gff-version 3\n" + _gff_row("chrB", 145001, 145600, "HIGH", "GOI_X|g_m2"))
    rows, _ = _run(tmp_path, _three_clusters(), gff=gff)
    chr_b = next(r for r in rows if r["chrom"] == "chrB")
    assert chr_b["goi_model_confidence"] == "HIGH"
    assert float(chr_b["goi_hit_best_identity"]) == 80.0
    chr_c = next(r for r in rows if r["chrom"] == "chrC")
    assert chr_c["goi_model_confidence"] == "" and chr_c["goi_hit_best_identity"] == ""


def test_best_goi_model_prefers_a_counted_model():
    models = [{"chrom": "c", "start": 100, "end": 200, "confidence": "MEDIUM",
               "goi_class": "genome_wide_strong_hit", "counted": False},
              {"chrom": "c", "start": 300, "end": 400, "confidence": "MEDIUM",
               "goi_class": "probable_goi", "counted": True},
              {"chrom": "d", "start": 100, "end": 200, "confidence": "HIGH",
               "goi_class": "confident_goi", "counted": True}]
    assert cg.best_goi_model(models, "c", 0, 1000)["goi_class"] == "probable_goi"
    assert cg.best_goi_model(models, "c", 0, 250)["goi_class"] == "genome_wide_strong_hit"
    assert cg.best_goi_model(models, "c", 600, 700) is None
    assert cg.best_goi_model(models, "c", 600, 700, padding_bp=250)["goi_class"] == "probable_goi"
