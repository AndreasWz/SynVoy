"""Region scores must not depend on the line order of the hit table.

The search tools write their tables in a thread-dependent order, so two runs on
identical inputs hand cluster_grs.py the same hits in a different order. Hits were
sorted by (chrom, start) only; hits sharing a start stayed in input order, the
order consistency (a longest-collinear-run over hits in target order) changed with
it, and so did the region score and the region NAME that embeds it. Measured on a
real run: five line orders of one file gave consistency 0.875-0.917 and five
different scores.
"""
from __future__ import annotations

import itertools
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import cluster_grs as cg  # noqa: E402

GENES = ["f0", "f1", "f2", "GOI_x", "f3", "f4", "f5", "f6"]


def _fixture(tmp: Path):
    bed = tmp / "syn.bed"
    bed.write_text("".join(f"home\t{i * 2000}\t{i * 2000 + 900}\t{g}\t0\t+\n"
                           for i, g in enumerate(GENES)))
    rows = []
    # A collinear neighbourhood on chr1 in which several hits share a start
    # coordinate -- exon-level queries of two genes landing on the same window --
    # plus a second, scrambled cluster on chr2.
    pos = 100_000
    for g in ["f0", "f1", "f2", "GOI_x", "f3", "f4", "f5", "f6"]:
        rows.append((g, "chr1", 80.0, pos, pos + 600, 1e-40))
        rows.append((f"{g}|exon_1", "chr1", 75.0, pos, pos + 300, 1e-20))      # same start
        pos += 20_000
    # ties between DIFFERENT genes at one start: order is undefined by position alone
    rows.append(("f5", "chr1", 60.0, 100_000, 100_450, 1e-9))
    rows.append(("f6|exon_2", "chr1", 55.0, 120_000, 120_200, 1e-6))
    rows.append(("f1", "chr1", 50.0, 220_000, 220_500, 1e-5))
    for g, p in zip(["f4", "f1", "f6", "f0", "f3"], itertools.count(500_000, 15_000)):
        rows.append((g, "chr2", 45.0, p, p + 500, 1e-8))
        rows.append((g, "chr2", 44.0, p, p + 480, 1e-7))
    return bed, rows


def _run(tmp: Path, bed: Path, rows, tag: str):
    m8 = tmp / f"{tag}.m8"
    m8.write_text("".join(f"{q}\t{t}\t{pid}\t100\t0\t0\t1\t100\t{s}\t{e}\t{ev}\t120\n"
                          for q, t, pid, s, e, ev in rows))
    out, scores = tmp / f"{tag}.bed", tmp / f"{tag}.scores.tsv"
    subprocess.run([sys.executable, str(ROOT / "bin" / "cluster_grs.py"),
                    "--hits", str(m8), "--synteny_bed", str(bed),
                    "--output", str(out), "--scores_output", str(scores),
                    "--flanking_count", "7", "--min_score", "0.1"],
                   check=True, capture_output=True, text=True)
    return out.read_text(), scores.read_text()


def test_scores_are_identical_for_any_line_order(tmp_path):
    bed, rows = _fixture(tmp_path)
    reference = _run(tmp_path, bed, rows, "ref")
    assert reference[1].count("\n") >= 2, "fixture produced no region to compare"
    orders = {"reversed": rows[::-1], "sorted": sorted(rows)}
    for seed in range(8):
        shuffled = rows[:]
        random.Random(seed).shuffle(shuffled)
        orders[f"shuffle{seed}"] = shuffled
    for tag, order in orders.items():
        got = _run(tmp_path, bed, order, tag)
        assert got == reference, f"output changed for line order '{tag}'"


def test_hit_sort_key_is_a_total_order():
    """Two hits that differ in any field must never compare equal."""
    a = {"chrom": "c", "start": 5, "end": 9, "query": "f1", "strand": "+",
         "identity": 80.0, "evalue": 1e-9}
    variants = [dict(a, end=10), dict(a, query="f2"), dict(a, strand="-"),
                dict(a, identity=70.0), dict(a, evalue=1e-3), dict(a, start=6),
                dict(a, chrom="d")]
    keys = {cg.hit_sort_key(h) for h in [a] + variants}
    assert len(keys) == len(variants) + 1


def test_clustering_sorts_ties_the_same_way_whatever_the_input_order():
    hits = [{"query": q, "chrom": "c", "start": 100, "end": 400, "strand": "+",
             "identity": 50.0, "evalue": 1e-5} for q in ("f3", "f1", "f2")]
    first = [h["query"] for h in cg.cluster_hits_proximity(list(hits), {}, 1000)[0]]
    for perm in itertools.permutations(hits):
        got = [h["query"] for h in cg.cluster_hits_proximity(list(perm), {}, 1000)[0]]
        assert got == first == ["f1", "f2", "f3"]
