"""§1m fumble fix — GOI synteny-hull rescue (docs/TODO.md §1m).

Pure-logic tests for the parts that were actually buggy in development:
  - collect_features: flanking/GOI come as mRNA OR gene features (dedup by ID).
  - compute_hull: must SPAN a flanking gap (where the GOI sits) while excluding a far outlier.
  - _identity_pct: miniprot --gff Identity is a FRACTION (0.96), not a percent.
The end-to-end miniprot path is validated on real data (cow decorin chr5:21 Mb → 90% HIGH).
"""
import os
import sys

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import rescue_goi_hull as h  # noqa: E402


def _write_gff(tmp_path):
    """Region GFF: flanking in TWO clusters bracketing a gap (GOI sits in the gap) +
    a far rearranged outlier; no HIGH GOI present. Flanking emitted as mRNA (the real
    flanking_miniprot form); a MEDIUM goi present to confirm it does NOT block rescue."""
    p = tmp_path / "region.gff"
    rows = [
        "##gff-version 3",
        # cluster A: 1,000,000-1,500,000 (3 flanking, HIGH, as mRNA)
        "chr1\tflanking_annotation\tmRNA\t1000000\t1010000\t.\t+\t.\tID=gene-A1;SynVoyRole=flanking;Confidence=HIGH",
        "chr1\tflanking_annotation\tmRNA\t1200000\t1210000\t.\t+\t.\tID=gene-A2;SynVoyRole=flanking;Confidence=HIGH",
        "chr1\tflanking_annotation\tmRNA\t1490000\t1500000\t.\t+\t.\tID=gene-A3;SynVoyRole=flanking;Confidence=HIGH",
        # cluster B: 5,500,000-6,000,000 (3 flanking) — 4 Mb gap (> max_gap) from A; the
        # widening (not clustering) must bridge it, holding the GOI in the gap.
        "chr1\tflanking_annotation\tmRNA\t5500000\t5510000\t.\t+\t.\tID=gene-B1;SynVoyRole=flanking;Confidence=HIGH",
        "chr1\tflanking_annotation\tmRNA\t5800000\t5810000\t.\t+\t.\tID=gene-B2;SynVoyRole=flanking;Confidence=HIGH",
        "chr1\tflanking_annotation\tmRNA\t5990000\t6000000\t.\t+\t.\tID=gene-B3;SynVoyRole=flanking;Confidence=HIGH",
        # far outlier (90 Mb away) — must NOT blow up the hull
        "chr1\tflanking_annotation\tmRNA\t95000000\t95010000\t.\t+\t.\tID=gene-OUT;SynVoyRole=flanking;Confidence=HIGH",
        # a MEDIUM goi inside the hull (low-identity fragment) — should not block rescue
        "chr1\tx\tgene\t3500000\t3501000\t.\t+\t.\tID=GOI_x;SynVoyRole=goi;Confidence=MEDIUM",
    ]
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def test_collect_features_reads_mrna_and_dedups(tmp_path):
    flank, high_goi = h.collect_features(_write_gff(tmp_path))
    assert len(flank["chr1"]) == 7          # all 7 HIGH flanking (mRNA features)
    assert high_goi == {} or "chr1" not in high_goi  # the goi is MEDIUM, not HIGH


def test_collect_features_dedups_gene_plus_mrna(tmp_path):
    p = tmp_path / "dup.gff"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t100\t200\t.\t+\t.\tID=g1;SynVoyRole=flanking;Confidence=HIGH\n"
        "chr1\tx\tmRNA\t100\t200\t.\t+\t.\tID=g1.mRNA;SynVoyRole=flanking;Confidence=HIGH\n"
    )
    flank, _ = h.collect_features(str(p))
    assert len(flank["chr1"]) == 1          # gene + its mRNA counted once


def test_compute_hull_spans_gap_excludes_outlier(tmp_path):
    flank, _ = h.collect_features(_write_gff(tmp_path))
    hs, he, n = h.compute_hull(flank["chr1"], cluster_max_gap=3_000_000,
                               max_window=20_000_000)
    # Hull spans BOTH clusters (1.0-6.0 Mb), bracketing the 4 Mb gap via widening...
    assert hs <= 1_000_000 and he >= 6_000_000
    # ...and the GOI at 3.5 Mb (in the gap) is inside the hull.
    assert hs <= 3_500_000 <= he
    # ...but the 95 Mb outlier is excluded.
    assert he < 95_000_000
    assert n == 6                            # 6 neighbourhood flanking, not the outlier


def test_compute_hull_respects_max_window(tmp_path):
    # With a tiny max_window the widening can't reach the far cluster (B at 5.5 Mb), so the
    # hull stays on the dominant cluster (A, ~1.0-1.5 Mb).
    flank, _ = h.collect_features(_write_gff(tmp_path))
    hs, he, n = h.compute_hull(flank["chr1"], cluster_max_gap=3_000_000,
                               max_window=1_000_000)
    assert he - hs <= 1_500_000              # bounded around the dominant cluster
    assert n == 3


def test_identity_pct_fraction_to_percent():
    assert h._identity_pct("Identity=0.962;Rank=1") == 962.0 / 10  # 96.2
    assert h._identity_pct("Identity=0.40") == 40.0
    # Already a percent (defensive): values > 1 pass through.
    assert h._identity_pct("Identity=87.0") == 87.0


def test_cds_coverage():
    # 3 CDS of 100 bp each = 300 bp = 100 aa; query 200 aa -> 0.5 coverage.
    cds = [["c", "s", "CDS", "1", "100"], ["c", "s", "CDS", "201", "300"],
           ["c", "s", "CDS", "401", "500"]]
    assert abs(h._cds_coverage(cds, 200) - 0.5) < 0.01
    assert h._cds_coverage(cds, 0) == 0.0


def test_high_goi_in_hull_blocks_rescue(tmp_path):
    """If a HIGH GOI already sits in the hull, the rescue must be a no-op for that chrom."""
    p = tmp_path / "withgoi.gff"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tx\tmRNA\t1000000\t1010000\t.\t+\t.\tID=f1;SynVoyRole=flanking;Confidence=HIGH\n"
        "chr1\tx\tmRNA\t1200000\t1210000\t.\t+\t.\tID=f2;SynVoyRole=flanking;Confidence=HIGH\n"
        "chr1\tx\tmRNA\t1400000\t1410000\t.\t+\t.\tID=f3;SynVoyRole=flanking;Confidence=HIGH\n"
        "chr1\tx\tmRNA\t1600000\t1610000\t.\t+\t.\tID=f4;SynVoyRole=flanking;Confidence=HIGH\n"
        "chr1\tx\tgene\t1300000\t1305000\t.\t+\t.\tID=GOI_good;SynVoyRole=goi;Confidence=HIGH\n"
    )
    flank, high_goi = h.collect_features(str(p))
    hs, he, _ = h.compute_hull(flank["chr1"], 3_000_000, 20_000_000)
    blocked = any(h._interval_overlaps(hs, he, gs, ge) for gs, ge in high_goi["chr1"])
    assert blocked is True
