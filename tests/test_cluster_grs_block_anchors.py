#!/usr/bin/env python3
"""Plot region selection must reuse the SEARCH's blocks, not re-cluster raw hits.

Background (dcn_18 cow): cluster_grs scored each 50 kb proximity fragment on its own
flanking, so the true decorin block — whose flanking are spread across a §18-bridged
rearrangement gap (decorin at 21 Mb, flanking at 19 and 23 Mb) — scored ~0 flanking and
the chrX biglycan paralog (compact 3-flanking block) was surfaced instead (the "50%, not
90%" plot bug). cluster_grs now groups GOI + flanking by their search block id (`_b<N>_`),
so a GOI region inherits its block's real HIGH-flanking support and full neighbourhood
span — decorin (more flanking) outranks biglycan on synteny, not identity.
"""
import os
import sys

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

import cluster_grs as c  # noqa: E402


def _goi(chrom, start, end, block, ident, conf="HIGH"):
    attrs = (f"ID=GOI_DCN|tgt_fna_b{block}_l1_exon_ann;Name=GOI_DCN;SynVoy_Parent=GOI_DCN;"
             f"SynVoyRole=goi;Confidence={conf};Identity={ident};BlockFlankingSupport=99")
    return "\t".join([chrom, "exon_annotation", "mRNA", str(start), str(end), str(ident), "-", ".", attrs])


def _flank(chrom, start, end, block, gene, conf="HIGH"):
    attrs = (f"ID={gene}|tgt_fna_b{block}_fl1_flank_ann;Name={gene};SynVoy_Parent={gene};"
             f"SynVoyRole=flanking;Confidence={conf};Identity=95.0")
    return "\t".join([chrom, "flanking_annotation", "mRNA", str(start), str(end), "95.0", "+", ".", attrs])


def _write(tmp_path, rows):
    p = tmp_path / "target.gff"
    p.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return str(p)


def _decorin_vs_biglycan_gff(tmp_path):
    rows = [
        # block b0 = true decorin: GOI 90% at chr5:21 Mb, flanking spread 19-23 Mb (bridged gap)
        _goi("chr5", 21_000_000, 21_040_000, 0, 90.0),
        _flank("chr5", 19_200_000, 19_250_000, 0, "gene-ATP2B1"),
        _flank("chr5", 19_300_000, 19_320_000, 0, "gene-DUSP6"),
        _flank("chr5", 22_000_000, 22_010_000, 0, "gene-BTG1"),
        _flank("chr5", 23_200_000, 23_210_000, 0, "gene-NUDT4"),
        _flank("chr5", 23_400_000, 23_410_000, 0, "gene-SOCS2"),
        # block b11 = biglycan paralog: GOI 55% on chrX, only 2 flanking nearby
        _goi("chrX", 36_736_000, 36_741_000, 11, 55.0),
        _flank("chrX", 36_765_000, 36_806_000, 11, "gene-ATP2B1"),
        _flank("chrX", 36_600_000, 36_610_000, 11, "gene-DUSP6"),
        # a LOW GOI in b11 must NOT make b11 anchor on it
        _goi("chrX", 36_900_000, 36_901_000, 11, 25.0, conf="LOW"),
    ]
    return _write(tmp_path, rows)


def test_block_anchor_counts_flanking_and_spans_gap(tmp_path):
    anchors = c.load_goi_block_anchors_from_gff(_decorin_vs_biglycan_gff(tmp_path), padding_bp=0)
    by_chrom = {a["chrom"]: a for a in anchors}
    assert set(by_chrom) == {"chr5", "chrX"}
    dec = by_chrom["chr5"]
    big = by_chrom["chrX"]
    # decorin block: 5 HIGH flanking, span covers the WHOLE bridged neighbourhood (19-23.4 Mb)
    assert dec["flanking_support"] == 5
    assert dec["start"] <= 19_200_000 and dec["end"] >= 23_410_000
    assert dec["start"] <= 21_000_000 <= dec["end"]   # decorin (in the gap) is inside the span
    assert dec["identity"] == 90.0 and dec["confidence"] == "HIGH"
    # biglycan block: only 2 HIGH flanking
    assert big["flanking_support"] == 2
    assert big["identity"] == 55.0   # the LOW 25% GOI did not become the anchor identity


def test_decorin_block_outranks_biglycan_block(tmp_path):
    anchors = c.load_goi_block_anchors_from_gff(_decorin_vs_biglycan_gff(tmp_path), padding_bp=0)
    # The decorin neighbourhood has more synteny support -> it must be the better anchor.
    dec = c._best_block_anchor("chr5", 21_000_000, 21_040_000, anchors)
    big = c._best_block_anchor("chrX", 36_736_000, 36_741_000, anchors)
    assert dec["flanking_support"] > big["flanking_support"]


def test_no_block_tags_is_noop(tmp_path):
    # Legacy GFF without _b<N>_ ids -> no anchors (feature degrades to legacy behaviour).
    rows = [
        "chr5\texon_annotation\tmRNA\t21000000\t21040000\t90.0\t-\t.\t"
        "ID=GOI_DCN|legacy_model;Name=GOI_DCN;SynVoyRole=goi;Confidence=HIGH;Identity=90.0",
    ]
    assert c.load_goi_block_anchors_from_gff(_write(tmp_path, rows)) == []


def test_block_without_goi_is_skipped(tmp_path):
    # A block with only flanking (no GOI model) is not an anchor.
    rows = [
        _flank("chr3", 1_000, 2_000, 7, "gene-A"),
        _flank("chr3", 3_000, 4_000, 7, "gene-B"),
    ]
    assert c.load_goi_block_anchors_from_gff(_write(tmp_path, rows)) == []
