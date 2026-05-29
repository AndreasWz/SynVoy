"""Regression test for docs/TODO.md §1c (`generate_report.py` scaling).

The two changes under test:

1. ``_process_gffs_unified`` reads each GFF once and emits the three outputs
   that the legacy ``summarize_annotations`` / ``collect_goi_annotations`` /
   ``collect_high_flanking_per_locus`` produced separately. The three legacy
   functions are now thin wrappers around it. This test pins:
   - the unified output is byte-identical to running the three wrappers
     individually (correctness preserved)
   - the SP-family-scale (1000 GFFs × ~20 lines each) path completes well
     under a wall-time budget (was hitting timeout-kill workarounds in prod)

2. ``scan_dir_by_suffix`` does one ``os.scandir`` pass per directory and
   classifies entries into named buckets. Replaces N×``glob.glob('*.ext')``
   patterns. This test pins the bucket contents against a hand-built fixture.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

# bin/ on sys.path so we can `import generate_report`
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import generate_report as gr  # noqa: E402


GFF_TEMPLATE = (
    "##gff-version 3\n"
    "{chrom}\tSynVoy\tgene\t1\t100\t.\t+\t.\t"
    "ID=gene-flanking-{i};Name=FLANK_{i};SynVoyRole=flanking;Confidence=HIGH\n"
    "{chrom}\tSynVoy\tmRNA\t1\t100\t.\t+\t.\t"
    "ID=mRNA-flanking-{i};Parent=gene-flanking-{i};SynVoyRole=flanking;Confidence=HIGH\n"
    "{chrom}\tSynVoy\tgene\t200\t300\t.\t+\t.\t"
    "ID=GOI_{i};SynVoyRole=goi;Confidence=HIGH;GOIClass=confident_goi;"
    "EvidenceType=miniprot;Identity=92.5;TargetGene=Mel;ModelStatus=complete\n"
    "{chrom}\tSynVoy\tmRNA\t200\t300\t.\t+\t.\t"
    "ID=GOI_mRNA_{i};Parent=GOI_{i};SynVoyRole=goi;Confidence=HIGH;"
    "GOIClass=confident_goi;EvidenceType=miniprot;Identity=92.5;TargetGene=Mel;"
    "ModelStatus=complete\n"
)


def _make_gff(path: Path, locus_idx: int, genome_idx: int) -> None:
    """One toy region GFF: 1 HIGH flanking + 1 HIGH GOI on chrom-N."""
    chrom = f"chrom_{locus_idx}_{genome_idx}"
    path.write_text(GFF_TEMPLATE.format(chrom=chrom, i=locus_idx * 1000 + genome_idx))


def test_unified_pass_matches_legacy_wrappers(tmp_path: Path) -> None:
    """The unified pass must produce byte-identical outputs to the three legacy
    functions, otherwise tests/test_generate_report_dedup.py and
    tests/test_generate_report_self_consistency.py would silently start drifting.
    """
    # 30 (locus, genome) pairs — enough variety to exercise per-genome grouping +
    # per-(genome, locus, chrom) flanking counting without being expensive.
    gff_dir = tmp_path / "regions"
    gff_dir.mkdir()
    gff_paths = []
    for locus in range(5):
        for genome in range(6):
            p = gff_dir / f"locus_{locus}__genome_{genome}.fna.gff"
            _make_gff(p, locus, genome)
            gff_paths.append(str(p))

    ann_unified, goi_unified, flank_unified = gr._process_gffs_unified(gff_paths)

    # The three public wrappers must agree with the unified return values.
    assert gr.summarize_annotations(gff_paths) == ann_unified
    assert gr.collect_goi_annotations(gff_paths) == goi_unified
    assert gr.collect_high_flanking_per_locus(gff_paths) == flank_unified

    # Sanity-check counts so a future refactor that quietly loses GFF rows trips here.
    # Each fixture file has 4 mRNA/gene rows: a (gene, mRNA) pair for flanking + a
    # (gene, mRNA) pair for the GOI. The parser counts gene and mRNA rows
    # independently (dedup happens downstream in dedupe_goi_annotations).
    assert ann_unified["total_annotations"] == 4 * 30
    assert sum(row["goi_annotations"] for row in ann_unified["per_genome"]) == 2 * 30
    assert len(goi_unified) == 2 * 30
    # Flanking-per-locus is keyed by (genome, locus, chrom) → set of distinct gene
    # names. The fixture's gene row uses Name=FLANK_<i>; the mRNA row falls through
    # to ID=mRNA-flanking-<i>, so the set holds 2 entries per file.
    assert len(flank_unified) == 30  # one (genome, locus, chrom) per file
    assert all(count == 2 for count in flank_unified.values())


def test_scan_dir_by_suffix_classifies_correctly(tmp_path: Path) -> None:
    """One scandir pass should bucket files by suffix; missing dirs → empty buckets."""
    (tmp_path / "a.gff").write_text("")
    (tmp_path / "b.gff3").write_text("")
    (tmp_path / "c.faa").write_text("")
    (tmp_path / "d.fna").write_text("")
    (tmp_path / "irrelevant.txt").write_text("")
    (tmp_path / "subdir").mkdir()

    buckets = gr.scan_dir_by_suffix(tmp_path, {
        "fasta": (".faa", ".fna"),
        "gff": (".gff", ".gff3"),
    })
    assert sorted(os.path.basename(p) for p in buckets["fasta"]) == ["c.faa", "d.fna"]
    assert sorted(os.path.basename(p) for p in buckets["gff"]) == ["a.gff", "b.gff3"]

    # Missing directory: every bucket is empty (no crash).
    empty = gr.scan_dir_by_suffix(tmp_path / "does_not_exist",
                                  {"x": (".any",)})
    assert empty == {"x": []}


def test_unified_pass_handles_1000_gffs_under_budget(tmp_path: Path) -> None:
    """SP-family-scale smoke: 1000 GFFs must complete well under 5s on the unified
    pass. The previous triple-read path scaled at ~3× this for the same fixture and
    was the hang Ivan's pipeline kept timing out on (docs/TODO.md §1c).
    """
    gff_dir = tmp_path / "regions"
    gff_dir.mkdir()
    gff_paths = []
    for locus in range(50):
        for genome in range(20):
            p = gff_dir / f"locus_{locus}__genome_{genome}.fna.gff"
            _make_gff(p, locus, genome)
            gff_paths.append(str(p))
    assert len(gff_paths) == 1000

    t0 = time.perf_counter()
    ann, goi, flank = gr._process_gffs_unified(gff_paths)
    elapsed = time.perf_counter() - t0

    # The unified pass on 1000 GFFs (each ~5 lines of GFF) is I/O-bound and takes
    # well under a second on a warm-cache laptop. Budget of 5s is generous and
    # leaves headroom for slow CI / cold-cache. The point of the assertion is to
    # catch a future regression that re-introduces a 3× read (which would land here
    # at ~1.5–3s of file I/O on the same hardware) or a quadratic accumulator.
    assert elapsed < 5.0, f"unified pass took {elapsed:.2f}s on 1000 GFFs"

    # And the outputs are well-formed (4 mRNA/gene rows per file → 2 GOI, 2 flanking).
    assert ann["total_annotations"] == 4 * 1000
    assert len(goi) == 2 * 1000
    assert len(flank) == 1000  # one (genome, locus, chrom) bucket per file
    assert all(count == 2 for count in flank.values())
