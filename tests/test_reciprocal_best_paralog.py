"""Tests for docs/TODO.md §1j Phase B — reciprocal-best paralog check.

Covers both halves of the implementation:

1. ``bin/reciprocal_best_paralog_check.py`` — the per-(locus, genome) script:
   the single-paralog no-op, the GFF/.faa join, and a smoke that a multi-paralog
   query produces well-formed rows with non-trivial bitscore gaps.

2. ``generate_report.build_paralog_confusion_flags`` — aggregation across (locus,
   genome) cells: modal-paralog inference, ``paralog_confusion`` flag when a row
   diverges from the modal best with a sufficient bitscore gap.
"""
from __future__ import annotations

import csv
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import reciprocal_best_paralog_check as rbp  # noqa: E402
import generate_report as gr  # noqa: E402


# ──────────────────────────────────────────────────────────────────────
# 1) Per-(locus, genome) script
# ──────────────────────────────────────────────────────────────────────

def _write(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip("\n"))


def test_single_paralog_query_is_a_noop(tmp_path: Path, monkeypatch) -> None:
    """Single-sequence home query → script emits header-only TSV; aligner is
    never called."""
    _write(tmp_path / "home.faa", """
        >TP53
        MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQK
    """)
    _write(tmp_path / "target.faa", """
        >GOI_target
        MEEPQSDPSVEPPLSQETFSDLW
    """)
    _write(tmp_path / "target.gff", """
        ##gff-version 3
        chrX\tSynVoy\tmRNA\t100\t200\t.\t+\t.\tID=GOI_target;SynVoyRole=goi;Confidence=HIGH
    """)
    out = tmp_path / "out.tsv"

    called = {"n": 0}
    def _trap(*a, **kw):
        called["n"] += 1
        return 0.0, 0.0
    monkeypatch.setattr(rbp, "_align", _trap)
    monkeypatch.setattr(sys, "argv", [
        "reciprocal_best_paralog_check.py",
        "--home_query", str(tmp_path / "home.faa"),
        "--target_faa", str(tmp_path / "target.faa"),
        "--target_gff", str(tmp_path / "target.gff"),
        "--locus_id", "locus_1",
        "--genome_name", "TargetA",
        "--output", str(out),
    ])
    rbp.main()
    assert called["n"] == 0, "must not invoke aligner on single-paralog input"
    lines = out.read_text().splitlines()
    assert lines[0] == "\t".join(rbp.HEADER)
    assert len(lines) == 1


def test_multi_paralog_writes_rows_and_picks_best(tmp_path: Path) -> None:
    """Multi-paralog query: identical-to-TP53 target should align best to TP53
    with a positive bitscore gap to TP63. Smoke uses real parasail."""
    # Two synthetic paralogs that share a short conserved region but differ
    # elsewhere — TP53 and TP63 are loosely modeled.
    _write(tmp_path / "home.faa", """
        >TP53
        MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQK
        >TP63
        MNFETSRCATLQYCPDPYIQRFVETPAHFSWKESYYRSTMSQSTQTNEFLSPEVFQHIWDFLEQPICSVQPIDLNFVDEPSEDGATNKIEISMDCIRMQDSDLSDPMWPQYTNLGLLNSMDQQIQNGSSSTSPYNTDHAQNSVTAPSPYAQPSSTFDALSPSPAIPSNT
    """)
    # Target = the TP53 sequence verbatim — best match must be TP53.
    _write(tmp_path / "target.faa", """
        >GOI_target_1
        MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQK
    """)
    _write(tmp_path / "target.gff", """
        ##gff-version 3
        chrX\tSynVoy\tmRNA\t100\t500\t.\t+\t.\tID=GOI_target_1;SynVoyRole=goi;Confidence=HIGH
    """)
    out = tmp_path / "out.tsv"

    import argparse
    args = argparse.Namespace(
        home_query=str(tmp_path / "home.faa"),
        target_faa=str(tmp_path / "target.faa"),
        target_gff=str(tmp_path / "target.gff"),
        locus_id="locus_1", genome_name="TargetA", output=str(out),
    )
    # Easier than monkeypatching argv: invoke main() with patched sys.argv.
    saved = sys.argv
    sys.argv = ["reciprocal_best_paralog_check.py"] + [
        "--home_query", args.home_query, "--target_faa", args.target_faa,
        "--target_gff", args.target_gff, "--locus_id", args.locus_id,
        "--genome_name", args.genome_name, "--output", args.output,
    ]
    try:
        rbp.main()
    finally:
        sys.argv = saved

    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert len(rows) == 1
    r = rows[0]
    assert r["genome"] == "TargetA"
    assert r["locus"] == "locus_1"
    assert r["best_paralog"] == "TP53"
    assert float(r["best_bit"]) > float(r["second_bit"])
    assert float(r["bitscore_gap"]) > 0


# ──────────────────────────────────────────────────────────────────────
# 2) Aggregation in generate_report.build_paralog_confusion_flags
# ──────────────────────────────────────────────────────────────────────

def _row(genome, locus, chrom, start, end, best, best_bit, second_bit):
    return {
        "genome": genome, "locus": locus, "chrom": chrom,
        "start": str(start), "end": str(end),
        "mrna_id": f"GOI_{chrom}_{start}",
        "best_paralog": best, "best_bit": str(best_bit),
        "second_paralog": "OTHER", "second_bit": str(second_bit),
        "bitscore_gap": str(best_bit - second_bit),
        # post-processed fields the loader adds:
        "best_bit_f": best_bit, "second_bit_f": second_bit,
        "gap_f": best_bit - second_bit,
        "start_i": start, "end_i": end,
    }


def test_no_paralog_rows_means_no_flags() -> None:
    flags, per_call, modal = gr.build_paralog_confusion_flags([], min_gap=5.0)
    assert flags == [] and per_call == [] and modal == {}


def test_modal_paralog_per_locus_emits_confusion_for_outliers() -> None:
    """Locus has 5 calls best-matching TP63 and 1 call best-matching TP73 with
    a 10-bit gap → 1 paralog_confusion flag; ties don't fire."""
    rows = [
        _row("GenomeA", "locus_3", "chr3", 100, 200, "TP63", 900, 850),
        _row("GenomeA", "locus_3", "chr3", 300, 400, "TP63", 880, 800),
        _row("GenomeA", "locus_3", "chr3", 500, 600, "TP63", 870, 810),
        _row("GenomeA", "locus_3", "chr3", 700, 800, "TP63", 860, 820),
        _row("GenomeA", "locus_3", "chr3", 900, 1000, "TP63", 855, 825),
        _row("GenomeA", "locus_3", "chr3", 1100, 1200, "TP73", 800, 790),  # gap=10, FLAG
    ]
    flags, per_call, modal = gr.build_paralog_confusion_flags(rows, min_gap=5.0)
    assert modal[("GenomeA", "locus_3")] == "TP63"
    assert len(flags) == 1
    f = flags[0]
    assert f["type"] == "paralog_confusion"
    assert f["best_match_paralog"] == "TP73"
    assert f["assigned_paralog"] == "TP63"
    assert f["bitscore_gap"] == 10
    # All 6 calls present in per_call summary, 5 marked not-confused.
    assert sum(1 for c in per_call if c["confused"]) == 1
    assert sum(1 for c in per_call if not c["confused"]) == 5


def test_min_gap_suppresses_borderline_calls() -> None:
    """Same outlier as above but with a 2-bit gap — should NOT flag at min_gap=5."""
    rows = [
        _row("GenomeA", "locus_3", "chr3", 100, 200, "TP63", 900, 850),
        _row("GenomeA", "locus_3", "chr3", 300, 400, "TP63", 880, 800),
        _row("GenomeA", "locus_3", "chr3", 500, 600, "TP63", 870, 810),
        _row("GenomeA", "locus_3", "chr3", 700, 800, "TP73", 800, 798),  # gap=2
    ]
    flags, _, _ = gr.build_paralog_confusion_flags(rows, min_gap=5.0)
    assert flags == []


def test_loader_skips_header_only_files(tmp_path: Path) -> None:
    """Single-paralog runs leave a header-only TSV; loader must yield zero rows."""
    pdir = tmp_path / "paralog_check"
    pdir.mkdir()
    header = "\t".join(rbp.HEADER) + "\n"
    (pdir / "GenomeA.paralog_check.tsv").write_text(header)
    (pdir / "GenomeB.paralog_check.tsv").write_text(
        header + "\t".join([
            "GenomeB", "locus_1", "chr1", "100", "200", "GOI_m1", "50",
            "TP53", "800", "92.0", "TP63", "750", "50.0", "2",
        ]) + "\n"
    )
    rows = gr.load_paralog_check_rows(str(pdir))
    assert len(rows) == 1
    assert rows[0]["best_paralog"] == "TP53"
    assert rows[0]["genome"] == "GenomeB"


def test_self_consistency_threads_paralog_flags_through() -> None:
    """When paralog flags are provided, build_self_consistency surfaces them in
    flags + summary, and lists paralog_confusion as a check_performed."""
    paralog_flags = [{
        "type": "paralog_confusion",
        "genome": "G", "locus": "locus_3", "chrom": "c", "start": 1, "end": 2,
        "assigned_paralog": "TP63", "best_match_paralog": "TP73",
        "bitscore_gap": 8.0, "mrna_id": "GOI_x", "advice": "...",
    }]
    paralog_per_call = [{"genome": "G", "locus": "locus_3", "confused": True,
                         "best_paralog": "TP73", "bitscore_gap": 8.0,
                         "chrom": "c", "start": 1, "end": 2,
                         "mrna_id": "GOI_x", "best_bit": 800,
                         "second_paralog": "TP63", "second_bit": 792,
                         "locus_modal_paralog": "TP63"}]
    res = gr.build_self_consistency(
        goi_annotations=[], goi_dedup={"records": []}, flanking_per_locus={},
        paralog_flags=paralog_flags,
        paralog_per_call=paralog_per_call,
        paralog_modal_per_locus={("G", "locus_3"): "TP63"},
        paralog_min_gap=5.0,
    )
    assert "paralog_confusion" in res["checks_performed"]
    assert res["summary"]["n_paralog_confusion"] == 1
    assert res["summary"]["total_flags"] == 1
    assert res["flags"][0]["type"] == "paralog_confusion"
    assert res["modal_paralog_per_locus"] == [
        {"genome": "G", "locus": "locus_3", "paralog": "TP63"}
    ]
