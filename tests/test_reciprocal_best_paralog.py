"""Tests for the reciprocal-best check (bin/reciprocal_best_paralog_check.py), which locus
ownership (§1m, ASSIGN_LOCUS_OWNERSHIP) runs with the home-paralog panel.

1. The per-(locus, genome) script: the single-sequence no-op, the GFF/.faa join, and a smoke
   that a multi-sequence query produces well-formed rows with non-trivial bitscore gaps.

2. The report reads its tables for locus ownership only: the reciprocal-best paralog check of
   the query file (RECIPROCAL_BEST_PARALOG) and its ``paralog_confusion`` flag were removed on
   2026-11-16.
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
# 2) The report reads the tables of locus ownership only
# ──────────────────────────────────────────────────────────────────────

def test_loader_skips_header_only_files(tmp_path: Path) -> None:
    """A one-sequence panel leaves a header-only TSV; the loader yields zero rows from it."""
    pdir = tmp_path / "locus_ownership"
    pdir.mkdir()
    header = "\t".join(rbp.HEADER) + "\n"
    (pdir / "locus_1__GenomeA.locus_ownership.tsv").write_text(header)
    (pdir / "locus_1__GenomeB.locus_ownership.tsv").write_text(
        header + "\t".join([
            "GenomeB", "locus_1", "chr1", "100", "200", "GOI_m1", "50",
            "locus_1|gene-a", "800", "92.0", "locus_2|gene-b", "750", "50.0", "2",
        ]) + "\n"
    )
    rows = gr.load_paralog_check_rows(str(pdir))
    assert len(rows) == 1
    assert rows[0]["best_paralog"] == "locus_1|gene-a"
    assert rows[0]["genome"] == "GenomeB"


def test_the_paralog_confusion_path_is_gone(tmp_path: Path) -> None:
    """Review fix (2026-11-17): RECIPROCAL_BEST_PARALOG went on 2026-11-16 (it always received one
    query sequence and could flag nothing), but the report still took --paralog_confusion_min_gap,
    read paralog_check/ and wrote the check's keys. They go with it."""
    import subprocess
    assert not hasattr(gr, "build_paralog_confusion_flags")
    assert "RECIPROCAL_BEST_PARALOG" not in gr._STEP_FAILURE_SENTENCE
    res = gr.build_self_consistency(goi_annotations=[], goi_dedup={"records": []}, flanking_per_locus={})
    assert "paralog_confusion" not in res["checks_performed"] + res["deferred_checks"]
    assert "paralog_confusion_min_gap" not in res["thresholds"]
    assert not {"paralog_alignment_table", "modal_paralog_per_locus"} & set(res)
    assert "n_paralog_confusion" not in res["summary"]
    for sub in ("regions", "hits", "scores"):
        (tmp_path / sub).mkdir()
    report = gr.build_report(str(tmp_path))
    assert "paralog_confusion_flags" not in report["summary"]
    out = subprocess.run([sys.executable, str(ROOT / "bin" / "generate_report.py"), "--results_dir",
                          str(tmp_path), "--output", str(tmp_path / "r.json"),
                          "--paralog_confusion_min_gap", "5"], capture_output=True, text=True)
    assert out.returncode == 2 and "unrecognized arguments" in out.stderr
