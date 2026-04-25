#!/usr/bin/env python3
"""Regression tests for bin/generate_report.py staging diagnostics and fail-loud behavior.

Covers two scenarios that the v2/v3/v4 melittin runs hit in production:
  1. Multi-locus staging with many GFFs + scores across several genomes.
     The report must accurately reflect totals and include staging_diagnostics
     with pattern_match_count > 0 for regions_gff and scores.
  2. Empty staged_results (plumbing failure upstream).
     Running bin/generate_report.py must exit non-zero with a diagnostic
     message that names the empty directories. The written report must
     carry staging_diagnostics.empty == True.

The --allow-empty flag should suppress the non-zero exit for legitimate
zero-hit runs.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest


BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

from generate_report import build_report  # noqa: E402


REPORT_SCRIPT = os.path.abspath(os.path.join(BIN_DIR, "generate_report.py"))


def _gff_lines(genome, locus_idx, n_annotations, base_pos=1000, role_mix=("goi", "flanking")):
    """Emit tab-separated GFF mRNA lines for a locus."""
    lines = ["##gff-version 3"]
    for i in range(n_annotations):
        role = role_mix[i % len(role_mix)]
        ident_prefix = "GOI_gene" if role == "goi" else f"gene{i}"
        start = base_pos + i * 200
        end = start + 150
        attrs = (
            f"ID={ident_prefix}|l{locus_idx}_{i};"
            f"Name={ident_prefix};"
            f"SynVoy_Parent={ident_prefix};"
            f"SynVoyRole={role};"
            f"Confidence=HIGH;"
            + ("GOIClass=confident_goi;" if role == "goi" else "")
            + "EvidenceType=exon_annotation;Identity=90.0;Exons=3"
        )
        lines.append(
            "\t".join(
                [f"chr{locus_idx}", "exon_annotation", "mRNA",
                 str(start), str(end), "90.0", "+", ".", attrs]
            )
        )
    return "\n".join(lines) + "\n"


def _scores_tsv(locus_idx, n_rows):
    header = (
        "region_rank\tregion_name\tchrom\tstart\tend\tstrand\tscore\t"
        "quality_score\tcoverage_score\tunique_genes\ttotal_genes_expected\t"
        "consistency\tstrand_consistency\tp_value\tgoi_overlap\tis_goi_anchor\t"
        "confidence\tselection_reason"
    )
    rows = [header]
    for i in range(n_rows):
        rows.append(
            "\t".join([
                str(i + 1),
                f"Reg{locus_idx}_{i}",
                f"chr{locus_idx}",
                str(i * 1000),
                str(i * 1000 + 900),
                "+",
                "0.85",
                "0.80",
                "0.90",
                "5",
                "7",
                "0.7",
                "1.0",
                "0.01",
                "true",
                "true",
                "HIGH",
                "best_score",
            ])
        )
    return "\n".join(rows) + "\n"


class TestMultiLocusStagingDiagnostics(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        for sub in ("regions", "hits", "scores"):
            os.makedirs(os.path.join(self.test_dir, sub), exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def _write(self, rel_path, content):
        path = os.path.join(self.test_dir, rel_path)
        with open(path, "w") as fh:
            fh.write(content)
        return path

    def test_multi_locus_report_has_diagnostics_and_non_zero_totals(self):
        genomes = ["GCF_000001.1.fna", "GCF_000002.1.fna", "GCF_000003.1.fna"]
        for g_idx, genome in enumerate(genomes, start=1):
            for locus_idx in range(1, 4):
                self._write(
                    f"regions/{genome}.locus{locus_idx}.gff",
                    _gff_lines(genome, locus_idx, n_annotations=4),
                )
                self._write(
                    f"scores/{genome}.locus{locus_idx}.scores.tsv",
                    _scores_tsv(locus_idx, n_rows=2),
                )
            self._write(
                f"regions/{genome}.faa",
                ">p1\nMKTA\n>p2\nMKTB\n>p3\nMKTC\n",
            )
            self._write(
                f"hits/synteny_block_locus_1_{genome}.m8",
                "q1\tchr1\t90\t100\t0\t0\t1\t100\t10\t110\t1e-40\t200\n",
            )

        report = build_report(self.test_dir)

        self.assertIn("staging_diagnostics", report)
        diag = report["staging_diagnostics"]
        self.assertFalse(diag["empty"])
        self.assertEqual(diag["match_counts"]["gff_files"], 9)
        self.assertEqual(diag["match_counts"]["score_files"], 9)
        self.assertEqual(diag["match_counts"]["fasta_files"], 3)
        self.assertEqual(diag["match_counts"]["hit_files"], 3)
        self.assertTrue(diag["dirs"]["regions_gff"]["exists"])
        self.assertGreaterEqual(diag["dirs"]["regions_gff"]["entry_count"], 9)
        self.assertGreater(len(diag["dirs"]["regions_gff"]["sample_matches"]), 0)

        self.assertGreater(report["summary"]["total_annotations"], 0)
        self.assertGreater(report["summary"]["total_goi_annotations"], 0)
        self.assertFalse(report["summary"]["staging_empty"])
        self.assertEqual(report["regions"]["total_regions"], 18)

        goi_class_counts = report["annotations"]["goi_class_counts"]
        self.assertIn("confident_goi", goi_class_counts)
        self.assertGreater(goi_class_counts["confident_goi"], 0)


class TestEmptyStagingFailsLoud(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        for sub in ("regions", "hits", "scores"):
            os.makedirs(os.path.join(self.test_dir, sub), exist_ok=True)
        self.output = os.path.join(self.test_dir, "synvoy_report.json")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_empty_staging_emits_diagnostics_and_exits_nonzero(self):
        proc = subprocess.run(
            [sys.executable, REPORT_SCRIPT,
             "--results_dir", self.test_dir,
             "--output", self.output],
            capture_output=True, text=True,
        )

        self.assertNotEqual(proc.returncode, 0,
                            msg=f"Expected non-zero exit. stderr={proc.stderr}")
        self.assertIn("zero annotation and zero region files", proc.stderr)
        self.assertIn("channel wiring is broken", proc.stderr)

        self.assertTrue(os.path.exists(self.output))
        with open(self.output) as fh:
            report = json.load(fh)
        self.assertTrue(report["staging_diagnostics"]["empty"])
        self.assertTrue(report["summary"]["staging_empty"])
        self.assertEqual(report["summary"]["total_annotations"], 0)

    def test_allow_empty_flag_suppresses_nonzero_exit(self):
        proc = subprocess.run(
            [sys.executable, REPORT_SCRIPT,
             "--results_dir", self.test_dir,
             "--output", self.output,
             "--allow-empty"],
            capture_output=True, text=True,
        )

        self.assertEqual(proc.returncode, 0,
                         msg=f"Expected zero exit with --allow-empty. stderr={proc.stderr}")
        self.assertIn("Warning (--allow-empty set)", proc.stderr)

        with open(self.output) as fh:
            report = json.load(fh)
        self.assertTrue(report["staging_diagnostics"]["empty"])


if __name__ == "__main__":
    unittest.main()
