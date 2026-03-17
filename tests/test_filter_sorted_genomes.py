#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPT = os.path.join(ROOT, "bin", "filter_sorted_genomes.py")


class TestFilterSortedGenomes(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _run(self, sorted_lines, qc_records, policy="drop"):
        sorted_path = os.path.join(self.test_dir, "sorted.tsv")
        qc_path = os.path.join(self.test_dir, "qc.json")
        output_path = os.path.join(self.test_dir, "filtered.tsv")

        with open(sorted_path, "w") as fh:
            fh.write("\n".join(sorted_lines) + "\n")
        with open(qc_path, "w") as fh:
            json.dump(qc_records, fh)

        proc = subprocess.run(
            [
                sys.executable,
                SCRIPT,
                "--sorted", sorted_path,
                "--qc_json", qc_path,
                "--output", output_path,
                "--policy", policy,
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        with open(output_path) as fh:
            output_lines = [line.strip() for line in fh if line.strip()]
        return output_lines, proc.stderr

    def test_drop_failed_and_keep_unknown(self):
        output_lines, stderr = self._run(
            [
                "GCF_000001.1.fna\t0.10",
                "GCF_000002.1.fna\t0.20",
                "GCF_000003.1.fna\t0.30",
            ],
            [
                {"genome": "GCF_000001.1.fna", "status": "PASS", "msg": "OK"},
                {"genome": "GCF_000002.1.fna", "status": "FAIL", "msg": "too fragmented"},
            ],
            policy="drop",
        )

        self.assertEqual(
            output_lines,
            ["GCF_000001.1.fna\t0.10", "GCF_000003.1.fna\t0.30"],
        )
        self.assertIn("kept=2 dropped=1 unknown=1", stderr)
        self.assertIn("without QC record", stderr)

    def test_keep_policy_preserves_failed_genomes(self):
        output_lines, stderr = self._run(
            [
                "GCF_000001.1.fna\t0.10",
                "GCF_000002.1.fna\t0.20",
            ],
            [
                {"genome": "GCF_000001.1.fna", "status": "PASS", "msg": "OK"},
                {"genome": "GCF_000002.1.fna", "status": "FAIL", "msg": "too fragmented"},
            ],
            policy="keep",
        )

        self.assertEqual(
            output_lines,
            ["GCF_000001.1.fna\t0.10", "GCF_000002.1.fna\t0.20"],
        )
        self.assertIn("policy=keep", stderr)

    def test_all_fail_falls_back_to_original_order(self):
        sorted_lines = [
            "GCF_000001.1.fna\t0.10",
            "GCF_000002.1.fna\t0.20",
        ]
        output_lines, stderr = self._run(
            sorted_lines,
            [
                {"genome": "GCF_000001.1.fna", "status": "FAIL", "msg": "bad"},
                {"genome": "GCF_000002.1.fna", "status": "FAIL", "msg": "bad"},
            ],
            policy="drop",
        )

        self.assertEqual(output_lines, sorted_lines)
        self.assertIn("All genomes failed QC", stderr)


if __name__ == "__main__":
    unittest.main()
