#!/usr/bin/env python3
"""
A0 reproducibility regression tests (TODO_JUN §A0).

The iterative search was non-deterministic run-to-run: mmseqs --threads>1
jitters marginal hits and the per-wave expansion DB was assembled in
as_completed() order, so a divergent GOI flipped confidence between identical
runs and cascaded through the wave seed set. These tests pin the two cheap,
always-on software fixes:

  * `parse_hits` imposes a total order on its output (task 2b), so the
    run-dependent .m8 row order can no longer flip downstream stable-sort
    tie-breaks.
  * the `--deterministic_goi_search` flag exists and defaults to True (task 2a),
    forcing the per-region augmented GOI mmseqs search to --threads 1.

The wave-assembly reordering (pinned source #1) lives inside main()'s wave loop
and is exercised by the end-to-end melittin rerun; here we cover the unit-level
pieces.
"""

import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from iterative_search_runner import parse_hits

SCRIPT = os.path.join(os.path.dirname(__file__), '..', 'bin', 'iterative_search_runner.py')

# query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
_ROWS = [
    "query1\tchr1\t90.0\t100\t10\t0\t1\t100\t1000\t1100\t1e-50\t200",
    "query5\tchr2\t85.0\t150\t22\t0\t1\t150\t5000\t5150\t1e-60\t250",
    "query3\tchr1\t88.0\t120\t14\t0\t1\t120\t3000\t3120\t1e-40\t180",
    "query2\tchr1\t92.0\t110\t9\t0\t1\t110\t2000\t2110\t1e-55\t210",
]


def _write(rows, path):
    with open(path, 'w') as fh:
        fh.write("\n".join(rows) + "\n")


class TestParseHitsDeterminism(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.dir, ignore_errors=True)

    def _parse(self, rows):
        f = os.path.join(self.dir, "h.m8")
        _write(rows, f)
        return parse_hits(f, min_identity=40.0, min_length=50, evalue_thresh=1e-5)

    def test_row_order_does_not_change_output_order(self):
        """Same hits in two different .m8 row orders → identical parsed order."""
        forward = self._parse(_ROWS)
        reversed_rows = self._parse(list(reversed(_ROWS)))
        shuffled = self._parse([_ROWS[2], _ROWS[0], _ROWS[3], _ROWS[1]])

        key = lambda hs: [(h['query'], h['chrom'], h['start']) for h in hs]
        self.assertEqual(key(forward), key(reversed_rows))
        self.assertEqual(key(forward), key(shuffled))

    def test_total_order_is_query_then_position(self):
        """Output is sorted by (query, chrom, start, ...) — fully deterministic."""
        hits = self._parse(list(reversed(_ROWS)))
        queries = [h['query'] for h in hits]
        self.assertEqual(queries, ['query1', 'query2', 'query3', 'query5'])

    def test_existing_two_hit_expectation_preserved(self):
        """The legacy test_core_functions expectation (query1 before query5) holds."""
        hits = self._parse([_ROWS[0], _ROWS[1]])
        self.assertEqual(hits[0]['query'], 'query1')
        self.assertEqual(hits[1]['query'], 'query5')


class TestDeterministicFlag(unittest.TestCase):
    def test_flag_present_in_help(self):
        """--deterministic_goi_search is a real CLI flag (and --help exits 0)."""
        out = subprocess.run(
            [sys.executable, SCRIPT, "--help"],
            capture_output=True, text=True, check=True,
        ).stdout
        self.assertIn("--deterministic_goi_search", out)

    def test_flag_defaults_true(self):
        """The flag defaults to True so determinism is on without extra config."""
        with open(SCRIPT) as fh:
            src = fh.read()
        idx = src.index('"--deterministic_goi_search"')
        decl = src[idx:idx + 200]
        self.assertIn("default=True", decl)


if __name__ == '__main__':
    unittest.main()
