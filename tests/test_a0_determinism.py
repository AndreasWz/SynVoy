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


# The real LRZ job 5778368 inputs: the home melittin and the A. florea model the
# wavefront adds in wave 2 are BOTH 70 aa and share parent `GOI_Melt`.
_HOME = {'id': 'GOI_Melt', 'seq': 'M' * 70}
_FLOREA = {'id': 'GOI_Melt|Apis_florea_fna_b0_l1_exon_ann', 'seq': 'K' * 70}


class TestGoiRepresentativeCollapse(unittest.TestCase):
    """Hash-order source #2: which GOI query represents a parent (LRZ job 5778368).

    The parent collapse kept the first-seen of two equal-length representatives, and
    "first" came from a set of strings, rebuilt per search block. With PYTHONHASHSEED
    unpinned, every block searched after the wave-2 expansion flipped between the two
    melittins: 17 of 20 GFFs differed between identical runs -- exactly the
    post-expansion genomes; the 3 identical files were A. cerana, A. florea (waves 1-2)
    and the home. LRZ R1 was a mixed draw (Tetramorium b0 florea, b1 home).
    """

    def setUp(self):
        from iterative_search_runner import collapse_goi_queries_by_parent
        self.collapse = collapse_goi_queries_by_parent

    def test_home_query_wins_a_length_tie_in_either_order(self):
        for order in ([_HOME, _FLOREA], [_FLOREA, _HOME]):
            self.assertIs(self.collapse(order)['GOI_Melt'], _HOME)

    def test_longer_expansion_model_still_wins(self):
        """The tie-break must not change the documented longest-wins rule."""
        longer = {'id': 'GOI_Melt|Apis_florea_fna_b0_l1_exon_ann', 'seq': 'K' * 71}
        self.assertIs(self.collapse([_HOME, longer])['GOI_Melt'], longer)

    def test_tie_between_expansion_models_takes_smallest_id(self):
        a = {'id': 'GOI_Melt|Bombus_b1_l1', 'seq': 'A' * 70}
        b = {'id': 'GOI_Melt|Apis_b0_l1', 'seq': 'B' * 70}
        self.assertIs(self.collapse([a, b])['GOI_Melt'], b)
        self.assertIs(self.collapse([b, a])['GOI_Melt'], b)

    def test_parents_come_out_in_a_fixed_order(self):
        x = {'id': 'GOI_X', 'seq': 'A' * 10}
        y = {'id': 'GOI_Y', 'seq': 'A' * 10}
        self.assertEqual(list(self.collapse([y, x])), ['GOI_X', 'GOI_Y'])


_CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
from iterative_search_runner import collapse_goi_queries_by_parent
ids = ['GOI_Melt', 'GOI_Melt|Apis_florea_fna_b0_l1_exon_ann', 'GOI_Melt|exon_1',
       'gene-LOC726866', 'gene-LOC409662', 'gene-LOC726827']
seqs = {i: ('M' if i == 'GOI_Melt' else 'K') * (70 if '|exon_' not in i else 39) for i in ids}
unordered = list(set(ids))          # hash-seed dependent, like `unique_queries`
goi = [{'id': i, 'seq': seqs[i]} for i in unordered if i.startswith('GOI_') and '|exon_' not in i]
print(','.join(unordered), collapse_goi_queries_by_parent(goi)['GOI_Melt']['id'])
"""


class TestHashSeedInvariance(unittest.TestCase):
    """The collapse must give one answer under every PYTHONHASHSEED.

    Runs in fresh interpreters because the hash seed is fixed at startup. The control
    assertion proves the test has teeth: the raw set order DOES change across these
    seeds, so a constant chosen id means the code is order-independent, not lucky.
    """

    def test_representative_is_independent_of_hash_seed(self):
        bindir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin')
        orders, chosen = set(), set()
        for seed in range(12):
            env = dict(os.environ, PYTHONHASHSEED=str(seed))
            out = subprocess.run([sys.executable, '-c', _CHILD, bindir], env=env,
                                 capture_output=True, text=True, check=True).stdout.split()
            orders.add(out[0])
            chosen.add(out[1])
        self.assertGreater(len(orders), 1, "control: set order should vary with the seed")
        self.assertEqual(chosen, {'GOI_Melt'})


class TestPipelinePinsHashSeed(unittest.TestCase):
    def test_nextflow_config_sets_pythonhashseed(self):
        cfg = os.path.join(os.path.dirname(__file__), '..', 'nextflow.config')
        with open(cfg) as fh:
            self.assertRegex(fh.read(), r"(?m)^\s*PYTHONHASHSEED\s*=\s*'0'")


if __name__ == '__main__':
    unittest.main()
