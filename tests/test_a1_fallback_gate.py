#!/usr/bin/env python3
"""
A1 "the fumble" regression tests (TODO_JUN §A1).

The GOI fallback gate used to veto on query COVERAGE, which is structurally low
for a short precursor whose conserved mature peptide is only a fraction of the
query. A real high-bitscore divergent hit (Xylocopa melittin b0: 50% / bits≈91)
was rejected at qcov≈0.18 while a spurious high-coverage/low-bitscore artifact in
another block (b1: 28% / qcov 0.40) passed — so the GOI was emitted in the WRONG
block. `is_valid_fallback` now gates short queries on aligned length + bitscore.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from iterative_search_runner import is_valid_fallback

# Defaults mirror nextflow.config / argparse.
SHORT_LEN = 150
MIN_ALN = 15
MIN_BITS = 30.0


def gate(ordered_hit_count, qcov, single_aln, aln_total, best_bits, query_len,
         short_query_len=SHORT_LEN, short_min_aln_aa=MIN_ALN, short_min_bits=MIN_BITS):
    return is_valid_fallback(
        ordered_hit_count=ordered_hit_count, qcov=qcov, single_aln=single_aln,
        aln_total=aln_total, best_bits=best_bits, query_len=query_len,
        short_query_len=short_query_len, short_min_aln_aa=short_min_aln_aa,
        short_min_bits=short_min_bits,
    )


class TestShortQueryFallback(unittest.TestCase):
    def test_xylopin_b0_recovered(self):
        """The real Xylocopa melittin b0 (70-aa query, 26-aa hit, bits 91) now passes."""
        self.assertTrue(gate(1, qcov=0.18, single_aln=26, aln_total=26, best_bits=91, query_len=70))

    def test_legacy_gate_would_have_rejected_b0(self):
        """With the kill switch (short_query_len=0) b0 is rejected as before — proves the fix matters."""
        self.assertFalse(
            gate(1, qcov=0.18, single_aln=26, aln_total=26, best_bits=91, query_len=70,
                 short_query_len=0)
        )

    def test_micro_window_rejected(self):
        """A 7-aa micro-window (the Polistes-style overcall) is still rejected (below aln floor)."""
        self.assertFalse(gate(1, qcov=0.10, single_aln=7, aln_total=7, best_bits=20, query_len=70))

    def test_low_bitscore_rejected(self):
        """A decently-long but weak (low-bitscore) short-query hit is rejected."""
        self.assertFalse(gate(1, qcov=0.40, single_aln=30, aln_total=30, best_bits=20, query_len=70))

    def test_tetramorium_multi_hit_unchanged(self):
        """Tetramorium OV788322 (76% cov, multi-exon, strong bits) passes — must not regress."""
        self.assertTrue(gate(2, qcov=0.76, single_aln=0, aln_total=53, best_bits=120, query_len=70))

    def test_multi_hit_short_aln_total_floor(self):
        """A short multi-hit fallback below the aligned-length floor is rejected."""
        self.assertFalse(gate(3, qcov=0.30, single_aln=0, aln_total=10, best_bits=50, query_len=70))


class TestLongQueryFallbackUnchanged(unittest.TestCase):
    def test_tp53_low_coverage_single_hit_rejected(self):
        """Long query (TP53, 393 aa): a low-coverage single hit is still vetoed by qcov."""
        self.assertFalse(gate(1, qcov=0.10, single_aln=40, aln_total=40, best_bits=90, query_len=393))

    def test_long_query_high_coverage_passes(self):
        """Long query with high coverage passes the legacy single-hit gate."""
        self.assertTrue(gate(1, qcov=0.80, single_aln=300, aln_total=300, best_bits=400, query_len=393))

    def test_long_query_short_single_hit_rejected(self):
        """Long query, single 24-aa hit (<25 aln floor) rejected exactly as before."""
        self.assertFalse(gate(1, qcov=0.30, single_aln=24, aln_total=24, best_bits=90, query_len=393))


if __name__ == '__main__':
    unittest.main()
