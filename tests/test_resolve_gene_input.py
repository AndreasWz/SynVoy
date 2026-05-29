#!/usr/bin/env python3
"""Tests for input-type detection in bin/resolve_gene_input.py.

Regression: NCBI RefSeq/GenBank accessions without a version suffix (e.g.
XP_031329057, the exact form a user pastes from an NCBI protein page) were
misclassified as gene symbols because NCBI_REFSEQ_RE required '.\\d+' and matched
the raw input, while the UniProt branch already strips the version. efetch resolves
both forms, so the version must be optional.
"""

import os
import sys
import unittest

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

from resolve_gene_input import detect_input_type  # noqa: E402


class TestDetectInputType(unittest.TestCase):
    def test_refseq_with_and_without_version(self):
        for acc in ("XP_031329057", "XP_031329057.1", "NP_001005484",
                    "NP_001005484.1", "WP_010895915", "YP_009724390.1"):
            kind, value = detect_input_type(acc)
            self.assertEqual(kind, "ncbi", msg=f"{acc} should be ncbi, got {kind}")
            self.assertEqual(value, acc, msg="raw accession (incl. any version) is preserved for efetch")

    def test_genbank_with_and_without_version(self):
        for acc in ("AAB12345", "AAB12345.1", "KAF1234567", "KAF1234567.1"):
            kind, _ = detect_input_type(acc)
            self.assertEqual(kind, "ncbi", msg=f"{acc} should be ncbi, got {kind}")

    def test_uniprot_still_detected(self):
        for acc in ("P01501", "P01501.2", "Q16553", "A0A6M3Z554"):
            kind, _ = detect_input_type(acc)
            self.assertEqual(kind, "uniprot", msg=f"{acc} should be uniprot, got {kind}")

    def test_symbols_fall_through(self):
        for sym in ("melittin", "LY6E", "TP53", "luciferase"):
            kind, _ = detect_input_type(sym)
            self.assertEqual(kind, "symbol", msg=f"{sym} should be symbol, got {kind}")

    def test_existing_file_wins(self):
        kind, value = detect_input_type(__file__)
        self.assertEqual(kind, "file")
        self.assertEqual(value, __file__)


if __name__ == "__main__":
    unittest.main()
