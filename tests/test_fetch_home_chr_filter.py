#!/usr/bin/env python3
"""Regression tests for filter_chromosomes_only() (fetch_home_genome.py).

Bug (found 2026-07-02 on the oskar easy-mode run): the filter kept only
NC_/CM_-prefixed sequences and discarded everything else as "unplaced
scaffolds". In the Drosophila melanogaster RefSeq assembly the main chromosome
arms 2L/2R/3L/3R are NT_ records (20-32 Mb), so they were silently stripped —
taking oskar (osk, on 3R = NT_033777.3) with them. LOCATE_GENE then found no
hits and the pipeline aborted at SPLIT_LOCI.

The naive fix (add NT_ to the keep-list) regresses human GRCh38, where NT_
records are sub-Mb *unlocalized scaffolds* that should stay dropped. So the fix
keeps NC_/CM_ by prefix, always drops NW_/NZ_, and additionally keeps any large
(>= min_chrom_len) sequence — length disambiguates the fly arm from the human
scaffold. These tests exercise the pure file-rewrite logic (no network).
"""
import os
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))
sys.path.insert(0, BIN_DIR)

from fetch_home_genome import filter_chromosomes_only  # noqa: E402


def _write_fasta(path: Path, records):
    """records: list of (seq_id, length_bp). Writes a FASTA of that shape."""
    with open(path, "w") as fh:
        for seq_id, length in records:
            fh.write(f">{seq_id} some description here\n")
            fh.write("A" * length + "\n")


def _seq_ids(path: Path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


class TestChrFilter(unittest.TestCase):
    # small threshold so we don't allocate megabytes in a unit test
    THRESH = 1000

    def _run(self, records):
        with TemporaryDirectory() as d:
            p = Path(d) / "home_genome.fna"
            _write_fasta(p, records)
            filter_chromosomes_only(p, min_chrom_len=self.THRESH)
            return set(_seq_ids(p))

    def test_fly_keeps_large_NT_arms(self):
        """Drosophila: NC_ molecules + large NT_ arm kept; small NW_ dropped."""
        kept = self._run([
            ("NC_004354.4", 800),    # chr X (small but NC_ -> kept by prefix)
            ("NT_033777.3", 2000),   # chr 3R arm (large NT_ -> kept by length; carries oskar)
            ("NW_001844851.1", 500), # unplaced scaffold -> dropped
        ])
        self.assertIn("NT_033777.3", kept)   # the bug: this used to be dropped
        self.assertIn("NC_004354.4", kept)
        self.assertNotIn("NW_001844851.1", kept)

    def test_human_drops_small_NT_scaffolds(self):
        """Human GRCh38: NC_ chromosome kept; small NT_/NW_ scaffolds dropped
        (identical to the old behaviour -> no benchmark regression)."""
        kept = self._run([
            ("NC_000012.12", 2000),   # chr 12 -> kept
            ("NT_187361.1", 500),     # unlocalized scaffold (< threshold) -> dropped
            ("NW_009646201.1", 800),  # alt locus -> dropped
        ])
        self.assertEqual(kept, {"NC_000012.12"})

    def test_large_NZ_scaffold_still_dropped(self):
        """NZ_ is always a scaffold prefix: dropped even when large."""
        kept = self._run([
            ("NC_000001.11", 2000),
            ("NZ_ABCD01000001.1", 5000),  # large WGS scaffold -> still dropped
        ])
        self.assertEqual(kept, {"NC_000001.11"})

    def test_pure_scaffold_assembly_untouched(self):
        """No chromosome-scale seqs -> keep everything (never blank the file)."""
        recs = [("NW_1", 300), ("NW_2", 400)]
        kept = self._run(recs)
        self.assertEqual(kept, {"NW_1", "NW_2"})

    def test_all_chromosomes_no_removal(self):
        """Only NC_ chromosomes present -> file unchanged."""
        kept = self._run([("NC_1", 800), ("NC_2", 900)])
        self.assertEqual(kept, {"NC_1", "NC_2"})


if __name__ == "__main__":
    unittest.main()
