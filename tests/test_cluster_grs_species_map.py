#!/usr/bin/env python3
"""Tests for cluster_grs.py species-map BED prefixing (plan item 3c)."""

import csv
import os
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")
sys.path.insert(0, BIN_DIR)

from cluster_grs import resolve_species, _strip_fasta_exts  # noqa: E402


CLUSTER_SCRIPT = os.path.join(BIN_DIR, "cluster_grs.py")


class TestResolveSpeciesHelper(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.map_path = os.path.join(self.tmp, "species_mapping.tsv")

    def _write_map(self, rows):
        with open(self.map_path, "w") as fh:
            for row in rows:
                fh.write("\t".join(row) + "\n")

    def test_returns_none_for_missing_file(self):
        self.assertIsNone(resolve_species(None, "Colletes_gigas"))
        self.assertIsNone(resolve_species("/does/not/exist", "Colletes_gigas"))

    def test_returns_none_for_sentinel_filename(self):
        sentinel = os.path.join(self.tmp, "NO_SPECIES_MAP")
        open(sentinel, "w").close()
        self.assertIsNone(resolve_species(sentinel, "Colletes_gigas"))

    def test_exact_match_pro_mode(self):
        self._write_map([
            ("Colletes_gigas", "Colletes gigas", "genome"),
            ("Euglossa_dilemma", "Euglossa dilemma", "genome"),
        ])
        self.assertEqual(
            resolve_species(self.map_path, "Colletes_gigas"),
            "Colletes_gigas",  # spaces normalized to underscores
        )

    def test_fasta_extension_stripped_for_lookup(self):
        self._write_map([("Colletes_gigas", "Colletes gigas", "genome")])
        self.assertEqual(
            resolve_species(self.map_path, "Colletes_gigas.fa"),
            "Colletes_gigas",
        )
        self.assertEqual(
            resolve_species(self.map_path, "Colletes_gigas.fna.gz"),
            "Colletes_gigas",
        )

    def test_version_suffix_fallback(self):
        self._write_map([("GCF_003672135", "Apis mellifera", "species")])
        self.assertEqual(
            resolve_species(self.map_path, "GCF_003672135.1.fna"),
            "Apis_mellifera",
        )

    def test_versioned_accession_exact_preferred_over_bare(self):
        self._write_map([
            ("GCF_003672135", "Wrong match", "species"),
            ("GCF_003672135.1", "Apis mellifera", "species"),
        ])
        # With .fna extension, first candidate (with version) should match
        self.assertEqual(
            resolve_species(self.map_path, "GCF_003672135.1.fna"),
            "Apis_mellifera",
        )

    def test_multi_dot_assembly_name_resolves(self):
        # Regression: AnoCar2.0 / GRCh38.p14 / Amel_HAv3.1 etc. embed dots in
        # the assembly name. The cascade must keep stripping until it reaches
        # the bare accession that STAGE_GENOMES truncates to (%%.* in bash).
        self._write_map([("GCF_000090745", "Anolis carolinensis", "genome")])
        self.assertEqual(
            resolve_species(self.map_path, "GCF_000090745.1_AnoCar2.0_genomic.fna.gz"),
            "Anolis_carolinensis",
        )

    def test_deeply_nested_dots_still_resolves(self):
        # Even more pathological: GRCh38.p14_assembly.fna.gz
        self._write_map([("GCF_000001405", "Homo sapiens", "genome")])
        self.assertEqual(
            resolve_species(self.map_path, "GCF_000001405.40_GRCh38.p14_assembly.fna.gz"),
            "Homo_sapiens",
        )

    def test_unknown_key_returns_none(self):
        self._write_map([("Colletes_gigas", "Colletes gigas", "genome")])
        self.assertIsNone(resolve_species(self.map_path, "not_in_map"))

    def test_rejects_fasta_header_leak_as_species_value(self):
        # Regression: STAGE_GENOMES used to write '>WUUM01000001.1' into the
        # species column when the FASTA header had no binomial. Defensive
        # sanitization should treat that as 'no mapping' so the BED name
        # stays bare instead of carrying '>WUUM01000001.1|Reg1_...'.
        self._write_map([
            ("Colletes_gigas", ">WUUM01000001.1", "genome"),
            ("Other_sp", "scaffold-1", "genome"),
        ])
        self.assertIsNone(resolve_species(self.map_path, "Colletes_gigas"))
        self.assertIsNone(resolve_species(self.map_path, "Other_sp"))

    def test_rejects_pipe_or_tab_in_species_value(self):
        self._write_map([("X_y", "Bad|name", "genome")])
        self.assertIsNone(resolve_species(self.map_path, "X_y"))

    def test_accepts_normal_binomial_with_underscore_or_spaces(self):
        self._write_map([
            ("a", "Apis mellifera", "genome"),
            ("b", "Apis_mellifera", "genome"),
        ])
        self.assertEqual(resolve_species(self.map_path, "a"), "Apis_mellifera")
        self.assertEqual(resolve_species(self.map_path, "b"), "Apis_mellifera")

    def test_strip_fasta_exts(self):
        self.assertEqual(_strip_fasta_exts("a.fa"), "a")
        self.assertEqual(_strip_fasta_exts("a.fna"), "a")
        self.assertEqual(_strip_fasta_exts("a.fasta"), "a")
        self.assertEqual(_strip_fasta_exts("a.fna.gz"), "a")
        self.assertEqual(_strip_fasta_exts("a.fa.gz"), "a")
        self.assertEqual(_strip_fasta_exts("a"), "a")


class TestClusterGrsBedPrefixing(unittest.TestCase):
    """End-to-end: run cluster_grs.py with/without --species_map and check output."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

        # Minimal synteny BED: 3 flanking genes on chr1
        self.synteny_bed = os.path.join(self.tmp, "synteny.bed")
        with open(self.synteny_bed, "w") as fh:
            fh.write("chr1\t100\t200\tgeneA\t.\t+\n")
            fh.write("chr1\t300\t400\tgeneB\t.\t+\n")
            fh.write("chr1\t500\t600\tgeneC\t.\t+\n")

        # Minimal hits m8: 3 hits on target chr7 clustered together
        self.hits = os.path.join(self.tmp, "hits.m8")
        with open(self.hits, "w") as fh:
            fh.write("geneA\tchr7\t95.0\t100\t0\t0\t1\t100\t1000\t1100\t1e-40\t200\n")
            fh.write("geneB\tchr7\t92.0\t100\t0\t0\t1\t100\t1200\t1300\t1e-35\t190\n")
            fh.write("geneC\tchr7\t90.0\t100\t0\t0\t1\t100\t1400\t1500\t1e-30\t180\n")

        # Tiny fake target genome (so --genome doesn't crash)
        self.genome = os.path.join(self.tmp, "Colletes_gigas.fa")
        with open(self.genome, "w") as fh:
            fh.write(">chr7\n" + "A" * 2000 + "\n")

        # Species map
        self.species_map = os.path.join(self.tmp, "species_mapping.tsv")
        with open(self.species_map, "w") as fh:
            fh.write("Colletes_gigas\tColletes gigas\tgenome\n")

        self.bed_out = os.path.join(self.tmp, "Colletes_gigas.regions.bed")
        self.scores_out = os.path.join(self.tmp, "Colletes_gigas.scores.tsv")

    def _run(self, extra_args):
        cmd = [
            sys.executable, CLUSTER_SCRIPT,
            "--hits", self.hits,
            "--synteny_bed", self.synteny_bed,
            "--genome", self.genome,
            "--output", self.bed_out,
            "--scores_output", self.scores_out,
            "--flanking_count", "3",
            "--cluster_distance", "1000",
            "--min_score", "0.3",
        ] + list(extra_args)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0,
                         msg=f"cluster_grs.py failed:\nstdout={proc.stdout}\nstderr={proc.stderr}")
        return proc

    def _read_bed_row(self):
        with open(self.bed_out) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    return line.split("\t")
        return None

    def _read_scores_row(self):
        with open(self.scores_out) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
            return rows[0] if rows else None, (reader.fieldnames or [])

    def test_without_species_map_bed_name_is_bare(self):
        self._run([])
        cols = self._read_bed_row()
        self.assertIsNotNone(cols)
        name = cols[3]
        self.assertTrue(name.startswith("Reg"),
                        msg=f"Expected bare 'Reg...' name, got: {name}")
        self.assertNotIn("|", name,
                         msg=f"Expected no species prefix, got: {name}")

        scores_row, fieldnames = self._read_scores_row()
        self.assertIn("species", fieldnames)
        self.assertEqual(scores_row["species"], "")

    def test_with_species_map_bed_name_has_prefix(self):
        self._run(["--species_map", self.species_map])
        cols = self._read_bed_row()
        self.assertIsNotNone(cols)
        name = cols[3]
        self.assertTrue(name.startswith("Colletes_gigas|Reg"),
                        msg=f"Expected 'Colletes_gigas|Reg...' name, got: {name}")

        scores_row, fieldnames = self._read_scores_row()
        self.assertIn("species", fieldnames)
        self.assertEqual(scores_row["species"], "Colletes_gigas")

    def test_explicit_genome_name_overrides_basename(self):
        # Use a genome file named differently, but pass --genome_name to force
        # the lookup key.
        alt_genome = os.path.join(self.tmp, "some_download.fna")
        with open(alt_genome, "w") as fh:
            fh.write(">chr7\n" + "A" * 2000 + "\n")
        cmd = [
            sys.executable, CLUSTER_SCRIPT,
            "--hits", self.hits,
            "--synteny_bed", self.synteny_bed,
            "--genome", alt_genome,
            "--genome_name", "Colletes_gigas",
            "--species_map", self.species_map,
            "--output", self.bed_out,
            "--scores_output", self.scores_out,
            "--flanking_count", "3",
            "--cluster_distance", "1000",
            "--min_score", "0.3",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        cols = self._read_bed_row()
        self.assertTrue(cols[3].startswith("Colletes_gigas|Reg"))

    def test_unknown_genome_falls_back_to_bare_name(self):
        # Species map doesn't contain this genome → no prefix, no error
        cmd = [
            sys.executable, CLUSTER_SCRIPT,
            "--hits", self.hits,
            "--synteny_bed", self.synteny_bed,
            "--genome", self.genome,
            "--genome_name", "Unknown_sp",
            "--species_map", self.species_map,
            "--output", self.bed_out,
            "--scores_output", self.scores_out,
            "--flanking_count", "3",
            "--cluster_distance", "1000",
            "--min_score", "0.3",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        cols = self._read_bed_row()
        self.assertFalse(cols[3].startswith("Unknown_sp|"))
        self.assertTrue(cols[3].startswith("Reg"))


if __name__ == "__main__":
    unittest.main()
