#!/usr/bin/env python3
"""Canary test for the melettin_gt_v3 regression fixture.

Asserts the fixture under `tests/ground_truth_test/melettin_gt_v3/` is
well-formed and has the expected shape. This does NOT run the pipeline —
it only checks that the pinned reference artifacts are intact.

See `tests/ground_truth_test/melettin_gt_v3/README.md` for the full
regeneration procedure.
"""

import csv
import os
import re
import unittest


FIXTURE_DIR = os.path.join(
    os.path.dirname(__file__), "ground_truth_test", "melettin_gt_v3"
)
REGIONS_DIR = os.path.join(FIXTURE_DIR, "regions")
TREE_PATH = os.path.join(FIXTURE_DIR, "locus_1_tree.nwk")
FLANKING_PATH = os.path.join(FIXTURE_DIR, "flanking_parents.tsv")

EXPECTED_SPECIES = {
    "Colletes_gigas",
    "Euglossa_dilemma",
    "Melipona_beecheii",
    "Tetragonula_carbonaria",
    "Xylocopa_violacea",
}

EXPECTED_TREE_LEAVES = {
    "GOI_Melt",
    "GOI_Melt|Tetragonula_carbonaria_fa_b0_l1_fallback",
    "GOI_Melt|Euglossa_dilemma_fa_b0_l2_fallback",
    "GOI_Melt|Xylocopa_violacea_fa_b0_l1_fallback",
    "GOI_Melt|Colletes_gigas_fa_b0_l2_fallback",
}


def _newick_leaf_labels(nwk_text):
    """Extract leaf labels from a Newick string (no external deps)."""
    stripped = re.sub(r":\d+(\.\d+)?([eE][+-]?\d+)?", "", nwk_text)
    stripped = stripped.rstrip(";").strip()
    tokens = re.split(r"[(),]", stripped)
    leaves = []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        if re.fullmatch(r"\d+(\.\d+)?", tok):
            continue
        leaves.append(tok)
    return set(leaves)


class TestMelettinGtV3FixtureIntegrity(unittest.TestCase):
    def test_fixture_directory_exists(self):
        self.assertTrue(os.path.isdir(FIXTURE_DIR),
                        f"Fixture directory missing: {FIXTURE_DIR}")
        self.assertTrue(os.path.isdir(REGIONS_DIR),
                        f"regions/ missing: {REGIONS_DIR}")
        self.assertTrue(os.path.exists(TREE_PATH),
                        f"Tree file missing: {TREE_PATH}")
        self.assertTrue(os.path.exists(FLANKING_PATH),
                        f"Flanking TSV missing: {FLANKING_PATH}")

    def test_bed_files_present_for_all_species(self):
        bed_species = {
            os.path.basename(p).replace(".fa.regions.bed", "")
            for p in os.listdir(REGIONS_DIR)
            if p.endswith(".regions.bed")
        }
        self.assertEqual(bed_species, EXPECTED_SPECIES)

    def test_bed_files_have_one_row_six_columns_and_valid_coords(self):
        for fname in sorted(os.listdir(REGIONS_DIR)):
            if not fname.endswith(".regions.bed"):
                continue
            path = os.path.join(REGIONS_DIR, fname)
            with open(path) as fh:
                rows = [l for l in fh.read().splitlines() if l.strip()]
            self.assertEqual(len(rows), 1,
                             msg=f"{fname}: expected exactly 1 BED row, got {len(rows)}")
            cols = rows[0].split("\t")
            self.assertEqual(len(cols), 6,
                             msg=f"{fname}: expected 6 BED columns, got {len(cols)}: {cols}")
            chrom, start, end, name, score, strand = cols
            self.assertTrue(chrom, msg=f"{fname}: empty chrom")
            self.assertTrue(start.isdigit() and end.isdigit(),
                            msg=f"{fname}: non-integer coords start={start} end={end}")
            self.assertLess(int(start), int(end),
                            msg=f"{fname}: start >= end ({start} >= {end})")
            # Name is either bare (e.g. "Reg1_G7_...") or species-prefixed
            # ("Colletes_gigas|Reg1_G7_..."). The fixture was captured before
            # species prefixing landed; regenerated fixtures will include it.
            self.assertRegex(name, r"(^|\|)Reg\d+_",
                             msg=f"{fname}: unexpected region name: {name}")
            self.assertIn(strand, {"+", "-"},
                          msg=f"{fname}: bad strand {strand}")
            self.assertGreater(float(score), 0.0,
                               msg=f"{fname}: score must be > 0, got {score}")

    def test_score_tsv_files_present_with_one_data_row(self):
        REQUIRED_COLS = {
            "region_rank", "region_name", "chrom", "start", "end", "strand",
            "score", "quality_score", "coverage_score", "unique_genes",
            "confidence", "selection_reason",
        }
        seen_species = set()
        for fname in sorted(os.listdir(REGIONS_DIR)):
            if not fname.endswith(".scores.tsv"):
                continue
            seen_species.add(fname.replace(".fa.scores.tsv", ""))
            path = os.path.join(REGIONS_DIR, fname)
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                rows = list(reader)
                header = set(reader.fieldnames or [])
            self.assertTrue(REQUIRED_COLS.issubset(header),
                            msg=f"{fname}: missing required columns {REQUIRED_COLS - header}")
            self.assertEqual(len(rows), 1,
                             msg=f"{fname}: expected exactly 1 data row, got {len(rows)}")
            row = rows[0]
            self.assertIn(row["confidence"], {"HIGH", "MEDIUM", "LOW"},
                          msg=f"{fname}: unexpected confidence {row['confidence']}")
            self.assertGreater(float(row["score"]), 0.0)
            self.assertGreaterEqual(int(row["unique_genes"]), 1)
        self.assertEqual(seen_species, EXPECTED_SPECIES)

    def test_tree_has_expected_leaf_set(self):
        with open(TREE_PATH) as fh:
            nwk = fh.read().strip()
        self.assertTrue(nwk.endswith(";"), "Tree must end with ';'")
        leaves = _newick_leaf_labels(nwk)
        self.assertEqual(
            leaves, EXPECTED_TREE_LEAVES,
            msg=(
                f"Leaf set changed. A 5th leaf (e.g. Melipona_beecheii) "
                f"appearing here means tree-building behavior shifted — "
                f"regenerate fixture intentionally, do not silently update."
            ),
        )

    def test_flanking_parents_tsv_has_ten_rows_and_expected_columns(self):
        REQUIRED_COLS = {"name", "chrom", "start", "end", "strand",
                         "position", "distance_to_melittin"}
        with open(FLANKING_PATH) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
            header = set(reader.fieldnames or [])
        self.assertTrue(REQUIRED_COLS.issubset(header),
                        msg=f"flanking_parents.tsv: missing {REQUIRED_COLS - header}")
        self.assertEqual(len(rows), 10,
                         msg=f"flanking_parents.tsv: expected 10 rows, got {len(rows)}")
        positions = [r["position"] for r in rows]
        self.assertEqual(positions.count("upstream"), 5)
        self.assertEqual(positions.count("downstream"), 5)


if __name__ == "__main__":
    unittest.main()
