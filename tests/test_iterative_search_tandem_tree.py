#!/usr/bin/env python3
"""Regression test for the tandem-copy GOI tree-input fix.

Background
----------
`iterative_search_runner.py` historically wrote a single `expanded_db.faa`
that served two purposes:
  (a) seed the next wave's MMseqs2 query DB,
  (b) feed MAFFT/IQ-TREE in COMPUTE_TREE.

To prevent contamination of (a), the seeding filter only kept GOI hits with
goi_class in {confident_goi, probable_goi}. Tandem-copy hits (goi_class =
tandem_goi_copy) were silently dropped from BOTH purposes, which meant
species like Apis cerana — whose only Melt evidence is a 100%-identity
tandem duplicate — never appeared in the tree.

The fix splits the two concerns:
  - `_classify_goi_for_seed_and_tree` returns (seed, tree_extra, suppressed_n)
  - the seed list is unchanged
  - tree_extra picks up MEDIUM/HIGH GOI hits that didn't qualify for seeding
  - main() writes goi_for_tree.faa = expanded_db.faa + tree_extra records

This test pins both the helper's behaviour and the suppression accounting.
"""

import os
import sys
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

from iterative_search_runner import _classify_goi_for_seed_and_tree  # noqa: E402


def _gene(gid, seq="MKT"):
    return {"id": gid, "seq": seq}


def _is_goi_query_id(gid):
    return gid.startswith("GOI_")


class TestSeedTreeClassifier(unittest.TestCase):

    def test_tandem_copy_goes_to_tree_not_seed(self):
        """Apis cerana scenario: only tandem_goi_copy GOI evidence."""
        all_genes = [
            _gene("GOI_copy_2|Apis_cerana_fna_b0_l1"),
            _gene("GOI_copy_1|Apis_cerana_fna_b0_l1"),
        ]
        meta = {
            "GOI_copy_2|Apis_cerana_fna_b0_l1": {
                "role": "goi", "confidence": "MEDIUM",
                "goi_class": "tandem_goi_copy",
            },
            "GOI_copy_1|Apis_cerana_fna_b0_l1": {
                "role": "goi", "confidence": "MEDIUM",
                "goi_class": "tandem_goi_copy",
            },
        }
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            all_genes, meta, _is_goi_query_id
        )
        self.assertEqual(seed, [],
                         msg="tandem_goi_copy must not seed the next wave")
        self.assertEqual(len(tree_extra), 2,
                         msg="MEDIUM tandem copies must reach the tree FASTA")
        self.assertEqual(suppressed, 2)

    def test_confident_goi_seeds_normally(self):
        all_genes = [_gene("GOI_Melt|Apis_florea_fna_b0_l1_exon_ann")]
        meta = {
            "GOI_Melt|Apis_florea_fna_b0_l1_exon_ann": {
                "role": "goi", "confidence": "HIGH",
                "goi_class": "confident_goi",
            },
        }
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            all_genes, meta, _is_goi_query_id
        )
        self.assertEqual(len(seed), 1)
        self.assertEqual(tree_extra, [])
        self.assertEqual(suppressed, 0)

    def test_low_confidence_dropped_from_both(self):
        """LOW-confidence hits feed neither seed nor tree (too noisy)."""
        all_genes = [_gene("GOI_copy_1|low_quality")]
        meta = {
            "GOI_copy_1|low_quality": {
                "role": "goi", "confidence": "LOW",
                "goi_class": "tandem_goi_copy",
            },
        }
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            all_genes, meta, _is_goi_query_id
        )
        self.assertEqual(seed, [])
        self.assertEqual(tree_extra, [])
        self.assertEqual(suppressed, 1,
                         msg="LOW-confidence GOI still counts as suppressed")

    def test_flanking_role_never_seeds_or_trees(self):
        all_genes = [_gene("gene-LOC102655466|flank_ann")]
        meta = {
            "gene-LOC102655466|flank_ann": {
                "role": "flanking", "confidence": "HIGH",
                "goi_class": "",
            },
        }
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            all_genes, meta, _is_goi_query_id
        )
        self.assertEqual(seed, [])
        self.assertEqual(tree_extra, [])
        self.assertEqual(suppressed, 0,
                         msg="non-GOI roles are not 'suppressed', they're skipped")

    def test_mixed_real_world_apis_cerana_case(self):
        """Full Apis cerana TSV row pattern: 2 tandem GOI copies + 9 flanking."""
        all_genes = (
            [_gene(f"GOI_copy_{i}|Apis_cerana_fna_b0_l1") for i in (1, 2)]
            + [_gene(f"gene-LOC{j}|Apis_cerana_fna_b0_fl1_flank_ann")
               for j in (102655466, 409662, 412898, 726866, 409659, 726827,
                        100578368, 726817, 726855)]
        )
        meta = {}
        for i in (1, 2):
            meta[f"GOI_copy_{i}|Apis_cerana_fna_b0_l1"] = {
                "role": "goi", "confidence": "MEDIUM",
                "goi_class": "tandem_goi_copy",
            }
        for j in (102655466, 409662, 412898, 726866, 409659, 726827,
                  100578368, 726817, 726855):
            meta[f"gene-LOC{j}|Apis_cerana_fna_b0_fl1_flank_ann"] = {
                "role": "flanking", "confidence": "HIGH", "goi_class": "",
            }
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            all_genes, meta, _is_goi_query_id
        )
        self.assertEqual(len(seed), 0)
        self.assertEqual(len(tree_extra), 2,
                         msg="Apis cerana must contribute 2 tree-only GOI hits")
        self.assertEqual(suppressed, 2)

    def test_unknown_meta_falls_back_to_id_heuristic(self):
        """If feature_meta is missing for a gene, role is inferred from the id."""
        all_genes = [_gene("GOI_Melt|some_evidence")]
        # No metadata entry — function uses is_goi_query_id_fn fallback,
        # but without confidence/goi_class it can't seed.
        seed, tree_extra, suppressed = _classify_goi_for_seed_and_tree(
            all_genes, feature_meta={}, is_goi_query_id_fn=_is_goi_query_id
        )
        self.assertEqual(seed, [])
        self.assertEqual(tree_extra, [],
                         msg="empty confidence is not in {HIGH, MEDIUM} so no tree-extra")
        self.assertEqual(suppressed, 1)


if __name__ == "__main__":
    unittest.main()
