#!/usr/bin/env python3

import json
import os
import sys
import tempfile
import time
import unittest
from unittest import mock
from pathlib import Path


sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import build_llm_context as blc


class TestFetchSpeciesInfo(unittest.TestCase):
    def test_fetch_species_info_handles_integer_lineage_and_camel_case_fields(self):
        tax_response = {
            "taxonomy_nodes": [
                {
                    "taxonomy": {
                        "taxId": 7460,
                        "rank": "SPECIES",
                        "currentScientificName": {"name": "Apis mellifera"},
                        "classification": {
                            "domain": {"name": "Eukaryota", "id": 2759},
                            "kingdom": {"name": "Metazoa", "id": 33208},
                            "phylum": {"name": "Arthropoda", "id": 6656},
                            "class": {"name": "Insecta", "id": 50557},
                            "order": {"name": "Hymenoptera", "id": 7399},
                            "family": {"name": "Apidae", "id": 7458},
                            "genus": {"name": "Apis", "id": 7459},
                            "species": {"name": "Apis mellifera", "id": 7460},
                        },
                        # This is the shape that crashed the Apr 12 run.
                        "lineage": [1, 131567, 2759, 33208, 6656, 50557, 7399, 7458, 7459],
                    }
                }
            ]
        }
        genome_response = {
            "reports": [
                {
                    "assemblyStats": {
                        "totalSequenceLength": "225000000",
                        "scaffoldN50": 1250000,
                    },
                    "annotationInfo": {
                        "stats": {
                            "geneCounts": {
                                "total": 12345,
                            }
                        }
                    },
                }
            ]
        }

        def fake_fetch(url, timeout=20):
            if url.startswith(blc.NCBI_TAXONOMY_BASE):
                return tax_response
            if url.startswith(blc.NCBI_DATASETS_BASE):
                return genome_response
            self.fail(f"Unexpected URL fetched: {url}")

        with mock.patch.object(blc, "_read_cache", return_value=None), mock.patch.object(
            blc, "_write_cache"
        ), mock.patch.object(blc, "_fetch_json", side_effect=fake_fetch):
            info = blc.fetch_species_info("Apis mellifera")

        self.assertEqual(info["kingdom"], "Animalia")
        self.assertEqual(info["phylum"], "Arthropoda")
        self.assertEqual(info["genome_size_mb"], 225.0)
        self.assertEqual(info["scaffold_n50"], 1250000)
        self.assertEqual(info["gene_count"], 12345)
        self.assertTrue(info["lineage"])
        self.assertTrue(all(isinstance(entry, dict) for entry in info["lineage"]))
        self.assertEqual(info["lineage"][0]["rank"], "superkingdom")
        self.assertEqual(info["lineage"][0]["name"], "Eukaryota")

    def test_fetch_species_info_does_not_cache_degraded_fallback_results(self):
        with mock.patch.object(blc, "_read_cache", return_value=None), mock.patch.object(
            blc, "_write_cache"
        ) as write_cache, mock.patch.object(blc, "_fetch_json", return_value=None):
            info = blc.fetch_species_info("Unknown species")

        self.assertEqual(info["kingdom"], "Unknown")
        self.assertEqual(info["genome_size_mb"], 0)
        self.assertEqual(info["gene_count"], 0)
        write_cache.assert_not_called()


class TestSpeciesCache(unittest.TestCase):
    def test_read_cache_ignores_degraded_entries(self):
        degraded = {
            "name": "Apis mellifera",
            "kingdom": "Unknown",
            "phylum": "Unknown",
            "genome_size_mb": 0,
            "gene_count": 0,
            "scaffold_n50": 0,
            "lineage": [{"name": "Apis mellifera", "rank": "species"}],
            "_cached_at": time.time(),
        }

        with tempfile.TemporaryDirectory() as tmpdir, mock.patch.object(
            blc, "CACHE_DIR", Path(tmpdir)
        ):
            cache_path = Path(tmpdir) / f"{blc._cache_key('Apis mellifera')}.json"
            cache_path.write_text(json.dumps(degraded))

            cached = blc._read_cache("Apis mellifera")

        self.assertIsNone(cached)


if __name__ == "__main__":
    unittest.main()
