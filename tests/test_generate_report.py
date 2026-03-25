#!/usr/bin/env python3

import json
import os
import shutil
import sys
import tempfile
import unittest


sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from generate_report import build_report


class TestGenerateReport(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.test_dir, "regions"), exist_ok=True)
        os.makedirs(os.path.join(self.test_dir, "hits"), exist_ok=True)
        os.makedirs(os.path.join(self.test_dir, "scores"), exist_ok=True)
        self.qc_json = os.path.join(self.test_dir, "qc.json")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _write(self, rel_path, content):
        path = os.path.join(self.test_dir, rel_path)
        with open(path, "w") as fh:
            fh.write(content)
        return path

    def test_build_report_summarizes_evidence_and_qc_propagation(self):
        with open(self.qc_json, "w") as fh:
            json.dump(
                [
                    {
                        "genome": "GCF_000001.1.fna",
                        "status": "PASS",
                        "msg": "OK",
                        "thresholds": {"min_n50": 5000, "max_contigs": 500000},
                    },
                    {
                        "genome": "GCF_000002.1.fna",
                        "status": "FAIL",
                        "msg": "too many contigs",
                    },
                ],
                fh,
            )

        self._write(
            "regions/GCF_000001.1.fna.gff",
            "\n".join(
                [
                    "##gff-version 3",
                    "chr1\texon_annotation\tmRNA\t100\t300\t92.0\t+\t.\tID=GOI_gene|pass;Name=GOI_gene;SynVoy_Parent=GOI_gene;SynVoyRole=goi;Confidence=HIGH;GOIClass=confident_goi;EvidenceType=exon_annotation;Identity=92.0;Exons=3",
                    "chr1\tflanking_annotation\tmRNA\t400\t600\t81.0\t+\t.\tID=geneA|pass;Name=geneA;SynVoy_Parent=geneA;SynVoyRole=flanking;Confidence=HIGH;EvidenceType=flanking_miniprot;Identity=81.0;Exons=2",
                ]
            ) + "\n",
        )
        self._write(
            "regions/GCF_000002.1.fna.gff",
            "\n".join(
                [
                    "##gff-version 3",
                    "chr7\tfallback_hits\tmRNA\t1000\t1200\t55.0\t-\t.\tID=GOI_gene|ambig;Name=GOI_gene;SynVoy_Parent=GOI_gene;SynVoyRole=goi;Confidence=LOW;GOIClass=ambiguous_goi_family_member;EvidenceType=fallback_hit_span;Identity=55.0;Exons=1",
                ]
            ) + "\n",
        )

        self._write("regions/GCF_000001.1.fna.faa", ">a\nMPEP\n>b\nMPEP\n")
        self._write("regions/GCF_000002.1.fna.faa", ">c\nMPEP\n")
        self._write(
            "hits/synteny_block_locus_1_GCF_000001.1.fna.m8",
            "q1\tchr1\t90\t100\t0\t0\t1\t100\t10\t110\t1e-40\t200\n",
        )
        self._write(
            "hits/synteny_block_locus_1_GCF_000002.1.fna.m8",
            "q2\tchr7\t55\t60\t0\t0\t1\t60\t1000\t1180\t1e-6\t70\n",
        )

        self._write(
            "scores/GCF_000001.1.fna.scores.tsv",
            "\n".join(
                [
                    "region_rank\tregion_name\tchrom\tstart\tend\tstrand\tscore\tquality_score\tcoverage_score\tunique_genes\ttotal_genes_expected\tconsistency\tstrand_consistency\tp_value\tgoi_overlap\tis_goi_anchor\tconfidence\tselection_reason",
                    "1\tReg1\tchr1\t0\t1000\t+\t0.92\t0.90\t0.94\t5\t6\t0.8\t1.0\t0.01\ttrue\ttrue\tHIGH\tbest_score",
                ]
            ) + "\n",
        )
        self._write(
            "scores/GCF_000002.1.fna.scores.tsv",
            "\n".join(
                [
                    "region_rank\tregion_name\tchrom\tstart\tend\tstrand\tscore\tquality_score\tcoverage_score\tunique_genes\ttotal_genes_expected\tconsistency\tstrand_consistency\tp_value\tgoi_overlap\tis_goi_anchor\tconfidence\tselection_reason",
                    "1\tReg1\tchr7\t900\t1400\t-\t0.31\t0.25\t0.36\t1\t6\t0.2\t1.0\t0.20\tfalse\tfalse\tLOW\tfallback_only",
                ]
            ) + "\n",
        )

        report = build_report(self.test_dir, qc_json=self.qc_json, qc_policy="drop")

        self.assertEqual(report["qc_summary"]["pass"], 1)
        self.assertEqual(report["qc_summary"]["fail"], 1)
        self.assertEqual(report["qc_summary"]["qc_fail_policy"], "drop")
        self.assertEqual(report["summary"]["total_goi_annotations"], 2)
        self.assertEqual(report["summary"]["ambiguous_goi_annotations"], 1)
        self.assertEqual(report["summary"]["fallback_goi_annotations"], 1)
        self.assertEqual(report["summary"]["low_confidence_regions"], 1)
        self.assertIn("GCF_000002.1", report["summary"]["goi_ambiguous_only_genomes"])
        self.assertIn("GCF_000002.1", report["summary"]["failed_qc_genomes_with_downstream_results"])
        self.assertEqual(report["regions"]["confidence_counts"]["HIGH"], 1)
        self.assertEqual(report["regions"]["confidence_counts"]["LOW"], 1)
        self.assertEqual(report["annotations"]["goi_class_counts"]["confident_goi"], 1)
        self.assertEqual(report["annotations"]["goi_class_counts"]["ambiguous_goi_family_member"], 1)


if __name__ == "__main__":
    unittest.main()
