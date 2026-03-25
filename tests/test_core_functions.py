#!/usr/bin/env python3
"""
Basic test suite for SynVoy pipeline core functions.
Tests critical functions for correctness.
"""

import unittest
import sys
import os
import tempfile
import shutil

# Add bin directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from iterative_search_runner import (
    normalize_coordinates,
    extract_base_gene_id,
    parse_hits,
    identify_best_synteny_block,
    merge_synteny_blocks
)
from sequence_utils import parse_fasta, write_fasta, reverse_complement, translate, extract_id


class TestCoordinateNormalization(unittest.TestCase):
    """Test coordinate handling functions."""
    
    def test_normalize_coordinates_forward(self):
        """Test that forward strand coordinates stay in order."""
        start, end = normalize_coordinates(100, 200)
        self.assertEqual(start, 100)
        self.assertEqual(end, 200)
    
    def test_normalize_coordinates_reverse(self):
        """Test that reverse strand coordinates get swapped."""
        start, end = normalize_coordinates(200, 100)
        self.assertEqual(start, 100)
        self.assertEqual(end, 200)
    
    def test_normalize_coordinates_equal(self):
        """Test edge case where start == end."""
        start, end = normalize_coordinates(100, 100)
        self.assertEqual(start, 100)
        self.assertEqual(end, 100)


class TestGeneIDExtraction(unittest.TestCase):
    """Test gene ID parsing and extraction."""
    
    def test_extract_base_simple(self):
        """Test simple gene ID without suffixes."""
        result = extract_base_gene_id("gene-LOC726866")
        self.assertEqual(result, "gene-LOC726866")
    
    def test_extract_base_with_exon(self):
        """Test gene ID with exon suffix."""
        result = extract_base_gene_id("gene-LOC726866|exon_1")
        self.assertEqual(result, "gene-LOC726866")
    
    def test_extract_base_with_variant(self):
        """Test gene ID with variant suffix."""
        result = extract_base_gene_id("gene-LOC726866|var1")
        self.assertEqual(result, "gene-LOC726866")
    
    def test_extract_base_with_genome(self):
        """Test gene ID with genome suffix."""
        result = extract_base_gene_id("gene-LOC726866|GCA_000001234_MP000001")
        self.assertEqual(result, "gene-LOC726866")
    
    def test_extract_goi_prefix(self):
        """Test GOI-prefixed gene ID."""
        result = extract_base_gene_id("GOI_gene-LOC726866")
        self.assertEqual(result, "GOI_gene-LOC726866")


class TestFastaIO(unittest.TestCase):
    """Test FASTA reading and writing."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.test_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.test_dir)
    
    def test_write_and_read_fasta(self):
        """Test round-trip FASTA writing and reading."""
        test_file = os.path.join(self.test_dir, "test.fasta")
        test_seqs = [
            ("seq1", "ATCGATCG"),
            ("seq2", "GCTAGCTA"),
            ("seq3_with_description", "AAAATTTTCCCCGGGG")
        ]
        
        write_fasta(test_seqs, test_file)
        self.assertTrue(os.path.exists(test_file))
        
        read_seqs = list(parse_fasta(test_file))
        self.assertEqual(len(read_seqs), 3)
        # parse_fasta returns (raw_header, clean_id, sequence)
        self.assertEqual(read_seqs[0][1], "seq1")  # clean_id
        self.assertEqual(read_seqs[0][2], "ATCGATCG")  # sequence
        self.assertEqual(read_seqs[1][1], "seq2")
        self.assertEqual(read_seqs[1][2], "GCTAGCTA")
    
    def test_fasta_with_line_wrapping(self):
        """Test FASTA with line-wrapped sequences."""
        test_file = os.path.join(self.test_dir, "wrapped.fasta")
        
        # Write FASTA with wrapped lines
        with open(test_file, 'w') as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCG\n")
            f.write("ATCGATCGATCG\n")  # Continuation
            f.write(">seq2\n")
            f.write("GCTA\n")
        
        seqs = list(parse_fasta(test_file))
        self.assertEqual(len(seqs), 2)
        # parse_fasta returns (raw_header, clean_id, sequence)
        self.assertEqual(seqs[0][2], "ATCGATCGATCGATCGATCGATCG")  # sequence
        self.assertEqual(seqs[1][2], "GCTA")

    def test_extract_id_preserves_long_loc_gene_ids(self):
        """Ensure long LOC-style IDs are not truncated by GenBank accession regex."""
        self.assertEqual(extract_id("LOC143834063"), "LOC143834063")
        self.assertEqual(extract_id("LOC143834063 some description"), "LOC143834063")

    def test_extract_id_keeps_genbank_accession_behavior(self):
        """Ensure canonical GenBank protein accessions are still extracted."""
        self.assertEqual(extract_id("AAA12345.1 hypothetical protein"), "AAA12345.1")


class TestSequenceOperations(unittest.TestCase):
    """Test sequence manipulation functions."""
    
    def test_reverse_complement(self):
        """Test reverse complement operation."""
        seq = "ATCGATCG"
        result = reverse_complement(seq)
        self.assertEqual(result, "CGATCGAT")
    
    def test_reverse_complement_lowercase(self):
        """Test reverse complement with lowercase input."""
        seq = "atcgatcg"
        result = reverse_complement(seq)
        self.assertEqual(result.upper(), "CGATCGAT")
    
    def test_translate_simple(self):
        """Test translation of simple codon."""
        seq = "ATG"  # Methionine
        result = translate(seq)
        self.assertEqual(result, "M")
    
    def test_translate_stop_codon(self):
        """Test translation with stop codon."""
        seq = "ATGTAA"  # M + STOP
        result = translate(seq)
        self.assertIn("M", result)
        # Stop codon should be translated to *
        self.assertIn("*", result)


class TestHitsFiltering(unittest.TestCase):
    """Test hit parsing and filtering."""
    
    def setUp(self):
        """Create temporary directory and test hits file."""
        self.test_dir = tempfile.mkdtemp()
        self.hits_file = os.path.join(self.test_dir, "test_hits.m8")
        
        # Write test hits (MMseqs2 format)
        # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
        with open(self.hits_file, 'w') as f:
            # Good hit
            f.write("query1\tchr1\t90.0\t100\t10\t0\t1\t100\t1000\t1100\t1e-50\t200\n")
            # Low identity (should be filtered)
            f.write("query2\tchr1\t30.0\t100\t70\t0\t1\t100\t2000\t2100\t1e-10\t50\n")
            # Short alignment (should be filtered)
            f.write("query3\tchr1\t90.0\t10\t1\t0\t1\t10\t3000\t3010\t1e-5\t20\n")
            # High evalue (should be filtered)
            f.write("query4\tchr1\t90.0\t100\t10\t0\t1\t100\t4000\t4100\t1.0\t10\n")
            # Another good hit
            f.write("query5\tchr2\t85.0\t150\t22\t0\t1\t150\t5000\t5150\t1e-60\t250\n")
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.test_dir)
    
    def test_parse_hits_filtering(self):
        """Test that hits are correctly filtered."""
        hits = parse_hits(self.hits_file, min_identity=40.0, min_length=50, evalue_thresh=1e-5)
        
        # Should have 2 good hits (query1 and query5)
        self.assertEqual(len(hits), 2)
        self.assertEqual(hits[0]['query'], 'query1')
        self.assertEqual(hits[1]['query'], 'query5')
    
    def test_parse_hits_strict_filtering(self):
        """Test stricter filtering parameters."""
        hits = parse_hits(self.hits_file, min_identity=87.0, min_length=100, evalue_thresh=1e-20)
        
        # Should have only 1 hit (query1)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]['query'], 'query1')
    
    def test_parse_nonexistent_file(self):
        """Test handling of nonexistent file."""
        hits = parse_hits("nonexistent.m8", min_identity=40.0, min_length=50, evalue_thresh=1e-5)
        self.assertEqual(len(hits), 0)


class TestSyntenyIdentification(unittest.TestCase):
    """Test synteny block identification."""
    
    def test_single_cluster(self):
        """Test identification with hits forming single cluster."""
        hits = [
            {'query': 'gene1', 'chrom': 'chr1', 'start': 1000, 'end': 1100},
            {'query': 'gene2', 'chrom': 'chr1', 'start': 1200, 'end': 1300},
            {'query': 'gene3', 'chrom': 'chr1', 'start': 1400, 'end': 1500},
        ]
        
        region = identify_best_synteny_block(hits, cluster_distance=500)
        
        self.assertIsNotNone(region)
        self.assertEqual(region['chrom'], 'chr1')
        self.assertEqual(region['genes_count'], 3)
        self.assertEqual(region['start'], 1000)
        self.assertEqual(region['end'], 1500)
    
    def test_multiple_clusters_choose_best(self):
        """Test choosing best cluster when multiple exist."""
        hits = [
            # Cluster 1: 2 genes
            {'query': 'gene1', 'chrom': 'chr1', 'start': 1000, 'end': 1100},
            {'query': 'gene2', 'chrom': 'chr1', 'start': 1200, 'end': 1300},
            # Cluster 2: 3 genes (should win)
            {'query': 'gene3', 'chrom': 'chr1', 'start': 10000, 'end': 10100},
            {'query': 'gene4', 'chrom': 'chr1', 'start': 10200, 'end': 10300},
            {'query': 'gene5', 'chrom': 'chr1', 'start': 10400, 'end': 10500},
        ]
        
        region = identify_best_synteny_block(hits, cluster_distance=500)
        
        self.assertIsNotNone(region)
        self.assertEqual(region['genes_count'], 3)
        self.assertEqual(region['start'], 10000)
        self.assertEqual(region['end'], 10500)
    
    def test_no_hits(self):
        """Test handling of empty hits list."""
        region = identify_best_synteny_block([], cluster_distance=500)
        self.assertIsNone(region)
    
    def test_different_chromosomes(self):
        """Test hits on different chromosomes."""
        hits = [
            {'query': 'gene1', 'chrom': 'chr1', 'start': 1000, 'end': 1100},
            {'query': 'gene2', 'chrom': 'chr2', 'start': 2000, 'end': 2100},
            {'query': 'gene3', 'chrom': 'chr1', 'start': 1200, 'end': 1300},
        ]
        
        region = identify_best_synteny_block(hits, cluster_distance=500)
        
        # Should pick chr1 cluster (2 genes)
        self.assertIsNotNone(region)
        self.assertEqual(region['chrom'], 'chr1')
        self.assertEqual(region['genes_count'], 2)


class TestRegionMerging(unittest.TestCase):
    """Test merging of overlapping synteny blocks."""
    
    def test_merge_overlapping_blocks(self):
        """Test merging of two blocks with overlapping search windows."""
        blocks = [
            {'chrom': 'chr1', 'start': 1000, 'end': 2000, 'score': 50, 'genes_count': 2},
            # Start 2100. Padding 200.
            # Block1 Search End = 2000 + 200 = 2200
            # Block2 Search Start = 2100 - 200 = 1900
            # 1900 <= 2200 -> Merge
            {'chrom': 'chr1', 'start': 2100, 'end': 3000, 'score': 60, 'genes_count': 3},
        ]
        merged = merge_synteny_blocks(blocks, padding=200)
        
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]['start'], 1000)
        self.assertEqual(merged[0]['end'], 3000)
        self.assertEqual(merged[0]['score'], 60) # Max score
        self.assertEqual(merged[0]['genes_count'], 5) # Sum counts

    def test_disjoint_blocks(self):
        """Test blocks that are too far apart to merge."""
        blocks = [
            {'chrom': 'chr1', 'start': 1000, 'end': 2000},
            # Start 3000. Padding 100.
            # B1 End = 2100
            # B2 Start = 2900
            # 2900 > 2100 -> distinct
            {'chrom': 'chr1', 'start': 3000, 'end': 4000},
        ]
        merged = merge_synteny_blocks(blocks, padding=100)
        self.assertEqual(len(merged), 2)

    def test_different_chromosomes(self):
        """Test blocks on different chromosomes never merge."""
        blocks = [
            {'chrom': 'chr1', 'start': 1000, 'end': 5000},
            {'chrom': 'chr2', 'start': 1000, 'end': 5000},
        ]
        merged = merge_synteny_blocks(blocks, padding=10000)
        self.assertEqual(len(merged), 2)

    def test_contained_block(self):
        """Test that a block fully contained in another search region is merged."""
        blocks = [
            {'chrom': 'chr1', 'start': 1000, 'end': 5000},
            {'chrom': 'chr1', 'start': 2000, 'end': 3000},
        ]
        merged = merge_synteny_blocks(blocks, padding=100)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]['start'], 1000)
        self.assertEqual(merged[0]['end'], 5000)


def run_tests():
    """Run all tests and return results."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestCoordinateNormalization))
    suite.addTests(loader.loadTestsFromTestCase(TestGeneIDExtraction))
    suite.addTests(loader.loadTestsFromTestCase(TestFastaIO))
    suite.addTests(loader.loadTestsFromTestCase(TestSequenceOperations))
    suite.addTests(loader.loadTestsFromTestCase(TestHitsFiltering))
    suite.addTests(loader.loadTestsFromTestCase(TestSyntenyIdentification))
    suite.addTests(loader.loadTestsFromTestCase(TestRegionMerging))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == '__main__':
    result = run_tests()
    sys.exit(0 if result.wasSuccessful() else 1)
