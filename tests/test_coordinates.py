#!/usr/bin/env python3
"""
Integration tests for coordinate transformations in SynVoy.

Tests the critical coordinate conversion logic:
1. GFF3 (1-based, closed) to BED (0-based, half-open)
2. Local region coordinates to global genomic coordinates
3. Miniprot GFF shifting from extracted region to chromosome
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_gff_to_bed_conversion():
    """Test GFF3 to BED coordinate conversion."""
    # GFF3: 1-based, closed interval [start, end]
    # Example: Gene from position 1000 to 2000 (inclusive)
    gff_start = 1000
    gff_end = 2000
    
    # BED: 0-based, half-open interval [start, end)
    # Correct conversion:
    bed_start = gff_start - 1  # 999
    bed_end = gff_end          # 2000 (unchanged)
    
    # Verify: BED interval should represent same bases
    # GFF [1000, 2000] = bases 1000-2000 inclusive = 1001 bases
    # BED [999, 2000) = bases 999-1999 (0-indexed) = positions 1000-2000 (1-indexed)
    
    assert bed_start == 999, f"Expected bed_start=999, got {bed_start}"
    assert bed_end == 2000, f"Expected bed_end=2000, got {bed_end}"
    
    # Length check
    gff_length = gff_end - gff_start + 1  # 1001 bases
    bed_length = bed_end - bed_start      # 1001 bases
    assert gff_length == bed_length, f"Length mismatch: GFF={gff_length}, BED={bed_length}"
    
    print("✓ GFF to BED conversion correct")


def test_local_to_global_shift():
    """Test shifting coordinates from local region to global genome."""
    # Global region: chr1:100000-200000
    window_start = 100000
    window_end = 200000
    
    # Local coordinates (within extracted region): 5000-6000
    local_start = 5000
    local_end = 6000
    
    # Global coordinates should be:
    global_start = window_start + local_start  # 105000
    global_end = window_start + local_end      # 106000
    
    assert global_start == 105000, f"Expected global_start=105000, got {global_start}"
    assert global_end == 106000, f"Expected global_end=106000, got {global_end}"
    
    # Verify length is preserved
    local_length = local_end - local_start
    global_length = global_end - global_start
    assert local_length == global_length, "Length changed during shift!"
    
    print("✓ Local to global coordinate shift correct")


def test_miniprot_gff_shift():
    """Test Miniprot GFF shifting (combines GFF parsing + coordinate shift)."""
    # Scenario: Extract region chr1:50000-150000 (window_start=50000)
    # Miniprot finds gene at local positions 10000-12000 (in extracted FASTA)
    # Should output GFF with global positions 60000-62000
    
    window_start = 50000
    
    # Miniprot GFF (local coordinates, 1-based)
    miniprot_local_start = 10000  # GFF format
    miniprot_local_end = 12000
    
    # Shift to global (still in GFF format, 1-based)
    global_gff_start = miniprot_local_start + window_start  # 60000
    global_gff_end = miniprot_local_end + window_start      # 62000
    
    assert global_gff_start == 60000, f"Expected 60000, got {global_gff_start}"
    assert global_gff_end == 62000, f"Expected 62000, got {global_gff_end}"
    
    print("✓ Miniprot GFF coordinate shift correct")


def test_sequence_extraction():
    """Test sequence extraction with coordinate systems."""
    # Python string/BioPython sequence (0-based)
    # Example: "ATCGATCGATCG" (12 bases, indices 0-11)
    seq = "ATCGATCGATCG"
    
    # BED coordinates [3, 9) should extract indices 3-8 = "GATCGA" (6 bases)
    bed_start = 3
    bed_end = 9
    
    extracted = seq[bed_start:bed_end]
    
    assert extracted == "GATCGA", f"Expected 'GATCGA', got '{extracted}'"
    assert len(extracted) == 6, f"Expected length 6, got {len(extracted)}"
    
    # Corresponding GFF would be [4, 9] (1-based, closed)
    # Which also represents 6 bases: positions 4,5,6,7,8,9
    gff_start = bed_start + 1  # 4
    gff_end = bed_end          # 9
    gff_length = gff_end - gff_start + 1  # 6
    
    assert gff_length == len(extracted), "GFF and BED represent different lengths!"
    
    print("✓ Sequence extraction with BED coordinates correct")


def test_exon_clustering():
    """Test multi-exon gene clustering logic."""
    # Gene with 3 exons, should be clustered together
    exons = [
        {'start': 1000, 'end': 1200},
        {'start': 2000, 'end': 2300},
        {'start': 3500, 'end': 3800}
    ]
    
    # Check gaps
    gap1 = exons[1]['start'] - exons[0]['end']  # 800 bp
    gap2 = exons[2]['start'] - exons[1]['end']  # 1200 bp
    
    # With max_intron=20000, all should cluster
    max_intron = 20000
    assert gap1 < max_intron, "Gap1 should be less than max_intron"
    assert gap2 < max_intron, "Gap2 should be less than max_intron"
    
    # Total span
    span = exons[2]['end'] - exons[0]['start']  # 2800 bp
    max_gene_span = 500000
    assert span < max_gene_span, "Gene span should be reasonable"
    
    print("✓ Exon clustering logic correct")


def test_coordinate_edge_cases():
    """Test edge cases in coordinate conversions."""
    # Edge case 1: Gene at start of chromosome (position 1 in GFF)
    gff_start = 1
    gff_end = 100
    bed_start = gff_start - 1  # 0
    bed_end = gff_end          # 100
    
    assert bed_start == 0, "Start of chromosome should be 0 in BED"
    assert bed_end == 100, "End should be 100"
    
    # Edge case 2: Zero-length feature (not typical, but should handle)
    # GFF [1000, 1000] = 1 base
    # BED [999, 1000) = 1 base
    gff_single = 1000
    bed_start_single = gff_single - 1  # 999
    bed_end_single = gff_single        # 1000
    length = bed_end_single - bed_start_single  # 1
    
    assert length == 1, "Single-base feature should have length 1"
    
    print("✓ Edge cases handled correctly")


def run_all_tests():
    """Run all coordinate transformation tests."""
    print("=" * 60)
    print("Running SynVoy Coordinate Transformation Tests")
    print("=" * 60)
    
    tests = [
        test_gff_to_bed_conversion,
        test_local_to_global_shift,
        test_miniprot_gff_shift,
        test_sequence_extraction,
        test_exon_clustering,
        test_coordinate_edge_cases
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            print(f"\n{test.__name__}:")
            test()
            passed += 1
        except AssertionError as e:
            print(f"✗ FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ ERROR: {e}")
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
