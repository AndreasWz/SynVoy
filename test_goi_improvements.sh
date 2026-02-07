#!/bin/bash
# Quick test of GOI annotation improvements

set -e

echo "==================================="
echo "Testing GOI Annotation Improvements"
echo "==================================="

# Test 1: Check if new functions are present
echo ""
echo "[Test 1] Checking for new functions in iterative_search_runner.py..."
if grep -q "def orf_based_annotation" bin/iterative_search_runner.py; then
    echo "  ✓ orf_based_annotation() found"
else
    echo "  ✗ orf_based_annotation() NOT found"
    exit 1
fi

if grep -q "def find_orfs_six_frame" bin/iterative_search_runner.py; then
    echo "  ✓ find_orfs_six_frame() found"
else
    echo "  ✗ find_orfs_six_frame() NOT found"
    exit 1
fi

if grep -q "sensitive: bool = True" bin/iterative_search_runner.py; then
    echo "  ✓ Sensitive miniprot parameter found"
else
    echo "  ✗ Sensitive parameter NOT found"
    exit 1
fi

# Test 2: Check if integration script exists
echo ""
echo "[Test 2] Checking for augmented GOI integration script..."
if [ -f "bin/integrate_augmented_goi.py" ]; then
    echo "  ✓ integrate_augmented_goi.py exists"
    if [ -x "bin/integrate_augmented_goi.py" ]; then
        echo "  ✓ Script is executable"
    else
        echo "  ✗ Script is not executable"
        exit 1
    fi
else
    echo "  ✗ Script NOT found"
    exit 1
fi

# Test 3: Check documentation
echo ""
echo "[Test 3] Checking documentation..."
if [ -f "docs/GOI_VS_FLANKING_ARCHITECTURE.md" ]; then
    echo "  ✓ Architecture documentation exists"
else
    echo "  ✗ Architecture doc NOT found"
    exit 1
fi

if [ -f "docs/GOI_ANNOTATION_IMPROVEMENTS.md" ]; then
    echo "  ✓ Improvements documentation exists"
else
    echo "  ✗ Improvements doc NOT found"
    exit 1
fi

# Test 4: Syntax check Python scripts
echo ""
echo "[Test 4] Checking Python syntax..."
if python3 -m py_compile bin/iterative_search_runner.py 2>/dev/null; then
    echo "  ✓ iterative_search_runner.py syntax OK"
else
    echo "  ✗ Syntax errors in iterative_search_runner.py"
    exit 1
fi

if python3 -m py_compile bin/integrate_augmented_goi.py 2>/dev/null; then
    echo "  ✓ integrate_augmented_goi.py syntax OK"
else
    echo "  ✗ Syntax errors in integrate_augmented_goi.py"
    exit 1
fi

# Test 5: Check for key improvements
echo ""
echo "[Test 5] Checking for key code improvements..."

if grep -q "ORF-based fallback" bin/iterative_search_runner.py; then
    echo "  ✓ ORF fallback integration found"
else
    echo "  ✗ ORF fallback NOT integrated"
    exit 1
fi

if grep -q "full_length_goi" bin/iterative_search_runner.py; then
    echo "  ✓ Full-length GOI preference found"
else
    echo "  ✗ Full-length preference NOT found"
    exit 1
fi

if grep -q "\-p.*0\\.4" bin/iterative_search_runner.py; then
    echo "  ✓ Sensitive miniprot parameters found"
else
    echo "  ✗ Sensitive parameters NOT found"
    exit 1
fi

# Summary
echo ""
echo "==================================="
echo "All tests passed! ✓"
echo "==================================="
echo ""
echo "Next steps:"
echo "  1. Re-run melettin test: nextflow run main.nf -profile test_melettin,conda -resume"
echo "  2. Check logs for 'ORF-based fallback' messages"
echo "  3. Verify database expansion in iterative_results/"
echo "  4. Check if augmented candidates can be annotated"
echo ""
