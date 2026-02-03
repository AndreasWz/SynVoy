#!/bin/bash
set -e

# SynTerra Integration Test
# Run with: ./integration_test.sh [QUERY_FILE] [TARGET_DIR] [HOME_GENOME]

QUERY=${1:-"test_data/query.fasta"}
TARGET_DIR=${2:-"test_data/targets/"}
HOME_GENOME=${3:-"test_data/home.fasta"}
OUTDIR="results_test"

echo "=== SynTerra Integration Test ==="
echo "Query: $QUERY"
echo "Targets: $TARGET_DIR"
echo "Output: $OUTDIR"

# Cleanup previous run
rm -rf $OUTDIR .nextflow* work

# Run Pipeline
# Use 'stub' or 'test' profile if available, else run normal
# Assuming we run normal execution on small test data
echo "Running Nextflow..."
nextflow run main.nf \
    --gene "$QUERY" \
    --target_genomes "$TARGET_DIR/*.fna" \
    --home_genome "$HOME_GENOME" \
    --mode pro \
    --outdir "$OUTDIR" \
    -resume

echo "=== Verification ==="

# 1. Check Region Files
echo "Checking regions..."
REGIONS_FOUND=$(find "$OUTDIR" -name "*.regions.bed" | wc -l)
if [ "$REGIONS_FOUND" -gt 0 ]; then
    echo "  [PASS] Found $REGIONS_FOUND region BED files."
else
    echo "  [FAIL] No region BED files found."
    exit 1
fi

# 2. Check for Empty/Missing Augmented Results
# Check if we have augmented outputs for at least one genome
AUG_BEDS=$(find "$OUTDIR/augmented" -name "*.candidates.bed" | wc -l)
if [ "$AUG_BEDS" -gt 0 ]; then
    echo "  [PASS] Found $AUG_BEDS augmented search BED files."
else
    echo "  [FAIL] No augmented search BED files found."
    exit 1
fi

# 3. Check Report
REPORT="$OUTDIR/synterra_report.html"
if [ -f "$REPORT" ]; then
    echo "  [PASS] Report generated."
    
    # Check for gene count > 0 in report
    # Simple grep for now, assuming JSON embedded or HTML table
    # Just checking file size > 0
    if [ -s "$REPORT" ]; then
        echo "  [PASS] Report is not empty."
    else
        echo "  [FAIL] Report is empty."
        exit 1
    fi
else
    echo "  [FAIL] Report not found."
    exit 1
fi

# 4. Check for Critical Errors in Log
if [ -f ".nextflow.log" ]; then
    ERRORS=$(grep -i "ERROR" .nextflow.log | grep -v "ErrorStrategy" | grep -v "ignore" | wc -l)
    if [ "$ERRORS" -eq 0 ]; then
        echo "  [PASS] No critical errors found in log."
    else
        echo "  [WARN] Found $ERRORS errors in log (check manually)."
        grep -i "ERROR" .nextflow.log | grep -v "ErrorStrategy" | head -n 5
    fi
fi

echo "=== Integration Test Complete ==="
exit 0
