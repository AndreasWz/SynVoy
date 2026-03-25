#!/bin/bash
# SynVoy Installation Script
# Sets up the Conda environment and verifies dependencies.
# Run from the SynVoy project root: ./install.sh

set -e

echo "==========================================="
echo "  SynVoy Installation Script"
echo "==========================================="

ENV_NAME="synvoy_env"

# ── Check for Conda ──────────────────────────────────────────
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda not found!"
    echo "Install Miniforge (recommended): https://github.com/conda-forge/miniforge"
    echo "  or Miniconda:                  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
echo "[ok] conda found"

# ── Check for Java (required by Nextflow) ────────────────────
if command -v java &> /dev/null; then
    echo "[ok] java found ($(java -version 2>&1 | head -1))"
else
    echo "[!!] Java not found — Nextflow requires Java >=11."
    echo "     It will be installed inside the Conda environment."
fi

# ── Check for Nextflow ───────────────────────────────────────
if command -v nextflow &> /dev/null; then
    echo "[ok] nextflow found ($(nextflow -version 2>&1 | grep 'version' | head -1 | awk '{print $NF}'))"
else
    echo "[!!] Nextflow not found. It will be installed inside the Conda environment."
fi

# ── Create Conda environment ────────────────────────────────
echo ""
echo "Creating Conda environment '${ENV_NAME}'..."
if conda env list | grep -qw "$ENV_NAME"; then
    echo "[!!] Environment '${ENV_NAME}' already exists."
    read -p "     Remove and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n "$ENV_NAME" -y
        conda env create -f environment.yml
    else
        echo "     Keeping existing environment."
    fi
else
    conda env create -f environment.yml
fi

echo ""
echo "[ok] Environment ready"

# ── Verify key tools ────────────────────────────────────────
echo ""
echo "Verifying dependencies inside '${ENV_NAME}'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

PASS=0
FAIL=0
check_tool() {
    if command -v "$1" &> /dev/null; then
        echo "  [ok] $1"
        PASS=$((PASS + 1))
    else
        echo "  [FAIL] $1 not found"
        FAIL=$((FAIL + 1))
    fi
}

check_tool nextflow
check_tool mmseqs
check_tool tblastn
check_tool makeblastdb
check_tool prodigal
check_tool miniprot
check_tool mafft
check_tool iqtree2

# Python packages
if python -c "import Bio, plotly, ete3, taxopy, parasail" 2>/dev/null; then
    echo "  [ok] Python packages (biopython, plotly, ete3, taxopy, parasail)"
    PASS=$((PASS + 1))
else
    echo "  [FAIL] One or more Python packages missing"
    FAIL=$((FAIL + 1))
fi

echo ""
if [ "$FAIL" -eq 0 ]; then
    echo "==========================================="
    echo "  All $PASS checks passed!"
    echo "==========================================="
else
    echo "==========================================="
    echo "  $PASS passed, $FAIL failed"
    echo "  Try: conda env remove -n ${ENV_NAME} && conda env create -f environment.yml"
    echo "==========================================="
    exit 1
fi

echo ""
echo "To get started:"
echo ""
echo "  conda activate ${ENV_NAME}"
echo ""
echo "  # Easy Mode (auto-fetch genomes):"
echo "  nextflow run main.nf --mode easy --query_id Q16553 --max_genomes 5 --outdir results -profile standard"
echo ""
echo "  # Pro Mode (local files):"
echo "  nextflow run main.nf --mode pro --query query.faa --home_genome genome.fna --target_genomes 'targets/*.fna' --outdir results -profile standard"
echo ""
echo "  # See USAGE.md for the full parameter reference."
echo ""
