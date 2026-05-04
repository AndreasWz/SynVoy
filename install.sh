#!/bin/bash
# SynVoy Installation Script
# Sets up the Conda environment and verifies dependencies.
# Run from the SynVoy project root: ./install.sh

set -euo pipefail

echo "==========================================="
echo "  SynVoy Installation Script"
echo "==========================================="

ENV_NAME="synvoy_env"

if command -v mamba &> /dev/null; then
    CONDA_FRONTEND="mamba"
else
    CONDA_FRONTEND="conda"
fi

# ── Check for Conda ──────────────────────────────────────────
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda not found!"
    echo "Install Miniforge (recommended): https://github.com/conda-forge/miniforge"
    echo "  or Miniconda:                  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
echo "[ok] conda found"
echo "[ok] using ${CONDA_FRONTEND} for environment operations"

# ── Check for Java (required by Nextflow) ────────────────────
if command -v java &> /dev/null; then
    JAVA_LINE="$(java -version 2>&1 | head -1)"
    JAVA_MAJOR="$(echo "${JAVA_LINE}" | sed -E 's/.*version "([0-9]+).*/\1/' || true)"
    if [[ "${JAVA_MAJOR}" =~ ^[0-9]+$ ]] && [ "${JAVA_MAJOR}" -ge 17 ]; then
        echo "[ok] java found (${JAVA_LINE})"
    else
        echo "[!!] Java found but version is too old for modern Nextflow: ${JAVA_LINE}"
        echo "     Java >=17 will be installed inside the Conda environment."
    fi
else
    echo "[!!] Java not found — modern Nextflow requires Java >=17."
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
        "${CONDA_FRONTEND}" env remove -n "$ENV_NAME" -y
        "${CONDA_FRONTEND}" env create -f environment.yml
    else
        echo "     Keeping existing environment."
    fi
else
    "${CONDA_FRONTEND}" env create -f environment.yml
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
check_tool java
check_tool mmseqs
check_tool tblastn
check_tool makeblastdb
check_tool prodigal
check_tool augustus
check_tool miniprot
check_tool mafft
if command -v iqtree2 &> /dev/null || command -v iqtree &> /dev/null; then
    echo "  [ok] IQ-TREE ($(command -v iqtree2 || command -v iqtree))"
    PASS=$((PASS + 1))
else
    echo "  [FAIL] IQ-TREE not found (expected iqtree2 or iqtree)"
    FAIL=$((FAIL + 1))
fi
check_tool datasets
check_tool esearch
check_tool efetch
check_tool xtract

# Version checks that catch common fresh-solver traps.
if python - <<'PY' 2>/dev/null
import sys
raise SystemExit(0 if (3, 10) <= sys.version_info < (3, 13) else 1)
PY
then
    echo "  [ok] Python version ($(python -V 2>&1))"
    PASS=$((PASS + 1))
else
    echo "  [FAIL] Python must be >=3.10,<3.13 (ete3 is not Python 3.13-ready)"
    FAIL=$((FAIL + 1))
fi

if python - <<'PY' 2>/dev/null
import subprocess, re
line = subprocess.run(["java", "-version"], capture_output=True, text=True).stderr.splitlines()[0]
m = re.search(r'version "(\d+)', line)
raise SystemExit(0 if m and int(m.group(1)) >= 17 else 1)
PY
then
    echo "  [ok] Java version >=17"
    PASS=$((PASS + 1))
else
    echo "  [FAIL] Java >=17 is required by current Nextflow"
    FAIL=$((FAIL + 1))
fi

# Python packages
if python -c "import Bio, plotly, ete3, taxopy, parasail, psutil, numpy" 2>/dev/null; then
    echo "  [ok] Python packages (biopython, plotly, ete3, taxopy, parasail, psutil, numpy)"
    PASS=$((PASS + 1))
else
    echo "  [FAIL] One or more Python packages missing"
    FAIL=$((FAIL + 1))
fi

if nextflow config -profile standard >/dev/null; then
    echo "  [ok] Nextflow config parses"
    PASS=$((PASS + 1))
else
    echo "  [FAIL] Nextflow config parsing failed"
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
