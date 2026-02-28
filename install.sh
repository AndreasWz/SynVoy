#!/bin/bash
# Installation script for SynTerra
# Run this to set up your environment

set -e

echo "==========================================="
echo "  SynTerra Installation Script"
echo "==========================================="

# Check for conda
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: Conda not found!"
    echo "Please install Miniconda or Anaconda first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✓ Found conda"

# Check for nextflow
if ! command -v nextflow &> /dev/null; then
    echo "⚠️  Nextflow not found. Installing..."
    conda install -c bioconda nextflow -y
else
    echo "✓ Found nextflow"
fi

# Create conda environment
echo ""
echo "Creating conda environment 'syntenyfinder'..."
if conda env list | grep -q syntenyfinder; then
    echo "⚠️  Environment 'syntenyfinder' already exists"
    read -p "Remove and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n syntenyfinder -y
        conda env create -f environment.yml
    fi
else
    conda env create -f environment.yml
fi

echo ""
echo "✓ Environment created successfully"

# Test installation
echo ""
echo "Testing installation..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate syntenyfinder

# Check key tools
echo "Checking dependencies..."
for tool in mmseqs blast prodigal miniprot plotly; do
    if command -v $tool &> /dev/null || python -c "import $tool" 2>/dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ⚠️  $tool not found (may not be critical)"
    fi
done

echo ""
echo "==========================================="
echo "  Installation Complete!"
echo "==========================================="
echo ""
echo "To use SynTerra:"
echo "  1. Activate environment:"
echo "     conda activate syntenyfinder"
echo ""
echo "  2. Test with example data:"
echo "     nextflow run main.nf -profile test"
echo ""
echo "  3. Run with your data:"
echo "     nextflow run main.nf --mode pro --query query.fasta --home_genome genome.fna --target_genomes 'targets/*.fna'"
echo ""
echo "For help:"
echo "  - See README.md for overview"
echo "  - See USAGE.md for detailed usage"
echo "  - Run: nextflow run main.nf --help"
echo ""
