# Installation Guide

Full setup instructions for SynVoy. For a 5-line quick install, see the [README](../README.md#quick-install).

---

## Prerequisites

| Requirement | Notes |
|---|---|
| **OS** | Linux (tested) or macOS |
| **Java** | 17 or newer (required by current Nextflow; the Conda env includes OpenJDK) |
| **Conda or Mamba** | [Miniforge](https://github.com/conda-forge/miniforge) recommended. Miniconda/Anaconda also work. |
| **Git** | To clone the repository |
| **Internet** | Easy Mode needs access to NCBI/UniProt for genome downloads |

---

## 1. Clone the Repository

```bash
git clone https://github.com/AndreasWz/SynVoy.git
cd SynVoy
```

## 2. Install Nextflow

Check if Nextflow is already installed:

```bash
nextflow -version
```

If the command is not found, install it:

```bash
# Option A: using Conda (simplest — it will be included in the env in step 3)
# Skip this step; the environment.yml already lists nextflow.

# Option B: standalone install into ~/bin
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin && mv nextflow ~/bin/
# Make sure ~/bin is on your PATH:
export PATH="$HOME/bin:$PATH"
# (add the line above to ~/.bashrc to make it permanent)
```

Current Nextflow requires **Java ≥17**. Verify with `java -version`. If missing, let Conda pull OpenJDK in with the environment below, or install Java 17+ via your system package manager.

## 3. Set Up the Conda Environment

The environment bundles Nextflow, all bioinformatics tools (MMseqs2, BLAST, Prodigal, miniprot, MAFFT, IQ-TREE), genome-fetching CLIs (NCBI datasets, Entrez Direct), and all Python dependencies.

```bash
# Create the environment (mamba is faster if available)
mamba env create -f environment.yml
# or: conda env create -f environment.yml

# Activate it
conda activate synvoy_env
```

> The environment is named `synvoy_env` (defined in `environment.yml`). You must activate it every time you open a new terminal before running the pipeline.

## 4. Verify the Installation

```bash
# All of these should print version info without errors:
nextflow -version
java -version
mmseqs version
tblastn -version
miniprot --version
prodigal -v
mafft --version
augustus --version
iqtree2 --version || iqtree --version
datasets version
esearch -version
python -c "import Bio; import plotly; import ete3; import taxopy; import parasail; import psutil; print('Python deps OK')"
```

If any tool is missing, re-create the environment:

```bash
conda env remove -n synvoy_env
conda env create -f environment.yml
```

---

## Alternative: Docker

If you prefer containers over Conda:

```bash
# Build the image (all tools are baked in)
docker build -t synvoy-local:latest .

# Run with the docker profile (no Conda needed)
nextflow run main.nf -profile docker --mode easy --query_id Q16553 --outdir results
```

For Singularity (common on HPC):

```bash
nextflow run main.nf -profile singularity --mode easy --query_id Q16553 --outdir results
```

---

## Troubleshooting

See [USAGE.md](USAGE.md#troubleshooting) for common installation and runtime issues.
