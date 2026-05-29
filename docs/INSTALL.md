# Installation Guide

Full setup instructions for SynVoy. For a 5-line quick install, see the [README](../README.md#quick-install).

---

## Prerequisites

| Requirement | Notes |
|---|---|
| **OS** | Linux (tested) or macOS |
| **Conda or Mamba** | [Miniforge](https://github.com/conda-forge/miniforge) recommended. Miniconda/Anaconda also work. |
| **Git** | To clone the repository |
| **Internet** | Easy Mode needs access to NCBI/UniProt for genome downloads |

> **Note on Java and Nextflow:** Both are provided by the `synvoy_env`
> conda environment (`openjdk >=17` and `nextflow >=25.10`). You do **not**
> need to install them separately. See the
> [HPC / standalone Nextflow](#hpc--standalone-nextflow-no-conda) section
> below if your environment forbids conda.

---

## 1. Clone the Repository

```bash
git clone https://github.com/AndreasWz/SynVoy.git
cd SynVoy
```

## 2. Set Up the Conda Environment

The environment bundles **Nextflow**, **OpenJDK 17**, all bioinformatics tools (MMseqs2, BLAST, Prodigal, Augustus, miniprot, MAFFT, IQ-TREE), genome-fetching CLIs (NCBI datasets, Entrez Direct), and all Python dependencies. This is the single canonical install path.

```bash
# Create the environment (mamba is faster if available)
mamba env create -f environment.yml
# or: conda env create -f environment.yml

# Activate it (must be done in every new terminal before running the pipeline)
conda activate synvoy_env
```

> The environment is named `synvoy_env` (defined in `environment.yml`).

## 3. Verify the Installation

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

## HPC / standalone Nextflow (no conda)

If your environment forbids Conda (some HPC systems load Nextflow as a
module), you can run SynVoy against a system Nextflow as long as you also
provide Java ≥17 and let Nextflow use the Docker or Singularity profile
to pick up all the bioinformatics tools:

```bash
# system / module-loaded Nextflow
nextflow -version    # must report ≥25.10 and Java ≥17
nextflow run main.nf -profile singularity --mode easy --query_id Q16553 --outdir results
```

In that path you do **not** need `environment.yml`. The conda env is
optional only when you go through `-profile docker` or `-profile singularity`.
With `-profile standard`, the conda env is mandatory.

---

## Troubleshooting

See [USAGE.md](USAGE.md#troubleshooting) for common installation and runtime issues.
