# SynVoy - Synteny Voyager

<p align="center">
  <img src="assets/logo.png" alt="SynVoy logo" width="220"/>
</p>

*Navigating synteny. Discovering orthology.*  
*Mapping evolutionary pathways through syntenic navigation.*

SynVoy is a Nextflow pipeline for finding orthologous genes across evolutionary distances using genomic synteny.

Standard sequence-similarity searches often fail when orthologs are highly divergent or consist of short, complex micro-exons. SynVoy addresses this by leveraging the conservation of gene order (macro-synteny): it identifies the conserved flanking genes around a Gene of Interest (GOI) in a reference species, locates the homologous genomic neighborhood in target species, and then runs a localized sequence search to find the GOI candidate.

> **Status:** Early development. Expect breaking changes between versions.

## Table of Contents

- [How It Works](#how-it-works)
- [Setup from Scratch](#setup-from-scratch)
  - [Prerequisites](#prerequisites)
  - [1. Clone the Repository](#1-clone-the-repository)
  - [2. Install Nextflow](#2-install-nextflow)
  - [3. Set Up the Conda Environment](#3-set-up-the-conda-environment)
  - [4. Verify the Installation](#4-verify-the-installation)
  - [Alternative: Docker](#alternative-docker)
- [Quick Start](#quick-start)
  - [Easy Mode](#easy-mode-automated-genome-retrieval)
  - [Pro Mode](#pro-mode-local-files)
- [Output](#output)
- [Further Reading](#further-reading)
- [License](#license)
- [Citation and Support](#citation-and-support)

---

## How It Works

1. **Input Resolution** — Accepts a UniProt/NCBI accession, a local FASTA, or an inline FASTA sequence (Easy Mode) and resolves it to a protein query.
2. **Genome Staging** — In Easy Mode, automatically fetches the reference ("home") genome and related target assemblies from NCBI. In Pro Mode, the user supplies local files.
3. **Gene Localization** — Maps the GOI onto the home genome with tblastn + MMseqs2 and annotates its exon structure (from GFF or *de novo* via Prodigal).
4. **Flanking Gene Extraction** — Extracts the *n* genes immediately upstream and downstream of the GOI locus. Genes that are similar to the GOI (e.g. tandem duplicates) are filtered out of the flanking set to avoid inflating synteny scores. Optionally, those GOI-similar neighbors are emitted as additional GOI queries (`--expand_goi_similar`), so that paralogs in other genomes are also discovered.
5. **Phylogenetic Ordering** — Sorts target genomes by evolutionary distance to the reference so that the iterative search proceeds from closest to most distant relatives.
6. **Iterative Synteny Search** — For each target genome, maps flanking genes with MMseqs2, clusters hits into candidate syntenic blocks, and runs localized tblastn + miniprot + Smith-Waterman searches inside those blocks to find the GOI (and any GOI-similar neighbor queries).
7. **Region Clustering & Scoring** — Filters and ranks candidate blocks by synteny score (fraction of conserved flanking genes).
8. **Phylogenetic Tree & Visualization** — Aligns all discovered GOI sequences across all genomes (MAFFT) and infers a phylogenetic tree (IQ-TREE with ultrafast bootstrap). When `--expand_goi_similar` is enabled, the tree includes paralogs and orthologs together, enabling resolution of duplication vs. speciation events. An interactive HTML synteny plot (Plotly) is generated alongside the tree.

---

## Setup from Scratch

### Prerequisites

| Requirement | Notes |
|---|---|
| **OS** | Linux (tested) or macOS |
| **Java** | 11 or 17 (required by Nextflow) |
| **Conda or Mamba** | [Miniforge](https://github.com/conda-forge/miniforge) recommended. Miniconda/Anaconda also work. |
| **Git** | To clone the repository |
| **Internet** | Easy Mode needs access to NCBI/UniProt for genome downloads |

### 1. Clone the Repository

```bash
git clone https://github.com/AndreasWz/SynVoy.git
cd SynVoy
```

### 2. Install Nextflow

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

Nextflow requires **Java ≥11**. Verify with `java -version`. If missing, install via your system package manager (`sudo apt install default-jdk` on Debian/Ubuntu) or let Conda pull it in with the environment below.

### 3. Set Up the Conda Environment

The environment bundles Nextflow, all bioinformatics tools (MMseqs2, BLAST, Prodigal, miniprot, MAFFT, IQ-TREE), genome-fetching CLIs (NCBI datasets, Entrez Direct), and all Python dependencies.

```bash
# Create the environment (use mamba for speed if available)
conda env create -f environment.yml
# or: mamba env create -f environment.yml

# Activate it
conda activate synvoy_env
```

> The environment is named `synvoy_env` (defined in `environment.yml`). You must activate it every time you open a new terminal before running the pipeline.

### 4. Verify the Installation

```bash
# All of these should print version info without errors:
nextflow -version
mmseqs version
tblastn -version
miniprot --version
prodigal -v
mafft --version
iqtree2 --version
python -c "import Bio; import plotly; import ete3; import taxopy; import parasail; print('Python deps OK')"
```

If any tool is missing, re-create the environment:

```bash
conda env remove -n synvoy_env
conda env create -f environment.yml
```

### Alternative: Docker

If you prefer containers over Conda:

```bash
# Build the image (all tools are baked in)
docker build -t synvoy-local:latest .

# Run with the docker profile (no Conda needed)
nextflow run main.nf -profile docker --query_id Q16553 --outdir results
```

For Singularity (common on HPC):

```bash
nextflow run main.nf -profile singularity --query_id Q16553 --outdir results
```

---

## Quick Start

### Easy Mode (automated genome retrieval)

Provide a UniProt/NCBI protein accession, a local FASTA (`--query`), or an inline sequence (`--query_seq`). SynVoy fetches the reference genome, downloads related target assemblies, and runs the full analysis:

```bash
nextflow run main.nf \
  --mode easy \
  --query_id Q16553 \
  --max_genomes 5 \
  --outdir results/ly6e_easy \
  -profile standard
```

Optional flags:

- `--home_species "Homo sapiens"` — override auto-detected species.
- `--target_species "Gallus gallus,Mus musculus"` — specify target species instead of auto-selecting.
- `--query_seq "MKT..."` — inline protein sequence; requires `--home_species`.

### Pro Mode (local files)

Supply your own query FASTA, reference genome, and target genomes:

```bash
nextflow run main.nf \
  --mode pro \
  --query queries/melittin.faa \
  --home_genome /path/to/apis_mellifera.fna \
  --home_gff /path/to/apis_mellifera.gff \
  --target_genomes "/path/to/targets/*.fna" \
  --outdir results/melittin_pro \
  -profile standard
```

> `--home_gff` is optional but strongly recommended — it provides much better flanking-gene extraction than Prodigal prediction alone.

> Use `-resume` to restart from the last successful step after a crash or parameter tweak.

---

## Output

Results are written to the directory specified by `--outdir`:

| File | Description |
|---|---|
| `*_synteny_plot.html` | Interactive HTML visualization of syntenic blocks across species |
| `*_tree.nwk` | Newick phylogenetic tree of all discovered GOI and GOI-similar sequences across genomes |
| `regions/*.regions.bed` | BED files with genomic coordinates of candidate syntenic blocks |
| `synvoy_report.json` | Machine-readable run summary (parameters, genome QC, exit codes) |
| `intermediate/` | Per-phase artifacts (flanking genes, MMseqs2 hits, GFFs, etc.) |

---

## Further Reading

- **[USAGE.md](USAGE.md)** — Full parameter reference, execution profiles, algorithm details, HPC/SLURM setup, and troubleshooting.


## License

SynVoy is distributed under the [GNU AGPLv3](LICENSE) License.

---

## Citation and Support

If SynVoy contributes to your research, please cite the software repository.

Recommended software citation:

> Weitz, F. A. SynVoy: Synteny-guided orthology discovery [Computer software]. GitHub. https://github.com/AndreasWz/SynVoy

You can also use this BibTeX template:

```bibtex
@software{synvoy,
  author  = {Weitz, Frank Andreas},
  title   = {SynVoy: Synteny-guided orthology discovery},
  year    = {2026},
  url     = {https://github.com/AndreasWz/SynVoy},
  note    = {GitHub repository. Accessed 2026}
}
```

If SynVoy is useful for your work, please consider starring the repository: https://github.com/AndreasWz/SynVoy
