# SynTerra

**Phylogenetically-informed syntenic ortholog discovery across divergent genomes**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

SynTerra discovers orthologous genes in newly sequenced genomes by exploiting **synteny conservation** — the principle that gene order is preserved across evolution. Unlike traditional BLAST-based approaches, SynTerra:

- **Uses genomic context**: Finds genes by their conserved neighborhood, not just sequence similarity
- **Iterative phylogenetic search**: Uses evolutionarily close species as "stepping stones" to find distant orthologs
- **High sensitivity**: Combines MMseqs2 with Smith-Waterman alignment for rigorous gene detection
- **Handles edge cases**: Finds pseudogenes, frameshifted genes, and highly divergent orthologs

## Key Innovation

The **iterative phylogenetic search** progressively builds a search database:

```
Query gene → Close species (90% ID) → Medium species (70% ID) → Distant species (50% ID)
                    ↓ add to DB             ↓ add to DB              ↓
              Direct search FAILS ────────────────────────────────────┘
              Iterative search SUCCEEDS ✓
```

Each discovered ortholog is added to the database, enabling detection of increasingly divergent genes.

## Quick Start

### Easy Mode (Recommended)

Just provide a UniProt ID and species name — SynTerra fetches everything automatically:

```bash
# Find LY6E orthologs across primates
nextflow run main.nf \
  --query_id Q16553 \
  --home_species "Homo sapiens" \
  --max_genomes 5 \
  --outdir results/LY6E

# Find melettin orthologs across bee species  
nextflow run main.nf \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/melettin
```

### Pro Mode (Custom Genomes)

Provide your own genome files:

```bash
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome reference.fna \
  --home_gff reference.gff \
  --target_genomes "genomes/*.fna" \
  --mode pro \
  --outdir results
```

## Installation

### Requirements
- **Nextflow** ≥ 22.10.1
- **Java** 11+ 
- **Conda/Mamba** or **Singularity/Docker**

### Setup

```bash
# Clone repository
git clone https://github.com/yourusername/SynTerra.git
cd SynTerra

# Install Nextflow (if needed)
curl -s https://get.nextflow.io | bash

# Run with Conda (creates environment automatically)
nextflow run main.nf -profile conda --query_id P01501 --home_species "Apis mellifera"

# Or with Singularity
nextflow run main.nf -profile singularity --query_id P01501 --home_species "Apis mellifera"
```

### Easy Mode Requirements

For automatic genome fetching, install NCBI tools:
```bash
conda install -c bioconda entrez-direct
conda install -c conda-forge ncbi-datasets-cli
```

## Parameters

### Core Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--query_id` | UniProt ID (e.g., P01501, Q16553) | - |
| `--gene` | Query gene FASTA (alternative to query_id) | - |
| `--home_species` | Species name (Easy mode) | - |
| `--max_genomes` | Number of related genomes to analyze | 10 |
| `--mode` | `easy` or `pro` | easy |
| `--outdir` | Output directory | results |

### Search Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--n_flanking_genes` | Flanking genes for synteny context | 10 |
| `--mmseqs_sensitivity` | MMseqs2 sensitivity (1-8.5) | 8.5 |
| `--enable_smith_waterman` | Use Smith-Waterman for GOI | true |
| `--sw_min_identity` | Minimum SW alignment identity | 20.0 |

## Output

```
results/
├── query/
│   └── Q16553.fasta              # Query sequence
├── genomes/
│   └── *.fna                     # Downloaded genomes
├── flanking/
│   └── flanking_proteins.faa     # Flanking gene proteins
├── search/
│   ├── *_hits.gff                # Gene predictions per genome
│   └── *_regions.bed             # Syntenic regions
├── trees/
│   └── gene_tree.nwk             # Phylogenetic tree (GOI only)
├── plots/
│   └── synteny_plot.html         # Interactive visualization
└── report.json                   # Pipeline summary
```

## Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 1: Input Preparation                                     │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ FETCH    │  │ FETCH_HOME    │  │ FETCH_RELATED           │  │
│  │ QUERY    │→ │ GENOME        │→ │ GENOMES (NCBI)          │  │
│  └──────────┘  └───────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 2: Synteny Context                                        │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ LOCATE   │  │ SPLIT         │  │ EXTRACT                 │  │
│  │ GENE     │→ │ LOCI          │→ │ FLANKING                │  │
│  └──────────┘  └───────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 3: Iterative Search                                       │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ PHYLO    │  │ ITERATIVE     │  │ CLUSTER                 │  │
│  │ SORT     │→ │ SEARCH        │→ │ REGIONS                 │  │
│  └──────────┘  │ (MMseqs2+SW)  │  └─────────────────────────┘  │
│                └───────────────┘                                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 4: Phylogenetics & Visualization                          │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ COMPUTE  │  │ PLOT          │  │ GENERATE                │  │
│  │ TREE     │→ │ SYNTENY       │→ │ REPORT                  │  │
│  └──────────┘  └───────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

## Running on HPC (SLURM)

### Quick Start

```bash
# Clone to cluster
git clone https://github.com/yourusername/SynTerra.git
cd SynTerra

# Submit job
sbatch slurm_submit.sh
```

### Manual Execution

```bash
# Load modules
module load java/17
module load singularity

# Run with SLURM + Singularity
nextflow run main.nf \
  -profile hpc_singularity \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/my_run \
  -work-dir $SCRATCH/work
```

### Available Profiles

| Profile | Description |
|---------|-------------|
| `hpc_singularity` | SLURM + Singularity (recommended for HPC) |
| `hpc_conda` | SLURM + Conda |
| `conda` | Local + Conda |
| `singularity` | Local + Singularity |
| `docker` | Local + Docker |

## Examples

### Insect Venom Peptides
```bash
# Melettin in bees
nextflow run main.nf --query_id P01501 --home_species "Apis mellifera" --max_genomes 10

# Apamin in bees  
nextflow run main.nf --query_id P01500 --home_species "Apis mellifera" --max_genomes 10
```

### Primate Immune Genes
```bash
# LY6E (viral restriction factor)
nextflow run main.nf --query_id Q16553 --home_species "Homo sapiens" --max_genomes 5
```

### Ant Venoms
```bash
# Using custom fasta
nextflow run main.nf \
  --gene test_data/tetramorium_query.fasta \
  --home_species "Tetramorium bicarinatum" \
  --max_genomes 5
```

## Citation

```bibtex
@software{synterra2026,
  title = {SynTerra: Phylogenetically-informed syntenic ortholog discovery},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/yourusername/SynTerra}
}
```

## License

MIT License - see [LICENSE](LICENSE)

## Acknowledgments

Built with [Nextflow](https://www.nextflow.io/), [MMseqs2](https://github.com/soedinglab/MMseqs2), [MAFFT](https://mafft.cbrc.jp/alignment/software/), [FastTree](http://www.microbesonline.org/fasttree/), and [Plotly](https://plotly.com/).
