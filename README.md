# SynTerra

**Phylogenetically-informed syntenic ortholog discovery across divergent genomes**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

SynTerra discovers orthologous genes in newly sequenced genomes by exploiting **synteny conservation** — the principle that gene order is preserved across evolution. Unlike traditional BLAST-based approaches, SynTerra:

- **Uses genomic context**: Finds genes by their conserved neighborhood, not just sequence similarity
- **Iterative phylogenetic search**: Uses evolutionarily close species as "stepping stones" to find distant orthologs
- **Exon-aware annotation**: Annotates GOI exons via GFF matching or hit-based splice-site detection
- **High sensitivity**: Combines MMseqs2 with Smith-Waterman alignment for rigorous gene detection
- **Always protein → DNA**: Uses `tblastn` / MMseqs2 `--search-type 2` throughout — never DNA → DNA

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
# Find melittin orthologs across bee species
nextflow run main.nf \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/melittin

# Find LY6E orthologs across primates
nextflow run main.nf \
  --query_id Q16553 \
  --home_species "Homo sapiens" \
  --max_genomes 5 \
  --outdir results/LY6E
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
| `--gene` | Query gene FASTA (alternative to `query_id`) | - |
| `--home_species` | Species name (Easy mode) | - |
| `--max_genomes` | Number of related genomes to analyze | 10 |
| `--mode` | `easy` or `pro` | easy |
| `--outdir` | Output directory | results |

### Search Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--n_flanking_genes` | Flanking genes for synteny context | 10 |
| `--mmseqs_sensitivity` | MMseqs2 sensitivity (1-8.5) | 8.5 |
| `--min_synteny_score` | Minimum synteny conservation (0-1) | 0.6 |
| `--enable_smith_waterman` | Use Smith-Waterman for GOI | true |
| `--sw_min_identity` | Minimum SW alignment identity | 20.0 |

### Pro Mode Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--home_genome` | Path to home genome FASTA | - |
| `--home_gff` | Home genome GFF3 (optional, improves accuracy) | - |
| `--target_genomes` | Glob pattern for targets (e.g., `"genomes/*.fna"`) | - |

## Output

```
results/
├── query/
│   └── P01501.fasta                # Query sequence (fetched from UniProt)
├── home_genome/
│   ├── home_genome.fna             # Home genome (fetched from NCBI)
│   └── home_genome.gff             # Home annotation (if available)
├── downloaded_genomes/
│   └── easy_mode_genomes/          # Related genomes (easy mode)
├── qc/
│   └── genome_qc_summary.json     # Assembly quality metrics (N50, L50)
├── iterative_results/
│   ├── expanded_db.faa             # All discovered orthologs
│   ├── hits/                       # Per-genome search hits
│   └── regions/
│       ├── *.gff                   # Gene predictions per genome
│       ├── *.faa                   # Translated proteins
│       └── *.homology.tsv          # Homology mappings
├── *_synteny_plot.html             # Interactive synteny visualization
└── synterra_report.json            # Pipeline summary
```

## Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 1: Gene Localization in Home Genome                      │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ LOCATE   │  │ ANNOTATE      │  │ SPLIT     EXTRACT       │  │
│  │ GENE     │→ │ GOI (exons)   │→ │ LOCI  →   FLANKING     │  │
│  └──────────┘  └───────────────┘  └─────────────────────────┘  │
│        ↓                                      ↓                 │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ PREPARE_INITIAL_DB (flanking genes + GOI exons)          │   │
│  └──────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 2: Iterative Phylogenetic Search                         │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ PHYLO    │  │ ITERATIVE     │  │ CLUSTER                 │  │
│  │ SORT     │→ │ SEARCH        │→ │ REGIONS                 │  │
│  └──────────┘  │ (MMseqs2+SW)  │  └─────────────────────────┘  │
│                └───────────────┘                                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 3: Phylogenetics & Visualization                         │
│  ┌──────────┐  ┌───────────────┐  ┌─────────────────────────┐  │
│  │ COMPUTE  │  │ PLOT          │  │ GENERATE                │  │
│  │ TREE     │→ │ SYNTENY       │→ │ REPORT                  │  │
│  └──────────┘  └───────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

### Pipeline Steps

1. **LOCATE_GENE** — `tblastn` + MMseqs2 protein→DNA search to find the GOI in the home genome
2. **ANNOTATE_GOI** — Annotate individual exons of the GOI:
   - *With GFF*: Match GOI to annotated gene by name (UniProt lookup), extract CDS/exon proteins
   - *Without GFF*: Detect exon boundaries from tblastn hits (splice sites, start/stop codons)
3. **SPLIT_LOCI** — Cluster hits into distinct loci (handles multi-copy genes)
4. **EXTRACT_FLANKING** — Extract flanking genes from GFF or Prodigal prediction (exon-level)
5. **PREPARE_INITIAL_DB** — Build search database: flanking proteins + GOI exons (fallback: arbitrary fragments)
6. **PHYLO_SORT** — Order target genomes by phylogenetic distance (NCBI taxonomy or alphabetical)
7. **ITERATIVE_SEARCH** — Wavefront-parallel search across genomes, expanding the database with each discovery. Combines MMseqs2, Smith-Waterman, and ORF-based annotation.
8. **CLUSTER_REGIONS** — Score and rank syntenic regions by flanking gene conservation
9. **COMPUTE_TREE** — MAFFT alignment + FastTree phylogeny of discovered GOI orthologs
10. **PLOT_SYNTENY** — Interactive Plotly visualization with gene arrows, homology links, and tree coloring
11. **GENERATE_REPORT** — JSON summary of all discoveries and quality metrics

## Running on HPC (SLURM)

```bash
# Submit with SLURM + Singularity
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
| `hpc_conda` | SLURM + Conda/Mamba |
| `conda` | Local + Conda |
| `singularity` | Local + Singularity |
| `docker` | Local + Docker |

## Examples

### Insect Venom Peptides
```bash
# Melittin in bees
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
# Using custom FASTA
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

MIT License — see [LICENSE](LICENSE)

## Acknowledgments

Built with [Nextflow](https://www.nextflow.io/), [MMseqs2](https://github.com/soedinglab/MMseqs2), [MAFFT](https://mafft.cbrc.jp/alignment/software/), [FastTree](http://www.microbesonline.org/fasttree/), and [Plotly](https://plotly.com/).
