# SynTerra: Mapping the syntenic landscape of divergent genomes.

**Novel synteny-guided gene finding across genomes using iterative phylogenetic search**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

SynTerra is a scientific tool for discovering genes in newly sequenced genomes by exploiting synteny conservation. Unlike traditional homology search, SynTerra:

1. **Uses synteny as evidence**: Searches for conserved genomic neighborhoods, not just sequence similarity
2. **Iterative phylogenetic approach**: Progressively searches through evolutionarily ordered genomes, using closer species as "stepping stones" to find distant orthologs
3. **Data augmentation**: Generates sequence variants to handle divergent genes, frameshifts, and domain shuffling
4. **Handles difficult cases**: Finds pseudogenes, partial genes, and highly divergent orthologs missed by standard BLAST

## Key Innovation

The **iterative phylogenetic search** (Step 3) is what makes SynTerra powerful:

```
Home genome → Close species (90% ID) → Medium species (70% ID) → Distant species (50% ID)
                    ↓ add to DB             ↓ add to DB              ↓
              Direct search would FAIL ───────────────────────────────┘
              But iterative search SUCCEEDS ✓
```

## Quick Start

### Requirements

- **Nextflow** ≥ 22.10.1
- **Conda** or **Docker** (recommended)
- Dependencies: MMseqs2, BLAST+, Prodigal, AUGUSTUS, Plotly

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/SynTerra.git
cd SynTerra

# Option 1: Conda (recommended for development)
conda env create -f environment.yml
conda activate syntenyfinder

# Option 2: Docker (recommended for reproducibility)
docker build -t synterra .

# Option 3: Singularity (for HPC)
singularity build synterra.sif docker://synterra:latest
```

### Running the Pipeline

```bash
# With Conda
nextflow run main.nf -profile conda --gene query.faa ...

# With Docker
nextflow run main.nf -profile docker --gene query.faa ...

# With Singularity (HPC)
nextflow run main.nf -profile singularity --gene query.faa ...
```

### Test Run

```bash
# Run with test data (honeybee melettin gene)
nextflow run main.nf -profile test_melettin,conda

# View results
open results/test_melettin/synteny_plot.html
```

## Usage

### Easy Mode (Auto-fetch related genomes)

```bash
# SynTerra automatically downloads related genomes from NCBI
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome home_genome.fna \
  --home_gff home_genome.gff \
  --mode easy \
  --easy_species "Apis mellifera" \
  --easy_max_genomes 10 \
  --outdir results
```

**Requirements for Easy Mode**:
- NCBI E-utilities: `conda install -c bioconda entrez-direct`
- (Optional) NCBI Datasets CLI: `conda install -c conda-forge ncbi-datasets-cli`

### Pro Mode (Provide your own genomes)

```bash
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome home_genome.fna \
  --home_gff home_genome.gff \
  --target_genomes "genomes/*.fna" \
  --mode pro \
  --outdir results
```

### Input Requirements

| Parameter | Description | Format | Required |
|-----------|-------------|--------|----------|
| `--gene` | Query gene sequence | FASTA (DNA or protein) | ✓ |
| `--home_genome` | Genome where gene is known | FASTA | ✓ |
| `--home_gff` | Annotation for home genome | GFF3 | Optional |
| `--mode` | Pipeline mode | 'easy' or 'pro' | ✓ (default: pro) |
| `--easy_species` | Species name for easy mode | String (e.g., "Apis mellifera") | Required for easy mode |
| `--easy_max_genomes` | Max genomes to fetch | Integer | Optional (default: 10) |
| `--target_genomes` | Target genomes (pro mode) | Glob pattern or list | Required for pro mode |

### Key Parameters

```groovy
// Synteny parameters
--n_flanking_genes 10           // Number of flanking genes (default: 10)
--min_synteny_score 0.6         // Minimum synteny conservation (default: 60%)

// Search sensitivity
--mmseqs_sensitivity 8.5        // MMseqs2 sensitivity (default: 8.5, very high)
--min_gene_identity 30          // Minimum sequence identity % (default: 30%)

// Augmentation (for divergent genes)
--mutation_rate 0.05            // Simulate 5% divergence (default)
--enable_frameshifts true       // Search for frameshifted genes
```

## Output Files

```
results/
├── synteny_plot.html           # Interactive synteny visualization
├── synteny_block.bed           # Home genome synteny block
├── flanking_proteins.faa       # Flanking gene proteins
├── {genome}_regions.bed        # Identified syntenic regions per genome
├── {genome}_gene_candidates.*  # Found genes per genome
└── {genome}_genes_annotated.gff # Final gene structures
```

## Workflow Steps

1. **STAGE_GENOME** - Prepare and validate input genomes
2. **LOCATE_GENE** - Find query in home genome
3. **EXTRACT_FLANKING** - Get conserved genomic environment
4. **PHYLO_SORT** - Order target genomes by distance
5. **ITERATIVE_SEARCH** - Cumulative search across ordered genomes
6. **ASSESS_QUALITY** - Quality assessment (N50, L50) of targets
7. **CLUSTER_REGIONS** - Identify syntenic blocks
8. **AUGMENTED_SEARCH** - High-sensitivity discovery in regions
9. **ANNOTATE_STRUCTURE** - Ab initio gene prediction
10. **HOMOLOGY_SEARCH** - Orthology validation
11. **PLOT_SYNTENY** - Interactive visualization
12. **GENERATE_REPORT** - Final pipeline summary (JSON)

## Scientific Considerations

### When to use SynTerra

✅ Finding genes in newly sequenced genomes  
✅ Distant orthologs (>40% evolutionary divergence)  
✅ Pseudogenes and partial genes  
✅ Gene families with synteny conservation  

### When NOT to use

❌ Recently duplicated genes (synteny unclear)  
❌ Highly mobile genes (transposons)  
❌ Genes with no synteny conservation  
❌ Single-exon genes with no context  

### Interpreting Results

- **Synteny score ≥ 0.7**: High confidence, strong synteny
- **Synteny score 0.5-0.7**: Moderate confidence, some rearrangement
- **Synteny score < 0.5**: Low confidence, verify manually

## Algorithm Details

See [instructions.md](instructions.md) for complete algorithm description and implementation details.

## Citation

If you use SynTerra, please cite:

```
Your Name et al. (2026). SynTerra: Iterative phylogenetic synteny-guided 
gene finding across genomes. Journal, Volume, Pages.
```

## License

MIT License

## Contact

- **Author**: Your Name
- **Email**: your.email@institution.edu
- **Issues**: https://github.com/yourusername/SynTerra/issues

## Acknowledgments

Built with Nextflow, MMseqs2, and Plotly.
