# SynTerra Usage Guide

## Table of Contents

1. [Quick Start](#quick-start)
2. [Input Preparation](#input-preparation)
3. [Running the Pipeline](#running-the-pipeline)
4. [Understanding Parameters](#understanding-parameters)
5. [Interpreting Results](#interpreting-results)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Usage](#advanced-usage)

## Quick Start

### Minimal Example

```bash
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome home_genome.fna \
  --target_genomes "target_genomes/*.fna"
```

This will:
1. Find your gene in the home genome
2. Extract flanking genes (synteny block)
3. Search for this synteny block in target genomes
4. Locate and annotate the gene in each target genome
5. Generate an interactive synteny plot

## Input Preparation

### 1. Query Gene (`--gene`)

**Format**: FASTA file (DNA or protein)

**Example**:
```fasta
>my_gene_of_interest
ATGGCTAGCTAGCTAGCTAGCTAG...
```

**Requirements**:
- Single sequence or multi-exon gene fragments
- DNA or amino acid sequence
- Minimum length: ~50 bp (DNA) or ~15 aa (protein)

**Tips**:
- If gene has multiple exons, provide all exons in one file
- SynTerra will detect them if they hit close together

### 2. Home Genome (`--home_genome`)

**Format**: FASTA file (genomic sequences)

**Example**:
```bash
--home_genome my_species.genome.fna
```

**Requirements**:
- Assembled genome (scaffolds/chromosomes)
- Can be draft assembly
- FASTA format

**Tips**:
- Use the highest quality assembly available
- Compress with gzip is OK (`.fna.gz`)

### 3. Home Annotation (`--home_gff`) - Optional

**Format**: GFF3 file

**Example**:
```bash
--home_gff my_species.genes.gff
```

**Requirements**:
- GFF3 format with gene features
- Must match chromosome names in genome file

**If not provided**:
- SynTerra will use ab initio gene prediction (Prodigal)
- Still works, but less accurate flanking gene identification

### 4. Target Genomes (`--target_genomes`)

**Format**: Multiple FASTA files

**Example**:
```bash
# Glob pattern
--target_genomes "genomes/*.fna"

# Or explicit list (in nextflow)
--target_genomes "genome1.fna,genome2.fna,genome3.fna"
```

**Requirements**:
- Assembled genomes (can be draft)
- FASTA format
- No annotation required

**Tips**:
- Include 3-10 genomes for best results
- Mix close and distant relatives
- More genomes = better phylogenetic coverage

## Running the Pipeline

### Easy Mode (Auto-fetch genomes)

**Best for**: First-time users, exploring a gene family, quick analyses

```bash
# Automatically download related genomes from NCBI
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome home.fna \
  --home_gff home.gff \
  --mode easy \
  --easy_species "Apis mellifera" \
  --easy_max_genomes 10 \
  --outdir results
```

**How it works**:
1. Searches NCBI for assemblies related to your species
2. Downloads top N related genomes (default: 10)
3. Runs synteny analysis on downloaded genomes
4. Cleans up temporary files

**Requirements**:
- NCBI E-utilities: `conda install -c bioconda entrez-direct`
- Internet connection

**Tips**:
- Use scientific name in quotes: `"Drosophila melanogaster"`
- Start with 5-10 genomes to test
- Downloaded genomes saved in `results/downloaded_genomes/`

### Pro Mode (Your own genomes)

**Best for**: Custom genome sets, offline analysis, specific comparisons

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome home.fna \
  --home_gff home.gff \
  --target_genomes "targets/*.fna" \
  --mode pro \
  --outdir results
```

### With Custom Parameters

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome home.fna \
  --target_genomes "targets/*.fna" \
  --n_flanking_genes 15 \
  --min_synteny_score 0.5 \
  --mmseqs_sensitivity 8.5 \
  --outdir my_results
```

### Using Test Data

```bash
# Quick test with provided test data
nextflow run main.nf -profile test

# Results will be in test_results/
```

### Resume Failed Run

```bash
# If pipeline failed or was interrupted
nextflow run main.nf -resume \
  --gene gene.fasta \
  --home_genome home.fna \
  --target_genomes "targets/*.fna"
```

## Understanding Parameters

### Core Parameters

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `--gene` | Required | Query gene sequence | - |
| `--home_genome` | Required | Home genome FASTA | - |
| `--home_gff` | Optional | Home genome annotation | Provide if available |
| `--mode` | 'pro' | Pipeline mode ('easy' or 'pro') | Use 'easy' for auto-fetch |
| `--easy_species` | Only for easy mode | Species name for NCBI | Required in easy mode |
| `--easy_max_genomes` | 10 | Max genomes to fetch | Adjust based on needs |
| `--target_genomes` | Only for pro mode | Target genome files | Required in pro mode |

### Synteny Parameters

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `--n_flanking_genes` | 10 | Number of flanking genes | Increase for large genomes (15-20) |
| `--min_synteny_score` | 0.6 | Minimum synteny conservation (0-1) | Lower (0.4-0.5) for distant species |
| `--cluster_distance` | 50000 | Max distance to cluster hits (bp) | Increase for gene-poor genomes |

### Search Sensitivity

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `--mmseqs_sensitivity` | 8.5 | MMseqs2 sensitivity (1-9) | Lower (7.5) for fast search |
| `--min_gene_identity` | 30 | Min % identity for gene hits | Increase (40-50) for close species |
| `--mutation_rate` | 0.05 | Simulation divergence rate | Increase (0.10) for very divergent |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | results | Output directory |

## Interpreting Results

### Output Files

```
results/
├── synteny_plot.html              # 🎨 Main result: Interactive plot
├── synterra_report.json           # 📊 Final pipeline summary & discovery stats
├── qc/
│   └── genome_qc_summary.json     # 🧪 Genome quality metrics (N50, L50)
├── synteny_block.bed              # Flanking genes in home genome
├── flanking_proteins.faa          # Protein sequences of flanking genes
├── {genome}_regions.bed           # Identified syntenic regions
├── {genome}_gene_candidates.bed   # Gene candidates per genome
├── {genome}_genes_annotated.gff   # Final annotated genes
└── {genome}_proteins.faa          # Translated proteins
```

### Interactive Synteny Plot (`synteny_plot.html`)

**How to use**:
1. Open in web browser
2. Each horizontal track = one genome
3. Colored blocks = genes
4. **Red blocks** = your query gene (target)
5. Gray/blue blocks = flanking genes
6. Lines connect orthologous genes

**Hover** over genes to see:
- Gene name
- Coordinates
- Synteny score

**What to look for**:
- **Strong synteny**: Flanking genes in same order across genomes
- **Gene location**: Query gene (red) in syntenic region
- **Rearrangements**: Breaks in synteny = genomic changes

### Synteny Scores

| Score | Interpretation | Confidence |
|-------|----------------|------------|
| ≥ 0.7 | Excellent synteny | High ✅ |
| 0.5-0.7 | Good synteny (some rearrangement) | Medium ⚠️ |
| 0.3-0.5 | Weak synteny (major rearrangement) | Low ⚠️ |
| < 0.3 | Very weak/no synteny | Very Low ❌ |

**Action**:
- Score ≥ 0.6: Trust the result
- Score 0.4-0.6: Manually verify with BLAST
- Score < 0.4: May be false positive

## Troubleshooting

### Common Issues

#### 1. "No hits found in home genome"

**Cause**: Gene sequence doesn't match home genome

**Solutions**:
- Check that gene sequence is from this species
- Try both DNA and protein versions
- Check for typos in gene sequence

#### 2. "No flanking genes found"

**Cause**: No annotation provided or region is gene-poor

**Solutions**:
- Provide `--home_gff` annotation file
- Check that gene is in gene-rich region
- Try with ab initio prediction (don't provide GFF)

#### 3. "Low synteny scores everywhere"

**Cause**: Gene region not conserved, or wrong parameters

**Solutions**:
- Lower `--min_synteny_score` to 0.4-0.5
- Increase `--n_flanking_genes` to 15-20
- Check if gene is in conserved region (core genes work best)

#### 4. "Pipeline runs forever"

**Cause**: Too many genomes or very large genomes

**Solutions**:
- Use `-resume` to continue
- Reduce `--mmseqs_sensitivity` to 7.5
- Test with 2-3 genomes first
- Use HPC cluster with `-profile cluster`

#### 5. "Multiple gene copies found"

**Cause**: Gene duplication (paralogs)

**Solutions**:
- This is expected for multi-copy genes
- Check synteny score to find ortholog (highest score)
- All copies are reported in output

### Getting Help

```bash
# Show all parameters
nextflow run main.nf --help

# Check pipeline info
nextflow info main.nf

# See detailed logs
cat .nextflow.log
```

## Advanced Usage

### Custom Phylogenetic Sorting

To enable phylogenetic ordering of genomes:

```bash
# 1. Download NCBI taxonomy
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C taxdb/

# 2. Set environment variable
export TAXDB=$PWD/taxdb

# 3. Run pipeline
nextflow run main.nf --gene gene.fasta ...
```

### Custom Species for Gene Prediction

```bash
# If Augustus supports your species
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome genome.fna \
  --target_genomes "targets/*.fna" \
  --augustus_species "drosophila"  # or human, arabidopsis, etc.
```

### Very Divergent Genes

For very divergent genes (>60% divergence):

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome genome.fna \
  --target_genomes "targets/*.fna" \
  --min_synteny_score 0.4 \
  --min_gene_identity 25 \
  --mutation_rate 0.10 \
  --mmseqs_sensitivity 9.0
```

### Gene-Poor Genomes

For genomes with sparse genes (e.g., large plant genomes):

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome genome.fna \
  --target_genomes "targets/*.fna" \
  --n_flanking_genes 5 \
  --cluster_distance 100000 \
  --prefer_large_genes true
```

### Keeping Intermediate Files

For debugging or analysis:

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome genome.fna \
  --target_genomes "targets/*.fna" \
  --keep_intermediate true
```

## Citation

If you use SynTerra in your research, please cite:

```
[Your citation here]
```

## Support

- **Documentation**: See README.md and instructions.md
- **Issues**: https://github.com/yourusername/SynTerra/issues
- **Email**: your.email@institution.edu
