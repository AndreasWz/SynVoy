# SynTerra Usage Guide

## Table of Contents

1. [Quick Start](#quick-start)
2. [Input Preparation](#input-preparation)
3. [Running the Pipeline](#running-the-pipeline)
4. [Understanding Parameters](#understanding-parameters)
5. [Interpreting Results](#interpreting-results)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Usage](#advanced-usage)

---

## Quick Start

### Minimal Example (Easy Mode)

```bash
nextflow run main.nf \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --outdir results/melittin
```

This will:
1. Fetch query protein from UniProt
2. Download the reference genome and GFF from NCBI
3. Download related genomes from the same taxonomic family
4. Locate the GOI in the home genome (tblastn + MMseqs2)
5. Annotate GOI exons (GFF match or hit-based splice-site detection)
6. Extract flanking genes (synteny block)
7. Iteratively search target genomes in phylogenetic order
8. Generate an interactive synteny plot and report

### Minimal Example (Pro Mode)

```bash
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome home_genome.fna \
  --target_genomes "target_genomes/*.fna" \
  --mode pro
```

---

## Input Preparation

### 1. Query Gene

You can provide the query in two ways:

**Option A — UniProt ID** (recommended for Easy mode):
```bash
--query_id P01501     # Melittin
--query_id Q16553     # LY6E
```
SynTerra fetches the protein sequence from the UniProt REST API automatically.

**Option B — FASTA file**:
```bash
--gene my_query.fasta
```

The query should be a **protein** sequence. SynTerra always searches protein → DNA (never DNA → DNA). If you provide a DNA sequence, it will be auto-detected and handled, but protein input is preferred.

**Requirements**:
- Single sequence recommended (multi-exon gene in one file is fine)
- Minimum length: ~15 amino acids
- Standard FASTA format

### 2. Home Genome (`--home_genome`)

The genome of the species where your query gene originates.

**In Easy mode**: Automatically fetched from NCBI based on `--home_species`.

**In Pro mode**:
```bash
--home_genome my_species.genome.fna
```

**Requirements**:
- Assembled genome (scaffolds or chromosomes)
- FASTA format (`.fna`, `.fasta`, `.fa`)
- Can be draft assembly — SynTerra handles fragmented assemblies

### 3. Home Annotation (`--home_gff`) — Optional but recommended

**In Easy mode**: Fetched automatically with the genome (if NCBI provides one).

**In Pro mode**:
```bash
--home_gff my_species.genes.gff
```

**Why it matters**: With a GFF, SynTerra can:
- Match the GOI to an annotated gene by name and extract proper CDS/exon sequences
- Extract accurate flanking genes with correct protein translations
- Without GFF, SynTerra falls back to Prodigal gene prediction and hit-based exon annotation

### 4. Target Genomes (`--target_genomes`)

**In Easy mode**: Automatically downloaded from NCBI — related species from the same genus/family.

**In Pro mode**:
```bash
--target_genomes "genomes/*.fna"
```

**Tips**:
- 5–15 genomes works best
- Mix close and distant relatives for best phylogenetic coverage
- No annotation required for targets — SynTerra annotates genes during search

---

## Running the Pipeline

### Easy Mode (Recommended for most users)

```bash
# Everything fetched automatically
nextflow run main.nf \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/melittin

# With a local query file instead of UniProt ID
nextflow run main.nf \
  --gene my_query.fasta \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results
```

**Requirements**:
- Internet connection (for NCBI/UniProt)
- NCBI datasets CLI: `conda install -c conda-forge ncbi-datasets-cli`
- NCBI E-utilities: `conda install -c bioconda entrez-direct`

### Pro Mode (Custom genomes)

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome home.fna \
  --home_gff home.gff \
  --target_genomes "targets/*.fna" \
  --mode pro \
  --outdir results
```

### Using Test Data

```bash
# Melittin test
nextflow run main.nf -profile test_melettin

# Tetramorium test
nextflow run main.nf -profile test_tetramorium

# Generic test
nextflow run main.nf -profile test
```

### Resume a Failed Run

```bash
nextflow run main.nf -resume \
  --query_id P01501 \
  --home_species "Apis mellifera"
```

### HPC / SLURM

```bash
# SLURM + Singularity
nextflow run main.nf \
  -profile hpc_singularity \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/my_run \
  -work-dir $SCRATCH/work

# SLURM + Conda
nextflow run main.nf \
  -profile hpc_conda \
  --query_id P01501 \
  --home_species "Apis mellifera"
```

---

## Understanding Parameters

### Mode Selection

| Mode | Use Case | Required Parameters |
|------|----------|---------------------|
| `easy` (default) | Auto-fetch from NCBI | `--query_id` or `--gene`, `--home_species` |
| `pro` | Custom genome files | `--gene`, `--home_genome`, `--target_genomes` |

### Easy Mode Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--query_id` | - | UniProt accession (e.g., `P01501`) |
| `--gene` | - | OR: Path to query protein FASTA |
| `--home_species` | Required | Species name (e.g., `"Apis mellifera"`) |
| `--max_genomes` | 10 | Number of related genomes to download |

### Pro Mode Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gene` | Required | Path to query gene FASTA |
| `--home_genome` | Required | Path to home genome FASTA |
| `--home_gff` | Optional | Path to home genome GFF3 annotation |
| `--target_genomes` | Required | Glob pattern for targets (e.g., `"genomes/*.fna"`) |

### Synteny Parameters

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `--n_flanking_genes` | 10 | Number of flanking genes per side | Increase to 15–20 for large genomes |
| `--min_synteny_score` | 0.6 | Minimum synteny conservation (0–1) | Lower to 0.4–0.5 for distant species |
| `--cluster_distance` | 50000 | Max distance to cluster hits (bp) | Increase for gene-poor genomes |
| `--prefer_large_genes` | true | Prefer longer flanking genes | Set false for compact genomes |
| `--min_flanking_size` | 500 | Minimum flanking gene size (bp) | Lower for small gene-dense genomes |

### Search Sensitivity

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `--mmseqs_sensitivity` | 8.5 | MMseqs2 sensitivity (1–9) | Lower to 7.5 for faster search |
| `--min_gene_identity` | 30 | Min % identity for gene hits | Raise to 40–50 for close species |
| `--enable_smith_waterman` | true | Use Smith-Waterman for GOI | Disable for speed |
| `--sw_min_identity` | 20.0 | Min SW alignment identity (%) | Lower for very divergent genes |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | results | Output directory |
| `--keep_intermediate` | false | Keep work dir files for debugging |

---

## Interpreting Results

### Output Files

```
results/
├── query/
│   └── P01501.fasta                # Query sequence
├── home_genome/                    # Home genome + GFF (easy mode)
├── downloaded_genomes/             # Target genomes (easy mode)
├── qc/
│   └── genome_qc_summary.json     # Assembly quality (N50, L50, contig count)
├── iterative_results/
│   ├── expanded_db.faa             # All discovered orthologs (growing DB)
│   ├── hits/                       # Per-genome MMseqs2 hit files
│   └── regions/
│       ├── {genome}.gff            # Gene annotations per genome
│       ├── {genome}.faa            # Translated proteins per genome
│       └── {genome}.homology.tsv   # Gene-to-home mappings
├── *_synteny_plot.html             # Interactive synteny visualization
└── synterra_report.json            # Pipeline summary
```

### Interactive Synteny Plot (`*_synteny_plot.html`)

Open in a web browser. Each horizontal track represents one genome.

- **Colored blocks** = genes (arrows show strand direction)
- **Red blocks** = your query gene / GOI
- **Gray/blue blocks** = flanking genes
- **Lines** connect orthologous genes across genomes
- **Hover** over genes for name, coordinates, and synteny score

**What to look for**:
- **Strong synteny**: Flanking genes in same order across genomes → high confidence
- **Gene present**: GOI (red) appears in the syntenic region → ortholog found
- **Rearrangements**: Breaks in gene order → genomic changes (inversions, translocations)

### Synteny Scores

| Score | Interpretation | Confidence |
|-------|----------------|------------|
| ≥ 0.7 | Excellent synteny | High |
| 0.5–0.7 | Good synteny (some rearrangement) | Medium |
| 0.3–0.5 | Weak synteny (major rearrangement) | Low |
| < 0.3 | Very weak / no synteny | Very Low |

**Guidelines**:
- Score ≥ 0.6: Trust the result
- Score 0.4–0.6: Manually verify with BLAST
- Score < 0.4: Likely false positive

---

## Troubleshooting

### Common Issues

#### "No hits found in home genome"

- Check that the query sequence is from the home species
- Try protein input instead of DNA
- Ensure query is long enough (≥15 aa)

#### "No flanking genes found"

- Provide `--home_gff` for accurate gene extraction
- Check that the gene is in a gene-rich region
- Try increasing `--n_flanking_genes`

#### "Low synteny scores everywhere"

- Lower `--min_synteny_score` to 0.4
- Increase `--n_flanking_genes` to 15–20
- The gene region may not be syntenic (common for lineage-specific genes)

#### "Pipeline runs forever"

- Reduce `--max_genomes` or `--mmseqs_sensitivity`
- Use `-resume` to avoid re-running completed steps
- Run on HPC with `-profile hpc_singularity`

#### "Multiple gene copies found"

- Expected for duplicated genes
- The copy with the highest synteny score is the ortholog
- All copies are reported in output

### Logs

```bash
# Full pipeline log
cat .nextflow.log

# See what ran
nextflow log
```

---

## Advanced Usage

### Phylogenetic Sorting with NCBI Taxonomy

To sort target genomes by phylogenetic distance (improves iterative search):

```bash
# 1. Download NCBI taxonomy dump
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p taxdb && tar -xzf taxdump.tar.gz -C taxdb/

# 2. Set environment variable
export TAXDB=$PWD/taxdb

# 3. Run — genomes will be searched in phylogenetic order
nextflow run main.nf --query_id P01501 --home_species "Apis mellifera"
```

Without `TAXDB`, genomes are searched in alphabetical order (still works, just less optimal).

### Very Divergent Genes

For genes with >60% divergence:

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome genome.fna \
  --target_genomes "targets/*.fna" \
  --mode pro \
  --min_synteny_score 0.4 \
  --min_gene_identity 25 \
  --sw_min_identity 15.0 \
  --mmseqs_sensitivity 8.5
```

### Gene-Poor Genomes

For genomes with sparse genes (e.g., large plant genomes):

```bash
nextflow run main.nf \
  --gene gene.fasta \
  --home_genome genome.fna \
  --target_genomes "targets/*.fna" \
  --mode pro \
  --n_flanking_genes 5 \
  --cluster_distance 100000 \
  --prefer_large_genes true
```

### Profiles Reference

| Profile | Executor | Container | Use Case |
|---------|----------|-----------|----------|
| `standard` | local | Conda | Default local |
| `conda` | local | Conda | Explicit Conda |
| `docker` | local | Docker | Docker users |
| `singularity` | local | Singularity | Singularity users |
| `slurm` | SLURM | (none) | Basic SLURM |
| `hpc_singularity` | SLURM | Singularity | HPC (recommended) |
| `hpc_conda` | SLURM | Conda/Mamba | HPC + Conda |
| `test` | local | Conda | Quick test |
| `test_melettin` | local | Conda | Melittin test |
| `test_tetramorium` | local | Conda | Tetramorium test |

---

## Support

- **Full parameter list**: See `nextflow.config`
- **Pipeline details**: See `README.md`
- **Issues**: https://github.com/yourusername/SynTerra/issues
