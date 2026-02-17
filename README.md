# SynTerra

SynTerra is a Nextflow pipeline for synteny-guided ortholog discovery in related genomes.
It is designed to work when target genomes are:
- well annotated,
- partially annotated,
- custom annotated,
- or completely unannotated.

## What SynTerra Does

1. Resolves your GOI query (UniProt/NCBI ID or FASTA).
2. Locates GOI loci in a home genome (MMseqs + BLAST).
3. Annotates GOI exon structure (from GFF when available, otherwise from protein-to-genome evidence).
4. Extracts flanking genes around each home locus.
5. Iteratively searches phylogenetically ordered target genomes.
6. Annotates GOI and flanking models in target regions.
7. Produces synteny plots, trees, and summary outputs.

## Key Logic Choices

- Search operates in protein space against genomic DNA (translated search), then refines exon/intron structures.
- GOI and flanking handling are separate by design.
- Flanking model collapsing is per-locus, not global:
  - one best model per parent ID per genomic locus,
  - multiple loci for the same parent ID are kept if they are distinct loci.
- Iterative DB expansion uses GOI-derived models; flanking models are kept for context/annotation and plotting.
- Parent IDs are handled generically (not tied to `XP_` accessions).

## Input Robustness

SynTerra is built to handle:
- NCBI GFF/GFF3,
- Ensembl-like GFF/GTF attribute styles,
- custom GFFs,
- missing GFF (`NO_GFF` fallback path with Prodigal + borrowed annotations).

## Setup (Choose One Backend)

You can run SynTerra with Conda (recommended), Docker, or Singularity/Apptainer.

### Option A: Conda/Mamba (recommended for most users)

1. Install Miniconda or Mambaforge.
2. Create and activate the environment.
3. Run SynTerra.

```bash
cd /path/to/SynTerra
conda env create -f environment.yml
conda activate syntenyfinder
nextflow -version
```

If `nextflow` is not on your `PATH`, use the bundled launcher in this repo:

```bash
./nextflow -version
```

Use the modern live terminal UI launcher:

```bash
./synterra --help
```

### Option B: Docker

Prerequisites:
- Docker installed and running.
- Nextflow installed on host.

Run with Docker profile:

```bash
nextflow run main.nf -profile docker --gene P01501 --mode easy --outdir results
```

### Option C: Singularity/Apptainer

Prerequisites:
- Singularity or Apptainer installed.
- Nextflow installed on host.

Recommended cache setup:

```bash
mkdir -p "$HOME/.singularity/cache"
export NXF_SINGULARITY_CACHEDIR="$HOME/.singularity/cache"
```

Run with Singularity profile:

```bash
nextflow run main.nf -profile singularity --gene P01501 --mode easy --outdir results
```

Optional quick smoke test (bundled test data):

```bash
nextflow run main.nf -profile test
```

## Quick Start

Easy mode (ID-based, automatic genome retrieval):

```bash
./synterra \
  --gene P01501 \
  --mode easy \
  --outdir results
```

Note: easy mode requires internet access to fetch genomes/metadata.
Assembly selection prefers reference/representative genomes, then ranks fallback assemblies by quality.  
If only low-quality assemblies exist, default policy is `--bad_quality_policy ask` with `--bad_quality_timeout 300` seconds (auto-NO on timeout).

Pro mode (user-provided files):

```bash
./synterra \
  --mode pro \
  --gene input/query.fasta \
  --home_genome input/home.fna \
  --home_gff input/home.gff \
  --target_genomes "input/targets/*.fna" \
  --outdir results
```

Resume an interrupted run (same command and same `--outdir`):

```bash
./synterra --gene P01501 --mode easy --outdir results -resume
```

Modern single-screen terminal UI is provided by `./synterra` and hides repetitive `executor` table reprints.
If you want classic Nextflow console output instead:

```bash
nextflow run main.nf --gene P01501 --mode easy --outdir results -resume
```

## Main Outputs

Under `--outdir`:
- `*_synteny_plot.html`: interactive synteny plots
- `*_tree.nwk`: GOI tree per locus
- `synterra_report.json`: run summary
- `regions/*.regions.bed`: clustered candidate regions
- `intermediate/`: phase-level intermediate artifacts

## Documentation

- Detailed setup, runtime profiles, and parameters: `USAGE.md`
- Detailed algorithm and data flow: `PIPELINE_DETAILED.md`
