# SynTerra

SynTerra is a Nextflow pipeline for synteny-guided ortholog discovery across related genomes.
It supports annotated genomes, partially annotated genomes, and no-annotation fallback paths.

## What SynTerra Does

1. Resolves GOI input from UniProt/NCBI ID or FASTA.
2. Retrieves or loads the home genome.
3. Locates GOI locus/loci in the home genome (handles home-genome tandem duplications natively).
4. Builds GOI exon-aware representation.
5. Extracts flanking genes per locus (strictly enforcing CDS-containing models).
6. Iteratively searches target genomes in phylogenetic order.
7. Annotates GOI + flanking context in target loci.
8. Performs **Unconstrained Chromosome-Scale Mapping**, gracefully handling distant translocations without losing block data.
9. Generates region clusters, synteny plots, trees, and summary reports.

## Current Runtime Model

- Wrapper command `./synterra` runs `nextflow run main.nf` with a compact live console UI.
- `./synterra --raw ...` disables the custom UI and prints raw Nextflow output.
- Raw logs are always stored in `.synterra_logs/run_YYYYMMDD_HHMMSS.log`.

## Setup

### Conda (recommended)

```bash
cd /path/to/SynTerra
conda env create -f environment.yml
conda activate syntenyfinder
nextflow -version
./synterra --help
```

If `nextflow` is not on your `PATH`, use `./nextflow`.

### Docker

```bash
docker build -t synterra-local:latest .
./synterra -profile docker --gene P01501 --mode easy --outdir results
```

To use another image:

```bash
./synterra -profile docker --docker_container your/image:tag --gene P01501 --mode easy --outdir results
```

### Singularity / Apptainer

```bash
mkdir -p "$HOME/.singularity/cache"
export NXF_SINGULARITY_CACHEDIR="$HOME/.singularity/cache"
./synterra -profile singularity --gene P01501 --mode easy --outdir results
```

## Quick Start

### Easy mode (automatic genome retrieval)

```bash
./synterra \
  --mode easy \
  --gene P01501 \
  --outdir results
```

### Pro mode (local genomes)

```bash
./synterra \
  --mode pro \
  --gene input/query.fasta \
  --home_genome input/home.fna \
  --home_gff input/home.gff \
  --target_genomes "input/targets/*.fna" \
  --outdir results
```

### Resume a stopped run

```bash
./synterra --mode easy --gene P01501 --outdir results -resume
```

## Practical 3FTx Example (3 snake species)

This pattern is tuned for desktop stability and lower crash risk:

```bash
conda activate syntenyfinder

NXF_OPTS='-Xms512m -Xmx2g' ./synterra \
  -profile docker \
  --mode easy \
  --gene P60615 \
  --home_species "Naja naja" \
  --target_species "Ophiophagus hannah,Bungarus multicinctus" \
  --max_genomes 2 \
  --n_flanking_genes 30 \
  --bad_quality_policy keep \
  --iterative_search_cpus 2 \
  --iterative_search_memory '6 GB' \
  --iterative_search_max_forks 1 \
  --mmseqs_split_memory_limit 2G \
  --mmseqs_verbosity 0 \
  --iterative_quiet_subtools true \
  --outdir results_3snake_3ftx
```

## Input Rules

- `--gene` accepts UniProt IDs (for example `P60615`), NCBI protein accessions, or local FASTA paths.
- `--query_id` is legacy and still supported.
- In easy mode, `--home_species` is optional when GOI input is a resolvable ID.
- In easy mode, when `--gene` is a local FASTA file, provide `--home_species`.
- If a FASTA path does not exist, input resolution fails early (`RESOLVE_GENE_INPUT` or input validation).

## Annotation Behavior

- If home GFF is available, SynTerra uses annotation-based extraction first.
- If no usable home GFF exists, SynTerra falls back to local prediction (Prodigal/borrowed annotations path).
- Flanking gene extraction strictly filters for valid CDS (Coding Sequence) parts, bypassing purely non-coding entries (e.g., lncRNAs) to guarantee exactly `N` fully protein-coding flanking genes.
- Flanking genes keep stable IDs and now also carry display labels derived from annotation names when available.
- Iterative expansion is GOI-driven; flanking models are retained for context and plotted with Unconstrained Chromosome-Scale mapping logic.

## Main Outputs

Inside `--outdir`:

- `synterra_report.json`: run summary
- `*_synteny_plot.html`: interactive synteny plots
- `*_tree.nwk`: GOI tree per locus
- `regions/*.regions.bed`: clustered regions
- `downloaded_genomes/easy_mode_genomes/assembly_quality.tsv`: easy-mode assembly quality report
- `intermediate/`: phase-level intermediate files

## Troubleshooting

### Long runtime in fetch phase

- Easy mode downloads large genomes and metadata; this can take hours on slow links.
- Check live raw logs in `.synterra_logs/`.
- Use explicit `--target_species` and small `--max_genomes` to reduce search/download scope.

### `RESOLVE_GENE_INPUT` failed

- Verify `--gene` path exists when passing a local FASTA.
- For ID input, verify the ID is valid and internet access is available.
- If using local FASTA in easy mode, provide `--home_species`.

### Desktop crashes / resource pressure

- Use Docker profile.
- Constrain Java with `NXF_OPTS='-Xms512m -Xmx2g'`.
- Keep `--iterative_search_cpus 2` and `--iterative_search_max_forks 1` for workstation runs.
- Lower `--mmseqs_split_memory_limit` if needed.

### Clean old runs

```bash
rm -rf results_* work
```

## Documentation

- Full setup, profiles, complete parameter reference, and recipes: `USAGE.md`
- Architecture and deep details: `PIPELINE_DETAILED.md`
