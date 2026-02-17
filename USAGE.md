# SynTerra Usage

## 1) Requirements

- Linux/macOS, or Windows with WSL2
- Nextflow `>=22.10.1`
- One runtime backend:
  - Conda/Mamba (recommended for most users)
  - Docker
  - Singularity or Apptainer

SynTerra profiles available in this repo:
- `standard` (default local profile): local executor + Conda (`environment.yml`)
- `conda`: explicit Conda profile
- `docker`: run all processes in container
- `singularity`: run all processes in Singularity/Apptainer container
- `hpc_singularity`: Slurm + Singularity
- `hpc_conda`: Slurm + Conda

## 2) Setup

### Option A: Conda/Mamba (recommended)

1. Create environment:

```bash
conda env create -f environment.yml
```

2. Activate:

```bash
conda activate syntenyfinder
```

3. Verify:

```bash
nextflow -version
python --version
```

If `nextflow` is not on your `PATH`, use the bundled launcher:

```bash
./nextflow -version
```

Optional helper script:

```bash
bash install.sh
```

Modern live logging launcher (recommended):

```bash
./synterra --help
```

### Option B: Docker

Prerequisites:
- Docker installed and running
- Nextflow installed on host

Run with:

```bash
nextflow run main.nf -profile docker --gene P01501 --mode easy --outdir results
```

### Option C: Singularity/Apptainer

Prerequisites:
- Singularity or Apptainer installed
- Nextflow installed on host

Recommended cache setup:

```bash
mkdir -p "$HOME/.singularity/cache"
export NXF_SINGULARITY_CACHEDIR="$HOME/.singularity/cache"
```

Run with:

```bash
nextflow run main.nf -profile singularity --gene P01501 --mode easy --outdir results
```

### HPC profiles (Slurm)

Conda-based Slurm:

```bash
nextflow run main.nf -profile hpc_conda --gene P01501 --mode easy --outdir results
```

Singularity-based Slurm:

```bash
nextflow run main.nf -profile hpc_singularity --gene P01501 --mode easy --outdir results
```

### Optional smoke test (bundled local test data)

```bash
nextflow run main.nf -profile test
```

## 3) Quick Run Modes

### Easy Mode

Use an ID or FASTA, and let SynTerra fetch home/related genomes.

```bash
./synterra \
  --gene P01501 \
  --mode easy \
  --outdir results
```

Easy mode requires internet access for genome/metadata retrieval.

Genome selection policy in easy mode:
- Prefer RefSeq `reference genome` / `representative genome` if available.
- Otherwise rank assemblies by quality (`--assembly-ranking`):
  - `hybrid` (default): reference/assembly-level + contiguity + N50/N80
  - `counts`: prioritize fewer chromosomes/scaffolds/contigs
  - `nstats`: prioritize higher N50/N80
- If only low-quality assemblies are available, policy is controlled by `--bad_quality_policy`:
  - `ask` (default): prompt and wait up to `--bad_quality_timeout` seconds (default `300`), then default to NO
  - `drop`: exclude low-quality assemblies
  - `keep`: keep them anyway

### Pro Mode

Provide your own home genome and target genomes.

```bash
./synterra \
  --mode pro \
  --gene input/query.fasta \
  --home_genome input/home.fna \
  --home_gff input/home.gff \
  --target_genomes "input/targets/*.fna" \
  --outdir results
```

### Resume

Use the same command and same `--outdir`, then add `-resume`. If you change parameters (or output folder), Nextflow may rerun download steps.

```bash
./synterra --gene P01501 --mode easy --outdir results -resume
```

### Cleaner terminal output

Use `./synterra` for a controlled single-line live UI (spinner + in-place status updates, no repeated `executor > local` tables).
If you want classic Nextflow console output:

```bash
nextflow run main.nf --gene P01501 --mode easy --outdir results -resume
```

## 4) Input Behavior

`--gene` accepts:
- UniProt ID
- NCBI protein accession
- local FASTA path

`--query_id` is a legacy alternative (still supported).

If query is nucleotide FASTA, SynTerra normalizes to protein space before downstream steps.

## 5) Home/Target Annotation Behavior

- If home GFF is present: use it.
- If home GFF is missing: run regional Prodigal prediction.
- Borrowed annotations from annotated targets are merged into fallback home annotations.
- Target genomes do not need GFF.
- Flanking annotations are retained for synteny context, but iterative DB expansion is GOI-driven.

## 6) Core Parameters

### Mode and I/O

- `--mode` (default: `easy`): `easy` or `pro`
- `--gene` (default: `null`): query ID or FASTA path
- `--query_id` (default: `null`): legacy ID input
- `--home_species` (default: `null`): required in some easy-mode cases
- `--home_genome` (default: `null`): required in pro mode
- `--home_gff` (default: `null`): optional in pro mode
- `--target_genomes` (default: `null`): target FASTA glob in pro mode
- `--target_species` (default: `null`): optional easy-mode species override list
- `--max_genomes` (default: `0`): easy-mode target count (`0` = auto policy)
- `--assembly_ranking` (default: `hybrid`): easy-mode assembly ranking mode (`hybrid`, `counts`, `nstats`)
- `--bad_quality_policy` (default: `ask`): low-quality assembly handling (`ask`, `drop`, `keep`)
- `--bad_quality_timeout` (default: `300`): seconds to wait for answer when policy=`ask`
- `--bad_max_contigs` (default: `100000`): low-quality threshold for contig count
- `--bad_max_scaffolds` (default: `50000`): low-quality threshold for scaffold count
- `--bad_min_n50` (default: `20000`): low-quality threshold for best available N50
- `--outdir` (default: `results`)

### Synteny and Search

- `--n_flanking_genes` (default: `10`)
- `--min_flanking_size` (default: `500`)
- `--prefer_large_genes` (default: `true`)
- `--exon_level_search` (default: `true`)
- `--cluster_distance` (default: `50000`)
- `--min_synteny_score` (default: `0.6`)
- `--min_hit_identity` (default: `40`)
- `--min_hit_length` (default: `100`)
- `--search_evalue` (default: `1e-5`)
- `--mmseqs_sensitivity` (default: `8.5`)
- `--max_intron` (default: `20000`)
- `--max_blocks_per_genome` (default: `80`)
- `--min_block_genes` (default: `2`)
- `--max_consecutive_empty_blocks` (default: `25`)

### GOI/Exon Refinement

- `--gff_search_window` (default: `100000`)
- `--gap_search_window` (default: `50000`)
- `--gap_min_size` (default: `10`)
- `--gap_evalue` (default: `10`)
- `--gap_min_identity` (default: `25.0`)
- `--gap_min_alnlen` (default: `10`)
- `--gap_max_hits` (default: `5`)
- `--min_exon_query_cov` (default: `0.25`)
- `--min_exon_alnlen` (default: `30`)

### Smith-Waterman Augmentation

- `--enable_smith_waterman` (default: `true`)
- `--sw_method` (default: `auto`)
- `--sw_min_score` (default: `50`)
- `--sw_min_identity` (default: `20.0`)
- `--sw_timeout_seconds` (default: `300`)

### Region Padding / Relaxed Augmented Search

- `--region_padding` (default: `150000`)
- `--padding_min` (default: `50000`)
- `--padding_max` (default: `200000`)
- `--aug_relaxed_evalue_mult` (default: `1000`)
- `--aug_relaxed_evalue_cap` (default: `10.0`)
- `--aug_relaxed_parse_evalue_mult` (default: `10`)
- `--aug_relaxed_identity_factor` (default: `0.6`)
- `--aug_relaxed_identity_min` (default: `25.0`)
- `--aug_relaxed_length_div` (default: `2`)
- `--aug_relaxed_length_min` (default: `15`)
- `--aug_dedup_bin_bp` (default: `100`)

### Prodigal Fallback

- `--pred_flank_window` (default: `50000`)
- `--pred_keep_pct` (default: `0.10`)
- `--prodigal_full_genome_fallback` (default: `false`)

### Synteny Scoring Weights

- `--synteny_weight_base` (default: `0.4`)
- `--synteny_weight_consistency` (default: `0.3`)
- `--synteny_weight_strand` (default: `0.3`)

### Reserved / Future Parameters

These are defined but currently not central to the main decision path:
- `--augustus_species`
- `--expand_db_threshold`
- `--diamond_sensitivity`
- `--enable_splice_variants`
- `--enable_frameshifts`
- `--mutation_rate`
- `--num_mutant_variants`
- `--keep_intermediate`
- `--max_retries`

## 7) Outputs

Main outputs:
- `synterra_report.json`
- `*_synteny_plot.html`
- `*_tree.nwk`
- `regions/*.regions.bed`

Plot behavior:
- Target genes are filtered to candidate region BED intervals before rendering, so off-locus annotations are not drawn.

Intermediate outputs:
- `intermediate/locate_gene/`
- `intermediate/annotate_goi/`
- `intermediate/flanking/`
- `intermediate/initial_db/`
- `intermediate/phylo_sort/`
- `intermediate/query/`

## 8) Standalone Local Recheck (ground-truth test harness)

Uses the project test script, useful during method tuning.

```bash
conda run --no-capture-output -n syntenyfinder \
  python scripts/reproduce_annotation.py \
  --outdir tests/ground_truth_test/output_recheck
```
