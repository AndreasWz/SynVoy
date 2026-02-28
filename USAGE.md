# SynTerra Usage

This document reflects the current pipeline behavior in this repository.
Use it as the primary runtime reference for setup, execution, parameters, outputs, and troubleshooting.

## 1) Requirements

- Linux or macOS (Windows via WSL2 is fine).
- Nextflow `>=22.10.1`.
- One execution backend: Conda/Mamba (recommended), Docker, or Singularity/Apptainer.
- Internet access for `--mode easy` (query resolution + genome downloads).

## 2) Setup

### Option A: Conda/Mamba (recommended)

```bash
cd /path/to/SynTerra
conda env create -f environment.yml
conda activate syntenyfinder
nextflow -version
nextflow run main.nf --help
```

If `nextflow` is not on your `PATH`, use `./nextflow`.

### Option B: Docker

```bash
docker build -t synterra-local:latest .
nextflow run main.nf -profile docker --mode easy --gene P01501 --outdir results
```

Override image:

```bash
nextflow run main.nf -profile docker --docker_container your/image:tag --mode easy --query_id P01501 --outdir results
```

### Option C: Singularity/Apptainer

```bash
mkdir -p "$HOME/.singularity/cache"
export NXF_SINGULARITY_CACHEDIR="$HOME/.singularity/cache"
nextflow run main.nf -profile singularity --mode easy --query_id P01501 --outdir results
```

### Optional smoke test

```bash
nextflow run main.nf -profile test
```

## 3) How to Run

`nextflow run main.nf` is the recommended launcher.
It forwards pipeline arguments to Nextflow and adds a compact live UI.

Examples:

```bash
nextflow run main.nf --mode easy --gene P01501 --outdir results
nextflow run main.nf --raw --mode easy --gene P01501 --outdir results
```

- `--raw` prints raw Nextflow output.
- Raw logs are written to `.synterra_logs/run_YYYYMMDD_HHMMSS.log`.

## 4) Run Modes

### Easy mode

SynTerra resolves query input, fetches the home genome, and fetches target genomes.

```bash
nextflow run main.nf \
  --mode easy \
  --query_id P01501 \
  --outdir results
```

### Pro mode

You provide local home and target genomes.

```bash
nextflow run main.nf \
  --mode pro \
  --query input/query.fasta \
  --home_genome input/home.fna \
  --home_gff input/home.gff \
  --target_genomes "input/targets/*.fna" \
  --outdir results
```

### Resume

Use the same command and same `--outdir`, then add `-resume`.

```bash
nextflow run main.nf --mode easy --gene P01501 --outdir results -resume
```

## 5) Query Input Behavior

SynTerra strictly separates input parameters by `--mode`:

- `--mode easy`: Requires `--query_id` (a UniProt ID or NCBI accession). Do not pass a file.
- `--mode pro`: Requires `--query` (a local FASTA file path).

Rules:
- In easy mode, `--home_species` is auto-detected from the UniProt/NCBI annotation. You can provide `--home_species` to override it.
- In pro mode, a FASTA path must be supplied for both `--query` and `--home_genome`.
- DNA query FASTA is normalized to protein space before search/annotation.

## 6) Easy-Mode Genome Selection

Home genome selection:

- Try reference/representative assembly first.
- Fallback to quality-ranked assembly for the requested species.
- Final fallback can use closest taxonomic relative if no species assembly exists.

Target genome selection:

- If `--target_species` is provided, SynTerra fetches those taxa directly.
- If not provided, SynTerra auto-searches related taxa (genus/family/order/class) up to `--max_genomes`.

Assembly ranking (`--assembly_ranking`):

- `hybrid` (default): contiguity counts + N-stats + assembly level priority.
- `counts`: prioritize fewer contigs/scaffolds/chromosomes.
- `nstats`: prioritize higher N50/N80.

Low-quality handling:

- `--bad_quality_policy ask|drop|keep` (default `ask`).
- `ask` waits up to `--bad_quality_timeout` seconds (default `300`) then defaults to NO.
- Threshold knobs: `--bad_max_contigs`, `--bad_max_scaffolds`, `--bad_min_n50`.

Quality report file:

- `downloaded_genomes/easy_mode_genomes/assembly_quality.tsv`

## 7) Annotation and Flanking Logic

- If home GFF exists and is usable, SynTerra prioritizes annotation-based extraction.
- Without usable home GFF, fallback paths predict genes locally.
- Borrowed annotation support is available for weak/no-annotation home genomes.
- Flanking genes are extracted per locus and retained for context.
- Flanking records keep stable IDs and include display labels from annotation names when available.
- Iterative region expansion remains GOI-driven.

## 8) Desktop-Stable Docker Recipe (3-snake 3FTx)

```bash
conda activate syntenyfinder

NXF_OPTS='-Xms512m -Xmx2g' nextflow run main.nf \
  -profile docker \
  --mode easy \
  --query_id P60615 \
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

## 9) Full Parameter Reference

Defaults are taken from `nextflow.config`.

### Mode and I/O

| Parameter | Default | Description |
|---|---|---|
| `--mode` | `easy` | `easy` or `pro`. |
| `--query` | `null` | Pro-mode path to local query FASTA file. |
| `--query_id` | `null` | Easy-mode query UniProt/NCBI Accession ID. |
| `--home_species` | `null` | Easy-mode home species name. |
| `--home_genome` | `null` | Pro-mode home genome FASTA path. |
| `--home_gff` | `null` | Pro-mode home GFF path (optional). |
| `--target_genomes` | `null` | Pro-mode target genome glob. |
| `--target_species` | `null` | Easy-mode comma-separated species list. |
| `--max_genomes` | `0` | Easy-mode target count (`0` = auto strategy). |
| `--outdir` | `results` | Output directory. |
| `--docker_container` | `synterra-local:latest` | Container image for docker profile. |

### Easy-mode assembly filtering and ranking

| Parameter | Default | Description |
|---|---|---|
| `--assembly_ranking` | `hybrid` | Assembly ranking mode: `hybrid`, `counts`, `nstats`. |
| `--bad_quality_policy` | `ask` | Low-quality policy: `ask`, `drop`, `keep`. |
| `--bad_quality_timeout` | `300` | Seconds to wait for `ask` policy. |
| `--bad_max_contigs` | `100000` | Assemblies above this are flagged low-quality. |
| `--bad_max_scaffolds` | `50000` | Assemblies above this are flagged low-quality. |
| `--bad_min_n50` | `20000` | Assemblies below this best-N50 are flagged low-quality. |

### Synteny and search controls

| Parameter | Default | Description |
|---|---|---|
| `--n_flanking_genes` | `10` | Number of flanking genes to extract per locus side. |
| `--prefer_large_genes` | `true` | Prefer larger gene models where relevant. |
| `--min_flanking_size` | `500` | Minimum flanking gene span (bp). |
| `--exon_level_search` | `true` | Enables exon-oriented flanking behavior. |
| `--cluster_distance` | `50000` | Region clustering distance (bp). |
| `--min_synteny_score` | `0.6` | Minimum synteny score threshold. |
| `--min_hit_identity` | `40` | Minimum identity for hits in iterative search. |
| `--min_hit_length` | `100` | Minimum hit length. |
| `--search_evalue` | `1e-5` | Core MMseqs e-value threshold. |
| `--max_intron` | `20000` | Max intron size used in model assembly. |
| `--mmseqs_sensitivity` | `8.5` | MMseqs sensitivity setting. |
| `--mmseqs_split_memory_limit` | `3G` | MMseqs internal split memory limit. |
| `--mmseqs_verbosity` | `1` | MMseqs verbosity. |
| `--min_gene_identity` | `30` | Minimum identity for gene-level consolidation steps. |

### GOI refinement and gap search

| Parameter | Default | Description |
|---|---|---|
| `--gff_search_window` | `100000` | Window around hits for GFF match recovery. |
| `--gap_search_window` | `50000` | GOI gap-search window. |
| `--gap_min_size` | `10` | Minimum gap size for gap search. |
| `--gap_evalue` | `10` | Gap search e-value threshold. |
| `--gap_min_identity` | `25.0` | Gap search minimum identity. |
| `--gap_min_alnlen` | `10` | Gap search minimum alignment length. |
| `--gap_max_hits` | `5` | Max gap hits retained. |
| `--min_exon_query_cov` | `0.25` | Minimum exon query coverage. |
| `--min_exon_alnlen` | `30` | Minimum exon alignment length. |

### Smith-Waterman augmentation

| Parameter | Default | Description |
|---|---|---|
| `--enable_smith_waterman` | `true` | Enables Smith-Waterman local refinement. |
| `--sw_method` | `auto` | SW engine selection (`auto`, `parasail`, `ssearch36`). |
| `--sw_min_score` | `50` | Minimum SW score. |
| `--sw_min_identity` | `20.0` | Minimum SW identity. |
| `--sw_timeout_seconds` | `300` | SW timeout per task. |

### Region padding and relaxed augmented search

| Parameter | Default | Description |
|---|---|---|
| `--region_padding` | `150000` | Base region padding around loci. |
| `--padding_min` | `50000` | Minimum adaptive padding. |
| `--padding_max` | `200000` | Maximum adaptive padding. |
| `--aug_relaxed_evalue_mult` | `1000` | Multiplier for relaxed e-value pass. |
| `--aug_relaxed_evalue_cap` | `10.0` | Cap for relaxed e-value. |
| `--aug_relaxed_parse_evalue_mult` | `10` | Relaxed parse threshold multiplier. |
| `--aug_relaxed_identity_factor` | `0.6` | Identity relaxation factor. |
| `--aug_relaxed_identity_min` | `25.0` | Floor for relaxed identity. |
| `--aug_relaxed_length_div` | `2` | Length relaxation divisor. |
| `--aug_relaxed_length_min` | `15` | Minimum relaxed length. |
| `--aug_dedup_bin_bp` | `100` | Deduplication bin size (bp). |
| `--max_blocks_per_genome` | `80` | Cap on candidate blocks per target genome. |
| `--min_block_genes` | `2` | Minimum genes per block. |
| `--max_consecutive_empty_blocks` | `25` | Early-stop guard for empty block streaks. |

### Resource guards

| Parameter | Default | Description |
|---|---|---|
| `--iterative_search_cpus` | `2` | CPUs for iterative search process. |
| `--iterative_search_memory` | `6 GB` | Memory for iterative search process. |
| `--iterative_search_max_forks` | `1` | Max parallel forks for iterative search. |
| `--iterative_quiet_subtools` | `true` | Quieter subtool logs in iterative phase. |
| `--locate_gene_cpus` | `1` | CPUs for home GOI localization process. |
| `--locate_gene_memory` | `3 GB` | Memory for home GOI localization process. |

### Fallback prediction and scoring

| Parameter | Default | Description |
|---|---|---|
| `--pred_flank_window` | `50000` | Window for fallback Prodigal extraction. |
| `--pred_keep_pct` | `0.10` | Fraction of fallback predictions retained. |
| `--prodigal_full_genome_fallback` | `false` | Enables full-genome Prodigal fallback path. |
| `--synteny_weight_base` | `0.4` | Base term weight in synteny scoring. |
| `--synteny_weight_consistency` | `0.3` | Consistency term weight in synteny scoring. |
| `--synteny_weight_strand` | `0.3` | Strand term weight in synteny scoring. |

### Advanced and reserved

| Parameter | Default | Description |
|---|---|---|
| `--keep_intermediate` | `false` | Keep extra intermediate artifacts. |
| `--max_retries` | `3` | Pipeline retry-related knob (advanced). |
| `--expand_db_threshold` | `1e-10` | Reserved for future wiring. |
| `--diamond_sensitivity` | `very-sensitive` | Reserved for future wiring. |
| `--enable_splice_variants` | `true` | Reserved for future wiring. |
| `--enable_frameshifts` | `true` | Reserved for future wiring. |
| `--mutation_rate` | `0.05` | Reserved for future wiring. |
| `--num_mutant_variants` | `10` | Reserved for future wiring. |

## 10) Output Structure

Typical outputs in `--outdir`:

- `synterra_report.json`
- `*_synteny_plot.html`
- `*_tree.nwk`
- `regions/*.regions.bed`
- `downloaded_genomes/easy_mode_genomes/genomes_manifest.txt` (easy mode)
- `downloaded_genomes/easy_mode_genomes/species_mapping.tsv` (easy mode)
- `downloaded_genomes/easy_mode_genomes/assembly_quality.tsv` (easy mode)
- `intermediate/locate_gene/`
- `intermediate/annotate_goi/`
- `intermediate/flanking/`
- `intermediate/initial_db/`
- `intermediate/phylo_sort/`
- `intermediate/query/`

## 11) Logs, Resume, and Cleanup

Last live log:

```bash
ls -1t .synterra_logs/run_*.log | head -n 1
```

Follow most recent log:

```bash
tail -f "$(ls -1t .synterra_logs/run_*.log | head -n 1)"
```

Resume run:

```bash
nextflow run main.nf --mode easy --query_id P01501 --outdir results -resume
```

Remove old outputs and work directory:

```bash
rm -rf results_* work
```

## 12) Troubleshooting

### `RESOLVE_GENE_INPUT` failed

- Check whether local FASTA path exists.
- For ID inputs, verify identifier correctness and internet connectivity.
- In easy mode with local FASTA, set `--home_species`.

### Pipeline appears stuck at fetch stage

- Easy mode can spend long time in NCBI metadata + genome download.
- Check `.synterra_logs` to confirm progress.
- Reduce search/download scope using explicit `--target_species` and small `--max_genomes`.

### Very long runtime or desktop instability

- Use docker profile.
- Limit JVM heap with `NXF_OPTS='-Xms512m -Xmx2g'`.
- Keep iterative search conservative: `--iterative_search_cpus 2`, `--iterative_search_max_forks 1`, `--iterative_search_memory '6 GB'`.
- Reduce `--mmseqs_split_memory_limit` when memory pressure persists.

### No loci found (`SPLIT_LOCI done: identified 0`)

- Increase search sensitivity or relax thresholds.
- Try higher `--mmseqs_sensitivity`.
- Try looser `--search_evalue`.
- Re-check query correctness (protein sequence, ID, species context).

### Synteny plot labels look generic

- Prefer assemblies with usable GFF annotations.
- Verify home/target GFF availability in run outputs.
- Check whether the selected assemblies in `assembly_quality.tsv` are scaffold-heavy or weakly annotated.
