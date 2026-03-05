# SynTerra — Usage & Reference Manual

Detailed reference for running and configuring SynTerra.  
For initial setup instructions, see the [README](README.md).

---

## Table of Contents

1. [Execution Modes](#1-execution-modes)
2. [Execution Profiles](#2-execution-profiles)
3. [Algorithm Overview](#3-algorithm-overview)
4. [Full Parameter Reference](#4-full-parameter-reference)
5. [Running on HPC / SLURM](#5-running-on-hpc--slurm)
6. [Output Files](#6-output-files)
7. [Resuming & Caching](#7-resuming--caching)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Execution Modes

SynTerra has two modes: **Easy** (automated genome retrieval) and **Pro** (local files).

### Easy Mode

Provide a UniProt or NCBI protein accession. SynTerra resolves the query, fetches the reference genome and related target assemblies from NCBI, and runs the full analysis.

```bash
nextflow run main.nf \
  --mode easy \
  --query_id Q16553 \
  --max_genomes 5 \
  --outdir results/my_run \
  -profile standard
```

**Required:**

| Flag | Description |
|---|---|
| `--mode easy` | Select Easy Mode |
| `--query_id` | UniProt accession (e.g. `Q16553`) or NCBI protein ID |

**Optional (Easy Mode only):**

| Flag | Default | Description |
|---|---|---|
| `--home_species` | auto-detected | Override the reference species (e.g. `"Homo sapiens"`) |
| `--target_species` | auto (taxonomic) | Comma-separated species list instead of auto-selection (e.g. `"Gallus gallus,Mus musculus"`) |
| `--max_genomes` | `0` (auto) | Number of related genomes to fetch. `0` = auto-detect (~3 per taxonomic level) |
| `--assembly_ranking` | `hybrid` | How to pick among multiple assemblies for one species: `hybrid`, `counts`, `nstats` |
| `--bad_quality_policy` | `drop` | What to do with low-quality assemblies: `drop`, `keep`, or `ask` (interactive prompt) |
| `--bad_quality_timeout` | `300` | Seconds to wait for user input when policy is `ask` |
| `--bad_max_contigs` | `500000` | Assemblies with more contigs are flagged as low quality |
| `--bad_max_scaffolds` | `500000` | Assemblies with more scaffolds are flagged as low quality |
| `--bad_min_n50` | `5000` | Assemblies with N50 below this are flagged as low quality |

### Pro Mode

Supply your own query FASTA, reference genome, and target genome files. Works offline.

```bash
nextflow run main.nf \
  --mode pro \
  --query queries/melittin.faa \
  --home_genome /data/apis_mellifera.fna \
  --home_gff /data/apis_mellifera.gff \
  --target_genomes "/data/targets/*.fna" \
  --outdir results/melittin \
  -profile standard
```

**Required:**

| Flag | Description |
|---|---|
| `--mode pro` | Select Pro Mode |
| `--query` | Path to query protein FASTA (DNA sequences are auto-translated to the best ORF) |
| `--home_genome` | Path to reference genome FASTA (`.fna` or `.fna.gz`) |
| `--target_genomes` | Glob pattern or comma-separated list matching target genome FASTAs |

**Optional:**

| Flag | Description |
|---|---|
| `--home_gff` | GFF annotation for the home genome. Highly recommended — provides much better flanking-gene extraction than Prodigal fallback. |

> **Tip:** `--target_genomes` accepts globs (`"genomes/*.fna"`), comma-separated paths (`"a.fna,b.fna"`), or Nextflow list syntax.

---

## 2. Execution Profiles

Append a profile with `-profile <name>` to control how resources are allocated. Combine with a comma for multi-profile runs (e.g. `-profile docker,laptop_safe` is **not** supported — pick one).

| Profile | Executor | Environment | Description |
|---|---|---|---|
| `standard` | local | Conda | Default. 2 CPUs / 6 GB RAM per iterative-search task, single-fork. Good baseline for workstations. |
| `conda` | local | Conda | Same as `standard` but explicitly disables Docker/Singularity. |
| `laptop_safe` | local | Conda | Conservative. 1 CPU, single task at a time, high memory ceiling (12 GB) but strict fork limits. Prevents system freezes on machines with limited RAM. |
| `docker` | local | Docker | Runs all processes inside the `synterra-local:latest` container. Build it first with `docker build -t synterra-local:latest .` |
| `docker_max` | local | Docker | Auto-detects all host CPUs and RAM. Allocates nearly everything to the heaviest tasks (MMseqs2, ITERATIVE_SEARCH). Single-fork to avoid OOM. Ideal for dedicated machines. |
| `singularity` | local | Singularity | Like `docker` but uses Singularity with auto-mounts. |
| `slurm` | SLURM | (none) | Submits tasks to a SLURM scheduler. Edit `nextflow.config` to set your partition and account. |
| `hpc_singularity` | SLURM | Singularity | SLURM + Singularity containers. Caches images in `~/.singularity/cache`. |
| `hpc_conda` | SLURM | Conda+Mamba | SLURM + Conda (uses Mamba for faster env creation). |
| `test` | local | Conda | Loads `conf/test.config` with small test data and relaxed thresholds for CI. |

---

## 3. Algorithm Overview

The pipeline proceeds through five phases:

### Phase 1 — Gene Localization

1. **Normalize Query:** If the input is DNA, the best ORF is translated to protein.
2. **Locate GOI:** The query protein is aligned against the home genome using tblastn and MMseqs2 to establish coordinates.
3. **Annotate GOI Exons:** If a GFF is available, the GOI is matched to an annotated gene and individual CDS/exons are extracted. Otherwise, exon boundaries are inferred from alignment hits.
4. **Split Loci:** If the GOI maps to multiple genomic locations (e.g. tandem duplicates), each locus is processed independently.
5. **Extract Flanking Genes:** The *n* genes upstream and downstream of each locus are identified from GFF or Prodigal prediction.
6. **Borrow Annotations:** When the home genome lacks a GFF, annotations can be borrowed from annotated target genomes via reciprocal best hits.

### Phase 2 — Phylogenetic Ordering & Iterative Search

7. **Phylo Sort:** Target genomes are ordered by evolutionary distance to the reference.
8. **Genome Quality Assessment:** Target assemblies are evaluated for contiguity (N50, scaffold count).
9. **Iterative Search:** For each target genome (nearest-first), flanking genes are mapped with MMseqs2. Hits are clustered into candidate syntenic blocks. Within each block, localized tblastn, miniprot, and Smith-Waterman searches attempt to find the GOI. Discovered genes are added to the search database, improving sensitivity for more distant species.

### Phase 3 — Region Clustering

10. **Cluster Regions:** Candidate blocks across all targets are filtered by synteny score and ranked.

### Phase 4 — Phylogenetics & Visualization

11. **Compute Tree:** Discovered GOI sequences are aligned (MAFFT) and a phylogenetic tree is inferred (IQ-TREE).
12. **Plot Synteny:** An interactive HTML plot shows the syntenic context of each hit, colored by homology, with the phylogenetic tree alongside.

### Phase 5 — Reporting

13. **Generate Report:** A JSON summary file is produced with run parameters, genome QC results, and per-target outcomes.

---

## 4. Full Parameter Reference

All parameters can be set on the command line (`--param value`) or in a custom config file (`-c my.config`).

### Synteny & Search

| Parameter | Default | Description |
|---|---|---|
| `--n_flanking_genes` | `10` | Number of flanking genes to extract on each side of the GOI |
| `--prefer_large_genes` | `true` | Prefer larger flanking genes (more informative for homology search) |
| `--min_flanking_size` | `500` | Minimum size (bp) for a flanking gene to be included |
| `--exon_level_search` | `true` | Search at exon level for better divergent-gene detection |
| `--cluster_distance` | `150000` | Max gap (bp) between flanking-gene hits to merge into one syntenic block |
| `--min_synteny_score` | `0.6` | Fraction of flanking genes that must map to a target to trigger local search |
| `--min_hit_identity` | `10` | Minimum alignment identity (%) for an individual hit |
| `--min_hit_length` | `10` | Minimum alignment length for an individual hit |
| `--search_evalue` | `0.01` | E-value threshold for tblastn/MMseqs2 searches |
| `--max_intron` | `20000` | Maximum intron length (bp) for miniprot gene models |
| `--region_padding` | `150000` | Extra flanking sequence (bp) appended to each side of a candidate block |
| `--padding_min` | `50000` | Minimum padding (bp) |
| `--padding_max` | `200000` | Maximum padding (bp) |
| `--max_blocks_per_genome` | `80` | Safety cap on candidate blocks per target genome |
| `--min_block_genes` | `2` | Minimum flanking-gene hits in a block to keep it |
| `--max_consecutive_empty_blocks` | `25` | Stop expanding after this many consecutive empty blocks |

### Smith-Waterman Local Search

| Parameter | Default | Description |
|---|---|---|
| `--enable_smith_waterman` | `true` | Use rigorous Smith-Waterman alignment (parasail) for GOI search |
| `--sw_method` | `auto` | Implementation: `auto`, `parasail`, or `ssearch36` |
| `--sw_min_score` | `50` | Minimum SW alignment score to report a hit |
| `--sw_min_identity` | `20.0` | Minimum identity (%) for SW hits |
| `--sw_timeout_seconds` | `300` | Timeout per SW search invocation |

### Relaxed / Augmented Search

Controls the increasingly permissive search passes used for highly divergent targets.

| Parameter | Default | Description |
|---|---|---|
| `--aug_relaxed_evalue_mult` | `1000` | Multiply base e-value by this factor in relaxed passes |
| `--aug_relaxed_evalue_cap` | `10.0` | Maximum e-value allowed even in relaxed mode |
| `--aug_relaxed_identity_factor` | `0.6` | Multiply normal identity threshold by this in relaxed mode |
| `--aug_relaxed_identity_min` | `25.0` | Absolute minimum identity (%) in relaxed mode |
| `--aug_relaxed_length_div` | `2` | Divide normal length threshold by this in relaxed mode |
| `--aug_relaxed_length_min` | `15` | Absolute minimum alignment length in relaxed mode |
| `--aug_dedup_bin_bp` | `100` | Bin size (bp) for deduplicating overlapping relaxed hits |

### MMseqs2

| Parameter | Default | Description |
|---|---|---|
| `--mmseqs_sensitivity` | `9.5` | MMseqs2 sensitivity (1–10+). Higher = slower but more sensitive |
| `--mmseqs_split_memory_limit` | `3G` | MMseqs2 memory limit for database splitting |
| `--mmseqs_verbosity` | `1` | MMseqs2 log verbosity (0 = silent) |
| `--min_gene_identity` | `30` | Minimum identity (%) for flanking-gene MMseqs2 matches |

### Annotation & Prodigal

| Parameter | Default | Description |
|---|---|---|
| `--gff_search_window` | `100000` | Window (bp) around GOI to search in GFF for flanking genes |
| `--gap_search_window` | `50000` | Window for gap-filling searches |
| `--gap_min_size` | `10` | Minimum gap size (bp) to attempt fill |
| `--gap_evalue` | `10` | E-value for gap search |
| `--gap_min_identity` | `25.0` | Minimum identity (%) for gap hits |
| `--gap_min_alnlen` | `10` | Minimum alignment length for gap hits |
| `--gap_max_hits` | `5` | Max gap hits to report |
| `--min_exon_query_cov` | `0.25` | Minimum query coverage fraction for exon annotation |
| `--min_exon_alnlen` | `30` | Minimum exon alignment length |
| `--pred_flank_window` | `50000` | Prodigal prediction window around locus |
| `--pred_keep_pct` | `0.10` | Fraction of Prodigal predictions to keep |
| `--prodigal_full_genome_fallback` | `false` | Run Prodigal on entire genome if windowed prediction fails |

### Synteny Scoring Weights

| Parameter | Default | Description |
|---|---|---|
| `--synteny_weight_base` | `0.4` | Weight for base synteny score |
| `--synteny_weight_consistency` | `0.3` | Weight for gene-order consistency |
| `--synteny_weight_strand` | `0.3` | Weight for strand conservation |
| `--synteny_goi_overlap_bonus` | `0.15` | Bonus for blocks that overlap a GOI annotation |
| `--max_regions` | `0` | Max regions to emit per locus. `0` = adaptive (all above threshold, capped at 6) |

### Visualization

| Parameter | Default | Description |
|---|---|---|
| `--plot_width` | `1500` | Width of the output HTML plot (px) |
| `--gap_threshold` | `50000` | Gaps larger than this (bp) are visually compressed |
| `--gap_visual_size` | `20000` | Size (bp) used to represent compressed gaps |
| `--flank_fallback_bp` | `1000000` | Maximum window (bp) rendered around distal targets |
| `--scale_bar_len` | `10000` | Scale bar size (bp) |

### Resource Tuning

These control per-process resource allocation. Override them for your hardware.

| Parameter | Default | Description |
|---|---|---|
| `--iterative_search_cpus` | `2` | CPUs for ITERATIVE_SEARCH tasks |
| `--iterative_search_memory` | `6 GB` | RAM for ITERATIVE_SEARCH tasks |
| `--iterative_search_max_forks` | `1` | Max parallel ITERATIVE_SEARCH tasks |
| `--locate_gene_cpus` | `1` | CPUs for LOCATE_GENE |
| `--locate_gene_memory` | `3 GB` | RAM for LOCATE_GENE |

### Output

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `results` | Directory for pipeline output |
| `--keep_intermediate` | `false` | Keep intermediate files (useful for debugging) |

---

## 5. Running on HPC / SLURM

A ready-made submission script is provided in `slurm_submit.sh`. Edit the variables at the top and submit:

```bash
# Edit slurm_submit.sh to set your query, species, partition, and account
sbatch slurm_submit.sh
```

The script submits the Nextflow **controller** as a SLURM job. Nextflow then submits individual pipeline tasks as separate SLURM jobs via the `hpc_singularity` or `hpc_conda` profile.

**Key environment variables for HPC:**

```bash
export NXF_WORK="${SCRATCH}/work"          # fast scratch for intermediates
export NXF_SINGULARITY_CACHEDIR="${HOME}/.singularity/cache"
```

**Manual SLURM example:**

```bash
nextflow run main.nf \
  -profile hpc_conda \
  --mode easy \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/melittin_hpc \
  -work-dir "${SCRATCH}/work" \
  -resume
```

> Edit the `slurm` profile in `nextflow.config` to change the default partition (`normal`) or add `--account`.

---

## 6. Output Files

All output goes into `--outdir` (default: `results/`):

| Path | Description |
|---|---|
| `*_synteny_plot.html` | Interactive HTML visualization. Open in a browser — shows syntenic blocks, gene arrows, homology links, and a phylogenetic tree. |
| `*_tree.nwk` | Newick-format phylogenetic tree of discovered GOI sequences. |
| `regions/*.regions.bed` | BED files with genomic coordinates of identified candidate syntenic blocks on each target genome. |
| `synterra_report.json` | Machine-readable JSON report: input parameters, genome QC metrics, per-target results, internal exit codes. |
| `intermediate/` | Per-phase artifacts — flanking gene FASTAs, MMseqs2 hit tables, per-target GFFs, miniprot alignments, etc. Only kept if `--keep_intermediate true`. |
| `downloaded_genomes/` | (Easy Mode only) Downloaded genome assemblies and `assembly_quality.tsv` with contiguity stats. |

---

## 7. Resuming & Caching

Nextflow caches completed tasks in the `work/` directory. To resume after a crash or parameter change:

```bash
nextflow run main.nf [your params] -resume
```

Only tasks whose inputs changed will be re-executed. This is especially useful for:
- Adding more target genomes to an existing run
- Tweaking visualization parameters (only `PLOT_SYNTENY` will re-run)
- Recovering from transient network failures in Easy Mode

**Cleaning up:**

```bash
# Remove work directory (frees disk space, loses cache)
rm -rf work/

# Or use Nextflow's built-in cleanup
nextflow clean -f
```

---

## 8. Troubleshooting

### Pipeline crashes with SIGKILL (exit code 137) during LOCATE_GENE or ITERATIVE_SEARCH

**Cause:** Out of memory. MMseqs2 database indexing and tblastn can be RAM-intensive.

**Fix:**
- Reduce MMseqs2 memory with `--mmseqs_split_memory_limit '1G'`
- Use `-profile laptop_safe` to constrain parallelism
- Increase available memory or switch to an HPC profile
- Check that `/tmp` has sufficient free space (MMseqs2 uses it for temporary files)

### ITERATIVE_SEARCH runs for a very long time (>40 min per genome)

**Cause:** The target genome contains large tandem duplication arrays or many candidate blocks, causing exhaustive local searches.

**What to do:**
- This can be normal for complex genomes. Check progress:
  ```bash
  # Find the task's work directory
  ls -lt work/*/*/.command.log | head -5
  # Follow the log
  tail -f work/<hash>/<hash>/.command.log
  ```
- To speed things up, reduce `--max_blocks_per_genome` (default 80) or increase `--min_synteny_score` (default 0.6)

### "Query FASTA not found" or "Home genome not found"

**Cause:** Paths are relative to the Nextflow launch directory, not to the script.

**Fix:** Use absolute paths or ensure you run `nextflow` from the SynTerra project root.

### Conda environment creation times out

**Cause:** Conda solver is slow.

**Fix:**
- Install [Mamba](https://github.com/mamba-org/mamba) — Nextflow will use it automatically (`conda.useMamba = true` is set in `nextflow.config`)
- Or increase the timeout: the config already sets `conda.createTimeout = '1 h'`

### Easy Mode fails to download genomes

**Cause:** NCBI API rate limiting or network issues.

**Fix:**
- Set an NCBI API key: `export NCBI_API_KEY=your_key_here`
- Re-run with `-resume` — completed downloads will be cached
- Check NCBI service status at https://www.ncbi.nlm.nih.gov/

### No synteny plot is generated

**Cause:** No candidate regions passed the synteny score threshold.

**Fix:**
- Lower `--min_synteny_score` (e.g. `0.3`)
- Increase `--n_flanking_genes` (e.g. `15`) to capture more genomic context
- Increase `--region_padding` to widen the search window
- Check `regions/*.bed` files — if they are empty, the flanking-gene mapping step did not find hits

### "Java not found" or Nextflow fails to start

**Fix:**
```bash
# Check Java version
java -version

# If missing, install via Conda
conda install -c conda-forge openjdk=17

# Or via system package manager
sudo apt install default-jdk   # Debian/Ubuntu
```
