# SynTerra Execution Manual

SynTerra is an experimental Nextflow pipeline developed for the identification of orthologous genes across related species using macro-syntenic mapping and local iterative alignment. 

This document details the configuration parameters, execution modes, output artifacts, and algorithm architecture.

---

## 1. Execution Environments

SynTerra requires Nextflow version `>=22.10.1`. Resource requirements vary significantly depending on the profile used and the parameters defined in `nextflow.config`.

### Recommended Execution (Conda)
For headless environments without Docker support:
```bash
conda env create -f environment.yml
conda activate syntenyfinder
nextflow run main.nf [parameters]
```

### Containerized Execution (Docker / Singularity)
The `-profile docker` or `-profile singularity` arguments bypass local environment setups.
```bash
docker build -t synterra-local:latest .
```

### Resource Profiles
The profile appended to the run command dictates how the pipeline scales on the host system:
*   `-profile standard`: Default executor. Utilizes `6 GB` RAM and `2` CPUs per iterative search task.
*   `-profile laptop_safe`: Strict resource restriction. Constrains execution to a single task at a time and strictly manages memory JVM heaps to prevent system stalls.
*   `-profile docker_max`: Hardware-aware autoscaling. Attempts to allocate the absolute maximum physical cores and available memory on the host system to the heaviest MMseqs2 tasks.
*   `-profile slurm`: Native configuration for HPC batch submission.

---

## 2. Core Execution Modes

The pipeline operates in two distinct operational modes.

### A. Easy Mode 
Automatically fetches taxonomy data, assemblies, and the query sequence from NCBI/UniProt databases given an accession identifier.

```bash
nextflow run main.nf \
  --mode easy \
  --query_id Q16553 \
  --max_genomes 5 \
  --outdir results
```
**Required Arguments:**
*   `--mode easy`
*   `--query_id`: A valid UniProt or NCBI identifier.

### B. Pro Mode
Expects locally provided target genome assemblies (`.fna` or `.fna.gz`), a local reference genome, and a local query FASTA file.

```bash
nextflow run main.nf \
  --mode pro \
  --query local_gene.fasta \
  --home_genome reference_assembly.fna \
  --target_genomes "target_assemblies/*.fna" \
  --outdir results
```
**Required Arguments:**
*   `--mode pro`
*   `--query`: Path to local sequence FASTA.
*   `--home_genome`: Path to local reference assembly.
*   `--target_genomes`: Glob pattern matching target genomes.

---

## 3. Algorithm Methodology

The default search behavior (`ITERATIVE_SEARCH`) is divided into several deterministic phases:
1.  **Reference Localization:** The pipeline aligns the user-provided query sequence to the `--home_genome` to define the baseline coordinate locus.
2.  **Context Extraction:** The *n* genes immediately upstream and downstream of the localized query are identified via GFF annotation or local Prodigal prediction.
3.  **Target Mappings:** The extracted flanking genes are queried against the `--target_genomes` dataset using `MMseqs2`.
4.  **Clustering:** Adjacent mappings on the target sequence are grouped into candidate syntenic blocks based on the `--cluster_distance` threshold.
5.  **Target Search:** The candidate blocks undergo an iterative sequence search, parsing alignment scores from increasingly relaxed sequence matching logic (`miniprot`, `tblastn`) to identify sequences morphologically matching the initial query.

---

## 4. Parameter Reference

The following parameters can be overridden in the command line invocation.

### Search and Synteny Boundaries

| Parameter | Default | Definition |
| :--- | :--- | :--- |
| `--n_flanking_genes` | `10` | The number of flanking genes parsed per boundary (upstream and downstream) from the reference locus. |
| `--cluster_distance` | `50000` | The minimum allowable gap (in base pairs) between target mappings to establish a unified syntenic candidate block. |
| `--min_synteny_score` | `0.6` | The fractional threshold of the retrieved flanking genes required to map onto the target to trigger sequence scanning. |
| `--region_padding` | `150000` | Genomic sequence (bp) artificially concatenated to the left and right borders of a synthesized candidate block prior to local search. |
| `--max_blocks_per_genome` | `80` | Safety limit on the maximum number of structural blocks evaluated per target sequence. |

### Quality Filtration and Fetch Constraints

Parameters evaluating automatically retrieved assemblies via Easy Mode.

| Parameter | Default | Definition |
| :--- | :--- | :--- |
| `--assembly_ranking` | `hybrid` | Scoring metric used when multiple genomes resolve to a requested target species. |
| `--bad_quality_policy` | `ask` | Action invoked when retrieved genomic metadata flags poor contiguity (`ask`, `drop`, `keep`). |
| `--bad_quality_timeout` | `300` | Duration (seconds) the runtime pauses awaiting user standard input during the `ask` policy. |
| `--bad_max_scaffolds` | `50000` | Threshold of maximum fragmentation permitted. Genomes breaching this limit are flagged. |
| `--bad_min_n50` | `20000` | Minimum N50 (contiguity depth) allowed. Assemblies beneath this are flagged. |

### Inference Constraints

| Parameter | Default | Definition |
| :--- | :--- | :--- |
| `--search_evalue` | `0.01` | The baseline e-value required for standard iteration logic to map a protein segment. |
| `--aug_relaxed_identity_min` | `25.0` | In extreme sequence divergence passes, the absolute identity (%) minimum the parser will allow. |
| `--mmseqs_sensitivity` | `9.5` | Operational sensitivity setting for underlying MMseqs2 calls. Values >10 incur massive computational penalties. |

### Visualization Constraints

| Parameter | Default | Definition |
| :--- | :--- | :--- |
| `--gap_threshold` | `50000` | Unannotated sequences mapping larger than this value (bp) are compressed in the layout. |
| `--gap_visual_size` | `3000` | The absolute pixel coordinate assigned to represent instances of compressed gaps. |
| `--flank_fallback_bp` | `1000000` | Maximum bounding boundary limit (bp) evaluated when rendering distal targets on the HTML canvas. |
| `--scale_bar_len` | `10000` | The physical representation metric evaluated by the plot scale bar (bp). |
| `--plot_width` | `1500` | Constraint width bounding the output HTML matrix (pixels). |

---

## 5. Directory Artifacts

Execution generates several output files within `--outdir`:

1.  `_synteny_plot.html`: Dynamic HTML representation mapping evolutionary relationships and boundary thresholds alongside compressed sequence gaps.
2.  `regions/*.bed`: Standardized tab-delimited text files denoting candidate bounding coordinates on target assemblies.
3.  `synterra_report.json`: Formatted trace detailing pipeline input arguments, genomic targets handled, memory limits utilized, and internal exit codes.
4.  `_tree.nwk`: Distance-matrix representation mapping output loci clustering via phylogenetic divergence.
5.  `downloaded_genomes/easy_mode_genomes/assembly_quality.tsv`: Generated log file reporting contiguity indices for remote assemblies pulled dynamically.

---

## 6. Known Execution States (Troubleshooting)

### Memory Halts During MMseqs2 Profile Building

**Symptom:** The pipeline crashes reporting a SIGKILL (137) during `LOCATE_GENE`.
**Cause:** Host system ran out of RAM or standard /tmp space mapping the index database.
**Resolution:** Supply a high-memory custom profile or restrict internal MMseqs chunk sizing using `--mmseqs_split_memory_limit '1G'`. Ensure `/tmp` is not constrained. 

### Prolonged Duration at ITERATIVE_SEARCH

**Symptom:** Nextflow interface lists `ITERATIVE_SEARCH` as processing for >40 minutes without state change on specific instances.
**Cause:** SynTerra has evaluated a target sequence containing massive tandem duplication arrays with complex inter-block spacing and is executing exhaustive localized `tblastn` matrices sequentially over all generated blocks.
**Resolution:** This is an expected mathematical calculation block. Await resolution unless explicitly violating the runtime timeout limit. You can dynamically read block state progress by performing a `tail -f` command against the target `.command.log` output inside the Nextflow `work/` directory hash.
