# SynTerra: The Complete Manual

SynTerra is a pipeline designed for extreme evolutionary distances. When you are looking for highly divergent orthologs (e.g., fast-mutating immunity genes, complex multi-exonic micro-peptides, or deeply ancient paralogs), standard aligners lose their statistical signal. 

**SynTerra assumes that while a gene sequence may mutate wildly, its position relative to its neighboring genes (its syntenic architecture) remains highly conserved.** 

This manual is exhaustively structured to guide you through every feature, parameter, runtime profile, and execution mode SynTerra offers.

---

## Table of Contents
1. [Installation & Requirements](#1-installation--requirements)
2. [Runtime Profiles (HPC vs. Desktop)](#2-runtime-profiles-hpc-vs-desktop)
3. [Execution Modes](#3-execution-modes)
   - [Easy Mode (Auto-Pilot)](#easy-mode-auto-pilot)
   - [Pro Mode (Bring-Your-Own-Data)](#pro-mode-bring-your-own-data)
4. [The SynTerra Search Algorithm](#4-the-synterra-search-algorithm)
5. [Complete Parameter Reference](#5-complete-parameter-reference)
   - [Core I/O Parameters](#core-io-parameters)
   - [Synteny & Region Tuning](#synteny--region-tuning)
   - [Iterative Search Stringency](#iterative-search-stringency)
   - [Visual Plotting Configuration](#visual-plotting-configuration)
   - [Quality Control & Assembly Filtering](#quality-control--assembly-filtering)
6. [Output Artifacts](#6-output-artifacts)
7. [Advanced Troubleshooting](#7-advanced-troubleshooting)

---

## 1. Installation & Requirements

SynTerra requires a Linux/macOS environment (or WSL2 on Windows) strictly running **Nextflow >=22.10.1**. Everything else is completely containerized.

### Option A: Conda/Mamba (Recommended for Headless Servers)
The supplied `environment.yml` contains the exact lock-versions of `MMseqs2`, `miniprot`, `tblastn`, `Prodigal`, and tree-builders needed.

```bash
git clone https://github.com/AndreasWz/SynTerra.git
cd SynTerra
conda env create -f environment.yml
conda activate syntenyfinder
```

### Option B: Docker (Recommended for Desktops & Workstations)
Containerized execution ensures zero host-system contamination regarding binary paths.

```bash
docker build -t synterra-local:latest .
```
You simply append `-profile docker` to any of your run commands.

### Option C: Singularity/Apptainer (Recommended for HPC/Slurm)
HPC clusters typically forbid Docker daemons. Nextflow integrates perfectly with Singularity.

```bash
mkdir -p "$HOME/.singularity/cache"
export NXF_SINGULARITY_CACHEDIR="$HOME/.singularity/cache"
```
Append `-profile singularity` to your run commands.

---

## 2. Runtime Profiles (HPC vs. Desktop)

SynTerra executes mathematically heavy all-vs-all sequence translations on raw genomic bins. It will happily consume 128 CPU cores and 1TB of RAM if you let it.

### `-profile standard` (Default)
Runs entirely locally using the host's native `PATH`. Limits itself to standard memory allocations (`6GB` for search, `2GB` for prep).

### `-profile laptop_safe` (The Safest Route)
**Highly recommended for desktop users.** Restricts the pipeline to sequential, single-fork processing using strict RAM limitations (`-Xmx`). Prevents the pipeline from instantly thrashing a 16GB laptop when analyzing massive, 3Gb mammalian genomes.

### `-profile docker_max` (Saturate The Machine)
Automatically detects your host's maximum logical cores and available RAM. It will dynamically rewrite the internal allocation matrices to consume ~90% of your system resources across heavily threaded `MMseqs2` forks. **Only use this on dedicated cloud instances.**

### `-profile slurm`
Configures standard Nextflow executor arrays to dispatch every single genome alignment and annotation chunk as an independent Slurm batch job.

---

## 3. Execution Modes

SynTerra operates in two distinct philosophies:

### Easy Mode (Auto-Pilot)
You provide a single Protein ID (`--query_id Q16553`). You do not provide any fasta files. You do not provide any target genome paths.

1. SynTerra queries the UniProt/NCBI REST APIs.
2. It identifies the origin species (e.g., *Homo sapiens*).
3. It downloads the highest-quality Reference Genome and GFF annotation directly from NCBI servers.
4. It navigates the taxonomic tree to identify related species.
5. It mass-downloads their assemblies (`.fna.gz`), running quality-control filtering to drop shattered/scaffold-only genomes.
6. It runs the full syntenic search pipeline.

**Example Command:**
```bash
nextflow run main.nf \
  --mode easy \
  --query_id Q16553 \
  --max_genomes 5 \
  --outdir results_easy \
  -profile laptop_safe
```

### Pro Mode (Bring-Your-Own-Data)
You operate in a completely offline, isolated environment. You provide all `.fasta` genomes via glob paths. 

*Note: In Pro mode, your `--query` must be a local FASTA file, not an ID.*

**Example Command:**
```bash
nextflow run main.nf \
  --mode pro \
  --query input/my_gene.fasta \
  --home_genome input/Homo_sapiens.fna.gz \
  --home_gff input/Homo_sapiens_primary.gff \
  --target_genomes "input/raw_vertebrates/*.fna" \
  --outdir results_pro \
  -profile docker
```

---

## 4. The SynTerra Search Algorithm

The magic of SynTerra lies in its `ITERATIVE_SEARCH` state machine. It does not blindly BLAST your GOI against a 3Gb genome.

1. **Anchoring:** SynTerra finds the exact coordinates of your GOI sequence inside your `home_genome`.
2. **Flanking:** It sweeps outwards `N` genes to the left and right (`--n_flanking_genes 10`), extracting their protein sequences.
3. **Target Sweeping:** It maps these 20 flanking genes against the massive target genomes using high-speed `MMseqs2` profiling.
4. **Clustering:** It identifies regions on the target genome where a subset of these flanking genes successfully bound within close proximity (`--cluster_distance 500000`). This is a **Candidate Block**.
5. **Micro-Scanning:** SynTerra completely disregards the rest of the 3Gb target genome. It zooms in exclusively on that Candidate Block. It drops the stringency barriers, switches to exhaustive translation matrices (`tblastn`, `miniprot`), and deeply sequence-sweeps the empty genomic space looking for the highly mutated remnants of your GOI.

---

## 5. Complete Parameter Reference

SynTerra's behavior can be radically altered by overriding the default parameters located inside `nextflow.config`. You can append any of these directly to your `nextflow run main.nf` invocation as `--parameter_name value`.

### Core I/O Parameters

| Flag | Description | Mode Requirement |
| :--- | :--- | :--- |
| `--mode` | Define execution variant (`easy` or `pro`). | Required. |
| `--query` | Path to your local protein FASTA file. | **Pro** Mode. |
| `--query_id` | UniProt/NCBI ID (e.g., `P60615`). | **Easy** Mode. |
| `--home_genome` | Path to your local reference assembly (`.fna` or `.fna.gz`). | **Pro** Mode. |
| `--home_gff` | Path to the known annotation (`.gff`) for the home genome. Optional but massively improves performance. | **Pro** Mode. |
| `--target_genomes` | Standard glob-string matching local files (e.g., `"data/*.fna"`). Wrap in quotes. | **Pro** Mode. |
| `--target_species` | Comma-separated list of exact species to download (e.g., `"Felis catus,Canis lupus"`). | **Easy** Mode. |
| `--max_genomes` | Limits the maximum amount of auto-taxonomic genomes downloaded. Specify `0` for unlimited. | **Easy** Mode. |
| `--outdir` | Path to place the final plots, trees, and reports. Default: `results`. | Universal. |

### Synteny & Region Tuning

These parameters govern how large the "net" is when hunting for candidate blocks.

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--n_flanking_genes` | `10` | The number of genes extracted to the immediate Left and Right of your GOI on the home genome. Increase this (e.g., `30`) for highly fragmented target genomes where synteny blocks may be shattered. |
| `--cluster_distance` | `50000` | (50kb). The maximum allowable distance between two mapped flanking hits for them to be considered part of the same syntenic block. |
| `--region_padding` | `150000` | (150kb). The absolute physical space appended to the left and right boundaries of a discovered candidate block before deep-scanning. Increase to ensure you don't miss GOIs lying on the periphery. |
| `--min_synteny_score` | `0.6` | Requires at least 60% of the originally pulled `--n_flanking_genes` to be successfully mapped in the target cluster to warrant a deep GOI sweep. Drop to `0.3` for extremely ancient lineages. |

### Iterative Search Stringency

Control the ferocity of the inner sequence-aligners when they are zooming into a candidate block.

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--search_evalue` | `0.01` | The baseline e-value required for the algorithm to trust a matching hit. |
| `--aug_relaxed_identity_min` | `25.0` | In the deepest sweep matrices, what is the absolute floor percentage identity SynTerra will accept as a valid mutant GOI hit? |
| `--mmseqs_sensitivity` | `9.5` | Top-end `k-mer` seeding depth. Values above `10.0` provide highly marginal returns at massive RAM cost. |
| `--gap_min_alnlen` | `10` | The absolute minimum sequence length (in amino acids) recovered from an empty gap to count as a micro-exon. |

### Visual Plotting Configuration

SynTerra automatically draws massive, interactive, browser-based HTML ribbons demonstrating the exact evolutionary conservation mapping. You can tune the visual output layout perfectly for screenshots or publication.

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--plot_width` | `1500` | Global width (in pixels) of the dynamically rendered HTML canvas. |
| `--gap_threshold` | `50000` | (50kb). If a genomic gap contains absolutely no mapped genes and exceeds this size, it will be visually squished/compressed to save screen space. |
| `--gap_visual_size` | `3000` | The physical pixel-coordinate width drawn to represent a compressed gap. |
| `--flank_fallback_bp` | `1000000` | (1Mb). The absolute maximum distance outward from the centered GOI to draw flanking genes in the plot. Extremely distally mapped genes are discarded from the visual to prevent zoom-out collapse. |
| `--scale_bar_len` | `10000` | The biological distance (bp) represented by the physical scale bar drawn at the bottom of the plot. |

### Quality Control & Assembly Filtering

Specifically for `Easy Mode` where vast NCBI downloads are automated.

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--bad_quality_policy` | `ask` | How to handle garbage-tier assemblies. Options: `ask` (pauses pipeline for user input), `drop` (silently ignores them), `keep` (forces synteny alignment on shattering contigs). |
| `--bad_max_scaffolds` | `50000` | If an assembly has more than 50k disconnected fragments, it triggers the bad quality flag. |
| `--bad_min_n50` | `20000` | If the N50 contiguity is lower than 20kb, synteny mapping is mathematically nearly impossible. It triggers the bad quality flag. |

---

## 6. Output Artifacts

If you successfully run the pipeline, target your `--outdir`.

1. **`*synteny_plot.html`**: Open this in Chrome/Firefox/Safari. You can interact with it, hover over genes to see their precise chromosomal coordinates and `% identity` matches, and take screenshots for publications.
2. **`regions/` subdirectory**: If you prefer to construct your own plots in R (e.g., `gggenomes`) or Python (`Matplotlib`), you can parse these `.bed` files directly.
3. **`synterra_report.json`**: An incredibly dense, programmatic JSON payload summarizing exactly which genomes were probed, how many candidate blocks were constructed on each, how many blocks actually contained the GOI, and the pipeline execution duration.
4. **`_tree.nwk`**: A raw Newick-format topological branching matrix defining the evolutionary proximity between the identified genomic blocks.

---

## 7. Advanced Troubleshooting

### "Error executing process > LOCATE_GENE" / "No Space Left on Device"
SynTerra extracts massive temporary databases during MMseqs2 generation. If your local `/tmp` environment is small, the pipeline will crash. Change Nextflow's working directory to a large storage volume:
```bash
export NXF_WORK=/path/to/massive/storage/work
```

### The pipeline froze/hung at "ITERATIVE_SEARCH"!
**It is not frozen, it is doing its job.** 
When SynTerra encounters a heavily duplicated locus (like a massive 3FTx snake-venom toxin array), it generates a candidate block containing dozens of empty sequence gaps. It spawns internal `tblastn` and `miniprot` evaluators against *every single gap* to ensure it doesn't miss a micro-exon. 
*Do not cancel the run!* If you inspect `.synterra_logs/run_....log`, you will see it actively grinding through the genomic windows. 

### My Plot Labels Look Bad / "Gene-123445"
SynTerra inherits the naming conventions of the GFF annotations provided to it. 
If an NCBI assembly does not contain recognizable labels (e.g., relying entirely on uninformative generic locus tags like `LOC1234`), SynTerra will display exactly that.
For the absolute highest quality visualization, supply your own `--home_gff` where you have manually corrected the `Name=` or `gene=` field tags.
