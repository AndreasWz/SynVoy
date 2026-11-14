# Running SynVoy

How to start SynVoy and how to fit a run to your data and your machine. A first run, step by step: [QUICKSTART.md](QUICKSTART.md). Installing: [INSTALL.md](INSTALL.md). Every parameter: [PARAMETERS.md](PARAMETERS.md). Every output file: [OUTPUT.md](OUTPUT.md). What SynVoy computes: [ALGORITHM.md](ALGORITHM.md).

## Contents

- [Starting a run](#starting-a-run)
1. [Modes](#1-modes)
2. [Execution profiles](#2-execution-profiles)
3. [Settings for your gene](#3-settings-for-your-gene)
4. [Parameters](#4-parameters)
5. [Running on HPC / SLURM](#5-running-on-hpc--slurm)
6. [Output files](#6-output-files)
7. [Resuming and caching](#7-resuming-and-caching)
8. [Troubleshooting](#8-troubleshooting)

---

## Starting a run

Start every run with the launcher, from the `SynVoy` folder:

```bash
./run_synvoy.sh [Nextflow options] [SynVoy parameters]
```

For example:

```bash
./run_synvoy.sh -profile standard --mode easy --query_id P01501 --max_genomes 5 --outdir results/melittin
```

Before it starts Nextflow, the launcher:

- prints the SynVoy version and commit (for example `SynVoy v1.0.1` and `main@1281e47`) and, when it can reach GitHub, warns when your copy is behind;
- activates the conda environment `synvoy_env` and checks that Nextflow, Java 17 or newer and the Smith-Waterman library parasail are there. When one is missing, it stops and says what is missing, for most checks with the command that fixes it;
- warns when a Python virtualenv is active, because a virtualenv can hide the environment's parasail from the pipeline's steps;
- adds `-profile auto,low_mem` when you give no `-profile` (see [§ 2](#2-execution-profiles)).

Then it runs `nextflow run main.nf` with everything you passed. `./run_synvoy.sh --check-only` runs the checks and nothing else.

Two things to know:

- **Paths are read from the `SynVoy` folder.** The launcher changes into its own folder before it starts Nextflow, so a relative path (`--query`, `--home_genome`, `--target_genomes`, `--outdir`, `-work-dir`, `-params-file`, `-c`) is read from there, not from the folder you are in. For data elsewhere, give absolute paths.
- **A `-profile` you give replaces the default completely.** To add a preset to the default, name the whole set, for example `-profile auto,low_mem,preset_short_peptide`.

You can also call `nextflow run main.nf` yourself, from the `SynVoy` folder. Then the launcher's checks do not run, and Nextflow uses the profile `standard` when you give none. You need this only in a few cases: a container setup without conda ([INSTALL.md](INSTALL.md#docker)), a SLURM controller job ([§ 5](#5-running-on-hpc--slurm)), and a run without parasail, because the launcher stops whenever parasail cannot be loaded. A run without parasail needs more than one switch: see [Smith-Waterman cannot run](#smith-waterman-cannot-run). Nextflow's own options that belong before `run`, such as `-log`, cannot be passed through the launcher either.

---

## 1. Modes

SynVoy has two modes. **Easy Mode** (`--mode easy`, the default) downloads the genomes from NCBI. **Pro Mode** (`--mode pro`) uses files you give it.

### Easy Mode

```bash
./run_synvoy.sh --mode easy --query_id P01501 --max_genomes 5 --outdir results/melittin
```

Give the query in one of three ways:

| Parameter | What it is |
|---|---|
| `--query_id` | A protein accession, from UniProt (`P01501`) or NCBI (`XP_006565763.1`). SynVoy downloads the sequence and reads the species from the record. A gene symbol also works together with `--home_species` (`--query_id Melt --home_species "Apis mellifera"`): SynVoy then takes the UniProt entry of that gene in that species, a reviewed one when there is one. |
| `--query` | A protein FASTA file of your own. Needs `--home_species`. |
| `--query_seq` | The protein sequence itself, as text. Needs `--home_species`. |

Then SynVoy downloads the reference genome of the home species and the genomes of related species, with their annotations:

| Parameter | Default | What it does |
|---|---|---|
| `--home_species` | read from the accession | The home species, e.g. `"Apis mellifera"`. Overrides the species of the accession. |
| `--max_genomes` | `0` | How many related genomes to search. `0` takes three for each taxonomic level from genus to class, so usually 12. |
| `--target_species` | none | The species to search, comma-separated (`"Bombus terrestris,Nomia melanderi"`), instead of species chosen by taxonomy. |
| `--assembly_ranking` | `hybrid` | Which assembly to take when a species has several. Every rule prefers NCBI's reference or representative assembly, then a RefSeq assembly, then the higher assembly level. After that, `hybrid` prefers fewer contigs and scaffolds and a higher N50, `counts` fewer chromosomes, scaffolds and contigs, and `nstats` a higher N80 and N50. |

Downloaded assemblies of low quality are left out or flagged:

| Parameter | Default | What it does |
|---|---|---|
| `--bad_quality_policy` | `auto` | What happens to a low-quality assembly. `auto` leaves out one that SynVoy chose and keeps one you chose (the home assembly, or a species named in `--target_species`), with a warning in the console, the report and the figures. `drop`, `keep` and `ask` apply one rule to all. |
| `--bad_min_n50` | `5000` | An assembly whose N50 is below this is of low quality. |
| `--bad_max_contigs`, `--bad_max_scaffolds` | `500000` | An assembly with more contigs or scaffolds is of low quality. These apply only when NCBI's record gives the count. |

Chromosome-level and complete assemblies always pass. Every target genome, downloaded or your own, is checked again on its FASTA file (N50 below `--bad_min_n50`, more sequences than `--bad_max_contigs`); a gzipped target is not checked. `--qc_fail_policy auto` keeps a genome that fails this check when you chose it (Pro Mode, or `--target_species`), flagged, and drops it when SynVoy chose it; `drop` and `keep` apply one rule to all.

A free NCBI API key, from your NCBI account, raises NCBI's limit from 3 to 10 requests per second. Set it before you start: `export NCBI_API_KEY=your_key`. Under the `docker` profile it must be passed into the container ([INSTALL.md](INSTALL.md#docker)).

### Pro Mode

```bash
./run_synvoy.sh --mode pro \
  --query /data/melittin.faa \
  --home_species "Apis mellifera" \
  --home_genome /data/Apis_mellifera.fna \
  --home_gff /data/Apis_mellifera.gff \
  --target_genomes /data/targets/ \
  --outdir results/melittin
```

| Parameter | What it is |
|---|---|
| `--query` | **Required.** Your gene's protein, as a FASTA file. |
| `--home_genome` | **Required.** The genome your gene comes from, as a FASTA file (`.fna`, `.fa` or `.fasta`). |
| `--home_species` | The home species, e.g. `"Apis mellifera"`. SynVoy uses it to order the target genomes by relatedness. Without it, SynVoy takes the name of the `--home_genome` file, without its extension, as the species (`Apis_mellifera.fna` is read as *Apis mellifera*), with a warning. It stops when that name is a placeholder such as `genome` or `reference`, or an assembly accession. `--allow_unknown_species true` runs without a species; the targets are then ordered without the home species' taxonomy. |
| `--home_gff` | *Recommended.* The home genome's gene annotation (GFF3). Without it, SynVoy predicts the home genome's genes (Augustus; Prodigal for prokaryotes) and transfers gene models from an annotated target genome, which is less reliable than a curated annotation. |
| `--target_genomes` | The genomes to search: a folder (every `.fna`, `.fa` or `.fasta` file in it, plain or gzipped; other files are ignored), a quoted glob (`"genomes/*.fna"`) or a comma-separated list (`"a.fna,b.fna"`). Without it, SynVoy analyses the home genome only. |
| `--target_gffs` | *Optional.* Annotations of the target genomes; see below. |

Pro Mode downloads no genomes, but it orders the targets by taxonomy: from ETE3's copy of the NCBI taxonomy (`~/.etetoolkit/taxa.sqlite`, downloaded on first use), or by looking up the lineages of assembly accessions in the file names at NCBI. When neither works, it orders them by their similarity to the home genome. The report records which way was used (`genome_order`).

**File names.** The name of a genome file, without its extension, becomes the genome's name in sequence IDs, GFF attributes and output file names. When target genomes are given, SynVoy stops at the start if the name of a query, home or target file contains whitespace or one of ``| ; = , ' " ` $ & ( ) < > * ? ! [ ] { } # % \``, or if two targets differ only in their extension (`Apis.fa` and `Apis.fna`). Use letters, digits, `.`, `_` and `-`, for example `Apis_mellifera.fna`.

**The query is one protein.** It may be protein or DNA; DNA is translated to its longest open reading frame. Lower case is read as upper case, and whitespace, digits, alignment gaps (`-` and `.`) and a trailing `*` are removed. An internal `*` is removed with a warning, and any other character that is not an amino-acid code stops the run. When the file holds several sequences, only the first is used, and the log says so. A query shorter than 20 aa stops the run (`--min_query_length`, `0` switches the check off), and one shorter than 40 aa runs with a warning.

**Target annotations.** `--target_genomes` takes FASTA files only, so annotations of the targets are given separately with `--target_gffs`: a folder, a quoted glob or a comma-separated list of `.gff`, `.gff3` or `.gtf` files, plain or gzipped. Each annotation is matched to its genome by file name (`cow.fna` and `cow.gff`) or by assembly accession (`GCF_002263795.1.fna` and `GCF_002263795.1_ARS-UCD1.2_genomic.gff`). Annotating only some targets is fine. An annotation that matches no genome stops the run; `--allow_unmatched_target_gffs true` turns that into a warning. Easy Mode downloads the annotations itself and ignores `--target_gffs`.

With target annotations, a result also names the annotated gene it lies on: `target_gene` and `target_product` in the report, `TargetGene`, `TargetProduct` and `TargetID` in the gene models. With a home GFF this is information only: the annotations of the targets change no block, score or label. Without a home GFF, gene models of the closest annotated target are transferred to the home genome and can become its neighbouring genes. Only a target whose genome file name ends in `.fna` (not gzipped) is used for this, so give an annotated target that name.

---

## 2. Execution profiles

A profile says where the tasks run and where their tools come from. Choose one with `-profile <name>`.

| Profile | Tasks run | Tools from | Use it for |
|---|---|---|---|
| `standard` | locally | conda | A workstation or server. Nextflow uses it when you call it without the launcher and without `-profile`. |
| `auto` | locally | conda | The same as `standard`. It also restates the defaults of SynVoy's automatic settings ([§ 3](#3-settings-for-your-gene)), so a run computes the same with either. |
| `conda` | locally | conda | The same as `standard`. |
| `laptop_safe`, `low_mem` | locally | conda | Laptops: memory tiers, see [below](#memory-tiers). |
| `docker` | locally | the Docker image `synvoy-local:latest` | A machine with Docker ([INSTALL.md](INSTALL.md#docker)). |
| `docker_max` | locally | the Docker image | A dedicated Linux machine: every task gets all CPUs and all memory but 1 GB (the search all but 2 GB), one task at a time. The memory is read from `/proc/meminfo`; without it, as on macOS, the profile assumes 8 GB. |
| `singularity` | locally | a Singularity image | A machine with Singularity or Apptainer. Give the image with `--docker_container /path/to/synvoy.sif` ([INSTALL.md](INSTALL.md#singularity-and-apptainer)). |
| `slurm` | as SLURM jobs | the compute nodes' `PATH` | A cluster where the tools are installed. Set its partition and account first ([§ 5](#5-running-on-hpc--slurm)). |
| `hpc_conda` | as SLURM jobs | conda | A cluster with conda and mamba ([§ 5](#5-running-on-hpc--slurm)). |
| `hpc_singularity` | as SLURM jobs | a Singularity image | A cluster with Singularity ([§ 5](#5-running-on-hpc--slurm)). |

Two more kinds of profile are added to one of these:

- `legacy_distances` sets the search distances to fixed numbers of base pairs, the same for every genome, instead of following the spacing of the genes ([CHANGELOG.md](../CHANGELOG.md)). Use it only to compare with runs made that way: `-profile standard,legacy_distances`.
- The presets `preset_short_peptide`, `preset_tandem_family`, `preset_paralog_discrimination` and `preset_single_copy` set a few values for one kind of gene ([§ 3](#3-settings-for-your-gene)).

The profiles `lrz_ai` and `lrz_ai_container` are examples for one GPU cluster ([§ 5](#5-running-on-hpc--slurm)). The profiles `test`, `test_melettin` and `test_tetramorium` are developer fixtures that need data the repository does not hold.

When you combine profiles, a later one overrides what an earlier one sets. So combine one profile from the table with, at most, `legacy_distances` and one preset.

### Memory tiers

`laptop_safe` and `low_mem` lower the memory of the heaviest steps so that SynVoy fits on a laptop. Through the launcher, a run without `-profile` uses `auto,low_mem`.

| Profile | For a machine with | Search task memory | MMseqs2 split limit | MMseqs2 sensitivity | Gene tree |
|---|---|---|---|---|---|
| `standard` or `auto` | 16 GB or more | `auto`: at least 10 GB, more for large genomes | 8G | 9.5 | yes |
| `laptop_safe` | about 16 GB | 8 GB | 4G | 8.0 | yes |
| `low_mem` | about 8 GB | 4 GB | 2G | 7.0 | no |

- Both tiers also give each search task 1 CPU, give the gene search in the home genome (`LOCATE_GENE`) 1 CPU and 2 GB, and run one task of each step at a time. `low_mem` builds no gene tree (`--skip_tree true`); add `--skip_tree false` to build it anyway.
- A lower sensitivity finds a few hits fewer for very divergent genes, so results of different tiers can differ.
- The split limit (`--mmseqs_split_memory_limit`) bounds how much of a target database MMseqs2 holds in memory at once. It is the value to lower when MMseqs2 says `cannot fit database into ...` or `not enough memory`.
- Under `standard`, the search task asks for at least 10 GB, and more when the genomes searched at the same time are large. A retry after a memory kill asks for twice, then three times as much, but on a local machine never more than the machine has. The console prints the request (`ITERATIVE_SEARCH memory: ...`).
- The tiers are local conda profiles themselves: they set the local executor and conda, and switch Docker and Singularity off. Combine them only with `standard`, `auto` or `conda`. On Docker, Singularity or SLURM, set the values as parameters instead, for example `--iterative_search_memory '8 GB' --mmseqs_split_memory_limit 4G`.
- One step that orders the genomes (`PHYLO_SORT`) asks for 8 GB in every profile, and Nextflow does not start a task that asks for more memory than the machine has. On a machine that reports less than 8 GB, see [Process requirement exceeds available memory](#process-requirement-exceeds-available-memory).

---

## 3. Settings for your gene

After SynVoy has found your gene in the home genome, it adapts some settings to it. It looks at the query (its length, cysteines, low-complexity stretches, a signal anchor), at the gene's exons and at its copies in the home genome, and applies every **preset module** whose trigger fires: for example `short_query` for a query shorter than 100 aa, or `local_family` when the gene has more than one copy at its home locus. Each module sets a few values by a fixed rule, except `shared_domain`, `cys_framework` and `home_context`, which are only recorded. One rule labels every call ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)) and no module changes it, but a module can change the evidence it reads: `signal_anchor` also measures completeness on the core of the protein, without its signal peptide or GPI tail. The modules, their triggers and what they set: [PARAMETERS.md](PARAMETERS.md#modules).

**What was applied.** The console names the modules (`[INFO] modules (auto), applied: ...`) and the values their rules set. The file `intermediate/locate_gene/effective_params.json` records the whole result: `modules_applied`, the final `settings`, for every value its `source` (`default`, `module:<id>`, `preset:<name>` or `user_set`) and, for every value that is not the default, what set it (`derivation`). The report keeps the same record under `run_settings`.

**Your own values win.** A value you set yourself, on the command line, in a params file, in a `-c` config or in a profile, is kept even where a module would set another one; `user_set_kept` in `effective_params.json` lists where that happened. SynVoy recognises your value by its difference from the shipped default, so a value you set equal to the default can still be changed by a module. To keep every default, use `--auto_apply_preset false`.

**Choosing the modules yourself.**

| Setting | What it does |
|---|---|
| `--modules auto` | The default: every module whose trigger fires. |
| `--modules 'auto,-dispersed_family'` | The same, without one module (or with one more: `'auto,+cys_framework'`). Write the value in quotes. |
| `--modules 'short_query,local_family'` | Exactly these modules, whether their trigger fires or not. |
| `--modules none` | No module. |
| `--modules legacy` | A chooser picks one of the four presets below from the gene's hits in the home genome and applies it; the modules' triggers are only recorded. |
| `--preset_override preset_short_peptide` | Applies this preset on top of the modules. |
| `--auto_apply_preset false` | No module and no preset: the defaults, and the values you set. |

**Presets.** The four presets hold only the values that have a reason for that kind of gene:

| Preset | For | Sets |
|---|---|---|
| `preset_short_peptide` | queries shorter than 100 aa | the coverage bands of a model's status: complete from 0.5 of the query, fragment below 0.3 |
| `preset_tandem_family` | several copies at the home locus | weaker relatives of the gene are not used as neighbouring genes either (`goi_relative_evalue 1e-3`) |
| `preset_paralog_discrimination` | the query aligns well at several places of the home genome | fewer candidate regions per genome (`adaptive_score_floor_frac 0.45`, `adaptive_score_floor_abs 0.08`, `adaptive_max_regions 4`) |
| `preset_single_copy` | everything else | nothing: the defaults |

A preset can also be given as a profile (`-profile auto,low_mem,preset_short_peptide`). Then its values count as values you set, so they win over the modules.

The length of the query is checked before any of this. To search a peptide shorter than 20 aa, set `--min_query_length 0` (or a lower length) at the start; no module or preset changes it.

---

## 4. Parameters

Give a SynVoy parameter with two dashes (`--max_genomes 10`) and a Nextflow option with one (`-profile`, `-resume`, `-work-dir`). Write switches as words: `--skip_tree false`. Every parameter, its default and when to change it: [PARAMETERS.md](PARAMETERS.md).

For a run with many parameters, write them into a params file and pass it with `-params-file`:

```yaml
# melittin.yaml
mode: pro
query: /data/melittin.faa
home_species: Apis mellifera
home_genome: /data/Apis_mellifera.fna
home_gff: /data/Apis_mellifera.gff
target_genomes: /data/targets/
outdir: results/melittin
```

```bash
./run_synvoy.sh -profile standard -params-file melittin.yaml
```

A configuration file passed with `-c` changes Nextflow's settings, for example the memory or time of one step:

```groovy
// my.config
process {
    withName: 'COMPUTE_TREE' { time = '12h' }
}
```

```bash
./run_synvoy.sh -profile standard -c my.config --mode easy --query_id P01501
```

What a run used is recorded with its results: the values that depend on the gene in `intermediate/locate_gene/effective_params.json` ([§ 3](#3-settings-for-your-gene)), the tool versions and resources in the report under `run_settings`.

---

## 5. Running on HPC / SLURM

On a cluster, Nextflow runs as one job, the controller, and submits every task as a job of its own. Choose how the tasks get their tools:

- **`hpc_conda`**: Nextflow builds the conda environment from `environment.yml` on the first run, with mamba, and reuses it ([INSTALL.md](INSTALL.md#4-the-first-run-builds-the-tool-environment)). On a slow shared file system this can take an hour. Set `NXF_CONDA_CACHEDIR` to keep the environment in one place for all your runs.
- **`hpc_singularity`**: every task runs in a Singularity image. Build `synvoy.sif` as in [INSTALL.md](INSTALL.md#singularity-and-apptainer) and give it with `--docker_container /path/to/synvoy.sif`.
- **`slurm`**: the tasks use the tools on the compute nodes' `PATH`. Before you use it, edit its partition (`normal`) and its placeholder account (`--account=your_account`) in `nextflow.config`. The `hpc_*` profiles set neither and use your cluster's defaults.

The controller needs Nextflow 25.10 or newer and Java 17 or newer, from your cluster's modules or from `synvoy_env` (`conda activate synvoy_env`).

**The submission script.** `scripts/slurm_submit.sh` submits the controller. Edit the lines at its top (query, species, number of genomes, profile) and its `#SBATCH` header (partition, account, time), then submit it:

```bash
sbatch scripts/slurm_submit.sh
```

- The script uses `hpc_singularity` by default. Add `--docker_container /path/to/synvoy.sif` to its `nextflow run` line, or set `PROFILE="hpc_conda"`.
- The controller must run as long as the whole run. Give it a time that covers the run: one search task may take up to 48 hours.
- The script calls `nextflow run main.nf` directly, so the launcher's checks do not run.
- It puts Nextflow's work folder on fast scratch (`$SCRATCH/work`, or `work/` in the folder you submit from when `SCRATCH` is not set) and the Singularity images in `~/.singularity/cache`.

**By hand.** Inside a job, or wherever your cluster lets a long-running controller run:

```bash
./run_synvoy.sh \
  -profile hpc_conda \
  --mode easy \
  --query_id P01501 \
  --max_genomes 10 \
  --outdir results/melittin_hpc \
  -work-dir "${SCRATCH}/work" \
  -resume
```

**GPU example profiles.** The search steps that use a protein language model or predicted structures (`--enable_plm_search`, `--enable_structural_search`) are experimental, off by default, and need the extra packages of `environment-ml.yml`. Two profiles show a setup for them on one GPU cluster; adapt them to your site.

- `lrz_ai_container` runs the whole run inside one container job, with the local executor. Its image must hold the packages of both `environment.yml` and `environment-ml.yml`; the repository's Dockerfile installs only the first.
- `lrz_ai` uses SLURM with conda; give your GPU partition with `--gpu_partition`, otherwise no GPU is requested. As under every conda profile, the tasks run in the environment Nextflow builds from `environment.yml`, which lacks these packages, and installing `environment-ml.yml` into `synvoy_env` does not change that environment. One way around it is to point `process.conda` (in a `-c` config) at the folder of an environment that holds both, such as `synvoy_env` after that install (`conda env list` shows its folder).

---

## 6. Output files

Everything is written to the folder given with `--outdir` (default `results`). Start with `synvoy_report.json`, whose `summary.headline` is the one-line answer, and with the figure `*_anchor_grid.html`. What every file contains, and how to read the report: [OUTPUT.md](OUTPUT.md).

At the end of a run, also a failed one, the task logs are collected in `<outdir>/logs/`: one folder per task that wrote any output, with `stdout.log`, `stderr.log` and `script.sh`. They are collected from the whole work folder, so the tasks of earlier runs in the same work folder are added too.

---

## 7. Resuming and caching

Nextflow keeps every task in its work folder (`work/` in the `SynVoy` folder, unless you give `-work-dir`). With `-resume`, a run reuses the tasks of the last run started in the `SynVoy` folder that have not changed:

```bash
./run_synvoy.sh [your parameters] -resume
```

- A changed input or parameter runs the tasks that read it again, and every task after them.
- A changed figure parameter runs only the figure step again. Exception: `--max_goi_per_genome` is part of the run's settings ([§ 3](#3-settings-for-your-gene)), so changing it runs the search again.
- One search task searches all target genomes of a home locus. Adding a genome reuses the steps on the home genome and runs the search again for every genome.
- After a failure, such as a network error in Easy Mode, the steps that had finished, such as the downloads, are reused.
- To continue an earlier run instead of the last one, list the runs with `nextflow log`, called in the `SynVoy` folder (the launcher starts every run there) with `synvoy_env` active, and give its name: `-resume <run name>`.

**Redrawing figures without rerunning.** `scripts/replot_run.sh <outdir>` draws every figure of a finished run again with your copy's figure code and writes them to `<outdir>/figures_redrawn_<date>/`; the calls stay those of the run. `scripts/replot.sh <outdir>/plot_inputs_synteny_block_locus_<N> <folder> [options]` redraws one home locus, with figure options of your choice. Run both with `synvoy_env` active.

**Cleaning up.** The work folder grows with every run. Delete it when you no longer need to resume: `rm -rf work/`. With conda profiles, the tool environment lives in `work/conda/` unless you set `NXF_CONDA_CACHEDIR`, so it is built again after that.

---

## 8. Troubleshooting

### The launcher stops before the run starts

Read its message: it names what is missing and, for most checks, the command that fixes it. The parasail message also offers `--allow_missing_smith_waterman true`; ignore that line, because the launcher stops on a missing parasail whatever you pass ([Smith-Waterman cannot run](#smith-waterman-cannot-run)).

- **`conda env 'synvoy_env' not found`:** run `./install.sh`.
- **`Java >=17 required` or `Java not found`:** Java ships in `synvoy_env`. Rebuild the environment with `./install.sh`.
- **`parasail is NOT importable`:** usually a Python virtualenv (`.venv`) shadows the conda environment; VS Code activates one in its terminal by itself. Run `deactivate` and start again. If parasail is really missing, rebuild the environment with `./install.sh`.

Problems with creating the environment, `mamba: command not found`, macOS and `No module named 'cgi'`: [INSTALL.md](INSTALL.md#troubleshooting).

### Smith-Waterman cannot run

The search stops when it should run Smith-Waterman and cannot: with `--sw_method auto` (the default) or `parasail` when `parasail` cannot be loaded in the task, and with `--sw_method ssearch36` when neither parasail nor the `ssearch36` program is installed (`environment.yml` does not install `ssearch36`). Smith-Waterman is what finds the most divergent orthologs, so SynVoy does not go on without it. Fix the environment (previous entry).

A run without parasail is a lesser run: its calls get no specificity terms, so the labels no longer test whether another home gene explains a call as well ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)), and the hull rescue and the paralog checks, which need parasail too, give no result. To run without it anyway, call `nextflow run main.nf` instead of the launcher and pass:

- `--enable_smith_waterman false`, or `--allow_missing_smith_waterman true` (which acts only with `--sw_method auto`); the search then finds fewer divergent orthologs. `--sw_method ssearch36` keeps Smith-Waterman when the `ssearch36` program of the FASTA package is installed;
- `--require_paralog_panel false`, because the step that builds the home-paralog panel (`BUILD_HOME_PARALOG_PANEL`) stops the run without parasail. Its failure is then recorded, and the report marks every call `ownership_not_checked`. `--disable_locus_ownership true` skips that check instead.

### Process requirement exceeds available memory

Nextflow does not start a task that asks for more memory than the machine has. The step that orders the genomes, `PHYLO_SORT`, asks for 8 GB in every profile, and a laptop sold with 8 GB usually reports a little less. Lower the request with a configuration file:

```bash
echo "process { withName: 'PHYLO_SORT' { memory = '6 GB' } }" > small_machine.config
./run_synvoy.sh -c small_machine.config --mode easy ...
```

The same message for another step: lower that step's memory parameter ([PARAMETERS.md](PARAMETERS.md)), or use the `low_mem` tier.

### Out of memory: exit code 137

A task was killed for memory, usually the search (`ITERATIVE_SEARCH`) or the gene search in the home genome (`LOCATE_GENE`) on a large genome. MMseqs2 builds a translated database of every genome, which for a genome of 1 to 3 Gb takes 3 to 9 GB.

- A killed task is retried, at most twice (a genome download three times); the step that scores the candidate regions (`CLUSTER_REGIONS`) is not retried. The search asks for more memory with each try, but on a local machine never more than the machine has.
- Through the launcher without `-profile`, the run already uses `low_mem`, the smallest tier. `laptop_safe` allows more memory, not less.
- When MMseqs2 says `cannot fit database into ...` or `not enough memory to keep dbreader/write in memory`, lower the split limit further: `--mmseqs_split_memory_limit 1G`.
- Otherwise search fewer or smaller genomes at once, or use a machine with more memory or a cluster ([§ 5](#5-running-on-hpc--slurm)).
- MMseqs2 also writes temporary files; check that the disk of the work folder has room.

### "process hasn't exited"

A task ran past its time limit on the local executor. Nextflow reports a time kill there with this message and does not retry it; the end summary names the step. Give the step more time with a configuration file ([§ 4](#4-parameters)), and run again with `-resume`. Most steps have a limit of 3 or 4 hours (`--default_time` sets the 4 hours), the search and the gene search in the home genome 48 hours.

### The search runs for a very long time

A search task can take hours for a large genome or for a gene with many relatives. To follow it, open the newest task log: `ls -lt work/*/*/.command.log | head -5`, then `tail -f` on the one you want.

- **Many candidate blocks.** In a genome rich in gene families (plants, genomes with old whole-genome duplications), every neighbouring gene hits many places, and two chance hits close together make a block. The log then says `Capping blocks: N -> 80`, and a large part of the genome is searched. `--anchor_max_places 5` builds blocks only at the five best places of each neighbouring gene. With `--modules` other than `legacy`, SynVoy does this by itself in a genome where the cap binds; the log says `genome axis: the block cap (80) bound at N blocks; rebuilding ...`. Lowering `--max_blocks_per_genome` (default 80) or raising `--min_block_genes` (default 2) also speeds it up.
- **Genes far apart.** The search distances follow the spacing of the genes around the home locus, so where genes lie far apart, as in many vertebrate genomes, every block and every search window is large. The log line `distances: home pitch ... -> block gap ..., window margin ... (DNA searched: ... per side)` shows the distances of each genome. If a run is too slow in such a genome, use `--anchor_max_places 5`, lower `--block_gap_genes`, or restore the fixed distances with `-profile <your profile>,legacy_distances`.

### "Query FASTA not found" or "Home genome not found"

The path is wrong as seen from the `SynVoy` folder: the launcher reads relative paths from there ([Starting a run](#starting-a-run)). Give an absolute path, or one relative to the `SynVoy` folder.

### "Could not determine home species from the query"

Easy Mode reads the species from the record of an accession. A FASTA file, a sequence or a record without an organism carries none, so give it: `--home_species "Genus species"`.

### Easy Mode fails to download genomes

NCBI limits how often it may be asked, and downloads can break.

- Set an NCBI API key: `export NCBI_API_KEY=your_key` ([Easy Mode](#easy-mode)).
- Run again with `-resume`: the steps that finished are reused.
- Check that NCBI is reachable from your machine, and see NCBI's status page.
- When SynVoy cannot find related genomes for your species, name the species to search with `--target_species`.

### "File name(s) with whitespace or special characters"

Rename the files as described under [Pro Mode](#pro-mode), then start again.

### "two records have the name ..." or "lies on the sequence ..., which this genome's FASTA does not hold"

A genome's sequence name is the first word of its FASTA header, whole: `>NODE_1_length_5000_cov_1.0 a description` is `NODE_1_length_5000_cov_1.0`. Every tool of the search names a sequence that way.

- **Two records have the name:** the task stops and names both headers and the file, because no tool can tell the two records apart. Rename one of them.
- **A hit or block lies on a sequence the FASTA does not hold:** a tool and the FASTA name the sequence differently. This check is a backstop and should not happen; please [report it](https://github.com/AndreasWz/SynVoy/issues) with the genome's first header lines.

### "the home sequence name ...", "the gene id ...", "share no sequence name", "The query id ..." or "The GOI name ..."

SynVoy stops on names that would otherwise lose results without a word. Each message names the id, the file and what to change.

- **A `|` in a flanking gene id** (an `ID=` of the home GFF, such as `ctg|1`): `EXTRACT_FLANKING` stops. SynVoy reads a gene id up to its first `|`, so `ctg|1` and `ctg|2` would be one gene. Rename the ids.
- **A `|` in a home sequence name, without a home GFF:** `ANNOTATE_GOI` stops. Predicted genes are named after their sequence, so every gene on `000000F|arrow` would read as one. Rename the sequences, for example with `sed -E '/^>/s/\|/_/g'`. With a home GFF, such names work.
- **A gene id with whitespace** (`ID=ORF 1`): the step that writes the gene's record stops, because a record is read by its first word, so `ORF 1` and `ORF 2` would be one gene.
- **The home FASTA and its GFF share no sequence name:** `ANNOTATE_GOI` stops. The GFF's first column must name the FASTA's records by the first word of their header (`NC_1.1`, not `NC_1.1|foo`). When only some sequences of the GFF are missing from the FASTA, their genes are left out, and the report says so (`home_gff_unmatched`, and a sentence in `warnings`).
- **A query id holding `|exon_`, `|frag_` or `|slide_`:** `NORMALIZE_QUERY` stops. SynVoy names the pieces of a protein that way, so the search would have no full-length query.
- **A query id that cannot name the gene:** `NORMALIZE_QUERY` stops on a `|` or on a character that GFF3 reserves (`,` `;` `=` `&` `%`) in the id of the query's header, and `ANNOTATE_GOI` on such characters, whitespace or the form `copy_<n>` in an Easy Mode `--query_id`. The id becomes the gene's name in every gene model. The common tagged headers are read as their accession: `sp|X|Y`, `tr|X|Y`, `ref|X|`, `gb|X|`, `emb|X|`, `dbj|X|` and `lcl|X` as `X`, `ENA|X|X.1` as `X.1`, and NCBI's CDS and protein downloads (`lcl|NC_1.1_cds_XP_1.1_1`, `lcl|NC_1.1_prot_XP_1.1_1`) as the protein's accession. Rename any other query in its header, for example with `sed -E '1s/[|,;=&%]/_/g' query.faa > renamed.faa`.

### No call in a genome where you expected one

An empty cell in the anchor grid means that no ortholog was placed there, not that the gene is absent. Check:

1. **The report's `rejected_candidates`:** a strong match outside the conserved neighbourhood is listed there as found but not syntenic ([README → Scope](../README.md#scope)).
2. **The genome's regions:** `regions/<locus>/<genome>.regions.bed`. When it is empty, the neighbouring genes were not found together in that genome. More neighbouring genes (`--n_flanking_genes 15`; the default is 10 per side) or a larger allowed gap between them (`--block_gap_genes`, default 7 gene spacings) can help where the region is broken up.
3. **The genome's quality:** a fragmented assembly can split the region or miss the gene. The report lists the low-quality genomes that were searched (`low_quality_genomes_searched`).
4. **The calls that are not counted:** an `AMBIGUOUS` call, a call the paralog check gave to another home gene (`paralog_not_goi`) and one with a high identity over a short part of the query (`identity_coverage_decoupled`) are in `goi_dedup.records` but not counted as orthologs. `LOW` models are not in the records: the headline counts them as "no significant match", and the search's LOW models are listed in `plot_inputs_synteny_block_locus_<N>/<genome>.homology.tsv` ([OUTPUT.md](OUTPUT.md#the-per-genome-table)).

### The report step stops: "zero annotation and zero region files"

`GENERATE_REPORT` found no gene models and no region files of the search. Either the search found nothing in any genome, or its results did not reach the report step.

1. The report is written before the step fails, but the files of a failed task are not copied to `--outdir`. Open it in the failed task's work folder (named in the error and in `.nextflow.log`): `jq '.staging_diagnostics' work/<hash>/synvoy_report.json`. It lists, per folder, how many files the step received.
2. If every count is zero, the search found nothing. Its log, `logs/ITERATIVE_SEARCH*/stdout.log` in `--outdir`, gives the hits per genome.
3. If files arrived but none matched (`*.gff`, `*.scores.tsv`), please [report it](https://github.com/AndreasWz/SynVoy/issues).

When a search with no result is what you expect, run the step's command again in its work folder with `--allow-empty` added to the `generate_report.py` call. The command is in the folder's `.command.sh`; run it with `synvoy_env` active.

### A step failed but the run finished: `failed_steps` in the report

The rescue, paralog-check and ownership steps may fail without stopping the run, because one failed rescue should not cost a whole result. Two kinds of failure still end it: a step that stops on one of SynVoy's own checks (exit code 3, for example a bad name or a failed tool), and a failed home-paralog panel (`BUILD_HOME_PARALOG_PANEL`), unless you set `--require_paralog_panel false`. Any other failure is recorded where the step's result would be, and the report names it: a sentence in `warnings` and a row in `failed_steps` (step, genome, locus, exit code, the last line the step wrote to its error output). A failed ownership check marks every call it would have checked `ownership_not_checked`: such a call was not compared with the other home genes and may belong to one of them. A failed rescue means that a gene only that rescue would have modelled is missing there. A step killed for memory or time is retried like any other, and after its retries it ends the run.

### A parameter I set seems to be ignored

Check these in order:

1. **Did a module or preset set it?** A value you set wins ([§ 3](#3-settings-for-your-gene)). Look it up in `intermediate/locate_gene/effective_params.json`: `source` says where each value came from, `derivation` what set it.
2. **Was it removed?** A removed parameter prints a warning and is ignored; a replaced one stops the run at the start and names its replacement. The list: [PARAMETERS.md](PARAMETERS.md#20-removed-parameters).
3. **Is it spelled right?** Check the parameter's name in [PARAMETERS.md](PARAMETERS.md).
4. **Is it an option of the figure script?** Figure options without a parameter of their own, such as `--print_grid_numbers`, go through `--plot_extra_args '...'` ([PARAMETERS.md](PARAMETERS.md#plot_extra_args)).

### "The output folder ... holds an earlier run's results"

A run overwrites the files it writes and leaves every other file in the output folder. When the folder already holds a report and the run is not resumed, the run says so at the start, and the end summary counts the files older than the run and does not list them as this run's results. Use an empty folder for a new run.

### `-resume` runs every step again

Nextflow reuses a task only when its inputs are the same. Any of these starts the task again, and every task after it:

1. An input file's path changed, even when its content did not.
2. A parameter changed, including a default in `nextflow.config`.
3. You moved the work folder or deleted parts of it, or started from another folder or with another `-work-dir`.
4. You switched profiles (`low_mem` and `standard` give different parameters).
5. `-resume` continued another run: it takes the last run started in the `SynVoy` folder. Name the run you mean ([§ 7](#7-resuming-and-caching)).

`./run_synvoy.sh ... -resume -dump-hashes` prints the hashes of every task, so you can see which input changed.

### GPU out of memory in the structural search

The experimental structural search (`--enable_structural_search`) predicts structures with ESMFold, whose memory grows with the square of the protein length. It folds every protein in overlapping pieces of 380 residues, and of 300 on a GPU with less than 10 GB or 150 below 6 GB. If it still runs out of memory, switch the step off (`--enable_structural_search false`, the default). `--structural_max_length` has no effect in this version.

### Finding out why a step failed

The error names the task's work folder (`Work dir: work/a7/8f3...`). In it:

```bash
cat work/a7/8f3*/.command.sh     # the command that ran
cat work/a7/8f3*/.command.err    # its error output
cat work/a7/8f3*/.command.out    # its normal output
cat work/a7/8f3*/.command.log    # both
```

The same logs, of every task that wrote any, are collected in `<outdir>/logs/` at the end of the run.

### Still stuck?

1. Run `./run_synvoy.sh --check-only`. It checks the environment without starting a run and prints your version.
2. Work through the [QUICKSTART](QUICKSTART.md), which is known to work.
3. [Open an issue](https://github.com/AndreasWz/SynVoy/issues) with the version line the launcher prints, your command, the end of `.nextflow.log` and the log of the failed step.
