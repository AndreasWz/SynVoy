# Installing SynVoy

The short version is in the [README](../README.md#install). This page gives the details: what is installed, how to check it, macOS, containers and clusters.

## What you need

| | |
|---|---|
| **Operating system** | Linux on x86-64, where SynVoy is developed and tested. For a Mac, see [macOS](#macos). On ARM Linux the environment cannot be created, because BLAST+ 2.17.0 has no build for it. |
| **conda** | [Miniforge](https://github.com/conda-forge/miniforge) is recommended: it comes with mamba, which builds SynVoy's tool environment faster than conda. With Miniconda or Anaconda, conda builds it, or install mamba too (`conda install -n base -c conda-forge mamba`). |
| **Memory** | 8 GB for the launcher's default profile, 16 GB or more for `-profile standard`, and more for vertebrate genomes. |
| **Disk** | Room for two conda environments built from `environment.yml` (`synvoy_env`, and the one Nextflow builds in `work/conda/`, see [step 4](#4-the-first-run-builds-the-tool-environment)), and for the genomes of each run: in Easy Mode the downloaded genomes stay in `work/` and are copied to the results folder. The quickstart's run needs about 5 GB besides the environments. |
| **git** | To download SynVoy and to update it. |
| **Network** | For the installation and the first run; in Easy Mode for downloading the query and the genomes; in both modes for NCBI's taxonomy and for gene names in the figures. Without it, the genomes are ordered by sequence similarity instead. |

You do not need to install Java, Nextflow or Python yourself: the conda environment brings them.

## 1. Download SynVoy

```bash
git clone https://github.com/AndreasWz/SynVoy.git
cd SynVoy
```

## 2. Create the environment

```bash
./install.sh
```

`install.sh` creates the conda environment `synvoy_env` from `environment.yml`, with mamba if you have it. Then it checks the tools: Nextflow, Java 17 or newer, MMseqs2, BLAST, Prodigal, Augustus, miniprot, MAFFT, IQ-TREE, NCBI datasets, Entrez Direct, Python 3.10 to 3.12 with ete3 and numpy, and a test alignment with parasail. When `synvoy_env` exists already, it asks whether to rebuild it.

At the end, `install.sh` prints example commands that call `nextflow run main.nf ... -profile standard` directly, for the human gene LY6E. Do not start with them: check the installation as in step 3, and start runs with `./run_synvoy.sh`, which runs its checks and uses the memory-safe default profile. A vertebrate gene such as LY6E needs a machine with at least 32 GB of memory ([QUICKSTART.md](QUICKSTART.md#going-further)).

To create the environment by hand instead:

```bash
mamba env create -f environment.yml     # or: conda env create -f environment.yml
```

The environment pins the versions SynVoy is tested with:

| Package | Version |
|---|---|
| Nextflow | 26.04.3 |
| Java (OpenJDK) | 17 to 25 |
| Python | 3.12.13 |
| MMseqs2 | 18.8cc5c |
| BLAST+ | 2.17.0 |
| miniprot | 0.18 |
| Augustus | 3.5.0 |
| Prodigal | 2.6.3 |
| MAFFT | 7.526 |
| IQ-TREE | 3.1.2 |
| samtools | 1.23 |
| NCBI datasets | 18.30.1 |
| parasail-python | 1.3.4 |

Entrez Direct, ete3 and numpy are not pinned. Every run records the versions of its tools in the report, under `run_settings.tool_versions`.

The experimental search steps that use a protein language model or predicted structures (`--enable_plm_search`, `--enable_structural_search`, both off by default) need more packages. `environment-ml.yml` holds the packages of `environment.yml` and these. Under the conda profiles Nextflow builds the environment of the search (`ITERATIVE_SEARCH`, the pipeline step they run in) from it instead of `environment.yml` when one of the two steps is on (step 4), so nothing has to be installed by hand. Neither `install.sh` nor the Docker image installs them. How to run these steps on a GPU cluster: the GPU example profiles in [USAGE.md](USAGE.md#5-running-on-hpc--slurm).

## 3. Check the installation

```bash
./run_synvoy.sh --check-only
```

The launcher activates `synvoy_env` and checks that conda or mamba, Nextflow and Java 17 or newer are there, that SynVoy's configuration can be read, and that parasail can be loaded. When something is missing, it stops and says what is missing, for most checks with the command that fixes it. With `--check-only` it runs nothing else.

## 4. The first run builds the tool environment

The pipeline's steps do not run in `synvoy_env`. With the conda profiles, which include the launcher's default `auto,low_mem` and `standard`, Nextflow builds a second environment from the same `environment.yml` when the first step starts, and runs every step in it. This needs the network and takes from a few minutes to an hour on a slow shared file system (Nextflow waits up to 3 hours). When a language-model or structure search step is on, the search (`ITERATIVE_SEARCH`) runs in an environment built from `environment-ml.yml` instead, and every other step keeps the first one. Switching one of these steps on or off under `-resume` therefore reruns the search and the steps after it, and the steps before it stay cached.

Nextflow keeps this environment in `work/conda/` and reuses it, so later runs start at once. When you delete `work/`, the environment is built again. To keep it in one place for all your runs, set `NXF_CONDA_CACHEDIR` to a folder of your choice before you start SynVoy.

## Updating

```bash
cd SynVoy
git pull
```

When `environment.yml` has changed, run `./install.sh` again and let it rebuild `synvoy_env`. Nextflow notices the change by itself and builds a new tool environment on the next run; it names the environment by the file's whole text, so any change of the file (a comment too) gives a new environment, and a `-resume` across it reruns every step. `./run_synvoy.sh` prints the version and the commit it runs and, when it can reach GitHub, warns when your copy is behind. With Docker or Singularity, build the image again after an update.

## macOS

SynVoy is developed and tested on Linux.

**Intel Macs:** the environment can be created for Intel macOS, so you install SynVoy as described above.

**Apple Silicon (M1 and later):** the environment cannot be created for arm64, because there the pinned Augustus and samtools need different versions of the htslib library. Two ways around this:

- **Create the Intel environment,** which macOS runs through Rosetta 2. Install Rosetta 2 once (`softwareupdate --install-rosetta`), then set `CONDA_SUBDIR` in every terminal where you install or run SynVoy, so that Nextflow builds its tool environment the same way:

  ```bash
  export CONDA_SUBDIR=osx-64
  ./install.sh
  ./run_synvoy.sh --check-only
  ```

- **Use Docker** with an image built for Intel: `docker build --platform linux/amd64 -t synvoy-local:latest .` (see [Docker](#docker)).

SynVoy's own tests run on Linux only, so none of these macOS routes is tested.

## Docker

The image holds every tool of `environment.yml` (not the extra packages of `environment-ml.yml`), so the machine that runs it needs Docker, Nextflow 25.10 or newer and Java 17 or newer, but no conda environment for the tools.

```bash
docker build -t synvoy-local:latest .
```

Build it in a fresh copy of SynVoy or before your first run: Docker sends the whole folder to the build, and `work/` and `results/` grow to many GB.

The `docker` profile runs every step in the image `synvoy-local:latest`; use `--docker_container` to give another name. It has the full memory settings of `standard`, so the search asks for at least 10 GB, and the memory tiers cannot be combined with it ([USAGE.md](USAGE.md#memory-tiers)). If you have `synvoy_env` as well, start through the launcher, which keeps its checks:

```bash
./run_synvoy.sh -profile docker --mode easy --query_id P01501 --max_genomes 5 --outdir results/melittin_docker
```

Without conda, the launcher stops, so call Nextflow yourself from the `SynVoy` folder:

```bash
nextflow run main.nf -profile docker --mode easy --query_id P01501 --max_genomes 5 --outdir results/melittin_docker
```

An `NCBI_API_KEY` set in your shell does not reach the tasks inside Docker. To pass it in, give a configuration file with the line `docker.envWhitelist = 'NCBI_API_KEY'` with `-c`.

## Singularity and Apptainer

Singularity needs the image as a file. Build the Docker image on a machine with Docker, as above, and save it:

```bash
docker save synvoy-local:latest -o synvoy.tar
```

Then convert it. The image sets the `PATH` to its tools, and the folder of Augustus' gene models, as image variables, which Singularity keeps (it does not run the Docker entrypoint):

```bash
singularity build synvoy.sif docker-archive://synvoy.tar    # or: apptainer build ...
```

Copy `synvoy.sif` to the machine that runs SynVoy and give its path with `--docker_container`:

```bash
nextflow run main.nf -profile singularity --docker_container /path/to/synvoy.sif \
    --mode easy --query_id P01501 --max_genomes 5 --outdir results/melittin_singularity
```

The Singularity route has not been tested end to end yet. If a step says that a tool is not found, check the image with `singularity exec synvoy.sif mmseqs version`. An image built from an earlier `Dockerfile` lacks these variables: build it again.

## Clusters (SLURM)

On a cluster you have two choices for the tools:

- **conda:** create `synvoy_env` with `./install.sh` and use the profile `hpc_conda`. Nextflow builds the tool environment on the first run, as above, which can take long on a shared file system.
- **Singularity:** build `synvoy.sif` as above and use the profile `hpc_singularity` with `--docker_container /path/to/synvoy.sif`. Then you need only Nextflow and Java 17 or newer on the cluster, for example from its module system.

A Nextflow of your own must be version 25.10 or newer; the environment pins 26.04.3. How to submit the run, and the script `scripts/slurm_submit.sh`, are described in [USAGE.md](USAGE.md#5-running-on-hpc--slurm).

## Troubleshooting

- **`conda env 'synvoy_env' not found`:** the launcher needs the environment; run `./install.sh`.
- **The environment does not solve on a Mac:** see [macOS](#macos).
- **`parasail` cannot be imported:** a Python virtualenv probably shadows the conda environment (VS Code activates one by itself). Run `deactivate` and start again. If parasail is really missing, rebuild the environment with `./install.sh`.
- **`No module named 'cgi'` when ete3 is imported:** the environment has Python 3.13, which ete3 does not support. Rebuild it from the current `environment.yml`, which pins Python 3.12.
- **Creating the environment takes very long:** use mamba (Miniforge has it). Nextflow builds the environment with mamba when a `mamba` is on the `PATH` of the shell that starts SynVoy (`command -v mamba` prints a path), and with conda otherwise (`conda.useMamba` in `nextflow.config`). Nextflow gives up after 3 hours.

More problems and their fixes: [USAGE.md § 8](USAGE.md#8-troubleshooting).
