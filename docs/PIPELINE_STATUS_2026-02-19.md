# SynTerra Status (2026-02-19)

## Run checked
- Run log: `.synterra_logs/run_20260219_192226.log`
- Command profile: `laptop_safe`, `--gene P60615`, home `Naja naja`, targets `Ophiophagus hannah,Bungarus multicinctus`
- State: pipeline completed (`[OK] PIPELINE`), output in `results_3snake_3ftx_kernfix_test/`

## Quick result check
- Plot exists: `results_3snake_3ftx_kernfix_test/synteny_block_locus_1_synteny_plot.html`
- Region BEDs exist:
  - `results_3snake_3ftx_kernfix_test/regions/GCA_000516915.1.fna.regions.bed`
  - `results_3snake_3ftx_kernfix_test/regions/GCA_023653725.1.fna.regions.bed`
- Tree exists: `results_3snake_3ftx_kernfix_test/locus_1_tree.nwk`
- GOI annotations exist in target GFFs:
  - `results_3snake_3ftx_kernfix_test/plot_inputs_synteny_block_locus_1/GCA_000516915.1.fna.gff`
  - `results_3snake_3ftx_kernfix_test/plot_inputs_synteny_block_locus_1/GCA_023653725.1.fna.gff`

## Current pipeline problems

### 1) Ground-truth anchor (TOP1MT) is not recovered in this run
- No `TOP1MT` string in home annotation:
  - `results_3snake_3ftx_kernfix_test/home_genome/home_genome/home_genome.gff`
- No `TOP1MT` string in run outputs:
  - `results_3snake_3ftx_kernfix_test/**`
- Practical impact: current synteny products cannot validate the expected `TOP1MT -> Ly6 -> 3FTx` story directly.

### 2) Flanking naming is still mostly technical IDs, not biology-friendly symbols
- Flanking BED uses mostly `gene-E2320_*` and many `hypothetical protein` entries:
  - `results_3snake_3ftx_kernfix_test/intermediate/flanking/synteny_block_locus_1.bed`
  - `results_3snake_3ftx_kernfix_test/home_genome/home_genome/home_genome.gff`
- Practical impact: difficult interpretation for end users and weak direct comparison to literature ground truth.

### 3) Region confidence remains weak (all medium, no high-confidence block)
- Region scores are around `0.43 - 0.48` only:
  - `results_3snake_3ftx_kernfix_test/regions/GCA_000516915.1.fna.regions.bed`
  - `results_3snake_3ftx_kernfix_test/regions/GCA_023653725.1.fna.regions.bed`
- Practical impact: region calls are still low/medium confidence for publication-quality inference.

### 4) GOI-overlap prioritization is over-inclusive in current output
- CLUSTER logs tag all top regions as `[GOI]`, including widely separated loci:
  - `work/2f/407ae973077b4cf32602d13c4de1cf/.command.err`
  - `work/31/1561c41d171db62ba68d5fa6eaebc0/.command.err`
- Practical impact: GOI overlap signal may be inflated by hit-derived intervals and can promote noisy regions.

### 5) Iterative search still has high noise and strong resource pressure
- Very high anchor hit counts and many seeded blocks:
  - `Parsed 445 hits`, `157 discrete search regions` (GCA_000516915.1)
  - `Parsed 420 hits`, `45 discrete search regions` (GCA_023653725.1)
- Repeated MMseqs prefilter failures, then low-memory retries:
  - `Error: Prefilter died` repeated many times in `work/2b/1942e88c87bfefdc95f34f485fda31/.command.err`
- Practical impact: unstable runtime behavior, high compute overhead, and elevated false-positive risk.

### 6) Report generation is not running (missing promised output)
- `synterra_report.json` is not present in output dir:
  - `results_3snake_3ftx_kernfix_test/`
- No `generate_report.py` task directory found in `work/` for this run.
- Practical impact: final summary artifact advertised in run footer is missing.

### 7) QC inconsistency between fetch-stage quality and runtime QC summary
- Fetch-stage marks both genomes as bad quality:
  - `results_3snake_3ftx_kernfix_test/downloaded_genomes/easy_mode_genomes/assembly_quality.tsv`
- Runtime QC marks both as PASS:
  - `results_3snake_3ftx_kernfix_test/qc/genome_qc_summary.json`
- Practical impact: contradictory quality signals make automatic decision logic unreliable.

### 8) Persistent Nextflow channel/wiring warnings remain
- Repeated warning: `The operator first is useless when applied to a value channel`
  - visible in `.synterra_logs/run_20260219_192226.log`
- Practical impact: not fatal, but indicates avoidable channel misuse and potential future wiring regressions.

## Net assessment
- Engineering status: pipeline now completes with current channel fixes.
- Scientific status: still not publication-ready for the stated TOP1MT/Ly6/3FTx biological question.
