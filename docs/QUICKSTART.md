# SynVoy — Quickstart (melittin example)

A complete end-to-end run of SynVoy on a small, well-understood problem:
finding **melittin** orthologs across related bee genomes. This is the
canonical onboarding example.

No local test data is required — SynVoy fetches the query and all
genomes from UniProt and NCBI automatically (Easy Mode).

**Runtime:** ~20–30 min on a 16 GB laptop (CPU-only, no GPU required).
Most of that is NCBI genome downloads on first run; subsequent runs
with `-resume` finish in a few minutes.
**Disk:** ~5 GB for the downloaded genomes + Nextflow `work/` cache.
**Network:** required (UniProt + NCBI).

## What you need

Before starting:

1. **A working SynVoy installation.** See the
   [README](../README.md) for conda setup, and come back here.
2. **The SynVoy conda environment activated:**
   ```
   conda activate synvoy_env
   ```
   The env provides Nextflow (≥25.10; 26.04 supported), OpenJDK (17 or newer) and
   Python 3.10–3.12 — you should not need to install any of them separately. Sanity check after activation:
   ```
   nextflow -version   # provided by the env
   java -version       # provided by the env
   ```
3. **(Optional but recommended) an NCBI API key** to raise the download
   rate limit from 3/s to 10/s:
   ```
   export NCBI_API_KEY=your_key_here
   ```
   Get one free at https://www.ncbi.nlm.nih.gov/account/.

## Run the pipeline (Easy Mode)

From the repository root:

```bash
./run_synvoy.sh \
    --mode easy \
    --query_id P01501 \
    --max_genomes 5 \
    --outdir results/quickstart_melittin
```

`./run_synvoy.sh` is the single entry point: it checks your version and
environment, then launches with a memory-safe profile by default — so this
command works on a laptop with no extra tuning.

Notes on the flags:

- `--query_id P01501` — the UniProt accession for *Apis mellifera*
  melittin (70 aa preprotein). SynVoy pulls the FASTA, reads the
  species from the entry, and fetches the matching *Apis mellifera*
  reference genome from NCBI.
- `--max_genomes 5` — cap the NCBI taxonomy walk at 5 related bee
  assemblies. Fewer genomes → faster run, less informative tree;
  more genomes → the opposite.

On a powerful machine, add `-profile standard` for full speed and the
phylogenetic tree: `./run_synvoy.sh -profile standard --mode easy ...`.
For the full parameter and profile reference (LLM advisor, presets, HPC),
see [USAGE.md](USAGE.md) and [PARAMETERS.md](PARAMETERS.md).

## What you should see

The pipeline logs one line per stage. Successful completion ends with a
summary like this (numbers will differ):

```
[OK  ] PIPELINE                 Pipeline completed successfully
Run Summary
Duration:          ~25m
Tasks Completed:   60
Results Directory: results/quickstart_melittin
Generated Outputs
  ✓ synvoy_report.json          (analysis summary)
    → 1 high-confidence + 2 medium-confidence GOI ortholog annotation(s) (+…) across 4 genome(s).
  ✓ synteny_block_locus_1_synteny_plot.html (interactive visualization)
  ✓ regions/                      (5 region BED file(s))
  ✓ logs/                         (task logs for … process(es))
```

A line `Tasks FAILED (ignored): N` in that summary means a rescue or paralog-check
step failed for some genome; the main results are still complete (see
[USAGE.md § 8](USAGE.md#the-run-ends-with-tasks-failed-ignored-n)).

Under `results/quickstart_melittin/` you should find roughly:

```
results/quickstart_melittin/
├── synvoy_report.json                        # the adjudicated result (see below)
├── synteny_block_locus_1_anchor_grid.html    # main figure (species × gene grid) + .svg
├── synteny_block_locus_1_synteny_plot.html   # interactive ribbon plot
├── synteny_block_locus_1_*_with_fragments.*  # the same figures with fragment hits drawn
├── synteny_block_locus_1_gene_positions.html # where each gene really sits
├── plot_inputs_synteny_block_locus_1/        # per-species gene models (.gff) + evidence table (.homology.tsv)
├── rescue/locus_1/                           # gene models found by the rescue passes (if any)
├── locus_1_tree.nwk                          # ortholog tree — PLACEHOLDER on the default profile (see note)
├── regions/
│   ├── <species_1>.fna.regions.bed
│   ├── <species_1>.fna.scores.tsv
│   └── ...
├── qc/  intermediate/  logs/
└── query/  home_genome/  downloaded_genomes/ # Easy Mode: what was fetched
```

> **Note — the phylogenetic tree is skipped by default.** The launcher's default
> profile is `auto,low_mem`, and `low_mem` sets `--skip_tree` to save memory, so
> `locus_1_tree.nwk` is a placeholder (`(GOI_placeholder:0.0);`) and
> `synteny_block_locus_1_tree.html` may be absent. Re-run with `-profile standard`
> for the real MAFFT + IQ-TREE phylogeny.

The exact species set depends on what NCBI has available when you run
(the taxonomy walk picks the best-quality assemblies at run time), so
filenames and scaffold IDs will vary. What to verify qualitatively:

- Several (typically 3–5) `*.regions.bed` files under `regions/`, one
  per target species returned.
- Each BED file has at least one row with a score > 0.3.
- *(only with `-profile standard`)* `locus_1_tree.nwk` has one leaf per
  target species that produced a candidate, plus the query. On the default
  laptop profile this file is a placeholder (the tree step is skipped).
- `synvoy_report.json` has a `summary.headline` that names at least one
  ortholog annotation, and `staging_diagnostics.empty` is `false`.
- The anchor grid shows the melittin column (red) filled for the close
  relatives (*Apis cerana*, *Apis florea*). Further out, expect `AMBIGUOUS`
  cells: a conserved neighbourhood whose sequence is not shown to be melittin.
  That is the honest answer — bumblebees, for example, have lost the gene.

If everything is zero, something went wrong — jump to
[Troubleshooting](#troubleshooting) below.

## Inspect the report

```bash
# The one-line answer, then the adjudicated calls
jq -r '.summary.headline' results/quickstart_melittin/synvoy_report.json
jq -r '.goi_dedup.records[] | [.genome, .chrom, .start, .end, .confidence, .identity] | @tsv' \
    results/quickstart_melittin/synvoy_report.json

# Top-level summary
jq '.summary' results/quickstart_melittin/synvoy_report.json

# Per-genome annotation counts
jq '.annotations.per_genome' results/quickstart_melittin/synvoy_report.json

# Diagnostic staging counts (always populated, useful on failure too)
jq '.staging_diagnostics.match_counts' results/quickstart_melittin/synvoy_report.json
```

## Open the interactive plots

The plots are standalone HTML — open them in a browser:

```bash
# Linux
xdg-open results/quickstart_melittin/synteny_block_locus_1_synteny_plot.html

# macOS
open results/quickstart_melittin/synteny_block_locus_1_synteny_plot.html
```

You should see one row per target species, with the melittin locus
and its flanking genes arranged in a (mostly) conserved order. Hover
tooltips show gene names and identities. For the overview across all
species, open `synteny_block_locus_1_anchor_grid.html` instead; what every
file is for is in [OUTPUT.md](OUTPUT.md).

## Going further

- **Try a different gene.** Pass any UniProt accession via `--query_id`
  (e.g. `Q16553` for human LY6E).
- **Switch to Pro Mode** to supply your own genomes and a home GFF —
  see [USAGE.md § 1](USAGE.md#pro-mode). Pro Mode is reproducible
  (no taxonomy walk) and is what the paper's benchmarks use.
- **Enable the LLM parameter advisor** for harder queries by providing
  a cloud API key:
  ```bash
  # Google Gemini (free tier available at aistudio.google.com)
  export GOOGLE_API_KEY=your_key
  ./run_synvoy.sh ... --auto_params true

  # OpenAI
  export OPENAI_API_KEY=your_key
  ./run_synvoy.sh ... --auto_params true --llm_provider openai
  ```
  Without an API key, `--auto_params true` falls back to built-in
  heuristics (still useful, just not LLM-quality).

  On current Nextflow the estimate is **advisory**: SynVoy cannot change its
  own parameters mid-run, so it prints the suggested flags (and saves them
  to `intermediate/estimate_params/estimated_params.json`) for you to re-run
  with.

## Troubleshooting

If the run fails, read the error message — SynVoy's CLI tools emit
"what broke / why / try this" guidance. Common issues and fixes are
documented in [USAGE.md § 8](USAGE.md#8-troubleshooting), especially:

- NCBI download hangs / rate-limited →
  [USAGE.md § Easy Mode fails to download genomes](USAGE.md#easy-mode-fails-to-download-genomes).
- Pipeline finishes with `0 annotations` →
  [USAGE.md § Pipeline finishes with ...](USAGE.md#pipeline-finishes-with-synvoy_reportjson-showing-0-annotations--0-regions).
- **OOM / `cannot fit database into ... not enough memory to keep dbreader/write in memory`**
  → you're running large (vertebrate-scale) targets on a small machine. Add
  `-profile standard,laptop_safe` (16 GB RAM) or `-profile standard,low_mem`
  (8 GB RAM). See [USAGE.md § Memory tiers](USAGE.md#memory-tiers-combine-with-an-execution-backend).
- `parasail` import error → the run **aborts**, by design: Smith-Waterman is
  load-bearing for divergent genes, so SynVoy refuses to search without it rather
  than quietly returning fewer hits. The usual cause is a `.venv` shadowing the
  conda env (VS Code activates one automatically) — run `deactivate` and relaunch.
  Otherwise re-create `synvoy_env`; Python must be `<3.13` because `ete3` is not
  Python 3.13-ready yet.
- `-resume` reruns everything → check that paths and params didn't change.
