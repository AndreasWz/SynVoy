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
2. **Java 17+** and **Nextflow 23.04+** on PATH:
   ```
   java -version   # expect 17 or higher
   nextflow -v     # expect 23.04 or higher
   ```
3. **The SynVoy conda environment activated:**
   ```
   conda activate synvoy
   ```
4. **(Optional but recommended) an NCBI API key** to raise the download
   rate limit from 3/s to 10/s:
   ```
   export NCBI_API_KEY=your_key_here
   ```
   Get one free at https://www.ncbi.nlm.nih.gov/account/.

## Run the pipeline (Easy Mode)

From the repository root:

```bash
nextflow run main.nf \
    -profile standard \
    --mode easy \
    --query_id P01501 \
    --max_genomes 5 \
    --outdir results/quickstart_melittin \
    --auto_params false --multi_profile false
```

Notes on the flags:

- `--query_id P01501` — the UniProt accession for *Apis mellifera*
  melittin (70 aa preprotein). SynVoy pulls the FASTA, reads the
  species from the entry, and fetches the matching *Apis mellifera*
  reference genome from NCBI.
- `--max_genomes 5` — cap the NCBI taxonomy walk at 5 related bee
  assemblies. Fewer genomes → faster run, less informative tree;
  more genomes → the opposite.
- `--auto_params false --multi_profile false` — disables the Ollama
  LLM parameter advisor for reproducibility. You can enable it later
  (see [Going further](#going-further)).

If your laptop has <16 GB RAM, add `-profile laptop_safe` instead of
`-profile standard`.

## What you should see

The pipeline will log progress for each stage. Successful completion
looks like:

```
executor > local (~40)
[...]
Completed at: ...
Duration    : ~25m
Succeeded   : 40
```

Under `results/quickstart_melittin/` you should find roughly:

```
results/quickstart_melittin/
├── synvoy_report.json          # structured summary (see below)
├── locus_1_tree.nwk            # phylogenetic tree of the orthologs
├── synteny_block_locus_1_synteny_plot.html   # interactive synteny plot
├── synteny_block_locus_1_tree.html           # interactive tree
├── regions/
│   ├── <species_1>.fna.regions.bed
│   ├── <species_2>.fna.regions.bed
│   └── ...
└── qc/
```

The exact species set depends on what NCBI has available when you run
(the taxonomy walk picks the best-quality assemblies at run time), so
filenames and scaffold IDs will vary. What to verify qualitatively:

- Several (typically 3–5) `*.regions.bed` files under `regions/`, one
  per target species returned.
- Each BED file has at least one row with a score > 0.3.
- `locus_1_tree.nwk` has one leaf per target species that produced a
  candidate, plus the query.
- `synvoy_report.json` summary shows `total_annotations > 0` and
  `total_goi_annotations > 0`.
- `staging_diagnostics.empty` is `false`.

If everything is zero, something went wrong — jump to
[Troubleshooting](#troubleshooting) below.

## Inspect the report

```bash
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
tooltips show gene names and identities.

## Going further

- **Try a different gene.** Pass any UniProt accession via `--query_id`
  (e.g. `Q16553` for human LY6E — the README's default example).
- **Switch to Pro Mode** to supply your own genomes and a home GFF —
  see [USAGE.md § 1](USAGE.md#pro-mode). Pro Mode is reproducible
  (no taxonomy walk) and is what the paper's benchmarks use.
- **Enable the LLM parameter advisor** for harder queries:
  `--auto_params true` (requires Ollama + `gemma3:4b` pulled locally,
  or `GOOGLE_API_KEY` set for the hosted Gemma endpoint).

## Troubleshooting

If the run fails, read the error message — SynVoy's CLI tools emit
"what broke / why / try this" guidance. Common issues and fixes are
documented in [USAGE.md § 8](USAGE.md#8-troubleshooting), especially:

- NCBI download hangs / rate-limited →
  [USAGE.md § Easy Mode fails to download genomes](USAGE.md#easy-mode-fails-to-download-genomes).
- Pipeline finishes with `0 annotations` →
  [USAGE.md § Pipeline finishes with ...](USAGE.md#pipeline-finishes-with-synvoy_reportjson-showing-0-annotations--0-regions).
- `parasail` import error → install parasail or use ssearch36.
- `-resume` reruns everything → check that paths and params didn't change.
