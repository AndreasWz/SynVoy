# SynVoy quickstart: melittin

This page walks you through one complete run: finding the orthologs of honeybee **melittin**, a short venom peptide, in the genomes of related bees. SynVoy downloads the query and all the genomes itself (Easy Mode), so you need no data of your own.

| | |
|---|---|
| **Time** | About 20 to 30 minutes on a 16 GB laptop, most of it downloading genomes. The very first run of SynVoy also builds its tool environment, once (a few minutes to an hour). A rerun with `-resume` takes a few minutes. |
| **Disk** | About 5 GB for the genomes and Nextflow's `work/` folder, plus the conda environments. |
| **Network** | Needed, for UniProt and NCBI. |

## Before you start

1. **Install SynVoy** as described in the [README](../README.md#install), then check the installation from the `SynVoy` folder:

   ```bash
   ./run_synvoy.sh --check-only
   ```

   It ends with `--check-only: environment is healthy. You're ready to run.` You do not need to activate the conda environment yourself: the launcher does it.

2. *(Optional)* **An NCBI API key** raises NCBI's limit from 3 to 10 requests per second, which helps when downloads are slow or fail. Get one free in your NCBI account (https://www.ncbi.nlm.nih.gov/account/) and set it before the run:

   ```bash
   export NCBI_API_KEY=your_key_here
   ```

## Run SynVoy

From the `SynVoy` folder:

```bash
./run_synvoy.sh --mode easy --query_id P01501 --max_genomes 5 --outdir results/quickstart_melittin
```

- `--query_id P01501` is the UniProt accession of *Apis mellifera* melittin (the 70-aa precursor). SynVoy downloads the protein, reads its species from the entry and downloads the *Apis mellifera* reference genome from NCBI.
- `--max_genomes 5` downloads up to 5 related genomes to search, chosen by taxonomy. Fewer genomes run faster.
- `--outdir` is the results folder. A relative path like this one lies inside the `SynVoy` folder.

The launcher prints the SynVoy version, checks the installation and starts Nextflow. Because no `-profile` is given, it uses the memory-safe profile `auto,low_mem`: one task of each step at a time, 4 GB for the search, MMseqs2 at sensitivity 7.0 and no gene tree. On a machine with more memory, add `-profile standard` for full sensitivity (9.5) and the gene tree. The profiles are described in [USAGE.md](USAGE.md#2-execution-profiles).

## What you should see

While the run goes on, Nextflow shows a table of its steps. At the end, SynVoy prints a summary like this one (your numbers will differ):

```
[OK  ] PIPELINE                 Pipeline completed successfully
Run Summary
Duration:          24m 37s
Tasks Completed:   58
Results Directory: results/quickstart_melittin
Generated Outputs
  ✓ synvoy_report.json          (analysis summary)
    → 2 ortholog(s) (1 high, 1 medium confidence); 2 need curation (2 gene order only); across 4 genome(s).
    → no GOI call placed in 1 genome(s)
  ✓ synteny_block_locus_1_anchor_grid.html (main figure)
  ✓ synteny_block_locus_1_synteny_plot.html (ribbon plot)
  – locus_1_tree.nwk             (placeholder, no tree inferred: --skip_tree, which the low_mem profile sets)
  ✓ regions/                      (4 region BED file(s) in 1 home-locus folder(s))
  ✓ logs/                         (task logs of 41 task(s))
```

The first `→` line is the report's headline: how many orthologs SynVoy found, how many candidates need curation, and in how many genomes it placed any model of the gene (`LOW` ones included). The lines after it count the genomes in which SynVoy placed no model of the gene, not even a `LOW` one: first those without a strong match of it, then those with a strong match that no synteny rule accepted. The report names them (`summary.goi_absent_genomes`, `summary.goi_found_but_not_syntenic`). On the default profile the tree step is skipped, and the summary names `locus_1_tree.nwk` a placeholder.

If a rescue or ownership (paralog check) step failed for a genome, the run still finishes, and the report says so: a sentence in `warnings` and a row in `failed_steps` (see [USAGE.md](USAGE.md#a-step-failed-but-the-run-finished-failed_steps-in-the-report)).

The results folder holds, among other things:

```
results/quickstart_melittin/
├── synvoy_report.json                          the result (see below)
├── synteny_block_locus_1_anchor_grid.html      main figure, interactive (and .svg)
├── synteny_block_locus_1_anchor_grid_print.svg the same at print size (legend text and cell table: .legend.txt, .source.tsv)
├── synteny_block_locus_1_synteny_plot.html     ribbon plot of the neighbourhoods
├── synteny_block_locus_1_gene_positions.html   where each gene lies in its genome
├── synteny_block_locus_1_*_with_fragments.*    ribbon plot and grid with fragment models (only if there are any)
├── locus_1_tree.nwk                            gene tree (a placeholder on the default profile)
├── plot_inputs_synteny_block_locus_1/          per genome: gene models (.gff) and evidence table (.homology.tsv)
├── regions/locus_1/                            per genome: candidate regions (.regions.bed) and their scores (.scores.tsv)
├── query/  home_genome/  downloaded_genomes/   what Easy Mode downloaded
└── qc/  intermediate/  logs/                   genome quality checks, intermediate files, task logs
```

`locus_1` is the place of melittin in the home genome. When a gene matches well in several places of the home genome, each place becomes a locus of its own (at most 5). The downloaded genomes are named by their assembly accession. What every file contains: [OUTPUT.md](OUTPUT.md).

## Read the report

The report is a JSON file. [jq](https://jqlang.github.io/jq/) makes it easy to read (`sudo apt install jq` or `brew install jq`):

```bash
cd results/quickstart_melittin

# The one-line answer
jq -r '.summary.headline' synvoy_report.json

# The HIGH, MEDIUM and AMBIGUOUS calls, one record per target gene: genome, scaffold, start, end, label, class, % identity
jq -r '.goi_dedup.records[] | [.genome, .chrom, .start, .end, .confidence, .goi_class, .identity] | @tsv' synvoy_report.json
```

A gene found from several places of the home genome is one record, and a gene-order-only call in a genome that has a counted call is left out ([One row per gene](OUTPUT.md#one-row-per-gene)). A `HIGH` or `MEDIUM` record is an ortholog, unless its `goi_class` is `paralog_not_goi` (it belongs to another gene of the home genome) or `identity_coverage_decoupled` (a high identity over only a short part of the query). An `AMBIGUOUS` record needs curation, unless its `goi_class` is `paralog_not_goi`: then it belongs to another home gene, and the headline counts it as `of another home gene`. `LOW` models, with no significant match, are not listed here; the headline only counts them. Its first number counts the orthologs alone. The labels are explained in the [README](../README.md#results) and, in full, in [OUTPUT.md](OUTPUT.md#confidence-levels).

What to expect: orthologs in the closest relatives, usually other *Apis* species. Further out, the neighbourhood of melittin is often still there while the gene is not found; such a genome has no call, or an `AMBIGUOUS` candidate that is not counted. The species SynVoy downloads depend on what NCBI holds when you run, so your genomes and numbers can differ from anyone else's.

## Look at the figures

The figures are standalone HTML files; open them in a browser:

```bash
xdg-open synteny_block_locus_1_anchor_grid.html     # Linux
open synteny_block_locus_1_anchor_grid.html         # macOS
```

The anchor grid has one row per genome and one column per gene: melittin in red in the middle, its neighbours on either side. The outline of each gene shows its confidence, and hovering shows the details. The ribbon plot (`*_synteny_plot.html`) draws the same neighbourhoods with every gene at its real size, and ribbons between matching genes. The README shows [an example](../README.md#what-it-does) of the print version.

## Going further

The commands below start from the `SynVoy` folder again (`cd ../..`).

- **Try another gene:** give another UniProt or NCBI protein accession to `--query_id`. Vertebrate genomes are about ten times larger than insect genomes: for a gene such as human LY6E (`Q16553`), use `-profile standard --max_genomes 3` on a machine with at least 32 GB of memory.
- **Use your own genomes** in Pro Mode, with a home annotation if you have one: see [USAGE.md](USAGE.md#pro-mode). Pro Mode does not depend on what NCBI holds, so a rerun searches the same genomes.
- **Rerun after a change:** add `-resume`, and Nextflow reuses every step whose inputs did not change.

## If something goes wrong

Read the error message first: SynVoy's messages say what broke, why, and what to try. The common problems are described in [USAGE.md § 8](USAGE.md#8-troubleshooting), in particular:

- **Downloads hang or fail:** set an NCBI API key (above) and rerun with `-resume`. See [Easy Mode fails to download genomes](USAGE.md#easy-mode-fails-to-download-genomes).
- **No ortholog in a genome where you expected one:** see [No call in a genome where you expected one](USAGE.md#no-call-in-a-genome-where-you-expected-one).
- **The report step stops with `zero annotation and zero region files`:** see [The report step stops](USAGE.md#the-report-step-stops-zero-annotation-and-zero-region-files).
- **Out of memory** (exit code 137): the default profile already uses the smallest memory tier. Lower MMseqs2's memory further (`--mmseqs_split_memory_limit 1500M`), search fewer or smaller genomes, or use a machine with more memory. When MMseqs2 says `Cannot fit databases into ...`, the split limit is too small for that genome: raise it instead. See [Out of memory](USAGE.md#out-of-memory-exit-code-137).
- **`Process requirement exceeds available memory`:** a step asks for more memory than your machine reports. See [USAGE.md](USAGE.md#process-requirement-exceeds-available-memory).
- **`parasail` cannot be imported:** the launcher stops before the run, on purpose, because every search needs parasail: the Smith-Waterman search, which finds divergent genes, and the label of each call. The usual cause is a Python virtualenv that shadows the conda environment (VS Code activates one by itself): run `deactivate` and start again. Otherwise rebuild the environment with `./install.sh`.
- **`-resume` reruns everything:** `-resume` continues the last run started in the `SynVoy` folder, and a changed path or parameter reruns the steps that use it. `nextflow log`, called in the `SynVoy` folder with `synvoy_env` active, lists the runs, and `-resume <run name>` continues a particular one.
