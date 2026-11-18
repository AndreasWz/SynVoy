# SynVoy live demo - budding yeasts

A pro-mode run small enough to execute in front of an audience, on real
chromosome-level genomes with real annotation.

## The genomes

| Role | Species | Assembly | Size | Seqs | Why |
|---|---|---|---|---|---|
| **home** | *Saccharomyces cerevisiae* S288C | `GCF_000146045.2` (R64) | 12.1 Mb | 16 | the best-annotated eukaryotic genome |
| target | *Naumovozyma castellii* | `GCF_000237345.1` | 11.2 Mb | 10 | post-WGD, closest of the three |
| target | *Kluyveromyces lactis* | `GCF_000002515.2` | 10.7 Mb | 7 | pre-WGD |
| target | *Lachancea thermotolerans* | `GCF_000142805.1` | 10.4 Mb | 8 | pre-WGD, farthest |

All four are **chromosome-level and annotated**, which matters for a demo:

- few contigs → the synteny signal is clean, not fragmented across scaffolds
- home GFF present → real flanking genes, no Augustus/Prodigal prediction step
- target GFFs present → recovered models carry real gene names (`TargetGene`),
  so the plot and report are readable rather than a wall of `LOC…` ids
- the three targets straddle the *Saccharomyces* whole-genome duplication, so
  divergence increases across them instead of being uniform

Total input ≈ 76 MB.

## The query

| File | Gene | Length | Character |
|---|---|---|---|
| `STE2.faa` | **STE2** (`D6VTK4`), α-factor pheromone receptor | 431 aa | 7-TM GPCR, single copy, fast-evolving - the interesting case |
| `PGK1.faa` | **PGK1** (`P00560`), phosphoglycerate kinase | 416 aa | highly conserved, single copy - the control that should always work |

STE2 is the default because a fast-evolving receptor is where a sequence-only
search starts to struggle across genera, which is the case SynVoy exists for.
PGK1 is there as the safe fallback if a live run needs to be guaranteed.

## Running it

```bash
demo/fetch_data.sh          # once - downloads into local_data/demo/ (gitignored)
demo/run_demo.sh            # default query STE2  -> results/demo_ste2
demo/run_demo.sh PGK1       # the control        -> results/demo_pgk1
demo/run_demo.sh STE2 -resume
```

Anything after the query name is passed straight through to Nextflow.

The script strips an active virtualenv from the environment before launching. A
`.venv` that VS Code activates shadows the conda environment, and without `parasail`
the launcher stops rather than search without Smith-Waterman.

## The equivalent raw command

```bash
./run_synvoy.sh -profile low_mem -c demo/demo.config --mode pro \
    --query          local_data/demo/query/STE2.faa \
    --home_genome    local_data/demo/home/Saccharomyces_cerevisiae.fna \
    --home_gff       local_data/demo/home/Saccharomyces_cerevisiae.gff \
    --home_species   "Saccharomyces cerevisiae" \
    --target_genomes 'local_data/demo/genomes/*.fna' \
    --target_gffs    'local_data/demo/target_gffs/*.gff' \
    --outdir         results/demo_ste2
```

## What to show

1. `results/demo_ste2/*_anchor_grid.html` (or its `.svg`): the payoff figure, with
   rows = species, columns = genes, and the STE2 column filled across all four.
2. `results/demo_ste2/synvoy_report.json`: read `summary.headline` and the records
   under `goi_dedup.records`; `target_gene` names the gene of the target's own
   annotation that each call lies on. Do not read `total_raw_search_hits`: it is always
   0 in this version, because the report does not receive the raw search hits.
3. If the search bridged a gap between flanking genes in some genome, its log says so:
   `grep -r "collinearity: bridged" results/demo_ste2/logs/`. That is the
   gap-bridging slide happening live.

## Run time

A few minutes on a workstation (about 4 minutes on 12 cores), so it fits a 10-minute
slot. Time it once on the machine you present from. `--anchor_max_places 5` makes the
search visit fewer windows and runs faster.

`low_mem` is deliberate: 4 GB per search task and one task at a time, so the demo
cannot run out of memory in front of an audience. On 11 Mb genomes that limit never
binds.

## Why the overrides live in `demo.config`

`low_mem` has one setting that is wrong for a *demo*: it sets `skip_tree = true`.
`demo/demo.config` turns the gene tree back on. It does **not** raise
`mmseqs_sensitivity` from 7.0 to 9.5: on genomes this small, with neighbours this well
conserved, the higher sensitivity costs minutes and gave the same STE2 calls.

The override is in a config file to keep the command short. Passing it as a flag
(`--skip_tree false`) works too.

## What to expect

One STE2 ortholog in each of the three targets, each on the gene that the target's own
annotation calls STE2 (`target_gene` in the record). The labels come from one rule for
every call: the E-value of the home protein in the place the neighbours predict, whether
another *S. cerevisiae* protein explains the model as well, and how much of the query the
model covers ([ALGORITHM.md § 7.8](../docs/ALGORITHM.md#78-confidence)). Identity alone
decides no label. Run the demo once before you present and note the headline it prints,
because labels and counts can change between versions.

## Things to know before presenting

- **The figures draw the search's own models,** not those of the rescue passes. A call
  that exists only as a rescue model is in the report and in
  `results/<run>/rescue/locus_N/`, but not in the anchor grid.
- **The LOW items in the report are hit fragments and weak models**, not candidates
  worth discussing. The default figures hide the fragment models; when a run has some,
  the `*_with_fragments` figures show them as well.
- **The search works in small neighbourhoods.** Distances are counted in gene spacings
  of the home locus (about 2.4 kb in yeast), so a block is a neighbourhood of some tens
  of kb. The first lines of each genome's GFF in `plot_inputs_*/` give the distances
  used (`##synvoy-distances`).

## Other demos

- `demo/run_demo_insect.sh` runs the same kind of pro-mode demo on *Drosophila*: home
  *D. melanogaster*, targets *D. simulans*, *D. yakuba* and *D. pseudoobscura*, default
  query Defensin (`Adh` as the conserved control). Fetch the data once with
  `demo/fetch_data_insect.sh`.
- `demo/run_uniprot.sh <UniProt ID>` lets the audience pick a *D. melanogaster*
  protein: it downloads the protein from UniProt and searches it in *D. pseudoobscura*
  and *Anopheles gambiae*, with the insect demo's data (`demo/fetch_data_insect.sh`
  downloads both). The script's header lists IDs that are known to work.
- `demo/demo_yeastscale.config` fixes the distances in base pairs, an experiment kept
  only to reproduce it. By default SynVoy derives the distances for each home locus and
  target genome.
