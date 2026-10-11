# SynVoy live demo — budding yeasts

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
| `STE2.faa` | **STE2** (`D6VTK4`), α-factor pheromone receptor | 431 aa | 7-TM GPCR, single copy, fast-evolving — the interesting case |
| `PGK1.faa` | **PGK1** (`P00560`), phosphoglycerate kinase | 416 aa | highly conserved, single copy — the control that should always work |

STE2 is the default because a fast-evolving receptor is where a sequence-only
search starts to struggle across genera, which is the case SynVoy exists for.
PGK1 is there as the safe fallback if a live run needs to be guaranteed.

## Running it

```bash
demo/fetch_data.sh          # once — downloads into local_data/demo/ (gitignored)
demo/run_demo.sh            # default query STE2  -> results/demo_ste2
demo/run_demo.sh PGK1       # the control        -> results/demo_pgk1
demo/run_demo.sh STE2 -resume
```

Anything after the query name is passed straight through to Nextflow.

The script strips an active virtualenv from the environment before launching. A
VS Code-activated `.venv` shadows the conda env inside Nextflow tasks and makes
`parasail` disappear, which disables Smith–Waterman — the launcher fails loud on
this rather than running degraded.

## The equivalent raw command

```bash
./run_synvoy.sh -profile low_mem --mode pro \
    --query          local_data/demo/query/STE2.faa \
    --home_genome    local_data/demo/home/Saccharomyces_cerevisiae.fna \
    --home_gff       local_data/demo/home/Saccharomyces_cerevisiae.gff \
    --home_species   "Saccharomyces cerevisiae" \
    --target_genomes 'local_data/demo/genomes/*.fna' \
    --target_gffs    'local_data/demo/target_gffs/*.gff' \
    --outdir         results/demo_ste2
```

## What to show

1. `results/demo_ste2/*_anchor_grid.svg` — the payoff figure: rows = species,
   columns = genes, the GOI column filled across all four.
2. `results/demo_ste2/synvoy_report.json` — read `summary.headline` and the three
   records under `goi_dedup.records`, not `total_raw_search_hits` (which is a
   diagnostic and is routinely 0 on a good run).
3. The console `[INFO] §1m collinearity: bridged …` lines, if any fire — that is
   the gap-bridging slide happening live.

## Measured runtimes

Workstation: 12 cores, 15 GB RAM (~4 GB free), disk at 97 %.

| Run | Config | Wall time |
|---|---|---|
| cold, no cache | `low_mem` as shipped | **3m49s** (50 tasks, 2026-07) |
| cold, no cache | `low_mem`, 2026-10-05 code | 3m56s (55 tasks incl. the parameter advisor) |
| `-resume` | `low_mem` | 3m37s (10 tasks cached) |
| cold, 2026-10-11 code (distances in gene units) | `low_mem`, test harness without the launcher | 3m59s; the code before the gene units in the same harness: 2m30s to 3m03s |

The newest code is about a third slower on this demo: it searches 101 small windows where the code before searched 28 large ones (a tenth to a third of each 11 Mb genome instead of two thirds and more). With `--anchor_max_places 5` it searches 18 windows and takes 1m10s, with the same three calls.

Comfortably inside a 10-minute slot, with room to spare. If you need to fill a
longer slot, add a 4th target genome (~+1 min) rather than slowing anything down.

`low_mem` is deliberate: 4 GB per search task and serialised forks, so the demo
cannot OOM in front of an audience. On 11 Mb genomes the ceiling never binds.

## Why the overrides live in `demo.config`

`low_mem` ships one setting that is wrong for a *demo*: it sets `skip_tree = true`.
`demo/demo.config` restores the gene tree. It deliberately does **not** raise
`mmseqs_sensitivity` back from 7.0 to 9.5: measured on this demo, 9.5 costs +8 min
for identical STE2 calls.

The override is in a config file to keep the command short. Passing it as a flag
(`--skip_tree false`) works too. (Until 2026-08-21 numeric flags such as
`--mmseqs_sensitivity 9.5` aborted the run with `Cannot compare java.lang.String …`;
if you still see that, your checkout is old — `git pull`.)

## What the result looks like (checked 2026-10-11)

```
2 high-confidence + 1 medium-confidence GOI ortholog annotation(s) (+17 low-confidence) across 3 genome(s).
```

| Species | Call | Identity | Target's own gene name |
|---|---|---|---|
| *N. castellii* | HIGH | 66.0 % | STE2 |
| *K. lactis* | HIGH | 47.8 % | STE2 |
| *L. thermotolerans* | MEDIUM | 41.9 % | STE2 |

Three orthologs, one per species, each on the gene that the target's own annotation calls
STE2. Two things decide the labels. STE2 has one exon at home, so a one-exon model can be
HIGH (until 2026-10-09 the main path asked for two exons, and *N. castellii* was HIGH only
through the hull rescue). And the HIGH threshold of 50 % is lowered for a genome whose
flanking genes are themselves diverged: in *K. lactis* they match at 55 % (median), which
puts the threshold at 45 %.

(Before 2026-07-26 the same run reported "1 high + 32 medium": weak hits at 21–28 %
identity were promoted by their neighbourhood alone. An identity floor ended that, and the
figures now hide fragment hits by default.)

## Things to know before presenting

- **Report and figures agree here.** All three calls come from the search itself, so the
  anchor grid shows the same confidence as `synvoy_report.json`. In general the figures
  draw the search's own models and not those of the rescue passes; a call that exists only
  as a rescue model is in the report and in `results/<run>/rescue/locus_N/`.
- **The 17 low-confidence items are hit fragments and weak models elsewhere in the
  genome**, not candidates worth discussing. They are in the `*_with_fragments` figures
  if someone asks.
- **The search visits more, smaller places than before.** Distances are counted in gene
  spacings of the home locus (2.4 kb here), so a block is a neighbourhood of some tens of
  kb (the largest: 130 to 630 kb) instead of one to two megabases. The log prints the
  distances used, one line per genome.
- **`--auto_params true` is advisory.** On yeast its built-in rules suggest
  `--n_flanking_genes 8` and print it as a flag to re-run with; it cannot change the
  running pipeline's parameters.
