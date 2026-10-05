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

## What the result looks like (checked 2026-10-05)

```
1 high-confidence + 2 medium-confidence GOI ortholog annotation(s) (+16 low-confidence) across 3 genome(s).
```

| Species | Call | Identity | Target's own gene name |
|---|---|---|---|
| *N. castellii* | HIGH | 66.0 % | — (found by the hull rescue) |
| *K. lactis* | MEDIUM | 47.8 % | STE2 |
| *L. thermotolerans* | MEDIUM | 41.9 % | STE2 |

Three orthologs, one per species, each matching the home STE2 gene best in the paralog
check. (Before 2026-07-26 the same run reported "1 high + 32 medium": weak hits at
21–28 % identity were promoted by their neighbourhood alone. An identity floor ended
that, and the figures now hide fragment hits by default.)

## Things to know before presenting

- **The report carries the final confidence, the figure the search-time one.** The
  *N. castellii* call is HIGH in `synvoy_report.json` because the hull rescue built a
  full-length model; the anchor grid draws the search's own models. Read the headline
  from the report. The rescued gene model itself is in
  `results/demo_ste2/rescue/locus_1/`.
- **The 16 low-confidence items are hit fragments**, not candidates worth discussing.
  They are in the `*_with_fragments` figures if someone asks.
- **`--auto_params true` is advisory.** On yeast its heuristic suggests fungal-scale
  values (`--max_intron 500 --cluster_distance 40000 …`) and prints them as flags to
  re-run with; it cannot change the running pipeline's parameters.
