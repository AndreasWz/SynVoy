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
2. `results/demo_ste2/synvoy_report.json` — read `headline` and
   `high_confidence_goi`, not `total_raw_search_hits` (which is a diagnostic and
   is routinely 0 on a good run).
3. The console `[INFO] §1m collinearity: bridged …` lines, if any fire — that is
   the gap-bridging slide happening live.

## Measured runtimes

Workstation: 12 cores, 15 GB RAM (~4 GB free), disk at 97 %.

| Run | Config | Wall time |
|---|---|---|
| cold, no cache | `low_mem` as shipped | **3m49s** (50 tasks) |
| `-resume` | `low_mem` | 3m37s (10 tasks cached) |

Comfortably inside a 10-minute slot, with room to spare. If you need to fill a
longer slot, add a 4th target genome (~+1 min) rather than slowing anything down.

`low_mem` is deliberate: 4 GB per search task and serialised forks, so the demo
cannot OOM in front of an audience. On 11 Mb genomes the ceiling never binds.

## Why the overrides live in `demo.config`, not on the command line

`low_mem` ships two settings that are wrong for a *demo*: it drops
`mmseqs_sensitivity` to 7.0 and sets `skip_tree = true`. `demo/demo.config`
restores 9.5 and the gene tree.

They are in a **config file** because passing them as `--flags` crashes the run:

```
$ ./run_synvoy.sh ... --mmseqs_sensitivity 9.5
ERROR ~ Cannot compare java.lang.String with value '9.5' and java.lang.Integer with value '1'
```

Nextflow hands every `--param value` in as a String and `main.nf`'s validation
compares numeric params with `<`/`>`. This affects `--max_intron`,
`--cluster_distance`, `--sw_timeout_seconds` and ~20 others — see `docs/TODO.md`
§1s. Config-file values keep their numeric type and validate cleanly.

## Known rough edges to avoid on stage

- **The GOI column in the anchor grid is noisy.** Each target row carries several
  red arrows, because `fallback_hit_span` / `rescued_exon` calls at 21–28 %
  identity are drawn next to the real ortholog. The report's headline reads
  "1 high + 32 medium" when the defensible answer is "3 orthologs at 66/48/42 %".
  There is no plot-side confidence filter today. Talk to the *columns* (anchors
  recovered across all four species) and to the three identities, not to the
  medium count.
