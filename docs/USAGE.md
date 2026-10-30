# SynVoy — Usage & Reference Manual

Complete reference for running and configuring SynVoy. For setup, see the [README](../README.md).

> **Run SynVoy with `./run_synvoy.sh`** ([README → Run](../README.md#run)). Every runnable
> example in this manual uses it. The launcher forwards **everything** you pass —
> `--params`, `-profile`, `-resume`, `-work-dir`, all of it — straight to
> `nextflow run main.nf`, so there is no flag it cannot take.
>
> It is not a convenience wrapper you can skip: before launching it activates
> `synvoy_env`, refuses a checkout that predates the multi-locus OOM cap, verifies
> Java/Nextflow, and — the one that actually bites — catches a `.venv` shadowing the
> conda env, which silently drops `parasail`/Smith-Waterman out of the search. Calling
> `nextflow run main.nf` yourself skips all of that and also loses the laptop-safe
> `-profile auto,low_mem` default. Use the bare form only when you already know why you
> want it (a controller job inside a SLURM script is the usual reason — see §5).

---

## Table of Contents

1. [Execution Modes](#1-execution-modes)
2. [Execution Profiles](#2-execution-profiles)
3. [Algorithm Overview](#3-algorithm-overview)
4. [Full Parameter Reference](#4-full-parameter-reference)
5. [Running on HPC / SLURM](#5-running-on-hpc--slurm)
6. [Output Files](#6-output-files)
7. [Resuming & Caching](#7-resuming--caching)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Execution Modes

SynVoy has two modes: **Easy** (automated genome retrieval) and **Pro** (local files).

### Easy Mode

Provide a UniProt/NCBI protein accession, a local FASTA (`--query`), or an inline sequence (`--query_seq`). SynVoy resolves the query, fetches the reference genome and related target assemblies from NCBI, and runs the full analysis.

```bash
./run_synvoy.sh \
  --mode easy \
  --query_id Q16553 \
  --max_genomes 5 \
  --outdir results/my_run \
  -profile standard
```

**Required (one of the query options):**

| Flag | Description |
|---|---|
| `--mode easy` | Select Easy Mode |
| `--query_id` | UniProt accession (e.g. `Q16553`) or NCBI protein ID |
| `--query` | Path to local FASTA (works in Easy Mode too) |
| `--query_seq` | Inline protein sequence or FASTA text (requires `--home_species`) |

**Optional (Easy Mode only):**

| Flag | Default | Description |
|---|---|---|
| `--home_species` | auto-detected | Override the reference species (e.g. `"Homo sapiens"`) |
| `--target_species` | auto (taxonomic) | Comma-separated species list instead of auto-selection (e.g. `"Gallus gallus,Mus musculus"`) |
| `--max_genomes` | `0` (auto) | Number of related genomes to fetch. `0` = auto-detect (~3 per taxonomic level) |
| `--assembly_ranking` | `hybrid` | How to pick among multiple assemblies for one species: `hybrid`, `counts`, `nstats` |
| `--bad_quality_policy` | `auto` | Low-quality assemblies: `auto` drops one the pipeline picked and keeps, flagged, one you asked for (`--target_species`, the home assembly); `drop`, `keep`, `ask` apply one rule to all |
| `--bad_quality_timeout` | `300` | Seconds to wait for user input when policy is `ask` |
| `--bad_max_contigs` | `500000` | Assemblies with more contigs are flagged as low quality |
| `--bad_max_scaffolds` | `500000` | Assemblies with more scaffolds are flagged as low quality |
| `--bad_min_n50` | `5000` | Assemblies with N50 below this are flagged as low quality |

> **How the quality check decides.** Chromosome-level and complete assemblies always pass.
> For the others the N50 decides in practice: NCBI's assembly records carry N50 values but
> no contig or scaffold counts, so the two count limits apply only when a count is known.

> **Note:** When using `--query_seq`, you must also provide `--home_species`.

### Pro Mode

Supply your own query FASTA, reference genome, and target genome files. Works offline.

```bash
./run_synvoy.sh \
  --mode pro \
  --query queries/melittin.faa \
  --home_genome /data/apis_mellifera.fna \
  --home_gff /data/apis_mellifera.gff \
  --home_species "Apis mellifera" \
  --target_genomes "/data/targets/*.fna" \
  --outdir results/melittin \
  -profile standard
```

**Required:**

| Flag | Description |
|---|---|
| `--mode pro` | Select Pro Mode |
| `--query` | Path to query protein FASTA (DNA sequences are auto-translated to the best ORF) |
| `--home_genome` | Path to reference genome FASTA (`.fna` or `.fna.gz`) |
| `--target_genomes` | A folder of genome FASTAs, a glob pattern, or a comma-separated list of target genome FASTAs |

**Optional:**

| Flag | Description |
|---|---|
| `--home_species` | Reference species (e.g. `"Apis mellifera"`), used for taxonomy ordering. Inferred from the `--home_genome` filename when it *is* a species name (e.g. `apis_mellifera.fna`); SynVoy **errors out with a clear message if it can't infer one** (a name like `reference.fna`). Pass it explicitly, or set `--allow_unknown_species true` to proceed with `kingdom=Unknown` / phylo distance 999 (mis-tunes per-target stringency). |
| `--home_gff` | GFF annotation for the home genome. Highly recommended — provides much better flanking-gene extraction than Prodigal fallback. |
| `--target_gffs` | Annotations for the **target** genomes (folder, glob, or comma-list). Lets SynVoy read each target's existing gene names/products instead of judging hits by sequence alone. Easy Mode downloads these automatically. **Added 2026-07-21** — `git pull` if your checkout is older (see the version note below). |
| `--allow_unmatched_target_gffs` | `true` to warn instead of failing when a supplied GFF matches no genome. Default `false`. |

> **Tip:** `--target_genomes` accepts a **folder** of genomes (`genomes/` — every `.fna`/`.fa`/`.fasta` inside is used; GFF/TSV/other files are ignored), a glob (`"genomes/*.fna"` — match your extension: `*.fa`, `*.fasta`), a comma-separated list (`"a.fna,b.fna"`), or Nextflow list syntax. SynVoy errors clearly if no genomes are found.

> **File names matter.** The genome file name (minus extension) becomes the genome's
> name inside sequence IDs, GFF attributes and output file names. SynVoy stops at launch
> if a query, home or target file name contains whitespace or any of
> `| ; = , ' " ( ) [ ] { } < > * ? ! # % & $`, or if two targets differ only by
> extension (`Apis.fa` + `Apis.fna`). Use letters, digits, `.`, `_` and `-` —
> e.g. `Apis_mellifera.fna`.

> **The query is one protein.** `--query` may be protein or DNA (DNA is translated to its
> longest ORF). Case, whitespace, digits, alignment gaps and a trailing `*` are cleaned
> away; anything else that is not an amino-acid code is an error. If the file holds
> several sequences, **only the first is used** (the task log says so).

> **Tip:** Because `--target_genomes` is FASTA-only, target annotations are passed **separately** via `--target_gffs` — keep genomes and GFFs in two folders:
>
> ```bash
> ./run_synvoy.sh --mode pro \
>   --query queries/decorin.faa \
>   --home_genome /data/human.fna --home_gff /data/human.gff \
>   --home_species "Homo sapiens" \
>   --target_genomes /data/targets/genomes \
>   --target_gffs /data/targets/gffs
> ```
>
> Each GFF is matched to its genome by **filename stem** (`cow.fna` ↔ `cow.gff`) or by **assembly accession** (`GCF_002263795.1.fna` ↔ `GCF_002263795.1_ARS-UCD1.2_genomic.gff`). Annotating only *some* targets is fine. A GFF that matches no genome is a **hard error** by default — pass `--allow_unmatched_target_gffs true` to downgrade it to a warning. `.gff`, `.gff3`, `.gtf` and their `.gz` forms are accepted.
>
> **Version note.** `--target_gffs` landed on `main` on **2026-07-21**. On an older
> checkout Nextflow accepts the flag and silently ignores it (unknown `--params` are not
> an error), so you get a run with no target annotations rather than a message. If you
> instead see `Unknown option: --target_gffs`, that is Nextflow's *CLI* parser, not
> SynVoy: the flag reached it before the `run` sub-command. Keep pipeline flags after
> the script — `./run_synvoy.sh` already does this for everything you pass it. Check
> what you are on with `./run_synvoy.sh --check-only` (it prints the version banner and
> warns when you are behind the remote), then `git pull`.
>
> **What it buys you:** matched target models gain `TargetGene` / `TargetProduct` / `TargetID` attributes in the output GFF and the `target_gene` / `target_product` report columns — i.e. you can see *which annotated gene* a call landed on. This is read-out evidence, not a scoring input: annotations do not move synteny blocks, region scores, or confidence. With `--goi_family_tokens` (your gene's family names), every GOI row also records whether its target gene name carries one of them (`GoiFamilyConsistent`); that is annotation only since 2026-10-15, when the strict gate that demoted calls by name (`--strict_goi_family`) was removed.

---

## 2. Execution Profiles

Append a profile with `-profile <name>` to control how resources are allocated. **Profiles are combinable** as long as you don't pick two from the same axis — pick **one execution profile** (`auto` / `standard` / `docker` / `singularity` / `slurm` / …), optionally **one memory tier** (`laptop_safe` / `low_mem`), and optionally **one biology preset** (`preset_single_copy` / `preset_short_peptide` / …). Example:

```bash
-profile standard,low_mem,preset_single_copy
```

> **Not sure how to tune your run? Use `-profile auto`.** It runs locally
> through the conda env like `standard`, and guarantees the full auto-tuning
> stack is on: the query preset is auto-picked from the hit distribution (§1f),
> the home-locus count is capped (§1d), strong-synteny rescue and the
> reciprocal-best paralog check run (§1e/§1j). It's the recommended starting
> point when you don't already know whether your query is short / divergent /
> single-copy / paralog-rich. On a small machine, add a memory tier:
> `-profile auto,low_mem`. (Docker / Singularity / HPC users already get this
> behaviour for free — the stack is on by default — so keep using your
> execution profile.)

Don't combine two execution backends (e.g. `docker,singularity`) or two memory tiers — later values silently override earlier ones, which is usually not what you want.

### Execution backends

| Profile | Executor | Environment | Description |
|---|---|---|---|
| `auto` | local | Conda | **Recommended zero-config default.** Like `standard`, plus it pins the full auto-tuning stack on (auto preset selection, locus cap, strong-synteny rescue, paralog check). Best when you don't want to hand-pick a preset. |
| `standard` | local | Conda | Default. 2 CPUs / at least 10 GB RAM per iterative-search task (more for large target genomes and on a retry, see `iterative_search_memory`), single-fork. Good baseline for workstations. |
| `conda` | local | Conda | Same as `standard` but explicitly disables Docker/Singularity. |
| `docker` | local | Docker | Runs all processes inside the `synvoy-local:latest` container. Build it first with `docker build -t synvoy-local:latest .` |
| `docker_max` | local | Docker | Auto-detects all host CPUs and RAM. Allocates nearly everything to the heaviest tasks (MMseqs2, ITERATIVE_SEARCH). Single-fork to avoid OOM. Ideal for dedicated machines. |
| `singularity` | local | Singularity | Like `docker` but uses Singularity with auto-mounts. |
| `slurm` | SLURM | (none) | Submits tasks to a SLURM scheduler. Edit `nextflow.config` to set your partition and account. |
| `hpc_singularity` | SLURM | Singularity | SLURM + Singularity containers. Caches images in `~/.singularity/cache`. |
| `hpc_conda` | SLURM | Conda+Mamba | SLURM + Conda (uses Mamba for faster env creation). |
| `lrz_ai` | SLURM | Conda | LRZ AI Systems (GPU) preset — sets `gpu_partition` / `gpu_cluster_options` and caps queue size. Conda-based, so it needs a system where conda works; the AI Systems themselves are container-only (see `lrz_ai_container`). |
| `lrz_ai_container` | SLURM | Enroot/pyxis | LRZ AI Systems container path — runs the pipeline from a `.sqsh` image instead of conda. |
| `test` | local | Conda | Loads `conf/test.config` with small test data and relaxed thresholds for CI. |
| `test_melettin` | local | Conda | `conf/test_melettin.config` — the melittin ground-truth regression fixture. |
| `test_tetramorium` | local | Conda | `conf/test_tetramorium.config` — the *Tetramorium* myrmicitoxin regression fixture. |

### Memory tiers (combine with an execution backend)

The defaults assume a workstation with at least ~12 GB free RAM. On laptops, lower the per-process ceilings — the single most important knob is `mmseqs_split_memory_limit`, which bounds the per-split target-DB size that MMseqs2 keeps resident during search. Pick the tier that matches your free RAM, not your total RAM (subtract whatever the OS + browser are using).

| Profile | Target machine | `iterative_search_memory` | `mmseqs_split_memory_limit` | `mmseqs_sensitivity` |
|---|---|---|---|---|
| (none — uses defaults) | ≥16 GB free RAM | `auto` (at least 10 GB; more for large genomes, x the attempt) | 8G | 9.5 |
| `laptop_safe` | ~16 GB RAM laptop | 8 GB | 4G | 8.0 |
| `low_mem` | ~8 GB RAM laptop | 4 GB | 2G | 7.0 |

Trade-off: tighter tiers are slower (more MMseqs splits, less parallelism) and slightly less sensitive on highly divergent queries. They're the right pick when you're seeing `cannot fit database into ... / not enough memory to keep dbreader/write in memory` errors from MMseqs2 — that means your `mmseqs_split_memory_limit` is too high for the genome you're searching.

### Biology presets (combine with an execution backend + memory tier)

| Profile | Use for |
|---|---|
| `preset_short_peptide` | Queries shorter than 100 aa: coverage bands of the model status 0.5 / 0.3. |
| `preset_paralog_discrimination` | The query aligns well at two or more places of the home genome (a family on several chromosomes): fewer candidate regions per genome in the figures. |
| `preset_tandem_family` | The home locus itself holds more than one copy of the gene: distant family members are not used as anchors either. |
| `preset_single_copy` | Everything else. Equal to the defaults. |

Each preset holds only the values that have a reason for that kind of query; the table of values and the measurements behind them are in [PARAMETERS.md](PARAMETERS.md#quick-start-pick-a-preset). No preset moves the label of a GOI call: one judge sets it for every query ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)).

> **Auto-apply (since 2026-05-29, default ON):** if you don't pick a preset, SynVoy picks one for you after `LOCATE_GENE` from the hit distribution and applies it at runtime to `EXTRACT_FLANKING` / `ITERATIVE_SEARCH` / `CLUSTER_REGIONS` via a value channel — no restart needed. The chosen preset is logged as `[INFO] §1f preset applied: …` and dumped to `intermediate/locate_gene/effective_params.json`. Pin one with `--preset_override preset_X`, opt out with `--auto_apply_preset false`. **A value you set yourself wins over the preset**: any preset-covered parameter that differs from the shipped default — set on the command line, in a params file, a `-c` config or a `-profile` — is kept, and `effective_params.json` lists it under `user_set_kept` (before 2026-10-05 the preset silently overwrote it). The one case this cannot see is passing the shipped default explicitly; use `--auto_apply_preset false` to pin that. The other caveat: `NORMALIZE_QUERY.min_query_length` is resolved before LOCATE_GENE runs, so short queries (< 30 aa) still need `--min_query_length 20` or `-profile preset_short_peptide` set at launch. See [docs/PARAMETERS.md §19](PARAMETERS.md#19-auto-apply-preset-self-consistency--rescue) for the full toggle list.

---

## 3. Algorithm Overview

The pipeline proceeds through five phases:

### Phase 1 — Gene Localization

1. **Normalize Query:** If the input is DNA, the best ORF is translated to protein.
2. **Locate GOI:** The query protein is aligned against the home genome using tblastn and MMseqs2 to establish coordinates.
3. **Annotate GOI Exons:** If a GFF is available, the GOI is matched to an annotated gene and individual CDS/exons are extracted. Otherwise, exon boundaries are inferred from alignment hits.
4. **Split Loci:** The annotated home gene that matches the query is the primary locus. Other places where the query aligns at least half as well (relatives elsewhere in the genome) become further loci, each processed independently. Places with at most `--locus_gap_genes` annotated genes between them are one locus.
5. **Extract Flanking Genes:** The *n* genes upstream and downstream of each locus are identified from GFF or Augustus/Prodigal prediction. Neighbours that are relatives of the GOI (pairwise alignment E-value at or below `--goi_relative_evalue`) are skipped: their hits would fall on the GOI's own relatives. A `--max_flanking_distance` cap can prevent walking into distant gene deserts.
6. **Expand GOI-Similar Neighbors** *(optional, off by default)*: When `--expand_goi_similar` is enabled, genes near the GOI that resemble it (e.g. tandem duplicates like MRJPs near Yellow-e3) are emitted as additional GOI queries with a `GOI_NEIGHBOR_` prefix. These are searched in all target genomes alongside the original GOI, and included in the phylogenetic tree — enabling resolution of paralogs vs. orthologs.
7. **Borrow Annotations:** When the home genome lacks a GFF, annotations can be borrowed from annotated target genomes via reciprocal best hits.

### Phase 2 — Phylogenetic Ordering & Iterative Search

7. **Phylo Sort:** Target genomes are ordered by evolutionary distance to the reference.
8. **Genome Quality Assessment:** Target assemblies are evaluated for contiguity (N50, scaffold count).
9. **Iterative Search:** For each target genome (nearest-first), flanking genes are mapped with MMseqs2. Hits are clustered into candidate syntenic blocks; two clusters separated by a rearrangement gap are bridged into one block when the flanking genes continue the home gene order. Within each block, localized MMseqs2, tblastn and Smith-Waterman searches (BLOSUM62 11/1 with real E-values) locate the GOI, and miniprot builds the gene model. The model's protein is miniprot's own translation, and its two ends are extended to the start and stop codon where the reading frame allows. Each call gets a confidence: `HIGH`, `MEDIUM`, `AMBIGUOUS` (in a conserved neighbourhood, but nothing shows the sequence is the ortholog) or `LOW`. Discovered genes are added to the search database for the next, more distant wave.

### Phase 3 — Region Clustering

10. **Cluster Regions:** Candidate blocks across all targets are filtered by synteny score and ranked.

### Phase 4 — Phylogenetics & Visualization

11. **Compute Tree:** All discovered GOI and GOI-similar sequences across all genomes are aligned (MAFFT) and a phylogenetic tree is inferred (IQ-TREE with automatic model selection and ultrafast bootstrap). Multiple hits per genome are preserved, so the tree can resolve paralogs from orthologs.
12. **Plot Synteny:** An interactive HTML plot shows the syntenic context of each hit, colored by homology, with the phylogenetic tree alongside.

### Phase 5 — Adjudication & Reporting

13. **Rescue passes:** where a neighbourhood is clearly the right one but holds no GOI model, a relaxed miniprot pass searches the block (strong-synteny rescue) or the whole flanking neighbourhood including the gaps between blocks (hull rescue).
14. **Adjudication:** calls are de-duplicated across home loci, checked against the home genome's paralogs (a call that matches a paralog better than the GOI becomes `paralog_not_goi`), demoted when a high identity rests on a short slice of the query, and tested for phylogenetic placement (advisory).
15. **Generate Report:** `synvoy_report.json` carries the adjudicated calls (`goi_dedup.records`), the headline counts, the self-consistency flags and the candidates that were rejected by a synteny gate.

---

## 4. Full Parameter Reference

All parameters can be set on the command line (`--param value`) or in a custom config file (`-c my.config`).

### Distances in gene units

"How far apart" means something else in every genome: ten neighbouring genes span about 20 kb in yeast, 150 kb in the honeybee and 1 Mb in human. The distances of the search (largest gap inside a block, extent of the neighbourhood around a block, largest gap a block may bridge, and two more) are therefore given in **gene pitches**, the mean step between neighbouring anchors of the home locus, and converted to base pairs for each target genome with the **expansion** measured there (how much larger or smaller the same neighbourhood is in that genome). You do not set them. The expansion is measured only where the gene order is conserved (at least five anchors in a row in home order); in a genome where it is not, the ratio of the genome sizes is used. The values used are printed per genome:

```
[Naumovozyma_castellii.fna] distances: home pitch 2.4 kb x expansion 1.30 (measured over 15 anchors in a row) -> block gap 22 kb, window margin 22 kb (DNA searched: 22 kb per side), GOI margin 6 kb, bridge up to 63 kb.
```

Each distance still has a parameter in base pairs (`auto` by default); a number there is used as it is for every genome. `-profile <execution>,legacy_distances` sets all of them to the fixed values used before 2026-10-09: for comparisons, and for a genome where a derived window turns out too large (see Troubleshooting). Definitions, defaults and the reasons for them: [PARAMETERS.md § 5](PARAMETERS.md#distances-in-gene-units).

### Synteny & Search

| Parameter | Default | Description |
|---|---|---|
| `--n_flanking_genes` | `10` | Number of non-GOI-similar flanking genes to extract on each side of the GOI. Four counts follow it (`auto`): the bridge limit, the strong-synteny and hull-rescue flanking counts and the gene gap between home loci |
| `--min_flanking_size` | `500` | Minimum size (bp) for a flanking gene to be included |
| `--goi_relative_evalue` | `1e-5` | A neighbour whose protein aligns to the GOI with an E-value at or below this (pairwise Smith–Waterman) is a relative of the GOI and is not used as an anchor. `0` = off. Replaces `--max_flanking_goi_similarity`, which removed about half of all genes of the GOI's size whatever their sequence. |
| `--max_flanking_distance` | `0` | Max distance (bp) from GOI center to walk for flanking genes. `0` = unlimited. Useful when the GOI neighbours a large tandem array. |
| `--expand_goi_similar` | `false` | Emit GOI-similar flanking genes as additional GOI queries (`GOI_NEIGHBOR_` prefix). Enables paralog discovery and phylogenetic resolution across genomes, but can flood paralog-rich searches; opt in deliberately. |
| `--expand_goi_similar_distance` | `300000` | Max distance (bp) from GOI to search for GOI-similar neighbor genes |
| `--exon_level_search` | `true` | Search at exon level for better divergent-gene detection |
| `--cluster_distance` | `auto` | Largest gap (bp) between two flanking-gene hits of one block. `auto` = `--block_gap_genes` gene pitches of the home locus × the expansion of the target genome (see "Distances in gene units" above); a number fixes it for every genome (150000 before 2026-10-09). |
| `--block_gap_genes` | `7` | The same gap in gene pitches (mean step between neighbouring anchors at home). |
| `--min_synteny_score` | `0.6` | Region score from which a candidate region is labelled HIGH (MEDIUM from half of it). A label only: it does not decide what is searched |
| `--min_hit_identity` | `10` | Minimum alignment identity (%) for an individual hit |
| `--min_hit_length` | `10` | Minimum alignment length for an individual hit |
| `--min_query_length` | `20` | Reject queries shorter than this (aa) during `NORMALIZE_QUERY`; `0` = no check. A query of 20–39 aa runs with a warning: its alignment alone is not significant in another species, so a call rests on the neighbourhood. |
| `--search_evalue` | `0.01` | E-value threshold for tblastn/MMseqs2 searches |
| `--max_intron` | `20000` | Maximum intron length (bp) for miniprot gene models |
| `--region_padding` | `auto` | Extent of the neighbourhood (bp) on each side of a block: blocks whose extents overlap are merged, anchors inside count as support. `auto` = `--window_margin_genes` pitches × expansion; the DNA cut out and searched is then at most this and not more than GOI margin + the room for one gene model. A number restores the former behaviour: an adaptive margin between `--padding_min` and `--padding_max`, searched in full. |
| `--window_margin_genes` | `7` | The same extent in gene pitches. |
| `--goi_margin_genes` / `--goi_margin` | `2` / `auto` | How far outside the block of flanking genes a GOI hit may lie and still be modelled (pitches / bp; 50000 before; `local_family`: the home array's reach + 1, at most `window_margin_genes`). |
| `--padding_min` | `50000` | Minimum padding (bp); only read when `--region_padding` is a number |
| `--padding_max` | `200000` | Maximum padding (bp) |
| `--max_blocks_per_genome` | `80` | Safety cap on candidate blocks per target genome |
| `--min_block_genes` | `2` | Minimum flanking-gene hits in a block to keep it |
| `--disable_synteny_collinearity` | `false` | Disable collinearity-aware block scoring + gap-bridging (revert to legacy gene-count clustering). When on (default), blocks are ranked by their longest run of home-ordered flanking genes and a GOI sitting in a rearrangement gap between two collinear clusters is kept inside one searched block. |
| `--synteny_bridge_max_gap` | `auto` | Largest gap (bp) between two same-chromosome flanking clusters that is bridged into one block when they continue the home gene order. `auto` = `--synteny_bridge_max_gap_genes` pitches × expansion (6000000 before). `0` disables bridging. |
| `--synteny_bridge_max_gap_genes` | `auto` | The same gap in gene pitches. `auto` = 2 × `--n_flanking_genes` (20 at the default): one home neighbourhood per side. `0` disables bridging. |
| `--anchor_max_places` | `auto` | Each flanking gene seeds blocks only at its N best places in a target genome (0 = every place). Cuts the blocks that two anchors form by chance: with `5`, the same orthologs from 18 windows instead of 101 in yeast and from 9 instead of 238 in four plant genomes. `auto` = every place; with any `--modules` other than `legacy`, 5 in a genome where the block cap binds (the genome axis). A number is used for every genome; see PARAMETERS.md |
| `--disable_anchor_resolution` | `false` | Where hits of different flanking genes cover the same target DNA only the best one counts as an anchor (a target gene is the ortholog of one home gene at most). `true` counts them all, as before. |
| `--disable_relative_query` | `false` | A window is searched with the home protein and, in addition, the gene model of the genome searched most recently; identity, coverage and confidence always refer to the home protein. `true` = home protein only. |
| `--decoy_mode` | `false` | Calibration: also search the reversed home protein and two shuffles of it in every window (and of the relative's query, where a window has one); their calls go to `<outdir>/decoys/` (never seeded or counted), the report gives the decoy calls of the searched cases. About four times the window search. |
| `--goi_relative_statistic` | `sw` | Statistic of the relative test of the flank walk: `sw` or `blastp_cbs` (composition-based; `low_complexity` module). |
| `--report_low_complexity` | `false` | `LowComplexity=` on GOI models (`low_complexity` module). |
| `--core_completeness` | `false` | HIGH's completeness also on the protein without signal peptide and GPI tail (`CoreCoverage`; `signal_anchor` module). |
| `--goi_max_models_per_locus` | `10` | Gene models kept per candidate locus: every miniprot model that does not lie on a better one (the copies of a tandem array). `1` = best model only, as before (`local_family`: home copies + 1, at most 60). |
| `--goi_miniprot_timeout` | `3600` | Guard of a miniprot run that models the GOI (search windows, the home gene model, the hull rescue): a run that does not finish stops the task (until 2026-10-28: 120 s, in the hull rescue 600 s, then no model). |
| `--flanking_miniprot_timeout` | `3600` | Guard of the other miniprot runs (flanking genes in a window, the transfer of target genes to a home without a GFF): a run that does not finish stops the task (until 2026-10-28: 120 s, then a hit-chain model, or the target gene was not transferred). |
| `--goi_candidate_regions_per_query` | `2` | GOI loci kept per query in a window (the `local_family` module raises it to the home copies + 1, at most 40). |
| `--report_copies_in_slot` | `false` | `CopiesInSlot` / `SlotRank` on GOI models and home : target copies per genome in the report (on with `local_family`). |
| `--window_tblastn_word_size` / `--window_tblastn_comp_stats` | `3` / `2` | tblastn options of the window search (finding only; the `short_query` module sets 2 / 0). |
| `--window_mmseqs_min_length` | `30` | MMseqs2 ORF floor (codons) of the window search (`short_query`: 10). |
| `--miniprot_short_orf_max_len` | `60` | Queries below this length get miniprot `-L 10` in the windows (`short_query`: 100). |
| `--synteny_bridge_max_rank_gap` | `5` | Max home-rank jump still treated as a collinear continuation when bridging. |
| `--synteny_bridge_min_anchors` | `2` | Min anchored flanking genes a block needs before it may bridge a gap (3 until 2026-10-14). |

### Smith-Waterman Local Search

| Parameter | Default | Description |
|---|---|---|
| `--enable_smith_waterman` | `true` | Use rigorous Smith-Waterman alignment (parasail) for GOI search |
| `--sw_method` | `auto` | Implementation: `auto`, `parasail`, or `ssearch36` |
| `--sw_min_score` | `20` | Minimum SW alignment score to report a hit |
| `--sw_min_identity` | `10.0` | Minimum identity (%) for SW hits |
| `--sw_timeout_seconds` | `3600` | Guard of the Smith-Waterman search of one window: a search that does not finish stops the task (until 2026-10-28: 300 s, then the window lost its Smith-Waterman hits) |
| `--sw_max_alignments_per_piece` | `1` | SW alignments per piece of a frame; above 1, further exons or copies in the same frame while E ≤ 1 (`local_family` / `exon_structure`: 3) |

### Window search thresholds

One search per window (the names date from a design with several relaxed rounds). A window hit of MMseqs2 or tblastn is kept at E ≤ 0.1, identity ≥ 15 %, length ≥ 15 aa; details in [PARAMETERS.md](PARAMETERS.md).

| Parameter | Default | Description |
|---|---|---|
| `--aug_relaxed_parse_evalue_mult` | `10` | Window hits are kept at E ≤ `search_evalue` × this (0.1) |
| `--aug_relaxed_evalue_mult` | `1000` | With the cap: the E-value handed to the search programs (10); does not change what is kept unless set below the line above |
| `--aug_relaxed_evalue_cap` | `10.0` | See above |
| `--aug_relaxed_identity_min` | `15.0` | Identity floor (%) of a window hit |
| `--aug_relaxed_identity_factor` | `0.6` | Floor becomes `min_hit_identity` × this when that is larger than the line above |
| `--aug_relaxed_length_min` | `15` | Length floor (aa) of a window hit (`short_query`: 10) |
| `--aug_relaxed_length_div` | `2` | Floor becomes `min_hit_length` / this when that is larger than the line above |
| `--fallback_short_query_len` | `150` | Queries shorter than this (aa) gate their GOI fallback on aligned length + bitscore instead of query coverage — a short mature peptide in a longer precursor (e.g. melittin) is structurally low-coverage even on a perfect hit. `0` = legacy coverage gate for all queries. |
| `--fallback_short_min_aln_aa` | `15` | Minimum aligned length (aa) for a short-query GOI fallback. Rejects micro-window overcalls while admitting short mature peptides (`short_query`: 10). |
| `--fallback_short_min_bits` | `30` | Minimum bitscore for a short-query GOI fallback. Lets a high-bitscore divergent hit pass where low coverage would have vetoed it. |

### MMseqs2

| Parameter | Default | Description |
|---|---|---|
| `--mmseqs_sensitivity` | `9.5` | MMseqs2 sensitivity (1–10+). Higher = slower but more sensitive |
| `--mmseqs_split_memory_limit` | `8G` | MMseqs2 memory limit for database splitting. Override to `3G` or `1G` on memory-constrained machines. |
| `--mmseqs_verbosity` | `1` | MMseqs2 log verbosity (0 = silent) |
| `--deterministic_goi_search` | `true` | Force the per-region augmented GOI MMseqs2 search to `--threads 1` so marginal/divergent GOI calls are reproducible run-to-run (regions are tiny, so the speed cost is small). Set `false` for the legacy multithreaded (non-deterministic) path. |

### Annotation & Gene Prediction

| Parameter | Default | Description |
|---|---|---|
| `--gff_search_window` | `100000` | Distance (bp) from the query's hits within which annotated home genes are compared with the query to find the home gene |
| `--gap_search_window` | `50000` | DNA (bp) added around a candidate's hits before miniprot builds its model; never less than the GOI's own span at home × expansion |
| `--gap_min_identity` | `15.0` | Lower bound of the identity for last-resort `raw_hit` rows; the bound used is `max(value, 90)` |
| `--gap_max_hits` | `5` | Candidate hits looked at for the last-resort `raw_hit` row of a locus (one row is written) |
| `--min_exon_alnlen` | `30` | Shortest alignment (aa) for a last-resort `raw_hit` row; never less than 60 % of the query |
| `--pred_flank_window` | `50000` | Prodigal prediction window around locus |
| `--pred_keep_pct` | `0.10` | Fraction of Prodigal predictions to keep |
| `--prodigal_full_genome_fallback` | `false` | Run Prodigal on entire genome if windowed prediction fails |
| `--gene_predictor` | `auto` | Gene predictor for a home genome without a GFF (and the optional ML steps): `auto`, `augustus`, or `prodigal`. With any `--modules` other than `legacy`, `auto` stops the run when a eukaryotic home would fall back to Prodigal; a failed Augustus run stops the prediction (no fallback to Prodigal since 2026-10-29) |
| `--gene_prediction_timeout` | `3600` | Guard of one Augustus or Prodigal run: a run that does not finish stops the task (until 2026-10-28: 600 s / 300 s) |
| `--augustus_species` | `auto` | Augustus species model. `auto`: `fly` under `--modules legacy`; otherwise, for a home without a GFF, the installed model nearest the home species (e.g. `candida_albicans`, `honeybee1`, `human`), recorded in `effective_params.json`. A model name is used as given; one Augustus cannot load stops the gene prediction with the nearest installed models |

### Gene Model Classification

The confidence label (HIGH / MEDIUM / AMBIGUOUS / LOW, with `Basis`) of a GOI call has no parameters: one judge sets it from the E-value in the expected place, specificity against the home proteome and completeness ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)); the identity thresholds of the earlier label were removed on 2026-10-15 ([PARAMETERS.md § 20](PARAMETERS.md#20-removed-parameters)). The parameters below set the model status labels (complete/partial/fragment) and the report's coverage rule.

| Parameter | Default | Description |
|---|---|---|
| `--goi_family_tokens` | _(empty)_ | Comma-separated family names. Every GOI row then records whether the target's own annotation of the gene carries one (`GoiFamilyConsistent`, `GoiFamilyReason`). Annotation only: no label changes. |
| `--identity_decoupled_min_identity` | `50.0` | With the next row: a call at or above this identity … |
| `--identity_decoupled_max_qcov` | `0.35` | … whose *recorded* query coverage is below this is relabelled `identity_coverage_decoupled` and dropped from the headline counts (a 10 aa window at 100 % is not an ortholog). |
| `--disable_coverage_demotion` | `false` | `true` keeps that rule advisory only (a flag in `self_consistency`, no demotion). |
| `--strong_synteny_no_goi_min_flanking` | `auto` | The report flags a genome, home locus and scaffold with at least this many HIGH flanking models and no HIGH GOI call (`strong_synteny_no_goi` in `self_consistency`; advisory). `auto` = 2 × `--n_flanking_genes` (20 at the default). |
| `--classify_fragment_max_qcov` | `0.4` | Query coverage below this marks a gene model as `fragment` in the ModelStatus field |
| `--classify_complete_min_qcov` | `0.7` | Query coverage above this (with multi-exon evidence) marks a model as `complete` |

**ModelStatus** is a GFF attribute independent of confidence that labels the completeness of a gene model:
- `complete` — query coverage >= 0.7 and multi-exon (or tandem copy)
- `partial` — between fragment and complete thresholds
- `fragment` — query coverage < 0.4, or evidence from rescued_exon / raw_hit only

### Synteny Scoring Weights

| Parameter | Default | Description |
|---|---|---|
| `--synteny_weight_base` | `0.4` | Weight for base synteny score |
| `--synteny_weight_consistency` | `0.3` | Weight for gene-order consistency |
| `--synteny_weight_strand` | `0.3` | Weight for strand conservation |
| `--synteny_goi_overlap_bonus` | `0.15` | Bonus for blocks that overlap a GOI annotation |
| `--max_regions` | `0` | Max regions to emit per locus. `0` = adaptive (all above threshold, capped at `adaptive_max_regions`) |
| `--adaptive_score_floor_frac` | `0.30` | Adaptive mode: fraction of best_score used as score floor. Raise to tighten (e.g. `0.45` for paralog discrimination) |
| `--adaptive_score_floor_abs` | `0.03` | Adaptive mode: absolute score floor. Permissive default preserves weak-synteny discovery for toxins/venoms/micro-exon peptides |
| `--adaptive_max_regions` | `6` | Adaptive mode: hard cap on emitted regions |
| `--adaptive_unique_gene_floor` | `3` | Adaptive mode: clusters with >= this many unique flanking hits are kept even below the score floor |
| `--disable_distant_synteny_rescue` | `false` | Skip the distant-macrosynteny rescue in region ranking, which promotes a strong GOI hit that is isolated from its local flanking genes but sits on a chromosome that still carries an anchor in a rearranged genome. |
| `--legacy_strand_score` | `false` | Reproduce pre-2026-07-26 region scores (strand term computed without the home strand). For comparing against old runs only. |

The region score is **coverage-dominated**: `synteny_score = Q × C` with `Q = 0.4·C + 0.3·K + 0.3·S`, where `C` is the fraction of flanking genes found, `K` the gene-order consistency and `S` the strand agreement with the home genome. Order and strand modulate the score by at most `0.3·C`; a uniformly inverted neighbourhood scores as conserved. The p-value is a permutation test on `synteny_score` only (the GOI bonus is not part of it).

### Visualization

| Parameter | Default | Description |
|---|---|---|
| `--plot_width` | `1500` | Width of the output HTML plot (px) |
| `--gap_threshold` | `50000` | Gaps larger than this (bp) are visually compressed |
| `--gap_visual_size` | `20000` | Size (bp) used to represent compressed gaps |
| `--flank_fallback_bp` | `1000000` | Maximum window (bp) rendered around distal targets |
| `--scale_bar_len` | `auto` | Scale bar of the ribbon plot: `auto` = a round length about 120 px wide at the scale of the drawing, or a length in bp |
| `--hide_goi_absent_tracks` | `false` | Leave genomes without any GOI call out of the figures. By default every searched genome is drawn |
| `--plot_block_selection` | `false` | A target row shows the search's blocks chosen by the label's scoring (the blocks of HIGH, MEDIUM and AMBIGUOUS GOI calls, then the significant blocks, up to 6) instead of the candidate regions; see [PARAMETERS.md](PARAMETERS.md#plot_block_selection) |
| `--max_goi_per_genome` | `10` | How many GOI models of a genome steer the choice of the neighbourhood its row shows (every GOI model inside the chosen view is drawn regardless). Settings key: the `local_family` module (A1) sets max(10, home copies + 1), at most 60. Until 2026-10-26 a `--plot_extra_args` option, which still works and wins. |
| `--plot_extra_args` | _(empty)_ | Extra options passed verbatim to `bin/plot_synteny.py`, for the figure controls below. Quote the whole value: `--plot_extra_args '--print_grid_numbers --no_fragment_variant'`. |

Figure controls available through `--plot_extra_args`, or after the run with `scripts/replot.sh <outdir>/plot_inputs_synteny_block_locus_1 <new_folder> [options]`, which re-draws the figures of a finished run without re-running the search:

| Option | Default | Effect |
|---|---|---|
| `--print_grid_numbers` | off | Print version of the anchor grid (`*_anchor_grid_print.svg`): write the identity into every flanking arrow. The full grid always has it. |
| `--print_width_mm` | `183` | Width of the print grid in mm (183 = two journal columns, 89 = one). Columns with the fewest orthologs are left out when the grid does not fit, and reported. |
| `--print_height_mm` | `170` | Largest height of the print grid; rows are packed tighter before it is exceeded. |
| `--no_print_grid` | off | Do not write the print grid. |
| `--grid_variants` | _(none)_ | Also write try-out versions of the print grid, comma-separated or `all`: `compact` (the flanking genes of each side as one grey stack with how many were placed; every GOI copy's exon model in the GOI column), `compact_panel` (the same with the exon panel on the right), `clades_muted` (the print grid with the flanking genes in muted colours and the GOI copies in GOI-tree clade colours). See [OUTPUT.md](OUTPUT.md). |
| `--grid_max_models` | `10` | Most GOI gene models drawn for one genome in the anchor grid; a longer tandem array shows `+N more` (every copy keeps its arrow). A drawing setting, set by no preset module. |
| `--no_fragment_variant` | off | Do not write the `*_with_fragments` figures. By default `ModelStatus=fragment` GOI models are hidden in the main figures and drawn in these siblings. |
| `--no_orient_to_home` | off | Do not flip scaffolds that read in the opposite direction to the home genome. |
| `--no_anchor_grid` | off | Skip the anchor grid. |
| `--gene_names_tsv` | _(none)_ | Names table that sets the label of the home genes it lists (`gene_id`, `label`, optional `full_name`, `source`). Start from the `<locus>_gene_names.tsv` a run writes. In `--plot_extra_args` give an absolute path. See [OUTPUT.md § Gene names in the figures](OUTPUT.md#gene-names-in-the-figures). |
| `--no_gene_legend` | off | Do not draw the "Gene names" table under the figures. |
| `--no_network` | off | Skip the NCBI lookups of the plot step: species common names, and the current gene symbols of home genes the GFF does not name. |
| `--max_goi_per_genome` | the run's | Wins over the parameter `--max_goi_per_genome` above; a re-draw uses the run's value from `plot_inputs_*/effective_params.json`, else 10. |
| `--goi_zoom`, `--goi_min_px`, `--caption_file` | — | Tandem-array figure controls. |

`--grid_goi_style` and `--grid_goi_max_models` were removed on 2026-10-07 and are now rejected; see [PARAMETERS.md § `plot_extra_args`](PARAMETERS.md#plot_extra_args).

### Resource Tuning

These control per-process resource allocation. Override them for your hardware.

| Parameter | Default | Description |
|---|---|---|
| `--iterative_search_cpus` | `2` | CPUs for ITERATIVE_SEARCH tasks |
| `--iterative_search_memory` | `auto` | RAM for ITERATIVE_SEARCH tasks: `auto` = max(10 GB, what the largest target genomes searched at once need), x the attempt, at most the local executor's memory; a size is used as it is (8 GB / 4 GB in the `laptop_safe` / `low_mem` tiers; `docker_max` sizes it from host RAM) |
| `--iterative_search_max_forks` | `1` | Max parallel ITERATIVE_SEARCH tasks |
| `--locate_gene_cpus` | `1` | CPUs for LOCATE_GENE |
| `--locate_gene_memory` | `3 GB` | RAM for LOCATE_GENE |
| `--skip_tree` | `false` | Skip the per-locus MAFFT + IQ-TREE phylogeny (`COMPUTE_TREE`) to save time/RAM on weak machines. Emits a placeholder newick (`(GOI_placeholder:0.0);`) so plotting still works. Forced on by the `low_mem` profile (and thus by the launcher's default `auto,low_mem`). |
| `--compute_tree_cpus` | `2` | CPUs for COMPUTE_TREE (MAFFT + IQ-TREE). |
| `--compute_tree_memory` | `4 GB` | RAM for COMPUTE_TREE. |
| `--default_cpus` | `1` | CPUs for every process without its own override. |
| `--default_memory` | `2 GB` | RAM for every process without its own override. |
| `--default_time` | `4h` | Wall-clock limit for every process without its own override. Every executor enforces it, the local one too: a task is killed at its limit (and the 3-4 h limits of the search steps apply to every run), so give a long run a `-c` config with its own `process { time = ... }`. |
| `--gpu_partition` | _(empty)_ | SLURM partition for the GPU-requiring steps (PLM / structural search). Only used by the SLURM profiles. |
| `--gpu_cluster_options` | `--gres=gpu:1` | Extra `sbatch` options for those steps. A GPU job without a `--gres` request sits permanently pending. |
| `--allow_missing_smith_waterman` | `false` | Proceed with Smith-Waterman disabled when parasail is unavailable, instead of failing loud. Degrades divergent-GOI recall — see [§8 parasail](#parasail-import-error-on-startup). |
| `--rank_wave_binning` | `true` | Bin search waves by phylo-distance **rank** rather than absolute distance, so a tight clade whose distances all saturate near 1.0 still grades close→far. `false` restores the pre-2026-07-26 binning. |
| `--disable_wavefront` | `false` | Measurement control, not a performance option: search every genome against the initial database only (one parallel wave), to test whether the expanding wavefront changes a result. |
| `--synteny_bridge_two_sided` | `false` | **(opt-in, unmeasured)** Weigh the flanking genes on both sides of a gap when deciding to bridge it, instead of only the left side. Bigger blocks raise flanking support, so this can promote more calls. |
| `--synteny_bridge_max_per_block` | `2` | Max gaps one block may chain when bridging (`0` = unlimited). |
| `--report_nonsyntenic_candidates` | `true` | Record GOI hits that clear the quality bar but are refused by a synteny gate; they are listed in `synvoy_report.json` under `rejected_candidates`. |
| `--quiet` | `true` | Compact console output. `false` prints the full per-task Nextflow log. |
| `--goi_fallback_intron_margin` | `5.0` | GOI fallback hit-chaining uses the home GOI's *own* largest intron × this margin as the max gap, so two far-apart spurious hits aren't stitched into one gene. |
| `--goi_fallback_intron_floor` | `10000` | Floor (bp) for that derived gap, covering compact / single-exon genes. Falls back to `--max_intron` when the home GOI structure is unavailable. |

### Home-Locus Selection & Multi-Locus Safety

How many places in the *home* genome SynVoy treats as "the gene". Every extra locus
multiplies the whole downstream search, so these are the run-time and OOM knobs.

| Parameter | Default | Description |
|---|---|---|
| `--max_loci` | `5` | Hard cap on distinct home loci to search. Firefly luciferase (an AMP-binding family member) produced 49 loci and a ~55 h run before this existed; the two real loci ranked 1st and 2nd. Raising it is the fastest way to make a run unbounded. |
| `--locus_min_bit_ratio` | `0.5` | Keep a secondary locus only if it scores ≥ this × the primary locus. The score is that of one copy of the gene (all its exons), not of the best single hit: a relative sharing one conserved domain no longer ranks like the gene itself. |
| `--family_warning_threshold` | `10` | Warn when the *pre-filter* locus count exceeds this — a signal that your query is a large family and the results are candidates, not orthologs. |
| `--home_goi_gene` | _(empty)_ | Name of the GOI gene in the home annotation, when the name from the query header or UniProt is not the one used there. |
| `--home_locus_mode` | `auto` | `auto`: the annotated home gene that matches the query is the primary locus, further loci are other places where the query aligns at least half as well. `gene_only`: that locus alone (replaces `--enable_name_locus`). `alignment`: loci from the alignments alone, as before 2026-10-08. |
| `--home_chain_max_intron` | `auto` | Largest gap (bp) between two hits of one copy of the gene at home. `auto`: the largest intron of the annotated genes the hits overlap × 5, within 10 kb to 1.1 Mb (10 kb without an annotation). |
| `--locus_gap_genes` / `--locus_gap` | `auto` / `auto` | Two places where the query aligns are one home locus when at most this many annotated genes lie between them (`auto`: half of `--n_flanking_genes`, 5 at the default). `--locus_gap <bp>` uses a distance instead (50000 when the home genome has no annotation). |

### Rescue Passes

Extra attempts to model the GOI when the main search comes up empty in a neighbourhood
that clearly *is* the right one. All are gated so they are no-ops when the GOI was
already recovered.

| Parameter | Default | Description |
|---|---|---|
| `--strong_synteny_min_flanking` | `auto` | A block with ≥ this many HIGH-confidence flanking genes but no GOI model is surfaced as `goi_missing_but_strong_synteny` (promoted to MEDIUM) instead of being dropped. `auto` = half of `--n_flanking_genes`, rounded up (5 at the default); `0` = off. |
| `--disable_strong_synteny_rescue` | `false` | Turn off the relaxed-miniprot pass inside those blocks. It emits a LOW-confidence model tagged `EvidenceType=relaxed_miniprot_rescue` — *some* model to inspect rather than a silent absence. |
| `--strong_synteny_rescue_miniprot_outc` | `0.05` | Min query coverage (fraction) for the relaxed pass to emit a model. |
| `--strong_synteny_rescue_window_pad` | `2000` | bp added each side of the block coordinates before the rescue search. |
| `--strong_synteny_rescue_miniprot_timeout` | `3600` | Guard of the rescue's miniprot run of one block: a run that does not finish ends the run (until 2026-10-28: 120 s, then the block was skipped). |
| `--disable_goi_hull_rescue` | `false` | Turn off the GOI **hull** rescue — miniprot of the GOI across the whole flanking neighbourhood of a target chromosome, which catches a gene sitting in a rearrangement gap *between* seeded blocks. Largely superseded by collinearity bridging (`--synteny_bridge_max_gap`), which covers the gap up front; kept as a backstop for when home gene order is unavailable. |
| `--goi_hull_min_flanking` | `auto` | Min HIGH-confidence flanking genes in the hull before a rescue is attempted. `auto` = round(0.4 × `--n_flanking_genes`) (4 at the default). |
| `--goi_hull_cluster_max_gap` | `auto` | bp gap that splits flanking genes into separate clusters when computing the hull. `auto` = the bridge limit of the search for that genome (3000000 before). |
| `--goi_hull_window_pad` | `auto` | bp padding around the hull window. `auto` = the GOI margin of the search (100000 before). |
| `--goi_hull_max_window` | `auto` | Skip hulls longer than this (bp). `auto` = twice the home neighbourhood × expansion (20000000 before). |
| `--goi_hull_min_identity` | `40.0` | Min miniprot % identity for a rescue model to be judged (a cost filter). |
| `--goi_hull_min_coverage` | `0.5` | Min CDS aa / query length for a rescue model to be judged (a cost filter; the label reads the query span of the alignment). Settings key: the `short_query` module (B1) sets 0.3, its fragment bound. |
| `--enable_dispersed_goi_rescue` | `true` | Seed a block from a strong GOI hit that landed *outside* every flanking-seeded block, provided the chromosome carries enough flanking anchors to define an envelope around it. |
| `--dispersed_goi_min_identity` | `40.0` | Min % identity for a dispersed GOI hit to qualify. |
| `--dispersed_goi_max_evalue` | `1e-10` | Max e-value for a dispersed GOI hit. |
| `--dispersed_goi_min_alnlen` | `50` | Min alignment length for a dispersed GOI hit (`short_query` / `exon_structure`: 0.6 × the longest exon or the query, 10-50). |
| `--dispersed_min_chrom_anchors` | `3` | The chromosome must carry at least this many flanking genes before an envelope is drawn around the dispersed hit. |
| `--dispersed_envelope_margin` | `auto` | bp padding each side of the flanking envelope. `auto` = `--dispersed_envelope_margin_genes` (9) pitches × expansion (200000 before). |

### Locus Ownership & Paralog Discrimination

Post-hoc checks that a recovered gene really is *your* gene and not its closest
paralog — the decorin-vs-biglycan class of error.

| Parameter | Default | Description |
|---|---|---|
| `--disable_locus_ownership` | `false` | Turn off reciprocal-best-hit assignment of each recovered gene to a home locus. With it on, a call whose best home match is a *family paralog* rather than the GOI is relabelled `paralog_not_goi` and excluded from the headline HIGH/MEDIUM counts. |
| `--require_paralog_panel` | `true` | Fail loud if the home-paralog panel cannot be built (e.g. parasail missing) instead of silently skipping ownership. Set `false` to tolerate a missing panel. |
| `--paralog_panel_min_identity` | `28.0` | Min % identity (query vs a home protein, with a coverage guard) for that protein to enter the paralog panel. The panel is built from the whole home proteome, not just the split-loci spans — otherwise the real paralogs are missing and only spurious domain cross-hits are present. |
| `--paralog_panel_max_family` | `15` | Cap on homology-derived family entries (closest paralogs by identity). Settings key: the `dispersed_family` module (A2) sets max(15, home loci of the family + 1), at most 50. |
| `--locus_ownership_pad` | `2000` | bp padding when intersecting genes with a locus span. |
| `--locus_ownership_max_genes_per_locus` | `5` | Cap panel entries per locus (the annotated genes overlapping the locus most). Settings key: the `local_family` module (A1) sets max(5, home copies + 1), at most 60, so every copy of a home array is in the panel. |
| `--locus_ownership_tiebreak_gap` | `10.0` | SW-score margin below which the RBH result counts as ambiguous and the synteny tiebreak decides. |
| `--disable_paralog_check` | `false` | Turn off the reciprocal-best paralog check. Automatically a no-op for single-sequence queries. |
| `--paralog_confusion_min_gap` | `5` | Min bitscore gap before a call whose best-matching paralog differs from the per-cell modal paralog is flagged `paralog_confusion` in `self_consistency`. |

### Phylogenetic Placement (advisory)

Asks whether a call's sequence divergence from the query fits its species' distance. A
lineage-specific paralog sits further away than its species should. The verdict
(`phylo_concordant` / `phylo_discordant` / `insufficient_data`) is written on each record
in the report; it does **not** change a confidence unless you ask it to. On the 10-family
benchmark it carried no usable signal (every discordant record was a real ortholog), so
treat it as a hint.

| Parameter | Default | Description |
|---|---|---|
| `--disable_phylo_placement` | `false` | Skip the check. |
| `--phylo_placement_min_calls` | `4` | Species with calls needed before the fit is attempted. |
| `--phylo_placement_z_threshold` | `2.0` | Residual z-score at or above which a call is `phylo_discordant`. |
| `--phylo_placement_promote` | `false` | Promote an `AMBIGUOUS` call that is `phylo_concordant` to `MEDIUM`. Not recommended yet. |

### Structural Discovery (experimental, GPU)

Predicts ORFs inside candidate regions, folds them, and matches by structure — for
orthologs that are invisible to sequence search. Requires the ML overlay
(`environment-ml.yml`) and a GPU to be practical. Off by default.

| Parameter | Default | Description |
|---|---|---|
| `--enable_structural_discovery` | `false` | Fold region ORFs to discover sequence-invisible orthologs and materialise them as GOI models. |
| `--discovery_min_tm` | `0.7` | Min Foldseek TM-score to materialise a discovery model. Deliberately conservative. |
| `--discovery_min_len_ratio` | `0.5` | Guard 1 — reject an ORF shorter than this × GOI length. |
| `--discovery_max_len_ratio` | `2.0` | Guard 1 — reject an ORF longer than this × GOI length. |
| `--discovery_min_goi_coverage` | `0.5` | Guard 2 — the SW alignment must span at least this fraction of the GOI. |
| `--discovery_paralog_margin` | `1.0` | Guard 3 — the GOI score must beat the best home-paralog score by this factor. |
| `--discovery_min_block_flanking` | `2` | Guard 4a — at least this many flanking anchors near the ORF window. |
| `--discovery_allow_offblock` | `false` | Guard 4b — allow a rescue off the sequence-GOI chromosome. |

### Output

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `results` | Directory for pipeline output |

### Automatic Parameter Estimation (LLM): removed

The LLM parameter advisor (`--auto_params`, `--force_disable_advanced_search` and the four
`--llm_*` options) was removed on 2026-10-22. On Nextflow 25.10 and later it could not
change a parameter after launch, so it only printed suggested flags. The settings that depend
on the gene come from the preset modules (`--modules`), the distances from the run itself
([Distances in gene units](#distances-in-gene-units)). A run that still sets one of the six
prints a warning and goes on ([PARAMETERS.md § 20](PARAMETERS.md#20-removed-parameters)).

---

## 5. Running on HPC / SLURM

A ready-made submission script is provided in `scripts/slurm_submit.sh`. Edit the variables at the top and submit:

```bash
# Edit scripts/slurm_submit.sh to set your query, species, partition, and account
sbatch scripts/slurm_submit.sh
```

The script submits the Nextflow **controller** as a SLURM job. Nextflow then submits individual pipeline tasks as separate SLURM jobs via the `hpc_singularity` or `hpc_conda` profile.

**Key environment variables for HPC:**

```bash
export NXF_WORK="${SCRATCH}/work"          # fast scratch for intermediates
export NXF_SINGULARITY_CACHEDIR="${HOME}/.singularity/cache"
```

**Manual SLURM example:**

```bash
./run_synvoy.sh \
  -profile hpc_conda \
  --mode easy \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --max_genomes 10 \
  --outdir results/melittin_hpc \
  -work-dir "${SCRATCH}/work" \
  -resume
```

> `scripts/slurm_submit.sh` calls `nextflow run main.nf` directly rather than the
> launcher, because the batch script already pins the environment itself. If you write
> your own submission script, prefer `./run_synvoy.sh` — its pre-flight checks (notably
> the parasail probe) are worth more, not less, on a shared cluster where a stray
> module load can shadow the env.

> Edit the `slurm` profile in `nextflow.config` to change the default partition (`normal`) or add `--account`.

---

## 6. Output Files

All output goes into `--outdir` (default: `results/`):

| Path | Description |
|---|---|
| `*_anchor_grid.html` | **Primary figure.** Species × gene grid (rows = species with the phylo tree at left, columns = home genes; the GOI is the red column). Always emitted. (`bin/plot_synteny.py` has a `--no_anchor_grid` switch, but `modules/plot_synteny.nf` does not forward it — there is no pipeline-level flag to turn the grid off.) |
| `*_synteny_plot.html` | Interactive HTML visualization. Open in a browser — shows syntenic blocks, gene arrows, homology links, and a phylogenetic tree. |
| `*_gene_names.tsv` | The label every home gene carries in the figures, with its full name, gene ID and the source of the label (home GFF, NCBI Gene, or an abbreviation of the product name). Editable: pass it back with `--gene_names_tsv`. |
| `*_tree.nwk` | Newick-format phylogenetic tree of all discovered GOI and GOI-similar sequences across genomes (multiple per genome when paralogs are found). |
| `regions/<locus>/*.regions.bed` | BED files with genomic coordinates of identified candidate syntenic blocks on each target genome, one folder per home locus. |
| `synvoy_report.json` | Machine-readable JSON report: input parameters, genome QC metrics, per-target results, internal exit codes. |
| `intermediate/` | Per-phase artifacts — flanking gene FASTAs, MMseqs2 hit tables, per-target GFFs, miniprot alignments, etc. **Always published** (each stage `publishDir`s into it unconditionally). |
| `downloaded_genomes/` | (Easy Mode only) Downloaded genome assemblies and `assembly_quality.tsv` with contiguity stats. |

---

## 7. Resuming & Caching

Nextflow caches completed tasks in the `work/` directory. To resume after a crash or parameter change:

```bash
./run_synvoy.sh [your params] -resume
```

Only tasks whose inputs changed will be re-executed. This is especially useful for:
- Adding more target genomes to an existing run
- Tweaking visualization parameters (only `PLOT_SYNTENY` will re-run)
- Recovering from transient network failures in Easy Mode

**Cleaning up:**

```bash
# Remove work directory (frees disk space, loses cache)
rm -rf work/

# Or use Nextflow's built-in cleanup
nextflow clean -f
```

---

## 8. Troubleshooting

### Pipeline crashes with SIGKILL (exit code 137) during LOCATE_GENE or ITERATIVE_SEARCH

**Cause:** Out of memory. MMseqs2's translated target database and tBLASTN's index can be RAM-intensive — especially on vertebrate-scale genomes (1-3 GB FASTA → 3-9 GB 6-frame translated DB).

**Fix:**
- A retry of ITERATIVE_SEARCH asks for more memory each time (`iterative_search_memory` `auto`; the console prints `ITERATIVE_SEARCH memory: ...` per attempt). On a laptop the request stops at the machine's memory, so there the next points matter.
- Use the right memory tier — `-profile standard,laptop_safe` (16 GB RAM) or `-profile standard,low_mem` (8 GB RAM). See [§2 — Memory tiers](#memory-tiers-combine-with-an-execution-backend).
- If the error message is `cannot fit database into ... / not enough memory to keep dbreader/write in memory`, drop `--mmseqs_split_memory_limit` further (try `1G`).
- Increase available memory or switch to an HPC profile.
- Check that `/tmp` has sufficient free space (MMseqs2 uses it for temporary files).

### ITERATIVE_SEARCH runs for a very long time (>40 min per genome)

**Cause:** The target genome contains large tandem duplication arrays or many candidate blocks, causing exhaustive local searches.

**What to do:**
- This can be normal for complex genomes. Check progress:
  ```bash
  # Find the task's work directory
  ls -lt work/*/*/.command.log | head -5
  # Follow the log
  tail -f work/<hash>/<hash>/.command.log
  ```
- To speed things up, reduce `--max_blocks_per_genome` (default 80) or raise `--min_block_genes` (default 2)
- In a genome rich in gene families (plants, genomes with old whole-genome duplications) each flanking gene hits dozens of places, and two chance hits near each other make a block: the log then says `Capping blocks: N -> 80`, and a large part of the genome is searched (a fifth of *Arabidopsis lyrata* for one locus of *A. thaliana*). `--anchor_max_places 5` seeds blocks only at the five best places of each flanking gene; in that case it searched one window and found the same gene. With any `--modules` other than `legacy` this happens on its own in each genome where the cap binds (the log says `genome axis: the block cap (80) bound at N blocks; rebuilding ...`).
- In a genome where genes lie far apart (a vertebrate home locus), the derived distances make every block and every search window large: the block gap is 7 gene spacings of the home locus, which is about 2 Mb where the spacing is 300 kb. The log line `distances: home pitch ... -> block gap ..., window margin ... (DNA searched: ... per side)` shows what was derived for each genome. The derived distances have been tested on yeast, fly, bee, ant, beetle and plant genomes (up to 0.6 Gb) and not yet on vertebrate genomes. If a run is too slow there, use `--anchor_max_places 5` (previous point); `-profile <execution>,legacy_distances` restores the former fixed distances, or lower `--block_gap_genes`.

### "Query FASTA not found" or "Home genome not found"

**Cause:** Paths are relative to the Nextflow launch directory, not to the script.

**Fix:** Use absolute paths or ensure you run `nextflow` from the SynVoy project root.

### Conda environment creation times out

**Cause:** Conda solver is slow.

**Fix:**
- Install [Mamba](https://github.com/mamba-org/mamba) and create the environment with `mamba env create -f environment.yml`
- Or increase the timeout: the config already sets `conda.createTimeout = '3 h'`
- Re-create the environment if Python 3.13 was selected by an older environment file; SynVoy pins Python `<3.13` because `ete3` still imports the removed `cgi` module.

### Easy Mode fails to download genomes

**Cause:** NCBI API rate limiting or network issues.

**Fix:**
- Set an NCBI API key: `export NCBI_API_KEY=your_key_here`
- Re-run with `-resume` — completed downloads will be cached
- Check NCBI service status at https://www.ncbi.nlm.nih.gov/

### No synteny plot is generated

**Cause:** No candidate regions passed the synteny score threshold.

**Fix:**
- Increase `--n_flanking_genes` (e.g. `15`) to capture more genomic context
- Increase `--block_gap_genes` (default 7 gene pitches) if the flanking genes lie further apart in the target than at home
- Check `regions/*.bed` files — if they are empty, the flanking-gene mapping step did not find hits

### "Java not found" or Nextflow fails to start

**Fix:**
```bash
# Check Java version
java -version

# If missing, install via Conda
mamba install -n synvoy_env -c conda-forge 'openjdk>=17'

# Or via system package manager
sudo apt install default-jdk   # Debian/Ubuntu
```

Current Nextflow releases require Java 17 or newer. The SynVoy Conda
environment includes OpenJDK, so `mamba run -n synvoy_env nextflow ...`
is a good fallback when your login shell does not have Java on `PATH`.

### `ete3` / `cgi` import error

**Symptom:** Python fails with `ModuleNotFoundError: No module named 'cgi'`
while importing `ete3`.

**Cause:** The environment was solved with Python 3.13. `cgi` was removed
from the Python standard library, and current `ete3` still imports it.

**Fix:**
```bash
mamba env remove -n synvoy_env
mamba env create -f environment.yml
conda activate synvoy_env
python -V   # should be >=3.10 and <3.13
```

### `iqtree2` not found

Some current IQ-TREE packages install the binary as `iqtree` instead of
`iqtree2`. SynVoy tries `iqtree2` first and falls back to `iqtree`, so this
is only a problem if neither command exists.

```bash
iqtree2 --version || iqtree --version
```

### Pipeline finishes with `synvoy_report.json` showing 0 annotations / 0 regions

**Symptom:** `GENERATE_REPORT` exits non-zero with a message like
"zero annotation and zero region files under staged_results".

**Cause:** Either `ITERATIVE_SEARCH` genuinely produced no hits, or the
Nextflow channel wiring did not stage the expected files into
`staged_results/`.

**How to diagnose:**
1. Open `synvoy_report.json` (it is still written even on failure). The
   `staging_diagnostics` block lists per-directory entry counts and sample
   filenames:
   ```bash
   jq '.staging_diagnostics' results/<run>/synvoy_report.json
   ```
2. If `match_counts` shows zero across the board **and** the sample
   entries are empty or only contain sentinels like `NO_REGIONS`, the
   upstream search truly found nothing. Inspect
   `logs/ITERATIVE_SEARCH__*/stdout.log` for the per-genome hit counts.
3. If the dirs contain files but none match the expected patterns
   (`*.gff`, `*.scores.tsv`), the module's channel wiring is wrong —
   check `modules/generate_report.nf` stageAs directives.

**If zero-hit is genuinely expected** (e.g. you're deliberately testing
a query with no orthologs in your targets), the Nextflow driver will
still fail the `GENERATE_REPORT` process. Re-run the reporter manually
against the staged output:
```bash
python3 bin/generate_report.py \
    --results_dir work/<generate_report_hash>/staged_results \
    --output synvoy_report.json \
    --allow-empty
```

### `parasail` import error on startup

**Symptom:** `ModuleNotFoundError: No module named 'parasail'` during
`ITERATIVE_SEARCH` or `smith_waterman_search.py`.

**Most common cause:** a Python virtualenv (`.venv`) is active and shadows the
conda env that Nextflow's tasks inherit. VS Code auto-activates one in its
integrated terminal. Run `deactivate`, then relaunch — `./run_synvoy.sh` probes
for exactly this before it starts anything.

**Fix:** Activate the SynVoy conda environment (it ships parasail via
`environment.yml`). If running outside conda:
```bash
pip install parasail
```
If `pip install parasail` fails to build, fall back to ssearch36:
```bash
conda install -c bioconda fasta3
./run_synvoy.sh ... --sw_method ssearch36
```

> **This is a hard failure, not a warning.** Smith-Waterman is the load-bearing
> tier for divergent GOI recovery (it is what re-finds a 26–40 % myrmicitoxin from
> a melittin query), so `ITERATIVE_SEARCH` **aborts** when SW is requested in
> `auto` mode and parasail is missing. Earlier versions silently disabled SW and
> produced quietly degraded results; that behaviour was removed deliberately. To
> proceed anyway — accepting worse recall on divergent queries — pass
> `--allow_missing_smith_waterman true`.

### GPU out-of-memory (OOM) during STRUCTURAL_SEARCH / ESMFold

**Symptom:** `CUDA out of memory. Tried to allocate X GiB`.

**Cause:** Your GPU's VRAM is smaller than what ESMFold needs for the
current query length. ESMFold scales quadratically in sequence length.

**Fix:**
- Lower `--structural_max_length` to 200 or 150 for <6 GB GPUs.
- For the 4 GB GTX 1650 class, the safe ceiling is ~150 aa.
- If your query is longer than the cap, disable structural search
  entirely — synteny + Smith-Waterman is usually enough:
  ```
  --enable_structural_search false
  ```

### `-resume` reruns every process instead of caching

**Symptom:** You expected `-resume` to skip completed processes, but
every task runs again from scratch.

**Cause:** Nextflow caches by content hash. Any of these invalidates
the cache for a given process and all downstream processes:
1. An input file's path changed, even if its content is identical.
2. A parameter changed (including defaults, if you touched
   `nextflow.config`).
3. You moved `work/` or deleted specific subdirectories.
4. You switched profiles (`-profile laptop_safe` vs. `-profile standard`).

**How to diagnose:**
```bash
# Show the cache-lookup result for each task
./run_synvoy.sh ... -resume -dump-hashes
```
Look for `CACHE FOUND` vs. `not found`. The first `not found` tells
you which process's inputs changed; everything downstream is forced
to re-run.

**Fix:** Revert the path/parameter change, or accept the re-run. Never
delete `work/` if you need `-resume` to work — delete it only when you
genuinely want a clean slate.

### Conda env build fails on macOS (Augustus)

**Symptom:** `conda env create -f environment.yml` fails while building
Augustus on macOS.

**Cause:** Augustus does not have a maintained conda recipe for Apple
Silicon / recent macOS.

**Fix:** Use the Docker container instead:
```bash
docker build -t synvoy .
nextflow run main.nf ... -profile docker
```
Or use `-profile singularity` on Linux systems where Docker is not
available. This is the one case where you call Nextflow directly:
`./run_synvoy.sh` aborts when `synvoy_env` is absent, which is exactly
the situation you are working around here.

### A parameter I set seems to be ignored

Check these in order:

1. **Was it replaced by a preset?** Open `intermediate/locate_gene/effective_params.json`.
   `source` says, per parameter, whether the value came from the default, a preset, a
   preset module (`module:A1`) or you (`user_set`), and `derivation` says what set each
   value that is not the default (for a module rule, from which inputs). Since 2026-10-05 a value you set
   wins; on an older checkout the auto-selected preset overwrote it — `git pull`.
2. **Was the parameter removed?** Parameters that no step read were removed on
   2026-10-09 (`--keep_intermediate`, `--max_retries`, `--multi_profile`, …): the run
   prints a warning and ignores them. A replaced parameter (`--max_flanking_goi_similarity`,
   `--enable_name_locus`) stops the run at start and names its replacement. The list is in
   [PARAMETERS.md § 20](PARAMETERS.md#20-removed-parameters).
3. **Is it a script-level option?** Figure controls such as `--print_grid_numbers` belong to
   `bin/plot_synteny.py`; pass them through `--plot_extra_args '…'`. (`--max_goi_per_genome`
   is a parameter of its own since 2026-10-26; given in `--plot_extra_args` it still wins.)
4. **Did you pass `false` to switch something off?** That works since 2026-10-05. On an
   older checkout `--require_paralog_panel false`, `--disable_goi_hull_rescue false` and
   a few others were read as *true*, because Nextflow hands every command-line value over
   as text.

### `SyntaxError: f-string expression part cannot include a backslash` in PLOT_SYNTENY

Your environment has Python 3.10 or 3.11 and your checkout dates from before 2026-10-05,
when one line of `bin/plot_synteny.py` only parsed on Python 3.12. `git pull`.

### `UnboundLocalError: cannot access local variable 'strong_synteny'` in CLUSTER_REGIONS

Fixed on 2026-08-21. Your checkout is older — `git pull` (and check which branch you are
on: the fix reached `main` later than `dev`).

### The run ends with "Tasks FAILED (ignored): N"

The rescue, paralog-check and ownership steps are allowed to fail without stopping the
run, because one failed rescue should not cost you a whole result. When any did fail the
summary says so. The main search results are complete; for the genomes concerned the
rescue or the paralog check is missing. Look in `logs/` for the step named
`RESCUE_GOI_HULL`, `RESCUE_STRONG_SYNTENY`, `RECIPROCAL_BEST_PARALOG` or
`ASSIGN_LOCUS_OWNERSHIP`.

### "File name(s) with whitespace or special characters"

See the file-name note under [Pro Mode](#pro-mode): rename the files, then relaunch.

### Tracking down a specific process failure

Nextflow keeps every task's working directory under `work/<hash>/`.
To find the directory for a failed task:
```bash
# From the error output, note the hash (e.g. "Work dir: /path/work/a7/8f3...")
ls /path/work/a7/8f3*/
# Inspect:
cat .command.log    # stderr from the process
cat .command.sh     # the exact command that was run
cat .command.out    # stdout from the process
```
This is usually the fastest way to understand why a single genome's
task failed without tripping the whole pipeline.

### Still stuck?

1. Check your environment first — `./run_synvoy.sh --check-only` verifies conda,
   Java, Nextflow, the config, and parasail without launching a run, and prints
   the version you are on. (There is no `--help` flag: use §4 above, or
   [PARAMETERS.md](PARAMETERS.md) for the annotated reference.)
2. See [QUICKSTART.md](QUICKSTART.md) for an end-to-end worked example on
   bee genomes (~20–30 min, mostly downloads).
3. Open an issue at https://github.com/AndreasWz/SynVoy/issues with
   your `nextflow run` command, the contents of
   `synvoy_report.json`'s `staging_diagnostics` block, and the
   relevant `work/*/.command.log`.
