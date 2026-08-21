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
| `--bad_quality_policy` | `drop` | What to do with low-quality assemblies: `drop`, `keep`, or `ask` (interactive prompt) |
| `--bad_quality_timeout` | `300` | Seconds to wait for user input when policy is `ask` |
| `--bad_max_contigs` | `500000` | Assemblies with more contigs are flagged as low quality |
| `--bad_max_scaffolds` | `500000` | Assemblies with more scaffolds are flagged as low quality |
| `--bad_min_n50` | `5000` | Assemblies with N50 below this are flagged as low quality |

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
> **What it buys you:** matched target models gain `TargetGene` / `TargetProduct` / `TargetID` attributes in the output GFF and the `target_gene` / `target_product` report columns — i.e. you can see *which annotated gene* a call landed on. This is read-out evidence, not a scoring input: annotations do not move synteny blocks, region scores, or confidence. The only mechanism that acts on them is the family-consistency gate (`--strict_goi_family`, default `false`), which can **demote** a weak call whose target gene name disagrees with `--goi_family_tokens`. If you enable it, set `--goi_family_tokens` for your gene explicitly — the automatic derivation does not work (see PARAMETERS.md).

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
| `auto` | local | Conda | **Recommended zero-config default.** Like `standard`, plus it pins the full auto-tuning stack on (auto preset selection, locus cap, strong-synteny rescue, paralog check). Best when you don't want to hand-pick a preset. Does *not* enable the LLM advisor (`--auto_params true` is still opt-in). |
| `standard` | local | Conda | Default. 2 CPUs / 10 GB RAM per iterative-search task, single-fork. Good baseline for workstations. |
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
| (none — uses defaults) | ≥16 GB free RAM | 10 GB | 8G | 9.5 |
| `laptop_safe` | ~16 GB RAM laptop | 8 GB | 4G | 8.0 |
| `low_mem` | ~8 GB RAM laptop | 4 GB | 2G | 7.0 |

Trade-off: tighter tiers are slower (more MMseqs splits, less parallelism) and slightly less sensitive on highly divergent queries. They're the right pick when you're seeing `cannot fit database into ... / not enough memory to keep dbreader/write in memory` errors from MMseqs2 — that means your `mmseqs_split_memory_limit` is too high for the genome you're searching.

### Biology presets (combine with an execution backend + memory tier)

| Profile | Use for |
|---|---|
| `preset_single_copy` | Housekeeping / single-copy conserved orthologs (e.g. SLRP genes DCN/FMOD/ASPN across vertebrates, ACTB, ribosomal proteins). |
| `preset_short_peptide` | Short, divergent peptides (e.g. mature melittin, defensins, signal-peptide-only queries). |
| `preset_tandem_family` | Tandem-duplicated families with many close paralogs (e.g. luciferases, opsins, cytochrome P450 clusters). |
| `preset_paralog_discrimination` | When the goal is to tell paralogs apart, not just find any family member (e.g. distinguishing TP53 from TP63/TP73). |

> **Auto-apply (since 2026-05-29, default ON):** if you don't pick a preset, SynVoy picks one for you after `LOCATE_GENE` from the hit distribution and applies it at runtime to `EXTRACT_FLANKING` / `ITERATIVE_SEARCH` / `CLUSTER_REGIONS` via a value channel — no restart needed. The chosen preset is logged as `[INFO] §1f preset applied: …` and dumped to `intermediate/locate_gene/effective_params.json`. Pin one with `--preset_override preset_X`, opt out with `--auto_apply_preset false`. The one caveat: `NORMALIZE_QUERY.min_query_length` is resolved before LOCATE_GENE runs, so short queries (< 30 aa) still need `--min_query_length 20` or `-profile preset_short_peptide` set at launch. See [docs/PARAMETERS.md §19](PARAMETERS.md#19-auto-apply-preset-self-consistency--rescue) for the full toggle list.

---

## 3. Algorithm Overview

The pipeline proceeds through five phases:

### Phase 1 — Gene Localization

1. **Normalize Query:** If the input is DNA, the best ORF is translated to protein.
2. **Locate GOI:** The query protein is aligned against the home genome using tblastn and MMseqs2 to establish coordinates.
3. **Annotate GOI Exons:** If a GFF is available, the GOI is matched to an annotated gene and individual CDS/exons are extracted. Otherwise, exon boundaries are inferred from alignment hits.
4. **Split Loci:** If the GOI maps to multiple genomic locations (e.g. tandem duplicates), each locus is processed independently.
5. **Extract Flanking Genes:** The *n* genes upstream and downstream of each locus are identified from GFF or Augustus/Prodigal prediction. Flanking candidates that are similar to the GOI (above `--max_flanking_goi_similarity`) are excluded to avoid inflating synteny scores. A `--max_flanking_distance` cap can prevent walking into distant gene deserts.
6. **Expand GOI-Similar Neighbors** *(optional, off by default)*: When `--expand_goi_similar` is enabled, genes near the GOI that resemble it (e.g. tandem duplicates like MRJPs near Yellow-e3) are emitted as additional GOI queries with a `GOI_NEIGHBOR_` prefix. These are searched in all target genomes alongside the original GOI, and included in the phylogenetic tree — enabling resolution of paralogs vs. orthologs.
7. **Borrow Annotations:** When the home genome lacks a GFF, annotations can be borrowed from annotated target genomes via reciprocal best hits.

### Phase 2 — Phylogenetic Ordering & Iterative Search

7. **Phylo Sort:** Target genomes are ordered by evolutionary distance to the reference.
8. **Genome Quality Assessment:** Target assemblies are evaluated for contiguity (N50, scaffold count).
9. **Iterative Search:** For each target genome (nearest-first), flanking genes are mapped with MMseqs2. Hits are clustered into candidate syntenic blocks. Within each block, localized tblastn, miniprot, and Smith-Waterman searches attempt to find the GOI. Discovered genes are added to the search database, improving sensitivity for more distant species.

### Phase 3 — Region Clustering

10. **Cluster Regions:** Candidate blocks across all targets are filtered by synteny score and ranked.

### Phase 4 — Phylogenetics & Visualization

11. **Compute Tree:** All discovered GOI and GOI-similar sequences across all genomes are aligned (MAFFT) and a phylogenetic tree is inferred (IQ-TREE with automatic model selection and ultrafast bootstrap). Multiple hits per genome are preserved, so the tree can resolve paralogs from orthologs.
12. **Plot Synteny:** An interactive HTML plot shows the syntenic context of each hit, colored by homology, with the phylogenetic tree alongside.

### Phase 5 — Reporting

13. **Generate Report:** A JSON summary file is produced with run parameters, genome QC results, and per-target outcomes.

---

## 4. Full Parameter Reference

All parameters can be set on the command line (`--param value`) or in a custom config file (`-c my.config`).

### Synteny & Search

| Parameter | Default | Description |
|---|---|---|
| `--n_flanking_genes` | `10` | Number of non-GOI-similar flanking genes to extract on each side of the GOI |
| `--prefer_large_genes` | `true` | Prefer larger flanking genes (more informative for homology search) |
| `--min_flanking_size` | `500` | Minimum size (bp) for a flanking gene to be included |
| `--max_flanking_goi_similarity` | `35.0` | Exclude flanking genes with k-mer similarity (%) to the GOI above this threshold. Prevents GOI paralogs (e.g. tandem duplicates) from being used as synteny anchors. Set to `100` to disable. |
| `--max_flanking_distance` | `0` | Max distance (bp) from GOI center to walk for flanking genes. `0` = unlimited. Useful when the GOI neighbours a large tandem array. |
| `--expand_goi_similar` | `false` | Emit GOI-similar flanking genes as additional GOI queries (`GOI_NEIGHBOR_` prefix). Enables paralog discovery and phylogenetic resolution across genomes, but can flood paralog-rich searches; opt in deliberately. |
| `--expand_goi_similar_distance` | `300000` | Max distance (bp) from GOI to search for GOI-similar neighbor genes |
| `--exon_level_search` | `true` | Search at exon level for better divergent-gene detection |
| `--cluster_distance` | `150000` | Max gap (bp) between flanking-gene hits to merge into one syntenic block |
| `--min_synteny_score` | `0.6` | Fraction of flanking genes that must map to a target to trigger local search |
| `--min_hit_identity` | `10` | Minimum alignment identity (%) for an individual hit |
| `--min_hit_length` | `10` | Minimum alignment length for an individual hit |
| `--min_query_length` | `30` | Reject queries shorter than this (aa) during `NORMALIZE_QUERY`. Set `0` to disable (e.g. searching a short motif or micro-exon). |
| `--search_evalue` | `0.01` | E-value threshold for tblastn/MMseqs2 searches |
| `--max_intron` | `20000` | Maximum intron length (bp) for miniprot gene models |
| `--region_padding` | `150000` | Extra flanking sequence (bp) appended to each side of a candidate block |
| `--padding_min` | `50000` | Minimum padding (bp) |
| `--padding_max` | `200000` | Maximum padding (bp) |
| `--max_blocks_per_genome` | `80` | Safety cap on candidate blocks per target genome |
| `--min_block_genes` | `2` | Minimum flanking-gene hits in a block to keep it |
| `--max_consecutive_empty_blocks` | `25` | Stop expanding after this many consecutive empty blocks |
| `--disable_synteny_collinearity` | `false` | Disable collinearity-aware block scoring + gap-bridging (revert to legacy gene-count clustering). When on (default), blocks are ranked by their longest run of home-ordered flanking genes and a GOI sitting in a rearrangement gap between two collinear clusters is kept inside one searched block. |
| `--synteny_bridge_max_gap` | `6000000` | Max bp gap between two same-chromosome flanking clusters to bridge into one block when they collinearly continue the home gene order. `0` disables bridging. |
| `--synteny_bridge_max_rank_gap` | `5` | Max home-rank jump still treated as a collinear continuation when bridging. |
| `--synteny_bridge_min_anchors` | `3` | Min anchored flanking genes a block needs before it may bridge a gap. |
| `--seed_on_flanking_support` | `false` | **(opt-in)** Also seed the expanding query DB (and thus later waves) with a HIGH/MEDIUM GOI feature regardless of its `goi_class` (e.g. a `tandem_goi_copy`) when flanking support and query coverage clear the floors below. Lets a close relative's perfect/divergent ortholog bridge to a more divergent species in the next wave. |
| `--seed_flanking_min_count` | `2` | Min flanking-gene support required for `--seed_on_flanking_support`. |
| `--seed_flanking_min_qcov` | `0.5` | Min query coverage required for `--seed_on_flanking_support`. |

### Smith-Waterman Local Search

| Parameter | Default | Description |
|---|---|---|
| `--enable_smith_waterman` | `true` | Use rigorous Smith-Waterman alignment (parasail) for GOI search |
| `--sw_method` | `auto` | Implementation: `auto`, `parasail`, or `ssearch36` |
| `--sw_min_score` | `20` | Minimum SW alignment score to report a hit |
| `--sw_min_identity` | `10.0` | Minimum identity (%) for SW hits |
| `--sw_timeout_seconds` | `300` | Timeout per SW search invocation |

### Relaxed / Augmented Search

Controls the increasingly permissive search passes used for highly divergent targets.

| Parameter | Default | Description |
|---|---|---|
| `--aug_relaxed_evalue_mult` | `1000` | Multiply base e-value by this factor in relaxed passes |
| `--aug_relaxed_evalue_cap` | `10.0` | Maximum e-value allowed even in relaxed mode |
| `--aug_relaxed_identity_factor` | `0.6` | Multiply normal identity threshold by this in relaxed mode |
| `--aug_relaxed_identity_min` | `15.0` | Absolute minimum identity (%) in relaxed mode |
| `--aug_relaxed_parse_evalue_mult` | `10` | Secondary e-value multiplier used when parsing relaxed-pass hits |
| `--aug_relaxed_length_div` | `2` | Divide normal length threshold by this in relaxed mode |
| `--aug_relaxed_length_min` | `15` | Absolute minimum alignment length in relaxed mode |
| `--aug_dedup_bin_bp` | `100` | Bin size (bp) for deduplicating overlapping relaxed hits |
| `--fallback_short_query_len` | `150` | Queries shorter than this (aa) gate their GOI fallback on aligned length + bitscore instead of query coverage — a short mature peptide in a longer precursor (e.g. melittin) is structurally low-coverage even on a perfect hit. `0` = legacy coverage gate for all queries. |
| `--fallback_short_min_aln_aa` | `15` | Minimum aligned length (aa) for a short-query GOI fallback. Rejects micro-window overcalls while admitting short mature peptides. |
| `--fallback_short_min_bits` | `30` | Minimum bitscore for a short-query GOI fallback. Lets a high-bitscore divergent hit pass where low coverage would have vetoed it. |

### MMseqs2

| Parameter | Default | Description |
|---|---|---|
| `--mmseqs_sensitivity` | `9.5` | MMseqs2 sensitivity (1–10+). Higher = slower but more sensitive |
| `--mmseqs_split_memory_limit` | `8G` | MMseqs2 memory limit for database splitting. Override to `3G` or `1G` on memory-constrained machines. |
| `--mmseqs_verbosity` | `1` | MMseqs2 log verbosity (0 = silent) |
| `--min_gene_identity` | `30` | Minimum identity (%) for flanking-gene MMseqs2 matches |
| `--deterministic_goi_search` | `true` | Force the per-region augmented GOI MMseqs2 search to `--threads 1` so marginal/divergent GOI calls are reproducible run-to-run (regions are tiny, so the speed cost is small). Set `false` for the legacy multithreaded (non-deterministic) path. |

### Annotation & Gene Prediction

| Parameter | Default | Description |
|---|---|---|
| `--gff_search_window` | `100000` | Window (bp) around GOI to search in GFF for flanking genes |
| `--gap_search_window` | `50000` | Window for gap-filling searches |
| `--gap_min_size` | `10` | Minimum gap size (bp) to attempt fill |
| `--gap_evalue` | `10` | E-value for gap search |
| `--gap_min_identity` | `15.0` | Minimum identity (%) for gap hits |
| `--gap_min_alnlen` | `10` | Minimum alignment length for gap hits |
| `--gap_max_hits` | `5` | Max gap hits to report |
| `--min_exon_query_cov` | `0.25` | Minimum query coverage fraction for exon annotation |
| `--min_exon_alnlen` | `30` | Minimum exon alignment length |
| `--pred_flank_window` | `50000` | Prodigal prediction window around locus |
| `--pred_keep_pct` | `0.10` | Fraction of Prodigal predictions to keep |
| `--prodigal_full_genome_fallback` | `false` | Run Prodigal on entire genome if windowed prediction fails |
| `--gene_predictor` | `auto` | Gene predictor for unannotated genomes: `auto`, `augustus`, or `prodigal` |
| `--augustus_species` | `fly` | Augustus species model used when Augustus is selected |

### Gene Model Classification

Controls the confidence labels (HIGH/MEDIUM/LOW) and model status labels (complete/partial/fragment) assigned to gene models in GFF output.

| Parameter | Default | Description |
|---|---|---|
| `--classify_high_min_identity` | `50.0` | Min identity (%) for HIGH-confidence exon_annotation models (lowered from 60 for cross-vertebrate realism) |
| `--classify_medium_min_identity` | `35.0` | Min identity (%) for MEDIUM-confidence exon_annotation models (lowered from 45 for divergent orthologs) |
| `--classify_high_min_collinear` | `3` | Once a block has ≥3 flanking genes, HIGH additionally requires a collinear run of home-ordered flanking at least this long — a scrambled paralog neighbourhood is capped at MEDIUM even at high identity. `0` restores the legacy count-only behaviour. |
| `--strict_goi_family` | `false` | Downgrade fallback/rescued_exon/raw_hit GOI calls whose annotated `TargetGene`/`TargetProduct` does not contain a family token. Useful for multi-paralog queries (e.g. TP53 family). |
| `--goi_family_tokens` | _(auto)_ | Comma-separated family name tokens for `--strict_goi_family`. If empty, auto-derived from query FASTA header (`GN=`, UniProt entry name). |
| `--classify_tandem_min_identity` | `40.0` | Min identity (%) for MEDIUM-confidence tandem copies. Below this, tandem copies are labeled LOW. |
| `--classify_tandem_min_qcov` | `0.35` | Min query coverage for a MEDIUM tandem copy. A high-identity but short local window (its true low coverage now recorded) is demoted to LOW — stops a 7–20 aa micro-window from masquerading as a full tandem copy. |
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

### Visualization

| Parameter | Default | Description |
|---|---|---|
| `--plot_width` | `1500` | Width of the output HTML plot (px) |
| `--gap_threshold` | `50000` | Gaps larger than this (bp) are visually compressed |
| `--gap_visual_size` | `20000` | Size (bp) used to represent compressed gaps |
| `--flank_fallback_bp` | `1000000` | Maximum window (bp) rendered around distal targets |
| `--scale_bar_len` | `10000` | Scale bar size (bp) |
| `--hide_goi_absent_tracks` | `true` | Omit target tracks from plots when no GOI hit was found in that genome |
| `--pub_svg` | `false` | Also emit the narrow publication SVG (`*_synteny_plot.svg`, vertical/condensed for a two-column journal figure). The `*_synteny_plot_view.svg` mirror of the HTML is emitted either way. |
| `--pub_width_mm` | `183` | Publication SVG width in mm (183 = Nature double column). |
| `--pub_palette` | `okabe_ito` | Colour palette for the publication SVG. Default is the colourblind-safe Okabe–Ito set. |
| `--enable_matrix_plot` | `false` | Also render the phylogeny-anchored presence/absence matrix (`*_synteny_matrix.svg`). Off because the anchor grid covers the same ground more compactly. |

### Resource Tuning

These control per-process resource allocation. Override them for your hardware.

| Parameter | Default | Description |
|---|---|---|
| `--iterative_search_cpus` | `2` | CPUs for ITERATIVE_SEARCH tasks |
| `--iterative_search_memory` | `10 GB` | RAM for ITERATIVE_SEARCH tasks (lowered to 8 GB / 4 GB by the `laptop_safe` / `low_mem` tiers; `docker_max` sizes it from host RAM) |
| `--iterative_search_max_forks` | `1` | Max parallel ITERATIVE_SEARCH tasks |
| `--locate_gene_cpus` | `1` | CPUs for LOCATE_GENE |
| `--locate_gene_memory` | `3 GB` | RAM for LOCATE_GENE |
| `--skip_tree` | `false` | Skip the per-locus MAFFT + IQ-TREE phylogeny (`COMPUTE_TREE`) to save time/RAM on weak machines. Emits a placeholder newick (`(GOI_placeholder:0.0);`) so plotting still works. Forced on by the `low_mem` profile (and thus by the launcher's default `auto,low_mem`). |
| `--compute_tree_cpus` | `2` | CPUs for COMPUTE_TREE (MAFFT + IQ-TREE). |
| `--compute_tree_memory` | `4 GB` | RAM for COMPUTE_TREE. |
| `--default_cpus` | `1` | CPUs for every process without its own override. |
| `--default_memory` | `2 GB` | RAM for every process without its own override. |
| `--default_time` | `4h` | Wall-clock limit for every process without its own override (matters on SLURM). |
| `--gpu_partition` | _(empty)_ | SLURM partition for the GPU-requiring steps (PLM / structural search). Only used by the SLURM profiles. |
| `--gpu_cluster_options` | `--gres=gpu:1` | Extra `sbatch` options for those steps. A GPU job without a `--gres` request sits permanently pending. |
| `--allow_missing_smith_waterman` | `false` | Proceed with Smith-Waterman disabled when parasail is unavailable, instead of failing loud. Degrades divergent-GOI recall — see [§8 parasail](#parasail-import-error-on-startup). |
| `--force_disable_advanced_search` | `false` | Prevent `--auto_params` from switching on the optional PLM / structural search layers, whatever the estimated evolutionary distance. |
| `--rank_wave_binning` | `false` | Bin search waves by phylo-distance **rank** rather than absolute distance, so a tight clade whose distances all saturate near 1.0 still grades close→far. |
| `--goi_fallback_intron_margin` | `5.0` | GOI fallback hit-chaining uses the home GOI's *own* largest intron × this margin as the max gap, so two far-apart spurious hits aren't stitched into one gene. |
| `--goi_fallback_intron_floor` | `10000` | Floor (bp) for that derived gap, covering compact / single-exon genes. Falls back to `--max_intron` when the home GOI structure is unavailable. |

### Home-Locus Selection & Multi-Locus Safety

How many places in the *home* genome SynVoy treats as "the gene". Every extra locus
multiplies the whole downstream search, so these are the run-time and OOM knobs.

| Parameter | Default | Description |
|---|---|---|
| `--max_loci` | `5` | Hard cap on distinct home loci to search. Firefly luciferase (an AMP-binding family member) produced 49 loci and a ~55 h run before this existed; the two real loci ranked 1st and 2nd. Raising it is the fastest way to make a run unbounded. |
| `--locus_min_bit_ratio` | `0.5` | Keep a secondary locus only if its best bitscore is ≥ this × the primary locus's. Ranking by bitscore ratio rather than e-value is deliberate — e-values saturate at 0.0 for a strong query and stop discriminating. |
| `--family_warning_threshold` | `10` | Warn when the *pre-filter* locus count exceeds this — a signal that your query is a large family and the results are candidates, not orthologs. |
| `--enable_name_locus` | `false` | **Opt-in.** Establish the home locus by looking the GOI's gene symbol up in `--home_gff` instead of re-discovering it by self-alignment. For an annotated, named gene in a family this avoids anchoring the search on a paralog's neighbourhood (the failure that filed asporin as decorin). Falls back loudly to alignment-locate if the query↔gene consistency check fails. |
| `--home_goi_gene` | _(empty)_ | Force the GOI gene symbol for the lookup above; implies `--enable_name_locus`. Otherwise the symbol resolved from UniProt is used. |
| `--name_locus_min_identity` | `60.0` | Min % identity (query vs the looked-up gene's protein) to accept the name match. Guards against a typo, the wrong species, or a paralog symbol. |
| `--name_locus_pad` | `0` | bp padding either side of the gene span in the emitted locus BED. |

### Rescue Passes

Extra attempts to model the GOI when the main search comes up empty in a neighbourhood
that clearly *is* the right one. All are gated so they are no-ops when the GOI was
already recovered.

| Parameter | Default | Description |
|---|---|---|
| `--strong_synteny_min_flanking` | `5` | A block with ≥ this many HIGH-confidence flanking genes but no GOI model is surfaced as `goi_missing_but_strong_synteny` (promoted to MEDIUM) instead of being dropped. |
| `--disable_strong_synteny_rescue` | `false` | Turn off the relaxed-miniprot pass inside those blocks. It emits a LOW-confidence model tagged `EvidenceType=relaxed_miniprot_rescue` — *some* model to inspect rather than a silent absence. |
| `--strong_synteny_rescue_miniprot_outc` | `0.05` | Min query coverage (fraction) for the relaxed pass to emit a model. |
| `--strong_synteny_rescue_window_pad` | `2000` | bp added each side of the block coordinates before the rescue search. |
| `--strong_synteny_rescue_miniprot_timeout` | `120` | Seconds per block for the rescue miniprot call. |
| `--disable_goi_hull_rescue` | `false` | Turn off the GOI **hull** rescue — miniprot of the GOI across the whole flanking neighbourhood of a target chromosome, which catches a gene sitting in a rearrangement gap *between* seeded blocks. Largely superseded by collinearity bridging (`--synteny_bridge_max_gap`), which covers the gap up front; kept as a backstop for when home gene order is unavailable. |
| `--goi_hull_min_flanking` | `4` | Min HIGH-confidence flanking genes in the hull before a rescue is attempted. |
| `--goi_hull_cluster_max_gap` | `3000000` | bp gap that splits flanking genes into separate clusters when computing the hull. |
| `--goi_hull_window_pad` | `100000` | bp padding around the hull window. |
| `--goi_hull_max_window` | `20000000` | Skip hull windows larger than this (bp) — stops a far rearranged outlier from ballooning the search. |
| `--goi_hull_min_identity` | `40.0` | Min miniprot % identity to emit a rescued model. |
| `--goi_hull_min_coverage` | `0.5` | Min query coverage (CDS aa / query length) for a rescued model. |
| `--enable_dispersed_goi_rescue` | `true` | Seed a block from a strong GOI hit that landed *outside* every flanking-seeded block, provided the chromosome carries enough flanking anchors to define an envelope around it. |
| `--dispersed_goi_min_identity` | `40.0` | Min % identity for a dispersed GOI hit to qualify. |
| `--dispersed_goi_max_evalue` | `1e-10` | Max e-value for a dispersed GOI hit. |
| `--dispersed_goi_min_alnlen` | `50` | Min alignment length for a dispersed GOI hit. |
| `--dispersed_min_chrom_anchors` | `3` | The chromosome must carry at least this many flanking genes before an envelope is drawn around the dispersed hit. |
| `--dispersed_envelope_margin` | `200000` | bp padding each side of the flanking envelope. |

### Locus Ownership & Paralog Discrimination

Post-hoc checks that a recovered gene really is *your* gene and not its closest
paralog — the decorin-vs-biglycan class of error.

| Parameter | Default | Description |
|---|---|---|
| `--disable_locus_ownership` | `false` | Turn off reciprocal-best-hit assignment of each recovered gene to a home locus. With it on, a call whose best home match is a *family paralog* rather than the GOI is relabelled `paralog_not_goi` and excluded from the headline HIGH/MEDIUM counts. |
| `--require_paralog_panel` | `true` | Fail loud if the home-paralog panel cannot be built (e.g. parasail missing) instead of silently skipping ownership. Set `false` to tolerate a missing panel. |
| `--paralog_panel_min_identity` | `28.0` | Min % identity (query vs a home protein, with a coverage guard) for that protein to enter the paralog panel. The panel is built from the whole home proteome, not just the split-loci spans — otherwise the real paralogs are missing and only spurious domain cross-hits are present. |
| `--paralog_panel_max_family` | `15` | Cap on homology-derived family entries (closest paralogs by identity). |
| `--locus_ownership_pad` | `2000` | bp padding when intersecting genes with a locus span. |
| `--locus_ownership_max_genes_per_locus` | `5` | Cap panel entries per locus (keeps tandem clusters from dominating). |
| `--locus_ownership_tiebreak_gap` | `10.0` | SW-score margin below which the RBH result counts as ambiguous and the synteny tiebreak decides. |
| `--locus_ownership_synteny_window` | `200000` | bp window for that flanking co-location tiebreak. |
| `--disable_paralog_check` | `false` | Turn off the reciprocal-best paralog check. Automatically a no-op for single-sequence queries. |
| `--paralog_confusion_min_gap` | `5` | Min bitscore gap before a call whose best-matching paralog differs from the per-cell modal paralog is flagged `paralog_confusion` in `self_consistency`. |

### Distance Auto-Tuning

Relaxes the HIGH/MEDIUM identity bars for distant targets, using the observed flanking
identity as the divergence estimate.

| Parameter | Default | Description |
|---|---|---|
| `--disable_distance_autotune` | `false` | Keep one global stringency for every target instead. |
| `--distance_autotune_close_pct` | `70.0` | Median flanking identity ≥ this → no relaxation. |
| `--distance_autotune_far_pct` | `40.0` | Median flanking identity ≤ this → full relaxation, and calls are tagged `manual_review`. |
| `--distance_autotune_max_relax` | `10.0` | Max identity-points subtracted from the HIGH/MEDIUM bars. |
| `--distance_autotune_min_flanking` | `3` | Minimum flanking genes needed before the estimate is trusted. |

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
| `--keep_intermediate` | `false` | **Currently inert.** Declared in `nextflow.config` but read by no process — `intermediate/` is published either way. Kept so existing command lines don't break. |

### Automatic Parameter Estimation (LLM)

SynVoy can automatically estimate optimal search parameters based on the biological context of your search — query gene size, home species genome architecture, and evolutionary distance to the targets. Without an API key it falls back to deterministic heuristics encoding the same biological rules.

**Disabled by default** (requires an API key to be useful). Enable with `--auto_params true`.

| Parameter | Default | Description |
|---|---|---|
| `--auto_params` | `false` | Enable automatic parameter estimation. When on, SynVoy analyzes your query gene, home species genome architecture, and target species distances to set optimal values for ~25 search parameters. |
| `--llm_provider` | `google` | LLM provider: `google` (Google Gemini) or `openai` (OpenAI or any OpenAI-compatible API). |
| `--llm_api_key` | _(empty)_ | API key for the chosen provider. Also read from `LLM_API_KEY`, `GOOGLE_API_KEY` (for google), or `OPENAI_API_KEY` (for openai) env vars. |
| `--llm_api_base_url` | _(empty)_ | Custom API base URL for OpenAI-compatible providers (Together, Groq, LM Studio, etc.). Ignored for google. |
| `--llm_model` | _(empty)_ | Model override. Defaults: `gemini-2.5-flash-lite` (google), `gpt-4o-mini` (openai). |
| `--multi_profile` | `false` | For small searches, run with multiple parameter profiles (sensitive/balanced/stringent) and automatically select the best result. |
| `--multi_profile_max_jobs` | `30` | Max total jobs (`loci × targets × 3`) allowed for multi-profile. If exceeded, only the LLM-estimated profile runs. |

**Setup:**

```bash
# Option A: Google Gemini (free tier at aistudio.google.com/apikey)
export GOOGLE_API_KEY=your_key_here
./run_synvoy.sh ... --auto_params true

# Option B: OpenAI
export OPENAI_API_KEY=your_key_here
./run_synvoy.sh ... --auto_params true --llm_provider openai

# Option C: OpenAI-compatible endpoint (Together, Groq, LM Studio, etc.)
export LLM_API_KEY=your_key_here
./run_synvoy.sh ... --auto_params true --llm_provider openai \
    --llm_api_base_url https://api.together.xyz \
    --llm_model meta-llama/Llama-3.1-8B-Instruct-Turbo
```

> **No API key?** `--auto_params true` still runs but uses built-in heuristic rules.
> These cover kingdom-specific intron lengths, distance-adaptive sensitivity, and
> query-size thresholds. Solid for most cases; the LLM adds nuance for edge cases.

**What gets estimated (~25 parameters):**
- **Genome architecture**: `max_intron`, `cluster_distance`, `region_padding` — adapted for plants (↑), bacteria (↓), vertebrates, fungi
- **Search sensitivity**: `mmseqs_sensitivity`, `search_evalue`, `min_hit_identity` — relaxed for distant searches, tightened for close species
- **Query-size tuning**: `sw_min_score`, `min_hit_length` — lowered for small peptides, raised for large proteins
- **Gene family handling**: `max_flanking_goi_similarity`, `expand_goi_similar` — tuned for tandem arrays and gene families
- **Advanced search**: `enable_plm_search`, `enable_structural_search` — auto-enabled for extreme evolutionary distances (>400 Mya)

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
| `*_tree.nwk` | Newick-format phylogenetic tree of all discovered GOI and GOI-similar sequences across genomes (multiple per genome when paralogs are found). |
| `regions/*.regions.bed` | BED files with genomic coordinates of identified candidate syntenic blocks on each target genome. |
| `synvoy_report.json` | Machine-readable JSON report: input parameters, genome QC metrics, per-target results, internal exit codes. |
| `intermediate/` | Per-phase artifacts — flanking gene FASTAs, MMseqs2 hit tables, per-target GFFs, miniprot alignments, etc. **Always published** (each stage `publishDir`s into it unconditionally); `--keep_intermediate` does not gate this. |
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
- To speed things up, reduce `--max_blocks_per_genome` (default 80) or increase `--min_synteny_score` (default 0.6)

### "Query FASTA not found" or "Home genome not found"

**Cause:** Paths are relative to the Nextflow launch directory, not to the script.

**Fix:** Use absolute paths or ensure you run `nextflow` from the SynVoy project root.

### Conda environment creation times out

**Cause:** Conda solver is slow.

**Fix:**
- Install [Mamba](https://github.com/mamba-org/mamba) and create the environment with `mamba env create -f environment.yml`
- Or increase the timeout: the config already sets `conda.createTimeout = '1 h'`
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
- Lower `--min_synteny_score` (e.g. `0.3`)
- Increase `--n_flanking_genes` (e.g. `15`) to capture more genomic context
- Increase `--region_padding` to widen the search window
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
   `logs/iterative_search/*.log` for the per-genome hit counts.
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

### LLM parameter advisor: API key error / no key set

**Symptom:** `LLM_PARAM_ADVISOR` logs "No LLM API key provided" and uses heuristics,
or fails with an HTTP 401 / 403 error.

**Cause:** `--auto_params true` is set but no API key was provided, or the key is invalid.

**Fix (pick one):**
- **Recommended for students / reproducibility:** disable LLM auto-params entirely:
  ```
  --auto_params false --multi_profile false
  ```
  The heuristic fallback is solid and covers most common cases.
- **Google Gemini** (free tier available at [aistudio.google.com/apikey](https://aistudio.google.com/apikey)):
  ```
  export GOOGLE_API_KEY=your_api_key
  ./run_synvoy.sh ... --auto_params true
  ```
- **OpenAI:**
  ```
  export OPENAI_API_KEY=your_api_key
  ./run_synvoy.sh ... --auto_params true --llm_provider openai
  ```

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
2. See [QUICKSTART.md](QUICKSTART.md) for a <15 min end-to-end worked
   example on small bee genomes.
3. Open an issue at https://github.com/AndreasWz/SynVoy/issues with
   your `nextflow run` command, the contents of
   `synvoy_report.json`'s `staging_diagnostics` block, and the
   relevant `work/*/.command.log`.
