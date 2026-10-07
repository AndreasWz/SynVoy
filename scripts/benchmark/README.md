# Benchmark scripts

Research code for benchmarking SynVoy. Not part of the pipeline — nothing here is
invoked by `main.nf`. There are two benchmarks:

1. **Melittin vs. other tools** — one gene, 19 Hymenoptera genomes, SynVoy against
   tblastn / MMseqs2 / miniprot / MCScanX. Scored twice: by *species presence* and by
   *coordinate overlap* with curated gene models (see below — the second exists
   because the first let two overclaims through).
2. **Family recovery** — ten curated venom-gene families across 33 ant genomes from
   one honeybee seed each, with the reachable set fixed in advance. Design and
   results: [`docs/NEXT_SESSION_FAMILY_BENCHMARK.md`](../../docs/NEXT_SESSION_FAMILY_BENCHMARK.md).

The curated coordinates (`tests/benchmark_truth/`) are unpublished and git-ignored,
so both scorers need the truth files on your machine; the scoring *logic* is tested
in CI against synthetic data (`tests/test_coordinate_benchmark.py`,
`tests/test_family_recovery.py`).

**v2 architecture (current):** genome-based competitors only. The proteome-based
tools (OrthoFinder, SonicParanoid, MMseqs2-RBH, GENESPACE) were dropped because
NCBI annotation coverage is too patchy across the target set to compare fairly —
their wrappers are kept for reference under `_archived_proteome_tools/`. Excluding
a tool because half the targets have no proteome is a methodological choice, not a
gap: the same annotation gap is why MCScanX is *kept*, as a synteny competitor that
visibly fails on unannotated genomes.

## Layout

```
scripts/benchmark/
├── README.md                   # this file
├── competitors.yml             # conda env for the competitor tools
├── fetch_benchmark_genomes.sh  # downloads the target genomes via NCBI Datasets
├── prep_proteomes.sh           # stages genomes/proteomes, reports which targets have a GFF
├── verify_gffs.py              # sanity-checks the staged annotations
├── run_all_sequential.sh       # the entry point — runs every tool, then scores
├── run_tblastn.sh              # universal baseline
├── run_mmseqs_genome.sh        # MMseqs2 sens=9.5, the paper's methodology
├── run_miniprot.sh             # splice-aware protein-to-genome
├── run_synvoy.sh               # SynVoy Pro Mode
├── run_mcscanx.sh              # synteny competitor (marks unannotated targets ABSENT)
├── run_toga_tier3.sh           # TOGA1 — separate, needs toga_setup.sh + 22 GB+ RAM
├── run_toga2_tier3.sh          # TOGA2 — reuses TOGA1's chains
├── toga_setup.sh               # heavier TOGA install
├── run_tp53.sh                 # the TP53 paralog-discrimination side benchmark
├── fetch_tp53_truth.py         # builds the TP53 truth table
├── score_benchmark.py          # species-presence scoring: tool outputs vs truth, confusion matrices
├── score_coordinates.py        # coordinate scoring of a SynVoy run against curated gene models
├── build_melittin_loci.py      # builds the coordinate truth table from the curated GFFs
├── family_recovery.py          # family benchmark harness: plan / stage / seedcheck / score / callcheck
├── run_family.sh               # one family end to end: stage -> SynVoy -> score -> callcheck
├── score_runs.sh               # re-score finished family runs with the current harness
└── _archived_proteome_tools/   # dropped v1 wrappers (OrthoFinder, SonicParanoid,
                                #   MMseqs2-RBH, GENESPACE) — kept for reference

tests/benchmark_truth/          # git-ignored: unpublished curation
├── melittin_orthologs.tsv      # presence truth, one row per species
├── melittin_orthologs_gr2.tsv  # GR2 (flanking-neighbourhood) variant
├── melittin_loci.tsv           # coordinate truth, one row per curated gene model
├── family_preregistration*.tsv # family benchmark: reachable set per family, hashed before the runs
└── tp53_orthologs.tsv          # TP53 paralog truth

local_data/benchmark/
└── genomes/                    # populated by fetch_benchmark_genomes.sh
    ├── home/Apis_mellifera/
    ├── tier1/Apis_cerana/
    ├── tier1/Apis_florea/
    ├── ...
    ├── tier2/Megachile_rotundata/
    └── tier3/Vollenhovia_emeryi/

benchmark_results/              # populated by the tool wrappers
├── _logs/                      # per-step logs from run_all_sequential.sh
├── synvoy/calls.tsv
├── tblastn/calls.tsv
├── ...
├── master_table.tsv            # written by score_benchmark.py
├── confusion_per_tool.tsv
└── confusion_per_tier.tsv
```

## Normalized tool output contract

Every `run_<tool>.sh` must emit `benchmark_results/<tool>/calls.tsv` with these
columns (TSV, header required):

| column | type | meaning |
|---|---|---|
| `species` | string | snake_case species name matching truth file |
| `accession` | string | NCBI accession or `local` |
| `called_status` | enum | `PRESENT`, `ABSENT`, or `AMBIGUOUS` |
| `locus_chrom` | string | chrom/scaffold ID, or `-` |
| `locus_start` | int | 1-based, or `-` |
| `locus_end` | int | 1-based, or `-` |
| `strand` | enum | `+`, `-`, or `-` (placeholder) |
| `confidence` | string | tool-native (e.g. `HIGH`, `1e-30`, `0.95`) |
| `extra` | string | free text notes |

`score_benchmark.py` joins these against `melittin_orthologs.tsv` and emits
TP/FP/FN/TN per tool and per phylogenetic tier.

## Running the full benchmark

```bash
# Phase A — fetch data
bash scripts/benchmark/fetch_benchmark_genomes.sh

# Phase B — install the competitor tools
mamba env create -f scripts/benchmark/competitors.yml

# Phase C+D — run every tool, then score. Long; run it detached.
tmux new -s bench
mamba activate synvoy_benchmark
bash scripts/benchmark/run_all_sequential.sh
# detach: Ctrl-b d   reattach: tmux attach -t bench
```

`run_all_sequential.sh` skips any tool whose `calls.tsv` already exists — delete
that file to force a re-run — and finishes by calling the scorer itself:

```bash
python scripts/benchmark/score_benchmark.py \
    --truth tests/benchmark_truth/melittin_orthologs.tsv \
    --calls-glob 'benchmark_results/*/calls.tsv' \
    --outdir benchmark_results/
```

TOGA is **not** part of that sequence — it needs the heavier install and 22 GB+ RAM:

```bash
bash scripts/benchmark/toga_setup.sh
bash scripts/benchmark/run_toga_tier3.sh    # TOGA1 (generates chains)
bash scripts/benchmark/run_toga2_tier3.sh   # TOGA2 (reuses those chains)
```

## Coordinate scoring (melittin)

Species-presence scoring asks "did the tool say PRESENT where the truth says
PRESENT". It recorded the ant melittin recovery as a clean true positive when the
call covered 132 bp of an 855 bp gene, and did not notice when a later run lost that
gene entirely, because another call on the same scaffold kept the species flag up.
`score_coordinates.py` scores a call by **where it landed** and reports the fraction
of each curated model that was recovered:

```bash
python scripts/benchmark/score_coordinates.py \
    --truth tests/benchmark_truth/melittin_loci.tsv \
    --run results/melittin_benchmark \
    --out benchmark_results/coordinate_per_model.tsv \
    --out-calls benchmark_results/coordinate_per_call.tsv
```

`run_all_sequential.sh` runs both scorers; report the two numbers together. A
curated model on a different assembly than the run used is reported as
`assembly_mismatch` and excluded, not counted as a miss.

## Family-recovery benchmark

```bash
# 1. Before any run: fix what counts as reachable, and commit the hash of the output
python scripts/benchmark/family_recovery.py --data <Ant_Venoms>/data plan

# 2. One family, end to end (stage -> SynVoy with default flags -> score -> callcheck)
SYNVOY_FAMILY_DATA=<Ant_Venoms>/data \
    scripts/benchmark/run_family.sh FAM_APYR <SEED> local_runs/family_bench

# 3. Re-score finished runs after a harness change
scripts/benchmark/score_runs.sh local_runs/family_bench/FAM_*__*
```

What the harness reports per family: recall at locus, address and family level over
the genes the seed can reach; per-call precision; and for every miss the pipeline
layer that lost it (A: no search block on the locus, B: block but no model, C:
model but confidence too low, D: right gene, wrong label). Genes at loci the home
genome does not have are `no_home_anchor` and are never counted — a hit there is
reported as `found_without_anchor`, not as a win. Calls are matched to truth genes
by coding-sequence overlap.

`run_family.sh` records anything passed after the output folder in
`<run>/run_args.txt`, so a tuned run cannot be mistaken for a default one. The
Slurm wrappers used on the cluster are site-specific and not in the repository.
