# Benchmark scripts

Research code for benchmarking SynVoy. Not part of the pipeline - nothing here is
invoked by `main.nf`. There are two benchmarks:

1. **Melittin vs. other tools** - one gene, 19 Hymenoptera genomes, SynVoy against
   tblastn / MMseqs2 / miniprot / MCScanX. Scored twice: by *species presence* and by
   *coordinate overlap* with curated gene models (see below).
2. **Family recovery** - curated gene families of one clade, each searched from one
   seed, with the reachable set fixed in advance. The curated data are not in the
   repository.

The curated coordinates (`tests/benchmark_truth/`) are unpublished and git-ignored,
so both scorers need the truth files on your machine; the scoring *logic* is tested
in CI against synthetic data (`tests/test_coordinate_benchmark.py`,
`tests/test_family_recovery.py`).

**Competitors:** genome-based tools only. The proteome-based tools (OrthoFinder,
SonicParanoid, MMseqs2-RBH, GENESPACE) are not compared, because NCBI annotation
coverage is too patchy across the target set to compare fairly; their wrappers are
kept for reference under `_archived_proteome_tools/`. Excluding
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
├── run_all_sequential.sh       # the entry point - runs every tool, then scores
├── run_tblastn.sh              # universal baseline
├── run_mmseqs_genome.sh        # MMseqs2 genome-wide, sensitivity 9.5, E <= 1e-3
├── run_miniprot.sh             # splice-aware protein-to-genome
├── run_synvoy.sh               # SynVoy Pro Mode
├── run_mcscanx.sh              # synteny competitor (marks unannotated targets ABSENT)
├── run_toga_tier3.sh           # TOGA1 - separate, needs toga_setup.sh + 22 GB+ RAM
├── run_toga2_tier3.sh          # TOGA2 - reuses TOGA1's chains
├── toga_setup.sh               # heavier TOGA install
├── run_tp53.sh                 # the TP53 paralog-discrimination side benchmark
├── fetch_tp53_truth.py         # builds the TP53 truth table
├── score_benchmark.py          # species-presence scoring: tool outputs vs truth, confusion matrices
├── score_coordinates.py        # coordinate scoring of a SynVoy run against curated gene models
├── build_melittin_loci.py      # builds the coordinate truth table from the curated GFFs
├── family_recovery.py          # family benchmark harness: plan / stage / seedcheck / score / callcheck
├── run_family.sh               # one family end to end: stage -> SynVoy -> score -> callcheck
├── score_runs.sh               # re-score finished family runs with the current harness
├── validate_melettin_gt_v3.py  # a melittin run vs the pinned fixture in tests/ground_truth_test/melettin_gt_v3/
├── run_ly6_ground_truth_check.py # LY6 ground-truth check of a finished run (--gt-dir; calls validate_ly6_results.py)
├── validate_ly6_results.py     # GOI calls vs a ground-truth GFF: CDS tracked, false positives
└── _archived_proteome_tools/   # wrappers of the proteome-based tools (OrthoFinder,
                                #   SonicParanoid, MMseqs2-RBH, GENESPACE), kept for reference

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

SynVoy's `calls.tsv` is written by `synvoy_calls.py` (called by `run_synvoy.sh` and
`run_tp53.sh`) from the report's records (`goi_dedup.records`). `PRESENT`: a record labelled
HIGH or MEDIUM that the report counts (its `goi_class` is not `paralog_not_goi` or
`identity_coverage_decoupled`, and it is no answer of the genome-wide path). `AMBIGUOUS`: else a
record labelled AMBIGUOUS that needs curation (not `paralog_not_goi`). `ABSENT`: neither. LOW
models (class `ambiguous_goi_family_member`) and regions set no status. The locus is the first
such record's span, 1-based; `extra` gives the counts and the top region for information. A
genome is named as in the report (its genome id: a GCA_/GCF_ accession, else the file name without
its suffixes), for the region BEDs too, so a `.fasta` target or an NCBI file name gives one row. A
`calls.tsv` written before 16 November 2026 read the per-genome class counts and the region names
(LOW rows gave `AMBIGUOUS`, uncounted rows and region labels gave `PRESENT`, and the locus was the
top region's 0-based start), so it is not comparable with one written now.

## Running the full benchmark

```bash
# Phase A - fetch data
bash scripts/benchmark/fetch_benchmark_genomes.sh

# Phase B - install the competitor tools
mamba env create -f scripts/benchmark/competitors.yml

# Phase C+D - run every tool, then score. Long; run it detached.
tmux new -s bench
mamba activate synvoy_benchmark
bash scripts/benchmark/run_all_sequential.sh
# detach: Ctrl-b d   reattach: tmux attach -t bench
```

`run_all_sequential.sh` skips any tool whose `calls.tsv` already exists - delete
that file to force a re-run - and finishes by calling the scorer itself:

```bash
python scripts/benchmark/score_benchmark.py \
    --truth tests/benchmark_truth/melittin_orthologs.tsv \
    --calls-glob 'benchmark_results/*/calls.tsv' \
    --outdir benchmark_results/
```

TOGA is **not** part of that sequence - it needs the heavier install and 22 GB+ RAM:

```bash
bash scripts/benchmark/toga_setup.sh
bash scripts/benchmark/run_toga_tier3.sh    # TOGA1 (generates chains)
bash scripts/benchmark/run_toga2_tier3.sh   # TOGA2 (reuses those chains)
```

## Coordinate scoring (melittin)

Species-presence scoring asks "did the tool say PRESENT where the truth says
PRESENT". It counts a call that covers a small part of a gene as a full recovery,
and a call elsewhere on the same scaffold keeps the species flag up when the gene
itself is missed. `score_coordinates.py` scores a call by **where it landed** and reports the fraction
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

`build_melittin_loci.py` writes the truth table `tests/benchmark_truth/melittin_loci.tsv` from
`melittin_orthologs.tsv` and the curated ant GFFs. `--ant-dir` names the folder of those GFFs
(default: `$SYNVOY_ANT_ANNOTATIONS`, else `../Ant_Venoms/data/annotations` next to this
repository); the script stops when the folder does not exist.

## LY6 ground-truth check

`run_ly6_ground_truth_check.py` compares the GOI calls of a finished LY6 run with curated
LY6/3FTx gene models. The curated models are not part of the repository; `--gt-dir` names
the folder that holds them (`Ly6/*.gff`, or the `Ly6` folder itself):

```bash
python scripts/benchmark/run_ly6_ground_truth_check.py \
    --gt-dir <folder with Ly6/*.gff> --outdir <SynVoy run folder>
```

Each curated file is matched to a genome's GFF in `<run>/plot_inputs_synteny_block_locus_1/`
by its sequence name and scored by `validate_ly6_results.py` (CDS bases tracked, GOI calls
outside every curated gene). `ly6_validation_report.json` and `.md` go to `--report-dir`
(default: the run folder), never into the truth folder. A `--gt-dir` without curated GFFs, or
a run folder without `plot_inputs_synteny_block_locus_1/`, stops the script with exit 2.

## Family-recovery benchmark

```bash
# 1. Before any run: fix what counts as reachable, and commit the hash of the output
python scripts/benchmark/family_recovery.py --data <DATA> plan

# 2. One family, end to end (stage -> SynVoy with default flags -> score -> callcheck)
SYNVOY_FAMILY_DATA=<DATA> \
    scripts/benchmark/run_family.sh <FAMILY> <SEED> local_runs/family_bench

# 3. Re-score finished runs after a harness change
scripts/benchmark/score_runs.sh local_runs/family_bench/FAM_*__*
```

`<DATA>` is the curated data folder, with `final_dataset/` and `genomes/`; `<FAMILY>`
is a family id from it, and `<SEED>` the home gene it is searched from.

What the harness reports per family: recall at locus, address and family level over
the genes the seed can reach; per-call precision; and for every miss the pipeline
layer that lost it (A: no search block on the locus, B: block but no model of its own,
C: model but confidence too low, D: right gene, wrong label). Genes at loci the home
genome does not have are `no_home_anchor` and are never counted - a hit there is
reported as `found_without_anchor`, not as a win. Calls are matched to truth genes
by coding-sequence overlap, and a call recovers only the gene it overlaps best
(`best_curated_gene`): one hit chain over a tandem array recovers one gene.

`run_family.sh` records anything passed after the output folder in
`<run>/run_args.txt`, so a tuned run cannot be mistaken for a default one. The
Slurm wrappers used on the cluster are site-specific and not in the repository.
