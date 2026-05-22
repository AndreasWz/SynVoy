# Benchmark scripts

One-shot research code for the SynVoy melittin benchmark. Not part of the
pipeline. See [`/BENCHMARK_PLAN.md`](/BENCHMARK_PLAN.md) for the full plan.

## Layout

```
scripts/benchmark/
├── README.md                       # this file
├── fetch_benchmark_genomes.sh      # downloads 13 new genomes via NCBI Datasets
├── score_benchmark.py              # joins tool outputs against truth, emits confusion matrices
├── run_synvoy.sh                   # (TODO) wraps SynVoy Pro Mode run on the 18 targets
├── run_orthofinder.sh              # (TODO) wraps OrthoFinder v3 on extracted proteomes
├── run_sonicparanoid.sh            # (TODO) wraps SonicParanoid2
├── run_mmseqs_rbh.sh               # (TODO) MMseqs2 reciprocal best hits baseline
├── run_genespace.sh                # (TODO) GENESPACE on annotated genomes
└── run_toga2.sh                    # (OPTIONAL) TOGA2, gated on install success

tests/benchmark_truth/
└── melittin_orthologs.tsv          # curated ground truth (18 species)

local_data/benchmark/
└── genomes/                        # populated by fetch_benchmark_genomes.sh
    ├── home/Apis_mellifera/
    ├── tier1/Apis_cerana/
    ├── tier1/Apis_florea/
    ├── ...
    ├── tier2/Megachile_rotundata/
    └── tier3/Vollenhovia_emeryi/

benchmark_results/                  # populated by tool wrappers
├── synvoy/calls.tsv
├── orthofinder/calls.tsv
├── ...
├── master_table.tsv                # written by score_benchmark.py
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
# Phase A — fetch data + verify
bash scripts/benchmark/fetch_benchmark_genomes.sh
python scripts/benchmark/verify_gffs.py   # TODO

# Phase B — install competitors
conda env create -f scripts/benchmark/competitors.yml   # TODO

# Phase C — run all tools
bash scripts/benchmark/run_synvoy.sh
bash scripts/benchmark/run_orthofinder.sh
bash scripts/benchmark/run_sonicparanoid.sh
bash scripts/benchmark/run_mmseqs_rbh.sh
bash scripts/benchmark/run_genespace.sh

# Phase D — score
python scripts/benchmark/score_benchmark.py \
    --truth tests/benchmark_truth/melittin_orthologs.tsv \
    --calls-glob 'benchmark_results/*/calls.tsv' \
    --outdir benchmark_results/
```
