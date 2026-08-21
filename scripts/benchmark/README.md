# Benchmark scripts

One-shot research code for the SynVoy melittin benchmark. Not part of the
pipeline — nothing here is invoked by `main.nf`.

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
├── score_benchmark.py          # joins tool outputs against truth, emits confusion matrices
└── _archived_proteome_tools/   # dropped v1 wrappers (OrthoFinder, SonicParanoid,
                                #   MMseqs2-RBH, GENESPACE) — kept for reference

tests/benchmark_truth/
├── melittin_orthologs.tsv      # curated ground truth
├── melittin_orthologs_gr2.tsv  # GR2 (flanking-neighbourhood) variant
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
