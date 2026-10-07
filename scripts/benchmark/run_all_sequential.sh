#!/usr/bin/env bash
# Sequential benchmark runner. v2 architecture: genome-based competitor tools
# only (proteome-based tools were dropped because of NCBI annotation gaps —
# see _archived_proteome_tools/).
#
# Usage:
#   tmux new -s bench
#   mamba activate synvoy_benchmark
#   bash scripts/benchmark/run_all_sequential.sh
#   # detach: Ctrl-b d ; reattach: tmux attach -t bench
#
# Skips any tool whose calls.tsv already exists.

set -uo pipefail

LOG_DIR="benchmark_results/_logs"
mkdir -p "${LOG_DIR}"

run_step() {
    local name="$1"
    local cmd="$2"
    local log="${LOG_DIR}/${name}.log"
    local out="benchmark_results/${name}/calls.tsv"

    if [[ -s "${out}" ]]; then
        echo "[${name}] calls.tsv exists, skipping. Delete to force re-run."
        return 0
    fi

    echo
    echo "=== [${name}] starting at $(date -Iseconds) ==="
    echo "    log: ${log}"
    if bash -c "${cmd}" 2>&1 | tee "${log}"; then
        echo "=== [${name}] DONE at $(date -Iseconds) ==="
    else
        echo "=== [${name}] FAILED — see ${log} ===" >&2
        return 1
    fi
}

# 0. Make sure all 21 target genomes are on disk and 3 new ant proteomes are
#    optional (only Cardiocondyla has a GFF among the new ants — others are FASTA-only).
run_step prep_proteomes "bash scripts/benchmark/prep_proteomes.sh" || true

# 1. tblastn — universal baseline (fast, ~5-15 min)
run_step tblastn "bash scripts/benchmark/run_tblastn.sh" || true

# 2. MMseqs2 sens=9.5 — paper's exact methodology (fast, ~10-20 min)
run_step mmseqs_genome "bash scripts/benchmark/run_mmseqs_genome.sh" || true

# 3. miniprot — splice-aware protein-to-genome (fast, ~5-10 min)
run_step miniprot "bash scripts/benchmark/run_miniprot.sh" || true

# 4. SynVoy Pro Mode — needs synvoy_env for nextflow (long, ~30-60 min for 3 new ants
#    via -resume; full re-run would be 6-12h but cached state should kick in).
run_step synvoy "bash scripts/benchmark/run_synvoy.sh" || true

# 5. MCScanX — synteny-based competitor. Will mark targets without GFF (e.g.
#    Tetramorium, Formica) as ABSENT with `mcscanx_no_annotation` — that's
#    the point of including it.
run_step mcscanx "bash scripts/benchmark/run_mcscanx.sh" || true

# 6a. Score by SPECIES PRESENCE (the historical number, all tools)
echo
echo "=== scoring [species presence] at $(date -Iseconds) ==="
python scripts/benchmark/score_benchmark.py \
    --truth tests/benchmark_truth/melittin_orthologs.tsv \
    --calls-glob 'benchmark_results/*/calls.tsv' \
    --outdir benchmark_results/ 2>&1 | tee "${LOG_DIR}/score.log"

# 6b. Score by COORDINATE OVERLAP (docs/TODO.md §1y) — SynVoy only, because it is
#     the only tool here that emits full gene models. Presence scoring recorded the
#     ant melittin "recovery" as a clean TP when it was 132 bp of an 855 bp gene,
#     and did not notice when a later run lost it entirely. Both numbers get
#     published; neither is reported alone.
echo
echo "=== scoring [coordinate overlap, SynVoy] at $(date -Iseconds) ==="
SYNVOY_RUN="${RESULTS_DIR:-results/melittin_benchmark}"
if [[ -d "${SYNVOY_RUN}" ]]; then
    python scripts/benchmark/score_coordinates.py \
        --truth tests/benchmark_truth/melittin_loci.tsv \
        --run "${SYNVOY_RUN}" \
        --out benchmark_results/coordinate_per_model.tsv \
        --out-calls benchmark_results/coordinate_per_call.tsv \
        2>&1 | tee "${LOG_DIR}/score_coordinates.log"
else
    echo "  SKIPPED: no SynVoy run dir at ${SYNVOY_RUN}" | tee "${LOG_DIR}/score_coordinates.log"
fi

echo
echo "=== all steps complete at $(date -Iseconds) ==="
echo
echo "Confusion per tool [species presence]:"
column -t -s $'\t' benchmark_results/confusion_per_tool.tsv 2>/dev/null
echo
echo "Coordinate recall [SynVoy] — read this next to the table above:"
grep -A6 "READ THIS NUMBER" "${LOG_DIR}/score_coordinates.log" 2>/dev/null
echo
echo "TOGA runs are separate (need the heavier toga_setup.sh install + 22GB+ RAM):"
echo "  bash scripts/benchmark/run_toga_tier3.sh    # TOGA1 (generates chains)"
echo "  bash scripts/benchmark/run_toga2_tier3.sh   # TOGA2 (reuses chains)"
