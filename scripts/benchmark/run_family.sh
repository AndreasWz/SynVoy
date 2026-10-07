#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# One family of the family-recovery benchmark, end to end: stage -> SynVoy -> score
# -> callcheck. Design: docs/NEXT_SESSION_FAMILY_BENCHMARK.md. Harness: family_recovery.py.
#
#   run_family.sh <FAMILY> <SEED> <OUT_ROOT> [extra nextflow args ...]
#
# SEED is a curated home gene model (GM_...), or for an address seed (Amendment 2) a
# home-genome gene id, or the name of an unannotated home ORF given by SEED_FAA +
# SEED_COORDS. SEED_LABEL: the curated label the seed stands for, from the address audit.
# e.g. (local, conda off, tools from synvoy_env):
#   PATH=/home/faw/miniforge3/envs/synvoy_env/bin:$PATH \
#   scripts/benchmark/run_family.sh FAM_APYR <SEED> local_runs/family_bench \
#       -c local_data/dcn_rerun/noconda.config --iterative_search_cpus 4
#
# Seeds (gitignored, unpublished): tests/benchmark_truth/family_preregistration*.tsv
# Env: SYNVOY_FAMILY_DATA (Ant_Venoms data/ dir), WORK_ROOT (Nextflow work parent,
# default <OUT_ROOT>/work), PYTHON (default python3), SEED_LABEL / SEED_FAA /
# SEED_COORDS (see above).
#
# Every run uses the SAME pipeline flags: the default configuration is what is being
# benchmarked. Anything passed after OUT_ROOT is recorded in <run>/run_args.txt, so a
# tuned run can never be mistaken for a default one.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
[ $# -ge 3 ] || { sed -n '2,24p' "$0"; exit 2; }
FAM=$1; SEED=$2; ROOT=$3; shift 3

REPO=$(cd "$(dirname "$0")/../.." && pwd)
PY=${PYTHON:-python3}
DATA=${SYNVOY_FAMILY_DATA:-$REPO/../Ant_Venoms_git/data}
HARNESS="$REPO/scripts/benchmark/family_recovery.py"
mkdir -p "$ROOT"
ROOT=$(cd "$ROOT" && pwd)
IN=$ROOT/inputs/$SEED
OUT=$ROOT/${FAM}__${SEED}
WORK=${WORK_ROOT:-$ROOT/work}/${FAM}__${SEED}

"$PY" "$HARNESS" --data "$DATA" stage --seed "$SEED" --family "$FAM" \
    ${SEED_LABEL:+--label "$SEED_LABEL"} ${SEED_FAA:+--seed-faa "$SEED_FAA"} \
    ${SEED_COORDS:+--seed-coords "$SEED_COORDS"} --out "$IN"

mkdir -p "$OUT"
cp "$IN/seed.json" "$OUT/seed.json"
{
  echo "family=$FAM seed=$SEED date=$(date -Iseconds)"
  echo "repo=$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)" \
       "dirty=$(git -C "$REPO" status --porcelain 2>/dev/null | wc -l)"
  echo "extra_args=$*"
} > "$OUT/run_args.txt"

nextflow -log "$OUT/nextflow.log" run "$REPO/main.nf" \
    --mode pro --query "$IN/query.faa" \
    --home_genome "$IN/home/Apis_mellifera.fna" --home_gff "$IN/home/Apis_mellifera.gff" \
    --home_species 'Apis mellifera' --target_genomes "$IN/targets/*" \
    --keep_intermediate true --outdir "$OUT" -w "$WORK" "$@"

"$PY" "$HARNESS" --data "$DATA" score \
    --run "$OUT" --family "$FAM" --out-prefix "$OUT/family_score" | tee "$OUT/family_score.txt"

# Diagnostic only: a callcheck failure must not fail a scored run.
"$PY" "$HARNESS" --data "$DATA" callcheck --run "$OUT" --targets "$IN/targets" \
    | tee "$OUT/family_score.callcheck.txt" || echo "** callcheck failed (run is scored)"
