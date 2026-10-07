#!/bin/bash
# score_runs.sh RUN_DIR... : rescore finished family-benchmark runs with the harness in this
# working tree (the run's own score files come from whatever harness version produced it), and
# print the key lines. Family name is taken from the run dir: <FAM>__<seed>. $EXTRA is passed
# through to `family_recovery.py score` (e.g. EXTRA="--label VA_CAP_loc").
PY=${PY:-/home/faw/miniforge3/envs/synvoy_env/bin/python3}
HERE=$(cd "$(dirname "$0")" && pwd)
for d in "$@"; do
  [ -d "$d" ] || { echo "no such run dir: $d" >&2; continue; }
  b=$(basename "$d"); fam=${b%%__*}
  $PY "$HERE/family_recovery.py" score --run "$d" --family "$fam" \
      --out-prefix "$d/family_score.local" $EXTRA > "$d/family_score.local.txt" 2>&1
  printf "%-44s " "$b"
  grep -E "^  (locus_recall |address_recall|reachable_genes|locus_precision|call_classes)" \
      "$d/family_score.local.txt" | sed 's/^  //; s/  */ /g' | tr '\n' '|'
  echo
done
