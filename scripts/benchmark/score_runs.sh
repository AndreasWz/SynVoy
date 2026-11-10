#!/bin/bash
# score_runs.sh RUN_DIR... : rescore finished family-benchmark runs with the harness in this
# working tree (the run's own score files come from whatever harness version produced it), and
# print the key lines. The family is the run's seed.json's, else the run dir's <FAM> of
# <FAM>__<seed>; a dir of neither form is skipped with a message. $EXTRA is passed
# through to `family_recovery.py score` (e.g. EXTRA="--min-confidence HIGH"). A run the scorer
# refuses (exit 2: unknown family, another truth) prints its ERROR line and the script ends with 1.
PY=${PY:-/home/faw/miniforge3/envs/synvoy_env/bin/python3}
HERE=$(cd "$(dirname "$0")" && pwd)
rc=0
for d in "$@"; do
  [ -d "$d" ] || { echo "no such run dir: $d" >&2; continue; }
  b=$(basename "$d"); fam=""
  if [ -s "$d/seed.json" ]; then
    fam=$($PY -c 'import json, sys; print(json.load(open(sys.argv[1])).get("family") or "")' "$d/seed.json")
  fi
  if [ -z "$fam" ]; then
    case "$b" in *__*) fam=${b%%__*} ;; esac
  fi
  [ -n "$fam" ] || { echo "$d: no seed.json with a family and no <FAM>__<seed> name: skipped" >&2; continue; }
  if ! $PY "$HERE/family_recovery.py" score --run "$d" --family "$fam" \
      --out-prefix "$d/family_score.local" $EXTRA > "$d/family_score.local.txt" 2>&1; then
    why=$(grep -m1 '^ERROR' "$d/family_score.local.txt" || tail -n 1 "$d/family_score.local.txt")
    echo "$b: not scored: $why" >&2
    rc=1; continue
  fi
  printf "%-44s " "$b"
  grep -E "^  (locus_recall |address_recall|reachable_genes|locus_precision|call_classes)" \
      "$d/family_score.local.txt" | sed 's/^  //; s/  */ /g' | tr '\n' '|'
  echo
done
exit $rc
