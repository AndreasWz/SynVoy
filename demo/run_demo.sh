#!/bin/bash
# =============================================================================
#  run_demo.sh — the SynVoy live demo (pro mode, four budding-yeast genomes)
# =============================================================================
#  Home:    Saccharomyces cerevisiae S288C  (12.1 Mb, 16 chromosomes, annotated)
#  Targets: Naumovozyma castellii           (11.2 Mb, post-WGD, closest)
#           Kluyveromyces lactis            (10.7 Mb, pre-WGD)
#           Lachancea thermotolerans        (10.4 Mb, pre-WGD, farthest)
#
#  All four are chromosome-level and annotated, so the synteny signal is clean
#  and target models carry real gene names. Total input: 76 MB.
#
#  Usage:
#     demo/run_demo.sh                 # default query: STE2
#     demo/run_demo.sh PGK1            # the conserved control
#     demo/run_demo.sh STE2 -resume    # anything extra is passed to Nextflow
#
#  Data lives in local_data/demo/ (gitignored). Re-fetch it with demo/fetch_data.sh.
# =============================================================================
set -o pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

# A VS Code-activated .venv shadows the conda env inside Nextflow tasks and makes
# parasail (Smith-Waterman) disappear — the launcher fails loud on this. Strip it
# here so the demo runs the same way from a terminal or from the IDE.
if [ -n "${VIRTUAL_ENV:-}" ]; then
    echo "[demo] dropping active virtualenv $VIRTUAL_ENV from the environment"
    PATH="$(echo "$PATH" | tr ':' '\n' | grep -v "^${VIRTUAL_ENV}/bin$" | paste -sd:)"
    unset VIRTUAL_ENV
    export PATH
fi

QUERY_NAME="${1:-STE2}"
shift 2>/dev/null

DEMO=local_data/demo
QUERY="$DEMO/query/${QUERY_NAME}.faa"

if [ ! -f "$QUERY" ]; then
    echo "ERROR: no query at $QUERY"
    echo "Available: $(ls "$DEMO/query" 2>/dev/null | tr '\n' ' ')"
    echo "Run demo/fetch_data.sh first if local_data/demo/ is missing."
    exit 1
fi

OUTDIR="results/demo_${QUERY_NAME,,}"

echo "=============================================================="
echo " SynVoy demo — query ${QUERY_NAME}, 3 target yeast genomes"
echo " output: ${OUTDIR}"
echo "=============================================================="

./run_synvoy.sh \
    -profile low_mem \
    -c demo/demo.config \
    --mode pro \
    --query          "$QUERY" \
    --home_genome    "$DEMO/home/Saccharomyces_cerevisiae.fna" \
    --home_gff       "$DEMO/home/Saccharomyces_cerevisiae.gff" \
    --home_species   "Saccharomyces cerevisiae" \
    --target_genomes "$DEMO/genomes/*.fna" \
    --target_gffs    "$DEMO/target_gffs/*.gff" \
    --outdir         "$OUTDIR" \
    "$@"
