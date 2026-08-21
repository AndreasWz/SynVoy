#!/bin/bash
# =============================================================================
#  run_demo_insect.sh — SynVoy on Drosophila (pro mode, four genomes, ~587 Mb)
# =============================================================================
#  Home:    Drosophila melanogaster   (144 Mb, gold-standard annotation)
#  Targets: Drosophila simulans       (132 Mb,  97 seqs, ~5 Mya)
#           Drosophila yakuba         (148 Mb,  64 seqs, ~10 Mya)
#           Drosophila pseudoobscura  (163 Mb,  71 seqs, ~25 Mya)
#
#  All four are chromosome-level and RefSeq-annotated, so target models carry
#  real gene names and the calls can be checked against an annotation SynVoy
#  did not produce. 13x more sequence than the yeast demo (44 Mb -> 587 Mb).
#
#  Usage:
#     demo/run_demo_insect.sh                  # default query: Defensin
#     demo/run_demo_insect.sh Adh              # the conserved control
#     demo/run_demo_insect.sh Defensin -resume
#
#  Data lives in local_data/demo_insect/ (gitignored). Fetch with
#  demo/fetch_data_insect.sh.
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

QUERY_NAME="${1:-Defensin}"
shift 2>/dev/null

DEMO=local_data/demo_insect
QUERY="$DEMO/query/${QUERY_NAME}.faa"

if [ ! -f "$QUERY" ]; then
    echo "ERROR: no query at $QUERY"
    echo "Available: $(ls "$DEMO/query" 2>/dev/null | tr '\n' ' ')"
    echo "Run demo/fetch_data_insect.sh first if local_data/demo_insect/ is missing."
    exit 1
fi

OUTDIR="results/demo_insect_${QUERY_NAME,,}"

echo "=============================================================="
echo " SynVoy insect demo — query ${QUERY_NAME}, 3 Drosophila targets"
echo " output: ${OUTDIR}"
echo "=============================================================="

./run_synvoy.sh \
    -profile low_mem \
    -c demo/demo_insect.config \
    --mode pro \
    --query          "$QUERY" \
    --home_genome    "$DEMO/home/Drosophila_melanogaster.fna" \
    --home_gff       "$DEMO/home/Drosophila_melanogaster.gff" \
    --home_species   "Drosophila melanogaster" \
    --target_genomes "$DEMO/genomes/*.fna" \
    --target_gffs    "$DEMO/target_gffs/*.gff" \
    --outdir         "$OUTDIR" \
    "$@"
