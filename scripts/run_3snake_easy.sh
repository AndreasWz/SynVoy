#!/usr/bin/env bash
set -euo pipefail

# Query can be a local FASTA path or an accession ID resolvable by SynTerra.
# Default is UniProt P60615 (alpha-bungarotoxin; 3FTx).
QUERY="${1:-P60615}"
OUTDIR="${2:-results_3snake_3ftx}"

shift $(( $# > 0 ? 1 : 0 ))
shift $(( $# > 0 ? 1 : 0 ))

./synterra \
  --mode easy \
  --gene "${QUERY}" \
  --home_species "Naja naja" \
  --target_species "Ophiophagus hannah,Bungarus multicinctus" \
  --max_genomes 2 \
  --bad_quality_policy keep \
  --outdir "${OUTDIR}" \
  "$@"
