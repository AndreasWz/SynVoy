#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEMO_FILE="${ROOT_DIR}/docs/DEMO_LOG_RICH.txt"

if [[ ! -f "${DEMO_FILE}" ]]; then
  echo "Demo log file not found: ${DEMO_FILE}" >&2
  exit 1
fi

mode="${1:-instant}"

case "${mode}" in
  instant)
    cat "${DEMO_FILE}"
    ;;
  live)
    frames=(
      "| 00:00:02  P1 Gene Localization  |  LOCATE_GENE  -  Locating GOI in home genome"
      "/ 00:00:07  P1 Gene Localization  |  ANNOTATE_GOI  -  Annotating GOI exons"
      "- 00:00:15  P1 Gene Localization  |  BORROW_ANNOTATIONS  -  Evaluating annotation borrowing"
      "\\ 00:00:24  P2 Iterative Search  |  ITERATIVE_SEARCH [2/12]  -  Processing GCA_046270085.1.fna"
      "| 00:01:02  P2 Iterative Search  |  ITERATIVE_SEARCH [7/12]  -  Processing GCF_900518725.1.fna"
      "/ 00:01:39  P3 Region Clustering  |  CLUSTER_REGIONS [10/12]  -  Building loci clusters"
      "- 00:01:56  P4 Trees and Plots  |  COMPUTE_TREE  -  Inferring GOI phylogeny"
      "\\ 00:02:06  P4 Trees and Plots  |  PLOT_SYNTENY  -  Rendering interactive figure"
      "| 00:02:12  P5 Reporting  |  GENERATE_REPORT  -  Writing synterra_report.json"
    )
    for frame in "${frames[@]}"; do
      printf '\r%-140s' "${frame}"
      sleep 0.35
    done
    printf '\r%-140s\n' "[OK] Pipeline finished successfully"
    printf 'Duration: 2m 12s\n'
    printf 'Tasks Completed: 26\n'
    printf 'Raw log: .synterra_logs/run_YYYYmmdd_HHMMSS.log\n'
    ;;
  *)
    echo "Usage: $0 [instant|live]" >&2
    exit 1
    ;;
esac
