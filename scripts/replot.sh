#!/bin/bash
# =============================================================================
#  replot.sh — re-draw the SynVoy figures of a finished run, without re-running
#  the search. Useful for trying figure options.
# =============================================================================
#  Usage:
#     scripts/replot.sh <plot_inputs_dir> <out_dir> [plot_synteny.py options]
#
#  <plot_inputs_dir> is the plot_inputs_synteny_block_locus_<N> folder inside a
#  run's --outdir; it holds everything the plot needs. The figures are written to
#  <out_dir> (created if missing), so the originals are not overwritten.
#
#  Examples:
#     scripts/replot.sh results/my_run/plot_inputs_synteny_block_locus_1 results/my_run/replot
#     scripts/replot.sh results/my_run/plot_inputs_synteny_block_locus_1 /tmp/fig \
#         --home_species "Apis mellifera" --grid_goi_style genomic --no_fragment_variant
#
#  --home_species labels the home row (the pipeline passes it; without it the row
#  reads "Home genome").
#
#  Run it with the SynVoy environment active (conda activate synvoy_env).
#  All options: python3 bin/plot_synteny.py --help
# =============================================================================
set -euo pipefail

if [ $# -lt 2 ]; then
    sed -n '2,23p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 2
fi

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[ -d "$1" ] || { echo "ERROR: no such folder: $1" >&2; exit 1; }
IN="$(cd "$1" && pwd)"
OUT="$2"
shift 2
mkdir -p "$OUT"

shopt -s nullglob
home_beds=("$IN"/synteny_block_*.bed)
query_beds=("$IN"/locus_*.bed)
trees=("$IN"/*tree.nwk)
homology=("$IN"/*.homology.tsv)
regions=("$IN"/*.regions.bed)
target_gffs=()
home_gff=""
for g in "$IN"/*.gff; do
    case "$g" in
        # A target model file is named after its genome FASTA: <genome>.fna.gff
        *.fna.gff|*.fa.gff|*.fasta.gff|*.fna.gz.gff|*.fa.gz.gff|*.fasta.gz.gff) target_gffs+=("$g") ;;
        *) [ -z "$home_gff" ] && home_gff="$g" ;;
    esac
done
shopt -u nullglob

[ ${#home_beds[@]} -ge 1 ] || { echo "ERROR: $IN has no synteny_block_*.bed — is it a plot_inputs_* folder?" >&2; exit 1; }
[ ${#target_gffs[@]} -ge 1 ] || { echo "ERROR: $IN has no <genome>.fna.gff target model files" >&2; exit 1; }

home_bed="${home_beds[0]}"
name="$(basename "$home_bed" .bed)"

args=(--home_bed "$home_bed" --target_gffs "${target_gffs[@]}")
[ ${#query_beds[@]} -ge 1 ] && args+=(--query_bed "${query_beds[0]}")
[ -n "$home_gff" ]          && args+=(--home_gff "$home_gff")
[ ${#homology[@]} -ge 1 ]   && args+=(--homology_tsvs "${homology[@]}")
[ ${#regions[@]} -ge 1 ]    && args+=(--candidate_beds "${regions[@]}")
[ ${#trees[@]} -ge 1 ]      && args+=(--tree "${trees[0]}")
[ -f "$IN/species_mapping.tsv" ] && args+=(--species_map "$IN/species_mapping.tsv")

# Same defaults as the pipeline (modules/plot_synteny.nf); later options win, so
# anything passed on the command line overrides these.
PYTHONHASHSEED=0 python3 "$REPO/bin/plot_synteny.py" \
    "${args[@]}" \
    --hide_goi_absent \
    --gap_threshold 50000 --gap_visual_size 20000 \
    --flank_fallback_bp 1000000 --scale_bar_len 10000 --plot_width 1500 \
    --no_network \
    --output "$OUT/${name}_synteny_plot.html" \
    "$@"

echo "Figures written to $OUT/"
ls "$OUT" | sed 's/^/  /'
