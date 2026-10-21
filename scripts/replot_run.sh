#!/bin/bash
# =============================================================================
#  replot_run.sh — re-draw every figure of a finished SynVoy run with this copy's
#  figure code, without re-running the search. One command per run.
# =============================================================================
#  Usage:
#     scripts/replot_run.sh <run_dir> [plot_synteny.py options]
#
#  <run_dir> is the run's --outdir (the folder that holds plot_inputs_synteny_block_locus_*).
#  The figures of every locus are written to <run_dir>/figures_redrawn_<date>/; the
#  original figures are not touched. Options after <run_dir> go to every locus, e.g.
#     scripts/replot_run.sh results/my_run --home_species "Photinus pyralis"
#
#  Only the drawing is new: the calls, labels and genomes are those of the original
#  run. REDRAWN.txt in the output folder records both versions.
#  Run it with the SynVoy environment active (conda activate synvoy_env). Plain python3
#  also draws every figure, but without ete3 and its NCBI taxonomy file
#  (~/.etetoolkit/taxa.sqlite) the rows follow the GOI gene tree, not the species tree;
#  the console says which tree ordered the rows. No network, except a gene-name lookup
#  at NCBI Gene for a run without a gene-names table.
# =============================================================================
set -euo pipefail

if [ $# -lt 1 ] || [ "${1:0:2}" = "--" ]; then
    sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 2
fi

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[ -d "$1" ] || { echo "ERROR: no such folder: $1" >&2; exit 1; }
RUN="$(cd "$1" && pwd)"
shift
OUT="$RUN/figures_redrawn_$(date +%Y-%m-%d)"

shopt -s nullglob
inputs=("$RUN"/plot_inputs_*/)
shopt -u nullglob
[ ${#inputs[@]} -ge 1 ] || { echo "ERROR: $RUN has no plot_inputs_* folder — is it a SynVoy --outdir?" >&2; exit 1; }

mkdir -p "$OUT"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

failed=0
older=""
gene_tree=""
old_labels=""
for dir in "${inputs[@]}"; do
    dir="${dir%/}"
    locus="$(basename "$dir")"
    in="$dir"
    # Runs from before the rename (SynTerra) or with the old query-BED name: draw from
    # a renamed copy, never change the run itself. -L: a linked file is copied, so the
    # renaming cannot reach the file it links to.
    if grep -qs -m1 'SynTerra' "$dir"/*.gff || { [ -f "$dir/home_gene_location.bed" ] && ! ls "$dir"/locus_*.bed >/dev/null 2>&1; }; then
        in="$TMP/$locus"
        cp -RL "$dir" "$in"
        for g in "$in"/*.gff; do sed -i.bak 's/SynTerra_Parent=/SynVoy_Parent=/g; s/SynTerraRole=/SynVoyRole=/g' "$g" && rm -f "$g.bak"; done
        n="${locus##*_locus_}"
        [ -f "$in/home_gene_location.bed" ] && [ ! -f "$in/locus_${n}.bed" ] && cp "$in/home_gene_location.bed" "$in/locus_${n}.bed"
        echo "[$locus] older run format: drawing from a renamed copy"
        older="$older $locus"
    fi
    # Labels from before the label rework of 2026-10-15 carry no Basis= (the judge's field).
    targets=("$dir"/*.fna.gff "$dir"/*.fa.gff "$dir"/*.fasta.gff "$dir"/*.fna.gz.gff "$dir"/*.fa.gz.gff "$dir"/*.fasta.gz.gff)
    if grep -qs -m1 'Role=goi' "${targets[@]}" && ! grep -qs -m1 'Basis=' "${targets[@]}"; then
        echo "[$locus] its labels come from code before the label rework of 2026-10-15; they are drawn as the run made them"
        old_labels="$old_labels $locus"
    fi
    echo "[$locus] drawing ..."
    if ! "$REPO/scripts/replot.sh" "$in" "$OUT" "$@" > "$OUT/${locus}.log" 2>&1; then
        echo "[$locus] FAILED, see $OUT/${locus}.log" >&2
        failed=$((failed + 1))
    elif grep -q 'GOI gene tree' "$OUT/${locus}.log"; then
        echo "[$locus] rows ordered by the GOI gene tree: no species tree was built (it needs ete3 and ~/.etetoolkit/taxa.sqlite, as in synvoy_env, and two species NCBI knows: pass --home_species)"
        gene_tree="$gene_tree $locus"
    fi
done

code_version="$(cat "$REPO/VERSION" 2>/dev/null || echo unknown)"
code_commit="$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)"
if [ "$code_commit" != unknown ] && ! git -C "$REPO" diff --quiet HEAD -- 2>/dev/null; then
    code_commit="$code_commit (with local changes)"
fi
{
    echo "Figures re-drawn on $(date '+%Y-%m-%d %H:%M') from $RUN"
    echo "Figure code: SynVoy $code_version, commit $code_commit"
    echo "Calls, labels and genomes: those of the original run, as saved in its"
    echo "plot_inputs_* folders; only the drawing is new."
    [ -n "$old_labels" ] && echo "Labels from code before the label rework of 2026-10-15, not today's labels:$old_labels"
    [ -n "$older" ] && echo "Older run format, drawn from a renamed copy:$older"
    [ -n "$gene_tree" ] && echo "Rows ordered by the GOI gene tree, not the species tree:$gene_tree"
    echo "Loci: ${#inputs[@]}, failed: $failed"
} > "$OUT/REDRAWN.txt"

echo "Figures written to $OUT/ (${#inputs[@]} loci, $failed failed)"
[ "$failed" -eq 0 ]
