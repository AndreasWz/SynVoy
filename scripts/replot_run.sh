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
#  Without --home_species the home species is read from the run (its task logs, or the
#  file name of its home GFF); the console says which one is used. At the end the
#  console names the figure for a slide of every locus (<locus>_anchor_grid_print.svg).
#
#  Only the drawing is new: the calls, labels and genomes are those of the original
#  run. REDRAWN.txt in the output folder records both versions.
#  Run it with the SynVoy environment active (conda activate synvoy_env). Plain python3
#  also draws every figure, but without ete3 and its NCBI taxonomy file
#  (~/.etetoolkit/taxa.sqlite) the rows follow the GOI gene tree, not the species tree;
#  the console says so before drawing, with the command that gets the file. No network,
#  except a gene-name lookup at NCBI Gene for a run without a gene-names table.
# =============================================================================
set -euo pipefail

if [ $# -lt 1 ] || [ "${1:0:2}" = "--" ]; then
    sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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

if ! mkdir -p "$OUT" 2>/dev/null || [ ! -w "$OUT" ]; then
    echo "ERROR: cannot write into $RUN: the folder is read-only (for example unpacked from a zip archive)." >&2
    echo "Make it writable, then run this command again:  chmod -R u+w \"$RUN\"" >&2
    exit 1
fi
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# The home species labels the home row and is a leaf of the species tree (with one target
# genome, the only other one). Not given: read from the run. Its task logs hold the value the
# pipeline used (logs/PLOT_SYNTENY*/script.sh, logs/PHYLO_SORT*/script.sh), but a work folder
# shared by several runs leaves the logs of all of them, so only a "Genus species" name counts,
# several are told apart by the file name of the home GFF (pro mode: Genus_species.gff), and
# that name alone is used when the run kept no logs.
genus_species() { sed -nE 's/^([A-Z][a-z]+ [a-z]+).*/\1/p'; }
same_species() { [ "$(printf '%s\n' "$1" | tr '_' ' ' | awk '{print tolower($1" "$2)}')" = "$(printf '%s\n' "$2" | tr '_' ' ' | awk '{print tolower($1" "$2)}')" ]; }
logged="$( { grep -hos -- '--home_species "[^"]*"' "$RUN"/logs/PLOT_SYNTENY*/script.sh || true
             grep -hos '^HOME_QUERY="[^"]*"' "$RUN"/logs/PHYLO_SORT*/script.sh || true; } \
    | sed -E 's/^(--home_species |HOME_QUERY=)"(.*)"$/\2/; s#^.*/([^/]+\.(fna|fa|fasta)(\.gz)?)$#\1#; s/\.(fna|fa|fasta)(\.gz)?$//; s/_/ /g; s/^ +//; s/ +$//' \
    | grep -E '^[A-Z][a-z]+ [a-z]+' | sort -u || true)"
gff_species="$(for d in "${inputs[@]}"; do
        for g in "$d"*.gff; do
            case "$g" in *.fna.gff|*.fa.gff|*.fasta.gff|*.fna.gz.gff|*.fa.gz.gff|*.fasta.gz.gff) continue ;; esac
            [ -f "$g" ] && basename "$g" .gff
            break
        done
    done | tr '_' ' ' | genus_species | sort -u || true)"
n_logged=$(printf '%s' "$logged" | grep -c . || true)
species="" species_from=""
if [ "$n_logged" -eq 1 ]; then
    species="$logged" species_from="its task logs"
elif [ "$n_logged" -gt 1 ]; then
    picked="$(printf '%s\n' "$logged" | while IFS= read -r v; do
        printf '%s\n' "$gff_species" | grep -qxF "$(printf '%s\n' "$v" | genus_species)" && printf '%s\n' "$v"; done || true)"
    [ "$(printf '%s' "$picked" | grep -c . || true)" -eq 1 ] && species="$picked" species_from="its task logs and the file name of its home GFF"
elif [ "$(printf '%s' "$gff_species" | grep -c . || true)" -eq 1 ]; then
    species="$gff_species" species_from="the file name of its home GFF"
fi
given="" given_set=0 prev=""
for a in "$@"; do
    [ "$prev" = "--home_species" ] && given="$a" given_set=1
    case "$a" in --home_species=*) given="${a#--home_species=}" given_set=1 ;; esac
    prev="$a"
done
if [ "$given_set" = 1 ]; then
    echo "Home species: $given (as given)"
    if [ -n "$species" ] && ! same_species "$given" "$species"; then
        echo "WARNING: the run names its home species $species ($species_from), not \"$given\"; the home row will read \"$given\". Leave --home_species out to use $species."
    fi
    species="$given" species_from="as given"
elif [ -n "$species" ]; then
    echo "Home species: $species (read from the run: $species_from; to use another, pass --home_species \"Genus species\")"
    set -- --home_species "$species" "$@"
elif [ "$n_logged" -gt 1 ]; then
    echo "Home species: not given, and the run's task logs name several ($(printf '%s' "$logged" | paste -sd ';' - | sed 's/;/, /g')): the home row reads \"Home genome\". Pass --home_species \"Genus species\"."
else
    echo "Home species: not given and not found in the run: the home row reads \"Home genome\", and a run with one target genome gets no species tree. Pass --home_species \"Genus species\"."
fi

# The species tree needs ete3 (imported as plot_synteny.py imports it) and the NCBI taxonomy
# file; without them every locus falls back to the GOI gene tree, so say it once, up front.
tree_missing=""
if ! python3 -c '
import sys, types
try:
    import cgi
except ModuleNotFoundError:
    m = types.ModuleType("cgi"); m.escape = lambda s, quote=True: s; sys.modules["cgi"] = m
import ete3' >/dev/null 2>&1; then
    tree_missing="no ete3"
    echo "NOTE: $(command -v python3 || echo python3) cannot import ete3, so no species tree can be built and the rows follow the GOI gene tree. Activate the SynVoy environment first: conda activate synvoy_env"
elif [ ! -e "${HOME:-}/.etetoolkit/taxa.sqlite" ]; then
    tree_missing="no taxonomy file"
    echo "NOTE: the NCBI taxonomy file ~/.etetoolkit/taxa.sqlite is missing, so no species tree can be built and the rows follow the GOI gene tree. Get it once (2-3 minutes; 760 MB on disk, 3 GB of memory while it is built), then run this command again:"
    echo "    (cd ~ && python3 -c 'from ete3 import NCBITaxa; NCBITaxa()' && rm -f taxdump.tar.gz)"
fi

failed=0
older=""
gene_tree=""
old_labels=""
slides=""
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
        chmod -R u+w "$in"     # a copy of a read-only run is read-only too
        for g in "$in"/*.gff; do sed -i.bak 's/SynTerra_Parent=/SynVoy_Parent=/g; s/SynTerraRole=/SynVoyRole=/g' "$g" && rm -f "$g.bak"; done
        n="${locus##*_locus_}"
        [ -f "$in/home_gene_location.bed" ] && [ ! -f "$in/locus_${n}.bed" ] && cp "$in/home_gene_location.bed" "$in/locus_${n}.bed"
        echo "[$locus] older run format: drawing from a renamed copy"
        older="$older $locus"
    fi
    # Labels from before the label rework of 2026-10-15 carry no Basis= (the judge's field).
    targets=("$dir"/*.fna.gff "$dir"/*.fa.gff "$dir"/*.fasta.gff "$dir"/*.fna.gz.gff "$dir"/*.fa.gz.gff "$dir"/*.fasta.gz.gff)
    old=0
    if grep -qs -m1 'Role=goi' "${targets[@]}" && ! grep -qs -m1 'Basis=' "${targets[@]}"; then
        old_labels="$old_labels $locus"
        old=1
    fi
    echo "[$locus] drawing ..."
    if ! "$REPO/scripts/replot.sh" "$in" "$OUT" "$@" > "$OUT/${locus}.log" 2>&1; then
        echo "[$locus] FAILED, see $OUT/${locus}.log" >&2
        failed=$((failed + 1))
    elif grep -q 'GOI gene tree' "$OUT/${locus}.log"; then
        if [ -n "$tree_missing" ]; then why="$tree_missing, see the NOTE above"
        elif [ -z "$species" ]; then why="fewer than two of its species are in NCBI's taxonomy; --home_species \"Genus species\" adds the home species"
        else why="fewer than two of its species are in NCBI's taxonomy"; fi
        echo "[$locus] rows ordered by the GOI gene tree, not the species tree ($why)"
        gene_tree="$gene_tree $locus"
    fi
    # The figure for a slide: the anchor grid at print size, named after the folder's home BED.
    bed="$(ls "$in"/synteny_block_*.bed 2>/dev/null | head -1 || true)"
    svg="$OUT/$(basename "${bed:-none}" .bed)_anchor_grid_print.svg"
    if [ -n "$bed" ] && [ -f "$svg" ]; then line="  $svg"
    else line="  [$locus] no print grid was drawn, see $OUT/${locus}.log"; fi
    slides="$slides${locus##*_locus_}"$'\t'"$line"$'\t'"$old"$'\n'
done

code_version="$(cat "$REPO/VERSION" 2>/dev/null || echo unknown)"
# Only this copy's own git history counts: a copy unpacked from an archive into another
# repository would otherwise report that repository's commit.
top="$(git -C "$REPO" rev-parse --show-toplevel 2>/dev/null || true)"
if [ -n "$top" ] && [ "$(cd "$top" && pwd -P)" = "$(cd "$REPO" && pwd -P)" ]; then
    code_commit="$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    if [ "$code_commit" != unknown ] && ! git -C "$REPO" diff --quiet HEAD -- 2>/dev/null; then
        code_commit="$code_commit (with local changes)"
    fi
else
    code_commit="unknown (not a git checkout: $(basename "$REPO"))"
fi
{
    echo "Figures re-drawn on $(date '+%Y-%m-%d %H:%M') from $RUN"
    echo "Figure code: SynVoy $code_version, commit $code_commit"
    echo "Calls, labels and genomes: those of the original run, as saved in its"
    echo "plot_inputs_* folders; only the drawing is new."
    if [ -n "$species" ]; then echo "Home species: $species ($species_from)"
    else echo "Home species: none (the home row reads \"Home genome\")"; fi
    [ -n "$old_labels" ] && echo "Labels from code before the label rework of 2026-10-15, not today's labels:$old_labels"
    [ -n "$older" ] && echo "Older run format, drawn from a renamed copy:$older"
    [ -n "$gene_tree" ] && echo "Rows ordered by the GOI gene tree, not the species tree:$gene_tree"
    echo "Loci: ${#inputs[@]}, failed: $failed"
} > "$OUT/REDRAWN.txt"

echo "Figures written to $OUT/ (${#inputs[@]} loci, $failed failed)"
echo "For a slide, the anchor grid at print size of each locus:"
check="check numbered copies before counting them"
n_old=$(echo $old_labels | wc -w | tr -d " ")
printf '%s' "$slides" | sort -t $'\t' -k1,1n | while IFS=$'\t' read -r _ line old; do
    if [ "$old" = 1 ] && [ "$n_old" -lt ${#inputs[@]} ]; then echo "$line  (labels from before 15 Oct 2026: $check)"
    else echo "$line"; fi
done
if [ "$n_old" -gt 0 ] && [ "$n_old" -eq ${#inputs[@]} ]; then
    echo "The labels of every locus come from code before 15 Oct 2026: $check."
fi
[ "$failed" -eq 0 ]
