#!/usr/bin/env bash
# GENESPACE on the 18 benchmark genomes.
# STATUS: SKELETON — needs hands-on R debugging on first run.
# GENESPACE is R-based and notoriously sensitive to R version + Bioconductor
# state. Allocate a full day buffer on first attempt.
#
# Usage: bash scripts/benchmark/run_genespace.sh

set -euo pipefail

PROTEOMES_DIR="${PROTEOMES_DIR:-local_data/benchmark/proteomes}"
GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
EXISTING_TARGETS="${EXISTING_TARGETS:-local_data/ground_truth/melettin/targets}"
OUTDIR="${OUTDIR:-benchmark_results/genespace}"
HOME_GFF="${HOME_GFF:-local_data/ground_truth/melettin/home/Apis_mellifera.gff}"

mkdir -p "${OUTDIR}/work/peptide" "${OUTDIR}/work/bed"

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found. conda env create -f scripts/benchmark/competitors.yml" >&2
    exit 1
fi

# GENESPACE expects:
#   <work>/peptide/<species>.fa     (proteomes — primary transcript only)
#   <work>/bed/<species>.bed         (gene coordinates: chr start end gene_id)
#
# We have proteomes from prep_proteomes.sh. We need to derive BEDs from each GFF.

derive_bed_from_gff() {
    local species="$1"
    local gff="$2"
    local out="${OUTDIR}/work/bed/${species}.bed"
    if [[ -s "${out}" ]]; then return 0; fi
    awk -F'\t' '$3 == "gene" {
        match($9, /ID=[^;]+/);
        id = substr($9, RSTART+3, RLENGTH-3);
        print $1"\t"($4-1)"\t"$5"\t"id;
    }' "${gff}" > "${out}"
    if [[ ! -s "${out}" ]]; then
        echo "  WARNING: empty BED for ${species}" >&2
    fi
}

# Stage proteomes
for f in "${PROTEOMES_DIR}"/*.faa; do
    sp=$(basename "${f}" .faa)
    cp -n "${f}" "${OUTDIR}/work/peptide/${sp}.fa"
done

# Stage BEDs from GFFs
# Newly fetched (with NCBI GFFs)
for tier in home tier1 tier2 tier3; do
    for species_dir in "${GENOMES_DIR}/${tier}"/*/; do
        sp=$(basename "${species_dir}")
        gff=$(ls "${species_dir}"*.gff 2>/dev/null | head -1)
        [[ -n "${gff}" ]] && derive_bed_from_gff "${sp}" "${gff}"
    done
done
# Existing 5 (with Koludarov GFFs — different attribute style; may need post-processing)
for sp in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    gff="${EXISTING_TARGETS}/${sp}.gff3"
    [[ -e "${gff}" ]] && derive_bed_from_gff "${sp}" "${gff}"
done

# Run GENESPACE via inline R — TODO: install GENESPACE in the env first
cat > "${OUTDIR}/work/run.R" <<'R'
# TODO: this is a SKELETON. Real run needs:
#   BiocManager::install("GENESPACE")
#   plus DIAMOND, MCScanX, OrthoFinder on PATH.
library(GENESPACE)
wd <- Sys.getenv("OUTDIR_R", "benchmark_results/genespace/work")
gpar <- init_genespace(wd = wd, path2mcscanx = Sys.which("MCScanX"))
out <- run_genespace(gpar, overwrite = FALSE)
# Extract the orthogroup containing melittin (gene-Melt on NC_037641.1 in Apis mellifera)
# and write benchmark_results/genespace/calls.tsv.
# TODO: parse out$pangenome to find the row containing the home-genome melittin gene
# and emit calls.tsv per the contract.
R

echo "Skeleton ready. To complete: install GENESPACE in the conda env and"
echo "fill in the R parsing logic in: ${OUTDIR}/work/run.R"
echo "Then: OUTDIR_R=${OUTDIR}/work Rscript ${OUTDIR}/work/run.R"
