#!/usr/bin/env bash
# Extract protein FASTAs from each benchmark genome + GFF using SynVoy's
# gff_to_faa.py. Produces proteomes used by OrthoFinder, SonicParanoid2,
# and MMseqs2 RBH wrappers.
#
# Usage: bash scripts/benchmark/prep_proteomes.sh

set -euo pipefail

GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
PROTEOMES_DIR="${PROTEOMES_DIR:-local_data/benchmark/proteomes}"
GFF_TO_FAA="${GFF_TO_FAA:-bin/gff_to_faa.py}"
EXISTING_DIR="${EXISTING_DIR:-local_data/ground_truth/melettin}"

if [[ ! -f "${GFF_TO_FAA}" ]]; then
    echo "ERROR: ${GFF_TO_FAA} not found. Run from repo root." >&2
    exit 1
fi

mkdir -p "${PROTEOMES_DIR}"

extract_one() {
    local species="$1"
    local fa="$2"
    local gff="$3"
    local out="${PROTEOMES_DIR}/${species}.faa"

    if [[ -s "${out}" ]]; then
        echo "[${species}] proteome exists, skipping"
        return 0
    fi
    if [[ -z "${fa}" || -z "${gff}" || ! -s "${fa}" || ! -s "${gff}" ]]; then
        echo "[${species}] SKIP: missing fa or gff (fa='${fa}', gff='${gff}')" >&2
        return 0   # do NOT abort the loop — set -e would kill us
    fi

    echo "[${species}] extracting proteome -> ${out}"
    if ! python "${GFF_TO_FAA}" --genome "${fa}" --gff "${gff}" --output "${out}"; then
        echo "  ERROR: gff_to_faa.py failed for ${species}" >&2
        rm -f "${out}"
        return 0
    fi
    local n
    n=$(grep -c '^>' "${out}" || true)
    echo "  -> ${n} proteins"
    if [[ "${n}" -lt 100 ]]; then
        echo "  WARNING: only ${n} proteins — GFF parse may have failed" >&2
    fi
}

# Fetched genomes: local_data/benchmark/genomes/{home,tier1,tier2,tier3}/<species>/
for tier_dir in "${GENOMES_DIR}"/{home,tier1,tier2,tier3}; do
    [[ -d "${tier_dir}" ]] || continue
    for species_dir in "${tier_dir}"/*/; do
        species=$(basename "${species_dir}")
        fa=$(ls "${species_dir}"*.fna 2>/dev/null | head -1)
        gff=$(ls "${species_dir}"*.gff 2>/dev/null | head -1)
        extract_one "${species}" "${fa}" "${gff}"
    done
done

# Existing 5 Koludarov-curated genomes — different layout
for species in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    fa="${EXISTING_DIR}/targets/${species}.fa"
    gff="${EXISTING_DIR}/targets/${species}.gff3"
    extract_one "${species}" "${fa}" "${gff}"
done

# Apis mellifera home (also under existing layout)
extract_one "Apis_mellifera_home" \
    "${EXISTING_DIR}/home/Apis_mellifera.fa" \
    "${EXISTING_DIR}/home/Apis_mellifera.gff"

echo
echo "Done. Proteomes in: ${PROTEOMES_DIR}/"
ls -la "${PROTEOMES_DIR}/" | tail -n +2
