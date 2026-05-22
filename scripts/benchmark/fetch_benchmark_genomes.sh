#!/usr/bin/env bash
# Fetch the 13 new genomes for the melittin SynVoy benchmark.
# Existing 5 (Colletes, Euglossa, Melipona, Tetragonula, Xylocopa) are already
# in local_data/ground_truth/melettin/ with hand-curated Koludarov 2023 GFFs.
#
# Requires: NCBI Datasets CLI v15+ (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
# Install: conda install -c conda-forge ncbi-datasets-cli   OR
#          curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' && chmod +x datasets
#
# Usage: bash scripts/benchmark/fetch_benchmark_genomes.sh [OUTDIR]
# Default OUTDIR: local_data/benchmark/genomes/

set -euo pipefail

OUTDIR="${1:-local_data/benchmark/genomes}"
mkdir -p "${OUTDIR}"

DATASETS="${DATASETS:-datasets}"
if ! command -v "${DATASETS}" >/dev/null 2>&1; then
    echo "ERROR: 'datasets' CLI not found. Install with:" >&2
    echo "  conda install -c conda-forge ncbi-datasets-cli" >&2
    exit 1
fi

# RefSeq + Ensembl accessions for the 13 species we need to add.
# Format: ACCESSION|TIER|SPECIES_TAG (used to rename output files)
ACCESSIONS=(
    "GCF_003254395.2|home|Apis_mellifera"
    "GCF_029169275.1|tier1|Apis_cerana"
    "GCF_000184785.3|tier1|Apis_florea"
    "GCF_910591885.1|tier1|Bombus_terrestris"
    "GCF_000188095.3|tier1|Bombus_impatiens"
    "GCF_000220905.1|tier2|Megachile_rotundata"
    "GCF_907164935.1|tier2|Osmia_bicornis"
    "GCF_003710045.2|tier2|Nomia_melanderi"
    "GCF_912470025.1|tier3|Vespa_velutina"
    "GCF_001313835.1|tier3|Polistes_canadensis"
    "GCF_016802725.1|tier3|Solenopsis_invicta"
    "GCF_000949405.1|tier3|Vollenhovia_emeryi"
    # v2 additions: ants with paper-validated GR2 (melittin) homologs.
    # Tetramorium and Formica are GenBank-only (no GFF). That's fine for the
    # genome-based competitor tools (tblastn, MMseqs-9.5, miniprot, TOGA) —
    # they don't need a target annotation.
    "GCA_928718305.1|tier3|Tetramorium_bicarinatum"   # 3 GR2 copies, U11/U13-MYRTX-Tb1a
    "GCF_019399895.1|tier3|Cardiocondyla_obscurior"   # 5 Clade_H tandem at GR2
    "GCA_009859135.1|tier3|Formica_selysi"            # 5 GR2 copies (Clade_M, N)
)

TMPDIR="$(mktemp -d)"
trap 'rm -rf "${TMPDIR}"' EXIT

for entry in "${ACCESSIONS[@]}"; do
    IFS='|' read -r acc tier species <<<"${entry}"
    species_dir="${OUTDIR}/${tier}/${species}"
    fa_path="${species_dir}/${species}.fna"
    gff_path="${species_dir}/${species}.gff"

    if [[ -s "${fa_path}" && -s "${gff_path}" ]]; then
        echo "[${species}] already present, skipping"
        continue
    fi

    echo "[${species}] fetching ${acc} (${tier})"
    mkdir -p "${species_dir}"
    zip_path="${TMPDIR}/${acc}.zip"

    "${DATASETS}" download genome accession "${acc}" \
        --include genome,gff3 \
        --filename "${zip_path}" \
        --no-progressbar

    unzip -q -o "${zip_path}" -d "${TMPDIR}/${acc}"

    src_fna=$(find "${TMPDIR}/${acc}/ncbi_dataset/data/${acc}" -maxdepth 1 -name '*.fna' | head -1)
    src_gff=$(find "${TMPDIR}/${acc}/ncbi_dataset/data/${acc}" -maxdepth 1 -name '*.gff' | head -1)

    if [[ -z "${src_fna}" ]]; then
        echo "  ERROR: ${acc} missing genome.fna" >&2
        ls "${TMPDIR}/${acc}/ncbi_dataset/data/${acc}/" >&2 || true
        rmdir "${species_dir}" 2>/dev/null || true
        continue
    fi
    cp "${src_fna}" "${fa_path}"

    # GFF is optional. Genome-based tools (SynVoy, tblastn, MMseqs-9.5, miniprot,
    # TOGA chain-based) work fine without one. Proteome-based tools cannot use
    # this species — that's expected for GenBank-only assemblies.
    if [[ -n "${src_gff}" ]]; then
        cp "${src_gff}" "${gff_path}"
        echo "  -> ${fa_path} (+ GFF)"
    else
        echo "  -> ${fa_path} (FASTA only — no annotation in NCBI Datasets)"
    fi
    rm -rf "${TMPDIR}/${acc}" "${zip_path}"
    echo "  -> ${fa_path}"
done

echo
echo "Done. Genomes in: ${OUTDIR}/{tier1,tier2,tier3,home}/"
echo "Total disk usage:"
du -sh "${OUTDIR}"
