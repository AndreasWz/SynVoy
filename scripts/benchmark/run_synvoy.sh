#!/usr/bin/env bash
# Run SynVoy in Pro Mode on the 18 benchmark genomes and emit calls.tsv.
# Reuses the existing ground-truth Apis mellifera home + the 5 Koludarov
# target genomes plus the 13 newly-fetched targets.
#
# Usage: bash scripts/benchmark/run_synvoy.sh

set -euo pipefail

QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
# Use the freshly-fetched RefSeq home (the original symlinks under
# local_data/ground_truth/melettin/home/ point to a deleted local_runs path).
HOME_FA="${HOME_FA:-local_data/benchmark/genomes/home/Apis_mellifera/Apis_mellifera.fna}"
HOME_GFF="${HOME_GFF:-local_data/benchmark/genomes/home/Apis_mellifera/Apis_mellifera.gff}"
GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
EXISTING_TARGETS="${EXISTING_TARGETS:-local_data/ground_truth/melettin/targets}"
OUTDIR="${OUTDIR:-benchmark_results/synvoy}"
RESULTS_DIR="${RESULTS_DIR:-results/melittin_benchmark}"

mkdir -p "${OUTDIR}"

# Stage all 17 target genomes (excluding Apis mellifera home) into a single dir
TARGETS_STAGED="${OUTDIR}/staged_targets"
mkdir -p "${TARGETS_STAGED}"

# Existing 5
for sp in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    [[ -e "${EXISTING_TARGETS}/${sp}.fa" ]] && cp -n "${EXISTING_TARGETS}/${sp}.fa" "${TARGETS_STAGED}/"
done

# Newly fetched 12 (skip home)
for tier in tier1 tier2 tier3; do
    for species_dir in "${GENOMES_DIR}/${tier}"/*/; do
        sp=$(basename "${species_dir}")
        fa=$(ls "${species_dir}"*.fna 2>/dev/null | head -1)
        if [[ -n "${fa}" ]]; then
            cp -n "${fa}" "${TARGETS_STAGED}/${sp}.fna"
        fi
    done
done

n_targets=$(ls "${TARGETS_STAGED}/" | wc -l)
echo "Staged ${n_targets} target genomes"

if [[ "${n_targets}" -lt 17 ]]; then
    echo "WARNING: fewer than 17 targets staged. Continuing anyway." >&2
fi

# Run SynVoy. Nextflow needs Java; both live in synvoy_env. Calling the
# nextflow binary directly inherits the calling shell's PATH (which usually
# lacks java). `mamba run -n synvoy_env` properly activates the env so java
# is on PATH for nextflow's startup.
NEXTFLOW_ENV="${NEXTFLOW_ENV:-synvoy_env}"
if ! mamba run -n "${NEXTFLOW_ENV}" which nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow not found in conda env '${NEXTFLOW_ENV}'." >&2
    echo "  install: mamba install -n ${NEXTFLOW_ENV} -c bioconda nextflow openjdk" >&2
    echo "  or set NEXTFLOW_ENV to a different env name." >&2
    exit 1
fi
echo "Using nextflow from conda env: ${NEXTFLOW_ENV}"

echo "Running SynVoy (Pro Mode, this can take 6-12 hours)..."
mamba run -n "${NEXTFLOW_ENV}" nextflow run main.nf -profile standard --mode pro \
    --query "${QUERY}" \
    --home_genome "${HOME_FA}" \
    --home_gff "${HOME_GFF}" \
    --target_genomes "${TARGETS_STAGED}/*" \
    --n_flanking_genes 5 \
    --outdir "${RESULTS_DIR}" \
    -resume 2>&1 | tail -40 || true

# Parse the calls from synvoy_report.json's records (scripts/benchmark/synvoy_calls.py: PRESENT
# = a counted HIGH or MEDIUM record, AMBIGUOUS = a record that needs curation; regions and LOW
# models set no status).
python3 scripts/benchmark/synvoy_calls.py "${RESULTS_DIR}" "${OUTDIR}/calls.tsv"

echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${OUTDIR}/calls.tsv" | sort | uniq -c
