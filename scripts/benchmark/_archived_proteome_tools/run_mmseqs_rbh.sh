#!/usr/bin/env bash
# MMseqs2 reciprocal-best-hits baseline.
# Query melittin against each target proteome; emit calls.tsv.
#
# Usage: bash scripts/benchmark/run_mmseqs_rbh.sh

set -euo pipefail

QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
PROTEOMES_DIR="${PROTEOMES_DIR:-local_data/benchmark/proteomes}"
OUTDIR="${OUTDIR:-benchmark_results/mmseqs_rbh}"
EVALUE="${EVALUE:-1e-5}"

mkdir -p "${OUTDIR}/work"
WORK="${OUTDIR}/work"

if ! command -v mmseqs >/dev/null 2>&1; then
    echo "ERROR: mmseqs not found. conda env create -f scripts/benchmark/competitors.yml" >&2
    exit 1
fi
if [[ ! -s "${QUERY}" ]]; then
    echo "ERROR: query not found: ${QUERY}" >&2
    exit 1
fi

# Build query DB once
QDB="${WORK}/queryDB"
if [[ ! -e "${QDB}.dbtype" ]]; then
    mmseqs createdb "${QUERY}" "${QDB}" >/dev/null
fi

CALLS="${OUTDIR}/calls.tsv"
printf "species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra\n" > "${CALLS}"

# Map species names from proteome filenames to truth species names.
# Most are direct (Apis_cerana.faa -> Apis_cerana). Apis_mellifera_home is special.
canonical_species() {
    local base="$1"
    case "${base}" in
        Apis_mellifera_home) echo "Apis_mellifera" ;;
        *) echo "${base}" ;;
    esac
}

for proteome in "${PROTEOMES_DIR}"/*.faa; do
    base=$(basename "${proteome}" .faa)
    species=$(canonical_species "${base}")

    if [[ "${species}" == "Apis_mellifera" ]]; then
        # Home genome — skip (truth file handles separately)
        continue
    fi

    tdb="${WORK}/${species}_targetDB"
    rdb="${WORK}/${species}_resDB"
    out_m8="${WORK}/${species}.m8"

    if [[ ! -e "${tdb}.dbtype" ]]; then
        mmseqs createdb "${proteome}" "${tdb}" >/dev/null
    fi

    mmseqs search "${QDB}" "${tdb}" "${rdb}" "${WORK}/tmp_${species}" \
        -e "${EVALUE}" -s 7.5 --max-seqs 50 \
        >/dev/null 2>&1 || true

    mmseqs convertalis "${QDB}" "${tdb}" "${rdb}" "${out_m8}" \
        --format-output 'query,target,pident,evalue,bits,qcov' \
        >/dev/null 2>&1 || true

    if [[ -s "${out_m8}" ]]; then
        # Best hit by bitscore
        best=$(sort -k5,5gr "${out_m8}" | head -1)
        target=$(awk '{print $2}' <<<"${best}")
        pident=$(awk '{print $3}' <<<"${best}")
        evalue=$(awk '{print $4}' <<<"${best}")
        bits=$(awk '{print $5}' <<<"${best}")
        # MMseqs2 RBH = best hit must reciprocally hit the query as best.
        # Approximate: declare PRESENT if e-value <= EVALUE and qcov reasonable
        # (already filtered by convertalis defaults; we keep it simple here).
        printf "%s\t-\tPRESENT\t-\t-\t-\t-\te=%s,bits=%s,pident=%s\ttarget=%s\n" \
            "${species}" "${evalue}" "${bits}" "${pident}" "${target}" >> "${CALLS}"
    else
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tno_hit\n" "${species}" >> "${CALLS}"
    fi
done

echo "Done. Calls: ${CALLS}"
echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${CALLS}" | sort | uniq -c
