#!/usr/bin/env bash
# MMseqs2 translated search at sensitivity 9.5 against each target genome.
# This replicates the EXACT methodology of Weitz 2026 Methods section:
#   "MMseqs2 was run with high sensitivity to maximize recovery of divergent
#    toxin homologs (sensitivity = 9.5) and an e-value cutoff of 1e-3 for
#    initial hit detection."
#
# Operates purely on genome FASTA — no annotation needed.
#
# Usage: bash scripts/benchmark/run_mmseqs_genome.sh

set -uo pipefail

QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
EXISTING_TARGETS="${EXISTING_TARGETS:-local_data/ground_truth/melettin/targets}"
OUTDIR="${OUTDIR:-benchmark_results/mmseqs_genome}"
EVALUE="${EVALUE:-1e-3}"
SENSITIVITY="${SENSITIVITY:-9.5}"
MIN_PIDENT="${MIN_PIDENT:-20}"
THREADS="${THREADS:-4}"

mkdir -p "${OUTDIR}/work"
WORK="${OUTDIR}/work"

if ! command -v mmseqs >/dev/null 2>&1; then
    echo "ERROR: mmseqs not found." >&2
    exit 1
fi
[[ -s "${QUERY}" ]] || { echo "ERROR: query missing: ${QUERY}" >&2; exit 1; }

CALLS="${OUTDIR}/calls.tsv"
printf "species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra\n" > "${CALLS}"

declare -A SPECIES_FASTA
for tier_dir in "${GENOMES_DIR}"/{tier1,tier2,tier3}; do
    [[ -d "${tier_dir}" ]] || continue
    for sd in "${tier_dir}"/*/; do
        sp=$(basename "${sd}")
        fa=$(ls "${sd}"*.fna 2>/dev/null | head -1)
        [[ -n "${fa}" ]] && SPECIES_FASTA[$sp]="$fa"
    done
done
for sp in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    fa="${EXISTING_TARGETS}/${sp}.fa"
    [[ -s "${fa}" ]] && SPECIES_FASTA[$sp]="$fa"
done

n_total=${#SPECIES_FASTA[@]}
i=0
for sp in $(echo "${!SPECIES_FASTA[@]}" | tr ' ' '\n' | sort); do
    i=$((i + 1))
    fa="${SPECIES_FASTA[$sp]}"
    out_m8="${WORK}/${sp}.m8"
    echo "[${i}/${n_total}] ${sp} ($(basename ${fa}))"

    if [[ ! -s "${out_m8}" ]]; then
        # mmseqs easy-search with translated nucleotide search-type 2:
        #   query is protein, target is nucleotide, MMseqs translates target.
        mmseqs easy-search "${QUERY}" "${fa}" "${out_m8}" "${WORK}/tmp_${sp}" \
            --search-type 2 \
            -s "${SENSITIVITY}" \
            -e "${EVALUE}" \
            --threads "${THREADS}" \
            --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov' \
            >/dev/null 2>&1 \
            || { echo "  mmseqs search failed for ${sp}" >&2; continue; }
        rm -rf "${WORK}/tmp_${sp}"
    fi

    if [[ -s "${out_m8}" ]]; then
        best=$(awk -v m="${MIN_PIDENT}" '$3 >= m' "${out_m8}" | sort -k12,12gr | head -1)
        if [[ -n "${best}" ]]; then
            chrom=$(awk '{print $2}' <<<"${best}")
            pident=$(awk '{print $3}' <<<"${best}")
            tstart=$(awk '{print $9}' <<<"${best}")
            tend=$(awk '{print $10}' <<<"${best}")
            evalue=$(awk '{print $11}' <<<"${best}")
            bits=$(awk '{print $12}' <<<"${best}")
            qcov=$(awk '{print $13}' <<<"${best}")
            if [[ "${tstart}" -gt "${tend}" ]]; then
                tmp=${tstart}; tstart=${tend}; tend=${tmp}
                strand="-"
            else
                strand="+"
            fi
            printf "%s\t-\tPRESENT\t%s\t%s\t%s\t%s\te=%s,bits=%s,pident=%s,qcov=%s\tsens=%s\n" \
                "${sp}" "${chrom}" "${tstart}" "${tend}" "${strand}" \
                "${evalue}" "${bits}" "${pident}" "${qcov}" "${SENSITIVITY}" >> "${CALLS}"
        else
            printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tmmseqs_below_pident_floor=%s\n" "${sp}" "${MIN_PIDENT}" >> "${CALLS}"
        fi
    else
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tmmseqs_no_hit_at_sens%s_e%s\n" "${sp}" "${SENSITIVITY}" "${EVALUE}" >> "${CALLS}"
    fi
done

echo
echo "Done. Calls: ${CALLS}"
echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${CALLS}" | sort | uniq -c
