#!/usr/bin/env bash
# Plain tblastn baseline: query melittin protein against each target genome
# (translated 6-frame). This is the simplest, most-cited "did you check?"
# baseline that anyone would do as a first pass. Operates purely on genome
# FASTA — no annotation needed, fair to compare against SynVoy.
#
# Usage: bash scripts/benchmark/run_tblastn.sh

set -uo pipefail

QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
EXISTING_TARGETS="${EXISTING_TARGETS:-local_data/ground_truth/melettin/targets}"
OUTDIR="${OUTDIR:-benchmark_results/tblastn}"
EVALUE="${EVALUE:-1e-3}"
MIN_PIDENT="${MIN_PIDENT:-20}"  # PRESENT requires identity >=20%, anything lower is junk
THREADS="${THREADS:-4}"

mkdir -p "${OUTDIR}/work"
WORK="${OUTDIR}/work"

if ! command -v tblastn >/dev/null 2>&1; then
    echo "ERROR: tblastn not found. Install: mamba install -n synvoy_benchmark -c bioconda blast" >&2
    exit 1
fi
[[ -s "${QUERY}" ]] || { echo "ERROR: query missing: ${QUERY}" >&2; exit 1; }

CALLS="${OUTDIR}/calls.tsv"
printf "species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra\n" > "${CALLS}"

# Helper: collect all (species, fasta_path) pairs from both NCBI tier dirs
# AND the existing Koludarov targets.
declare -A SPECIES_FASTA
for tier_dir in "${GENOMES_DIR}"/{tier1,tier2,tier3}; do
    [[ -d "${tier_dir}" ]] || continue
    for sd in "${tier_dir}"/*/; do
        sp=$(basename "${sd}")
        fa=$(ls "${sd}"*.fna 2>/dev/null | head -1)
        [[ -n "${fa}" ]] && SPECIES_FASTA[$sp]="$fa"
    done
done
# Existing 5 (Koludarov)
for sp in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    fa="${EXISTING_TARGETS}/${sp}.fa"
    [[ -s "${fa}" ]] && SPECIES_FASTA[$sp]="$fa"
done

n_total=${#SPECIES_FASTA[@]}
i=0
for sp in $(echo "${!SPECIES_FASTA[@]}" | tr ' ' '\n' | sort); do
    i=$((i + 1))
    fa="${SPECIES_FASTA[$sp]}"
    db="${WORK}/${sp}_db"
    out_m8="${WORK}/${sp}.m8"
    echo "[${i}/${n_total}] ${sp} ($(basename ${fa}))"

    # Build nucleotide BLAST DB
    if [[ ! -s "${db}.nhr" && ! -s "${db}.nsq" ]]; then
        makeblastdb -in "${fa}" -dbtype nucl -out "${db}" 2>/dev/null >/dev/null \
            || { echo "  makeblastdb failed for ${sp}" >&2; continue; }
    fi

    # tblastn
    if [[ ! -s "${out_m8}" ]]; then
        tblastn -query "${QUERY}" -db "${db}" \
            -evalue "${EVALUE}" \
            -num_threads "${THREADS}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs sstrand" \
            -out "${out_m8}" 2>/dev/null \
            || { echo "  tblastn failed for ${sp}" >&2; continue; }
    fi

    # Pick best hit (max bitscore) over MIN_PIDENT
    if [[ -s "${out_m8}" ]]; then
        best=$(awk -v m="${MIN_PIDENT}" '$3 >= m' "${out_m8}" | sort -k12,12gr | head -1)
        if [[ -n "${best}" ]]; then
            chrom=$(awk '{print $2}' <<<"${best}")
            pident=$(awk '{print $3}' <<<"${best}")
            sstart=$(awk '{print $9}' <<<"${best}")
            send=$(awk '{print $10}' <<<"${best}")
            evalue=$(awk '{print $11}' <<<"${best}")
            bits=$(awk '{print $12}' <<<"${best}")
            qcov=$(awk '{print $13}' <<<"${best}")
            sstrand=$(awk '{print $14}' <<<"${best}")
            # Normalize coords (sstart > send means minus strand)
            if [[ "${sstart}" -gt "${send}" ]]; then
                tmp=${sstart}; sstart=${send}; send=${tmp}
                strand="-"
            else
                strand="+"
            fi
            [[ "${sstrand}" == "minus" ]] && strand="-"
            printf "%s\t-\tPRESENT\t%s\t%s\t%s\t%s\te=%s,bits=%s,pident=%s,qcov=%s\ttblastn_min_pident=%s\n" \
                "${sp}" "${chrom}" "${sstart}" "${send}" "${strand}" \
                "${evalue}" "${bits}" "${pident}" "${qcov}" "${MIN_PIDENT}" >> "${CALLS}"
        else
            # Hits exist but all below MIN_PIDENT
            printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\ttblastn_below_pident_floor=%s\n" "${sp}" "${MIN_PIDENT}" >> "${CALLS}"
        fi
    else
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\ttblastn_no_hit_at_e=%s\n" "${sp}" "${EVALUE}" >> "${CALLS}"
    fi
done

echo
echo "Done. Calls: ${CALLS}"
echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${CALLS}" | sort | uniq -c
