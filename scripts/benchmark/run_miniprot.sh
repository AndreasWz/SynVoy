#!/usr/bin/env bash
# miniprot: splice-aware protein-to-genome aligner.
# Built for exactly this problem — aligns a query protein across exon/intron
# boundaries against an unannotated genome. Modern standard for protein->genome.
#
# Usage: bash scripts/benchmark/run_miniprot.sh

set -uo pipefail

QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
EXISTING_TARGETS="${EXISTING_TARGETS:-local_data/ground_truth/melettin/targets}"
OUTDIR="${OUTDIR:-benchmark_results/miniprot}"
THREADS="${THREADS:-4}"
MIN_IDENTITY="${MIN_IDENTITY:-0.20}"  # 20% identity, fraction (miniprot uses fractions)

mkdir -p "${OUTDIR}/work"
WORK="${OUTDIR}/work"

# Try both env locations: synvoy_env (where SynVoy uses miniprot) and synvoy_benchmark
MINIPROT="${MINIPROT:-}"
if [[ -z "${MINIPROT}" ]]; then
    if command -v miniprot >/dev/null 2>&1; then
        MINIPROT=$(command -v miniprot)
    elif [[ -x /home/faw/miniforge3/envs/synvoy_env/bin/miniprot ]]; then
        MINIPROT=/home/faw/miniforge3/envs/synvoy_env/bin/miniprot
    else
        echo "ERROR: miniprot not found. Install: mamba install -n synvoy_benchmark -c bioconda miniprot" >&2
        exit 1
    fi
fi
[[ -s "${QUERY}" ]] || { echo "ERROR: query missing: ${QUERY}" >&2; exit 1; }
echo "Using miniprot: ${MINIPROT}"

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
    gff_out="${WORK}/${sp}.gff"
    echo "[${i}/${n_total}] ${sp} ($(basename ${fa}))"

    if [[ ! -s "${gff_out}" ]]; then
        # miniprot outputs paf (default) or gff with --gff
        # -t threads, --gff for GFF3, no -k since we want best alignment
        "${MINIPROT}" -t "${THREADS}" --gff "${fa}" "${QUERY}" > "${gff_out}" 2>/dev/null \
            || { echo "  miniprot failed for ${sp}" >&2; continue; }
    fi

    # Parse the best mRNA from the GFF
    # miniprot GFF has lines like:
    #   ##PAF gene-Melt 70 0 70 + NC_037641.1 ... AS:i:280 ms:i:280 ...
    #   chrom miniprot mRNA start end identity strand . ID=MP000001;Identity=0.985;Positive=0.992;...
    if [[ -s "${gff_out}" ]]; then
        best=$(awk -F'\t' '$3=="mRNA"' "${gff_out}" | sort -k6,6gr | head -1)
        if [[ -n "${best}" ]]; then
            chrom=$(awk -F'\t' '{print $1}' <<<"${best}")
            start=$(awk -F'\t' '{print $4}' <<<"${best}")
            end=$(awk -F'\t' '{print $5}' <<<"${best}")
            score=$(awk -F'\t' '{print $6}' <<<"${best}")
            strand=$(awk -F'\t' '{print $7}' <<<"${best}")
            attrs=$(awk -F'\t' '{print $9}' <<<"${best}")
            ident=$(grep -oE 'Identity=[0-9.]+' <<<"${attrs}" | head -1 | sed 's/Identity=//')
            posit=$(grep -oE 'Positive=[0-9.]+' <<<"${attrs}" | head -1 | sed 's/Positive=//')
            ident="${ident:-0.0}"
            # Use awk for float comparison (bash can't do floats natively)
            passes=$(awk -v i="${ident}" -v t="${MIN_IDENTITY}" 'BEGIN{print (i+0 >= t+0) ? 1 : 0}')
            if [[ "${passes}" == "1" ]]; then
                printf "%s\t-\tPRESENT\t%s\t%s\t%s\t%s\tident=%s,positive=%s,score=%s\tminiprot\n" \
                    "${sp}" "${chrom}" "${start}" "${end}" "${strand}" \
                    "${ident}" "${posit}" "${score}" >> "${CALLS}"
            else
                printf "%s\t-\tABSENT\t-\t-\t-\t-\tident=%s\tminiprot_below_min_ident=%s\n" \
                    "${sp}" "${ident}" "${MIN_IDENTITY}" >> "${CALLS}"
            fi
        else
            printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tminiprot_no_mRNA\n" "${sp}" >> "${CALLS}"
        fi
    else
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tminiprot_no_output\n" "${sp}" >> "${CALLS}"
    fi
done

echo
echo "Done. Calls: ${CALLS}"
echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${CALLS}" | sort | uniq -c
