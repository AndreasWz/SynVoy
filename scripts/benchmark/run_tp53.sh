#!/usr/bin/env bash
# Stage and run the TP53 SynVoy benchmark across 20 vertebrate targets.
#
# Produces two paired runs sharing fetched genomes via -resume:
#   A) sensitivity      — default params + wider flanking window. Tests
#                         cross-vertebrate recall over ~500 MY of divergence.
#   B) paralog_discrim  — preset_paralog_discrimination + TP53/63/73 tokens.
#                         Tests whether SynVoy separates the three paralogs
#                         cleanly at loci where multiple are present (elephant
#                         shark, all tetrapods).
#
# Intentionally NOT auto-launched. Review, then run by hand:
#     bash scripts/benchmark/run_tp53.sh fetch       # download genomes only
#     bash scripts/benchmark/run_tp53.sh sensitivity # Run A
#     bash scripts/benchmark/run_tp53.sh paralog     # Run B
#     bash scripts/benchmark/run_tp53.sh score       # rescore both calls.tsv
#                                                      against tp53_orthologs.tsv
#
# Wall-time rough estimate on a 16 GB PC: 8-24 h per run (20 targets, many
# large mammalian genomes). Comment out lungfish if disk is tight (40 GB).

set -euo pipefail

ACTION="${1:-help}"

QUERY="${QUERY:-local_data/queries/P04637_TP53_human.faa}"
HOME_ACC="GCF_000001405.40"   # Homo sapiens GRCh38.p14
GENOMES_ROOT="${GENOMES_ROOT:-local_data/benchmark_tp53}"
NEXTFLOW_ENV="${NEXTFLOW_ENV:-synvoy_env}"

# (snake_case, scientific_name, GCF accession, tier)
TARGETS=(
    "Pan_troglodytes;Pan troglodytes;GCF_028858775.2;tier1"
    "Macaca_mulatta;Macaca mulatta;GCF_049350105.2;tier1"
    "Mus_musculus;Mus musculus;GCF_000001635.27;tier1"
    "Rattus_norvegicus;Rattus norvegicus;GCF_036323735.1;tier1"
    "Heterocephalus_glaber;Heterocephalus glaber;GCF_000247695.1;tier1"
    "Canis_lupus_familiaris;Canis lupus familiaris;GCF_011100685.1;tier1"
    "Bos_taurus;Bos taurus;GCF_002263795.3;tier1"
    "Loxodonta_africana;Loxodonta africana;GCF_030014295.1;tier1"
    "Monodelphis_domestica;Monodelphis domestica;GCF_027887165.1;tier1"
    "Gallus_gallus;Gallus gallus;GCF_016699485.2;tier2"
    "Anolis_carolinensis;Anolis carolinensis;GCF_035594765.1;tier2"
    "Crocodylus_porosus;Crocodylus porosus;GCF_001723895.1;tier2"
    "Xenopus_tropicalis;Xenopus tropicalis;GCF_000004195.4;tier2"
    "Latimeria_chalumnae;Latimeria chalumnae;GCF_037176945.1;tier3"
    # Comment the next line if disk is tight — Protopterus annectens is ~40 Gb.
    "Protopterus_annectens;Protopterus annectens;GCF_019279795.1;tier3"
    "Danio_rerio;Danio rerio;GCF_049306965.1;tier3"
    "Oryzias_latipes;Oryzias latipes;GCF_053564925.1;tier3"
    "Takifugu_rubripes;Takifugu rubripes;GCF_901000725.3;tier3"
    "Callorhinchus_milii;Callorhinchus milii;GCF_018977255.1;tier4"
    "Petromyzon_marinus;Petromyzon marinus;GCF_010993605.1;tier4"
)


need_datasets() {
    if ! mamba run -n "${NEXTFLOW_ENV}" which datasets >/dev/null 2>&1; then
        echo "ERROR: 'datasets' CLI not found in env '${NEXTFLOW_ENV}'." >&2
        exit 1
    fi
}

fetch_one() {
    local snake="$1" sci="$2" acc="$3"
    local outdir="${GENOMES_ROOT}/genomes/${snake}"
    if [[ -s "${outdir}/${snake}.fna" && -s "${outdir}/${snake}.gff" ]]; then
        echo "  [skip] ${snake} already present"
        return 0
    fi
    mkdir -p "${outdir}"
    local zip="${outdir}/${acc}.zip"
    echo "  [fetch] ${snake} ${acc}"
    mamba run -n "${NEXTFLOW_ENV}" datasets download genome accession "${acc}" \
        --include genome,gff3 --filename "${zip}" >/dev/null 2>&1 || {
            echo "    WARN: download failed for ${acc}"
            return 0
        }
    (cd "${outdir}" && unzip -o -q "${zip}" >/dev/null)
    local fna gff
    fna="$(find "${outdir}/ncbi_dataset" -name '*.fna' | head -1)"
    gff="$(find "${outdir}/ncbi_dataset" -name '*.gff' | head -1)"
    if [[ -n "${fna}" ]]; then mv "${fna}" "${outdir}/${snake}.fna"; fi
    if [[ -n "${gff}" ]]; then mv "${gff}" "${outdir}/${snake}.gff"; fi
    rm -rf "${outdir}/ncbi_dataset" "${zip}" "${outdir}/README.md" "${outdir}/md5sum.txt" 2>/dev/null || true
}

fetch_all() {
    need_datasets
    mkdir -p "${GENOMES_ROOT}/genomes"
    echo "Fetching home genome (Homo sapiens)..."
    fetch_one "Homo_sapiens" "Homo sapiens" "${HOME_ACC}"
    echo "Fetching ${#TARGETS[@]} target genomes..."
    for entry in "${TARGETS[@]}"; do
        IFS=';' read -r snake sci acc tier <<< "${entry}"
        fetch_one "${snake}" "${sci}" "${acc}"
    done
    echo "Done. Genomes under ${GENOMES_ROOT}/genomes/"
}

prepare_query() {
    if [[ -s "${QUERY}" ]]; then return 0; fi
    mkdir -p "$(dirname "${QUERY}")"
    cat > "${QUERY}" <<'EOF'
>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
EOF
}

stage_targets() {
    local dir="${GENOMES_ROOT}/staged_targets"
    mkdir -p "${dir}"
    for entry in "${TARGETS[@]}"; do
        IFS=';' read -r snake sci acc tier <<< "${entry}"
        local fna="${GENOMES_ROOT}/genomes/${snake}/${snake}.fna"
        if [[ -s "${fna}" ]]; then
            ln -sfn "$(realpath "${fna}")" "${dir}/${snake}.fna"
        else
            echo "  WARN: ${snake} fna missing — skipping"
        fi
    done
    echo "${dir}"
}

run_nextflow() {
    local outdir="$1"; shift
    local label="$1"; shift
    local staged
    staged="$(stage_targets)"
    local home_fna="${GENOMES_ROOT}/genomes/Homo_sapiens/Homo_sapiens.fna"
    local home_gff="${GENOMES_ROOT}/genomes/Homo_sapiens/Homo_sapiens.gff"
    if [[ ! -s "${home_fna}" || ! -s "${home_gff}" ]]; then
        echo "ERROR: home genome missing at ${home_fna}. Run 'fetch' first." >&2
        exit 1
    fi
    echo "Running SynVoy ${label} (outdir: ${outdir})..."
    mamba run -n "${NEXTFLOW_ENV}" nextflow run main.nf -profile standard --mode pro \
        --query "${QUERY}" \
        --home_genome "${home_fna}" \
        --home_gff "${home_gff}" \
        --target_genomes "${staged}/*" \
        --outdir "${outdir}" \
        "$@" \
        -resume
}

run_sensitivity() {
    prepare_query
    run_nextflow "results/tp53_benchmark_sensitivity" "SENSITIVITY" \
        --n_flanking_genes 7
}

run_paralog() {
    prepare_query
    run_nextflow "results/tp53_benchmark_paralog" "PARALOG_DISCRIM" \
        -profile preset_paralog_discrimination \
        --goi_family_tokens 'TP53,TP63,TP73,TRP53,TRP63,TRP73,P53,P63,P73'
}

parse_calls() {
    # Parse a results dir into benchmark_results/<tool>/calls.tsv from the report's records
    # (scripts/benchmark/synvoy_calls.py, as scripts/benchmark/run_synvoy.sh does).
    local results_dir="$1"
    local calls_path="$2"
    python3 scripts/benchmark/synvoy_calls.py "${results_dir}" "${calls_path}"
}

score_runs() {
    parse_calls "results/tp53_benchmark_sensitivity" \
                "benchmark_results/synvoy_tp53_sensitivity/calls.tsv"
    parse_calls "results/tp53_benchmark_paralog" \
                "benchmark_results/synvoy_tp53_paralog/calls.tsv"
    python3 scripts/benchmark/score_benchmark.py \
        --truth tests/benchmark_truth/tp53_orthologs.tsv \
        --calls-glob 'benchmark_results/synvoy_tp53_*/calls.tsv' \
        --outdir benchmark_results/_tp53_score
    echo
    echo "=== confusion_per_tool.tsv ==="
    cat benchmark_results/_tp53_score/confusion_per_tool.tsv
}

usage() {
    cat <<EOF
Usage: $0 <action>
Actions:
  fetch       Download home + 20 target genomes (~150 GB total without lungfish, ~200 GB with)
  sensitivity Run A — default params + n_flanking=7
  paralog     Run B — preset_paralog_discrimination with TP53/63/73 tokens
  score       Parse both run outputs and score against tp53_orthologs.tsv
  help        This message
EOF
}

case "${ACTION}" in
    fetch)       fetch_all ;;
    sensitivity) run_sensitivity ;;
    paralog)     run_paralog ;;
    score)       score_runs ;;
    help|*)      usage ;;
esac
