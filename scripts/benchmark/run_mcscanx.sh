#!/usr/bin/env bash
# MCScanX: collinear-block detection between Apis mellifera and each target.
# This is the legitimate synteny-based competitor for SynVoy in the
# annotation-dependent regime. The point of including it (per user) is to
# show *exactly* what happens when annotation is incomplete: targets that
# have no NCBI GFF (Formica_selysi, Tetramorium_bicarinatum) cannot be
# scored at all by MCScanX, while SynVoy still finds melittin in them via
# the genome-anchored flanking-gene search.
#
# Method per target species:
#   1. Require both .faa (built by prep_proteomes.sh) and .gff/.gff3.
#      If either is missing → ABSENT with mcscanx_no_annotation note.
#   2. DIAMOND blastp: Apis ∪ target  vs  Apis ∪ target  (-evalue 1e-10,
#      tabular outfmt 6, top 5 hits/query).
#   3. Build MCScanX simplified GFF (chr_id, gene_id, start, end) with a
#      2-letter species prefix on chrom names so MCScanX can tell sides apart.
#   4. Run MCScanX → prefix.collinearity, prefix.html (we don't need html).
#   5. Parse .collinearity: PRESENT if gene-Melt (Apis) appears in any
#      inter-species block; ABSENT otherwise.
#
# Outputs:
#   benchmark_results/mcscanx/calls.tsv
#   benchmark_results/mcscanx/work/<species>/Apis_vs_<species>.collinearity
#
# Usage:
#   mamba activate synvoy_benchmark
#   bash scripts/benchmark/prep_proteomes.sh    # if not done yet
#   bash scripts/benchmark/run_mcscanx.sh

set -uo pipefail

REPO="${REPO:-/home/faw/dev/projects/SynVoy}"
PROTEOMES_DIR="${PROTEOMES_DIR:-${REPO}/local_data/benchmark/proteomes}"
GENOMES_DIR="${GENOMES_DIR:-${REPO}/local_data/benchmark/genomes}"
EXISTING_DIR="${EXISTING_DIR:-${REPO}/local_data/ground_truth/melettin/targets}"
HOME_GFF="${HOME_GFF:-${GENOMES_DIR}/home/Apis_mellifera/Apis_mellifera.gff}"
APIS_FAA="${APIS_FAA:-${PROTEOMES_DIR}/Apis_mellifera.faa}"
OUTDIR="${OUTDIR:-${REPO}/benchmark_results/mcscanx}"
WORK="${OUTDIR}/work"
THREADS="${THREADS:-4}"
APIS_PREFIX="am"  # 2-letter species prefix used in simplified GFF

mkdir -p "${WORK}"

# ── pre-flight ──────────────────────────────────────────────────────────
if ! command -v diamond >/dev/null 2>&1; then
    echo "ERROR: diamond not on PATH. Install: mamba install -n synvoy_benchmark -c bioconda diamond" >&2
    exit 1
fi
if ! command -v MCScanX >/dev/null 2>&1; then
    echo "ERROR: MCScanX not on PATH. Install: mamba install -n synvoy_benchmark -c bioconda mcscanx" >&2
    exit 1
fi
[[ -s "${APIS_FAA}" ]] || { echo "ERROR: Apis proteome missing: ${APIS_FAA}" >&2;
    echo "  Run scripts/benchmark/prep_proteomes.sh first." >&2; exit 1; }
[[ -s "${HOME_GFF}" ]] || { echo "ERROR: Apis GFF missing: ${HOME_GFF}" >&2; exit 1; }

CALLS="${OUTDIR}/calls.tsv"
printf "species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra\n" > "${CALLS}"

# ── species enumeration ────────────────────────────────────────────────
# Build a map species -> (faa_path, gff_path). gff may be empty if missing.
declare -A SPECIES_FAA
declare -A SPECIES_GFF

# Tier 1/2/3 from local_data/benchmark/genomes/<tier>/<sp>/
for tier_dir in "${GENOMES_DIR}"/{tier1,tier2,tier3}; do
    [[ -d "${tier_dir}" ]] || continue
    for sd in "${tier_dir}"/*/; do
        sp=$(basename "${sd}")
        faa="${PROTEOMES_DIR}/${sp}.faa"
        gff=$(ls "${sd}"*.gff 2>/dev/null | head -1)
        SPECIES_FAA[$sp]="${faa}"
        SPECIES_GFF[$sp]="${gff:-}"
    done
done
# Existing 5 Koludarov species
for sp in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    faa="${PROTEOMES_DIR}/${sp}.faa"
    gff="${EXISTING_DIR}/${sp}.gff3"
    SPECIES_FAA[$sp]="${faa}"
    SPECIES_GFF[$sp]="${gff}"
done

# ── helpers ──────────────────────────────────────────────────────────────
species_prefix() {
    # Stable 2-letter species prefix: first letter of genus + first letter of
    # species (lowercase). Apis_mellifera → "am", Vespa_velutina → "vv".
    # Collisions resolved by appending a digit per call (we don't expect any
    # in this 18-species set).
    local sp="$1"
    local first_letter
    first_letter=$(echo "${sp}" | cut -c1 | tr '[:upper:]' '[:lower:]')
    local genus_species
    genus_species=$(echo "${sp}" | awk -F'_' '{print tolower(substr($1,1,1)) tolower(substr($2,1,1))}')
    echo "${genus_species}"
}

# Build MCScanX simplified GFF from a real GFF. Output rows:
#   <prefix>_<chrom>\t<sanitized_gene_id>\t<start>\t<end>
# Sanitization: '-' → '_' in gene IDs (MCScanX is strict about non-alnum),
# and chrom names lose non-alnum chars too.
build_simple_gff() {
    local gff_in="$1"
    local prefix="$2"
    local out="$3"
    awk -v PREFIX="${prefix}" -F'\t' '
        $1 ~ /^#/ {next}
        NF >= 9 && $3 == "gene" {
            split($9, a, ";")
            id=""
            for (k in a) {
                if (a[k] ~ /^ID=/) { sub(/^ID=/, "", a[k]); id=a[k]; break }
            }
            if (id != "") {
                gsub(/-/, "_", id)
                chrom=$1
                gsub(/[^A-Za-z0-9_]/, "_", chrom)
                print PREFIX "_" chrom "\t" id "\t" $4 "\t" $5
            }
        }
    ' "${gff_in}" > "${out}"
}

run_one() {
    local sp="$1"
    local sp_faa="${SPECIES_FAA[$sp]:-}"
    local sp_gff="${SPECIES_GFF[$sp]:-}"
    local sp_prefix
    sp_prefix=$(species_prefix "${sp}")
    if [[ "${sp_prefix}" == "${APIS_PREFIX}" ]]; then
        # Collision with Apis — bump query letter.
        sp_prefix="${sp_prefix}2"
    fi

    if [[ ! -s "${sp_faa}" ]]; then
        echo "  [${sp}] proteome missing (${sp_faa})"
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tmcscanx_no_proteome\n" "${sp}" >> "${CALLS}"
        return 0
    fi
    if [[ -z "${sp_gff}" || ! -s "${sp_gff}" ]]; then
        echo "  [${sp}] NO GFF (${sp_gff:-none}) — MCScanX cannot score, marking ABSENT"
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tmcscanx_no_annotation\n" "${sp}" >> "${CALLS}"
        return 0
    fi

    local sp_work="${WORK}/${sp}"
    mkdir -p "${sp_work}"
    local prefix_base="${sp_work}/Apis_vs_${sp}"
    local combined_faa="${prefix_base}.faa"
    local combined_gff="${prefix_base}.gff"
    local blast_out="${prefix_base}.blast"
    local collinearity="${prefix_base}.collinearity"

    # 1. Combined proteome (Apis + target). MCScanX matches gene IDs from the
    #    simplified GFF to the BLAST results, and only allows alphanumerics
    #    and underscores in gene IDs — so we sanitize `-` → `_` everywhere.
    #    `gene-Melt` becomes `gene_Melt` throughout (proteome headers and
    #    simplified GFF gene_id column).
    if [[ ! -s "${combined_faa}" ]]; then
        # Sanitize header IDs by replacing '-' with '_' in everything after '>'.
        cat "${APIS_FAA}" "${sp_faa}" | sed -E '/^>/ s/-/_/g' > "${combined_faa}"
    fi

    # 2. Simplified GFF for both sides. build_simple_gff already sanitizes IDs.
    if [[ ! -s "${combined_gff}" ]]; then
        local apis_simple="${sp_work}/_apis.simple.gff"
        local sp_simple="${sp_work}/_${sp}.simple.gff"
        build_simple_gff "${HOME_GFF}" "${APIS_PREFIX}" "${apis_simple}"
        build_simple_gff "${sp_gff}" "${sp_prefix}" "${sp_simple}"
        cat "${apis_simple}" "${sp_simple}" > "${combined_gff}"
        rm -f "${apis_simple}" "${sp_simple}"
    fi

    # 3. DIAMOND blastp all-vs-all (Apis ∪ target).
    if [[ ! -s "${blast_out}" ]]; then
        echo "  [${sp}] diamond blastp (this is the slow step, ~2-5 min)"
        local db="${prefix_base}.dmnd"
        diamond makedb --in "${combined_faa}" --db "${db}" --threads "${THREADS}" \
            > "${sp_work}/diamond_makedb.log" 2>&1 \
            || { echo "    diamond makedb FAILED"; return 1; }
        diamond blastp --query "${combined_faa}" --db "${db}" \
            --outfmt 6 -e 1e-10 -k 5 --threads "${THREADS}" \
            -o "${blast_out}" \
            > "${sp_work}/diamond_blastp.log" 2>&1 \
            || { echo "    diamond blastp FAILED"; return 1; }
        rm -f "${db}"
    fi

    # 4. Run MCScanX. It writes <prefix>.collinearity next to the input files.
    #    `-b 2` restricts output to INTER-species blocks (excludes Apis-Apis
    #    paralog blocks that would otherwise false-positive melittin via its
    #    own paralogs). Other params left at defaults.
    if [[ ! -s "${collinearity}" ]]; then
        echo "  [${sp}] MCScanX (inter-species collinear blocks, -b 2)"
        MCScanX -b 2 "${prefix_base}" > "${sp_work}/mcscanx.log" 2>&1 \
            || { echo "    MCScanX FAILED — see ${sp_work}/mcscanx.log"; return 1; }
    fi

    # 5. Parse: does any inter-species block contain gene_Melt (sanitized)?
    #    Anchor rows: "  N-  gene_A  gene_B  evalue"
    if [[ ! -s "${collinearity}" ]]; then
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tmcscanx_no_collinearity\n" "${sp}" >> "${CALLS}"
        return 0
    fi

    local melt_rows
    melt_rows=$(grep -E '[[:space:]]gene_Melt[[:space:]]' "${collinearity}" 2>/dev/null)
    if [[ -z "${melt_rows}" ]]; then
        printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\tmcscanx_melt_not_in_block\n" "${sp}" >> "${CALLS}"
        return 0
    fi

    # Pick the best-evalue row; partner is whichever of cols 2/3 isn't gene_Melt.
    local best partner ev
    best=$(echo "${melt_rows}" | awk '{print $0}' | sort -k4,4g | head -1)
    partner=$(awk '{
        for (i=2; i<=3; i++) {
            if ($i != "gene_Melt" && $i ~ /[A-Za-z]/) { print $i; exit }
        }
    }' <<<"${best}")
    ev=$(awk '{print $NF}' <<<"${best}")
    if [[ -z "${partner}" ]]; then
        printf "%s\t-\tAMBIGUOUS\t-\t-\t-\t-\t%s\tmcscanx_partner_unresolved\n" "${sp}" "${ev}" >> "${CALLS}"
        return 0
    fi
    printf "%s\t-\tPRESENT\t-\t-\t-\t-\t%s\tmcscanx_anchor=%s\n" \
        "${sp}" "${ev}" "${partner}" >> "${CALLS}"
}

# ── orchestration ───────────────────────────────────────────────────────
n_total=${#SPECIES_FAA[@]}
i=0
for sp in $(echo "${!SPECIES_FAA[@]}" | tr ' ' '\n' | sort); do
    i=$((i + 1))
    echo "[${i}/${n_total}] ${sp}"
    run_one "${sp}" || echo "  [${sp}] failed (continuing)"
done

echo
echo "Done. Calls: ${CALLS}"
echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${CALLS}" | sort | uniq -c
echo
echo "ABSENT-with-mcscanx_no_annotation rows (the comparison-of-interest):"
awk -F'\t' '$3=="ABSENT" && $9 ~ /no_annotation/ {print "  " $1}' "${CALLS}"
