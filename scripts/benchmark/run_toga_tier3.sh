#!/usr/bin/env bash
# TOGA1 on Tier 3 (deep aculeate homology). The legitimate competitor for
# SynVoy in the synteny-aware-genome-projection space.
#
# We use TOGA1 (Kirilenko et al., Science 2023) rather than TOGA2. TOGA2 is
# early-access / unpublished; for a benchmark in a peer-reviewed paper, the
# Science 2023 implementation is what reviewers can look up and verify.
#
# WHERE TO RUN: this script is designed for a host with ≥22 GB RAM
# (your other PC, NOT the 16 GB laptop). LASTZ chain generation on
# whole-genome pairs blows out memory on smaller hosts.
#
# RESOURCES:
#   - ~5-15 hours of compute for 7 Tier 3 species
#   - ~30 GB scratch disk (chain files + 2bit + intermediate)
#   - 4-8 CPU cores (per make_lastz_chains recommendation)
#
# EXPECTED OUTPUT: benchmark_results/toga/calls.tsv with the same contract
# as the other competitor wrappers, comparable via score_benchmark.py.
#
# RESUMABILITY: every stage skips work that's already done. Re-running is
# safe and cheap.
#
# Usage on the target PC:
#   git pull   # pull latest scripts
#   tmux new -s toga
#   bash scripts/benchmark/run_toga_tier3.sh
#   # detach: Ctrl-b d ; reattach: tmux attach -t toga

set -uo pipefail

# ── paths ────────────────────────────────────────────────────────────────
REPO="${REPO:-/home/faw/dev/projects/SynVoy}"
GENOMES_DIR="${GENOMES_DIR:-${REPO}/local_data/benchmark/genomes}"
OUTDIR="${OUTDIR:-${REPO}/benchmark_results/toga}"
WORK="${OUTDIR}/work"
HOME_FNA="${HOME_FNA:-${REPO}/local_data/benchmark/genomes/home/Apis_mellifera/Apis_mellifera.fna}"
HOME_GFF="${HOME_GFF:-${REPO}/local_data/benchmark/genomes/home/Apis_mellifera/Apis_mellifera.gff}"
QUERY="${QUERY:-${REPO}/local_data/ground_truth/melettin/query_melittin.faa}"

# Tier 3 species — paths under genomes_dir/tier3/<species>/<species>.fna
TIER3_SPECIES=(
    Tetramorium_bicarinatum
    Cardiocondyla_obscurior
    Formica_selysi
    Vespa_velutina
    Polistes_canadensis
    Solenopsis_invicta
    Vollenhovia_emeryi
)

THREADS="${THREADS:-4}"
LASTZ_JOB_MEMORY_GB="${LASTZ_JOB_MEMORY_GB:-8}"
LASTZ_JOB_TIME_HOURS="${LASTZ_JOB_TIME_HOURS:-12}"
MLC_MAX_RETRIES="${MLC_MAX_RETRIES:-0}"

# LASTZ params tuned for distant-homology comparisons (Apis vs ants ≈ 155 Mya).
# Defaults in make_lastz_chains are vertebrate-tuned (K=2400 L=3000 H=2000
# Y=9400) and produce ZERO alignments for Hymenoptera family-level pairs.
# These values follow Hiller lab's distant-homology recommendations
# (Kirilenko et al. 2023 Supp. for >100My; HoxD55-style permissive).
# Override per-run via env vars if you want to tighten back up.
LASTZ_K="${LASTZ_K:-1500}"        # seed score threshold (lower = more sensitive)
LASTZ_L="${LASTZ_L:-2500}"        # gap-free extension threshold
LASTZ_H="${LASTZ_H:-0}"           # HSP score filter (0 = disable)
LASTZ_Y="${LASTZ_Y:-3400}"        # chain gap penalty cap
# Target (Apis) chunk size for make_lastz_chains. The make_lastz_chains default
# is 175 Mb — i.e. the WHOLE Apis genome goes into a single lastz job; combined
# with sensitive K that produces 100M+ HSPs and overflows lastz's 32-bit segment
# table. Smaller chunks bound each job's segment count (belt-and-suspenders on
# top of the soft-masking fix below) and parallelize better. Override via env.
LASTZ_SEQ1_CHUNK="${LASTZ_SEQ1_CHUNK:-30000000}"
# These default to the env + paths created by scripts/benchmark/toga_setup.sh.
# Run toga_setup.sh ONCE before this script. Override to reuse a different env.
TOGA_ENV="${TOGA_ENV:-synvoy_toga}"
LASTZ_ENV="${LASTZ_ENV:-synvoy_toga}"
TOOLS_DIR="${TOOLS_DIR:-${HOME}/dev/tools}"
MAKE_CHAINS_PY="${MAKE_CHAINS_PY:-${TOOLS_DIR}/make_lastz_chains/make_chains.py}"
TOGA_PY="${TOGA_PY:-${TOOLS_DIR}/TOGA/toga.py}"

# To run a subset (e.g. a single-species smoke test on a 16 GB laptop):
#   ONLY=Cardiocondyla_obscurior bash run_toga_tier3.sh
#   ONLY="Cardiocondyla_obscurior Vespa_velutina" bash run_toga_tier3.sh
ONLY="${ONLY:-}"

mkdir -p "${WORK}"/{2bit,chains,toga_runs,bed,logs}

patch_make_lastz_chains() {
    # make_lastz_chains 2.0.8 ignores the memory/time arguments passed to
    # NextflowConfig and instead hard-codes 16 GB and 1 hour. On a 16 GB laptop
    # this prevents retry submission, and the 1h cap kills legitimate LASTZ
    # jobs. Patch the local checkout in-place, controlled by env vars below.
    local wrapper
    wrapper="$(dirname "${MAKE_CHAINS_PY}")/parallelization/nextflow_wrapper.py"
    if [[ ! -f "${wrapper}" ]]; then
        echo "  ERROR: nextflow wrapper not found at ${wrapper}" >&2
        return 1
    fi

    python3 - "${wrapper}" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text()
original = text

if "import os\n" not in text:
    text = text.replace(
        '"""Module to manage Nextflow processes."""\n\n',
        '"""Module to manage Nextflow processes."""\n\nimport os\n',
        1,
    )

text = re.sub(
    r"(?m)^(\s*)self\.memory = .*$",
    r'\1self.memory = int(os.environ.get("MLC_JOB_MEMORY_GB", memory))',
    text,
    count=1,
)
text = re.sub(
    r"(?m)^(?P<indent>\s*)(?:raw_time = os\.environ\.get\(\"MLC_JOB_TIME_HOURS\", time\)\n(?P=indent))?self\.time = .*$",
    lambda match: (
        f'{match.group("indent")}raw_time = os.environ.get("MLC_JOB_TIME_HOURS", time)\n'
        f'{match.group("indent")}self.time = int(str(raw_time).rstrip("h"))'
    ),
    text,
    count=1,
)
text = re.sub(
    r"(?m)^(\s*)self\.maxRetries = .*$",
    r'\1self.maxRetries = int(os.environ.get("MLC_MAX_RETRIES", "0"))',
    text,
    count=1,
)
text = re.sub(
    r'f\.write\(f"    memory = .*?\\n"\)',
    lambda _: 'f.write(f"    memory = {self.memory}.GB\\n")',
    text,
    count=1,
)
text = re.sub(
    r'f\.write\(f"    time = .*?\\n"\)',
    lambda _: 'f.write(f"    time = {self.time}.hour\\n")',
    text,
    count=1,
)

if text != original:
    path.write_text(text)
PY

    export MLC_JOB_MEMORY_GB="${LASTZ_JOB_MEMORY_GB}"
    export MLC_JOB_TIME_HOURS="${LASTZ_JOB_TIME_HOURS}"
    export MLC_MAX_RETRIES
    echo "  patched make_lastz_chains Nextflow wrapper:"
    echo "    memory=${MLC_JOB_MEMORY_GB} GB, time=${MLC_JOB_TIME_HOURS} h, maxRetries=${MLC_MAX_RETRIES}"

    # ── soft-masking fix ──────────────────────────────────────────────────
    # make_lastz_chains 2.0.8 (upstream PR #110 "fix/py2bit-large-genome")
    # extracts each lastz partition with py2bit's tb.sequence(), which silently
    # STRIPS soft-masking (returns all-uppercase; only hard-N blocks survive).
    # The genomes ARE soft-masked (~27-42% lowercase), so lastz ends up seeding
    # across all repeats -> ~100M HSPs -> overflow of lastz's 32-bit segment
    # table ("table size ... exceeds allocation limit of 4,294,967,279") at the
    # sensitive K we need for Apis-vs-ant (~155 My) homology. twoBitToFa
    # preserves lowercase soft-masking by default AND reads v1/-long 2bits, so
    # it keeps PR #110's large-genome support. Patch the extractor in place.
    local run_lastz_py
    run_lastz_py="$(dirname "${MAKE_CHAINS_PY}")/standalone_scripts/run_lastz.py"
    if [[ ! -f "${run_lastz_py}" ]]; then
        echo "  ERROR: run_lastz.py not found at ${run_lastz_py}" >&2
        return 1
    fi
    python3 - "${run_lastz_py}" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text()
if "SYNVOY_MASK_FIX" in text:
    sys.exit(0)  # idempotent: already patched

new_func = '''def extract_twobit_partition(two_bit_path, chrom, start, end, tmp_dir):
    """Extract one partition from a .2bit to temp FASTA via twoBitToFa.

    SYNVOY_MASK_FIX: upstream PR #110 switched this to py2bit, whose
    .sequence() silently drops soft-masking (returns uppercase). Unmasked
    input makes lastz seed across repeats and overflow its 32-bit segment
    table at sensitive K. twoBitToFa preserves lowercase soft-masking and
    still reads v1/-long 2bits, so PR #110's large-genome support is kept.
    0-based half-open coords match the partition strings.
    """
    import subprocess as _subprocess
    fasta_path = os.path.join(tmp_dir, f"{_gen_random_string(8)}_partition.fa")
    raw_path = os.path.join(tmp_dir, f"{_gen_random_string(8)}_raw.fa")
    _subprocess.run(
        ["twoBitToFa", f"-seq={chrom}", f"-start={start}", f"-end={end}",
         two_bit_path, raw_path],
        check=True,
    )
    # twoBitToFa emits ">chrom:start-end"; rewrite to bare ">chrom" so the
    # sequence name matches what the rest of the pipeline expects.
    with open(raw_path) as _fin, open(fasta_path, "w") as _fout:
        for _i, _line in enumerate(_fin):
            if _i == 0:
                print(">" + chrom, file=_fout)  # bare name (twoBitToFa wrote ">chrom:start-end")
            else:
                _fout.write(_line)
    os.remove(raw_path)
    return fasta_path
'''

text, n = re.subn(
    r"def extract_twobit_partition\(.*?\n    return fasta_path\n",
    new_func,
    text,
    count=1,
    flags=re.DOTALL,
)
if n != 1:
    sys.stderr.write("ERROR: could not locate extract_twobit_partition to patch\n")
    sys.exit(2)
path.write_text(text)
PY
    local mask_rc=$?
    if (( mask_rc != 0 )); then
        echo "  ERROR: run_lastz.py soft-masking patch failed (rc=${mask_rc})" >&2
        return 1
    fi
    echo "  patched run_lastz.py: twoBitToFa partition extraction (preserves soft-masking)"
}

# ── stage 1: verify env exists (created separately by toga_setup.sh) ─────
setup_envs() {
    echo "=== [stage 1] verify TOGA env (created by toga_setup.sh) ==="
    # Modern mamba env list has a 2-space leading indent on env names; match
    # both the indented table format and any older single-column format.
    if ! mamba env list 2>/dev/null | grep -qE "^\s*${TOGA_ENV}(\s|$)"; then
        echo "  ERROR: conda env '${TOGA_ENV}' not found." >&2
        echo "  Run scripts/benchmark/toga_setup.sh first to install it." >&2
        return 1
    fi
    if [[ ! -x "${MAKE_CHAINS_PY}" && ! -f "${MAKE_CHAINS_PY}" ]]; then
        echo "  ERROR: make_chains.py not at ${MAKE_CHAINS_PY}" >&2
        echo "  Run scripts/benchmark/toga_setup.sh first." >&2
        return 1
    fi
    if [[ ! -x "${TOGA_PY}" && ! -f "${TOGA_PY}" ]]; then
        echo "  ERROR: TOGA1 toga.py not at ${TOGA_PY}" >&2
        echo "  Run scripts/benchmark/toga_setup.sh first." >&2
        return 1
    fi
    echo "  env: ${TOGA_ENV}"
    echo "  make_chains.py: ${MAKE_CHAINS_PY}"
    echo "  TOGA1 toga.py: ${TOGA_PY}"
    patch_make_lastz_chains || return 1
}

# ── stage 2: reference prep ──────────────────────────────────────────────
prep_reference() {
    echo "=== [stage 2] prep Apis (reference) inputs ==="
    if [[ ! -s "${HOME_FNA}" || ! -s "${HOME_GFF}" ]]; then
        echo "  ERROR: Apis mellifera home inputs missing." >&2
        echo "    HOME_FNA: ${HOME_FNA}" >&2
        echo "    HOME_GFF: ${HOME_GFF}" >&2
        echo "  Fetch them first with:" >&2
        echo "    bash scripts/benchmark/fetch_benchmark_genomes.sh" >&2
        return 1
    fi
    local missing_tier3=()
    for sp in "${TIER3_SPECIES[@]}"; do
        [[ -s "${GENOMES_DIR}/tier3/${sp}/${sp}.fna" ]] || missing_tier3+=("${sp}")
    done
    if (( ${#missing_tier3[@]} > 0 )); then
        echo "  WARNING: missing Tier 3 genomes (will be skipped at stage 3): ${missing_tier3[*]}" >&2
        echo "  Fetch with: bash scripts/benchmark/fetch_benchmark_genomes.sh" >&2
    fi

    local apis_2bit="${WORK}/2bit/Apis_mellifera.2bit"
    local apis_bed="${WORK}/bed/Apis_mellifera.melittin.bed12"

    if [[ ! -s "${apis_2bit}" ]]; then
        echo "  converting Apis genome -> 2bit (stripping NCBI .N version suffixes)"
        # make_lastz_chains rejects dot-containing chrom names and rewrites the
        # 2bit with `-long`, producing a v1 file that its own twobitreader can't
        # read (known incompatibility). Side-step by stripping NCBI version
        # suffixes (.1, .2, ...) before conversion. Names become NC_037638
        # instead of NC_037638.1 — coords stay identical.
        local fa_clean="${WORK}/2bit/Apis_mellifera.clean.fa"
        if [[ ! -s "${fa_clean}" ]]; then
            sed -E 's/^(>[^[:space:].]+)\.[0-9]+(.*)/\1\2/' "${HOME_FNA}" > "${fa_clean}"
        fi
        mamba run -n "${LASTZ_ENV}" faToTwoBit "${fa_clean}" "${apis_2bit}" \
            2> "${WORK}/logs/apis_2bit.log" || return 1
        rm -f "${fa_clean}"  # save disk: ~250 MB of fasta we don't need
    fi
    if [[ ! -s "${apis_bed}" ]]; then
        echo "  extracting Apis melittin BED12 from GFF"
        # GFF hierarchy: gene → mRNA(Parent=gene) → CDS(Parent=mRNA).
        # We need to follow the chain via IDs, not by attribute substring.
        python3 - "${HOME_GFF}" "${apis_bed}" <<'PY' 2> "${WORK}/logs/apis_bed.log"
import re, sys
gff_path, bed_path = sys.argv[1:]

def attrs_dict(s):
    return {k: v for k, _, v in (a.partition("=") for a in s.rstrip(";").split(";"))}

mrna_ids_for_melt = set()  # mRNA IDs whose Parent is gene-Melt
gene_row = None
mrnas = {}  # mrna_id -> (chrom, start, end, strand)
cds_by_mrna = {}  # mrna_id -> [(start, end), ...]

with open(gff_path) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        a = attrs_dict(f[8])
        if f[2] == "gene" and a.get("ID") == "gene-Melt":
            gene_row = (f[0], int(f[3]) - 1, int(f[4]), f[6])
        elif f[2] == "mRNA" and a.get("Parent") == "gene-Melt":
            mid = a.get("ID")
            mrna_ids_for_melt.add(mid)
            mrnas[mid] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
        elif f[2] == "CDS":
            parent = a.get("Parent", "")
            for pid in parent.split(","):
                if pid in mrna_ids_for_melt:
                    cds_by_mrna.setdefault(pid, []).append((int(f[3]) - 1, int(f[4])))

if gene_row is None:
    sys.exit("ERROR: gene-Melt not found in GFF")
if not mrnas:
    sys.exit("ERROR: no mRNA with Parent=gene-Melt found")

# Pick the mRNA with the most CDS exons (canonical longest isoform).
best_mid = max(mrnas, key=lambda m: len(cds_by_mrna.get(m, [])))
chrom, start, end, strand = mrnas[best_mid]
# Strip NCBI .N version suffix from chrom name to match the cleaned 2bit
chrom = re.sub(r"\.[0-9]+$", "", chrom)
cds = sorted(cds_by_mrna.get(best_mid, []))
if not cds:
    sys.exit(f"ERROR: no CDS rows for {best_mid}; "
             f"GFF parse failed? Found mRNAs: {list(mrnas.keys())}")

block_sizes = [str(b - a) for a, b in cds]
block_starts = [str(a - start) for a, _ in cds]
# BED12 order: chrom start end name score strand thickStart thickEnd rgb
#              blockCount blockSizes blockStarts
bed = "\t".join([
    chrom, str(start), str(end), "gene-Melt", "0", strand,
    str(cds[0][0]), str(cds[-1][1]), "0",
    str(len(cds)),
    ",".join(block_sizes) + ",",
    ",".join(block_starts) + ",",
])
with open(bed_path, "w") as out:
    out.write(bed + "\n")
print(f"wrote {bed_path}: {bed}")
print(f"  mRNA: {best_mid}, {len(cds)} CDS blocks", file=sys.stderr)
PY
    fi
    echo "  reference: ${apis_2bit}"
    echo "  bed12:     ${apis_bed}"
}

# ── stage 3: per-target chain generation + TOGA run ──────────────────────
run_target() {
    local sp="$1"
    local target_fna="${GENOMES_DIR}/tier3/${sp}/${sp}.fna"
    if [[ ! -s "${target_fna}" ]]; then
        echo "  [${sp}] missing genome FASTA, skipping"
        return 0
    fi

    local target_2bit="${WORK}/2bit/${sp}.2bit"
    local chain_file="${WORK}/chains/Apis_mellifera.${sp}.final.chain.gz"
    local toga_out="${WORK}/toga_runs/${sp}"
    mkdir -p "${toga_out}"

    # 3a. 2bit conversion (~1 min). Strip NCBI version suffixes from chrom
    # names to keep make_lastz_chains' twobitreader happy (see prep_reference).
    if [[ ! -s "${target_2bit}" ]]; then
        echo "  [${sp}] convert -> 2bit (stripping version suffixes)"
        local fa_clean="${WORK}/2bit/${sp}.clean.fa"
        sed -E 's/^(>[^[:space:].]+)\.[0-9]+(.*)/\1\2/' "${target_fna}" > "${fa_clean}"
        mamba run -n "${LASTZ_ENV}" faToTwoBit "${fa_clean}" "${target_2bit}" \
            2> "${WORK}/logs/${sp}_2bit.log" || { echo "    2bit FAILED"; rm -f "${fa_clean}"; return 1; }
        rm -f "${fa_clean}"
    fi

    # 3b. chain generation (the slow step: ~30min - several hours per pair)
    if [[ ! -s "${chain_file}" ]]; then
        echo "  [${sp}] generating chain (LASTZ — this is the slow step, expect 30min-3h)"
        local chain_proj="${WORK}/chains/_run_${sp}"
        # IMPORTANT: do NOT pre-create chain_proj — make_chains.py creates it
        # itself and errors if it already exists. Use --force to overwrite if
        # a previous failed run left a partial dir.
        rm -rf "${chain_proj}"
        # make_chains.py CLI: positional (target_name, query_name, target.2bit,
        # query.2bit) + flags. --pd (project dir), --cluster_executor (executor
        # type), --force (overwrite existing project dir).
        mamba run -n "${LASTZ_ENV}" \
            python3 "${MAKE_CHAINS_PY}" \
                "Apis_mellifera" "${sp}" \
                "${WORK}/2bit/Apis_mellifera.2bit" "${target_2bit}" \
                --pd "${chain_proj}" \
                --cluster_executor local \
                --num_fill_jobs "${THREADS}" \
                --lastz_k "${LASTZ_K}" \
                --lastz_l "${LASTZ_L}" \
                --lastz_h "${LASTZ_H}" \
                --lastz_y "${LASTZ_Y}" \
                --seq1_chunk "${LASTZ_SEQ1_CHUNK}" \
                --force \
            > "${WORK}/logs/${sp}_chain.log" 2>&1 \
            || { echo "    chain generation FAILED — see ${WORK}/logs/${sp}_chain.log"; return 1; }
        # The pipeline puts the result at: ${chain_proj}/Apis_mellifera.${sp}.final.chain.gz
        if [[ -s "${chain_proj}/Apis_mellifera.${sp}.final.chain.gz" ]]; then
            mv "${chain_proj}/Apis_mellifera.${sp}.final.chain.gz" "${chain_file}"
        else
            echo "    chain file not found at expected path; check ${chain_proj}"
            return 1
        fi
    fi

    # 3c. run TOGA1 (~10-30 min per pair)
    # TOGA1 CLI: positional args (chain, bed12, ref2bit, query2bit) + flags
    #   --chn N    chain-extract jobs (parallelism)
    #   --cjn N    CESAR alignment jobs (parallelism)
    #   --pn / --project_name  output dir
    if [[ ! -s "${toga_out}/orthology_classification.tsv" ]]; then
        echo "  [${sp}] running TOGA1 (orthology projection)"
        mamba run -n "${TOGA_ENV}" python3 "${TOGA_PY}" \
            "${chain_file}" \
            "${WORK}/bed/Apis_mellifera.melittin.bed12" \
            "${WORK}/2bit/Apis_mellifera.2bit" \
            "${target_2bit}" \
            --project_name "${sp}" \
            --project_dir "${toga_out}" \
            --chain_jobs_num "${THREADS}" \
            --cesar_jobs_num "${THREADS}" \
            > "${WORK}/logs/${sp}_toga.log" 2>&1 \
            || { echo "    TOGA FAILED — see ${WORK}/logs/${sp}_toga.log"; return 1; }
    fi

    echo "  [${sp}] DONE"
}

# ── stage 4: parse TOGA outputs into normalized calls.tsv ────────────────
parse_outputs() {
    echo "=== [stage 4] parse TOGA outputs into calls.tsv ==="
    local calls="${OUTDIR}/calls.tsv"
    printf "species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra\n" > "${calls}"

    for sp in "${TIER3_SPECIES[@]}"; do
        local toga_out="${WORK}/toga_runs/${sp}"
        # TOGA1 emits orthology_classification.tsv with classes:
        #   one2one, one2many, many2one, many2many, one2zero
        # one2zero = no ortholog projected → ABSENT
        # everything else (one2one, etc.) = ortholog present → PRESENT
        local class_tsv="${toga_out}/orthology_classification.tsv"
        if [[ ! -s "${class_tsv}" ]]; then
            printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\ttoga_no_output\n" "${sp}" >> "${calls}"
            continue
        fi
        # Look for the gene-Melt row
        local row
        row=$(grep -E "gene-Melt|Melt" "${class_tsv}" | head -1)
        if [[ -z "${row}" ]]; then
            printf "%s\t-\tABSENT\t-\t-\t-\t-\t-\ttoga_no_melittin_in_classification\n" "${sp}" >> "${calls}"
            continue
        fi
        local class_label
        class_label=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i~/^one2(one|many|zero)$|^many2(one|many)$/){print $i; exit}}' <<<"${row}")
        local status
        case "${class_label}" in
            one2zero|"") status="ABSENT" ;;
            *) status="PRESENT" ;;
        esac
        printf "%s\t-\t%s\t-\t-\t-\t-\t%s\ttoga1_class=%s\n" \
            "${sp}" "${status}" "${class_label:--}" "${class_label:--}" >> "${calls}"
    done

    echo "  wrote ${calls}"
    awk -F'\t' 'NR>1 {print $3}' "${calls}" | sort | uniq -c
}

# ── orchestration ────────────────────────────────────────────────────────
echo "=== TOGA-Tier3 benchmark on $(hostname) at $(date -Iseconds) ==="
echo "Memory:  $(free -h | awk '/^Mem:/{print $2}') total, $(nproc) CPUs"
echo

setup_envs        || { echo "STAGE 1 FAILED"; exit 1; }
prep_reference    || { echo "STAGE 2 FAILED"; exit 1; }

echo "=== [stage 3] per-target chain + TOGA runs ==="
# Apply ONLY filter if set
if [[ -n "${ONLY}" ]]; then
    declare -A keep
    for sp in ${ONLY}; do keep[$sp]=1; done
    SPECIES_TO_RUN=()
    for sp in "${TIER3_SPECIES[@]}"; do
        [[ -n "${keep[$sp]:-}" ]] && SPECIES_TO_RUN+=("$sp")
    done
    echo "  ONLY filter active — running: ${SPECIES_TO_RUN[*]}"
else
    SPECIES_TO_RUN=("${TIER3_SPECIES[@]}")
fi

for sp in "${SPECIES_TO_RUN[@]}"; do
    run_target "${sp}" || echo "  [${sp}] failed (continuing)"
done

parse_outputs
echo
echo "Done. calls.tsv: ${OUTDIR}/calls.tsv"
echo "After this script finishes on the 22GB PC, copy benchmark_results/toga/calls.tsv"
echo "back to the laptop and re-score:"
echo "  python scripts/benchmark/score_benchmark.py \\"
echo "    --truth tests/benchmark_truth/melittin_orthologs_gr2.tsv \\"
echo "    --calls-glob 'benchmark_results/*/calls.tsv' \\"
echo "    --outdir benchmark_results/_gr2"
