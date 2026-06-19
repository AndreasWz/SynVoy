#!/bin/bash
# =============================================================================
#  run_synvoy.sh — the one safe way to launch SynVoy
# =============================================================================
#  Designed for non-expert users (the Bachelor cohort) and weak laptops. It:
#    1. prints the exact SynVoy version you are about to run (no more guessing);
#    2. REFUSES to run a dangerously old checkout (the stale-code → 49-locus
#       runaway → OOM that has been crashing everyone);
#    3. runs a pre-flight check of every tool + Python package and FAILS CLEARLY
#       with the one command to fix it (it never silently rebuilds your env);
#    4. launches Nextflow with a memory-safe profile by default (auto,low_mem).
#
#  Usage (easy mode, auto-fetch genomes — the common case):
#     ./run_synvoy.sh --mode easy --query_id Q16553 --max_genomes 5 \
#                     --outdir results/my_run
#
#  Pro mode (local files):
#     ./run_synvoy.sh --mode pro --query q.faa --home_genome g.fna \
#                     --home_gff g.gff --target_genomes 'targets/*.fna' \
#                     --outdir results/my_run
#
#  Just check my setup, don't run anything:
#     ./run_synvoy.sh --check-only
#
#  Override the profile (e.g. a beefy workstation):
#     ./run_synvoy.sh -profile standard --mode easy --query_id Q16553 ...
#
#  Everything after the script name is passed straight through to
#  `nextflow run main.nf`, so any normal SynVoy flag works.
# =============================================================================
# NB: intentionally NOT using `set -u` — `conda activate` references unbound
# shell vars and would abort the launcher. We default-initialise our own vars.
set -o pipefail

ENV_NAME="synvoy_env"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || { echo "ERROR: cannot cd to repo dir $SCRIPT_DIR"; exit 1; }

# ── colours (skip if not a TTY) ──────────────────────────────────────────────
if [ -t 1 ]; then BOLD=$'\033[1m'; RED=$'\033[31m'; GRN=$'\033[32m'; YEL=$'\033[33m'; RST=$'\033[0m'
else BOLD=""; RED=""; GRN=""; YEL=""; RST=""; fi
say()  { echo "${GRN}[ok]${RST}   $*"; }
warn() { echo "${YEL}[warn]${RST} $*"; }
die()  { echo "${RED}${BOLD}[FAIL]${RST} $*" >&2; exit 1; }

CHECK_ONLY=0
for a in "$@"; do [ "$a" = "--check-only" ] && CHECK_ONLY=1; done

# ── 1. version banner ────────────────────────────────────────────────────────
VERSION="$(cat VERSION 2>/dev/null || echo 'unknown')"
GIT_SHA="$(git rev-parse --short HEAD 2>/dev/null || echo 'no-git')"
GIT_BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
GIT_DIRTY=""
git diff --quiet 2>/dev/null || GIT_DIRTY=" ${YEL}(local edits)${RST}"
echo "${BOLD}========================================================${RST}"
echo "${BOLD}  SynVoy v${VERSION}${RST}  —  ${GIT_BRANCH}@${GIT_SHA}${GIT_DIRTY}"
echo "${BOLD}========================================================${RST}"

# ── 2. stale-code guard (network-free, the anti-OOM gate) ────────────────────
# These two markers are the OOM-critical fixes. If the checkout predates them,
# a gene-family query (luciferase, LY6, MRJP) fans out to dozens of loci and
# OOMs the machine. Refuse to run rather than crash hours later.
STALE=0
grep -q -- '--max_loci' modules/split_loci.nf 2>/dev/null || STALE=1
grep -q 'skip_tree'      nextflow.config       2>/dev/null || STALE=1
if [ "$STALE" -eq 1 ]; then
    echo
    echo "${RED}${BOLD}  ###############################################################${RST}"
    echo "${RED}${BOLD}  ##  WARNING: this SynVoy checkout looks DANGEROUSLY OLD       ##${RST}"
    echo "${RED}${BOLD}  ###############################################################${RST}"
    echo "${YEL}  It is missing the multi-locus safety cap and/or the low-memory${RST}"
    echo "${YEL}  tree switch. On a gene-family query (luciferase, LY6, MRJP) this${RST}"
    echo "${YEL}  is the #1 cause of out-of-memory crashes. Strongly recommended:${RST}"
    echo
    echo "      ${BOLD}cd $SCRIPT_DIR && git pull${RST}"
    echo "      ${BOLD}(if 'git pull' complains about local changes: git stash, then git pull)${RST}"
    echo
    echo "${YEL}  Continuing anyway in 5s (Ctrl-C to abort)...${RST}"
    sleep 5
else
    say "version guard passed (OOM-critical fixes present)"
fi

# Best-effort: warn if behind the remote. Short timeout, never fatal, never blocks.
if command -v git >/dev/null 2>&1 && [ -d .git ]; then
    if timeout 8 git fetch --quiet 2>/dev/null; then
        BEHIND="$(git rev-list --count HEAD..@{u} 2>/dev/null || echo 0)"
        if [ "${BEHIND:-0}" -gt 0 ]; then
            warn "you are ${BEHIND} commit(s) behind the remote — consider: git pull"
        fi
    fi
fi

# ── 3. activate the conda env ────────────────────────────────────────────────
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    if conda env list 2>/dev/null | grep -qw "$ENV_NAME"; then
        conda activate "$ENV_NAME" 2>/dev/null && say "activated conda env '$ENV_NAME'"
    else
        die "conda env '$ENV_NAME' not found. Create it once with:

           ${BOLD}cd $SCRIPT_DIR && ./install.sh${RST}
        or  ${BOLD}mamba env create -f environment.yml${RST}"
    fi
else
    warn "conda not on PATH — using whatever tools are already on PATH"
fi

# ── 4. pre-flight check (fails clearly, never auto-rebuilds) ──────────────────
# IMPORTANT: the bioinformatics tools (mmseqs, samtools, augustus, miniprot, …)
# and Python packages are provisioned by Nextflow itself from environment.yml
# (process.conda) — they do NOT need to be on your PATH. So here we only HARD-
# CHECK what *you* must supply to launch: conda/mamba, Java 17+, Nextflow, and a
# parseable config. Tool presence is reported as info, never fatal. (The full
# tool audit lives in ./install.sh, run once at setup.)

# 4a. Launch essentials (fatal).
command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1 \
    || die "conda/mamba not found. Install Miniforge: https://github.com/conda-forge/miniforge"
command -v nextflow >/dev/null 2>&1 \
    || die "nextflow not found. It ships inside '$ENV_NAME' — run ./install.sh, or 'conda activate $ENV_NAME'."
if command -v java >/dev/null 2>&1; then
    JV="$(java -version 2>&1 | head -1 | sed -E 's/.*version "([0-9]+).*/\1/')"
    [ "${JV:-0}" -ge 17 ] 2>/dev/null \
        || die "Java >=17 required by current Nextflow (found: $(java -version 2>&1 | head -1))."
else
    die "Java not found. Java >=17 ships inside '$ENV_NAME' — run ./install.sh."
fi
say "launch essentials present (conda, Nextflow, Java >=$JV)"

# 4b. Config must parse (catches a broken nextflow.config before any compute).
nextflow config -profile auto,low_mem >/dev/null 2>&1 \
    || die "nextflow.config failed to parse for 'auto,low_mem'. Run 'git pull' or report this output:
$(nextflow config -profile auto,low_mem 2>&1 | tail -5)"
say "config parses (auto,low_mem)"

# 4c. Tool/package probe — INFO ONLY. Missing here is fine: Nextflow builds them
# from environment.yml on first run. We surface it so you know what the first
# run will install (and to spot a genuinely broken setup early).
ABSENT=""
for b in mmseqs tblastn makeblastdb prodigal augustus miniprot mafft samtools datasets; do
    command -v "$b" >/dev/null 2>&1 || ABSENT="$ABSENT $b"
done
if [ -n "$ABSENT" ]; then
    warn "tools not on PATH:$ABSENT"
    warn "  ^ this is normal — Nextflow installs them from environment.yml on the first run"
    warn "    (first run therefore needs internet + a few minutes to build the tool env)"
fi

if [ "$CHECK_ONLY" -eq 1 ]; then
    echo
    say "${BOLD}--check-only: environment is healthy. You're ready to run.${RST}"
    exit 0
fi

# ── 5. launch ────────────────────────────────────────────────────────────────
# Default to the memory-safe profile unless the user supplied their own -profile.
PROFILE_ARGS=()
case " $* " in
    *" -profile "*|*" --profile "*) : ;;                       # user chose a profile
    *) PROFILE_ARGS=(-profile auto,low_mem)
       warn "no -profile given → using memory-safe default: ${BOLD}auto,low_mem${RST}"
       warn "  (on a powerful machine add: -profile standard  for full speed/tree)"
       ;;
esac

echo
echo "${BOLD}Launching:${RST} nextflow run main.nf ${PROFILE_ARGS[*]} $*"
echo
exec nextflow run main.nf "${PROFILE_ARGS[@]}" "$@"
