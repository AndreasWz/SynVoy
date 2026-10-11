#!/bin/bash
# run_job.sh <state_dir> <command> [args...]
#
# Runs one search for the web worker and writes the command's exit code to
# <state_dir>/exit_code when it ends. The worker learns the outcome from that file,
# not from being the parent process, so a search survives a restart of the worker.
#
# With SYNVOY_WEB_CONDA_ENV set, that conda environment is activated first (it holds
# Nextflow and Java; the pipeline builds the environment of its tasks itself).
#
# The PATH a search starts with decides which python3 its tasks run: a task switches
# to its own conda environment by replacing the active one's entry in PATH where it
# stands. So a Python virtualenv, or anything else with a python3, that sits in front
# of the conda environment stays in front inside every task and hides the pipeline's
# packages (parasail goes missing; see run_synvoy.sh, section 4d). Two measures:
# every virtualenv is taken out of PATH, whether VIRTUAL_ENV names it or not, and the
# activated conda environment is moved to the very front.
#
# Exit code 96 = the environment could not be set up (the search never started).
set -o pipefail

state_dir="$1"
shift

finish() {
    printf '%s\n' "$1" > "$state_dir/exit_code.tmp" && mv -f "$state_dir/exit_code.tmp" "$state_dir/exit_code"
    exit "$1"
}

# PATH without the entries that are the bin folder of a Python virtualenv.
clean_path=""
IFS=':' read -r -a path_entries <<< "$PATH"
for entry in "${path_entries[@]}"; do
    [ -z "$entry" ] && continue
    [ -f "$entry/../pyvenv.cfg" ] && continue
    [ -n "${VIRTUAL_ENV:-}" ] && [ "$entry" = "${VIRTUAL_ENV}/bin" ] && continue
    clean_path="${clean_path:+$clean_path:}$entry"
done
PATH="$clean_path"
export PATH
unset VIRTUAL_ENV

if [ -n "${SYNVOY_WEB_CONDA_ENV:-}" ]; then
    base="${SYNVOY_WEB_CONDA_BASE:-}"
    if [ -z "$base" ] && command -v conda >/dev/null 2>&1; then
        base="$(conda info --base 2>/dev/null)"
    fi
    if [ -z "$base" ] || [ ! -f "$base/etc/profile.d/conda.sh" ]; then
        echo "[synvoy-web] conda not found (conda_base='${base}'). Set conda_base in the settings."
        finish 96
    fi
    # shellcheck disable=SC1091
    source "$base/etc/profile.d/conda.sh"
    if ! conda activate "$SYNVOY_WEB_CONDA_ENV"; then
        echo "[synvoy-web] cannot activate the conda environment '$SYNVOY_WEB_CONDA_ENV'."
        finish 96
    fi
    if [ -z "${CONDA_PREFIX:-}" ] || [ ! -d "$CONDA_PREFIX/bin" ]; then
        echo "[synvoy-web] the conda environment '$SYNVOY_WEB_CONDA_ENV' has no bin folder."
        finish 96
    fi
    # conda puts the environment where the previous one stood, which can be behind
    # other folders. Move it to the front.
    rest=""
    IFS=':' read -r -a path_entries <<< "$PATH"
    for entry in "${path_entries[@]}"; do
        [ "$entry" = "$CONDA_PREFIX/bin" ] && continue
        rest="${rest:+$rest:}$entry"
    done
    PATH="$CONDA_PREFIX/bin${rest:+:$rest}"
    export PATH
fi

"$@"
finish $?
