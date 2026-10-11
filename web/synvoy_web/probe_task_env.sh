#!/bin/bash
# probe_task_env.sh <conda environment folder>
#
# Switches to a conda environment in exactly the way a Nextflow task of the pipeline
# does, and checks that python3 is then that environment's and can load parasail (the
# package whose absence silently weakens the search; the pipeline stops on it).
# Prints the python3 in use. Used by the worker's start-up check.
env_dir="$1"
base="$(conda info --json 2>/dev/null | awk '/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }')"
if [ -z "$base" ] || [ ! -f "$base/bin/activate" ]; then
    echo "conda does not report its base folder"
    exit 1
fi
# shellcheck disable=SC1091
source "$base/bin/activate" "$env_dir" || { echo "cannot switch to $env_dir"; exit 1; }
python3 -c "import parasail, sys; print(sys.executable)"
