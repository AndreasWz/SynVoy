#!/usr/bin/env bash
# Install TOGA1 + TOGA2 + the shared make_lastz_chains pipeline.
# Each TOGA version gets its own conda env (different Python pins, different
# Nextflow versions). Both versions consume the same chain files produced by
# make_lastz_chains, so chain generation in run_toga_tier3.sh is paid once and
# both TOGA1 and TOGA2 reuse the result.
#
# Compute on a 22 GB PC: 5-15 h for 7 Tier 3 chain pairs (lastz dominates).
# Install size estimate:
#   synvoy_toga  (TOGA1, py3.11)         ~2 GB
#   synvoy_toga2 (TOGA2, py3.12+ML deps) ~6 GB (tensorflow, xgboost, spliceai)
#
# Why both:
#   TOGA1 (Kirilenko et al., Science 2023, DOI: 10.1126/science.abn3107) is
#   citable and reviewer-verifiable. TOGA2 is the upstream-current version
#   under early access; we run it for completeness because the user asked,
#   even though the paper's headline TOGA comparison is the published TOGA1.
#
# Usage:
#   bash scripts/benchmark/toga_setup.sh
#     # installs both versions
#   SKIP_TOGA2=1 bash scripts/benchmark/toga_setup.sh
#     # only TOGA1 (fast, citable)
#   SKIP_TOGA1=1 bash scripts/benchmark/toga_setup.sh
#     # only TOGA2 (heavy, requires Py3.12)
#
# After install, paths are:
#   $HOME/dev/tools/make_lastz_chains/make_chains.py   (shared)
#   $HOME/dev/tools/TOGA/toga.py                       (TOGA1)
#   $HOME/dev/tools/TOGA2/toga2.py                     (TOGA2)

set -euo pipefail

ENV_TOGA1="${ENV_TOGA1:-${ENV:-synvoy_toga}}"
ENV_TOGA2="${ENV_TOGA2:-synvoy_toga2}"
TOOLS_DIR="${TOOLS_DIR:-${HOME}/dev/tools}"
SKIP_TOGA1="${SKIP_TOGA1:-0}"
SKIP_TOGA2="${SKIP_TOGA2:-0}"
mkdir -p "${TOOLS_DIR}"

if ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: 'mamba' not found on PATH. Install miniforge first:" >&2
    echo "  https://github.com/conda-forge/miniforge#install" >&2
    exit 1
fi

env_exists() {
    mamba env list 2>/dev/null | grep -qE "^\s*$1(\s|$)"
}

# ─────────────────────────────────────────────────────────────────────────
# TOGA1 (Kirilenko 2023) + make_lastz_chains
# ─────────────────────────────────────────────────────────────────────────
if [[ "${SKIP_TOGA1}" != "1" ]]; then
    echo "=== [TOGA1] creating conda env: ${ENV_TOGA1} ==="
    if env_exists "${ENV_TOGA1}"; then
        echo "  env ${ENV_TOGA1} already exists, skipping create"
    else
        # rust meta-package on conda-forge includes cargo (no separate `cargo` pkg).
        # Drop ucsc-chainantirepeat (not in any channel under that name).
        mamba create -n "${ENV_TOGA1}" -c bioconda -c conda-forge -y \
            python=3.11 \
            nextflow openjdk=17 \
            lastz \
            ucsc-fatotwobit ucsc-twobittofa \
            ucsc-axtchain ucsc-chainmergesort \
            ucsc-chaincleaner ucsc-chainsort ucsc-chainscore ucsc-chainnet \
            ucsc-chainfilter ucsc-pslsortacc ucsc-axttopsl \
            rust \
            2>&1 | tail -10
    fi

    echo
    echo "=== [TOGA1] installing uv into ${ENV_TOGA1} ==="
    mamba run -n "${ENV_TOGA1}" pip install --quiet uv 2>&1 | tail -3 || true

    echo
    echo "=== [shared] cloning make_lastz_chains ==="
    if [[ ! -d "${TOOLS_DIR}/make_lastz_chains" ]]; then
        git clone --depth 1 https://github.com/hillerlab/make_lastz_chains.git "${TOOLS_DIR}/make_lastz_chains"
    fi

    (cd "${TOOLS_DIR}/make_lastz_chains" && \
        mamba run -n "${ENV_TOGA1}" bash -c "uv venv && source .venv/bin/activate && uv pip install ."  2>&1 | tail -5 || true)

    echo
    echo "=== [TOGA1] cloning TOGA (Kirilenko et al., Science 2023) ==="
    if [[ ! -d "${TOOLS_DIR}/TOGA" ]]; then
        git clone --depth 1 https://github.com/hillerlab/TOGA.git "${TOOLS_DIR}/TOGA"
    fi

    echo "  installing TOGA1 Python deps into ${ENV_TOGA1}..."
    (cd "${TOOLS_DIR}/TOGA" && \
        mamba run -n "${ENV_TOGA1}" python3 -m pip install -r requirements.txt 2>&1 | tail -10 || true)

    # TOGA1 ships configure.sh which downloads CESAR2.0 + builds C utilities.
    if [[ -x "${TOOLS_DIR}/TOGA/configure.sh" ]]; then
        echo "  running TOGA1 configure.sh (downloads CESAR2.0, builds C utils)..."
        (cd "${TOOLS_DIR}/TOGA" && \
            mamba run -n "${ENV_TOGA1}" bash configure.sh 2>&1 | tail -20) || \
            echo "  (configure.sh non-fatal errors; verify with toga.py --help below)"
    fi
fi

# ─────────────────────────────────────────────────────────────────────────
# TOGA2 (upstream early access)
# ─────────────────────────────────────────────────────────────────────────
if [[ "${SKIP_TOGA2}" != "1" ]]; then
    echo
    echo "=== [TOGA2] creating conda env: ${ENV_TOGA2} ==="
    if env_exists "${ENV_TOGA2}"; then
        echo "  env ${ENV_TOGA2} already exists, skipping create"
    fi

    echo
    echo "=== [TOGA2] cloning hillerlab/TOGA2 ==="
    if [[ ! -d "${TOOLS_DIR}/TOGA2" ]]; then
        git clone --recurse-submodules https://github.com/hillerlab/TOGA2.git "${TOOLS_DIR}/TOGA2"
    else
        echo "  TOGA2 dir exists; updating submodules"
        (cd "${TOOLS_DIR}/TOGA2" && git submodule update --init --recursive 2>&1 | tail -5 || true)
    fi

    # TOGA2 supplies its own conda env spec; use it rather than rolling our own.
    if [[ ! -f "${TOOLS_DIR}/TOGA2/conda.yaml" ]]; then
        echo "  ERROR: TOGA2/conda.yaml missing; clone may be incomplete" >&2
    else
        if ! env_exists "${ENV_TOGA2}"; then
            # Use upstream's conda.yaml but force our env name (it's hardcoded
            # to 'toga2' in the file). Pipe the yaml through sed to swap names.
            tmp_yaml=$(mktemp --suffix=.yaml)
            sed "s/^name:.*/name: ${ENV_TOGA2}/" "${TOOLS_DIR}/TOGA2/conda.yaml" > "${tmp_yaml}"
            mamba env create -f "${tmp_yaml}" 2>&1 | tail -10
            rm -f "${tmp_yaml}"
        fi
    fi

    # TOGA2 uses `make` to build C/Rust/Cython parts, download CESAR, install
    # the Python deps from requirements.txt, and train/download ML models.
    echo "  running TOGA2 'make' (builds C/Cython/Rust, pip-installs deps,"
    echo "  downloads training models — expect 5-20 minutes)"
    (cd "${TOOLS_DIR}/TOGA2" && \
        mamba run -n "${ENV_TOGA2}" make 2>&1 | tail -30) || \
        echo "  (make exited non-zero; verify with toga2.py --help below)"
fi

# ─────────────────────────────────────────────────────────────────────────
# Verification
# ─────────────────────────────────────────────────────────────────────────
echo
echo "=== verification ==="
if [[ "${SKIP_TOGA1}" != "1" ]]; then
    mamba run -n "${ENV_TOGA1}" lastz --help 2>&1 | head -1 || echo "  lastz: ✗"
    mamba run -n "${ENV_TOGA1}" faToTwoBit 2>&1 | head -1 || echo "  faToTwoBit: ✓ (no-arg help)"
    mamba run -n "${ENV_TOGA1}" axtChain 2>&1 | head -1 || echo "  axtChain: ✓"
    echo "  make_chains.py:  $(ls ${TOOLS_DIR}/make_lastz_chains/make_chains.py 2>/dev/null && echo ✓ || echo ✗)"
    echo "  TOGA1 toga.py:   $(ls ${TOOLS_DIR}/TOGA/toga.py 2>/dev/null && echo ✓ || echo ✗)"
fi
if [[ "${SKIP_TOGA2}" != "1" ]]; then
    echo "  TOGA2 toga2.py:  $(ls ${TOOLS_DIR}/TOGA2/toga2.py 2>/dev/null && echo ✓ || echo ✗)"
    mamba run -n "${ENV_TOGA2}" python3 "${TOOLS_DIR}/TOGA2/toga2.py" --help 2>&1 | head -3 || \
        echo "  toga2.py --help failed; check 'make' output above"
fi

echo
echo "Done."
[[ "${SKIP_TOGA1}" != "1" ]] && echo "  TOGA1: mamba activate ${ENV_TOGA1}"
[[ "${SKIP_TOGA2}" != "1" ]] && echo "  TOGA2: mamba activate ${ENV_TOGA2}"
echo "Tools at: ${TOOLS_DIR}/{make_lastz_chains,TOGA,TOGA2}"
