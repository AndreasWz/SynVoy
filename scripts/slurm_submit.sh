#!/bin/bash
#SBATCH --job-name=synvoy
#SBATCH --output=synvoy_%j.log
#SBATCH --error=synvoy_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal
# Optional: #SBATCH --account=your_account

# ============================================================
# SynVoy SLURM Submission Script
# ============================================================
# This script submits the Nextflow controller job.
# Nextflow will then submit individual tasks to SLURM.
# ============================================================

# Load required modules (adjust for your cluster)
module load java/17 2>/dev/null || module load java/11 2>/dev/null || true
module load singularity 2>/dev/null || true
module load conda 2>/dev/null || true

# Ensure Nextflow is available
export PATH="${HOME}/bin:${PATH}"
export NXF_HOME="${HOME}/.nextflow"

# Work directory (use fast scratch if available)
export NXF_WORK="${SCRATCH:-${PWD}}/work"
mkdir -p "${NXF_WORK}"

# Singularity cache (for container images)
export NXF_SINGULARITY_CACHEDIR="${HOME}/.singularity/cache"
mkdir -p "${NXF_SINGULARITY_CACHEDIR}"

# ============================================================
# EDIT THESE PARAMETERS FOR YOUR RUN
# ============================================================

QUERY_ID="P01501"                    # UniProt ID (or switch to --mode pro with --query for local FASTA)
HOME_SPECIES="Apis mellifera"        # Species name for Easy Mode
MAX_GENOMES=10                       # Number of related genomes to fetch
OUTDIR="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

# Choose profile: hpc_singularity (recommended) or hpc_conda
PROFILE="hpc_singularity"

# ============================================================
# RUN PIPELINE
# ============================================================

echo "Starting SynVoy pipeline..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Query: ${QUERY_ID}"
echo "Species: ${HOME_SPECIES}"
echo "Output: ${OUTDIR}"
echo "Profile: ${PROFILE}"
echo "Work dir: ${NXF_WORK}"
echo "=============================================="

nextflow run main.nf \
    -profile ${PROFILE} \
    --query_id "${QUERY_ID}" \
    --home_species "${HOME_SPECIES}" \
    --max_genomes ${MAX_GENOMES} \
    --outdir "${OUTDIR}" \
    -work-dir "${NXF_WORK}" \
    -resume

EXIT_CODE=$?

echo "=============================================="
echo "Pipeline finished with exit code: ${EXIT_CODE}"
echo "Results: ${OUTDIR}"

# Optional: cleanup work directory on success
# if [ ${EXIT_CODE} -eq 0 ]; then
#     rm -rf "${NXF_WORK}"
# fi

exit ${EXIT_CODE}
