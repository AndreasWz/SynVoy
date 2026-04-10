#!/bin/bash
CONDA_BASE=/home/faw/miniforge3
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate syntenyfinder

JAVA_CMD=/home/faw/miniforge3/envs/syntenyfinder/bin/java \
nextflow run main.nf \
  --query_id "Q16553" \
  --home_species "Homo sapiens" \
  --target_species "Gallus gallus,Naja naja,Mus musculus,Canis lupus familiaris" \
  --outdir "local_runs/results_high_quality_test" \
  -profile laptop_safe \
  -resume
