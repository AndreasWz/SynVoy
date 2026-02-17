#!/bin/bash
set -e

# Snake Test Run
# Query: PLA2G2A (Phospholipase A2 Group IIA)
# Home: Homo sapiens (Reference)
# Targets: Naja naja, Pseudonaja textilis, Protobothrops mucrosquamatus, Python bivittatus, Bungarus multicinctus, Anolis carolinensis

# Ensure reference directory exists
mkdir -p ref

# Move Human genome if it's in genomes/ (post-download step)
if [ -f "genomes/Homo_sapiens.fna.gz" ]; then
    echo "Moving Homo_sapiens to ref/..."
    mv genomes/Homo_sapiens.fna.gz ref/
fi

echo "Running SynTerra Snake Test..."

# Run Nextflow
# Note: Using "pro" mode with explicit file paths
./nextflow run main.nf \
    --query_id "PLA2G2A" \
    --home_genome "ref/Homo_sapiens.fna.gz" \
    --home_species "Homo sapiens" \
    --target_genomes "genomes/*.fna.gz" \
    --mode pro \
    --outdir "results/snake_test" \
    --params.max_genomes 10 \
    -resume \
    -profile standard

echo "Snake Test Completed."
