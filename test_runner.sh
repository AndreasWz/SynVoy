#!/bin/bash
set -e

# Create dummy input
echo ">gene1" > test_db.faa
echo "MKAILIVGLSLWKS" >> test_db.faa

echo ">chrom1" > test_genome.fna
# Create random sequence
# Sequence generation loop follows
# Generate random sequence using python
python3 -c "import random; print(''.join(random.choices('ACGT', k=2000)) + 'ATGAAGGCCATCCTGATCGTGGGCCTGAGCCTGTGGAAGAGC' + ''.join(random.choices('ACGT', k=2000)))" >> test_genome.fna

echo "$(pwd)/test_genome.fna" > sorted.txt

# Create dummy home db
mkdir -p test_home_db
mmseqs createdb test_db.faa test_home_db/db

# Run
python3 bin/iterative_search_runner.py \
    --initial_db test_db.faa \
    --sorted_genomes sorted.txt \
    --output_dir test_out \
    --home_db_dir test_home_db \
    --min_length 10 \
    --evalue 10 \
    --threads 1

echo "Run finished"
if [ -f test_out/hits/test_genome.fna.m8 ]; then
    echo "Hits file created (probably empty)"
else
    echo "Hits file missing"
fi
