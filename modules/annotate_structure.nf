process ANNOTATE_STRUCTURE {
    tag "$genome_name"
    
    input:
    tuple val(unique_id), val(genome_name), path(regions_bed), path(genomes_dir)
    val augustus_species

    output:
    tuple val(unique_id), path("annotations/${genome_name}.augustus.gff"), emit: gff
    tuple val(unique_id), path("annotations/${genome_name}.augustus.aa"), emit: proteins

    script:
    """
    mkdir -p annotations
    
    # Resolve target genome
    target_genome=\$(find -L $genomes_dir -name "${genome_name}*" -type f | head -n 1)
    
    if [ -z "\$target_genome" ]; then
         echo "Error: Could not find genome for ${genome_name}"
         exit 1
    fi
    
    # Extract region sequences first (with padding already in previous step, but let's be safe)
    # Actually, previous step (AUGMENTED_SEARCH) output candidates, but we want the full region for structure prediction.
    # The input here is regions_bed and genome.
    
    # We need to extract the FASTA for the regions to run Augustus on them efficiently,
    # OR run Augustus on the whole genome (slow) limited to regions.
    # Identifying regions first is better.
    
    # Extract regions using Python (biopython)
    python3 -c "
import sys
import os
from Bio import SeqIO

if os.stat('$regions_bed').st_size == 0:
    print('Empty regions file')
    sys.exit(0)

genome = SeqIO.to_dict(SeqIO.parse('\$target_genome', 'fasta'))
with open('$regions_bed') as f:
    for line in f:
        if not line.strip(): continue
        parts = line.strip().split('\t')
        if len(parts) < 4: continue
        chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
        if chrom in genome:
            seq = genome[chrom].seq[start:end]
            print(f'>{name}')
            print(seq)
    " > regions.fasta
    
    if [ -s regions.fasta ]; then
        # Run Augustus
        augustus --species=${augustus_species} --gff3=on regions.fasta > ${genome_name}.augustus.gff
        
        # Extract protein sequences
        if command -v getAnnoFasta.pl &> /dev/null; then
            getAnnoFasta.pl ${genome_name}.augustus.gff
            mv ${genome_name}.augustus.aa annotations/${genome_name}.augustus.aa || true
        else
           # Fallback
           grep -A 20 "^# protein sequence" ${genome_name}.augustus.gff | sed 's/# protein sequence = \\[//' | sed 's/\\]//' | sed 's/#//g' | grep -v "^--" > annotations/${genome_name}.augustus.aa
        fi
        
        # Move GFF
        mv ${genome_name}.augustus.gff annotations/${genome_name}.augustus.gff
    else
        echo "No regions to annotate."
        touch annotations/${genome_name}.augustus.gff
        touch annotations/${genome_name}.augustus.aa
    fi
    
    # Final safeguard
    touch annotations/${genome_name}.augustus.gff
    touch annotations/${genome_name}.augustus.aa
    """
}
