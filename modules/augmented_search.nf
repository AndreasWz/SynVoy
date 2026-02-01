process AUGMENTED_SEARCH {
    tag "augmented_search"
    label 'process_high'

    input:
    tuple val(unique_id), val(genome_name), path(regions_bed), path(genomes_dir)
    path query_gene
    val padding

    output:
    tuple val(unique_id), path("augmented/${genome_name}.candidates.fna"), emit: proteins
    tuple val(unique_id), path("augmented/${genome_name}.candidates.bed"), emit: bed

    script:
    """
    mkdir -p augmented
    
    # Resolve target genome from directory
    # Assume genome_name matches filename prefix or exact match
    target_genome=\$(find -L $genomes_dir -name "${genome_name}*" -type f | head -n 1)
    
    if [ -z "\$target_genome" ]; then
        echo "Error: Could not find genome for ${genome_name} in $genomes_dir"
        exit 1
    fi
    
    echo "Using target genome: \$target_genome"
    
    augmented_search_runner.py \\
        --regions_bed $regions_bed \\
        --target_genome \$target_genome \\
        --query_gene $query_gene \\
        --output_base augmented/${genome_name}.candidates \\
        --padding $padding \\
        --mmseqs_sens "8.5"
    """
}
