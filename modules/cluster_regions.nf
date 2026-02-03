process CLUSTER_REGIONS {
    tag "cluster_regions"
    label 'process_medium'
    errorStrategy 'ignore'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(genome_name), val(payload), path(hits_file), path(genomes_dir)
    tuple val(locus_id), path(synteny_bed)
    val flanking_count
    val min_score

    output:
    tuple val(genome_name), val(payload), val(locus_id), path("regions/${genome_name}.regions.bed"), emit: bed
    tuple val(genome_name), path("regions/${genome_name}.scores.tsv"), emit: scores

    script:
    """
    mkdir -p regions
    
    # Resolve genome file
    target_genome=\$(find -L $genomes_dir -name "${genome_name}*" -type f | head -n 1)
    
    cluster_grs.py \\
        --hits $hits_file \\
        --synteny_bed $synteny_bed \\
        --flanking_count $flanking_count \\
        --genome "\$target_genome" \\
        --output regions/${genome_name}.regions.bed \\
        --min_score $min_score
        
    # Create simple scores output
    grep "Region" regions/${genome_name}.regions.bed > regions/${genome_name}.scores.tsv || touch regions/${genome_name}.scores.tsv
    """
}
