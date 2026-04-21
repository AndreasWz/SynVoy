process CLUSTER_REGIONS {
    tag "cluster_regions"
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(genome_name), val(payload), path(hits_file), path(genomes_dir), path(target_gff)
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
        --target_gff $target_gff \\
        --synteny_bed $synteny_bed \\
        --flanking_count $flanking_count \\
        --genome "\$target_genome" \\
        --output regions/${genome_name}.regions.bed \\
        --scores_output regions/${genome_name}.scores.tsv \\
        --min_score $min_score \\
        --cluster_distance ${params.cluster_distance} \\
        --weight_base ${params.synteny_weight_base} \\
        --weight_consistency ${params.synteny_weight_consistency} \\
        --weight_strand ${params.synteny_weight_strand} \\
        --goi_overlap_bonus ${params.synteny_goi_overlap_bonus} \\
        --max_regions ${params.max_regions} \\
        --adaptive_score_floor_frac ${params.adaptive_score_floor_frac} \\
        --adaptive_score_floor_abs ${params.adaptive_score_floor_abs} \\
        --adaptive_max_regions ${params.adaptive_max_regions} \\
        --adaptive_unique_gene_floor ${params.adaptive_unique_gene_floor}
    """
}
