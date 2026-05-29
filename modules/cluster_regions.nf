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
    path species_map
    // §1f: preset-affected params resolved at runtime (see RESOLVE_EFFECTIVE_PARAMS).
    val settings

    output:
    tuple val(genome_name), val(payload), val(locus_id), path("regions/${genome_name}.regions.bed"), emit: bed
    // Carry locus_id on scores too so GENERATE_REPORT can stage per-locus copies under
    // distinct prefixes (docs/TODO.md §1h follow-up; previously the genome-only basename
    // collided across loci and under-counted region scores in multi-locus runs).
    tuple val(locus_id), val(genome_name), path("regions/${genome_name}.scores.tsv"), emit: scores

    script:
    def species_arg = species_map.name != 'NO_SPECIES_MAP' ? "--species_map ${species_map}" : ""
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
        --genome_name "${genome_name}" \\
        --output regions/${genome_name}.regions.bed \\
        --scores_output regions/${genome_name}.scores.tsv \\
        --min_score $min_score \\
        --cluster_distance ${params.cluster_distance} \\
        --weight_base ${params.synteny_weight_base} \\
        --weight_consistency ${params.synteny_weight_consistency} \\
        --weight_strand ${params.synteny_weight_strand} \\
        --goi_overlap_bonus ${params.synteny_goi_overlap_bonus} \\
        --max_regions ${params.max_regions} \\
        --adaptive_score_floor_frac ${settings.adaptive_score_floor_frac} \\
        --adaptive_score_floor_abs ${settings.adaptive_score_floor_abs} \\
        --adaptive_max_regions ${settings.adaptive_max_regions} \\
        --adaptive_unique_gene_floor ${settings.adaptive_unique_gene_floor} \\
        --strong_synteny_min_flanking ${params.strong_synteny_min_flanking} \\
        ${species_arg}
    """
}
