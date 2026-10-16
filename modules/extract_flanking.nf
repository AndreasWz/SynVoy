process EXTRACT_FLANKING {
    tag "$bed"
    label 'process_medium'
    publishDir "${params.outdir}/intermediate/flanking", mode: 'copy'

    input:
    tuple val(locus_id), path(bed)
    path gff
    path genome
    val n_flank
    val min_size
    path goi_faa
    // §1f: preset-affected params resolved at runtime (see RESOLVE_EFFECTIVE_PARAMS).
    val settings

    output:
    tuple val(locus_id), path("synteny_block_${locus_id}.bed"), emit: bed
    tuple val(locus_id), path("flanking_proteins_${locus_id}.faa"), emit: faa

    script:
    def goi_arg = (goi_faa && goi_faa.name != 'NO_GOI') ? "--goi_faa ${goi_faa} --goi_relative_evalue ${settings.goi_relative_evalue} --goi_relative_statistic ${settings.goi_relative_statistic}" : ""
    def dist_arg = params.max_flanking_distance > 0 ? "--max_flanking_distance ${params.max_flanking_distance}" : ""
    def expand_arg = settings.expand_goi_similar.toString().toBoolean() ? "--expand_goi_similar true --expand_goi_similar_distance ${params.expand_goi_similar_distance}" : ""
    def gc_arg = GeneticCodes.homeArg(params)   // the home genome's code; '' for table 1
    """
    # Relatives of the GOI are skipped as anchors (alignment E-value); the proteins carry
    # the anchors' home coordinates, from which the search derives its distances.
    extract_flanking_genes.py \\
        --bed $bed \\
        --gff $gff \\
        --genome $genome \\
        --n_flank $n_flank \\
        --min_size $min_size \\
        --exon_mode ${params.exon_level_search} \\
        --pred_flank_window ${params.pred_flank_window} \\
        --pred_keep_pct ${params.pred_keep_pct} \\
        --gene_predictor ${params.gene_predictor} \\
        --augustus_species ${params.augustus_species} \\
        $goi_arg \\
        $dist_arg \\
        $expand_arg \\
        $gc_arg \\
        --out_bed synteny_block_${locus_id}.bed \\
        --out_faa flanking_proteins_${locus_id}.faa
    """
}
