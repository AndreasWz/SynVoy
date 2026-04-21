process ITERATIVE_SEARCH {
    tag "iterative"
    label 'process_long'

    input:
    tuple val(locus_id), path(initial_db)
    path sorted_genomes_file
    path genomes // List or Dir
    path home_db // MMseqs DB directory
    val flanking_genes
    val min_score
    val sensitivity

    output:
    tuple val(locus_id), path("iterative_results/expanded_db.faa"), emit: expanded_db
    tuple val(locus_id), path("iterative_results/hits"), emit: hits, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.faa"), emit: region_genes, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.gff"), emit: gff, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.homology.tsv"), emit: homology, optional: true


    script:
    """
    ${projectDir}/bin/iterative_search_runner.py \\
        --initial_db $initial_db \\
        --sorted_genomes $sorted_genomes_file \\
        --genomes_dir $genomes \\
        --home_db_dir $home_db \\
        --output_dir iterative_results \\
        --threads ${task.cpus} \\
        --mmseqs_sens ${sensitivity} \\
        --mmseqs_split_memory_limit ${params.mmseqs_split_memory_limit} \\
        --mmseqs_verbosity ${params.mmseqs_verbosity} \\
        --evalue ${params.search_evalue} \\
        --min_identity ${params.min_hit_identity} \\
        --min_length ${params.min_hit_length} \\
        --max_intron ${params.max_intron} \\
        --cluster_distance ${params.cluster_distance} \\
        --min_gene_identity ${params.min_gene_identity} \\
        --region_padding ${params.region_padding} \\
        --padding_min ${params.padding_min} \\
        --padding_max ${params.padding_max} \\
        --enable_smith_waterman ${params.enable_smith_waterman} \\
        --sw_method ${params.sw_method} \\
        --sw_min_score ${params.sw_min_score} \\
        --sw_min_identity ${params.sw_min_identity} \\
        --sw_timeout_seconds ${params.sw_timeout_seconds} \\
        --aug_relaxed_evalue_mult ${params.aug_relaxed_evalue_mult} \\
        --aug_relaxed_evalue_cap ${params.aug_relaxed_evalue_cap} \\
        --aug_relaxed_parse_evalue_mult ${params.aug_relaxed_parse_evalue_mult} \\
        --aug_relaxed_identity_factor ${params.aug_relaxed_identity_factor} \\
        --aug_relaxed_identity_min ${params.aug_relaxed_identity_min} \\
        --aug_relaxed_length_div ${params.aug_relaxed_length_div} \\
        --aug_relaxed_length_min ${params.aug_relaxed_length_min} \\
        --aug_dedup_bin_bp ${params.aug_dedup_bin_bp} \\
        --gap_search_window ${params.gap_search_window} \\
        --gap_min_size ${params.gap_min_size} \\
        --gap_evalue ${params.gap_evalue} \\
        --gap_min_identity ${params.gap_min_identity} \\
        --gap_min_alnlen ${params.gap_min_alnlen} \\
        --gap_max_hits ${params.gap_max_hits} \\
        --min_exon_query_cov ${params.min_exon_query_cov} \\
        --min_exon_alnlen ${params.min_exon_alnlen} \\
        --max_blocks_per_genome ${params.max_blocks_per_genome} \\
        --min_block_genes ${params.min_block_genes} \\
        --max_consecutive_empty_blocks ${params.max_consecutive_empty_blocks} \\
        --quiet_subtools ${params.iterative_quiet_subtools} \\
        --classify_high_min_identity ${params.classify_high_min_identity} \\
        --classify_medium_min_identity ${params.classify_medium_min_identity} \\
        --classify_tandem_min_identity ${params.classify_tandem_min_identity} \\
        --classify_fragment_max_qcov ${params.classify_fragment_max_qcov} \\
        --classify_complete_min_qcov ${params.classify_complete_min_qcov} \\
        --strict_goi_family ${params.strict_goi_family} \\
        --goi_family_tokens "${params.goi_family_tokens}" \\
        --gene_predictor ${params.gene_predictor} \\
        --augustus_species ${params.augustus_species} \\
        --enable_plm_search ${params.enable_plm_search} \\
        --plm_device ${params.plm_device} \\
        --plm_similarity_threshold ${params.plm_similarity_threshold} \\
        --plm_medium_threshold ${params.plm_medium_threshold} \\
        --plm_high_threshold ${params.plm_high_threshold} \\
        --enable_structural_search ${params.enable_structural_search} \\
        --structural_device ${params.structural_device} \\
        --structural_tm_threshold ${params.structural_tm_threshold} \\
        --structural_medium_threshold ${params.structural_medium_threshold} \\
        --structural_high_threshold ${params.structural_high_threshold} \\
        --structural_max_length ${params.structural_max_length} \\
        --prefix ${locus_id} \\
        --resume
    """
}
