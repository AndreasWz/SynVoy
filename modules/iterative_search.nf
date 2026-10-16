process ITERATIVE_SEARCH {
    tag "iterative"
    label 'process_long'
    // Decoy calls (--decoy_mode) are published as they are written; nothing else reads them.
    publishDir path: { "${params.outdir}/decoys/${locus_id}" }, mode: 'copy',
        pattern: 'iterative_results/decoys/*', saveAs: { fn -> file(fn).name }

    input:
    tuple val(locus_id), path(initial_db)
    path sorted_genomes_file
    path genomes // List or Dir
    path home_db // MMseqs DB directory
    val sensitivity
    // §1f: preset-affected params resolved at runtime (see RESOLVE_EFFECTIVE_PARAMS).
    // All other params still read directly from `params.X`.
    val settings
    // Home GOI structure (goi_info.json) → data-driven GOI fallback max-intron.
    path goi_info

    output:
    tuple val(locus_id), path("iterative_results/expanded_db.faa"), emit: expanded_db
    tuple val(locus_id), path("iterative_results/goi_for_tree.faa"), emit: goi_for_tree, optional: true
    tuple val(locus_id), path("iterative_results/hits"), emit: hits, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.faa"), emit: region_genes, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.gff"), emit: gff, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.homology.tsv"), emit: homology, optional: true
    tuple val(locus_id), path("iterative_results/decoys/*"), emit: decoys, optional: true


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
        --deterministic_goi_search ${params.deterministic_goi_search} \\
        --evalue ${params.search_evalue} \\
        --min_identity ${settings.min_hit_identity} \\
        --min_length ${settings.min_hit_length} \\
        --max_intron ${params.max_intron} \\
        --goi_info $goi_info \\
        --goi_fallback_intron_margin ${params.goi_fallback_intron_margin} \\
        --goi_fallback_intron_floor ${params.goi_fallback_intron_floor} \\
        --block_gap_genes ${params.block_gap_genes} \\
        --cluster_distance ${params.cluster_distance} \\
        --window_margin_genes ${params.window_margin_genes} \\
        --region_padding ${params.region_padding} \\
        --goi_margin_genes ${params.goi_margin_genes} \\
        --goi_margin ${params.goi_margin} \\
        --padding_min ${params.padding_min} \\
        --padding_max ${params.padding_max} \\
        --enable_smith_waterman ${params.enable_smith_waterman} \\
        --sw_method ${params.sw_method} \\
        --allow_missing_smith_waterman ${params.allow_missing_smith_waterman} \\
        --sw_min_score ${settings.sw_min_score} \\
        --sw_min_identity ${settings.sw_min_identity} \\
        --sw_timeout_seconds ${params.sw_timeout_seconds} \\
        --aug_relaxed_evalue_mult ${params.aug_relaxed_evalue_mult} \\
        --aug_relaxed_evalue_cap ${params.aug_relaxed_evalue_cap} \\
        --aug_relaxed_parse_evalue_mult ${params.aug_relaxed_parse_evalue_mult} \\
        --aug_relaxed_identity_factor ${params.aug_relaxed_identity_factor} \\
        --aug_relaxed_identity_min ${params.aug_relaxed_identity_min} \\
        --aug_relaxed_length_div ${params.aug_relaxed_length_div} \\
        --aug_relaxed_length_min ${params.aug_relaxed_length_min} \\
        --gap_search_window ${params.gap_search_window} \\
        --gap_min_identity ${params.gap_min_identity} \\
        --gap_max_hits ${params.gap_max_hits} \\
        --min_exon_alnlen ${params.min_exon_alnlen} \\
        --max_blocks_per_genome ${params.max_blocks_per_genome} \\
        --min_block_genes ${params.min_block_genes} \\
        --disable_synteny_collinearity ${params.disable_synteny_collinearity} \\
        --synteny_bridge_max_gap_genes ${params.synteny_bridge_max_gap_genes} \\
        --synteny_bridge_max_gap ${params.synteny_bridge_max_gap} \\
        --disable_anchor_resolution ${params.disable_anchor_resolution} \\
        --anchor_max_places ${params.anchor_max_places} \\
        --disable_relative_query ${params.disable_relative_query} \\
        --goi_max_models_per_locus ${params.goi_max_models_per_locus} \\
        --synteny_bridge_max_rank_gap ${params.synteny_bridge_max_rank_gap} \\
        --synteny_bridge_min_anchors ${params.synteny_bridge_min_anchors} \\
        --synteny_bridge_two_sided ${params.synteny_bridge_two_sided} \\
        --synteny_bridge_max_per_block ${params.synteny_bridge_max_per_block} \\
        --fallback_short_query_len ${params.fallback_short_query_len} \\
        --fallback_short_min_aln_aa ${params.fallback_short_min_aln_aa} \\
        --fallback_short_min_bits ${params.fallback_short_min_bits} \\
        --quiet_subtools ${params.iterative_quiet_subtools} \\
        --decoy_mode ${params.decoy_mode} \\
        --report_nonsyntenic_candidates ${params.report_nonsyntenic_candidates} \\
        --rank_wave_binning ${params.rank_wave_binning} \\
        --disable_wavefront ${params.disable_wavefront} \\
        --classify_fragment_max_qcov ${settings.classify_fragment_max_qcov} \\
        --classify_complete_min_qcov ${settings.classify_complete_min_qcov} \\
        --goi_family_tokens "${settings.goi_family_tokens}" \\
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
        --enable_structural_discovery ${params.enable_structural_discovery} \\
        --discovery_min_tm ${params.discovery_min_tm} \\
        --discovery_min_len_ratio ${params.discovery_min_len_ratio} \\
        --discovery_max_len_ratio ${params.discovery_max_len_ratio} \\
        --discovery_min_goi_coverage ${params.discovery_min_goi_coverage} \\
        --discovery_paralog_margin ${params.discovery_paralog_margin} \\
        --discovery_min_block_flanking ${params.discovery_min_block_flanking} \\
        --discovery_allow_offblock ${params.discovery_allow_offblock} \\
        --enable_dispersed_goi_rescue ${params.enable_dispersed_goi_rescue} \\
        --dispersed_goi_min_identity ${params.dispersed_goi_min_identity} \\
        --dispersed_goi_max_evalue ${params.dispersed_goi_max_evalue} \\
        --dispersed_goi_min_alnlen ${params.dispersed_goi_min_alnlen} \\
        --dispersed_min_chrom_anchors ${params.dispersed_min_chrom_anchors} \\
        --dispersed_envelope_margin_genes ${params.dispersed_envelope_margin_genes} \\
        --dispersed_envelope_margin ${params.dispersed_envelope_margin} \\
        --prefix ${locus_id} \\
        --resume
    """
}
