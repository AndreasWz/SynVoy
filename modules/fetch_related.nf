process FETCH_RELATED_GENOMES {
    tag "easy_mode"
    publishDir "${params.outdir}/downloaded_genomes", mode: 'copy'
    
    input:
    val home_species
    val max_genomes
    val target_species
    
    output:
    path "easy_mode_genomes", emit: genomes_dir
    path "easy_mode_genomes/genomes_manifest.txt", emit: manifest
    path "easy_mode_genomes/species_mapping.tsv", emit: species_map
    
    script:
    def target_arg = target_species ? "--target-species \"${target_species}\"" : ""
    def ranking_arg = "--assembly-ranking \"${params.assembly_ranking}\""
    def bad_policy_arg = "--bad-quality-policy \"${params.bad_quality_policy}\""
    def bad_timeout_arg = "--bad-quality-timeout ${params.bad_quality_timeout}"
    def bad_contigs_arg = "--bad-max-contigs ${params.bad_max_contigs}"
    def bad_scaffolds_arg = "--bad-max-scaffolds ${params.bad_max_scaffolds}"
    def bad_n50_arg = "--bad-min-n50 ${params.bad_min_n50}"
    """
    fetch_related_genomes.py \\
        --home-species "${home_species}" \\
        --max ${max_genomes} \\
        --outdir easy_mode_genomes \\
        ${ranking_arg} \\
        ${bad_policy_arg} \\
        ${bad_timeout_arg} \\
        ${bad_contigs_arg} \\
        ${bad_scaffolds_arg} \\
        ${bad_n50_arg} \\
        ${target_arg}
    """
}
