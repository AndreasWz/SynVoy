process FETCH_HOME_GENOME {
    tag "${species}"
    publishDir "${params.outdir}/home_genome", mode: 'copy'
    
    input:
    val species
    
    output:
    path "home_genome/home_genome.fna", emit: genome
    path "home_genome/home_genome.gff", emit: gff, optional: true
    
    script:
    def ranking_arg = "--assembly-ranking \"${params.assembly_ranking}\""
    def bad_policy_arg = "--bad-quality-policy \"${params.bad_quality_policy}\""
    def bad_timeout_arg = "--bad-quality-timeout ${params.bad_quality_timeout}"
    def bad_contigs_arg = "--bad-max-contigs ${params.bad_max_contigs}"
    def bad_scaffolds_arg = "--bad-max-scaffolds ${params.bad_max_scaffolds}"
    def bad_n50_arg = "--bad-min-n50 ${params.bad_min_n50}"
    """
    fetch_home_genome.py \\
        --species "${species}" \\
        --outdir home_genome \\
        ${ranking_arg} \\
        ${bad_policy_arg} \\
        ${bad_timeout_arg} \\
        ${bad_contigs_arg} \\
        ${bad_scaffolds_arg} \\
        ${bad_n50_arg}
    """
}
