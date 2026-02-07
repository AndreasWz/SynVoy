process FETCH_HOME_GENOME {
    tag "${species}"
    publishDir "${params.outdir}/home_genome", mode: 'copy'
    
    input:
    val species
    
    output:
    path "home_genome/home_genome.fna", emit: genome
    path "home_genome/home_genome.gff", emit: gff, optional: true
    
    script:
    """
    fetch_home_genome.py \\
        --species "${species}" \\
        --outdir home_genome
    """
}
