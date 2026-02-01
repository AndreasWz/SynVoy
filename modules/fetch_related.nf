process FETCH_RELATED_GENOMES {
    tag "easy_mode"
    publishDir "${params.outdir}/downloaded_genomes", mode: 'copy'
    
    input:
    val species_name
    val max_genomes
    
    output:
    path "easy_mode_genomes", emit: genomes_dir
    path "easy_mode_genomes/genomes_manifest.txt", emit: manifest
    
    script:
    """
    fetch_related_genomes.py \\
        --species "${species_name}" \\
        --max ${max_genomes} \\
        --outdir easy_mode_genomes
    """
}
