process FETCH_RELATED_GENOMES {
    tag "easy_mode"
    publishDir "${params.outdir}/downloaded_genomes", mode: 'copy'
    
    input:
    val home_species
    val max_genomes
    
    output:
    path "easy_mode_genomes", emit: genomes_dir
    path "easy_mode_genomes/genomes_manifest.txt", emit: manifest
    path "easy_mode_genomes/species_mapping.tsv", emit: species_map
    
    script:
    """
    fetch_related_genomes.py \\
        --home-species "${home_species}" \\
        --max ${max_genomes} \\
        --outdir easy_mode_genomes
    """
}
