process STAGE_GENOMES {
    tag "staging"
    
    input:
    path genomes
    
    output:
    path "staged_genomes", emit: dir
    
    script:
    """
    mkdir staged_genomes
    # Copy files into directory
    # Use cp -r to handle potential directories or multiple files
    cp -f -L -r $genomes staged_genomes/
    """
}
