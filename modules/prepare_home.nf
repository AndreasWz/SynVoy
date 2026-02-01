process PREPARE_HOME_PROTEOME {
    tag "home_prot"
    
    input:
    path home_genome
    
    output:
    path "home_proteome.faa", emit: faa
    path "home_proteome_db", emit: db
    
    script:
    """
    # Predict proteins
    prodigal -i $home_genome -a home_proteome.faa -p meta -q
    
    # Create MMseqs DB for fast searching
    mkdir home_proteome_db
    mmseqs createdb home_proteome.faa home_proteome_db/db
    """
}
