process HOMOLOGY_SEARCH {
    tag "homology_search"
    // publishDir "${params.outdir}/homology", mode: 'copy'
    
    input:
    tuple val(genome_name), path(target_proteins)
    path(home_proteins)
    
    output:
    tuple val(genome_name), path("${genome_name}.homology.tsv"), emit: tsv
    
    script:
    """
    # Search target proteins (query) against home proteins (target/db)
    # We want to find the best match for each target protein in the home set.
    
    # Validating input
    if [ ! -s "${target_proteins}" ]; then
        echo "Input target proteins file is empty. No homology search performed."
        touch ${genome_name}.homology.tsv
        exit 0
    fi

    mmseqs easy-search \\
        ${target_proteins} \\
        ${home_proteins} \\
        ${genome_name}.homology.m8 \\
        tmp_homology \\
        --format-output "query,target,evalue,pident" \\
        -e 1e-3
        
    # Convert to simple TSV: TargetGene \t HomeGene
    # Take best hit only (sort by evalue (col 3), then distinct query (col 1))
    
    if [ -s "${genome_name}.homology.m8" ]; then
        sort -k1,1 -k3,3g ${genome_name}.homology.m8 | sort -u -k1,1 > best_hits.m8
        # Extract Gene Name from Transcript ID (remove .t*)
        # Input: g1.t1 \t P12345
        # Output: g1 \t P12345
        cut -f1,2 best_hits.m8 | sed 's/\\.t[0-9]*\\t/\\t/' > ${genome_name}.homology.tsv
    else
        touch ${genome_name}.homology.tsv
    fi
    """
}
