process PHYLO_SORT {
    tag "sorting"
    publishDir "${params.outdir}/intermediate/phylo_sort", mode: 'copy'
    
    input:
    tuple val(locus_id), path(home_genome)
    path genomes_dir // Directory instead of list
    
    output:
    tuple val(locus_id), path("sorted_genomes.txt"), emit: sorted_list
    
    script:
    """
    # Try phylogenetic sorting if taxonomy database available
    # Otherwise fall back to simple listing
    
    TAXDB=\${TAXDB:-NO_TAXDB}
    HOME_NAME="${params.home_species ?: ''}"
    if [ -n "\$HOME_NAME" ]; then
        HOME_QUERY="\$HOME_NAME"
    else
        HOME_QUERY="$home_genome"
    fi
    
    if [ -d "\$TAXDB" ] && [ -f "\$TAXDB/nodes.dmp" ]; then
        echo "Attempting phylogenetic sorting using TaxDB: \$TAXDB"
        phylo_sort.py \\
            --home "\$HOME_QUERY" \\
            --targets_dir $genomes_dir \\
            --taxdb \$TAXDB \\
            --output sorted_genomes.txt || {
                echo "Phylogenetic sorting failed, falling back to alphabetical"
                phylo_sort.py --home "\$HOME_QUERY" --targets_dir $genomes_dir --taxdb "NO_DB" --output sorted_genomes.txt
            }
    else
        echo "No TaxDB found (set TAXDB environment variable)"
        echo "Using alphabetical order (no phylogenetic sorting)"
        # Use phylo_sort's fallback logic
        phylo_sort.py --home "\$HOME_QUERY" --targets_dir $genomes_dir --taxdb "NO_DB" --output sorted_genomes.txt
    fi
    """
}
