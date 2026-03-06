process PHYLO_SORT {
    tag "sorting"
    publishDir "${params.outdir}/intermediate/phylo_sort", mode: 'copy'
    
    input:
    tuple val(locus_id), path(home_genome)
    path genomes_dir // Directory instead of list
    val home_species_name // Species name or genome path for home identification
    
    output:
    tuple val(locus_id), path("sorted_genomes.txt"), emit: sorted_list
    
    script:
    """
    # Try phylogenetic sorting if taxonomy database available
    # Otherwise fall back to simple listing
    
    TAXDB=\${TAXDB:-NO_TAXDB}
    # Determine home query: prefer explicit species name, then fall back to genome file
    HOME_QUERY="${home_species_name}"
    if [ -z "\$HOME_QUERY" ] || [ "\$HOME_QUERY" = "null" ]; then
        HOME_QUERY="$home_genome"
    fi
    # If still a generic name, try to extract accession from FASTA headers
    if ! echo "\$HOME_QUERY" | grep -qE 'GC[AF]_[0-9]+'; then
        if [ -f "$home_genome" ]; then
            FASTA_ACC=\$(head -20 "$home_genome" | grep -oE 'GC[AF]_[0-9]+\\.[0-9]+' | head -1 || true)
            if [ -n "\$FASTA_ACC" ]; then
                HOME_QUERY="\${FASTA_ACC}.fna"
                echo "Extracted home accession from FASTA headers: \$HOME_QUERY"
            fi
        fi
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
