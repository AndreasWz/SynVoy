process FETCH_QUERY_FROM_ID {
    tag "uniprot"
    publishDir "${params.outdir}/query", mode: 'copy'
    
    input:
    val uniprot_id
    
    output:
    path "${uniprot_id}.fasta", emit: fasta
    
    script:
    """
    echo "Fetching sequence for ${uniprot_id} from UniProt..."
    curl -f -L "https://rest.uniprot.org/uniprotkb/${uniprot_id}.fasta" -o ${uniprot_id}.fasta || {
        echo "Failed to fetch from UniProt. Checking if ID exists..."
        exit 1
    }
    """
}
