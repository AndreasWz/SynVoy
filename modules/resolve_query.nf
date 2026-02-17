process RESOLVE_GENE_INPUT {
    tag "resolve"
    publishDir "${params.outdir}/query", mode: 'copy'
    
    input:
    val gene_input
    val species_override
    
    output:
    path "resolved_query/*.fasta", emit: fasta
    path "resolved_query/resolved_input.json", emit: metadata
    path "resolved_query/resolved_species.txt", emit: species
    
    script:
    def species_arg = species_override ? "--species \"${species_override}\"" : ""
    """
    resolve_gene_input.py \\
        --input "${gene_input}" \\
        --outdir resolved_query \\
        ${species_arg}
    """
}
