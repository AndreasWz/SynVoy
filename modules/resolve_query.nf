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
    // If --gene is a local file, pass an absolute path so the resolver can
    // still find it from the process work directory.
    def resolved_input = gene_input
    try {
        def maybe_file = file(gene_input.toString())
        if (maybe_file.exists()) {
            resolved_input = maybe_file.toAbsolutePath().toString()
        }
    } catch (Exception ignored) {
        resolved_input = gene_input
    }

    def species_arg = species_override ? "--species \"${species_override}\"" : ""
    """
    resolve_gene_input.py \\
        --input "${resolved_input}" \\
        --outdir resolved_query \\
        ${species_arg}
    """
}
