process RESOLVE_GENE_INPUT {
    tag "resolve"
    publishDir "${params.outdir}/query", mode: 'copy'
    
    input:
    val gene_input
    val species_override
    val inline_input
    
    output:
    path "resolved_query/*.fasta", emit: fasta
    path "resolved_query/resolved_input.json", emit: metadata
    path "resolved_query/resolved_species.txt", emit: species
    
    script:
    def inline_flag = inline_input ? "true" : "false"
    def inline_has_header = inline_input && gene_input.toString().trim().startsWith('>')
    def inline_header_flag = inline_has_header ? "true" : "false"

    // If --gene is a local file, pass an absolute path so the resolver can
    // still find it from the process work directory.
    def resolved_input = gene_input
    if (!inline_input) {
        try {
            def maybe_file = file(gene_input.toString())
            if (maybe_file.exists()) {
                resolved_input = maybe_file.toAbsolutePath().toString()
            }
        } catch (Exception ignored) {
            resolved_input = gene_input
        }
    }

    def species_arg = species_override ? "--species \"${species_override}\"" : ""
    """
    resolved_input_path="${resolved_input}"
    if [ "${inline_flag}" = "true" ]; then
        if [ "${inline_header_flag}" = "true" ]; then
            cat <<'EOF' > inline_query.fasta
${gene_input}
EOF
        else
            cat <<'EOF' > inline_query.fasta
>inline_query
${gene_input}
EOF
        fi
        resolved_input_path="inline_query.fasta"
    fi

    resolve_gene_input.py \\
        --input "\${resolved_input_path}" \\
        --outdir resolved_query \\
        ${species_arg}
    """
}
