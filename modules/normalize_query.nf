process NORMALIZE_QUERY {
    tag "normalize_query"
    label 'process_low'
    publishDir "${params.outdir}/intermediate/query", mode: 'copy'

    input:
    path gene

    output:
    path "normalized_query.faa", emit: fasta

    script:
    // The home genome's genetic code, only when it is not the standard one (table 1), so a
    // default run's command line is unchanged (lib/GeneticCodes.groovy).
    def gc_arg = GeneticCodes.homeArg(params)
    """
    normalize_query.py --input $gene --output tmp_query_out.faa --min_length ${params.min_query_length} ${gc_arg}
    mv tmp_query_out.faa normalized_query.faa
    """
}
