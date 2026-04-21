process NORMALIZE_QUERY {
    tag "normalize_query"
    label 'process_low'
    publishDir "${params.outdir}/intermediate/query", mode: 'copy'

    input:
    path gene

    output:
    path "normalized_query.faa", emit: fasta

    script:
    """
    normalize_query.py --input $gene --output tmp_query_out.faa --min_length ${params.min_query_length}
    mv tmp_query_out.faa normalized_query.faa
    """
}
