process FILTER_SORTED_GENOMES {
    publishDir "${params.outdir}/intermediate/qc", mode: 'copy'

    input:
    tuple val(locus_id), path(sorted_list)
    path qc_json
    val qc_policy

    output:
    tuple val(locus_id), path("filtered_sorted_genomes.txt"), emit: sorted_list

    script:
    """
    filter_sorted_genomes.py \\
        --sorted $sorted_list \\
        --qc_json $qc_json \\
        --policy $qc_policy \\
        --output filtered_sorted_genomes.txt
    """
}
