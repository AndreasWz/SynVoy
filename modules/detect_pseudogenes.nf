process DETECT_PSEUDOGENES {
    tag "${genome_name}"
    label 'process_low'
    publishDir "${params.outdir}/pseudogenes", mode: 'copy'

    input:
    tuple val(genome_name), path(gff), path(genome)
    path reference_faa

    output:
    tuple val(genome_name), path("${genome_name}.pseudogenes.tsv"), emit: report

    script:
    """
    detect_pseudogenes.py \\
        --gff $gff \\
        --reference $reference_faa \\
        --genome $genome \\
        --output ${genome_name}.pseudogenes.tsv \\
        --min_coverage 0.5 \\
        --min_identity 30.0
    """
}
