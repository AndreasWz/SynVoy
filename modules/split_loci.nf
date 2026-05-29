process SPLIT_LOCI {
    tag "split"
    publishDir "${params.outdir}/intermediate/split_loci", mode: 'copy'
    
    input:
    path bed
    
    output:
    path "locus_*.bed", emit: beds
    
    script:
    """
    split_loci.py --bed $bed --output_prefix locus \\
        --max_loci ${params.max_loci} \\
        --min_bit_ratio ${params.locus_min_bit_ratio} \\
        --family_warning_threshold ${params.family_warning_threshold}
    """
}
