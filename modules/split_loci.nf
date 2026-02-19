process SPLIT_LOCI {
    tag "split"
    publishDir "${params.outdir}/intermediate/split_loci", mode: 'copy'
    
    input:
    path bed
    
    output:
    path "locus_*.bed", emit: beds
    
    script:
    """
    split_loci.py --bed $bed --output_prefix locus
    """
}
