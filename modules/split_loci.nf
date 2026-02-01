process SPLIT_LOCI {
    tag "split"
    
    input:
    path bed
    
    output:
    path "locus_*.bed", emit: beds
    
    script:
    """
    split_loci.py --bed $bed --output_prefix locus
    """
}
