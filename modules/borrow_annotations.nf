process BORROW_ANNOTATIONS {
    tag "borrow_annotations"
    label 'process_medium'
    
    input:
    path home_genome
    path home_proteins  // Prodigal-predicted proteins
    path genomes_dir    // Directory with target genomes + optional GFFs
    path goi_bed        // GOI location from LOCATE_GENE
    val n_flanking
    
    output:
    path "borrowed_annotations.gff", emit: gff
    path "borrowed_proteins.faa", emit: proteins
    
    script:
    """
    borrow_annotations.py \\
        --home_genome $home_genome \\
        --home_proteins $home_proteins \\
        --genomes_dir $genomes_dir \\
        --goi_bed $goi_bed \\
        --output_gff borrowed_annotations.gff \\
        --output_proteins borrowed_proteins.faa \\
        --n_flanking $n_flanking
    """
}
