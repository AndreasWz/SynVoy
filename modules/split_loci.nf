process SPLIT_LOCI {
    tag "split"
    publishDir "${params.outdir}/intermediate/split_loci", mode: 'copy'

    input:
    path bed
    path gff             // home annotation, or the NO_GFF placeholder
    path home_gene_bed   // goi_annotation.bed of ANNOTATE_GOI
    path home_gene_info  // goi_info.json of ANNOTATE_GOI
    path blast_hits
    path mmseqs_hits
    val home_locus_mode  // auto | gene_only | alignment

    output:
    path "locus_*.bed", emit: beds

    script:
    def gff_arg = gff.name != 'NO_GFF' ? "--gff ${gff}" : ''
    """
    split_loci.py --bed $bed --output_prefix locus \\
        ${gff_arg} \\
        --locus_gap_genes ${LocusGap.genes(params)} \\
        --locus_gap ${params.locus_gap} \\
        --home_gene_bed ${home_gene_bed} \\
        --home_gene_info ${home_gene_info} \\
        --home_locus_mode ${home_locus_mode} \\
        --hits ${blast_hits} \\
        --hits ${mmseqs_hits} \\
        --max_evalue ${params.search_evalue} \\
        --max_loci ${params.max_loci} \\
        --min_bit_ratio ${params.locus_min_bit_ratio} \\
        --family_warning_threshold ${params.family_warning_threshold}
    """
}
