process ANNOTATE_GOI {
    tag "annotate_goi"
    label 'process_medium'
    publishDir "${params.outdir}/intermediate/annotate_goi", mode: 'copy'
    
    input:
    path query_gene
    path home_genome
    path home_gff
    path blast_hits
    path mmseqs_hits
    val query_id
    val gene_name   // --home_goi_gene: the gene's name in the home annotation, or ''
    
    output:
    path "goi_exons.faa", emit: exons
    path "goi_annotation.bed", emit: bed
    path "goi_info.json", emit: info
    
    script:
    def gff_arg = home_gff.name != 'NO_GFF' ? home_gff : 'NO_GFF'
    def blast_arg = blast_hits.name != 'NO_BLAST_HITS' ? "--blast_hits ${blast_hits}" : ""
    def mmseqs_arg = mmseqs_hits.name != 'NO_MMSEQS_HITS' ? "--mmseqs_hits ${mmseqs_hits}" : ""
    def qid_arg = query_id ? "--query_id ${query_id}" : ""
    def name_arg = gene_name ? "--gene_name '${gene_name}'" : ""
    def gc_arg = GeneticCodes.homeArg(params)   // the home genome's code; '' for table 1
    """
    annotate_goi_exons.py \\
        --query ${query_gene} \\
        --genome ${home_genome} \\
        --gff ${gff_arg} \\
        ${blast_arg} \\
        ${mmseqs_arg} \\
        ${qid_arg} \\
        ${name_arg} \\
        ${gc_arg} \\
        --max_evalue ${params.search_evalue} \\
        --gff_search_window ${params.gff_search_window} \\
        --gap_search_window ${params.gap_search_window} \\
        --output_exons goi_exons.faa \\
        --output_bed goi_annotation.bed \\
        --output_info goi_info.json
    """
}
