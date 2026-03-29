process EXTRACT_FLANKING {
    tag "$bed"
    label 'process_medium'
    publishDir "${params.outdir}/intermediate/flanking", mode: 'copy'

    input:
    tuple val(locus_id), path(bed)
    path gff
    path genome
    val n_flank
    val min_size
    val prefer_large
    path goi_faa

    output:
    tuple val(locus_id), path("synteny_block_${locus_id}.bed"), emit: bed
    tuple val(locus_id), path("flanking_proteins_${locus_id}.faa"), emit: faa

    script:
    def goi_arg = (goi_faa && goi_faa.name != 'NO_GOI') ? "--goi_faa ${goi_faa} --max_goi_similarity ${params.max_flanking_goi_similarity}" : ""
    """
    # v4: GOI-similarity filter + expanded window
    extract_flanking_genes.py \\
        --bed $bed \\
        --gff $gff \\
        --genome $genome \\
        --n_flank $n_flank \\
        --min_size $min_size \\
        --prefer_large $prefer_large \\
        --exon_mode ${params.exon_level_search} \\
        --pred_flank_window ${params.pred_flank_window} \\
        --pred_keep_pct ${params.pred_keep_pct} \\
        $goi_arg \\
        --out_bed synteny_block_${locus_id}.bed \\
        --out_faa flanking_proteins_${locus_id}.faa
    """
}
