process EXTRACT_FLANKING {
    tag "$bed"
    label 'process_medium'

    input:
    tuple val(locus_id), path(bed)
    path gff
    path genome
    val n_flank
    val min_size
    val prefer_large

    output:
    tuple val(locus_id), path("synteny_block.bed"), emit: bed
    tuple val(locus_id), path("flanking_proteins.faa"), emit: faa

    script:
    """
    # v3: GOI overlap re-injection + fallback pseudo-gene (force re-run)
    extract_flanking_genes.py \\
        --bed $bed \\
        --gff $gff \\
        --genome $genome \\
        --n_flank $n_flank \\
        --min_size $min_size \\
        --prefer_large $prefer_large \\
        --exon_mode true \\
        --out_bed synteny_block.bed \\
        --out_faa flanking_proteins.faa
    """
}
