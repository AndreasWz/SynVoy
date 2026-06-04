// §1n — Two-path home-locus establishment (name-lookup path).
// Resolves the GOI gene symbol against the home GFF to a SINGLE annotated locus
// (+ query↔gene consistency check), bypassing the noisy locate_gene/split_loci
// fan-out for annotated genes. Emits an EMPTY bed on any miss so main.nf falls
// back to the alignment-locate path. Gated in main.nf (enable_name_locus / home_goi_gene).
process RESOLVE_HOME_LOCUS {
    tag "name_locus"
    publishDir "${params.outdir}/intermediate/resolve_home_locus", mode: 'copy'

    input:
    path query
    val  gene_symbol
    path home_gff
    path home_genome

    output:
    // Named locus_1.bed so locus_id = "locus_1" matches split_loci's single-locus
    // convention downstream (EXTRACT_FLANKING / ITERATIVE_SEARCH / report locus_id).
    path "locus_1.bed",            emit: bed
    path "name_locus_status.json", emit: status

    script:
    def sym = (gene_symbol ?: '').toString().replaceAll('"', '').trim()
    """
    resolve_home_locus.py \\
        --query ${query} \\
        --gene_symbol "${sym}" \\
        --home_gff ${home_gff} \\
        --home_genome ${home_genome} \\
        --out_bed locus_1.bed \\
        --out_status name_locus_status.json \\
        --min_consistency_identity ${params.name_locus_min_identity} \\
        --pad ${params.name_locus_pad}
    """
}
