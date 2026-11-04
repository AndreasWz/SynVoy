// The home gene's own neighbourhood and scaffold (bin/home_context.py), read by the
// preset modules (bin/select_modules.py): family members of the query among the genes
// next to it (the label's family test), and how far the gene lies from a scaffold end.
// Runs next to SPLIT_LOCI; only RESOLVE_EFFECTIVE_PARAMS waits for it. A failure stops the
// task (fix-phase d5): until 2026-11-02 it left no file and one line in the log, the resolver
// read a sentinel and the modules that need these fields fell back (A1 on the copy count of
// the home hits), so the settings depended on whether this step worked.
process HOME_CONTEXT {
    tag "home_context"
    label 'process_low'
    publishDir "${params.outdir}/intermediate/locate_gene", mode: 'copy'

    input:
    path query
    path home_genome
    path home_gff
    path goi_bed

    output:
    path "home_context.json", emit: json, optional: true

    script:
    """
    home_context.py \\
        --query ${query} \\
        --genome ${home_genome} \\
        --gff ${home_gff} \\
        --goi_bed ${goi_bed} \\
        --n_flanking_genes ${params.n_flanking_genes} \\
        --output home_context.json
    """
}
