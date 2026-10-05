// Phylogenetic placement check (docs/STATE_OF_THE_PROJECT.md F9).
//
// Asks whether each recovered GOI call's sequence divergence tracks its species'
// divergence. A lineage-specific paralog retained where the true ortholog was lost
// (bombolitin in Bombus) is anomalously diverged for how close its species is -- a
// signal that only exists ACROSS calls, which is why the per-call classifier cannot
// see it and why three earlier per-call fixes failed.
//
// Advisory: it emits verdicts, GENERATE_REPORT decides what to do with them.
process PHYLO_PLACEMENT_CHECK {
    tag "phylo_placement_${locus_id}"
    label 'process_low'
    publishDir "${params.outdir}/intermediate/phylo_placement", mode: 'copy', overwrite: true

    input:
    tuple val(locus_id), path(goi_faa)
    path query_faa
    path sorted_genomes

    output:
    tuple val(locus_id), path("${locus_id}.phylo_placement.tsv"), emit: tsv

    when:
    !params.disable_phylo_placement.toString().toBoolean()

    script:
    """
    phylo_placement_check.py \\
        --goi_faa ${goi_faa} \\
        --query ${query_faa} \\
        --sorted_genomes ${sorted_genomes} \\
        --min_calls ${params.phylo_placement_min_calls} \\
        --z_threshold ${params.phylo_placement_z_threshold} \\
        --disable_phylo_placement ${params.disable_phylo_placement} \\
        --output ${locus_id}.phylo_placement.tsv
    """
}
