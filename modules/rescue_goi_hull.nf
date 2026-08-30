// GOI synteny-hull rescue (docs/TODO.md §1m fumble fix).
//
// SynVoy only models the GOI inside flanking-SEEDED sub-blocks. When the GOI sits in
// a GAP between flanking blocks (its immediate neighbours rearranged in the target),
// no search window covers it and the true ortholog is missed — even at near-identity.
// Real case: cow decorin at chr5:21.0 Mb fell between the seeded blocks (17.8-19.5 and
// 22.5-24 Mb), so the run recovered only paralog fragments (55 %) while miniprot finds
// the true 90 %+ decorin there. This pass runs miniprot of the GOI across the whole
// flanking neighbourhood (hull) per target chromosome, but only when the hull has strong
// flanking support and no HIGH GOI model already inside it. Output GFF (Confidence by
// identity, EvidenceType=synteny_hull_rescue) mixes into the per-locus region GFFs, so the
// report machinery picks it up via canonical_genome_id() without changes.
process RESCUE_GOI_HULL {
    tag "${locus_id}/${genome_name}"
    errorStrategy 'ignore'  // never block the pipeline on a single failed rescue

    input:
    // Same 4-tuple as paralog_inputs_ch (region_faa unused here) so main.nf needn't
    // build an extra channel — keeps the main workflow under the JVM method-size limit.
    tuple val(locus_id), val(genome_name), path(region_faa), path(region_gff)
    path(genomes_dir)
    path(query_faa)

    output:
    tuple val(locus_id), path("${genome_name}.hull_rescue.gff"), emit: gff
    // Protein of the rescued model, so the §1m ownership check can RBH it against the
    // home paralog panel. Rescue models used to bypass ownership entirely (they are
    // mixed in at STAGE_REGION_GFF, downstream of it), which is how a mislabelled
    // rescue call reached the headline unchecked.
    tuple val(locus_id), val(genome_name), path("${genome_name}.hull_rescue.faa"), emit: faa

    when:
    !params.disable_goi_hull_rescue

    script:
    """
    target_genome=\$(find -L ${genomes_dir} -name "${genome_name}*" -type f \\
        \\( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fna.gz" \\) \\
        | head -n 1)
    if [[ -z "\$target_genome" ]]; then
        echo "##gff-version 3" > ${genome_name}.hull_rescue.gff
        echo "# RESCUE_GOI_HULL: target genome ${genome_name} not found in ${genomes_dir}" \\
            >> ${genome_name}.hull_rescue.gff
        : > ${genome_name}.hull_rescue.faa
        exit 0
    fi

    rescue_goi_hull.py \\
        --region_gff ${region_gff} \\
        --target_genome "\$target_genome" \\
        --query ${query_faa} \\
        --genome_name ${genome_name} \\
        --output ${genome_name}.hull_rescue.gff \\
        --output_faa ${genome_name}.hull_rescue.faa \\
        --min_flanking ${params.goi_hull_min_flanking} \\
        --cluster_max_gap ${params.goi_hull_cluster_max_gap} \\
        --window_pad ${params.goi_hull_window_pad} \\
        --max_window ${params.goi_hull_max_window} \\
        --min_identity ${params.goi_hull_min_identity} \\
        --min_coverage ${params.goi_hull_min_coverage} \\
        --classify_high_min_identity ${params.classify_high_min_identity} \\
        --classify_medium_min_identity ${params.classify_medium_min_identity}
    """
}
