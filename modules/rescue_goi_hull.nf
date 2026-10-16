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
// the search's judge since 2026-10-14, EvidenceType=synteny_hull_rescue) mixes into the per-locus region GFFs, so the
// report machinery picks it up via canonical_genome_id() without changes.
process RESCUE_GOI_HULL {
    tag "${locus_id}/${genome_name}"
    // Never block the pipeline on a single failed rescue, but retry a task the scheduler
    // killed for memory or time: 'ignore' alone overrode the config's retry for those.
    errorStrategy { task.exitStatus in [137, 140, 143] && task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    // Publish the rescued models. Until 2026-10-05 nothing did: a hull-rescued ortholog
    // appeared in synvoy_report.json as a span only, with its exon coordinates and
    // protein left behind in the work directory. One folder per home locus (the file
    // name carries only the genome, so a flat folder would let loci overwrite each
    // other); header-only files from a rescue that did not fire are not published.
    publishDir path: { "${params.outdir}/rescue/${locus_id}" }, mode: 'copy',
               saveAs: { fn -> task.workDir.resolve(fn).toFile().readLines().any { l -> l && !l.startsWith('#') } ? fn : null }

    input:
    // Same 4-tuple as paralog_inputs_ch (region_faa unused here) so main.nf needn't
    // build an extra channel — keeps the main workflow under the JVM method-size limit.
    tuple val(locus_id), val(genome_name), path(region_faa), path(region_gff)
    path(genomes_dir)
    path(query_faa)
    // The home protein database folder (home_proteome.faa + home_genes.tsv): a rescue model
    // is judged like a candidate of the main search (bin/synvoy_label.py), with its
    // specificity against the other home proteins and the GOI's home neighbourhood.
    path(home_db)
    // goi_info.json of ANNOTATE_GOI: the home gene the home step took as the query's is its
    // own gene in the specificity terms, as in the main search (review fix RF-2b).
    path(goi_info)

    output:
    tuple val(locus_id), path("${genome_name}.hull_rescue.gff"), emit: gff
    // Protein of the rescued model, so the §1m ownership check can RBH it against the
    // home paralog panel. Rescue models used to bypass ownership entirely (they are
    // mixed in at STAGE_REGION_GFF, downstream of it), which is how a mislabelled
    // rescue call reached the headline unchecked.
    tuple val(locus_id), val(genome_name), path("${genome_name}.hull_rescue.faa"), emit: faa

    when:
    !params.disable_goi_hull_rescue.toString().toBoolean()

    script:
    """
    # The exact file of this genome: matching on the name as a prefix also found
    # Apis_cerana_japonica.fna for Apis_cerana (bug hunt 2026-10-13, A11).
    target_genome=\$(find -L ${genomes_dir} -type f \\( -name "${genome_name}" -o -name "${genome_name}.fa" \\
        -o -name "${genome_name}.fna" -o -name "${genome_name}.fasta" -o -name "${genome_name}.fa.gz" \\
        -o -name "${genome_name}.fna.gz" -o -name "${genome_name}.fasta.gz" \\) | sort | head -n 1)
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
        --n_flanking_genes ${params.n_flanking_genes} \\
        --home_db_dir ${home_db} \\
        --goi_info ${goi_info} ${GeneticCodes.targetArgs(params)}
    """
}
