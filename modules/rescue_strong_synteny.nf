// Relaxed-miniprot rescue pass for goi_missing_but_strong_synteny blocks
// (docs/TODO.md §1e follow-up).
//
// CLUSTER_REGIONS classifies syntenic blocks; when one has strong flanking-gene
// support but no GOI model, cluster_grs.py emits region_class=goi_missing_but_strong_synteny
// in scores.tsv. This process runs miniprot inside those specific blocks with
// permissive thresholds (--outs 0.2 --outc 5%) so the user gets some model to
// look at instead of only an empty BED row. Output is staged into the per-locus
// region GFFs (Confidence=LOW, EvidenceType=relaxed_miniprot_rescue), so the
// existing report machinery picks it up via canonical_genome_id() without changes.
process RESCUE_STRONG_SYNTENY {
    tag "${locus_id}/${genome_name}"
    errorStrategy 'ignore'  // never block the pipeline on a single failed rescue

    input:
    tuple val(locus_id), val(genome_name), path(scores_tsv)
    path(genomes_dir)
    path(query_faa)

    output:
    tuple val(locus_id), path("${genome_name}.rescue.gff"), emit: gff

    when:
    !params.disable_strong_synteny_rescue

    script:
    """
    target_genome=\$(find -L ${genomes_dir} -name "${genome_name}*" -type f \\
        \\( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fna.gz" \\) \\
        | head -n 1)
    if [[ -z "\$target_genome" ]]; then
        echo "##gff-version 3" > ${genome_name}.rescue.gff
        echo "# RESCUE_STRONG_SYNTENY: target genome ${genome_name} not found in ${genomes_dir}" \\
            >> ${genome_name}.rescue.gff
        exit 0
    fi

    rescue_strong_synteny.py \\
        --scores ${scores_tsv} \\
        --target_genome "\$target_genome" \\
        --query ${query_faa} \\
        --genome_name "${genome_name}" \\
        --output ${genome_name}.rescue.gff \\
        --miniprot_outc ${params.strong_synteny_rescue_miniprot_outc} \\
        --window_pad ${params.strong_synteny_rescue_window_pad} \\
        --miniprot_timeout ${params.strong_synteny_rescue_miniprot_timeout}
    """
}
