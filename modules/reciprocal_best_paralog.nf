// Per-(locus, genome) reciprocal-best paralog check (docs/TODO.md §1j Phase B).
//
// Aligns each recovered GOI target protein against every home paralog in
// params.query. Emits a TSV that generate_report.py joins onto the GOI annotation
// list to detect ``paralog_confusion`` cases — calls whose best home-paralog
// match differs from the modal best at the same locus, indicating the MEDIUM-
// confidence iterative search dragged in a sibling paralog.
//
// Skipped when ``params.query`` is single-sequence (nothing to be reciprocal
// against) or when the kill switch ``--disable_paralog_check`` is set.
process RECIPROCAL_BEST_PARALOG {
    tag "${locus_id}/${genome_name}"
    errorStrategy 'ignore'  // one bad alignment must not block the whole report

    input:
    tuple val(locus_id), val(genome_name), path(target_faa), path(target_gff)
    path home_query

    output:
    tuple val(locus_id), path("${genome_name}.paralog_check.tsv"), emit: tsv

    when:
    !params.disable_paralog_check

    script:
    """
    reciprocal_best_paralog_check.py \\
        --home_query ${home_query} \\
        --target_faa ${target_faa} \\
        --target_gff ${target_gff} \\
        --locus_id "${locus_id}" \\
        --genome_name "${genome_name}" \\
        --output ${genome_name}.paralog_check.tsv
    """
}
