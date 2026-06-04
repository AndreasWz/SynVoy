// Locus ownership via home-paralog reciprocal-best (docs/TODO.md §1m).
//
// Counterpart to the §1m diagnosis: home loci search independently, and because
// flanking genes cross-match between sister paralogs, the WRONG home locus can
// model a gene (the chr9 OMD locus modelled real decorin on cow BTA5, while the
// chr12 decorin locus produced nothing there). This module assigns each recovered
// GOI target to the home locus it TRULY belongs to:
//
//   BUILD_HOME_PARALOG_PANEL — one protein per home GOI-paralog, labelled
//       ``>locus_<id>|<gene>`` (the home gene overlapping each split locus).
//   ASSIGN_LOCUS_OWNERSHIP   — reuses reciprocal_best_paralog_check.py to SW-align
//       each recovered GOI target against the panel; best panel entry = owning locus.
//
// generate_report.py joins the ownership TSVs onto the deduped GOI records to
// re-attribute mis-filed orthologs and relabel paralogs carried under the GOI name.
// Disable with --disable_locus_ownership.

process BUILD_HOME_PARALOG_PANEL {
    tag "home_paralog_panel"
    errorStrategy 'ignore'  // a missing/odd home GFF must not block the run

    input:
    path home_gff
    path home_proteome
    path locus_beds
    path query

    output:
    path "home_paralog_panel.faa", emit: faa
    path "home_paralog_panel.faa.meta.tsv", emit: meta

    when:
    !params.disable_locus_ownership

    script:
    def query_arg = query.name != 'NO_QUERY' ? "--query ${query}" : ""
    """
    build_home_paralog_panel.py \\
        --home_gff ${home_gff} \\
        --home_proteome ${home_proteome} \\
        --locus_beds ${locus_beds} \\
        ${query_arg} \\
        --pad ${params.locus_ownership_pad} \\
        --max_genes_per_locus ${params.locus_ownership_max_genes_per_locus} \\
        --paralog_min_identity ${params.paralog_panel_min_identity} \\
        --max_family_paralogs ${params.paralog_panel_max_family} \\
        --output home_paralog_panel.faa

    # Guarantee both outputs exist even if the panel came out empty, so the
    # collect() downstream never hangs waiting on an optional file.
    touch home_paralog_panel.faa home_paralog_panel.faa.meta.tsv
    """
}

process ASSIGN_LOCUS_OWNERSHIP {
    tag "${locus_id}/${genome_name}"
    errorStrategy 'ignore'  // one bad alignment must not block the whole report

    input:
    tuple val(locus_id), val(genome_name), path(target_faa), path(target_gff)
    // NOTE: the caller MUST pass the panel as a value channel (`.first()`). The single
    // panel has to broadcast to EVERY (locus, genome) target; as a plain queue channel it
    // pairs positionally with the per-genome tuple and only the first genome runs, so the
    // other genomes' paralog calls are never RBH-checked / demoted (the §16 fan-out bug).
    path panel

    output:
    // Prefix with locus_id: the same genome is processed once per home locus, so
    // a bare ${genome_name}.locus_ownership.tsv would collide when GENERATE_REPORT
    // flat-copies all loci into one staging dir (same family as docs/TODO.md §1h).
    // generate_report.py keys ownership rows by the in-file `locus`/`genome`
    // columns, so the prefix only guarantees filename uniqueness.
    tuple val(locus_id), path("${locus_id}__${genome_name}.locus_ownership.tsv"), emit: tsv

    when:
    !params.disable_locus_ownership

    script:
    """
    # Reuse the §1j reciprocal-best machinery, but the "paralog panel" is the
    # home-locus panel (>locus_<id>|<gene>), so best_paralog names the owning locus.
    reciprocal_best_paralog_check.py \\
        --home_query ${panel} \\
        --target_faa ${target_faa} \\
        --target_gff ${target_gff} \\
        --locus_id "${locus_id}" \\
        --genome_name "${genome_name}" \\
        --output ${locus_id}__${genome_name}.locus_ownership.tsv
    """
}
