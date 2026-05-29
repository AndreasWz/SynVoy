// Prefix a per-locus region file with its home-locus id before it is collected
// for GENERATE_REPORT. Invoked once per file type via aliased imports in main.nf
// (STAGE_REGION_GFF / STAGE_REGION_FAA / STAGE_REGION_HOMOLOGY / STAGE_REGION_SCORES).
//
// Why: ITERATIVE_SEARCH / CLUSTER_REGIONS emit one region file per (locus, genome)
// but name them by genome only (e.g. GCA_036250285.1.fna.gff), so the same basename
// recurs across loci. GENERATE_REPORT's staging step flat-cps all collected files into
// one directory, where same-named files from different loci overwrite each other —
// the report then silently saw only one locus per genome (HIGH-confidence GOI under-
// counted; the cross-locus duplicates of docs/TODO.md §1h hidden before dedup; and the
// §1h follow-up `total_new_genes` + per-locus region-score counts under-reported).
//
// Tagging with "<locus_id>__" makes every staged file unique and lets generate_report.py
// strip the prefix via canonical_genome_id() and dedupe per genome. Cluster/plot
// consumers continue to read the unprefixed ITERATIVE_SEARCH/CLUSTER_REGIONS outputs
// directly and are unaffected.
process STAGE_REGION_FOR_REPORT {
    tag "${locus_id}"

    input:
    tuple val(locus_id), path(region_file)

    output:
    path "${locus_id}__${region_file.name}"

    script:
    """
    cp -L "${region_file}" "${locus_id}__${region_file.name}"
    """
}
