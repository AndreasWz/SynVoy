// Prefix a per-locus region file with its home-locus id before it is collected
// for GENERATE_REPORT.
//
// Why: ITERATIVE_SEARCH emits one region GFF per (locus, genome) but names them by
// genome only (e.g. GCA_036250285.1.fna.gff), so the same basename recurs across loci.
// The GENERATE_REPORT staging step flat-copies all collected files into one directory,
// where same-named files from different loci overwrite each other — the report then
// silently saw only one locus per genome (HIGH-confidence GOI undercounted, and the
// cross-locus duplicates of docs/TODO.md §1h hidden before dedup could run).
//
// Tagging with "<locus_id>__" makes every staged file unique and lets
// generate_report.py attribute each GOI hit to its source home locus for dedup.
// Cluster/plot consumers read ITERATIVE_SEARCH.out.gff directly and are unaffected.
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
