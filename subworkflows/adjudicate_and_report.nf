// Adjudication + staging + report — extracted from main.nf's `workflow {}` body.
//
// WHY THIS FILE EXISTS: the main workflow body sits near the JVM 64 kB method-size
// limit (a 4-line comment has previously re-triggered `UTF8 string too large`), which
// blocked every new process call — including routing rescue models through the locus
// ownership check. A named sub-workflow gets its own method budget, so this is where
// new adjudication/staging steps belong from now on. Keep main.nf's body lean.
//
// Everything here runs AFTER the per-locus search: the two rescue passes, the
// reciprocal-best paralog check, locus ownership, the per-locus staging that keeps
// same-named files from different loci from overwriting each other (§1h), and the
// final report.

include { GENERATE_REPORT } from '../modules/generate_report.nf'
include { STAGE_REGION_FOR_REPORT as STAGE_REGION_GFF       } from '../modules/stage_for_report.nf'
include { STAGE_REGION_FOR_REPORT as STAGE_REGION_FAA       } from '../modules/stage_for_report.nf'
include { STAGE_REGION_FOR_REPORT as STAGE_REGION_HOMOLOGY  } from '../modules/stage_for_report.nf'
include { STAGE_REGION_FOR_REPORT as STAGE_REGION_SCORES    } from '../modules/stage_for_report.nf'
include { RESCUE_STRONG_SYNTENY } from '../modules/rescue_strong_synteny.nf'
include { RESCUE_GOI_HULL } from '../modules/rescue_goi_hull.nf'
include { RECIPROCAL_BEST_PARALOG } from '../modules/reciprocal_best_paralog.nf'
include { BUILD_HOME_PARALOG_PANEL; ASSIGN_LOCUS_OWNERSHIP } from '../modules/assign_locus_ownership.nf'

workflow ADJUDICATE_AND_REPORT {

    take:
    iter_gff            // ITERATIVE_SEARCH.out.gff          (locus_id, [gff...])
    iter_region_genes   // ITERATIVE_SEARCH.out.region_genes (locus_id, [faa...])
    iter_homology       // ITERATIVE_SEARCH.out.homology     (locus_id, [tsv...])
    iter_hits           // ITERATIVE_SEARCH.out.hits         (locus_id, m8)
    cluster_scores      // CLUSTER_REGIONS.out.scores        (locus_id, genome, scores.tsv)
    genomes_dir         // value channel: dir holding the target FASTAs
    query_faa           // value channel: the NORMALIZED GOI protein (see F8)
    home_gff            // value channel: home GFF
    home_proteome       // value channel: home proteome FASTA
    locus_beds          // collected per-locus BEDs
    qc_summary          // genome QC summary JSON

    main:

    // Use sentinel files so that collect() always yields a valid path list
    // that Nextflow can stage into the process work directory.
    def no_regions_sentinel = file("${projectDir}/assets/sentinels/NO_REGIONS")
    def no_gffs_sentinel    = file("${projectDir}/assets/sentinels/NO_GFFS")
    def no_homology_sentinel = file("${projectDir}/assets/sentinels/NO_HOMOLOGY")
    def no_hits_sentinel    = file("${projectDir}/assets/sentinels/NO_HITS")
    def no_augmented_sentinel = file("${projectDir}/assets/sentinels/NO_AUGMENTED")
    def no_scores_sentinel = file("${projectDir}/assets/sentinels/NO_SCORES")
    
    // Tag every per-locus region file (.gff / .faa / .homology.tsv / .scores.tsv)
    // with its home-locus id before GENERATE_REPORT's flat-cp staging, so same-named
    // files from different loci no longer overwrite each other. The GFF case landed
    // first (docs/TODO.md §1h); this also covers .faa / .homology.tsv / .scores.tsv
    // (the §1h follow-up that previously under-reported `total_new_genes` and the
    // per-locus region score counts in multi-locus runs).
    // Python side: generate_report.py's canonical_genome_id() strips the locus prefix
    // before grouping by genome, so no parser change is needed.

    // §1e follow-up: rescue pass against blocks cluster_grs classified as
    // goi_missing_but_strong_synteny. Inputs: scores.tsv per (locus, genome),
    // the genomes dir to resolve the target FASTA, and the GOI query FASTA.
    // Output GFFs are mixed with iterative_search's per-locus GFFs into the
    // staging channel so the report parser sees them via the same path.
    rescue_query_ch = query_faa

    RESCUE_STRONG_SYNTENY(
        cluster_scores,
        genomes_dir,
        rescue_query_ch
    )

    // §1j Phase B: reciprocal-best paralog check. Per (locus, genome), align
    // each recovered GOI target protein against every paralog in the home
    // query FASTA. generate_report.py uses the result to flag MEDIUM-confidence
    // calls whose best home paralog differs from the modal best at the same
    // locus — the TP63 ↔ TP73 bleed observed in docs/TP53_PARALOG_ASSESSMENT.md.
    // No-op when params.query has a single sequence (nothing reciprocal).
    paralog_faa_ch = iter_region_genes.transpose()
        .map { locus_id, faa ->
            def stem = faa.name.replaceFirst(/\.faa$/, '')
            tuple("${locus_id}::${stem}", locus_id, stem, faa)
        }
    paralog_gff_ch = iter_gff.transpose()
        .map { locus_id, gff ->
            def stem = gff.name.replaceFirst(/\.gff3?$/, '')
            tuple("${locus_id}::${stem}", gff)
        }
    paralog_inputs_ch = paralog_faa_ch
        .map { key, locus_id, stem, faa -> tuple(key, locus_id, stem, faa) }
        .join(paralog_gff_ch)
        .map { key, locus_id, stem, faa, gff ->
            def genome_name = stem.replaceFirst(/\.(fna|fa|fasta)$/, '')
            tuple(locus_id, genome_name, faa, gff)
        }

    RECIPROCAL_BEST_PARALOG(
        paralog_inputs_ch,
        rescue_query_ch   // same home query as the rescue pass — multi-FASTA = paralog panel
    )

    // §1m fumble fix: GOI synteny-hull rescue (reuses paralog_inputs_ch's region GFF;
    // models a GOI sitting in a flanking GAP, e.g. cow decorin chr5:21 Mb). Mixes into
    // STAGE_REGION_GFF like the strong-synteny rescue.
    RESCUE_GOI_HULL(paralog_inputs_ch, genomes_dir, rescue_query_ch)

    // §1m: locus ownership. Build a panel of the home GOI-paralog at each home
    // locus (>locus_<id>|<gene>), then SW-align every recovered GOI target
    // against it so the report can re-attribute orthologs filed under the wrong
    // home locus (decorin recovered by the chr9 OMD locus belongs to the chr12
    // DCN locus) and relabel paralogs carried under the GOI name (asporin at
    // 49.5% labelled GOI_DCN). No-op when --disable_locus_ownership.
    BUILD_HOME_PARALOG_PANEL(
        home_gff,
        home_proteome,
        locus_beds.flatten().collect(),
        query_faa
    )

    // Rescue models are mixed into the region GFFs at STAGE_REGION_GFF, which is
    // DOWNSTREAM of ownership — so until now they were never RBH-checked against the
    // home panel, i.e. the guard against a mislabelled GOI never saw the calls most
    // likely to be one. Feed them in as their own (locus, genome) ownership task,
    // tagged `.hull_rescue` so the output filename cannot collide with the main task.
    // canonical_genome_id() already strips that suffix, so the report folds the rows
    // back onto the real genome (the §18 phantom-genome fix).
    hull_rescue_ownership_ch = RESCUE_GOI_HULL.out.faa
        .map { locus_id, genome_name, faa -> tuple("${locus_id}::${genome_name}", faa) }
        .join(
            RESCUE_GOI_HULL.out.gff.map { locus_id, gff ->
                tuple("${locus_id}::${gff.name.replaceFirst(/\.hull_rescue\.gff$/, '')}", gff)
            }
        )
        .map { key, faa, gff ->
            def parts = key.split("::")
            tuple(parts[0], "${parts[1]}.hull_rescue", faa, gff)
        }
        .filter { _l, _g, faa, _gff -> faa.size() > 0 }

    ASSIGN_LOCUS_OWNERSHIP(
        paralog_inputs_ch.mix(hull_rescue_ownership_ch),
        BUILD_HOME_PARALOG_PANEL.out.faa.first()  // value channel: broadcast panel to every genome (see module)
    )

    STAGE_REGION_GFF(
        iter_gff.transpose()
            .mix(RESCUE_STRONG_SYNTENY.out.gff)
            .mix(RESCUE_GOI_HULL.out.gff)
    )
    STAGE_REGION_GFF.out
        .ifEmpty(no_gffs_sentinel)
        .collect()
        .set { collected_region_gffs }

    STAGE_REGION_FAA(
        iter_region_genes.transpose()
    )
    STAGE_REGION_FAA.out
        .ifEmpty(no_regions_sentinel)
        .collect()
        .set { collected_regions }

    STAGE_REGION_HOMOLOGY(
        iter_homology.transpose()
    )
    STAGE_REGION_HOMOLOGY.out
        .ifEmpty(no_homology_sentinel)
        .collect()
        .set { collected_homology }

    iter_hits
        .map { rec -> rec[1] }
        .ifEmpty(no_hits_sentinel)
        .collect()
        .set { collected_hits }

    // cluster_scores tuple is (locus_id, genome_name, score_file).
    // Drop genome_name and feed (locus_id, score_file) into STAGE_REGION_SCORES.
    STAGE_REGION_SCORES(
        cluster_scores.map { rec -> tuple(rec[0], rec[2]) }
    )
    STAGE_REGION_SCORES.out
        .ifEmpty(no_scores_sentinel)
        .collect()
        .set { collected_scores }

    // §1j Phase B paralog-check TSVs (one per (locus, genome) when params.query
    // has >=2 paralogs; the script no-ops to a header-only TSV when it has 1).
    def no_paralog_sentinel = file("${projectDir}/assets/sentinels/NO_PARALOG_CHECK")
    RECIPROCAL_BEST_PARALOG.out.tsv
        .map { rec -> rec[1] }
        .ifEmpty(no_paralog_sentinel)
        .collect()
        .set { collected_paralog_check }

    // §1m locus-ownership TSVs (one per (locus, genome)) + the panel meta sidecar.
    def no_ownership_sentinel = file("${projectDir}/assets/sentinels/NO_LOCUS_OWNERSHIP")
    ASSIGN_LOCUS_OWNERSHIP.out.tsv
        .map { rec -> rec[1] }
        .ifEmpty(no_ownership_sentinel)
        .collect()
        .set { collected_locus_ownership }
    collected_panel_meta = BUILD_HOME_PARALOG_PANEL.out.meta.ifEmpty(no_ownership_sentinel)

    // No standalone augmented proteins - pass sentinel file
    collected_augmented = channel.value(no_augmented_sentinel)

    GENERATE_REPORT(
        collected_regions,
        collected_region_gffs,
        collected_homology,
        collected_hits,
        collected_augmented,
        qc_summary,
        collected_scores,
        collected_paralog_check,
        params.paralog_confusion_min_gap,
        params.qc_fail_policy,
        collected_locus_ownership,
        collected_panel_meta
    )

    emit:
    report = GENERATE_REPORT.out.report
}
