process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    cache 'deep'
    
    input:
    // stageAs '?/*' places each collected file under a unique numbered subdir
    // (e.g. 1/file.faa, 2/file.faa ...) to avoid filename collision across loci.
    path region_files,   stageAs: 'region_genes/?/*'
    path region_gffs,    stageAs: 'region_gffs/?/*'
    path homology_files, stageAs: 'homology/?/*'
    path hits_dirs,      stageAs: 'hits/?/*'
    path augmented       // unused sentinel – single file, no collision
    path qc_json
    path score_files,    stageAs: 'scores/?/*'
    path paralog_check_files, stageAs: 'paralog_check/?/*'  // §1j Phase B
    val paralog_confusion_min_gap
    val qc_policy
    path locus_ownership_files, stageAs: 'locus_ownership/?/*'  // §1m
    path panel_meta            // §1m panel sidecar — single file or sentinel
    path phylo_placement_files, stageAs: 'phylo_placement/?/*'  // F9
    path effective_params      // RESOLVE_EFFECTIVE_PARAMS output, or a sentinel -> run_settings

    output:
    path "synvoy_report.json", emit: report

    script:
    """
    mkdir -p staged_results/regions staged_results/hits staged_results/scores staged_results/paralog_check staged_results/locus_ownership staged_results/phylo_placement

    # Stage F9 phylo-placement verdicts from numbered subdirs.
    for f in \$(find -L phylo_placement -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_PHYLO_PLACEMENT|NO_PARALOG_CHECK|NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_GFFS|NO_HOMOLOGY|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/phylo_placement/ 2>/dev/null || true
    done

    # Stage paralog-check TSVs from numbered subdirs (§1j Phase B).
    for f in \$(find -L paralog_check -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_PARALOG_CHECK|NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_GFFS|NO_HOMOLOGY|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/paralog_check/ 2>/dev/null || true
    done

    # Stage region outputs (.faa / .gff / .homology.tsv) from numbered subdirs
    for f in \$(find -L region_genes region_gffs homology -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_GFFS|NO_HOMOLOGY|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/regions/ 2>/dev/null || true
    done

    # Stage the per-genome tables the search writes into its hits/ folder: the §1x
    # non-syntenic candidates and the genomes it could not search. Each hits/ folder
    # arrives as hits/<n>/hits/<file>, one level deeper than the old loops looked, so none
    # of them ever reached the report (bug hunt 2026-10-13). Their names carry the locus,
    # so one flat folder is safe. The raw .m8 hits are not staged: they are large and the
    # old copies of them were never read.
    find -L hits -type f \\( -name '*.nonsyntenic.tsv' -o -name '*.not_searched.tsv' \\) \\
        -exec cp {} staged_results/hits/ \\; 2>/dev/null || true

    # Stage structured region scores
    for f in \$(find -L scores -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/scores/ 2>/dev/null || true
    done

    # Stage §1m locus-ownership TSVs from numbered subdirs.
    for f in \$(find -L locus_ownership -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_LOCUS_OWNERSHIP|NO_PARALOG_CHECK|NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_GFFS|NO_HOMOLOGY|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/locus_ownership/ 2>/dev/null || true
    done

    # §1m panel meta sidecar (single file; sentinel when ownership disabled/empty).
    # generate_report.py expects it INSIDE the ownership dir under its build-time name:
    #   <results_dir>/locus_ownership/home_paralog_panel.faa.meta.tsv
    case "${panel_meta}" in
        NO_LOCUS_OWNERSHIP) : ;;
        *) if [ -s "${panel_meta}" ]; then cp "${panel_meta}" staged_results/locus_ownership/home_paralog_panel.faa.meta.tsv; fi ;;
    esac

    ${projectDir}/bin/generate_report.py \\
        --results_dir staged_results \\
        --qc_json "${qc_json}" \\
        --qc_policy "${qc_policy}" \\
        --paralog_confusion_min_gap ${paralog_confusion_min_gap} \\
        --locus_ownership_tiebreak_gap ${params.locus_ownership_tiebreak_gap} \\
        --identity_decoupled_min_identity ${params.identity_decoupled_min_identity} \\
        --identity_decoupled_max_qcov ${params.identity_decoupled_max_qcov} \\
        ${params.disable_coverage_demotion.toString().toBoolean() ? '--disable_coverage_demotion' : ''} \\
        --phylo_placement_dir staged_results/phylo_placement \\
        ${params.phylo_placement_promote.toString().toBoolean() ? '--phylo_placement_promote' : ''} \\
        --effective_params ${effective_params} \\
        --output synvoy_report.json
    """
}
