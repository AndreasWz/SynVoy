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

    output:
    path "synvoy_report.json", emit: report

    script:
    """
    mkdir -p staged_results/regions staged_results/hits staged_results/scores staged_results/paralog_check staged_results/locus_ownership

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

    # Stage hits directories (containing .m8 files)
    for d in \$(find -L hits -mindepth 1 -maxdepth 1 -type d 2>/dev/null); do
        dname=\$(basename "\$d")
        case "\$dname" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        if [ -d "\$d" ]; then
            cp -r "\$d"/* staged_results/hits/ 2>/dev/null || true
        fi
    done
    # Also pick up any plain .m8 files staged directly under hits/
    for f in \$(find -L hits -maxdepth 2 -name '*.m8' -type f 2>/dev/null); do
        cp "\$f" staged_results/hits/ 2>/dev/null || true
    done

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
        --locus_ownership_tiebreak_gap ${params.locus_ownership_tiebreak_gap} \
        --identity_decoupled_min_identity ${params.identity_decoupled_min_identity} \
        --identity_decoupled_max_qcov ${params.identity_decoupled_max_qcov} \
        ${params.disable_coverage_demotion ? '--disable_coverage_demotion' : ''} \\
        --output synvoy_report.json
    """
}
