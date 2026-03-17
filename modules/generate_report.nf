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
    val qc_policy
    
    output:
    path "synterra_report.json", emit: report
    
    script:
    """
    mkdir -p staged_results/regions staged_results/hits staged_results/scores

    # Stage region outputs (.faa / .gff / .homology.tsv) from numbered subdirs
    for f in \$(find region_genes region_gffs homology -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_GFFS|NO_HOMOLOGY|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/regions/ 2>/dev/null || true
    done

    # Stage hits directories (containing .m8 files)
    for d in \$(find hits -mindepth 1 -maxdepth 1 -type d 2>/dev/null); do
        dname=\$(basename "\$d")
        case "\$dname" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        if [ -d "\$d" ]; then
            cp -r "\$d"/* staged_results/hits/ 2>/dev/null || true
        fi
    done
    # Also pick up any plain .m8 files staged directly under hits/
    for f in \$(find hits -maxdepth 2 -name '*.m8' -type f 2>/dev/null); do
        cp "\$f" staged_results/hits/ 2>/dev/null || true
    done

    # Stage structured region scores
    for f in \$(find scores -type f 2>/dev/null); do
        fname=\$(basename "\$f")
        case "\$fname" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_SPECIES_MAP|NO_SCORES) continue ;;
        esac
        cp "\$f" staged_results/scores/ 2>/dev/null || true
    done

    ${projectDir}/bin/generate_report.py \\
        --results_dir staged_results \\
        --qc_json "${qc_json}" \\
        --qc_policy "${qc_policy}" \\
        --output synterra_report.json
    """
}
