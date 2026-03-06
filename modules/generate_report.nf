process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path region_files  // collected region .faa files (may be sentinel placeholder)
    path hits_dirs     // collected hits directories (may be sentinel placeholder)
    path augmented     // unused, kept for signature compatibility
    path qc_json
    
    output:
    path "synterra_report.json", emit: report
    
    script:
    """
    mkdir -p staged_results/regions staged_results/hits

    # Stage region files (individual .faa files), skipping sentinel placeholders
    for f in ${region_files}; do
        case "\$f" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_SPECIES_MAP) continue ;;
        esac
        if [ -f "\$f" ]; then
            cp "\$f" staged_results/regions/ 2>/dev/null || true
        fi
    done

    # Stage hits (directories containing .m8 files), skipping sentinel placeholders
    for d in ${hits_dirs}; do
        case "\$d" in
            NO_REGIONS|NO_HITS|NO_AUGMENTED|NO_GFF|NO_SPECIES_MAP) continue ;;
        esac
        if [ -d "\$d" ]; then
            cp -r "\$d"/* staged_results/hits/ 2>/dev/null || true
        elif [ -f "\$d" ] && echo "\$d" | grep -q '\\.m8\$'; then
            cp "\$d" staged_results/hits/ 2>/dev/null || true
        fi
    done

    ${projectDir}/bin/generate_report.py \\
        --results_dir staged_results \\
        --qc_json "${qc_json}" \\
        --output synterra_report.json
    """
}
