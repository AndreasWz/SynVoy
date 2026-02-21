process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path region_files  // collected region .faa files (may be empty placeholder)
    path hits_dirs     // collected hits directories (may be empty placeholder)
    path augmented     // unused, kept for signature compatibility
    path qc_json
    
    output:
    path "synterra_report.json", emit: report
    
    script:
    """
    mkdir -p staged_results/regions staged_results/hits

    # Stage region files (individual .faa files)
    for f in ${region_files}; do
        if [ -f "\$f" ]; then
            cp "\$f" staged_results/regions/ 2>/dev/null || true
        fi
    done

    # Stage hits (directories containing .m8 files)
    for d in ${hits_dirs}; do
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
