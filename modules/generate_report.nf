process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path region_files
    path hits_dirs
    path augmented_genes
    path qc_json
    
    output:
    path "synterra_report.json", emit: report
    
    script:
    """
    mkdir -p results/regions results/hits results/augmented
    
    # Check if inputs exist
    for d in ${hits_dirs}; do
        if [ -d "\$d" ]; then
            cp -r \$d/* results/hits/ 2>/dev/null || true
        fi
    done
    
    for f in ${augmented_genes}; do
       if [ -f "\$f" ]; then
           cp \$f results/augmented/ 2>/dev/null || true
       fi
    done
    
    for f in ${region_files}; do
       if [ -f "\$f" ]; then
           cp \$f results/regions/ 2>/dev/null || true
       fi
    done
    
    generate_report.py --results_dir results --qc_json ${qc_json} --output synterra_report.json
    """
}
