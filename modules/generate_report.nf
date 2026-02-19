process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    val region_files
    val hits_dirs
    val augmented_genes
    path qc_json
    
    output:
    path "synterra_report.json", emit: report
    
    script:
    def flattenInput = { obj ->
        if (obj == null) {
            return []
        }
        if (obj instanceof Collection) {
            def out = []
            obj.each { item ->
                out.addAll(flattenInput(item))
            }
            return out
        }
        return [obj]
    }
    def regionList = flattenInput(region_files).collect { it.toString() }.join('\n')
    def hitsList = flattenInput(hits_dirs).collect { it.toString() }.join('\n')
    def augmentedList = flattenInput(augmented_genes).collect { it.toString() }.join('\n')
    """
    mkdir -p results/regions results/hits results/augmented

    cat > .hits.list <<'EOF'
${hitsList}
EOF
    cat > .regions.list <<'EOF'
${regionList}
EOF
    cat > .augmented.list <<'EOF'
${augmentedList}
EOF
    
    # Check if inputs exist
    while IFS= read -r d; do
        [ -z "\$d" ] && continue
        if [ -d "\$d" ]; then
            cp -r "\$d"/* results/hits/ 2>/dev/null || true
        fi
    done < .hits.list
    
    while IFS= read -r f; do
       [ -z "\$f" ] && continue
       if [ -f "\$f" ]; then
           cp "\$f" results/augmented/ 2>/dev/null || true
       fi
    done < .augmented.list
    
    while IFS= read -r f; do
       [ -z "\$f" ] && continue
       if [ -f "\$f" ]; then
           cp "\$f" results/regions/ 2>/dev/null || true
       fi
    done < .regions.list
    
    ${projectDir}/bin/generate_report.py --results_dir results --qc_json "${qc_json}" --output synterra_report.json
    """
}
