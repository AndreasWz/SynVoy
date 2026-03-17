process ASSESS_GENOME_QUALITY {
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    path genomes_dir
    
    output:
    path "genome_qc_summary.json", emit: json
    
    script:
    """
    echo "[" > genome_qc_summary.json
    first=true
    # Find files safely - follow symlinks
    find -L ${genomes_dir} -maxdepth 1 -type f \\( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \\) > genome_files.txt
    while read f; do
        if [ "\$first" = "true" ]; then first=false; else echo "," >> genome_qc_summary.json; fi
        assess_genome_quality.py \\
            --genome "\$f" \\
            --output tmp.json \\
            --min_n50 ${params.bad_min_n50} \\
            --max_contigs ${params.bad_max_contigs}
        cat tmp.json >> genome_qc_summary.json
    done < genome_files.txt
    echo "]" >> genome_qc_summary.json
    """
}
