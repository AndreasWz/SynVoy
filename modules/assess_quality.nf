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
    find -L ${genomes_dir} -maxdepth 1 -type f \\( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \\) | while read f; do
        if [ "\$first" = "true" ]; then first=false; else echo "," >> genome_qc_summary.json; fi
        assess_genome_quality.py --genome "\$f" --output tmp.json
        cat tmp.json >> genome_qc_summary.json
    done
    echo "]" >> genome_qc_summary.json
    """
}
