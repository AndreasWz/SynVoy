process LOCATE_GENE {
    tag "$home_genome"
    label 'process_medium'
    publishDir "${params.outdir}/intermediate/locate_gene", mode: 'copy'

    input:
    path gene
    path home_genome

    output:
    path "home_gene_location.bed", emit: bed
    path "hits_blast.txt", emit: blast_hits
    path "hits_mmseqs.m8", emit: mmseqs_hits
    path "hit_profile.json", emit: hit_profile, optional: true
    path "auto_preset.json", emit: auto_preset, optional: true

    script:
    def preset_override_arg = params.preset_override ? "--preset_override ${params.preset_override}" : ''
    """
    # 1. MMseqs2 Search
    # Check query type (Protein vs DNA) - look for amino acid-specific chars
    # Exclude E (valid IUPAC nucleotide) to avoid false positives
    is_prot=\$(grep -v "^>" $gene | head -c 200 | grep -q "[DFHIKLMPQRSVWY]" && echo "true" || echo "false")
    
    if [ "\$is_prot" = "true" ]; then
        echo "Detected Protein Query"
        SEARCH_TYPE=2 # translated
        BLAST_CMD="tblastn"
    else
        echo "Detected Nucleotide Query"
        SEARCH_TYPE=3 # nucleotide
        BLAST_CMD="blastn"
    fi

    # MMSEQS
    # Keep an empty placeholder so downstream steps can proceed on BLAST-only
    # fallback if MMseqs is killed on very large/fragmented genomes.
    touch hits_mmseqs.m8
    mmseqs easy-search $gene $home_genome hits_mmseqs.m8 tmp \\
        --search-type \$SEARCH_TYPE \\
        --threads ${task.cpus} \\
        --split-memory-limit ${params.mmseqs_split_memory_limit} \\
        -v ${params.mmseqs_verbosity} \\
        -e ${params.search_evalue} \\
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits" \\
        -s ${params.mmseqs_sensitivity} || echo "MMSeqs search failed or found no hits"

    # 2. BLAST Search
    makeblastdb -in $home_genome -dbtype nucl
    
    if [ "\$BLAST_CMD" = "tblastn" ]; then
        tblastn -query $gene -db $home_genome -out hits_blast.txt -outfmt 6
    else
        blastn -query $gene -db $home_genome -out hits_blast.txt -outfmt 6
    fi

    # 3. Combine and convert to BED for merging
    
    # Process MMseqs results
    touch mmseqs.bed
    if [ -s hits_mmseqs.m8 ]; then
        awk -v OFS="\\t" '{
            s=\$9; e=\$10;
            if(s>e){t=s; s=e; e=t; str="-"} else {str="+"}
            print \$2, s-1, e, "mmseqs", \$11, str, \$12
        }' hits_mmseqs.m8 > mmseqs.bed
    fi

    # Process BLAST results
    touch blast.bed
    if [ -s hits_blast.txt ]; then
        awk -v OFS="\\t" '{
            s=\$9; e=\$10;
            if(s>e){t=s; s=e; e=t; str="-"} else {str="+"}
            print \$2, s-1, e, "blast", \$11, str, \$12
        }' hits_blast.txt > blast.bed
    fi

    # Run Python Merge Script
    merge_hits.py \\
        --mmseqs mmseqs.bed \\
        --blast blast.bed \\
        --max_evalue ${params.search_evalue} \\
        --output home_gene_location.bed

    # §1f: profile the hit distribution and recommend a query-type preset (advisory).
    # The recommendation is printed to this process log + saved to intermediate/locate_gene/
    # so a user who didn't pick a preset can re-run with the suggested -profile preset_*.
    # Best-effort: never fail LOCATE_GENE if profiling hiccups.
    profile_hits.py \\
        --blast hits_blast.txt \\
        --mmseqs hits_mmseqs.m8 \\
        --query $gene \\
        --output hit_profile.json || echo "hit profiling skipped"
    if [ -s hit_profile.json ]; then
        auto_select_preset.py \\
            --hit_profile hit_profile.json \\
            --output auto_preset.json \\
            ${preset_override_arg} \\
            || echo "preset advisor skipped"
    fi
    """
}
