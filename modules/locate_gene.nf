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

    script:
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
    # Removed --mask-mode 0 (invalid)
    mmseqs easy-search $gene $home_genome hits_mmseqs.m8 tmp \\
        --search-type \$SEARCH_TYPE \\
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
            print \$2, s-1, e, "mmseqs", \$11, str 
        }' hits_mmseqs.m8 > mmseqs.bed
    fi

    # Process BLAST results
    touch blast.bed
    if [ -s hits_blast.txt ]; then
        awk -v OFS="\\t" '{ 
            s=\$9; e=\$10; 
            if(s>e){t=s; s=e; e=t; str="-"} else {str="+"} 
            print \$2, s-1, e, "blast", \$11, str 
        }' hits_blast.txt > blast.bed
    fi

    # Run Python Merge Script
    merge_hits.py --mmseqs mmseqs.bed --blast blast.bed --output home_gene_location.bed
    """
}
