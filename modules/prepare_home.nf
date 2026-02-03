process PREPARE_HOME_PROTEOME {
    tag "home_prot"
    
    input:
    path home_genome
    path home_gff
    
    output:
    path "home_proteome.faa", emit: faa
    path "home_proteome_db", emit: db
    path "home_predicted.gff", emit: gff, optional: true
    
    script:
    """
    if [ "$home_gff" != "NO_GFF" ]; then
        echo "Extracting proteins from Home GFF..."
        gff_to_faa.py --gff $home_gff --genome $home_genome --output home_proteome.faa
        # No predicted GFF needed - use provided one
        touch home_predicted.gff.skip
    else
        echo "No Home GFF provided. Predicting proteins with Prodigal..."
        
        # Split genome into individual sequences to avoid Prodigal 32MB limit
        mkdir -p split_seqs
        awk '/^>/{if(fn) close(fn); fn="split_seqs/seq_"++c".fna"} {print > fn}' $home_genome
        
        # Run Prodigal on each sequence, collect both proteins and GFF
        > home_proteome.faa
        > home_predicted.gff
        for seq in split_seqs/*.fna; do
            prodigal -i "\$seq" -a tmp_proteins.faa -f gff -o tmp.gff -p meta -q 2>/dev/null || true
            if [ -f tmp_proteins.faa ]; then
                cat tmp_proteins.faa >> home_proteome.faa
                rm tmp_proteins.faa
            fi
            if [ -f tmp.gff ]; then
                # Skip header lines for all but first file
                grep -v "^#" tmp.gff >> home_predicted.gff 2>/dev/null || true
                rm tmp.gff
            fi
        done
        rm -rf split_seqs
        
        # Add GFF header
        if [ -s home_predicted.gff ]; then
            mv home_predicted.gff tmp_body.gff
            echo "##gff-version 3" > home_predicted.gff
            cat tmp_body.gff >> home_predicted.gff
            rm tmp_body.gff
        fi
    fi
    
    # Create MMseqs DB for fast searching
    mkdir home_proteome_db
    mmseqs createdb home_proteome.faa home_proteome_db/db
    """
}
