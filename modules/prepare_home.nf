process PREPARE_HOME_PROTEOME {
    tag "home_prot"
    
    input:
    path home_genome
    path home_gff
    path goi_bed
    
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
        echo "No Home GFF provided. Predicting proteins with Prodigal (GOI regions only)..."
        prodigal_on_regions.py \\
            --genome $home_genome \\
            --goi_bed $goi_bed \\
            --window ${params.pred_flank_window} \\
            --output_faa home_proteome.faa \\
            --output_gff home_predicted.gff \\
            --fallback_full_genome ${params.prodigal_full_genome_fallback}
    fi
    
    # Create MMseqs DB for fast searching
    mkdir home_proteome_db
    mmseqs createdb home_proteome.faa home_proteome_db/db
    """
}
