process PREPARE_HOME_PROTEOME {
    tag "home_prot"
    
    input:
    path home_genome
    path home_gff
    path goi_bed
    val settings    // augustus_species of a home without a GFF (lib/HomeGenes.groovy, stage 5c)
    
    output:
    path "home_proteome.faa", emit: faa
    path "home_proteome_db", emit: db
    path "home_predicted.gff", emit: gff, optional: true
    
    script:
    def gc_arg = GeneticCodes.homeArg(params)   // the home genome's code; '' for table 1
    def augustus = HomeGenes.augustusSpecies(params, settings)
    """
    if [ "$home_gff" != "NO_GFF" ]; then
        echo "Extracting proteins from Home GFF..."
        gff_to_faa.py --gff $home_gff --genome $home_genome --output home_proteome.faa \\
            --coords_tsv home_genes.tsv ${gc_arg}
        # No predicted GFF needed - use provided one
        touch home_predicted.gff.skip
    else
        echo "No Home GFF provided. Predicting proteins (GOI regions only)..."
        prodigal_on_regions.py \\
            --genome $home_genome \\
            --goi_bed $goi_bed \\
            --window ${params.pred_flank_window} \\
            --output_faa home_proteome.faa \\
            --output_gff home_predicted.gff \\
            --fallback_full_genome ${params.prodigal_full_genome_fallback} \\
            --gene_predictor ${params.gene_predictor} \\
            --augustus_species ${augustus} ${gc_arg}
    fi

    # Guard: if GFF extraction produced an empty file, fall back to gene prediction
    if [ ! -s home_proteome.faa ]; then
        echo "WARNING: home_proteome.faa is empty after GFF extraction; falling back to gene prediction..."
        prodigal_on_regions.py \\
            --genome $home_genome \\
            --goi_bed $goi_bed \\
            --window ${params.pred_flank_window} \\
            --output_faa home_proteome.faa \\
            --output_gff home_predicted.gff \\
            --fallback_full_genome true \\
            --gene_predictor ${params.gene_predictor} \\
            --augustus_species ${augustus} ${gc_arg}
    fi

    # Create MMseqs DB for fast searching
    mkdir home_proteome_db
    mmseqs createdb home_proteome.faa home_proteome_db/db
    # The search's specificity check (bin/synvoy_specificity.py) aligns each GOI model
    # to every home protein and needs gene positions: both travel in the folder the
    # search already receives (--home_db_dir), so ITERATIVE_SEARCH needs no new input.
    cp home_proteome.faa home_proteome_db/home_proteome.faa
    if [ -s home_genes.tsv ]; then cp home_genes.tsv home_proteome_db/home_genes.tsv; fi
    """
}
