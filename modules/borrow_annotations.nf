process BORROW_ANNOTATIONS {
    tag "borrow_annotations"
    label 'process_medium'

    input:
    path home_genome
    path predicted_gff   // the home gene predictions (PREPARE_HOME_PROTEOME, no home GFF)
    path genomes_dir     // Directory with target genomes + optional GFFs
    path goi_bed         // the primary home locus BED (SPLIT_LOCI, by index: lib/HomeLoci.groovy)
    path sorted_genomes  // PHYLO_SORT's order: the closest annotated target is borrowed from
    path goi_annotation  // ANNOTATE_GOI: genes over the query's home model stay out of the proteome
    val n_flanking
    val run_token        // Emits only when home GFF is missing

    output:
    // One merged home annotation (predictions + borrowed gene models, batch HG-5b) and the
    // home proteome built from it, in the folder layout of PREPARE_HOME_PROTEOME's.
    path "merged_home_annotations.gff", emit: gff
    path "borrowed_proteins.faa", emit: proteins
    path "home_proteome_db", emit: db

    script:
    """
    borrow_annotations.py \\
        --home_genome $home_genome \\
        --predicted_gff $predicted_gff \\
        --genomes_dir $genomes_dir \\
        --goi_bed $goi_bed \\
        --sorted_genomes $sorted_genomes \\
        --goi_annotation $goi_annotation \\
        --n_flanking $n_flanking \\
        --window ${params.pred_flank_window} \\
        --target_flank ${params.borrow_target_flank} \\
        --output_gff merged_home_annotations.gff \\
        --output_proteins borrowed_proteins.faa \\
        --output_db_dir home_proteome_db ${GeneticCodes.homeArg(params)} ${GeneticCodes.targetArgs(params)}
    if [ -s home_proteome_db/home_proteome.faa ]; then
        mmseqs createdb home_proteome_db/home_proteome.faa home_proteome_db/db
    fi
    """
}
