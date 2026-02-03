process ITERATIVE_SEARCH {
    tag "iterative"
    label 'process_long'

    input:
    tuple val(locus_id), path(initial_db)
    path sorted_genomes_file
    path genomes // List or Dir
    path home_db // MMseqs DB directory
    val flanking_genes
    val min_score
    val sensitivity

    output:
    tuple val(locus_id), path("iterative_results/expanded_db.faa"), emit: expanded_db
    tuple val(locus_id), path("iterative_results/hits"), emit: hits, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.faa"), emit: region_genes, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.gff"), emit: gff, optional: true
    tuple val(locus_id), path("iterative_results/regions/*.homology.tsv"), emit: homology, optional: true


    script:
    """
    iterative_search_runner.py \\
        --initial_db $initial_db \\
        --sorted_genomes $sorted_genomes_file \\
        --genomes_dir $genomes \\
        --home_db_dir $home_db \\
        --output_dir iterative_results \\
        --threads ${task.cpus} \\
        --prefix ${locus_id}
    """
}
