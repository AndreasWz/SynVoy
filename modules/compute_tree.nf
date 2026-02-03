process COMPUTE_TREE {
    tag "tree"
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(locus_id), path(fasta_files)

    output:
    tuple val(locus_id), path("*.nwk"), emit: tree

    script:
    """
    # Concatenate all fasta files (flanking + new hits)
    cat ${fasta_files} > all_sequences.faa
    
    compute_tree.py \\
        --input all_sequences.faa \\
        --output ${locus_id}_tree.nwk \\
        --threads ${task.cpus}
    """
}
