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
    # Concatenate all fasta files
    cat ${fasta_files} > all_sequences.faa
    
    # CRITICAL: Filter to only GOI sequences for phylogenetic tree
    # This excludes flanking genes - tree should only show GOI homologs
    # Match headers containing "GOI_" or the query ID pattern
    python3 -c "
import sys
keep = False
seen_names = {}
for line in open('all_sequences.faa'):
    if line.startswith('>'):
        name = line.strip()[1:]
        is_goi = 'GOI_' in name or 'GOI|' in name
        if is_goi:
            if name in seen_names:
                seen_names[name] += 1
                unique = name + '_dup' + str(seen_names[name])
            else:
                seen_names[name] = 1
                unique = name
            sys.stdout.write('>' + unique + chr(10))
            keep = True
        else:
            keep = False
    elif keep:
        sys.stdout.write(line)
" > goi_only.faa
    
    # Check if we have sequences
    count=\$(grep -c '^>' goi_only.faa 2>/dev/null || echo 0)
    echo "Filtered to \$count GOI sequences for tree"
    
    if [ "\$count" -lt 3 ]; then
        echo "Not enough GOI sequences (<3) for tree, creating placeholder"
        echo "(GOI_placeholder:0.0);" > ${locus_id}_tree.nwk
    else
        compute_tree.py \\
            --input goi_only.faa \\
            --output ${locus_id}_tree.nwk \\
            --threads ${task.cpus}
    fi
    """
}
