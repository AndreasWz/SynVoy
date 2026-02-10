process PLOT_SYNTENY {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path home_bed
    path query_bed
    path home_gff
    path target_gffs
    val target_names
    path candidate_beds
    path homology_tsvs
    path tree
    path species_map

    output:
    path "*_synteny_plot.html", emit: plot
    path "*_tree.html", emit: tree, optional: true
    path "plot_inputs_*", emit: inputs, optional: true

    script:
    // Handle empty collections gracefully
    def gffs_str = target_gffs ? target_gffs.join(' ') : ''
    def names_str = target_names ? target_names.join(' ') : ''
    def cands_str = candidate_beds ? candidate_beds.join(' ') : ''
    def homo_str = homology_tsvs ? homology_tsvs.join(' ') : ''
    
    // Only include arguments if they have values
    def gffs_arg = gffs_str ? "--target_gffs ${gffs_str}" : ""
    def names_arg = names_str ? "--target_names ${names_str}" : ""
    def cands_arg = cands_str ? "--candidate_beds ${cands_str}" : ""
    def homo_arg = homo_str ? "--homology_tsvs ${homo_str}" : ""
    def species_arg = species_map.name != 'NO_SPECIES_MAP' ? "--species_map ${species_map}" : ""
    
    """
    inputs_dir="plot_inputs_${home_bed.baseName}"
    mkdir -p "\$inputs_dir"
    cp $home_bed "\$inputs_dir/" || true
    cp $query_bed "\$inputs_dir/" || true
    if [ "$home_gff" != "NO_GFF" ]; then
        cp $home_gff "\$inputs_dir/" || true
    fi
    if [ -n "${gffs_str}" ]; then
        cp ${gffs_str} "\$inputs_dir/" || true
    fi
    if [ -n "${cands_str}" ]; then
        cp ${cands_str} "\$inputs_dir/" || true
    fi
    if [ -n "${homo_str}" ]; then
        cp ${homo_str} "\$inputs_dir/" || true
    fi
    if [ "$tree" != "NO_TREE" ]; then
        cp $tree "\$inputs_dir/" || true
    fi
    if [ "$species_map" != "NO_SPECIES_MAP" ]; then
        cp $species_map "\$inputs_dir/" || true
    fi

    plot_synteny.py \\
        --query_bed $query_bed \\
        --home_bed $home_bed \\
        --home_gff $home_gff \\
        $gffs_arg \\
        $names_arg \\
        $cands_arg \\
        $homo_arg \\
        --tree $tree \\
        $species_arg \\
        --output ${home_bed.baseName}_synteny_plot.html
    """
}
