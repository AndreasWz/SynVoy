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

    output:
    path "*_synteny_plot.html", emit: plot

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
    
    """
    plot_synteny.py \\
        --query_bed $query_bed \\
        --home_bed $home_bed \\
        --home_gff $home_gff \\
        $gffs_arg \\
        $names_arg \\
        $cands_arg \\
        $homo_arg \\
        --tree $tree \\
        --output ${home_bed.baseName}_synteny_plot.html
    """
}
