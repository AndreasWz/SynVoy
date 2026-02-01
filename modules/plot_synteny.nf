process PLOT_SYNTENY {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path home_bed
    path target_gffs
    val target_names
    path candidate_beds
    path homology_tsvs

    output:
    path "*_synteny_plot.html", emit: plot

    script:
    // target_names is a list like ['genome1', 'genome2']
    // We need to pass them as space-separated string
    def names_str = target_names.join(' ')
    def cands_str = candidate_beds.join(' ')
    def homo_str = homology_tsvs.join(' ')
    def out_name = "${home_bed.baseName}_synteny_plot.html"
    
    """
    plot_synteny.py \\
        --home_bed $home_bed \\
        --target_gffs $target_gffs \\
        --target_names $names_str \\
        --candidate_beds $cands_str \\
        --homology_tsvs $homo_str \\
        --output $out_name
    """
}
