// Resolves preset-affected params at workflow runtime — docs/TODO.md §1f, and since
// 2026-10-15 the preset modules (bin/select_modules.py).
//
// LOCATE_GENE emits an advisory auto_preset.json and the hit profile (with the query's
// own features), ANNOTATE_GOI the home gene model, HOME_CONTEXT the home neighbourhood.
// This process merges them with the defaults from nextflow.config (passed as a JSON
// value), the user's preset_override / disable-auto knobs and --modules, then emits a
// single effective_params.json that downstream processes consume as a Map via
// `val settings` (see modules/iterative_search.nf, cluster_regions.nf,
// extract_flanking.nf), and that the report reads as its run_settings. It also holds the
// profile the module rules read, the derivation of every key a module rule, the preset or
// the user set, and the constraint checks; a violation that only the rules made exits 2
// (the run stops before the search). With any --modules other than legacy and a home genome
// without a GFF it also chooses the Augustus model of the home's gene prediction
// (augustus_species, stage 5c), and stops the run when a eukaryotic home would fall back to
// Prodigal.
//
// We intentionally keep this process tiny + cheap (cache 'deep'): the
// resolution logic lives in bin/resolve_effective_params.py for testability.
// Each optional input arrives as its own sentinel file when missing (NO_AUTO_PRESET,
// NO_HIT_PROFILE, NO_HOME_CONTEXT): two inputs with one name would collide in the task
// folder, and an empty channel would leave this process, and every step waiting for
// `settings`, never started.

process RESOLVE_EFFECTIVE_PARAMS {
    tag "effective_params"
    label 'process_low'
    cache 'deep'
    publishDir "${params.outdir}/intermediate/locate_gene", mode: 'copy'

    input:
    path auto_preset
    path hit_profile
    path goi_info
    path home_context
    val defaults_json
    val preset_override
    val disable_auto_apply
    val modules_spec
    val home_species        // for the Augustus model of a home without a GFF (stage 5c)
    val home_gff_present

    output:
    path "effective_params.json", emit: json

    script:
    def override_arg = preset_override ? "--preset_override ${preset_override}" : ''
    // The run's resources as it used them (fix-phase d6): a profile (laptop_safe, low_mem,
    // docker_max) changes these and effective_params.json did not say so.
    def resources = [nextflow_profile: workflow.profile, executor: task.executor,
                     nextflow_version: workflow.nextflow.version,
                     mmseqs_sensitivity: params.mmseqs_sensitivity,
                     mmseqs_split_memory_limit: params.mmseqs_split_memory_limit,
                     iterative_search_cpus: params.iterative_search_cpus,
                     iterative_search_memory: params.iterative_search_memory,
                     iterative_search_max_forks: params.iterative_search_max_forks,
                     locate_gene_cpus: params.locate_gene_cpus, locate_gene_memory: params.locate_gene_memory,
                     compute_tree_cpus: params.compute_tree_cpus, skip_tree: params.skip_tree.toString().toBoolean()]
    def resources_arg = resources.collect { k, v -> '--run_resource ' + ShellArgs.q("${k}=${v}") }.join(' ')
    def disable_arg  = disable_auto_apply ? '--disable_auto_apply' : ''
    // lib/HomeGenes.groovy: '' under --modules legacy, so that task is the one it was before.
    def home_args = HomeGenes.resolverArgs(params, home_species, home_gff_present)
    """
    # The defaults map comes in as a Groovy/JSON string. Write it to a sidecar
    # so we don't have to worry about quoting in the python invocation.
    cat > defaults.json <<'__SYNVOY_DEFAULTS__'
${defaults_json}
__SYNVOY_DEFAULTS__

    ${projectDir}/bin/resolve_effective_params.py \\
        --auto_preset ${auto_preset} \\
        --hit_profile ${hit_profile} \\
        --goi_info ${goi_info} \\
        --home_context ${home_context} \\
        --modules '${modules_spec}' \\
        --n_flanking_genes ${params.n_flanking_genes} \\
        --window_margin_genes ${params.window_margin_genes}${home_args} \\
        --defaults defaults.json \\
        --shipped_config ${projectDir}/nextflow.config \\
        ${resources_arg} \\
        --tool_versions \\
        ${override_arg} \\
        ${disable_arg} \\
        --output effective_params.json
    """
}
