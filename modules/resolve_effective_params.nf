// Resolves preset-affected params at workflow runtime — docs/TODO.md §1f, and since
// 2026-10-15 the preset modules (bin/select_modules.py).
//
// LOCATE_GENE emits an advisory auto_preset.json and the hit profile (with the query's
// own features), ANNOTATE_GOI the home gene model, HOME_CONTEXT the home neighbourhood.
// This process merges them with the defaults from nextflow.config (passed as a JSON
// value), the user's preset_override / disable-auto knobs and --modules, then emits a
// single effective_params.json that downstream processes consume as a Map via
// `val settings` (see modules/iterative_search.nf, cluster_regions.nf,
// extract_flanking.nf), and that the report reads as its run_settings.
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

    output:
    path "effective_params.json", emit: json

    script:
    def override_arg = preset_override ? "--preset_override ${preset_override}" : ''
    def disable_arg  = disable_auto_apply ? '--disable_auto_apply' : ''
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
        --defaults defaults.json \\
        --shipped_config ${projectDir}/nextflow.config \\
        ${override_arg} \\
        ${disable_arg} \\
        --output effective_params.json
    """
}
