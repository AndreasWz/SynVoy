// Resolves preset-affected params at workflow runtime — docs/TODO.md §1f.
//
// LOCATE_GENE emits an advisory auto_preset.json. This process merges that
// advice with the defaults from nextflow.config (passed as a JSON value) and
// the user's preset_override / disable-auto knobs, then emits a single
// effective_params.json that downstream processes consume as a Map via
// `val settings` (see modules/iterative_search.nf, cluster_regions.nf,
// extract_flanking.nf).
//
// We intentionally keep this process tiny + cheap (cache 'deep'): the
// resolution logic lives in bin/resolve_effective_params.py for testability.

process RESOLVE_EFFECTIVE_PARAMS {
    tag "effective_params"
    label 'process_low'
    cache 'deep'
    publishDir "${params.outdir}/intermediate/locate_gene", mode: 'copy'

    input:
    path auto_preset
    val defaults_json
    val preset_override
    val disable_auto_apply

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
        --defaults defaults.json \\
        ${override_arg} \\
        ${disable_arg} \\
        --output effective_params.json
    """
}
