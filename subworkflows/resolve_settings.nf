// The settings of a run: the old preset or the preset modules, merged with the defaults
// and the values the user set (bin/resolve_effective_params.py, bin/select_modules.py).
// A sub-workflow of its own because the main `workflow {}` body is at the JVM method-size
// limit. Emits the parsed settings map (a value channel: every EXTRACT_FLANKING /
// ITERATIVE_SEARCH / CLUSTER_REGIONS task reads it) and effective_params.json itself,
// which the report reads as its run_settings.

include { HOME_CONTEXT } from '../modules/home_context.nf'
include { RESOLVE_EFFECTIVE_PARAMS } from '../modules/resolve_effective_params.nf'

workflow RESOLVE_SETTINGS {

    take:
    query           // value channel: the normalized query FASTA
    home_genome     // home genome FASTA
    home_gff        // value channel: home GFF or the NO_GFF sentinel
    goi_bed         // ANNOTATE_GOI.out.bed
    goi_info        // ANNOTATE_GOI.out.info
    hit_profile     // LOCATE_GENE.out.hit_profile (optional)
    auto_preset     // LOCATE_GENE.out.auto_preset (optional)
    defaults_json   // val: presetDefaultKeys() of params, as JSON
    preset_override // val
    disable_auto    // val: auto_apply_preset = false
    modules_spec    // val: --modules

    main:
    HOME_CONTEXT(query, home_genome, home_gff, goi_bed)

    RESOLVE_EFFECTIVE_PARAMS(
        auto_preset.ifEmpty(file("${projectDir}/assets/sentinels/NO_AUTO_PRESET")),
        hit_profile.ifEmpty(file("${projectDir}/assets/sentinels/NO_HIT_PROFILE")),
        goi_info,
        HOME_CONTEXT.out.json.ifEmpty(file("${projectDir}/assets/sentinels/NO_HOME_CONTEXT")),
        defaults_json,
        preset_override,
        disable_auto,
        modules_spec
    )

    settings = RESOLVE_EFFECTIVE_PARAMS.out.json
        .map { f ->
            def parsed = new groovy.json.JsonSlurper().parse(f)
            // Surface the chosen preset and modules in the console banner for auditability.
            log.info "[INFO] §1f preset applied: ${parsed.preset_applied ?: 'none'} (${parsed.apply_reason})"
            def fmt = { list -> list ? list.collect { m -> "${m.id} ${m.name} (${m.reason})" }.join(', ') : 'none' }
            def shown = parsed.modules_mode == 'legacy' ? 'triggers firing (recorded only)' : 'applied'
            def mods = parsed.modules_mode == 'legacy' ? parsed.modules_would_fire : parsed.modules_applied
            log.info "[INFO] modules (${parsed.modules_mode}), ${shown}: ${fmt.call(mods)}"
            // What each module rule set and from what (select_modules.RULES); none under legacy.
            def ruled = (parsed.derivation ?: [:]).findAll { k, d -> (parsed.source ?: [:])[k]?.toString()?.startsWith('module:') }
            if (ruled) {
                log.info "[INFO] settings from the module rules: ${ruled.collect { k, d -> d.line }.sort().join('; ')}"
            }
            (parsed.constraints ?: []).findAll { con -> con.applies && !con.ok }.each { con ->
                log.warn "constraint ${con.id} (${con.doc}): ${con.detail}; kept, a value you set or a default is involved"
            }
            parsed.settings as Map
        }
        .first()

    emit:
    settings = settings
    json = RESOLVE_EFFECTIVE_PARAMS.out.json.first()
}
