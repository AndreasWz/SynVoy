// Gene prediction for a home genome without a GFF (profile stage 5c). augustus_species 'auto' (the
// default): under --modules legacy Augustus's 'fly' model for every home, as before; with any other
// --modules (GenomeAxis.on) the installed model nearest the home species, which
// bin/resolve_effective_params.py chooses for a home without a GFF and records in the settings and
// the derivation of effective_params.json. A model the user names is used as given.
class HomeGenes {

    static Object augustusSpecies(params, settings) {
        def chosen = settings instanceof Map ? settings.augustus_species : null
        if (chosen) {
            return chosen
        }
        return params.augustus_species.toString().trim() == 'auto' ? 'fly' : params.augustus_species
    }

    // RESOLVE_EFFECTIVE_PARAMS's arguments for the choice, appended to an existing argument line:
    // nothing under legacy, so that task stays as it was
    static String resolverArgs(params, homeSpecies, homeGffPresent) {
        if (!GenomeAxis.on(params)) {
            return ''
        }
        def species = homeSpecies.toString().replace("'", "'\"'\"'")
        return " --home_gff_present ${homeGffPresent} --home_species '${species}'" +
               " --augustus_species '${params.augustus_species}' --gene_predictor ${params.gene_predictor}"
    }
}
