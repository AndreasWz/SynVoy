// The genome axis (profile stage 5). The preset modules run once per query and see only the home
// side, so a setting that should differ per target genome is set inside the search, where
// bin/iterative_search_runner.py derives each genome's distances (--genome_axis there). On with
// any --modules other than legacy. Under legacy the search gets the arguments it had before: an
// 'auto' value is passed as its legacy value and the switch is left out, so the task, its hash
// and every file stay as they were. A value the user set is passed as given and wins for every
// genome.
//   anchor_max_places  'auto': every place (0), unless the block cap binds in a genome, which is
//                      then rebuilt from the 5 best places of each flanking gene (axis on only)
class GenomeAxis {

    // bin/select_modules.parse_spec: an empty --modules is legacy, else its first word is the mode
    static boolean on(params) {
        def text = (params.modules == null ? '' : params.modules).toString().trim()
        def first = text.split(',').collect { p -> p.trim() }.find { p -> p } ?: 'legacy'
        return first.toLowerCase() != 'legacy'
    }

    static Object anchorMaxPlaces(params) {
        if (params.anchor_max_places.toString().trim() != 'auto') {
            return params.anchor_max_places
        }
        return on(params) ? 'auto' : 0
    }

    // the runner's switch, appended to an existing argument line: nothing under legacy
    static String args(params) {
        return on(params) ? ' --genome_axis true' : ''
    }
}
