// Genetic codes of a run (preset plan of 2026-10-15, PR-3). The home genome has one NCBI
// table (--home_genetic_code); each target genome has the table named for it in
// --target_genetic_codes (a TSV: genome<TAB>table, the genome as its file name with or
// without the extension) or else --default_target_genetic_code. Every one is 1 by default,
// and then no script or engine gets a genetic-code option at all: a default run's commands
// do not change. Shared by the modules that translate; kept here (lib/ is on the
// pipeline's class path) so the TSV is read the same way everywhere.
class GeneticCodes {

    // The tables every engine of the search translates with (BLAST+ 2.17, MMseqs2 18,
    // miniprot 0.18); the same list as sequence_utils.SUPPORTED_GENETIC_CODES.
    static final List<Integer> SUPPORTED = (1..6) + (9..16) + (21..31)

    // [[genome, table], ...] of the TSV; null or '' gives []. Blank lines, lines starting
    // with # and a header whose second column is not a number are skipped.
    static List<List> parse(Object tsv) {
        if (!tsv || !tsv.toString().trim()) {
            return []
        }
        def rows = []
        new File(tsv.toString()).readLines().each { line ->
            def text = line.trim()
            if (!text || text.startsWith('#')) {
                return
            }
            def f = text.split(/\s+/)
            if (f.size() >= 2 && f[1] ==~ /\d+/) {
                rows << [f[0], f[1].toInteger()]
            }
        }
        return rows
    }

    // The options of a script that translates target genomes; '' when every target uses
    // table 1. The scripts pick each genome's table from 'genome=table,...'.
    static String targetArgs(Map params) {
        def spec = parse(params.target_genetic_codes).collect { r -> "${r[0]}=${r[1]}" }.join(',')
        def dflt = params.default_target_genetic_code.toString().toInteger()
        if (!spec && dflt == 1) {
            return ''
        }
        return "--target_genetic_codes '${spec}' --default_target_genetic_code ${dflt}"
    }

    // The option of a script that translates the home genome; '' for table 1.
    static String homeArg(Map params) {
        def n = params.home_genetic_code.toString().toInteger()
        return n != 1 ? "--genetic_code ${n}" : ''
    }

    // What is wrong with the three parameters (checked before any task starts).
    static List<String> errors(Map params) {
        def out = []
        ['home_genetic_code', 'default_target_genetic_code'].each { name ->
            def text = (params[name] == null ? '1' : params[name]).toString().trim()
            if (!(text ==~ /\d+/) || !(text.toInteger() in SUPPORTED)) {
                out << "Invalid --${name} '${text}': an NCBI genetic code that every search engine supports (${SUPPORTED.join(', ')})."
            }
        }
        def tsv = params.target_genetic_codes
        if (tsv && tsv.toString().trim()) {
            if (!new File(tsv.toString()).exists()) {
                out << "--target_genetic_codes: file not found: ${tsv}"
            } else {
                parse(tsv).each { r ->
                    if (!(r[1] in SUPPORTED)) {
                        out << "--target_genetic_codes: genome ${r[0]} has table ${r[1]}, which not every search engine supports (${SUPPORTED.join(', ')})."
                    }
                    if (r[0] ==~ /.*[=,'"].*/) {
                        out << "--target_genetic_codes: genome name '${r[0]}' contains one of = , ' \" ."
                    }
                }
            }
        }
        return out
    }
}
