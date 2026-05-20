#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// Import Modules
include { LOCATE_GENE } from './modules/locate_gene.nf'
include { SPLIT_LOCI } from './modules/split_loci.nf'
include { EXTRACT_FLANKING } from './modules/extract_flanking.nf'
include { PREPARE_INITIAL_DB } from './modules/prepare_initial_db.nf'
include { ITERATIVE_SEARCH } from './modules/iterative_search.nf'
include { CLUSTER_REGIONS } from './modules/cluster_regions.nf'
include { PREPARE_HOME_PROTEOME } from './modules/prepare_home.nf'
include { PLOT_SYNTENY } from './modules/plot_synteny.nf'
include { COMPUTE_TREE } from './modules/compute_tree.nf'
include { ANNOTATE_GOI } from './modules/annotate_goi.nf'

// New Modules
include { STAGE_GENOMES } from './modules/stage_genomes.nf'
include { ASSESS_GENOME_QUALITY } from './modules/assess_quality.nf'
include { FETCH_QUERY_FROM_ID } from './modules/fetch_query.nf'
include { RESOLVE_GENE_INPUT } from './modules/resolve_query.nf'
include { PHYLO_SORT } from './modules/phylo_sort.nf'
include { FETCH_RELATED_GENOMES } from './modules/fetch_related.nf'
include { FETCH_HOME_GENOME } from './modules/fetch_home.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'
include { BORROW_ANNOTATIONS } from './modules/borrow_annotations.nf'
include { NORMALIZE_QUERY } from './modules/normalize_query.nf'
include { FILTER_SORTED_GENOMES } from './modules/filter_targets.nf'
include { ESTIMATE_PARAMS } from './modules/estimate_params.nf'

// ==============================================================================
// Console UI helpers
// ==============================================================================

// ANSI color codes as a function (NF 26.x forbids top-level variable declarations)
def colors() {
    return [
        reset:  "\033[0m",
        bold:   "\033[1m",
        dim:    "\033[2m",
        black:  "\033[0;30m",
        red:    "\033[0;31m",
        green:  "\033[0;32m",
        yellow: "\033[0;33m",
        blue:   "\033[0;34m",
        purple: "\033[0;35m",
        cyan:   "\033[0;36m",
        white:  "\033[0;37m"
    ]
}

def uiRule() {
    def c = colors()
    return "${c.blue}${'═' * 63}${c.reset}"
}

def uiStatus(String level, String task, String detail = '') {
    def c = colors()
    def levelColors = [
        'RUN ': c.blue,
        'OK  ': c.green,
        'INFO': c.cyan,
        'WARN': c.yellow,
        'SKIP': c.dim,
        'FAIL': c.red
    ]
    def key = (level ?: 'INFO').padRight(4).substring(0, 4)
    def levelColor = levelColors.get(key, c.white)
    def taskCol = task ? "${c.white}${task.padRight(24)}${c.reset}" : ''
    def detailCol = detail ?: ''
    def prefix = "${levelColor}[${key}]${c.reset}"
    if (taskCol && detailCol) {
        log.info "${prefix} ${taskCol} ${detailCol}"
    } else if (taskCol) {
        log.info "${prefix} ${taskCol}"
    } else if (detailCol) {
        log.info "${prefix} ${detailCol}"
    } else {
        log.info prefix
    }
}

def uiPhase(int idx, String title) {
    def c = colors()
    log.info ""
    log.info uiRule()
    log.info "${c.white}Phase ${idx}${c.reset} ${c.dim}|${c.reset} ${c.cyan}${title}${c.reset}"
    log.info uiRule()
}

def printHeader() {
    def c = colors()
    log.info ""
    log.info uiRule()
    log.info "${c.cyan}${c.bold}SynVoy${c.reset} ${c.dim}v2.0${c.reset}"
    log.info "${c.dim}Phylogenetically-informed syntenic ortholog discovery${c.reset}"
    log.info uiRule()
}

def looksLikeInlineSequence(value) {
    if (!value) {
        return false
    }
    def text = value.toString().trim()
    if (!text) {
        return false
    }
    if (text.contains('\n') || text.contains('\r')) {
        return true
    }
    if (text.startsWith('>')) {
        return true
    }
    if (text.size() >= 30 && (text ==~ /[A-Za-z\\*\\-]+/)) {
        return true
    }
    return false
}

def paramBool(value) {
    if (value instanceof Boolean) {
        return value
    }
    if (value == null) {
        return false
    }
    def text = value.toString().trim().toLowerCase()
    return text in ['true', '1', 'yes', 'y', 'on']
}

def flattenNestedList(value) {
    if (!(value instanceof List)) {
        return [value]
    }
    def out = []
    value.each { item ->
        out.addAll(flattenNestedList(item))
    }
    return out
}

def normalizeCombineRecord(record, int leftWidth) {
    if (!(record instanceof List)) {
        return null
    }
    // Nextflow has changed combine tuple flattening behavior across versions.
    // Accept both [left_tuple, right] and [left_tuple..., right] forms.
    if (record.size() == 2 && record[0] instanceof List) {
        def left = record[0]
        if (left.size() < leftWidth) {
            return null
        }
        return flattenNestedList(left) + [record[1]]
    }
    if (record.size() >= leftWidth + 1) {
        return record
    }
    return null
}

def printParams() {
    def c = colors()
    def inline_query = params.query_seq ?: (looksLikeInlineSequence(params.query_id) ? params.query_id : null)
    def query_display = 'N/A'
    if (params.mode == 'easy') {
        if (inline_query) {
            query_display = 'inline_sequence'
        } else if (params.query) {
            query_display = new File(params.query).name
        } else {
            query_display = params.query_id ?: 'N/A'
        }
    } else {
        query_display = params.query ? new File(params.query).name : 'N/A'
    }
    def home_display = params.mode == 'easy' ? params.home_species : (params.home_genome ? new File(params.home_genome).name : 'N/A')
    def target_display = params.target_species ?: 'auto (taxonomic search)'
    
    log.info uiRule()
    log.info "${c.white}Run Configuration${c.reset}"
    log.info uiRule()
    log.info "${c.dim}Query Gene     ${c.reset} ${c.green}${query_display}${c.reset}"
    log.info "${c.dim}Home Genome    ${c.reset} ${c.green}${home_display}${c.reset}"
    log.info "${c.dim}Mode           ${c.reset} ${c.yellow}${params.mode}${c.reset}"
    log.info "${c.dim}Target Species ${c.reset} ${c.cyan}${target_display}${c.reset}"
    log.info "${c.dim}Flanking Genes ${c.reset} ${params.n_flanking_genes}"
    log.info "${c.dim}MMseqs Sens.   ${c.reset} ${params.mmseqs_sensitivity}"
    if (params.mode == 'easy') {
        log.info "${c.dim}Asm Ranking    ${c.reset} ${params.assembly_ranking}"
        log.info "${c.dim}LowQ Policy    ${c.reset} ${params.bad_quality_policy} (timeout=${params.bad_quality_timeout}s)"
    }
    log.info "${c.dim}Output Dir     ${c.reset} ${params.outdir}"
    log.info uiRule()
}

workflow {
    // Color codes for logging (defined once, used throughout workflow)
    def c = colors()

    // Print header and params at workflow start
    printHeader()
    printParams()

    // Stable sentinel files for optional path inputs.
    def no_gff_file = file("${projectDir}/assets/sentinels/NO_GFF")
    def no_species_map_file = file("${projectDir}/assets/sentinels/NO_SPECIES_MAP")
    def inline_query_mode = false

    log.info ""
    
    // ========== INPUT VALIDATION ==========
    // Collect ALL errors so the user sees every problem at once, not one at a time.
    def validationErrors = []
    def validationWarnings = []

    // --- Mode validation ---
    if (!(params.mode in ['easy', 'pro'])) {
        validationErrors << "Invalid --mode '${params.mode}'. Must be 'easy' or 'pro'."
    }

    // --- Mode-specific input validation ---
    if (params.mode == 'easy') {
        def query_id_is_inline = looksLikeInlineSequence(params.query_id)
        def has_query = params.query_id || params.query || params.query_seq
        if (!has_query) {
            validationErrors << "Easy mode requires at least one of: --query_id, --query, or --query_seq"
        }
        // All three supplied → ambiguous intent, refuse rather than silently pick one.
        // Allow --query_id when it's inline-sequence syntax (legacy path) paired with --query_seq.
        if (params.query_id && params.query && params.query_seq && !query_id_is_inline) {
            validationErrors << "Easy mode: --query_id, --query, and --query_seq are all set — choose one. Priority order if you insist on two: --query_seq > --query > --query_id."
        }
        if ((params.query_seq || query_id_is_inline) && !params.home_species) {
            validationErrors << "Easy mode with inline sequence requires --home_species (cannot auto-detect from raw sequence)"
        }
        if (params.query_seq && params.query) {
            validationWarnings << "Both --query_seq and --query provided; using --query_seq."
        } else if (params.query && params.query_id && !params.query_seq && !query_id_is_inline) {
            validationWarnings << "Both --query and --query_id provided; using --query."
        } else if (params.query_id && query_id_is_inline) {
            validationWarnings << "Inline sequence detected in --query_id; treating as inline FASTA input."
        }

        // Easy-mode enum checks
        if (!(params.assembly_ranking in ['hybrid', 'counts', 'nstats'])) {
            validationErrors << "Invalid --assembly_ranking '${params.assembly_ranking}'. Must be: hybrid, counts, or nstats."
        }
        if (!(params.bad_quality_policy in ['ask', 'drop', 'keep'])) {
            validationErrors << "Invalid --bad_quality_policy '${params.bad_quality_policy}'. Must be: ask, drop, or keep."
        }

    } else if (params.mode == 'pro') {
        if (!params.query) {
            validationErrors << "Pro mode requires --query (path to query FASTA file)"
        } else if (!file(params.query).exists()) {
            validationErrors << "Query FASTA not found: ${params.query}"
        } else {
            // Pre-flight query length check — normalize_query.py enforces this
            // at runtime, but surfacing it here saves a full pipeline startup.
            def q_seq_chars = 0
            try {
                file(params.query).readLines().each { line ->
                    if (!line.startsWith('>')) {
                        q_seq_chars += line.replaceAll('\\s', '').length()
                    }
                }
            } catch (Exception e) {
                validationWarnings << "Could not pre-read --query for length check: ${e.message}"
            }
            if (q_seq_chars > 0 && q_seq_chars < 30) {
                // Warn (not error): the file may be DNA that normalize_query.py
                // will translate, turning a 30 nt file into a ~10 aa protein, or
                // vice versa. normalize_query.py issues the hard stop.
                validationWarnings << "Query FASTA '${params.query}' contains only ${q_seq_chars} non-header chars. If protein: below SynVoy's 30 aa minimum — searches will be noisy. If DNA: fine (will be translated). normalize_query.py will enforce the final check."
            }
        }
        if (!params.home_genome) {
            validationErrors << "Pro mode requires --home_genome (path to home genome FASTA)"
        } else if (!file(params.home_genome).exists()) {
            validationErrors << "Home genome not found: ${params.home_genome}"
        }
        if (params.home_gff && !file(params.home_gff).exists()) {
            validationErrors << "Home GFF not found: ${params.home_gff}"
        } else if (!params.home_gff && params.home_genome) {
            validationWarnings << "No --home_gff provided. Flanking-gene extraction will fall back to Prodigal gene prediction on the home genome, which is substantially less reliable than a curated GFF. Supply --home_gff if you have one."
        }
        if (!params.target_genomes) {
            validationWarnings << "No --target_genomes provided; will run home-genome-only analysis (no iterative search)."
        } else {
            // Resolve the glob / list / comma-separated string up front so
            // students don't see a silent "0 target genomes" downstream.
            def tg = params.target_genomes
            def matches = []
            try {
                if (tg instanceof List) {
                    tg.each { p -> if (file(p).exists()) matches << p }
                } else if (tg.toString().contains(',')) {
                    tg.toString().split(',').collect { item -> item.trim() }.each { p ->
                        if (file(p).exists()) matches << p
                    }
                } else {
                    def resolved = file(tg.toString())
                    if (resolved instanceof List) {
                        matches = resolved
                    } else if (resolved.exists()) {
                        matches = [resolved]
                    }
                }
            } catch (Exception e) {
                validationWarnings << "Could not pre-resolve --target_genomes pattern '${tg}': ${e.message}"
            }
            if (matches.isEmpty()) {
                validationErrors << "Target genomes pattern '${tg}' matched zero files. Check the path (relative paths are resolved against the Nextflow launch dir), the glob (quote it to prevent shell expansion: --target_genomes \"path/to/*.fa\"), and that the files exist."
            }
        }
    }

    // --- Universal parameter range validation ---
    if (!(params.qc_fail_policy in ['drop', 'keep'])) {
        validationErrors << "Invalid --qc_fail_policy '${params.qc_fail_policy}'. Must be 'drop' or 'keep'."
    }
    if (params.n_flanking_genes < 1) {
        validationErrors << "Invalid --n_flanking_genes (${params.n_flanking_genes}). Must be >= 1."
    }
    if (params.min_synteny_score < 0 || params.min_synteny_score > 1) {
        validationErrors << "Invalid --min_synteny_score (${params.min_synteny_score}). Must be between 0 and 1."
    }
    if (params.min_hit_identity < 0 || params.min_hit_identity > 100) {
        validationErrors << "Invalid --min_hit_identity (${params.min_hit_identity}). Must be between 0 and 100."
    }
    if (params.min_hit_length < 1) {
        validationErrors << "Invalid --min_hit_length (${params.min_hit_length}). Must be >= 1."
    }
    if (params.search_evalue <= 0) {
        validationErrors << "Invalid --search_evalue (${params.search_evalue}). Must be > 0."
    }
    if (params.max_intron < 0) {
        validationErrors << "Invalid --max_intron (${params.max_intron}). Must be >= 0."
    }
    if (params.cluster_distance < 0 && params.cluster_distance != -1) {
        validationErrors << "Invalid --cluster_distance (${params.cluster_distance}). Must be >= 0 or -1 (auto)."
    }
    if (params.mmseqs_sensitivity < 1 || params.mmseqs_sensitivity > 12) {
        validationWarnings << "Unusual --mmseqs_sensitivity (${params.mmseqs_sensitivity}). Typical range is 1-9.5."
    }
    if (params.sw_timeout_seconds < 1) {
        validationErrors << "Invalid --sw_timeout_seconds (${params.sw_timeout_seconds}). Must be >= 1."
    }
    if (!(params.sw_method in ['auto', 'parasail', 'ssearch36'])) {
        validationErrors << "Invalid --sw_method '${params.sw_method}'. Must be: auto, parasail, or ssearch36."
    }
    if (params.sw_min_identity < 0 || params.sw_min_identity > 100) {
        validationErrors << "Invalid --sw_min_identity (${params.sw_min_identity}). Must be between 0 and 100."
    }
    if (params.region_padding < 0 || params.padding_min < 0 || params.padding_max < 0) {
        validationErrors << "Padding values (--region_padding, --padding_min, --padding_max) must be >= 0."
    }
    if (params.padding_max < params.padding_min) {
        validationErrors << "Invalid padding: --padding_max (${params.padding_max}) must be >= --padding_min (${params.padding_min})."
    }
    if (params.max_blocks_per_genome < 0) {
        validationErrors << "Invalid --max_blocks_per_genome (${params.max_blocks_per_genome}). Must be >= 0."
    }
    if (params.min_block_genes < 1) {
        validationErrors << "Invalid --min_block_genes (${params.min_block_genes}). Must be >= 1."
    }
    if (params.max_flanking_goi_similarity < 0 || params.max_flanking_goi_similarity > 100) {
        validationErrors << "Invalid --max_flanking_goi_similarity (${params.max_flanking_goi_similarity}). Must be between 0 and 100."
    }
    if (params.min_flanking_size < 0) {
        validationErrors << "Invalid --min_flanking_size (${params.min_flanking_size}). Must be >= 0."
    }
    if (params.bad_quality_timeout < 0) {
        validationErrors << "Invalid --bad_quality_timeout (${params.bad_quality_timeout}). Must be >= 0."
    }

    // Classification thresholds
    if (params.classify_high_min_identity < 0 || params.classify_high_min_identity > 100) {
        validationErrors << "Invalid --classify_high_min_identity (${params.classify_high_min_identity}). Must be between 0 and 100."
    }
    if (params.classify_medium_min_identity < 0 || params.classify_medium_min_identity > 100) {
        validationErrors << "Invalid --classify_medium_min_identity (${params.classify_medium_min_identity}). Must be between 0 and 100."
    }
    if (params.classify_tandem_min_identity < 0 || params.classify_tandem_min_identity > 100) {
        validationErrors << "Invalid --classify_tandem_min_identity (${params.classify_tandem_min_identity}). Must be between 0 and 100."
    }
    if (params.classify_fragment_max_qcov < 0 || params.classify_fragment_max_qcov > 1) {
        validationErrors << "Invalid --classify_fragment_max_qcov (${params.classify_fragment_max_qcov}). Must be between 0 and 1."
    }
    if (params.classify_complete_min_qcov < 0 || params.classify_complete_min_qcov > 1) {
        validationErrors << "Invalid --classify_complete_min_qcov (${params.classify_complete_min_qcov}). Must be between 0 and 1."
    }
    if (params.classify_fragment_max_qcov >= params.classify_complete_min_qcov) {
        validationWarnings << "Unusual thresholds: --classify_fragment_max_qcov (${params.classify_fragment_max_qcov}) >= --classify_complete_min_qcov (${params.classify_complete_min_qcov}). Fragment and complete ranges overlap."
    }

    // Gene predictor validation
    if (!(params.gene_predictor in ['auto', 'augustus', 'prodigal'])) {
        validationErrors << "Invalid --gene_predictor (${params.gene_predictor}). Must be 'auto', 'augustus', or 'prodigal'."
    }

    // PLM embedding search thresholds
    if (paramBool(params.enable_plm_search)) {
        if (!(params.plm_device in ['cpu', 'cuda'])) {
            validationErrors << "Invalid --plm_device (${params.plm_device}). Must be 'cpu' or 'cuda'."
        }
        if (params.plm_similarity_threshold < 0 || params.plm_similarity_threshold > 1) {
            validationErrors << "Invalid --plm_similarity_threshold (${params.plm_similarity_threshold}). Must be between 0 and 1."
        }
        if (params.plm_medium_threshold < 0 || params.plm_medium_threshold > 1) {
            validationErrors << "Invalid --plm_medium_threshold (${params.plm_medium_threshold}). Must be between 0 and 1."
        }
        if (params.plm_high_threshold < 0 || params.plm_high_threshold > 1) {
            validationErrors << "Invalid --plm_high_threshold (${params.plm_high_threshold}). Must be between 0 and 1."
        }
        if (params.plm_medium_threshold >= params.plm_high_threshold) {
            validationWarnings << "PLM thresholds: --plm_medium_threshold (${params.plm_medium_threshold}) >= --plm_high_threshold (${params.plm_high_threshold})."
        }
    }

    // Structural search (ESMFold + Foldseek) thresholds
    if (paramBool(params.enable_structural_search)) {
        if (!(params.structural_device in ['cpu', 'cuda'])) {
            validationErrors << "Invalid --structural_device (${params.structural_device}). Must be 'cpu' or 'cuda'."
        }
        if (params.structural_tm_threshold < 0 || params.structural_tm_threshold > 1) {
            validationErrors << "Invalid --structural_tm_threshold (${params.structural_tm_threshold}). Must be between 0 and 1."
        }
        if (params.structural_medium_threshold < 0 || params.structural_medium_threshold > 1) {
            validationErrors << "Invalid --structural_medium_threshold (${params.structural_medium_threshold}). Must be between 0 and 1."
        }
        if (params.structural_high_threshold < 0 || params.structural_high_threshold > 1) {
            validationErrors << "Invalid --structural_high_threshold (${params.structural_high_threshold}). Must be between 0 and 1."
        }
        if (params.structural_medium_threshold >= params.structural_high_threshold) {
            validationWarnings << "Structural thresholds: --structural_medium_threshold (${params.structural_medium_threshold}) >= --structural_high_threshold (${params.structural_high_threshold})."
        }
        if (params.structural_max_length < 10) {
            validationErrors << "Invalid --structural_max_length (${params.structural_max_length}). Must be >= 10."
        }
    }

    // Synteny scoring weights should sum to ~1
    def weightSum = (params.synteny_weight_base ?: 0) + (params.synteny_weight_consistency ?: 0) + (params.synteny_weight_strand ?: 0)
    if (Math.abs(weightSum - 1.0) > 0.01) {
        validationWarnings << "Synteny score weights sum to ${weightSum} (expected ~1.0). Scoring may behave unexpectedly."
    }

    // --- Print all warnings ---
    validationWarnings.each { msg ->
        uiStatus('WARN', 'INPUT', msg)
    }

    // --- Print all errors and abort ---
    if (validationErrors) {
        log.info ""
        log.info "${c.red}${c.bold}Input validation failed with ${validationErrors.size()} error(s):${c.reset}"
        validationErrors.eachWithIndex { msg, idx ->
            log.info "${c.red}  ${idx + 1}. ${msg}${c.reset}"
        }
        log.info ""
        log.info "${c.dim}Run with --help or see USAGE.md for parameter documentation.${c.reset}"
        exit 1
    }

    uiStatus('OK', 'INPUT', 'All parameters validated')
    
    // Channel setup — depends on mode
    if (params.mode == 'easy') {
        // ID/symbol/file/inline mode: resolve input via resolver process.
        def query_id_is_inline = looksLikeInlineSequence(params.query_id)
        def gene_input = params.query_id
        inline_query_mode = false
        if (params.query_seq) {
            gene_input = params.query_seq
            inline_query_mode = true
        } else if (params.query) {
            gene_input = params.query
        } else if (query_id_is_inline) {
            gene_input = params.query_id
            inline_query_mode = true
        }
        def species_override = params.home_species ?: ''
        
        def resolve_label = inline_query_mode ? 'Easy mode: resolving inline FASTA input' : (params.query ? 'Easy mode: resolving local FASTA' : 'Easy mode: resolving query_id input')
        uiStatus('RUN ', 'RESOLVE_QUERY', resolve_label)
        RESOLVE_GENE_INPUT(gene_input, species_override, inline_query_mode)
        
        // Use resolved FASTA as query
        raw_gene_ch = RESOLVE_GENE_INPUT.out.fasta
        
        // Get resolved species (auto-detected from ID, or user-provided)
        resolved_species_ch = RESOLVE_GENE_INPUT.out.species.map { species_file -> species_file.text.trim() }
        
        // Determine home species: user-provided takes priority, else auto-detected
        home_species_ch = resolved_species_ch.map { resolved ->
            def species = params.home_species ?: resolved
            if (!species) {
                log.error "${c.red}Could not detect species. Please provide --home_species${c.reset}"
                exit 1
            }
            return species
        }
        
        // Fetch home genome automatically for easy mode
        uiStatus('RUN ', 'FETCH_HOME', 'Downloading home genome from NCBI')
        FETCH_HOME_GENOME(home_species_ch)
        home_genome_ch = FETCH_HOME_GENOME.out.genome
        FETCH_HOME_GENOME.out.genome.view { genome ->
            def sizeMb = String.format("%.1f", genome.size() / (1024.0 * 1024.0))
            "${c.green}[OK  ]${c.reset} ${c.white}${'FETCH_HOME'.padRight(24)}${c.reset} home genome ready (${sizeMb} MB)"
        }
        // Use GFF if available, otherwise mark as missing
        home_gff_ch = FETCH_HOME_GENOME.out.gff.ifEmpty(no_gff_file)
        
        // Fetch related genomes for easy mode
        def max_genomes = (params.max_genomes == null ? 10 : params.max_genomes as Integer)
        if (max_genomes < 3) {
            log.warn("max_genomes=${max_genomes}: synteny scoring derives signal from consensus across species; with <3 target genomes, fallback GOI calls tend to be classified as 'ambiguous' (no multi-genome conservation evidence). Consider raising max_genomes to >=3.")
        }
        def target_species = params.target_species ?: ''
        uiStatus('RUN ', 'FETCH_RELATED', "Fetching related genomes${target_species ? ' (user-specified species)' : ' (auto-detect from taxonomy)'}")
        FETCH_RELATED_GENOMES(home_species_ch, max_genomes, target_species)
        genomes_dir_ch = FETCH_RELATED_GENOMES.out.genomes_dir
        FETCH_RELATED_GENOMES.out.genomes_dir.view { dir ->
            def count = new File(dir.toString()).listFiles()?.findAll { genome_file -> genome_file.name.endsWith('.fna') || genome_file.name.endsWith('.fna.gz') || genome_file.name.endsWith('.fa') || genome_file.name.endsWith('.fasta') }?.size() ?: 0
            "${c.green}[OK  ]${c.reset} ${c.white}${'FETCH_RELATED'.padRight(24)}${c.reset} downloaded ${count} target genome(s)"
        }
        species_map_ch = FETCH_RELATED_GENOMES.out.species_map.first()
        // Species name for phylogenetic sorting
        home_species_for_sort_ch = home_species_ch
        
    } else {
        // --- Pro mode: User provides files directly ---
        
        // 1. Query Setup
        raw_gene_ch = channel.fromPath(params.query)
        
        // 2. Home Genome Setup
        home_genome_ch = channel.fromPath(params.home_genome)
        
        if (params.home_gff) {
            home_gff_ch = channel.value(file(params.home_gff, checkIfExists: true))
        } else {
            home_gff_ch = channel.value(no_gff_file)
        }
        
        // 3. Target Genomes Setup
        if (params.target_genomes) {
            uiStatus('RUN ', 'STAGE_GENOMES', 'Loading target genomes list')
            // Support both glob patterns ("genomes/*.fna") and comma-separated
            // lists ("a.fna,b.fna,c.fna") as well as Nextflow list syntax.
            def tg = params.target_genomes
            if (tg instanceof List) {
                target_genomes_list = channel.fromPath(tg).collect()
            } else if (tg.toString().contains(',')) {
                target_genomes_list = channel
                    .fromPath(tg.toString().split(',').collect { target_path -> target_path.trim() })
                    .collect()
            } else {
                target_genomes_list = channel.fromPath(tg).collect()
            }
            
            // Show count
            target_genomes_list.view { genomes ->
                "${c.green}[OK  ]${c.reset} ${c.white}${'STAGE_GENOMES'.padRight(24)}${c.reset} staged ${genomes.size()} target genomes"
            }
            
            STAGE_GENOMES(target_genomes_list)
            genomes_dir_ch = STAGE_GENOMES.out.dir
            species_map_ch = STAGE_GENOMES.out.species_map.first()
            
        } else {
            uiStatus('WARN', 'STAGE_GENOMES', 'No target genomes provided; running home-genome-only analysis')
            genomes_dir_ch = channel.empty()
            species_map_ch = channel.value(no_species_map_file)
        }
        // Species name for phylogenetic sorting + matrix home-row label.
        // Pro mode default: take the home_genome filename and strip
        // FASTA suffixes. The full path used to leak through as the
        // species name (so the matrix row label became the file path).
        def _home_stem = ''
        if (params.home_genome) {
            _home_stem = new File(params.home_genome).name
                .replaceFirst(/\.gz$/, '')
                .replaceFirst(/\.(fna|fa|fasta)$/, '')
        }
        home_species_for_sort_ch = channel.value(params.home_species ?: _home_stem)
    }

    // Normalize query to protein space (DNA queries are translated to best ORF)
    // to keep downstream search/annotation behavior consistent.
    NORMALIZE_QUERY(raw_gene_ch)
    normalized_gene_ch = NORMALIZE_QUERY.out.fasta

    // ========== LLM PARAMETER ESTIMATION (Phase 0.5) ==========
    // When auto_params is enabled, analyze the query/species context and
    // estimate optimal pipeline parameters via Gemma 4 (or heuristic fallback).
    if (paramBool(params.auto_params)) {
        uiPhase(0, 'Automatic Parameter Estimation')
        uiStatus('RUN ', 'ESTIMATE_PARAMS', 'Estimating optimal parameters for this search')

        // Determine inputs for the estimator
        def est_home_species = params.mode == 'easy' ? home_species_ch : channel.value(params.home_species ?: '')
        def est_target_species = channel.value(params.target_species ?: '')

        // Resolved metadata JSON is available in easy mode; create a stub for pro mode
        if (params.mode == 'easy') {
            est_metadata_ch = RESOLVE_GENE_INPUT.out.metadata
        } else {
            // Create a minimal resolved_input.json for pro mode
            est_metadata_ch = normalized_gene_ch.map { fasta ->
                def meta = file("${workDir}/resolved_input_stub.json")
                meta.text = "{\"source\": \"file\", \"fasta_path\": \"${fasta}\", \"species\": \"${params.home_species ?: ''}\"}"
                return meta
            }
        }

        ESTIMATE_PARAMS(
            est_metadata_ch,
            est_home_species,
            est_target_species,
            normalized_gene_ch
        )

        // Apply estimated parameters synchronously via .map{} (not .view{})
        // .view{} is asynchronous and can race with downstream processes.
        params_applied_ch = ESTIMATE_PARAMS.out.params_json.map { json_file ->
            try {
                def est = new groovy.json.JsonSlurper().parse(json_file)
                def overrides = est.get('parameters', [:])
                def backend = est.get('backend', 'unknown')
                def count = overrides.size()
                def summary = est.get('context_summary', [:])

                // Log what was estimated
                def msg = "${count} parameter(s) estimated via ${backend}"
                if (summary) {
                    msg += " (kingdom=${summary.get('kingdom','?')}, genome=${summary.get('genome_size_mb',0)}Mb, query=${summary.get('query_length_aa',0)}aa)"
                }

                // Apply overrides to params
                // Only allow known Tier 1+2 parameters to be overridden
                def allowed = [
                    'max_intron', 'cluster_distance', 'n_flanking_genes', 'min_synteny_score',
                    'region_padding', 'padding_min', 'padding_max', 'search_evalue',
                    'min_hit_identity', 'min_hit_length', 'mmseqs_sensitivity',
                    'max_flanking_goi_similarity', 'max_flanking_distance',
                    'expand_goi_similar', 'expand_goi_similar_distance',
                    'min_gene_identity', 'enable_smith_waterman', 'sw_min_score',
                    'sw_min_identity', 'enable_plm_search', 'enable_structural_search',
                    'max_blocks_per_genome', 'min_block_genes', 'max_consecutive_empty_blocks',
                    'aug_relaxed_evalue_mult', 'gap_search_window',
                    'prefer_large_genes', 'min_flanking_size', 'exon_level_search'
                ] as Set

                // Safety lock: prevent LLM from auto-enabling advanced ML structural/PLM searches
                if (paramBool(params.force_disable_advanced_search)) {
                    allowed.remove('enable_plm_search')
                    allowed.remove('enable_structural_search')
                }

                overrides.each { key, value ->
                    if (allowed.contains(key)) {
                        def old_val = params.get(key)
                        params.put(key, value)
                        log.info "${c.dim}  [auto] ${key}: ${old_val} → ${value}${c.reset}"
                    }
                }

                // Log any warnings/issues
                def warnings = est.get('warnings', [])
                warnings.each { w -> log.info "${c.yellow}  [auto-warn] ${w}${c.reset}" }
                def issues = est.get('issues', [])
                issues.each { iss -> log.info "${c.red}  [auto-issue] ${iss}${c.reset}" }

                log.info "${c.green}[OK  ]${c.reset} ${c.white}${'ESTIMATE_PARAMS'.padRight(24)}${c.reset} ${msg}"
            } catch (Exception e) {
                log.warn "${c.yellow}[WARN]${c.reset} ${c.white}${'ESTIMATE_PARAMS'.padRight(24)}${c.reset} Could not apply estimated params: ${e.message}"
            }
            return true  // gate signal
        }
    } else {
        uiStatus('SKIP', 'ESTIMATE_PARAMS', 'Auto parameter estimation disabled (--auto_params false)')
        params_applied_ch = channel.value(true)
    }

    normalized_gene_ready_ch = normalized_gene_ch
        .combine(params_applied_ch)
        .map { fasta, _ready -> fasta }

    // PHASE 1: Core Localization
    uiPhase(1, 'Gene Localization in Home Genome')
    
    uiStatus('RUN ', 'LOCATE_GENE', 'Locating GOI in home genome')
    LOCATE_GENE(normalized_gene_ready_ch, home_genome_ch)

    LOCATE_GENE.out.bed.view { bed ->
        def lines = bed.readLines().findAll { it.trim() && !it.startsWith('#') }
        if (lines) {
            def loci = lines.collect { l -> def f = l.split('\t'); "${f[0]}:${f[1]}-${f[2]}" }.join(', ')
            "${c.green}[OK  ]${c.reset} ${c.white}${'LOCATE_GENE'.padRight(24)}${c.reset} ${lines.size()} hit(s) → ${loci}"
        } else {
            "${c.red}[WARN]${c.reset} ${c.white}${'LOCATE_GENE'.padRight(24)}${c.reset} No hits found — gene absent or below e-value threshold (${params.search_evalue}). Pipeline will abort with a diagnostic."
        }
    }

    // 4b. ANNOTATE GOI EXONS
    // Uses hits from LOCATE_GENE to annotate individual exons of the GOI
    // If GFF available: matches GOI to annotated gene and extracts CDS/exons
    // If no GFF: uses tblastn hits to detect exon boundaries (splice sites, start/stop codons)
    uiStatus('RUN ', 'ANNOTATE_GOI', 'Annotating GOI exons')
    
    // Determine query_id for name-based GFF matching
    def effective_query_id = inline_query_mode ? '' : (params.query_id ?: '')
    
    ANNOTATE_GOI(
        normalized_gene_ready_ch.first(),
        home_genome_ch,
        home_gff_ch,
        LOCATE_GENE.out.blast_hits,
        LOCATE_GENE.out.mmseqs_hits,
        effective_query_id
    )
    
    ANNOTATE_GOI.out.info.view { info ->
        "${c.green}[OK  ]${c.reset} ${c.white}${'ANNOTATE_GOI'.padRight(24)}${c.reset} GOI exon annotation complete"
    }
    
    // 5. SPLIT LOCI
    SPLIT_LOCI(LOCATE_GENE.out.bed)
    
    SPLIT_LOCI.out.beds.flatten().count().view { count ->
        "${c.green}[OK  ]${c.reset} ${c.white}${'SPLIT_LOCI'.padRight(24)}${c.reset} identified ${count} locus/loci"
    }
    
    distinct_loci_ch = SPLIT_LOCI.out.beds.flatten()
        .map { file -> tuple(file.baseName, file) }
    
    // Prepare effective home GFF (borrowed annotations + predictions) when targets exist
    def has_targets = params.target_genomes || params.mode == 'easy'
    if (has_targets) {
        home_gff_ch.branch { gff ->
            real: gff.name != 'NO_GFF'
            missing: true
        }.set { gff_status }

        // Prepare Home Proteome (Run once, used for RBH + borrowing)
        uiStatus('RUN ', 'PREPARE_HOME', 'Preparing home proteome database')
        PREPARE_HOME_PROTEOME(home_genome_ch, home_gff_ch, LOCATE_GENE.out.bed.first())
        home_proteome_db_ch = PREPARE_HOME_PROTEOME.out.db
        PREPARE_HOME_PROTEOME.out.db.view { db ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'PREPARE_HOME'.padRight(24)}${c.reset} home proteome ready"
        }

        // Borrow only when home genome has no usable GFF.
        gff_status.real.view { gff ->
            "${c.dim}[SKIP]${c.reset} ${c.white}${'BORROW_ANNOT'.padRight(24)}${c.reset} home GFF found (${gff.name})"
        }
        uiStatus('RUN ', 'BORROW_ANNOT', 'Checking targets for annotation borrowing when home GFF is missing')
        BORROW_ANNOTATIONS(
            home_genome_ch,
            PREPARE_HOME_PROTEOME.out.faa,
            genomes_dir_ch,
            LOCATE_GENE.out.bed.first(),
            params.n_flanking_genes,
            gff_status.missing.map { true }
        )
        BORROW_ANNOTATIONS.out.gff.view { gff ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'BORROW_ANNOT'.padRight(24)}${c.reset} borrowed annotations generated"
        }

        // Build effective GFF from all available sources:
        // 1. User-provided / NCBI GFF (if present and real)
        // 2. Prodigal-predicted GFF (when no annotation was available)
        // 3. Borrowed annotations from annotated target genomes
        fallback_gff_ch = PREPARE_HOME_PROTEOME.out.gff
            .mix(BORROW_ANNOTATIONS.out.gff)
            .collectFile(name: 'merged_home_annotations.gff')
            .ifEmpty(no_gff_file)
        
        effective_home_gff_ch = gff_status.real
            .concat(
                gff_status.missing.combine(fallback_gff_ch).map { pair -> pair[1] }
            )
            .first()
    } else {
        effective_home_gff_ch = home_gff_ch
    }

    // 6. Extract Flanking Genes
    uiStatus('RUN ', 'EXTRACT_FLANKING', "Extracting flanking genes (n=${params.n_flanking_genes})")
    
    EXTRACT_FLANKING(
        distinct_loci_ch,
        effective_home_gff_ch,
        home_genome_ch.first(),
        params.n_flanking_genes,
        params.min_flanking_size,
        params.prefer_large_genes,
        normalized_gene_ready_ch.first()  // GOI protein for similarity-based flanking filter
    )
    
    // 6b. CRITICAL FIX: Prepare Initial Database with GOI included
    // Combine flanking genes with query gene for iterative search
    uiStatus('RUN ', 'PREPARE_DB', 'Preparing initial database with query gene')
    
    PREPARE_INITIAL_DB(
        EXTRACT_FLANKING.out.faa,
        ANNOTATE_GOI.out.exons.first()  // .first() → value channel so it pairs with ALL loci
    )
        // Only run if we have targets
    if (params.target_genomes || params.mode == 'easy') {
        
        uiPhase(2, 'Phylogenetic Ordering and Iterative Search')
        
        EXTRACT_FLANKING.out.bed
            .map { locus_id, bed -> locus_id }
            .combine(home_genome_ch)
            .map { rec ->
                def norm = normalizeCombineRecord(rec, 1)
                norm ? tuple(norm[0], norm[1]) : null
            }
            .filter { rec -> rec != null }
            .set { phylo_sort_inputs } // [locus_id, home_genome]

        uiStatus('RUN ', 'PHYLO_SORT', 'Sorting genomes by phylogenetic distance')
        
        PHYLO_SORT(
            phylo_sort_inputs,
            genomes_dir_ch,
            home_species_for_sort_ch
        )
        
        PHYLO_SORT.out.sorted_list.view { locus, sorted ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'PHYLO_SORT'.padRight(24)}${c.reset} ordering complete for ${locus}"
        }
        
        // QC
        uiStatus('RUN ', 'GENOME_QC', "Assessing target genome quality (policy: ${params.qc_fail_policy})")

        ASSESS_GENOME_QUALITY(genomes_dir_ch)
        qc_summary_ch = ASSESS_GENOME_QUALITY.out.json

        ASSESS_GENOME_QUALITY.out.json.view { qc_json ->
            try {
                def qc = new groovy.json.JsonSlurper().parse(qc_json)
                def pass_count = qc.count { genome_qc -> genome_qc.qc_status == 'PASS' }
                def fail_count = qc.count { genome_qc -> genome_qc.qc_status == 'FAIL' }
                def total = qc.size()
                def msg = "QC complete: ${pass_count}/${total} passed"
                if (fail_count > 0) {
                    def failed_names = qc.findAll { genome_qc -> genome_qc.qc_status == 'FAIL' }.collect { genome_qc -> genome_qc.genome_id ?: genome_qc.genome ?: 'unknown' }.take(3).join(', ')
                    def suffix = fail_count > 3 ? " (+${fail_count - 3} more)" : ""
                    msg += ", ${fail_count} failed [${failed_names}${suffix}]"
                    if (params.qc_fail_policy == 'drop') msg += " (will be dropped)"
                }
                return "${c.green}[OK  ]${c.reset} ${c.white}${'GENOME_QC'.padRight(24)}${c.reset} ${msg}"
            } catch (Exception e) {
                return "${c.green}[OK  ]${c.reset} ${c.white}${'GENOME_QC'.padRight(24)}${c.reset} QC assessment complete"
            }
        }

        FILTER_SORTED_GENOMES(
            PHYLO_SORT.out.sorted_list,
            qc_summary_ch,
            params.qc_fail_policy
        )

        FILTER_SORTED_GENOMES.out.sorted_list.view { locus, sorted ->
            def count = sorted.readLines().findAll { line -> line.trim() }.size()
            "${c.green}[OK  ]${c.reset} ${c.white}${'QC_FILTER'.padRight(24)}${c.reset} ${count} target genome(s) passed QC filter for ${locus}"
        }

        // 8. Iterative Search (FOR EACH LOCUS) - Using FIXED database with GOI
        PREPARE_INITIAL_DB.out.db
            .join(FILTER_SORTED_GENOMES.out.sorted_list)
            .set { iterative_search_inputs_partial } // [locus_id, initial_db, sorted_list]
            
        iterative_search_inputs_partial
            .combine(genomes_dir_ch)
            .map { rec ->
                def norm = normalizeCombineRecord(rec, 3)
                norm ? tuple(norm[0], norm[1], norm[2], norm[3]) : null
            }
            .filter { rec -> rec != null }
            .set { iterative_search_inputs } // [locus_id, faa, sorted_list, genomes_dir]

        iterative_search_inputs
            .combine(home_proteome_db_ch)
            .map { rec ->
                def norm = normalizeCombineRecord(rec, 4)
                norm ? tuple(norm[0], norm[1], norm[2], norm[3], norm[4]) : null
            }
            .filter { rec -> rec != null }
            .set { iterative_search_final_inputs } // [locus_id, faa, sorted_list, genomes_dir, home_db]

        uiStatus('RUN ', 'ITERATIVE_SEARCH', 'Running iterative phylogenetic search')
        
        ITERATIVE_SEARCH(
            iterative_search_final_inputs.map { rec -> tuple(rec[0], rec[1]) }, // [locus, faa]
            iterative_search_final_inputs.map { rec -> rec[2] }, // sorted_list
            iterative_search_final_inputs.map { rec -> rec[3] }, // genomes_dir
            iterative_search_final_inputs.map { rec -> rec[4] }, // home_db
            params.n_flanking_genes,
            params.min_synteny_score,
            params.mmseqs_sensitivity
        )
        
        ITERATIVE_SEARCH.out.expanded_db.view { locus, db ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'ITERATIVE_SEARCH'.padRight(24)}${c.reset} complete for ${locus}"
        }

        // PHASE 3: Region Identification & Augmented Search
        uiPhase(3, 'Region Clustering')
        uiStatus('RUN ', 'CLUSTER_REGIONS', 'Clustering genomic regions by synteny')
        
        ITERATIVE_SEARCH.out.hits
            .join(EXTRACT_FLANKING.out.faa)
            .join(EXTRACT_FLANKING.out.bed)
            .flatMap { locus_id, hits_dir, faa_file, bed_file ->
                // Explode hits directory
                def dir_file = new File(hits_dir.toString())
                if (dir_file.exists() && dir_file.isDirectory()) {
                    (dir_file.listFiles() ?: [])
                        .findAll { hit_candidate -> hit_candidate.name.endsWith(".m8") }
                        .collect { hit_file ->
                            def genome_name = hit_file.name.replace(".m8", "").replace("${locus_id}_", "")
                            tuple(genome_name, faa_file, hit_file.toPath(), locus_id, bed_file)
                        }
                } else {
                     [] 
                }
            }
            .combine(genomes_dir_ch) // Combine with genomes_dir -> [..., genomes_dir]
            // Normalize combine output across Nextflow tuple-shape variants.
            .map { rec ->
                def norm = normalizeCombineRecord(rec, 5)
                norm ? tuple(norm[0], norm[1], norm[2], norm[3], norm[4], norm[5]) : null
            }
            .filter { rec ->
                rec != null && rec[2] != null && rec[4] != null && rec[5] != null
            }
            .set { clustering_inputs_raw }

        // Key iterative target GFFs by locus+genome so CLUSTER_REGIONS can
        // prioritize regions that actually overlap final GOI annotations.
        gff_keyed_ch = ITERATIVE_SEARCH.out.gff
            .transpose()
            .map { locus_id, gff ->
                def genome_name = gff.baseName
                tuple("${locus_id}::${genome_name}", gff)
            }

        clustering_inputs_raw
            .map { rec ->
                def key = "${rec[3]}::${rec[0]}"
                tuple(key, rec)
            }
            .join(gff_keyed_ch, remainder: true)
            // remainder:true can emit right-only tuples with rec=null;
            // keep only left-origin records for clustering inputs.
            .filter { key, rec, gff -> rec != null }
            .map { key, rec, gff ->
                def effective_target_gff = gff ?: no_gff_file
                tuple(rec[0], rec[1], rec[2], rec[3], rec[4], rec[5], effective_target_gff)
            }
            .filter { rec ->
                rec[2] != null && rec[4] != null && rec[5] != null && rec[6] != null
            }
            .set { clustering_inputs }

        CLUSTER_REGIONS(
            clustering_inputs.map { rec -> tuple(rec[0], rec[1], rec[2], rec[5], rec[6]) }, // [genome, payload, hit, genomes_dir, target_gff]
            clustering_inputs.map { rec -> tuple(rec[3], rec[4]) }, // [locus_id, synteny_bed]
            params.n_flanking_genes,
            params.min_synteny_score,
            species_map_ch
        )
        CLUSTER_REGIONS.out.bed.count().view { count ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'CLUSTER_REGIONS'.padRight(24)}${c.reset} generated ${count} clustered region set(s)"
        }

        // --- PHYLOGENY & PLOTTING ---
        // The tree FASTA is goi_for_tree.faa (= expanded_db.faa + tree-only
        // GOI hits like tandem_copy that are deliberately withheld from wave
        // seeding). Falls back to expanded_db if the new file isn't present
        // (older runs / partial reruns).

        uiPhase(4, 'Phylogenetics and Visualization')
        uiStatus('RUN ', 'COMPUTE_TREE', 'Computing GOI phylogenetic trees')

        tree_input_ch = ITERATIVE_SEARCH.out.goi_for_tree
            .ifEmpty( ITERATIVE_SEARCH.out.expanded_db )

        COMPUTE_TREE(
            tree_input_ch
        )
        
        COMPUTE_TREE.out.tree.view { locus, tree ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'COMPUTE_TREE'.padRight(24)}${c.reset} tree computed for ${locus}"
        }
        
        // PHASE 4: Miniprot-based Annotation (Replacing Augustus/HomologySearch)
        // Collect Data for Plotting (ALL loci)
        // Group GFFs and TSVs by locus_id using direct join (not collect+combine)
        gffs_by_locus_ch = ITERATIVE_SEARCH.out.gff
            .transpose()
            .map { locus_id, gff -> tuple(locus_id, gff) }
            .groupTuple()  // [locus_id, [gff1, gff2, ...]]

        tsvs_by_locus_ch = ITERATIVE_SEARCH.out.homology
            .transpose()
            .map { locus_id, tsv -> tuple(locus_id, tsv) }
            .groupTuple()  // [locus_id, [tsv1, tsv2, ...]]

        cluster_by_locus_ch = CLUSTER_REGIONS.out.bed
            .map { genome_name, payload, locus_id, bed -> tuple(locus_id, genome_name, bed) }
            .groupTuple()  // produces [locus_id, [names], [beds]]

        home_bed_by_locus_ch = EXTRACT_FLANKING.out.bed  // [locus_id, bed]
        tree_by_locus_ch = COMPUTE_TREE.out.tree         // [locus_id, tree]
        // Per-locus GOI-only BED from SPLIT_LOCI (used as plot query_bed).
        // Previously LOCATE_GENE.out.bed.first() was used for all loci, which
        // mis-anchored multi-paralog plots (e.g. TP53/TP63/TP73 → same query
        // shown on all three locus plots). Now each locus gets its own slice.
        query_bed_by_locus_ch = distinct_loci_ch        // [locus_id, bed]

        plot_inputs = home_bed_by_locus_ch
            .join(query_bed_by_locus_ch)                  // [locus_id, home_bed, query_bed]
            .join(cluster_by_locus_ch)                    // [locus_id, home_bed, query_bed, names, beds]
            .join(tree_by_locus_ch)                       // [..., tree]
            .join(gffs_by_locus_ch, remainder: true)      // [..., tree, gffs_or_null]
            .join(tsvs_by_locus_ch, remainder: true)      // [..., tree, gffs_or_null, tsvs_or_null]
            .map { entry ->
                def locus_id = entry[0]
                def home_bed = entry[1]
                def query_bed = entry[2]
                def names = entry[3]
                def beds = entry[4]
                def tree = entry[5]
                def gffs = entry[6] ?: []
                def tsvs = entry[7] ?: []
                tuple(home_bed, query_bed, names, beds, gffs, tsvs, tree)
            }

        plot_inputs.multiMap { item ->
            home_bed: item[0]
            query_bed: item[1]
            target_names: item[2]
            candidate_beds: item[3]
            target_gffs: item[4]
            homology_tsvs: item[5]
            tree: item[6]
        }.set { plot_inputs_split }

        uiStatus('RUN ', 'PLOT_SYNTENY', 'Generating synteny visualizations')

        PLOT_SYNTENY(
            plot_inputs_split.home_bed,       // home_bed
            plot_inputs_split.query_bed,      // query_bed (per-locus, was global)
            effective_home_gff_ch,            // home_gff (user-provided or Prodigal-predicted)
            plot_inputs_split.target_gffs,    // target_gffs
            plot_inputs_split.target_names,   // target_names
            plot_inputs_split.candidate_beds, // candidate_beds
            plot_inputs_split.homology_tsvs,  // homology_tsvs
            plot_inputs_split.tree,           // tree
            species_map_ch,                   // species_mapping.tsv (already a value channel)
            home_species_for_sort_ch          // home species label for the matrix-view home row
        )
        
        PLOT_SYNTENY.out.plot.view { plot ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'PLOT_SYNTENY'.padRight(24)}${c.reset} synteny visualization complete"
        }
        
        // Final Reporting
        uiPhase(5, 'Report Generation')
        uiStatus('RUN ', 'GENERATE_REPORT', 'Generating comprehensive report')
        
        // Use sentinel files so that collect() always yields a valid path list
        // that Nextflow can stage into the process work directory.
        def no_regions_sentinel = file("${projectDir}/assets/sentinels/NO_REGIONS")
        def no_gffs_sentinel    = file("${projectDir}/assets/sentinels/NO_GFFS")
        def no_homology_sentinel = file("${projectDir}/assets/sentinels/NO_HOMOLOGY")
        def no_hits_sentinel    = file("${projectDir}/assets/sentinels/NO_HITS")
        def no_augmented_sentinel = file("${projectDir}/assets/sentinels/NO_AUGMENTED")
        def no_scores_sentinel = file("${projectDir}/assets/sentinels/NO_SCORES")
        
        ITERATIVE_SEARCH.out.region_genes
            .map { rec -> rec[1] }
            .flatten()
            .ifEmpty(no_regions_sentinel)
            .collect()
            .set { collected_regions }

        ITERATIVE_SEARCH.out.gff
            .map { rec -> rec[1] }
            .flatten()
            .ifEmpty(no_gffs_sentinel)
            .collect()
            .set { collected_region_gffs }

        ITERATIVE_SEARCH.out.homology
            .map { rec -> rec[1] }
            .flatten()
            .ifEmpty(no_homology_sentinel)
            .collect()
            .set { collected_homology }
            
        ITERATIVE_SEARCH.out.hits
            .map { rec -> rec[1] }
            .ifEmpty(no_hits_sentinel)
            .collect()
            .set { collected_hits }

        CLUSTER_REGIONS.out.scores
            .map { rec -> rec[1] }
            .ifEmpty(no_scores_sentinel)
            .collect()
            .set { collected_scores }
            
        // No standalone augmented proteins - pass sentinel file
        collected_augmented = channel.value(no_augmented_sentinel)
        
        GENERATE_REPORT(
            collected_regions,
            collected_region_gffs,
            collected_homology,
            collected_hits,
            collected_augmented,
            qc_summary_ch,
            collected_scores,
            params.qc_fail_policy
        )
        
        GENERATE_REPORT.out.report.view { report ->
            "${c.green}[OK  ]${c.reset} ${c.white}${'GENERATE_REPORT'.padRight(24)}${c.reset} analysis report generated"
        }
    }

    // --- Workflow completion handler (must be inside workflow block for NF 26.x) ---
    workflow.onComplete = {
        def done_c = colors()
        // Collect task logs into results/logs/ regardless of success/failure
        try {
            collectTaskLogs(params.outdir)
        } catch (Exception e) {
            log.info "${done_c.yellow}[WARN] Could not collect task logs: ${e.message}${done_c.reset}"
        }

        log.info ""
        log.info uiRule()
        if (workflow.success) {
            uiStatus('OK  ', 'PIPELINE', 'Pipeline completed successfully')
            log.info uiRule()
            log.info "${done_c.white}Run Summary${done_c.reset}"
            log.info uiRule()
            log.info "${done_c.dim}Duration:         ${done_c.reset} ${workflow.duration}"
            log.info "${done_c.dim}Tasks Completed:  ${done_c.reset} ${workflow.stats.succeedCount}"
            if (workflow.stats.cachedCount > 0) {
                log.info "${done_c.dim}Tasks Cached:     ${done_c.reset} ${workflow.stats.cachedCount} (reused from previous run)"
            }
            log.info "${done_c.dim}Results Directory: ${done_c.reset} ${done_c.cyan}${params.outdir}${done_c.reset}"
            log.info uiRule()

            // Scan for actual output files and report what was generated
            log.info "${done_c.white}Generated Outputs${done_c.reset}"
            def outdir = new File(params.outdir.toString())
            def found_outputs = false

            // Report
            def report_file = new File(outdir, 'synvoy_report.json')
            if (report_file.exists()) {
                log.info "${done_c.green}  ✓${done_c.reset} synvoy_report.json          ${done_c.dim}(analysis summary)${done_c.reset}"
                // Try to extract key stats from report
                try {
                    def report = new groovy.json.JsonSlurper().parse(report_file)
                    def summary = report.summary
                    if (summary) {
                        def goi_count = summary.total_goi_annotations ?: 0
                        def genomes_hit = summary.genomes_with_annotations ?: 0
                        def absent = summary.goi_absent_genomes?.size() ?: 0
                        log.info "${done_c.dim}    → GOI found in ${genomes_hit} genome(s) (${goi_count} annotation(s) total)${done_c.reset}"
                        if (absent > 0) {
                            log.info "${done_c.yellow}    → GOI absent in ${absent} genome(s)${done_c.reset}"
                        }
                    }
                } catch (Exception ignored) {}
                found_outputs = true
            } else {
                log.info "${done_c.yellow}  ✗${done_c.reset} synvoy_report.json          ${done_c.yellow}(not generated)${done_c.reset}"
            }

            // Plots
            def plots = outdir.listFiles()?.findAll { output_file -> output_file.name.endsWith('_synteny_plot.html') } ?: []
            if (plots) {
                plots.each { p ->
                    log.info "${done_c.green}  ✓${done_c.reset} ${p.name.padRight(28)} ${done_c.dim}(interactive visualization)${done_c.reset}"
                }
                found_outputs = true
            }

            // Trees
            def trees = outdir.listFiles()?.findAll { output_file -> output_file.name.endsWith('_tree.nwk') } ?: []
            if (trees) {
                trees.each { t ->
                    log.info "${done_c.green}  ✓${done_c.reset} ${t.name.padRight(28)} ${done_c.dim}(GOI phylogeny)${done_c.reset}"
                }
                found_outputs = true
            }

            // Regions
            def regions_dir = new File(outdir, 'regions')
            if (regions_dir.exists()) {
                def beds = regions_dir.listFiles()?.findAll { output_file -> output_file.name.endsWith('.regions.bed') } ?: []
                if (beds) {
                    log.info "${done_c.green}  ✓${done_c.reset} regions/                      ${done_c.dim}(${beds.size()} region BED file(s))${done_c.reset}"
                    found_outputs = true
                }
            }

            if (!found_outputs) {
                log.info "${done_c.yellow}  No output files found in ${params.outdir}${done_c.reset}"
            }

            // Logs
            def logsDir = new File(params.outdir.toString(), 'logs')
            if (logsDir.exists()) {
                def logDirs = logsDir.listFiles()?.findAll { output_file -> output_file.isDirectory() } ?: []
                if (logDirs) {
                    log.info "${done_c.green}  ✓${done_c.reset} logs/                         ${done_c.dim}(task logs for ${logDirs.size()} process(es))${done_c.reset}"
                }
            }

            log.info uiRule()
        } else {
            uiStatus('FAIL', 'PIPELINE', 'Pipeline execution failed')
            log.info uiRule()
            log.info "${done_c.dim}Duration:         ${done_c.reset} ${workflow.duration}"
            log.info "${done_c.dim}Tasks Completed:  ${done_c.reset} ${workflow.stats.succeedCount}"
            log.info "${done_c.dim}Tasks Failed:     ${done_c.reset} ${done_c.red}${workflow.stats.failedCount}${done_c.reset}"
            log.info ""
            log.info "${done_c.red}Error: ${workflow.errorMessage}${done_c.reset}"
            log.info ""
            log.info "${done_c.dim}Troubleshooting tips:${done_c.reset}"
            log.info "${done_c.dim}  • Check task logs:            ${params.outdir}/logs/ (collected per process)${done_c.reset}"
            log.info "${done_c.dim}  • Check the Nextflow log:     .nextflow.log${done_c.reset}"
            log.info "${done_c.dim}  • Re-run with -resume to pick up from the last successful step${done_c.reset}"
            if (workflow.stats.failedCount > 0) {
                log.info "${done_c.dim}  • Common causes: missing tools (tblastn, mmseqs, miniprot), OOM, network timeout${done_c.reset}"
            }
            log.info uiRule()
        }
    }
}

// Collect task logs from work/ into results/logs/ for easy debugging.
// Structure: logs/<PROCESS_NAME>/stdout.log, stderr.log, script.sh
// When a process runs multiple times (e.g. per-locus), the tag from .command.run
// is used to disambiguate (e.g. logs/CLUSTER_REGIONS__GCA_001234/...).
def collectTaskLogs(outdir) {
    def logsDir = new File(outdir.toString(), 'logs')
    logsDir.mkdirs()
    def workDir = new File(workflow.workDir.toString())
    if (!workDir.exists()) return

    def seenNames = [:] as Map  // track duplicates
    workDir.eachDirRecurse { taskDir ->
        def cmdRun = new File(taskDir, '.command.run')
        def cmdLog = new File(taskDir, '.command.log')
        def cmdErr = new File(taskDir, '.command.err')
        def cmdSh  = new File(taskDir, '.command.sh')
        if (!cmdRun.exists()) return  // not a task directory

        // Skip if no log content at all
        if ((!cmdLog.exists() || cmdLog.length() == 0) && (!cmdErr.exists() || cmdErr.length() == 0)) return

        // Extract process name and tag from .command.run header
        def processName = 'unknown'
        def tag = ''
        try {
            def lines = cmdRun.readLines().take(10)
            lines.each { line ->
                def m = (line =~ /name:\s*'([^']+)'/)
                if (m.find()) {
                    processName = m.group(1).replaceAll(/[^a-zA-Z0-9_.-]/, '_')
                }
                // Also extract Nextflow outputs for tag hints
                def mOut = (line =~ /- '([^']+)'/)
                if (mOut.find() && !tag) {
                    tag = mOut.group(1).replaceAll(/.*\//, '').replaceAll(/[^a-zA-Z0-9_.-]/, '_')
                }
            }
        } catch (Exception ignored) {}

        // Build unique directory name
        def dirName = processName
        if (tag) {
            dirName = "${processName}__${tag}".take(120)
        }
        // Handle duplicates by appending counter
        def count = seenNames.getOrDefault(dirName, 0) + 1
        seenNames[dirName] = count
        if (count > 1) {
            dirName = "${dirName}__${count}"
        }

        def destDir = new File(logsDir, dirName)
        destDir.mkdirs()

        // Copy files with actual content
        [[cmdLog, 'stdout.log'], [cmdErr, 'stderr.log'], [cmdSh, 'script.sh']].each { entry ->
            def src = entry[0]
            def destName = entry[1]
            if (src.exists() && src.length() > 0) {
                try {
                    new File(destDir, destName).text = src.text
                } catch (Exception ignored) {}
            }
        }
    }
}
