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

// ==============================================================================
// ANSI Color Codes (Script-level variables)
// ==============================================================================

c_reset = "\033[0m"
c_bold = "\033[1m"
c_dim = "\033[2m"
c_black = "\033[0;30m"
c_red = "\033[0;31m"
c_green = "\033[0;32m"
c_yellow = "\033[0;33m"
c_blue = "\033[0;34m"
c_purple = "\033[0;35m"
c_cyan = "\033[0;36m"
c_white = "\033[0;37m"

// ==============================================================================
// Console UI helpers
// ==============================================================================

def uiRule() {
    return "${c_blue}${'═' * 63}${c_reset}"
}

def uiStatus(String level, String task, String detail = '') {
    def levelColors = [
        'RUN ': c_blue,
        'OK  ': c_green,
        'INFO': c_cyan,
        'WARN': c_yellow,
        'SKIP': c_dim,
        'FAIL': c_red
    ]
    def key = (level ?: 'INFO').padRight(4).substring(0, 4)
    def levelColor = levelColors.get(key, c_white)
    def taskCol = task ? "${c_white}${task.padRight(24)}${c_reset}" : ''
    def detailCol = detail ?: ''
    def prefix = "${levelColor}[${key}]${c_reset}"
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
    log.info ""
    log.info uiRule()
    log.info "${c_white}Phase ${idx}${c_reset} ${c_dim}|${c_reset} ${c_cyan}${title}${c_reset}"
    log.info uiRule()
}

def printHeader() {
    log.info ""
    log.info uiRule()
    log.info "${c_cyan}${c_bold}SynTerra${c_reset} ${c_dim}v2.0${c_reset}"
    log.info "${c_dim}Phylogenetically-informed syntenic ortholog discovery${c_reset}"
    log.info uiRule()
}

def printParams() {
    def query_display = params.mode == 'easy' ? params.query_id : (params.query ? new File(params.query).name : 'N/A')
    def home_display = params.mode == 'easy' ? params.home_species : (params.home_genome ? new File(params.home_genome).name : 'N/A')
    def target_display = params.target_species ?: 'auto (taxonomic search)'
    
    log.info uiRule()
    log.info "${c_white}Run Configuration${c_reset}"
    log.info uiRule()
    log.info "${c_dim}Query Gene     ${c_reset} ${c_green}${query_display}${c_reset}"
    log.info "${c_dim}Home Genome    ${c_reset} ${c_green}${home_display}${c_reset}"
    log.info "${c_dim}Mode           ${c_reset} ${c_yellow}${params.mode}${c_reset}"
    log.info "${c_dim}Target Species ${c_reset} ${c_cyan}${target_display}${c_reset}"
    log.info "${c_dim}Flanking Genes ${c_reset} ${params.n_flanking_genes}"
    log.info "${c_dim}MMseqs Sens.   ${c_reset} ${params.mmseqs_sensitivity}"
    if (params.mode == 'easy') {
        log.info "${c_dim}Asm Ranking    ${c_reset} ${params.assembly_ranking}"
        log.info "${c_dim}LowQ Policy    ${c_reset} ${params.bad_quality_policy} (timeout=${params.bad_quality_timeout}s)"
    }
    log.info "${c_dim}Output Dir     ${c_reset} ${params.outdir}"
    log.info uiRule()
}

printHeader()
printParams()

// Stable sentinel files for optional path inputs.
def no_gff_file = file("${projectDir}/assets/sentinels/NO_GFF")
def no_species_map_file = file("${projectDir}/assets/sentinels/NO_SPECIES_MAP")

workflow {
    log.info ""
    
    // ========== INPUT VALIDATION ==========
    
    // Mode-specific validation
    if (params.mode == 'easy') {
        if (!params.query_id) {
            log.error "${c_red}Easy mode requires --query_id (UniProt or NCBI ID)${c_reset}"
            exit 1
        }
        if (params.query) {
            uiStatus('WARN', 'INPUT', 'Pro mode flag --query was provided but ignored in easy mode; using --query_id.')
        }
    } else if (params.mode == 'pro') {
        if (!params.query) {
            log.error "${c_red}Pro mode requires --query (path to FASTA file)${c_reset}"
            exit 1
        }
        if (!params.home_genome) { 
            log.error "${c_red}Pro mode requires --home_genome (path to home FASTA)${c_reset}"
            exit 1
        }
        if (!file(params.home_genome).exists()) {
            log.error "${c_red}Home genome not found: ${params.home_genome}${c_reset}"
            exit 1
        }
        if (!file(params.query).exists()) {
            log.error "${c_red}Query FASTA not found: ${params.query}${c_reset}"
            exit 1
        }
    } else {
        log.error "${c_red}Invalid mode: ${params.mode}. Expected 'easy' or 'pro'${c_reset}"
        exit 1
    }
    
    uiStatus('OK', 'INPUT', 'Validation passed')
    
    // Channel setup — depends on mode
    if (params.mode == 'easy') {
        // ID/symbol mode: resolve input via resolver process.
        def gene_input = params.query_id
        def species_override = params.home_species ?: ''
        
        uiStatus('RUN ', 'RESOLVE_QUERY', 'Easy mode: resolving query_id input')
        RESOLVE_GENE_INPUT(gene_input, species_override)
        
        // Use resolved FASTA as query
        raw_gene_ch = RESOLVE_GENE_INPUT.out.fasta
        
        // Get resolved species (auto-detected from ID, or user-provided)
        resolved_species_ch = RESOLVE_GENE_INPUT.out.species.map { it.text.trim() }
        
        // Determine home species: user-provided takes priority, else auto-detected
        home_species_ch = resolved_species_ch.map { resolved ->
            def species = params.home_species ?: resolved
            if (!species) {
                log.error "${c_red}Could not detect species. Please provide --home_species${c_reset}"
                exit 1
            }
            return species
        }
        
        // Fetch home genome automatically for easy mode
        FETCH_HOME_GENOME(home_species_ch)
        home_genome_ch = FETCH_HOME_GENOME.out.genome
        // Use GFF if available, otherwise mark as missing
        home_gff_ch = FETCH_HOME_GENOME.out.gff.ifEmpty(no_gff_file)
        
        // Fetch related genomes for easy mode
        def max_genomes = (params.max_genomes == null ? 10 : params.max_genomes as Integer)
        def target_species = params.target_species ?: ''
        FETCH_RELATED_GENOMES(home_species_ch, max_genomes, target_species)
        genomes_dir_ch = FETCH_RELATED_GENOMES.out.genomes_dir
        species_map_ch = FETCH_RELATED_GENOMES.out.species_map.first()
        // Species name for phylogenetic sorting
        home_species_for_sort_ch = home_species_ch
        
    } else {
        // --- Pro mode: User provides files directly ---
        
        // 1. Query Setup
        raw_gene_ch = Channel.fromPath(params.query)
        
        // 2. Home Genome Setup
        home_genome_ch = Channel.fromPath(params.home_genome)
        
        if (params.home_gff) {
            home_gff_ch = Channel.value(file(params.home_gff, checkIfExists: true))
        } else {
            home_gff_ch = Channel.value(no_gff_file)
        }
        
        // 3. Target Genomes Setup
        if (params.target_genomes) {
            uiStatus('RUN ', 'STAGE_GENOMES', 'Loading target genomes list')
            // Support both glob patterns ("genomes/*.fna") and comma-separated
            // lists ("a.fna,b.fna,c.fna") as well as Nextflow list syntax.
            def tg = params.target_genomes
            if (tg instanceof List) {
                target_genomes_list = Channel.fromPath(tg).collect()
            } else if (tg.toString().contains(',')) {
                target_genomes_list = Channel
                    .fromPath(tg.toString().split(',').collect { it.trim() })
                    .collect()
            } else {
                target_genomes_list = Channel.fromPath(tg).collect()
            }
            
            // Show count
            target_genomes_list.view { genomes ->
                "${c_green}[OK  ]${c_reset} ${c_white}${'STAGE_GENOMES'.padRight(24)}${c_reset} staged ${genomes.size()} target genomes"
            }
            
            STAGE_GENOMES(target_genomes_list)
            genomes_dir_ch = STAGE_GENOMES.out.dir
            species_map_ch = STAGE_GENOMES.out.species_map.first()
            
        } else {
            uiStatus('WARN', 'STAGE_GENOMES', 'No target genomes provided; running home-genome-only analysis')
            genomes_dir_ch = Channel.empty()
            species_map_ch = Channel.value(no_species_map_file)
        }
        // Species name for phylogenetic sorting (pro mode: use param or extract from filename)
        home_species_for_sort_ch = Channel.value(params.home_species ?: params.home_genome)
    }

    // Normalize query to protein space (DNA queries are translated to best ORF)
    // to keep downstream search/annotation behavior consistent.
    NORMALIZE_QUERY(raw_gene_ch)
    normalized_gene_ch = NORMALIZE_QUERY.out.fasta

    // PHASE 1: Core Localization
    uiPhase(1, 'Gene Localization in Home Genome')
    
    uiStatus('RUN ', 'LOCATE_GENE', 'Locating GOI in home genome')
    LOCATE_GENE(normalized_gene_ch, home_genome_ch)
    
    // 4b. ANNOTATE GOI EXONS
    // Uses hits from LOCATE_GENE to annotate individual exons of the GOI
    // If GFF available: matches GOI to annotated gene and extracts CDS/exons
    // If no GFF: uses tblastn hits to detect exon boundaries (splice sites, start/stop codons)
    uiStatus('RUN ', 'ANNOTATE_GOI', 'Annotating GOI exons')
    
    // Determine query_id for name-based GFF matching
    def effective_query_id = params.query_id ?: ''
    
    ANNOTATE_GOI(
        normalized_gene_ch.first(),
        home_genome_ch,
        home_gff_ch,
        LOCATE_GENE.out.blast_hits,
        LOCATE_GENE.out.mmseqs_hits,
        effective_query_id
    )
    
    ANNOTATE_GOI.out.info.view { info ->
        "${c_green}[OK  ]${c_reset} ${c_white}${'ANNOTATE_GOI'.padRight(24)}${c_reset} GOI exon annotation complete"
    }
    
    // 5. SPLIT LOCI
    SPLIT_LOCI(LOCATE_GENE.out.bed)
    
    SPLIT_LOCI.out.beds.flatten().count().view { count ->
        "${c_green}[OK  ]${c_reset} ${c_white}${'SPLIT_LOCI'.padRight(24)}${c_reset} identified ${count} locus/loci"
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
            "${c_green}[OK  ]${c_reset} ${c_white}${'PREPARE_HOME'.padRight(24)}${c_reset} home proteome ready"
        }

        // Borrow only when home genome has no usable GFF.
        gff_status.real.view { gff ->
            "${c_dim}[SKIP]${c_reset} ${c_white}${'BORROW_ANNOT'.padRight(24)}${c_reset} home GFF found (${gff.name})"
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
            "${c_green}[OK  ]${c_reset} ${c_white}${'BORROW_ANNOT'.padRight(24)}${c_reset} borrowed annotations generated"
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
                gff_status.missing.combine(fallback_gff_ch).map { it[1] }
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
        params.prefer_large_genes
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
            .set { phylo_sort_inputs } // [locus_id, home_genome]

        uiStatus('RUN ', 'PHYLO_SORT', 'Sorting genomes by phylogenetic distance')
        
        PHYLO_SORT(
            phylo_sort_inputs,
            genomes_dir_ch,
            home_species_for_sort_ch
        )
        
        PHYLO_SORT.out.sorted_list.view { locus, sorted ->
            "${c_green}[OK  ]${c_reset} ${c_white}${'PHYLO_SORT'.padRight(24)}${c_reset} ordering complete for ${locus}"
        }
        
        // QC
        uiStatus('RUN ', 'GENOME_QC', 'Assessing target genome quality')
        
        ASSESS_GENOME_QUALITY(genomes_dir_ch)
        qc_summary_ch = ASSESS_GENOME_QUALITY.out.json

        FILTER_SORTED_GENOMES(
            PHYLO_SORT.out.sorted_list,
            qc_summary_ch,
            params.qc_fail_policy
        )

        FILTER_SORTED_GENOMES.out.sorted_list.view { locus, sorted ->
            "${c_green}[OK  ]${c_reset} ${c_white}${'QC_FILTER'.padRight(24)}${c_reset} filtered target list for ${locus}"
        }

        // 8. Iterative Search (FOR EACH LOCUS) - Using FIXED database with GOI
        PREPARE_INITIAL_DB.out.db
            .join(FILTER_SORTED_GENOMES.out.sorted_list)
            .set { iterative_search_inputs_partial } // [locus_id, initial_db, sorted_list]
            
        iterative_search_inputs_partial
            .combine(genomes_dir_ch)
            .set { iterative_search_inputs } // [locus_id, faa, sorted_list, genomes_dir]

        iterative_search_inputs
            .combine(home_proteome_db_ch)
            .set { iterative_search_final_inputs } // [locus_id, faa, sorted_list, genomes_dir, home_db]

        uiStatus('RUN ', 'ITERATIVE_SEARCH', 'Running iterative phylogenetic search')
        
        ITERATIVE_SEARCH(
            iterative_search_final_inputs.map { tuple(it[0], it[1]) }, // [locus, faa]
            iterative_search_final_inputs.map { it[2] }, // sorted_list
            iterative_search_final_inputs.map { it[3] }, // genomes_dir
            iterative_search_final_inputs.map { it[4] }, // home_db
            params.n_flanking_genes,
            params.min_synteny_score,
            params.mmseqs_sensitivity
        )
        
        ITERATIVE_SEARCH.out.expanded_db.view { locus, db ->
            "${c_green}[OK  ]${c_reset} ${c_white}${'ITERATIVE_SEARCH'.padRight(24)}${c_reset} complete for ${locus}"
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
                        .findAll { it.name.endsWith(".m8") }
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
                if (!(rec instanceof List)) {
                    return null
                }
                // Variant A: [left_tuple, genomes_dir]
                if (rec.size() == 2 && rec[0] instanceof List) {
                    def left = rec[0]
                    if (left.size() < 5) return null
                    return tuple(left[0], left[1], left[2], left[3], left[4], rec[1])
                }
                // Variant B: already flattened [genome, faa, hits, locus, bed, genomes_dir]
                if (rec.size() >= 6) {
                    return tuple(rec[0], rec[1], rec[2], rec[3], rec[4], rec[5])
                }
                return null
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
            clustering_inputs.map { tuple(it[0], it[1], it[2], it[5], it[6]) }, // [genome, payload, hit, genomes_dir, target_gff]
            clustering_inputs.map { tuple(it[3], it[4]) }, // [locus_id, synteny_bed]
            params.n_flanking_genes,
            params.min_synteny_score
        )
        CLUSTER_REGIONS.out.bed.count().view { count ->
            "${c_green}[OK  ]${c_reset} ${c_white}${'CLUSTER_REGIONS'.padRight(24)}${c_reset} generated ${count} clustered region set(s)"
        }

        // --- PHYLOGENY & PLOTTING ---
        // Collect all proteins for tree: Flanking Genes + Discovered Genes (from Expanded DB or Regions?)
        // Iterative Search output might be best source of all gene sequences found.
        // expanded_db contains everything found so far.
        // Let's use expanded_db per locus.
        
        uiPhase(4, 'Phylogenetics and Visualization')
        uiStatus('RUN ', 'COMPUTE_TREE', 'Computing GOI phylogenetic trees')
        
        COMPUTE_TREE(
            ITERATIVE_SEARCH.out.expanded_db
        )
        
        COMPUTE_TREE.out.tree.view { locus, tree ->
            "${c_green}[OK  ]${c_reset} ${c_white}${'COMPUTE_TREE'.padRight(24)}${c_reset} tree computed for ${locus}"
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

        plot_inputs = home_bed_by_locus_ch
            .join(cluster_by_locus_ch)                    // [locus_id, home_bed, names, beds]
            .join(tree_by_locus_ch)                       // [locus_id, home_bed, names, beds, tree]
            .join(gffs_by_locus_ch, remainder: true)      // [locus_id, home_bed, names, beds, tree, gffs_or_null]
            .join(tsvs_by_locus_ch, remainder: true)      // [locus_id, home_bed, names, beds, tree, gffs_or_null, tsvs_or_null]
            .map { entry ->
                def locus_id = entry[0]
                def home_bed = entry[1]
                def names = entry[2]
                def beds = entry[3]
                def tree = entry[4]
                def gffs = entry[5] ?: []
                def tsvs = entry[6] ?: []
                tuple(home_bed, names, beds, gffs, tsvs, tree)
            }

        plot_inputs.multiMap { item ->
            home_bed: item[0]
            target_names: item[1]
            candidate_beds: item[2]
            target_gffs: item[3]
            homology_tsvs: item[4]
            tree: item[5]
        }.set { plot_inputs_split }

        uiStatus('RUN ', 'PLOT_SYNTENY', 'Generating synteny visualizations')
            
        PLOT_SYNTENY(
            plot_inputs_split.home_bed,   // home_bed
            LOCATE_GENE.out.bed.first(),  // query_bed
            effective_home_gff_ch,        // home_gff (user-provided or Prodigal-predicted)
            plot_inputs_split.target_gffs,    // target_gffs
            plot_inputs_split.target_names,   // target_names
            plot_inputs_split.candidate_beds, // candidate_beds
            plot_inputs_split.homology_tsvs,  // homology_tsvs
            plot_inputs_split.tree,           // tree
            species_map_ch                    // species_mapping.tsv (already a value channel)
        )
        
        PLOT_SYNTENY.out.plot.view { plot ->
            "${c_green}[OK  ]${c_reset} ${c_white}${'PLOT_SYNTENY'.padRight(24)}${c_reset} synteny visualization complete"
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
            .map { it[1] } 
            .flatten()
            .collect()
            .ifEmpty(no_regions_sentinel)
            .set { collected_regions }

        ITERATIVE_SEARCH.out.gff
            .map { it[1] }
            .flatten()
            .collect()
            .ifEmpty(no_gffs_sentinel)
            .set { collected_region_gffs }

        ITERATIVE_SEARCH.out.homology
            .map { it[1] }
            .flatten()
            .collect()
            .ifEmpty(no_homology_sentinel)
            .set { collected_homology }
            
        ITERATIVE_SEARCH.out.hits
            .map { it[1] } 
            .collect()
            .ifEmpty(no_hits_sentinel)
            .set { collected_hits }

        CLUSTER_REGIONS.out.scores
            .map { it[1] }
            .collect()
            .ifEmpty(no_scores_sentinel)
            .set { collected_scores }
            
        // No standalone augmented proteins - pass sentinel file
        collected_augmented = Channel.value(no_augmented_sentinel)
        
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
            "${c_green}[OK  ]${c_reset} ${c_white}${'GENERATE_REPORT'.padRight(24)}${c_reset} analysis report generated"
        }
    }
}

workflow.onComplete {
    log.info ""
    log.info uiRule()
    if (workflow.success) {
        uiStatus('OK  ', 'PIPELINE', 'Pipeline completed successfully')
        log.info "${c_dim}Results Directory:${c_reset} ${c_cyan}${params.outdir}${c_reset}"
        log.info "${c_dim}Duration:         ${c_reset} ${workflow.duration}"
        log.info "${c_dim}Tasks Completed:  ${c_reset} ${workflow.stats.succeedCount}"
        log.info uiRule()
        log.info "${c_dim}Key outputs:${c_reset}"

        // Check report existence before listing it
        def report_file = file("${params.outdir}/synterra_report.json")
        if (report_file.exists()) {
            log.info "${c_dim}- synterra_report.json          (analysis summary)${c_reset}"
        } else {
            log.info "${c_yellow}- synterra_report.json          (NOT GENERATED)${c_reset}"
        }

        log.info "${c_dim}- *_synteny_plot.html           (interactive visualization)${c_reset}"
        log.info "${c_dim}- *_tree.nwk                    (GOI phylogeny)${c_reset}"
        log.info "${c_dim}- regions/*.regions.bed         (candidate regions)${c_reset}"
        log.info "${c_dim}- intermediate/                 (per-phase artifacts)${c_reset}"
        log.info uiRule()
    } else {
        uiStatus('FAIL', 'PIPELINE', 'Pipeline execution failed')
        log.info "${c_dim}Duration:${c_reset} ${workflow.duration}"
        log.info "${c_dim}Error:   ${c_reset} ${c_red}${workflow.errorMessage}${c_reset}"
        log.info uiRule()
    }
}
