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
include { PHYLO_SORT } from './modules/phylo_sort.nf'
include { FETCH_RELATED_GENOMES } from './modules/fetch_related.nf'
include { FETCH_HOME_GENOME } from './modules/fetch_home.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'
include { BORROW_ANNOTATIONS } from './modules/borrow_annotations.nf'
include { NORMALIZE_QUERY } from './modules/normalize_query.nf'

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
// ASCII Banner & Pipeline Info
// ==============================================================================

def printHeader() {
    log.info """
    ${c_blue}═══════════════════════════════════════════════════════════════${c_reset}
    ${c_cyan}S Y N T E R R A${c_reset}   ${c_dim}v2.0${c_reset}
    ${c_dim}Phylogenetically-informed syntenic ortholog discovery${c_reset}
    ${c_blue}═══════════════════════════════════════════════════════════════${c_reset}
    """.stripIndent()
}

def printParams() {
    def query_display = params.gene ? new File(params.gene).name : params.query_id
    def home_display = params.mode == 'easy' ? params.home_species : (params.home_genome ? new File(params.home_genome).name : 'N/A')
    def target_display = params.target_species ?: 'auto (taxonomic search)'
    
    log.info """
    ${c_blue}═══════════════════════════════════════════════════════════════
    ${c_white}RUN CONFIGURATION${c_reset}
    ${c_blue}═══════════════════════════════════════════════════════════════${c_reset}
    ${c_dim}Query Gene      :${c_reset} ${c_green}${query_display}${c_reset}
    ${c_dim}Home Genome     :${c_reset} ${c_green}${home_display}${c_reset}
    ${c_dim}Mode            :${c_reset} ${c_yellow}${params.mode}${c_reset}
    ${c_dim}Target Species  :${c_reset} ${c_cyan}${target_display}${c_reset}
    ${c_dim}Flanking Genes  :${c_reset} ${params.n_flanking_genes}
    ${c_dim}MMseqs Sens.    :${c_reset} ${params.mmseqs_sensitivity}
    ${c_dim}Output Dir      :${c_reset} ${params.outdir}
    ${c_blue}═══════════════════════════════════════════════════════════════${c_reset}
    """.stripIndent()
}

printHeader()
printParams()

workflow {
    log.info ""
    
    // ========== INPUT VALIDATION ==========
    
    // Check if we have EITHER gene file OR query_id
    if (!params.gene && !params.query_id) { 
        log.error """
        ${c_red}╔════════════════════════════════════════╗
        ║  ERROR: No query provided!         ║
        ╚════════════════════════════════════════╝${c_reset}
        Please provide either --gene or --query_id
        """
        exit 1
    }
    
    if (params.gene && params.query_id) {
        log.warn "${c_yellow}WARNING: Both --gene and --query_id provided. Using --gene${c_reset}"
    }
    
    // Easy mode: require home_species, Pro mode: require home_genome
    if (params.mode == 'easy') {
        if (!params.home_species) {
            log.error "${c_red}ERROR: Easy mode requires --home_species parameter!${c_reset}"
            log.error "${c_red}Example: --home_species 'Apis mellifera'${c_reset}"
            exit 1
        }
    } else {
        if (!params.home_genome) { 
            log.error "${c_red}ERROR: Pro mode requires --home_genome parameter!${c_reset}"
            exit 1
        }
        
        if (!file(params.home_genome).exists()) {
            log.error "${c_red}ERROR: Home genome not found: ${params.home_genome}${c_reset}"
            exit 1
        }
    }
    
    if (params.gene && !file(params.gene).exists()) {
        log.error "${c_red}ERROR: Gene file not found: ${params.gene}${c_reset}"
        exit 1
    }
    
    log.info "${c_green}Input validation passed${c_reset}"
    
    // Channel setup
    if (params.gene) {
        raw_gene_ch = Channel.fromPath(params.gene)
    } else {
        FETCH_QUERY_FROM_ID(params.query_id)
        raw_gene_ch = FETCH_QUERY_FROM_ID.out.fasta
    }
    
    // Normalize query: translate nucleotide queries to protein
    NORMALIZE_QUERY(raw_gene_ch)
    normalized_gene_ch = NORMALIZE_QUERY.out.fasta
    
    normalized_gene_ch.multiMap { it ->
        loc: it
        aug: it
    }.set { gene_inputs }
    
    query_gene_source_ch = gene_inputs.loc
    aug_query_gene_ch = gene_inputs.aug.first()

    // Handle Easy vs Pro mode for home genome and targets
    if (params.mode == 'easy') {
        if (!params.home_species) {
            log.error "${c_red}ERROR: Easy mode requires --home_species parameter!${c_reset}"
            log.error "${c_red}Example: --home_species 'Apis mellifera'${c_reset}"
            exit 1
        }
        
        log.info "${c_cyan}Easy mode: Fetching genomes for ${c_white}${params.home_species}${c_reset}"
        
        // Fetch home genome from NCBI
        log.info "${c_cyan}  - Downloading reference genome...${c_reset}"
        FETCH_HOME_GENOME(params.home_species)
        home_genome_ch = FETCH_HOME_GENOME.out.genome
        home_gff_ch = FETCH_HOME_GENOME.out.gff.ifEmpty(file("NO_GFF"))
        
        // Fetch related genomes
        log.info "${c_cyan}  - Downloading related genomes...${c_reset}"
        def target_species_val = params.target_species ?: ''
        FETCH_RELATED_GENOMES(params.home_species, params.max_genomes, target_species_val)
        genomes_dir_ch = FETCH_RELATED_GENOMES.out.genomes_dir
        species_map_ch = FETCH_RELATED_GENOMES.out.species_map
        
        // Count genomes found
        genomes_dir_ch.view { dir ->
            def genome_count = new File(dir.toString()).listFiles().findAll { 
                it.name.endsWith('.fna') || it.name.endsWith('.fasta') || it.name.endsWith('.fa')
            }.size()
            "${c_green}Downloaded ${genome_count} related genome(s)${c_reset}"
        }
        
    } else {
        // Pro mode - user provides files
        if (!params.home_genome) {
            log.error "${c_red}ERROR: Pro mode requires --home_genome parameter!${c_reset}"
            exit 1
        }
        
        home_genome_ch = Channel.fromPath(params.home_genome)
        
        if (params.home_gff) {
            home_gff_ch = Channel.fromPath(params.home_gff).first()
        } else {
            home_gff_ch = Channel.value(file("NO_GFF"))
        }
        
        if (params.target_genomes) {
            log.info "${c_cyan}Loading target genomes list...${c_reset}"
            target_genomes_list = Channel.fromPath(params.target_genomes).collect()
            
            // Show count
            target_genomes_list.view { genomes ->
                "${c_green}[STAGE] Staged ${genomes.size()} target genomes${c_reset}"
            }
            
            STAGE_GENOMES(target_genomes_list)
            genomes_dir_ch = STAGE_GENOMES.out.dir
            species_map_ch = Channel.fromPath("NO_SPECIES_MAP")
            
        } else {
            log.warn "${c_yellow}WARNING: No target genomes provided - running home genome analysis only${c_reset}"
            genomes_dir_ch = Channel.empty()
            species_map_ch = Channel.fromPath("NO_SPECIES_MAP")
        }
    }

    // PHASE 1: Core Localization
    log.info "\n${c_blue}═══════════════════════════════════════════════════════════════"
    log.info "${c_white}PHASE 1: Gene Localization in Home Genome${c_reset}"
    log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
    
    LOCATE_GENE(query_gene_source_ch, home_genome_ch)
    
    // 4b. ANNOTATE GOI EXONS
    // Uses hits from LOCATE_GENE to annotate individual exons of the GOI
    // If GFF available: matches GOI to annotated gene and extracts CDS/exons
    // If no GFF: uses tblastn hits to detect exon boundaries (splice sites, start/stop codons)
    log.info "${c_cyan}Annotating GOI exons...${c_reset}"
    
    // Determine query_id for name-based GFF matching
    def effective_query_id = params.query_id ?: ''
    
    ANNOTATE_GOI(
        query_gene_source_ch.first(),
        home_genome_ch,
        home_gff_ch,
        LOCATE_GENE.out.blast_hits,
        LOCATE_GENE.out.mmseqs_hits,
        effective_query_id
    )
    
    ANNOTATE_GOI.out.info.view { info ->
        "${c_green}GOI exon annotation complete${c_reset}"
    }
    
    // 5. SPLIT LOCI
    SPLIT_LOCI(LOCATE_GENE.out.bed)
    
    SPLIT_LOCI.out.beds.flatten().count().view { count ->
        "${c_green}Identified ${count} distinct locus/loci${c_reset}"
    }
    
    distinct_loci_ch = SPLIT_LOCI.out.beds.flatten()
        .map { file -> tuple(file.baseName, file) }
    
    // Prepare effective home GFF (borrowed annotations + predictions) when targets exist
    def has_targets = params.target_genomes || params.mode == 'easy'
    if (has_targets) {
        // Prepare Home Proteome (Run once, used for RBH + borrowing)
        log.info "${c_cyan}[PREPARE] Preparing home proteome database...${c_reset}"
        PREPARE_HOME_PROTEOME(home_genome_ch, home_gff_ch, LOCATE_GENE.out.bed.first())
        home_proteome_db_ch = PREPARE_HOME_PROTEOME.out.db
        
        // Always borrow annotations - valuable when home genome lacks GFF
        log.info "${c_cyan}[BORROW] Checking for annotated target genomes to borrow gene models...${c_reset}"
        BORROW_ANNOTATIONS(
            home_genome_ch,
            PREPARE_HOME_PROTEOME.out.faa,
            genomes_dir_ch,
            LOCATE_GENE.out.bed.first(),
            params.n_flanking_genes
        )
        BORROW_ANNOTATIONS.out.gff.view { gff ->
            "${c_green}[BORROW] Borrowed annotations generated${c_reset}"
        }

        // Build effective GFF from all available sources:
        // 1. User-provided / NCBI GFF (if present and real)
        // 2. Prodigal-predicted GFF (when no annotation was available)
        // 3. Borrowed annotations from annotated target genomes
        home_gff_ch.branch { gff ->
            real: gff.name != 'NO_GFF'
            missing: true
        }.set { gff_status }
        
        fallback_gff_ch = PREPARE_HOME_PROTEOME.out.gff
            .mix(BORROW_ANNOTATIONS.out.gff)
            .collectFile(name: 'merged_home_annotations.gff')
            .ifEmpty(file("NO_GFF"))
        
        effective_home_gff_ch = gff_status.real
            .concat(
                gff_status.missing.combine(fallback_gff_ch).map { it[1] }
            )
            .first()
    } else {
        effective_home_gff_ch = home_gff_ch
    }

    // 6. Extract Flanking Genes
    log.info "${c_cyan}Extracting flanking genes (n=${params.n_flanking_genes})...${c_reset}"
    
    EXTRACT_FLANKING(
        distinct_loci_ch, 
        effective_home_gff_ch, 
        home_genome_ch,
        params.n_flanking_genes,
        params.min_flanking_size,
        params.prefer_large_genes
    )
    
    // 6b. CRITICAL FIX: Prepare Initial Database with GOI included
    // Combine flanking genes with query gene for iterative search
    log.info "${c_cyan}Preparing initial database with query gene...${c_reset}"
    
    PREPARE_INITIAL_DB(
        EXTRACT_FLANKING.out.faa,
        ANNOTATE_GOI.out.exons.first()  // .first() → value channel so it pairs with ALL loci
    )
        // Only run if we have targets
    if (params.target_genomes || params.mode == 'easy') {
        
        log.info "\n${c_blue}═══════════════════════════════════════════════════════════════"
        log.info "${c_white}PHASE 2: Phylogenetic Ordering & Iterative Search${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        
        EXTRACT_FLANKING.out.bed
            .map { locus_id, bed -> locus_id }
            .combine(home_genome_ch)
            .set { phylo_sort_inputs } // [locus_id, home_genome]

        log.info "${c_cyan}[PHYLO] Sorting genomes by phylogenetic distance...${c_reset}"
        
        PHYLO_SORT(
            phylo_sort_inputs,
            genomes_dir_ch
        )
        
        PHYLO_SORT.out.sorted_list.view { locus, sorted ->
            "${c_green}[PHYLO] Phylogenetic ordering complete for ${locus}${c_reset}"
        }
        
        // 8. Iterative Search (FOR EACH LOCUS) - Using FIXED database with GOI
        PREPARE_INITIAL_DB.out.db
            .join(PHYLO_SORT.out.sorted_list) 
            .set { iterative_search_inputs_partial } // [locus_id, initial_db, sorted_list]
            
        iterative_search_inputs_partial
            .combine(genomes_dir_ch)
            .set { iterative_search_inputs } // [locus_id, faa, sorted_list, genomes_dir]

        // QC
        log.info "${c_cyan}[QC] Assessing genome quality...${c_reset}"
        
        ASSESS_GENOME_QUALITY(genomes_dir_ch)
        qc_summary_ch = ASSESS_GENOME_QUALITY.out.json

        iterative_search_inputs
            .combine(home_proteome_db_ch)
            .set { iterative_search_final_inputs } // [locus_id, faa, sorted_list, genomes_dir, home_db]

        log.info "${c_cyan}[SEARCH] Running iterative phylogenetic search...${c_reset}"
        
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
            "${c_green}[SEARCH] Iterative search complete: ${locus}${c_reset}"
        }

        // PHASE 3: Region Identification & Augmented Search
        log.info "\n${c_blue}═══════════════════════════════════════════════════════════════"
        log.info "${c_white}PHASE 3: Region Clustering & Augmented Search${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        
        log.info "${c_cyan}[CLUSTER] Clustering genomic regions by synteny...${c_reset}"
        
        ITERATIVE_SEARCH.out.hits
            .join(EXTRACT_FLANKING.out.faa)
            .join(EXTRACT_FLANKING.out.bed)
            .map { locus_id, hits_dir, faa_file, bed_file ->
                // Explode hits directory
                def dir_file = new File(hits_dir.toString())
                if (dir_file.exists() && dir_file.isDirectory()) {
                    dir_file.listFiles()
                        .findAll { it.name.endsWith(".m8") }
                        .collect { hit_file ->
                            def genome_name = hit_file.name.replace(".m8", "").replace("${locus_id}_", "")
                            tuple(genome_name, faa_file, hit_file.toPath(), locus_id, bed_file)
                        }
                } else {
                     [] 
                }
            }
            .flatten()
            .collate(5) // [genome_name, faa_file, hit_file, locus_id, bed_file]
            .combine(genomes_dir_ch) // Combine with genomes_dir -> [..., genomes_dir]
            .set { clustering_inputs }

        CLUSTER_REGIONS(
            clustering_inputs.map { tuple(it[0], it[1], it[2], it[5]) }, // [genome, payload, hit, genomes_dir]
            clustering_inputs.map { tuple(it[3], it[4]) }, // [locus_id, synteny_bed]
            params.n_flanking_genes,
            params.min_synteny_score
        )
        
        def clustered_regions_ch = CLUSTER_REGIONS.out.bed
            .map { genome_name, payload_faa, locus_id, region_bed ->
                def unique_id = "${genome_name}_${locus_id}"
                tuple(unique_id, locus_id, region_bed, payload_faa, genome_name)
            }
            .combine(genomes_dir_ch) // [unique, locus, bed, faa, gname, genomes_dir]
            .set { joined_ch }
            
        CLUSTER_REGIONS.out.bed.view { genome, payload, locus, bed ->
            "${c_green}[CLUSTER] Clustered regions for ${locus} in ${genome}${c_reset}"
        }

        // --- PHYLOGENY & PLOTTING ---
        // Collect all proteins for tree: Flanking Genes + Discovered Genes (from Expanded DB or Regions?)
        // Iterative Search output might be best source of all gene sequences found.
        // expanded_db contains everything found so far.
        // Let's use expanded_db per locus.
        
        log.info "\n${c_blue}═══════════════════════════════════════════════════════════════"
        log.info "${c_white}PHASE 4: Phylogenetics & Visualization${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        
        log.info "${c_cyan}[TREE] Computing phylogenetic trees...${c_reset}"
        
        COMPUTE_TREE(
            ITERATIVE_SEARCH.out.expanded_db
        )
        
        COMPUTE_TREE.out.tree.view { locus, tree ->
            "${c_green}[TREE] Phylogenetic tree computed: ${locus}${c_reset}"
        }
        
        // PHASE 4: Miniprot-based Annotation (Replacing Augustus/HomologySearch)
        // Collect Data for Plotting (ALL loci)
        miniprot_gffs_ch = ITERATIVE_SEARCH.out.gff
            .transpose()
            .map { locus_id, gff -> tuple(locus_id, gff) }

        miniprot_tsvs_ch = ITERATIVE_SEARCH.out.homology
            .transpose()
            .map { locus_id, tsv -> tuple(locus_id, tsv) }

        gffs_by_locus_map_ch = miniprot_gffs_ch
            .groupTuple()
            .collect()
            .map { items ->
                def m = [:]
                items.each { entry ->
                    def locus_id = entry[0]
                    def gffs = entry[1]
                    m[locus_id] = gffs
                }
                m
            }

        tsvs_by_locus_map_ch = miniprot_tsvs_ch
            .groupTuple()
            .collect()
            .map { items ->
                def m = [:]
                items.each { entry ->
                    def locus_id = entry[0]
                    def tsvs = entry[1]
                    m[locus_id] = tsvs
                }
                m
            }

        cluster_by_locus_ch = CLUSTER_REGIONS.out.bed
            .map { genome_name, payload, locus_id, bed -> tuple(locus_id, genome_name, bed) }
            .groupTuple()
            .map { locus_id, items ->
                def names = items.collect { it[0] }
                def beds = items.collect { it[1] }
                tuple(locus_id, names, beds)
            }

        home_bed_by_locus_ch = EXTRACT_FLANKING.out.bed  // [locus_id, bed]
        tree_by_locus_ch = COMPUTE_TREE.out.tree         // [locus_id, tree]

        plot_inputs = home_bed_by_locus_ch
            .join(cluster_by_locus_ch)   // [locus_id, home_bed, names, beds]
            .join(tree_by_locus_ch)      // [locus_id, home_bed, names, beds, tree]
            .combine(gffs_by_locus_map_ch)
            .combine(tsvs_by_locus_map_ch)
            .map { locus_id, home_bed, names, beds, tree, gff_map, tsv_map ->
                def gffs = gff_map.get(locus_id, [])
                def tsvs = tsv_map.get(locus_id, [])
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

        log.info "${c_cyan}[PLOT] Generating synteny visualizations...${c_reset}"
            
        PLOT_SYNTENY(
            plot_inputs_split.home_bed,   // home_bed
            LOCATE_GENE.out.bed.first(),  // query_bed
            effective_home_gff_ch,        // home_gff (user-provided or Prodigal-predicted)
            plot_inputs_split.target_gffs,    // target_gffs
            plot_inputs_split.target_names,   // target_names
            plot_inputs_split.candidate_beds, // candidate_beds
            plot_inputs_split.homology_tsvs,  // homology_tsvs
            plot_inputs_split.tree,           // tree
            species_map_ch                    // species_mapping.tsv
        )
        
        PLOT_SYNTENY.out.plot.view { plot ->
            "${c_green}[PLOT] Synteny visualization complete${c_reset}"
        }
        
        // Final Reporting
        log.info "\n${c_blue}═══════════════════════════════════════════════════════════════"
        log.info "${c_white}PHASE 5: Report Generation${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        
        log.info "${c_cyan}[REPORT] Generating comprehensive report...${c_reset}"
        
        ITERATIVE_SEARCH.out.region_genes
            .map { it[1] } 
            .flatten()
            .collect()
            .ifEmpty([])
            .set { collected_regions }
            
        ITERATIVE_SEARCH.out.hits
            .map { it[1] } 
            .collect()
            .ifEmpty([])
            .set { collected_hits }
            
        // No standalone augmented proteins - integrated into iterative search
        collected_augmented = Channel.of([]).collect()
        
        GENERATE_REPORT(collected_regions, collected_hits, collected_augmented, qc_summary_ch)
        
        GENERATE_REPORT.out.report.view { report ->
            "${c_green}[REPORT] Analysis report generated successfully${c_reset}"
        }
    }
}

workflow.onComplete {
    log.info "\n${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
    if (workflow.success) {
        log.info "${c_green}Pipeline completed successfully!${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        log.info "${c_white}  Results Directory:${c_reset} ${c_cyan}${params.outdir}${c_reset}"
        log.info "${c_white}  Duration:${c_reset}          ${c_dim}${workflow.duration}${c_reset}"
        log.info "${c_white}  Tasks Completed:${c_reset}   ${c_dim}${workflow.stats.succeedCount}${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        log.info "${c_dim}  Key outputs:${c_reset}"
        log.info "${c_dim}    - synterra_report.json    (Analysis summary)${c_reset}"
        log.info "${c_dim}    - synteny_plot.pdf        (Visualization)${c_reset}"
        log.info "${c_dim}    - expanded_databases/     (Ortholog databases)${c_reset}"
        log.info "${c_dim}    - augmented_regions/      (Identified regions)${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}\n"
    } else {
        log.info "${c_red}Pipeline execution failed${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
        log.info "${c_white}  Duration:${c_reset}    ${c_dim}${workflow.duration}${c_reset}"
        log.info "${c_white}  Error:${c_reset}       ${c_red}${workflow.errorMessage}${c_reset}"
        log.info "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}\n"
    }
}
