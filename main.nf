#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// Import Modules
include { LOCATE_GENE } from './modules/locate_gene.nf'
include { SPLIT_LOCI } from './modules/split_loci.nf'
include { EXTRACT_FLANKING } from './modules/extract_flanking.nf'
include { PREPARE_INITIAL_DB } from './modules/prepare_initial_db.nf'
include { ITERATIVE_SEARCH } from './modules/iterative_search.nf'
include { CLUSTER_REGIONS } from './modules/cluster_regions.nf'
include { AUGMENTED_SEARCH } from './modules/augmented_search.nf'
include { ANNOTATE_STRUCTURE } from './modules/annotate_structure.nf'
include { HOMOLOGY_SEARCH } from './modules/homology_search.nf'
include { PREPARE_HOME_PROTEOME } from './modules/prepare_home.nf'
include { PLOT_SYNTENY } from './modules/plot_synteny.nf'
include { COMPUTE_TREE } from './modules/compute_tree.nf'

// New Modules
include { STAGE_GENOMES } from './modules/stage_genomes.nf'
include { ASSESS_GENOME_QUALITY } from './modules/assess_quality.nf'
include { FETCH_QUERY_FROM_ID } from './modules/fetch_query.nf'
include { PHYLO_SORT } from './modules/phylo_sort.nf'
include { FETCH_RELATED_GENOMES } from './modules/fetch_related.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

// Log parameters
log.info """
    S Y N T E R R A
    ===========================
    Gene            : ${params.gene ?: params.query_id}
    Home Genome     : ${params.home_genome}
    Target Genomes  : ${params.target_genomes}
    Mode            : ${params.mode}
    """

workflow {
    log.info "SynTerra pipeline started..."
    
    // ========== INPUT VALIDATION ==========
    
    // Check if we have EITHER gene file OR query_id
    if (!params.gene && !params.query_id) { 
        error """
        ❌ ERROR: No query provided!
        """
    }
    
    if (params.gene && params.query_id) {
        log.warn "Both --gene and --query_id provided. Using --gene (${params.gene}) and ignoring ID."
    }
    
    if (!params.home_genome) { 
        error "❌ ERROR: No home genome provided!"
    }
    
    if (params.gene && !file(params.gene).exists()) {
        error "❌ ERROR: Gene file not found: ${params.gene}"
    }
    
    if (!file(params.home_genome).exists()) {
        error "❌ ERROR: Home genome file not found: ${params.home_genome}"
    }
    
    // Channel setup
    if (params.gene) {
        raw_gene_ch = Channel.fromPath(params.gene)
    } else {
        FETCH_QUERY_FROM_ID(params.query_id)
        raw_gene_ch = FETCH_QUERY_FROM_ID.out.fasta
    }
    
    raw_gene_ch.multiMap { it ->
        loc: it
        aug: it
    }.set { gene_inputs }
    
    query_gene_source_ch = gene_inputs.loc
    aug_query_gene_ch = gene_inputs.aug.first()

    home_genome_ch = Channel.fromPath(params.home_genome)
    
    // Track if user provided GFF or not
    user_provided_gff = params.home_gff ? true : false
    
    if (params.home_gff) {
        home_gff_ch = Channel.fromPath(params.home_gff).first()
    } else {
        home_gff_ch = Channel.value(file("NO_GFF"))
    }

    // Targets - Handle Easy vs Pro mode
    if (params.mode == 'easy') {
        if (!params.easy_species) {
            error "❌ ERROR: Easy mode requires --easy_species parameter!"
        }
        
        FETCH_RELATED_GENOMES(params.easy_species, params.easy_max_genomes)
        genomes_dir_ch = FETCH_RELATED_GENOMES.out.genomes_dir
        
    } else {
        if (params.target_genomes) {
            target_genomes_list = Channel.fromPath(params.target_genomes).collect()
            STAGE_GENOMES(target_genomes_list)
            genomes_dir_ch = STAGE_GENOMES.out.dir
            
        } else {
            genomes_dir_ch = Channel.empty()
        }
    }

    // PHASE 1: Core Localization
    LOCATE_GENE(query_gene_source_ch, home_genome_ch)
    
    // 5. SPLIT LOCI
    SPLIT_LOCI(LOCATE_GENE.out.bed)
    
    distinct_loci_ch = SPLIT_LOCI.out.beds.flatten()
        .map { file -> tuple(file.name, file) }
    
    // 6. Extract Flanking Genes
    EXTRACT_FLANKING(
        distinct_loci_ch, 
        home_gff_ch, 
        home_genome_ch,
        params.n_flanking_genes,
        params.min_flanking_size,
        params.prefer_large_genes
    )
        # 6b. CRITICAL FIX: Prepare Initial Database with GOI included
    // Combine flanking genes with query gene for iterative search
    PREPARE_INITIAL_DB(
        EXTRACT_FLANKING.out.faa,
        query_gene_source_ch.first()  // Use first instance of query gene
    )
        // Only run if we have targets
    if (params.target_genomes || params.mode == 'easy') {
        
        EXTRACT_FLANKING.out.faa
            .combine(genomes_dir_ch)
            .set { phylo_sort_inputs } // [locus, faa, genomes_dir]

        PHYLO_SORT(
            phylo_sort_inputs.map { tuple(it[0], it[1]) },
            phylo_sort_inputs.map { it[2] }  // genomes_dir
        )
        
        // 8. Iterative Search (FOR EACH LOCUS) - Using FIXED database with GOI
        PREPARE_INITIAL_DB.out.db
            .join(PHYLO_SORT.out.sorted_list) 
            .set { iterative_search_inputs_partial } // [locus_id, initial_db, sorted_list]
            
        iterative_search_inputs_partial
            .combine(genomes_dir_ch)
            .set { iterative_search_inputs } // [locus_id, faa, sorted_list, genomes_dir]

        // QC
        ASSESS_GENOME_QUALITY(genomes_dir_ch)
        qc_summary_ch = ASSESS_GENOME_QUALITY.out.json

        // Prepare Home Proteome (Run once)
        PREPARE_HOME_PROTEOME(home_genome_ch, home_gff_ch)
        home_proteome_db_ch = PREPARE_HOME_PROTEOME.out.db
        
        // Use predicted GFF if no user GFF was provided
        // If user provided GFF, use that; otherwise use Prodigal-generated GFF
        effective_home_gff_ch = user_provided_gff 
            ? home_gff_ch 
            : PREPARE_HOME_PROTEOME.out.gff.ifEmpty(file("NO_GFF"))

        iterative_search_inputs
            .combine(home_proteome_db_ch)
            .set { iterative_search_final_inputs } // [locus_id, faa, sorted_list, genomes_dir, home_db]

        ITERATIVE_SEARCH(
            iterative_search_final_inputs.map { tuple(it[0], it[1]) }, // [locus, faa]
            iterative_search_final_inputs.map { it[2] }, // sorted_list
            iterative_search_final_inputs.map { it[3] }, // genomes_dir
            iterative_search_final_inputs.map { it[4] }, // home_db
            params.n_flanking_genes,
            params.min_synteny_score,
            params.mmseqs_sensitivity
        )
        
        ITERATIVE_SEARCH.out.expanded_db.view { "Final Expanded DB: $it" }

        // PHASE 3: Region Identification & Augmented Search
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
            
        // Tracker for plotting
        locus_tracker_ch = joined_ch.map { unique, locus, rb, pf, gname, gd -> 
            tuple(unique, locus) 
        }

        joined_ch.multiMap { unique, locus, rb, pf, gname, gd ->
            // Inputs for processes - Pass unique_id as first val for key propagation
            aug: tuple(unique, gname, rb, gd) 
            // anno: tuple(unique, gname, rb, gd)
            // homology: tuple(unique, pf) 
        }.set { phase3_inputs }
        
        AUGMENTED_SEARCH(
            phase3_inputs.aug,
            aug_query_gene_ch, // Reusable Value Channel
            params.region_padding
        )
        
        /*
        ANNOTATE_STRUCTURE(
            phase3_inputs.anno,
            params.augustus_species
        )
        
        // Homology Search
        ANNOTATE_STRUCTURE.out.proteins
            .join(phase3_inputs.homology) 
            .map { unique_id, proteins_file, flanking_faa ->
                 tuple(unique_id, proteins_file, flanking_faa)
            }
            .set { homology_inputs }

        HOMOLOGY_SEARCH(
            homology_inputs.map { tuple(it[0], it[1]) }, 
            homology_inputs.map { it[2] } 
        )
        */
        
        // --- PHYLOGENY & PLOTTING ---
        // Collect all proteins for tree: Flanking Genes + Discovered Genes (from Expanded DB or Regions?)
        // Iterative Search output might be best source of all gene sequences found.
        // expanded_db contains everything found so far.
        // Let's use expanded_db per locus.
        
        COMPUTE_TREE(
            ITERATIVE_SEARCH.out.expanded_db
        )
        
        // PHASE 4: Miniprot-based Annotation (Replacing Augustus/HomologySearch)
        miniprot_gffs_ch = ITERATIVE_SEARCH.out.gff
            .transpose()
            .map { locus_id, gff ->
                def gname = gff.name.replace(".gff", "")
                def unique_id = "${gname}_${locus_id}"
                tuple(unique_id, gff)
            }

        miniprot_tsvs_ch = ITERATIVE_SEARCH.out.homology
            .transpose()
            .map { locus_id, tsv ->
                def gname = tsv.name.replace(".homology.tsv", "")
                def unique_id = "${gname}_${locus_id}"
                tuple(unique_id, tsv)
            }

        // Collect Data for Plotting
        plot_data_by_unique = miniprot_gffs_ch
            .join(AUGMENTED_SEARCH.out.bed)
            .join(miniprot_tsvs_ch)
            
        plot_data_with_locus = plot_data_by_unique
            .join(locus_tracker_ch) 
            
        grouped_plot_data = plot_data_with_locus
            .map { unique, gff, bed, tsv, locus ->
                tuple(locus, gff, bed, tsv, unique)
            }
            .groupTuple(by: 0) 
        
        // Simplified plot input preparation
        // Collect all files instead of complex joins
        all_gffs = miniprot_gffs_ch.map { it[1] }.collect().ifEmpty([])
        all_beds = AUGMENTED_SEARCH.out.bed.map { it[1] }.collect().ifEmpty([])
        all_tsvs = miniprot_tsvs_ch.map { it[1] }.collect().ifEmpty([])
        all_names = miniprot_gffs_ch.map { it[0] }.collect().ifEmpty([])
        
        // Get first home_bed and tree (should be same for all loci in single-locus case)
        home_bed_ch = EXTRACT_FLANKING.out.bed.map { it[1] }.first()
        tree_ch = COMPUTE_TREE.out.tree.map { it[1] }.first()
            
        PLOT_SYNTENY(
            home_bed_ch,                     // home_bed
            LOCATE_GENE.out.bed.first(),     // query_bed
            effective_home_gff_ch,           // home_gff (user-provided or Prodigal-predicted)
            all_gffs,                        // target_gffs
            all_names,                       // target_names
            all_beds,                        // candidate_beds
            all_tsvs,                        // homology_tsvs
            tree_ch                          // tree
        )
        
        // Final Reporting
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
            
        AUGMENTED_SEARCH.out.proteins
            .map { it[1] }
            .collect()
            .ifEmpty([])
            .set { collected_augmented }
        
        GENERATE_REPORT(collected_regions, collected_hits, collected_augmented, qc_summary_ch)
    }
}
