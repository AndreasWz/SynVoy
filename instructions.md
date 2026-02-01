okay das is zu naiv implementiert. schau dir das gesamte script nochmal an. schau dir auch die @organized_results an, ich denke hier passiert ein ersteer fehler. was wir wollen:
wir suchen eine GR von unserem home genome oder (nach x iterartionen) eben diese GR allgemein. das ist eine GR. jetzt wollen wir die finden, also eine list an genen in einer bestimmten reihenfolge.

dazu wollen wir folgendes beachten: es geht nicht darum, im target genome hits zu clustern und dann 10 GR zu finden, sondern eben die eine GR zu finden. dazu wollen wir für jedes gen schauen: wo haben wir die meisten CDS hits? sind das unterschiedliche CDD von dem gen (also die einzelenn CDS vom exon) (WIR WOLLEN DIE CDS VON JEDEM EXON IN ALLEN FALNKING GENES EINZUEN SUCHEN) 
und anhand von den infos wollen wir dann das gen einmalig localisieren. wenn wir also im target genome unsere flanking genes gefunden haben, schauen wir uns an, ob die GR konserviert ist, oder eben nicht. wenn es also einen cluster unserer flanking genes gibt, dann wollen wir dort unsere GR definieren.


vergleiche das, was ich geschrieben habe, mit dem aktuellen Wissenschaftlichen Stand in der Genomik. 

dann schaue die nextflow pipeline an und passe alles dementsprechend an. ! :)

kürze nichts, priorisiere nichts und überspringe den rest. nimm dir zeit. wenn es zu viel auf einmal ist, schreibe eine lsite mit allem, was getan werden muss, und arbeite erstmal nur ein paar schritte davon ab, sodass danach noch weitergearbeitet werden kann.



# Task: Implement SynTerra Pipeline in Nextflow DSL2

## Background
We are building a NEW tool for synteny-guided gene finding across genomes. 
This is NOT a migration of existing code - it's a fresh implementation of a 
novel algorithm.

## Core Algorithm (SynTerra)

### Input
- **Gene sequence** (DNA or amino acid, FASTA format)
- **Home genome** (where gene is known to exist, NCBI accession OR FASTA file)
- **Target genomes** (list of accessions OR directory of FASTA files)
- **Mode**: "easy" (auto-fetch related genomes) or "pro" (user-provided genomes)

### Output
- Interactive synteny plot (HTML) showing gene and flanking genes across all genomes
- Table of all found genes with coordinates, exon structure, synteny score
- MSA of found genes (FASTA)
- GFF3 files with annotated gene structures

### Workflow Steps

#### STEP 1: Locate Gene in Home Genome
**Input:** 
- `gene.fasta` (query gene sequence)
- `home_genome.fna` (genome FASTA)
- `home_genome.gff` (optional, genome annotation)

**Process:**
1. Use BLAST or MMseqs2 to find gene location(s) in home genome
2. If GFF provided, check if gene is already annotated
3. Handle multi-exon genes (multiple hits close together = exons)
4. Output: BED file with gene coordinates

**Output:** `home_gene_location.bed`

---

#### STEP 2: Extract Flanking Genes (Define Synteny Block)
**Input:**
- `home_gene_location.bed` (from Step 1)
- `home_genome.gff` (genome annotation)

**Process:**
1. Get genomic region (GR) = gene ± flanking region
2. Extract N genes upstream and N genes downstream (default N=10)
3. **Prioritize large genes** (>2kb) as they're better synteny markers
4. If large gene spans position -15 to +5, include it even though it's "far"
5. Extract protein sequences for all flanking genes
6. Store as "Synteny Block Database" (SBD)

**Parameters:**
- `n_flanking_genes`: 10 (default)
- `min_flanking_gene_size`: 500 bp (minimum)
- `prefer_large_genes`: true (prioritize >2kb genes)

**Output:** 
- `synteny_block.bed` (coordinates of all genes in block)
- `flanking_proteins.faa` (protein sequences for flanking genes)

---

#### STEP 3: Iterative Genome Search (Phylogenetic Order)
**This is the KEY innovation of SynTerra**

**Input:**
- `flanking_proteins.faa` (initial search database)
- List of target genomes (sorted by phylogenetic distance)

**Process:**
1. **Sort genomes by evolutionary distance** from home genome
   - Use NCBI taxonomy API to get phylogenetic tree
   - OR use AAI (Average Amino acid Identity) if genomes provided
   - Closest species first, most distant last

2. **For each genome (in order):**
   a. Search flanking proteins using DIAMOND (faster than MMseqs2 for proteins)
```bash
      diamond blastp --query flanking_proteins.faa \
                     --db target_genome.dmnd \
                     --outfmt 6 \
                     --evalue 1e-5 \
                     --max-target-seqs 50
```
   
   b. Find hits and extract hit sequences from genome
   
   c. **ADD high-quality hits to search database** for next iteration
      - Filter: E-value < 1e-10, alignment length > 100 aa, identity > 40%
      - This allows finding distant orthologs in later species
   
   d. Update `flanking_proteins.faa` with new sequences

3. **Why iterative?**
   - Direct search in distant species fails (too divergent)
   - Using intermediate species as "stepping stones" improves sensitivity
   - Example: Home → Species1 (90% ID) → Species2 (80% ID) → Species3 (60% ID)
   - Direct Home → Species3 would fail, but iterative succeeds

**Parameters:**
- `diamond_sensitivity`: "very-sensitive"
- `min_hit_identity`: 40%
- `min_hit_length`: 100 aa
- `expand_db_threshold`: 1e-10 (E-value cutoff for adding to DB)

**Output:**
- `{genome}_flanking_hits.tsv` (DIAMOND results for each genome)
- `expanded_flanking_db.faa` (accumulated database after all iterations)

---

#### STEP 4: Identify Genomic Regions (GRs) via Clustering
**Input:**
- `{genome}_flanking_hits.tsv` (from Step 3)
- Target genome FASTA

**Process:**
1. For each genome, cluster DIAMOND hits by genomic proximity
   - Hits within 50kb of each other = same Genomic Region (GR)
   - Use DBSCAN or simple distance-based clustering

2. Score each GR by synteny conservation:
```
   synteny_score = (# of conserved flanking genes) / (total flanking genes)
   
   Example:
   Home genome: 10 upstream + gene + 10 downstream = 21 genes
   Target GR: Found 7 upstream + ??? + 8 downstream = 15/20 flanking genes
   Synteny score = 15/20 = 0.75 (good!)
```

3. Filter GRs by minimum synteny score (default 0.6 = 60% conservation)

4. **Handle duplications:**
   - Multiple GRs per genome possible (gene duplications)
   - Report all GRs with synteny_score > threshold

**Parameters:**
- `cluster_distance`: 50000 bp (50kb)
- `min_synteny_score`: 0.6 (60%)

**Output:**
- `{genome}_genomic_regions.bed` (coordinates of all candidate GRs)
- `{genome}_synteny_scores.tsv` (synteny score for each GR)

---

#### STEP 5: Augmented Gene Search in GRs
**This is where we find the actual gene (not just flanking genes)**

**Input:**
- `gene.fasta` (original query gene)
- `{genome}_genomic_regions.bed` (from Step 4)
- Target genome FASTA

**Process:**
1. For each GR, extract genomic sequence (GR start - 10kb to GR end + 10kb)

2. **Generate augmented gene variants** for sensitive search:
   
   a. **Splice variants** (if gene has multiple exons):
      - Try all possible exon combinations
      - Example: 3 exons → search [E1], [E2], [E3], [E1+E2], [E2+E3], [E1+E2+E3]
   
   b. **Frameshift tolerance**:
      - Search in all 6 reading frames
      - Allow +1/-1 frameshifts (common in pseudogenes)
   
   c. **Evolutionary divergence simulation**:
      - Generate 10 variants with 5% random mutations
      - Simulate codon substitutions (not just any mutation)
      - Use codon substitution matrix (e.g., allow synonymous changes more)

3. **High-sensitivity search** in GR:
```bash
   # Use MMseqs2 translated search with very high sensitivity
   mmseqs easy-search gene_variants.faa \
                      GR_sequence.fna \
                      hits.m8 \
                      tmp \
                      --search-type 2 \        # translated
                      -s 8.5 \                 # very high sensitivity
                      --min-seq-id 0.3 \       # accept 30% identity
                      -e 0.01                  # relaxed E-value
```

4. Report ALL candidates (even weak hits) with scores

**Parameters:**
- `enable_splice_variants`: true
- `enable_frameshifts`: true
- `mutation_rate`: 0.05 (5%)
- `num_mutant_variants`: 10
- `mmseqs_sensitivity`: 8.5
- `min_identity`: 30%

**Output:**
- `{genome}_gene_candidates.bed` (all gene hits in all GRs)
- `{genome}_gene_sequences.faa` (extracted sequences)

---

#### STEP 6: Annotate Gene Structure (Exon/Intron Boundaries)
**Input:**
- `{genome}_gene_candidates.bed` (from Step 5)
- Target genome FASTA
- Target genome GFF (if available)

**Process:**
1. **If gene already annotated in target genome:**
   - Search GFF for gene name (use synonym matching)
   - Use existing annotation as candidate

2. **If not annotated:**
   - Use SpliceAI to predict exon/intron boundaries
   - Or use AUGUSTUS/BRAKER3 with gene as evidence
   
3. **Validate predicted structure:**
   - Check for start codon (ATG)
   - Check for stop codon in correct frame
   - Check for canonical splice sites (GT-AG)
   - Flag pseudogenes (premature stop, frameshift)

4. **Extract final protein sequence**

**Output:**
- `{genome}_genes_annotated.gff` (GFF3 with exon structure)
- `{genome}_proteins.faa` (translated proteins)
- `{genome}_genes_status.tsv` (functional/pseudogene status)

---

#### STEP 7: Visualize Synteny Across All Genomes
**Input:**
- All `{genome}_genes_annotated.gff` files
- All `{genome}_genomic_regions.bed` files
- All synteny scores

**Process:**
1. Create interactive HTML plot with Plotly:
   - Each genome = horizontal track
   - Genes = colored boxes (color by orthology)
   - Flanking genes = gray boxes
   - Query gene = red box
   - Lines connecting orthologous genes across genomes

2. **Features:**
   - Hover: show gene name, product, coordinates, synteny score
   - Click: highlight ortholog group
   - Zoom: focus on specific region
   - Export: PDF/SVG for publication

3. **Layout:**
   - Sort genomes by phylogenetic distance
   - Align by query gene position
   - Show synteny breaks (rearrangements)

**Output:**
- `synteny_plot.html` (interactive plot)
- `synteny_plot.pdf` (static plot for paper)

---

#### STEP 8: Generate Summary Outputs
**Process:**
1. Align all found genes with MAFFT
2. Build phylogenetic tree (optional)
3. Create summary table

**Outputs:**
- `all_genes_alignment.fasta` (MSA)
- `all_genes_tree.nwk` (phylogeny)
- `summary_table.csv`:
```
   genome,gene_id,start,end,strand,exons,synteny_score,status
   conus_consors,gene_001,123456,125678,+,3,0.85,functional
   conus_betulinus,gene_002,234567,236789,-,3,0.75,functional
   pomacea,gene_003,345678,347890,+,2,0.65,pseudogene
```

---

## Implementation Requirements

### Nextflow Structure
```
syntenyfinder/
├── main.nf                      # Entry point
├── nextflow.config              # Configuration
├── modules/
│   ├── locate_gene.nf          # Step 1
│   ├── extract_flanking.nf     # Step 2
│   ├── iterative_search.nf     # Step 3 (KEY)
│   ├── cluster_regions.nf      # Step 4
│   ├── augmented_search.nf     # Step 5 (KEY)
│   ├── annotate_structure.nf   # Step 6
│   └── plot_synteny.nf         # Step 7
├── bin/
│   ├── phylo_sort.py           # Sort genomes by distance
│   ├── generate_variants.py    # Create augmented gene variants
│   ├── cluster_grs.py          # Cluster genomic regions
│   └── plot_synteny.py         # Create interactive plot
└── conf/
    ├── base.config
    └── test.config
```

### Key Parameters (nextflow.config)
```groovy
params {
    // Input
    gene = null
    home_genome = null
    target_genomes = null
    mode = 'pro'  // 'easy' or 'pro'
    
    // Synteny parameters
    n_flanking_genes = 10
    prefer_large_genes = true
    min_flanking_size = 500
    
    // Search parameters
    cluster_distance = 50000
    min_synteny_score = 0.6
    
    // Augmentation
    enable_splice_variants = true
    enable_frameshifts = true
    mutation_rate = 0.05
    
    // Sensitivity
    mmseqs_sensitivity = 8.5
    min_gene_identity = 30
}
```

### Error Handling
- Fail gracefully if gene not found in home genome
- Warn if synteny score < threshold (but still report)
- Handle genomes with no annotation (GFF missing)
- Handle multi-copy genes (duplications)
- Report pseudogenes separately

### Testing
Include test dataset:
- 1 query gene (conotoxin, ~100bp)
- 1 home genome (Conus geographus)
- 3 target genomes (2 Conus + 1 outgroup)
- Expected: Find gene in 2 Conus with synteny >0.7

---

## CRITICAL: What NOT to do

❌ DO NOT migrate existing Snakemake code
❌ DO NOT keep broken download/search/synteny rules
❌ DO NOT use the fake ML classifier
❌ DO NOT use Clinker (we're making custom plots)

✅ DO implement the algorithm described above
✅ DO write clean, new Nextflow processes
✅ DO focus on the iterative search (Step 3)
✅ DO focus on augmented search (Step 5)

---

## Questions for Clarification
1. Should easy mode auto-download from NCBI or just fail if not provided?
2. What if gene has no annotated orthologs in any target genome?
3. Should we report ALL GRs or just the best one per genome?
4. Output format preference for plots (just HTML or also PDF)?

---

Please implement this pipeline following the steps above. Start with 
Steps 1-2 (locate gene + extract flanking), then implement the critical 
Step 3 (iterative search), then Steps 4-5 (clustering + augmented search).

The visualization (Steps 6-7) can be added last.
```

---

## 🎯 WHY THIS PROMPT IS CORRECT

### ✅ **It describes the ALGORITHM, not code migration**
- Explains WHY each step exists
- Gives biological motivation
- Specifies parameters with reasoning

### ✅ **It's complete and specific**
- AI knows exactly what each process should do
- No ambiguity about inputs/outputs
- Clear success criteria

### ✅ **It prevents common mistakes**
- Explicitly says "DON'T migrate old code"
- Focuses on novel parts (iterative search, augmentation)
- Avoids broken approaches (Clinker, fake ML)

### ✅ **It's implementable**
- Realistic timeline (each step is 1-2 days)
- Uses standard tools (DIAMOND, MMseqs2, SpliceAI)
- Testable at each stage


early sketch:

┌─────────────────────────────────────────────────────────┐
│ INPUT                                                   │
├─────────────────────────────────────────────────────────┤
│ • Gene sequence (DNA/AA)                               │
│ • Home genome (accession or FASTA)                     │
│ • Target genomes (list of accessions OR folder)        │
│ • Outgroup (optional, auto-selected if Easy Mode)      │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 1: Locate Gene in Home Genome                     │
├─────────────────────────────────────────────────────────┤
│ • BLAST/MMseqs2: Find gene position                    │
│ • Check GFF annotation (if available)                  │
│ • Handle multi-exon genes                              │
│ • Extract Genomic Region (GR): gene ± flanking genes   │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 2: Extract Flanking Genes (Synteny Markers)       │
├─────────────────────────────────────────────────────────┤
│ • Get 10 upstream + 10 downstream genes                │
│ • Priority: Large genes (>2kb) are better markers      │
│ • Store as "Synteny Block Database" (SBD)              │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 3: Iterative Genome Search (Phylogenetic Order)   │
├─────────────────────────────────────────────────────────┤
│ for genome in sorted_by_distance(target_genomes):      │
│   1. Search flanking genes (DIAMOND/MMseqs2)           │
│   2. Cluster hits → identify Genomic Regions (GRs)     │
│   3. Add found genes to SBD (for next iteration)       │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 4: Sensitive Gene Search in GRs                   │
├─────────────────────────────────────────────────────────┤
│ • High sensitivity search (MMseqs2 -s 8.5)             │
│ • Data augmentation:                                   │
│   - Splice variants (try all exon combinations)       │
│   - Allow frameshifts (pseudogenes)                    │
│   - Codon-aware mutations (simulate divergence)       │
│ • Report all candidates (even weak hits)               │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ STEP 5: Annotation & Validation                        │
├─────────────────────────────────────────────────────────┤
│ • Predict exon/intron boundaries (SpliceAI/BRAKER)    │
│ • Validate ORFs (check start/stop codons)             │
│ • Calculate synteny score (% flanking genes conserved) │
│ • Flag duplications (multiple GRs per genome)          │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ OUTPUT                                                  │
├─────────────────────────────────────────────────────────┤
│ • Interactive synteny plot (Clinker-style)             │
│ • MSA of found genes (MAFFT)                           │
│ • Phylogenetic tree (optional)                         │
│ • Summary table (CSV)                                  │
└─────────────────────────────────────────────────────────┘


also check if everyting makes sense and is best practise



basic idea:

input pro mode: gene sequence (DNA, or amino) and home genome of the gene, genome ncbi accesion ids or folder with all genome fasta files in there
input easy mode: gene sequence (DNA or amino) and home genome of the gene. maybe in future like this: just name of uniprot id for gene and genome name, it does the rest

easy mode searches on ncbi for close genomes, also takes an outgroup
pro mode searches in genomes given

1. searches for gene in home genome, finds it (if multi exon u will find multiple parts of it near each other, but diff seq it, u can annotate exon boundaries).
2. if genome annotated takes flanking genes (~10 or so, if a big gene is at position -15 until +15 still take it, as bigger genes are better markers in Syntheny) (we call the genomic region GR, that is flanking genes with the gene)
3. search the flanking genes in the other genomes, its good to do it iteravly, start with the clostest species, do everyting there, than add the found stuff to search db and go on, otherwise it gets hard to find anything in very evolutionary far species
4. identify correct position of the GR in new genome, then search for gene in that GR (very sensitivly). use data augmentation for the gene (cropping, splicing, random mutations, etc)   (! multiple GRs per genome possible!)
5. annotate all fragments and genes found (exon intron boundaries!)
6. make a plot that shows the GR across all genomes showing the flanking genes and the located genes.


if multiple genes given, run this workflow for each gene, but check if u can merge GRs at the end (e.g. if they have same flanking genes etc)

hints: if gene is already annotated in home genome, search for name in other genomes and take that as candidate Location.

## Implementation Status (Jan 2026)

### Phase 4: Annotation & Visualization
- **Annotation**: Uses `AUGUSTUS` (species: `honeybee1` by default) to predict gene structures in identified regions.
- **Visualization**: `results/synteny_plot.html` is an interactive Plotly HTML file.
  - **Synteny Anchors**: Displayed as colored blocks (blue/orange/green) based on genome track.
  - **Candidate Genes**: The specific target gene found by the **Augmented Search** (Step 5) is highlighted in **RED** on the plot. This allows users to pinpoint the divergent ortholog within the syntenic neighborhood.
- **Pipeline Cleanup**: Intermediate files stored in `work/` can be cleaned up; `results/` contains the final outputs.