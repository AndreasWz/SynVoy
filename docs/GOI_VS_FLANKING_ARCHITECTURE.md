# SynTerra Architecture: GOI vs Flanking Genes

## Core Concept

SynTerra uses **two types of genes** with **different roles** in the pipeline:

### 1. **Flanking Genes** (Context Markers)
- **Purpose**: Identify conserved genomic regions (synteny)
- **Characteristics**: 
  - Well-conserved across species
  - Used ONLY for region detection
  - Do NOT need to be annotated precisely in target genomes
  - Extracted from home genome annotation (or predicted with Prodigal)

### 2. **Gene of Interest (GOI)** (Discovery Target)
- **Purpose**: The actual gene we want to find and annotate across genomes
- **Characteristics**:
  - May be highly divergent
  - May be short (like melettin: 70 amino acids)
  - Requires precise annotation (exons, splice sites)
  - Needs data augmentation (variants) for divergent orthologs
  - **Gets iteratively expanded**: newly found GOIs become queries for next genome

## Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────┐
│ Phase 1: Home Genome Analysis                                   │
│                                                                  │
│  Input: Query Gene (GOI) + Home Genome                          │
│                                                                  │
│  1. LOCATE_GENE: Find GOI position in home genome               │
│  2. EXTRACT_FLANKING: Get neighboring genes (context markers)   │
│  3. PREPARE_INITIAL_DB:                                          │
│     - Add flanking genes (for synteny detection)                │
│     - Add GOI with "GOI_" prefix (for tracking)                 │
│     - Generate GOI fragments (for sensitivity)                  │
│                                                                  │
│  Output: Initial database with flanking genes + GOI             │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ Phase 2: Iterative Search (Per-Genome Wavefront)                │
│                                                                  │
│  For each target genome (phylogenetically ordered):             │
│                                                                  │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │ Step 1: SYNTENY DETECTION (Using Flanking Genes)        │   │
│  │  - MMseqs2 search with flanking genes                   │   │
│  │  - Identify conserved genomic regions                   │   │
│  │  - Cluster hits into syntenic blocks                    │   │
│  └────────────────────────┬────────────────────────────────┘   │
│                           ↓                                     │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │ Step 2: GOI ANNOTATION (In Syntenic Regions)            │   │
│  │                                                          │   │
│  │  A. Miniprot Alignment (Sensitive Mode)                 │   │
│  │     - Align GOI proteins to genomic region              │   │
│  │     - Find precise exon boundaries                      │   │
│  │     - Extract CDS sequences                             │   │
│  │                                                          │   │
│  │  B. Augmented Search (If Miniprot Fails)                │   │
│  │     - Generate GOI variants (mutations, indels)         │   │
│  │     - MMseqs2 search with variants                      │   │
│  │     - Find divergent orthologs                          │   │
│  │                                                          │   │
│  │  C. ORF-Based Fallback (If Both Fail)                   │   │
│  │     - 6-frame translation of region                     │   │
│  │     - Find best-matching ORFs                           │   │
│  │     - Simple gene prediction                            │   │
│  └────────────────────────┬────────────────────────────────┘   │
│                           ↓                                     │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │ Step 3: VALIDATION & EXPANSION                           │   │
│  │  - RBH check (validate orthology)                       │   │
│  │  - Add validated GOIs to database                       │   │
│  │  - These become queries for NEXT genome                 │   │
│  └─────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ───────────────── Repeat for Next Genome ─────────────────    │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ Phase 3: Post-Processing                                         │
│  - Phylogenetic tree construction                               │
│  - Synteny visualization                                        │
│  - Comprehensive report                                         │
└─────────────────────────────────────────────────────────────────┘
```

## Key Design Decisions

### Why Separate GOI from Flanking Genes?

1. **Different Conservation Levels**
   - Flanking genes: Housekeeping, well-conserved (e.g., ribosomal proteins)
   - GOI: May be rapidly evolving (e.g., venom peptides, immune genes)

2. **Different Annotation Requirements**
   - Flanking genes: Approximate location is enough for synteny
   - GOI: Needs precise exon/intron structure for downstream analysis

3. **Different Search Strategies**
   - Flanking genes: Standard MMseqs2 works well
   - GOI: Needs augmentation + sensitive alignment + fallbacks

### GOI-Specific Enhancements

#### 1. Data Augmentation
```python
# Generate variants of GOI for divergent search
GOI sequence:  MKFLVNVALVFMVVYISY...
Variants:      MKFLVNVSLVFMVVYISY...  (mutation)
               MKFLVNVALVF---YISY...  (deletion)
               MKFLVNVALVFAAVYISY...  (insertion)
```

#### 2. Fragment Generation
```python
# For very short proteins, create overlapping fragments
Full:     MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG
Frag 1:   MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADP  (1-35)
Frag 2:                              EPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG  (36-70)
Frag 3:   MKFLVNVALVFMVVYISY  (1-23)
...
```

#### 3. Sensitive Miniprot Parameters
```bash
miniprot \
  -I              # Output introns
  --gff           # GFF3 format
  -p 0.4          # Lower identity threshold (default 0.75)
  --aln           # Output alignment details
  -G 100000       # Max intron size
  genome.fna query.faa
```

#### 4. ORF-Based Fallback
When miniprot fails (e.g., very short or divergent GOI):
1. Extract syntenic region identified by flanking genes
2. Perform 6-frame translation
3. Find ORFs matching GOI profile
4. Use as candidates for database expansion

### Integration of Augmented Search

**Problem**: In melettin test run, augmented search found 45 candidates, but these weren't used!

**Solution**: 
1. Augmented search runs on syntenic regions identified by flanking genes
2. Candidates are annotated (simple translation or Prodigal)
3. Validated candidates (RBH check) are added to expanding database
4. These GOI sequences become queries for the next genome

**Workflow Integration**:
```groovy
// After CLUSTER_REGIONS identifies syntenic regions
AUGMENTED_SEARCH(regions, query_gene, padding)

// Annotate augmented search results
INTEGRATE_AUGMENTED_GOI(
    AUGMENTED_SEARCH.out.bed,  // Candidate regions
    genomes_dir,                // Target genome
    query_gene                  // Original GOI
)

// Add to expanding database for next genome
// (This happens inside ITERATIVE_SEARCH now)
```

## Current Implementation Status

### ✅ Implemented
- [x] GOI marking with "GOI_" prefix
- [x] GOI fragment generation
- [x] Miniprot integration in iterative search
- [x] Augmented search for divergent GOIs
- [x] RBH validation

### ⚠️ Improvements Added
- [x] Sensitive miniprot parameters
- [x] ORF-based fallback when miniprot fails
- [x] Full-length GOI preference over fragments
- [x] Integration script for augmented results

### 🔄 Needs Integration
- [ ] Feed augmented search results back to iterative search
- [ ] Use augmented GOI annotations in next genome iteration
- [ ] Compare miniprot vs ORF-based annotations
- [ ] Log which method succeeded for each genome

## Melettin Test Case Analysis

**Melettin characteristics**:
- Length: 70 amino acids (very short!)
- Type: Venom peptide (rapidly evolving)
- Structure: Signal peptide + mature peptide

**Results**:
- Flanking genes: Successfully identified syntenic regions in all 4 genomes
- Miniprot: Found 0 hits (too short/divergent)
- Augmented search: Found 45 candidates across genomes ✅
- **Issue**: Augmented results not used for annotation

**Why this happened**:
1. Iterative search runs first, uses miniprot
2. Miniprot fails on short/divergent melettin
3. Augmented search runs later, finds candidates
4. But augmented results don't feed back to database expansion
5. Result: No GOI annotations, no database growth

**Fix needed**:
- Integrate augmented results into iterative search
- When miniprot fails, use augmented candidates
- Annotate those candidates (ORF or Prodigal)
- Add to database for next genome

## Testing Recommendations

### Test with Different Gene Types

1. **Short peptide** (like melettin, 50-100aa)
   - Should trigger augmented + ORF fallback
   
2. **Medium protein** (500-1000aa)
   - Should work well with standard miniprot
   
3. **Large multi-domain** (>2000aa)
   - Should use exon-level search
   
4. **Highly divergent** (immune genes, venom)
   - Should benefit from augmentation

### Validation Steps

1. Check if augmented candidates are added to database
2. Compare miniprot vs ORF-based annotation accuracy
3. Verify GOI sequences are iteratively expanded
4. Ensure flanking genes are NOT over-annotated

## Future Enhancements

1. **Hybrid annotation**: Combine miniprot + Augustus + Prodigal
2. **Domain-aware search**: Use conserved domains for validation
3. **Structure-based validation**: Predict protein structure, compare
4. **ML-based gene prediction**: Train on known orthologs
5. **Interactive refinement**: Allow manual curation of borderline cases

---

**Last Updated**: February 5, 2026
**Authors**: SynTerra Development Team
