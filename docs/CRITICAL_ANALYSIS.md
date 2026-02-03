# SynTerra: Critical Analysis & Architectural Review

**Date**: February 3, 2026  
**Reviewer**: Technical Analysis  
**Status**: 🔴 Major Logic Flaws Identified - Requires Fundamental Redesign

---

## Executive Summary

SynTerra is an ambitious bioinformatics pipeline for synteny-guided gene discovery. While the **core concept is scientifically sound** and addresses a real need (finding divergent orthologs using genomic context), the **implementation has critical architectural flaws** that fundamentally undermine its primary objective.

### The Fatal Flaw

**The pipeline searches for flanking genes but fails to properly search for the gene of interest (GOI) during the iterative search phase.** This is the biological equivalent of building a treasure map to find landmarks around the treasure, but forgetting to include the treasure itself on the map.

### Severity Classification

- 🔴 **Critical (Pipeline-Breaking)**: 3 issues
- 🟠 **Major (Logic Errors)**: 7 issues  
- 🟡 **Moderate (Design Weaknesses)**: 5 issues
- 🔵 **Minor (Code Quality)**: 8 issues

---

## Table of Contents

1. [Critical Issues (Pipeline-Breaking)](#1-critical-issues-pipeline-breaking)
2. [Major Logic Errors](#2-major-logic-errors)
3. [Moderate Design Weaknesses](#3-moderate-design-weaknesses)
4. [Minor Code Quality Issues](#4-minor-code-quality-issues)
5. [Positive Aspects](#5-positive-aspects)
6. [Recommended Architecture](#6-recommended-architecture)
7. [Implementation Priority](#7-implementation-priority)

---

## 1. Critical Issues (Pipeline-Breaking)

### 🔴 **CRITICAL 1: The Query Gene is NOT Searched Iteratively**

**Location**: `bin/iterative_search_runner.py:530-545`

**The Problem**:
```python
# Lines 530-545 in process_single_genome()
relevant_hits = [h for h in hits if h['chrom'] == chrom]
unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)

# Extract sequences from DB
for query_id in unique_queries:
    if query_id in db_sequences:
        found_queries.append(db_sequences[query_id])
```

**What Actually Happens**:
1. Search for flanking genes in target genome ✅
2. Identify synteny region based on flanking gene hits ✅
3. Extract the synteny region ✅
4. **Extract ONLY the flanking genes that had hits** ❌
5. Run Miniprot with ONLY flanking genes ❌
6. **The original query gene is NEVER included in the Miniprot search** ❌

**Why This Breaks Everything**:
- The entire point of the pipeline is to find the query gene in divergent genomes
- Iterative search should progressively use newly found query genes as templates for the next genome
- Currently, it only finds flanking genes iteratively, not the actual target gene
- The pipeline relies 100% on "augmented search" to find the query gene, but this only works if regions are correctly identified

**Visual Representation**:
```
HOME GENOME: [FlankL2][FlankL1][QUERY_GENE][FlankR1][FlankR2]
                                    ↓
INITIAL DB:  [FlankL2][FlankL1] [MISSING!] [FlankR1][FlankR2]
                                    ↓
TARGET SEARCH: MMseqs → Find flanking genes → Define region
                                    ↓
MINIPROT SEARCH: [FlankL2][FlankL1] vs Region  ← Wrong queries!
                 Should be: [FlankL2][FlankL1][QUERY_GENE] vs Region
```

**Impact**: **CATASTROPHIC** - Pipeline cannot fulfill its primary function.

---

### 🔴 **CRITICAL 2: Miniprot Step in Iterative Search is Redundant**

**Location**: `bin/iterative_search_runner.py:550-580`

**The Problem**:
After finding flanking genes with MMseqs2 and identifying synteny regions, the code runs Miniprot to re-annotate genes in those regions. However:

```python
# Prepare queries for Miniprot
found_queries = []  # Only contains flanking genes that had MMseqs hits

# Run Miniprot
miniprot_hits = run_miniprot(temp_fa, query_mini_fa, miniprot_paf)
```

**What's Wrong**:
1. Miniprot searches the region with **the same flanking genes that MMseqs already found**
2. This just re-finds the same genes with slightly better exon boundaries
3. **The query gene is not included**, so Miniprot can't find it either
4. The step is computationally expensive but adds minimal value

**What It Should Do**:
1. Use Miniprot to search for **all** genes (including the query gene)
2. Use MMseqs hits only to identify the synteny region
3. Then use Miniprot for comprehensive annotation of everything in that region

**Impact**: **Wasted computation** and **missed opportunity** to find the query gene with splice-aware alignment.

---

### 🔴 **CRITICAL 3: Empty Region Problem - Silent Failure Mode**

**Location**: `bin/cluster_grs.py:280-300`, `bin/augmented_search_runner.py:65-75`

**The Problem**:
If synteny scoring is too stringent, no regions are passed to augmented search, causing silent failure:

```python
# cluster_grs.py
passes_score = best['score'] >= (args.min_score * 0.5)

if passes_score:
    f_out.write(f"{best['chrom']}\t{best['start']}...")
# else: Write NOTHING → Empty BED file

# augmented_search_runner.py  
regions = extract_regions(args.regions_bed, genome_file, padding)

if not regions:
    print("No regions extracted. Exiting.")
    # Create empty outputs
    open(f"{args.output_base}.bed", 'w').close()
    return  # SILENT FAILURE - no error raised!
```

**Why This is Critical**:
1. If only 2 out of 10 flanking genes are conserved → synteny score too low → no region
2. Empty BED → Augmented search has nothing to search
3. Pipeline continues normally, no error thrown
4. User gets empty results with no explanation why

**Evidence from Your Own Config**:
```python
# conf/test_tetramorium.config
params {
    min_synteny_score = 0.01  # You lowered to 1% because you hit this problem!
}
```

**Impact**: **Pipeline fails silently** for moderately divergent genomes (50-70% flanking gene conservation).

---

## 2. Major Logic Errors

### 🟠 **MAJOR 1: Insufficient Search Region Padding**

**Location**: `bin/iterative_search_runner.py:510-515`, `bin/augmented_search_runner.py:25`

**The Problem**:
```python
# iterative_search_runner.py - padding for Miniprot region
w_start = max(0, best_region['start'] - 20000)  # Only ±20kb!
w_end = min(slen, best_region['end'] + 20000)

# But augmented_search_runner.py uses
parser.add_argument("--padding", type=int, default=10000)  # 10kb default
# Overridden in config to 50kb, but still inconsistent
```

**Why This is Wrong**:
1. Gene spacing varies widely (5kb to >100kb between genes)
2. If query gene is 30kb upstream of first flanking gene → outside search window
3. **The adaptive padding function exists but uses wrong default**:
```python
def calculate_adaptive_padding(hits, best_region, default=100000):  # Says 100kb
    # But is called with no default override, uses 20kb in actual extraction!
```

**Impact**: Query gene missed if it's outside the narrow window around flanking genes.

---

### 🟠 **MAJOR 2: Exon Mode Implementation is Incomplete**

**Location**: `bin/extract_flanking_genes.py:230-280`

**The Problem**:
The pipeline has exon-level fragmentation code for flanking genes, but:

```python
parser.add_argument("--exon_mode", type=str, default="false")  # Disabled by default
```

And more critically:
- Exon mode works for flanking genes extraction ✓
- But iterative search doesn't properly handle exon-level queries ✗
- Fragment generation for query gene (planned feature) is not integrated ✗

**From Your Own Documentation** (IMPLEMENTATION_PLAN.md):
> "Desired Implementation: Output EACH EXON as separate CDS sequence"  
> **Status**: Partially implemented but not tested/integrated

**Why This Matters**:
- Exon-level search is crucial for finding divergent genes (as you noted)
- If exon 3 is conserved but exons 1-2 are lost → gene can still be detected
- Current implementation would miss such cases

**Impact**: Reduced sensitivity for detecting partial genes, pseudogenes, and frameshifted orthologs.

---

### 🟠 **MAJOR 3: Query Gene Fragmentation Not Integrated**

**Location**: `bin/generate_variants.py:120-150`, workflow missing

**The Problem**:
Your plan mentions fragmenting the query gene:
> "Search WHOLE query sequence, HALVES, THIRDS, QUARTERS"

The code exists (`bin/fragment_query.py`), but:
1. **Not called anywhere in the workflow** ❌
2. `generate_variants.py` has fragment generation code but it's not invoked in the pipeline
3. Only used for data augmentation (random mutations), not systematic fragmentation

**Current Reality**:
```python
# bin/augmented_search_runner.py
run_command([
    "python3", gen_var_script,
    "--query", args.query_gene,
    "--output", variants_faa,
    "--mutation_rate", "0.05",
    "--num_variants", "10"  # Only 10 mutated variants, no fragments!
])
```

**What's Missing**:
- No systematic halves/thirds/quarters generation
- Fragment hits are not consolidated back to full gene predictions
- No metadata tracking which fragment found which hit

**Impact**: Can't detect partial genes, cleaved propeptides, or domain-shuffled variants.

---

### 🟠 **MAJOR 4: Homology Validation (RBH) Logic is Flawed**

**Location**: `bin/iterative_search_runner.py:420-445`

**The Problem**:
```python
# RBH check
parent = cand_map[cand_id]
parent_base = extract_base_gene_id(parent).strip()
target_base = extract_base_gene_id(target_id).strip()

if parent_base == target_base or parent == target_id or target_id == parent:
    valid_ids.add(cand_id)
```

**Issues**:
1. **String matching for gene names is fragile**: "gene-1" vs "gene-11" problem (partially fixed but still risky)
2. **RBH back to home genome only**: Doesn't validate against other confirmed orthologs in phylogenetic order
3. **No domain validation**: Should use InterProScan/Pfam to confirm conserved domains
4. **E-value threshold too relaxed**: `evalue=1e-5` allows weak hits to pass

**From Your Documentation** (ARCHITECTURE_PROPOSAL.md):
> "3. VALIDATE: RBH Check + Domain Annotation  
>    - InterProScan/Pfam für Domänen  
>    - Nutze Domänen als zusätzlichen Filter"

**Status**: NOT implemented!

**Impact**: False positives from paralogous genes or unrelated genes with similar names.

---

### 🟠 **MAJOR 5: Phylogenetic Ordering is Disconnected from Iterative Expansion**

**Location**: `main.nf:149-158`, `bin/phylo_sort.py`

**The Problem**:
```groovy
// main.nf - PHYLO_SORT happens BEFORE iterative search
PHYLO_SORT(...)
ITERATIVE_SEARCH(...)  // Uses sorted genomes
```

But then:
```python
# iterative_search_runner.py
# Defines waves based on distance, processes in parallel WITHIN waves
waves = []
# ...
for i, wave in enumerate(waves):
    # All genomes in a wave are processed in parallel
```

**Why This is Suboptimal**:
1. **Parallelization breaks strict phylogenetic order**: Genomes at distance 0.15 and 0.151 are in the same wave
2. **No database updates between close genomes**: Should update DB after each genome (strict serial) for maximum sensitivity
3. **Wave definition is arbitrary**: Distance threshold 0.001 has no biological justification

**Better Approach**:
- Process strictly in order for closest genomes (distance < 0.1)
- Only parallelize for distant genomes where order doesn't matter
- Update database after EACH genome in close-species mode

**Impact**: Reduced effectiveness of iterative search for closely related species.

---

### 🟠 **MAJOR 6: Cluster Distance Parameter is Auto-Detected but Unreliable**

**Location**: `bin/iterative_search_runner.py:362-375`

**The Problem**:
```python
def estimate_cluster_dist(genome_file, default_dist=50000):
    size = os.path.getsize(genome_file)
    if size < 10_000_000: # < 10MB
        return 20000
    elif size > 100_000_000: # > 100MB
        return 100000
    return default_dist
```

**Why This is Wrong**:
1. **File size ≠ gene density**: A 200 MB genome could be:
   - Mammal (sparse, ~20kb/gene) → need 100kb clustering
   - Highly AT-rich plant (dense) → need 30kb clustering
2. **Ignoring GFF information**: If GFF provided, could calculate actual gene spacing
3. **No per-chromosome adjustment**: Some chromosomes are gene-dense, others sparse

**Better Approach**:
- Parse actual gene coordinates from GFF/predictions
- Calculate median inter-gene distance
- Use 2-3× median as cluster distance

**Impact**: Sub-optimal synteny block detection, especially for unusual genomes (polyploid plants, highly repetitive).

---

### 🟠 **MAJOR 7: No Pseudogene Detection Despite Documentation Claims**

**Location**: Throughout pipeline, especially `README.md:14`

**Documentation Claims**:
> "4. Handles difficult cases: Finds **pseudogenes**, partial genes, and highly divergent orthologs"

**Reality Check**:
❌ No explicit pseudogene detection logic  
❌ No frameshift detection beyond Miniprot's basic reporting  
❌ No stop codon analysis  
❌ No comparison of exon structure to reference  
❌ No reporting of truncated genes  

**What Would Be Needed**:
1. Detect premature stop codons in alignments
2. Identify frameshifts (Miniprot reports these but pipeline doesn't parse)
3. Compare exon counts: if reference has 5 exons but target has 2 → pseudogene
4. Check for transposon insertions (common pseudogenization mechanism)
5. Flag genes with low coverage (<50% of reference length)

**Impact**: **Misleading documentation** - feature is claimed but not implemented.

---

## 3. Moderate Design Weaknesses

### 🟡 **MODERATE 1: Wasteful Recomputation of GFF Parsing**

**Location**: Multiple Python scripts parse GFF independently

**The Problem**:
```python
# extract_flanking_genes.py has its own parse_gff()
# plot_synteny.py has its own parse_gff()  
# sequence_utils.py has parse_gff() but not always used
# Each slightly different!
```

**Impact**: Code duplication, inconsistent behavior, maintenance burden.

**Fix**: Centralize in `sequence_utils.py` and import everywhere.

---

### 🟡 **MODERATE 2: Miniprot GFF Parsing is Fragile**

**Location**: `bin/iterative_search_runner.py:260-335`

**The Problem**:
```python
# Lots of try-except blocks because parsing is unreliable
try:
    hit = {
        'start': int(parts[3]),
        'identity': float(info.get('Identity', 0)) * 100,
        # ...
    }
except (ValueError, IndexError) as e:
    print(f"Warning: Failed to parse mRNA line: {e}", file=sys.stderr)
    continue
```

**Why This Happens**:
- Miniprot GFF format slightly non-standard
- Custom attributes (`SynTerra_Parent`, `SynTerra_ID`) are added but not consistently parsed
- No validation of GFF before processing

**Better Approach**:
- Use a proper GFF parsing library (`gffutils` or `bcbio-gff`)
- Validate GFF structure before parsing
- Standardize attribute naming conventions

---

### 🟡 **MODERATE 3: No Resume Capability for Long Pipelines**

**Location**: Workflow orchestration (Nextflow handles this, but...)

**The Problem**:
- Long pipelines (50+ genomes) can take days
- If failure occurs at genome #47, restart from scratch
- Nextflow caching helps but doesn't persist well across major changes

**Impact**: Wasted compute time during development/debugging.

**Better Approach**:
- Checkpoint after each genome in iterative search
- Save intermediate databases
- Allow explicit restart from checkpoint

---

### 🟡 **MODERATE 4: MMseqs2 Sensitivity Hardcoded in Multiple Places**

**Location**: Various scripts have different defaults

**The Problem**:
```python
# iterative_search_runner.py:483
"-s", "7.5",  # Hardcoded!

# augmented_search_runner.py:105
"-s", args.mmseqs_sens,  # Uses parameter
```

**Inconsistency**:
- Iterative search: 7.5 (hardcoded)
- Augmented search: 8.5 (from config)
- Should be same or configurable

---

### 🟡 **MODERATE 5: Tree Computation is Disconnected from Results**

**Location**: `modules/compute_tree.nf`, `main.nf:277`

**The Problem**:
```groovy
COMPUTE_TREE(
    ITERATIVE_SEARCH.out.expanded_db
)
```

Tree is computed from expanded database (all found genes), but:
1. **Not used for any downstream analysis** beyond plotting colors
2. Could be used for phylogenetic validation (e.g., flag genes on wrong branches)
3. Could inform search strategy (prioritize clades with better gene conservation)

**Impact**: Underutilized valuable information.

---

## 4. Minor Code Quality Issues

### 🔵 **MINOR 1: Inconsistent Coordinate Systems**

**Location**: Throughout codebase

**The Problem**:
```python
# GFF: 1-based, closed interval [start, end]
# BED: 0-based, half-open interval [start, end)
# Python slices: 0-based, half-open
# MMseqs output: Who knows? (depends on format)
```

**Code Comments Show Confusion**:
```python
# extract_flanking_genes.py:48-51
# "GFF3 Standard: 1-based, closed interval [start, end].
#  Python/BED Standard: 0-based, half-open interval [start, end).
#  Conversion: BED Start = GFF Start - 1, BED End = GFF End"
```

But then conversions are inconsistent across files.

**Impact**: Off-by-one errors in coordinate transformations (latent bugs).

---

### 🔵 **MINOR 2: Magic Numbers Throughout Code**

**The Problem**:
```python
# bin/iterative_search_runner.py
max_intron=20000  # Why 20kb?
cluster_dist=50000  # Why 50kb?
padding=20000  # Why 20kb?
MAX_GENE_SPAN = 500000  # Why 500kb?

# bin/cluster_grs.py
BASE_QUALITY_WEIGHT = 0.4  # Why 0.4?
CONSISTENCY_WEIGHT = 0.3  # Why 0.3?
```

No documentation explaining these thresholds.

**Better**: Move to config file with documentation.

---

### 🔵 **MINOR 3: Verbose Debug Output Not Cleaned Up**

**Location**: Multiple scripts still have print statements

**Example**:
```python
# bin/iterative_search_runner.py:555
print(f"[{genome_name}] CMD: miniprot -I --gff {temp_fa} {query_mini_fa}", flush=True)
```

**Impact**: Cluttered logs, hard to parse for real errors.

**Fix**: Use proper logging module with levels (DEBUG, INFO, WARNING, ERROR).

---

### 🔵 **MINOR 4: No Input Validation in Many Scripts**

**Example**:
```python
# Many scripts don't check:
# - Does file exist?
# - Is file empty?
# - Is file format correct?
# - Are parameters in valid ranges?
```

**Impact**: Cryptic error messages when input is malformed.

---

### 🔵 **MINOR 5: Inconsistent Error Handling**

**The Problem**:
```python
# Some functions use try-except
# Others use subprocess.run(..., check=True)
# Others silently fail and return empty lists
# No consistent error propagation strategy
```

---

### 🔵 **MINOR 6: Unused Imports and Functions**

**Location**: Multiple files

**Example**:
```python
# bin/iterative_search_runner.py:18-28
try:
    from fragment_query import ...
    FRAGMENT_SUPPORT = True
except ImportError:
    FRAGMENT_SUPPORT = False

# But then FRAGMENT_SUPPORT is never checked!
```

---

### 🔵 **MINOR 7: No Type Hints**

**Impact**: Hard to understand function signatures, easier to make mistakes.

**Modern Best Practice**: Use type hints for all function arguments and returns.

---

### 🔵 **MINOR 8: Test Data is Minimal**

**Location**: `test_data/` directory

**Problem**:
- Only two test configs (melettin, tetramorium)
- No unit tests for individual Python functions
- No integration tests with known-good outputs
- No regression tests

**Impact**: Hard to verify fixes don't break existing functionality.

---

## 5. Positive Aspects

Despite the critical issues, SynTerra has strong foundations:

### ✅ **Excellent Core Concept**
The idea of using synteny as evidence for gene prediction is scientifically sound and innovative.

### ✅ **Comprehensive Documentation**
- Multiple detailed markdown files (README, ARCHITECTURE_PROPOSAL, IMPLEMENTATION_PLAN)
- Clear documentation of intended features
- Good user-facing documentation (USAGE.md)

### ✅ **Well-Structured Codebase**
- Modular Nextflow workflow
- Separated Python scripts for each task
- Clean directory structure

### ✅ **No External Dependencies Beyond Standard Bioinformatics Tools**
- Self-contained FASTA/GFF parsing (no BioPython dependency)
- Uses standard tools (MMseqs2, Miniprot, Prodigal)

### ✅ **Robust Sequence Utilities**
`bin/sequence_utils.py` is well-written with:
- Multiple ID extraction patterns
- Fallback handling
- Good documentation

### ✅ **Adaptive Features**
- Adaptive padding calculation (even if not properly used)
- Auto-detection of genome type (prokaryote vs eukaryote)
- Flexible input handling (with/without GFF)

### ✅ **Phylogenetic Awareness**
- PHYLO_SORT module for ordering genomes
- Tree computation for visualization
- Wavefront parallel processing concept

---

## 6. Recommended Architecture

### Core Principles

1. **Always include the query gene in iterative search database**
2. **Search systematically**: Find synteny → Use synteny to search for query gene
3. **Multi-level search strategy**: Whole gene → Exons → Fragments
4. **Validate orthology**: Use RBH + domain annotation + phylogenetic position
5. **Never fail silently**: Always output regions, even if score is low

---

### Proposed Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                     PHASE 1: INITIALIZATION                     │
└─────────────────────────────────────────────────────────────────┘
   │
   ├─→ Locate GOI in home genome
   ├─→ Extract flanking genes (use GFF or Prodigal)
   ├─→ Build initial database: [FLANKING_GENES] + [GOI] + [GOI_EXONS] + [GOI_FRAGMENTS]
   ├─→ Sort target genomes by phylogenetic distance
   │
┌─────────────────────────────────────────────────────────────────┐
│              PHASE 2: ITERATIVE SEARCH (Per Genome)             │
└─────────────────────────────────────────────────────────────────┘
   │
   FOR EACH genome IN phylo_sorted_order:
   │
   ├─→ STEP 1: FLANKING GENE SEARCH (MMseqs2)
   │   └─→ Search current_db (flanking genes) vs genome
   │   └─→ Identify synteny blocks
   │   └─→ Score blocks, select best (ALWAYS output best, even if score is low)
   │
   ├─→ STEP 2: QUERY GENE SEARCH (Miniprot + MMseqs2)
   │   ├─→ Extract synteny region (generous padding: 100-200kb)
   │   ├─→ Search GOI variants + fragments vs region (MMseqs2 tblastn)
   │   ├─→ Search GOI full/exons vs region (Miniprot splice-aware)
   │   └─→ Merge hits, deduplicate
   │
   ├─→ STEP 3: ANNOTATION & VALIDATION
   │   ├─→ Run Miniprot on all hits for precise exon boundaries
   │   ├─→ RBH validation (vs home genome)
   │   ├─→ Domain validation (InterProScan/Pfam) [OPTIONAL]
   │   ├─→ Flag pseudogenes (stop codons, frameshifts, truncations)
   │   └─→ Assign confidence scores
   │
   ├─→ STEP 4: DATABASE UPDATE
   │   ├─→ Add validated GOI orthologs to current_db
   │   ├─→ Add their exons as separate queries
   │   └─→ (Optional) Add newly discovered flanking genes
   │
   └─→ NEXT genome (using updated database)

┌─────────────────────────────────────────────────────────────────┐
│                  PHASE 3: FINAL ANALYSIS                        │
└─────────────────────────────────────────────────────────────────┘
   │
   ├─→ Build phylogenetic tree of GOI orthologs
   ├─→ Analyze synteny conservation
   ├─→ Detect gene duplications
   ├─→ Identify pseudogenes and losses
   ├─→ Generate interactive synteny plot
   └─→ Generate comprehensive report
```

---

### Key Changes from Current Implementation

| Aspect | Current | Proposed |
|--------|---------|----------|
| **Query gene in iterative DB** | ❌ Not included | ✅ Always included |
| **Search strategy** | Flanking genes only | Flanking + GOI + Exons + Fragments |
| **Miniprot usage** | Re-finds flanking genes | Searches for GOI with splice awareness |
| **Empty regions** | Silent failure | Always output best region |
| **Padding** | 20kb (too small) | 100-200kb adaptive |
| **Exon-level search** | Partial, not integrated | Fully integrated |
| **Fragmentation** | Code exists, not used | Systematic halves/thirds/quarters |
| **Pseudogene detection** | Not implemented | Explicit detection + reporting |
| **Orthology validation** | RBH only | RBH + domains + phylogeny |
| **Database updates** | After wave (parallel) | After each genome (serial for close species) |

---

## 7. Implementation Priority

### 🔴 **CRITICAL (Do First)**

**Week 1-2**: Fix the fatal flaw
1. **Add GOI to initial database** (`extract_flanking_genes.py` or new step)
2. **Include GOI in Miniprot search** (`iterative_search_runner.py:530`)
3. **Remove region filtering** (`cluster_grs.py` always output best region)
4. **Increase padding** to 100kb (configurable)

**Validation**: Run test pipeline and verify GOI is found in at least 1 target genome.

---

### 🟠 **MAJOR (Week 3-4)**

5. **Integrate query fragmentation** into workflow
6. **Implement systematic exon-level search** for GOI
7. **Fix RBH validation logic** (exact matching, domain validation)
8. **Add pseudogene detection** (stop codons, frameshifts, coverage)
9. **Improve phylogenetic ordering** (strict serial for close species)

**Validation**: Run on real dataset with 20+ genomes, verify increased sensitivity.

---

### 🟡 **MODERATE (Week 5-6)**

10. **Centralize GFF parsing** (use `sequence_utils.py`)
11. **Add proper logging** (replace print statements)
12. **Implement checkpointing** for resume capability
13. **Standardize coordinate systems** (document and enforce)
14. **Make sensitivity parameters configurable** consistently

**Validation**: Code review, refactor tests to use centralized functions.

---

### 🔵 **MINOR (Week 7+)**

15. **Add type hints** to all functions
16. **Create unit tests** for individual modules
17. **Add integration tests** with known-good outputs
18. **Move magic numbers** to config files
19. **Clean up unused code** and imports
20. **Improve error messages** and input validation

**Validation**: Achieve 80%+ test coverage, pass all integration tests.

---

## 8. Specific Code Fixes

### Fix 1: Add GOI to Initial Database

**File**: `modules/extract_flanking.nf` or new module

**Add after EXTRACT_FLANKING**:
```groovy
process PREPARE_INITIAL_DB {
    input:
    tuple val(locus_id), path(flanking_faa)
    path query_gene  // Original GOI
    
    output:
    tuple val(locus_id), path("initial_db.faa"), emit: db
    
    script:
    """
    # Combine flanking genes + query gene
    cat $flanking_faa > initial_db.faa
    
    # Add query gene with special ID marking
    cat $query_gene | sed 's/>/>GOI_/' >> initial_db.faa
    
    # Generate query fragments
    fragment_query.py --query $query_gene --output fragments.faa
    cat fragments.faa >> initial_db.faa
    
    # Generate query exons (if multi-exon)
    # [Add exon extraction logic here]
    """
}
```

---

### Fix 2: Include GOI in Miniprot Search

**File**: `bin/iterative_search_runner.py`

**Location**: Lines 530-545

**Replace**:
```python
# OLD CODE:
unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)

# NEW CODE:
unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)

# CRITICAL FIX: Always include queries marked with GOI prefix
all_queries_in_db = set(extract_base_gene_id(clean_id) for _, clean_id, _ in parse_fasta(db_path))
goi_queries = {q for q in all_queries_in_db if q.startswith('GOI_')}
unique_queries.update(goi_queries)  # Force include GOI and its variants
```

---

### Fix 3: Remove Region Filtering

**File**: `bin/cluster_grs.py`

**Location**: Lines 290-300

**Replace**:
```python
# OLD CODE:
if passes_score:
    f_out.write(f"{best['chrom']}\t...")
# else: write nothing!

# NEW CODE:
# ALWAYS output best region, even if score is low
# Add quality flag instead of filtering
quality = "HIGH" if passes_score else "LOW"
f_out.write(f"{best['chrom']}\t{best['start']}\t{best['end']}\t"
            f"{name}\t{best['score']:.2f}\t{region_strand}\t{quality}\n")

# Log warning if score is low
if not passes_score:
    print(f"Warning: Low synteny score ({best['score']:.2f}) for {name}. "
          f"Found {len(best['genes'])} / {total_expected} flanking genes. "
          f"Will still search for query gene in this region.", 
          file=sys.stderr)
```

---

### Fix 4: Increase Padding

**File**: `bin/iterative_search_runner.py`

**Location**: Lines 510-515

**Replace**:
```python
# OLD CODE:
w_start = max(0, best_region['start'] - 20000)
w_end = min(slen, best_region['end'] + 20000)

# NEW CODE:
# Use adaptive padding (already calculated) or large default
padding = calculate_adaptive_padding(hits, best_region, default=150000)  # 150kb default
w_start = max(0, best_region['start'] - padding)
w_end = min(slen, best_region['end'] + padding)

print(f"[{genome_name}] Using {padding/1000:.0f}kb padding for region extraction.", flush=True)
```

---

## 9. Testing Strategy

### Unit Tests Needed

1. **Coordinate conversions** (GFF ↔ BED ↔ Python)
2. **ID extraction** (all patterns in `sequence_utils.py`)
3. **Synteny scoring** (`cluster_grs.py`)
4. **Fragment generation** (`fragment_query.py`)
5. **RBH validation** (`iterative_search_runner.py:batch_rbh_check`)

### Integration Tests Needed

1. **End-to-end with synthetic data**:
   - Create mock genomes with known gene positions
   - Verify correct gene recovery
   
2. **Progressive divergence test**:
   - Home genome + 5 targets at increasing distances
   - Verify iterative search improves over direct search
   
3. **Edge cases**:
   - Gene loss (pseudogene should be reported)
   - Gene duplication (both copies found)
   - Synteny break (still finds gene if regions passed)
   - No flanking genes conserved (falls back to augmented search)

### Regression Tests

Save outputs from current version as baseline:
- Run full pipeline on test datasets
- Save all intermediate files
- After fixes, compare outputs
- Verify improvements (more genes found, correct annotation)

---

## 10. Documentation Updates Needed

### README.md

**Add Limitations Section**:
```markdown
## Limitations

- Requires at least 2-3 conserved flanking genes for synteny detection
- Works best with high-quality genome assemblies
- Performance degrades for genomes with extensive chromosomal rearrangements
- Pseudogene detection is experimental (requires manual curation)
```

### USAGE.md

**Add Troubleshooting Section**:
```markdown
## Troubleshooting

### No genes found in any target genome

Possible causes:
1. Check `results/*/regions/*.regions.bed` - are files empty?
   → Lower `--min_synteny_score` parameter
2. Check `results/*/augmented/*.candidates.bed` - are files empty?
   → Gene may be too divergent; try increasing `--region_padding`
3. Check logs for "No synteny region found"
   → Flanking genes not conserved; consider using fewer flanking genes
```

---

## 11. Final Recommendations

### For Immediate Action

1. **Stop claiming the pipeline works for pseudogene detection** until implemented
2. **Fix the critical GOI search bug** (highest priority)
3. **Add prominent warning in README** about current limitations
4. **Create GitHub issues** for all critical/major bugs identified here

### For Long-Term Success

1. **Hire a bioinformatics software engineer** to refactor the codebase
2. **Publish a preprint** only AFTER critical bugs are fixed and validated on real data
3. **Create comprehensive test suite** before adding new features
4. **Set up continuous integration** (GitHub Actions) to prevent regressions
5. **Consider benchmark comparison** with existing tools (Satsuma, GECKO, etc.)

### Architecture Decision

**Should you refactor or rewrite?**

**Refactor** (recommended):
- Core logic is sound, just misimplemented
- Nextflow workflow structure is good
- Most Python scripts are salvageable
- Estimated time: 4-6 weeks

**Rewrite** (if you have time):
- Opportunity to use modern best practices from the start
- Could switch to Snakemake (better Python integration)
- Estimated time: 3-4 months

---

## 12. Conclusion

SynTerra is a **promising tool with a fatal implementation flaw**. The core scientific idea is excellent, the documentation is thorough, and the modular structure is well-designed. However, the pipeline currently **fails to achieve its primary objective** (iteratively finding the query gene) due to a critical logic error.

### The Good News

All identified issues are **fixable with focused engineering effort**. The codebase is clean enough that refactoring is viable. With 4-6 weeks of dedicated work, this could become a highly valuable tool for the comparative genomics community.

### The Bad News

**The pipeline in its current state should not be used for production research.** Results will be incomplete or incorrect. Extensive testing is needed after fixes are implemented.

### Severity Summary

- **3 Critical Bugs**: Pipeline-breaking, must fix immediately
- **7 Major Issues**: Seriously compromise accuracy/completeness
- **5 Moderate Issues**: Reduce robustness and usability
- **8 Minor Issues**: Code quality and maintainability

### Overall Assessment

**Scientific Value**: ⭐⭐⭐⭐⭐ (5/5) - Excellent concept  
**Current Implementation**: ⭐⭐☆☆☆ (2/5) - Fundamentally flawed  
**Code Quality**: ⭐⭐⭐☆☆ (3/5) - Good structure, but issues  
**Documentation**: ⭐⭐⭐⭐☆ (4/5) - Thorough but overclaims  
**Potential After Fixes**: ⭐⭐⭐⭐⭐ (5/5) - Could be groundbreaking

---

## Appendix: Quick Reference

### Files with Critical Issues

1. `bin/iterative_search_runner.py` - Lines 530-580 (GOI not searched)
2. `bin/cluster_grs.py` - Lines 290-300 (silent failure on empty regions)
3. `bin/augmented_search_runner.py` - Lines 65-75 (handles empty input)
4. `modules/extract_flanking.nf` - Missing GOI inclusion step

### Config Parameters to Adjust Immediately

```groovy
params {
    min_synteny_score = 0.01  // Lower from 0.6 until bug is fixed
    region_padding = 150000    // Increase from 50000
    n_flanking_genes = 5      // Lower from 10 for divergent genomes
}
```

### Priority Fixes (In Order)

1. Add GOI to initial database ← **DO THIS FIRST**
2. Include GOI in Miniprot search
3. Remove region filtering (always output)
4. Increase padding to 100-150kb
5. Add validation tests

---

**End of Critical Analysis**

*This analysis was conducted with the goal of helping you build a robust, scientifically valuable tool. While the critique is harsh, please know that the core idea is excellent and the project is absolutely worth saving. With focused effort on the critical issues, SynTerra can become a reference tool in comparative genomics.*

*If you need help prioritizing fixes or want to discuss specific implementations, I'm here to help.*

---

**Addendum: Response to Your Comments**

Based on your statement:
> "we want to also search the gene of interest iteratively, as this is the hack to find even diverged seqs"

**You are 100% correct**, and this confirms that the critical bug identified (GOI not being searched iteratively) completely undermines your intended design. The fixes outlined in Section 8 will restore your original vision.

Your goals are all achievable with the proposed architecture:
- ✅ Iterative GOI search (fix in progress)
- ✅ Synteny mapping (already works)
- ✅ Exon-intron structure (Miniprot provides this)
- ✅ Gene duplications (detection logic needed)
- ✅ Pseudogenes (detection logic needed)

The fundamental architecture is correct; the implementation just deviated from the plan. Priority: **Fix Critical bugs 1-4**, then add pseudogene/duplication detection logic.
