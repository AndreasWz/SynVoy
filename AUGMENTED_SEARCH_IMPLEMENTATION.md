# Implementation Summary: Critical Augmented Search Integration

## What Was Fixed

### The Core Problem
The pipeline had a **fatal architectural flaw**: augmented search (MMseqs2 with variants to find divergent GOI homologs) ran AFTER all iterative searches completed, so its results were never integrated back into the expanding database.

### The Solution
**Moved augmented search INSIDE the iterative genome-by-genome loop** so that:
1. For each genome, after miniprot runs
2. Check which GOI queries miniprot failed to find
3. For missing GOIs, run augmented search (MMseqs2 with variants) on the SAME syntenic region
4. Annotate augmented hits using exon-aware ORF fallback
5. **Add augmented results to database BEFORE moving to next genome**

## Code Changes

### 1. New Function: `run_augmented_search()`
**Location:** `bin/iterative_search_runner.py` (lines ~828-895)

**Purpose:** Run MMseqs2 with sequence variants for divergent GOI discovery

**Logic:**
```python
def run_augmented_search(region_fasta, goi_queries, genome_name, args, unique_id, threads):
    # 1. Generate variants (fragments) for each GOI query
    if FRAGMENT_SUPPORT:
        fragments = generate_fragments(query['seq'], query['id'], min_size=15)
    
    # 2. Run MMseqs2 with relaxed thresholds
    #    - Identity: 60% of normal (e.g., 24% instead of 40%)
    #    - Length: 50% of normal (e.g., 25aa instead of 50aa)
    augmented_hits = parse_hits(aug_hits_file, relaxed_identity, relaxed_length, evalue * 10)
    
    # 3. Return hits for use as ORF guides
    return augmented_hits
```

**Key Features:**
- Uses existing `generate_fragments()` from fragment_query.py
- Relaxed thresholds to catch divergent sequences
- Returns MMseqs2 hits in same format as flanking gene search
- Clean temp file management

### 2. Modified Logic in `process_single_genome()`
**Location:** `bin/iterative_search_runner.py` (lines ~1060-1120)

**Before:**
```python
# Run miniprot
miniprot_hits = run_miniprot(region_fasta, queries, ...)

# If NO hits at all, try ORF fallback
if len(miniprot_hits) == 0:
    miniprot_hits = orf_based_annotation(...)
```

**After:**
```python
# 1. Run miniprot
miniprot_hits = run_miniprot(region_fasta, queries, ...)

# 2. Check which GOI queries were NOT found
miniprot_found_genes = set(extract_base_gene_id(hit['parent_query']) for hit in miniprot_hits)
missing_goi = goi_queries - miniprot_found_genes

# 3. For missing GOIs, run augmented search
if missing_goi:
    augmented_hits_mmseqs = run_augmented_search(region_fasta, missing_goi_queries, ...)
    
    if augmented_hits_mmseqs:
        # Annotate with exon-aware ORF fallback using augmented hits as guides
        augmented_annotations = orf_based_annotation(..., mmseqs_hits=augmented_hits_mmseqs)
        
        # CRITICAL: Merge augmented results with miniprot results
        miniprot_hits.extend(augmented_annotations)

# 4. Final fallback if still no hits
if len(miniprot_hits) == 0:
    miniprot_hits = orf_based_annotation(...)  # Naive mode
```

**Key Changes:**
1. **Selective augmentation:** Only run for GOIs that miniprot missed
2. **Same region:** Augmented search uses the same syntenic region identified by flanking genes
3. **Exon guidance:** Augmented MMseqs2 hits guide ORF extraction (multi-exon aware)
4. **Database integration:** All results (miniprot + augmented) returned as `new_genes`

### 3. Database Expansion
**Location:** `bin/iterative_search_runner.py` main() function (lines ~1445-1470)

**Unchanged but critical:** The existing wave completion logic now receives augmented results:
```python
# After each wave completes
if wave_results:  # wave_results now includes miniprot + augmented genes
    new_genes_fasta = write_fasta(wave_results)
    
    # Append to database
    next_db = concatenate(latest_db, new_genes_fasta)
    
    latest_db = next_db  # Next wave uses expanded DB
```

## Data Flow

### Before Fix (BROKEN):
```
Home Genome
    ↓
[Initial DB: GOI + Flanking]
    ↓
Genome 1:
  - Flanking search → syntenic region ✓
  - Miniprot (GOI) → 0 hits ✗
  - ORF fallback → some hits ✓
  - Add to DB: flanking + ORF hits
    ↓
[DB: GOI + Flanking + some ORFs]
    ↓
Genome 2:
  - Flanking search → syntenic region ✓
  - Miniprot (GOI) → 0 hits ✗
  - ORF fallback → some hits ✓
  - Add to DB: flanking + ORF hits
    ↓
... (NO IMPROVEMENT across genomes)
    ↓
AUGMENTED_SEARCH (runs AFTER all iterations):
  - Finds 45 candidates in Genome 1 ✓
  - Finds 30 candidates in Genome 2 ✓
  - Results → REPORT ONLY ✗
  - NOT added to DB ✗
```

### After Fix (CORRECT):
```
Home Genome
    ↓
[Initial DB: GOI + Flanking]
    ↓
Genome 1:
  - Flanking search → syntenic region ✓
  - Miniprot (GOI) → 0 hits ✗
  - AUGMENTED SEARCH (integrated):
      - Generate variants (fragments) ✓
      - MMseqs2 with relaxed thresholds ✓
      - Find 45 candidates ✓
      - ORF annotate with exon guidance ✓
  - Add to DB: flanking + augmented GOIs ✓
    ↓
[DB: GOI + Flanking + 45 augmented GOIs] ← EXPANDED!
    ↓
Genome 2:
  - Flanking search → syntenic region ✓
  - Miniprot (GOI) → Using 46 GOI queries (1 original + 45 from Genome 1) ✓
  - Miniprot finds 12 hits ✓ (benefits from expanded DB!)
  - AUGMENTED SEARCH for remaining 34 missing GOIs:
      - Find 18 more candidates ✓
      - ORF annotate ✓
  - Add to DB: flanking + miniprot + augmented ✓
    ↓
[DB: GOI + Flanking + 45 + 12 + 18 = 75 GOI variants] ← GROWING!
    ↓
Genome 3:
  - Now searches with 76 GOI queries ✓
  - Better chance of finding divergent homologs ✓
  - Database continues to expand ✓
```

## Testing Plan

### 1. Syntax Validation
```bash
python3 -m py_compile bin/iterative_search_runner.py
```
**Status:** ✅ PASSED

### 2. Melettin Test
```bash
nextflow run main.nf -profile test_melettin -resume
```

**Expected Outcomes:**
- Wave 1 log should show: "Augmented search found N MMseqs2 hits"
- Wave 1 log should show: "Added N GOI candidates from augmented search"
- `db_iter_2.faa` should be larger than `db_iter_1.faa`
- `iter_1_new_genes.faa` should contain augmented candidates
- Genome 2+ should find more hits due to expanded database

### 3. Database Growth Verification
```bash
# Check database sizes across iterations
for db in work/*/*/db_iter_*.faa; do
    echo "$db: $(grep -c '^>' $db) sequences"
done

# Check for augmented candidates
grep -r "augmented" work/*/*/*.log

# Verify new genes files contain GOI variants
grep "GOI_" work/*/*/iter_*_new_genes.faa | wc -l
```

## Cleanup Needed

### AUGMENTED_SEARCH Module (Redundant)
**Files to update:**
- `main.nf`: Remove AUGMENTED_SEARCH process call (lines ~347-353)
- `main.nf`: Remove AUGMENTED_SEARCH from report collection (lines ~476-480)
- `modules/augmented_search.nf`: Mark as deprecated or remove
- `bin/augmented_search_runner.py`: Mark as deprecated or remove

**Rationale:** Augmented search is now integrated into iterative_search_runner.py, so the standalone module is redundant.

**Note:** Keeping it temporarily for backward compatibility, but it should be removed after validating the new integrated approach works correctly.

## Risk Assessment

### Low Risk:
- ✅ Syntax validated (compiles successfully)
- ✅ Uses existing functions (generate_fragments, parse_hits, orf_based_annotation)
- ✅ Fallback logic preserved (if augmented search fails, continue with miniprot results)
- ✅ No changes to database format or file structure

### Medium Risk:
- ⚠️ Performance impact: Each genome now runs additional MMseqs2 search
  - Mitigation: Only runs for missing GOIs, not all queries
  - Mitigation: Uses same threads allocation as main search
- ⚠️ Database size growth: More sequences added per genome
  - Mitigation: Expected behavior - this is the fix!
  - Mitigation: RBH validation still filters low-quality hits

### Testing Required:
- ⏳ Melettin test with full pipeline
- ⏳ Verify database expansion is reasonable (not exponential)
- ⏳ Check memory usage doesn't spike
- ⏳ Validate final report quality improves

## Success Criteria

1. ✅ Code compiles without errors
2. ⏳ Melettin test finds GOI homologs in divergent genomes
3. ⏳ Database grows across iterations (can verify with `grep -c '^>' db_iter_*.faa`)
4. ⏳ Augmented candidates are RBH-validated before adding to database
5. ⏳ Final report shows improved GOI coverage compared to old pipeline

## Next Steps

1. **Run melettin test** with updated pipeline
2. **Monitor logs** for augmented search activity
3. **Validate database expansion** is working correctly
4. **Compare results** to previous run (should find more GOI homologs)
5. **Remove redundant AUGMENTED_SEARCH module** from main.nf if integrated version works
6. **Update documentation** with new architecture

## Date
2026-02-05

## Related Documentation
- [CRITICAL_FIX_AUGMENTED_INTEGRATION.md](CRITICAL_FIX_AUGMENTED_INTEGRATION.md) - Problem statement and solution
- [ORF_FALLBACK_ANALYSIS.md](ORF_FALLBACK_ANALYSIS.md) - Exon-aware ORF annotation
- [GOI_ANNOTATION_IMPROVEMENTS.md](GOI_ANNOTATION_IMPROVEMENTS.md) - Multi-exon handling
