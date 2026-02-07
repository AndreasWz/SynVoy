# CRITICAL FIX: Augmented Search Integration

## The Problem (LIFE-THREATENING BUG)

### What SHOULD Happen:
1. Find GOI in home genome ✓
2. Take flanking genes from GOI region ✓
3. Search flanking genes in first target genome ✓
4. Identify genomic region (GR) in first target genome ✓
5. **Search GOI in the GR** (augmented with variants) ✓
6. **Annotate found GOI hits deeply** (multi-exon, fragments, ORFs) ✓
7. **ADD FOUND GOIs AND FLANKING GENES TO DATABASE** ✗ **THIS WAS BROKEN!**
8. Search flanking genes in next target genome (should use expanded DB) ✗ **WAS NOT HAPPENING!**

### What WAS Happening (BROKEN):

```
Phase 2: ITERATIVE_SEARCH
├── For each genome wave:
│   ├── Search flanking genes (MMseqs2) ✓
│   ├── Identify syntenic region ✓
│   ├── Search GOI with miniprot ✓
│   └── Add miniprot results to DB ✓
│
└── (All waves complete) ← DATABASE IS NOW FROZEN
    ↓
Phase 3: AUGMENTED_SEARCH (TOO LATE!)
├── For each genome:
│   ├── Run MMseqs2 with variants ✓
│   ├── Find divergent GOI candidates ✓
│   └── Send to REPORT only ✗ NEVER ADDED TO DATABASE!
```

**Critical Failure:** Augmented search results were never integrated back into the database, so each genome searched with the same limited query set instead of benefiting from previous discoveries.

### Melettin Example:
- **Genome 1:** Miniprot finds 0 hits (70aa peptide too divergent), augmented search finds 45 candidates
- **Genome 2:** Still searches with only original 70aa query (45 candidates were never added to DB!)
- **Result:** Pipeline fails to propagate GOI discoveries across genomes

## The Fix

### New Architecture (CORRECT):

```
For each genome wave:
├── 1. Search flanking genes (MMseqs2)
├── 2. Identify syntenic region
├── 3. Search GOI with miniprot (sensitive mode)
├── 4. Check which GOI queries were NOT found
├── 5. For missing GOIs:
│   ├── Generate sequence variants (fragments)
│   ├── Run augmented MMseqs2 search on region
│   ├── Use augmented hits as guides for ORF annotation
│   └── Merge augmented results with miniprot results
├── 6. Add ALL found genes to database:
│   ├── Flanking genes from MMseqs2
│   ├── GOI genes from miniprot
│   └── GOI genes from augmented search ← NEW!
└── 7. Next genome uses EXPANDED database ✓
```

### Implementation Details:

**Location:** `bin/iterative_search_runner.py`

**New Function:** `run_augmented_search()`
- Takes: Region FASTA, missing GOI queries, genome name, args
- Generates: Sequence fragments (halves, thirds, quarters)
- Runs: MMseqs2 with relaxed thresholds (60% of normal identity, half normal length)
- Returns: MMseqs2 hits for use as exon guides in ORF annotation

**Modified Logic in `process_single_genome()`:**
```python
# After miniprot runs:
miniprot_found_genes = set(extract_base_gene_id(hit['parent_query']) for hit in miniprot_hits)
missing_goi = goi_queries - miniprot_found_genes

if missing_goi:
    # Run augmented search on the SAME region
    augmented_hits_mmseqs = run_augmented_search(region_fasta, missing_goi_queries, ...)
    
    # Annotate with exon-aware ORF fallback
    augmented_annotations = orf_based_annotation(..., mmseqs_hits=augmented_hits_mmseqs)
    
    # CRITICAL: Merge with miniprot results
    miniprot_hits.extend(augmented_annotations)
```

**Result Flow:**
- All annotated genes (miniprot + augmented) → `new_genes`
- `new_genes` → added to `wave_results`
- `wave_results` → appended to database after wave completes
- Next wave searches with **expanded database**

## Impact

### Before Fix:
- Database grew only from miniprot hits
- Divergent GOI sequences were found but discarded
- Each genome searched independently with same queries
- Melettin test: 0 propagation across genomes

### After Fix:
- Database grows from miniprot + augmented search
- Divergent GOI sequences are annotated and propagated
- Each genome benefits from previous discoveries
- Melettin test: Should find homologs in increasingly divergent genomes

## Testing

To validate the fix works:

```bash
# Run melettin test with updated pipeline
nextflow run main.nf -profile test_melettin -resume

# Check database expansion:
ls -lh work/*/*/iter_*_new_genes.faa
ls -lh work/*/*/db_iter_*.faa

# Verify augmented candidates are added:
grep "augmented" work/*/*/iter_*.log
```

Expected output:
- `iter_1_new_genes.faa` should contain miniprot + augmented results
- `db_iter_2.faa` should be larger than `db_iter_1.faa`
- Logs should show "Added N GOI candidates from augmented search"

## Files Modified

1. `bin/iterative_search_runner.py`:
   - Added `run_augmented_search()` function
   - Modified `process_single_genome()` to call augmented search for missing GOIs
   - Integrated augmented results into database expansion

## Related Documentation

- [ORF_FALLBACK_ANALYSIS.md](ORF_FALLBACK_ANALYSIS.md) - Exon-aware ORF annotation
- [GOI_ANNOTATION_IMPROVEMENTS.md](GOI_ANNOTATION_IMPROVEMENTS.md) - Multi-exon handling
- [GOI_VS_FLANKING_ARCHITECTURE.md](GOI_VS_FLANKING_ARCHITECTURE.md) - Query selection logic

## Status

- ✅ Code implemented
- ✅ Syntax validated (py_compile)
- ⏳ Needs testing with melettin dataset
- ⏳ Needs validation of database expansion

## Date

2026-02-05
