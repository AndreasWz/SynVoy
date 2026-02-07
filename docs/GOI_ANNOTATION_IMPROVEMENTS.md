# SynTerra GOI Annotation Improvements - February 5, 2026

## Summary

Enhanced the pipeline to properly handle Gene of Interest (GOI) annotation in the iterative search, with special focus on short and divergent proteins like melettin.

## Problem Identified

In the melettin test run:
- ✅ Flanking genes successfully identified syntenic regions (all 4 genomes)
- ❌ Miniprot found **0 hits** in all genomes (too short/divergent)
- ✅ Augmented search found **45 candidates** across genomes
- ❌ Augmented results were **not used** for GOI annotation
- **Result**: No GOI sequences added to expanding database → iterative search stalled

## Root Cause

The pipeline had a **disconnect between phases**:
1. **Iterative Search** (Phase 2): Used miniprot for GOI annotation → failed
2. **Augmented Search** (Phase 3): Found GOI candidates → not integrated
3. No fallback mechanism when miniprot failed
4. Augmented results not fed back into database expansion

## Changes Made

### 1. Enhanced Miniprot Sensitivity (`iterative_search_runner.py`)

**File**: `bin/iterative_search_runner.py`

**Changes**:
```python
# OLD: Basic miniprot command
miniprot -I --gff genome.fna query.faa

# NEW: Sensitive parameters for short/divergent proteins
miniprot \
  -I              # Output introns
  --gff           # GFF3 format  
  -p 0.4          # Lower identity (was 0.75 default)
  --aln           # Alignment details
  -G 100000       # Max intron size
  genome.fna query.faa
```

**Rationale**: Melettin (70aa) and fragments (23-36aa) need relaxed thresholds.

### 2. ORF-Based Fallback Annotation

**File**: `bin/iterative_search_runner.py`

**New Functions**:
- `orf_based_annotation()`: 6-frame translation + ORF finding
- `find_orfs_six_frame()`: Find ORFs in all reading frames
- `simple_protein_similarity()`: Calculate percent identity

**Workflow**:
```python
# Try miniprot first (sensitive mode)
miniprot_hits = run_miniprot(region, queries, sensitive=True)

# If miniprot fails, use ORF fallback
if len(miniprot_hits) == 0:
    print("Miniprot found 0 hits. Trying ORF-based fallback...")
    miniprot_hits = orf_based_annotation(region, queries, genome_seq, offset)
```

**Rationale**: Provides backup when miniprot can't align (too divergent/short).

### 3. Improved GOI Query Selection

**File**: `bin/iterative_search_runner.py`

**Change**:
```python
# OLD: Use all GOI queries (including tiny fragments)
unique_queries.update(goi_queries)

# NEW: Prefer full-length GOI over fragments
full_length_goi = [goi for goi in goi_queries if '|frag_' not in goi]
if full_length_goi:
    unique_queries.update(full_length_goi)  # Use only full-length
else:
    unique_queries.update(goi_queries)      # Fallback to fragments
```

**Rationale**: 
- Fragments (23-36aa) are too short for reliable alignment
- Full-length sequence (70aa) has better chance
- Fragments still available as fallback

### 4. Augmented Search Integration Script

**New File**: `bin/integrate_augmented_goi.py`

**Purpose**: Convert augmented search candidates into annotated GOI sequences

**Features**:
- Simple annotation: 6-frame translation, best ORF selection
- Prodigal annotation: More sophisticated ORF prediction
- Validates minimum length (20aa)
- Adds "GOI_" prefix for tracking
- Outputs FASTA for database expansion

**Usage**:
```bash
integrate_augmented_goi.py \
  --candidate_bed augmented_candidates.bed \
  --genome target_genome.fna \
  --query original_query.faa \
  --output annotated_goi.faa \
  --method simple  # or 'prodigal'
```

### 5. Architecture Documentation

**New File**: `docs/GOI_VS_FLANKING_ARCHITECTURE.md`

**Contents**:
- Clear explanation of GOI vs Flanking genes roles
- Detailed pipeline flow diagrams
- Design rationale for each component
- Melettin test case analysis
- Integration recommendations
- Future enhancement roadmap

## Impact on Melettin Run

### Before Changes
```
Wave 1: Find syntenic region → Run miniprot → 0 hits → No genes added
Wave 2: Find syntenic region → Run miniprot → 0 hits → No genes added
Wave 3: Find syntenic region → Run miniprot → 0 hits → No genes added
Wave 4: Find syntenic region → Run miniprot → 0 hits → No genes added

Augmented Search: Find 45 candidates → Not used for annotation
```

### After Changes
```
Wave 1: Find region → Miniprot (sensitive) → If fail → ORF fallback → Add genes
Wave 2: Find region → Search with new GOIs → Annotate → Add genes
Wave 3: Find region → Search with expanded DB → Annotate → Add genes
Wave 4: Find region → Search with expanded DB → Annotate → Add genes

Augmented Search: Candidates annotated → Integrated into database
```

## Next Steps

### Immediate Testing
1. Re-run melettin test with new parameters
2. Check if ORF fallback finds candidates
3. Verify augmented results can be integrated
4. Compare miniprot vs ORF annotation quality

### Integration TODO
1. **Modify iterative search** to call `integrate_augmented_goi.py` when miniprot fails
2. **Update main workflow** to pass augmented results to next genome iteration
3. **Add logging** to track which method succeeded (miniprot/ORF/augmented)
4. **Create validation** to compare annotation methods

### Code Changes Needed in `main.nf`

```groovy
// After AUGMENTED_SEARCH
INTEGRATE_AUGMENTED_GOI(
    AUGMENTED_SEARCH.out.bed,
    genomes_dir,
    query_gene
)

// Feed results back to ITERATIVE_SEARCH for next genome
// (requires restructuring the wavefront approach)
```

### Testing Plan

1. **Short proteins** (50-100aa): melettin, other venom peptides
2. **Divergent proteins**: immune genes, rapidly evolving families
3. **Standard proteins**: housekeeping genes (positive control)
4. **Multi-exon genes**: test exon boundary detection

## Files Modified

1. **bin/iterative_search_runner.py**
   - Enhanced `run_miniprot()` with sensitive parameters
   - Added `orf_based_annotation()` fallback
   - Added `find_orfs_six_frame()` helper
   - Added `simple_protein_similarity()` helper
   - Improved GOI query filtering (full-length preference)

2. **bin/integrate_augmented_goi.py** (NEW)
   - Converts augmented candidates to GOI annotations
   - Two methods: simple (6-frame) and prodigal
   - Ready to integrate into workflow

3. **docs/GOI_VS_FLANKING_ARCHITECTURE.md** (NEW)
   - Comprehensive architecture documentation
   - Clarifies GOI vs Flanking gene roles
   - Explains design decisions

## Performance Considerations

### Computational Cost
- **Miniprot sensitive mode**: ~same as default (just different params)
- **ORF fallback**: Fast (6-frame translation is cheap)
- **Augmented integration**: Negligible (just annotation step)

### Memory Usage
- No significant increase (ORF finding is in-memory)
- Augmented results already computed, just reusing

### Accuracy Trade-offs
- **Miniprot**: High precision, may miss divergent
- **ORF-based**: Lower precision, catches divergent
- **Best approach**: Try both, validate with RBH

## Validation Strategy

For each annotated GOI:
1. **RBH check**: Reciprocal best hit against home proteome
2. **Length check**: Within 50-150% of query length
3. **Domain check**: (Future) Conserved domains present
4. **Synteny check**: In expected genomic context

## Success Metrics

Pipeline is working correctly when:
- [ ] Miniprot finds hits in most genomes
- [ ] ORF fallback activates only when needed
- [ ] Augmented candidates are validated and integrated
- [ ] Database grows with each genome iteration
- [ ] Final report shows GOI annotations in all target genomes

## Known Limitations

1. **ORF-based annotation**: Less accurate than miniprot (no splice awareness)
2. **Augmented integration**: Not yet fully automatic (requires workflow update)
3. **Fragment queries**: Still used as fallback (may cause noise)
4. **No exon refinement**: ORF method finds continuous ORFs, not exons

## Future Work

1. **Combine methods**: Use miniprot + ORF + augmented consensus
2. **Machine learning**: Train on known orthologs for better prediction
3. **Domain-guided**: Use Pfam/InterPro domains as anchors
4. **Structure prediction**: Validate with AlphaFold-like tools
5. **Interactive mode**: Allow manual curation of uncertain cases

---

**Status**: Code changes complete, testing needed
**Priority**: High (melettin pipeline currently stalled)
**Next Action**: Test on melettin dataset, measure improvement
