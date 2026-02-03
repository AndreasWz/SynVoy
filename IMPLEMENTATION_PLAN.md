# SynTerra Redesign: Exon-Level Search Strategy

## ✅ Implementation Status

### Completed Changes

| Feature | Status | Files Modified |
|---------|--------|----------------|
| Plot color mapping fix | ✅ Done | `bin/plot_synteny.py` |
| SynTerra_Parent attribute parsing | ✅ Done | `bin/plot_synteny.py` |
| Clade-based tree coloring | ✅ Done | `bin/plot_synteny.py` |
| Exon-level CDS extraction | ✅ Done | `bin/extract_flanking_genes.py` |
| Exon ID handling | ✅ Done | `bin/iterative_search_runner.py` |
| Query fragmentation | ✅ Done | `bin/generate_variants.py` |
| Fragment utility module | ✅ Done | `bin/fragment_query.py` (new) |

---

## 🎯 Summary of Requested Changes

Based on the detailed German description, three major architectural changes are requested:

---

## 1. Flanking Gene Search: Exon-Level CDS Approach

### Current Implementation
```
extract_flanking_genes.py:
  1. Parse GFF → Get genes with CDS parts
  2. For each gene: CONCATENATE all CDS → Single protein sequence
  3. Search this whole protein in target genomes
  4. Match by gene ID
```

### Desired Implementation
```
NEW Approach:
  1. Parse GFF → Get genes with individual CDS exons
  2. For each gene: Output EACH EXON as separate CDS sequence
     - ID format: gene_id|exon_1, gene_id|exon_2, etc.
  3. Search each exon CDS in target genomes
  4. In targets: Even if only exons 1,2,5 are found (3,4 missing) → Gene is present
  5. Later: Annotate exon-intron boundaries with Miniprot
```

### Key Changes Required

**A. `extract_flanking_genes.py` Modifications:**
- New output mode: `--exon_mode` (default: True)
- Output individual exon CDS sequences
- ID format: `{gene_id}|exon_{n}` where n is 1-indexed position
- Keep exon coordinates in output for later mapping

**B. `iterative_search_runner.py` Modifications:**
- Parse exon IDs: extract `gene_id` and `exon_n` from hit query
- Group hits by `gene_id` (not just substring matching)
- Gene is "found" if ANY exon hits (not requiring all)
- Consolidate exon hits into gene loci correctly
- Report which exons were found/missing per gene

**C. BED/FAA Output Changes:**
- FAA: Individual exon protein sequences
- BED: Could include exon metadata
- New column or naming convention for exon tracking

---

## 2. Gene of Interest (GOI) Search: Fragmentation Approach

### Current Implementation
```
iterative_search.nf / iterative_search_runner.py:
  1. Search whole query protein once
  2. Find hits
  3. Use Miniprot to align/annotate
```

### Desired Implementation
```
NEW Fragmentation Approach:
  1. Search WHOLE query sequence
  2. Search HALVES (1-50%, 51-100%)
  3. Search THIRDS (1-33%, 34-66%, 67-100%)
  4. Search QUARTERS (1-25%, 26-50%, 51-75%, 76-100%)
  5. Stop when fragments get too small (< min_fragment_size, e.g., 20 AA)
  6. Combine all hits
  7. Use Miniprot to annotate:
     - Cleavage sites (look for V, K, R, etc.)
     - Signal peptides
     - Mature peptide boundaries
     - Frame shifts / Stop codons
```

### Key Changes Required

**A. New script: `fragment_query.py`:**
```python
def generate_fragments(seq, min_size=20):
    fragments = [seq]  # Whole
    length = len(seq)
    
    # Halves
    if length // 2 >= min_size:
        fragments.append(seq[:length//2])
        fragments.append(seq[length//2:])
    
    # Thirds
    if length // 3 >= min_size:
        third = length // 3
        fragments.append(seq[:third])
        fragments.append(seq[third:2*third])
        fragments.append(seq[2*third:])
    
    # Quarters
    if length // 4 >= min_size:
        quarter = length // 4
        fragments.append(seq[:quarter])
        fragments.append(seq[quarter:2*quarter])
        fragments.append(seq[2*quarter:3*quarter])
        fragments.append(seq[3*quarter:])
    
    return fragments  # With position metadata
```

**B. Modify `iterative_search_runner.py`:**
- Call fragment generation
- Search each fragment
- Merge overlapping hits
- Annotate fragment boundaries in output

**C. Post-processing annotation:**
- Use hit boundaries to infer cleavage sites
- Look for characteristic cleavage motifs (V, KR, RR)
- Compare fragment coverage across species

---

## 3. Synteny Plot Coloring

### Current Issue
"aktuell is alles grau in den target genomes" - Everything is gray in target genomes.

### Desired Behavior

**A. Flanking Genes:**
- All flanking genes should have the SAME color (consistent across genomes)
- Color per gene ID, not per genome
- If gene A in home is blue → gene A in all targets is blue

**B. Genes of Interest (GOI):**
- Build phylogenetic tree from all GOI sequences (all genomes)
- Assign colors by clade
- Closely related GOIs get similar colors
- Distantly related GOIs get different colors

### Current Color Logic Analysis
```python
# plot_synteny.py current logic:
# - Parses home_bed for gene names
# - Tries to match target genes to home genes
# - Falls back to gray if no match
```

### Problem
- ID matching fails because:
  - Target GFF has different ID format than home BED
  - Exon-level IDs don't match gene-level IDs
  - No proper ortholog mapping passed to plotter

### Solution

**A. Pass ortholog/homology mapping:**
```python
# homology_tsv format (already exists but not used properly):
# home_gene_id | target_gene_id | identity | target_genome
```

**B. Modify `plot_synteny.py`:**
- Parse homology mapping
- Build ID translation table: target_id → home_gene_id
- Use home_gene_id for color assignment
- Separate color logic for flanking vs GOI

**C. GOI Clade Coloring:**
- Input: Tree from COMPUTE_TREE process
- Parse tree with ete3
- Assign clades by distance cutoff or bootstrap
- Map clade → color gradient

---

## 📁 Files to Modify

| File | Changes |
|------|---------|
| `bin/extract_flanking_genes.py` | Add exon-level output mode |
| `bin/iterative_search_runner.py` | Fragment query, handle exon IDs, consolidate |
| `bin/plot_synteny.py` | Fix color mapping, add clade coloring |
| `bin/generate_variants.py` | May need updates for fragment approach |
| `modules/extract_flanking.nf` | Pass new parameters |
| `modules/iterative_search.nf` | Handle fragment outputs |

---

## 🔧 Implementation Priority

1. **HIGH: Fix Plot Coloring** (Quick win, visible issue)
   - Fix ID matching in plot_synteny.py
   - Use homology TSV for mapping

2. **HIGH: Exon-Level CDS Extraction**
   - Modify extract_flanking_genes.py
   - Update iterative_search_runner.py for exon handling

3. **MEDIUM: Query Fragmentation**
   - New fragment_query.py script
   - Integrate into search pipeline

4. **LOW: Cleavage Site Annotation**
   - Post-processing script
   - Requires biological knowledge input

---

## ⚠️ Potential Issues

1. **MMseqs2 sensitivity for small fragments**
   - May need relaxed E-value for fragments
   - Consider using shorter fragment mode

2. **Exon boundary accuracy**
   - GFF quality varies
   - Phase information crucial for translation

3. **Memory usage**
   - More queries = more search time
   - Consider batching fragments

4. **Tree-based coloring**
   - Requires ete3 (optional dependency)
   - Graceful fallback if tree missing

---

## 🚀 Next Steps

1. Start with **plot coloring fix** (immediate visual improvement)
2. Implement **exon-level extraction** (core functionality)
3. Add **fragmentation search** (enhanced discovery)
4. Polish with **annotation features** (cleavage sites, etc.)

Shall I proceed with implementation?
