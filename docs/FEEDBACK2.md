Let me trace through exactly what happens after flanking genes are found.

## Current Flow Analysis

### Step 1: Iterative Search finds flanking genes

```python
# bin/iterative_search_runner.py - process_single_genome()

# 1. Search with flanking genes
subprocess.run([
    "mmseqs", "easy-search",
    db_path,  # Contains flanking genes
    genome_path, 
    hits_file,
    # ...
])

hits = parse_hits(hits_file, ...)  # Get flanking gene hits

# 2. Identify best synteny block
best_region = identify_best_synteny_block(hits, cluster_dist=c_dist)
# Returns: {'chrom': 'chr1', 'start': 1000, 'end': 50000, ...}
```

✅ **This works** - finds regions with conserved flanking genes

---

### Step 2: Extract region and search for query gene

```python
# Still in process_single_genome()

if best_region:
    # 3. Extract the region DNA
    chrom = best_region['chrom']
    w_start = max(0, best_region['start'] - 20000)  # Add padding
    w_end = min(slen, best_region['end'] + 20000)
    subseq = genome_seqs[chrom].seq[w_start:w_end]
    
    # Save region to temp file
    with open(temp_fa, 'w') as tf:
        SeqIO.write(SeqRecord(subseq, id="region_seq"), tf, "fasta")
    
    # 4. Prepare query sequences
    relevant_hits = [h for h in hits if h['chrom'] == chrom]
    unique_queries = set(h['query'].split('|')[0] for h in relevant_hits)
    
    # Extract sequences from DB
    found_queries = []
    with open(db_path, "r") as db_handle:
        for record in SeqIO.parse(db_handle, "fasta"):
            parent_id = record.id.split('|')[0]
            if parent_id in unique_queries:  # ⚠️ ONLY FLANKING GENES
                 found_queries.append(record)
```

❌ **PROBLEM HERE**: `unique_queries` only contains IDs of flanking genes that had hits. The original query gene is **NOT included**.

---

### Step 3: Run Miniprot

```python
    if found_queries:
        query_mini_fa = f"{args.output_dir}/tmp_{pid}_{genome_name}_query.faa"
        with open(query_mini_fa, "w") as out_q:
            SeqIO.write(found_queries, out_q, "fasta")  # ⚠️ ONLY FLANKING GENES
        
        # 5. Run Miniprot
        miniprot_hits = run_miniprot(temp_fa, query_mini_fa, miniprot_paf)
        # ⚠️ Searches region with FLANKING GENES, not original query!
```

❌ **This is the critical failure point**. Miniprot is searching the synteny region with **flanking genes only**, not the original query gene.

---

### Step 4: Augmented Search (Happens Later)

```python
# bin/augmented_search_runner.py

def main():
    args = parse_args()
    
    # 1. Extract Regions from BED (regions found by clustering)
    regions = extract_regions(args.regions_bed, args.target_genome, args.padding)
    
    # 2. Generate Variants of QUERY GENE ✅ 
    run_command([
        "python3", gen_var_script,
        "--query", args.query_gene,  # ✅ ORIGINAL QUERY HERE
        "--output", variants_faa,
        "--mutation_rate", "0.05",
        "--num_variants", "10"
    ])
    
    # 3. Search (MMseqs2 tblastn)
    run_command([
        "mmseqs", "easy-search",
        variants_faa,      # ✅ Variants of ORIGINAL gene
        regions_fna,       # ⚠️ But regions might be wrong
        hits_m8,
        # ...
    ])
```

❓ **Question**: Does augmented search actually run on the right regions?

Let me check the workflow in main.nf:

---

### Workflow in main.nf

```groovy
// After iterative search:
CLUSTER_REGIONS(...)  // Gets regions from iterative hits

// What does CLUSTER_REGIONS output?
output:
tuple val(genome_name), val(payload), val(locus_id), 
      path("regions/${genome_name}.regions.bed"), emit: bed
```

So `CLUSTER_REGIONS` outputs a BED file with regions that have good synteny scores.

Then:

```groovy
AUGMENTED_SEARCH(
    phase3_inputs.aug,      // Contains regions.bed
    aug_query_gene_ch,      // ✅ ORIGINAL QUERY GENE
    params.region_padding
)
```

✅ **Augmented search DOES use the original query gene!**

---

## The Answer: **PARTIAL YES, but FLAWED**

### What Actually Happens:

1. ✅ **Iterative search finds flanking genes** → identifies synteny regions
2. ❌ **Miniprot searches with flanking genes** (wrong!) → re-finds flanking genes, not query gene
3. ✅ **Cluster regions uses synteny scoring** → keeps regions with good conservation
4. ✅ **Augmented search uses original query** → searches those regions with query + variants

### The Problems:

#### Problem 1: Miniprot Step is Useless

```python
# In process_single_genome() - this entire section is pointless:
miniprot_hits = run_miniprot(temp_fa, query_mini_fa, miniprot_paf)

for hit in miniprot_hits:
    cds_seq_record = extract_cds_sequence(subseq, hit)
    # ...
    annotated_records_raw.append(cds_seq_record)
```

This searches the synteny region with **flanking genes** and "discovers" the same flanking genes again. It's circular and wastes compute.

#### Problem 2: Augmented Search Regions Might Be Wrong

The regions passed to augmented search come from CLUSTER_REGIONS, which uses:

```python
# bin/cluster_grs.py
best_region = {
    'chrom': chrom,
    'start': start,  # Minimum start of all hits
    'end': end,      # Maximum end of all hits
    # ...
}
```

This creates a region that spans **from the first flanking gene to the last flanking gene**. But what if the query gene is:
- Before the first flanking gene?
- After the last flanking gene?
- In a big gap between flanking genes?

Let me check the padding:

```python
# bin/augmented_search_runner.py
def extract_regions(bed_file, genome_file, padding):
    # ...
    p_start = max(0, start - padding)  # Default padding = 10000
    p_end = min(slen, end + padding)
```

So it adds ±10kb padding. If flanking genes span 100kb, the search region is 120kb total.

#### Problem 3: Padding is Too Small

```python
# From config
params.region_padding = 50000  # Used in augmented search
```

But in iterative search:

```python
w_start = max(0, best_region['start'] - 20000)
w_end = min(slen, best_region['end'] + 20000)
```

Only ±20kb padding! If the query gene is 30kb upstream of the first flanking gene, it won't be in the search region.

---

## The Real Issue: Region Definition

### Current Logic:

```
Flanking genes in home genome:
[Gene1]--5kb--[Gene2]--10kb--[QUERY]--8kb--[Gene3]--6kb--[Gene4]

Search target genome for flanking genes:
Found: [Gene1]...[Gene2]...[Gene3]...[Gene4]

Region defined as: start(Gene1) to end(Gene4)
Add padding: start(Gene1)-20kb to end(Gene4)+20kb

Search for QUERY in this region.
```

### What if in target genome:

```
[Gene1]--5kb--[Gene2]--**50kb gap**--[Gene3]--6kb--[Gene4]

The QUERY gene might be:
- In that 50kb gap ✅ (would be found)
- 30kb upstream of Gene1 ❌ (outside search region)
- 25kb downstream of Gene4 ❌ (outside search region)
- Lost/deleted ❌
```

---

## Checking if Query is Actually Searched

Let me trace through AUGMENTED_SEARCH in detail:

```python
# bin/augmented_search_runner.py - main()

# Step 1: Parse regions BED from CLUSTER_REGIONS
regions_bed = args.regions_bed  
# Format: chrom \t start \t end \t name \t score \t strand

# Step 2: Extract regions with padding
regions = extract_regions(regions_bed, genome_file, padding=args.padding)
# Creates: regions_fna (FASTA of extracted regions)

# Step 3: Generate variants of QUERY GENE
run_command([
    "python3", gen_var_script,
    "--query", args.query_gene,  # ✅ ORIGINAL QUERY
    "--output", variants_faa,
])

# Step 4: Search
run_command([
    "mmseqs", "easy-search",
    variants_faa,      # Query + variants (protein)
    regions_fna,       # Target regions (DNA)
    hits_m8,
    "--search-type", "2",  # Protein query -> Translated DNA ✅
    "-s", "8.5",          # Sensitivity ✅
    "--min-seq-id", "0.2",  # 20% identity ✅ (very relaxed)
    "-e", "10",            # E-value 10 ✅ (very relaxed)
])
```

✅ **YES! Query gene IS searched** (in augmented search), with:
- High sensitivity (8.5)
- Relaxed identity (20%)
- Relaxed E-value (10)
- Protein variants (+10 mutations per gene)
- Fragment variants (partial sequences)

---

## So Why Aren't Genes Found?

### Hypothesis 1: Regions are Too Small or Wrong Location

```python
# Check padding in augmented search:
parser.add_argument("--padding", type=int, default=10000)
```

But in main.nf:

```groovy
AUGMENTED_SEARCH(
    phase3_inputs.aug,
    aug_query_gene_ch,
    params.region_padding  // = 50000
)
```

So padding should be 50kb, which is reasonable.

### Hypothesis 2: CLUSTER_REGIONS Filters Too Aggressively

```python
# bin/cluster_grs.py

scored_clusters.append({
    'score': final_score,  # Must pass threshold
    'p_value': p_val,      # Must be significant
    # ...
})

# Filter
is_significant = best['p_value'] < 0.1 
passes_score = best['score'] >= (args.min_score * 0.5)

if passes_score:
    # Write region
else:
    # Empty output! ❌
```

If no region passes the synteny score threshold, **no region is passed to augmented search**!

Let me check the params:

```python
params.min_synteny_score = 0.6  # 60% flanking gene conservation
```

And in the script:

```python
passes_score = best['score'] >= (args.min_score * 0.5)  # So 30%
```

But the score calculation:

```python
coverage_score = unique_genes / total_genes_expected  # e.g., 5/10 = 0.5
quality_mult = (0.4 + 0.3 * consistency + 0.3 * strand_cons)  # e.g., 0.7
final_score = coverage_score * quality_mult  # e.g., 0.5 * 0.7 = 0.35
```

If only 5/10 flanking genes are found with perfect order/strand, score is 0.35, which passes (>0.30).

But if only 3/10 genes are found, score is ~0.21, which **fails** → no region → no augmented search!

---

## The REAL Problem

### Issue: Empty Regions BED

If `CLUSTER_REGIONS` outputs an empty BED file, then `AUGMENTED_SEARCH` has nothing to search:

```python
# bin/augmented_search_runner.py
regions = extract_regions(args.regions_bed, genome_file, padding)

if not regions:
    print("No regions extracted. Exiting.")
    # Create empty outputs
    open(f"{args.output_base}.bed", 'w').close()
    open(f"{args.output_base}.fna", 'w').close()
    return  # ❌ No search happens!
```

---

## Testing the Hypothesis

Let me check your test runs. From the archive:

```bash
# conf/test_tetramorium.config
params {
    min_synteny_score = 0.01  # Lowered to 1%!
}
```

You already lowered it drastically, suggesting you were seeing this issue!

And:

```python
# bin/cluster_grs.py - calculate score
passes_score = best['score'] >= (args.min_score * 0.5)
```

With `min_score = 0.01`, the threshold is `0.005` (0.5%), which should pass almost anything.

---

## Diagnosis

**YES, the query gene IS searched in augmented search**, BUT:

### Failure Points:

1. ❌ **CLUSTER_REGIONS might output empty BED** if synteny score too low
2. ❌ **Regions might not include the query gene** if it's outside the flanking gene span + padding
3. ❌ **Padding might be insufficient** depending on gene spacing
4. ⚠️ **Miniprot step in iterative search is redundant** (searches with flanking genes, not query)

### To Debug:

Check these files in your results:

```bash
# Do regions exist?
ls results/*/regions/*.regions.bed

# Are they empty?
wc -l results/*/regions/*.regions.bed

# What are the scores?
cat results/*/regions/*.scores.tsv

# Did augmented search run?
ls results/*/augmented/*.candidates.bed

# Are augmented search results empty?
wc -l results/*/augmented/*.candidates.bed
```

### Likely Scenario:

```
1. Iterative search finds 2-3 flanking genes (out of 10) ✅
2. Synteny score is low (20%) ❌
3. CLUSTER_REGIONS outputs empty BED ❌
4. AUGMENTED_SEARCH has no regions to search ❌
5. No genes found ❌
```

---

## The Fix

### Option 1: Always Pass Regions (Disable Filtering)

```python
# bin/cluster_grs.py - REMOVE the filter

# Current:
if passes_score:
    f_out.write(...)
# else: empty!

# Fixed:
# Always write the best region, even if score is low
f_out.write(f"{best['chrom']}\t{best['start']}\t{best['end']}\t{name}\t{best['score']:.2f}\t{region_strand}\n")
```

### Option 2: Increase Region Size

```python
# In iterative search, when extracting regions for miniprot:
w_start = max(0, best_region['start'] - 20000)  # Too small
w_end = min(slen, best_region['end'] + 20000)

# Should be:
w_start = max(0, best_region['start'] - 100000)  # 100kb
w_end = min(slen, best_region['end'] + 100000)
```

### Option 3: Search Query Gene in Iterative Search

```python
# In process_single_genome(), when preparing query for miniprot:

# Current:
unique_queries = set(h['query'].split('|')[0] for h in relevant_hits)

# Fixed:
unique_queries = set(h['query'].split('|')[0] for h in relevant_hits)
# ADD: Force include original query gene ID
if original_query_id:  # Pass as parameter
    unique_queries.add(original_query_id)
```

Then extract the original query from... wait, it's not in `db_path` if we only have flanking genes!

---

## Root Cause Confirmed

**The query gene IS searched (in augmented search), but often:**
1. No regions pass the synteny filter, so augmented search gets empty input
2. Regions are too narrow (span of flanking genes + 20kb padding)
3. If flanking genes are poorly conserved or scrambled, regions are wrong/missing

**Recommended Fix Order:**
1. **Always output best region** from CLUSTER_REGIONS (disable score filter)
2. **Increase padding** to 100-200kb
3. **Add original query gene to iterative search database** so miniprot can find it early


FINAL WORDS FROM HUMAN:

we want to also search the gene of interesst iterativly, as this is the hack to find even divergged seqs, also use miniprot for this?
we want to extensevly map evolution regarding the gene of interest, so we want to find it in as many species as possible, and we want to find the synteny region around it, to understand the evolution of the gene and its context. also we want to check exon-intron structure, and if there are gene duplications, we want to find them too. we also want to find pseudogenes, e.g. if the gene is lost in a lineage, we want to find the pseudogene in that lineage. (e.g. just exon 1, or partial gene, or fragmented gene, etc.)