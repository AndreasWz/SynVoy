# SynVoy Feature Analysis: From Sequence-Synteny to Structure- & Semantics-Synteny

**Date:** 2026-04-10  
**Scope:** Technical feasibility, integration strategy, critical evaluation, and cost analysis for the four feature proposals in `feature_ideas.txt`.

---

## Table of Contents

1. [Feature 1: Structural Search (Foldseek) & PLM Embeddings (ESM-2/ProTrek)](#1-the-ai-leap-structural--semantic-search)
   - [1A: Foldseek Integration](#1a-foldseek-integration---on-the-fly-structural-search)
   - [1B: PLM-Based Embedding Search](#1b-plm-based-embedding-search-esm-2--protrek)
   - [1C: Hosting & Compute Cost Estimates](#1c-hosting--compute-cost-estimates)
2. [Feature 2: Pangenomics — Graph Genomes](#2-pangenomics--graph-genomes)
3. [Feature 3: Non-Coding Elements & Operons (Deferred)](#3-non-coding-elements--operons-deferred)
4. [Feature 4: nf-core Standardization & Community](#4-nf-core-standardization--community)
5. [Priority Ranking & Roadmap](#5-priority-ranking--roadmap)
6. [Sources](#sources)

---

## 1. The AI Leap: Structural & Semantic Search

### The Core Thesis

SynVoy's unique value proposition is that it reduces the search space dramatically through synteny — it narrows down "where to look" using conserved gene neighborhoods. Currently, the "how to look" step within that narrowed block relies on sequence-based methods (MMseqs2 → tblastn → Smith-Waterman). The thesis is: if we replace or augment the "how to look" step with structure- or embedding-based methods, we can find orthologs that have diverged beyond sequence recognition (<10-15% identity), because 3D structure and biological function are conserved far longer than amino acid sequence.

**This thesis is fundamentally sound.** It's one of the most exciting directions for SynVoy. But the devil is in the details.

---

### 1A: Foldseek Integration — On-the-Fly Structural Search

#### What Foldseek Does

Foldseek uses a "3Di" structural alphabet that encodes the 3D geometric relationship between each residue and its nearest spatial neighbor into a 20-letter alphabet. This converts a 3D structure comparison into something like a sequence alignment — enabling MMseqs2-like speed (4-5 orders of magnitude faster than TM-align/Dali) while retaining 86-133% of their sensitivity.

#### How It Would Integrate Into SynVoy

The integration point is inside `run_augmented_search()` in `iterative_search_runner.py` (lines 1571-1753). Currently:

```
Region FASTA → MMseqs2 (fast) → tblastn (sensitive) → Smith-Waterman (very sensitive) → combined hits
```

With Foldseek, this would become:

```
Region FASTA → MMseqs2 → tblastn → Smith-Waterman → [Foldseek structural] → combined hits
```

**Concrete steps:**

1. **Predict ORFs in the syntenic region** — Run Prodigal or miniprot on the extracted region to get candidate proteins. SynVoy already does this partially via miniprot (line 2307+), but it happens *after* hits are found. For Foldseek, we'd need ORF predictions *before* or in parallel with the sequence search.

2. **Fold the candidate proteins** — Use ESMFold to predict structures for the region's ORFs. This is the expensive step.

3. **Fold the GOI query** (once, cached) — Predict or retrieve the reference structure for the gene of interest.

4. **Run Foldseek** — Compare the GOI structure against all region ORF structures. Foldseek itself is extremely fast (milliseconds per comparison); the bottleneck is step 2.

5. **Convert hits to the existing m8 format** — Foldseek outputs alignment scores (TM-score, 3Di identity) that would need mapping to a format compatible with the existing hit pipeline.

6. **Feed into classification** — Add a new `evidence_type="structural_match"` to `_classify_goi_evidence()` (line 596+).

#### Critical Evaluation

**The good:**
- The synteny-first approach makes this *uniquely* feasible. SynVoy doesn't need to fold entire proteomes — only the ~5-50 ORFs within a single syntenic block per genome. This is a massive advantage over a naive "fold everything" approach.
- For a typical syntenic region of ~200kb, Prodigal would predict maybe 20-80 ORFs (depending on gene density). At ~14 seconds per ESMFold prediction on a V100 GPU, that's 5-19 minutes per block. For a pipeline that already runs hours, this is acceptable.
- Structure diverges ~10x slower than sequence. For deeply divergent orthologs (e.g., insect vs. mammal homologs of fast-evolving genes), this would catch what tblastn misses.

**The problematic:**
- **ESMFold accuracy drops for remote homologs.** Benchmarks show ESMFold achieves median TM-score 0.95 and pLDDT 87.4 overall — but for proteins with few homologs (exactly the case we care about!), accuracy drops. AlphaFold2 does better because it uses MSAs, but is much slower (~40 min/protein including MSA search). ESMFold's speed advantage (14s/protein) comes with a quality trade-off precisely where we need it most. However, even imperfect structures are useful for Foldseek — a rough fold is often enough to detect the right topology.
- **GPU dependency.** ESMFold requires a GPU with ≥8GB VRAM for proteins up to ~400 residues, ≥16GB for ~700 aa, ≥40GB for ~2000 aa. This means SynVoy would need a GPU profile or cloud offloading, which breaks the current "run on a laptop with conda" simplicity.
- **The ORF prediction chicken-and-egg problem.** Foldseek compares protein structures. But the syntenic region is genomic DNA. We need to first translate ORFs before we can fold them. For well-annotated genomes (with GFF), this is trivial. For unannotated genomes (SynVoy's bread and butter), we rely on ab initio prediction (Prodigal), which itself can miss highly divergent genes — exactly the ones we want Foldseek to find. However: Prodigal is good at finding ORFs even when it can't predict function; the structure step adds the functional identification layer on top.

**Verdict: HIGH value, MEDIUM-HIGH complexity. The synteny-guided approach makes this far more feasible than it would be for a general tool. The ~20-80 ORFs per block keep compute manageable. The main challenge is the GPU requirement, which changes SynVoy's deployment model.**

#### Recommended Implementation Strategy

**Phase 1 — AlphaFold DB lookup (no GPU needed):**
Before predicting structures on-the-fly, check if the GOI or close homologs already have structures in the AlphaFold Protein Structure Database (>200M predicted structures). For many well-studied proteins, the reference structure is already available. This is free and instant.

**Phase 2 — ESMFold on-the-fly (GPU mode):**
Add an optional `--enable_structural_search` flag. When enabled:
- Require GPU availability (auto-detect CUDA)
- Run Prodigal/miniprot to predict ORFs in each syntenic block
- Fold each ORF with ESMFold
- Compare against GOI structure via Foldseek
- Merge structural hits with sequence hits before classification

**Phase 3 — ColabFold/batch service:**
For users without local GPUs, support offloading structure predictions to a remote API or pre-compute structures during genome preparation.

---

### 1B: PLM-Based Embedding Search (ESM-2 / ProTrek)

#### What This Means

Instead of comparing amino acid sequences with substitution matrices (BLOSUM62), extract high-dimensional "meaning" vectors (embeddings) from protein language models. Proteins with similar function/structure will be close in embedding space even with <10% sequence identity.

Two concrete options:

**ESM-2 (Meta/FAIR):** Generates per-residue and per-protein embeddings. The 650M-parameter model produces 1280-dimensional vectors. Cosine similarity between embeddings correlates with structural and functional similarity.

**ProTrek (Westlake University):** A trimodal model that jointly embeds protein sequence, structure, AND natural language function descriptions. ProTrek surpasses Foldseek and MMseqs2 in speed and accuracy for identifying functionally related proteins. It has been validated experimentally and supports a database of >5 billion pre-computed embeddings.

#### How It Would Integrate Into SynVoy

Two integration approaches, each with different trade-offs:

**Approach A — Embedding-based hit rescue (post-search filter):**

After the existing search pipeline finds candidate ORFs in a syntenic block, compute ESM-2 embeddings for all candidates AND the GOI. Rank by cosine similarity. This could "rescue" true orthologs that scored too low on sequence identity but are structurally/functionally correct.

- Integration point: after `run_augmented_search()` returns hits, before classification
- Minimal pipeline disruption
- Only embeds proteins that were already found (doesn't find new ones)

**Approach B — Embedding-based primary search (replace/augment sequence search):**

Predict all ORFs in the syntenic block (Prodigal), embed them with ESM-2, embed the GOI, then find nearest neighbors in embedding space. This is a fundamentally different search modality.

- Integration point: inside `run_augmented_search()`, as a 4th search method alongside MMseqs2/tblastn/SW
- Requires ORF prediction first (same as Foldseek)
- Can find completely novel hits invisible to sequence methods

**Approach C — ProTrek function-guided search:**

ProTrek's unique feature: you can search by natural language description. E.g., "antimicrobial peptide with amphipathic helix" would find functionally similar proteins regardless of sequence. This could be used as a verification step — after finding a candidate in a syntenic block, ask "does this protein's embedding match the functional description of my GOI?"

- Extremely powerful for annotation and confidence scoring
- Requires text description of GOI function (which is usually available from UniProt)

#### Critical Evaluation

**The good:**
- **ESM-2 is fast on CPU.** Unlike Foldseek (which needs prior structure prediction), ESM-2 embedding inference runs on CPU in ~0.5-2 seconds per protein (for the 650M model). No GPU required. This preserves SynVoy's "run anywhere" philosophy.
- **ProTrek is a game-changer for confidence classification.** Instead of relying purely on identity/coverage thresholds (the current `_classify_goi_evidence()`), an embedding-based similarity score adds a completely orthogonal signal. A protein at 12% sequence identity but 0.95 embedding similarity to the GOI is almost certainly a true ortholog.
- **Pre-computed databases exist.** ProTrek has >5 billion pre-computed embeddings. For standard UniProt proteins, you don't even need to compute — just look up.
- **Complementary to Foldseek.** Embeddings capture function/fold-family membership; Foldseek captures detailed structural geometry. Together they're very powerful.

**The problematic:**
- **ESM-2 model size.** The 650M model is ~2.5 GB. The 3B model is ~12 GB. These need to be downloaded and loaded into memory. For a bioinformatics pipeline that users install via conda, bundling a 2.5 GB neural network is a significant addition.
- **The ORF prediction bottleneck (again).** For Approach B, we need all possible ORFs in the region before we can embed them. Same chicken-and-egg as Foldseek.
- **Embedding similarity ≠ orthology.** Embeddings capture fold-level similarity. Two proteins from the same superfamily will have high embedding similarity even if they're paralogs, not orthologs. Synteny helps here (if it's in the right genomic context, it's more likely the ortholog), but it doesn't eliminate the problem. This is especially relevant for gene families like LY6/3FTx where SynVoy already struggles with GOI vs. paralog discrimination.
- **Validation gap.** There are no published benchmarks of "PLM embedding search within syntenic blocks" — SynVoy would be pioneering this combination. That's exciting for a paper, but risky for reliability.

**Verdict: VERY HIGH value, MEDIUM complexity. ESM-2 embeddings (Approach A — post-search rescue/re-ranking) can be integrated with minimal disruption and no GPU requirement. ProTrek adds a powerful function-aware dimension. This is arguably the highest-impact, lowest-barrier feature of all four proposals.**

#### Recommended Implementation Strategy

**Phase 1 — ESM-2 embedding re-ranking (CPU, no new dependencies beyond torch):**
- After augmented search finds GOI candidates in a block, compute ESM-2 embeddings for candidates + GOI
- Add `embedding_similarity` as a new signal to the classification system
- Use it to upgrade LOW→MEDIUM confidence when embedding similarity is high despite low sequence identity
- Flag for users: `--enable_embedding_rerank` (default off until validated)

**Phase 2 — Embedding-based primary search:**
- Predict all ORFs with Prodigal in each block
- Embed everything with ESM-2
- k-NN search against GOI embedding
- Merge with sequence-based hits

**Phase 3 — ProTrek integration:**
- Accept GOI function description as input (or auto-fetch from UniProt)
- Use ProTrek's trimodal space for candidate verification
- Add natural language annotations to output reports

---

### 1C: Hosting & Compute Cost Estimates

#### The User's Question: "What would it cost to host SynVoy on a strong server with compute to use Foldseek on the fly?"

This depends heavily on usage patterns. Let me break it down:

#### Scenario 1: Single-User Research Server

A dedicated GPU workstation for one lab group running SynVoy.

| Component | Specification | Cost |
|---|---|---|
| GPU | NVIDIA RTX 4090 (24GB VRAM) | ~$1,600 (one-time) |
| CPU | AMD Ryzen 9 / Intel i9 (16 cores) | ~$500 |
| RAM | 128 GB DDR5 | ~$400 |
| Storage | 2 TB NVMe SSD | ~$150 |
| **Total hardware** | | **~$2,650 one-time** |
| Power + cooling (yearly) | ~300W average | ~$400/year |

- ESMFold on RTX 4090: ~5-8 seconds per protein (faster than V100 benchmarks)
- Throughput: ~50-100 structures per syntenic block = ~8-13 minutes per block
- For a typical run (15 genomes × 3 blocks × 50 ORFs): ~2,250 predictions = ~5 hours of folding
- This is in the same ballpark as the existing runtime — acceptable.

#### Scenario 2: Cloud GPU (On-Demand)

Pay-per-use on AWS, GCP, or budget providers.

| Instance Type | GPU | Price/Hour (2026) | Per SynVoy Run (~5h folding) |
|---|---|---|---|
| AWS p3.2xlarge | 1× V100 (16GB) | ~$1.50/h (spot) | ~$7.50 |
| AWS p4d.xlarge | 1× A100 (40GB) | ~$1.00/h (spot) | ~$5.00 |
| RunPod A100 | 1× A100 (80GB) | ~$0.80/h | ~$4.00 |
| Thunder Compute H100 | 1× H100 (80GB) | ~$1.50/h | ~$7.50 |
| Google Cloud A100 | 1× A100 (40GB) | ~$1.20/h (spot) | ~$6.00 |

**Key insight:** GPU prices have dropped dramatically. A100s are now sub-$1/hour on the spot market. The cost per SynVoy run is in the $4-8 range — negligible for a research setting.

#### Scenario 3: Hosted SynVoy Service (Multi-User)

If SynVoy were hosted as a web service for the community:

| Monthly Usage | Est. GPU Hours | Cloud Cost/Month | Notes |
|---|---|---|---|
| 10 runs/month (small lab) | ~50 h | ~$50-80 | Spot instances |
| 100 runs/month (institute) | ~500 h | ~$400-600 | Reserved instances better |
| 1000 runs/month (public service) | ~5000 h | ~$3,000-5,000 | Needs autoscaling, queue management |

**For a single lab:** $50-80/month is extremely affordable. A university cloud allocation easily covers this.

**For a public service:** $3,000-5,000/month is substantial but comparable to other bioinformatics SaaS tools. Could be grant-funded or sustained through institutional subscriptions.

#### Scenario 4: ESM-2 Embeddings Only (No GPU Needed!)

If you implement only PLM embeddings (Feature 1B, Phase 1) without Foldseek:

| Component | Cost |
|---|---|
| Additional compute | **$0** — runs on existing CPU |
| Storage for ESM-2 model | ~2.5 GB disk (650M model) |
| Runtime overhead per run | ~5-15 minutes additional |

**This is the most cost-effective path.** ESM-2 embedding re-ranking adds AI-powered remote homolog detection to SynVoy with zero infrastructure cost increase.

---

## 2. Pangenomics — Graph Genomes

### The Proposal

Replace linear reference genomes (.fna) with graph genomes (e.g., from minigraph-cactus) that represent all known structural variants (inversions, deletions, insertions) of a population.

### How SynVoy Currently Handles Genomes

SynVoy treats each genome as a collection of linear contigs (chromosomes/scaffolds). The core data structure is:

```python
genome_seqs = {contig_name: sequence_string, ...}  # dict of linear sequences
```

Key operations that assume linearity:
- **Region extraction** (`process_region_block()`, line 2086): `subseq = genome_seqs[chrom][w_start:w_end]` — simple string slicing
- **Hit coordinates**: All hits use `(chrom, start, end, strand)` on a linear axis
- **Synteny blocks**: Defined as contiguous intervals on a single chromosome
- **Adaptive padding**: Assumes a continuous genomic region around the block
- **Cluster distance**: Estimated from inter-gene spacing on a linear genome

### What Would Need to Change

A graph genome represents sequences as nodes connected by edges, where alternative paths through the graph represent structural variants. The fundamental data model changes from:

```
Linear: chr1 ──────────────────────────────────────────────────>
                    [Gene A]  [Gene B]  [GOI]  [Gene C]

Graph:                        ╭──[Gene B']──╮
        chr1 ──[Gene A]──────┤              ├──[GOI]──[Gene C]──>
                              ╰──[Gene B]───╯
                                  ↕ (inversion)
```

**Changes required:**

1. **Genome loading** — Replace FASTA dict with GFA (Graphical Fragment Assembly) format parser. Store a graph structure (nodes + edges + paths) instead of flat sequences.

2. **Region extraction** — Cannot use string slicing. Need to enumerate all paths through the relevant subgraph and extract each as a candidate region. A single syntenic block might correspond to multiple haplotype-specific sequences.

3. **Hit coordinate system** — Replace `(chrom, start, end)` with path-aware coordinates `(path, node_offset, start, end)`. The existing GFF output format would need modification or extension.

4. **Synteny scoring** — Currently counts shared flanking genes. With a graph, the same gene might appear on different paths. Need to define what "syntenic" means when the block structure itself varies between haplotypes.

5. **MMseqs2/tblastn** — These tools work on linear FASTA. Would need to extract all paths as separate sequences before running search tools. This is a linearization step that partially defeats the purpose of using graphs.

6. **Output and visualization** — All downstream tools (plotting, GFF, tree building) assume linear coordinates.

### Critical Evaluation

**The good:**
- The scientific case is strong. Pangenomes genuinely capture variation that single references miss. For species with high structural variation (plants, humans), a pangenome reference finds genes that are simply absent from any single haplotype.
- Minigraph-cactus is mature, well-supported software (Nature Biotechnology 2023). The output formats (GFA, VCF) are standardized.
- SynVoy's iterative cross-species search is conceptually related to pangenomics — both acknowledge that a single reference is insufficient.

**The problematic:**
- **SynVoy operates across species, not within species.** Pangenomes represent intra-species variation. SynVoy searches for orthologs across species separated by millions of years of evolution. The structural variation captured by pangenomes (a few % sequence difference) is dwarfed by the inter-species divergence SynVoy already handles. The additional information from a pangenome is marginal for SynVoy's core use case.
- **Massive engineering effort.** Every core data structure and algorithm in the 4,400-line search engine assumes linear coordinates. Refactoring to graph-aware processing is essentially a rewrite of the coordinate system, region extraction, and hit mapping — conservatively 2-3 months of full-time work.
- **Tool ecosystem not ready.** MMseqs2, tblastn, miniprot, Prodigal — none of SynVoy's search tools natively support graph genomes. Each syntenic block would need to be "linearized" before running any search tool, which means extracting all paths and searching them separately. This multiplies compute by the number of haplotypes.
- **User base mismatch.** Pangenome data is currently available only for a handful of species (human, rice, sorghum, tomato). The vast majority of SynVoy's target users work with single-reference assemblies. Building pangenome support for a niche within a niche is hard to justify.
- **The real gain is Presence/Absence Variation (PAV).** If a gene is deleted in the reference haplotype but present in another, a pangenome would let SynVoy find it. But SynVoy already mitigates this by searching multiple assemblies per species (easy mode's `max_genomes` parameter). Searching 3-5 different assemblies of the same species catches most PAV without any graph infrastructure.

**Verdict: LOW-MEDIUM value for current use case, VERY HIGH complexity. The engineering cost is enormous, the tool ecosystem isn't ready, and SynVoy's multi-assembly approach already captures most of the benefit. This is a "maybe in 2-3 years when the ecosystem matures" feature, not a near-term priority.**

### Recommended Alternative: Multi-Assembly Mode

Instead of graph genomes, enhance SynVoy's existing multi-genome approach:

1. **Per-species assembly pooling:** In easy mode, automatically download 2-3 assemblies per species and merge results. This captures PAV without any new data model.
2. **Consensus synteny blocks:** When the same gene appears in multiple assemblies, combine evidence across assemblies for higher-confidence block identification.
3. **Variant-aware annotation:** Flag when a GOI is found in assembly A but absent in assembly B of the same species, indicating possible PAV.

This achieves 80% of the pangenome benefit at 10% of the engineering cost.

---

## 3. Non-Coding Elements & Operons (Deferred)

*The user indicated this is not a priority, so this section provides a brief assessment for future reference.*

### Conserved Non-coding Elements (CNEs) / Enhancers / lncRNAs

**The idea:** Use synteny to find regulatory elements that have no sequence conservation across distant species.

**Assessment:** Scientifically fascinating, but fundamentally different from SynVoy's current architecture:
- SynVoy's entire search pipeline (MMseqs2, tblastn, miniprot, SW) works on protein sequences. Non-coding elements require DNA-level tools (e.g., BLAST for short conserved motifs, or epigenomic data integration).
- The classification system assumes protein-coding genes (exon counts, query coverage, amino acid identity).
- Would essentially require a parallel pipeline for non-coding search.
- Better as a separate companion tool (`SynVoy-NC`) than a bolt-on to the existing pipeline.

### Biosynthetic Gene Clusters (BGCs) in Metagenomes

**The idea:** Search for operons/gene clusters rather than single genes.

**Assessment:** More compatible with SynVoy's architecture (it already searches for gene neighborhoods), but:
- Metagenomes are fragmented, short-contig assemblies. Synteny blocks spanning multiple genes are harder to identify on 5-20kb contigs.
- Existing tools like antiSMASH are mature and well-established for BGC detection.
- SynVoy's value-add would be cross-metagenome synteny, which is an interesting niche but a different product.

---

## 4. nf-core Standardization & Community

### Current State

SynVoy is **not nf-core compliant.** Key gaps:

| Requirement | SynVoy Status | Effort |
|---|---|---|
| nf-core template structure | Custom layout | HIGH — restructure entire pipeline |
| `nextflow_schema.json` | Missing — uses `nextflow.config` only | MEDIUM — generate from existing params |
| Module format (meta maps, emit blocks) | Custom `.nf` modules | HIGH — rewrite all 21 modules |
| Container requirements (Docker + Singularity) | Docker only (partial) | MEDIUM |
| CI tests with nf-test | No tests | HIGH |
| MIT license | Not specified | LOW |
| Community ownership | Single developer | Governance change |
| Stub processes for dry-run | None | MEDIUM |
| MultiQC integration | None | LOW-MEDIUM |

### Critical Evaluation

**The good:**
- nf-core adoption would give SynVoy instant visibility among thousands of bioinformatics users at hospitals and institutes that standardize on nf-core.
- The quality standards (CI, tests, schema validation) would improve SynVoy's robustness significantly.
- MultiQC integration would make the output much more accessible.

**The problematic:**
- **The refactoring is massive.** nf-core requires a very specific project structure, module format, and workflow pattern. SynVoy's custom architecture (especially the Python-heavy `iterative_search_runner.py` which is essentially a pipeline-within-a-pipeline) doesn't fit the nf-core module paradigm where each process does one thing.
- **Community ownership requirement.** nf-core pipelines must be maintained by the community, not a single developer. For a bachelor thesis project, this governance change might be premature.
- **nf-core tooling overhead.** The template adds significant boilerplate (params validation, CI configs, container builds) that increases maintenance burden.
- **SynVoy is unique.** There's no existing nf-core pipeline in the synteny/ortholog space to build on. This would be a from-scratch nf-core pipeline, not an adaptation.

**Verdict: MEDIUM value, VERY HIGH effort. The visibility benefits are real, but the refactoring cost is enormous and the community governance requirement may not fit SynVoy's current development stage.**

### Recommended Alternative: Incremental Standards Adoption

Instead of full nf-core compliance, adopt the most valuable standards piece by piece:

1. **`nextflow_schema.json`** — Generate parameter schema for built-in validation and GUI support (Nextflow Tower, nf-core launch). This alone gives users a much better experience. Effort: 1-2 days.

2. **MultiQC module** — Write a custom MultiQC plugin for SynVoy's report JSON. This makes results visible in institutional MultiQC dashboards. Effort: 2-3 days.

3. **Singularity/Apptainer support** — Add a `singularity` profile. Many HPC clusters require this. Effort: 1 day.

4. **CI tests** — Add basic integration tests with nf-test or simple pytest. Effort: 3-5 days.

5. **Full nf-core port** — Defer until SynVoy has a stable user base and community contributors.

---

## 5. Priority Ranking & Roadmap

Based on impact, feasibility, and alignment with SynVoy's mission:

### Tier 1: High Impact, Achievable Now

| Feature | Why | Effort | Prerequisites |
|---|---|---|---|
| **ESM-2 embedding re-ranking** (1B, Phase 1) | Adds AI-powered remote homolog detection with zero infrastructure cost. CPU-only. Directly improves the weakest part of the pipeline (detecting divergent orthologs). | 1-2 weeks | PyTorch + ESM-2 model download |
| **`nextflow_schema.json`** (4, Step 1) | Instant UX improvement for all users. Enables Nextflow Tower integration. | 1-2 days | None |

### Tier 2: High Impact, Moderate Effort

| Feature | Why | Effort | Prerequisites |
|---|---|---|---|
| **Foldseek + ESMFold structural search** (1A) | Catches orthologs invisible to sequence methods. Synteny-first approach keeps compute tractable. | 3-4 weeks | GPU profile, ESMFold + Foldseek dependencies |
| **ProTrek trimodal verification** (1B, Phase 3) | Adds function-aware confidence scoring. Could revolutionize classification accuracy. | 2-3 weeks | ProTrek model, API integration |
| **MultiQC integration** (4, Step 2) | Institutional adoption enabler. | 2-3 days | MultiQC plugin structure |

### Tier 3: Strategic, Long-Term

| Feature | Why | Effort | Prerequisites |
|---|---|---|---|
| **Multi-assembly PAV mode** (2, alternative) | Captures pangenome benefits without graph infrastructure. | 1-2 weeks | Refactor easy-mode genome selection |
| **Full nf-core port** (4, Step 5) | Community adoption at scale. | 2-3 months | Stable API, community interest |
| **Graph genome support** (2) | Future-proofing for when ecosystem matures. | 3-4 months | GFA ecosystem maturation |

### Suggested Implementation Order

```
Now ──── nextflow_schema.json (2 days)
  │
  ├──── ESM-2 embedding re-ranking (2 weeks)
  │       └── Validate on ground truth (melittin, known divergent cases)
  │
  ├──── MultiQC plugin (3 days)
  │
  ├──── Foldseek/ESMFold structural search (4 weeks)
  │       └── GPU profile, optional flag
  │       └── Validate: does it find what sequence search misses?
  │
  ├──── ProTrek integration (3 weeks)
  │       └── Function-aware classification
  │
  ├──── Multi-assembly PAV mode (2 weeks)
  │
  └──── nf-core port (when community exists)
         Graph genomes (when ecosystem ready)
```

---

## Summary: What Makes Sense and What Doesn't

| Idea | Makes Sense? | Why / Why Not |
|---|---|---|
| Foldseek on-the-fly | **Yes, strongly** | SynVoy's synteny-first approach uniquely limits the search space, making on-the-fly folding tractable. Few other tools can do this. |
| ESM-2 embeddings | **Yes, strongly** | CPU-only, low barrier, orthogonal signal to sequence search. Best value-for-effort of all proposals. |
| ProTrek trimodal | **Yes** | Function-aware search is genuinely novel in this context. But depends on quality of GOI function descriptions. |
| Graph genomes | **Not yet** | Engineering cost is extreme, tool ecosystem isn't ready, and multi-assembly approach captures most benefit. Revisit in 2-3 years. |
| Non-coding elements | **Not for SynVoy** | Fundamentally different search paradigm. Better as a companion tool. |
| BGC/operon search | **Interesting niche** | Compatible architecture but crowded field (antiSMASH). Low priority. |
| nf-core full port | **Not yet** | Too much refactoring for current development stage. Adopt standards incrementally. |
| nextflow_schema.json | **Yes, immediately** | Low effort, high UX payoff. Do this first. |

The vision in `feature_ideas.txt` is exactly right: SynVoy's future lies in crossing its synteny logic with deep-learning search methods. The synteny-first architecture gives SynVoy a unique advantage — it reduces the AI compute problem from "fold the entire proteome" to "fold 50 proteins in a small genomic window" — and that's what makes these features feasible where they wouldn't be for other tools.

---

## Sources

- [Foldseek — GitHub (steineggerlab)](https://github.com/steineggerlab/foldseek)
- [Foldseek — Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01773-0)
- [ESMFold — BioLM](https://biolm.ai/models/esmfold/)
- [ESMFold — Science 2022 (Lin et al.)](https://www.science.org/doi/10.1126/science.ade2574)
- [ESM-2 — GitHub (facebookresearch)](https://github.com/facebookresearch/esm)
- [ESMFold vs AlphaFold accuracy — NAR Genomics 2026](https://academic.oup.com/nargab/article/8/1/lqag002/8427121)
- [ESMFold vs AlphaFold vs OmegaFold — Frontiers in Genetics 2025](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1715037/full)
- [OpenFold-TRT — bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.03.11.711233v1)
- [ProTrek — Nature Biotechnology 2025](https://www.nature.com/articles/s41587-025-02836-0)
- [ProTrek — GitHub (westlake-repl)](https://github.com/westlake-repl/ProTrek)
- [PLM remote homology review — PMC 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12168796/)
- [Dense Homolog Retriever — Nature Biotechnology 2024](https://www.nature.com/articles/s41587-024-02353-6)
- [Minigraph-Cactus — Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01793-w)
- [Pangenome graph evaluation — GigaScience 2025](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf121/8364989)
- [vg Giraffe — bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.09.29.678807v1)
- [H100 rental prices 2026 — IntuitionLabs](https://intuitionlabs.ai/articles/h100-rental-prices-cloud-comparison)
- [AWS GPU price cuts 2026 — DCD](https://www.datacenterdynamics.com/en/news/aws-cuts-costs-for-h100-h200-and-a100-instances-by-up-to-45/)
- [Cloud GPU pricing 2026 — SynpixCloud](https://www.synpixcloud.com/blog/cloud-gpu-pricing-comparison-2026)
- [NVIDIA protein structure acceleration](https://developer.nvidia.com/blog/accelerate-protein-structure-inference-over-100x-with-nvidia-rtx-pro-6000-blackwell-server-edition/)
- [nf-core pipeline guidelines](https://nf-co.re/docs/guidelines/pipelines/overview)
- [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules)
