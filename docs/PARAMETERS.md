# SynVoy — Comprehensive Parameter Reference

This document provides in-depth explanations for every configurable parameter in SynVoy. Each entry describes the parameter's purpose, its effect on the pipeline, practical guidance for tuning, and the biological reasoning behind its default value. Parameters are grouped by functional category.

---

## Table of Contents

1. [Input & Mode Selection](#1-input--mode-selection)
2. [Easy Mode](#2-easy-mode)
3. [Pro Mode](#3-pro-mode)
4. [Synteny & Flanking Gene Extraction](#4-synteny--flanking-gene-extraction)
5. [Search Parameters](#5-search-parameters)
6. [Smith-Waterman Local Search](#6-smith-waterman-local-search)
7. [Augmented / Relaxed Search](#7-augmented--relaxed-search)
8. [MMseqs2 Configuration](#8-mmseqs2-configuration)
9. [Annotation & Gap Filling](#9-annotation--gap-filling)
10. [Gene Prediction](#10-gene-prediction)
11. [Gene Model Classification](#11-gene-model-classification)
12. [Synteny Scoring Weights](#12-synteny-scoring-weights)
13. [PLM (Protein Language Model) Search](#13-plm-protein-language-model-search)
14. [Structural Search (ESMFold + Foldseek)](#14-structural-search-esmfold--foldseek)
15. [Visualization](#15-visualization)
16. [Resource Tuning](#16-resource-tuning)
17. [LLM Parameter Estimation](#17-llm-parameter-estimation)
18. [Advanced & Output](#18-advanced--output)
19. [Reserved / Unimplemented](#19-reserved--unimplemented)

---

## 1. Input & Mode Selection

### `mode`
**Type:** String | **Default:** `'easy'` | **Options:** `'easy'`, `'pro'`

Controls which execution path the pipeline follows. In **Easy Mode**, SynVoy handles genome retrieval automatically: it resolves the query protein via UniProt or NCBI, fetches the home genome and GFF annotation from NCBI Datasets, selects related target genomes based on taxonomic relationships, and downloads them. This is the recommended mode for most users who want a turnkey experience. In **Pro Mode**, the user supplies all input files directly — the query FASTA, the home genome, an optional GFF, and target genome files. Pro Mode works entirely offline and gives the user full control over which genomes are searched. The mode choice fundamentally determines which downstream processes are activated: Easy Mode triggers genome fetching, assembly QC, and taxonomic auto-selection processes, while Pro Mode skips all of these and feeds files directly into the analysis pipeline.

### `query`
**Type:** Path (String) | **Default:** `null`

Path to a local query protein FASTA file. This is the primary input for Pro Mode and an alternative input for Easy Mode. The file should contain a single protein sequence (or DNA sequence, which SynVoy will auto-translate to the best open reading frame). The query protein defines the Gene of Interest (GOI) — the gene you want to find orthologs of across target genomes. When provided alongside `--mode easy`, this file is used instead of fetching the protein from UniProt/NCBI. The pipeline extracts the protein name, estimates exon count, and checks for signal peptides from this sequence. If the FASTA contains multiple sequences, only the first is used as the GOI. The file path should be absolute or relative to the Nextflow launch directory.

### `query_id`
**Type:** String | **Default:** `null`

A UniProt accession (e.g., `P01501`, `Q16553`) or NCBI protein ID used in Easy Mode to automatically fetch the query protein sequence, its annotation, and associated metadata. When SynVoy receives a `query_id`, it contacts the UniProt API (or NCBI Protein if UniProt fails) to download the full protein sequence, identify the source organism (which becomes the home species unless overridden), and gather functional annotations including signal peptides, domain families, and gene family size estimates. This metadata feeds into the LLM parameter estimation system, allowing the pipeline to tailor search parameters to the specific biology of your query protein. Using `query_id` is the simplest way to run SynVoy — a single accession plus a target species list is enough to launch a complete analysis.

### `query_seq`
**Type:** String | **Default:** `null`

An inline protein sequence or FASTA-formatted text passed directly on the command line, without needing a file on disk. This is useful for quick ad-hoc searches where creating a FASTA file would be unnecessary overhead. The sequence can be a raw amino acid string (e.g., `MKFLVNVALVFMVVYISYIY...`) or a complete FASTA block with a header line. When using `query_seq`, the pipeline cannot automatically determine the source organism, so **`--home_species` is required**. SynVoy writes the inline sequence to a temporary FASTA file internally, then proceeds identically to the `--query` path. This parameter is mutually exclusive with `--query` and `--query_id`; providing more than one query source will cause a validation error at pipeline startup.

---

## 2. Easy Mode

### `home_species`
**Type:** String | **Default:** `null` (auto-detected)

Specifies the reference organism whose genome will serve as the "home" genome — the genome in which the GOI is first located and from which flanking synteny anchors are extracted. In Easy Mode, this is normally auto-detected from the query protein's source organism in UniProt/NCBI. However, you can override it when the auto-detected species is wrong (e.g., the protein was characterized in a model organism but you want to use a different reference assembly), or when using `--query_seq` where no organism metadata exists. The species name must match an NCBI Taxonomy entry (e.g., `"Homo sapiens"`, `"Apis mellifera"`). SynVoy uses this to fetch the best available genome assembly from NCBI Datasets, query species-level genomic statistics (genome size, gene count, N50), and establish the phylogenetic root for target genome ordering.

### `max_genomes`
**Type:** Integer | **Default:** `0` (auto)

Controls how many related target genomes SynVoy fetches in Easy Mode. When set to `0` (the default), SynVoy uses an adaptive taxonomic sampling strategy: it walks up the taxonomic tree from the home species (genus → family → order → class → phylum) and selects approximately 3 representative species at each level, choosing those with the best available genome assemblies. This produces a phylogenetically balanced set that samples both close relatives (for high-confidence synteny) and distant relatives (for evolutionary reach). Setting a specific number (e.g., `--max_genomes 5`) caps the total number of target genomes fetched. Lower values speed up the pipeline and reduce downloads; higher values provide broader phylogenetic coverage but increase runtime linearly. For quick tests, `--max_genomes 1` with a specific `--target_species` is recommended.

> **Warning — single-genome searches are under-powered.** Synteny scoring derives its signal from *consensus across species*: when 7 flanking genes appear in the same order across 3+ target genomes, fallback GOI hits can be classified as `probable_goi`. With a single target genome there is no cross-species conservation evidence, so fallback hits almost always end up labelled `ambiguous_goi_family_member` or `tandem_goi_copy`, even when the correct locus has been found. The pipeline logs a warning when `max_genomes < 3`. Prefer `max_genomes >= 3` for any analysis where orthology calls matter.

### `target_species`
**Type:** String (comma-separated) | **Default:** `null` (auto-selected)

A comma-separated list of species names to use as target genomes instead of automatic taxonomic sampling. When provided, SynVoy fetches the best available genome assembly for each named species from NCBI Datasets, bypassing the taxonomic tree walker entirely. This is useful when you have specific organisms of interest (e.g., `"Gallus gallus,Mus musculus,Danio rerio"`) or when the automatic selection picks unsuitable genomes. Each species name must match an NCBI Taxonomy entry. If a species cannot be found or has no genome assembly available, SynVoy logs a warning and continues with the remaining species. This parameter can be combined with `--max_genomes` — if both are provided, the explicit species list takes priority and `max_genomes` serves as a cap on how many of the listed species are actually processed.

### `assembly_ranking`
**Type:** String | **Default:** `'hybrid'` | **Options:** `'hybrid'`, `'counts'`, `'nstats'`

Determines how SynVoy ranks multiple genome assemblies for the same species to select the best one. When NCBI has several assemblies for a species (common for model organisms), this parameter controls the selection heuristic. **`hybrid`** (recommended) combines both contiguity statistics (N50, scaffold count) and annotation completeness (gene count, RefSeq status) into a composite score, balancing assembly quality with annotation richness. **`counts`** ranks purely by annotation completeness — preferring assemblies with more annotated genes and RefSeq status, which is useful when you need GFF annotations for flanking gene extraction. **`nstats`** ranks purely by contiguity — preferring assemblies with higher N50 and fewer scaffolds, ideal when assembly completeness matters more than annotation. The hybrid approach works best in most cases because SynVoy benefits from both good contiguity (for accurate synteny blocks) and rich annotations (for flanking gene extraction).

### `bad_quality_policy`
**Type:** String | **Default:** `'drop'` | **Options:** `'drop'`, `'keep'`, `'ask'`

Controls what happens when a downloaded genome assembly fails quality checks (excessive contigs, low N50, or too many scaffolds as defined by the `bad_max_*` and `bad_min_*` thresholds). **`drop`** silently removes low-quality assemblies from the analysis — they are downloaded but not processed further. This is the safest default because fragmented assemblies produce unreliable synteny blocks: flanking genes may land on different scaffolds, breaking the contiguity assumption that underpins synteny analysis. **`keep`** processes all assemblies regardless of quality, which is appropriate when working with non-model organisms where all available assemblies are draft-quality. **`ask`** prompts the user interactively (via terminal) for each flagged assembly, showing quality metrics and waiting up to `bad_quality_timeout` seconds for a response. If no response is received, the assembly is dropped.

### `bad_quality_timeout`
**Type:** Integer (seconds) | **Default:** `300`

The number of seconds SynVoy waits for a user response when `bad_quality_policy` is set to `'ask'`. During this window, the pipeline pauses and displays the assembly name, contig count, scaffold count, and N50 value, then asks whether to keep or drop the assembly. If the timeout expires without user input, the assembly is dropped (fail-safe behavior). This parameter only has an effect when `bad_quality_policy = 'ask'`; it is ignored for `'drop'` and `'keep'` policies. A 5-minute default gives users enough time to evaluate the metrics, but you can increase it for batch runs where you might not be watching the terminal continuously. Setting it very high (e.g., 86400) effectively makes the pipeline wait indefinitely for input.

### `bad_max_contigs`
**Type:** Integer | **Default:** `500000`

The maximum number of contigs an assembly can have before it is flagged as low quality. Assemblies with more contigs than this threshold are subject to the `bad_quality_policy`. Contig count is a proxy for assembly fragmentation — a highly fragmented assembly has hundreds of thousands of small contigs, meaning genes and their flanking regions are often split across multiple pieces. This breaks synteny analysis because the pipeline cannot reliably identify the genomic neighborhood of a gene when that neighborhood spans multiple unlinked contigs. The default of 500,000 is permissive; most reasonable genome assemblies have far fewer contigs. Lowering this threshold (e.g., to 50,000) enforces stricter quality requirements, which improves synteny reliability but may exclude draft-quality genomes for non-model organisms.

### `bad_max_scaffolds`
**Type:** Integer | **Default:** `500000`

The maximum number of scaffolds an assembly can have before it is flagged as low quality. Similar to `bad_max_contigs`, this measures assembly fragmentation but at the scaffold level (scaffolds are ordered/oriented contigs joined by gaps). Scaffold count matters because SynVoy identifies syntenic blocks by clustering flanking gene hits that are nearby on the same scaffold. When a genome has too many scaffolds, flanking genes may be scattered across different scaffolds even if they are physically adjacent in the real genome, producing artificially low synteny scores and missed candidate regions. The default of 500,000 is generous. For chromosome-level assemblies this value is typically under 100, while for draft assemblies it might be in the thousands. Tighten this value if you want to restrict analysis to higher-quality assemblies only.

### `bad_min_n50`
**Type:** Integer (bp) | **Default:** `5000`

The minimum N50 value (in base pairs) an assembly must have to pass quality filtering. N50 is the length such that 50% of the total assembly is contained in contigs/scaffolds of this length or longer — it measures assembly contiguity. A low N50 means the assembly is composed of many short fragments, which is problematic for synteny analysis because candidate syntenic blocks (typically 100–500 kb) may be larger than the longest scaffold. The default of 5,000 bp is extremely permissive, accepting essentially any assembly. For reliable synteny analysis, an N50 of at least 100,000 bp (100 kb) is recommended. Consider setting this to 50,000–100,000 for stricter quality control. Assemblies below this threshold are handled according to `bad_quality_policy`.

### `qc_fail_policy`
**Type:** String | **Default:** `'drop'` | **Options:** `'drop'`, `'keep'`

Controls how assemblies that fail downstream quality checks (after the initial `bad_quality_*` screening) are handled. This is a second-tier quality gate applied after genome downloading and initial QC. While `bad_quality_policy` catches assemblies with obviously poor contiguity statistics, `qc_fail_policy` handles assemblies that pass initial screening but fail more detailed quality assessments during analysis — for example, assemblies where flanking gene extraction fails entirely or where the home gene cannot be located at all. Setting this to `'drop'` removes such assemblies from results; `'keep'` retains them with appropriate warnings in the output report. This two-tier approach lets users be permissive during initial selection (to avoid discarding potentially useful genomes too early) while still filtering out assemblies that demonstrably fail during analysis.

---

## 3. Pro Mode

### `home_genome`
**Type:** Path (String) | **Default:** `null`

Path to the reference genome FASTA file (`.fna` or `.fna.gz`) in Pro Mode. This genome serves as the "home" where the GOI is located and flanking synteny anchors are extracted. The quality of this genome directly impacts the entire analysis: a fragmented home genome will produce incomplete flanking gene sets, reducing synteny sensitivity across all targets. Chromosome-level or scaffold-level assemblies are strongly preferred. The genome can be gzip-compressed (`.fna.gz`) — SynVoy handles decompression automatically. The file must contain nucleotide sequences (not protein). SynVoy locates the GOI within this genome using tblastn and MMseqs2, then uses the GFF annotation (if provided via `--home_gff`) or Prodigal/Augustus gene prediction to identify the flanking genes that form the synteny anchors.

### `home_gff`
**Type:** Path (String) | **Default:** `null`

Path to a GFF3 annotation file for the home genome. While optional, providing a GFF is **highly recommended** because it dramatically improves flanking gene extraction. With a GFF, SynVoy can identify the exact gene models surrounding the GOI, extract their protein sequences with correct exon-intron boundaries, and compute precise GOI-similarity scores for filtering. Without a GFF, the pipeline falls back to ab initio gene prediction using Augustus (eukaryotes) or Prodigal (prokaryotes), which is less accurate — predicted genes may have incorrect boundaries, miss small exons, or merge adjacent genes. The GFF must be in GFF3 format and correspond to the same assembly as `--home_genome`. Features of type `CDS`, `gene`, and `mRNA` are parsed. Annotations from NCBI, Ensembl, or other standard sources work without modification.

### `target_genomes`
**Type:** String (glob/paths) | **Default:** `null`

Specifies the target genome files to search for GOI orthologs in Pro Mode. This parameter accepts three formats: (1) a glob pattern like `"genomes/*.fna"` that matches multiple files, (2) a comma-separated list of paths like `"genome_a.fna,genome_b.fna"`, or (3) a Nextflow list/channel expression. Each target genome is an independent nucleotide FASTA file (`.fna` or `.fna.gz`) representing a different species or strain. SynVoy processes targets in phylogenetic order (nearest to the home species first), using evolutionary distance estimates from NCBI taxonomy. The iterative search strategy means that genes discovered in closer relatives are added to the search database, improving sensitivity for more distant targets. Providing targets ordered from close to distant relatives maximizes this iterative benefit, though SynVoy handles reordering automatically when taxonomic information is available.

---

## 4. Synteny & Flanking Gene Extraction

### `n_flanking_genes`
**Type:** Integer | **Default:** `10` | **Range:** 2–30

The number of non-GOI-similar flanking genes to extract on **each side** of the GOI in the home genome. These flanking genes serve as synteny anchors — their presence in target genomes indicates conserved gene order (synteny). With the default of 10, SynVoy extracts up to 10 upstream and 10 downstream genes, for up to 20 total anchors. The optimal value depends on genome architecture: gene-dense genomes (bacteria, compact insects) need fewer anchors (5–8) because genes are closely packed and synteny is well-preserved, while gene-sparse or heavily rearranged genomes (plants, large vertebrates) benefit from more anchors (12–20) to increase the chance that at least some are conserved. More anchors increase computational cost (each must be searched in every target) but improve synteny detection reliability. The `min_synteny_score` threshold works relative to this count — with 10 flanking genes and a 0.6 score, at least 6 must be found.

### `prefer_large_genes`
**Type:** Boolean | **Default:** `true`

When enabled, the flanking gene extraction process preferentially selects larger genes as synteny anchors. Larger genes produce longer protein sequences that are more informative for homology search — they have more sequence to align, yielding higher-confidence MMseqs2 hits in target genomes. Small genes (e.g., short peptides, tRNAs) are poor synteny anchors because their short sequences produce weak, ambiguous alignments that may not pass identity and length thresholds. With this option enabled, the flanking gene selector sorts candidates by CDS length and prioritizes longer ones, while still respecting the spatial constraint of extracting from the immediate neighborhood of the GOI. Disabling this may be useful in very gene-dense regions where all neighboring genes are small (e.g., tightly packed operons in bacteria), but for most eukaryotic genomes the default is strongly recommended.

### `min_flanking_size`
**Type:** Integer (bp) | **Default:** `500` | **Range:** 100–5000

The minimum nucleotide length for a gene to be considered as a flanking synteny anchor. Genes shorter than this threshold are skipped during flanking gene extraction. This filter complements `prefer_large_genes` by providing a hard minimum floor. The default of 500 bp (~167 amino acids) excludes very small genes, pseudogenes, and short non-coding features that would make unreliable anchors. For bacteria and archaea, where genes are typically shorter, lowering this to 200–300 bp is recommended (the heuristic estimator does this automatically). For large eukaryotic genomes with multi-exon genes, the default is appropriate. Setting this too high risks excluding legitimate flanking genes in gene-dense regions, while setting it too low allows short, poorly conserved genes that produce noisy search results and inflate false-positive synteny scores.

### `exon_level_search`
**Type:** Boolean | **Default:** `true`

When enabled, flanking gene protein sequences are decomposed into individual exon-level fragments for the initial MMseqs2 search in target genomes. This dramatically improves sensitivity for detecting divergent orthologs because individual exons may be conserved even when full-length protein similarity has dropped below detection thresholds. Multi-exon genes that have accumulated mutations unevenly across exons benefit most — a single well-conserved exon can anchor the synteny block even if other exons have diverged beyond recognition. The pipeline searches both the full-length protein and individual exon sequences, taking the best hit. Disabling this saves some computational overhead but significantly reduces sensitivity for cross-order or cross-class searches (>200 Mya divergence). It is strongly recommended to keep this enabled unless you are searching only within the same genus or family where full-length matches are expected.

### `max_flanking_goi_similarity`
**Type:** Float (%) | **Default:** `35.0` | **Range:** 10–100

The maximum k-mer-based sequence similarity (as a percentage) that a flanking gene candidate may have to the GOI before it is excluded as a synteny anchor. This filter is critical for gene families: if the GOI belongs to a tandem array (e.g., MRJP cluster, LY6 family, defensins), its immediate neighbors may be paralogous copies that share significant sequence similarity. Using such paralogs as synteny anchors is counterproductive — they would match any region containing a gene family member, not just the true orthologous locus. The default of 35% strikes a balance between excluding obvious paralogs and retaining distant but unrelated neighbors. For large gene families with many tandem duplicates, lower this to 20–25%. For unique genes with no close paralogs, you can relax this to 50–100% (set to 100 to disable filtering entirely). The similarity is computed using a fast k-mer comparison, not full alignment.

### `max_flanking_distance`
**Type:** Integer (bp) | **Default:** `0` (unlimited) | **Range:** 0–5,000,000

The maximum distance in base pairs from the GOI center that the flanking gene extractor will walk to find synteny anchors. When set to 0 (default), there is no distance limit — the extractor walks as far as needed along the chromosome to collect `n_flanking_genes` anchors per side. This can be problematic when the GOI sits adjacent to a large tandem gene array: one direction fills up entirely with GOI-similar genes (which are filtered out by `max_flanking_goi_similarity`), forcing the extractor to walk far into a distant genomic region to find enough non-similar anchors. These distant anchors may belong to a different synteny block entirely, confusing downstream analysis. Setting a distance cap (e.g., 500,000 bp) prevents this pathological case. If fewer than `n_flanking_genes` anchors are found within the distance window, the pipeline proceeds with whatever was found. This parameter is automatically adjusted by the LLM estimator for known gene family cases.

### `expand_goi_similar`
**Type:** Boolean | **Default:** `false` (was `true` until 2026-04-25)

When enabled, genes near the GOI that are similar to it (detected during flanking gene extraction) are not just excluded as anchors — they are emitted as additional GOI queries with a `GOI_NEIGHBOR_` prefix. These neighbor queries are searched in all target genomes alongside the original GOI and included in the phylogenetic tree.

**Why default changed to `false`:** for paralog-rich gene families (e.g. LY6 with ~15 chr8 paralogs) the NEIGHBOR queries flood the per-target m8 hits with paralog matches. Every cluster picks up a GOI-overlap bonus from at least one paralog, so the bonus loses its discriminating power. The true ortholog scaffold drops out of `adaptive_max_regions=6`. Empirical: LY6 ground-truth scaffold recall went from **8/9 (off)** to **0/9 (on)**. Melittin (single-copy gene, no nearby paralogs) is unaffected either way.

**When to turn it on:** you have a single-copy GOI sitting next to a small known tandem array, and you want the paralogs in the same tree to disambiguate. Examples: a Yellow gene with one Yellow-e3 paralog two genes away, a venom peptide flanked by 1–2 tandem duplicates. In these narrow cases the bonus does help.

**When to leave it off:** any gene-family query (LY6, TP53, MRJP, 3FTx, defensins) — the paralog flood will swamp the synteny signal. The phylogenetic tree downstream still resolves paralogs from the iterative search's existing GOI-similarity classification; you don't need this flag to get a paralog-aware tree.

The similarity threshold for detecting GOI-like neighbors reuses `max_flanking_goi_similarity`.

### `expand_goi_similar_distance`
**Type:** Integer (bp) | **Default:** `300000` | **Range:** 10,000–2,000,000

The maximum distance in base pairs from the GOI within which to search for GOI-similar neighbor genes for expansion (see `expand_goi_similar`). Genes beyond this distance that happen to resemble the GOI are not emitted as additional queries. The default of 300 kb covers most tandem gene arrays, which typically span 50–500 kb in eukaryotic genomes. For very large tandem arrays (e.g., olfactory receptor clusters in vertebrates, which can span several megabases), you may need to increase this value. For compact genomes (bacteria, fungi), a smaller value (50–100 kb) is more appropriate to avoid picking up unrelated genes that coincidentally share sequence similarity. Each additional GOI-similar neighbor adds a full search pass across all target genomes, so expanding this distance excessively can significantly increase runtime without improving results.

---

## 5. Search Parameters

### `cluster_distance`
**Type:** Integer (bp) | **Default:** `150000` | **Range:** 1,000–2,000,000

The maximum gap in base pairs between flanking gene hits on a target genome scaffold before they are considered part of separate candidate syntenic blocks. When flanking gene proteins are searched against a target genome, the resulting hits are clustered by genomic proximity: hits within `cluster_distance` of each other are merged into a single candidate block. This parameter must scale with genome architecture. In compact bacterial genomes where genes are tightly packed in operons (intergenic distances <1 kb), a small value (10,000–30,000 bp) is appropriate. In vertebrate genomes where intergenic regions span hundreds of kilobases, a larger value (200,000–500,000 bp) is needed to avoid splitting a single syntenic block into many fragments. In plants with very large genomes and huge intergenic/intronic regions, values up to 1,000,000 bp may be necessary. Setting this too small causes real syntenic blocks to be fragmented; too large causes unrelated hits on the same chromosome to merge into a single spurious block.

### `min_synteny_score`
**Type:** Float | **Default:** `0.6` | **Range:** 0.1–1.0

The minimum fraction of flanking gene anchors that must have hits in a candidate syntenic block for the block to be considered valid and trigger a local GOI search. With the default of 0.6 and 10 flanking genes per side (20 total), at least 12 anchors must map to the target region. This is the primary precision/recall knob for synteny analysis. Higher values (0.7–0.8) demand strong synteny conservation, giving high confidence but missing orthologous regions where chromosomal rearrangements have disrupted gene order — common in distantly related species or lineages with high rearrangement rates (e.g., Drosophila, plants). Lower values (0.2–0.4) accept regions with weaker synteny, casting a wider net that catches rearranged orthologs but also produces more false-positive candidate regions that must be evaluated. The LLM estimator adjusts this based on evolutionary distance: close species get strict scores, distant species get relaxed scores.

### `min_hit_identity`
**Type:** Float (%) | **Default:** `10` | **Range:** 0–100

The minimum percent identity for an individual alignment hit to be retained during the initial flanking gene and GOI search phases. This is a coarse filter applied early in the pipeline to remove obvious noise from BLAST/MMseqs2 results. The default of 10% is intentionally very permissive, allowing even highly divergent homologs through the initial filter. Subsequent scoring (synteny score, Smith-Waterman, classification) provides the actual quality control. Raising this value (e.g., to 20–30%) speeds up the pipeline by discarding more hits early, but risks losing legitimate divergent orthologs in distant species searches. Lowering it below 10% is rarely useful because hits at <10% identity are typically random noise even for protein alignments. The LLM estimator may lower this for small peptides that naturally produce lower identity scores due to their short length.

### `min_query_length`
**Type:** Integer (amino acids) | **Default:** `30` | **Range:** 0–∞

The minimum length (in amino acids) a query protein must have before the pipeline will proceed. Queries shorter than this threshold are rejected with a non-zero exit code during `NORMALIZE_QUERY`, with a message explaining why. This is a pre-flight safety check: MMseqs2 and tblastn produce unreliable noise on very short queries (<30 aa), and short peptides disproportionately confuse the iterative-search and synteny-scoring logic, which are both tuned around query-length heuristics. If you genuinely need to search with a short fragment (e.g. a conserved motif, or a micro-exon probe), override with `--min_query_length 0` to disable the check. The check runs once in `normalize_query.py`, *after* DNA-to-protein translation, so a 60-bp DNA query that translates to a 20-aa ORF is still caught. Typical threshold choices: 30 aa for conservative runs, 20 aa when searching defensins and small antimicrobial peptides, 0 to disable entirely.

### `min_hit_length`
**Type:** Integer (amino acids) | **Default:** `10` | **Range:** 5–500

The minimum alignment length in amino acid residues for an individual hit to be retained. Short alignments, even with high identity, are often spurious matches to low-complexity regions or common protein motifs rather than genuine homology. The default of 10 amino acids is permissive enough to detect even small conserved domains while filtering out trivially short alignments. For large proteins (>500 aa), raising this to 30+ reduces noise without losing real hits. For small peptides (<100 aa), this should be lowered to 8 — a 70-amino-acid peptide might only produce a 15-residue alignment against a divergent ortholog, and a 10-residue minimum is already close to the true signal boundary. This parameter interacts with `min_hit_identity`: together they define the minimum quality envelope for retaining initial search hits. The LLM estimator automatically lowers this for small query proteins.

### `search_evalue`
**Type:** Float | **Default:** `0.01` | **Range:** 1e-10–100

The E-value threshold for tblastn and MMseqs2 homology searches. The E-value represents the expected number of hits with equal or better score by chance in a database of the given size. Lower values (0.001, 1e-5) are more stringent, reporting only highly significant hits; higher values (0.1, 1.0, 10) are more permissive, allowing weaker matches through. The default of 0.01 provides a good balance for typical within-class or within-order searches. For very close species (same genus), tightening to 0.001 reduces noise. For distant cross-phylum searches (>400 Mya), relaxing to 0.1–1.0 is necessary because genuine orthologs may have diverged enough that their alignment scores are only marginally above background. The pipeline applies additional filters (synteny score, identity, length) downstream, so a permissive E-value here is safer than a restrictive one — it lets the synteny context disambiguate weak hits rather than discarding them prematurely.

### `max_intron`
**Type:** Integer (bp) | **Default:** `20000` | **Range:** 0–500,000

The maximum intron length in base pairs for gene model construction, used by miniprot when building gene models from protein-to-genome alignments and by the exon clustering logic. This parameter must match the biology of the organisms being searched. Bacterial genes have no introns (set to 0). Fungal introns are typically short (50–500 bp). Insect introns range from hundreds to tens of thousands of base pairs. Vertebrate introns can exceed 100,000 bp (the dystrophin gene has introns >400 kb). Plant introns are highly variable, from hundreds of base pairs in compact genes to over 100,000 bp in large genes. Setting this too small causes multi-exon genes to be fragmented into separate incomplete gene models, because miniprot cannot bridge exons across gaps larger than this value. Setting it too large risks merging exons from adjacent genes into a single spurious gene model. The LLM estimator sets kingdom-appropriate values automatically.

---

## 6. Smith-Waterman Local Search

### `enable_smith_waterman`
**Type:** Boolean | **Default:** `true`

Enables rigorous Smith-Waterman (SW) local alignment for GOI detection within candidate syntenic regions. After flanking gene clustering identifies candidate blocks, the pipeline extracts the genomic sequence of each block and searches for the GOI using three methods: tblastn, miniprot, and (when enabled) Smith-Waterman alignment. SW alignment via parasail or ssearch36 is the most sensitive local alignment method — it guarantees finding the optimal local alignment between the query protein and every possible six-frame translation of the region. This catches GOI copies that tblastn misses due to frameshift-tolerant scoring, compositional bias, or extreme divergence. The cost is runtime: SW alignment is O(m*n) in sequence lengths and can be slow for large regions. It is strongly recommended to keep this enabled; disabling it only makes sense if runtime is critical and you are searching very close relatives where tblastn alone suffices.

### `sw_method`
**Type:** String | **Default:** `'auto'` | **Options:** `'auto'`, `'parasail'`, `'ssearch36'`

Selects the Smith-Waterman implementation. **`auto`** (recommended) tries parasail first, falling back to ssearch36 (from the FASTA suite) if parasail is not installed. **`parasail`** is a fast SIMD-accelerated library that uses SSE/AVX instructions for hardware-accelerated alignment; it is typically 10–100x faster than ssearch36 and is the preferred choice. **`ssearch36`** is the classic reference implementation from the FASTA package — slower but universally available and well-validated. Both produce identical optimal alignments (Smith-Waterman is an exact algorithm), so the choice affects only performance. If parasail is available in your environment, the auto setting will always prefer it. Manual override to `ssearch36` might be useful if parasail produces unexpected behavior on a specific platform, but this has not been observed in practice.

### `sw_min_score`
**Type:** Float | **Default:** `20` | **Range:** 5–200

The minimum Smith-Waterman alignment score for a GOI candidate hit to be reported. SW scores depend on the scoring matrix (typically BLOSUM62), gap penalties, and the length and quality of the alignment. A score of 20 corresponds to roughly a 10-residue perfect match or a 20-residue match at ~50% identity — this is the minimum threshold where alignment scores are distinguishable from random noise for typical-length proteins. For small peptides (<100 aa), the maximum achievable SW score is inherently lower because there is less sequence to align, so this threshold must be lowered to 10–15 to avoid missing genuine hits entirely. For large proteins (>500 aa), raising this to 30–50 reduces false positives without sacrificing sensitivity. The LLM estimator automatically adjusts this based on query protein length, which is one of the most impactful automatic parameter tuning decisions.

### `sw_min_identity`
**Type:** Float (%) | **Default:** `10.0` | **Range:** 0–100

The minimum percent identity for Smith-Waterman alignment hits. This works alongside `sw_min_score` as a secondary quality filter. Even when an alignment score is above the threshold, if the overall identity is very low, the hit may be a false positive arising from a compositionally biased region rather than genuine homology. The default of 10% is very permissive, which is appropriate because SW alignment scores already incorporate substitution quality information — a high-scoring low-identity hit typically indicates a long alignment with many conservative substitutions, which can still represent genuine distant homology. Raising this threshold (e.g., to 20–30%) is appropriate for same-family or same-genus searches where you expect clear sequence similarity. Keeping it low is important for cross-class or cross-phylum searches where genuine orthologs may share less than 20% sequence identity.

### `sw_timeout_seconds`
**Type:** Integer | **Default:** `300` | **Range:** 1–3600

The maximum time in seconds allowed for a single Smith-Waterman search invocation. SW alignment is O(m*n) where m is the query length and n is the target region length; for large regions (>1 Mb) or large query proteins (>1000 aa), this can take several minutes even with SIMD-accelerated parasail. The timeout prevents a single pathological case (e.g., a highly repetitive region) from stalling the entire pipeline. When the timeout is hit, the SW search for that specific region is aborted and the pipeline proceeds with whatever tblastn and miniprot hits were found. The default of 300 seconds (5 minutes) is sufficient for most cases. Increase this for very large query proteins or when running on slow hardware; decrease it if you prioritize throughput over exhaustive sensitivity. The timeout applies per-region, not per-genome.

---

## 7. Augmented / Relaxed Search

### `region_padding`
**Type:** Integer (bp) | **Default:** `150000` | **Range:** 5,000–1,000,000

The amount of extra flanking sequence in base pairs appended to each side of a candidate syntenic block before gene prediction and GOI search are performed. After clustering flanking gene hits into a candidate block, the pipeline extends the block by `region_padding` on each end. This extension captures the GOI even when it lies just outside the flanking gene cluster (e.g., at the edge of a syntenic block that has been partially disrupted by a chromosomal rearrangement) and provides sufficient context for gene prediction tools to build complete gene models. The value should scale with genome architecture: use ~0.5–1x the `cluster_distance` value. Compact genomes (bacteria) need only 15,000–20,000 bp of padding, while large plant or vertebrate genomes may need 200,000–400,000 bp. Setting this too small risks missing the GOI when it is not tightly surrounded by the synteny anchors; too large wastes compute on gene prediction in irrelevant regions.

### `padding_min`
**Type:** Integer (bp) | **Default:** `50000` | **Range:** 5,000–500,000

The minimum padding in base pairs applied to candidate syntenic blocks, regardless of the adaptive padding algorithm's output. The pipeline uses adaptive padding that considers the block size, flanking gene density, and region characteristics. This parameter sets a floor to prevent the padding from becoming too small in gene-dense regions where the adaptive algorithm might calculate a very tight window. The minimum ensures there is always enough surrounding sequence for gene prediction tools to work with and for the GOI to be captured even if it lies slightly outside the core syntenic block. For bacteria, 5,000 bp is appropriate; for most eukaryotes, the default of 50,000 bp works well; for large genomes (vertebrates, plants), 100,000 bp is recommended. This value must always be less than or equal to `padding_max` — the pipeline will auto-fix violations but log a warning.

### `padding_max`
**Type:** Integer (bp) | **Default:** `200000` | **Range:** 10,000–1,000,000

The maximum padding in base pairs applied to candidate syntenic blocks, capping the adaptive padding algorithm. This prevents excessively large search regions that would slow down gene prediction and Smith-Waterman alignment without meaningfully improving sensitivity. Very large padding creates regions with thousands of predicted genes, most of which are irrelevant — this increases the false-positive rate for GOI candidates and multiplies downstream computation. The cap should scale with genome architecture: bacterial genomes need at most 30,000 bp; compact animal genomes (insects) work well with 200,000 bp; large vertebrate and plant genomes may need 350,000–500,000 bp. This value must always be greater than or equal to `padding_min`. The `detect_breaking_combos()` validator will auto-fix violations by setting `padding_max = padding_min + 50000`.

### `aug_relaxed_evalue_mult`
**Type:** Float | **Default:** `1000` | **Range:** 10–100,000

The multiplier applied to the base `search_evalue` during relaxed/augmented search passes. When the standard search fails to find the GOI in a candidate syntenic block, the pipeline enters progressively more relaxed search rounds. In the relaxed pass, the E-value threshold becomes `search_evalue * aug_relaxed_evalue_mult`, subject to `aug_relaxed_evalue_cap`. With the default E-value of 0.01 and this multiplier of 1000, the relaxed E-value becomes 10 (but is capped at `aug_relaxed_evalue_cap`). This aggressive relaxation is intentional — by the time the relaxed pass runs, the pipeline already has synteny evidence suggesting the GOI should be in this region, so a very permissive E-value is appropriate because false positives are constrained by the syntenic context. Lower this value if you want tighter control over relaxed-pass hits, but this may cause the pipeline to miss highly divergent orthologs.

### `aug_relaxed_evalue_cap`
**Type:** Float | **Default:** `10.0` | **Range:** 0.1–100

The absolute maximum E-value allowed even during the most relaxed search passes. Regardless of how high `search_evalue * aug_relaxed_evalue_mult` evaluates to, the actual E-value used in relaxed searches will never exceed this cap. The default of 10.0 is very permissive — an E-value of 10 means you expect 10 equally good hits by chance in the database, which would normally indicate noise. However, in the context of augmented search within a synteny-validated region, this permissiveness is acceptable because the syntenic context provides strong prior evidence that the GOI exists in the region. The cap prevents completely degenerate searches (E-value of 1000+) that would report essentially every possible alignment. Lowering this to 1.0 makes the relaxed pass more conservative; values below 0.1 effectively disable the augmented search benefit.

### `aug_relaxed_parse_evalue_mult`
**Type:** Float | **Default:** `10`

A secondary E-value multiplier used for parsing (post-filtering) relaxed search results. After the relaxed search with the wide E-value window produces hits, this multiplier is applied to the base E-value for initial hit parsing and filtering. This creates a two-stage filter: the search itself uses a very wide E-value to ensure all possible hits are found, then the parsing stage applies a tighter (but still relaxed) E-value to discard the worst noise before more expensive downstream processing. With `search_evalue=0.01` and this multiplier of 10, the parsing E-value is 0.1 — much tighter than the search E-value (10.0) but looser than the standard E-value (0.01). This intermediate filtering step reduces the number of hits that need to be processed by Smith-Waterman alignment and gene prediction, improving runtime without sacrificing sensitivity for borderline hits.

### `aug_relaxed_identity_factor`
**Type:** Float | **Default:** `0.6` | **Range:** 0.1–1.0

The factor by which `min_hit_identity` is multiplied during relaxed search passes. With the default identity threshold of 10% and this factor of 0.6, the relaxed identity threshold becomes 6%. This allows hits with lower sequence identity through the initial filter during augmented search, which is necessary for finding highly divergent orthologs in distant species. The rationale is the same as for E-value relaxation: when synteny evidence already supports the presence of the GOI in a region, accepting lower-identity hits is justified because the spatial context reduces the false-positive risk. Lowering this factor (e.g., to 0.3) makes the relaxed pass even more permissive; raising it (e.g., to 0.8) keeps the relaxed pass closer to the standard thresholds. The absolute minimum identity is bounded by `aug_relaxed_identity_min`.

### `aug_relaxed_identity_min`
**Type:** Float (%) | **Default:** `15.0` | **Range:** 5–50

The absolute minimum percent identity for hits in the relaxed search pass, regardless of how low `min_hit_identity * aug_relaxed_identity_factor` evaluates to. This hard floor prevents degenerate matches — even in the most permissive relaxed search, a hit must show at least this much sequence identity to be considered. The default of 15% is set to be above the "twilight zone" of protein sequence similarity (~20–35% for reliable homology inference) but below it, reflecting that within a synteny-validated region, lower-identity matches have higher prior probability of being real. Genuine distant orthologs at 15% identity are borderline but detectable with additional evidence (synteny, domain architecture, structural similarity). Setting this below 10% risks accepting random alignments; setting it above 25% defeats the purpose of the relaxed pass for distant searches.

### `aug_relaxed_length_div`
**Type:** Integer | **Default:** `2`

The divisor applied to `min_hit_length` during relaxed search passes. With the default minimum length of 10 and this divisor of 2, the relaxed minimum length becomes 5 amino acids. This allows shorter alignment fragments through the filter, which is important for finding divergent orthologs where only a small conserved core domain remains alignable. Short fragments are individually weak evidence, but when multiple short fragments map to the same genomic region that already has synteny support, the collective evidence can be compelling. The divisor approach ensures the relaxation scales proportionally with the base threshold — if the base is raised to 30 for large proteins, the relaxed threshold becomes 15, which is still meaningful. This is bounded by `aug_relaxed_length_min` to prevent trivially short matches.

### `aug_relaxed_length_min`
**Type:** Integer (aa) | **Default:** `15`

The absolute minimum alignment length (in amino acids) for hits in the relaxed search pass, regardless of how low `min_hit_length / aug_relaxed_length_div` evaluates to. This hard floor prevents accepting trivially short alignments that are almost certainly noise. An alignment of 15 amino acids provides enough residues for a meaningful substitution pattern — below this, even the distinction between homologous and chance similarity becomes unreliable. With BLOSUM62 scoring, a 15-residue alignment needs roughly 40% identity to achieve statistical significance, providing a natural quality floor even at the permissive E-values used in relaxed search. This parameter is particularly important when the base `min_hit_length` has been lowered for small peptides — the relaxed pass should not accept alignments shorter than 15 residues even for very small query proteins.

### `aug_dedup_bin_bp`
**Type:** Integer (bp) | **Default:** `100`

The bin size in base pairs for deduplicating overlapping hits from relaxed search passes. When the augmented search produces multiple overlapping alignments to the same genomic region (common when the E-value is relaxed), hits whose start positions fall within the same `aug_dedup_bin_bp`-sized bin are deduplicated, keeping only the best-scoring hit per bin. The default of 100 bp is fine-grained enough to distinguish truly separate hits (e.g., exons separated by introns) while collapsing redundant alignments to the same location. Increasing this value (e.g., to 500 bp) collapses more aggressively, which can help when many noisy overlapping hits slow downstream processing. Decreasing it (e.g., to 50 bp) preserves more hits but may retain redundant alignments. This parameter has minimal impact on final results and is primarily a performance tuning knob.

### `max_blocks_per_genome`
**Type:** Integer | **Default:** `80` | **Range:** 5–500

A safety cap on the total number of candidate syntenic blocks evaluated per target genome. In genomes with many repetitive regions or large gene families, flanking gene hits may cluster into hundreds of candidate blocks. Evaluating each block involves gene prediction, tblastn, SW alignment, and classification — so the cost is approximately linear in block count. The default of 80 is sufficient for even complex genomes with multiple GOI loci (e.g., a gene family with 10 tandem copies on 5 chromosomes = 50 blocks). For very large gene families or whole-genome searches, increase to 120–200. For quick exploratory runs, decrease to 20–30. When the cap is hit, blocks are prioritized by synteny score (highest first), so the most promising regions are always evaluated. The cap prevents runaway searches in pathological genomes from stalling the entire pipeline for hours.

### `min_block_genes`
**Type:** Integer | **Default:** `2` | **Range:** 1–10

The minimum number of flanking gene hits required in a candidate syntenic block for it to be kept and evaluated further. Blocks with fewer than this many hits are discarded as noise. With the default of 2, at least two flanking gene anchors must map to the same target region (within `cluster_distance`). Setting this to 1 is extremely permissive — a single flanking gene hit somewhere on a chromosome does not constitute synteny evidence and is almost certainly a standalone ortholog of that flanking gene rather than a conserved syntenic block. Setting this higher than 3 is very strict and may miss real syntenic regions in organisms with high chromosomal rearrangement. The combination of `min_block_genes` and `min_synteny_score` together determines which candidate regions survive for GOI search: `min_block_genes` is the absolute floor, while `min_synteny_score` is the proportional threshold.

### `max_consecutive_empty_blocks`
**Type:** Integer | **Default:** `25` | **Range:** 5–100

Controls the iterative search early-stopping behavior. When evaluating candidate syntenic blocks in order of decreasing synteny score, if this many consecutive blocks fail to produce a GOI hit (i.e., are "empty"), the pipeline stops evaluating remaining blocks for that target genome. This is a heuristic optimization: if the top 25 candidate blocks (by synteny score) all fail to contain the GOI, it is very unlikely that lower-scored blocks will succeed. The default of 25 is generous — even in complex cases with many false-positive syntenic blocks, the true positive is almost always found within the top 10–15 candidates. Reducing this to 10 speeds up processing but risks missing GOIs in unusual genomic contexts (e.g., highly rearranged genomes where the true orthologous region has a mediocre synteny score). Increasing it beyond 50 is rarely beneficial.

---

## 8. MMseqs2 Configuration

### `mmseqs_sensitivity`
**Type:** Float | **Default:** `9.5` | **Range:** 1–12

Controls the sensitivity of MMseqs2 prefilter and alignment stages, directly trading speed for detection power. MMseqs2 uses a k-mer-based prefilter to identify candidate target sequences before running full alignment; higher sensitivity values use more k-mers and longer seed matching, finding more distant homologs at the cost of slower execution. The scale runs from 1 (fastest, least sensitive — suitable for finding near-identical sequences) to 12 (slowest, most sensitive — comparable to BLAST with all optimizations). The default of 9.5 is already quite high, suitable for finding orthologs within the same class or order. For cross-phylum searches (>400 Mya), increase to 10–11. For very distant searches (>700 Mya, cross-kingdom), push to 11–12 and accept significant runtime increase. For within-genus searches, 7–8 is sufficient and much faster. This parameter has the largest single impact on runtime of any sensitivity control.

### `mmseqs_split_memory_limit`
**Type:** String | **Default:** `'8G'` (config default)

The memory limit for MMseqs2 database splitting. When the target database is larger than this limit, MMseqs2 splits it into chunks and processes each chunk sequentially, reducing peak memory usage at the cost of slightly increased runtime (due to repeated I/O). The default of 8 GB works well for most workstations. On memory-constrained systems (laptops with 8 GB RAM), reduce to `'3G'` or `'1G'` to prevent OOM kills — the `laptop_safe` profile already sets this to `'8G'`. On servers with ample RAM (32+ GB), you can increase to `'16G'` or higher to avoid splitting altogether, which slightly improves performance. If MMseqs2 processes crash with exit code 137 (SIGKILL from OOM killer), this is the first parameter to reduce. The value must include a unit suffix (`G` for gigabytes, `M` for megabytes).

### `mmseqs_verbosity`
**Type:** Integer | **Default:** `1` | **Options:** `0`, `1`, `2`, `3`

Controls the verbosity of MMseqs2 log output. Level 0 suppresses all output except errors; level 1 shows warnings and basic progress; higher levels show detailed progress information including per-sequence alignments. The default of 1 provides useful diagnostic information without flooding logs. The `laptop_safe` profile sets this to 0 to reduce log noise during single-threaded execution. For debugging failed searches, temporarily set to 2 or 3 to see which sequences are being searched and what hits are found. This parameter does not affect results, only logging. It is passed directly to MMseqs2's `--verbosity` flag. Note that even at verbosity 0, critical errors are still reported through Nextflow's process error handling.

### `min_gene_identity`
**Type:** Float (%) | **Default:** `30` | **Range:** 5–100

The minimum percent identity for flanking gene MMseqs2 hits to be considered valid synteny anchors. This is stricter than `min_hit_identity` because flanking genes serve a specific purpose: they must be confidently identified as orthologs in the target genome to provide reliable synteny evidence. At 30% identity, protein sequences are in the "safe zone" above the twilight zone of homology detection (20–35%), meaning the homology relationship is reliable. For close species (same genus or family), this can be increased to 40% to reduce noise. For distant searches (cross-order, >200 Mya), lowering to 15–20% is necessary because flanking genes diverge faster than well-known marker genes. Setting this too low allows non-homologous matches to masquerade as synteny anchors, inflating synteny scores; too high excludes real anchors in divergent species, reducing synteny sensitivity.

---

## 9. Annotation & Gap Filling

### `gff_search_window`
**Type:** Integer (bp) | **Default:** `100000`

The window size in base pairs around the GOI position within which SynVoy searches the GFF annotation file for flanking genes. When a GFF is provided (via `--home_gff`), the pipeline locates the GOI in the home genome by coordinate overlap, then extracts annotated gene features within this window on each side. The default of 100,000 bp (100 kb) typically captures 5–15 genes per side in eukaryotic genomes, which is more than enough for the default `n_flanking_genes=10`. For gene-dense genomes (bacteria: ~1 gene per kb), a smaller window (20,000 bp) suffices. For gene-sparse large genomes where flanking genes may be very far apart, increasing to 200,000–500,000 bp may be necessary. This parameter only affects flanking gene extraction from the home genome GFF — target genome analysis uses different windowing controlled by `cluster_distance` and `region_padding`.

### `gap_search_window`
**Type:** Integer (bp) | **Default:** `50000`

The window size for gap-filling searches within candidate syntenic regions. During GOI search within a candidate block, if the initial tblastn/MMseqs2 search finds partial hits (e.g., only some exons), the gap-filling module searches within this window around existing hits for additional exons or fragments that may have been missed. The default of 50,000 bp covers the typical span of a multi-exon gene in most eukaryotes. For organisms with very large introns (some vertebrate genes span >500 kb), increase this to capture all exons. For bacteria and fungi where genes are compact, 10,000–20,000 bp is sufficient. Gap filling is particularly important for multi-exon genes in distant species where some exons have diverged beyond initial detection thresholds — the gap filler uses relaxed parameters to rescue these missing exons.

### `gap_min_size`
**Type:** Integer (bp) | **Default:** `10`

The minimum gap size in base pairs between existing hits for the gap-filling module to attempt a search. Gaps smaller than this are considered too small to contain a missing exon and are skipped. The default of 10 bp is extremely permissive — almost any gap between hits triggers a gap-filling search. This is appropriate because even very small exons (microexons of 3–30 bp) exist in some gene families, and the gap-filling search cost is low since it only examines a small region. For most practical purposes, this parameter does not need adjustment. Setting it higher (e.g., 50–100 bp) would skip tiny gaps but might miss microexons in genes known to contain them (e.g., some ion channels, cadherins).

### `gap_evalue`
**Type:** Float | **Default:** `10`

The E-value threshold for gap-filling searches. This is deliberately very permissive — much more so than the main `search_evalue` (default 0.01). The rationale is that gap filling occurs within a region already validated by synteny evidence, so we expect the GOI to be there. The search target is small (only the gap region), which means even weak alignments to short exon fragments are worth capturing. An E-value of 10 in a small search space has a different significance than E-value of 10 in a whole-genome search; the effective false-positive rate is controlled by the restricted search region. This permissive threshold enables rescue of highly divergent exons that would be invisible under standard thresholds. Tightening this value risks losing the marginal exons that gap filling is specifically designed to recover.

### `gap_min_identity`
**Type:** Float (%) | **Default:** `15.0`

The minimum percent identity for gap-filling search hits. This is slightly more permissive than the main `min_hit_identity` (default 10%) to allow detection of highly divergent exon fragments. At 15% identity, alignments are in the deep twilight zone, but within a synteny-validated region, even such weak matches can represent genuine exon fragments. The gap-filling module combines identity with the gap context (position relative to other confirmed exons, frame consistency, splice site signals) to evaluate whether a weak hit is a genuine exon. Lowering this below 10% is not recommended as alignments below that threshold are indistinguishable from random for protein sequences. Setting it higher (e.g., 25%) would be appropriate for within-family searches where you expect clear sequence conservation even in exon fragments.

### `gap_min_alnlen`
**Type:** Integer (aa) | **Default:** `10`

The minimum alignment length in amino acids for gap-filling hits. Combined with `gap_min_identity`, this defines the quality floor for rescued exon fragments. An alignment of 10 amino acids represents approximately 30 nucleotides of coding sequence — this is short enough to capture microexons and highly diverged exon fragments but long enough to have some statistical significance. For multi-domain proteins where individual exons encode distinct domains, even short alignments to a single domain's fragment provide valuable positional information. For most use cases, the default is appropriate. Increase to 20–30 if you want to restrict gap filling to more confident hits, at the cost of potentially missing small exons in divergent species.

### `gap_max_hits`
**Type:** Integer | **Default:** `5`

The maximum number of gap-filling hits to report per gap. When the gap-filling search finds multiple candidate fragments in a single gap region, only the top N (by score) are retained. This prevents a single large gap in a repetitive region from generating hundreds of weak candidate fragments that would overwhelm downstream processing. The default of 5 is usually sufficient — in a genuine gap between exons, there is typically at most one real missing exon plus a few noise hits. Increasing this is rarely helpful because the additional hits beyond the top 5 are almost always noise. Decreasing to 1 forces the pipeline to commit to the single best hit per gap, which is fine for clean genomes but may miss alternative exon configurations in polyploid or duplicated regions.

### `min_exon_query_cov`
**Type:** Float | **Default:** `0.25` | **Range:** 0–1

The minimum query coverage fraction for an exon alignment to be included in the gene model annotation. Query coverage measures what fraction of the query protein (GOI) a single exon alignment covers. With a default of 0.25, an individual exon must cover at least 25% of the query protein length to be included. For multi-exon genes, individual exons may cover much less than 25% each — a 5-exon gene might have exons covering 15–25% each. In such cases, this threshold may exclude some smaller exons from the annotation. For single-exon genes, this effectively requires the alignment to cover at least a quarter of the protein. Lowering this allows smaller exon fragments to be annotated, increasing sensitivity but potentially including short spurious matches. This threshold affects annotation completeness, which in turn affects the gene model classification (complete/partial/fragment).

### `min_exon_alnlen`
**Type:** Integer (aa) | **Default:** `30`

The minimum alignment length for an individual exon-level alignment to be included in gene model annotation. This provides a hard floor independent of query coverage — even if an alignment covers 50% of a very small query, it must still span at least 30 amino acids to be included. This prevents trivially short alignments from being annotated as exons. The value of 30 amino acids (approximately 90 nucleotides) is chosen as the approximate minimum length for reliable protein sequence homology detection — alignments shorter than this have high false-positive rates even at high identity. For very small query proteins (<80 aa), this threshold may be too strict; consider lowering to 15–20 in those cases. For large proteins, this is rarely a limiting factor since individual exons are typically longer than 30 residues.

---

## 10. Gene Prediction

### `gene_predictor`
**Type:** String | **Default:** `'auto'` | **Options:** `'auto'`, `'augustus'`, `'prodigal'`

Selects the ab initio gene prediction tool used to identify open reading frames in candidate syntenic regions of target genomes. **`auto`** (recommended) automatically selects based on the organism: Augustus for eukaryotes, Prodigal for prokaryotes. The selection is driven by the `home_species` kingdom detected during context building. **`augustus`** is a eukaryotic gene predictor that models intron-exon structure using species-specific hidden Markov models — it can predict multi-exon genes with realistic splice sites. **`prodigal`** is a prokaryotic gene finder optimized for compact genomes without introns — it is fast and accurate for bacteria and archaea but cannot predict spliced genes. Using Prodigal for eukaryotes will miss all multi-exon genes, which is most eukaryotic genes. Using Augustus for prokaryotes works but is unnecessarily slow. The `auto` setting handles this correctly in virtually all cases.

### `augustus_species`
**Type:** String | **Default:** `'fly'`

The species model used by Augustus for gene prediction. Augustus uses pre-trained statistical models that capture species-specific gene structure patterns (splice site signals, intron length distributions, codon usage, etc.). The model must match the target organism as closely as possible for accurate predictions. The default `'fly'` (Drosophila melanogaster) is a reasonable general-purpose model for insects and other arthropods. For vertebrate genomes, use `'human'` or `'chicken'`. For nematodes, use `'caenorhabditis'`. For fungi, use `'aspergillus_nidulans'` or another fungal model. Run `augustus --species=help` to see all available models. Using an inappropriate model (e.g., `'fly'` for a mammalian genome) causes Augustus to predict incorrect exon-intron boundaries, which can fragment gene models or merge adjacent genes. When using `--auto_params`, the LLM estimator does not currently adjust this parameter — manual selection based on your organisms is recommended.

### `pred_flank_window`
**Type:** Integer (bp) | **Default:** `50000`

The window size around each candidate locus for gene prediction. When Augustus or Prodigal is run on a candidate syntenic region, this parameter defines how much flanking sequence is included around the region for gene prediction context. Gene prediction tools perform better with surrounding context because they can detect gene boundaries more accurately when they see the transitions from intergenic to genic regions. The default of 50,000 bp provides ample context for most eukaryotic genes. For bacteria, this could be reduced to 10,000–20,000 bp. For very large eukaryotic genes that span hundreds of kilobases, increasing this ensures the entire gene plus some context is captured. Note: this parameter is defined in nextflow.config but its wiring to the gene prediction processes may be incomplete in the current pipeline version.

### `pred_keep_pct`
**Type:** Float | **Default:** `0.10` | **Range:** 0–1

The fraction of gene predictions to retain from Augustus/Prodigal output. Gene prediction tools often produce many candidate gene models, including low-confidence predictions, alternative transcripts, and fragments. This parameter filters to keep only the top percentage by prediction score. The default of 0.10 (10%) keeps only the highest-confidence predictions, reducing noise in downstream analysis. For organisms where gene prediction is less reliable (non-model species with limited training data), increasing this to 0.20–0.30 ensures more candidate gene models are available for GOI matching. Setting this too high floods the analysis with low-quality predictions; too low may discard the gene model that actually contains the GOI in a divergent form. Note: this parameter is defined in nextflow.config but its full implementation may be pending.

### `prodigal_full_genome_fallback`
**Type:** Boolean | **Default:** `false`

When enabled, Prodigal runs on the entire target genome if the windowed prediction around the candidate region fails to produce any gene models. The windowed approach (default) extracts just the candidate syntenic region and runs Prodigal on that fragment, which is fast but can fail if the extracted region is too small for Prodigal to train its self-calibrating model (Prodigal requires a minimum of ~20 kb of sequence for reliable self-training). Full-genome fallback re-runs Prodigal on the complete genome and then filters to predictions overlapping the candidate region. This is more reliable but much slower (processing an entire bacterial genome takes seconds, but a large eukaryotic genome could take minutes). The fallback is primarily useful for very small candidate regions in draft-quality genomes. For most analyses, the default of `false` is appropriate.

---

## 11. Gene Model Classification

### `classify_high_min_identity`
**Type:** Float (%) | **Default:** `50.0`

Minimum alignment identity percentage for an exon_annotation gene model to receive HIGH confidence. Lowered from 60 → 50 in 2026-04 because the previous default systematically mis-classified genuine cross-vertebrate orthologs (e.g. human TP53 vs fish tp53 at ~52% identity) as MEDIUM. At 50% identity the hit is comfortably above the twilight zone (~20–35%) and, combined with the required multi-exon evidence and flanking support, represents a high-confidence ortholog call across Metazoa. Raise toward 60 for closely related clades; lower toward 40 for very deep phylogenies.

### `classify_medium_min_identity`
**Type:** Float (%) | **Default:** `35.0`

Minimum alignment identity for MEDIUM confidence exon_annotation models. Lowered from 45 → 35 in 2026-04 to avoid collapsing divergent but real orthologs (e.g. TP53/TP63/TP73 across fish) into LOW. 35% is the lower edge of reliable homology detection; below this, signal is indistinguishable from random family-member noise. Models between 35 and 50 are MEDIUM, models below 35 drop to LOW. This threshold is fixed and not adjusted by the LLM estimator.

### `strict_goi_family`
**Type:** Bool | **Default:** `false`

Enables the family-consistency gate. When `true`, every GOI call with evidence_type in {`fallback_hit_span`, `rescued_exon`, `raw_hit`} that does **not** have a `TargetGene`/`TargetProduct` containing one of the expected family tokens is downgraded to LOW confidence / `ambiguous_goi_family_member`. Prevents fallback-heavy output from masquerading as probable GOI when the annotated locus is clearly a different gene (e.g. DNAH2 labelled as "probable TP53"). All GOI features also gain `GoiFamilyConsistent=true/false` and `GoiFamilyReason=...` attributes regardless of strict mode, which makes post-hoc filtering possible. Exon-annotation models (miniprot-supported multi-exon predictions) are never downgraded — their gene model is independent evidence.

### `goi_family_tokens`
**Type:** Comma-separated string | **Default:** `''` (auto)

Family name tokens used by `--strict_goi_family`. When empty, SynVoy parses the query FASTA header: UniProt `GN=XYZ` → token `XYZ`; UniProt entry-name `sp|ACC|NAME_SPECIES` → token `NAME` (e.g. `P04637` → `{TP53, P53}`). For multi-paralog queries where the run should accept all paralog labels (e.g. running TP53 but accepting TP63/TP73 annotations), override explicitly: `--goi_family_tokens TP53,TP63,TP73,TRP53,TRP63,TRP73`. Matching is case-insensitive substring after normalization (strip non-alphanumeric). `P53` correctly matches `Trp53` (mouse ortholog) because `P53` is a substring of `TRP53` once normalized.

### `classify_tandem_min_identity`
**Type:** Float (%) | **Default:** `40.0`

The minimum alignment identity for tandem copy gene models to receive MEDIUM confidence. Tandem copies (genes adjacent to the primary GOI hit that also match the query) use a slightly lower identity threshold than standard exon-annotation models because their genomic context (physical proximity to a confirmed ortholog) provides additional supporting evidence for their homology. A tandem copy at 40% identity next to a 70% identity primary hit is more likely to be a genuine paralog than an isolated 40% hit elsewhere in the genome. Below this threshold, tandem copies are classified as LOW confidence. This parameter is important for gene family analysis: when studying tandem gene arrays (e.g., olfactory receptors, defensins, MRJPs), the more divergent members of the array may fall just below standard MEDIUM thresholds but are still biologically meaningful given their genomic context.

### `classify_fragment_max_qcov`
**Type:** Float | **Default:** `0.4` | **Range:** 0–1

The maximum query coverage fraction below which a gene model is labeled as a `fragment` in the ModelStatus annotation field. Query coverage measures what fraction of the query protein length is covered by the aligned portion of the predicted gene. A model covering less than 40% of the query is considered a fragment — it likely represents an incomplete gene prediction (truncated by scaffold boundaries, missing exons, or pseudogenization). Fragment status is independent of the confidence classification (HIGH/MEDIUM/LOW): a fragment can be HIGH-confidence if the aligned portion has high identity, but it is still incomplete. This distinction is valuable for downstream analysis: fragments may represent pseudogenes, assembly artifacts, or genes at scaffold edges where the rest of the sequence is missing. Setting this higher (e.g., 0.5) classifies more models as fragments; lower (e.g., 0.3) is more lenient.

### `classify_complete_min_qcov`
**Type:** Float | **Default:** `0.7` | **Range:** 0–1

The minimum query coverage fraction for a gene model to be labeled as `complete` in the ModelStatus field. Models must cover at least 70% of the query protein and have multi-exon evidence (or be identified as tandem copies) to earn the `complete` status. Models between the fragment threshold (40%) and the complete threshold (70%) are labeled `partial` — they have significant alignment to the query but likely miss one or more exons or terminal regions. The `complete` status indicates a gene model that is likely to represent a full-length or near-full-length ortholog suitable for functional annotation and phylogenetic analysis without qualification. These thresholds are intentionally conservative: 70% coverage allows for some divergence at the N- or C-terminus (signal peptides, disordered tails) while ensuring the core protein is well-represented in the gene model.

---

## 12. Synteny Scoring Weights

### `synteny_weight_base`
**Type:** Float | **Default:** `0.4`

The weight assigned to the base synteny component in the composite synteny score. The composite score combines three components: base synteny (fraction of flanking genes found), gene-order consistency (whether the relative order of flanking genes is preserved), and strand conservation (whether flanking genes maintain their orientation). With a weight of 0.4, the base synteny fraction accounts for 40% of the final score. This is the largest single component because the mere presence of flanking gene orthologs in a target region is the strongest indicator of conserved synteny. This parameter is defined in nextflow.config but its implementation in the scoring pipeline may be incomplete in the current version. The three weights should ideally sum to 1.0 (0.4 + 0.3 + 0.3 = 1.0).

### `synteny_weight_consistency`
**Type:** Float | **Default:** `0.3`

The weight assigned to gene-order consistency in the composite synteny score. Gene-order consistency measures whether the flanking genes found in the target region maintain the same relative order as in the home genome (upstream gene A before upstream gene B, etc.). Conserved gene order provides stronger synteny evidence than mere gene presence, because gene order is more likely to be disrupted by chromosomal rearrangements over time. With a weight of 0.3, this component contributes 30% of the final score. Regions where all flanking genes are found but in scrambled order (indicating micro-inversions or rearrangements) receive a lower composite score than regions with fully conserved order. This helps distinguish true orthologous syntenic blocks from regions that happen to contain the same gene set due to convergent evolution or ancestral duplication followed by differential loss.

### `synteny_weight_strand`
**Type:** Float | **Default:** `0.3`

The weight assigned to strand conservation in the composite synteny score. Strand conservation checks whether flanking genes in the target region maintain the same transcriptional direction (sense/antisense) as in the home genome. Inversions that flip a segment of the genome reverse the strand of all genes within that segment. With a weight of 0.3, strand conservation accounts for 30% of the composite score. Perfect strand conservation (all flanking genes on the same relative strand) strongly supports synteny because large-scale inversions are relatively rare between closely related species. This component is particularly diagnostic for detecting inversions: a block with all expected genes in the correct order but reversed strand likely represents a simple inversion event, which still preserves synteny in a meaningful biological sense. The scoring accounts for whole-block inversions differently from gene-by-gene strand discordance.

### `synteny_goi_overlap_bonus`
**Type:** Float | **Default:** `0.15`

A bonus score added to candidate syntenic blocks that physically overlap an annotated GOI position. When a target genome has GFF annotations and the candidate block overlaps a gene annotated with the same name or ortholog group as the query GOI, this bonus is applied. The bonus rewards regions where independent evidence (gene annotation) corroborates the synteny-based prediction, increasing the score of true-positive blocks relative to false positives. The default of 0.15 (15% bonus) is moderate — enough to differentiate between otherwise equal candidates but not so large that it overrides strong synteny evidence. This is most useful in well-annotated genomes (model organisms) where GFF annotations are reliable. In draft genomes without GFF, this bonus never applies and scoring relies entirely on the three base components.

### `max_regions`
**Type:** Integer | **Default:** `0` (adaptive)

The maximum number of candidate regions emitted per locus for downstream evaluation. When set to 0 (default), the pipeline uses an adaptive strategy: it emits all regions that pass the synteny score threshold, with a hard cap of 6 per locus. This adaptive approach handles both single-copy genes (typically 1–2 regions per target) and multi-copy gene families (which may have legitimate multiple loci). Setting a specific number (e.g., `--max_regions 3`) forces the pipeline to emit exactly that many top-scoring regions per locus, regardless of how many pass the threshold. This is useful when you know your GOI is single-copy and want to minimize false positive processing. Setting it higher (e.g., 10) is useful for highly duplicated gene families. Regions are ranked by composite synteny score, so only the best N are emitted.

---

## 13. PLM (Protein Language Model) Search

### `enable_plm_search`
**Type:** Boolean | **Default:** `false`

Enables protein language model (PLM) embedding-based search using ProtT5-XL-UniRef50. When enabled, the pipeline computes a high-dimensional embedding vector for the query protein and for all candidate ORFs predicted within syntenic regions, then identifies ORFs whose embedding cosine similarity exceeds the threshold. PLM embeddings capture structural and functional information that goes beyond raw sequence alignment — two proteins can have completely different amino acid sequences but produce similar embeddings if they fold similarly and perform analogous functions. This makes PLM search extremely powerful for detecting remote homologs where sequence identity has dropped below 15–20%, well below the threshold of reliable sequence-based detection. The tradeoff is computational cost: computing ProtT5 embeddings requires significant CPU or GPU resources and the model itself requires ~5 GB of memory. Enable this for cross-phylum searches (>400 Mya divergence) where sequence methods alone may miss genuine orthologs.

### `plm_device`
**Type:** String | **Default:** `'cpu'` | **Options:** `'cpu'`, `'cuda'`

The device on which ProtT5 embeddings are computed. **`cpu`** works on any machine but is slow (minutes per protein for long sequences). **`cuda`** uses NVIDIA GPU acceleration and is 10–50x faster, making PLM search practical for large-scale analyses. GPU acceleration requires PyTorch with CUDA support and an NVIDIA GPU with at least 6 GB VRAM (the ProtT5-XL model requires ~5 GB). If `cuda` is specified but no GPU is available, the pipeline will fail with a CUDA error. For batch processing of many target genomes, GPU acceleration is strongly recommended when PLM search is enabled. On CPU, consider enabling PLM search only for specific high-priority distant targets rather than the entire target set.

### `plm_similarity_threshold`
**Type:** Float | **Default:** `0.5` | **Range:** 0–1

The minimum cosine similarity between the query protein embedding and a candidate ORF embedding for the ORF to be reported as a PLM hit and considered as a potential GOI ortholog. Cosine similarity of 1.0 means identical embedding vectors (essentially the same protein); 0.0 means orthogonal (completely unrelated). The default of 0.5 is a moderate threshold that captures proteins with significant structural/functional similarity while excluding clearly unrelated proteins. For highly conserved gene families, raise this to 0.6–0.7 to reduce false positives. For extremely divergent searches where you want maximum sensitivity, lower to 0.3–0.4 at the cost of more false positives. Unlike sequence identity, there is no well-established "twilight zone" for PLM similarity — the threshold should be calibrated empirically for your specific gene family and divergence level.

**Short-query auto-adjustment:** when the longest GOI query is under 100 aa, the pipeline automatically raises this threshold to at least `0.75` per-genome. Mean-pooled ProtT5 embeddings are noisy for short peptides (few residues dominate the pooled vector), so the default of 0.5 lies below the noise floor and produces spurious hits for queries like melittin, defensins, and small venom peptides. The auto-raise is conservative and can be overridden by explicitly setting a higher `plm_similarity_threshold` on the command line.

### `plm_medium_threshold`
**Type:** Float | **Default:** `0.7` | **Range:** 0–1

The PLM cosine similarity threshold above which a LOW-confidence gene model can be boosted to MEDIUM confidence. This integrates PLM evidence into the confidence classification system: a gene model that has poor sequence-level evidence (LOW confidence based on identity/coverage) but strong embedding similarity (>0.7) to the query suggests that the protein is a genuine homolog whose sequence has diverged beyond reliable alignment-based detection. The boost recognizes that PLM embeddings provide orthogonal evidence to sequence alignment — high embedding similarity despite low sequence identity is a hallmark of remote homology. The default of 0.7 requires strong embedding similarity for the boost, which is appropriate because boosting a LOW model to MEDIUM has significant implications for the interpretation of results.

### `plm_high_threshold`
**Type:** Float | **Default:** `0.85` | **Range:** 0–1

The PLM cosine similarity threshold above which a MEDIUM-confidence gene model can be boosted to HIGH confidence. An embedding similarity of 0.85 indicates near-structural equivalence, providing very strong evidence for homology even when sequence-based evidence is only moderate. This is the most conservative PLM boost tier — only clearly related proteins at the structural/functional level receive the HIGH designation based on PLM evidence. This threshold should rarely be lowered because HIGH confidence is used for making definitive orthology claims, and PLM embeddings, while powerful, are not infallible. The gap between the medium threshold (0.7) and the high threshold (0.85) creates a meaningful distinction between "probably homologous" (MEDIUM boost) and "almost certainly homologous" (HIGH boost) based on PLM evidence.

---

## 14. Structural Search (ESMFold + Foldseek)

### `enable_structural_search`
**Type:** Boolean | **Default:** `false`

Enables 3D structure prediction and structural comparison for GOI candidate detection. When enabled, the pipeline predicts the 3D structure of the query protein and all candidate ORFs using ESMFold (Meta's efficient protein structure predictor), then compares structures using Foldseek's 3Di structural alphabet. Structural comparison is the ultimate method for detecting remote homologs: protein 3D structure is conserved approximately 10x longer than protein sequence during evolution. Proteins that share <15% sequence identity can still have near-identical folds. The tradeoff is extreme computational cost — ESMFold structure prediction requires significant GPU memory (4–16 GB depending on protein length) and takes seconds to minutes per protein. Enable this only for very distant cross-kingdom searches (>600 Mya) where both sequence alignment and PLM embeddings may fail. Requires both PyTorch and Foldseek to be installed.

**Short-query auto-skip:** when the longest GOI query is under 50 aa, structural search is automatically skipped for that genome (a one-line log message is emitted). ESMFold produces degenerate structures for tiny peptides and TM-score is not meaningful for a single short amphipathic helix (e.g. melittin, defensins). Keeping the feature enabled for a mixed run is safe — only the sub-50-aa cases are bypassed.

### `structural_device`
**Type:** String | **Default:** `'cpu'` | **Options:** `'cpu'`, `'cuda'`

The device for ESMFold structure prediction. **GPU (`cuda`) is strongly recommended** for structural search — ESMFold on CPU is impractically slow (minutes to hours per protein versus seconds on GPU). The ESMFold model requires approximately 4 GB of GPU VRAM for small proteins (<200 aa) and up to 16 GB for large proteins (~700 aa, the maximum set by `structural_max_length`). If you don't have a suitable GPU, consider disabling structural search entirely and relying on PLM embeddings (which are feasible on CPU). CUDA support requires PyTorch compiled with CUDA and an NVIDIA GPU. AMD GPUs via ROCm may work with the appropriate PyTorch build but are untested. The Foldseek comparison step itself is CPU-based and runs quickly regardless of this setting.

### `structural_tm_threshold`
**Type:** Float | **Default:** `0.3` | **Range:** 0–1

The minimum TM-score (Template Modeling score) from Foldseek structural comparison for a candidate ORF to be reported as a structural hit. TM-score ranges from 0 to 1, where 1.0 is an identical structure. TM-scores above 0.5 generally indicate the same fold; above 0.7 indicates the same protein family. The default of 0.3 is permissive, allowing detection of proteins that share a common structural core but may have significant structural differences in peripheral regions. At 0.3, proteins typically share the same superfold but may belong to different protein families. This low threshold is appropriate because structural search is designed for extreme divergence cases where any structural similarity is biologically meaningful. For more conservative structural matching (only proteins with clearly similar folds), raise to 0.5. The TM-score is length-independent, making it more reliable than RMSD for comparing proteins of different sizes.

### `structural_medium_threshold`
**Type:** Float | **Default:** `0.5` | **Range:** 0–1

The TM-score threshold above which structural evidence can boost a LOW-confidence gene model to MEDIUM confidence. A TM-score of 0.5 indicates that the candidate shares the same overall fold as the query — this is strong structural evidence for homology even when sequence-based evidence is weak. Structural homology at this level is rarely coincidental: convergent evolution to the same fold does occur but is relatively rare, and the syntenic context further reduces the chance of false positives. This boost is particularly valuable for ancient gene families where sequence has diverged beyond recognition but the 3D fold is preserved. The combination of synteny evidence plus structural similarity at TM>0.5 provides high-confidence orthology assignment even for cross-kingdom comparisons.

### `structural_high_threshold`
**Type:** Float | **Default:** `0.7` | **Range:** 0–1

The TM-score threshold above which structural evidence can boost a MEDIUM-confidence gene model to HIGH confidence. A TM-score of 0.7 indicates strong structural similarity — the proteins are almost certainly in the same protein family with a shared evolutionary origin. At this level, structural homology provides evidence comparable to moderate sequence identity (30–40%), justifying the HIGH confidence designation. Combined with synteny context, a TM-score >0.7 hit is a very reliable ortholog identification even in the complete absence of detectable sequence similarity. This is the state-of-the-art for detecting orthologs separated by >500 million years of evolution, where sequence methods have essentially zero sensitivity.

### `structural_max_length`
**Type:** Integer (aa) | **Default:** `700` | **Range:** 50–2000

The maximum protein length in residues for ESMFold structure prediction. Proteins longer than this threshold are skipped during structural search because ESMFold's memory requirements scale quadratically with sequence length — a 700-residue protein requires approximately 8–12 GB GPU VRAM, while a 1000-residue protein would require 20+ GB. The default of 700 provides a practical ceiling that works with most consumer and workstation GPUs (RTX 3090 with 24 GB, A100 with 40/80 GB). On GPUs with limited VRAM (6–8 GB), reduce to 400–500. On high-memory GPUs (A100 80 GB), increase to 1000+. Query proteins longer than this threshold still get structure-based analysis for their individual domains if domain boundaries are detected, but the full-length structural comparison is skipped.

**Auto VRAM guard:** when `structural_device = cuda`, the pipeline probes the GPU at fold time. If the device has less than 20 GB total VRAM and `structural_max_length > 400`, the effective length is automatically capped to 400 aa to avoid mid-fold OOM, and a warning is logged. Explicit user settings below 400 are respected as-is; the guard only lowers values, never raises them.

---

## 15. Visualization

### `plot_width`
**Type:** Integer (pixels) | **Default:** `1500`

The width of the output SVG synteny plot in pixels. The synteny plot is the primary visual output of SynVoy, showing gene arrows along chromosomal tracks with homology links between genomes and a phylogenetic tree alongside. The default of 1500 pixels provides a good balance between detail and readability on standard monitors. For presentations or publications, increase to 2000–3000 for higher resolution. For quick previews on small screens, reduce to 800–1000. The plot is rendered as an SVG (Scalable Vector Graphics), so it can be scaled without quality loss regardless of this setting — the width primarily controls the initial rendering layout and how much genomic context fits horizontally. Wider plots can show more flanking genes and larger genomic regions without overlapping labels.

### `gap_threshold`
**Type:** Integer (bp) | **Default:** `50000`

The minimum gap size in base pairs in the genomic sequence that triggers visual compression in the synteny plot. When the distance between adjacent genes in the plot exceeds this threshold, the gap is visually compressed to `gap_visual_size` to prevent large empty spaces from dominating the plot and making nearby genes too small to see. The default of 50,000 bp works well for most eukaryotic genomes where intergenic regions of 50+ kb are common but not the focus of the visualization. For compact genomes (bacteria), lower this to 5,000–10,000 bp. For large vertebrate genomes where 50 kb gaps are routine, increase to 100,000–200,000 bp. This parameter only affects visualization — it does not change any analytical results. The compression is indicated visually in the plot with a special gap symbol.

### `gap_visual_size`
**Type:** Integer (bp) | **Default:** `20000`

The visual size in base pairs used to represent compressed gaps in the synteny plot. When a genomic gap exceeds `gap_threshold`, it is drawn as if it were `gap_visual_size` base pairs long, with a visual indicator showing that the gap has been compressed. The default of 20,000 bp provides enough visual space to indicate a gap while keeping it compact relative to the gene arrows. The ratio between `gap_threshold` and `gap_visual_size` determines the compression factor — the default ratio of 50,000:20,000 means gaps are compressed to about 40% of the threshold size. Decrease this value for more aggressive compression (smaller gaps in the plot); increase it if you want compressed gaps to remain more visible. This is purely a cosmetic parameter with no analytical impact.

### `flank_fallback_bp`
**Type:** Integer (bp) | **Default:** `1000000`

The maximum genomic window in base pairs rendered around distal target hits when the normal flanking gene context is not available. In some cases, a GOI hit in a target genome may not have a well-defined syntenic block (e.g., only the GOI was found with minimal flanking gene support). In these cases, the plot renderer falls back to showing a fixed window of `flank_fallback_bp` around the hit, showing whatever genes are annotated in that window. The default of 1,000,000 bp (1 Mb) provides substantial genomic context. For compact genomes, this shows many genes; for large genomes, it may show only a few. This is primarily a visualization safety net — the vast majority of hits have proper syntenic block context. Reduce to 200,000–500,000 for more focused views; increase if you want maximum genomic context around poorly-anchored hits.

### `scale_bar_len`
**Type:** Integer (bp) | **Default:** `10000`

The length of the scale bar in the synteny plot, in base pairs. The scale bar provides a visual reference for interpreting distances in the plot. The default of 10,000 bp (10 kb) is appropriate for typical eukaryotic synteny plots where gene-to-gene distances are in the range of 1–100 kb. For bacterial genome plots where genes are much more tightly packed (1–5 kb intergenic distances), a 1,000–2,000 bp scale bar is more appropriate. For large vertebrate genome plots where distances span megabases, a 50,000–100,000 bp scale bar provides better visual reference. The scale bar is drawn at a fixed position in the plot and labeled with the length value. This is a purely cosmetic parameter that helps interpret the spatial relationships in the visualization.

### `hide_goi_absent_tracks`
**Type:** Boolean | **Default:** `true`

When enabled, genomic tracks (species rows) where no GOI candidate was found are hidden from the synteny plot. This keeps the visualization focused on species where the GOI was detected, avoiding visual clutter from empty tracks that provide no useful information. With the default of `true`, only species with at least one GOI hit appear in the plot. Setting this to `false` shows all target species including those where the GOI was not found — this can be informative for understanding the phylogenetic distribution of gene loss events (where did the GOI disappear?). For publications where you want to explicitly show which species lack the gene, set this to `false`. For routine analysis where you care primarily about the species where orthologs were found, the default `true` produces cleaner, more interpretable plots.

---

## 16. Resource Tuning

### `iterative_search_cpus`
**Type:** Integer | **Default:** `2`

The number of CPU cores allocated to each ITERATIVE_SEARCH process. This is the most computationally intensive task in the pipeline — it runs MMseqs2, tblastn, miniprot, and Smith-Waterman alignment for each target genome. More CPUs accelerate the MMseqs2 prefiltering and alignment stages, which are the bottleneck. The default of 2 is conservative, designed for shared workstations. On dedicated machines, increase to 4–8 for significantly faster processing. The `docker_max` profile automatically detects and uses all available CPUs. On HPC clusters, this should match your SLURM job CPU allocation. Note that giving one ITERATIVE_SEARCH task all CPUs while running multiple tasks in parallel (`iterative_search_max_forks > 1`) can cause oversubscription — the product of CPUs × forks should not exceed your total core count. Memory pressure typically limits parallelism before CPU count.

### `iterative_search_memory`
**Type:** String | **Default:** `'10 GB'`

The memory allocated to each ITERATIVE_SEARCH process. MMseqs2 database indexing and search are the primary memory consumers. The required amount depends on target genome size: bacterial genomes need 2–4 GB, insect genomes need 4–8 GB, and large vertebrate genomes can require 16–32 GB. The default of 10 GB works for most invertebrate and small vertebrate genomes. If processes fail with exit code 137 (OOM kill), increase this value. The `mmseqs_split_memory_limit` parameter provides a secondary control — if the database exceeds that limit, MMseqs2 splits it and processes chunks sequentially, reducing peak memory at the cost of speed. On memory-constrained systems, it is better to reduce `mmseqs_split_memory_limit` than to lower this value, because the process needs baseline memory for the main workflow beyond just MMseqs2.

### `iterative_search_max_forks`
**Type:** Integer | **Default:** `1`

The maximum number of ITERATIVE_SEARCH processes that can run simultaneously. Each fork processes a different target genome. The default of 1 (sequential processing) is the safest option — it prevents memory oversubscription and is appropriate for laptops and workstations. On machines with ample RAM (64+ GB) and many CPU cores, increase to 2–4 to parallelize across target genomes. The memory requirement multiplies with fork count: if each ITERATIVE_SEARCH needs 10 GB and you set max_forks to 4, you need at least 40 GB of available RAM. The `laptop_safe` profile keeps this at 1 to prevent system freezes. On HPC clusters, this is less relevant because each target genome is typically submitted as a separate SLURM job. Setting this too high on limited hardware is the most common cause of pipeline crashes due to OOM kills.

### `iterative_quiet_subtools`
**Type:** Boolean | **Default:** `true`

When enabled, suppresses verbose output from sub-tools (MMseqs2, tblastn, miniprot) called within the ITERATIVE_SEARCH process. This reduces log file size and I/O overhead, which can be significant when processing many target genomes with high MMseqs2 verbosity. The default of `true` keeps logs clean and focused on pipeline-level progress messages. Set to `false` when debugging a specific search failure — the full tool output can help diagnose why a particular target genome produced no hits or crashed. The `laptop_safe` profile enables this to minimize disk I/O during single-threaded execution. This parameter does not affect search behavior or results, only the amount of diagnostic output captured in the `.command.log` files within the Nextflow work directory.

### `locate_gene_cpus`
**Type:** Integer | **Default:** `1`

The number of CPU cores allocated to the LOCATE_GENE process, which finds the GOI in the home genome during the initial localization phase. This step runs tblastn and MMseqs2 against a single genome (the home genome), so it is typically less resource-intensive than ITERATIVE_SEARCH. The default of 1 CPU is usually sufficient because the home genome is a single genome searched once. Increasing to 2–4 CPUs can speed up the MMseqs2 search for large genomes (vertebrates with 1000+ Mb). This task is not a bottleneck in most pipeline runs — it completes quickly and only runs once. The `docker_max` profile increases this to use all available CPUs for maximum speed on dedicated machines.

### `locate_gene_memory`
**Type:** String | **Default:** `'3 GB'`

The memory allocated to the LOCATE_GENE process. This should accommodate the MMseqs2 database for the home genome plus working memory for tblastn. The default of 3 GB is sufficient for most genomes up to ~500 Mb. For large vertebrate genomes (1–3 GB), increase to 8–16 GB. If this process fails with exit code 137, it is running out of memory — increase this value or reduce `mmseqs_split_memory_limit` to force MMseqs2 to split the database. The `laptop_safe` profile sets this to 2 GB, which works for smaller genomes but may need manual increase for large vertebrate genomes. As with `iterative_search_memory`, the `mmseqs_split_memory_limit` provides a secondary control mechanism.

---

## 17. LLM Parameter Estimation

### `auto_params`
**Type:** Boolean | **Default:** `true`

Master switch for automatic parameter estimation. When enabled, SynVoy analyzes the biological context of your search — the query protein characteristics (size, exon count, signal peptide, gene family), the home species genome architecture (kingdom, genome size, gene count, intron lengths), and the target species evolutionary distances — to automatically set optimal values for approximately 25 search parameters. The estimation uses a three-tier system: (1) Ollama local LLM (Gemma 4), (2) Google Cloud Gemini API, (3) deterministic heuristic rules. All three backends encode the same biological reasoning (kingdom-specific intron sizes, distance-adaptive sensitivity, query-size thresholds), but the LLM backends add nuance for edge cases. Disable this with `--auto_params false` if you want full manual control over all parameters or if the estimation is producing suboptimal results for your specific use case. When disabled, all parameters use their nextflow.config defaults.

### `llm_model`
**Type:** String | **Default:** `'auto'` | **Options:** `'auto'`, `'gemma4:e4b'`, `'gemma4:26b'`, `'gemma4:31b'`

The Ollama model used for local LLM-powered parameter estimation. **`auto`** detects system resources (RAM, GPU VRAM) and selects the best model that will fit: `gemma4:e4b` (4B parameters, ~4 GB, laptops), `gemma4:26b` (26B MoE, ~16 GB, workstations), or `gemma4:31b` (31B dense, ~24 GB, servers). The larger models produce more nuanced parameter estimates, especially for unusual biological contexts (polyploid plants, highly rearranged insect lineages, extremophile bacteria). The 4B model is sufficient for common cases (standard insects, vertebrates, plants) but may struggle with edge cases. All models use the same system prompt and produce JSON output. If Ollama is not installed or the model is not pulled, the pipeline seamlessly falls back to the Google Cloud API or heuristic estimation. Specify a model explicitly to force a particular size regardless of detected resources.

### `ollama_url`
**Type:** String | **Default:** `'http://localhost:11434'`

The URL of the Ollama server for local LLM inference. The default assumes Ollama is running locally on the standard port. Change this if Ollama is running on a different machine (e.g., `'http://gpu-server:11434'`) or a non-standard port. The pipeline verifies server reachability with a version check before attempting model inference. If the server is unreachable (network error, Ollama not running), the pipeline falls back to Google Cloud or heuristic estimation without error — the Ollama backend is always optional. For containerized deployments, you may need to use the Docker host IP or a service name instead of localhost. The pipeline sends a single chat API request per run, so network latency is not a concern — even a remote server on a different continent would add only a few hundred milliseconds to the total estimation time.

### `ollama_timeout`
**Type:** Integer (seconds) | **Default:** `480`

The maximum time in seconds to wait for a response from the Ollama server. CPU inference with Gemma 4 models can be slow — the 4B model (`gemma4:e4b`) takes approximately 3–4 minutes on a modern laptop CPU, while the 26B and 31B models can take 10+ minutes without GPU acceleration. The default of 480 seconds (8 minutes) provides adequate headroom for CPU inference of the smallest model. If using GPU inference, responses typically return in 5–30 seconds, making the timeout effectively irrelevant. If the timeout is hit, the pipeline falls back to Google Cloud or heuristic estimation. Increase this value if running larger models on CPU. Decrease it if you want faster fallback to heuristic when Ollama is slow. The timeout applies to the total API response time including model loading, inference, and response generation.

### `google_api_key`
**Type:** String | **Default:** `''` (empty)

The Google Cloud Gemini API key for cloud-based LLM parameter estimation. This is the second-tier fallback after Ollama — if Ollama is not available or fails, and this key is provided, the pipeline sends the biological context to the Gemini 2.5 Flash Lite model via Google's generative AI API. The key can be provided via this parameter (`--google_api_key YOUR_KEY`), or more safely via the `GOOGLE_API_KEY` environment variable (recommended to avoid exposing the key in command line history or Nextflow logs). The API key is associated with your Google Cloud project billing — each call uses approximately 1,000–2,000 tokens, which is negligible cost even on the free tier. The Gemini API typically responds in 2–3 seconds. If the API returns an error (rate limiting, authentication failure), the pipeline falls back to heuristic estimation. **Never commit API keys to version control.**

### `multi_profile`
**Type:** Boolean | **Default:** `true`

Enables multi-profile mode, where SynVoy runs the same search with multiple parameter profiles (sensitive, balanced, stringent) and automatically selects the best result. When the LLM estimates parameters, it produces a single "balanced" profile. Multi-profile mode generates two additional variants: a "sensitive" profile with relaxed thresholds (lower identity, lower synteny score, higher sensitivity) and a "stringent" profile with tighter thresholds (higher identity, higher synteny score). The best result per target is selected based on GOI detection confidence and synteny score. This triples the computational work but significantly improves the chance of finding the optimal parameters for each target species — close species benefit from stringent parameters while distant species benefit from sensitive parameters. The mode is automatically disabled when the total job count (loci × targets × 3 profiles) would exceed `multi_profile_max_jobs`.

### `multi_profile_max_jobs`
**Type:** Integer | **Default:** `30`

The maximum total number of jobs (loci × targets × 3 profiles) allowed before multi-profile mode is automatically disabled. When the estimated job count exceeds this threshold, only the LLM-estimated balanced profile runs, saving computational resources. The default of 30 accommodates typical small-to-medium analyses (e.g., 1 locus × 10 targets × 3 profiles = 30 jobs). For larger analyses (many target genomes or multiple loci), multi-profile would create excessive parallelism — a 3-locus × 20-target analysis would need 180 profile runs, which is likely not worth the 3x computational cost. Increase this value on HPC clusters where computational resources are abundant. Decrease it on laptops or when optimizing for runtime. When multi-profile is disabled (either manually or by exceeding this cap), only the single best-estimated parameter set runs.

---

## 18. Advanced & Output

### `outdir`
**Type:** Path (String) | **Default:** `'results'`

The directory where all pipeline output files are written. This includes synteny plots (SVG/HTML), phylogenetic trees (Newick), BED region files, the JSON summary report, and optionally intermediate files. The directory is created automatically if it does not exist. Using descriptive directory names (e.g., `results/melittin_apis_vs_bombus`) helps organize multiple runs. When using `-resume`, Nextflow caches task results in the `work/` directory and regenerates outputs in `outdir` — so changing `outdir` and resuming will write results to the new location without re-running completed tasks. For HPC runs, consider pointing this to a permanent storage location rather than scratch space that may be automatically purged.

### `keep_intermediate`
**Type:** Boolean | **Default:** `false`

When enabled, intermediate files from each pipeline stage are preserved in the output directory. This includes per-target MMseqs2 hit tables, tblastn alignments, miniprot output, Augustus/Prodigal predictions, Smith-Waterman alignment results, flanking gene FASTAs, and per-target GFF annotations. These files are invaluable for debugging unexpected results — for example, examining the raw tblastn output for a target where no GOI was found, or inspecting the flanking gene sequences to understand why synteny scores are low. The default of `false` keeps only the final output files (plots, trees, reports) to minimize disk usage. Intermediate files can be large (hundreds of MB for multi-genome analyses), so enable this judiciously. Note that Nextflow's `work/` directory always contains intermediate files for cached tasks, regardless of this setting — `keep_intermediate` controls what is copied to `outdir`.

### `max_retries`
**Type:** Integer | **Default:** `3`

The maximum number of times a failed Nextflow process is retried before the pipeline aborts. Combined with the `errorStrategy` configuration (which retries on signal-based failures: SIGKILL, SIGTERM, SIGSEGV, etc.), this provides resilience against transient failures like OOM kills, network timeouts, and temporary file system errors. When a process fails with an OOM kill (exit code 137), Nextflow retries it — if the process has a dynamic memory directive that scales with retry count, the retry may succeed with more memory. The default of 3 retries provides robust fault tolerance without infinite loops on persistent failures. For HPC environments with occasional node failures, 3 retries is appropriate. For local runs where failures are more likely systematic (wrong parameters) than transient, you might reduce to 1–2 to fail faster. This parameter applies globally to all processes.

### `docker_container`
**Type:** String | **Default:** `'synvoy-local:latest'`

The Docker/Singularity container image used for running pipeline processes when a container profile is active. The default expects a locally built image tagged `synvoy-local:latest`. Build it from the project Dockerfile: `docker build -t synvoy-local:latest .`. The container bundles all dependencies (MMseqs2, tblastn, Augustus, Prodigal, miniprot, parasail, MAFFT, IQ-TREE, Python with BioPython, etc.) in a reproducible environment. Override this to use a pre-built remote image (e.g., from Docker Hub or a private registry). The `beforeScript` directive ensures workspace scripts in `bin/` override any scripts packaged inside the container, so local code changes take effect immediately without rebuilding. This parameter is only used when a Docker or Singularity profile is active — Conda profiles ignore it entirely.

---

## 19. Reserved / Unimplemented

The following parameters are defined in `nextflow.config` but are not currently wired into the pipeline. They are reserved for planned future features:

| Parameter | Default | Intended Purpose |
|---|---|---|
| `expand_db_threshold` | `1e-10` | E-value threshold for expanding the iterative search database with discovered homologs |
| `diamond_sensitivity` | `"very-sensitive"` | DIAMOND alignment sensitivity (planned DIAMOND integration as MMseqs2 alternative) |
| `enable_splice_variants` | `true` | Detect and report alternative splice variants at each locus |
| `enable_frameshifts` | `true` | Detect and report frameshift mutations in pseudogene candidates |
| `mutation_rate` | `0.05` | Expected mutation rate for evolutionary distance calibration |
| `num_mutant_variants` | `10` | Number of mutant variants to generate for sensitivity testing |

These parameters can be safely ignored. Setting them has no effect on current pipeline behavior.
