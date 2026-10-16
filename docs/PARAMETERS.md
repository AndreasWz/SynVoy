# SynVoy — Comprehensive Parameter Reference

This document provides in-depth explanations for every configurable parameter in SynVoy. Each entry describes the parameter's purpose, its effect on the pipeline, practical guidance for tuning, and the biological reasoning behind its default value. Parameters are grouped by functional category.

---

## Quick start: pick a preset

Most users don't need to tune individual parameters, and since 2026-10-09 not the distances either: how far apart two anchors may lie, how much DNA around a block is searched and similar distances are derived for each home locus and each target genome ([§ 5, Distances in gene units](#distances-in-gene-units)).

SynVoy ships four named presets. Each holds only the values that have a reason for that kind of query; everything else is the default.

| Preset | Chosen automatically when | What it changes |
|---|---|---|
| `preset_short_peptide` | the query is shorter than 100 aa | the coverage bands of the model status, 0.5 / 0.3 (one missed exon is a large share of a short protein) |
| `preset_paralog_discrimination` | the query aligns well at two or more places of the home genome, the second at least half as strong as the first | fewer candidate regions per genome in the figures (`adaptive_max_regions 4`, stricter score floors) |
| `preset_tandem_family` | the home locus itself holds more than one copy of the gene | a neighbour counts as a relative of the gene of interest, and is not used as an anchor, up to a pairwise E-value of 1e-3 instead of 1e-5 |
| `preset_single_copy` | none of the above | nothing |

No preset moves the label of a GOI call (since 2026-10-15). One judge sets HIGH / MEDIUM / AMBIGUOUS / LOW for every call, from the E-value in the expected place, specificity against the home proteome and a completeness bar of 0.7 query coverage ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)); its thresholds are significance levels, the same for every query. The identity thresholds that the presets used to move (HIGH 40 / MEDIUM 30 / tandem 40 for short peptides; `strict_goi_family` in three presets) were removed with the identity-based label ([§ 20](#20-removed-parameters)).

Until 2026-10-09 `preset_short_peptide` also lowered four hit floors (`sw_min_score 10`, `sw_min_identity 8`, `min_hit_identity 8`, `min_hit_length 8`). They are no longer part of it, because lowering them changes nothing: the E-value filters of the search are stricter. A Smith-Waterman hit must reach E ≤ 1, which needs a raw score of 31 even for a 20-aa query in a 2 kb window (the floor is 20); inside a window the identity and length floors are `max(15, …)` for either value; and no alignment of 8–9 residues or 8–10 % identity reaches E ≤ 0.01 in a genome-wide search.

"Aligns well" is counted in places, not in alignment rows: each exon and each search program adds rows, so a row count says little about a family (until 2026-10-08 the chooser used `n_hits > 15` and half of all queries, most with one strong locus, were given `preset_paralog_discrimination`). A place is scored as one copy of the gene (`one_copy_score` in `bin/sequence_utils.py`).

To pin one by hand: `-profile <execution>,preset_short_peptide` (at launch) or `--preset_override preset_short_peptide`.

`<execution>` is one of `auto`, `standard`, `conda`, `docker`, `docker_max`, `singularity`,
`slurm`, `hpc_singularity`, `hpc_conda`, `lrz_ai`, `lrz_ai_container` — see
[USAGE.md § 2](USAGE.md#2-execution-profiles) for the full table. A memory tier
(`laptop_safe` / `low_mem`) is a *separate* axis you may add on top; presets compose with
both via comma, e.g. `-profile auto,low_mem,preset_short_peptide`.

Equivalent invocation with `-c` (no profile composition):
```
./run_synvoy.sh -c conf/presets/<preset_name>.config ...
```

Preset source files live in [`conf/presets/`](../conf/presets/) — copy and edit them as a starting point for query-specific tweaks. The full reference below documents every individual parameter for users who need finer control.

### Preset auto-apply (default ON since 2026-05-29)

You don't actually have to pick a preset by hand. After `LOCATE_GENE` finishes, SynVoy profiles the BLAST/MMseqs2 hit distribution (`intermediate/locate_gene/hit_profile.json`), picks the matching preset (`intermediate/locate_gene/auto_preset.json`), and **applies it at runtime** to the downstream steps (`EXTRACT_FLANKING`, `ITERATIVE_SEARCH`, `CLUSTER_REGIONS`) via a value channel — no pipeline restart needed. The chosen preset is logged as:

```
[INFO] §1f preset applied: preset_paralog_discrimination (auto_preset.json chose preset_paralog_discrimination)
```

Knobs:

- **A value you set yourself wins over the preset** (since 2026-10-05). Any preset-covered parameter whose launch-time value differs from the default shipped in `nextflow.config` — whether you set it on the command line, in a params file, a `-c` config or through a `-profile` — is kept, and `intermediate/locate_gene/effective_params.json` records it (`source: user_set`, and `user_set_kept` when the preset wanted a different value). Before that date the preset silently overwrote it. The one case that cannot be detected is explicitly passing the *shipped default* (e.g. `--goi_relative_evalue 1e-5`): to pin a default against a preset, use the next option.
- `--auto_apply_preset false` keeps the recommendation advisory-only: no preset value is applied at all.
- `--preset_override preset_X` pins a specific preset, bypassing the hit-profile decision. Combine with `--auto_apply_preset true` (the default) to apply it at runtime, or use `-profile preset_X` to apply it at launch.

The minimum query length is not part of any preset: it is 20 aa for every run (`min_query_length`), and a query of 20–39 aa gets a warning. The preset-covered parameters are the Smith–Waterman and hit floors, the model-status coverage bounds, `expand_goi_similar`, the family names, the adaptive-region limits and `goi_relative_evalue`.

The Python `PRESET_OVERRIDES` table in `bin/resolve_effective_params.py` is the single source of truth; `tests/test_resolve_effective_params.py::test_preset_overrides_match_groovy_configs` keeps the Groovy `.config` files in lockstep so `-profile preset_X` and `--auto_apply_preset` produce identical effective values.

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
19. [Auto-Apply Preset, Self-Consistency & Rescue](#19-auto-apply-preset-self-consistency--rescue)
20. [Removed parameters](#20-removed-parameters)

> **Not everything is in this file.** Several newer parameter groups are documented as
> tables in [USAGE.md § 4](USAGE.md#4-full-parameter-reference) rather than as long-form
> entries here: the multi-locus cap (`--max_loci`), the rescue passes (strong-synteny, GOI hull, dispersed-GOI),
> locus ownership and paralog discrimination, distance auto-tuning and structural
> discovery.
> Between the two documents every parameter in `nextflow.config` is covered.

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
**Type:** String | **Default:** `'auto'` | **Options:** `'auto'`, `'drop'`, `'keep'`, `'ask'`

What happens to an assembly that fails the quality limits (`bad_min_n50`, `bad_max_contigs`, `bad_max_scaffolds`) when genomes are fetched in easy mode. Chromosome-level and complete assemblies always pass.

- **`auto`** (default) depends on who chose the genome. An assembly the pipeline picked on its own, while walking up the taxonomy, is dropped: another species takes its place. An assembly the user asked for is kept and flagged: a species named in `--target_species`, and the home assembly, which follows from the query. For these there is nothing to fall back to, and leaving the genome out silently would read as "no gene found there".
- **`drop`** removes every low-quality assembly, also a named species. A low-quality home assembly then stops the run.
- **`keep`** keeps every one.
- **`ask`** prompts on the terminal for each (see `bad_quality_timeout`).

A kept low-quality genome is marked in the console, in `synvoy_report.json` (`warnings`) and in the figures (a dagger after the species name). A gene that is not found in it may be missing from the assembly, not from the species.

### `bad_quality_timeout`
**Type:** Integer (seconds) | **Default:** `300`

The number of seconds SynVoy waits for a user response when `bad_quality_policy` is set to `'ask'`. During this window, the pipeline pauses and displays the assembly name, contig count, scaffold count, and N50 value, then asks whether to keep or drop the assembly. If the timeout expires without user input, the assembly is dropped (fail-safe behavior). This parameter only has an effect when `bad_quality_policy = 'ask'`; it is ignored for `'drop'` and `'keep'` policies. A 5-minute default gives users enough time to evaluate the metrics, but you can increase it for batch runs where you might not be watching the terminal continuously. Setting it very high (e.g., 86400) effectively makes the pipeline wait indefinitely for input.

### `bad_max_contigs`
**Type:** Integer | **Default:** `500000`

Applies only when the contig count of an assembly is known. NCBI's assembly records carry N50 values but no counts, so for most assemblies the decision rests on the assembly level (chromosome-level and complete assemblies always pass) and on `bad_min_n50`. The same holds for `bad_max_scaffolds`.

The maximum number of contigs an assembly can have before it is flagged as low quality. Assemblies with more contigs than this threshold are subject to the `bad_quality_policy`. Contig count is a proxy for assembly fragmentation — a highly fragmented assembly has hundreds of thousands of small contigs, meaning genes and their flanking regions are often split across multiple pieces. This breaks synteny analysis because the pipeline cannot reliably identify the genomic neighborhood of a gene when that neighborhood spans multiple unlinked contigs. The default of 500,000 is permissive; most reasonable genome assemblies have far fewer contigs. Lowering this threshold (e.g., to 50,000) enforces stricter quality requirements, which improves synteny reliability but may exclude draft-quality genomes for non-model organisms.

### `bad_max_scaffolds`
**Type:** Integer | **Default:** `500000`

The maximum number of scaffolds an assembly can have before it is flagged as low quality. Similar to `bad_max_contigs`, this measures assembly fragmentation but at the scaffold level (scaffolds are ordered/oriented contigs joined by gaps). Scaffold count matters because SynVoy identifies syntenic blocks by clustering flanking gene hits that are nearby on the same scaffold. When a genome has too many scaffolds, flanking genes may be scattered across different scaffolds even if they are physically adjacent in the real genome, producing artificially low synteny scores and missed candidate regions. The default of 500,000 is generous. For chromosome-level assemblies this value is typically under 100, while for draft assemblies it might be in the thousands. Tighten this value if you want to restrict analysis to higher-quality assemblies only.

### `bad_min_n50`
**Type:** Integer (bp) | **Default:** `5000`

The minimum N50 value (in base pairs) an assembly must have to pass quality filtering. N50 is the length such that 50% of the total assembly is contained in contigs/scaffolds of this length or longer — it measures assembly contiguity. A low N50 means the assembly is composed of many short fragments, which is problematic for synteny analysis because candidate syntenic blocks (typically 100–500 kb) may be larger than the longest scaffold. The default of 5,000 bp is extremely permissive, accepting essentially any assembly. For reliable synteny analysis, an N50 of at least 100,000 bp (100 kb) is recommended. Consider setting this to 50,000–100,000 for stricter quality control. Assemblies below this threshold are handled according to `bad_quality_policy`.

### `qc_fail_policy`
**Type:** String | **Default:** `'auto'` | **Options:** `'auto'`, `'drop'`, `'keep'`

What happens to a target genome that fails the quality check made on the FASTA file itself before the search: an N50 below `bad_min_n50` or more sequences than `bad_max_contigs`. The check runs in both modes.

- **`auto`** (default) keeps the genome, flagged, when the user chose the genomes: in pro mode, and in easy mode with `--target_species`. It drops the genome when the pipeline chose the targets (easy mode without `--target_species`).
- **`drop`** removes every failed genome before the search. If all fail, all are kept.
- **`keep`** searches every genome.

A kept genome that failed is named in a `[WARN]` line of the console and in the run summary, listed in `synvoy_report.json` under `warnings` and `summary.low_quality_genomes_searched`, its calls carry `low_quality_assembly: true`, and its row is marked with a dagger (†) in the figures, with one line in the legend. Before 2026-10-08 the default was `drop` in every mode, so a genome supplied in pro mode could disappear from the results without a row.

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

### `target_gffs`
**Type:** String (folder/glob/paths) | **Default:** `null` | **Pro Mode only**

Optional annotations for the target genomes. `--target_genomes` is FASTA-only by design, so annotations are supplied through this separate parameter — genomes and GFFs therefore live in two different folders (or globs, or comma-separated lists). Each annotation is matched to its genome by filename stem (`cow.fna` ↔ `cow.gff`) or, failing that, by assembly accession (`GCF_002263795.1.fna` ↔ `GCF_002263795.1_ARS-UCD1.2_genomic.gff`), then linked next to that genome under the one filename the search layer probes. `.gff`, `.gff3` and `.gtf` are accepted, plain or gzipped; gzipped input is decompressed during staging so every consumer can read it. Annotating a subset of targets is legitimate and does not fail the run, but an annotation that matches *no* genome is a hard error by default — this is deliberate, because the silent failure mode (a run that quietly produces no annotation-derived evidence at all) is far more expensive than a launch-time stop. Set `allow_unmatched_target_gffs = true` to downgrade that to a warning.

When a model overlaps an annotated feature, SynVoy stamps `TargetGene`, `TargetProduct` and `TargetID` onto the emitted GFF row and fills the `target_gene` / `target_product` columns in the report. Easy Mode needs no equivalent flag: `fetch_related_genomes.py` already downloads each target's annotation alongside its genome, which is why this gap only ever existed in Pro Mode.

**Scope — what this does and does not affect.** Target annotations are consulted *after* a model has been built. They do **not** influence where SynVoy searches (synteny blocks are placed from flanking-gene search hits), nor the region scoring in `cluster_grs.py`, nor the label of a call (the judge decides that from the E-value in the expected place, specificity and completeness). With `--goi_family_tokens`, every GOI row records whether `TargetGene`/`TargetProduct` carries one of the names (`GoiFamilyConsistent`, `GoiFamilyReason`); that is annotation only since 2026-10-15, when the strict gate that demoted calls by name (`strict_goi_family`) was removed. Treat `TargetGene`/`TargetProduct` as independent evidence *for you to read*, not as a scoring input.

### `allow_unmatched_target_gffs`
**Type:** Boolean | **Default:** `false`

Controls what happens when a file supplied via `--target_gffs` cannot be matched to any target genome (a typo, a leftover annotation for a genome that is not in this run, or a naming scheme the stem/accession matcher cannot resolve). The default fails the run at the staging step with a list of the offending files, on the principle that an unnoticed mismatch produces a run whose annotation-derived evidence is silently absent. Setting this to `true` reports the same list as a warning and proceeds with whatever annotations *did* match — appropriate when you knowingly point at a shared GFF directory that covers more genomes than the current run uses.

---

## 4. Synteny & Flanking Gene Extraction

### `n_flanking_genes`
**Type:** Integer | **Default:** `10` | **Range:** 2–30

The number of non-GOI-similar flanking genes to extract on **each side** of the GOI in the home genome. These flanking genes serve as synteny anchors — their presence in target genomes indicates conserved gene order (synteny). With the default of 10, SynVoy extracts up to 10 upstream and 10 downstream genes, for up to 20 total anchors. The optimal value depends on genome architecture: gene-dense genomes (bacteria, compact insects) need fewer anchors (5–8) because genes are closely packed and synteny is well-preserved, while gene-sparse or heavily rearranged genomes (plants, large vertebrates) benefit from more anchors (12–20) to increase the chance that at least some are conserved. More anchors increase computational cost (each must be searched in every target) but improve synteny detection reliability.


The flanking genes are also the unit of every distance of the search ([Distances in gene units](#distances-in-gene-units)): the pitch is their mean spacing, and the expansion of a target genome is measured on them. That measurement needs five of them in a row in the target, so with `n_flanking_genes` of 2 or less (at most four anchors) the ratio of the genome sizes is always used instead.
### `min_flanking_size`
**Type:** Integer (bp) | **Default:** `500` | **Range:** 100–5000

The minimum coding length for a gene to be written as a flanking protein. **It only acts when `exon_level_search` is `false`**, and on genes predicted without an annotation; with the default exon-level search every annotated neighbour is used whatever its length. A very short protein is a weak anchor (few residues to align), so when you switch the exon-level search off, keep this at a few hundred bp.

### `exon_level_search`
**Type:** Boolean | **Default:** `true`

When enabled, flanking gene protein sequences are decomposed into individual exon-level fragments for the initial MMseqs2 search in target genomes. This dramatically improves sensitivity for detecting divergent orthologs because individual exons may be conserved even when full-length protein similarity has dropped below detection thresholds. Multi-exon genes that have accumulated mutations unevenly across exons benefit most — a single well-conserved exon can anchor the synteny block even if other exons have diverged beyond recognition. The pipeline searches both the full-length protein and individual exon sequences, taking the best hit. Disabling this saves some computational overhead but significantly reduces sensitivity for cross-order or cross-class searches (>200 Mya divergence). It is strongly recommended to keep this enabled unless you are searching only within the same genus or family where full-length matches are expected.

### `goi_relative_evalue`
**Type:** Float (E-value) | **Default:** `1e-5` | `0` = off

A neighbour of the gene of interest (GOI) is not used as an anchor when it is a **relative** of the GOI: its hits in a target genome would fall on the GOI's own relatives and count the same sequence evidence a second time as "neighbourhood". Relative means: the neighbour's protein aligns to the GOI protein (Smith–Waterman, BLOSUM62, gap costs 11/1) with an E-value at or below this value. The E-value is computed for the two sequences alone (search space = product of their lengths).

It replaces `max_flanking_goi_similarity` (longest common subsequence divided by the longer length, threshold 35). That measure gives about 35 for two *unrelated* proteins of similar length, so it removed about half of all genes of the GOI's size whatever their sequence: on five benchmark loci 22 of 33 removed genes had no significant alignment to the GOI, and in yeast the two nearest neighbours of STE2 were dropped. With the E-value, unrelated proteins are not removed at any length (0 of 60 random pairs at 1e-5 in the tests).

`preset_tandem_family` uses `1e-3`, which also removes distant members of the family. `0` switches the test off for the neighbours; the gene at the GOI's own position is never an anchor.

### `max_flanking_distance`
**Type:** Integer (bp) | **Default:** `0` (unlimited) | **Range:** 0–5,000,000

The maximum distance in base pairs from the GOI center that the flanking gene extractor will walk to find synteny anchors. When set to 0 (default), there is no distance limit — the extractor walks as far as needed along the chromosome to collect `n_flanking_genes` anchors per side. This can be problematic when the GOI sits adjacent to a large tandem gene array: one direction fills up entirely with relatives of the GOI (which are skipped, see `goi_relative_evalue`), forcing the extractor to walk far into a distant genomic region to find enough non-similar anchors. These distant anchors may belong to a different synteny block entirely, confusing downstream analysis. Setting a distance cap (e.g., 500,000 bp) prevents this pathological case. If fewer than `n_flanking_genes` anchors are found within the distance window, the pipeline proceeds with whatever was found. This parameter is automatically adjusted by the LLM estimator for known gene family cases.

### `expand_goi_similar`
**Type:** Boolean | **Default:** `false` (was `true` until 2026-04-25)

When enabled, genes near the GOI that are similar to it (detected during flanking gene extraction) are not just excluded as anchors — they are emitted as additional GOI queries with a `GOI_NEIGHBOR_` prefix. These neighbor queries are searched in all target genomes alongside the original GOI and included in the phylogenetic tree.

**Why default changed to `false`:** for paralog-rich gene families (e.g. LY6 with ~15 chr8 paralogs) the NEIGHBOR queries flood the per-target m8 hits with paralog matches. Every cluster picks up a GOI-overlap bonus from at least one paralog, so the bonus loses its discriminating power. The true ortholog scaffold drops out of `adaptive_max_regions=6`. Empirical: LY6 ground-truth scaffold recall went from **8/9 (off)** to **0/9 (on)**. Melittin (single-copy gene, no nearby paralogs) is unaffected either way.

**When to turn it on:** you have a single-copy GOI sitting next to a small known tandem array, and you want the paralogs in the same tree to disambiguate. Examples: a Yellow gene with one Yellow-e3 paralog two genes away, a venom peptide flanked by 1–2 tandem duplicates. In these narrow cases the bonus does help.

**When to leave it off:** any gene-family query (LY6, TP53, MRJP, 3FTx, defensins) — the paralog flood will swamp the synteny signal. The phylogenetic tree downstream still resolves paralogs from the iterative search's existing GOI-similarity classification; you don't need this flag to get a paralog-aware tree.

A neighbour is "GOI-like" by the same test as for anchors (`goi_relative_evalue`).

### `expand_goi_similar_distance`
**Type:** Integer (bp) | **Default:** `300000` | **Range:** 10,000–2,000,000

The maximum distance in base pairs from the GOI within which to search for GOI-similar neighbor genes for expansion (see `expand_goi_similar`). Genes beyond this distance that happen to resemble the GOI are not emitted as additional queries. The default of 300 kb covers most tandem gene arrays, which typically span 50–500 kb in eukaryotic genomes. For very large tandem arrays (e.g., olfactory receptor clusters in vertebrates, which can span several megabases), you may need to increase this value. For compact genomes (bacteria, fungi), a smaller value (50–100 kb) is more appropriate to avoid picking up unrelated genes that coincidentally share sequence similarity. Each additional GOI-similar neighbor adds a full search pass across all target genomes, so expanding this distance excessively can significantly increase runtime without improving results.

### `locus_gap_genes` / `locus_gap`
**Type:** Integer (genes) / `auto` or Integer (bp) | **Defaults:** `5` / `auto`

The query usually aligns at several places of the home genome (exons, relatives). Two such places are **one home locus** when at most `locus_gap_genes` annotated genes lie between them. The unit is genes because the consequence is counted in genes: each locus gets its own walk for `n_flanking_genes` neighbours on each side, and two places five genes apart would share most of them. Protein-coding genes are counted where the annotation states a biotype, every `gene` row otherwise.

`locus_gap` gives the gap in base pairs instead. `auto` (default) counts genes when the home genome has an annotation and uses 50,000 bp when it has none (this was the fixed value for every genome until 2026-10-08; 50 kb are about 20 genes in yeast and about one gene or less in human).

### `home_locus_mode`
**Type:** `auto` / `gene_only` / `alignment` | **Default:** `auto`

Which home loci are searched.

- `auto`: when an annotated home gene matches the query (by name or by position, verified by alignment), **that gene is the primary locus**, as one interval with the gene's own span. Further loci are other places where the query aligns at least half as well as at the primary locus (`locus_min_bit_ratio`), up to `max_loci`. Without a matching annotated gene the loci come from the alignments alone.
- `gene_only`: only the annotated gene's locus. Use it when the further loci are relatives you do not want (replaces `--enable_name_locus`).
- `alignment`: loci from the alignments alone, as before 2026-10-08.

A locus is scored as **one copy of the gene**: the best set of its hits that covers the query at most once. Before, a locus was ranked by its single best hit, so a relative sharing one well-conserved domain scored like the gene itself (human decorin: four loci of other leucine-rich-repeat genes were searched in every genome and the locus of its real relatives was cut; now the loci are decorin, the asporin cluster and biglycan). The primary locus is never removed by `max_loci`.

---

## 5. Search Parameters

### Distances in gene units

A distance that describes a neighbourhood means something else in every genome: ten neighbouring genes span about 20 kb in yeast, 150 kb in the honeybee and 1 Mb in human. Until 2026-10-08 these distances were fixed numbers of base pairs, the same for every genome. They are now given in **gene pitches** and converted to base pairs for each home locus and each target genome:

- **pitch (home)** = span of the anchors of the home locus / (number of anchors − 1), i.e. the mean step from one anchor to the next.
- **expansion** of a target genome = the length of the longest run of anchors that lie in home order on one target sequence, divided by the length of the same run at home. One step of the run goes at most three anchors ahead and may be 0.2 to 8 times as long as at home; it is measured between two hits in the target and between the two query records of those hits at home (two exons with exon-level search, otherwise two genes). The measurement is used only where the gene order is conserved: the run must hold at least 5 anchors. Otherwise the expansion is the ratio of the genome sizes, kept between 0.2 and 8; without that it is 1. (Until 2026-10-11 the expansion was the median of the distance ratios of nearby anchor pairs. A neighbourhood mostly grows by a few large insertions, which a median does not see: *A. lyrata* at the FT locus of *A. thaliana* got 1.5 where the neighbourhood is 4.9 times as long.)
- **distance in bp** = number of pitches × pitch (home) × expansion.

| Distance | Parameter in pitches | Default | Parameter in bp | Value before |
|---|---|---|---|---|
| largest gap between two anchors of one block | `block_gap_genes` | 7 | `cluster_distance` | 150,000 |
| extent of the neighbourhood on each side of a block | `window_margin_genes` | 7 | `region_padding` | 150,000 |
| how far outside the block a GOI hit may lie | `goi_margin_genes` | 2 | `goi_margin` | 50,000 |
| largest gap a block may bridge | `synteny_bridge_max_gap_genes` | 20 | `synteny_bridge_max_gap` | 6,000,000 |
| margin of the flanking envelope for a dispersed GOI | `dispersed_envelope_margin_genes` | 9 | `dispersed_envelope_margin` | 200,000 |

Each bp parameter takes `auto` (default: derive it) or a number, which is then used as it is for every genome. The defaults in pitches are the former bp values at the locus where they were set: at the median honeybee locus of the family benchmark the pitch is 23 kb, so 150 kb are 7 pitches and 50 kb are 2; the 6 Mb bridge limit was set for cow decorin, where it is 18 pitches (20 pitches are a whole home neighbourhood of 10 genes per side).

Why "5 anchors in a row": in a genome where the gene order is lost, hits of two anchors still lie within eight times their home distance of each other by chance, and such pairs give a large factor. Measured on 23 genome–locus pairs, the longest run of anchors in home order was 6–19 where the neighbourhood is conserved (yeasts, flies, bees) and at most 3 where it is not (fruit fly against mosquito, firefly against fruit fly and mosquito). An earlier rule ("a quarter of the pairs that share a target sequence") let 11 chance pairs through for firefly against fruit fly: factor 3.7 where the genome sizes say 0.3, a block gap of 956 kb, and one search window that needed 10 GB of memory.

How much DNA is cut out around a block: the window margin is the extent of the neighbourhood (blocks whose extents overlap are merged into one, and the anchors inside it count as support for a call). Hits of the gene of interest are used only when they touch the block or lie within the GOI margin of it, and a gene model may use the DNA up to `gap_search_window` (or the gene's own span at home × expansion, if that is more) beyond its hits. DNA further out would be searched for nothing, so what is cut out and searched is the smaller of the window margin and GOI margin + that room for the model. At a yeast or honeybee locus this is the window margin; at a human locus it is about 0.7 Mb on each side where the window margin is 2.3 Mb. A window margin given in bp is searched in full, as before.

Two more rules follow the home locus. The block gap is never smaller than 1.25 × the largest gap between two neighbouring anchors at home, so a neighbourhood that is unchanged in the target comes out as one block. And hits of one flanking gene belong to one occurrence of that gene up to 1.5 × the gene's own span at home × expansion (never less than `max_intron`), instead of a fixed 20 kb for every gene (15 of 17 honeybee benchmark loci have a flanking gene with an intron longer than 20 kb).

What was used is written to the log, one line per genome, and into each result GFF as a `##synvoy-distances` line ([OUTPUT.md](OUTPUT.md)); the region scorer and the hull rescue read it from there.

`-profile <execution>,legacy_distances` sets these distances, the three of the hull rescue and the gap between home loci to their former fixed values. It is meant for comparisons (what do the derived distances change?) and for a genome where a derived window turns out too large. It does not bring back an older version: the limit for hits of one flanking gene still follows that gene's span at home, and the other changes of 2026-10 (queries, anchors, gene models, thresholds) are not distances.

### `cluster_distance` / `block_gap_genes`
**Type:** `auto` or Integer (bp) / Float (gene pitches) | **Defaults:** `auto` / `7`

The largest gap between two neighbouring flanking-gene hits inside one syntenic block. With `auto` it is `block_gap_genes` pitches of the home locus times the expansion of the target genome, and at least 1.25 × the largest gap between neighbouring anchors at home (see [Distances in gene units](#distances-in-gene-units)). A number of base pairs switches the derivation off for this distance. Measured on saved runs: with the fixed 150 kb the yeast genomes (pitch 2.4 kb) gave 7–11 blocks of up to 1.2 Mb that together covered most of an 11 Mb genome; derived, the block gap there is 20–35 kb and the block with the gene is 20–130 kb long.

### `window_margin_genes`
**Type:** Float (gene pitches) | **Default:** `7`

The extent of the neighbourhood on each side of a block, in pitches (in bp: `region_padding`). Blocks whose extents overlap are merged, and the anchors inside the extent count as support for a call. The DNA that is cut out and searched is at most this wide, and narrower when the usable part is narrower: GOI margin plus the room for one gene model (see [Distances in gene units](#distances-in-gene-units)). The log line of each genome prints both, for example `window margin 35 kb (DNA searched: 35 kb per side)` for a yeast genome, where the two are the same.

### `goi_margin_genes` / `goi_margin`
**Type:** Float (gene pitches) / `auto` or Integer (bp) | **Defaults:** `2` / `auto`

A hit of the gene of interest is used for a gene model only when it lies in the block of flanking genes or at most this far outside it. Until 2026-10-08 this was `gap_search_window` (50 kb), which also set the DNA given to miniprot; the two are separate now.

### `synteny_bridge_max_gap_genes` / `synteny_bridge_max_gap`
**Type:** Float (gene pitches) / `auto` or Integer (bp) | **Defaults:** `20` / `auto`

Two clusters of flanking genes on one sequence stay one block across a gap up to this size when the flanking genes continue the home order across it (a local rearrangement the gene may sit in). `0` in either parameter switches bridging off. The rules for *when* a gap is bridged are `synteny_bridge_max_rank_gap`, `synteny_bridge_min_anchors`, `synteny_bridge_two_sided` and `synteny_bridge_max_per_block`.

`synteny_bridge_min_anchors` (default `2`, `3` until 2026-10-14): the number of distinct flanking genes a block must hold before it may bridge a gap. With 2, a block of two genes bridges when the first gene after the gap continues their home order within `synteny_bridge_max_rank_gap` ranks, so `A B … D E F` is one block (two at `3`). Its direction then rests on those two genes. The same value is the minimum collinear run of the two-sided rule (off by default), where 2 is a weak bar.

### `dispersed_envelope_margin_genes` / `dispersed_envelope_margin`
**Type:** Float (gene pitches) / `auto` or Integer (bp) | **Defaults:** `9` / `auto`

A strong hit of the gene of interest outside every block is still accepted when it lies inside the envelope of the flanking genes on that sequence, widened by this margin; outside it the hit is listed under `rejected_candidates` in the report.

### `min_synteny_score`
**Type:** Float | **Default:** `0.6` | **Range:** 0.1–1.0

The score from which a candidate **region** is labelled HIGH in `regions/*.bed` and the region tables; half of it is the limit for MEDIUM (`cluster_grs.py`). It is a label of the region scorer only: it does not decide which blocks are searched, and it does not enter the confidence of any gene model. Which blocks are searched is decided by `min_block_genes` and `max_blocks_per_genome`. (Earlier versions of this page described it as the fraction of flanking genes needed to trigger a search; no such gate exists.)

### `min_hit_identity`
**Type:** Float (%) | **Default:** `10` | **Range:** 0–100

A floor on the identity of an alignment hit, applied in the genome-wide search for the anchors. It is a guard, not the filter that decides: hits are kept by E-value (`search_evalue`), and an alignment that reaches E ≤ 0.01 in a whole genome is far above 10 % identity. Inside a search window the floor in use is `max(aug_relaxed_identity_min, min_hit_identity × aug_relaxed_identity_factor)` = 15 % at the defaults. Lowering the value therefore changes nothing (which is why `preset_short_peptide` no longer sets it); raising it, for example to 30 for a search among close species, removes weak anchor hits and the chance blocks they create.

### `min_query_length`
**Type:** Integer (amino acids) | **Default:** `20` | `0` = no check

Queries shorter than this stop the run in `NORMALIZE_QUERY` (after DNA-to-protein translation). A query of 20–39 aa runs with a warning in the console and in the end summary.

The default was 30. Measured on the honeybee genome with peptides cut from annotated proteins: an exact peptide is found at its own place with E ≤ 0.01 in 98 % of cases at 18 aa and 89 % at 15 aa, so the home step works from about 20 aa. What gets weak below 40 aa is the search in other species. An ortholog at 60 % identity scores about 25 bits for a 20-aa query, 35 bits for 30 aa and 45 bits for 40 aa (simulated, BLOSUM62 11/1). A 300 Mb genome is expected to hold about 270 chance alignments as good as the first, 0.4 as good as the second and 0.0005 as good as the third. Below 40 aa a call therefore rests on the neighbourhood, not on the alignment; that is what the warning says.

### `min_hit_length`
**Type:** Integer (amino acids) | **Default:** `10` | **Range:** 5–500

A floor on the aligned length of a hit, applied in the genome-wide search for the anchors. Like `min_hit_identity` it is a guard behind the E-value filter: an alignment of fewer than 10 residues scores about 20 bits at best, and E ≤ 0.01 in a genome needs 35 bits or more. Inside a search window the floor in use is `max(aug_relaxed_length_min, min_hit_length / aug_relaxed_length_div)` = 15 aa at the defaults. Lowering the value changes nothing; raising it removes short anchor hits.

### `search_evalue`
**Type:** Float | **Default:** `0.01` | **Range:** 1e-10–100

The E-value threshold for tblastn and MMseqs2 homology searches. The E-value represents the expected number of hits with equal or better score by chance in a database of the given size. Lower values (0.001, 1e-5) are more stringent, reporting only highly significant hits; higher values (0.1, 1.0, 10) are more permissive, allowing weaker matches through. The default of 0.01 provides a good balance for typical within-class or within-order searches. For very close species (same genus), tightening to 0.001 reduces noise. For distant cross-phylum searches (>400 Mya), relaxing to 0.1–1.0 is necessary because genuine orthologs may have diverged enough that their alignment scores are only marginally above background. The pipeline applies additional filters (synteny score, identity, length) downstream, so a permissive E-value here is safer than a restrictive one — it lets the synteny context disambiguate weak hits rather than discarding them prematurely.

### `max_intron`
**Type:** Integer (bp) | **Default:** `20000` | **Range:** 0–500,000

The floor of the largest gap between two hits of one flanking gene that still count as one occurrence of that gene. The limit used is `max(max_intron, 1.5 × the gene's span at home × expansion)`, so a gene with long introns is not split into several "copies" and a compact gene keeps the floor. For the gene of interest the corresponding limit comes from its home model (`goi_fallback_intron_margin`, `goi_fallback_intron_floor`). miniprot is not bound by this value: it gets the whole stretch of DNA of a candidate and may place an intron anywhere in it.

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

A floor on the raw Smith-Waterman score (BLOSUM62, gap 11/1) of a hit in a search window. Since September 2026 every Smith-Waterman hit must also reach E ≤ 1 over the six reading frames of its window, and that filter is the stricter one for any value up to about 30: the lowest raw score with E ≤ 1 is 31 for a 20-aa query in a 2 kb window, 54 for a 70-aa query in 300 kb and 60 for a 300-aa query in 300 kb. The default, and the value of 10 that `preset_short_peptide` used to set, therefore never decide anything. Only a value above those numbers has an effect (it makes the Smith-Waterman search stricter than E ≤ 1).

### `sw_min_identity`
**Type:** Float (%) | **Default:** `10.0` | **Range:** 0–100

A floor on the identity of a Smith-Waterman hit, applied together with `sw_min_score` and before the E ≤ 1 filter. An alignment that reaches E ≤ 1 with less than 10 % identical residues does not occur in practice, so the default is a guard that does not bind; raising it (for example to 25) makes the Smith-Waterman search stricter.

### `sw_timeout_seconds`
**Type:** Integer | **Default:** `300` | **Range:** 1–3600

The maximum time in seconds allowed for a single Smith-Waterman search invocation. SW alignment is O(m*n) where m is the query length and n is the target region length; for large regions (>1 Mb) or large query proteins (>1000 aa), this can take several minutes even with SIMD-accelerated parasail. The timeout prevents a single pathological case (e.g., a highly repetitive region) from stalling the entire pipeline. When the timeout is hit, the SW search for that specific region is aborted and the pipeline proceeds with whatever tblastn and miniprot hits were found. The default of 300 seconds (5 minutes) is sufficient for most cases. Increase this for very large query proteins or when running on slow hardware; decrease it if you prioritize throughput over exhaustive sensitivity. The timeout applies per-region, not per-genome.

---

## 7. Augmented / Relaxed Search

### `region_padding`
**Type:** `auto` or Integer (bp) | **Default:** `auto`

The extent of the neighbourhood on each side of a block, in bp. Blocks whose extents overlap are merged before the search, and the anchors inside the extent count as support for a call. `auto`: `window_margin_genes` pitches of the home locus times the expansion of the target ([Distances in gene units](#distances-in-gene-units)); the DNA that is cut out and searched is then at most this wide and no wider than GOI margin + the room for one gene model. A number of bp restores the former behaviour: the margin is adapted to the spacing of the hits (twice the mean gap between consecutive hits, between `padding_min` and `padding_max`) and searched in full.

### `padding_min`
**Type:** Integer (bp) | **Default:** `50000`

Lower limit of the adaptive margin. Only read when `region_padding` is a number.

### `padding_max`
**Type:** Integer (bp) | **Default:** `200000`

Upper limit of the adaptive margin. Only read when `region_padding` is a number.

### The thresholds of the window search (`aug_relaxed_*`)

The names come from an earlier design with several rounds of ever more relaxed searches. There is one search per window now. MMseqs2 (E 10) and tblastn (E 1000) are run on the window with a wide E-value, and of what they return a hit is kept when it reaches

- E ≤ `search_evalue × aug_relaxed_parse_evalue_mult` (0.01 × 10 = **0.1**, counted over the window, not the genome),
- identity ≥ `max(aug_relaxed_identity_min, min_hit_identity × aug_relaxed_identity_factor)` (**15 %**), and
- aligned length ≥ `max(aug_relaxed_length_min, min_hit_length / aug_relaxed_length_div)` (**15 aa**).

Smith-Waterman hits have their own rule (E ≤ 1 over the six frames of the window, see `sw_min_score`).

### `aug_relaxed_parse_evalue_mult`
**Type:** Float | **Default:** `10`

The factor between the E-value for anchors in the whole genome (`search_evalue`) and the E-value at which a hit inside a search window is kept. This is the one of the seven that decides which window hits exist. An E-value is already scaled to the search space, so 0.1 in a window of a few hundred kb asks for fewer bits than 0.01 in a genome: that is the intended relaxation where anchors support the place.

### `aug_relaxed_evalue_mult`
**Type:** Float | **Default:** `1000` | **Range:** 10–100,000

### `aug_relaxed_evalue_cap`
**Type:** Float | **Default:** `10.0` | **Range:** 0.1–100

Together they set the E-value handed to MMseqs2 for the window: `min(aug_relaxed_evalue_cap, search_evalue × aug_relaxed_evalue_mult)` = 10. Hits between that value and the kept threshold of 0.1 are computed and then discarded, so neither parameter changes the result as long as their product stays above the kept threshold. Lowering them below it would tighten the window search; raising them does nothing.

tblastn does not use them (since 2026-10-14): it always reports to E 1000. For tblastn the cutoff is more than an output filter. It combines HSPs by sum statistics only among the ones it keeps, so with the old cutoff of 10 weak HSPs that would have joined a significant set were thrown away first, and the E-value of what it did report changed with the setting. On the Candida pilot a curated ortholog was absent from tblastn's output in SynVoy's own window at 10 and reported at E 0.005 at 1000. The kept threshold of 0.1 applies to tblastn as before.

### `aug_relaxed_identity_min`
**Type:** Float (%) | **Default:** `15.0` | **Range:** 5–50

### `aug_relaxed_identity_factor`
**Type:** Float | **Default:** `0.6` | **Range:** 0.1–1.0

The identity floor of a window hit is the larger of `aug_relaxed_identity_min` and `min_hit_identity × aug_relaxed_identity_factor`. At the defaults that is max(15, 6) = 15 %; the factor only matters when `min_hit_identity` is set above 25.

### `aug_relaxed_length_min`
**Type:** Integer (aa) | **Default:** `15`

### `aug_relaxed_length_div`
**Type:** Integer | **Default:** `2`

The length floor of a window hit is the larger of `aug_relaxed_length_min` and `min_hit_length / aug_relaxed_length_div`. At the defaults that is max(15, 5) = 15 aa; the divisor only matters when `min_hit_length` is set above 30. The 15 aa are a real limit for short queries: an exon of fewer than 15 residues cannot produce a window hit from MMseqs2 or tblastn (Smith-Waterman has no such floor beyond 10 aa).

### `max_blocks_per_genome`
**Type:** Integer | **Default:** `80` | **Range:** 5–500

A safety cap on the total number of candidate syntenic blocks evaluated per target genome. In genomes with many repetitive regions or large gene families, flanking gene hits may cluster into hundreds of candidate blocks. Evaluating each block involves gene prediction, tblastn, SW alignment, and classification — so the cost is approximately linear in block count. The default of 80 is sufficient for even complex genomes with multiple GOI loci (e.g., a gene family with 10 tandem copies on 5 chromosomes = 50 blocks). For very large gene families or whole-genome searches, increase to 120–200. For quick exploratory runs, decrease to 20–30. When the cap is hit, blocks are prioritized by synteny score (highest first), so the most promising regions are always evaluated. The cap prevents runaway searches in pathological genomes from stalling the entire pipeline for hours.

### `min_block_genes`
**Type:** Integer | **Default:** `2` | **Range:** 1–10

The minimum number of flanking gene hits required in a candidate syntenic block for it to be kept and evaluated further. Blocks with fewer than this many hits are discarded as noise. With the default of 2, at least two flanking gene anchors must map to the same target region (within `cluster_distance`). Setting this to 1 is extremely permissive — a single flanking gene hit somewhere on a chromosome does not constitute synteny evidence and is almost certainly a standalone ortholog of that flanking gene rather than a conserved syntenic block. Setting this higher than 3 is very strict and may miss real syntenic regions in organisms with high chromosomal rearrangement. The combination of `min_block_genes` and `min_synteny_score` together determines which candidate regions survive for GOI search: `min_block_genes` is the absolute floor, while `min_synteny_score` is the proportional threshold.

### `fallback_short_query_len` / `fallback_short_min_aln_aa` / `fallback_short_min_bits`
**Types:** Integer (aa) / Integer (aa) / Float (bits) | **Defaults:** `150` / `15` / `30`

These control the **short-query GOI fallback gate** — the "fallback" model emitted (from tblastn/MMseqs2/Smith-Waterman hits) when miniprot cannot build a spliced gene model in a flanking-seeded block. The legacy gate vetoed a fallback on low *query coverage*, which is structurally low for a short precursor whose conserved region is only a fraction of the query: a 26 aa mature melittin peptide is ~37% of the 70 aa preprotein, so a perfect hit still looks "low coverage" and was rejected — while a high-coverage but low-identity artifact in a *different* block passed, stamping the GOI on the wrong locus.

For a query shorter than `fallback_short_query_len` (aa), SynVoy instead gates the fallback on **absolute aligned length** (`fallback_short_min_aln_aa`) AND **bitscore** (`fallback_short_min_bits`) — coverage is no longer a veto — so the real, divergent hit is emitted in its own flanking-seeded block. A micro-window (a few aa) is still rejected by the aligned-length floor, and a statistically weak hit by the bitscore floor. Queries at or above `fallback_short_query_len` keep the coverage floor unchanged, so long single-copy or paralog-family genes (e.g. LY6, TP53) are unaffected, and a high-coverage divergent ortholog (e.g. the *Tetramorium* myrmicitoxin at ~40% identity × ~76% coverage) clears either path. Set `fallback_short_query_len 0` to disable the short-query path entirely (legacy coverage gate for all queries). Lower `fallback_short_min_aln_aa`/`fallback_short_min_bits` for extremely short or divergent peptide families; raise them if a short-query run produces spurious LOW fallbacks. This is the core of the §A1 "fumble" fix.

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

### `deterministic_goi_search`
**Type:** Boolean | **Default:** `true`

Makes the iterative per-target GOI search **reproducible run-to-run**. MMseqs2 with more than one thread can give slightly different results for marginal/divergent hits between identical runs (non-deterministic prefilter/alignment tie-breaks). On its own that jitter is small, but the iterative search amplifies it: a borderline GOI call flips between confidence classes, which flips whether it seeds the expanding query database for later waves, which cascades into wholesale-different calls. With this on (default), the per-region augmented GOI search runs at `--threads 1`; because each augmented search covers only a small flanking-seeded block (not the whole genome), the speed cost is minor. (tblastn already runs single-threaded; the genome-wide flanking MMseqs2 search is unaffected.) Combined with deterministic ordering of the expanding database and of hit parsing, this yields byte-identical GOI calls and seed sets across reruns — important for publication-grade reproducibility. Set to `false` for the legacy multithreaded path if you are willing to trade reproducibility for a small speed gain on the augmented passes.

---

## 9. Annotation & Gap Filling

### `gff_search_window`
**Type:** Integer (bp) | **Default:** `100000`

Used once, when the home gene model is built (`annotate_goi_exons.py`). If the GOI is not found in the home GFF by name, the proteins of the annotated genes within this distance of the query's hits are compared with the query, and the best match is taken as the home gene. It does not affect the choice of flanking genes (that is `n_flanking_genes` and `max_flanking_distance`) or anything in the target genomes. Raise it only if the query's hits lie far from the gene they belong to.

### `gap_search_window`
**Type:** Integer (bp) | **Default:** `50000`

DNA added on each side of the hits of one candidate before miniprot builds its gene model. It is a floor: the margin is never smaller than the span of the gene of interest at home times the expansion of the target, so a gene of 200 kb gets room for its introns. The home gene model uses the same value when it falls back to miniprot. There is no separate gap-filling search: missing exons are found by miniprot inside this stretch, or not at all. (Until 2026-10-08 the same number also decided how far outside a block a hit may lie; that is `goi_margin` now.)

### `goi_max_models_per_locus`
**Type:** Integer | **Default:** `10` | `1` = best model only

Gene models kept per candidate locus. miniprot reports every place where the query aligns in the stretch of DNA; until 2026-10-08 the first model was used and the rest dropped, so a tandem array gave one model, and which copy it was depended on a score difference of under one per cent. Now every model is kept whose coding sequence overlaps the coding sequence of a better model by less than half. Where hits of the window search remain on no model, the coding exons found so far are masked and miniprot runs again (with a long allowed intron it does not report a diverged neighbouring copy in the first run: of two copies 2 kb apart the one at 80 % identity is missing from 10 kb on). A further model must reach 30 % of the best score, the limit miniprot applies itself. Further models carry `m2`, `m3` … in their id (`…_b0_l1m2_exon_ann`) and are classified like any other model.

### `gap_min_identity`
**Type:** Float (%) | **Default:** `15.0`

Despite the name, this is only a lower bound of the identity required for a `raw_hit` GOI row, the last-resort output when neither miniprot nor a hit chain gives a model. The bound is `max(gap_min_identity, 90)`, so the default has no effect and only a value above 90 changes anything. A `raw_hit` is always LOW confidence.

### `gap_max_hits`
**Type:** Integer | **Default:** `5`

The number of candidate hits looked at for the `raw_hit` GOI row of a candidate locus (the last-resort output described under `gap_min_identity`): the hits are sorted by identity, the first one that passes is written, and the search stops there, so a locus gets one such row at most. Must be at least 1.

### `min_exon_alnlen`
**Type:** Integer (aa) | **Default:** `30`

Despite the name, this no longer filters exon hits (that filter was switched off in March 2026, because it removed the short exon hits the gene model needs). It is one of the three bounds of the last-resort `raw_hit` row (see `gap_min_identity`): the hit must align at least `max(min_exon_alnlen, 60 % of the query length)` residues. With the default it therefore only matters for queries shorter than 50 aa. None of the 7,470 result files of earlier runs contains a `raw_hit` row. With miniprot installed the row also needs two hits, 75 % query coverage and 60 bits.

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

The **label** of a GOI call (`Confidence` HIGH / MEDIUM / AMBIGUOUS / LOW, with `Basis`) has no parameters: one judge sets it from the E-value in the expected place, specificity against the home proteome and a completeness bar of 0.7 query coverage, all fixed significance levels ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)). The identity thresholds of the earlier label and its switches were removed on 2026-10-15 ([§ 20](#20-removed-parameters)). The parameters below set the **model status** (`ModelStatus`: fragment / partial / complete), the family-name annotation, and the report's coverage rule.

### `goi_family_tokens`
**Type:** Comma-separated string | **Default:** `''` (empty)

Family names. With them, every GOI row records whether the target genome's own annotation of the gene (`TargetGene` / `TargetProduct`) carries one of them: `GoiFamilyConsistent=true/false` and `GoiFamilyReason` (`matched_gene:<name>`, `matched_product:<name>`, `no_family_match`, `no_target_annotation`). Annotation only: it changes no label (the strict gate that demoted hit-chain calls without a matching name, `strict_goi_family`, was removed on 2026-10-15). Without names the check is off. Matching is case-insensitive substring after normalization (strip non-alphanumeric): `P53` matches `Trp53` (mouse ortholog) because `P53` is a substring of `TRP53` once normalized. For TP53: `--goi_family_tokens TP53,TP63,TP73,TRP53,TRP63,TRP73,P53,P63,P73` (or `-c conf/tp53_strict.config`). No preset supplies names: `preset_paralog_discrimination` used to ship the TP53 names and, being auto-applied to any large gene family, made the strict gate demote calls of unrelated queries (fixed 2026-10-05). (Until 2026-10-09 the code also tried to derive names from the first record of the search database. That record is a flanking exon written by the pipeline, so nothing was ever derived.)

### `classify_fragment_max_qcov`
**Type:** Float | **Default:** `0.4` | **Range:** 0–1

The maximum query coverage fraction below which a gene model is labeled as a `fragment` in the ModelStatus annotation field. Query coverage measures what fraction of the query protein length is covered by the aligned portion of the predicted gene. A model covering less than 40% of the query is considered a fragment — it likely represents an incomplete gene prediction (truncated by scaffold boundaries, missing exons, or pseudogenization). Fragment status is independent of the label (HIGH/MEDIUM/AMBIGUOUS/LOW): the label's completeness test is the judge's own fixed bar (0.7 query coverage for HIGH), and a fragment can be MEDIUM when its E-value in the expected place is significant. This distinction is valuable for downstream analysis: fragments may represent pseudogenes, assembly artifacts, or genes at scaffold edges where the rest of the sequence is missing. Setting this higher (e.g., 0.5) classifies more models as fragments; lower (e.g., 0.3) is more lenient.

### `classify_complete_min_qcov`
**Type:** Float | **Default:** `0.7` | **Range:** 0–1

The minimum query coverage fraction for a gene model to be labeled as `complete` in the ModelStatus field. Models must cover at least 70% of the query protein and have multi-exon evidence (or be identified as tandem copies) to earn the `complete` status. Models between the fragment threshold (40%) and the complete threshold (70%) are labeled `partial` — they have significant alignment to the query but likely miss one or more exons or terminal regions. The `complete` status indicates a gene model that is likely to represent a full-length or near-full-length ortholog suitable for functional annotation and phylogenetic analysis without qualification. These thresholds are intentionally conservative: 70% coverage allows for some divergence at the N- or C-terminus (signal peptides, disordered tails) while ensuring the core protein is well-represented in the gene model. This parameter sets the model status only: the judge's completeness bar for HIGH is a fixed 0.7 that `preset_short_peptide` (0.5 here) does not move.

### `identity_decoupled_min_identity` / `identity_decoupled_max_qcov`
**Type:** Float | **Defaults:** `50.0` (%) / `0.35`

A call at or above the identity bound whose *recorded* query coverage is below the coverage bound is a short local window, not an ortholog (ten residues at 100 %). The report relabels it `identity_coverage_decoupled` and leaves it out of the HIGH/MEDIUM headline counts. A call with no recorded coverage is never demoted by this rule.

### `disable_coverage_demotion`
**Type:** Boolean | **Default:** `false`

`true` keeps the rule above advisory: the flag still appears in `self_consistency`, the call keeps its confidence.

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

The weight of strand agreement in the quality term. For each flanking gene found in the target region SynVoy compares its strand with the strand of its home gene and takes the larger of "agree" and "disagree" over all of them, so a **uniformly inverted** neighbourhood scores 1.0 — an inversion keeps relative orientation and is conserved synteny — while a neighbourhood with mixed agreement scores lower. (Before 2026-07-26 this term never looked at the home genome and rewarded scrambled neighbourhoods; `--legacy_strand_score` reproduces those numbers.)

Note what the three weights do and do not mean. The score is `synteny_score = Q × C` with `Q = w_base·C + w_consistency·K + w_strand·S`, where `C` is the fraction of flanking genes found. Coverage therefore enters twice and dominates: order (`K`) and strand (`S`) can move the score by at most `0.3·C` each. Describe the score as coverage-dominated with order and strand as modifiers, not as a 40/30/30 split.

### `synteny_goi_overlap_bonus`
**Type:** Float | **Default:** `0.15`

A bonus added to the **ranking** score of a candidate region that overlaps a GOI call made by the search itself (it has nothing to do with the target genome's own annotation). It keeps the region that actually holds a GOI model ahead of an equally syntenic region that does not. The bonus is prior knowledge SynVoy injects, not evidence from the neighbourhood, so it is deliberately excluded from the permutation p-value, which tests `synteny_score` alone; both numbers are written to `regions/*.scores.tsv`.

### `disable_distant_synteny_rescue`
**Type:** Boolean | **Default:** `false`

In region ranking, a strong GOI hit (high identity, very low E-value) that is isolated from its local flanking genes but lies on a chromosome that still carries an anchor, in a genome that is visibly rearranged, is promoted to a floor score instead of being dropped — it is treated as a dispersed ortholog. A region with strong flanking support never takes this route. Set `true` to skip the rescue.

### `legacy_strand_score`
**Type:** Boolean | **Default:** `false`

Compute the strand term as before 2026-07-26 (`max(#plus, #minus) / k` over the target hits, ignoring the home strand). Only for reproducing or comparing against region scores from older runs.

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

The PLM cosine similarity threshold above which a LOW-confidence gene model can be boosted to MEDIUM confidence (in a block with at least one flanking gene). Since 2026-10-14 this acts only on a GOI row that the judge did not label (no `Basis`, which in practice no row of the main search lacks): a judged row keeps the judge's label, and the score is recorded on it (`EmbeddingSimilarity` / `StructuralSimilarity`). A gene model that has poor sequence-level evidence but strong embedding similarity (>0.7) to the query suggests that the protein is a genuine homolog whose sequence has diverged beyond reliable alignment-based detection. The boost recognizes that PLM embeddings provide orthogonal evidence to sequence alignment — high embedding similarity despite low sequence identity is a hallmark of remote homology. The default of 0.7 requires strong embedding similarity for the boost, which is appropriate because boosting a LOW model to MEDIUM has significant implications for the interpretation of results.

### `plm_high_threshold`
**Type:** Float | **Default:** `0.85` | **Range:** 0–1

The PLM cosine similarity threshold above which a MEDIUM-confidence gene model can be boosted to HIGH confidence (in a block with at least two flanking genes). Since 2026-10-14 this acts only on a GOI row that the judge did not label (no `Basis`, which in practice no row of the main search lacks): a judged row keeps the judge's label, and the score is recorded on it (`EmbeddingSimilarity` / `StructuralSimilarity`). An embedding similarity of 0.85 indicates near-structural equivalence, providing very strong evidence for homology even when sequence-based evidence is only moderate. This is the most conservative PLM boost tier — only clearly related proteins at the structural/functional level receive the HIGH designation based on PLM evidence. This threshold should rarely be lowered because HIGH confidence is used for making definitive orthology claims, and PLM embeddings, while powerful, are not infallible. The gap between the medium threshold (0.7) and the high threshold (0.85) creates a meaningful distinction between "probably homologous" (MEDIUM boost) and "almost certainly homologous" (HIGH boost) based on PLM evidence.

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

The TM-score threshold above which structural evidence can boost a LOW-confidence gene model to MEDIUM confidence (in a block with at least one flanking gene). Since 2026-10-14 this acts only on a GOI row that the judge did not label (no `Basis`, which in practice no row of the main search lacks): a judged row keeps the judge's label, and the score is recorded on it (`EmbeddingSimilarity` / `StructuralSimilarity`). A TM-score of 0.5 indicates that the candidate shares the same overall fold as the query — this is strong structural evidence for homology even when sequence-based evidence is weak. Structural homology at this level is rarely coincidental: convergent evolution to the same fold does occur but is relatively rare, and the syntenic context further reduces the chance of false positives. This boost is particularly valuable for ancient gene families where sequence has diverged beyond recognition but the 3D fold is preserved. The combination of synteny evidence plus structural similarity at TM>0.5 provides high-confidence orthology assignment even for cross-kingdom comparisons.

### `structural_high_threshold`
**Type:** Float | **Default:** `0.7` | **Range:** 0–1

The TM-score threshold above which structural evidence can boost a MEDIUM-confidence gene model to HIGH confidence (in a block with at least two flanking genes). Since 2026-10-14 this acts only on a GOI row that the judge did not label (no `Basis`, which in practice no row of the main search lacks): a judged row keeps the judge's label, and the score is recorded on it (`EmbeddingSimilarity` / `StructuralSimilarity`). A TM-score of 0.7 indicates strong structural similarity — the proteins are almost certainly in the same protein family with a shared evolutionary origin. At this level, structural homology provides evidence comparable to moderate sequence identity (30–40%), justifying the HIGH confidence designation. Combined with synteny context, a TM-score >0.7 hit is a very reliable ortholog identification even in the complete absence of detectable sequence similarity. This is the state-of-the-art for detecting orthologs separated by >500 million years of evolution, where sequence methods have essentially zero sensitivity.

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
**Type:** `'auto'` or Integer (bp) | **Default:** `'auto'` (10000 until 2026-10-09)

The length of the scale bar in the ribbon plot. `auto` picks a round length (1, 2 or 5 × 10^k bp) whose bar is about 120 px wide at the scale of the drawing. The former fixed 10 kb was 34 px wide in a yeast plot and 3 px in a human one (measured on the STE2 and the decorin run), so it could not be read where genes are far apart. A number fixes the length. The bar refers to the uncompressed parts of a track; gaps longer than `gap_threshold` are drawn shortened.

### `hide_goi_absent_tracks`
**Type:** Boolean | **Default:** `false` (`true` until 2026-10-09)

By default every searched genome is drawn in every figure, also one in which no GOI call was made: a searched genome without the gene is a result. In the anchor grid such a genome keeps its row with the flanking genes that were placed, its GOI cell carries a dash, and a row in which nothing was placed reads "no flanking gene placed". The figure legend (`*_anchor_grid_print.legend.txt`) names these genomes with the number of flanking genes placed in each. In the ribbon plot they are drawn below the genomes with a call and labelled "GOI not placed"; the gene-position map prints "no orthologs mapped" for an empty row. A missing model is not evidence that the gene is absent.

Set it to `true` to leave these genomes out of the figures (the former default; on 26 finished loci it removed 16 of 90 searched genomes from the figures without a word in them). They stay in `synvoy_report.json` either way.

### `plot_extra_args`
**Type:** String | **Default:** `''`

Extra options appended verbatim to the `bin/plot_synteny.py` call. Use it for the figure controls that have no parameter of their own, quoting the whole value:

```
--plot_extra_args '--print_grid_numbers --no_fragment_variant'
```

| Option | Default | Effect |
|---|---|---|
| `--print_grid_numbers` | off | Write the identity into every flanking arrow of the print anchor grid (`*_anchor_grid_print.svg`). Off by default: at 183 mm the number is 5 pt, and whether it helps depends on the size the figure is printed at. The full grid (`*_anchor_grid.svg`) always has it. A query coverage below 80 % is written in brackets under the arrow, which makes that row about 1 mm taller; a figure that would pass `--print_height_mm` that way marks the identity with `*` instead (the log says so). |
| `--print_width_mm` | `183` | Width of the print grid in mm: 183 for two journal columns, 89 for one. When the grid does not fit, the spacing shrinks first; then the columns with the fewest orthologs are left out and named in the log and in `_print.legend.txt`. Text stays at 5 to 7 pt. |
| `--print_height_mm` | `170` | Largest height of the print grid in mm. Rows are packed tighter before it is exceeded; a grid that still does not fit is written anyway and the log says so. |
| `--no_print_grid` | off | Do not write the print grid. |
| `--grid_max_models` | `10` | Most GOI gene models drawn for one genome in the anchor grid. A longer tandem array shows the first ten and `+N more`; every copy keeps its small arrow and its row in `_print.source.tsv`. |
| `--no_fragment_variant` | off | Skip the `*_with_fragments` figures. By default GOI models with `ModelStatus=fragment` are hidden in the main ribbon plot and grid (they are single-exon hits that otherwise stack into extra lanes) and drawn in these sibling files. |
| `--no_orient_to_home` | off | By default a scaffold whose flanking genes run opposite to the home genome is mirrored so the ribbons run straight. This switches the mirroring off. |
| `--no_anchor_grid` | off | Do not write the anchor grid. |
| `--gene_names_tsv` | _(none)_ | Names table (`gene_id`, `label`, optional `full_name`, `source`) that sets the label of the home genes it lists; give an absolute path. Without it a home gene is labelled with its symbol in the home GFF, else its current symbol at NCBI Gene, else an abbreviation of its product name (marked `*`), else its ID. Each run writes the table it used to `<locus>_gene_names.tsv`. |
| `--no_gene_legend` | off | Do not draw the "Gene names" table (label, full name, gene ID, source) under the figures. |
| `--no_network` | off | Skip the plot step's NCBI lookups: species common names, and gene symbols of home genes the GFF leaves unnamed. |
| `--max_goi_per_genome` | `10` | How many GOI models steer the choice of which neighbourhood to draw. Every GOI model inside the chosen view is drawn regardless. |
| `--goi_zoom`, `--goi_min_px`, `--max_legend_entries`, `--caption_file` | — | Controls for tandem-array figures; see `bin/plot_synteny.py --help`. |

Do not repeat an option the pipeline already sets (`--plot_width`, `--gap_threshold`, …); use the parameter instead. You can also re-render any figure by hand from a finished run's `plot_inputs_*` folder.

**Removed on 2026-10-07: `--grid_goi_style` and `--grid_goi_max_models`.** They belonged to the earlier anchor grid: the first chose one of four ways to draw a GOI cell (`cds`, `cds_aligned`, `genomic`, `notched`), the second set how many models a cell drew before it showed the best one plus `×N`. The grid now has one GOI layout ([OUTPUT.md § The GOI column of the anchor grid](OUTPUT.md#the-goi-column-of-the-anchor-grid)), so neither has anything left to control. `plot_synteny.py` stops with `unrecognized arguments` when it gets one, and a pipeline run then ends at the plot step: take the option out of `--plot_extra_args` and re-run with `-resume` (the search is not repeated). `--grid_max_models` is the nearest equivalent of the second.

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
**Type:** Boolean | **Default:** `false`

> **Advisory on Nextflow ≥ 25.** Current Nextflow does not allow a pipeline to change its parameters after launch, so the estimated values cannot be applied within the run that computed them. SynVoy prints a warning with the exact flags to re-run with and writes the estimate to `intermediate/estimate_params/estimated_params.json`. Until 2026-10-05 the log claimed every value was applied while every step kept its launch-time value. To use an estimate: run once with `--auto_params true`, copy the suggested flags, re-run with them.

Master switch for parameter estimation (default off). From the kingdom of the home species and the evolutionary distance to the targets, SynVoy suggests values for up to 10 search parameters: the number of flanking genes, search sensitivity and E-value, the block limits, and the optional search layers. With an API key (`llm_api_key` or the relevant environment variable) a cloud LLM makes the suggestion; without one, built-in rules do. Distances, hit floors and the family options are not estimated: distances are derived from the run (see "Distances in gene units"), and lowering a hit floor changes nothing (see `sw_min_score`).

### `llm_provider`
**Type:** String | **Default:** `'google'` | **Options:** `'google'`, `'openai'`

The LLM provider for parameter estimation. `'google'` uses Google Gemini (default model: `gemini-2.5-flash-lite`). `'openai'` uses OpenAI's chat completions API (default model: `gpt-4o-mini`) or any OpenAI-compatible endpoint when `llm_api_base_url` is set — this covers Together, Groq, Fireworks, LM Studio, and similar services. The API key is resolved from `llm_api_key`, then from the provider-specific env vars (`GOOGLE_API_KEY` for google, `OPENAI_API_KEY` for openai), then from the generic `LLM_API_KEY` env var.

### `llm_api_key`
**Type:** String | **Default:** `''` (empty)

The API key for the chosen LLM provider. Can be provided here or (recommended) via environment variable to avoid exposing the key in command history or Nextflow logs. The key is resolved in this order: `--llm_api_key` CLI flag → `LLM_API_KEY` env var → `GOOGLE_API_KEY` env var (for google provider) → `OPENAI_API_KEY` env var (for openai provider). Each call uses approximately 1,000–2,000 tokens, negligible cost on any provider's free or standard tier. **Never commit API keys to version control.**

### `llm_api_base_url`
**Type:** String | **Default:** `''` (empty)

Custom API base URL for OpenAI-compatible providers. Leave empty to use the default endpoint for the chosen provider (`https://api.openai.com` for openai, Google's endpoint for google). Set this to use Together (`https://api.together.xyz`), Groq (`https://api.groq.com/openai`), LM Studio (`http://localhost:1234`), or any other OpenAI-compatible service. Ignored for the google provider.

### `llm_model`
**Type:** String | **Default:** `''` (empty = provider default)

Model name override. When empty, uses the provider default: `gemini-2.5-flash-lite` for google, `gpt-4o-mini` for openai. Set this when using an OpenAI-compatible provider with a specific model name (e.g. `meta-llama/Llama-3.1-8B-Instruct-Turbo` on Together). The model receives a detailed system prompt with SynVoy-specific biological reasoning — any instruction-following model with JSON output support works.

---

## 18. Advanced & Output

### `outdir`
**Type:** Path (String) | **Default:** `'results'`

The directory where all pipeline output files are written. This includes synteny plots (SVG/HTML), phylogenetic trees (Newick), BED region files, the JSON summary report, and optionally intermediate files. The directory is created automatically if it does not exist. Using descriptive directory names (e.g., `results/melittin_apis_vs_bombus`) helps organize multiple runs. When using `-resume`, Nextflow caches task results in the `work/` directory and regenerates outputs in `outdir` — so changing `outdir` and resuming will write results to the new location without re-running completed tasks. For HPC runs, consider pointing this to a permanent storage location rather than scratch space that may be automatically purged.

### `docker_container`
**Type:** String | **Default:** `'synvoy-local:latest'`

The Docker/Singularity container image used for running pipeline processes when a container profile is active. The default expects a locally built image tagged `synvoy-local:latest`. Build it from the project Dockerfile: `docker build -t synvoy-local:latest .`. The container bundles all dependencies (MMseqs2, tblastn, Augustus, Prodigal, miniprot, samtools, parasail, MAFFT, IQ-TREE, and the Python stack) in a reproducible environment. Override this to use a pre-built remote image (e.g., from Docker Hub or a private registry). The `beforeScript` directive ensures workspace scripts in `bin/` override any scripts packaged inside the container, so local code changes take effect immediately without rebuilding. This parameter is only used when a Docker or Singularity profile is active — Conda profiles ignore it entirely.

---

## 19. Auto-Apply Preset, Self-Consistency & Rescue

A small group of orthogonal toggles that landed together in the 2026-05-27/29 batch. All are safe defaults — change them only if you understand what you give up.

### `auto_apply_preset`
**Type:** Boolean | **Default:** `true`

When true (default), the preset chosen from the LOCATE_GENE hit profile is applied at runtime to every downstream process via a value channel (see §1f auto-apply in the Quick start above). When false, the recommendation is kept advisory-only and the params used by `EXTRACT_FLANKING`, `ITERATIVE_SEARCH`, and `CLUSTER_REGIONS` come straight from `params.*` — so any CLI override you pass for a preset-affected param (e.g. `--goi_relative_evalue 1e-4`) is honoured end-to-end. Disable this when you're benchmarking parameter sensitivity, or when you want the static `-profile preset_X` startup path to be the only place presets take effect. NORMALIZE_QUERY's `min_query_length` is the one preset-affected param that this toggle can't touch — it runs before LOCATE_GENE, so for short queries you still need `--min_query_length 20` or `-profile preset_short_peptide` at launch.

### `preset_override`
**Type:** String | **Default:** `''` (empty — let auto-detection choose)

Pin a specific preset and skip the hit-profile decision. Valid values: `preset_short_peptide`, `preset_tandem_family`, `preset_single_copy`, `preset_paralog_discrimination`. With `auto_apply_preset = true`, the pinned preset's overrides are applied at runtime to the downstream steps. With `auto_apply_preset = false`, this only affects what gets written to `auto_preset.json` — the params themselves come from `-profile preset_X` (startup) or the CLI. Useful when you already know the query type and don't want SynVoy to second-guess based on a noisy hit distribution.

### `allow_unknown_species`
**Type:** Boolean | **Default:** `false`

By default the pipeline aborts at startup if `--home_species` is empty AND can't be inferred from the genome filename — that empty value used to cascade silently into `phylo_sort.py` and `auto_params`, producing kingdom=Unknown and broken phylo distances (the §1i incident on the luciferase run). Set to true to opt out of that gate when you genuinely don't know the species (e.g. a draft assembly with a meaningless stem). The pipeline then proceeds with a loud `WARN` and degraded phylogenetic ordering. Always prefer fixing `--home_species`; this flag exists for the rare case where you can't.

### `disable_strong_synteny_rescue`
**Type:** Boolean | **Default:** `false`

§1e rescue pass: when a syntenic block has ≥5 HIGH-confidence flanking genes but no GOI annotation (a "right neighbourhood, missed GOI" case), SynVoy runs a relaxed-stringency `miniprot` pass over the block window to emit a LOW-confidence rescue GOI model (`EvidenceType=relaxed_miniprot_rescue`). Set true to skip this — useful when you're already running many loci and don't want extra LOW-confidence noise in the report.

### `strong_synteny_rescue_miniprot_outc`
**Type:** Float | **Default:** `0.05` (5 % query coverage)

Minimum miniprot query-coverage fraction for a rescue model to be emitted. Lower = more permissive (rescues more sketchy fragments). Raise to `0.2` if you only want rescue models that cover at least a fifth of the query.

### `strong_synteny_rescue_window_pad`
**Type:** Integer (bp) | **Default:** `2000`

Bp padded each side of the syntenic block coordinates before extracting the window for the rescue miniprot run. Wider = catches GOI bodies that start just outside the strong-flanking span, at the cost of larger search windows.

### `strong_synteny_rescue_miniprot_timeout`
**Type:** Integer (seconds) | **Default:** `120`

Per-block timeout for the rescue miniprot invocation. Bumps to `300` are reasonable on long blocks; keep it short to avoid stalling the report on pathological blocks.

### `disable_paralog_check`
**Type:** Boolean | **Default:** `false`

§1j Phase B reciprocal-best paralog check: for paralog-family queries (multi-sequence home FASTA like TP53 + TP63 + TP73), each recovered GOI target protein is Smith–Waterman-aligned against every home paralog with `parasail`. The per-call best/runner-up bitscores feed `generate_report`'s `paralog_confusion` flag — when a call's best home paralog differs from the *modal* best at the same (genome, locus) with a sufficient bitscore gap, the report flags it. Single-paralog queries are automatic no-ops (the alignment step is skipped). Set true to skip even on multi-paralog runs (useful for short-peptide runs where the paralog list is irrelevant).

> **`disable_paralog_check` is a no-op today, either way.** The check compares each call against *several* home query sequences, but the query that reaches it is always the single normalised GOI protein (only the first record of `--query` is used). The job it was meant to do — "is this call really my gene, or its paralog?" — is done by locus ownership (`--disable_locus_ownership`), which builds its paralog panel from the home proteome.

### `paralog_confusion_min_gap`
**Type:** Float | **Default:** `5.0`

Minimum bitscore-gap (best vs second-best paralog) before a paralog-confusion flag fires. The default of 5 keeps borderline ties off the report. Raise to `15`/`20` for stricter flagging (only clear-cut paralog reassignments), lower to `2` if you want every disagreeing call surfaced.

### `disable_phylo_placement`
**Type:** Boolean | **Default:** `false`

Phylogenetic placement check (advisory). For every adjudicated GOI call SynVoy compares the call's sequence divergence from the query with its species' taxonomic distance from the home species, fits the trend across species, and labels a call that sits anomalously far for its species `phylo_discordant` (others `phylo_concordant`, or `insufficient_data` when too few species have calls). The idea: a paralog that arose before the speciation is more diverged than an ortholog in the same species. The verdict and its z-score are written on each record in `synvoy_report.json` (`goi_dedup.records[].phylo_verdict`, `phylo_z`) and summarised under `self_consistency.phylo_placement_summary`.

It does not change any confidence by default, and it should be read as a hint: divergence comes from a local alignment, species distance is a taxonomic rank difference rather than time, and the test assumes comparable rates across lineages — which fast-evolving genes violate. On the 10-family benchmark every `phylo_discordant` record was a real ortholog. Set `true` to skip the step.

### `phylo_placement_min_calls`
**Type:** Integer | **Default:** `4`

Number of species with a call needed before the trend is fitted. Below it every call is `insufficient_data`.

### `phylo_placement_z_threshold`
**Type:** Float | **Default:** `2.0`

Residual z-score at or above which a call is `phylo_discordant`.

### `phylo_placement_promote`
**Type:** Boolean | **Default:** `false`

Promote an `AMBIGUOUS` call that is `phylo_concordant` to `MEDIUM`. Leave it off unless a run on your own data has shown the verdict is trustworthy.

### `disable_wavefront`
**Type:** Boolean | **Default:** `false`

A measurement control, not a performance option. Normally targets are searched closest-first in waves, and each wave's recovered orthologs are added to the query database for the next. `true` puts every genome in one parallel wave so the database never grows. Run a set both ways to see whether the wavefront changes the result. On the melittin benchmark it did not; it has not yet been tested on a set where the farthest target is out of reach of the original query.

### `disable_relative_query`
**Type:** Boolean | **Default:** `false`

Inside a window the gene of interest is searched with the **home protein**, and in addition with the gene model of the genome searched most recently (the closest relative in search order; of its models the one that aligns best to the home protein). A model found only with the relative is measured again against the home protein: identity, coverage and confidence refer to the home protein, the id carries `_r<n>_` instead of `_l<n>_`, and `SearchQuery` / `SearchQueryIdentity` name the relative and the identity to it. `true` searches with the home protein only.

Until 2026-10-08 a window was searched with ONE sequence per gene, the longest of the home protein and all earlier models. A model three residues longer than the home protein then replaced it for the rest of the run, and identity and coverage of every later call referred to that model. Calls of later genomes can therefore show a lower identity than in older runs.

### `anchor_max_places`
**Type:** Integer | **Default:** `0` (= every place)

Each flanking gene seeds blocks only at its N best places in a target genome. A place is a run of hits of that gene within its own gap limit; the best places are those with the largest summed bit score. The ortholog of a flanking gene is one place; its further places are relatives and shared domains, and they are what puts two anchors near each other by chance, which is how most searched blocks arise. The hits that are left out still count inside a window (flanking support, flanking gene models); only the seeding of blocks changes.

Measured by replaying the saved genome-wide hits of 31 genome searches (yeasts, flies, bees, ants, plants): with N = 5 the windows searched fall from 512 to 121 and the DNA searched from 130 to 25 Mb, and of the 29 HIGH, MEDIUM or ambiguous calls that lie in a searched block only one is no longer covered, a relative of the gene that the report labels `paralog_not_goi` (N = 10: 166 windows; N = 3: 106). Five test cases were run with N = 5 and gave the same ortholog calls as without it:

| Case | Windows without / with N = 5 | Run time |
|---|---|---|
| yeast STE2 in three yeasts | 101 / 18 | 4.0 / 1.2 min |
| *Arabidopsis* FT in three Brassicaceae and papaya | 238 / 9 | 12.5 / 4.4 min |
| firefly luciferase in fruit fly and mosquito | 295 / 67 | 23 / 9 min |
| fly Adh in three *Drosophila* and a mosquito | 45 / 23 | 4.6 / 3.8 min |
| apyrase in seven ants and a hornet | 35 / 35 | unchanged |

The limit matters where flanking genes have many places: about 2 per gene in an ant genome (nothing changes), 7 to 9 in a yeast, and 92 hits per flanking gene in *Arabidopsis lyrata*, where the default search reaches the cap of 80 blocks per genome and cuts out a fifth of the genome for one home locus. **For plant genomes, and whenever the log says `Capping blocks`, use `--anchor_max_places 5`.**

One side effect. Where a target genome has no conserved neighbourhood, its blocks with two anchors are chance blocks, and the limit removes them all. The search then falls back to every block with one anchor, as it always does when no block passes `min_block_genes` (up to `max_blocks_per_genome`). Fruit fly against mosquito: 3 windows without the limit, 19 with it.

It is off by default because it has not been run on the large families of the benchmark.

### `disable_anchor_resolution`
**Type:** Boolean | **Default:** `false`

A target gene is the ortholog of at most one home gene. Where hits of different flanking genes cover the same target DNA (overlap of at least half of the shorter hit: two kinases, a tandem pair), only the best-scoring one counts as an anchor. Before, one target gene could count as two flanking genes, which let a place with a single gene pass the "at least two flanking genes" test of a block and raised the flanking support of calls. `true` counts every hit, as before 2026-10-08.

### `goi_hull_cluster_max_gap` / `goi_hull_window_pad` / `goi_hull_max_window`
**Type:** `auto` or Integer (bp) | **Default:** `auto`

Sizes of the hull rescue (a miniprot search of the gene across the whole flanking neighbourhood of a sequence, run when no HIGH model lies there). `auto` takes them from the distances of the search for that genome: flanking genes form one cluster up to the bridge limit, the hull is widened by the GOI margin, and a hull longer than twice the home neighbourhood (times the expansion) is skipped. The former fixed values (3 Mb, 100 kb, 20 Mb) were set on a mammalian chromosome and are used when a result GFF carries no distances. The rescue is skipped when a HIGH model already lies in the DNA it would search (before: in the hull without its margin, which modelled a gene just outside the outermost flanking gene a second time).

### `synteny_bridge_two_sided` / `synteny_bridge_max_per_block`
**Type:** Boolean / Integer | **Defaults:** `false` / `2`

The default gap-bridging rule counts the anchored flanking genes on the left of a gap only, so the same evidence bridges or not depending on where the gap falls. `synteny_bridge_two_sided true` weighs both sides together. It is **off pending measurement**: larger blocks raise the flanking support that some confidence rules read, so it can promote more calls. `synteny_bridge_max_per_block` caps how many gaps one block may chain (`0` = unlimited).

### `decoy_mode`
**Type:** Boolean | **Default:** `false`

Calibration of the labels without a truth set. In every search window the reversed home protein and two shuffles of it (seeded by the query id, so a rerun repeats them) are searched like the home protein, modelled, labelled and given the same E-values and specificity terms. They are chosen among themselves only (never against real candidates), never seed the next wave and never reach the result GFFs, the report or the figures; they go to `<outdir>/decoys/<locus>/`. Every decoy call is false, so the share of decoy cases that reach a label or an E-value threshold is that rule's false-call rate on this data. Costs about four times the window search.

### `report_nonsyntenic_candidates`
**Type:** Boolean | **Default:** `true`

Record GOI hits that clear the identity / E-value / length bar but are refused by a synteny gate, with the gate that refused them and the distance to the nearest flanking anchor. They are listed in `synvoy_report.json` under `rejected_candidates`. They are not orthology calls; they exist so that a genome holding a strong unplaceable hit is not read as a genome with no hit.

### `quiet`
**Type:** Boolean | **Default:** `true`

Compact console output: the per-step RUN/SKIP lines are dropped. `false` prints them.

---

## 20. Removed parameters

Removed on 2026-10-09 (the two figure switches on 2026-10-13, the seed and context-promotion options on 2026-10-14, the parameters of the identity-based GOI label on 2026-10-15). Passing one of the first group stops the run at start with the replacement named; one of the second group prints a warning and is ignored.

**Replaced**

| Parameter | Use instead |
|---|---|
| `max_flanking_goi_similarity` | `goi_relative_evalue` (an E-value, not a percentage) |
| `enable_name_locus`, `name_locus_min_identity`, `name_locus_pad` | `home_locus_mode` (`auto` already makes the annotated home gene the primary locus; `gene_only` = that locus alone). `--home_goi_gene SYMBOL` still names the gene |
| `pub_svg` | Nothing: the publication SVG of the ribbon plot was removed. The figure for print is `*_anchor_grid_print.svg`; `*_synteny_plot_view.svg` is the ribbon plot as a picture |
| `enable_matrix_plot` | Nothing: the matrix figure was removed. The anchor grid (`*_anchor_grid.svg`) shows the same and more |
| `seed_on_flanking_support` | Nothing (2026-10-14): seeds follow one rule, the HIGH calls of a query in a genome, else its MEDIUM calls, each when no other home protein competes for it (ALGORITHM.md §7.12) |
| `classify_context_min_collinear` | Nothing (2026-10-14): the context promotion was removed. A weak match in the slot of a block with significant gene order is labelled AMBIGUOUS with `Basis=gene_order_only` ("needs curation", never counted), always on. The promotion was on by default for one day; reversed and shuffled home proteins met its rule in 30 of 181 Candida decoy cases |
| `classify_high_min_identity`, `classify_medium_min_identity`, `classify_tandem_min_identity`, `classify_tandem_min_qcov` | Nothing (2026-10-15): one judge labels every GOI call, tandem copies included, from the E-value in the expected place, specificity against the home proteome and completeness ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)); identity does not decide a label. Since step 4 of the label rework (2026-10-14) these thresholds shaped only a provisional label that the judge overwrote |
| `classify_high_min_collinear` | Nothing (2026-10-15): the order of the anchors enters the label through the block's E_syn (a block is significant at E_syn ≤ 0.1, and only a significant block gives a call its slot or core as the place of its E-value) |
| `disable_ambiguous_tier` | Nothing (2026-10-15): AMBIGUOUS ("needs curation", never counted) is part of the one label, with `Basis=not_specific` or `gene_order_only`, and cannot be switched off |
| `disable_distance_autotune` | Nothing (2026-10-15): the per-genome relaxation of the identity thresholds (§1g) went with the thresholds. The judge computes each call's E-value in the expected place of its own genome. The median flanking identity it rested on is still logged per genome |

**Had no effect** (declared and accepted, read by nothing)

| Parameter | What it was meant for |
|---|---|
| `prefer_large_genes` | prefer long genes as anchors |
| `min_gene_identity` | identity floor for flanking hits (the floor in use is `min_hit_identity`) |
| `aug_dedup_bin_bp` | de-duplication of window hits |
| `max_consecutive_empty_blocks` | early stop of the search in a genome |
| `gap_min_size`, `gap_evalue`, `gap_min_alnlen` | a gap-filling search that was replaced by one miniprot run per candidate |
| `keep_intermediate` | `intermediate/` is always written |
| `max_retries` | retries are set in the `process {}` block of `nextflow.config` |
| `multi_profile`, `multi_profile_max_jobs` | a run with several parameter sets, never wired |
| `locus_ownership_synteny_window` | read by nothing |
| `min_exon_query_cov` | a coverage floor for exon hits; the line that applied it was commented out in March 2026 |
| `pub_width_mm`, `pub_palette` | read by nothing (the print figure has a fixed size, see OUTPUT.md) |
| `seed_flanking_min_count`, `seed_flanking_min_qcov` | floors of `seed_on_flanking_support` (2026-10-14); they acted only when that was switched on |
| `classify_context_min_identity`, `classify_context_min_qcov` | floors of the context promotion (2026-10-14), removed with it |
| `classify_fallback_strong_min_identity_floor` | identity floor of the strong-flanking route by which a hit chain reached MEDIUM (2026-10-15; that route went with the identity-based label) |
| `strict_goi_family` | demotion of hit-chain calls whose target annotation lacks a family name (2026-10-15). The family-name check is annotation only now: `--goi_family_tokens` writes `GoiFamilyConsistent` / `GoiFamilyReason` and changes no label |
| `distance_autotune_close_pct`, `distance_autotune_far_pct`, `distance_autotune_max_relax`, `distance_autotune_min_flanking`, `distance_autotune_high_floor`, `distance_autotune_medium_floor` | settings of the per-genome relaxation of the identity thresholds (2026-10-15), removed with it |

`plot_synteny.py` no longer accepts `--target_names`, `--max_legend_entries`, `--ribbon_alpha_dense`, `--pub_width`, `--pub_palette` (all were parsed and never read) and, since 2026-10-13, `--pub_svg` (its figure was removed). A run that passes one of them in `--plot_extra_args` ends at PLOT_SYNTENY with `unrecognized arguments`; drop the option and `-resume`. `plot_synteny_matrix.py` no longer exists.

Figures removed on 2026-10-13: the anchor-positions figure (`*_anchor_positions.html` / `.svg`; it repeated the anchor grid and the gene-position map), the publication SVG of the ribbon plot (`*_synteny_plot.svg`) and the matrix (`*_synteny_matrix.svg`). The figures of a run are now the anchor grid with its print version, the ribbon plot, the gene-position map and the tree.
