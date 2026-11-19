# SynVoy: the algorithm

This document describes what SynVoy 1.0.1 computes, step by step, with the rule and the
default value at every decision. It is written from the code. Code locations are given as
`file` → `function`, because line numbers drift.

It does not explain how to run the pipeline ([USAGE.md](USAGE.md)), what the output
files contain ([OUTPUT.md](OUTPUT.md)) or what every parameter does
([PARAMETERS.md](PARAMETERS.md)).

Contents: [1 Idea](#1-the-idea) · [2 Terms](#2-terms) · [3 Overview](#3-overview) ·
[4 Inputs](#4-inputs) · [5 Home locus](#5-the-home-locus) ·
[6 Target order and waves](#6-target-order-quality-check-and-waves) ·
[7 Search in one genome](#7-the-search-in-one-target-genome) ·
[8 Region scores](#8-region-scores) · [9 Rescue passes](#9-rescue-passes) ·
[10 Adjudication and report](#10-adjudication-and-the-report) · [11 Tree](#11-tree) ·
[12 Figures](#12-figures) · [13 Reproducibility](#13-reproducibility) ·
[14 Properties and limits](#14-properties-and-known-limits) ·
[15 Defaults](#15-default-values)

---

## 1. The idea

A sequence search over a whole genome needs strict thresholds, so it misses a gene
that is short or has diverged. The genes next to it are often better conserved.
SynVoy therefore:

1. finds the gene in one annotated genome and takes its neighbours,
2. finds where those neighbours lie, in the same order, in each other genome,
3. searches for the gene only inside that stretch, with thresholds that would be
   useless genome-wide,
4. labels each gene model it builds by how well the sequence and the neighbourhood
   support it.

Genomes are searched from the closest relative outwards. A gene found with the label HIGH
in a close relative (MEDIUM when there is no HIGH call) is added to the query set for the
more distant ones, as long as no other home protein competes for it
([§7.12](#712-what-is-added-to-the-query-set)).

## 2. Terms

| Term | Meaning |
|---|---|
| **GOI** | Gene of interest: the query. |
| **Home genome** | The annotated genome the query comes from. |
| **Target genome** | A genome in which the GOI is searched. |
| **Home locus** | One place in the home genome where the query aligns. A gene family gives several. Everything from [§5.7](#57-flanking-genes) on runs once per home locus. |
| **Flanking gene** | A protein-coding neighbour of the GOI in the home genome, used as an anchor. |
| **Flank walk** | The walk outwards from the home gene along its scaffold that picks the flanking genes ([§5.7](#57-flanking-genes)). |
| **Home rank** | The position of a flanking gene in home genomic order (0, 1, 2 …). |
| **Gene pitch** | The mean step from one flanking gene of a home locus to the next: the span of the flanking genes divided by their number minus one. The unit of the search's distances. |
| **Expansion** | How much larger or smaller the home neighbourhood is in one target genome ([Distances](#distances-in-gene-units)). |
| **Hit** | One local alignment of a query protein to genomic DNA. |
| **Hit interval** | The query's hits at home that form one copy of its gene: chained while each continues the query, within the intron bound and one annotated gene ([§5.2](#52-locating-the-query)). |
| **Block** | A run of flanking-gene hits on one target scaffold: a candidate neighbourhood. |
| **Collinear run** | The longest subsequence of a block's flanking genes whose home ranks only rise or only fall along the target scaffold. |
| **Model** | A gene structure (exons, protein) built in a target genome. |
| **Query coverage (qcov)** | The share of the query that a model covers. For a model built from the home gene's exons, and for a tandem copy, it is the span from the first to the last aligned query residue, divided by the query length. For a model built from a chain of hits, and for a locus of the hit profile, it is the union of the aligned query intervals. |
| **Label** | The confidence of a GOI model: HIGH, MEDIUM, AMBIGUOUS or LOW ([§7.8](#78-confidence)). |
| **Judge** | The one function that labels every GOI candidate, and the order in which candidates are ranked (`synvoy_label.judge`, `rank_key`; [§7.8](#78-confidence)). |
| **ModelStatus** | `complete`, `partial` or `fragment`: how much of the query a model covers. It is not the label and does not change it ([§7.8](#78-confidence)). |
| **Settings map** | The 32 parameters that travel from the settings step to the later steps: the 25 that a preset module may change and 7 fixed ones ([§5.3](#53-parameter-preset-and-preset-modules)). |
| **Wave** | A group of target genomes searched with the same query set. |

## 3. Overview

| # | Step | Nextflow process | Script | Result |
|---|---|---|---|---|
| 1 | Resolve the input | `RESOLVE_GENE_INPUT`, `FETCH_HOME_GENOME`, `FETCH_RELATED_GENOMES` (easy mode); `STAGE_GENOMES` (pro mode) | `resolve_gene_input.py`, `fetch_home_genome.py`, `fetch_related_genomes.py` | query, home genome + GFF, target genomes |
| 2 | Normalise the query | `NORMALIZE_QUERY` | `normalize_query.py` | one protein sequence |
| 3 | Locate the query at home | `LOCATE_GENE` | MMseqs2, tblastn, `merge_hits.py`, `profile_hits.py`, `auto_select_preset.py` | hit intervals with bit scores, hit profile |
| 4 | Choose the settings: the preset modules (the old preset under `--modules legacy`); runs after step 5, because it reads the home gene model | `HOME_CONTEXT`, `RESOLVE_EFFECTIVE_PARAMS` | `home_context.py`, `select_modules.py`, `resolve_effective_params.py` | the settings map, `effective_params.json` |
| 5 | Model the home gene | `ANNOTATE_GOI` | `annotate_goi_exons.py` | GOI protein, exons, largest intron |
| 6 | Split into home loci | `SPLIT_LOCI` | `split_loci.py` | up to 5 loci, the annotated home gene first |
| 7 | Build the home proteins and pick flanking genes | `PREPARE_HOME_PROTEOME`, `BORROW_ANNOTATIONS`, `EXTRACT_FLANKING`, `PREPARE_INITIAL_DB` | `gff_to_faa.py`, `prodigal_on_regions.py`, `borrow_annotations.py`, `extract_flanking_genes.py` | home proteome, initial query set |
| 8 | Order and check targets | `PHYLO_SORT`, `ASSESS_GENOME_QUALITY`, `FILTER_SORTED_GENOMES` | `phylo_sort.py`, `assess_genome_quality.py` | targets, closest first |
| 9 | Search | `ITERATIVE_SEARCH` | `iterative_search_runner.py` | gene models per genome |
| 10 | Score regions | `CLUSTER_REGIONS` | `cluster_grs.py` | ranked regions per genome |
| 11 | Rescue | `RESCUE_STRONG_SYNTENY`, `RESCUE_GOI_HULL` | `rescue_strong_synteny.py`, `rescue_goi_hull.py` | extra GOI models |
| 12 | Adjudicate and report | `BUILD_HOME_PARALOG_PANEL`, `ASSIGN_LOCUS_OWNERSHIP`, `PHYLO_PLACEMENT_CHECK`, `GENERATE_REPORT` | `build_home_paralog_panel.py`, `reciprocal_best_paralog_check.py`, `phylo_placement_check.py`, `generate_report.py` | `synvoy_report.json` |
| 13 | Tree and figures | `COMPUTE_TREE`, `PLOT_SYNTENY` | `compute_tree.py`, `plot_synteny.py` | tree (a placeholder when `skip_tree` is true, as under the `low_mem` profile that `run_synvoy.sh` adds when no profile is given), HTML/SVG figures |

The home proteins come from the home GFF. Without one, the home genes are predicted
(Augustus, or Prodigal) and the genes of an annotated target genome are transferred onto
the home genome ([§5.7](#57-flanking-genes)).

Nextflow runs a step as soon as its inputs exist, so steps without a dependency on each
other run at the same time. Step 10 and the tree read the gene models of step 9. Step 11
reads both: the hull rescue reads the models of step 9, and the strong-synteny rescue reads
the regions of step 10. Only the report (step 12) reads the rescue models of step 11, so
they appear in the report but not in the tree or the figures
([§14](#14-properties-and-known-limits)). The figures are drawn last, from the models of
step 9, the regions of step 10 and the tree; they read the report and leave out the calls
it drops.

## 4. Inputs

**Pro mode** takes local files: the query (a protein, or a DNA sequence that is translated,
[§5.1](#51-query-normalisation)), the home genome, optionally its GFF (without one, the home
genes are predicted and the genes of an annotated target are transferred,
[§5.7](#57-flanking-genes)), and the target genomes (optionally their GFFs).

**Easy mode** takes an accession or a sequence and fetches the rest.

1. `resolve_gene_input.py` → `detect_input_type` classifies the input as a UniProt
   accession, an NCBI protein accession (version optional), a FASTA file or a gene
   symbol, fetches the protein and records the species and the gene symbol.
2. `fetch_home_genome.py` downloads the reference genome of that species, else any
   assembly of it, else the genome of its closest relative (with a warning that a
   surrogate is used). It takes the annotation when the assembly has one.
3. `fetch_related_genomes.py` chooses the targets by walking up the taxonomy:
   - **Levels.** The genus of the home species (when it resolves), then its family,
     order and class.
   - **How many.** With `max_genomes = 0`, three per level, at most 20. Each level
     gets `max(2, max_genomes // levels)`; places a level cannot fill are carried to
     the next level.
   - **Candidates.** All assemblies under the level's taxon. Each species is used at
     most once, and the home species is excluded (compared without strain names).
   - **Ranking** (`assembly_rank_tuple`, mode `hybrid`): chromosome-level before
     complete before scaffold before contig; then reference genome before
     representative before other RefSeq before GenBank; then fewer contigs, fewer
     scaffolds, larger contig N50, larger scaffold N50. It chooses the assembly of a
     species and, in that order, the species that fill a level's places.
   - **Quality gate** (`is_bad_quality`): chromosome-level and complete assemblies
     always pass. Any other assembly is low quality when its best N50 is below 5 kb,
     when a known contig or scaffold count exceeds 500,000, or when it is contig-level
     and has no quality figures at all. With `bad_quality_policy = auto` such an
     assembly is dropped from the automatic search, and kept with a warning when the
     user named the species (`--target_species`) or it is the home assembly. The gate
     runs on each candidate as it is about to fill a place, so the next candidate of
     the level takes the place of a dropped one.

**Genetic codes.** Each genome is translated with its NCBI table: `--home_genetic_code`
for the home side, `--target_genetic_codes` (a TSV of genome and table) or
`--default_target_genetic_code` for each target. All default to table 1, and then no
engine gets a genetic-code option. The table is a property of the process
(`sequence_utils.set_genetic_code`): the search sets the target's table at the start of
each genome, and `translate`, the stop codons of the gene-model walks, the frame guard
(§7.8) and the options of MMseqs2, tblastn, miniprot and Smith-Waterman follow it.
Supported are tables 1–6, 9–16 and 21–31, the ones every engine translates with.

## 5. The home locus

### 5.1 Query normalisation

`normalize_query.py` upper-cases the sequence and removes whitespace, digits, gap
characters and a trailing `*`. Internal `*` are removed with a warning. Other
non-sequence characters are an error. A file with several sequences is searched with the
first one only, with a warning. A nucleotide query (DNA or RNA) is translated in six
frames and the longest open reading frame is taken. A protein shorter than
`min_query_length` (20) stops the run; one of 20 to 39 residues runs with a warning,
because an alignment that short is not significant on its own in another genome and a
call then rests on the neighbourhood.

The query's id is read once, here. The known forms of public databases (`sp|P12345|NAME`
and `P12345|NAME`) are read as the accession (`P12345`). An id that still holds a `|`, a
character that GFF reserves (`,;=&%`), or the form SynVoy gives the pieces of a protein
(`|exon_`, `|frag_`, `|slide_`) is an error. Every file SynVoy writes after this step is
read by the first word of each header ([§5.7](#57-flanking-genes)).

### 5.2 Locating the query

`LOCATE_GENE` searches the home genome twice and pools the result:

- MMseqs2 `easy-search`, translated search, sensitivity `mmseqs_sensitivity` (9.5; 7.0
  under the `low_mem` profile, which `run_synvoy.sh` adds when no `-profile` is given;
  8.0 under `laptop_safe`), E ≤ `search_evalue` (0.01);
- `tblastn` with default settings.

`merge_hits.py` chains the hits into **hit intervals**, one per copy of the gene. It
keeps hits with E ≤ `search_evalue` (0.01), one per stretch of DNA (of two hits that
overlap by more than half of the shorter one the higher bit score stays: both programs
report every exon), and reads them per scaffold and strand in coding order. A hit
continues a chain when

- it continues the query: it starts at most 20 residues before the chain's last hit
  ends, and ends after it. An exon's alignment often runs into the intron before it,
  whose translation matches the end of the previous exon's part of the query; a second
  copy of the gene repeats the query and starts a chain of its own, however close it
  lies;
- it lies at most `home_chain_max_intron` after that hit. `auto` is the largest intron
  of the annotated genes the hits overlap × `goi_fallback_intron_margin` (5), at least
  `goi_fallback_intron_floor` (10 kb) and at most 1.1 Mb; the floor without an
  annotation (the home gene model's intron, which the search uses, is not known yet);
- the two are not in two different annotated genes, and no annotated gene lies
  entirely between them (genes nested in an intron of the gene do not stop it).

A hit that several chains could take goes to the one whose last query end is closest
to its query start. A hit interval covers its hits and the other program's hits of the
same exons and carries their lowest E-value and highest bit score.

### 5.3 Parameter preset and preset modules

`profile_hits.py` summarises the hits with E ≤ `search_evalue`. Hits that cover the
same DNA are reduced to the best one (both programs report every exon). Hits with at
most `locus_gap_genes` (`auto` = ⌊`n_flanking_genes`/2⌋, 5) annotated genes between them
form one locus, as in §5.4. A locus is scored as **one copy of the gene**: the best total
bit score of hits that cover the query at most once (`one_copy_score`). A locus is
*strong* when its score is at least half of the top locus's; `bit_ratio` is the top score
divided by the second; `primary_copies` is the number of aligned query residues at the
top locus divided by the query length (about 1 for one gene with any number of exons,
about 2 for a tandem pair). Each locus also gets its **coverage**, the share of the query
covered by the union of its hits: `secondary_full_length` says that a strong locus other
than the top one covers at least half of the query (a full-length relative),
`shared_domain_only` that there are two or more strong loci and every secondary one covers
less (relatives that share a domain). The profile carries the query's own features too
(`query_features.py`): length, self bit score, the share of low-complexity residues
(12-residue windows of Shannon entropy ≤ 2.2 bits, an approximation of BLAST's SEG), the
share in repeated 4-mers, an N-terminal signal peptide (n-, h- and c-region with the
−1/−3 rule; a heuristic, not SignalP), a C-terminal GPI-anchor signal, and cysteines.

**The old chooser.** `auto_select_preset.py` → `select_preset` applies the first rule
that matches (its preset is applied under `--modules legacy` only; under the default
`auto` it is recorded as `legacy_preset`):

| Order | Condition | Preset |
|---|---|---|
| 1 | no hits | `preset_single_copy` |
| 2 | query shorter than 100 aa | `preset_short_peptide` |
| 3 | at least 2 strong loci and `bit_ratio` < 2 (relatives elsewhere in the genome) | `preset_paralog_discrimination` |
| 4 | `primary_copies` ≥ 1.5 (several copies at the home locus) | `preset_tandem_family` |
| 5 | otherwise | `preset_single_copy` |

Rule 3 is in effect "at least 2 strong loci": a strong locus scores at least half of the
top one, so the ratio of the top two is at most 2. When rule 3 or 4 decides, more than 5
strong loci or more than 10 loci add a large-family warning; the other outcomes and
`--preset_override` never warn. The number of alignment rows decides nothing: a row is one
exon in one search program.

What each preset changes (`resolve_effective_params.py` → `PRESET_OVERRIDES`):

| Preset | Changes from the defaults |
|---|---|
| `preset_short_peptide` | the ModelStatus of gene models, not their label: fragment below 0.3 query coverage (default 0.4), complete from 0.5 (default 0.7) |
| `preset_paralog_discrimination` | fewer candidate regions per genome (the regions that the figures draw and the strong-synteny rescue reads): a score floor of 0.45 × the best region's score (default 0.30) and of at least 0.08 (default 0.03), at most 4 regions (default 6) |
| `preset_tandem_family` | a neighbour is a relative of the GOI up to E = 1e-3 (default 1e-5; §5.7) |
| `preset_single_copy` | nothing |

**Preset modules** (`select_modules.py`). A gene is not one kind: it can be short,
cysteine-rich, secreted and in a tandem array at once, and the first-match rule gives it
one preset. A module is one kind of gene; each tests its own trigger, they stack, and none
comes first. The triggers read the hit profile with the query features, the home gene
model (`goi_info.json`) and the home neighbourhood (`home_context.py`): the family members
of the query among the ±`n_flanking_genes` protein-coding home genes next to it, tested
with the label's family test (blastp with composition-based statistics at E ≤ 1e-3; found
only with a home GFF and blastp), and the distance of the home gene to the scaffold end.

| Id | Module | Trigger | Parameters it sets (values in the rules table below) |
|---|---|---|---|
| A1 | `local_family` | `primary_copies` ≥ 1.5, or a family member among the ±`n_flanking_genes` home genes next to the home gene | `goi_relative_evalue`, `goi_candidate_regions_per_query`, `goi_max_models_per_locus`, `goi_margin_genes`, `sw_max_alignments_per_piece`, `report_copies_in_slot`, `locus_ownership_max_genes_per_locus`, `max_goi_per_genome` |
| A2 | `dispersed_family` | ≥ 2 strong loci and `secondary_full_length` | `adaptive_score_floor_frac`, `adaptive_score_floor_abs`, `adaptive_max_regions`, `paralog_panel_max_family` |
| A3 | `shared_domain` | `shared_domain_only` | none (a flag; A2 does not fire) |
| B1 | `short_query` | query < 100 aa (its headroom, the self bit score over the bits an alignment needs for E ≤ 1e-3 in a 100-kb window, is logged) | `classify_complete_min_qcov`, `classify_fragment_max_qcov`, `window_tblastn_word_size`, `window_tblastn_comp_stats`, `window_mmseqs_min_length`, `miniprot_short_orf_max_len`, `aug_relaxed_length_min`, `fallback_short_min_aln_aa`, `dispersed_goi_min_alnlen`, `goi_hull_min_coverage`; also the LOCATE_GENE retry (below) |
| B2 | `low_complexity` | low-complexity share ≥ 0.10 or repeated 4-mers ≥ 0.25 | `goi_relative_statistic`, `report_low_complexity` |
| B3 | `signal_anchor` | a signal peptide or a GPI tail | `core_completeness` |
| B4 | `cys_framework` | ≥ 6 Cys in a 90-residue window | none (advisory) |
| C1 | `exon_structure` | a home gene of ≥ 2 exons (tandem copies of the home model not counted) | `sw_max_alignments_per_piece`, `dispersed_goi_min_alnlen` |
| C2 | `home_context` | the home gene < 40 kb from the end of a scaffold of ≥ 1 Mb | none (advisory) |

The trigger thresholds come from the preset chooser (100 aa, 1.5 copies) or were set by
judgement; none has been calibrated on data.

**One rule per parameter** (`RULES`). Every parameter a module may change has exactly one
rule: its owners (the modules that switch it on), one formula of the run's design
quantities, a cap, and the evidence for it. The rule applies when any of its owners is
applied, so two modules can never disagree about a key and their order never matters. The
design quantities are computed once from the features (`design_profile`):

- `c`, the home copies: the larger of `primary_copies` (rounded half up) and 1 + the
  number of family members next to the home gene (as for A1);
- `R`, the reach of the home array in genes: the largest |rank| of such a family member,
  0 without one;
- `Fam`, the home loci of the family (the loci of the hit profile);
- `Lq`, the query length; `Lx`, the longest real exon of the home gene; `Lpiece` =
  min(`Lq`, `Lx`); `Lcore`, the query without its signal peptide and GPI tail;
- `nx`, the real exons of the home gene; `n` = `n_flanking_genes`; `W` =
  `window_margin_genes`.

| Parameter | What it sets | Default | Modules | Value | Cap |
|---|---|---|---|---|---|
| `goi_relative_evalue` | the E-value up to which a neighbour is a relative of the GOI and no anchor (§5.7) | 1e-5 | A1 | 1e-3 | |
| `goi_candidate_regions_per_query` | GOI loci kept per query in a search window | 2 | A1 | max(2, `c` + 1) | 40 |
| `goi_max_models_per_locus` | gene models per GOI locus (masked miniprot runs) | 10 | A1 | max(10, `c` + 1) | 60 |
| `goi_margin_genes` | the GOI margin, in pitches ([Distances](#distances-in-gene-units)) | 2 | A1 | min(`W`, max(2, `R` + 1)) | |
| `sw_max_alignments_per_piece` | Smith-Waterman alignments per window piece (a second copy or exon in the same frame) | 1 | A1, C1 | 3 if `c` ≥ 2 or `nx` ≥ 2, else 1 | |
| `report_copies_in_slot` | `CopiesInSlot=` and `SlotRank=` on GOI models (§7.13) | false | A1 | on | |
| `locus_ownership_max_genes_per_locus` | home genes per home locus in the ownership panel (§10) | 5 | A1 | max(5, `c` + 1) | 60 |
| `max_goi_per_genome` | GOI models of a genome that steer which neighbourhood a figure draws (§12) | 10 | A1 | max(10, `c` + 1) | 60 |
| `adaptive_score_floor_frac` | candidate regions: score floor as a share of the best region's score (§8) | 0.30 | A2 | 0.45 | |
| `adaptive_score_floor_abs` | candidate regions: absolute score floor (§8) | 0.03 | A2 | 0.08 | |
| `adaptive_max_regions` | candidate regions kept per genome (§8) | 6 | A2 | 4 | |
| `paralog_panel_max_family` | home relatives in the ownership panel (§10) | 15 | A2 | max(15, `Fam` + 1) | 50 |
| `classify_complete_min_qcov` | ModelStatus: complete from this coverage (not the label: HIGH needs 0.7 query coverage, §7.8) | 0.7 | B1 | 0.5 | |
| `classify_fragment_max_qcov` | ModelStatus: fragment below this coverage (the figures hide fragments; the report's piece rule, §10, reads only fragments) | 0.4 | B1 | 0.3 | |
| `window_tblastn_word_size` | tblastn word size in the windows | 3 | B1 | 2 | |
| `window_tblastn_comp_stats` | tblastn composition statistics in the windows (a run with tblastn's defaults still gives `EeffTblastn`) | 2 | B1 | 0 | |
| `window_mmseqs_min_length` | the shortest ORF, in codons, that MMseqs2 searches in the windows | 30 | B1 | 10 | |
| `miniprot_short_orf_max_len` | queries shorter than this (aa) get miniprot `-L 10` in the windows | 60 | B1 | 100 | |
| `aug_relaxed_length_min` | the shortest GOI hit kept in a window, in aligned query residues | 15 | B1 | 10 | |
| `fallback_short_min_aln_aa` | the shortest hit of the short-query fallback, in aligned query residues | 15 | B1 | 10 | |
| `dispersed_goi_min_alnlen` | the shortest hit of the dispersed rescue, in aligned query residues (§7.11) | 50 | B1, C1 | min(50, max(10, round(0.6 × `Lpiece`))) | |
| `goi_hull_min_coverage` | the hull rescue's coverage filter (§9) | 0.5 | B1 | 0.3 | |
| `goi_relative_statistic` | how the flank walk's relative test scores a pair (§5.7) | `sw` | B2 | `blastp_cbs` (blastp with composition-based statistics) | |
| `report_low_complexity` | `LowComplexity=` on GOI models | false | B2 | on | |
| `core_completeness` | completeness for HIGH also on the core of the home protein (`CoreCoverage=`; §7.8) | false | B3 | on | |

Rules that hold for every module:

- **No module changes the label rules.** The judge's significance levels and its
  completeness bar (§7.8) are constants, not parameters. A module can still change a
  call's label, because it changes the evidence the judge reads: it may widen the search,
  tighten the inputs of the specificity test, or measure completeness on the core of the
  protein (B3). A call can therefore be HIGH with B3 and MEDIUM without it. Report labels
  such as `LowComplexity=` change nothing.
- **Keys no rule may write** (`FORBIDDEN_FOR_MODULES`): a switch (a whole step on or off,
  such as `decoy_mode`, a run-level option, or `tandem_copy_parent`), `n_flanking_genes`,
  `expand_goi_similar`, the module mechanism itself, and a value of the genome axis, which
  the search sets per target genome (§7.1).
- **Constraints** (`check_constraints`, after the rules and the values the user set; the
  constraint ids are not module ids): the GOI margin lies within the window margin, where
  the DNA searched ends (constraint C1); with A1, the GOI loci kept per query and the
  models per GOI locus reach `c` + 1 within their caps (a cap that binds is named;
  constraint C3), as do the ownership panel's genes per home locus (C4) and the GOI models
  that steer a figure (C5); a per-hit length floor a rule sets lies between 10 and
  max(10, 0.6 × `Lpiece`) (C7); a value a rule set stays within its cap (C11); the
  ModelStatus fragment bound lies below the complete bound (C12). A violation that the
  rules alone made is a bug and stops the run before the search; one that involves a
  value the user set is a warning, and the user's value is kept. Under
  `core_completeness`, every completeness test of a gene model reads the larger of
  `QueryCoverage` and `CoreCoverage` (§7.8); one function computes this, so it is not
  checked at run time.
- **Precedence**, lowest to highest: shipped default < modules < preset (the chooser's
  under `legacy`, or `--preset_override` in any mode) < a value the user set at launch
  (a `-profile preset_X` included). A value counts as set by the user only when it
  differs from the shipped default in `nextflow.config`. To keep a shipped default
  against a module, remove the module (`--modules auto,-local_family`) or use
  `--modules none` or `--auto_apply_preset false`.

**Modes.** `--modules` selects the mode:

- `auto`, the default (a value with no word in it counts as `auto` too): every module
  whose trigger fires;
- a list (`short_query,local_family`): exactly those modules;
- deltas (`auto,-dispersed_family,+cys_framework`): `auto`, edited;
- `none`: no module;
- `legacy`: the old chooser's preset; the triggers that fire are only recorded.

Each preset corresponds to a module set: `preset_short_peptide` = {B1},
`preset_tandem_family` = {A1}, `preset_paralog_discrimination` = {A2},
`preset_single_copy` = {}. Each module sets its preset's keys to the same values. A1, A2
and B1 also set further keys (listed above), which `legacy` never applies. Under
`legacy`, `modules_applied` names the module that matches the chosen preset, but only the
preset's keys change.

Every mode other than `legacy`, `none` included, also sets `anchor_max_places` per target
genome (§7.1). For a home genome without a GFF, such a mode also chooses the installed
Augustus model nearest the home species, and stops the run when a eukaryotic home would
fall back to Prodigal (§5.7). Under `legacy` the search uses every anchor place, and
Augustus uses its `fly` model. All of this follows `--modules`, even with
`--auto_apply_preset false`.

LOCATE_GENE runs before the modules are resolved, so its retry follows the `--modules`
value itself. When a query under 100 aa has no home hit at E ≤ `search_evalue`, tblastn
runs once more without its low-complexity filter (`-seg no`, recorded as `locate_retry`)
if `--modules` is `auto` (unless it removes `short_query`) or a list that names
`short_query`; never under `legacy` or `none`. A run without home hits stops, so the
retry can only rescue a run.

`--auto_apply_preset false` keeps every setting of the map at its default or at the value
the user set. No module and no preset is applied, not even one named with
`--preset_override`. `--preset_override` applies a preset by hand, in every mode, on top
of the modules (`preset_single_copy` adds nothing).

**Where the values go.** The chosen values travel to the later steps as the settings map
of 32 parameters (`presetDefaultKeys()` in `main.nf`): the 25 a rule may set and 7 fixed
pass-through keys that no rule sets (`PASS_THROUGH_KEYS`: the hit floors `sw_min_score`,
`sw_min_identity`, `min_hit_identity`, `min_hit_length`, the forbidden
`expand_goi_similar`, the annotation input `goi_family_tokens` and the region floor
`adaptive_unique_gene_floor`). So `effective_params.json` (in `intermediate/locate_gene/`)
holds every value of these parameters that the run used; with any mode other than
`legacy` and a home genome without a GFF, it also holds the chosen `augustus_species`. The
flank walk, the search and the region scorer read the map directly; the hull rescue, the
ownership panel, the report (`core_completeness`) and the figures read their keys from
`effective_params.json` (`bin/synvoy_settings.py`, read-only; an option on the script's own
command line wins). The home-proteome step and the flank walk read `augustus_species` from
the map. `effective_params.json` records the mode, the modules applied and why, the
triggers that fire, the source of every value and the features, and the report copies
these as `run_settings`. The same file also holds the design quantities (`profile`), the
derivation of every value a rule, the preset or the user set (`derivation`; for a rule its
formula and inputs, e.g. `goi_candidate_regions_per_query 4 = max(2, c + 1), c = 3 (A1)`),
and the constraint checks (`constraints`); the report's `run_settings` carries these three
as well.

### 5.4 Splitting into home loci

`split_loci.py` runs after the home gene model (§5.6) is built:

1. Sort the hit intervals by position. Two neighbouring intervals on a scaffold
   belong to one locus when at most `locus_gap_genes` (`auto` = ⌊`n_flanking_genes`/2⌋,
   5 at the default) annotated genes lie entirely between them (protein-coding genes when
   the annotation states a biotype, every gene otherwise). On a scaffold without
   annotated genes, and for a home genome without annotation, the limit is 50 kb.
   `--locus_gap <bp>` uses a distance everywhere.
2. **Primary locus.** When the home gene model came from the annotation
   (`method: gff_annotation`), every group of intervals that overlaps that gene is
   replaced by one interval, the coding span of the home gene's chosen transcript (§5.6),
   and this is the primary locus (§5.5). Otherwise the locus with the highest score is
   the primary locus.
3. Each locus is scored as one copy of the gene (`one_copy_score` over the raw hits
   that lie in it, §5.3), not by its best single hit.
4. Keep a further locus only if its score is at least `locus_min_bit_ratio` (0.5) of
   the primary's.
5. Keep at most `max_loci` (5): the primary locus and the highest-scoring others.
6. **Rows of a kept locus.** Each hit interval that an annotated gene covers at least
   half of is replaced by that gene's coding span (the longest coding model, as the home
   proteome takes it); its hit span stays in columns 8 and 9. Intervals on one gene
   become one row. With several rows the row holding the best copy of the gene
   (`one_copy_score` of the hits in it) comes first and is named `best_interval`: the
   flank walk of the locus starts there (§5.7). The primary locus made from the annotated
   home gene stays its one row.

Scoring and filtering (steps 3 to 5) use the hit intervals as chained, before step 6.
More than 10 loci before filtering produce a large-family warning.

### 5.5 The annotated home gene as primary locus

`home_locus_mode` decides which loci are searched:

| Mode | Loci |
|---|---|
| `auto` (default) | the annotated gene that the query is, plus further loci that pass steps 4 and 5 of §5.4 |
| `gene_only` | the annotated gene's locus alone |
| `alignment` | the loci of the alignments alone (§5.4 without step 2) |

Naming the home gene with `--home_goi_gene` sets `gene_only`, whatever `home_locus_mode`
says: only that gene's locus is searched. Without a matching annotated gene, `auto` and
`gene_only` fall back to the alignments (`gene_only` with a warning). The gene is what the
query is; the alignments around it are that gene's exons and, often, hits in a relative
next door.

### 5.6 The home gene model

`annotate_goi_exons.py` builds the GOI's own exon structure, trying in this order:

1. **Annotation by name** (`match_goi_in_gff` → `_search_gff_by_name`). A name given
   with `--home_goi_gene` is looked up first; then the UniProt names of the query (gene
   names, synonyms, the protein name) and the query id, with case variants and the
   `gene-` prefix, are looked up in the home GFF. A name only points at a gene, so each
   match is then compared with the query: up to 50 matches are compared, best name match
   first, and the first gene whose protein agrees with the query is taken
   (Smith-Waterman, identity ≥ 60 % and a score of at least half the query's
   self-score). Of that gene, the transcript that agrees best with the query is used. The
   best other gene compared is recorded as the runner-up.
2. **Annotation by position and sequence** (`_search_gff_by_region_and_sequence`).
   The candidates are the annotated genes that overlap a hit of the locating step with
   E ≤ `search_evalue`, or lie within `gff_search_window` (100 kb) of one. Each of their
   transcripts is compared with the query by the same check as in step 1 (a transcript
   that shares almost no 3-residue words with the query is skipped before the
   alignment). Among the genes that pass, one that overlaps a hit comes before one that
   is only near a hit, then the higher score fraction. When another gene comes within
   0.05 of the winner's score fraction, the winner is still taken, with a warning and
   `home_gene_ambiguous: true`; a gene below 95 % identity is named in a warning (a query
   from another species, or a relative where the query's own gene is not annotated).
   Without parasail the step uses a coarser test and says so in the log: the best share
   of shared 3-residue words, scaled by the length ratio, above 30 %. That test needs
   about 75 to 80 % identity at equal length, so a query from another species or a short
   isoform falls through to miniprot. (A pipeline run without parasail stops at the start
   of the search, [§7.6](#76-search-inside-the-window); this test serves the step run on its own.)
3. **miniprot**, when no annotated gene matches: on each scaffold, the stretch from the
   first to the last hit of the locating step (at any E-value) ± `gap_search_window`
   (50 kb). Scaffolds are tried in the order of the hit tables, and the first that gives
   a model is used.
4. **Query protein only** (`found_by: none`), when miniprot builds no model from the
   hits. A query without any home hit stops the run when the hits are split into loci
   (§5.4).

Cases 1 and 2 take the exons from the GFF's CDS rows (`method: gff_annotation`). The
output is the GOI protein, one sequence per exon, and `goi_info.json` with the gene's
genomic span and its largest intron and with how the gene was found: `found_by`
(`name`, `position`, `miniprot` or `none`; `position_3mer` without parasail), `gene_id`,
`transcript_id`, `identity`, `score_fraction`, `overlaps_hit`, `runner_up_gene_id`,
`runner_up_score_fraction` and `home_gene_ambiguous`. The GOI protein written here, not
the query file, is what the search uses ([§5.8](#58-the-initial-query-set)): for an
annotated gene it is the protein of the chosen transcript.

### 5.7 Flanking genes

`extract_flanking_genes.py` walks once per home locus, from its best interval (the
first row of the locus BED, named `best_interval` when the locus has several rows,
[§5.4](#54-splitting-into-home-loci)). The locus's other intervals are copies or
relatives of the gene: a gene there is no anchor, and it is written into the home BED
as a `GOI_` row without a protein.

1. Take the annotated gene whose centre is closest to the centre of the best interval.
2. Walk outwards on each side. Collect genes until `n_flanking_genes` (10)
   protein-coding genes are found on that side. A side ends earlier only at the end of
   the scaffold, without a home GFF at the end of the stretch SynVoy annotated around the
   locus (below), at `max_flanking_distance` from the centre when that is set (off by
   default), or at the limit: the walk has stepped over more than 5 × `n_flanking_genes`
   genes (50; every gene counts) and has met at least `n_flanking_genes` + R
   protein-coding genes, where R is the number of protein-coding genes on that side of
   the scaffold that overlap a home hit of the query (relatives the walk will skip).
   Genes at the locus's other intervals count towards the 50, but neither towards
   `n_flanking_genes` + R nor towards R. Each side is logged (`[flank-walk]`: anchors,
   genes visited, relatives skipped, what ended it), with a warning when the limit ends a
   side with fewer than `n_flanking_genes` anchors.
3. Skip a gene that is a relative of the GOI: its protein aligns to the query protein
   (as normalised in §5.1) (Smith-Waterman, BLOSUM62, 11/1) with an E-value at or below
   `goi_relative_evalue` (1e-5), computed for the two sequences alone. Such a gene would
   match wherever the GOI's relatives are. The `local_family` module raises the limit to
   1e-3, and the `low_complexity` module scores the pair with blastp and
   composition-based statistics (`goi_relative_statistic blastp_cbs`; §5.3). `0`
   switches the test off for the neighbours; the gene at the GOI's own position is never
   taken.
4. Genes without CDS are kept for the figures but do not count towards the 10.

Five counts follow `n_flanking_genes` = n. Each is a parameter whose default `auto`
derives it from n (`lib/FlankingCounts.groovy`, `lib/LocusGap.groovy`), and a number sets
it directly: the bridge limit `synteny_bridge_max_gap_genes`, 2n pitches (§7.3); the HIGH
flanking models of a strong-synteny region, `strong_synteny_min_flanking`, n / 2 rounded
up (§8); those in the hull of the hull rescue, `goi_hull_min_flanking`, round(0.4 n)
(§9); those at which the report flags strong synteny without a HIGH GOI call,
`strong_synteny_no_goi_min_flanking`, 2n (§10); and the gene gap between home loci,
`locus_gap_genes`, ⌊n / 2⌋ (§5.4). At n = 10 they are 20, 5, 4, 20 and 5. The walk's
limit above, 5n, follows n too, but has no parameter of its own.

Without a home GFF, a substitute annotation is built first. Genes are predicted
around the query's hits (`prodigal_on_regions.py`, Augustus or Prodigal), and
`borrow_annotations.py` transfers the genes of one annotated target genome.

Augustus's model (`augustus_species`, default `auto`) is chosen by
`resolve_effective_params.py`. With any `--modules` other than `legacy` it is the
installed model nearest the home species in the NCBI taxonomy on disk
(`gene_predictor.nearest_augustus_model`): the model whose species shares the deepest
taxon with the home species, so the species' own model first. *Candida albicans* gets
`candida_albicans`; a firefly gets `tribolium2012` through the suborder Polyphaga. The
choice is recorded in `effective_params.json` (`settings`, `source` `home_lineage`,
`derivation`). Without the taxonomy database the model is `fly`, with a warning. Under
`--modules legacy`, `auto` is `fly`. A model the user names is used as given.

With any `--modules` other than `legacy`, a home that the taxonomy database on disk
places among the eukaryotes stops the run before the search when `gene_predictor` is
`auto` and Augustus is not installed, because Prodigal cannot predict genes with
introns. Without the database, or for a species it does not know, the run goes on, and
`auto` predicts with Prodigal when Augustus is missing. A failed Augustus run stops the
prediction (it does not fall back to Prodigal), and a run that does not finish within
`gene_prediction_timeout` stops the task. Target genomes are not gene-predicted in the
default search (their models come from miniprot). Only the optional structural discovery
(`--enable_structural_discovery`, which runs only with `--enable_plm_search` or
`--enable_structural_search` on) runs Augustus on target windows, with the same model; the
PLM and structural searches on their own run no Augustus.

The borrowing (`borrow_annotations.py`):

1. the home locus is the best interval of the primary home locus (§5.4), with a window of
   `pred_flank_window` (50 kb) on each side, the flank doubled until the window holds
   2 × `n_flanking_genes` predicted genes or spans 500 kb. The doubling that first reaches
   500 kb is kept, so a window that had to grow ends below 1 Mb (about 800 kb at the
   default for a short locus); the first window, the locus ± `pred_flank_window`, is kept
   as it is when it already spans 500 kb (a `pred_flank_window` of 250 kb or more gives a
   window of the locus plus twice that);
2. the target is the closest annotated target (at least 10 CDS rows) by the run's
   phylogenetic distances (PHYLO_SORT runs before this step). A target is any genome file
   the search reads (`.fna`, `.fa` or `.fasta`, plain or gzipped) with its annotation beside
   it under the genome's file name with `.gff` for its last extension (`cow.fna` → `cow.gff`,
   `cow.fna.gz` → `cow.fna.gff`), as easy mode and `--target_gffs` place it. Several at that
   distance, or no distances, are told apart by the median identity of the hits of the
   window's predicted proteins;
3. the window's predicted proteins are searched against it (MMseqs2, E ≤ 1e-3), their
   hits clustered per scaffold with the 50 kb gap, and the cluster with the most home
   proteins, ± `borrow_target_flank` (100 kb), is the target region;
4. each annotated protein of that region is aligned to the home window with miniprot;
   a model with identity ≥ 30 % and query coverage ≥ 50 % whose best protein in the
   target region is the one it came from becomes a gene model with CDS rows
   (`ID=borrowed_<target gene>`, the target gene's `Name` and, when it has one, its
   `product`, `transferred_from=<target>`) and replaces a prediction it overlaps by more
   than half. The figures name it after its target gene (its symbol, else its product).

The predictions and the borrowed models form one merged GFF. The home proteome of such a
run is built from it (the longest coding model per gene), without the genes over the
query's own home model, and the report flags it as partial (`home_annotation`).
The home genome holds annotation only in stretches: the prediction windows around the
query's home hits and, when gene models were transferred, the stretch they cover, from
the locus to the outermost transferred model on each side (not the whole window of step 1,
which can hold the prediction windows of other home hits). Their header lines
(`##synvoy-annotated-window <scaffold> <start> <end> predicted|borrowed`) travel in the
merged GFF, and the flank walk stays inside the stretch that holds the locus (windows that
overlap or touch are one stretch): beyond it the next gene the file lists was predicted
around another hit, possibly megabases away, and the genes in between are not annotated.
A prediction window that overlaps or touches the locus's stretch is part of it, so the walk
goes on into it there, over DNA that is annotated without a gap. A side that ends at the stretch's edge with fewer than `n_flanking_genes`
anchors is logged with a warning (raise `pred_flank_window`, or give the home GFF).
Predictions are written as gene, mRNA and CDS rows, and a CDS that belongs to no feature
of the merged file still counts as a gene; AUGUSTUS proteins carry their region and
genome coordinates like Prodigal's. In a run without target genomes nothing is predicted
or borrowed beforehand. `extract_flanking_genes.py` then predicts genes itself within
`pred_flank_window` (50 kb) of the locus and keeps the longest `pred_keep_pct` of them (a
tenth), at least 2 × `n_flanking_genes` + 1 (21).

Each flanking gene is written as its full protein and, with `exon_level_search` (on),
as one sequence per exon. The FASTA header of a protein carries the gene's home
coordinates and the length of the home genome (`home=start-end genome=bp`), the header
of an exon its own coordinates (`coords=`). The search derives its distances from
them ([§7, Distances](#distances-in-gene-units)). Every file SynVoy writes is read by the
first word of each header, so `FBgn0000001|exon_2` stays that record and no other.

The home BED also gets a marker row for the GOI (`GOI_<scaffold>_<start>`) when no
annotated gene is left at the best interval, and one row named the same way for the
annotated gene at each other interval of the locus. They are used by the region scoring
and the figures and are not queries. The `goi=` field of the anchors' headers spans the
hits of every interval of the locus (the hit span of a snapped row, columns 8-9).

### 5.8 The initial query set

`PREPARE_INITIAL_DB` writes one FASTA per home locus: the flanking sequences in home
order, then the GOI protein and its exons (for a tandem array that miniprot found, its
copies `GOI_copy_<n>`). If the GOI has neither exon sequences nor copies, the protein is
cut into halves, thirds and quarters (pieces of at least 20 aa) instead, and above 100 aa
into overlapping 50-aa windows as well. These pieces are never a GOI query of a window,
the GOI entry of the ownership panel or a leaf of the tree.

The order of the flanking genes in this file defines the home rank
(`iterative_search_runner.py` → `build_home_rank`).

## 6. Target order, quality check and waves

**Order.** `phylo_sort.py` gives each target a distance to the home species: the
number of taxonomy nodes that the two lineages do not share. It uses a local taxonomy
database when one is installed, otherwise the NCBI lineage of the assembly accession
in the file name, otherwise a MinHash distance between the genomes (k = 21). Targets
are sorted closest first. ETE3 reads the folder named by the environment variable
`TAXDB` only when that folder holds both `nodes.dmp` and `taxa.sqlite`. Otherwise it
opens its default database (`~/.etetoolkit/taxa.sqlite`), so two machines can order the
waves differently; `intermediate/phylo_sort/sort_provenance.json` records the strategy,
the database that was opened (path, size, date, md5) and the order.

**Quality.** `assess_genome_quality.py` fails a genome with N50 below 5 kb or more
than 500,000 sequences. With `qc_fail_policy = auto`, a failed genome is removed when
the pipeline chose the targets (easy mode without `--target_species`), and kept when
the user chose them (pro mode, or `--target_species`). A kept genome is flagged in the
console, in `synvoy_report.json` (`warnings`) and in the figures (†).

**Waves** (`define_waves` → `assign_waves_by_rank`). Genomes are grouped by their
position in the sorted list, not by the distance value. With `n` genomes, the genome
at index `i` has `q = i / (n − 1)`:

| q | Wave size |
|---|---|
| below 0.10 | 1 |
| 0.10 to below 0.35 | 2 |
| 0.35 to below 0.70 | 3 |
| 0.70 and above | 5 |

Five genomes give waves of 1, 2 and 2. The genomes of a wave are searched with the same
query set, up to `iterative_search_cpus` of them at a time (2; 1 under the `low_mem`
profile that `run_synvoy.sh` uses when no profile is given). After each wave, the GOI
models that qualify ([§7.12](#712-what-is-added-to-the-query-set)) are added to the query
set of the following waves.

## 7. The search in one target genome

`iterative_search_runner.py` → `process_single_genome`.

### 7.1 Genome-wide search

MMseqs2 `easy-search` aligns the whole query set (flanking genes, GOI, and GOI models
from earlier waves) to the genome. Hits are kept at identity ≥ 10 %, at least 10
aligned query residues and E ≤ 0.01. MMseqs2 keeps at most `mmseqs_max_seqs` targets
per query (300, its own default; a target is a stretch of one reading frame, about one
exon).

Every length floor of the search counts **aligned query residues**, qend − qstart + 1:
the query is a protein for MMseqs2, tblastn and Smith-Waterman alike, so this is the one
length all three report in the same unit. The engines' own alignment-length column is not
used for a floor: MMseqs2's translated search, run without a backtrace, writes the target
span in nucleotides there, tblastn and Smith-Waterman write alignment columns (query
residues plus gaps).

**One anchor per place** (`resolve_anchor_hits`). Where hits of different flanking
genes overlap by at least half of the shorter one, only the hit with the highest bit
score stays (ties: lower E-value, then gene id). A target gene is the ortholog of at
most one home gene; two flanking genes that are relatives of each other would
otherwise count one target gene as two anchors. `--disable_anchor_resolution` keeps
every hit.

**The best places of each anchor** (`limit_anchor_places`, `anchor_max_places`; for
`auto`, the default, see the genome axis below). A flanking gene usually hits more than
its ortholog: relatives and shared domains give further places, and those are what puts
two anchors near each other by chance. With `--anchor_max_places N` a gene seeds blocks
only at its N best places (a place = its hits within its own gap limit; best = largest
summed bit score); the GOI's own home gene among the anchors (a GOI proxy) keeps every
place. The other hits still count inside a window. Where a genome has no conserved
neighbourhood the limit can remove every block with two anchors (such blocks arise there
by chance); the search then falls back to the blocks with one anchor, as it does whenever
no block passes `min_block_genes`.

**The genome axis** (`genome_axis_blocks`). The preset modules run once per query and see
only the home side; a value that should differ per target genome is set here, in the
search, with any `--modules` other than `legacy`, so by default (`lib/GenomeAxis.groovy`
passes `--genome_axis`). Under `--modules legacy` the axis is off and
`anchor_max_places auto` keeps every place. The axis holds one value,
`anchor_max_places`: with `auto`, a genome whose blocks exceed the cap (§7.4) is built
again once from the 5 best places of each flanking gene, and the cap then cuts the
rebuild as before if it still binds; where the cap does not bind, every place counts, as
without the axis. The GOI's own home gene among the anchors (a GOI proxy) keeps every
place. The rescue of GOI hits far from their anchors (§7.11) is chosen again against the
rebuilt blocks: a hit that a dropped block's GOI zone held becomes a seed if it passes the
rescue's rules, and a hit that a rebuilt block's zone holds is left to it. The places
counted for E_syn (§7.13) are taken from all hits before the blocks and do not change. A
number the user gives is used for every genome. Each genome's value is logged
(`genome axis: anchor_max_places N (source; ...)`), written
to the search's table `hits/<locus>_<genome>.target_policy.tsv` (read into
`synvoy_report.json` → `target_policy`) with its source (`default`, `genome_axis`,
`user_set`) and the block counts before and after a rebuild, and to the result GFF
(`##synvoy-genome-axis`). The modules may not set an axis value (`GENOME_AXIS_KEYS`, part
of `FORBIDDEN_FOR_MODULES`).

Only hits of flanking genes are used to build blocks. GOI hits are not, with one
exception ([§7.11](#711-goi-hits-far-from-their-anchors)). If no flanking gene has a
hit, all hits are used: every cluster of GOI hits then seeds a block of its own, with no
neighbour test, and the search log warns. When the cap (§7.4) binds there, the blocks
holding the strongest hits (bits) are kept. Their GOI models carry `DispersedRescue=1`
with `RescueRule=no_flanking_hits`, and the report lists the genome under
`summary.genomes_without_flanking_hits`.

### Distances in gene units

Distances that describe a neighbourhood are given in *gene pitches* and converted to
base pairs for this home locus and this target genome (`synvoy_scale.py`,
`genome_scale_args`).

```
pitch      = span of the home anchors / (number of anchors − 1)
expansion  = length in the target / length at home
             of the longest run of anchors that lie in home order on one target scaffold
distance   = pitches × pitch × expansion
```

A run (`anchor_run`) is a series of anchor hits along one scaffold whose home ranks all
rise or all fall. One step goes at most 3 ranks ahead and may be 0.2 to 8 times as long
as the same step at home. Each step is measured in the target between two hits and at
home between the two query records those hits belong to: two exons when the queries are
exons, otherwise two genes. Both lengths are sums over the steps of the run.

The quotient of two lengths is used, not a median of the step ratios, because a
neighbourhood rarely grows evenly. Most neighbouring genes keep their distance and a
few large insertions make up the difference; a median does not see them. For example,
the same 17 genes span 78 kb at the FT locus of *A. thaliana* and 387 kb in *A. lyrata*
(4.9; the run gives 5.0), while the median of the ratios of nearby anchor pairs is 1.5.
An insertion that is more than 8 times its home step is crossed by a step of two or
three ranks when that brings the ratio back inside the bounds, and ends the run when it
does not.

The measured value is used only where the gene order is conserved: the run must hold at
least 5 anchors. In a genome where the order is lost, hits of two or three anchors still
lie near each other by chance and would give any factor up to the bound; in such genomes
runs stay short, while a conserved neighbourhood gives long runs. Otherwise the expansion
is the length of the target genome divided by the length of the home genome (kept
between 0.2 and 8), and 1 when that is unknown.

| Distance | Pitches | bp parameter | Fixed value (no home coordinates, or `-profile legacy_distances`) |
|---|---|---|---|
| block gap (§7.2) | `block_gap_genes` 7, and at least 1.25 × the largest gap between neighbouring anchors at home × expansion | `cluster_distance` | 150 kb |
| window margin (§7.4, §7.5) | `window_margin_genes` 7 | `region_padding` | 150 kb (adapted per block to 50–200 kb when `region_padding` is set as a number) |
| GOI margin (§7.5) | `goi_margin_genes` 2 | `goi_margin` | 50 kb |
| bridge limit (§7.3) | `synteny_bridge_max_gap_genes` 20 (2 × `n_flanking_genes`) | `synteny_bridge_max_gap` | 6 Mb |
| envelope margin (§7.11) | `dispersed_envelope_margin_genes` 9 | `dispersed_envelope_margin` | 200 kb |

A bp parameter set to a number replaces the pitches for every genome. For
`region_padding` the number is only a starting value: each block's margin is adapted to
the spacing of the hits on its scaffold (twice the average gap) and kept between
`padding_min` (50 kb) and `padding_max` (200 kb). `synteny_bridge_max_gap 0` switches
bridging off. When fewer than two anchors of the home locus carry home coordinates, every
distance not set in bp takes the value of the last column, and the window margin is
150 kb for every block (it is adapted to 50–200 kb only when `region_padding` is set as a
number, as `-profile legacy_distances` does). The values used are logged per genome and
written into the result GFF as a `##synvoy-distances` line, which the region scorer (§8)
and the hull rescue (§9) read.

### 7.2 Blocks

`identify_synteny_blocks`:

1. **Per gene.** Hits of the same flanking gene (its exons count as the gene) form one
   *gene locus* when the gap between them is below 1.5 × the gene's span at home ×
   expansion, and at least `max_intron` (20 kb). A larger gap starts another locus,
   which is a second copy or a paralog.
2. **Along the scaffold.** Gene loci are sorted by position. A locus joins the current
   block when it starts less than the block gap ([Distances](#distances-in-gene-units)) after the previous one ends;
   otherwise the gap is tested for bridging (§7.3), and if that fails a new block
   begins.
3. **Per block.** The collinear run is computed from one position per flanking gene
   (its first), ordered along the scaffold. Blocks are sorted by collinear run, then
   by the number of distinct flanking genes.

### 7.3 Bridging a gap

A rearrangement or an insertion can put a long stretch between two parts of one
neighbourhood. A block continues across a gap when all of these hold (`_can_bridge`):

- same scaffold, and the gap is below the bridge limit ([Distances](#distances-in-gene-units): 2 ×
  `n_flanking_genes` = 20 pitches, the length of a home neighbourhood of n genes per side);
- the block has bridged fewer than `synteny_bridge_max_per_block` (2; 0: no limit) gaps so far;
- the block already holds at least `synteny_bridge_min_anchors` (2) distinct flanking genes;
- the first gene after the gap continues the block's home order: its home rank lies
  beyond the block's outermost rank, in the block's own direction, by at most
  `synteny_bridge_max_rank_gap` (5).

A bridge limit of 0, in genes or in bp, switches bridging off.

The direction is taken from the block's collinear run, so an inverted neighbourhood
bridges like a forward one. A neighbourhood whose direction reverses at the gap does
not.

| Flanking genes on both sides of the gap | Result |
|---|---|
| `A B C` … `D E F` | one block |
| `F E D` … `C B A` | one block |
| `D E F` … `C B A` | two blocks |
| `A B` … `D E F` | one block (two blocks with `synteny_bridge_min_anchors 3`) |
| `A` … `B C D` | two blocks (1 anchor) |

With two anchors the direction of the block rests on those two genes alone. Under
`synteny_bridge_two_sided` the same number is also the shortest collinear run that the
genes on both sides of the gap must form.

`synteny_bridge_two_sided` (off) weighs the whole cluster behind the gap instead of
its first gene. `--disable_synteny_collinearity` turns bridging and the collinear
ranking off. Without home order no block has a collinear run, so no block is
significant (E_syn ≥ 1): every GOI call is judged in the genome place, with no slot, no
gene-order-only call and no positional tie-break. The blocks are then not sorted before
the cap: those that pass the gene-count filter come in position order, and the GOI-seeded
blocks kept below `min_block_genes` come after them; the cap still keeps the GOI-seeded
blocks first (§7.4).

### 7.4 Merging, filtering, capping

1. **Merge** (`merge_synteny_blocks`). Two blocks on one scaffold are merged when
   their spans, each widened by the window margin ([Distances](#distances-in-gene-units)), overlap. Blocks less than
   two margins apart therefore become one. The merged block holds the union of the parts'
   flanking genes, so a gene with two copies counts once.
2. **Filter.** Blocks with fewer than `min_block_genes` (2) flanking genes are
   dropped. A block is kept regardless when it was seeded by the GOI's own home gene,
   by a rescued GOI hit (§7.11) or by a GOI hit (blocks are seeded from GOI hits when no
   flanking gene has a hit). If no block is left, all blocks are used. When the rescue's
   seeds are added (§7.11), the blocks are built again, and every block of the build
   without them is kept as well, so a fallback to all blocks survives a seed elsewhere.
3. **Cap.** At most `max_blocks_per_genome` (80; 0: no cap) blocks are searched. Blocks
   are sorted GOI-seeded first, then by collinear run, then by gene count; without home
   order they are not sorted (those that pass the filter in position order, the GOI-seeded
   blocks kept below `min_block_genes` after them). When the cap binds, it keeps first the
   GOI-seeded blocks and the blocks whose GOI zone holds a strong GOI hit that the rescue
   left to them (§7.11), then the others in that order, and the kept blocks keep their
   order, so no block number moves where the cap does not bind. In a genome without any
   flanking hit it keeps the blocks holding the strongest GOI hits (§7.1). Each dropped
   block is logged and counted (`blocks_capped`). On the genome axis (§7.1) a genome where
   the cap binds is built again from the 5 best places of each flanking gene; the same
   steps (`build_genome_blocks`, the rescue) then filter, sort and cap the rebuild.

### 7.5 Windows

`process_region_block`, per block:

- **Neighbourhood** = block ± the window margin ([Distances](#distances-in-gene-units)).
- **Core window** = block ± the GOI margin. A hit that lies entirely outside
  it is discarded, so every GOI model starts from a hit within the GOI margin of the
  flanking genes (2 pitches; the `local_family` module: up to the farthest home family
  member next to the GOI plus one gene, never beyond the window margin).
- **Search window** = the sequence that is cut out and searched: block ± the smaller of
  the window margin and (GOI margin + model margin). The model margin is the DNA a gene
  model may use beyond its hits: `gap_search_window`, or the home gene's span × expansion
  if that is more. DNA further out could only give hits that are discarded. Where genes
  lie close together the search window is the whole neighbourhood; where they lie far
  apart it is much smaller than the neighbourhood. (With `region_padding` set to a number
  the whole margin is searched: twice the mean gap between consecutive hits on that
  scaffold, limited to `padding_min`–`padding_max`, or the number itself when the
  scaffold has fewer than two hits. The whole margin is also searched when the home locus
  has no coordinates.)

Two numbers describe the block, written on every model of it (`BlockFlankingSupport`,
`BlockCollinearSupport`):

- **Flanking support** `B`: the number of distinct flanking genes with a genome-wide
  hit in the neighbourhood.
- **Collinear support**: the collinear run of those genes.

The label reads neither directly: the collinear run enters only through the block's
E-value of gene order (E_syn, §7.13).

### 7.6 Search inside the window

Two queries are used per home GOI (`select_goi_queries`):

- the **home protein**, always; and
- the gene model of the **relative searched most recently**: of the models that
  earlier waves added, those of the genome added last, and among them the one that
  aligns best to the home protein. A model identical to the home protein adds nothing.

Both queries are searched, and the window's hits are thinned by place within each query:
of hits of one query that overlap by more than half, only the one with the higher bit score
is kept. The two queries are thinned apart, so where the relative's protein aligns
better, the home protein still keeps its own hit there. Everything reported about a gene
model refers to the home protein. A model that only the relative's query produced is
aligned to the home protein again (`rebase_relative_candidates`): its identity (over all
alignment columns, gaps included, as BLAST counts it), coverage and confidence come from
that alignment, its id carries `_r<n>_` instead of `_l<n>_`, and
`SearchQuery` / `SearchQueryIdentity` record the relative and the identity to it.
Where the home protein produced a model on the same DNA, the relative's model is
dropped.
`--disable_relative_query` searches with the home protein only.

| Engine | Query | Kept when |
|---|---|---|
| MMseqs2, one thread | the protein, its halves, thirds and quarters (pieces ≥ 15 aa) and, above 100 aa, overlapping 50-aa windows in steps of 25; a hit's query positions are read back in the whole protein | identity ≥ 15 %, ≥ 15 aligned query residues (for a piece shorter than 25 aa: 60 % of its length, at least 12), E ≤ 0.1 |
| tblastn (`-seg no`) | the protein | identity ≥ 15 %, ≥ 15 aligned query residues (for a query shorter than 25 aa: 60 % of its length, at least 12), E ≤ 0.1 |
| Smith-Waterman (parasail), six frames | the protein | score ≥ 20, identity ≥ 10 %, ≥ 10 aligned query residues, E ≤ 1 |

With the `short_query` module the window ORFs and the hits count from 10 residues
(`window_mmseqs_min_length`, `aug_relaxed_length_min` 10), and tblastn runs with word size
2 and without composition-based statistics.

Smith-Waterman uses BLOSUM62 with gap costs 11/1. Its E-value is the Karlin-Altschul
value for the six-frame length of the window. It returns the best alignment of each
reading frame. The traceback needs memory in proportion to query length × frame length,
so a frame beyond 200 million cells is aligned in overlapping pieces (overlap: twice the
query length) and the best alignment of each piece is returned; up to 1.1 Mb of window for a 550-aa query nothing changes. A query
longer than 10,000 aa is not searched with Smith-Waterman (a warning says so): the
memory of one piece grows with the square of the query length. Without
parasail the search stops at its start: the bit scores and the specificity terms of the label
([§7.8](#78-confidence)) need it as well.
Each piece is scored without a traceback first (parasail's score-only kernel); the
traceback runs only where the score can give a hit that is kept (E ≤ 1, read from the
printed value: one raw score below the first score with E ≤ 1), so the kept hits are the
same as with a traceback on every piece, and the work of a window, 2 × window bp × query
aa cells, takes a fraction of the time. The time decides nothing: a window whose
Smith-Waterman search has not finished after `sw_timeout_seconds` (3,600 s), or whose
search fails in any of the three engines, stops the genome's search and the task (when
every failure among the genomes searched at the same time is a kill for memory: exit 137,
which Nextflow retries with more memory).
With `sw_max_alignments_per_piece` above 1 (default 1; the `local_family` and
`exon_structure` modules set 3 when the home gene has two or more copies or two or more
exons), a piece is aligned again after each alignment, with the residues of that
alignment masked in the frame (`*`), and
further alignments are returned while their E-value is at most 1, up to that number per
piece: a second exon or a second copy of the gene in the same frame is otherwise hidden
behind the best alignment. An alignment that runs through a masked span ends the piece;
one that ends on the last residue of a piece is left to the next piece, which holds it
whole.

MMseqs2 runs with E ≤ 10 and tblastn with E ≤ 1000; the limit of 0.1 is applied when
their output is read. tblastn's cutoff is wide on purpose: through its composition-based
statistics (on by default) the cutoff also changes the E-values tblastn reports, in either
direction, so with a cutoff of 10 a weak true gene could vanish from its output. The hits
of the three engines are pooled in order of bit score, and a hit is dropped when more than
half of it overlaps a hit of the same query already kept. The genome-wide hits
that lie in the window are then added (GOI hits only for the two queries above).

### 7.7 From hits to gene models

For each GOI query:

1. **Loci.** Hits less than 40 kb apart form a locus (more when the home gene has
   long introns: 5 × its largest intron). The best `goi_candidate_regions_per_query`
   loci are kept (2; the `local_family` module raises it to the home copies + 1, at most
   40), ranked by query coverage, then best bit score, then number of hits. When loci are
   cut, one INFO line per window and query names them (position, strands, hits, coverage,
   bits; the first 10, strongest first) and the weakest locus kept; decoys are not logged.
2. **Gene models** (`annotate_models_from_hit_list` → `miniprot_gene_models`). miniprot
   aligns the query to the locus ± `gap_search_window` (50 kb; at least the home
   gene's span × expansion), with a maximum intron of that stretch (at least 20 kb)
   and sensitive settings (`-n 2 -p 0.3 -N 50 --outs=0.3`; below 60 aa also `-L 10`,
   below 100 aa with the `short_query` module). Queries below 30 aa are aligned without
   splicing. **Every model is kept** whose coding sequence overlaps the coding sequence
   of a better model by less than half (`exon_annotation`; ids `…_l1_…`, `…_l1m2_…`, …;
   at most `goi_max_models_per_locus`, 10; the `local_family` module: max(10, home
   copies + 1), at most 60). While a hit of the locus lies on no coding exon of a kept
   model, the coding exons found so far are masked and miniprot runs again. A model of a
   later round must reach 30 % of the first run's best score (miniprot's own `--outs`
   bound in that run), unless its coding sequence overlaps a hit of the locus that no kept
   model explains: the hit says the query aligns there, and the judge (§7.8) labels the
   model by its own E-value. When the cap stops the rounds with such a hit, or an
   alignment that would have been kept, still left, the report says so in `warnings`
   (`limits.goi_model_cap_stopped`). A model's protein is miniprot's own translation.
3. **Tandem copies** (`detect_tandem_duplications`). Hits within 5 kb of each other are
   one candidate copy. When two or more copies each cover at least 35 % of the query (or
   three or more hits mostly cover the same part of it), each such copy that lies on no
   model is written as a one-segment copy (`tandem_copy`). When miniprot returns no model,
   all of them are. Two copies of the gene within 5 kb of each other are therefore one
   candidate copy, and no tandem: a compact pair is found only as miniprot's models
   (step 2).
   With `tandem_query_restart` (off by default) copies are built the way the
   hit chain is: per strand, in coding order, a hit extends a copy when it continues the
   query (10 residues of overlap at most) within the GOI's intron bound. A copy counts at
   a union query coverage of ≥ 0.35 with a best hit of ≥ 20 aligned residues or ≥ 30 bits
   (both for a query under 150 aa); copies that share genome are one (a second program's
   hit of the same exon cannot fake a copy); it is a tandem when two counted copies overlap
   in the query by more than half of the shorter one (an inverted copy counts). Each copy
   is one mRNA row with one CDS per hit and the union coverage, and is judged on its hits
   with their stops. Without it (the default) two exons of one gene that each cover 35 % of
   the query and lie more than 5 kb apart count as two tandem copies.
4. **Start and stop** (`refine_model_termini`). miniprot reports the aligned part
   only. The first exon is extended upstream, in frame, to the nearest `ATG`, and the
   last exon downstream to the first stop codon. The walk covers at most the number
   of unaligned query residues plus 3 codons (10 codons when the query is aligned to
   its end; never more than 30 upstream, 20 downstream). It does not start at a
   splice site (`AG` before the first base, `GT` after the last), does not cross a
   stop codon, and changes nothing when it finds no start or stop.
   With `goi_first_exon_guard` (off by default) the walk is skipped where a
   short home first exon is missing: the home gene has ≥ 2 exons, its first is at most
   30 aa, and the query residues before the model's first CDS are within 5 aa of that
   length (for an exon of 8 aa or less: at least max(4, length − 2)). The nearest ATG then
   lies in the intron. The model keeps miniprot's start and carries
   `MissingQueryStart=1-<p−1>` and `FirstExon=missing`; the same at the stop end
   (`MissingQueryEnd`, `LastExon=missing`). Only for the home query and its decoys, not
   for the relative's query.
5. **Hit chain** (`fallback_hit_span`), when miniprot returns no model and no tandem copy
   is written. The hits on the strand with the larger summed bit score are chained in
   genomic order so that their query positions only rise (or only fall, on the minus
   strand). Two hits are not linked when they
   overlap by more than 30 bp on the genome, or when the gap between them exceeds
   5 × the home gene's largest intron (at least 10 kb). Lengths are aligned query
   residues, summed over the chain's hits. The chain is accepted when:
   - query below 150 aa: ≥ 15 aligned query residues and best bit score ≥ 30;
   - longer query, one hit: coverage ≥ 0.25 and ≥ 25 aligned query residues;
   - longer query, several hits: coverage ≥ 0.25, or ≥ 35 aligned query residues in all,
     or bit score ≥ 60;
   - and its genomic span is at most the largest of 100 kb, 150 bp per query residue and
     1.5 × the GOI's own span at home × the target's expansion (the rule of an anchor's
     gap limit).

   The first three tests apply only outside a block with E_syn ≤ 0.1 (§7.13). Inside such a
   block, which predicted the gene's place, a chain that fails them is written all the same,
   with `RefusalLifted=` naming the tests it failed (`short_aln`, `short_bits`, `long_qcov`,
   `long_aln`, `long_bits`), and judged by its E-value in its place like every candidate
   (§7.8): a small or hidden fragment of the gene where miniprot finds no model gets a model
   and the label it earns. Such a call neither seeds nor goes to the tree (§7.12). The span
   test applies everywhere.

   The identity of a chain is the mean over its hits, weighted by aligned query residues.
   The mRNA row lists the query interval of each hit (`ChainQuery=1-40,41-88,...`).
6. **Single hits.** A hit that no kept model covers (more than half) is written as
   `rescued_exon`. A locus with no model and no tandem copy can give one `raw_hit`:
   identity ≥ 90 % (`max(gap_min_identity, 90)`), at least max(`min_exon_alnlen` 30,
   0.6 × query length) aligned query residues and coverage ≥ 0.6; with miniprot
   installed only at a locus with ≥ 2 hits, union coverage ≥ 0.75 and a best hit of
   ≥ 60 bits. Both are always fragments (`ModelStatus=fragment`), get an E-value in their
   place and are judged like every GOI candidate (§7.8); a `raw_hit` can be HIGH when it
   covers at least 0.7 of the query.
   With `goi_join_leftover_hits` (off by default) a `rescued_exon` hit is first
   offered to a gene model of its own query that the window's judge calls HIGH or MEDIUM
   (the home query or its decoy, not the relative's query): the hit must lie outside the
   model on its strand and continue the query from the model's first or last CDS (it
   starts at most 10 residues before the other ends and ends after it, on the minus
   strand the mirror image; the residues it skips fit into the
   genomic gap, 3 bp each less 30 bp; a gap under 40 bp only between two pieces of one
   exon, in one reading frame; within the GOI's intron bound), with no CDS of another gene model or
   tandem copy in between. miniprot then runs again over the model and the hit (± 1 kb),
   and the new model, judged like every candidate, replaces the old one only when its CDS
   cover the hit and the judge's order (`rank_key`) puts it strictly first. It keeps the
   old model's id and lists the joined hits (`JoinedHits=`); the `rescued_exon` rows of
   the hits it covers are dropped. The join is tried after the window's judge, not inside
   the model loop, because the label needs the window's place terms.

Every candidate gets a bit score: the Smith-Waterman score against the home protein,
BLOSUM62 with gap costs 11/1, in bits (λ 0.267, K 0.041), written as `Bits=`, of the best
of the proteins the candidate is judged on (§7.13): its label protein in pieces between
stops (the model's translation through the frame guard, or a chain's, copy's or single
hit's hits translated again with their stops; two consecutive hits in one frame with no
stop between them lie in one open reading frame and are one piece, translated through the
gap) and its hit regions. The specificity uses that same protein.
The candidates of a window are judged first and then ranked in the judge's order
(`synvoy_label.rank_key`, §7.8: label, E-value, the block's E_syn, bits, the shorter span; a
chain over two copies scores what its better copy scores, and must not win over that copy);
one that overlaps an accepted candidate by more than half of its own length is dropped.
The strict fallback on single hits (`raw_hit`) ranks hits by bit score too.

### 7.8 Confidence

One function judges every GOI candidate (`bin/synvoy_label.py`, `judge`): in the choice
among overlapping candidates in a window, in the final label and in the hull rescue. The
report's dedup and the figures order calls by its sort key (`rank_key`). Two row types are
not judged (§14). The judge's inputs are the E-value of the home protein against the
model in the place SynVoy expected the gene (§7.13), the model's specificity against the
home proteome (§7.12) and its query coverage:

| Confidence | `Basis` | Rule |
|---|---|---|
| **HIGH** | `significant_specific` | Eeff ≤ 1e-6; margin ≥ 10 bits or no real competitor; coverage ≥ 0.7 |
| MEDIUM | `significant_specific` | Eeff ≤ 1e-3 and specific |
| AMBIGUOUS | `not_specific` | Eeff ≤ 1e-3, not specific: another home gene fits as well |
| AMBIGUOUS | `gene_order_only` | 1e-3 < Eeff ≤ 1, in the slot of a block with E_syn ≤ 0.1, and specific |
| LOW | `not_significant` | everything else |

- **Coverage** is `QueryCoverage`, the share of the home protein the model covers. With
  `core_completeness` (the `signal_anchor` module) a gene model also carries
  `CoreCoverage`, the share of the protein's core (without the signal peptide and the GPI
  tail found by `query_features.py`) covered by its exons' query intervals, and the judge
  reads the larger of the two; 0.7 stays the bar. So do the other completeness tests of a
  model (C6): the model status, the hull rescue's cost filter and the coverage its judge reads
  (the query span of its alignment; it writes `CoreCoverage` too), and the report's coverage
  demotion and `identity_coverage_decoupled` flag. A hit chain, a tandem copy or a rescued exon
  has no `CoreCoverage` and reads the whole protein, and so do the search's bars on single hits
  (`fallback_refusal`, the tandem detector) and the fused-model rule of the report's dedup,
  where a real signal-peptide exon behind a long intron has to count.
- **Coverage below 0.7** rules out HIGH, not MEDIUM. A model that covers less, such as a
  fragment (`ModelStatus=fragment`) or a hit chain written in a block with E_syn ≤ 0.1
  although the fallback gate refused it (`RefusalLifted`, §7.7), is MEDIUM when its E-value
  in its place is at most 1e-3 and it is specific, and is then counted unless the report's
  coverage rule relabels it (§10: identity ≥ 50 % at a coverage below 0.35). The report
  counts the counted fragments and lifted chains (`summary.counted_fragment_calls`,
  `summary.refusal_lifted_calls`).
- **Specific**: no real competitor, or a margin of at least 3.3 bits over the best other
  home protein. Within ±3.3 bits a model in the slot or core of a significant block is
  decided for the gene (positional tie-break), unless that protein lies in the gene's own
  home neighbourhood, or its position is unknown: in a tandem array position cannot tell
  the two apart.
- AMBIGUOUS calls ("needs curation") are kept, drawn and reported, never counted. A
  gene-order-only call is offered only in a genome where no home locus and no rescue pass
  produced a counted call of the gene; counted means after the ownership and coverage
  demotions of the report ([§10](#10-adjudication-and-the-report)). Where the report drops one, the per-genome
  GFF and the figures leave it out too (`apply_report_drops.py`; the figures are drawn
  after the report).
- **Slot.** Without the home locus span (`goi=` in the anchor headers) no block has a
  slot; the search then warns once and judges in the core, the window or the genome.
- **parasail** computes `Bits=` and the specificity terms, so the search stops at its
  start when parasail cannot be loaded, whatever the Smith-Waterman settings.
- **Precision.** The judge and the seed bar work on the full-precision values; the GFF
  holds E-values with 3 significant digits and `SpecMargin` at 0.1 bit, for reading only:
  a value that rounds to a bar does not pass it.
- **Ranking among candidates** (`rank_key`): label, then E-value (values below 1e-10
  count as equal), then the block's E_syn, then bits, then the shorter span, then
  position, then a tail without the GOI's name: the kind of model (gene model, tandem copy,
  hit chain, rescued exon, raw hit), strand, CDS rows and the id behind its `|`.
  In a window the label is provisional, on the E-value of the place before the
  genome's multiplier; the final label is set at the end of the genome on Eeff.
- **Frame guard**: a miniprot model's E-value is computed on the translation of its CDS
  (first exon's phase honoured). When that has a stop before its end or a length that is
  not a multiple of 3, another frame over the same exons is used if it is clean
  (`FrameFixed=1`), else only the model's hit regions are judged (`FrameBroken=1`).
- The hull rescue judges its models with the same function and the same inputs: the
  label protein through the frame guard (a broken frame is judged on the model's own
  translation in pieces between stops, since the rescue has no search hits), completeness
  as the query span of the alignment, specificity with the GOI's home neighbourhood (the
  anchors of the locus found in the genome). Its place follows the search's rule (§7.13)
  on the search's block table of the genome (`blocks.tsv`): the core or the window of a
  significant block that holds the model, the smallest of them; its E-value there is
  multiplied by `SigBlocks` (`Eeff`). Elsewhere, and without a block table, its place is the
  whole genome. The slot is not used (the table does not hold the placed neighbours), so a
  model in a slot is judged in its core. It writes no model that repeats a GOI model of the
  main search (same strand, each covering more than 80 % of the other) with an equal or
  better label: judged in a larger place than the main model, such a twin could read
  AMBIGUOUS where the main path's model reads MEDIUM.
- A model carries a provisional LOW (`InferenceReason=awaiting_label`) from the moment it
  is built; every GOI row with an E-value is relabelled at the end of its genome. A GOI
  row the judge cannot label (no Eeff) keeps the LOW, never counts and never seeds, and
  is counted in the log as `LOW/unjudged`.
- ML scores (`--enable_plm_search`, `--enable_structural_search`) are recorded on the
  rows (`EmbeddingSimilarity`, `StructuralSimilarity`); they lift a label by one step only
  on a row the judge did not label (`_ml_relabel`).
- **The genome-wide path** (§7.11, `goi_genome_wide_path`): a model of one of its windows is
  an answer only when the judge labels it HIGH in the genome place
  (`synvoy_label.genome_wide_answer`); it is then written as MEDIUM with
  `GOIClass=genome_wide_strong_hit` (a tandem copy too) and never counted (§10). A window of
  the path writes no other model.
- A model also gets a status that does not depend on confidence: `fragment` below
  0.4 coverage, `complete` from 0.7 coverage with at least 2 exons (a tandem copy: one exon
  is enough), else `partial`; a rescued exon or raw hit is always `fragment`. The
  `short_query` module sets 0.3 and 0.5.

### 7.9 No relaxation for distant targets

The judge of §7.8 has no identity threshold, and its E-value is computed in the expected
place of each genome, so nothing is relaxed for a genome whose flanking genes match
poorly. The median of the best identity of each flanking gene is logged per genome
(`Median flanking identity: … over N flanking gene(s)`); it changes nothing.

### 7.10 Flanking gene models

Flanking genes with hits in the core window are modelled the same way (miniprot, else a
hit chain) and labelled by `_classify_flanking_evidence`: a miniprot model is HIGH with
at least 2 exons or identity ≥ 55 %, otherwise MEDIUM; a hit chain is MEDIUM with
coverage ≥ 0.65 and identity ≥ 55 %, otherwise LOW. One model per flanking gene and
locus is kept across blocks.

A flanking gene that got no model in any block is looked for on the scaffolds that hold
no block, in a genome with at least one block. Only its strongest locus there is used, and
only when it has a hit of at least 40 bits; it is modelled and written as a rearranged
flanking gene, so the figures can show it. When that gene is the GOI's own home gene (a
GOI proxy), the locus gets a full GOI search instead.

### 7.11 GOI hits far from their anchors

`select_dispersed_goi_seeds` (on by default) takes a strong GOI hit that lies outside the
GOI zone of every block (block ± the GOI margin, the core of §7.5); it does not look at the
gene order. The blocks are built without it first (`blocks_with_rescue`), and their GOI zones
decide. A genome-wide GOI hit with identity ≥ 40 %, E ≤ 1e-10 and at least 50 aligned query
residues (qend − qstart + 1; the `short_query` and `exon_structure` modules: min(50, max(10,
round(0.6 × Lpiece))), Lpiece the shorter of the query and the home gene's longest exon):

- inside a GOI zone, is left to that block, whose own search looks for the gene there;
- outside every zone, becomes a block of its own when the best hits of at least 3 distinct
  flanking genes lie on its sequence (one per flanking gene: its hit with the highest bit
  score in the genome) and it lies between the outermost of them, widened by the envelope
  margin ([Distances](#distances-in-gene-units): 9 pitches) (rule `envelope`);
- outside every zone, on a sequence that carries fewer such best hits, becomes a block of its
  own when it lies within the GOI margin of one of them (rule `near_neighbour`): a single
  neighbour forms no block of its own (a block needs `min_block_genes` genes), so no other
  search looks beside it;
- otherwise is refused: outside the envelope of its sequence (`outside_flanking_envelope`), far
  from the best hits of a sequence that carries fewer than 3 (`chromosome_lacks_min_anchors`),
  or in a genome where no sequence carries 3 (`no_flanking_envelope_in_genome`).

Of the hits that pass, one per stretch of `region_padding` of a scaffold seeds the block: the
home protein's hit first, then those of its exons, pieces and earlier genomes' models, then other
GOI parents (seeded copies), decoys last; within one of these by the query id behind its `|` (a
copy's parent cut), then by position, strand and bits. With seeds, the blocks are built again from
the flanking hits and the seeds (§7.4), and a seed near a block merges into it as any seed hit
does. A seeded block is searched and judged like every other. Its GOI models that lie outside
every GOI zone of the blocks built without the rescue carry `DispersedRescue=1`, `RescueRule`
(`envelope` or `near_neighbour`), `RescueAnchors` (the best hits of flanking genes on the
sequence) and `RescueAnchorBp` (the distance to the nearest of them); the report counts the
counted calls among them (`summary.dispersed_rescue_goi`). A hit that a GOI zone held but whose
block the cap dropped is counted as `dead_zone` in the limit counts; the cap keeps such blocks
first (§7.4).

A refused GOI hit is listed in the report under `rejected_candidates`, with the rule it failed
and the nearest flanking gene, and merged with the other rows of its place, so the report can
say "found, not syntenic" (`summary.goi_found_but_not_syntenic`: genomes without a counted call
that hold a place no counted call of the run overlaps; LOW or AMBIGUOUS models there do not take a
genome out). A refused hit that the GOI zone of a block seeded by the rescue then holds (another
exon of a gene whose near hit seeded that block) was searched there and is no candidate. The hits
outside every searched GOI zone that fail the quality bar, and the dead-zone hits, are listed under
`rejected_candidates.dropped`, with the test each failed (the strongest 200 per genome search).
The searched GOI zones are the blocks'; a window of the genome-wide path (below) is not one of them.
Each flanking gene counts at its best
hit in the genome only. As a measure, each refused row also gives the flanking genes with any
hit on its sequence (`genes_on_sequence`), and `dispersed_refused_per_piece_pass` counts the
refused hits that counting each gene on every sequence where it has a hit would have seeded.

With `goi_genome_wide_path` (on by default; `false` switches it off), a GOI hit that lies outside every searched GOI
zone and outside every model the search wrote, where no block predicted the gene, gets an answer
of its own when it is strong: the strong hit without synteny. The step runs in every searched
genome that has the home proteome, with or without a block, after its windows:

- the GOI hits outside every searched GOI zone (block ± the GOI margin) and outside every
  written model are grouped per home GOI and sequence into
  places: hits within the larger of the model margin and the GOI's intron bound of each other,
  never across a zone. A hit between a GOI zone and the edge of its block's window (searched by
  the window, not kept by the zone) is such a hit; its place is near that block (`NearBlock`);
- a place is tried when the terms of the judge's HIGH hold for the home protein against its hits
  (in their frames, in pieces between stops), from one blastp run in the genome place
  (`-searchsp` = |q| × ⌊genome/3⌋): E ≤ 1e-6; the union of its HSPs with E ≤ 1 (at most one
  chance alignment that good in the genome) covers ≥ 0.7 of the query (with
  `core_completeness`, the larger of the query's and its core's coverage); and no other home
  protein is a real competitor within 10 bits (§7.12). The HSPs are not tested against the home
  gene's exon junctions: a local alignment of an ortholog with its introns can run 10 to 20
  residues past an exon end into the intron, so its query span says nothing about the target's
  structure;
- per home GOI the places are tried in the order E (equal below 1e-10), bits, coverage,
  position, until one gives an answer, at most `max_blocks_per_genome` windows per genome
  (`genome_wide_tries_cap`; the decoy pass below has a bound of its own): a window of the
  place ± the model margin, cut at the GOI zones and at the step's earlier windows (a place that
  an earlier window overlaps gets none; these windows are kept apart from the blocks' GOI zones, so
  the dead-zone count and `rejected_candidates.dropped` are as without the step), searched
  without the flanking pass, its
  candidates judged in the genome place (no multiplier);
- a candidate is an answer when the judge labels it HIGH there, no exon of its own gene model
  covers an exon junction of the home gene model with 10 residues on each side
  (`RetrocopyCandidate=1`: the target encodes both exons without an intron, as a retrocopy does),
  it is not frame-broken without an intron (`FrameBroken=1` with every gap between its CDS rows
  shorter than 40 bp: miniprot can split a frameshifted copy into two pieces a few bp apart), and
  it overlaps no model written before. The junctions are those of a home gene model of ≥ 2 exons
  (§5.6): the annotated gene, or the miniprot model when every exon-exon boundary has canonical
  splice sites (the donor of the exon before it GT or GC, the acceptor of the exon after it AG).
  The test reads the candidate's miniprot exons only, so a hit chain, a tandem copy and a
  relative's model rebased onto the home protein carry no `RetrocopyCandidate`. A candidate
  the test did not run on (no known junctions, or such a model) is no answer
  (`retrocopy_untested`), unless the home gene model has one exon: there is no junction to
  test. A window keeps its best answer and writes nothing else, so the step adds no LOW row
  and no call that needs curation;
- an exon across junctions whose CDS is at least 30 nt per junction longer than 3 × its aligned
  query span is taken to hold an intron at each, read as coding, and does not count as crossing
  them: miniprot aligns a short intron that is a multiple of 3 and free of stops as an insertion
  (on synthetic genes every 45- and 60-bp such intron, and some of 90 bp, gave one exon across
  the junction). The nt that the walk to the start and stop codon adds beyond the aligned span
  (§7.7, `refine_model_termini`) are not counted, and split codons and a stop codon add at most
  3 nt. Introns shorter than 30 nt are known from a few protists, ciliates such as *Paramecium*
  among them; there a read-through intron still reads as a retrocopy. The excess is counted,
  not located (the model's exons carry no alignment path): an exon across two junctions is
  excused by 60 nt of excess wherever in the exon they lie, so one read-through intron of 60 nt
  or more next to a junction whose intron was lost excuses it too, and the model gets
  `RetrocopyCandidate=0` although the target has no intron at that junction (an intron loss
  next to a kept intron, not a retrocopy of the whole gene); below 60 nt such an exon stays
  flagged. A target exon with an insertion of 10 or more residues per junction reads the same,
  so a retrocopy with one is not flagged.

An answer is written as MEDIUM (no synteny supports it), `GOIClass=genome_wide_strong_hit`,
`Basis=significant_specific`, with `GenomeWide=1`, `GenomeWideCoverage` (the place's union
coverage), `RetrocopyCandidate=0` (absent only with a one-exon home gene model), `NeighbourMedianIdentity` and, inside a
block's searched window, `NearBlock` (that block's number). `NeighbourMedianIdentity` is the
median identity of the flanking genes' best hits in the genome, recorded for the curator and never
a test: it has no stated false-call rate, flanking genes are the conserved anchors and so do not
measure how far the GOI has diverged, and it is missing in the genomes without flanking hits. An
answer never seeds and is not in the tree (§7.12, §11); the report keeps at most one per genome,
only where no call is counted (§10). Every place is listed in
`hits/<locus>_<genome>.genome_wide.tsv` and in the report's `genome_wide_path` block with its
terms and what became of it.

With `--decoy_mode` the genome place is calibrated too (§7.13): the reversed and shuffled home
proteins are searched in the genome in the same MMseqs2 run as the query set, their rows leave the
genome's hit table before anything reads it, and their places take the same steps against the GOI
zones as they were before the real pass (the decoy pass). Decoys calibrate chance and composition,
not a homologous non-ortholog: a paralog, a retrocopy that kept its frame or an assembly duplicate
can still be an answer, which is one reason the answer is never counted.

### 7.12 What is added to the query set

**Specificity** (`bin/synvoy_specificity.py`), computed for every GOI candidate with an
E-value of at most 1 in its place, when the home proteome and Smith-Waterman bits are
there. The model's protein x (the protein its `Bits=` come from) is aligned
(Smith-Waterman, BLOSUM62, gaps 11/1) to every protein of the home proteome (one per gene,
the longest isoform). The best one other than the GOI's own home gene is p*
(`BestOtherHome`); `SpecMargin` = bits(home protein, x) − bits(p*, x).
The GOI's own home gene is the gene the home step took as the query's (`gene_id` in
`goi_info.json`, found by name or by position, [§5.6](#56-the-home-gene-model)); a gene the
home step flagged as ambiguous (with a warning: it may be a relative of the query's own gene)
only when it is itself near-identical to the GOI's home protein. Without one, and always
for a decoy or a seeded copy, it is the first of the three home proteins that align best
to the GOI's home protein that is near-identical to it (≥ 95 % identity over ≥ 90 % of
it; the hull rescue compares the query). The home step's gene comes first because a home
annotation can differ from the query (an extra exon, a partial model).
A home protein in the GOI's own home locus (the own gene and its `n_flanking_genes` neighbours
on each side, and this locus's hull when it holds the own gene: the flank walk skips relatives, so
a hull spans an array that n genes do not) that passes the same identity rule against the GOI or its
own gene is no p* for a model either, unless the model favours it over the GOI and its own gene by 10
bits (the seed's bar). Without this, a home copy identical to the GOI would tie with it at
0 bits and make every copy in a target not specific. A target's copy of such a neighbour
favours it by more than 10 bits and stays not specific. A paralog below that identity
stays a competitor.
p* is a real competitor (`RealCompetitor=1`) when it aligns to x at E ≤ 0.01 over the
whole proteome, or when it is a home family member of the GOI (the two align at E ≤ 1e-3,
blastp with composition-based statistics and SEG, the GOI as the query: in a pairwise run
`-seg yes` masks the low-complexity stretches of the query, and the composition-based
statistics mask the subject's, so both sides are masked. The test is not symmetric, so the
GOI is always the query; a failed blastp call is logged and counts as no family).
`BestOtherInNeighbourhood=1` when p* lies between the outermost flanking genes of the home
locus (the span extended to the locus itself), or among the GOI's own home gene and its
`n_flanking_genes` neighbours on each side, also in a run from a secondary home locus:
there, position cannot tell the GOI from p*. The search is exhaustive, with no prefilter:
a prefilter can miss home paralogs. The proteome and its gene positions reach the search
inside the home MMseqs2 database folder (`home_proteome.faa`, `home_genes.tsv`). The
label (§7.8) and the seed rule below use these terms.

**Seeds** (`_classify_goi_for_seed_and_tree`), per query in a genome and home locus,
after all blocks are labelled: the HIGH calls, or the MEDIUM calls when there is no
HIGH call; of these, the ones with no real competitor or a margin of at least 10 bits.
They join the query set of the following waves, and the window search's relative query
is chosen among them. Every other HIGH or MEDIUM call is passed to the tree only;
AMBIGUOUS and LOW calls never seed. A hit chain written although the fallback gate refused
it (`RefusalLifted`, §7.7) neither seeds nor goes to the tree, and is left out before the
HIGH/MEDIUM pool is formed: a lifted HIGH call does not make the genome's MEDIUM calls
tree-only. So is an answer of the genome-wide path (`GenomeWide=1`, §7.11): later waves search
as they would without it, and a MEDIUM call of the main search in the same genome seeds by the
rule above. A seeded tandem copy is, by default
(`tandem_copy_parent seed`), recorded as `GOI_copy_<n>` and therefore a gene of its own in
the following waves, judged against the seed; a model found with it carries its identity,
coverage and bits against the home protein as well (`HomeIdentity`, `HomeQueryCoverage`,
`HomeBits`; no label reads them). With `--tandem_copy_parent home` it is a
model of the home gene and can be its relative query. Its record id ends with the query it
was found with (`GOI_copy_<n>|<genome>_b<k>_l<j>_by_<parent>`): `n` is counted per parent,
block and locus, so copies of two parents never share an id.

### 7.13 E-values in the expected place

Every kept GOI model carries the evidence the label uses (`bin/synvoy_evalue.py`):

- **`BlockEsyn`**: the expected number of windows as long as the block's neighbourhood
  that hold, by chance, a run of at least the block's collinear support
  (`BlockCollinearSupport`) of anchors in home order, either direction. Each home
  neighbour with m places in the genome lies in a random window of length W with
  probability 1 − exp(−mW/G). A place is one occurrence of the gene: its hits on one
  sequence, each at most the gene's own gap limit (1.5 × its home span × the expansion,
  at least `max_intron`) after the furthest end of the place so far, as for the anchor
  places of §7.1, so the exons of a gene with long introns are one place, not one each.
  The places are counted on every flanking hit of the genome-wide search, before one
  anchor per place. The number present is a
  Poisson binomial; for n present anchors in random order the chance of a monotone run
  ≥ k is computed exactly (Robinson–Schensted–Knuth). A block with E_syn ≤ 0.1 predicted
  the gene's place.
- **`Place`**: the smallest of slot (between the nearest found home neighbours on either
  side of the gene, widened by one gene pitch; the side of each neighbour comes from the
  home locus span, `goi=` in the anchor headers), core (block ± GOI margin) and window
  that holds the model, in a block with E_syn ≤ 0.1; else `genome`. A hull-rescue model
  (§9) takes the core or window of any significant block of its genome that holds it, never
  the slot.
- **`EPlace`**: blastp of the home protein against the model's protein and against each of
  its search-hit regions (cut at the aligned span, translated in the hit's frame, split
  at stop codons), `-comp_based_stats 2 -seg no -word_size 2 -threshold 11`, with the
  search space set explicitly to |q| × ⌊place/3⌋ (`-searchsp`), read at full precision
  (JSON). **`EPlaceTblastn`**: the best tblastn row of the home protein that reaches into
  the model (at least 15 bp of it in its own frame inside the model, the rule of the hit
  regions above), as tblastn reported it for the window (before the keep
  floor of 0.1 and the dedup against the other engines' hits), its E rescaled from the
  window to the place. The window tblastn runs once into an archive that `blast_formatter`
  writes twice: as the tabular output the search reads (byte-identical to a direct run) and
  as JSON, from which this E is read at full precision (tabular output prints E with 3
  decimals from 0.0009 to 0.1, so 0.00123 read as 0.001). The judge does not use it.
- **`Eeff`**, **`EeffTblastn`**: the two multiplied by `SigBlocks`, the number of
  significant blocks searched in this genome (at least 1; each is one predicted place); a
  model in the genome-wide place is not multiplied.

Every block searched, with or without a call, is written with its E_syn at full precision,
its neighbourhood and the final labels of its GOI models to
`blocks/locus_<N>/locus_<N>_<genome>.blocks.tsv` (it changes no call). With
`--plot_block_selection true` the figures choose from it (§8).

**Copies in the slot** (`--report_copies_in_slot`, off; the `local_family` module turns it
on). After the genome's final labels, every GOI model of a block that has a slot gets
`CopiesInSlot=`: the GOI models with `Place=slot` and a label of AMBIGUOUS or better in that
slot. Each of them gets `SlotRank=`, counted from the copy nearest the home neighbour on
the left (the slot's orientation comes from its two neighbours). The report pairs the home
copies (`home_copies` of `run_settings`) with the copies in the slot of each genome's best
call (`copies_in_slot`, e.g. `2:3`). Both decide nothing: equal rank is not positional
orthology (Dewey 2011), since each gain or loss in an array shifts every rank behind it.

**Decoy mode** (`--decoy_mode`, off): in every window the reversed home protein and two
shuffles of it (seeded by the decoy kind and the protein sequence, so a rerun repeats them and a
renamed gene gets the same decoys) are searched like the home protein, modelled,
labelled and judged with the same terms, and written to
`decoys/locus_<N>/locus_<N>_<genome>.decoy.gff` and `.decoy.tsv`. Each decoy case (home
gene × decoy kind) is chosen among its own candidates, so a call of one kind never
suppresses a call of another. A window that also searches the relative's query
(`<gene>~rel`) searches its reversed and shuffled copies too (`~rev~rel`, `~shuf1~rel`,
`~shuf2~rel`): like a real model of the relative's query, their models are measured against
the home decoy of the same kind, filed under it (`SearchQuery=` names the decoy query) and
chosen in its case, so the rates cover both paths of the search. Each decoy row also
carries `PlaceBp=` and its E-value in every nested place that holds it (`EPlaceSlot`,
`EPlaceCore`, `EPlaceWindow`, `EPlaceGenome`, before `SigBlocks`), so a rule can be
calibrated in any place. They never seed and never reach the
result GFF. Every decoy call is false, so their rate per label and place calibrates the
labels without a truth set. Each window is searched with three decoy queries besides
each real query, so the window search takes several times as long. The report sums the
decoy cases that were searched (home gene × target genome with a search window ×
decoy kind; the home genome itself is never searched) into `decoy_rate`: per label the
calls and the upper end of the exact 95 % interval. Unless `goi_genome_wide_path` is false, the decoys of
the home protein are also searched in the whole genome and take the genome-wide path's steps
(§7.11); their calls carry `GenomeWide=1`, are left out of these counts and are counted apart
(`decoy_gw_cases`, `decoy_gw_answers`: a case with a decoy answer and no window decoy call of
MEDIUM or better; `decoy_rate.genome_wide`). Decoy mode is a run-level option; no preset module
sets it.

## 8. Region scores

`cluster_grs.py` scores candidate regions per genome. They are used for the figures
and for the strong-synteny rescue ([§9](#9-rescue-passes)). They do not change the
confidence of a gene model. With `--plot_block_selection true` (off by default) a figure row
shows the search's blocks instead: the blocks of the HIGH, MEDIUM and AMBIGUOUS GOI calls
(an answer of the genome-wide path draws none), then the significant blocks without a counted
call by E_syn, up to 6 (fragment calls take part, hidden from view); a genome without its block
table keeps the regions.

**Clusters.** The genome-wide hits are put into one fixed order and clustered: a hit
joins the cluster when it starts less than the block gap of the search for this
genome after the furthest end of the cluster's hits so far (read from the result GFF;
150 kb for a GFF without the record). A cluster longer than twice the home neighbourhood
times the genome's expansion (`2 × home_span_bp × expansion` of the search's distance record,
the largest hull of §9) is split at its largest gaps until every piece fits; without the
record there is no such limit.

**GOI models and GOI hits.** A cluster is gene-bearing when it lies within 20 kb of a GOI
model of the search that is not LOW and no answer of the genome-wide path
([§7.11](#711-goi-hits-far-from-their-anchors)): an answer has no synteny and is never
counted, so the regions are those of a run without the path (the columns
`goi_model_confidence` and `goi_model_class` still name it where it lies in a region). Raw
GOI hits of the genome-wide search stand in for models only where the search wrote no GOI
model of any label other than an answer (the target GFF holds no other GOI row, or there is
none): beside a model, LOW included, they are not used, so a place whose only model is LOW
does not come back through its hits. "GOI model or GOI hit" below means this.

**Score.** With `N` flanking genes at home:

```
C = max(U, B) / N        coverage: U = distinct flanking genes in the cluster,
                          B = for a cluster within 20 kb of a GOI model (not LOW,
                          no answer) or a GOI hit: the HIGH flanking models of a
                          search block that holds such a GOI model and whose models,
                          ± 20 kb, the cluster touches (the block with most of
                          them); else 0
K = L / k                order: L = longest collinear run of the k ranked hits
S = max(a, d) / (a + d)  strand: a hits agree with the home strand, d disagree

quality       = 0.4·C + 0.3·K + 0.3·S
synteny_score = quality · C
score         = synteny_score + 0.15     if the cluster lies within 20 kb of a GOI
                                         model (not LOW, no answer) or of a GOI hit
                                         of the genome-wide search
                                         (synteny_goi_overlap_bonus, at most 0.15)
```

A cluster that takes `B` from a block also takes over the block's span. Identity is not
part of the score: the score describes the neighbourhood, the confidence of §7.8
describes the gene.

**p-value** (`estimate_pvalue`). The gene labels of the cluster's hits are redrawn 200
times, with replacement, from the labels of all hits in the genome (fixed seed 42),
and `synteny_score` is recomputed each time. `p = (count(null ≥ observed) + 1) / 201`.

**Order and selection.** Clusters are sorted by `B` (descending), then score, then p.
Clusters near a GOI model or GOI hit (as for the bonus) are placed first. When no cluster
lies near a GOI model or hit, each GOI model (not LOW, no answer, ± 20 kb) that no cluster
touches is put first as a region of its own (`goi_anchor`), with a score of 0.60 for a HIGH
model and 0.40 otherwise. With `max_regions = 0`, a region whose GOI model the search counted (the
best GOI model within 20 kb, as in `goi_model_confidence`, is HIGH or MEDIUM and no answer
of the genome-wide path) and every `goi_anchor` region are written whatever the floors and
the cap. Then every cluster with a score at or above the floor
(`max(0.03, 0.30 × the best score of all clusters)`) or with max(U, B) at least 3
(`adaptive_unique_gene_floor`) is added, in this order, until 6 regions are written
(`adaptive_max_regions`, counting the regions kept for their model; the `dispersed_family`
module: floors 0.45 and 0.08, at most 4); when nothing is written so, the first cluster is.
With `max_regions` N > 0 the first N clusters in this order are written.

**Strong synteny without a gene.** A cluster holding at least
`strong_synteny_min_flanking` (⌈`n_flanking_genes`/2⌉ = 5) HIGH flanking models, with no GOI
model (not LOW, no answer) within 20 kb of it and, where the search wrote no GOI model of any
label but an answer, no raw GOI hit there either, is written as `goi_missing_but_strong_synteny`. Beside any GOI model of the
search (a LOW one too) a raw GOI hit does not prevent this, so such a region can hold one
(`goi_hit_best_identity`). Up to two such regions are added beyond the cap,
unless they lie inside a selected GOI region.

**Dispersed hit.** A GOI-overlapping cluster with max(U, B) ≤ 1 that is not a
strong-synteny cluster gets a score of at least 0.45 when its GOI hits reach the dispersed
rescue's bar of §7.11 (identity ≥ `dispersed_goi_min_identity` 40 % and E ≤
`dispersed_goi_max_evalue` 1e-10, the best of each over its GOI hits), at least one flanking
gene lies on the same scaffold, and the genome's rearrangement score is at least 0.5. The
region step places no gene: such a region is `goi_dispersed_macrosynteny` only where the
search modelled the gene there (within 20 kb) as HIGH or MEDIUM, not as an answer of the
genome-wide path; elsewhere it is `goi_hit_no_call` (no model, or only an AMBIGUOUS or LOW
one). Both keep the score, so the regions chosen do not depend on the class. That score (`target_rearrangement_score`) is the mean of two fractions:
flanking genes off the main scaffold, and pairs of home-adjacent flanking genes that are
not within the search's bridge limit (§7.3; 2 Mb for a GFF without the record) of each
other in the target, each gene placed by its best hit (0 with fewer than two genes
placed).

Region confidence: HIGH from `min_synteny_score` (0.6), MEDIUM from half of it, else
LOW; a region with strong synteny and no GOI model or GOI hit near it
(`goi_missing_but_strong_synteny`) and a `goi_dispersed_macrosynteny` region are at least
MEDIUM; a `goi_hit_no_call` region takes the label of its score. A region's label rates the
neighbourhood, not a gene: the report lists every gene-bearing region of MEDIUM or better
(and every `goi_hit_no_call` region) without a counted call, with the reason
(`regions.goi_regions_without_call`).

## 9. Rescue passes

Both passes run miniprot with the normalised query on a window chosen from the
flanking genes alone. Both extend the start and stop as in §7.7, with the fixed limits
only (30 codons upstream, 20 downstream).

| | Strong-synteny rescue | Hull rescue |
|---|---|---|
| Script | `rescue_strong_synteny.py` | `rescue_goi_hull.py` |
| Runs when | a region is `goi_missing_but_strong_synteny` | a scaffold has ≥ `goi_hull_min_flanking` (round(0.4 × `n_flanking_genes`) = 4) HIGH flanking models in the hull and no HIGH GOI model lies in the DNA that would be searched |
| Window | the region ± 2 kb | the hull ± the GOI margin; skipped when the hull is longer than twice the home neighbourhood × expansion |
| Accepts | every miniprot model covering ≥ 5 % of the query | every model with identity ≥ 40 % and CDS ≥ `goi_hull_min_coverage` (0.5; the `short_query` module 0.3) of the query length, or, with `core_completeness`, an alignment over that share of the query's core (cost filters in front of the judge); the one the judge ranks first |
| Confidence | always LOW | the judge of the main search (§7.8): E-value in the place the search's rule gives the model (the core or window of a significant block of the genome that holds it, × `SigBlocks`; else the genome), specificity with the GOI's home neighbourhood, completeness = query span of the alignment |
| Evidence type | `relaxed_miniprot_rescue` | `synteny_hull_rescue` |

**The hull** (`compute_hull`). The HIGH flanking models of a scaffold are clustered
with the bridge limit of the search as gap limit. The centre of the largest cluster is
taken, and the hull is the span of all flanking models within one home neighbourhood
(× expansion) of it. It covers a gene that lies between two blocks. The three sizes
come from the distances the search recorded for the genome ([Distances](#distances-in-gene-units));
for a GFF without that record they are 3 Mb, 100 kb and 20 Mb.

The rescue models are published under `rescue/locus_<N>/` and reach the report.

## 10. Adjudication and the report

`generate_report.py` → `build_report` reads all loci and genomes once and applies, in
this order:

1. **Deduplication** (`dedupe_goi_records`, `dedupe_goi_annotations`). The answers of the
   genome-wide path (`GOIClass=genome_wide_strong_hit`, §7.11) are deduplicated among
   themselves, marked `genome_wide_strong_hit` (a field, since the class can become
   `cross_locus_duplicate`) and listed after the other records, which are exactly as without
   them; the steps below see both. HIGH, MEDIUM and AMBIGUOUS calls of
   one genome and scaffold are one gene when they overlap reciprocally by more than
   80 %, when one lies at least 90 % inside the other, or when two gene models on one
   strand meet only where the upstream model's last CDS row and the downstream model's
   first one overlap, partly and in one frame: the two halves of one reading frame,
   which miniprot can model apart where the query aligns
   twice inside it (`same_orf_merged`; a model whose terminal exon lies inside its neighbour's
   exon stays apart). The call the judge ranks first is kept (`_record_rank`: label, then
   `synvoy_label.rank_key`, i.e. E-value in its place, the block's E_syn, bits, span; bits
   or identity for a record without `Eeff` (a row the judge did not label, such as a
   `structural_discovery` row)), with the list of home loci that found it. The figures
   draw these records ([§12](#12-figures)). A HIGH or
   MEDIUM model contained in a wider one replaces it when the wider model is at least
   twice as long and adds less than 0.15 query coverage (both coverages recorded). Otherwise,
   when the wider model is a gene model or a rescued hit, AMBIGUOUS with `Basis=not_specific`,
   its verdict stands against a HIGH or MEDIUM fragment (`ModelStatus=fragment`) inside it, on its
   strand and coding bases and covering no more of the query: a piece of a gene that another home
   protein explains as well is no ortholog call (`overruled_counted_piece`). A hit chain or a tandem
   row can lie over several genes, and a hull rescue's verdict can come from its place alone, so
   none of them overrules on its own. When one of them is itself AMBIGUOUS with
   `Basis=not_specific` and the cluster holds a not-specific gene model or rescued hit of the
   same gene (each covering more than 0.8 of the other), the piece is tested against that
   twin. When the record later goes to another model, such a piece is clustered on its own
   again.
2. **Phylogenetic placement** (`phylo_placement_check.py`). Across the HIGH and MEDIUM
   models of a home locus (the tree's input), sequence divergence (1 − identity to the
   query) is fitted against the taxonomic distance of the species with a Theil-Sen line.
   A call whose residual, minus the median residual and divided by 1.4826 × the median
   absolute deviation, is 2 or more is flagged `phylo_discordant`: it has diverged more
   than its species should have, which is what a paralog looks like. The check needs at
   least `phylo_placement_min_calls` (4) calls with a species distance, at three or more
   distinct distances, and is advisory.
3. **Locus ownership** (`build_home_paralog_panel.py`, `build_locus_ownership`). A
   panel of home proteins is built: the annotated genes at each home locus (up to
   `locus_ownership_max_genes_per_locus`, 5; the `local_family` module max(5, `c` + 1)),
   plus up to `paralog_panel_max_family` (15; the `dispersed_family` module
   max(15, `Fam` + 1)) proteins of the home proteome that align to the query with
   identity ≥ 28 % and a Smith-Waterman score of at least a quarter of the query's
   self-score. The panel's GOI entry is the GOI's own home gene, found as the
   specificity check of §7.12 finds it: the gene the home step took (`goi_info.json`, by
   name or by position; a gene it flagged as ambiguous only when that gene is itself
   ≥ 95 % identical over ≥ 90 % of the GOI's home protein), else the first of the three
   home proteins that align best to the GOI's home protein that is ≥ 95 % identical over
   ≥ 90 % of it (the GOI's home protein is the GOI record of `goi_exons.faa`, §5.8). When
   no home gene is the GOI's own (a home genome without a GFF, whose proteome leaves out
   the gene over the GOI's home model; an unannotated ORF), the GOI's home protein itself
   is the GOI entry, at the primary home locus. A GOI gene that the gene cap leaves out of
   its locus is added there; one outside every locus span is entered at a locus of its own
   (`home_<gene>`), from which no search started, so a call it owns is marked as found from
   another home locus (`locus_reattributed`, and `goi_owner_search_fumble`). Each recovered
   protein is aligned to the
   panel. A panel gene that passes the own-gene rule against the GOI's home protein or its own gene
   (≥ 95 % identity over ≥ 90 % of it) counts as the GOI for a call that does not favour it over
   the GOI's own entry by 10 bits, as in the specificity (§7.12); the call then stays at a searched
   locus. A call is relabelled `paralog_not_goi` only when a paralog beats the GOI's best panel
   entry by at least 10 (`locus_ownership_tiebreak_gap`, SW score of BLOSUM62 at gaps 11/1, the
   specificity's scale); below that the call stays the GOI's (`ownership_near_tie`). If the
   two best scores differ by less than 10 and both or neither of them are the GOI, the locus
   with more HIGH flanking models decides.
4. **Coverage demotion** (`apply_coverage_demotion`). A HIGH or MEDIUM call with
   identity ≥ 50 % and a recorded coverage below 0.35 is relabelled
   `identity_coverage_decoupled`.
5. **Strong hits without synteny** (`drop_genome_wide_strong_hits`). An answer of the
   genome-wide path that ownership did not relabel is dropped in a genome that has a counted
   call (HIGH or MEDIUM, not relabelled in steps 3–4, no answer), where it overlaps a HIGH,
   MEDIUM or AMBIGUOUS record of any other source (another home locus, a rescue pass), and
   where another answer of its genome ranks first; the dropped ones are listed in
   `goi_dedup.genome_wide_strong_hits_dropped_records` with the reason, and the figures leave
   them out (`apply_report_drops.py`). One relabelled `paralog_not_goi` stays a visible record
   of that class and is no answer.
6. **Gene-order-only calls** (`drop_gene_order_only`). An AMBIGUOUS call with
   `Basis=gene_order_only` is dropped in a genome that has a counted call (HIGH or MEDIUM,
   not relabelled in steps 3–4, no answer of the genome-wide path); the dropped calls are
   listed in `goi_dedup.gene_order_only_dropped_records`. It runs after steps 3–4, so a
   genome whose only counted call is relabelled keeps its gene-order-only candidate.
7. **Headline.** HIGH and MEDIUM calls are counted without the two relabelled classes and
   without the answers of the genome-wide path. AMBIGUOUS calls are reported as needing
   curation, split by `Basis`; those that ownership gave to another home gene are counted
   apart. The kept answers are named apart, "strong hit(s) without synteny (not counted)",
   and those whose place lies inside a block's window "strong hit(s) near a block, outside its
   GOI zone (not counted)"; their genomes leave the lists of genomes without the gene
   (`strong_hit_without_synteny`, `strong_hit_near_block`). LOW calls are reported
   separately.
8. **Self-consistency flags**: strong synteny without a GOI (a genome, home locus and
   scaffold with at least 2 × `n_flanking_genes` HIGH flanking models and no HIGH GOI
   call; `strong_synteny_no_goi_min_flanking`), the same gene found from several home
   loci, and the flags of steps 2–4.

## 11. Tree

`compute_tree.py` aligns the GOI proteins of the search (the home GOI, the seeds and every
other HIGH or MEDIUM call, §7.12, but neither a lifted hit chain nor an answer of the
genome-wide path; no exon records or query pieces) with MAFFT (`--auto`,
the iterative refinement in one thread, `--threadit 0`, so the alignment does not depend
on the threads) and builds a tree with IQ-TREE (model selection `-m MFP`; 1000 ultrafast
bootstrap replicates from 4 sequences). With fewer than 3 sequences a placeholder tree
is written. With `skip_tree` (false; the `low_mem` profile sets true, and so does the
launcher's default `-profile auto,low_mem`) no tree is built and the placeholder tree is
written.

## 12. Figures

`plot_synteny.py` draws one set of figures per home locus from the home annotation,
the models of [§7](#7-the-search-in-one-target-genome) and the regions of
[§8](#8-region-scores). The figures are drawn after the report: `apply_report_drops.py`
copies the per-genome GFFs into `plot_inputs_*` without the calls the report dropped
([§10](#10-adjudication-and-the-report), steps 5 and 6; with the report's `figures` block,
exactly the models it marks `dropped`), writes on each GOI model of the home locus what the
report did with it (that block: `ReportOutcome`, `ReportRecord`, `ReportDrawnBy`,
`AlsoFoundFrom`, under the header line `##synvoy-report-marks`), and the figures are drawn
from those copies.

- **Labels.** A GOI call is drawn as resolved when its label is HIGH or MEDIUM, whatever
  its class; AMBIGUOUS calls are drawn in the needs-curation style. The best GOI
  call of a row, the overlap dedup and the GOI cap rank by label, then by the strength of
  the match (E in its place, the block's E_syn, bits). The GOI cap (`max_goi_per_genome`,
  10; the `local_family` module max(10, `c` + 1)) only steers which neighbourhood a row
  shows: every GOI model inside the chosen view is drawn. A counted call (HIGH or MEDIUM)
  that the view leaves out is listed in the anchor grid, placed by the rule of the other
  calls (at the locus when it lies within reach of the best call in view, else after `//`),
  and counted in the ribbon plot's label (`+N elsewhere`).
- **Records.** The GOI calls of a row are the report's records. A record is drawn once in
  each home locus whose search modelled it, by one model: its representative, else its
  model at the record's span, else the one that overlaps it most (the rule of the tree
  page; the rescue passes' models are no figure input), taken among the record's models
  that are no fragment where the locus has one, as the default figures hide fragments. The
  record's other models (`ReportOutcome=merged`, `ReportDrawnBy` the model that draws it)
  are not drawn; the anchor grid's source table lists them where the model that draws their
  record is drawn or listed in that figure. The models of a record the report dropped are
  not in the plot inputs. The overlap dedup does not remove a HIGH, MEDIUM or AMBIGUOUS call
  for another one, so two records that overlap are both drawn; it still removes a LOW call
  or a flanking model under a call. A call the report does not count (`paralog_not_goi`,
  `identity_coverage_decoupled`, an overruled piece) is drawn by its label, with the reason
  in its tooltip, the source table and the legend text; so is an answer of the genome-wide
  path that the ownership check gives to another home gene, which is no answer. An answer
  the report kept is drawn outlined (next point), with the report's outcome in its tooltip.
  A GFF without the header line is deduped by overlap alone (reciprocal overlap of at least
  0.5), as before.
- **A strong hit without synteny.** An answer of the genome-wide path that the report kept
  (`GOIClass=genome_wide_strong_hit`, MEDIUM, [§7.11](#711-goi-hits-far-from-their-anchors);
  `ReportOutcome=genome_wide_strong_hit` in a GFF with the report's marks) is no counted call
  in the figures: it is no copy, never a row's best call, in none of the counts of a row's
  label (`×N`, `+N weak`, `+N elsewhere`), and draws no block of the block choice. No model
  of the path places a locus: the regions ([§8](#8-region-scores)) and the stretch a row
  shows (below) are chosen without it, and a row's scaffold and position (the anchor grid's
  location column and the reach of a flanking copy, the gene-position map) are those of its
  other models. It is drawn outlined where it lies inside the row's stretch, with no GOI
  ribbon; the anchor grid lists it also outside that stretch (at the locus or after `//`),
  its source table gives it its own state, the legend text names the report's counts
  (`summary.genome_wide_strong_hits`, and `summary.near_block_strong_hits` for one inside a
  block's searched window), and the ribbon plot marks a row that holds one and no counted
  call ("Strong hit without synteny, not counted"). The gene-position map draws one GOI point
  per row: an answer only when it is the row's only GOI call in view and lies at the row's
  locus (the grid's test); otherwise the tooltip of the row's location names it, and one on
  another scaffold counts in `(+N)`.
- **GOI name.** The GOI column, the titles, the legends, the home row's label and the GOI's
  row of the Gene names table name the GOI after the annotated home gene on whose CDS the
  query hits fall. Without such a gene the name is the GOI's own, as its target models carry
  it, and `_gene_names.tsv` holds it as a row `GOI_<name>` (source `query_name`): a home gene
  that only overlaps the hits (an lncRNA, a gene with the GOI in an intron) is still drawn in
  the GOI's place on the home row, and the table names it there (`home row: …`), but it does
  not name the GOI.

- **Home row.** The GOI is drawn from the home annotation: the transcript whose CDS
  covers the hits best. Only coding exons are drawn.
- **What a target row shows** (without block tables). The models inside the genome's
  candidate regions ([§8](#8-region-scores)), after the GOI cap. When the regions that hold
  a GOI model hold at least half of these models (and at least 3), the row is narrowed to
  them and to the regions the scorer wrote as `goi_missing_but_strong_synteny`
  (`Reg<k>_FLANKonly<N>_…`): a neighbourhood without a GOI model, such as a homeolog that
  lost the gene, stays in the row. The plot log names every region left out that holds a
  model the kept regions do not also hold. When the regions hold no GOI model, the row
  shows the models within `flank_fallback_bp` of the GOI models on the scaffold that holds
  most of them. When the models lie on 20 sequences or fewer, a GOI scaffold without
  flanking genes is left out when another GOI scaffold has them, and the row keeps the GOI
  scaffolds, every scaffold with at least 3 flanking genes, and every scaffold that holds a
  flanking gene placed on none of those (most flanking genes first); when that would keep
  fewer than 30 % of the models on more than 3 scaffolds, it keeps them all. A model of the
  genome-wide path takes part in none of these steps, so a row shows what it shows in a run
  without the path; such a model is drawn where it lies inside that view (on a scaffold the
  view shows, between its first and last model), and the anchor grid lists the others.
- **Target rows.** A scaffold with at least 3 flanking genes is mirrored when at least
  70 % of them are on the opposite strand to home and their order does not contradict it,
  or when their order is clearly reversed and at least half of them are on the opposite
  strand. The counted calls on it that the view leaves out are mirrored with it, so the
  anchor grid draws every call of the scaffold in the same orientation.
- **Synteny plot.** Genes are stacked in at most 3 lanes; fragment models are hidden
  in the default figure and shown in a `_with_fragments` variant.
- **Anchor grid.** One row per species, one column per home gene. A flanking gene with
  more than one copy at the locus (on the scaffold of its best model, no farther from
  it than the neighbourhood is long) is drawn as a stacked arrow; the gene-positions
  figure (`*_gene_positions.html`) writes that count as `×N`. The GOI
  column has one small arrow per copy; a genome with more than ten copies shows the
  first ten models and `+N more` (`--grid_max_models`, a drawing setting that drops no
  copy, so no rule sets it). The number after a GOI model is its identity to the query;
  for a model found from a seeded tandem copy (`SynVoy_Parent=GOI_copy_<n>`,
  [§7.12](#712-what-is-added-to-the-query-set)) it is the identity to that copy, which the
  legend text and the source table (`identity_to`) say. Without the numbers, a flanking
  arrow of a HIGH ortholog is paler the lower its identity to the home gene; a MEDIUM or
  LOW one has a fixed pale fill.
- **Tree page.** The midpoint-rooted gene tree is split into `--clade_count` (4) clades from
  the root: the clade with the most leaves is split into its children, splits that leave no
  one-leaf clade first; branch lengths are not used.

## 13. Reproducibility

The same inputs give the same output files:

- the searches inside a window run on one thread: MMseqs2 because of
  `deterministic_goi_search`, tblastn and parasail always;
- the results of a wave are joined in genome order and sorted by id;
- sets are sorted before they are iterated, and `PYTHONHASHSEED` is fixed to 0;
- region hits are put into one total order before scoring;
- the permutation test uses a fixed seed;
- the gene tree: IQ-TREE has a fixed seed and runs in one thread (`compute_tree_cpus`, 1 by
  default), so the tree is the same every time. With more threads its branch lengths and
  bootstrap supports can differ slightly from run to run; the report then warns.

## 14. Properties and known limits

Consequences of the rules above that a reader of the results should know.

**Confidence**

- HIGH needs a query coverage of at least 0.7, whatever the number of exons.
- The hull rescue's models are judged like the main search's (§7.8). Its 40 % identity
  and 0.5 CDS-coverage gates (0.3 with the `short_query` module; with `core_completeness` the
  aligned share of the query's core counts as well) only decide which models
  are judged.
- Flanking support is a property of the window (every candidate in a window has the
  same value, also a copy at the edge of the window), and the label does not read it.
- Two GOI row types are not judged: `relaxed_miniprot_rescue` rows of the strong-synteny
  rescue (on by default; written LOW, with no `Eeff` and no `Basis`; LOW is never
  counted) and `structural_discovery` rows (opt-in, `--enable_structural_discovery`). A
  discovery model gets its label from its scores: HIGH from a TM-score of
  `structural_high_threshold` (0.7) or an embedding similarity of 0.85, MEDIUM from
  `structural_medium_threshold` (0.5) or 0.7, else LOW. It is kept only at a TM-score of
  at least `discovery_min_tm` (0.7), so at the defaults every discovery model is HIGH.
- A genome can hold several HIGH models for one home locus: the copies of a tandem
  array, or a relative inside the searched window. Each of them that the report keeps as
  its own record is counted in the headline's "ortholog(s)"; the pipeline does not say
  which one is the ortholog. Within a window, `l1` is the locus with the best coverage and
  score.
- Identity and coverage refer to the home protein, for every genome.

**Search**

- A stretch between two blocks that is neither bridged nor merged is not searched by
  the main path. The hull rescue covers it only when the scaffold has at least
  `goi_hull_min_flanking` (4) HIGH flanking models.
- Every block that passes the filter and the cap of §7.4 is searched; there is no
  early stop. Most windows hold no GOI model, and each can add LOW rows.
- A genome whose MMseqs2 search fails for lack of memory (a memory message of MMseqs2, or a
  step that dies without a message), after retries, is not searched: the report lists it,
  with the reason, under `summary.genomes_not_searched` and says in `warnings` that it is
  not GOI-absent. A search that fails with another error stops the run. When MMseqs2
  cannot fit the genome's database into its split limit (`Cannot fit databases into ...`)
  no retry is made, since every rung keeps that limit: the reason is
  `mmseqs_split_memory_limit_too_small`, and `warnings` names `--mmseqs_split_memory_limit`
  and `--iterative_search_memory`.
- MMseqs2 keeps at most `mmseqs_max_seqs` (300) targets per query; the report counts the
  queries that reached it (`limits`), a lower bound.
- Only the first sequence of a multi-sequence query is used (with a warning).
- There is no gap-filling search; missing exons are found by miniprot inside its
  window or not at all.
- miniprot does not report every copy of a gene in one run when it is allowed long
  introns (§7.7). The masked rounds find a copy only where the window search has a
  hit on it; a copy too diverged for any of the three engines stays unfound.
- Distances that describe a neighbourhood scale with the home locus and the target
  genome ([Distances](#distances-in-gene-units)). Three limits that describe one
  gene keep a floor in base pairs and follow the gene's own size above it:
  `max_intron` (20 kb), `gap_search_window` (50 kb) and the 40 kb that separate two
  GOI loci in a window (§7.7). In a compact genome these floors span several genes.
- The expansion is one number per target genome and home locus: the mean growth over
  the conserved run of anchors. A neighbourhood that grew in one place only (one large
  insertion next to the gene) has gaps several times that mean; the bridge rule (§7.3)
  is what covers such a gap. A run of only 5 or 6 anchors over a few tens of kb gives a
  rough number.
- The pitch is the mean step between anchors, so one very long gap at home (a gene
  desert next to the locus) raises every derived distance of that locus.
- The cost of a search grows with the pitch. A block needs two anchors within the
  block gap, and chance alone puts pairs of anchor hits that close; each such block is
  searched, and its window is as large as the pitch makes it. Where genes lie far apart
  (a pitch of several hundred kb, as at many vertebrate loci) blocks and windows are large
  and the search is slow; `-profile legacy_distances` gives fixed distances instead.
- The number of chance blocks grows with the number of places a flanking gene hits.
  Where flanking genes hit many places, the cap of 80 blocks per genome can be reached and
  the windows can cover much of the genome. `--anchor_max_places 5` (§7.1) lets each
  flanking gene seed blocks at its 5 best places only; on the genome axis (`--modules`
  other than `legacy`, the default) a genome where the cap binds gets that limit on its
  own.
- Where the gene order between home and target is lost (no five anchors in a row), the
  expansion is the ratio of the genome sizes, a coarse stand-in for the spacing of
  genes.

**Gene families**

The pipeline searches one gene at a time and reports the copies it meets: every gene
model of a window is kept (§7.7), and a home locus that holds several copies is recognised
(§5.3). A family is not treated as a unit. `--expand_goi_similar true` adds the relatives
next to the GOI as further queries (`GOI_NEIGHBOR_…`); it is off by default because in a
large family the extra queries fill the candidate regions.

**Scores**

- `synteny_score = 0.4·C² + 0.3·K·C + 0.3·S·C`. Coverage dominates; the weights are
  not shares of the score.
- The permutation test cannot produce the `max(U, B)` term, so p is too small when
  `B > U`. The p-value is reported and breaks ties in the region order; it does not
  select regions and does not enter any gene's confidence.
- Region selection adds every region that passes a floor, in the order of §8, until the
  cap. That order is not by score alone (regions near a GOI model or GOI hit first, then by
  the block support `B`), so when the cap binds, a region with a higher score can be left
  out behind them.

**Reporting**

- The figures and the tree are built from the search output. They do not include
  rescue models: a record made of rescue models alone is in no figure. The figures draw
  each record once, but a call that the report relabelled or does not count is drawn by
  its label (a kept answer of the genome-wide path outlined), with the reason in its
  tooltip (the tree page draws such leaves hollow).
- Locus ownership compares a call with proteins of the *home* genome. A paralog that
  exists only in the target lineage is not in that panel; the only check that can
  flag it is phylogenetic placement, which is advisory.

**Easy mode**

- Each taxonomy level requests the record of every assembly under the taxon. A level
  with several thousand assemblies can exceed the 600 s limit and is then skipped, so
  the run continues with fewer genomes.

## 15. Default values

The defaults are those of `nextflow.config`. `run_synvoy.sh` adds `-profile auto,low_mem`
when no profile is given; `low_mem` sets `mmseqs_sensitivity` 7.0 and `skip_tree` true.

| Parameter | Default | Used in |
|---|---|---|
| `search_evalue` | 0.01 | §5.2, §7.1 |
| `mmseqs_sensitivity` | 9.5 (`low_mem`: 7.0) | §5.2, §7.1 |
| `min_query_length` | 20 (warning below 40) | §5.1 |
| `max_loci` / `locus_min_bit_ratio` | 5 / 0.5 | §5.4 |
| `locus_gap_genes` / `locus_gap` | auto (⌊n_flanking_genes/2⌋ = 5) / auto | §5.3, §5.4 |
| `home_locus_mode` | auto | §5.5 |
| `n_flanking_genes` | 10 per side | §5.7 |
| `goi_relative_evalue` | 1e-5 (A1: 1e-3) | §5.7 |
| `rank_wave_binning` | true | §6 |
| `min_hit_identity` / `min_hit_length` | 10 / 10 | §7.1 |
| `max_intron` | 20,000 (floor of the per-gene limit) | §7.2 |
| `block_gap_genes` / `cluster_distance` | 7 pitches / auto | §7.2, §8 |
| `synteny_bridge_max_gap_genes` / `synteny_bridge_max_gap` | auto (2 × n_flanking_genes = 20 pitches) / auto | §7.3 |
| `synteny_bridge_min_anchors` / `_max_rank_gap` / `_max_per_block` | 2 / 5 / 2 | §7.3 |
| `window_margin_genes` / `region_padding` | 7 pitches / auto | §7.4, §7.5 |
| `goi_margin_genes` / `goi_margin` | 2 pitches (A1: min(W, max(2, R + 1))) / auto | §7.5 |
| `min_block_genes` / `max_blocks_per_genome` | 2 / 80 | §7.4 |
| `anchor_max_places` | `auto` (every place; on the genome axis 5 where the cap binds) | §7.1, §7.4 |
| `gap_search_window` | 50,000 (floor) | §7.7 |
| `goi_candidate_regions_per_query` | 2 (A1: max(2, c + 1), at most 40) | §7.7 |
| `goi_max_models_per_locus` | 10 (A1: max(10, c + 1), at most 60) | §7.7 |
| `sw_min_score` / `sw_min_identity` | 20 / 10 | §7.6 |
| `sw_max_alignments_per_piece` | 1 (A1, C1: 3 with ≥ 2 copies or ≥ 2 exons) | §7.6 |
| `goi_fallback_intron_margin` / `_floor` | 5 / 10,000 | §7.7 |
| `fallback_short_query_len` / `_min_aln_aa` / `_min_bits` | 150 / 15 (B1: 10) / 30 | §7.7 |
| `classify_fragment_max_qcov` / `classify_complete_min_qcov` | 0.4 / 0.7 (B1: 0.3 / 0.5) | §7.8 |
| `dispersed_goi_min_identity` / `_max_evalue` / `_min_alnlen` | 40 / 1e-10 / 50 (B1, C1: min(50, max(10, round(0.6 × Lpiece)))) | §7.11 |
| `dispersed_envelope_margin_genes` / `dispersed_envelope_margin` | 9 pitches / auto | §7.11 |
| `synteny_weight_base` / `_consistency` / `_strand` | 0.4 / 0.3 / 0.3 | §8 |
| `synteny_goi_overlap_bonus` | 0.15 | §8 |
| `adaptive_score_floor_frac` / `_abs` / `adaptive_max_regions` | 0.30 / 0.03 / 6 (A2: 0.45 / 0.08 / 4) | §8 |
| `strong_synteny_min_flanking` | auto (⌈n_flanking_genes/2⌉ = 5) | §8, §9 |
| `min_synteny_score` | 0.6 | §8 |
| `goi_hull_min_flanking` / `_min_identity` / `_min_coverage` | auto (round(0.4 × n_flanking_genes) = 4) / 40 / 0.5 (B1: 0.3) | §9 |
| `goi_hull_window_pad` / `_max_window` / `_cluster_max_gap` | auto / auto / auto | §9 |
| `paralog_panel_min_identity` / `_max_family` | 28 / 15 (A2: max(15, Fam + 1), at most 50) | §10 |
| `locus_ownership_max_genes_per_locus` | 5 (A1: max(5, c + 1), at most 60) | §10 |
| `locus_ownership_tiebreak_gap` | 10 | §10 |
| `strong_synteny_no_goi_min_flanking` | auto (2 × n_flanking_genes = 20) | §10 |
| `skip_tree` | false (`low_mem`: true) | §11 |
| `max_goi_per_genome` | 10 (A1: max(10, c + 1), at most 60) | §12 |
| `identity_decoupled_min_identity` / `_max_qcov` | 50 / 0.35 | §10 |
| `phylo_placement_min_calls` / `_z_threshold` | 4 / 2.0 | §10 |
