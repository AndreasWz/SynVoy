# SynVoy: the algorithm

This document describes what SynVoy computes, step by step, with the rule and the
default value at every decision. It is written from the code, not from earlier
descriptions: every statement was checked against the source of branch `dev` on
2026-10-07 (commit `f839bd8`). Code locations are given as `file` → `function`,
because line numbers drift.

It does not explain how to run the pipeline ([USAGE.md](USAGE.md)), what the output
files contain ([OUTPUT.md](OUTPUT.md)) or what every parameter does
([PARAMETERS.md](PARAMETERS.md)). Measurements behind several of the rules are in
[ALGORITHM_AUDIT_2026-09.md](ALGORITHM_AUDIT_2026-09.md).

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

Genomes are searched from the closest relative outwards. A gene found in a close
relative is added to the query set for the more distant ones.

## 2. Terms

| Term | Meaning |
|---|---|
| **GOI** | Gene of interest: the query. |
| **Home genome** | The annotated genome the query comes from. |
| **Target genome** | A genome in which the GOI is searched. |
| **Home locus** | One place in the home genome where the query aligns. A gene family gives several. Everything from [§5.7](#57-flanking-genes) on runs once per home locus. |
| **Flanking gene** | A protein-coding neighbour of the GOI in the home genome, used as an anchor. |
| **Home rank** | The position of a flanking gene in home genomic order (0, 1, 2 …). |
| **Gene pitch** | The mean step from one flanking gene of a home locus to the next: the span of the flanking genes divided by their number minus one. The unit of the search's distances. |
| **Expansion** | How much larger or smaller the home neighbourhood is in one target genome ([Distances](#distances-in-gene-units)). |
| **Hit** | One local alignment of a query protein to genomic DNA. |
| **Hit interval** | The query's hits at home that form one copy of its gene: chained while each continues the query, within the intron bound and one annotated gene ([§5.2](#52-locating-the-query)). |
| **Block** | A run of flanking-gene hits on one target scaffold: a candidate neighbourhood. |
| **Collinear run** | The longest subsequence of a block's flanking genes whose home ranks only rise or only fall along the target scaffold. |
| **Model** | A gene structure (exons, protein) built in a target genome. |
| **Query coverage (qcov)** | Fraction of the query's residues that are aligned, counted as the union of aligned intervals. |
| **Wave** | A group of target genomes searched in parallel with the same query set. |

## 3. Overview

| # | Step | Nextflow process | Script | Result |
|---|---|---|---|---|
| 1 | Resolve the input | `RESOLVE_GENE_INPUT`, `FETCH_HOME_GENOME`, `FETCH_RELATED_GENOMES` (easy mode); `STAGE_GENOMES` (pro mode) | `resolve_gene_input.py`, `fetch_home_genome.py`, `fetch_related_genomes.py` | query, home genome + GFF, target genomes |
| 2 | Normalise the query | `NORMALIZE_QUERY` | `normalize_query.py` | one protein sequence |
| 3 | Locate the query at home | `LOCATE_GENE` | MMseqs2, tblastn, `merge_hits.py` | hit intervals with bit scores |
| 4 | Choose a parameter preset | `RESOLVE_EFFECTIVE_PARAMS` | `profile_hits.py`, `auto_select_preset.py`, `resolve_effective_params.py` | the `settings` map |
| 5 | Model the home gene | `ANNOTATE_GOI` | `annotate_goi_exons.py` | GOI protein, exons, largest intron |
| 6 | Split into home loci | `SPLIT_LOCI` | `split_loci.py` | up to 5 loci, the annotated home gene first |
| 7 | Pick flanking genes | `EXTRACT_FLANKING`, `PREPARE_INITIAL_DB` | `extract_flanking_genes.py` | initial query set |
| 8 | Order and check targets | `PHYLO_SORT`, `ASSESS_GENOME_QUALITY`, `FILTER_SORTED_GENOMES` | `phylo_sort.py`, `assess_genome_quality.py` | targets, closest first |
| 9 | Search | `ITERATIVE_SEARCH` | `iterative_search_runner.py` | gene models per genome |
| 10 | Score regions | `CLUSTER_REGIONS` | `cluster_grs.py` | ranked regions per genome |
| 11 | Rescue | `RESCUE_STRONG_SYNTENY`, `RESCUE_GOI_HULL` | `rescue_strong_synteny.py`, `rescue_goi_hull.py` | extra GOI models |
| 12 | Adjudicate and report | `BUILD_HOME_PARALOG_PANEL`, `ASSIGN_LOCUS_OWNERSHIP`, `PHYLO_PLACEMENT_CHECK`, `GENERATE_REPORT` | `build_home_paralog_panel.py`, `reciprocal_best_paralog_check.py`, `phylo_placement_check.py`, `generate_report.py` | `synvoy_report.json` |
| 13 | Tree and figures | `COMPUTE_TREE`, `PLOT_SYNTENY` | `compute_tree.py`, `plot_synteny.py` | tree, HTML/SVG figures |

Nextflow runs a step as soon as its inputs exist, so steps without a dependency on
each other run at the same time. Steps 10, 11 and 13 all read the output of step 9;
only step 12 reads the output of step 11 ([§14](#14-properties-and-known-limits)).

## 4. Inputs

**Pro mode** takes local files: the query protein, the home genome and its GFF, and
the target genomes (optionally their GFFs).

**Easy mode** takes an accession or a sequence and fetches the rest.

1. `resolve_gene_input.py` → `detect_input_type` classifies the input as a UniProt
   accession, an NCBI protein accession (version optional), a FASTA file or a gene
   symbol, fetches the protein and records the species and the gene symbol.
2. `fetch_home_genome.py` downloads an assembly of that species with its annotation.
3. `fetch_related_genomes.py` chooses the targets by walking up the taxonomy:
   - **Levels.** The genus of the home species (when it resolves), then its family,
     order and class.
   - **How many.** With `max_genomes = 0`, three per level, at most 20. Each level
     gets `max(2, max_genomes // levels)`; places a level cannot fill are carried to
     the next level.
   - **Candidates.** All assemblies under the level's taxon. Each species is used at
     most once, and the home species is excluded (compared without strain names).
   - **Ranking** (`assembly_rank_tuple`, mode `hybrid`): reference genome before
     representative before other RefSeq before GenBank; then chromosome-level before
     complete before scaffold before contig; then fewer contigs, fewer scaffolds,
     larger contig N50, larger scaffold N50.
   - **Quality gate** (`is_bad_quality`): chromosome-level and complete assemblies
     always pass. Any other assembly is low quality when its best N50 is below 5 kb,
     when a known contig or scaffold count exceeds 500,000, or when it is contig-level
     and has no quality figures at all. With `bad_quality_policy = auto` such an
     assembly is dropped from the automatic search, and kept with a warning when the
     user named the species (`--target_species`) or it is the home assembly.

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
characters and a trailing `*`. Other non-sequence characters are an error. A
nucleotide query (DNA or RNA) is translated in six frames and the longest open
reading frame is taken. A protein shorter than `min_query_length` (20) stops the run;
one of 20 to 39 residues runs with a warning, because an alignment that short is not
significant on its own in another genome and a call then rests on the neighbourhood.

### 5.2 Locating the query

`LOCATE_GENE` searches the home genome twice and pools the result:

- MMseqs2 `easy-search`, translated search, sensitivity 9.5, E ≤ 0.01;
- `tblastn` with default settings.

`merge_hits.py` chains the hits into **hit intervals**, one per copy of the gene. It
keeps hits with E ≤ `search_evalue` (0.01), one per stretch of DNA (of two hits that
overlap by more than half of the shorter one the higher bit score stays: both programs
report every exon), and reads them per scaffold and strand in coding order. A hit
continues a chain when

- it continues the query: it starts at most 20 residues before the chain's last hit
  ends, and ends after it. An exon's alignment often runs into the intron before it,
  whose translation matches the end of the previous exon's part of the query (up to 19
  residues between hits inside one annotated gene in the test data); a second copy of
  the gene repeats the query (43 residues or more there) and starts a chain of its own,
  however close it lies;
- it lies at most `home_chain_max_intron` after that hit. `auto` is the largest intron
  of the annotated genes the hits overlap × `goi_fallback_intron_margin` (5), at least
  `goi_fallback_intron_floor` (10 kb) and at most 1.1 Mb; the floor without an
  annotation (the home gene model's intron, which the search uses, is not known yet);
- the two are not in two different annotated genes, and no annotated gene lies
  entirely between them (genes nested in an intron of the gene do not stop it).

A hit that several chains could take goes to the one whose last query end is closest
to its query start. A hit interval covers its hits and the other program's hits of the
same exons and carries their lowest E-value and highest bit score. Until 2026-10-16
(batch HG-1A) hits were merged when the next began less than 1 kb after the previous
one, whatever their query coordinates: two copies 500 bp apart became one interval
(fly *Adh* and *Adhr*), and the exons of a gene with introns over 1 kb stayed apart
(honeybee DPP4: 15 intervals).

### 5.3 Parameter preset and preset modules

`profile_hits.py` summarises the hits with E ≤ `search_evalue`. Hits that cover the
same DNA are reduced to the best one (both programs report every exon). Hits with at
most `locus_gap_genes` (`auto` = ⌊`n_flanking_genes`/2⌋, 5) annotated genes between them form one locus, as in §5.4. A
locus is scored as **one copy of the gene**: the best total bit score of hits that
cover the query at most once (`one_copy_score`). A locus is *strong* when its score is
at least half of the top locus's; `bit_ratio` is the top score divided by the second;
`primary_copies` is the number of aligned query residues at the top locus divided by
the query length (about 1 for one gene with any number of exons, about 2 for a tandem
pair). Each locus also gets its **coverage**, the share of the query covered by the
union of its hits: `secondary_full_length` says that a strong locus other than the top
one covers at least half of the query (a full-length relative), `shared_domain_only`
that there are two or more strong loci and every secondary one covers less (relatives
that share a domain). The profile carries the query's own features too
(`query_features.py`): length, self bit score, the share of low-complexity residues
(12-residue windows of Shannon entropy ≤ 2.2 bits, an approximation of BLAST's SEG),
the share in repeated 4-mers, an N-terminal signal peptide (n-, h- and c-region with the
−1/−3 rule; a heuristic, not SignalP), a C-terminal GPI-anchor signal, and cysteines.

**The old chooser.** `auto_select_preset.py` → `select_preset` applies the first rule
that matches:

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
| `preset_short_peptide` | model status: fragment below 0.3 coverage, complete from 0.5 |
| `preset_paralog_discrimination` | at most 4 regions per genome with a higher score floor (figures) |
| `preset_tandem_family` | a neighbour is a relative of the GOI up to E = 1e-3 (§5.7) |
| `preset_single_copy` | nothing |

No preset moves the label of a GOI call (§7.8). Until 2026-10-15 the short-peptide preset
also lowered the identity thresholds of the earlier label (HIGH 40, MEDIUM 30, tandem
copies 40) and three presets switched on a family-name gate (`strict_goi_family`); both
went with the identity-based label.

**Preset modules** (`select_modules.py`, since 2026-10-15). A gene is not one kind: it can
be short, cysteine-rich, secreted and in a tandem array at once, and the first-match rule
gives it one preset. A module is one kind of gene; each tests its own trigger, they
stack, and none comes first. The triggers read the hit profile with the query features,
the home gene model (`goi_info.json`) and the home neighbourhood (`home_context.py`:
family members of the query among the ±`n_flanking_genes` genes next to it, tested with
the label's family test, blastp with composition-based statistics at E ≤ 1e-3; and the
distance to the scaffold end).

| Id | Module | Trigger | What it changes (today) |
|---|---|---|---|
| A1 | `local_family` | `primary_copies` ≥ 1.5, or a family member within the flank walk | `goi_relative_evalue` 1e-3 (= `preset_tandem_family`); GOI loci kept per query in a window (`goi_candidate_regions_per_query`) 2 → home copies + 1, at most 40, where home copies = the larger of `primary_copies` (rounded half up) and 1 + the family members within the flank walk; gene models per GOI locus (`goi_max_models_per_locus`) 10 → max(10, home copies + 1), at most 60; the GOI margin (`goi_margin_genes`) 2 → min(`window_margin_genes`, max(2, R + 1)); Smith-Waterman alignments per piece 1 → 3; `report_copies_in_slot` on (§7.13) |
| A2 | `dispersed_family` | ≥ 2 strong loci and `secondary_full_length` | region floors 0.45 / 0.08 / 4 regions (= `preset_paralog_discrimination`) |
| A3 | `shared_domain` | `shared_domain_only` | nothing yet (a flag; A2 does not fire) |
| B1 | `short_query` | query < 100 aa (the headroom, self bits over the bits an alignment needs for E ≤ 1e-3 in a 100-kb window, is logged) | model status 0.5 / 0.3 (= `preset_short_peptide`); finding in the windows: tblastn word size 2 and no composition statistics (the run with tblastn's defaults still gives `EeffTblastn`), MMseqs2 ORFs from 10 codons (30), miniprot `-L 10` below 100 aa (60); the length floors of window hits and of the short-query fallback 15 → 10 aa, the dispersed rescue's 50 → 0.6 × `Lpiece` (10-50); LOCATE_GENE: tblastn once more with `-seg no` when a query < 100 aa has no home hit (recorded as `locate_retry`) |
| B2 | `low_complexity` | low-complexity share ≥ 0.10 or repeated 4-mers ≥ 0.25 | stricter only: the relative test of the flank walk with blastp and composition-based statistics (`goi_relative_statistic blastp_cbs`); `LowComplexity=` on GOI models |
| B3 | `signal_anchor` | a signal peptide or a GPI tail | completeness for HIGH also on the core of the home protein (`core_completeness`: `CoreCoverage=` on gene models, the judge reads the larger of the two coverages; §7.8) |
| B4 | `cys_framework` | ≥ 6 Cys in a 90-residue window | nothing (advisory) |
| C1 | `exon_structure` | a home gene of ≥ 2 exons (tandem copies of the home model not counted) | Smith-Waterman alignments per piece 1 → 3 (a second exon in the same frame); the dispersed rescue's length floor 50 → 0.6 × `Lpiece` (10-50) |
| C2 | `home_context` | the home gene < 40 kb from the end of a scaffold of ≥ 1 Mb | nothing (advisory) |

The thresholds marked as triggers are today's values or assumed, never tuned. Rules that
hold for every module:

- **A module never moves the label.** It may widen the search, tighten the inputs of the
  specificity test, change how completeness is measured, or add report labels; the
  judge's significance levels (§7.8) are constants, not parameters.
- **One rule per parameter** (`RULES`, since 2026-10-20). Every parameter a module may
  change has exactly one rule: its owners (the modules that switch it on), one formula of
  the run's design quantities, a cap, and the evidence for it. The rule applies when any
  of its owners is applied, so two modules can never disagree about a key and their order
  never matters. The design quantities are computed once from the features
  (`design_profile`): `c` the home copies (as defined for A1 above), `R` the reach of the
  home array in genes (the largest |rank| of a family member in the flank walk, 0 without
  one), `Fam` the home loci of the family, `Lq` the query length, `Lx` the longest real
  exon of the home gene, `Lpiece` = min(`Lq`, `Lx`), `Lcore` the query without its signal
  peptide and GPI tail, `nx` the real exons of the home gene, `n` =
  `n_flanking_genes`, and `W` = `window_margin_genes`. Today's rules give the values in the
  table above: max(2, `c` + 1) for the GOI loci per query (cap 40), max(10, `c` + 1) for the models per GOI locus (cap 60),
  min(`W`, max(2, `R` + 1)) for the GOI margin, 3 Smith-Waterman alignments per piece when
  `c` ≥ 2 or `nx` ≥ 2, min(50, max(10, round(0.6 × `Lpiece`))) for the dispersed rescue's
  floor; every other rule is a constant. No rule may write a switch (a whole step on or off,
  such as `decoy_mode`,
  a run-level option, or `tandem_copy_parent`), `n_flanking_genes`, `expand_goi_similar`,
  or the module mechanism itself (`FORBIDDEN_FOR_MODULES`).
- **Constraints** (`check_constraints`, after the rules and the values the user set): the
  GOI margin lies within the window margin, where the DNA searched ends (C1); with A1, the
  GOI loci kept per query and the models per GOI locus reach `c` + 1 within their caps (a
  cap that binds is named; C3); a per-hit length floor a rule sets lies between 10 and
  max(10, 0.6 × `Lpiece`) (C7); a value a rule set stays within its cap (C11); the
  ModelStatus fragment bound lies below the complete bound (C12). A violation that the rules
  alone made is a bug and stops the run before
  the search; one that involves a value the user set is a warning, and the user's value is
  kept.
- **Precedence**, lowest to highest: shipped default < modules < preset (the chooser's
  under `legacy`, or `--preset_override` in any mode) < a value the user set at launch
  (a `-profile preset_X` included) < `--user_override`.

`--modules` selects the mode: `legacy` (**the default** until each module is measured)
applies the chooser's preset exactly as before and only records which triggers fire;
`auto` applies every module whose trigger fires; `none` applies nothing; a list
(`short_query,local_family`) applies exactly those modules; deltas
(`auto,-dispersed_family,+cys_framework`) edit `auto`. The old preset names are the
module sets {B1}, {A1}, {A2} and {}, with the same values (A1 and B1 add the keys above,
which `legacy` never applies). LOCATE_GENE runs before the modules are resolved, so its
retry follows the `--modules` value itself: `auto` (unless it removes `short_query`) or a
list that names it; never `legacy` or `none`. A run without home hits stops today, so the
retry can only rescue a run.

The chosen values travel to the later steps as a `settings` map of 28 parameters
(`presetDefaultKeys()` in `main.nf`). A value the user changed at launch is kept over a
preset or a module. `--auto_apply_preset false` leaves the defaults in place;
`--preset_override` picks a preset by hand. `effective_params.json` records the mode,
the modules applied and why, the triggers that fire, the source of every value and the
features, and the report copies these as `run_settings`. The same file also holds the
design quantities (`profile`), the derivation of every value a rule, the preset or the
user set (`derivation`; for a rule its formula and inputs, e.g.
`goi_candidate_regions_per_query 4 = max(2, c + 1), c = 3 (A1)`), and the constraint
checks (`constraints`).

### 5.4 Splitting into home loci

`split_loci.py` runs after the home gene model (§5.6) is built:

1. Sort the hit intervals by position. Two neighbouring intervals on a scaffold
   belong to one locus when at most `locus_gap_genes` (`auto` = ⌊`n_flanking_genes`/2⌋,
   5 at the default) annotated genes lie entirely between them. On a scaffold without annotated genes, and for a home genome
   without annotation, the limit is 50 kb. `--locus_gap <bp>` uses a distance
   everywhere.
2. **Primary locus.** When the home gene model came from the annotation
   (`method: gff_annotation`), every group of intervals that overlaps that gene is
   replaced by one interval, the gene's own span, and this is the primary locus
   (§5.5). Otherwise the locus with the highest score is the primary locus.
3. Each locus is scored as one copy of the gene (`one_copy_score` over the raw hits
   that lie in it, §5.3), not by its best single hit.
4. Keep a further locus only if its score is at least `locus_min_bit_ratio` (0.5) of
   the primary's.
5. Keep at most `max_loci` (5): the primary locus and the highest-scoring others.
6. **Rows of a kept locus** (batch HG-1B). Each hit interval that an annotated gene
   covers at least half of is replaced by that gene's coding span (the longest coding
   model, as the home proteome takes it); its hit span stays in columns 8 and 9.
   Intervals on one gene become one row. With several rows the row holding the best
   copy of the gene (`one_copy_score` of the hits in it) comes first and is named
   `best_interval`: the flanking walk of the locus starts there (§5.7). The primary
   locus made from the annotated home gene stays its one row.

Scoring and filtering (steps 3 to 5) use the hit intervals as chained, before step 6.
More than 10 loci before filtering produce a large-family warning.

### 5.5 The annotated home gene as primary locus

`home_locus_mode` decides which loci are searched:

| Mode | Loci |
|---|---|
| `auto` (default) | the annotated gene that the query is, plus further loci that pass steps 4 and 5 of §5.4 |
| `gene_only` | the annotated gene's locus alone |
| `alignment` | the loci of the alignments alone (§5.4 without step 2) |

Without a matching annotated gene, `auto` and `gene_only` fall back to the
alignments (`gene_only` with a warning). The gene is what the query is; the
alignments around it are that gene's exons and, often, hits in a relative next door.

### 5.6 The home gene model

`annotate_goi_exons.py` builds the GOI's own exon structure, trying in this order:

1. **Annotation by name** (`match_goi_in_gff` → `_search_gff_by_name`). The UniProt
   names of the query (gene names, synonyms, the protein name) and the query id, with
   case variants and the `gene-` prefix, are looked up in the home GFF. A name only
   points at a gene, so each match is then compared with the query: up to 50 matches
   are compared, best name match first, and the first gene whose protein agrees with
   the query is taken (Smith-Waterman, identity ≥ 60 % and a score of at least half the
   query's self-score). Of that gene, the transcript that agrees best with the query
   is used. The best other gene compared is recorded as the runner-up.
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
   Without parasail the step uses the rule from before 2026-10-15 and says so in the
   log: the best share of shared 3-residue words, scaled by the length ratio, above 30 %.
   That rule needed about 75 to 80 % identity at equal length, so a query from another
   species or a short isoform fell through to miniprot. The old windows also came from
   every hit, including tblastn's hits up to E 10.
3. **miniprot** on the hits ± `gap_search_window` (50 kb), when no annotated gene matches.
4. **Query protein only**, when there are no hits to work from.

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
as a `GOI_` row without a protein. Until 2026-10-16 (batch HG-1B) every interval had a
walk and a GOI marker of its own, so a tandem array of three copies gave three walks.

1. Take the annotated gene whose centre is closest to the centre of the best interval.
2. Walk outwards on each side. Collect genes until `n_flanking_genes` (10)
   protein-coding genes are found on that side. A side stops early only when it has
   visited more than 5 × `n_flanking_genes` genes (50) **and** `n_flanking_genes` + R
   protein-coding genes, where R is the number of protein-coding genes on that side of
   the scaffold that overlap a home hit of the query (relatives the walk will skip;
   genes at the locus's own intervals are not counted, neither in R nor in the walk).
   Until 2026-10-16 (batch HG-4) every visited gene counted against the 50, so a family
   larger than that, or a run of tRNA genes, could end a side short of anchors; the
   rule only makes walks longer. Each side is logged (`[flank-walk]`: anchors, genes
   visited, relatives skipped, what ended it), with a warning when the limit ends a side
   with fewer than `n_flanking_genes` anchors.
3. Skip a gene that is a relative of the GOI: its protein aligns to the GOI protein
   (Smith-Waterman, BLOSUM62, 11/1) with an E-value at or below `goi_relative_evalue`
   (1e-5), computed for the two sequences alone. Such a gene would match wherever the
   GOI's relatives are. `0` switches the test off for the neighbours; the gene at the
   GOI's own position is never taken.
4. Genes without CDS are kept for the figures but do not count towards the 10.

Five counts follow `n_flanking_genes` = n (constraint C9; each is `auto` by default, resolved
in `lib/FlankingCounts.groovy` and `lib/LocusGap.groovy`): the bridge limit, 2n pitches
(§7.3); the HIGH flanking models of a strong-synteny region, n / 2 rounded up (§8); those
in the hull of the hull rescue, round(0.4 n) (§9); the gene gap between home loci,
⌊n / 2⌋ (§5.4); and the walk's limit above, 5n. At n = 10 they are 20, 5, 4, 5 and 50,
the fixed values they had (until 2026-10-23; the locus gap until 2026-10-15). A number
sets a count directly.

Without a home GFF, a substitute annotation is built first. Genes are predicted
around the query's hits (`prodigal_on_regions.py`, Augustus or Prodigal), and
`borrow_annotations.py` transfers the genes of one annotated target genome (batch HG-5b):

1. the home locus is the best interval of the primary home locus (§5.4), with a window of
   `pred_flank_window` (50 kb) around it, widened until it holds 2 × `n_flanking_genes`
   predicted genes or reaches 500 kb;
2. the target is the closest annotated target (at least 10 CDS rows) by the run's
   phylogenetic distances (PHYLO_SORT, which now runs before this step); several at that
   distance, or no distances, are told apart by the median identity of the hits of the
   window's predicted proteins;
3. the window's predicted proteins are searched against it (MMseqs2, E ≤ 1e-3), their
   hits clustered per scaffold with the 50 kb gap, and the cluster with the most home
   proteins, ± `borrow_target_flank` (100 kb), is the target region;
4. each annotated protein of that region is aligned to the home window with miniprot;
   a model with identity ≥ 30 % and query coverage ≥ 50 % whose best protein in the
   target region is the one it came from becomes a gene model with CDS rows
   (`ID=borrowed_<target gene>`, `Note=transferred from <target>`) and replaces a
   prediction it overlaps by more than half.

The predictions and the borrowed models form one merged GFF. The home proteome of such a
run is built from it (the longest coding model per gene), without the genes over the
query's own home model, and the report flags it as partial (`home_annotation`). Until
2026-10-16 the borrowed rows had no CDS (the walk read them as non-coding: drawn, never
anchors) and came from the best-annotated target, placed by the first row of the hit
table. If no annotation can be built at all, `extract_flanking_genes.py` predicts genes
within 50 kb of the locus and keeps the longest tenth, at least 21. Predictions are written as gene, mRNA and CDS
rows, and a CDS that belongs to no feature of the merged file still counts as a gene
(until 2026-10-15 one borrowed gene row made every Prodigal prediction vanish); AUGUSTUS
proteins carry their region and genome coordinates like Prodigal's.

Each flanking gene is written as its full protein and, with `exon_level_search` (on),
as one sequence per exon. The FASTA header of a protein carries the gene's home
coordinates and the length of the home genome (`home=start-end genome=bp`), the header
of an exon its own coordinates (`coords=`). The search derives its distances from
them ([§7, Distances](#distances-in-gene-units)).

The home BED also gets a marker row for the GOI (`GOI_<scaffold>_<start>`) when no
annotated gene is left at the best interval, and one row named the same way for the
annotated gene at each other interval of the locus. They are used by the region scoring
and the figures and are not queries. The `goi=` field of the anchors' headers spans the
hits of every interval of the locus (the hit span of a snapped row, columns 8-9).

### 5.8 The initial query set

`PREPARE_INITIAL_DB` writes one FASTA per home locus: the flanking sequences in home
order, then the GOI protein and its exons. If the GOI has no exon sequences, the
protein is cut into halves, thirds and quarters (pieces of at least 20 aa) instead.

The order of the flanking genes in this file defines the home rank
(`iterative_search_runner.py` → `build_home_rank`).

## 6. Target order, quality check and waves

**Order.** `phylo_sort.py` gives each target a distance to the home species: the
number of taxonomy nodes that the two lineages do not share. It uses a local taxonomy
database when one is installed, otherwise the NCBI lineage of the assembly accession
in the file name, otherwise a MinHash distance between the genomes (k = 21). Targets
are sorted closest first. Without a usable `TAXDB` (a folder with `taxa.sqlite`), ETE3
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

Five genomes give waves of 1, 2 and 2. The genomes of a wave are searched in parallel
with the same query set. After each wave, the GOI models that qualify
([§7.12](#712-what-is-added-to-the-query-set)) are added to the query set of the
following waves.

## 7. The search in one target genome

`iterative_search_runner.py` → `process_single_genome`.

### 7.1 Genome-wide search

MMseqs2 `easy-search` aligns the whole query set (flanking genes, GOI, and GOI models
from earlier waves) to the genome. Hits are kept at identity ≥ 10 %, alignment length
≥ 10 and E ≤ 0.01.

**One anchor per place** (`resolve_anchor_hits`). Where hits of different flanking
genes overlap by at least half of the shorter one, only the hit with the highest bit
score stays (ties: lower E-value, then gene id). A target gene is the ortholog of at
most one home gene; two flanking genes that are relatives of each other would
otherwise count one target gene as two anchors. `--disable_anchor_resolution` keeps
every hit.

**Optional: the best places of each anchor** (`limit_anchor_places`, `anchor_max_places`,
off by default). A flanking gene usually hits more than its ortholog: relatives and
shared domains give further places, and those are what puts two anchors near each other
by chance. With `--anchor_max_places N` a gene seeds blocks only at its N best places (a
place = its hits within its own gap limit; best = largest summed bit score). The other
hits still count inside a window. Replayed on the saved hits of 31 genome searches
(yeasts, flies, bees, ants, plants), N = 5 left 121 of 512 windows, and the one call that
was no longer in a searched block was a relative of the gene, not an ortholog. In five
test cases run with N = 5 the ortholog calls were the same (yeast: 18 windows instead of
101; *Arabidopsis* FT in four plants: 9 instead of 238). Where a genome has no conserved
neighbourhood the limit removes every block with two anchors (they were chance blocks),
and the search falls back to the blocks with one anchor, as it does whenever no block
passes `min_block_genes`: 19 windows instead of 3 for fruit fly against mosquito. It has
not been run on the large families of the benchmark, which is why it is not the default.

Only hits of flanking genes are used to build blocks. GOI hits are not, with one
exception ([§7.11](#711-goi-hits-far-from-their-anchors)). If no flanking gene has a
hit, all hits are used.

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
few large insertions make up the difference; a median does not see them. *A. lyrata* at
the FT locus of *A. thaliana*: the median of the ratios of nearby anchor pairs, the
measure until 2026-10-11, was 1.5, while the same 17 genes span 387 kb there and 78 kb
at home (4.9; the run gives 5.0). An insertion that is more than 8 times its home step
is crossed by a step of two or three ranks when that brings the ratio back inside the
bounds, and ends the run when it does not.

The measured value is used only where the gene order is conserved: the run must hold at
least 5 anchors. In a genome where the order is lost, hits of two or three anchors still
lie near each other by chance and would give any factor up to the bound; such genomes
have runs of 3 at most (fly against mosquito, firefly against fly), conserved
neighbourhoods runs of 6 to 19. Otherwise the expansion is the length of the target
genome divided by the length of the home genome (kept between 0.2 and 8), and 1 when
that is unknown.

Checked against a quantity the estimate does not use, the span of the flanking gene
models of the best block divided by the span of the same genes at home, in 27
genome–locus pairs with a conserved neighbourhood (ants, a wasp, flies, bees, yeasts,
three plants): the run quotient is within 5 % in 21 pairs and within 20 % in 24, with a
largest deviation of 33 % (a run of 6 anchors over 30 kb); the median was within 5 % in
14, within 20 % in 20, and 69 % too small in the worst case.

| Distance | Pitches | bp parameter | Fixed value before 2026-10-09 |
|---|---|---|---|
| block gap (§7.2) | `block_gap_genes` 7, and at least 1.25 × the largest gap between neighbouring anchors at home | `cluster_distance` | 150 kb |
| window margin (§7.4, §7.5) | `window_margin_genes` 7 | `region_padding` | 150 kb, adapted per block to 50–200 kb |
| GOI margin (§7.5) | `goi_margin_genes` 2 | `goi_margin` | 50 kb |
| bridge limit (§7.3) | `synteny_bridge_max_gap_genes` 20 (2 × `n_flanking_genes`) | `synteny_bridge_max_gap` | 6 Mb |
| envelope margin (§7.11) | `dispersed_envelope_margin_genes` 9 | `dispersed_envelope_margin` | 200 kb |

A bp parameter set to a number is used as it is. A query set without home coordinates
(from an older run) falls back to the fixed values. The values used are logged per
genome and written into the result GFF as a `##synvoy-distances` line, which the
region scorer (§8) and the hull rescue (§9) read.

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
- the block has bridged fewer than `synteny_bridge_max_per_block` (2) gaps so far;
- the block already holds at least `synteny_bridge_min_anchors` (2; 3 until
  2026-10-14) distinct flanking genes;
- the first gene after the gap continues the block's home order: its home rank lies
  beyond the block's outermost rank, in the block's own direction, by at most
  `synteny_bridge_max_rank_gap` (5).

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

With two anchors the direction of the block rests on those two genes alone. The floor
of 2 also is the minimum collinear run of the two-sided rule below, which is off; a
run of 2 is a weak bar there and would need its own floor if that rule is switched on.

`synteny_bridge_two_sided` (off) weighs the whole cluster behind the gap instead of
its first gene. `--disable_synteny_collinearity` turns bridging and the collinear
ranking off.

### 7.4 Merging, filtering, capping

1. **Merge** (`merge_synteny_blocks`). Two blocks on one scaffold are merged when
   their spans, each widened by the window margin ([Distances](#distances-in-gene-units)), overlap. Blocks less than
   two margins apart therefore become one. The merged block holds the union of the parts'
   flanking genes, so a gene with two copies counts once.
2. **Filter.** Blocks with fewer than `min_block_genes` (2) flanking genes are
   dropped. A block is kept regardless when it was seeded by the GOI's own home gene
   or by a rescued GOI hit (§7.11). If no block is left, all blocks are used.
3. **Cap.** At most `max_blocks_per_genome` (80) blocks are searched: GOI-seeded
   blocks first, then by collinear run, then by gene count.

### 7.5 Windows

`process_region_block`, per block:

- **Neighbourhood** = block ± the window margin ([Distances](#distances-in-gene-units)).
- **Core window** = block ± the GOI margin. A hit that lies entirely outside
  it is discarded, so every GOI model starts from a hit within 2 pitches of the
  flanking genes.
- **Search window** = the sequence that is cut out and searched: block ± the smaller of
  the window margin and (GOI margin + model margin). The model margin is the DNA a gene
  model may use beyond its hits: `gap_search_window`, or the home gene's span × expansion
  if that is more. DNA further out could only give hits that are discarded. At a yeast
  or honeybee locus the search window is the neighbourhood; at a human locus it is about
  0.7 Mb per side where the neighbourhood reaches 2.3 Mb. (With `region_padding` set to a
  number the whole margin is searched, and the margin is twice the mean gap between
  consecutive hits on that scaffold, limited to `padding_min`–`padding_max`.)

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
is kept. The two queries are thinned apart (since 2026-10-14), so where the relative's
protein aligns better, the home protein still keeps its own hit there. Everything reported
about a gene model
refers to the home protein. A model that only the relative's query produced is aligned
to the home protein again (`rebase_relative_candidates`): its identity (over all
alignment columns, gaps included, as BLAST counts it; until 2026-10-14 the gap columns
were left out and the identity read several points higher), coverage and confidence
come from that alignment, its id carries `_r<n>_` instead of `_l<n>_`, and
`SearchQuery` / `SearchQueryIdentity` record the relative and the identity to it.
Where the home protein produced a model on the same DNA, the relative's model is
dropped.
`--disable_relative_query` searches with the home protein only.

Until 2026-10-08 one query was used: the longest protein among the home GOI and all
earlier models. A model a few residues longer than the home protein replaced it for
the rest of the run, and identity and coverage of later calls referred to that model.

| Engine | Query | Kept when |
|---|---|---|
| MMseqs2, one thread | the protein and its halves, thirds and quarters (pieces ≥ 15 aa) | identity ≥ 15 %, length ≥ 15, E ≤ 0.1 |
| tblastn (`-seg no`) | the protein | identity ≥ 15 %, length ≥ 15, E ≤ 0.1 |
| Smith-Waterman (parasail), six frames | the protein | score ≥ 20, identity ≥ 10 %, length ≥ 10, E ≤ 1 |

Smith-Waterman uses BLOSUM62 with gap costs 11/1. Its E-value is the Karlin-Altschul
value for the six-frame length of the window. It returns the best alignment of each
reading frame. The traceback needs memory in proportion to query length × frame length
(4.2 bytes per cell, measured), so a frame beyond 200 million cells is aligned in
overlapping pieces (overlap: twice the query length) and the best alignment of each
piece is returned; up to 1.1 Mb of window for a 550-aa query nothing changes. A query
longer than 10,000 aa is not searched with Smith-Waterman (a warning says so): the
memory of one piece grows with the square of the query length. A missing
Smith-Waterman backend stops the run (`allow_missing_smith_waterman = false`).
With `sw_max_alignments_per_piece` above 1 (default 1), a piece is aligned again after
each alignment, with the residues of that alignment masked in the frame (`*`), and
further alignments are returned while their E-value is at most 1, up to that number per
piece: a second exon or a second copy of the gene in the same frame is otherwise hidden
behind the best alignment. An alignment that runs through a masked span ends the piece;
one that ends on the last residue of a piece is left to the next piece, which holds it
whole.

MMseqs2 runs with E ≤ 10 and tblastn with E ≤ 1000; the limit of 0.1 is applied when
their output is read. tblastn's cutoff is wide on purpose: through its composition-based
statistics (on by default) the cutoff also changes the E-values tblastn reports, in either
direction, so with a cutoff of 10 a weak true gene could vanish from its output. The hits of the three engines are pooled in order of bit score, and a hit is
dropped when more than half of it overlaps a hit already kept. The genome-wide hits
that lie in the window are then added (GOI hits only for the two queries above).

### 7.7 From hits to gene models

For each GOI query:

1. **Loci.** Hits less than 40 kb apart form a locus (more when the home gene has
   long introns: 5 × its largest intron). The best `goi_candidate_regions_per_query`
   loci are kept (2; the `local_family` module raises it to the home copies + 1, at most
   10), ranked by query coverage, then best bit score, then number of hits.
2. **Gene models** (`annotate_models_from_hit_list` → `miniprot_gene_models`). miniprot
   aligns the query to the locus ± `gap_search_window` (50 kb; at least the home
   gene's span × expansion), with a maximum intron of that stretch (at least 20 kb)
   and sensitive settings (`-n 2 -p 0.3 -N 50 --outs=0.3`). Queries below 30 aa are
   aligned without splicing. **Every model is kept** whose coding sequence overlaps
   the coding sequence of a better model by less than half (`exon_annotation`; ids
   `…_l1_…`, `…_l1m2_…`, …; at most `goi_max_models_per_locus`, 10). While a hit of
   the locus lies on no coding exon of a kept model, the coding exons found so far
   are masked and miniprot runs again; a model of a later round must reach 30 % of
   the best score. A model's protein is miniprot's own translation.
3. **Tandem copies** (`detect_tandem_duplications`). When two or more hits each cover
   at least 35 % of the query (or three or more hits mostly cover the same part of
   it), each such hit that lies on no model is written as a one-segment copy
   (`tandem_copy`). When miniprot returns no model, all of them are.
   With `tandem_query_restart` (off by default, batch G2) copies are built the way the
   hit chain is: per strand, in coding order, a hit extends a copy when it continues the
   query (10 residues of overlap at most) within the GOI's intron bound. A copy counts at
   a union query coverage of ≥ 0.35 with a best hit of ≥ 20 aligned residues or ≥ 30 bits
   (both for a query under 150 aa); copies that share genome are one (a second program's
   hit of the same exon cannot fake a copy); it is a tandem when two counted copies overlap
   in the query by more than half of the shorter one (an inverted copy counts). Each copy
   is one mRNA row with one CDS per hit and the union coverage, and is judged on its hits
   with their stops. Without the switch two exons of one gene that each cover 35 % of the
   query were a tandem.
4. **Start and stop** (`refine_model_termini`). miniprot reports the aligned part
   only. The first exon is extended upstream, in frame, to the nearest `ATG`, and the
   last exon downstream to the first stop codon. The walk covers at most the number
   of unaligned query residues plus 3 codons (10 codons when the query is aligned to
   its end; never more than 30 upstream, 20 downstream). It does not start at a
   splice site (`AG` before the first base, `GT` after the last), does not cross a
   stop codon, and changes nothing when it finds no start or stop.
   With `goi_first_exon_guard` (off by default, batch G3A) the walk is skipped where a
   short home first exon is missing: the home gene has ≥ 2 exons, its first is at most
   30 aa, and the query residues before the model's first CDS are within 5 aa of that
   length (for an exon of 8 aa or less: at least max(4, length − 2)). The nearest ATG then
   lies in the intron. The model keeps miniprot's start and carries
   `MissingQueryStart=1-<p−1>` and `FirstExon=missing`; the same at the stop end
   (`MissingQueryEnd`, `LastExon=missing`). Only for the home query and its decoys, not
   for the relative's query.
5. **Hit chain** (`fallback_hit_span`), when miniprot returns no model. The hits on
   the majority strand are chained in genomic order so that their query positions
   only rise (or only fall, on the minus strand). Two hits are not linked when they
   overlap by more than 30 bp on the genome, or when the gap between them exceeds
   5 × the home gene's largest intron (at least 10 kb). The chain is accepted when:
   - query below 150 aa: aligned length ≥ 15 aa and best bit score ≥ 30;
   - longer query, one hit: coverage ≥ 0.25 and length ≥ 25 aa;
   - longer query, several hits: coverage ≥ 0.25, or total length ≥ 35 aa, or bit
     score ≥ 60;
   - and its genomic span is at most 600 × the query length in bp (at least 100 kb)
     for queries below 150 aa, or 150 × (at least 5 kb) for longer ones.

   The identity of a chain is the mean over its hits, weighted by alignment length.
   The mRNA row lists the query interval of each hit (`ChainQuery=1-40,41-88,...`).

   With `goi_chain_query_aware` (off by default, batch G1) the chain follows the query
   instead: a hit continues the chain only when it starts at most 10 residues before the
   previous hit ends and ends after it (on the minus strand the mirror image); the residues
   it skips must fit into the genomic gap (3 bp each, less 30 bp); a gap under 40 bp joins
   only two pieces of one exon (same reading frame, the gap within 15 bp of the skipped
   codons), which are then written as one CDS row. The chain with the most bits wins (ties:
   more query covered, then fewer hits). The label is still computed on the hits
   themselves, and the chain still carries the bits of its best single hit.
6. **Single hits.** Remaining hits can be written as `rescued_exon` or `raw_hit`.
   Both get an E-value in their place and are judged like every GOI candidate (§7.8); a
   `raw_hit` needs 60 % query coverage to be written, so it can be HIGH.
   With `goi_join_leftover_hits` (off by default, batch G3B) a `rescued_exon` hit is first
   offered to a gene model of its own query that the window's judge calls HIGH or MEDIUM
   (the home query or its decoy, not the relative's query): the hit must lie outside the
   model on its strand and continue the query from the model's first or last CDS by the
   rule of `goi_chain_query_aware` (at most 10 residues of overlap, the skipped residues
   fit into the gap, within the GOI's intron bound), with no CDS of another gene model or
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
hit's hits translated again with their stops) and its hit regions. The specificity uses
that same protein (review fix RF-1e, 2026-10-15; before, both used the record protein,
for a chain the join of its hits with the stops deleted).
The candidates of a window are judged first and then ranked in the judge's order
(`synvoy_label.rank_key`, §7.8: label, E-value, the block's E_syn, bits, the shorter span; a
chain over two copies scores what its better copy scores, and must not win over that copy);
one that overlaps an accepted candidate by more than half of its own length is dropped.
Until 2026-10-14 they were ranked by identity, longer span first, and only hit chains had a
bit score (that of their best hit). The strict fallback on single hits ranks hits by bit
score too. (The report's dedup chose by bits too since 2026-10-14; on the Candida pilot that
picked the model on the curated gene in 29 of the 31 genome pairs with several models,
identity in 16.)

### 7.8 Confidence

Since 2026-10-14 (label rework, step 4) one function judges every GOI candidate, wherever
one is judged: the choice among overlapping candidates in a window, the label, the
report's dedup, the hull rescue, the figures (`bin/synvoy_label.py`, `judge`). Its inputs
are the E-value of the home protein against the model in the place SynVoy expected the
gene (§7.13), the model's specificity against the home proteome (§7.12) and its query
coverage:

| Confidence | `Basis` | Rule |
|---|---|---|
| **HIGH** | `significant_specific` | Eeff ≤ 1e-6; margin ≥ 10 bits or no real competitor; coverage ≥ 0.7 |
| MEDIUM | `significant_specific` | Eeff ≤ 1e-3 and specific |
| AMBIGUOUS | `not_specific` | Eeff ≤ 1e-3, not specific: a relative where the gene should be |
| AMBIGUOUS | `gene_order_only` | 1e-3 < Eeff ≤ 1, in the slot of a block with E_syn ≤ 0.1, and specific |
| LOW | `not_significant` | everything else |

- **Coverage** is `QueryCoverage`, the share of the home protein the model covers. With
  `core_completeness` (the `signal_anchor` module) a gene model also carries
  `CoreCoverage`, the share of the protein's core (without the signal peptide and the GPI
  tail found by `query_features.py`) covered by its exons' query intervals, and the judge
  reads the larger of the two; 0.7 stays the bar. The hull rescue reads the query span of
  its alignment, without the core.
- **Specific**: no real competitor, or a margin of at least 3.3 bits over the best other
  home protein. Within ±3.3 bits a model in the slot or core of a significant block is
  decided for the gene (positional tie-break), unless that protein lies in the gene's own
  home neighbourhood: in a tandem array position cannot tell the two apart.
- AMBIGUOUS calls ("needs curation") are kept, drawn and reported, never counted. A
  gene-order-only call is offered only in a genome where no home locus and no rescue pass
  produced a counted call of the gene; counted means after the ownership and coverage
  demotions of the report ([§10](#10-adjudication-and-the-report)). Where the report drops one, the per-genome
  GFF and the figures leave it out too (`apply_report_drops.py`; the figures are drawn
  after the report).
- **Slot.** Without the home locus span (`goi=` in the anchor headers) no block has a
  slot; the search then warns once and judges in the core, the window or the genome.
- **Without parasail** no `Bits=` is computed: every candidate type is ordered by the best
  bit score of its own search hits and gets no specificity terms (one warning). Until
  2026-10-15 hit chains kept the bits of their best hit and every other type had 0.
- **Precision.** The judge and the seed bar work on the full-precision values; the GFF
  holds E-values with 3 significant digits and `SpecMargin` at 0.1 bit, for reading only
  (review fix RF-1f, 2026-10-15: an E of 1.0004e-6 was read back as 1e-06 and judged HIGH
  at the bar of 1e-6, a margin of 9.96 bits as 10.0 and allowed to seed).
- **Ranking among candidates** (`rank_key`): label, then E-value (values below 1e-10
  count as equal), then the block's E_syn, then bits, then the shorter span, then
  position. In a window the label is provisional, on the E-value of the place before the
  genome's multiplier; the final label is set at the end of the genome on Eeff.
- **Frame guard**: a miniprot model's E-value is computed on the translation of its CDS
  (first exon's phase honoured). When that has a stop before its end or a length that is
  not a multiple of 3, another frame over the same exons is used if it is clean
  (`FrameFixed=1`), else only the model's hit regions are judged (`FrameBroken=1`).
- The hull rescue judges its models with the same function and the same inputs: the
  label protein through the frame guard (a broken frame is judged on the model's own
  translation in pieces between stops, since the rescue has no search hits), completeness
  as the query span of the alignment, specificity with the GOI's home neighbourhood (the
  anchors of the locus found in the genome). It has no block E-value, so its place is the
  whole genome.
- A model carries a provisional LOW (`InferenceReason=awaiting_label`) from the moment it
  is built; every GOI row with an E-value is relabelled at the end of its genome. A GOI
  row the judge cannot label (no Eeff) keeps the LOW, never counts and never seeds, and
  is counted in the log as `LOW/unjudged`. Until 2026-10-14 the label came from identity
  thresholds: HIGH at ≥ 50 % identity with two exons and two flanking genes in home
  order, MEDIUM at ≥ 40 %, AMBIGUOUS for hit chains in blocks of ≥ 5 flanking genes, and
  a context promotion (one day on by default, removed) lifted weak models on gene order
  alone. Its parameters (`classify_high_min_identity` and its kin) were removed on
  2026-10-15.
- ML scores (`--enable_plm_search`, `--enable_structural_search`) are recorded on the
  rows (`EmbeddingSimilarity`, `StructuralSimilarity`); they lift a label by one step only
  on a row the judge did not label (`_ml_relabel`).
- A model also gets a status that does not depend on confidence: `fragment` below
  0.4 coverage, `complete` from 0.7 coverage with at least 2 exons, else `partial`.

### 7.9 Thresholds relaxed for distant targets (removed)

Until 2026-10-15 the HIGH and MEDIUM identity thresholds of the earlier label were lowered
for a genome whose flanking genes matched poorly (the §1g distance auto-tune: up to 10
points below a median flanking identity of 70 %, down to floors of 45 and 30 %). The
judge of §7.8 has no identity threshold, and its E-value is computed in the expected place
of each genome, so the relaxation went with the thresholds. The median of the best
identity of each flanking gene is still computed and logged per genome
(`Median flanking identity: … over N flanking gene(s)`); it changes nothing.

### 7.10 Flanking gene models

Flanking genes found in the window are modelled the same way (miniprot, else a hit
chain) and labelled by `_classify_flanking_evidence`: a miniprot model is HIGH with
at least 2 exons or identity ≥ 55 %, otherwise MEDIUM; a hit chain is MEDIUM with
coverage ≥ 0.65 and identity ≥ 55 %, otherwise LOW. One model per flanking gene and
locus is kept across blocks.

A flanking gene that no block contains is searched on the other scaffolds and, if
found, written as a rearranged flanking gene, so the figures can show it.

### 7.11 GOI hits far from their anchors

`select_dispersed_goi_seeds` (on by default) handles a rearranged neighbourhood. A
genome-wide GOI hit with identity ≥ 40 %, E ≤ 1e-10 and length ≥ 50 becomes a block
of its own when:

- its scaffold carries at least 3 distinct flanking genes,
- it lies between the outermost of them, widened by the envelope margin
  ([Distances](#distances-in-gene-units): 9 pitches), and
- it is not already within the GOI margin of a flanking gene (a block keeps the GOI hits
  within the GOI margin of it; until 2026-10-22 the radius here was the window margin, so a
  hit between the two was used by neither path).

A strong GOI hit that fails these conditions is written to
`hits/<genome>.nonsyntenic.tsv`, so the report can say "found, not syntenic".

With `goi_genome_wide_path` (off by default; FINAL step 5, batch of 2026-10-15) a GOI hit
outside every searched window is modelled without these gates, in a genome that has at
least one block at E_syn ≤ 0.1:

- the GOI hits outside every window are grouped per chromosome into places (hits within
  the GOI's intron bound of each other);
- a place qualifies when the home protein's blastp E against its hits (in their frames, in
  pieces between stops) over the whole genome (`-searchsp` = |q| × genome / 3) is at most
  1e-3, the label's MEDIUM bar, and the specificity holds in the genome place: no real
  competitor, or a margin of at least 3.3 bits (§7.8). A place another home protein
  explains as well is not modelled, so it never becomes "needs curation";
- at most N places per home GOI and genome are modelled, best E first, N = the query's
  home family members + 1: the home proteins other than the query's own gene that the
  label's family test calls family members (blastp E ≤ 1e-3 with composition-based
  statistics and SEG; only proteins whose Smith-Waterman score gives a Karlin-Altschul
  E ≤ 0.1 are tested);
- each place is searched in a window of its hits ± the model margin, outside every other
  window, without decoys and without the flanking pass. Its candidates are judged in the
  genome place (no multiplier), carry `GenomeWide=1`, and one that another home protein
  explains as well is not kept.

Every place is listed in `hits/<genome>.genome_wide.tsv` with its E over the genome, its
specificity terms, its rank, the cap and what became of it (`modelled`, `cut_by_cap`,
`not_specific`, `not_significant`); the report's `genome_wide_path` block lists them. The
genome place has no decoy calibration: decoys are searched in block windows only.

With `goi_genome_wide_without_blocks` as well (PR-9) the path also runs in a genome without
a block at E_syn ≤ 0.1, and in one where no block is found at all (the search then goes on
with no block instead of stopping; the cross-chromosome flanking recovery is skipped there).
A genome-wide call carries `RetrocopyCandidate=1` when one of its CDS, or one hit of the home
query in it, aligns 10 residues on both sides of an exon junction of the home gene
(annotated, ≥ 2 exons): the target encodes both exons there without an intron, as a
retrocopy does (and an ortholog that lost the intron).

### 7.12 What is added to the query set

**Specificity** (`bin/synvoy_specificity.py`, written on every kept GOI model). The
model's protein x (the protein its `Bits=` come from) is aligned (Smith-Waterman, BLOSUM62, gaps 11/1) to every protein of
the home proteome (one per gene, the longest isoform). The best one other than the GOI's
own home gene is p* (`BestOtherHome`); `SpecMargin` = bits(home protein, x) − bits(p*, x).
The GOI's own home gene is the gene the home step took as the query's (`gene_id` in
`goi_info.json`, found by name or by position, [§5.6](#56-the-home-gene-model)). Without
one, or when the home step flagged its choice as ambiguous (with a warning: the gene it
took may be a relative of the query's own gene), it is the home protein near-identical to
the query (≥ 95 % identity over ≥ 90 % of it). Until 2026-10-16 (review fix RF-2b) only that
identity rule existed, so where the home annotation differs from the query (an extra exon,
a partial model) the GOI's own gene was p*.
p* is a real competitor (`RealCompetitor=1`) when it aligns to x at E ≤ 0.01 over the
whole proteome, or when it is a home family member of the GOI (the two align at E ≤ 1e-3,
blastp with composition-based statistics and SEG, the GOI as the query: in a pairwise run
`-seg yes` masks the low-complexity stretches of the query, and the composition-based
statistics mask the subject's, so both sides are masked (with `-seg no` only the subject's
were). The test is not symmetric: in 55 of 949 home pairs the two directions disagree at
1e-3 under `-seg no`, 34 under `-seg yes`; a failed blastp call is logged and counts as no
family).
`BestOtherInNeighbourhood=1` when p* lies between the outermost flanking genes of the home
locus (the span extended to the locus itself), or among the GOI's own home gene and its
`n_flanking_genes` neighbours on each side, also in a run from a secondary home locus:
there, position cannot tell the GOI from p*. Until 2026-10-15 (review fix RF-2) only the
first counted, so flanking genes on one side of the locus, or a secondary locus, let the
positional tie-break decide for a tandem neighbour of the GOI. The search is exhaustive: on the Candida pilot a prefilter (MMseqs2,
top 15) missed two home paralogs. The proteome and its gene positions reach the search
inside the home MMseqs2 database folder (`home_proteome.faa`, `home_genes.tsv`). The
label (§7.8) and the seed rule below use these terms.

**Seeds** (`_classify_goi_for_seed_and_tree`), per query in a genome and home locus,
after all blocks are labelled: the HIGH calls, or the MEDIUM calls when there is no
HIGH call; of these, the ones with no real competitor or a margin of at least 10 bits.
They join the query set of the following waves, and the window search's relative query
is chosen among them. Every other HIGH or MEDIUM call is passed to the tree only;
AMBIGUOUS, LOW and context-promoted calls never seed. Until 2026-10-14 a call seeded
when it was HIGH or MEDIUM of class `confident_goi` or `probable_goi`: a tandem copy
never did, and a relative's MEDIUM next to a HIGH call in the same genome did (apyrase:
three genomes, margins of −339 to −363 bits); the opt-in `--seed_on_flanking_support`
was removed. A seeded single-segment tandem copy is, by default (`tandem_copy_parent
seed`), recorded as `GOI_copy_<n>` and therefore a gene of its own in the following waves,
judged against the seed; with `--tandem_copy_parent home` it is a model of the home gene
and can be its relative query.

The window search keeps the hits of the home protein and of the relative's model apart
when it removes overlapping hits (since 2026-10-14): pooled, the relative's hit at a
relative's gene removed the home protein's own hit there.

### 7.13 E-values in the expected place

Since 2026-10-14 (label rework, step 3) every kept GOI model also carries the evidence
the new label will use (`bin/synvoy_evalue.py`):

- **`BlockEsyn`**: the expected number of windows of this length in the genome that hold,
  by chance, a run of the block's length (`BlockCollinearSupport`) of anchors in home
  order, either direction. Each home neighbour with m places in the genome lies in a
  random window of length W with probability 1 − exp(−mW/G); the number present is a
  Poisson binomial; for n present anchors in random order the chance of a monotone run
  ≥ k is computed exactly (Robinson–Schensted–Knuth). A block with E_syn ≤ 0.1 predicted
  the gene's place.
- **`Place`**: the smallest of slot (between the nearest found home neighbours on either
  side of the gene, widened by one gene pitch; the side of each neighbour comes from the
  home locus span, `goi=` in the anchor headers), core (block ± GOI margin) and window
  that holds the model, in a block with E_syn ≤ 0.1; else `genome`.
- **`EPlace`**: blastp of the home protein against the model's protein and against each of
  its search-hit regions (cut at the aligned span, translated in the hit's frame, split
  at stop codons), `-comp_based_stats 2 -seg no -word_size 2 -threshold 11`, with the
  search space set explicitly to |q| × ⌊place/3⌋ (`-searchsp`; `-dbsize` let cut-out hit
  regions look rarer than they are), read at full precision (JSON). **`EPlaceTblastn`**:
  the best tblastn row of the home protein that reaches into the model (at least 15 bp of
  it in its own frame inside the model, the rule of the hit regions above; until 2026-10-15
  the row had to lie wholly inside), as tblastn reported it for the window (before the keep
  floor of 0.1 and the dedup against the other engines' hits), its E rescaled from the
  window to the place. The window tblastn runs once into an archive that `blast_formatter`
  writes twice: as the tabular output the search reads (byte-identical to a direct run) and
  as JSON, from which this E is read at full precision (tabular output prints E with 3
  decimals from 0.0009 to 0.1, so 0.00123 read as 0.001). The judge does not use it.
- **`Eeff`**, **`EeffTblastn`**: the two multiplied by `SigBlocks`, the number of
  significant blocks searched in this genome (each is one predicted place); a model in
  the genome-wide place is not multiplied.

Every block searched, with or without a call, is written with its E_syn at full precision,
its neighbourhood and the final labels of its GOI models to `blocks/<locus>/<genome>.blocks.tsv`
(since 2026-10-19; it changes no call). With `--plot_block_selection true` the figures choose
from it (§8).

**Copies in the slot** (`--report_copies_in_slot`, off; the `local_family` module turns it
on). After the genome's final labels, every GOI model of a block that has a slot gets
`CopiesInSlot=`: the GOI models with `Place=slot` and a label of AMBIGUOUS or better in that
slot. Each of them gets `SlotRank=`, counted from the copy nearest the home neighbour on
the left (the slot's orientation comes from its two neighbours). The report pairs the home
copies (`home_copies` of `run_settings`) with the copies in the slot of each genome's best
call (`copies_in_slot`, e.g. `2:3`). Both decide nothing: equal rank is not positional
orthology (Dewey 2011), since each gain or loss in an array shifts every rank behind it.

**Decoy mode** (`--decoy_mode`, off): in every window the reversed home protein and two
shuffles of it (seeded by the query id) are searched like the home protein, modelled,
labelled and judged with the same terms, and written to
`decoys/<locus>/<genome>.decoy.gff` and `.decoy.tsv`. Each decoy case (home gene × decoy
kind) is chosen among its own candidates: pooled, a call of one kind suppressed an
overlapping call of another, so the rate was a lower bound (rates from before 2026-10-15
are not comparable). A window that also searches the relative's query (`<gene>~rel`)
searches its reversed and shuffled copies too (`~rev~rel`, `~shuf1~rel`, `~shuf2~rel`, since
2026-10-15): like a real model of the relative's query, their models are measured against
the home decoy of the same kind, filed under it (`SearchQuery=` names the decoy query) and
chosen in its case, so the rates cover both paths of the search. Each decoy row also
carries `PlaceBp=` and its E-value in every nested place that holds it (`EPlaceSlot`,
`EPlaceCore`, `EPlaceWindow`, `EPlaceGenome`, before `SigBlocks`), so a rule can be
calibrated in any place. They never seed and never reach the
result GFF. Every decoy call is false, so their rate per label and place calibrates the
labels without a truth set. The window search takes about four times as long. The report
sums the decoy cases that were searched (home gene × target genome with a search window ×
decoy kind; the home genome itself is never searched) into `decoy_rate`: per label the
calls and the upper end of the exact 95 % interval. Decoy mode is a run-level option; no
preset module sets it.

## 8. Region scores

`cluster_grs.py` scores candidate regions per genome. They are used for the figures
and for the strong-synteny rescue ([§9](#9-rescue-passes)). They do not change the
confidence of a gene model. With `--plot_block_selection true` (off by default) a figure row
shows the search's blocks instead: the blocks of the HIGH, MEDIUM and AMBIGUOUS GOI calls,
then the significant blocks without a counted call by E_syn, up to 6 (fragment calls take
part, hidden from view); a genome without its block table keeps the regions.

**Clusters.** The genome-wide hits are put into one fixed order and clustered: a hit
joins the cluster when it starts less than the block gap of the search for this
genome after the previous hit ends (read from the result GFF; 150 kb for a GFF
without the record).

**Score.** With `N` flanking genes at home:

```
C = max(U, B) / N        coverage: U = distinct flanking genes in the cluster,
                          B = flanking support of the search block of a GOI model
                          that the cluster overlaps (0 without one)
K = L / k                order: L = longest collinear run of the k ranked hits
S = max(a, d) / (a + d)  strand: a hits agree with the home strand, d disagree

quality       = 0.4·C + 0.3·K + 0.3·S
synteny_score = quality · C
score         = synteny_score + 0.15     if the cluster overlaps a GOI call
```

A cluster that overlaps a GOI model also takes over the span of that model's search
block. Identity is not part of the score: the score describes the neighbourhood, the
confidence of §7.8 describes the gene.

**p-value** (`estimate_pvalue`). The gene labels of the cluster's hits are redrawn 200
times, with replacement, from the labels of all hits in the genome (fixed seed 42),
and `synteny_score` is recomputed each time. `p = (count(null ≥ observed) + 1) / 201`.

**Order and selection.** Clusters are sorted by `B` (descending), then score, then p.
Clusters that overlap a GOI model are placed first. With `max_regions = 0`, regions
are written from the top until one has both a score below the floor
(`max(0.03, 0.30 × best score)`) and fewer than 3 distinct flanking genes; at least
one and at most 6 are written.

**Strong synteny without a gene.** A cluster holding at least
`strong_synteny_min_flanking` (⌈`n_flanking_genes`/2⌉ = 5) HIGH flanking models but no GOI call is written as
`goi_missing_but_strong_synteny`. Up to two such regions are added beyond the cap,
unless they lie inside a selected GOI region.

**Dispersed hit.** A GOI-overlapping cluster with at most one flanking gene gets a
score of at least 0.45 (`goi_dispersed_macrosynteny`) when its GOI hit has identity
≥ 40 % and E ≤ 1e-10, at least one flanking gene lies on the same scaffold, and the
genome's rearrangement score is at least 0.5. That score (`target_rearrangement_score`)
is the mean of two fractions: flanking genes off the main scaffold, and pairs of
home-adjacent flanking genes that are not within 2 Mb of each other in the target.

Region confidence: HIGH from `min_synteny_score` (0.6), MEDIUM from half of it, else
LOW.

## 9. Rescue passes

Both passes run miniprot with the normalised query on a window chosen from the
flanking genes alone. Both apply the start/stop refinement of §7.7.

| | Strong-synteny rescue | Hull rescue |
|---|---|---|
| Script | `rescue_strong_synteny.py` | `rescue_goi_hull.py` |
| Runs when | a region is `goi_missing_but_strong_synteny` | a scaffold has ≥ `goi_hull_min_flanking` (round(0.4 × `n_flanking_genes`) = 4) HIGH flanking models in the hull and no HIGH GOI model lies in the DNA that would be searched |
| Window | the region ± 2 kb | the hull ± the GOI margin; skipped when the hull is longer than twice the home neighbourhood × expansion |
| Accepts | every miniprot model covering ≥ 5 % of the query | every model with identity ≥ 40 % and CDS ≥ 0.5 of the query length (cost filters in front of the judge, kept by decision of 2026-10-15); the one the judge ranks first |
| Confidence | always LOW | the judge of the main search (§7.8): E-value in the genome-wide place, specificity with the GOI's home neighbourhood, completeness = query span of the alignment |
| Evidence type | `relaxed_miniprot_rescue` | `synteny_hull_rescue` |

**The hull** (`compute_hull`). The HIGH flanking models of a scaffold are clustered
with the bridge limit of the search as gap limit. The centre of the largest cluster is
taken, and the hull is the span of all flanking models within one home neighbourhood
(× expansion) of it. It covers a gene that lies between two blocks. The three sizes
come from the distances the search recorded for the genome ([Distances](#distances-in-gene-units));
for a GFF without that record they are 3 Mb, 100 kb and 20 Mb, the fixed values used
before.

The rescue models are published under `rescue/locus_N/` and reach the report.

## 10. Adjudication and the report

`generate_report.py` → `build_report` reads all loci and genomes once and applies, in
this order:

1. **Deduplication** (`dedupe_goi_annotations`). HIGH, MEDIUM and AMBIGUOUS calls of
   one genome and scaffold are one gene when they overlap reciprocally by more than
   80 %. The call the judge ranks first is kept (`_record_rank`: label, then
   `synvoy_label.rank_key`, i.e. E-value in its place, the block's E_syn, bits, span; bits
   or identity for a record without `Eeff`, from a run before 2026-10-14 or a rescue
   pass), with the list of home loci that found it. The figures choose their
   representative GOI model by label, then the same strength. A model contained in a wider one replaces it when the
   wider model is at least twice as long and adds less than 0.15 query coverage.
2. **Phylogenetic placement** (`phylo_placement_check.py`). Across all calls, sequence
   divergence (1 − identity to the query) is fitted against the taxonomic distance of
   the species with a Theil-Sen line. A call whose residual, scaled by the median
   absolute deviation, is 2 or more above the line is flagged `phylo_discordant`: it
   has diverged more than its species should have, which is what a paralog looks
   like. The check needs calls from at least 4 species and is advisory.
3. **Locus ownership** (`build_home_paralog_panel.py`, `build_locus_ownership`). A
   panel of home proteins is built: the annotated gene at each home locus, plus up to
   15 proteins of the home proteome that align to the query with identity ≥ 28 % and
   a Smith-Waterman score of at least a quarter of the query's self-score. Each
   recovered protein is aligned to the panel. If its best match is a paralog and not
   the GOI, the call is relabelled `paralog_not_goi`. If the two best scores differ
   by less than 10, the locus with more HIGH flanking models decides.
4. **Coverage demotion** (`apply_coverage_demotion`). A HIGH or MEDIUM call with
   identity ≥ 50 % and a recorded coverage below 0.35 is relabelled
   `identity_coverage_decoupled`.
5. **Gene-order-only calls** (`drop_gene_order_only`). An AMBIGUOUS call with
   `Basis=gene_order_only` is dropped in a genome that has a counted call (HIGH or MEDIUM,
   not relabelled in steps 3–4); the dropped calls are listed in
   `goi_dedup.gene_order_only_dropped_records`. Until 2026-10-15 this ran inside step 1, so
   a genome whose only call steps 3–4 then relabelled had lost its candidate as well.
6. **Headline.** HIGH and MEDIUM calls are counted without the two relabelled classes.
   AMBIGUOUS and LOW calls are reported separately.
7. **Self-consistency flags**: strong synteny without a GOI, the same gene found from
   several home loci, and the flags of steps 2–4.

## 11. Tree

`compute_tree.py` aligns the GOI proteins of §7.12 with MAFFT (`--auto`) and builds a
tree with IQ-TREE (model selection `-m MFP`; 1000 ultrafast bootstrap replicates from
4 sequences). With fewer than 3 sequences a placeholder tree is written.

## 12. Figures

`plot_synteny.py` draws one set of figures per home locus from the home annotation,
the models of [§7](#7-the-search-in-one-target-genome) and the regions of
[§8](#8-region-scores). The figures are drawn after the report: `apply_report_drops.py`
copies the per-genome GFFs into `plot_inputs_*` without the gene-order-only calls the
report dropped ([§10](#10-adjudication-and-the-report), step 5), and the figures are drawn
from those copies.

- **Labels.** A GOI call is drawn as resolved when its label is HIGH or MEDIUM, whatever
  its class (a HIGH tandem copy was drawn as unresolved, an AMBIGUOUS call as resolved
  until 2026-10-15); AMBIGUOUS calls are drawn in the needs-curation style. The best GOI
  call of a row, the overlap dedup and the GOI cap rank by label, then by the strength of
  the match (E in its place, the block's E_syn, bits).

- **Home row.** The GOI is drawn from the home annotation: the transcript whose CDS
  covers the hits best. Only coding exons are drawn.
- **Target rows.** A scaffold is mirrored when at least 70 % of its flanking genes
  are on the opposite strand to home and their order agrees.
- **Synteny plot.** Genes are stacked in at most 3 lanes; fragment models are hidden
  in the default figure and shown in a `_with_fragments` variant.
- **Anchor grid.** One row per species, one column per home gene. `×N` above a
  flanking arrow counts copies on the same scaffold within the neighbourhood; a GOI
  cell with more than 3 models shows the best one and `×N`.

## 13. Reproducibility

The same inputs give the same output files:

- the search inside a window runs single-threaded (`deterministic_goi_search`);
- the results of a wave are joined in genome order and sorted by id;
- sets are sorted before they are iterated, and `PYTHONHASHSEED` is fixed to 0;
- region hits are put into one total order before scoring;
- the permutation test uses a fixed seed.

## 14. Properties and known limits

Consequences of the rules above that a reader of the results should know.

**Confidence**

- HIGH needs a query coverage of at least 0.7, whatever the number of exons (until
  2026-10-14 a one-exon model could be HIGH only when the home gene had one exon).
- The hull rescue's models are judged like the main search's (§7.8). Its 40 % identity
  and 0.5 CDS-coverage gates only decide which models are judged.
- Flanking support is a property of the window (every candidate in a window has the
  same value, also a copy at the edge of the window), and the label does not read it.
- Two GOI row types are not judged: `relaxed_miniprot_rescue` rows of the strong-synteny
  rescue (on by default; written LOW, with no `Eeff` and no `Basis`; LOW is never
  counted) and `structural_discovery` rows (opt-in; HIGH or MEDIUM by TM-score and
  embedding similarity). Whether to judge them is open.
- `Place=` is the plan's `SearchSpace=`; it will be renamed (with `PlaceBp` and the
  `EPlace*` fields), the old names readable for one release.
- A genome can hold several HIGH models for one home locus: the copies of a tandem
  array, or a relative inside the searched window. All of them are counted in the
  headline as "GOI ortholog annotations"; the pipeline does not say which one is the
  ortholog. Within a window, `l1` is the locus with the best coverage and score.
- Identity and coverage refer to the home protein, for every genome. Runs made before
  2026-10-09 could report them against the model of an earlier genome, which gave
  higher values for late genomes.

**Search**

- A stretch between two blocks that is neither bridged nor merged is not searched by
  the main path. The hull rescue covers it only when the scaffold has at least
  `goi_hull_min_flanking` (4) HIGH flanking models.
- Every block that passes the filter and the cap of §7.4 is searched; there is no
  early stop. Most windows hold no GOI model, and each can add LOW rows.
- A genome whose MMseqs2 search fails for lack of memory, after retries, is skipped
  with a warning and reads as "no hits".
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
  rough number (33 % off in one case).
- The pitch is the mean step between anchors, so one very long gap at home (a gene
  desert next to the locus) raises every derived distance of that locus.
- The cost of a search grows with the pitch. A block needs two anchors within the
  block gap, and chance alone puts pairs of anchor hits that close; each such block is
  searched, and its window is as large as the pitch makes it. The derived distances
  have been run on yeast, fly, bee, ant, wasp, beetle and plant genomes (up to 0.6 Gb).
  At a vertebrate home locus (pitch of several hundred kb, block gap about 2 Mb) they
  have not been run; `-profile legacy_distances` gives the former fixed distances there.
- The number of chance blocks grows with the number of places a flanking gene hits. In
  an ant genome that is about 2 per gene; in *Arabidopsis lyrata* a flanking gene of
  *A. thaliana* has 92 hits, the cap of 80 blocks per genome is reached, and one home
  locus cuts out a fifth of the genome (the FT orthologs are still found).
  `--anchor_max_places 5` (§7.1) brings that to one window.
- Where the gene order between home and target is lost (no five anchors in a row), the
  expansion is the ratio of the genome sizes, a coarse stand-in for the spacing of
  genes.
- The effect of the growing query set ([§6](#6-target-order-quality-check-and-waves))
  has been measured on one test set, where switching it off gave the same result.

**Not built: a family mode**

For a gene family the pipeline searches one member at a time and reports the copies it
meets on the way: every gene model of a window is kept (§7.7), and a home locus that
holds several copies is recognised (§5.3). The family is not yet the unit. A family
mode would need three things:

1. every member of the family in the home genome as a query (the home-paralog panel
   of §10 already finds them);
2. each target model assigned to its closest home member (the ownership step of §10
   computes this, today only to demote a call that belongs to a paralog);
3. a result per family member and genome instead of one headline count.

An older switch goes part of the way: `--expand_goi_similar true` adds the relatives
next to the GOI as further queries (`GOI_NEIGHBOR_…`). It is off since 2026-04-25,
because in a large family the extra queries filled the candidate regions: on the test
set of that time the right scaffold was recovered in 0 of 9 genomes with it and in 8
of 9 without.

**Scores**

- `synteny_score = 0.4·C² + 0.3·K·C + 0.3·S·C`. Coverage dominates; the weights are
  not shares of the score.
- The permutation test cannot produce the `max(U, B)` term, so p is too small when
  `B > U`. The p-value is reported and breaks ties in the region order; it does not
  select regions and does not enter any gene's confidence.
- Region selection stops at the first region below both floors, but the list is not
  sorted by score alone.

**Reporting**

- The figures and the tree are built from the search output. They do not include
  rescue models, and they still show calls that the report relabelled.
- `RECIPROCAL_BEST_PARALOG` still runs but cannot flag anything: it needs several
  query sequences and always receives one. Paralogs are handled by locus ownership
  ([§10](#10-adjudication-and-the-report)), which is not affected.
- Locus ownership compares a call with proteins of the *home* genome. A paralog that
  exists only in the target lineage is not in that panel; the only check that can
  flag it is phylogenetic placement, which is advisory.

**Easy mode**

- Each taxonomy level requests the record of every assembly under the taxon. A level
  with several thousand assemblies can exceed the 600 s limit and is then skipped, so
  the run continues with fewer genomes.

## 15. Default values

| Parameter | Default | Used in |
|---|---|---|
| `search_evalue` | 0.01 | §5.2, §7.1 |
| `mmseqs_sensitivity` | 9.5 | §5.2, §7.1 |
| `min_query_length` | 20 (warning below 40) | §5.1 |
| `max_loci` / `locus_min_bit_ratio` | 5 / 0.5 | §5.4 |
| `locus_gap_genes` / `locus_gap` | auto (⌊n_flanking_genes/2⌋ = 5) / auto | §5.3, §5.4 |
| `home_locus_mode` | auto | §5.5 |
| `n_flanking_genes` | 10 per side | §5.7 |
| `goi_relative_evalue` | 1e-5 | §5.7 |
| `rank_wave_binning` | true | §6 |
| `min_hit_identity` / `min_hit_length` | 10 / 10 | §7.1 |
| `max_intron` | 20,000 (floor of the per-gene limit) | §7.2 |
| `block_gap_genes` / `cluster_distance` | 7 pitches / auto | §7.2, §8 |
| `synteny_bridge_max_gap_genes` / `synteny_bridge_max_gap` | auto (2 × n_flanking_genes = 20 pitches) / auto | §7.3 |
| `synteny_bridge_min_anchors` / `_max_rank_gap` / `_max_per_block` | 2 / 5 / 2 | §7.3 |
| `window_margin_genes` / `region_padding` | 7 pitches / auto | §7.4, §7.5 |
| `goi_margin_genes` / `goi_margin` | 2 pitches / auto | §7.5 |
| `min_block_genes` / `max_blocks_per_genome` | 2 / 80 | §7.4 |
| `gap_search_window` | 50,000 (floor) | §7.7 |
| `goi_max_models_per_locus` | 10 | §7.7 |
| `sw_min_score` / `sw_min_identity` | 20 / 10 | §7.6 |
| `goi_fallback_intron_margin` / `_floor` | 5 / 10,000 | §7.7 |
| `fallback_short_query_len` / `_min_aln_aa` / `_min_bits` | 150 / 15 / 30 | §7.7 |
| `classify_fragment_max_qcov` / `classify_complete_min_qcov` | 0.4 / 0.7 | §7.8 |
| `dispersed_goi_min_identity` / `_max_evalue` / `_min_alnlen` | 40 / 1e-10 / 50 | §7.11 |
| `dispersed_envelope_margin_genes` / `dispersed_envelope_margin` | 9 pitches / auto | §7.11 |
| `synteny_weight_base` / `_consistency` / `_strand` | 0.4 / 0.3 / 0.3 | §8 |
| `synteny_goi_overlap_bonus` | 0.15 | §8 |
| `adaptive_score_floor_frac` / `_abs` / `adaptive_max_regions` | 0.30 / 0.03 / 6 | §8 |
| `strong_synteny_min_flanking` | auto (⌈n_flanking_genes/2⌉ = 5) | §8, §9 |
| `min_synteny_score` | 0.6 | §8 |
| `goi_hull_min_flanking` / `_min_identity` / `_min_coverage` | auto (round(0.4 × n_flanking_genes) = 4) / 40 / 0.5 | §9 |
| `goi_hull_window_pad` / `_max_window` / `_cluster_max_gap` | auto / auto / auto | §9 |
| `paralog_panel_min_identity` / `_max_family` | 28 / 15 | §10 |
| `locus_ownership_tiebreak_gap` | 10 | §10 |
| `identity_decoupled_min_identity` / `_max_qcov` | 50 / 0.35 | §10 |
| `phylo_placement_min_calls` / `_z_threshold` | 4 / 2.0 | §10 |
