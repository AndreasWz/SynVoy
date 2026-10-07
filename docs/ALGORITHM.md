# SynVoy: the algorithm

This document describes what SynVoy computes, step by step, with the rule and the
default value at every decision. It is written from the code, not from earlier
descriptions: every statement was checked against the source of branch `dev` on
2026-10-07 (commit `18300a9` plus the block-merge fix of that day). Code locations
are given as `file` → `function`, because line numbers drift.

It does not explain how to run the pipeline ([USAGE.md](USAGE.md)), what the output
files contain ([OUTPUT.md](OUTPUT.md)) or what every parameter does
([PARAMETERS.md](PARAMETERS.md)). The history of each rule, with the measurements
behind it, is in [STATE_OF_THE_PROJECT.md](STATE_OF_THE_PROJECT.md) and
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
| **Hit** | One local alignment of a query protein to genomic DNA. |
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
| 5 | Split into home loci | `SPLIT_LOCI` (optionally `RESOLVE_HOME_LOCUS`) | `split_loci.py`, `resolve_home_locus.py` | up to 5 loci |
| 6 | Model the home gene | `ANNOTATE_GOI` | `annotate_goi_exons.py` | GOI protein, exons, largest intron |
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
     always pass. Any other assembly is dropped when its best N50 is below 5 kb, when
     a known contig or scaffold count exceeds 500,000, or when it is contig-level and
     has no quality figures at all (`bad_quality_policy = drop`).

## 5. The home locus

### 5.1 Query normalisation

`normalize_query.py` upper-cases the sequence and removes whitespace, digits, gap
characters and a trailing `*`. Other non-sequence characters are an error. A
nucleotide query (DNA or RNA) is translated in six frames and the longest open
reading frame is taken. A protein shorter than `min_query_length` (30) stops the run.

### 5.2 Locating the query

`LOCATE_GENE` searches the home genome twice and pools the result:

- MMseqs2 `easy-search`, translated search, sensitivity 9.5, E ≤ 0.01;
- `tblastn` with default settings.

`merge_hits.py` keeps hits with E ≤ `search_evalue` (0.01), sorts them by position
and merges a hit into the previous one when both are on the same scaffold and strand
and it starts less than 1 kb after the previous one ends. A merged interval carries
the lowest E-value and the highest bit score of its members.

### 5.3 Parameter preset

`profile_hits.py` summarises the raw hits: hits within 100 kb form one locus; a locus
is *strong* when its best bit score is at least half of the top locus's; `bit_ratio`
is the top locus's bit score divided by the second's. `auto_select_preset.py` →
`select_preset` then applies the first rule that matches:

| Order | Condition | Preset |
|---|---|---|
| 1 | no hits | `preset_single_copy` |
| 2 | query shorter than 100 aa | `preset_short_peptide` |
| 3 | more than 15 hits, or more than 5 strong loci | `preset_paralog_discrimination` (with a warning) |
| 4 | at least 2 strong loci and `bit_ratio` < 2 (or undefined) | `preset_tandem_family` |
| 5 | otherwise | `preset_single_copy` |

What each preset changes (`resolve_effective_params.py` → `PRESET_OVERRIDES`):

| Preset | Main changes from the defaults |
|---|---|
| `preset_short_peptide` | HIGH ≥ 40 % identity, MEDIUM ≥ 25 %; hits kept from 8 % identity and 8 aa; Smith-Waterman from score 10; fragment below 0.3 coverage, complete from 0.5 |
| `preset_single_copy` | HIGH ≥ 55 %, MEDIUM ≥ 40 %; complete from 0.75 coverage |
| `preset_tandem_family` | flanking genes rejected from 25 % similarity to the GOI; MEDIUM ≥ 30 %; tandem copies from 35 % |
| `preset_paralog_discrimination` | HIGH ≥ 55 %, MEDIUM ≥ 40 %; tandem copies from 50 %; at most 4 regions per genome with a higher score floor |

The chosen values travel to the later steps as a `settings` map of 17 parameters. A
value the user changed at launch is kept over the preset. `--auto_apply_preset false`
leaves the defaults in place; `--preset_override` picks a preset by hand.

### 5.4 Splitting into home loci

`split_loci.py`:

1. Sort the merged intervals by position. Start a new locus when the scaffold changes
   or the gap to the previous interval exceeds 50 kb.
2. Order the loci by their best E-value. The first is the primary locus.
3. Keep a further locus only if its best bit score is at least
   `locus_min_bit_ratio` (0.5) of the primary's.
4. Keep at most `max_loci` (5), the ones with the highest bit scores.

More than 10 loci before filtering produce a large-family warning.

### 5.5 Optional: locus by gene name

With `--enable_name_locus` or `--home_goi_gene SYMBOL` (default off),
`resolve_home_locus.py` looks the gene symbol up in the home GFF and replaces the loci
of §5.4 with that one annotated gene. The match is accepted only when the gene's
protein and the query are at least 60 % identical in a Smith-Waterman alignment;
otherwise the run continues with the loci of §5.4.

### 5.6 The home gene model

`annotate_goi_exons.py` builds the GOI's own exon structure, trying in this order:

1. **Annotation by name** (`match_goi_in_gff`). The UniProt gene names of the query
   and the query id, with case variants and the `gene-` prefix, are looked up in the
   home GFF.
2. **Annotation by sequence.** The proteins of the annotated genes within 100 kb of
   the hits are compared with the query; the best match is taken.
3. **miniprot** on the hits ± 50 kb, when no annotated gene matches.
4. **Query protein only**, when there are no hits to work from.

Cases 1 and 2 take the exons from the GFF's CDS rows (`method: gff_annotation`). The
output is the GOI protein, one sequence per exon, and `goi_info.json` with the gene's
genomic span and its largest intron.

### 5.7 Flanking genes

`extract_flanking_genes.py`, per home locus:

1. Take the annotated gene whose centre is closest to the locus centre.
2. Walk outwards on each side. Collect genes until `n_flanking_genes` (10)
   protein-coding genes are found on that side, looking at no more than 50 genes per
   side.
3. Skip a gene whose protein is too similar to the GOI: longest common subsequence,
   divided by the longer length, at or above `max_flanking_goi_similarity` (35 %).
   Such a gene is a family member and would match everywhere.
4. Genes without CDS are kept for the figures but do not count towards the 10.

Without a home GFF, a substitute annotation is built first: genes are predicted
around the locus (`prodigal_on_regions.py`, Augustus or Prodigal) and annotations are
transferred from annotated target genomes (`borrow_annotations.py`). If that yields
nothing either, `extract_flanking_genes.py` predicts genes within 50 kb of the locus
and keeps the longest tenth, at least 21.

Each flanking gene is written as its full protein and, with `exon_level_search` (on),
as one sequence per exon.

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
are sorted closest first.

**Quality.** `assess_genome_quality.py` fails a genome with N50 below 5 kb or more
than 500,000 sequences. Failed genomes are removed (`qc_fail_policy = drop`).

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

Only hits of flanking genes are used to build blocks. GOI hits are not, with one
exception ([§7.11](#711-goi-hits-far-from-their-anchors)). If no flanking gene has a
hit, all hits are used.

### 7.2 Blocks

`identify_synteny_blocks`:

1. **Per gene.** Hits of the same flanking gene (its exons count as the gene) that lie
   less than `max_intron` (20 kb) apart form one *gene locus*. A larger gap starts
   another locus, which is a second copy or a paralog.
2. **Along the scaffold.** Gene loci are sorted by position. A locus joins the current
   block when it starts less than `cluster_distance` (150 kb) after the previous one
   ends; otherwise the gap is tested for bridging (§7.3), and if that fails a new
   block begins.
3. **Per block.** The collinear run is computed from one position per flanking gene
   (its first), ordered along the scaffold. Blocks are sorted by collinear run, then
   by the number of distinct flanking genes.

### 7.3 Bridging a gap

A rearrangement or an insertion can put several megabases between two parts of one
neighbourhood. A block continues across a gap when all of these hold (`_can_bridge`):

- same scaffold, and the gap is below `synteny_bridge_max_gap` (6 Mb);
- the block has bridged fewer than `synteny_bridge_max_per_block` (2) gaps so far;
- the block already holds at least `synteny_bridge_min_anchors` (3) distinct flanking
  genes;
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
| `A B` … `D E F` | two blocks (2 anchors) |

`synteny_bridge_two_sided` (off) weighs the whole cluster behind the gap instead of
its first gene. `--disable_synteny_collinearity` turns bridging and the collinear
ranking off.

### 7.4 Merging, filtering, capping

1. **Merge** (`merge_synteny_blocks`). Two blocks on one scaffold are merged when
   their spans, each widened by `region_padding` (150 kb), overlap. Blocks less than
   300 kb apart therefore become one. The merged block holds the union of the parts'
   flanking genes, so a gene with two copies counts once.
2. **Filter.** Blocks with fewer than `min_block_genes` (2) flanking genes are
   dropped. A block is kept regardless when it was seeded by the GOI's own home gene
   or by a rescued GOI hit (§7.11). If no block is left, all blocks are used.
3. **Cap.** At most `max_blocks_per_genome` (80) blocks are searched: GOI-seeded
   blocks first, then by collinear run, then by gene count.

### 7.5 Windows

`process_region_block`, per block:

- **Search window** = block ± padding. The padding is twice the mean gap between
  consecutive hits on that scaffold, limited to 50–200 kb (150 kb when the scaffold
  has fewer than two hits). Only this sequence is searched.
- **Core window** = block ± `gap_search_window` (50 kb). A hit that lies entirely
  outside it is discarded, so every GOI model starts from a hit within 50 kb of the
  flanking genes, although the sequence searched reaches further.

Two numbers describe the block and go into every confidence decision:

- **Flanking support** `B`: the number of distinct flanking genes with a genome-wide
  hit in the search window.
- **Collinear support**: the collinear run of those genes.

### 7.6 Search inside the window

The queries are the full-length GOI proteins in the query set: the home GOI and the
models added by earlier waves, one per origin (the longest).

| Engine | Query | Kept when |
|---|---|---|
| MMseqs2, one thread | the protein and its halves, thirds and quarters (pieces ≥ 15 aa) | identity ≥ 15 %, length ≥ 15, E ≤ 0.1 |
| tblastn (`-seg no`) | the protein | identity ≥ 15 %, length ≥ 15, E ≤ 0.1 |
| Smith-Waterman (parasail), six frames | the protein | score ≥ 20, identity ≥ 10 %, length ≥ 10, E ≤ 1 |

Smith-Waterman uses BLOSUM62 with gap costs 11/1. Its E-value is the Karlin-Altschul
value for the six-frame length of the window. A missing Smith-Waterman backend stops
the run (`allow_missing_smith_waterman = false`).

MMseqs2 and tblastn run with E ≤ 10; the limit of 0.1 is applied when their output is
read. The hits of the three engines are pooled in order of bit score, and a hit is
dropped when more than half of it overlaps a hit already kept. The genome-wide hits
that lie in the window are then added.

### 7.7 From hits to gene models

For each GOI query:

1. **Loci.** Hits less than 40 kb apart form a locus. The two best loci are kept,
   ranked by query coverage, then best bit score, then number of hits.
2. **Tandem copies** (`detect_tandem_duplications`). When two or more hits each cover
   at least 35 % of the query, every such hit is one copy (`tandem_copy`), not an exon
   of one gene. The same holds for three or more hits when most pairs of them cover
   the same part of the query.
3. **Gene model** (`annotate_using_miniprot`). Otherwise miniprot aligns the query to
   the locus ± 50 kb, with a maximum intron of the window length (at least 20 kb) and
   sensitive settings (`-n 2 -p 0.3 -N 50 --outs=0.3`). Queries below 30 aa are aligned
   without splicing. The best model is taken (`exon_annotation`). Its protein is
   miniprot's own translation.
4. **Start and stop** (`refine_model_termini`). miniprot reports the aligned part
   only. The first exon is extended upstream, in frame, to the nearest `ATG`, and the
   last exon downstream to the first stop codon. The walk covers at most the number
   of unaligned query residues plus 3 codons (10 codons when the query is aligned to
   its end; never more than 30 upstream, 20 downstream). It does not start at a
   splice site (`AG` before the first base, `GT` after the last), does not cross a
   stop codon, and changes nothing when it finds no start or stop.
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
6. **Single hits.** Remaining hits can be written as `rescued_exon` or `raw_hit`.
   Both are always LOW.

All candidates of a window are ranked, and one that overlaps an accepted candidate by
more than half of its own length is dropped.

### 7.8 Confidence

`_classify_goi_evidence` assigns each GOI model a confidence and a class. `flank` is
the flanking support `B`; identities are in percent.

| Evidence | Condition | Confidence, class |
|---|---|---|
| `exon_annotation` | ≥ 2 exons, identity ≥ 50, flank ≥ 2, order conserved | **HIGH**, `confident_goi` |
| | the same, order not conserved | MEDIUM, `probable_goi` |
| | identity ≥ 35 and (flank ≥ 1 and the model is not a fragment, or coverage ≥ 0.65) | MEDIUM, `probable_goi` |
| | otherwise | LOW |
| `fallback_hit_span` | flank ≥ 2, coverage ≥ 0.75, identity ≥ 60 | MEDIUM, `probable_goi` |
| | flank ≥ 5, identity ≥ 30, and (coverage ≥ 0.25 or identity ≥ 35) | **AMBIGUOUS**, `syntenic_candidate_unconfirmed` |
| | otherwise | LOW |
| `tandem_copy` | identity ≥ 40, coverage ≥ 0.35, and (flank ≥ 5 or coverage ≥ 0.65) | MEDIUM, `tandem_goi_copy` |
| | otherwise | LOW, `tandem_goi_copy` |
| `rescued_exon`, `raw_hit` | always | LOW |

- **Order conserved** means: collinear support ≥ `classify_high_min_collinear` (3).
  The condition is waived when the block has fewer than 3 flanking genes or the home
  order is unknown.
- **AMBIGUOUS** marks a candidate that is supported by its neighbourhood only. A
  conserved neighbourhood looks the same whether the gene is there or was lost, so
  such a call is kept and drawn but not counted as an ortholog.
  `--disable_ambiguous_tier` labels it MEDIUM.
- A model also gets a status that does not depend on confidence: `fragment` below
  0.4 coverage, `complete` from 0.7 coverage with at least 2 exons, else `partial`.

### 7.9 Thresholds relaxed for distant targets

Before classification, `_apply_distance_adaptive_thresholds` takes the best identity
of each flanking gene in this genome and computes the median `m`. With at least 3
flanking genes, the HIGH and MEDIUM identity thresholds are lowered by

```
relax = 0                               if m ≥ 70
relax = 10 · (70 − m) / (70 − 40)       if 40 < m < 70
relax = 10                              if m ≤ 40
```

percentage points. Targets with `m ≤ 40` are marked for manual review.

The 70 and 40 are identities of the *flanking* genes, not limits for the GOI. With
the default thresholds, the lowest GOI identity that can still be HIGH is 40 %, and
MEDIUM 25 %.

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
- it lies between the outermost of them, widened by 200 kb, and
- it is not already within `region_padding` of a flanking gene.

A strong GOI hit that fails these conditions is written to
`hits/<genome>.nonsyntenic.tsv`, so the report can say "found, not syntenic".

### 7.12 What is added to the query set

`_classify_goi_for_seed_and_tree`:

- A GOI model with confidence HIGH or MEDIUM and class `confident_goi` or
  `probable_goi` is added to the query set of the following waves.
- A HIGH or MEDIUM model of another class (in practice a tandem copy) is not added,
  but is passed to the tree.

## 8. Region scores

`cluster_grs.py` scores candidate regions per genome. They are used for the figures
and for the strong-synteny rescue ([§9](#9-rescue-passes)). They do not change the
confidence of a gene model.

**Clusters.** The genome-wide hits are put into one fixed order and clustered: a hit
joins the cluster when it starts less than `cluster_distance` (150 kb) after the
previous hit ends.

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
`strong_synteny_min_flanking` (5) HIGH flanking models but no GOI call is written as
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
| Runs when | a region is `goi_missing_but_strong_synteny` | a scaffold has ≥ 4 HIGH flanking models in the hull and no HIGH GOI model overlaps it |
| Window | the region ± 2 kb | the hull ± 100 kb; skipped above 20 Mb |
| Accepts | every miniprot model covering ≥ 5 % of the query | the most identical model with identity ≥ 40 % and coverage ≥ 0.5 |
| Confidence | always LOW | by identity only: HIGH from 50 %, MEDIUM from 35 % (the shipped values; a preset's thresholds do not apply here) |
| Evidence type | `relaxed_miniprot_rescue` | `synteny_hull_rescue` |

**The hull** (`compute_hull`). The HIGH flanking models of a scaffold are clustered
with a gap limit of 3 Mb. The centre of the largest cluster is taken, and the hull is
the span of all flanking models within 10 Mb of it. It covers a gene that lies
between two blocks.

The rescue models are published under `rescue/locus_N/` and reach the report.

## 10. Adjudication and the report

`generate_report.py` → `build_report` reads all loci and genomes once and applies, in
this order:

1. **Deduplication** (`dedupe_goi_annotations`). HIGH, MEDIUM and AMBIGUOUS calls of
   one genome and scaffold are one gene when they overlap reciprocally by more than
   80 %. The call with the higher confidence, then identity, is kept, with the list
   of home loci that found it. A model contained in a wider one replaces it when the
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
5. **Headline.** HIGH and MEDIUM calls are counted without the two relabelled classes.
   AMBIGUOUS and LOW calls are reported separately.
6. **Self-consistency flags**: strong synteny without a GOI, the same gene found from
   several home loci, and the flags of steps 2–4.

## 11. Tree

`compute_tree.py` aligns the GOI proteins of §7.12 with MAFFT (`--auto`) and builds a
tree with IQ-TREE (model selection `-m MFP`; 1000 ultrafast bootstrap replicates from
4 sequences). With fewer than 3 sequences a placeholder tree is written.

## 12. Figures

`plot_synteny.py` draws one set of figures per home locus from the home annotation,
the models of [§7](#7-the-search-in-one-target-genome) and the regions of
[§8](#8-region-scores).

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

- A single-exon gene cannot be HIGH on the main path, because HIGH requires at least
  2 exons (fixed in the code, not a parameter). The rule does not look at how many
  exons the home gene has.
- The hull rescue assigns confidence from identity alone. Exon count and gene order
  are not checked there, and it uses the shipped thresholds even when a preset raised
  them for the main path. A gene between the two thresholds is MEDIUM on the main
  path and HIGH through the rescue.
- Together: an intronless gene is MEDIUM at best on the main path, which leaves "no
  HIGH GOI" on its scaffold, so the hull rescue runs when the scaffold has at least 4
  HIGH flanking models. Its HIGH calls then come from the path without the order
  check. With fewer flanking models the call stays MEDIUM.
- Flanking support is a property of the window. Every candidate in a window has the
  same value.

**Search**

- A stretch between two blocks that is neither bridged nor merged is not searched by
  the main path. The hull rescue covers it only when the scaffold has at least 4 HIGH
  flanking models.
- `max_consecutive_empty_blocks` has no effect: the code logs the streak and
  continues.
- The hull rescue checks for an existing HIGH GOI model inside the span of the
  flanking genes, but searches that span plus 100 kb on each side. A gene the main
  path already found just outside the span is therefore modelled a second time. The
  report merges the two; when they tie, the main-path model is the one reported.
- A genome whose MMseqs2 search fails for lack of memory, after retries, is skipped
  with a warning and reads as "no hits".
- Only the first sequence of a multi-sequence query is used (with a warning).
- `gap_min_size`, `gap_evalue` and `gap_min_alnlen` are accepted but not read. There
  is no gap-filling search; missing exons are found by miniprot inside its window or
  not at all.
- All distances (`cluster_distance`, the paddings, the bridge limit) are in base
  pairs and do not scale with gene density. In a compact genome a window holds many
  more genes than in a large one.
- The effect of the growing query set ([§6](#6-target-order-quality-check-and-waves))
  has been measured on one test set, where switching it off gave the same result
  ([STATE_OF_THE_PROJECT.md](STATE_OF_THE_PROJECT.md), F4).

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
| `min_query_length` | 30 | §5.1 |
| `max_loci` / `locus_min_bit_ratio` | 5 / 0.5 | §5.4 |
| `n_flanking_genes` | 10 per side | §5.7 |
| `max_flanking_goi_similarity` | 35 | §5.7 |
| `rank_wave_binning` | true | §6 |
| `min_hit_identity` / `min_hit_length` | 10 / 10 | §7.1 |
| `max_intron` | 20,000 | §7.2 |
| `cluster_distance` | 150,000 | §7.2, §8 |
| `synteny_bridge_max_gap` | 6,000,000 | §7.3 |
| `synteny_bridge_min_anchors` / `_max_rank_gap` / `_max_per_block` | 3 / 5 / 2 | §7.3 |
| `region_padding` / `padding_min` / `padding_max` | 150,000 / 50,000 / 200,000 | §7.4, §7.5 |
| `min_block_genes` / `max_blocks_per_genome` | 2 / 80 | §7.4 |
| `gap_search_window` | 50,000 | §7.5, §7.7 |
| `sw_min_score` / `sw_min_identity` | 20 / 10 | §7.6 |
| `goi_fallback_intron_margin` / `_floor` | 5 / 10,000 | §7.7 |
| `fallback_short_query_len` / `_min_aln_aa` / `_min_bits` | 150 / 15 / 30 | §7.7 |
| `classify_high_min_identity` / `classify_medium_min_identity` | 50 / 35 | §7.8 |
| `classify_high_min_collinear` | 3 | §7.8 |
| `classify_tandem_min_identity` / `_min_qcov` | 40 / 0.35 | §7.8 |
| `classify_fallback_strong_min_identity_floor` | 30 | §7.8 |
| `classify_fragment_max_qcov` / `classify_complete_min_qcov` | 0.4 / 0.7 | §7.8 |
| `distance_autotune_close_pct` / `_far_pct` / `_max_relax` | 70 / 40 / 10 | §7.9 |
| `dispersed_goi_min_identity` / `_max_evalue` / `_min_alnlen` | 40 / 1e-10 / 50 | §7.11 |
| `synteny_weight_base` / `_consistency` / `_strand` | 0.4 / 0.3 / 0.3 | §8 |
| `synteny_goi_overlap_bonus` | 0.15 | §8 |
| `adaptive_score_floor_frac` / `_abs` / `adaptive_max_regions` | 0.30 / 0.03 / 6 | §8 |
| `strong_synteny_min_flanking` | 5 | §8, §9 |
| `min_synteny_score` | 0.6 | §8 |
| `goi_hull_min_flanking` / `_min_identity` / `_min_coverage` | 4 / 40 / 0.5 | §9 |
| `goi_hull_window_pad` / `_max_window` / `_cluster_max_gap` | 100,000 / 20,000,000 / 3,000,000 | §9 |
| `paralog_panel_min_identity` / `_max_family` | 28 / 15 | §10 |
| `locus_ownership_tiebreak_gap` | 10 | §10 |
| `identity_decoupled_min_identity` / `_max_qcov` | 50 / 0.35 | §10 |
| `phylo_placement_min_calls` / `_z_threshold` | 4 / 2.0 | §10 |
