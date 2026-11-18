# SynVoy parameters

Every parameter of SynVoy: what it does, its default and when to change it. The defaults are
those in `nextflow.config`. How to start a run and choose a profile: [USAGE.md](USAGE.md).
What each step computes: [ALGORITHM.md](ALGORITHM.md).

**Setting a parameter.** On the command line with two dashes (`--max_genomes 10`), in a params
file (`-params-file run.yaml`) or in a config file (`-c my.config`, as
`params.max_genomes = 10`). Write a switch as a word (`--skip_tree false`), and a value that
starts with a minus sign in quotes (`--modules 'auto,-short_query'`).

**What else changes a default.**

- **The profile.** `./run_synvoy.sh` without `-profile` runs `-profile auto,low_mem`. The memory
  tiers `low_mem` and `laptop_safe` change these values
  ([USAGE.md § 2](USAGE.md#memory-tiers)):

  | Parameter | Default | `low_mem` | `laptop_safe` |
  |---|---|---|---|
  | `mmseqs_sensitivity` | 9.5 | 7.0 | 8.0 |
  | `mmseqs_split_memory_limit` | 8G | 2G | 4G |
  | `mmseqs_verbosity` | 1 | 0 | 0 |
  | `iterative_search_cpus` | 2 | 1 | 1 |
  | `iterative_search_memory` | auto | 4 GB | 8 GB |
  | `locate_gene_memory` | 3 GB | 2 GB | 2 GB |
  | `skip_tree` | false | true | false |

- **The preset modules.** After SynVoy has found the gene in the home genome, the modules set up
  to 25 parameters for that gene ([§ 1](#1-presets-and-preset-modules)). The entries below say
  which module can change them.
- **Your values win.** A value you set is kept over a module and a preset, as long as it
  differs from the default in `nextflow.config`.

What a run used is recorded with its results: the 32 settings that modules and presets decide,
with the source of each value, in `intermediate/locate_gene/effective_params.json`
([OUTPUT.md](OUTPUT.md#run_settings)). A parameter that SynVoy no longer has is listed under
[Removed parameters](#20-removed-parameters).

## Contents

1. [Presets and preset modules](#1-presets-and-preset-modules)
2. [Query and mode](#2-query-and-mode)
3. [Easy mode: target genomes](#3-easy-mode-target-genomes)
4. [Pro mode: your files](#4-pro-mode-your-files)
5. [The home locus](#5-the-home-locus)
6. [Flanking genes](#6-flanking-genes)
7. [Distances in gene units](#7-distances-in-gene-units)
8. [Genome-wide search and blocks](#8-genome-wide-search-and-blocks)
9. [Search inside a window](#9-search-inside-a-window)
10. [Gene models and labels](#10-gene-models-and-labels)
11. [Gene prediction for a home genome without a GFF](#11-gene-prediction-for-a-home-genome-without-a-gff)
12. [Region scores](#12-region-scores)
13. [Rescue passes](#13-rescue-passes)
14. [Ownership, paralog check and placement](#14-ownership-paralog-check-and-placement)
15. [Optional steps: language model, structure, discovery](#15-optional-steps-language-model-structure-discovery)
16. [Figures and tree](#16-figures-and-tree)
17. [Resources](#17-resources)
18. [Output and console](#18-output-and-console)
19. [Calibration and comparison switches](#19-calibration-and-comparison-switches)
20. [Removed parameters](#20-removed-parameters)

---

## 1. Presets and preset modules

After SynVoy has found the gene in the home genome, it sets up to 25 parameters for that
gene: by default with the **preset modules**, under `--modules legacy` with one of four
**presets**. A value you set yourself is kept. The rules, with every formula and cap:
[ALGORITHM.md § 5.3](ALGORITHM.md#53-parameter-preset-and-preset-modules).

### `modules`
**Default:** `auto`

Which preset modules are applied.

| Value | Applies |
|---|---|
| `auto` | every module whose trigger fires; an empty value counts as `auto` |
| `'auto,-dispersed_family,+cys_framework'` | `auto` without, or with, the modules named |
| `'short_query,local_family'` | exactly these modules, whether their trigger fires or not |
| `none` | no module |
| `legacy` | no module: the preset chooser picks one of the four presets and applies it ([below](#the-preset-chooser-under---modules-legacy)) |

Write a value with a minus sign in quotes: Nextflow reads a bare `-x` as `true`. A module can be
named by its id or its name (`A1` or `local_family`).

| Id | Module | Fires when | Sets |
|---|---|---|---|
| A1 | `local_family` | the query covers its home locus 1.5 times or more, or a family member of the gene lies among the `n_flanking_genes` home genes on either side (found only with a home GFF and blastp) | `goi_relative_evalue`, `goi_candidate_regions_per_query`, `goi_max_models_per_locus`, `goi_margin_genes`, `sw_max_alignments_per_piece`, `report_copies_in_slot`, `locus_ownership_max_genes_per_locus`, `max_goi_per_genome` |
| A2 | `dispersed_family` | the query aligns strongly at two or more home loci, and a locus other than the best covers at least half of it | `adaptive_score_floor_frac`, `adaptive_score_floor_abs`, `adaptive_max_regions`, `paralog_panel_max_family` |
| A3 | `shared_domain` | the query aligns strongly at two or more home loci, and every other locus covers less than half of it | nothing (recorded only) |
| B1 | `short_query` | the query is shorter than 100 aa | `classify_complete_min_qcov`, `classify_fragment_max_qcov`, `window_tblastn_word_size`, `window_tblastn_comp_stats`, `window_mmseqs_min_length`, `miniprot_short_orf_max_len`, `aug_relaxed_length_min`, `fallback_short_min_aln_aa`, `dispersed_goi_min_alnlen`, `goi_hull_min_coverage`; and a second home search without tblastn's low-complexity filter when the query has no hit in the home genome |
| B2 | `low_complexity` | 10 % or more of the query lies in low-complexity windows, or 25 % or more in repeated 4-mers | `goi_relative_statistic`, `report_low_complexity` |
| B3 | `signal_anchor` | the query has a signal peptide or a GPI-anchor signal | `core_completeness` |
| B4 | `cys_framework` | six or more cysteines lie within 90 residues | nothing (recorded only) |
| C1 | `exon_structure` | the home gene has two or more exons | `sw_max_alignments_per_piece`, `dispersed_goi_min_alnlen` |
| C2 | `home_context` | the home gene lies within 40 kb of the end of a scaffold of 1 Mb or more | nothing (recorded only) |

A strong locus scores at least half as much as the best one. Several modules can fire for one
gene. Each parameter has a single rule, so two modules never set it to different values. No
module changes the rules of the label; a module can still change a call's label, because it
changes the evidence the label reads (B3, for example, also measures completeness on the core of
the protein, [ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)).

Every mode except `legacy`, `none` included, also lets the search set `anchor_max_places` per
target genome ([§ 8](#anchor_max_places)) and, for a home genome without a GFF, chooses the
Augustus model nearest the home species ([§ 11](#augustus_species)). Under `legacy` neither
happens.

The console lists the applied modules in one line, each with the reason it fired
(`[INFO] modules (auto), applied: …`), then the values their rules set. `effective_params.json` records the modules, every value and its source
(`default`, `module:A1`, `preset:<name>` or `user_set`), and `user_set_kept` lists the values of
yours that were kept against a module or a preset.

### `auto_apply_preset`
**Default:** `true`

`false` applies no module and no preset, not even one named with `--preset_override`: every one
of the 32 settings keeps its default or the value you set. Use it to run with the plain
defaults. What `--modules` does outside these settings stays: `anchor_max_places` per genome,
the Augustus model, and the second home search of `short_query`.

### `preset_override`
**Default:** none

Applies one preset by name, in every mode: under `--modules legacy` instead of the chooser's
pick, otherwise on top of the modules, whose values it replaces where both set one. A value you
set wins over both. `preset_single_copy` adds nothing.

### Presets

| Preset | Chosen under `--modules legacy` when | Sets (default) |
|---|---|---|
| `preset_short_peptide` | the query is shorter than 100 aa | `classify_complete_min_qcov 0.5` (0.7), `classify_fragment_max_qcov 0.3` (0.4): the coverage bands of a model's status, not its label |
| `preset_paralog_discrimination` | the query aligns strongly at two or more home loci | `adaptive_score_floor_frac 0.45` (0.30), `adaptive_score_floor_abs 0.08` (0.03), `adaptive_max_regions 4` (6): fewer candidate regions per genome |
| `preset_tandem_family` | the query covers its home locus 1.5 times or more | `goi_relative_evalue 1e-3` (1e-5): weaker relatives of the gene are not used as anchors either |
| `preset_single_copy` | no hits, or none of the above | nothing |

The module of the same kind (B1, A2, A1, none) sets the same values, and B1, A1 and A2 set more.

A preset can also be given as a profile, after the execution profile:
`-profile auto,low_mem,preset_short_peptide`. Its values then count as values you set: they win
over the modules, and the modules still set everything else. `-profile preset_short_peptide`
alone replaces the launcher's default `auto,low_mem`. The profile files are in `conf/presets/`
and are named differently: `short_divergent_peptide.config`, `paralog_discrimination.config`,
`tandem_family.config` and `single_copy_conserved.config`; one can also be passed with `-c`.

### The preset chooser under `--modules legacy`

The chooser reads the query's hits in the home genome and takes the first rule that matches: no
hits, `preset_single_copy`; a query shorter than 100 aa, `preset_short_peptide`; two or more
strong loci, `preset_paralog_discrimination`; a home locus covered 1.5 times or more,
`preset_tandem_family`; otherwise `preset_single_copy`. Under the other modes its pick is only
recorded (`legacy_preset` in `effective_params.json`). The console names the preset:

```
[INFO] §1f preset applied: preset_paralog_discrimination (auto_preset.json chose preset_paralog_discrimination)
```

Under the default `--modules auto` the same line reads `preset applied: none (modules auto: …)`,
followed by the line of the modules.

### The 32 settings

Modules and presets act on 32 parameters, which travel together from the settings step to the
later steps: the 25 a module can set (listed above) and 7 that no module sets
(`sw_min_score`, `sw_min_identity`, `min_hit_identity`, `min_hit_length`, `expand_goi_similar`,
`goi_family_tokens`, `adaptive_unique_gene_floor`). The minimum query length is not among them:
it is checked before, and no module or preset changes it.

---

## 2. Query and mode

### `mode`
**Default:** `easy` (`easy` or `pro`)

- `easy`: SynVoy fetches the query protein (by accession, gene symbol or sequence), the home
  genome with its annotation and the target genomes from UniProt and NCBI, and chooses the
  targets by taxonomy ([§ 3](#3-easy-mode-target-genomes)).
- `pro`: you give the files ([§ 4](#4-pro-mode-your-files)). SynVoy downloads no query and no
  genome. It still uses the network when it can, to order the target genomes by NCBI taxonomy;
  without it, it orders them by a sequence distance. The first run also needs the network to
  build the tool environment.

Easy mode with `--home_genome`, `--home_gff` or `--target_genomes` stops at its start and asks for
`--mode pro`: easy mode would download other genomes and leave these files unused.
`--target_gffs` in easy mode gives a warning (easy mode fetches the annotations with the genomes).

Both modes check the target FASTA files before the search ([`qc_fail_policy`](#qc_fail_policy)).

### `query`
**Default:** none

A FASTA file with the query: a protein, or a DNA sequence that SynVoy translates (the longest
open reading frame in six frames). Only the first sequence of the file is used, with a warning.
Pro mode needs it. In easy mode it replaces the download, and then `--home_species` is required,
because a file names no organism.

Give an absolute path: `run_synvoy.sh` starts the pipeline from the SynVoy folder, so a relative
path is read from there, not from the folder you started the launcher in.

The first word of the header names the gene in every output. It may not contain `|exon_`,
`|frag_`, `|slide_`, a `|` or any of `,` `;` `=` `&` `%`. The usual database forms are read as
their accession: `sp|P12345|NAME` and `tr|…|…` (UniProt), `ref|XP_1.1|`, `gb|…`, `emb|…`,
`dbj|…`, `lcl|X`.

### `query_id`
**Default:** none (easy mode)

An accession or a gene symbol. SynVoy reads the source from the format: a UniProt accession
(`P01501`) is fetched from UniProt, an NCBI protein accession (`XP_006565763.1`, version
optional) from NCBI Protein. Anything else is taken as a gene symbol and looked up in UniProt
for the species given with `--home_species`. From the record SynVoy takes the sequence, the
source organism (the home species, unless you set `--home_species`), the protein name and, from
UniProt, the gene symbol.

### `query_seq`
**Default:** none (easy mode)

The query as text on the command line: a sequence or a whole FASTA record. `--home_species` is
required with it. Pro mode ignores it and needs `--query`.

Give one query source. With two, easy mode uses the first in the order `--query_seq`, `--query`,
`--query_id`, and warns when one of the two is `--query`; with all three, the run stops at the
start (unless `--query_id` holds a sequence).

### `home_species`
**Default:** none (easy mode: the organism of the query's record)

The species whose genome is the home genome, as an NCBI Taxonomy name (`"Apis mellifera"`).

- **Easy mode.** Set it when the query names no organism (always with `--query_seq` or a local
  `--query` file) or to use another species than the record's. SynVoy downloads that species'
  reference genome, else its best-ranked assembly, else, with a warning, the genome of its
  closest relative that has one.
- **Pro mode.** It orders the target genomes. Without it, SynVoy takes the home genome's file
  name as the species name (`Apis_mellifera.fna`), with a warning. A file name that is an
  assembly accession (`GCF_…`, `GCA_…`) or a generic word (`genome`, `reference`, `home`, …)
  stops the run at the start, unless you set `--allow_unknown_species true`.

### `allow_unknown_species`
**Default:** `false` (pro mode)

`true` lets a pro-mode run start when no home species is given and none can be read from the
home genome's file name. The targets are then ordered without the home species' taxonomy.

### `min_query_length`
**Default:** `20` (aa)

A shorter query stops the run before the search; `0` switches the check off. A query file with
fewer sequence characters than this stops the run at its start, before any step (a DNA query
translates to at most a third of its length). A query of 20 to
39 aa runs with a warning in the console: in another genome only a nearly identical ortholog of
so short a sequence can be told from a chance match, so calls rest on the neighbourhood. If it
is the mature part of a longer precursor, the precursor is the better query. No module or preset
changes this value.

---

## 3. Easy mode: target genomes

### `max_genomes`
**Default:** `0` (auto)

How many target genomes easy mode fetches. With `0` SynVoy walks up the NCBI taxonomy from the
home species (genus, family, order and class) and takes 3 species per level, at most 20 in all,
each with its best-ranked assembly ([`assembly_ranking`](#assembly_ranking)). A number caps the
total; each level gets at least 2, and places a level cannot fill go to the next one. It is not
used with `--target_species`.

With `--max_genomes` 1 or 2 the console warns. Each call is labelled from the evidence in its own
genome, so the labels mean the same with one target as with many. More targets help in other
ways: the confident calls of closer genomes become queries for the more distant ones, and the
phylogenetic placement check needs at least 4 calls ([§ 14](#phylo_placement_min_calls)).

### `target_species`
**Default:** none

A comma-separated list of species to use as targets instead of the taxonomy walk
(`"Gallus gallus,Mus musculus"`). SynVoy fetches each species' best-ranked assembly; a species
that NCBI does not know or that has no assembly is skipped with a warning. Every listed species
is fetched: `--max_genomes` is not used. For a quick test, name two or three species.

### `assembly_ranking`
**Default:** `hybrid` (`hybrid`, `counts` or `nstats`)

Which assembly of a species is used, and which species fill a taxonomy level's places. Every mode
first prefers the assembly level (chromosome, complete, scaffold, contig), because the search reads
the target's DNA and its contiguity decides the synteny it can find, and then the NCBI status
(reference genome, representative genome, other RefSeq, GenBank). Then:

- `hybrid`: fewer contigs, fewer scaffolds, then a larger contig N50 and scaffold N50;
- `counts`: fewer chromosomes, scaffolds and contigs, then a larger N50;
- `nstats`: a larger N80 and N50, then fewer contigs and scaffolds.

Counts are used only where NCBI's record gives them. No mode looks at the gene annotation. For
the home genome, the species' reference genome is taken when there is one, and otherwise the
status still comes before the level, because the home genome's annotation gives the flanking genes.

### `bad_quality_policy`
**Default:** `auto` (`auto`, `drop`, `keep` or `ask`)

What happens at download (easy mode) to an assembly of low quality: its larger N50 (contig or
scaffold) is below `bad_min_n50`, a contig or scaffold count that NCBI gives is above
`bad_max_contigs` or `bad_max_scaffolds`, or it is a contig-level assembly without any of these
figures. Chromosome-level and complete assemblies always pass.

- `auto`: an assembly the taxonomy walk picked is dropped. An assembly you asked for is kept with a
  warning: a species named in `--target_species`, and the home assembly.
- `drop`: every low-quality assembly is dropped. A low-quality home assembly stops the run.
- `keep`: every one is kept at download. Without `--target_species`, a kept target can still be
  dropped by the check of its FASTA file, unless you also set `--qc_fail_policy keep`.
- `ask`: SynVoy asks on the terminal ([`bad_quality_timeout`](#bad_quality_timeout)).

The quality check runs on each candidate of the taxonomy walk as it is about to fill a place, so
the place of a dropped pick (with `auto`, `drop`, or an answer of no) goes to the next species of
its taxonomy level, and a place that no species of the level can fill carries over to the next
level. The run has fewer targets only when the candidates run out. A species named in
`--target_species` that is dropped is not replaced.

A low-quality home assembly is not marked in the results: its warning appears only in the log of
the download step. Targets are marked by the FASTA check ([`qc_fail_policy`](#qc_fail_policy)).

### `bad_quality_timeout`
**Default:** `300` (seconds)

Only with `bad_quality_policy ask`: how long SynVoy waits for an answer. Without an answer the
answer is no: a target assembly is dropped, and a home assembly stops the run. Without an
interactive terminal SynVoy waits the full time and then answers no.

### `bad_min_n50`
**Default:** `5000` (bp)

The smallest N50 of an assembly that is not low quality. At download it is compared with the
larger of the contig and scaffold N50 in NCBI's record ([`bad_quality_policy`](#bad_quality_policy));
the check of the target FASTA files compares it with the N50 of the file's sequences
([`qc_fail_policy`](#qc_fail_policy)). The default passes almost every assembly. A gene whose
neighbourhood is split over several pieces is hard to place, so raise it to 50,000–100,000 to
keep only contiguous assemblies.

### `bad_max_contigs`
**Default:** `500000`

The largest number of contigs of an assembly that is not low quality. At download it acts only
when NCBI's record gives a contig count, which is rare; the check of the target FASTA files
counts the sequences of each file and compares them with it, at every assembly level.

### `bad_max_scaffolds`
**Default:** `500000`

The same for scaffolds, at download only and only when NCBI's record gives a scaffold count. To
keep only contiguous assemblies, raise `bad_min_n50` instead.

### `qc_fail_policy`
**Default:** `auto` (`auto`, `drop` or `keep`)

What happens to a target genome that fails the check of its FASTA file before the search: an N50
below `bad_min_n50`, or more sequences than `bad_max_contigs`. The check runs in both modes, on
plain `.fna`, `.fa` and `.fasta` files; a gzipped target is not checked.

- `auto`: the genome is kept, and marked, when you chose the targets (pro mode, or easy mode with
  `--target_species`); it is dropped when SynVoy chose them.
- `drop`: every failed genome is removed. If all fail, all are kept.
- `keep`: every genome is searched.

A kept genome that failed is named in a `[WARN]` line of the console and in the run summary,
listed in `synvoy_report.json` under `warnings` and `summary.low_quality_genomes_searched`, its
calls carry `low_quality_assembly: true`, and its row in the figures has a dagger (†). A gene
that is not found in such a genome may be missing from the assembly, not from the species.

---

## 4. Pro mode: your files

File names become genome names inside sequence ids and GFF attributes. Use only letters,
digits, `.`, `_` and `-`, and give each target a name of its own before the extension
(`Apis.fa` and `Apis.fna` collide); otherwise the run stops at the start. A sequence is named by
the first word of its FASTA header, in the home genome and in every target; database forms such
as `lcl|X` or `ref|X|` work. A step stops when two records share a name, when the id of a home
gene used as a flanking gene contains a `|`, or, for a home genome without a GFF, when a sequence
name contains a `|` (SynVoy reads a gene id up to its first `|`).

### `home_genome`
**Default:** none

The home genome as a nucleotide FASTA file, in pro mode (in easy mode the run stops and asks for
`--mode pro`). Give it uncompressed (`.fna`, `.fa` or `.fasta`): the home search builds its BLAST
database from the file as it is. A contiguous assembly gives complete sets of flanking genes. Give
`--home_species` with it ([§ 2](#home_species)).

### `home_gff`
**Default:** none (optional)

The home genome's annotation, as plain (uncompressed) GFF3 of the same assembly. Its first column
must name the FASTA's records by the first word of their headers; a GFF that shares no sequence
name with the FASTA stops the run. With it, SynVoy uses the annotated gene models around the GOI
and their proteins. Each neighbour is tested for being a relative of the GOI before it becomes
an anchor ([`goi_relative_evalue`](#goi_relative_evalue)).

Without it, SynVoy predicts the home genes near the query's hits (Augustus when it is installed,
Prodigal otherwise; [§ 11](#11-gene-prediction-for-a-home-genome-without-a-gff)) and adds the
gene models of the closest annotated target genome. Predicted genes can have wrong boundaries,
miss small exons or merge neighbouring genes, and the report marks the home annotation as
partial.

### `target_genomes`
**Default:** none

The target genomes: a folder (every `.fna`, `.fa` or `.fasta` file in it, plain or gzipped), a
quoted glob (`"genomes/*.fna"`), a comma-separated list (`"a.fna,b.fna"`), or a list in a config
file (`params.target_genomes = ['a.fna', 'b.fna']`). Each file is one species or strain. The
order in which you list them does not matter: SynVoy orders them by their distance to the home
species, by NCBI taxonomy (a local ETE3 database, or NCBI Datasets for files named by a `GCF_` or
`GCA_` accession) and, where that fails, by a MinHash distance to the home genome.
`intermediate/phylo_sort/sort_provenance.json` records the method. Without targets, SynVoy only
analyses the home genome. A listed file, folder or pattern that matches nothing stops the run at
its start, with its name.

### `target_gffs`
**Default:** none (pro mode)

Annotations of the target genomes, given apart from the FASTA files: a folder, a glob or a
comma-separated list of `.gff`, `.gff3` or `.gtf` files, plain or gzipped. Each annotation is
matched to a genome by file stem (`cow.fna` and `cow.gff`), then by assembly accession, then by
an unambiguous name prefix. A subset of the targets may be annotated. The run stops in the
staging step, before any search, when an annotation matches no genome or several, or when a
genome gets a second annotation; a `--target_gffs` that finds no file, or that lists a file,
folder or pattern that matches nothing, stops it at the start.
Easy mode fetches the annotations with the genomes.

Each gene model of the main search then takes `TargetGene`, `TargetProduct` and `TargetID` from
the target's annotated gene or transcript that overlaps it most, or, when none overlaps, from the
nearest one within 5 kb. The per-genome table `<genome>.homology.tsv` has them as `target_gene`
and `target_product`, and the report's records carry `target_gene`. On the models of the two
rescue passes, `TargetGene` holds the query's own name, not a target annotation.

When the home genome has its own GFF, target annotations are read only after a model is built:
they do not move the blocks, the region scores or the label of a call. Without `--home_gff`, the
closest annotated target (at least 10 CDS rows) also supplies gene models for the home genome,
which become anchors and part of the home proteome that the specificity check compares a call
with. With [`goi_family_tokens`](#goi_family_tokens), every GOI row records whether the target's
annotation names the gene's family.

### `allow_unmatched_target_gffs`
**Default:** `false`

`true` turns the stop for an annotation that matches no target genome into a warning, for
example when `--target_gffs` points at a shared folder with annotations of more genomes than
the run uses.

### `home_genetic_code` / `target_genetic_codes` / `default_target_genetic_code`
**Defaults:** `1` / none / `1`

The NCBI genetic code of each genome. `--home_genetic_code` is used for every translation on the
home side: a DNA query, the home search, the home gene model, the home proteome and the flanking
proteins. Each target genome uses its row of `--target_genetic_codes`, a TSV of
`genome<TAB>table` (the genome as its file name, with or without the extension; a header line and
`#` lines are skipped), else `--default_target_genetic_code`: in the genome-wide and window
searches, Smith-Waterman, miniprot, the label's own translations and the two rescue passes.
Supported are 1–6, 9–16 and 21–31, the tables every tool translates with. Examples: the CUG-Ser1
yeasts (*Candida albicans* and relatives) use table 12, *Pachysolen* 26, most ciliates 6.

With table 1 everywhere, no tool gets a genetic-code option. A table other than 1 is written as
`##synvoy-genetic-code N` into each target's result GFF, and the `genetic_code` column of
`target_policy` in the report gives every target's table. Augustus gets the home table; Prodigal
(`--gene_predictor prodigal`, or `auto` without Augustus) runs in metagenome mode, translates
with its own tables 11 or 4 and only warns.

---

## 5. The home locus

The query's hits in the home genome are chained into one hit interval per copy of the gene, and
the intervals are grouped into home loci. Each home locus is searched on its own
([ALGORITHM.md § 5.2–5.5](ALGORITHM.md#52-locating-the-query)).

### `home_locus_mode`
**Default:** `auto` (`auto`, `gene_only` or `alignment`)

Which home loci are searched.

- `auto`: when an annotated home gene matches the query (by name or by position, checked by
  alignment), that gene is the primary locus. Further loci are other places where the query
  aligns at least half as well ([`locus_min_bit_ratio`](#locus_min_bit_ratio)), up to
  [`max_loci`](#max_loci). Without a matching annotated gene the loci come from the alignments.
- `gene_only`: only the annotated gene's locus. Use it when the further loci are relatives you do
  not want. Without a matching annotated gene the loci come from the alignments, with a warning.
- `alignment`: the loci of the alignments alone.

Naming the home gene with [`home_goi_gene`](#home_goi_gene) sets `gene_only`, whatever this
parameter says.

### `home_goi_gene`
**Default:** none

The name of the home gene in the home GFF, when the query's own names are not in it (for
example a query from another species). It is looked up first, and it sets
`home_locus_mode gene_only`: only that gene's locus is searched.

### `gff_search_window`
**Default:** `100000` (bp)

When no annotated home gene carries the query's name, SynVoy looks for the gene by its position:
the annotated genes that overlap a significant hit of the query, or lie within this many bp of
one, are compared with the query, and the best match that passes is the home gene. A gene that
overlaps a hit comes first. Without such a gene, the home gene model is built with miniprot.

### `max_loci`
**Default:** `5`

The largest number of home loci searched: the primary locus and the best-scoring others. A locus
is scored as one copy of the gene, the best set of its hits that covers the query at most once,
so a relative that shares one conserved domain scores below the gene itself. The primary locus
is never removed. Each locus is a search of its own, so a large family costs time.

### `locus_min_bit_ratio`
**Default:** `0.5`

A further home locus is kept only when its score reaches this share of the primary locus's.

### `family_warning_threshold`
**Default:** `10`

More home loci than this, before `locus_min_bit_ratio` and `max_loci`, give a large-family
warning.

### `home_chain_max_intron`
**Default:** `auto` (or a number of bp)

How far apart two hits of the query in the home genome may lie and still be one copy of the
gene. A hit continues a chain when it continues the query (it starts at most 20 residues before
the chain's last hit ends, and ends after it), lies at most this far after that hit, the two are
not in two different annotated genes, and, unless both lie in one annotated gene, no annotated
gene lies entirely between them. `auto` takes the largest intron of the annotated genes the hits
overlap, times [`goi_fallback_intron_margin`](#goi_fallback_intron_margin--goi_fallback_intron_floor)
(5), at least `goi_fallback_intron_floor` (10 kb) and at most 1.1 Mb; without a home annotation,
the floor. A number sets the bound directly.

### `goi_fallback_intron_margin` / `goi_fallback_intron_floor`
**Defaults:** `5.0` / `10000` (bp)

The largest gap between two hits that are chained into one gene: the home gene's largest intron
times the margin, but at least the floor. They set `home_chain_max_intron auto` above and, in a
target genome, the largest gap of a GOI model built from a chain of hits, where the home gene
model's own intron is used; without it, [`max_intron`](#max_intron).

### `locus_gap_genes` / `locus_gap`
**Defaults:** `auto` / `auto`

Two places where the query aligns in the home genome are one locus when at most
`locus_gap_genes` annotated genes lie between them (protein-coding genes when the annotation
states a biotype, every gene otherwise). `auto` is half of
[`n_flanking_genes`](#n_flanking_genes), rounded down: 5 at the default of 10. Two places that
close would share most of their flanking genes.

`locus_gap` gives the gap in base pairs instead, for every scaffold. `auto` counts genes when the
home genome has an annotation, and uses 50,000 bp on a scaffold without annotated genes and for
a home genome without annotation. `-profile legacy_distances` sets 50,000.

---

## 6. Flanking genes

### `n_flanking_genes`
**Default:** `10` (an integer of 1 or more; 5–20 are usual)

How many protein-coding neighbours of the gene, on each side, are taken as anchors in the home
genome. Relatives of the gene are skipped and do not count. More anchors find a neighbourhood
that is partly rearranged and cost more search time; gene-dense, conserved genomes need fewer.

The same number sets:

- the gene's own home neighbourhood in the specificity check: the home gene and `n_flanking_genes`
  genes on each side ([ALGORITHM.md § 7.12](ALGORITHM.md#712-what-is-added-to-the-query-set));
- the home genes the module triggers look at for family members ([§ 1](#modules));
- five counts, each `auto` by default (a number sets it directly):
  `synteny_bridge_max_gap_genes` = 2 × n pitches (20 at n = 10), `strong_synteny_min_flanking` =
  n / 2 rounded up (5), `goi_hull_min_flanking` = round(0.4 × n) (4),
  `strong_synteny_no_goi_min_flanking` = 2 × n (20) and `locus_gap_genes` = n / 2 rounded down
  (5);
- the limit of the flank walk: a side stops before it has n anchors only after it has passed
  more than 5 × n genes and has met n + R protein-coding genes, where R is the number of
  protein-coding genes on that side of the scaffold that overlap a home hit of the query.

The flanking genes are also the unit of the search's distances
([§ 7](#7-distances-in-gene-units)). Measuring how much a neighbourhood grew in a target needs
five anchors in a row there, so with `n_flanking_genes` of 2 or less the ratio of the genome
sizes is always used instead.

### `min_flanking_size`
**Default:** `500` (bp of coding sequence)

The shortest coding length of a flanking protein. It acts only with
`--exon_level_search false` and on genes the flank walk predicts itself (a home genome without a
GFF in a run without target genomes). With the default exon-level search, every gene of the home
annotation is used whatever its length.

### `exon_level_search`
**Default:** `true`

Each flanking gene is searched as its full-length protein and as each of its exons. All these
hits count for the gene, so one conserved exon is enough to place an anchor where the whole
protein has diverged. `false` searches the full-length proteins only (with
[`min_flanking_size`](#min_flanking_size)). The GOI is always searched with its exons.

### `goi_relative_evalue`
**Default:** `1e-5` · module A1 sets `1e-3`

A neighbour is not used as an anchor when it is a relative of the gene: its protein aligns to the
query protein with an E-value at or below this value, computed for the two sequences alone
(Smith-Waterman, BLOSUM62, gap costs 11/1; see `goi_relative_statistic`). Its hits would fall on
the gene's relatives in a target and count the same evidence twice. Unrelated proteins are not
counted as relatives at any length.

The module `local_family` (A1) sets `1e-3`, which also removes distant members of the family, as
`preset_tandem_family` does under `--modules legacy`. `0` switches the test off for the
neighbours. A gene at the GOI's own position that aligns to the query (E ≤ 1e-5, or this value
when it is larger) is never an anchor; another gene there, such as one that holds the GOI in an
intron, can be.

### `goi_relative_statistic`
**Default:** `sw` · module B2 sets `blastp_cbs`

How the test of `goi_relative_evalue` scores a pair: `sw` is the Smith-Waterman E-value;
`blastp_cbs` is blastp with composition-based statistics and low-complexity masking, the test
the specificity check uses for family members. With `blastp_cbs` a low-complexity neighbour no
longer passes as a relative by its composition alone; the module `low_complexity` (B2) sets it.
The gene at the GOI's own position is always tested with `sw`.

### `max_flanking_distance`
**Default:** `0` (bp; 0 = no limit)

A side of the flank walk also stops this far from the gene. With 0, a side ends when it has its
anchors, at the end of the scaffold, or at the limit given under
[`n_flanking_genes`](#n_flanking_genes). A limit such as 500,000 keeps the walk out of a gene
desert next to a large tandem array, at the cost of fewer anchors on that side. A fraction is
rounded down to whole base pairs.

### `expand_goi_similar` / `expand_goi_similar_distance`
**Defaults:** `false` / `300000` (bp)

With `true`, every relative of the gene that the flank walk passes within
`expand_goi_similar_distance` of it becomes a query of its own (`GOI_NEIGHBOR_…`), searched like
the gene. Only relatives between the gene and the last anchor on each side are candidates. In a
gene family their hits draw many candidate regions, so the regions in the figures and the ones
offered to the strong-synteny rescue can be those of relatives instead of the ortholog; this is
why it is off. It can help when the targets are themselves nearby tandem copies of a single-copy
gene. You do not need it for paralogs in the tree: the specificity check already marks calls
that a relative explains better. No module sets either value.

---

## 7. Distances in gene units

Ten neighbouring genes span about 20 kb in a yeast, 150 kb in the honeybee and 1 Mb in human, so
a distance that describes a neighbourhood cannot be one number of base pairs. SynVoy gives these
distances in **gene pitches** and converts them to base pairs for each home locus and each target
genome:

- **pitch (home)**: the span of the anchors of the home locus divided by (number of anchors − 1),
  the mean step from one anchor to the next.
- **expansion** of a target genome: the length of the longest run of anchors that lie in home
  order on one target sequence, divided by the length of the same run at home. One step of the run
  goes at most three anchors ahead and may be 0.2 to 8 times as long as at home; it is measured
  between two hits in the target and between the two query records of those hits at home (two
  exons with the exon-level search, otherwise two genes). The run must hold at least 5 anchors,
  because in a genome where the gene order is lost, pairs of anchors still lie close together by
  chance. Without such a run the expansion is the ratio of the genome sizes (in bases), kept
  between 0.2 and 8, and without that 1.
- **distance in bp** = pitches × pitch (home) × expansion.

| Distance | In pitches | Default | In bp | Fixed value |
|---|---|---|---|---|
| largest gap between two anchors of one block | `block_gap_genes` | 7 | `cluster_distance` | 150,000 |
| extent of the neighbourhood on each side of a block | `window_margin_genes` | 7 | `region_padding` | 150,000 |
| how far outside a block a GOI hit may lie | `goi_margin_genes` | 2 | `goi_margin` | 50,000 |
| largest gap a block may bridge | `synteny_bridge_max_gap_genes` | `auto`: 2 × `n_flanking_genes` = 20 | `synteny_bridge_max_gap` | 6,000,000 |
| margin of the flanking envelope of the dispersed-GOI rescue | `dispersed_envelope_margin_genes` | 9 | `dispersed_envelope_margin` | 200,000 |

Each bp parameter takes `auto` (the default: derive it) or a number, which is then used for every
genome. In `cluster_distance`, `goi_margin` and `dispersed_envelope_margin`, `0` means `auto`; in
`synteny_bridge_max_gap`, `0` switches bridging off; for `region_padding` a number changes the rule
(see its entry). When the home locus gives no pitch (no coordinates), each distance takes its fixed
value. The defaults in pitches match the fixed values at a typical insect locus, with a pitch of
about 23 kb.

**How much DNA is searched.** The window margin is the extent of the neighbourhood: blocks whose
extents overlap are merged, and the anchors inside it count as support for a call. Hits of the
gene of interest are used only when they touch the block or lie within the GOI margin of it, and a
gene model may use DNA up to `gap_search_window` beyond its hits (more for a gene that is long at
home). DNA further out would be searched for nothing, so the DNA cut out and searched is the
smaller of the window margin and the GOI margin plus that room for a model. Where genes lie close
together this is the whole neighbourhood; where they lie far apart it is much less. A window margin
given in bp is searched in full, and so is the whole margin when the home locus has no
coordinates.

**Two more rules follow the home locus.** The block gap is never smaller than 1.25 × the largest
gap between two neighbouring anchors at home, times the expansion, so a neighbourhood that is
unchanged in the target comes out as one block. Hits of one flanking gene belong to one occurrence
of that gene up to 1.5 × the gene's own span at home × the expansion, and never less than
[`max_intron`](#max_intron).

What was used is written to the search log, one line per genome, and to each result GFF as a
`##synvoy-distances` line ([OUTPUT.md](OUTPUT.md)); the region scorer and the hull rescue read it
from there.

`-profile <execution>,legacy_distances` sets these distances, the three sizes of the hull rescue
and the gap between home loci to their fixed values. Use it to see what the derived distances
change, or for a genome where a derived window turns out too large. The limit for hits of one
flanking gene still follows that gene's span at home.

### `cluster_distance` / `block_gap_genes`
**Defaults:** `auto` / `7` (pitches)

The largest gap between two neighbouring flanking-gene hits inside one block. With `auto` it is
`block_gap_genes` pitches times the expansion, and at least 1.25 × the largest gap between
neighbouring anchors at home times the expansion. A number of bp is used for every genome.

### `window_margin_genes`
**Default:** `7` (pitches)

The extent of the neighbourhood on each side of a block (in bp: `region_padding`). Blocks whose
extents overlap are merged, and the anchors inside the extent count as support for a call. The DNA
that is searched is at most this wide, and narrower when the GOI margin plus the room for one gene
model is narrower. The search log prints both for each genome, for example
`window margin 35 kb (DNA searched: 35 kb per side)`.

### `region_padding` / `padding_min` / `padding_max`
**Defaults:** `auto` / `50000` / `200000` (bp)

`auto`: `window_margin_genes` pitches times the expansion. A number switches to a margin from the
spacing of the hits on each scaffold: twice the mean gap between consecutive hits, kept between
`padding_min` and `padding_max`, and searched in full; the number itself is used only on a scaffold
with fewer than two hits. `padding_min` and `padding_max` are read only then, and `padding_max`
must not be smaller than `padding_min`.

### `goi_margin_genes` / `goi_margin`
**Defaults:** `2` (pitches) / `auto` · module A1 sets `goi_margin_genes`

A hit of the gene of interest is used for a gene model only when it lies in the block of flanking
genes or at most this far outside it. A wider margin adds no DNA to the search, but it widens the
GOI core, the place in which a call's E-value is counted
([ALGORITHM.md § 7.13](ALGORITHM.md#713-e-values-in-the-expected-place)).

The module `local_family` (A1) sets `goi_margin_genes` to min(`window_margin_genes`, max(2, R + 1)),
where R is the largest distance, in genes, of a family member next to the home gene, so the copies
of a home array are inside the margin in the target too. It never passes `window_margin_genes`,
where the searched DNA ends. A `goi_margin` in bp that you set wins over the rule.

### `synteny_bridge_max_gap_genes` / `synteny_bridge_max_gap`
**Defaults:** `auto` / `auto`

Two clusters of flanking genes on one sequence stay one block across a gap up to this size, when
the flanking genes continue the home order across it (a local rearrangement the gene may sit in).
`auto` (in pitches) is 2 × [`n_flanking_genes`](#n_flanking_genes), one home neighbourhood: 20 at
the default of 10. `synteny_bridge_max_gap 0` switches bridging off; `synteny_bridge_max_gap_genes 0`
does so only while `synteny_bridge_max_gap` is `auto`. When a gap is bridged:
[`synteny_bridge_min_anchors`](#synteny_bridge_min_anchors) and the entries after it.

### `dispersed_envelope_margin_genes` / `dispersed_envelope_margin`
**Defaults:** `9` (pitches) / `auto`

How far the envelope of the flanking genes on a chromosome is widened for the dispersed-GOI rescue
([§ 13](#enable_dispersed_goi_rescue)).

---

## 8. Genome-wide search and blocks

Each target genome is searched once, with MMseqs2, for the flanking proteins (and their exons) and
for the query. The hits of the flanking genes are grouped into **blocks**: runs of anchors on one
sequence, in home order where possible. The blocks are filtered and capped, and the neighbourhood
of each is searched for the gene ([§ 9](#9-search-inside-a-window)). Genomes are searched from the
closest relative outwards, in waves, and the genes found in one wave join the queries of the next
([ALGORITHM.md § 6–7](ALGORITHM.md)).

### `search_evalue`
**Default:** `0.01` (greater than 0)

The E-value bound of the searches over a whole genome: LOCATE_GENE in the home genome, the flank
walk, and the search of each target genome for the anchors keep hits at E ≤ this value. The search
inside a window derives its own bounds from it (`aug_relaxed_*`, [§ 9](#9-search-inside-a-window)).
A higher value lets weaker anchors through; the blocks and the labels still decide what a hit is
worth. For very close species, 0.001 removes chance anchors.

### `min_hit_identity`
**Default:** `10` (%)

A floor on the identity of a hit of the genome-wide search. It is a guard, not the filter that
decides: hits are kept by E-value, and a hit that reaches E ≤ 0.01 in a whole genome is far above
10 % identity, so lowering it changes nothing. Raising it, for example to 30 among close species,
removes weak anchor hits and the chance blocks they make. Inside a window the floor is
max(`aug_relaxed_identity_min`, `min_hit_identity` × `aug_relaxed_identity_factor`), 15 % at the
defaults. One of the 32 settings that no module changes.

### `min_hit_length`
**Default:** `10` (aligned query residues)

A floor on the aligned length of a hit of the genome-wide search. Like `min_hit_identity` it is a
guard behind the E-value: fewer than 10 residues score about 20 bits at best, and E ≤ 0.01 in a
genome needs 35 bits or more. Inside a window the floor is max(`aug_relaxed_length_min`,
`min_hit_length` / `aug_relaxed_length_div`), 15 residues at the defaults. One of the 32 settings that no
module changes.

Every length floor of the search counts aligned query residues, qend − qstart + 1, the same for
MMseqs2, tblastn and Smith-Waterman. The engines' own alignment-length column differs: MMseqs2's
translated search writes the target span in nucleotides there, about three times the residues,
and tblastn and Smith-Waterman write alignment columns (residues plus gaps).

### `mmseqs_sensitivity`
**Default:** `9.5` · `low_mem` 7.0, `laptop_safe` 8.0

MMseqs2's sensitivity (`-s`) in LOCATE_GENE and the genome-wide search: higher values find more
distant hits and take longer. Values from about 7 (faster, for close species) to 9.5 are usual.
After a search of a genome is killed for memory, it is retried with one thread and a sensitivity
of at most 7.0, then 6.0, at the same split limit; the search of a window also tries a split limit
of 2G. The search log says so, and the report counts such searches (`limits`). Below 1 or above
12, the start and the search log warn.

### `mmseqs_split_memory_limit`
**Default:** `8G` · `low_mem` 2G, `laptop_safe` 4G, `docker_max` 70 % of the memory it gives the
search

MMseqs2's `--split-memory-limit` in LOCATE_GENE and the search: a target database larger than this
is searched in parts, which lowers the peak memory and costs some time. Give a unit (`G` or `M`).
When a search is killed for memory (exit 137), lower this first. `iterative_search_memory auto`
counts it in its request ([§ 17](#iterative_search_memory)).

The limit also has a floor: MMseqs2 needs a part of it whatever the number of parts, and stops
with `Cannot fit databases into ...` when the limit is too small for the genome. With MMseqs2
18-8cc5c, 1G is too small even for a genome of 200 kb, 4G was too small for a genome of 3.2 Gb
(6G worked), and a genome of 3.6 Gb needed 8G. Then raise this limit, and `iterative_search_memory`
with it when that is set to a size (as in `laptop_safe` and `low_mem`); a lower sensitivity or
fewer threads do not help. The search does not retry such a genome: it is recorded as not searched
with the reason `mmseqs_split_memory_limit_too_small`, and the report's `warnings` names these two
settings ([OUTPUT.md](OUTPUT.md#genomes-without-a-call)). A window whose search fails so stops the
search with the same advice, after its 2G rung when the limit is below 2G.

### `mmseqs_max_seqs`
**Default:** `300` (MMseqs2's own default)

MMseqs2's `--max-seqs` in LOCATE_GENE, the genome-wide search and the window searches: the most
targets per query that pass MMseqs2's prefilter. In a translated search a target is one reading
frame stretch, roughly one exon, so a query with many hits in a genome (a gene family, a repeated
domain) can lose hits beyond this number. The report's `limits` block counts the queries with at
least this many result rows (`mmseqs_maxseqs_genome`, `mmseqs_maxseqs_window`), and `warnings`
has a sentence when one is not 0. The count is a lower bound: a query whose prefilter kept this
many targets but fewer of them passed the E-value is not counted. A larger value can find more
hits and takes longer. A whole number of at least 1.

### `mmseqs_verbosity`
**Default:** `1` · `low_mem` and `laptop_safe` 0

Passed to MMseqs2 as `-v` (0 to 3; `mmseqs -h` gives the levels). It changes the messages only.
The searches (LOCATE_GENE, the genome-wide search and the windows) run at 1 or more, so that a
failure prints MMseqs2's reason; at 0 MMseqs2 says only `Prefilter died`.

### `deterministic_goi_search`
**Default:** `true`

Runs the MMseqs2 search of each window with one thread, so a rerun gives the same GOI calls.
MMseqs2 with several threads can order marginal hits differently between runs, and a borderline
call that flips also changes the queries of the later waves. tblastn and Smith-Waterman (parasail)
run in one thread anyway, and a window is small, so it costs little time. `false` lets the window
search use the task's threads.

### `max_blocks_per_genome`
**Default:** `80` (0 = no cap)

The most blocks searched in one target genome for one home locus. When the cap binds, it keeps
first the blocks seeded by the gene itself and the blocks whose GOI zone holds a strong hit of the
gene that the dispersed-GOI rescue left to them ([§ 13](#enable_dispersed_goi_rescue)), then the
blocks with the longer collinear run, then those with more genes; the kept blocks keep their order.
Without the home order (`--disable_synteny_collinearity`) the other blocks are not sorted: they come
in position order. In a genome without any hit of a flanking gene the cap keeps the blocks of the
gene's strongest hits. The log says `Capping blocks` and names each dropped block. With any
`--modules` other than `legacy`, a genome in which the cap binds is first built again from the 5
best places of each flanking gene ([`anchor_max_places`](#anchor_max_places)), and the cap applies
to that. Each block costs window searches and gene models, so the run time grows with the number
of blocks.

### `min_block_genes`
**Default:** `2`

Blocks with fewer distinct flanking genes than this are not searched. A block seeded by the gene
itself is kept whatever its size: by the gene's own home gene among the anchors, by a rescued GOI
hit, or by a GOI hit where no flanking gene has a hit at all. When no block passes, every block is
searched (up to `max_blocks_per_genome`), and these blocks stay when rescued GOI hits add blocks of
their own. `1` is very permissive: one flanking hit on a chromosome is no evidence of synteny.

### `anchor_max_places`
**Default:** `auto` (or a number; 0 = every place)

Each flanking gene seeds blocks only at its N best places in a target genome. A place is a run of
hits of that gene within its own gap limit, and the best places have the largest summed bit score.
The ortholog of a flanking gene is one place; its other places are relatives and shared domains,
and those are what put two anchors near each other by chance, which is how most searched blocks
arise. Hits that are left out still count inside a window (flanking support, flanking models); only
the seeding of blocks changes.

- `auto`: every place, except that with any `--modules` other than `legacy`, a genome in which
  `max_blocks_per_genome` binds is built again from the 5 best places of each flanking gene (the
  genome axis, [§ 1](#modules)); the cap cuts that again if it still binds.
- A number is used for every genome.

The gene's own home gene among the anchors keeps every place, and the dispersed-GOI rescue is
chosen again against the rebuilt blocks. It matters where flanking genes have many places, as in
plant genomes; under `--modules legacy`, use `--anchor_max_places 5` there and whenever the log
says `Capping blocks`. Where a genome has no conserved neighbourhood, its two-anchor blocks are
chance blocks and the limit removes them; the search then falls back to every block, as it does
whenever no block passes `min_block_genes`.

The value used for each genome, and its source (`default`, `genome_axis` or `user_set`), is in the
search log, in the report's `target_policy` block and in the `##synvoy-genome-axis` line of the
genome's result GFF.

### `disable_anchor_resolution`
**Default:** `false`

A target gene is the anchor of at most one home gene: where hits of different flanking genes cover
the same target DNA (overlapping by at least half of the shorter hit, as for two kinases or a
tandem pair), only the best-scoring hit counts as an anchor. `true` counts every hit, so one target
gene can count as two flanking genes.

### `disable_synteny_collinearity`
**Default:** `false`

`true` builds blocks from gene counts alone: no bridging of gaps and no ranking by collinear run.
Without the home order no block has a collinear run, so no block is significant: every GOI call is
then judged in the genome place, with no slot, no gene-order-only call and no positional tie-break
([ALGORITHM.md § 7.3](ALGORITHM.md#73-bridging-a-gap)). For comparisons.

### `synteny_bridge_min_anchors`
**Default:** `2`

How many distinct flanking genes a block must hold before it may bridge a gap
([`synteny_bridge_max_gap_genes`](#synteny_bridge_max_gap_genes--synteny_bridge_max_gap)). With 2, a
block of two genes bridges when the first gene after the gap continues their home order within
`synteny_bridge_max_rank_gap` ranks, so `A B … D E F` is one block; its direction then rests on
those two genes. Under `synteny_bridge_two_sided` the same number is the shortest collinear run the
genes on both sides of the gap must form.

### `synteny_bridge_max_rank_gap`
**Default:** `5`

The largest jump in home rank that still continues the home order across a gap.

### `synteny_bridge_two_sided` / `synteny_bridge_max_per_block`
**Defaults:** `false` / `2`

The default rule looks at the anchored flanking genes on one side of a gap, so the same evidence
bridges or not depending on where the gap falls. `synteny_bridge_two_sided true` weighs both sides
together; it is off because its effect on the calls has not been measured.
`synteny_bridge_max_per_block` caps how many gaps one block may bridge (`0` = no limit).

### `rank_wave_binning`
**Default:** `true`

How the target genomes, sorted from the closest to the farthest, are grouped into waves. `true`:
by rank, the closest tenth one genome at a time, then pairs up to 35 % of the list, threes up to
70 % and fives after that, so the first genes found reach the queries before the distant genomes
are searched. `false`: by phylogenetic distance (closer than 0.05 one at a time, up to 0.15 at most
three of similar distance, beyond that at most five), which gives few, large waves when the targets
are all about equally distant.

### `disable_wavefront`
**Default:** `false`

`true` searches every genome in one wave, so the genes found never join the queries. A control:
run a set both ways to see what the growing query set changes.

### `disable_relative_query`
**Default:** `false`

Inside a window the gene is searched with the **home protein** and also with the gene model of the
genome searched most recently (of its models, the one that aligns best to the home protein). A
model found only with that relative's protein is measured again against the home protein: its
identity, coverage and confidence refer to the home protein, its id carries `_r<n>_` instead of
`_l<n>_`, and `SearchQuery` / `SearchQueryIdentity` name the relative and the identity to it.
`true` searches with the home protein only.

---

## 9. Search inside a window

Each window (a block with its neighbourhood) is searched for the gene with MMseqs2, tblastn and,
by default, Smith-Waterman; miniprot then builds gene models from the hits
([ALGORITHM.md § 7.6–7.7](ALGORITHM.md#76-search-inside-the-window)). A hit of MMseqs2 or tblastn
is kept when it reaches E ≤ 0.1 counted over the window, 15 % identity and 15 aligned query residues (the
`aug_relaxed_*` values below); a Smith-Waterman hit when it reaches E ≤ 1 over the window's six
reading frames. The module `short_query` (B1) lets reading frames and hits count from 10 residues
and runs tblastn with word size 2 and without composition-based statistics.

### `enable_smith_waterman`
**Default:** `true`

Adds a Smith-Waterman search of the query against the six-frame translation of each window. It
finds hits that MMseqs2 and tblastn miss through biased composition or strong divergence. Its time
grows with query length × window length; queries longer than 10,000 aa skip it. Keep it on.

### `sw_method`
**Default:** `auto` (`auto`, `parasail` or `ssearch36`)

The Smith-Waterman backend.

- `auto` and `parasail`: parasail, a SIMD library: each piece of a translated frame gets its best
  alignment (more with `sw_max_alignments_per_piece`).
- `ssearch36`: another search: the best alignment over all six frames, masked, then the next,
  round after round, run in parasail (recorded as `parasail_best_hit_rounds`). Its hits differ
  from those of `auto`.

Every search needs parasail, whatever this setting and `enable_smith_waterman`: the bit score and
the specificity terms of each call, from which its label is computed
([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)), are parasail alignments. Without it the search
stops at its start.

The backend, the parasail version and the SIMD instruction set of a run are in the report's
`run_settings.smith_waterman`.

### `sw_min_score`
**Default:** `20` (whole bits)

A floor on the raw Smith-Waterman score (BLOSUM62, gap costs 11/1) of a window hit. Every
Smith-Waterman hit must also reach E ≤ 1 over the six frames of its window, and that is the
stricter test up to about 30: the lowest score with E ≤ 1 is 31 for a 20-aa query in a 2-kb window,
54 for a 70-aa query in 300 kb and 60 for a 300-aa query in 300 kb. Only a value above those has an
effect. One of the 32 settings that no module changes.

### `sw_min_identity`
**Default:** `10.0` (%)

A floor on the identity of a Smith-Waterman hit, applied with `sw_min_score` and before the E ≤ 1
filter. At the default it does not bind; a higher value, such as 25, makes the Smith-Waterman
search stricter. One of the 32 settings that no module changes.

### `sw_timeout_seconds`
**Default:** `3600` (whole seconds, 1 or more)

A guard, not a limit on the result. When the Smith-Waterman search of one window has not finished
after this many seconds, the genome's search stops, and the task with it; the message names the
genome and block (`<genome>_b<N>`), the window's length, its queries and its cells. A failure of any
of the three window searches stops it the same way. A search killed for memory ends the task with
exit 137 (when that is the only kind of failure among the genomes searched at the same time), and
Nextflow retries it with more memory ([`iterative_search_memory`](#iterative_search_memory)).

The work of a window is fixed by its inputs: every query against the window's six frames,
2 × window bp × query aa cells, which the search log prints for each window with its time
(`Smith-Waterman: N cells in S s, K hits kept.`). Raise the guard on a very slow or busy machine.

### `sw_max_alignments_per_piece`
**Default:** `1` · modules A1 and C1 set `3` when the home gene has two or more copies or two or
more exons

How many Smith-Waterman alignments are returned per piece of a translated frame (a frame is one
piece unless query length × frame length exceeds 200 million cells). Above 1, the residues of each
alignment are masked and the piece is aligned again; further alignments are returned while their
E-value is at most 1, up to this number. Otherwise two copies of the gene, or two exons, in one
frame of a window come out as one: the second is hidden behind the best alignment. Under the
default `--modules auto` most genes with several exons therefore run with 3. It works with
`sw_method` `auto` or `parasail`; `ssearch36` ignores it with a warning.

### The thresholds of the window hits (`aug_relaxed_*`)

MMseqs2 and tblastn search a window with a wide E-value, and of what they return a hit is kept when
it passes three bounds. Startup checks none of these values except that `aug_relaxed_length_div` is
above 0.

| Parameter | Default | What it sets |
|---|---|---|
| `aug_relaxed_parse_evalue_mult` | `10` | A window hit is kept at E ≤ `search_evalue` × this (0.1, counted over the window, not the genome). This is the one that decides which window hits exist: 0.1 in a window of a few hundred kb asks for fewer bits than 0.01 in a genome, the intended relaxation where anchors support the place. |
| `aug_relaxed_evalue_mult` / `aug_relaxed_evalue_cap` | `1000` / `10.0` | The E-value MMseqs2 searches the window with: min(cap, `search_evalue` × mult) = 10. Hits between it and the kept bound are computed and discarded, so neither changes the result while it stays above 0.1. tblastn always reports to E 1000 and is filtered at the kept bound too. |
| `aug_relaxed_identity_min` / `aug_relaxed_identity_factor` | `15.0` / `0.6` | Identity floor: max(min, `min_hit_identity` × factor) = 15 %. The factor matters only when `min_hit_identity` is above 25. |
| `aug_relaxed_length_min` / `aug_relaxed_length_div` | `15` / `2` · module B1 sets `aug_relaxed_length_min` 10 | Length floor: max(min, `min_hit_length` / div) = 15 aligned query residues (qend − qstart + 1, for MMseqs2 and tblastn alike). An exon of fewer residues gives no MMseqs2 or tblastn hit; Smith-Waterman has no such floor beyond 10 residues. |

### `window_tblastn_word_size` / `window_tblastn_comp_stats`
**Defaults:** `3` / `2` (tblastn's own) · module B1 sets `2` / `0`

tblastn's `-word_size` and `-comp_based_stats` in the windows. With values other than tblastn's
defaults, tblastn runs twice per window: the hits come from the run with these values, and
`EeffTblastn`, the label's second E-value, from the run with the defaults. A short query has few
words of 3 and loses score to the composition adjustment, so B1 sets 2 and 0, also for a short
query of low complexity: with the composition statistics on, short queries lost real window hits
in tests. The label's own E-value (blastp, `Eeff`) is not affected.

### `window_mmseqs_min_length`
**Default:** `30` (codons, MMseqs2's own) · module B1 sets `10`

MMseqs2's `--min-length` in the windows: open reading frames shorter than this are not searched.
B1's 10 lets an exon of a short peptide between two stop codons be found. The genome-wide search is
not affected.

### `goi_candidate_regions_per_query`
**Default:** `2` (1 or more) · module A1 sets max(2, home copies + 1), at most 40

GOI loci kept per query in a window. The hits of a query in a window are split into loci (hits
closer than max(5 kb, min(100 kb, 2 × `max_intron`), the gene's intron bound): 40 kb at the
defaults), ranked by the share of the query they cover, then the best bit score, then the number of
hits, and this many are modelled; the rest are dropped. A1's value lets a tandem array with more
than two copies get each copy modelled. The dropped loci are counted (`limits`, `regions_cap`) and
named in the search log, one INFO line per window and query (`kept 2 of 5 GOI loci ...; cut:
...; weakest kept: ...`, at most 10 named; decoys are not logged).

### `max_intron`
**Default:** `20000` (bp)

The floor of the largest gap between two hits of one flanking gene that still count as one
occurrence of that gene: the limit is max(`max_intron`, 1.5 × the gene's span at home × the
expansion), so a gene with long introns is not split into several copies. It also sets the gap at
which the hits of the gene of interest are split into loci (see `goi_candidate_regions_per_query`),
the gap between flanking loci, and the gene's own intron bound when its home gene structure is
unknown (otherwise [`goi_fallback_intron_margin` and `goi_fallback_intron_floor`](#goi_fallback_intron_margin--goi_fallback_intron_floor)).
miniprot is not bound by it: it gets the whole stretch of DNA of a candidate and may place an
intron anywhere in it.

### `gap_search_window`
**Default:** `50000` (bp)

DNA added on each side of the hits of one candidate before miniprot builds its gene model. It is a
floor: for the gene of interest the margin is at least its span at home × the expansion, for a
flanking gene 1.5 × that gene's span at home × the expansion, so a long gene gets room for its
introns. The home gene model uses the same value when it falls back to miniprot. Missing exons are
found by miniprot inside this stretch, or not at all.

### `goi_max_models_per_locus`
**Default:** `10` (1 = the best model only) · module A1 sets max(10, home copies + 1), at most 60

Gene models kept per GOI locus. miniprot reports every place where the query aligns in the DNA;
every model is kept whose coding sequence overlaps that of a better model by less than half. Where
hits of the window search remain on no model, the exons found so far are masked and miniprot runs
again, because with a long allowed intron it does not report a diverged neighbouring copy in the
first run. A further model must reach 30 % of the first run's best score, miniprot's own limit,
unless it lies on a window hit that no kept model explains (counted as `goi_models_kept_below_outs`).
Further models carry `m2`, `m3` … in their id (`…_b0_l1m2_exon_ann`) and are labelled like any
other. When the cap ends the copy rounds of a locus while such a hit, or an alignment that would
have been kept, is left, the report's `limits` block counts it (`goi_model_cap_stopped`) and
`warnings` says that a copy may be missing there.

### `miniprot_short_orf_max_len`
**Default:** `60` (aa) · module B1 sets `100`

Queries shorter than this are aligned by miniprot with `-L 10` (reading frames from 10 aa;
miniprot's default is 30) when the gene models of a window are built. Queries under 30 aa are also
aligned without splicing (`-S`, `-G 2000`). The home gene model is built before the modules and
always uses 60; the rescue passes use `-L 10` for every query.

### `goi_miniprot_timeout` / `flanking_miniprot_timeout`
**Defaults:** `3600` / `3600` (s)

Guards of a miniprot run: a run that has not finished after this many seconds stops the task, with
a message that names the run and the parameter (in the hull rescue the run ends with exit 3).
`goi_miniprot_timeout` covers every run that models the gene of interest: the copy rounds of a
window, the home gene model of ANNOTATE_GOI when no annotated gene matches, and the hull rescue.
`flanking_miniprot_timeout` covers the flanking genes of a window and, for a home genome without a
GFF, the transfer of a target's genes to the home (BORROW_ANNOTATIONS). The search does not start
without miniprot. The strong-synteny rescue has its own guard
([`strong_synteny_rescue_miniprot_timeout`](#strong_synteny_rescue_miniprot_timeout)).

### `tandem_query_restart`
**Default:** `false`

How tandem copies of the gene are found among its hits in a window. Off: hits within 5 kb of each
other are one candidate copy, two or more copies that each cover at least 35 % of the query make a
tandem, and each copy is the span of its hits; two exons of one gene that each cover 35 % of the
query and lie more than 5 kb apart count as two copies, and two copies less than 5 kb apart are one
candidate copy (miniprot's models find such a pair). On: copies are built per strand from hits that
continue the query (a hit starts at most 10 residues before the previous one ends, and ends after
it); a copy counts at a union query coverage of at least 0.35 with a best hit of at least 20
aligned residues or 30 bits (both for a query under 150 aa); copies that share genome are one; and
a tandem needs two counted copies that overlap in the query by more than half. Each copy is written
as one mRNA with one CDS per hit and its union coverage, which the label's completeness bar reads.
Two exons of one gene, or two programs' hits of one exon, then no longer make a tandem. Works with
both [`tandem_copy_parent`](#tandem_copy_parent) settings.

### `goi_first_exon_guard`
**Default:** `false`

A GOI model from miniprot is extended to its own start and stop codon. The window hits cannot see
an exon below about 10 to 15 aa, and miniprot rarely places a first exon of 5 to 16 aa, so a model
often begins at the second exon, and the walk upstream then reaches an ATG in the intron. With this
switch the walk is skipped when the home gene has two or more exons, its first exon is at most
30 aa, and the query residues missing before the model's first CDS are within 5 aa of that exon's
length; the model then carries `FirstExon=missing` and `MissingQueryStart=1-<n>`. The same holds at
the stop end (`LastExon=missing`, `MissingQueryEnd`). Labels do not read these attributes. Applies to
the home query and its decoys only.

### `goi_join_leftover_hits`
**Default:** `false`

A GOI hit that no gene model of its locus covers is written as a `rescued_exon` call, even when it
is an exon that miniprot left out of a model next to it. With this switch such a hit is first
offered to a gene model of its own query that the window's judge labels HIGH or MEDIUM. It must lie
outside the model on its strand; continue the model's query from its first or last CDS (it starts
at most 10 residues before the other ends and ends after it, the residues it skips fit into the
genomic gap, and a gap under 40 bp joins only two pieces of one exon in one reading frame); lie
within the gene's intron bound (see [`max_intron`](#max_intron)); and no CDS of another model or
tandem copy may lie in between. miniprot then models both again (± 1 kb), and the new model
replaces the old one only when it covers the hit and the label's order ranks it strictly first. It
keeps the old model's id and lists the joined hits in `JoinedHits=`. Models found with the
relative's query are not joined.

### `fallback_short_query_len` / `fallback_short_min_aln_aa` / `fallback_short_min_bits`
**Defaults:** `150` (aa) / `15` (aligned query residues) / `30` (bits) · module B1 sets `fallback_short_min_aln_aa` 10

When miniprot builds no model at a GOI locus, the hits are chained into a model of their own (a
hit chain, `fallback_hit_span`, [ALGORITHM.md § 7.7](ALGORITHM.md#77-from-hits-to-gene-models)).
For a query shorter than `fallback_short_query_len`, the chain is accepted on its aligned query
residues (`fallback_short_min_aln_aa`; qend − qstart + 1, summed over its hits) and its best bit
score (`fallback_short_min_bits`), not on query
coverage: the conserved part of a short precursor is a small share of it (the 26-aa mature melittin
is 37 % of its 70-aa precursor), so a perfect hit still covers little. A longer query is refused
with one hit below 25 % query coverage or below 25 aligned residues, and with several hits only
when all three are low (coverage below 25 %, fewer than 35 aligned residues, below 60 bits).
`fallback_short_query_len 0` applies the rule of long queries to every query.

These tests apply only outside a block with E_syn ≤ 0.1. In such a block, which predicted the
gene's place, a chain that fails them is written all the same, with `RefusalLifted=` naming the
tests it failed, and its E-value in its place decides its label: a fragment can be MEDIUM, and
never HIGH while it covers less than 0.7 of the query. Such a call never seeds and is not in the
tree; the search counts the chains as `fallback_lifted`, the report the counted ones as
`summary.refusal_lifted_calls`. The span limit of a chain applies everywhere.

### `gap_min_identity` / `gap_max_hits` / `min_exon_alnlen`
**Defaults:** `15.0` (%) / `5` / `30` (aligned query residues)

Despite their names, these set only the last-resort `raw_hit` row: a locus that got no model and no
tandem copy can give one, from a single hit with identity ≥ max(`gap_min_identity`, 90), an
alignment of at least max(`min_exon_alnlen`, 0.6 × the query length, 15) query residues and a hit
coverage of at least 0.6, at a locus with at least 2 hits, a union coverage of at least 0.75 and a
best hit of at least 60 bits. At most `gap_max_hits` hits are tried, sorted by bit score (then
identity, then aligned query residues), and the first that passes is written. So `gap_min_identity` changes
anything only above 90, and `min_exon_alnlen` only for queries under 50 aa. A `raw_hit` is always
`ModelStatus=fragment`; its label comes from the judge, like any call's.

---

## 10. Gene models and labels

Every GOI call carries a **label** (`Confidence`: HIGH, MEDIUM, AMBIGUOUS or LOW, with the reason
in `Basis`). The label has no parameters: one judge sets it from the E-value of the match in the
place where the gene is expected, its specificity against the home proteome and a completeness bar
of 0.7 query coverage, at fixed significance levels
([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)). The identity thresholds of earlier versions are
gone ([§ 20](#20-removed-parameters)). The parameters below set the **model status**
(`ModelStatus`), what a model reports beside its label, and the checks of the report.

### `classify_fragment_max_qcov` / `classify_complete_min_qcov`
**Defaults:** `0.4` / `0.7` (0 to 1) · module B1 sets `0.3` / `0.5`

How complete a gene model is, written as `ModelStatus`: `fragment` below the first bound,
`complete` from the second when the model has two or more exons or is a tandem copy, otherwise
`partial`. A single-exon model other than a tandem copy is therefore at most `partial`, and a
`raw_hit` or `rescued_exon` row is always `fragment`. Models of the hull and strong-synteny rescues
carry `ModelStatus=rescue` instead. Coverage is the share of the home protein
that the model's alignment covers; with [`core_completeness`](#core_completeness), the larger of
that and the coverage of the protein's core.

The status does not set the label: HIGH needs 0.7 query coverage whatever these bounds are (module
B1's 0.5 included), and a fragment can be MEDIUM when its match in the expected place is
significant. What reads the status: the figures hide fragment models and draw them in the
`*_with_fragments` files ([OUTPUT.md](OUTPUT.md#the-figures)); in the report only a fragment can be
the piece that a wider AMBIGUOUS model of the same gene overrules (`overruled_counted_piece`), and
`summary.counted_fragment_calls` counts the counted calls that are fragments.
The fragment bound must lie below the complete bound; the search warns when it does not.
`-c conf/tp53_strict.config` sets 0.35 and 0.75 (see [`goi_family_tokens`](#goi_family_tokens)).

### `core_completeness`
**Default:** `false` · module B3 turns it on

HIGH needs a call that covers at least 0.7 of the home protein. In a secreted or GPI-anchored
protein the signal peptide and the GPI tail are often the parts a model misses or that are least
conserved, so a model of the whole mature protein can fall below 0.7 of the precursor. With this
switch a gene model also carries `CoreCoverage` (the share of the protein without its signal
peptide and GPI tail that the model's exons cover), and every completeness test of a gene model
reads the larger of `QueryCoverage` and `CoreCoverage`: the label, the model status, the hull
rescue's cost filter and the label of its models, and the report's coverage check. MEDIUM needs no
completeness and is unaffected. A hit chain and a tandem copy have no `CoreCoverage` and are
measured on the whole protein, as are the search's bars on single hits. In conotoxins,
three-finger toxins and phospholipases A2 the signal peptide is the most conserved part of the
precursor.

### `identity_decoupled_min_identity` / `identity_decoupled_max_qcov`
**Defaults:** `50.0` (%) / `0.35`

The report's coverage check. A HIGH or MEDIUM call with at least this identity whose recorded query
coverage lies below the second value is a short local match, not an ortholog: ten residues at
100 % identity say little about a 70-residue protein. The report sets its `goi_class` to
`identity_coverage_decoupled` (with `coverage_demoted`) and leaves it out of the HIGH and MEDIUM
counts; the record stays. A call without a recorded coverage is not demoted. The same test, as an
advisory flag, is listed under `self_consistency` for every call.

### `disable_coverage_demotion`
**Default:** `false`

`true` keeps the coverage check advisory: the calls keep their class and are counted, and
`self_consistency` still lists them.

### `goi_family_tokens`
**Default:** none · one of the 32 settings that no module changes

Family names, separated by commas. With them, every GOI model records whether the target genome's
own annotation of the gene (`TargetGene`, `TargetProduct`) carries one of them:
`GoiFamilyConsistent=true` or `false`, and `GoiFamilyReason` (`matched_gene:<name>`,
`matched_product:<name>`, `no_family_match` or `no_target_annotation`). It changes no label. A name
matches as part of the annotation, ignoring case and every character that is not a letter or a
digit, so `P53` matches `Trp53`. For TP53:
`--goi_family_tokens TP53,TP63,TP73,TRP53,TRP63,TRP73,P53,P63,P73`. `-c conf/tp53_strict.config`
sets these names and also the model-status bounds 0.35 and 0.75 and the region settings of the
`dispersed_family` module (`adaptive_score_floor_frac 0.45`, `adaptive_score_floor_abs 0.08`,
`adaptive_max_regions 4`).

### `report_copies_in_slot`
**Default:** `false` · module A1 turns it on

Counts the copies of the gene in the slot of a block, the stretch between the nearest home
neighbours of the gene ([ALGORITHM.md § 7.13](ALGORITHM.md#713-e-values-in-the-expected-place)).
Every GOI model of a block with a slot gets `CopiesInSlot` (the GOI models in that slot labelled
AMBIGUOUS or better), and each of those copies gets `SlotRank`, counted from the copy nearest the
home neighbour on the left. The report's records carry both, and its `copies_in_slot` block gives,
per genome, the copies at home and in the slot of the genome's best call (`unit`, for example
`2:3`). Advisory: an equal rank is not evidence of positional orthology.

### `report_low_complexity`
**Default:** `false` · module B2 turns it on

Writes `LowComplexity` on the GOI models: the share of the model's protein that lies in
low-complexity windows, the measure the B2 trigger reads (close to SEG). No label reads it.

### `tandem_copy_parent`
**Default:** `seed` (`seed` or `home`)

Whose model a tandem copy is once it seeds the next wave. A tandem copy here is a single-segment
hit at a place where no gene model lies (`EvidenceType=tandem_copy`).

- `seed`: the copy is recorded as `GOI_copy_<n>|<genome>_b<block>_l<locus>_by_<parent>`, where the
  parent is the query it was found with. Once it seeds, `GOI_copy_<n>` is a gene of interest of its
  own in the later waves: the models it finds are judged against the copy, not against the home
  protein, and their identity and coverage (in the GFF, the report and the figures) are to the
  copy. Such a model also carries its identity, coverage and bits against the home protein
  (`HomeIdentity`, `HomeQueryCoverage`, `HomeBits`; the report's `copy_lineage` and `home_*`
  fields), which no label reads. The anchor grid's legend text says how many such models it draws
  with a number, and its source table names the reference of every value (`identity_to`).
- `home`: the copy is recorded under the home gene (`<parent>|<genome>_b<block>_l<locus>_tc<n>`),
  so it can become that gene's relative query, and every model it finds is reported and judged
  against the home protein.

No module changes this switch.

---

## 11. Gene prediction for a home genome without a GFF

Without `--home_gff`, SynVoy predicts the home genes near the query's locus (PREPARE_HOME_PROTEOME),
transfers the gene models of the closest annotated target genome where it can
(BORROW_ANNOTATIONS), and takes the flanking genes from these. The home genome then holds
annotation only around the query's home hits and in the stretch the transferred gene models
cover, from the locus to the outermost of them; the flank walk stays inside the stretch around
the locus and does not step over unannotated DNA to genes predicted around another hit (it goes
on into another hit's prediction window only where that window overlaps or touches the stretch). The same happens when a home GFF
yields no protein, except that its genes are then predicted with Augustus's `fly` model (for
`augustus_species auto`) and without the stop for a eukaryote home without Augustus. With a usable home GFF, only `gene_predictor`, `gene_prediction_timeout` and
`augustus_species` still matter, and only for the optional structural discovery, which predicts the
genes of target windows with Augustus ([§ 15](#15-optional-steps-language-model-structure-discovery)).

### `gene_predictor`
**Default:** `auto` (`auto`, `augustus` or `prodigal`)

The gene finder. `auto` takes Augustus when it is installed, else Prodigal with a warning. With
any `--modules` other than `legacy`, a run whose home species is a eukaryote (in the NCBI taxonomy
on disk) stops before the search instead, because Prodigal cannot predict genes with introns; the
message says how to go on: install Augustus, give the home GFF, or set `--gene_predictor prodigal`
to accept Prodigal's predictions. `augustus` stops the prediction when Augustus is missing, and an
Augustus run that fails stops it with Augustus's exit code and the end of its messages.
`prodigal` is meant for bacteria and archaea: it finds no gene with an intron, so in a eukaryote
it misses most genes. The structural discovery always predicts with Augustus and stops the search
at its start when Augustus is missing.

### `gene_prediction_timeout`
**Default:** `3600` (s, a whole number, at least 1)

Guard of one Augustus or Prodigal run: of the home prediction, of the flank walk, and of the
structural discovery's runs on target windows. A run that has not finished after this many
seconds stops the task.

### `augustus_species`
**Default:** `auto`

The Augustus model. `auto`, with any `--modules` other than `legacy`: the installed model nearest
the home species in the NCBI taxonomy on disk (`~/.etetoolkit/taxa.sqlite`, the database the
species tree reads; it is not downloaded for this). That is the species' own model when there is
one, else the model of the species that shares the deepest taxon with it: *Candida albicans* →
`candida_albicans`, *Apis mellifera* → `honeybee1`, *Homo sapiens* → `human`, an ant →
`camponotus_floridanus`, a firefly → `tribolium2012`. Without the database, for a species it does
not know, and under `--modules legacy`, the model is `fly` (*Drosophila melanogaster*), with a
warning outside `legacy`. The console names the choice (`[INFO] home genome without a GFF, gene
prediction: augustus_species …`), and `effective_params.json` records it with the shared taxon and
its rank.

A name you give is used as given, and a few common names are translated (`honeybee`, `bombus`,
`beetle` and others: `AUGUSTUS_SPECIES_MAP` in `bin/gene_predictor.py`). A name Augustus cannot
load stops the prediction with a message that names the value, the model it maps to and the
nearest installed models. `augustus --species=help` lists the models; the `Genus_species` folders
in Augustus's `config/species` are models too. The Augustus of `synvoy_env` has no silkworm
(`bombyx_mori`) and no mouse model. A model far from the organism places exon and intron borders
wrongly, which splits gene models or merges neighbouring genes.

### `pred_flank_window`
**Default:** `50000` (bp)

Genes are predicted within this many bp on each side of each region of the query's home hits. The
flank walk, when it predicts genes itself, uses the same window around the centre of each region,
and the annotation transfer starts from it and doubles it until the window holds
2 × `n_flanking_genes` predicted genes or spans 500 kb. The doubling that first reaches 500 kb is
kept, so a window that had to grow ends below 1 Mb (about 800 kb at the default 50 kb for a short
locus); with 250 kb or more the first window already spans 500 kb and is kept, the locus plus twice
this value (1.2 Mb at 600 kb). The flank walk takes its anchors only from the stretch the
predictions and the transferred models annotate around the locus. Where genes lie far apart, 50 kb on each side holds fewer genes than `n_flanking_genes`,
and a side of the walk can end at the stretch's edge with fewer anchors (a warning says so); a
larger window gives the flank walk more to choose from and takes longer to predict.

### `borrow_target_flank`
**Default:** `100000` (bp)

The annotation transfer (BORROW_ANNOTATIONS). The home genes predicted near the query's locus are
searched in the closest annotated target. The annotated genes of the hit cluster that holds most of
them, plus this many bp on each side, are aligned to the home genome with miniprot, and a model is
kept at an identity of at least 30 % and a query coverage of at least 50 % when it is the
reciprocal best match. The transferred models are anchors like annotated genes. A larger flank
transfers more neighbours; the flank walk still takes `n_flanking_genes` per side.

### `pred_keep_pct`
**Default:** `0.10` (above 0, at most 1; other values act as 0.10)

When the flank walk predicts the home genes itself, it keeps only the longest predictions: this
share of them, at least 2 × `n_flanking_genes` + 1, and always the genes that overlap the query's
home region. Gene finders predict many short reading frames, Prodigal in a eukaryote above all; a
larger share keeps more of them as possible flanking genes.

### `prodigal_full_genome_fallback`
**Default:** `false`

For the home genome only: predict genes on the whole genome when no home region of the query was
found, instead of predicting nothing. A home GFF that yields no protein turns it on. A whole
eukaryotic genome takes long to predict.

---

## 12. Region scores

After the search, the region scorer groups the hits of each target genome into candidate regions,
scores the neighbourhood of each and keeps the best of them, per genome and home locus
(`regions/locus_<N>/<genome>.scores.tsv`, [OUTPUT.md](OUTPUT.md#region-scores)). The regions are what
the figures draw by default and what the strong-synteny rescue searches. They do not change the
label of any call ([ALGORITHM.md § 8](ALGORITHM.md#8-region-scores)).

### `synteny_weight_base` / `synteny_weight_consistency` / `synteny_weight_strand`
**Defaults:** `0.4` / `0.3` / `0.3`

The weights of the region score. With C the share of the flanking genes found in the region (the
larger of its distinct flanking genes and the flanking support of the search block of a GOI model
it overlaps, divided by the number of flanking genes of the home locus, both sides, normally
2 × `n_flanking_genes`), K the longest run of its hits in home order divided
by their number, and S the strand agreement:

```
quality       = w_base·C + w_consistency·K + w_strand·S
synteny_score = quality · C
```

Coverage enters twice and dominates: order and strand can move the score by at most 0.3 × C each.
Read the score as coverage, with order and strand as modifiers, not as a split of 40, 30 and 30.
The strand term compares the strand of each flanking gene with that of its home gene and takes the
larger of "agree" and "disagree", so a neighbourhood inverted as a whole scores 1 and a mixed one
less.

### `legacy_strand_score`
**Default:** `false`

`true` computes the strand term as older SynVoy versions did: max(plus, minus) / k over the
region's hits, without looking at the home strands. Only to compare region scores with old runs.

### `synteny_goi_overlap_bonus`
**Default:** `0.15` (at most 0.15)

Added to the ranking score of a region that lies within 20 kb of a GOI model of the search that is
not LOW, or of any GOI hit of the genome-wide search, so the region that holds the gene stays ahead of an equally syntenic region without it.
Values above 0.15 act as 0.15. The bonus is prior knowledge, not evidence from the neighbourhood,
so the region's permutation p-value tests `synteny_score` without it; the region table gives both.

### `min_synteny_score`
**Default:** `0.6` (0 to 1)

The label of a region (`confidence` in the region table, `C<label>` in the region names): HIGH from
this score, MEDIUM from half of it, LOW below; a strong-synteny region without a GOI model
(`goi_missing_but_strong_synteny`) and a distant-macrosynteny region where the search modelled the
gene as HIGH or MEDIUM (`goi_dispersed_macrosynteny`) are at least MEDIUM. It decides no search, no
rescue and no call: a region's label rates the neighbourhood, and the report lists the gene-bearing
regions of MEDIUM or better that hold no counted call (`regions.goi_regions_without_call`).

### `max_regions`
**Default:** `0` (the adaptive selection)

How many candidate regions the scorer keeps per genome and home locus. The regions are ordered with
those near such a GOI model or GOI hit first ([`synteny_goi_overlap_bonus`](#synteny_goi_overlap_bonus)),
then by the flanking support of the search block, then by score. With `0` every region in this order
whose score reaches max(`adaptive_score_floor_abs`, the best score of all regions ×
`adaptive_score_floor_frac`), or that holds at least `adaptive_unique_gene_floor` distinct flanking
genes, is kept, at most `adaptive_max_regions`; when none does, the first is kept. A number N keeps the first N. Either way, up to two strong-synteny regions
without a GOI model are added beyond the cap ([`strong_synteny_min_flanking`](#strong_synteny_min_flanking)).

### `adaptive_score_floor_frac` / `adaptive_score_floor_abs` / `adaptive_max_regions` / `adaptive_unique_gene_floor`
**Defaults:** `0.30` / `0.03` / `6` / `3` · module A2 sets `0.45` / `0.08` / `4`

The adaptive selection of `max_regions 0` (above). A higher floor or a lower cap keeps fewer regions
per genome; the `dispersed_family` module and `preset_paralog_discrimination` set the values given.
`adaptive_unique_gene_floor` is one of the 32 settings that no module changes.

### `strong_synteny_min_flanking`
**Default:** `auto` (half of `n_flanking_genes`, rounded up: 5) · `0` switches it off

A region that holds at least this many HIGH flanking models, with no GOI model (not LOW) and no GOI
hit of the genome-wide search within 20 kb of it, is written as `goi_missing_but_strong_synteny`: the neighbourhood is there, but no model of the gene was built.
Up to two such regions are kept beyond the region cap, unless they lie inside a GOI region already
kept, and the strong-synteny rescue searches them ([§ 13](#disable_strong_synteny_rescue)). `0`
switches this region class off, and with it that rescue.

### `strong_synteny_no_goi_min_flanking`
**Default:** `auto` (2 × `n_flanking_genes`: 20)

The report's check `strong_synteny_no_goi` under `self_consistency`: a genome, home locus and
scaffold with at least this many HIGH flanking models and no HIGH GOI call is flagged. Advisory.
The report gives the value a run used (`self_consistency.thresholds.strong_flanking_min`).

### `disable_distant_synteny_rescue`
**Default:** `false`

A gene-bearing region that holds at most one flanking gene scores low. When its GOI hit reaches the
dispersed-GOI rescue's bar ([`dispersed_goi_min_identity` / `dispersed_goi_max_evalue`](#dispersed_goi_min_identity--dispersed_goi_max_evalue--dispersed_goi_min_alnlen--dispersed_min_chrom_anchors),
40 % and 1e-10), a flanking gene lies on the same scaffold, and the genome is visibly rearranged
(rearrangement score of at least 0.5, [ALGORITHM.md § 8](ALGORITHM.md#8-region-scores)), the region
gets a score of at least 0.45. Its class says what the search made of the gene there:
`goi_dispersed_macrosynteny` where it modelled it as HIGH or MEDIUM (the region is then at least
MEDIUM), `goi_hit_no_call` where it did not (no model, an AMBIGUOUS or LOW one, or an answer of the
genome-wide path; the label comes from the score). The rearrangement bar is fixed. A region with
strong flanking support never takes this route. `true` switches it off. It changes regions only,
not calls.

---

## 13. Rescue passes

Three passes look for the gene where no window reached it or no model was built there: the
dispersed-GOI rescue inside the search, and after the search the hull rescue and the strong-synteny
rescue. A fourth, the genome-wide path ([`goi_genome_wide_path`](#goi_genome_wide_path)), gives a
strong hit outside every searched GOI zone an answer of its own, which is never counted. The judge
labels the models of the rescues like those of the search, except the strong-synteny rescue's,
which stay LOW. The rules:
[ALGORITHM.md § 7.11](ALGORITHM.md#711-goi-hits-far-from-their-anchors) and
[§ 9](ALGORITHM.md#9-rescue-passes); the files: [OUTPUT.md](OUTPUT.md#rescue-models).

### `enable_dispersed_goi_rescue`
**Default:** `true`

A strong GOI hit of the genome-wide search that lies outside the GOI zone of every block (block ±
the GOI margin) is in no search. The blocks are built first; with this pass, a hit that reaches the
identity, E-value and length below and lies outside every GOI zone seeds a block of its own, which
is then searched and modelled like any other, when it lies on a sequence that holds the anchors
(the best hit of each flanking gene in the genome) of at least `dispersed_min_chrom_anchors`
distinct flanking genes, inside their envelope: from the first to the last of them, widened on each
side by [`dispersed_envelope_margin_genes`](#dispersed_envelope_margin_genes--dispersed_envelope_margin).
On a sequence with fewer anchors, it seeds a block when it lies within the GOI margin of one of
them: a single neighbour forms no block of its own. A hit inside a GOI zone is left to that block.
The pass does not look at the gene order. Where the neighbourhood is conserved, the gene lies in the
zone of its block and the pass adds nothing. A model found only through the pass carries
`DispersedRescue=1` and the rule that seeded it. When a block whose zone held such a hit is dropped
by the cap, the limit counts give the hit as `dead_zone`. A hit that passes the quality tests but
not the place is listed in the report under `rejected_candidates`, and a hit outside every searched
GOI zone that fails them under `rejected_candidates.dropped`
([`report_nonsyntenic_candidates`](#report_nonsyntenic_candidates)). `false` switches the pass
off.

### `dispersed_goi_min_identity` / `dispersed_goi_max_evalue` / `dispersed_goi_min_alnlen` / `dispersed_min_chrom_anchors`
**Defaults:** `40.0` (%) / `1e-10` / `50` (aligned query residues) / `3` · modules B1 and C1 set `dispersed_goi_min_alnlen`
to 0.6 × the shorter of the query and its longest home exon, between 10 and 50

The tests of the dispersed-GOI rescue (above). A lower identity or a weaker E-value lets relatives
of the gene seed blocks; the sequence and envelope tests keep a paralog elsewhere in the genome
out. A hit of a short query, or of a gene with short exons, cannot reach 50 aligned residues, hence
the module rule. The length is the hit's aligned query residues (qend − qstart + 1); the
`alnlen` that `rejected_candidates` lists is MMseqs2's own column, the target span in nucleotides.
The region step uses the same identity and E-value for its distant-macrosynteny regions
([`disable_distant_synteny_rescue`](#disable_distant_synteny_rescue)).

### `goi_genome_wide_path`
**Default:** `false`

The strong hit without synteny: a GOI hit of the genome-wide search that lies outside every
searched GOI zone (block ± the GOI margin) and outside every model the search wrote lies where no
block predicted the gene. With this switch, every searched genome that has the home proteome, with
or without a block, gets one more step after its windows. The hits are grouped into places (per home
gene and sequence, hits within the larger of the model margin and the gene's intron bound of each
other, never across a GOI zone), and a place is tried when the terms of the judge's HIGH hold for
its hits in the genome place ([ALGORITHM.md § 7.11](ALGORITHM.md#711-goi-hits-far-from-their-anchors)):

- the home protein's blastp E-value against them over the whole genome is at most 1e-6;
- the HSPs with an E-value of at most 1 cover at least 0.7 of the query, their union (with
  [`core_completeness`](#core_completeness), the larger of the query's and its core's coverage);
- no other home protein is a real competitor within 10 bits of it.

Per home gene the places are tried best first (E-value, then bits, then coverage) until one gives
an answer, at most [`max_blocks_per_genome`](#max_blocks_per_genome) windows per genome
(`genome_wide_tries_cap` in the limit counts). A try searches the place ± the model margin, cut at
the GOI zones and at the earlier tries' windows, without flanking genes, and judges its models in the genome place. Its best model
is the answer when the judge labels it HIGH there, no exon of its gene model runs across an exon
junction of the annotated home gene (≥ 2 exons) with 10 residues on each side (`RetrocopyCandidate=1`:
the target would encode both exons without an intron, as a retrocopy does; the model's exons are
tested, not the hits, whose ends can run into an intron), it is not frame-broken without an intron (a frameshifted copy split into
two pieces a few bp apart), and it overlaps no model the search wrote before.

An answer is written as MEDIUM, `GOIClass=genome_wide_strong_hit`, with `GenomeWide=1`,
`GenomeWideCoverage`, `NeighbourMedianIdentity` (the genome's median flanking identity, recorded
for the curator and never a test) and, when its place lies inside the searched window of a block
(outside its GOI zone), `NearBlock` (that block's number). It never joins the query set of the next
wave and is not in the gene tree. The report counts no answer: it keeps one per genome, only in a
genome without a counted call and where no other record lies at its place, and names the kept ones
in the headline as "(not counted)", those near a block apart
([OUTPUT.md](OUTPUT.md#the-headline)). With [`decoy_mode`](#decoy_mode) the reversed and shuffled
home proteins are searched in the genome as well and take the same steps, so the report gives the
decoy rate of the genome place (`decoy_rate.genome_wide`). Every place is listed under
`genome_wide_path`, with what became of it ([OUTPUT.md](OUTPUT.md#the-other-top-level-keys)). An
answer is a candidate for curation: the rule's false-call rate on genes it was not developed on
has not been measured.

### `disable_goi_hull_rescue`
**Default:** `false`

The hull rescue: after the search, miniprot searches the gene across the whole flanking
neighbourhood (the hull) of a target sequence that holds at least `goi_hull_min_flanking` HIGH
flanking models, when no HIGH GOI model already lies in the DNA it would search. It finds a gene
that lies between two blocks, where no window reached. Its models are judged like the search's and
can be counted (`goi_class` `synteny_hull_rescue`). `true` switches it off.

### `goi_hull_min_flanking`
**Default:** `auto` (0.4 × `n_flanking_genes`, rounded half up, at least 1: 4)

The HIGH flanking models a hull must hold for the hull rescue to search it.

### `goi_hull_min_identity` / `goi_hull_min_coverage`
**Defaults:** `40.0` (%) / `0.5` · module B1 sets `goi_hull_min_coverage` `0.3`

Cost filters in front of the judge: a miniprot model of the hull rescue with a lower identity, or
whose coding sequence covers less than this share of the query, is not judged. With
[`core_completeness`](#core_completeness) a model also passes when its alignment covers this share
of the query's core. B1 sets the coverage filter to its own fragment bound, so a hull model of a
short query that B1 would not call a fragment reaches the judge.

### `goi_hull_cluster_max_gap` / `goi_hull_window_pad` / `goi_hull_max_window`
**Defaults:** `auto` (bp)

The sizes of the hull. `auto` takes them from the distances of the search for that genome
([§ 7](#7-distances-in-gene-units)): the HIGH flanking models of a sequence form one cluster up to
the bridge limit, the hull is widened by the GOI margin, and a hull longer than twice the home
neighbourhood (times the expansion) is skipped. A number fixes a size for every genome. For a
result GFF without distances, and under `-profile legacy_distances`, they are 3 Mb, 100 kb and
20 Mb.

### `disable_strong_synteny_rescue`
**Default:** `false`

The strong-synteny rescue: miniprot, with relaxed settings, searches each
`goi_missing_but_strong_synteny` region ([`strong_synteny_min_flanking`](#strong_synteny_min_flanking))
for the gene. Its models are always LOW (`GOIClass=relaxed_rescue`,
`EvidenceType=relaxed_miniprot_rescue`): leads for a closer look, never counted calls. `true`
switches it off.

### `strong_synteny_rescue_miniprot_outc` / `strong_synteny_rescue_window_pad`
**Defaults:** `0.05` / `2000` (bp)

miniprot's `--outc` in the strong-synteny rescue: a model must cover at least this share of the
query (`0.2` keeps only models that cover a fifth of it). The window searched is the region plus
`strong_synteny_rescue_window_pad` on each side.

### `strong_synteny_rescue_miniprot_timeout`
**Default:** `3600` (s)

Guard of the strong-synteny rescue's miniprot run on one region: a run that has not finished after
this many seconds ends the run, and so do a missing miniprot and a third kill of the rescue task for
memory or time. Any other failure lets the run go on: the report names the genome and the step
under `warnings` and `failed_steps`, and a gene that only this rescue would have modelled is
missing there.

---

## 14. Ownership, paralog check and placement

The report checks each call against the home genome: whether it is the gene of interest or one of
its relatives (locus ownership), and whether its divergence fits its species (phylogenetic
placement). [ALGORITHM.md § 10](ALGORITHM.md#10-adjudication-and-the-report).

### `disable_locus_ownership`
**Default:** `false`

Locus ownership. SynVoy builds a panel of home proteins: the annotated genes of each home locus and
the home relatives of the gene (the two entries below). Every recovered GOI protein is aligned to
the panel (Smith-Waterman). A call whose best match is a home relative that beats the gene's own
best entry by at least `locus_ownership_tiebreak_gap` gets the `goi_class` `paralog_not_goi` and is
not counted. A call of the gene that was found from another home locus moves to the locus it
belongs to (`locus_reattributed`). When two entries of different home loci are nearly tied, the
locus whose HIGH flanking genes lie on the call's scaffold wins. `true` skips the step: the calls
keep the locus that found them, and a relative found under the gene's name is counted.

### `locus_ownership_max_genes_per_locus` / `paralog_panel_max_family`
**Defaults:** `5` / `15` · module A1 sets `locus_ownership_max_genes_per_locus` max(5, home copies
+ 1), at most 60; module A2 sets `paralog_panel_max_family` max(15, home loci of the family + 1),
at most 50

The size of the panel. `locus_ownership_max_genes_per_locus`: annotated home genes per home locus
(those that overlap the locus most), so with A1 every copy of a home array is in the panel.
`paralog_panel_max_family`: home proteins added for their similarity to the query (identity of at
least `paralog_panel_min_identity`, the closest first).

### `locus_ownership_pad` / `paralog_panel_min_identity` / `locus_ownership_tiebreak_gap`
**Defaults:** `2000` (bp) / `28.0` (%) / `10.0` (score)

`locus_ownership_pad`: added on each side of a home locus when its genes are collected for the
panel. `paralog_panel_min_identity`: the identity to the query from which a home protein joins the
panel as a relative (with a coverage guard). `locus_ownership_tiebreak_gap`: the Smith-Waterman
score by which a relative must beat the gene to take a call; below it the call stays the gene's
(`ownership_near_tie`).

### `require_paralog_panel`
**Default:** `true`

When the panel cannot be built, the run stops. `false` lets it go on: the failure is recorded, and
the report marks every call `ownership_not_checked`.

### `disable_phylo_placement`
**Default:** `false`

The phylogenetic placement check, advisory. For every GOI model that the search labelled HIGH or
MEDIUM, SynVoy compares the model's divergence from the home protein with its species' taxonomic
distance from the home species, fits the trend over all of them (Theil-Sen), and labels a model
`phylo_discordant` when it is much more diverged than its species predicts (see
`phylo_placement_z_threshold`); the others are `phylo_concordant`. The idea: a paralog that arose
before the speciation is more diverged than an ortholog in the same species. The verdict and its
z-score are on each record whose own model has one (`phylo_verdict`, `phylo_z`), and
`self_consistency.phylo_placement_summary` sums them up. It changes no label. Read it as a hint:
the divergence comes from a local alignment, the species distance is a difference of taxonomic
ranks rather than time, and a gene that evolves fast in one lineage breaks the assumption of equal
rates. In tests on real gene families every `phylo_discordant` model was a real ortholog. `true`
skips the step.

### `phylo_placement_min_calls`
**Default:** `4`

The trend is fitted only when at least this many models have a species distance and the distances
take at least three different values; otherwise every model is `insufficient_data`.

### `phylo_placement_z_threshold`
**Default:** `2.0`

A model is `phylo_discordant` when (residual − median residual) / (1.4826 × the median absolute
deviation) reaches this value. One-sided: only more divergence than expected counts.

---

## 15. Optional steps: language model, structure, discovery

Experimental and off by default. These steps need the packages of `environment-ml.yml` (PyTorch,
transformers, Foldseek; under a conda profile Nextflow builds the environment of the search
(`ITERATIVE_SEARCH`, the pipeline step they run in) from that file when one of these steps is on, and
it holds the packages of `environment.yml` too), model weights that are
downloaded from HuggingFace on first use (ProtT5 about 2.5 GB, ESMFold about 2.8 GB; set `HF_HOME`
to keep them out of your home folder) and, at a useful speed, an NVIDIA GPU. After the search they
score the GOI models of the search: with ProtT5 embeddings (cosine similarity to the home protein,
`EmbeddingSimilarity`) and with ESMFold structures compared by Foldseek (TM-score,
`StructuralSimilarity`). A step you switch on stops the run when its packages are missing, or the
GPU when its device is `cuda`.

The scores change the label only of a GOI row that the judge did not label: LOW becomes MEDIUM from
the medium threshold with at least one flanking gene in its block, and MEDIUM becomes HIGH from the
high threshold with an exon annotation and at least two. Every row of the search carries the
judge's label, so in practice the scores are recorded and change no label. The structural
discovery adds models of its own.

### `enable_plm_search` / `plm_device`
**Defaults:** `false` / `cpu` (`cpu` or `cuda`)

Scores the GOI models with ProtT5 embeddings. `plm_device cuda` runs the model on an NVIDIA GPU;
on a CPU it is slow. With `enable_structural_search` on, both models run on `structural_device`
and `plm_device` is ignored.

### `plm_medium_threshold` / `plm_high_threshold`
**Defaults:** `0.7` / `0.85`

The embedding similarity from which a GOI row that the judge did not label is lifted to MEDIUM or
to HIGH (rules above). They also label the models of the structural discovery: HIGH from the high
threshold, MEDIUM from the medium one, as the TM-score thresholds do.

### `plm_similarity_threshold`
**Default:** `0.5`

Only for the structural discovery: a predicted reading frame of a window becomes a discovery
candidate from this embedding similarity, or from a TM-score of `structural_tm_threshold`.

### `enable_structural_search` / `structural_device`
**Defaults:** `false` / `cpu` (`cpu` or `cuda`)

Scores the GOI models by structure: ESMFold predicts the structure of each model's protein and of
the home protein, and Foldseek compares them. Each protein is folded in overlapping pieces of 380
residues (300 on a GPU with less than 10 GB of memory, 150 below 6 GB). On a CPU this is very slow;
use `structural_device cuda`. When the GPU runs out of memory:
[USAGE.md](USAGE.md#gpu-out-of-memory-in-the-structural-search).

### `structural_medium_threshold` / `structural_high_threshold`
**Defaults:** `0.5` / `0.7`

The TM-score from which a GOI row that the judge did not label is lifted to MEDIUM or to HIGH
(rules above). They also label the models of the structural discovery: HIGH from the high
threshold, MEDIUM from the medium one.

### `structural_tm_threshold`
**Default:** `0.3`

Only for the structural discovery: a predicted reading frame becomes a discovery candidate from
this TM-score.

### `enable_structural_discovery`
**Default:** `false`

With `enable_structural_search` on: Augustus predicts the genes in each window, and a reading frame
whose structure matches the home protein's (TM-score of at least `discovery_min_tm`) and that passes
the guards below becomes a GOI model, labelled from its scores. It looks for a gene that the
sequence search did not find at all. It needs Augustus and stops the search at its start without
it.

| Parameter | Default | The reading frame is kept when |
|---|---|---|
| `discovery_min_tm` | `0.7` | its TM-score against the home protein reaches this value |
| `discovery_min_len_ratio` / `discovery_max_len_ratio` | `0.5` / `2.0` | its length lies between these multiples of the home protein's |
| `discovery_min_goi_coverage` | `0.5` | its Smith-Waterman alignment covers at least this share of the home protein |
| `discovery_min_block_flanking` | `2` | at least this many flanking anchors lie within 100 kb of it |
| `discovery_allow_offblock` | `false` | `false`: in a genome where the sequence search found the gene, it lies on a sequence that holds such a call |


### `gpu_partition` / `gpu_cluster_options`
**Defaults:** none / `--gres=gpu:1`

Under SLURM: the queue and the extra `sbatch` options of the search task when
`enable_plm_search` or `enable_structural_search` is on. Without a partition no GPU is requested.
The example profiles `lrz_ai` and `lrz_ai_container` show a setup for one cluster
([USAGE.md](USAGE.md#5-running-on-hpc--slurm)).

---

## 16. Figures and tree

For each home locus the plot step draws the anchor grid with its print version, the ribbon plot and
the gene-position map, and the tree step builds the gene tree
([OUTPUT.md](OUTPUT.md#the-figures)). These parameters change the drawing only, never a call. A
finished run can be drawn again from its `plot_inputs_*` folders without searching again
([OUTPUT.md](OUTPUT.md#drawing-the-figures-again)).

### `plot_width`
**Default:** `1500` (px)

The width of the ribbon plot (`*_synteny_plot.html`), at least 800. `0` picks a width from the
drawn span, between 1200 and 6000 px. It sets the scale, not which genes are drawn.

### `gap_threshold` / `gap_visual_size`
**Defaults:** `50000` / `20000` (bp)

In the ribbon plot, every gap between genes that is longer than `gap_threshold` is drawn
`gap_visual_size` long, whatever its real length, and labelled with its real length. Lower the
threshold for a compact genome; raise it where long gaps are common and should keep their scale.

### `flank_fallback_bp`
**Default:** `1000000` (bp)

When the candidate regions of a genome hold none of its GOI models, its row shows the genes within
this many bp on each side of the GOI models on the sequence that holds most of them.

### `scale_bar_len`
**Default:** `auto` (or a length in bp)

The scale bar of the ribbon plot. `auto` (or `0`) picks a round length (1, 2 or 5 × 10^k bp) that
is about 120 px long at the drawing's scale; a number fixes the length. The bar refers to the
parts of a track drawn to scale: gaps longer than `gap_threshold` are drawn shortened.

### `hide_goi_absent_tracks`
**Default:** `false`

By default every searched genome is drawn, also one without a GOI call: a searched genome without
the gene is a result. In the anchor grid such a genome keeps its row with the flanking genes that
were placed, its GOI cell shows a dash, and a row with nothing placed reads "no flanking gene
placed"; the legend of the print grid (`*_anchor_grid_print.legend.txt`) names these genomes with
the number of flanking genes placed in each. In the ribbon plot they are drawn below the genomes
with a call and labelled "GOI not placed", and the gene-position map prints "no orthologs mapped"
for an empty row. A missing model is not evidence that the gene is absent. `true` leaves these
genomes out of the figures when at least one genome has a call; the report keeps them either way.

### `plot_block_selection`
**Default:** `false`

What the row of a target genome shows. `false`: the candidate regions that the region scorer kept
([`max_regions`](#max_regions)), narrowed to those with a GOI model and the strong neighbourhoods
without one when those hold most of the row's models ([ALGORITHM.md § 12](ALGORITHM.md#12-figures)).
`true`: the search's own blocks, chosen as the label sees them: the block of every HIGH, MEDIUM and
AMBIGUOUS GOI call, then the significant blocks (E_syn ≤ 0.1, the predicted places) that hold no
counted call, by E_syn, up to 6 blocks in all. A LOW call draws no block of its own, and a block
that is neither significant nor holds a call is not drawn.
Fragment models take part in the choice and stay hidden, so the `*_with_fragments` figures show
the same blocks. A genome with nothing to draw keeps its row, empty. The cap of 6 and the
threshold are fixed in the figure code.

With `true` the block tables are copied into the `plot_inputs_*` folders, so a re-draw follows the
same rule. A run made with `false` has its tables in `blocks/` only; to draw it by its blocks:
`scripts/replot_run.sh <run> --block_tables <run>/blocks/*/*.blocks.tsv`.

### `max_goi_per_genome`
**Default:** `10` · module A1 sets max(10, home copies + 1), at most 60

How many GOI models of a genome steer the choice of the neighbourhood its row shows: the best
ones, by label, then by the strength of the match. Every GOI model inside the chosen view is
drawn, so the cap hides no copy there; it keeps many weak calls elsewhere in a genome from pulling
the view away. A counted call (HIGH or MEDIUM) outside the chosen view, for example in a second
neighbourhood whose models fell below the cap, is listed in the anchor grid (at the locus when it
lies there, else after `//`) and counted in the ribbon plot's label (`+N elsewhere`). With A1, a tandem array larger than ten steers the
view with every copy. The figures read the value from `effective_params.json`, which the plot step
copies into each `plot_inputs_*` folder, so a re-draw uses the run's value. `--max_goi_per_genome`
given in `plot_extra_args` wins over it.

### `plot_extra_args`
**Default:** none

Options added as they are to the call of `bin/plot_synteny.py`, for the figure controls that have
no parameter of their own. Quote the whole value:

```
--plot_extra_args '--print_grid_numbers --no_fragment_variant'
```

| Option | Default | Effect |
|---|---|---|
| `--print_grid_numbers` | off | Write the identity into every flanking arrow of the print grid (`*_anchor_grid_print.svg`). At 183 mm the number is 5 pt. The full grid (`*_anchor_grid.svg`) always has it. A query coverage below 80 % is written in brackets under the arrow, which makes that row about 1 mm taller; when the grid would then exceed `--print_height_mm`, the identity is marked with `*` instead, and the log says so. |
| `--print_width_mm` | `183` | Width of the print grid in mm: 183 for two journal columns, 89 for one. When the grid does not fit, the spacing shrinks first; then the columns with the fewest orthologs are left out and named in the log and in `_print.legend.txt`. Text stays at 5 to 7 pt. |
| `--print_height_mm` | `170` | Largest height of the print grid in mm. Rows are packed tighter first; a grid that still does not fit is written anyway, and the log says so. |
| `--no_print_grid` | off | Do not write the print grid. |
| `--grid_variants` | none | Also write try-out versions of the print grid, separated by commas, or `all`: `compact` (the flanking genes of each side as one grey stack with a count, and the exon model of every GOI copy in the GOI column), `compact_panel` (the same with the exon panel), `clades_muted` (the flanking genes in muted colours, the GOI copies in the clade colours of the gene tree). Files `*_anchor_grid_<variant>_print.svg`. |
| `--grid_max_models` | `10` | Most GOI models drawn for one genome in the anchor grid. A longer array shows the first ten and `+N more`; every copy keeps its small arrow and its row in `_print.source.tsv`. It drops no copy, so no module changes it. |
| `--no_fragment_variant` | off | Skip the `*_with_fragments` figures. GOI models with `ModelStatus=fragment` are hidden in the main ribbon plot and grid and drawn in these files. |
| `--no_orient_to_home` | off | Do not mirror a scaffold whose flanking genes run opposite to the home genome. By default it is mirrored, so the ribbons run straight. |
| `--no_anchor_grid` | off | Do not write the anchor grid. |
| `--clade_count` | `4` | Clades coloured in the tree figure and the clade grid: the midpoint-rooted gene tree is split from the top, the largest clade first, until there are this many; branch lengths are not used. |
| `--species_tree_newick` | none | A species tree in Newick format, with species names as leaves, for the cladogram and the row order of the grid, instead of the tree built from the NCBI taxonomy. Give an absolute path. |
| `--common_names` | `both` | Species labels: `both` (`Scientific (common)`), `common`, `scientific` or `off`. |
| `--common_names_tsv` | none | A table of two tab-separated columns (scientific name, common name) that overrides the common names from NCBI. Absolute path. |
| `--gene_names_tsv` | none | A names table (`gene_id`, `label`, optional `full_name`, `source`) that sets the labels of the home genes it lists. Absolute path. Without it a home gene is labelled with its symbol in the home GFF, else its current symbol at NCBI Gene, else an abbreviation of its product name (marked `*`), else its ID. Each run writes the table it used to `<locus>_gene_names.tsv`. |
| `--no_gene_legend` | off | Do not draw the table of gene names (label, full name, gene ID, source) under the figures. |
| `--no_network` | off | Skip the plot step's NCBI lookups: the common names of the species, and the symbols of home genes the GFF leaves unnamed. |
| `--plot_height` | automatic | Smallest height of the ribbon plot in px. |
| `--goi_zoom` | `1.0` | Magnify the stretch of the ribbon plot that holds the GOI models by this factor and leave the flanks at their scale, so the copies of a dense tandem array stop overlapping (about 8 to 12 for such an array). |
| `--goi_min_px` | `0` | Smallest width of a GOI gene in the ribbon plot, in px (about 30), so its exon structure stays visible at the scale of the locus; its exon lengths then become schematic. |
| `--caption_file` | none | A UTF-8 text file whose lines are written as a caption under the ribbon plot. Absolute path. |

Do not repeat an option that the pipeline already sets (`--plot_width`, `--gap_threshold`,
`--gap_visual_size`, `--flank_fallback_bp`, `--scale_bar_len`); use the parameter instead.
`--max_goi_per_genome` given here wins over the parameter and over a module's value. An option
that `plot_synteny.py` does not know ends the run at the plot step with `unrecognized arguments`:
take it out and run again with `-resume`, which does not repeat the search. Options that were
removed: [§ 20](#20-removed-parameters).

### `skip_tree`
**Default:** `false` · `low_mem` sets `true`

Skips the gene tree (MAFFT and IQ-TREE over the GOI models of each locus) and writes a placeholder
tree; the figures are drawn all the same. The launcher's default profile `auto,low_mem` sets it,
so a run started with `run_synvoy.sh` without `-profile` has no gene tree: add `--skip_tree false`, or use another
memory tier, to get one.

---

## 17. Resources

How much CPU, memory and time each step asks for. None of them changes a result. The memory tiers
and `docker_max` set several of them ([USAGE.md § 2](USAGE.md#memory-tiers)). A step killed for
memory or time is retried twice, a genome download three times; most steps also retry a crash with
exit code 104, 134 or 139, and the region scoring (CLUSTER_REGIONS) is not retried.

### `iterative_search_cpus`
**Default:** `2` · `low_mem` and `laptop_safe` 1, `docker_max` every CPU

CPUs of one search task. The genomes of one wave are searched at the same time by workers, at most
this many (a wave holds 1, 2, 3 or 5 genomes), and MMseqs2 uses the CPUs too. The memory a task
asks for grows with the number of genomes it searches at once (`iterative_search_memory auto`).

### `iterative_search_memory`
**Default:** `auto` · `low_mem` 4 GB, `laptop_safe` 8 GB, `docker_max` the machine's memory − 2 GB

The memory each search task asks for. One task searches every target genome of a home locus, and
the genomes of a wave are searched at the same time, each by a worker that runs its own MMseqs2
(up to `mmseqs_split_memory_limit`) and holds the whole genome in memory. `auto` asks for
max(10 GB, the need of the largest genomes searched at once): 1 GB plus, for each of the W largest
targets, min(`mmseqs_split_memory_limit`, 12 GB per Gb of genome) + the genome + 1 GB, where W is
the most genomes a wave searches at once (at most `iterative_search_cpus`). A size you give is used
as it is. Attempt n of a task asks for n times the request, so a retry after a kill for memory
asks for more; under the local executor at most the memory the executor has (`executor.memory`,
else the machine's), because Nextflow refuses a task that asks for more. The console prints each
request:

```
ITERATIVE_SEARCH memory: 25 GB (attempt 1 x 25 GB; max(10 GB, 25 GB for 1 genome(s) at once, the largest 15 GB))
```

Measured peaks, at a split limit of 16G on a cluster node: one genome of 1.0 Gb at a time,
13.9 GiB (`auto` asks for 15.4 GB); three genomes of 1.4 to 2.3 Gb at a time, 37 to 38 GiB
(58 GB). On a laptop at 2G no test went above 2.7 GiB. When a task is killed for memory (exit code
137), lower `mmseqs_split_memory_limit` first.

### `iterative_search_max_forks`
**Default:** `1`

How many search tasks run at once. Each task searches one home locus with all its target genomes,
so with one home locus a larger value changes nothing. With several loci and enough memory, 2 or
more run them side by side, and their memory adds up.

### `locate_gene_cpus` / `locate_gene_memory`
**Defaults:** `1` / `3 GB` · `low_mem` and `laptop_safe` 2 GB; `docker_max` every CPU and the
machine's memory − 1 GB

CPUs and memory of LOCATE_GENE, which finds the gene in the home genome (tblastn and MMseqs2, once
per run). A large home genome, such as a vertebrate's, needs 8 to 16 GB; a task that ends with exit
code 137 ran out of memory.

### `compute_tree_cpus` / `compute_tree_memory`
**Defaults:** `1` / `4 GB` · `docker_max` the machine's memory − 1 GB

CPUs and memory of the gene tree (MAFFT and IQ-TREE). IQ-TREE runs with a fixed seed, and in one
thread a rerun on the same alignment gives the same tree. With more threads the tree step is
faster (about 1.8 times with two), but a rerun can give other branch lengths and supports: the
tree's info (`trees` in the report) then records `reproducible: false`, and the report warns.

### `default_cpus` / `default_memory` / `default_time`
**Defaults:** `1` / `2 GB` / `4h` · `docker_max` every CPU and the machine's memory − 1 GB

CPUs, memory and time limit of every step without settings of its own. Every executor enforces the
time limit, the local one too, and a retry after a time kill gets n times `default_time` at
attempt n. Steps with their own limit keep it on every attempt: the search and LOCATE_GENE 48 h, the tree 4 h, the home preparation,
the home gene model and the flank walk 3 h, the genome downloads 3 h. To change one of these, give
a config file with `-c`, for example
`process { withName: 'ITERATIVE_SEARCH' { time = '96h' } }`.

### `docker_container`
**Default:** `synvoy-local:latest`

The image every step runs in under the Docker and Singularity profiles; the conda profiles ignore
it. Build the default from the repository's `Dockerfile` ([INSTALL.md](INSTALL.md#docker)) or give
another image. The scripts in `bin/` of your SynVoy folder are used ahead of those in the image, so
a change to them needs no new image.

---

## 18. Output and console

### `outdir`
**Default:** `results`

The folder all results are written to; it is created when missing. What it holds:
[OUTPUT.md](OUTPUT.md). Nextflow's work folder (`work/`) is separate, and `-resume` reuses it
([USAGE.md § 7](USAGE.md#7-resuming-and-caching)).

### `quiet`
**Default:** `true`

Compact console output: the RUN and SKIP lines of each step and the phase banners are left out, so
they do not interrupt Nextflow's progress table. The banner, the `[OK]` result lines, warnings,
errors and the final summary are always printed. `false` prints everything.

### `iterative_quiet_subtools`
**Default:** `true`

Hides the messages printed while the search builds gene models with miniprot. `false` shows them
in the task's log (`.command.log` in the work folder), for debugging.

### `report_nonsyntenic_candidates`
**Default:** `true`

Keeps the GOI hits that pass the identity, E-value and length tests of the dispersed-GOI rescue
([§ 13](#enable_dispersed_goi_rescue)) but fail its tests of place, with the test that refused each
and the nearest flanking gene. The report lists them under `rejected_candidates`, row by row and
merged into places, and a genome with such a place that no counted call of the run covers, and no
call, reads `goi_found_but_not_syntenic` instead of looking like a genome without a hit. They are
not calls. The same table keeps the GOI hits outside every searched GOI zone that failed the
rescue's quality tests (`rejected_candidates.dropped`, the strongest 200 per genome search, with the
test each failed). With `false` no table is written; `target_policy`'s `rescue_dropped_outside_zones`
still counts the dropped hits.

---

## 19. Calibration and comparison switches

These measure SynVoy rather than improve a result. Leave them at their defaults for a normal run.

### `decoy_mode`
**Default:** `false`

Calibrates the labels without a truth set. In every search window the reversed home protein and two
shuffles of it (seeded by the decoy kind and the protein sequence, so a rerun repeats them and a
renamed gene gets the same decoys), and those of the relative's query
where the window has one, are searched like the home protein: modelled, labelled, and given the
same E-values and specificity terms. They compete only with each other, never seed the next wave,
and never reach the result GFFs or the figures; they are written to `decoys/locus_<N>/`. Every
decoy call is false, so the share of decoy cases that reach a label is that label's false-call rate
on your data. The report gives it as `decoy_rate`, with the upper end of the exact (Clopper-Pearson)
95 % interval ([OUTPUT.md](OUTPUT.md#decoy-calls---decoy_mode)). With
[`goi_genome_wide_path`](#goi_genome_wide_path) the same decoys are also searched in the whole
genome, in the same MMseqs2 run as the query set, and take the genome-wide path's steps
(`decoy_rate.genome_wide`). It makes the window search about four times as long. No module sets
it.

### The comparison switches

Each of these turns off a part of SynVoy, or sets it as an older version did, so that its effect
on a set of genomes can be measured. Run the set both ways and compare.

| Switch | Entry | What the run then does |
|---|---|---|
| `--modules legacy`, `--modules none` | [§ 1](#modules) | the preset chooser instead of the modules, or neither |
| `--auto_apply_preset false` | [§ 1](#auto_apply_preset) | the plain defaults |
| `-profile <execution>,legacy_distances` | [§ 7](#7-distances-in-gene-units) | fixed distances in bp instead of gene pitches |
| `--anchor_max_places N` | [§ 8](#anchor_max_places) | a fixed number of places per flanking gene in every genome |
| `--disable_anchor_resolution` | [§ 8](#disable_anchor_resolution) | every hit counts as an anchor |
| `--disable_synteny_collinearity` | [§ 8](#disable_synteny_collinearity) | blocks without home order, so no block is significant |
| `--synteny_bridge_two_sided` | [§ 8](#synteny_bridge_two_sided--synteny_bridge_max_per_block) | gaps bridged on the evidence of both sides |
| `--rank_wave_binning false` | [§ 8](#rank_wave_binning) | waves by phylogenetic distance instead of rank |
| `--disable_wavefront` | [§ 8](#disable_wavefront) | one wave, so the genes found never join the queries |
| `--disable_relative_query` | [§ 8](#disable_relative_query) | windows searched with the home protein only |
| `--tandem_query_restart`, `--goi_first_exon_guard`, `--goi_join_leftover_hits` | [§ 9](#9-search-inside-a-window) | other rules for tandem copies, first exons and leftover hits |
| `--tandem_copy_parent home` | [§ 10](#tandem_copy_parent) | tandem copies judged against the home protein |
| `--disable_coverage_demotion` | [§ 10](#disable_coverage_demotion) | the report's coverage check advisory only |
| `--legacy_strand_score` | [§ 12](#legacy_strand_score) | the strand term of older versions in the region score |
| `--goi_genome_wide_path` | [§ 13](#goi_genome_wide_path) | a strong hit outside every searched GOI zone gets an answer of its own (never counted) |

Write a switch with its value (`--disable_wavefront true`).

---

## 20. Removed parameters

These parameters were removed before version 1.0.0, or after version 1.0.1 where the table says so. What a run does when you still give one:

- **Replaced**: the run stops at its start and names what to use instead.
- **Had no effect**, and **the LLM parameter advisor**: a warning, and the value is ignored.
- **`goi_chain_query_aware`**, **`allow_missing_smith_waterman`**: `true` stops the run at its start, `false` gives a warning.
- **`goi_genome_wide_without_blocks`**: `false` stops the run at its start, `true` gives a warning.

**Replaced**

| Parameter | Use instead |
|---|---|
| `max_flanking_goi_similarity` | `goi_relative_evalue` (an E-value, not a percentage) |
| `enable_name_locus`, `name_locus_min_identity`, `name_locus_pad` | `home_locus_mode` (`auto` already makes the annotated home gene the primary locus; `gene_only` searches that locus alone). `--home_goi_gene SYMBOL` still names the gene |
| `pub_svg` | Nothing: the publication SVG of the ribbon plot is gone. The figure for print is `*_anchor_grid_print.svg`; `*_synteny_plot_view.svg` is the ribbon plot as a picture |
| `enable_matrix_plot` | Nothing: the matrix figure is gone. The anchor grid (`*_anchor_grid.svg`) shows the same and more |
| `seed_on_flanking_support` | Nothing: the genes that join the queries follow one rule, the HIGH calls of a query in a genome, else its MEDIUM calls, each when no other home protein competes for it ([ALGORITHM.md § 7.12](ALGORITHM.md#712-what-is-added-to-the-query-set)) |
| `classify_context_min_collinear` | Nothing: the context promotion is gone. A weak match in the slot of a block with significant gene order is labelled AMBIGUOUS with `Basis=gene_order_only` ("needs curation", never counted), always. Reversed and shuffled home proteins met the promotion's rule in 30 of 181 decoy cases of a test set |
| `classify_high_min_identity`, `classify_medium_min_identity`, `classify_tandem_min_identity`, `classify_tandem_min_qcov` | Nothing: one judge labels every GOI call, tandem copies included, from the E-value in the expected place, specificity against the home proteome and completeness ([ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence)); identity decides no label |
| `classify_high_min_collinear` | Nothing: the order of the anchors enters the label through the block's E_syn (a block is significant at E_syn ≤ 0.1, and only a significant block gives a call its slot or core as the place of its E-value) |

**Had no effect** (accepted, but read by nothing)

| Parameter | What it was meant for |
|---|---|
| `prefer_large_genes` | prefer long genes as anchors |
| `min_gene_identity` | an identity floor for flanking hits (the floor in use is `min_hit_identity`) |
| `aug_dedup_bin_bp` | de-duplication of window hits |
| `max_consecutive_empty_blocks` | an early stop of the search in a genome |
| `gap_min_size`, `gap_evalue`, `gap_min_alnlen` | a gap-filling search, replaced by one miniprot run per candidate |
| `keep_intermediate` | `intermediate/` is always written |
| `max_retries` | retries are set in the `process {}` block of `nextflow.config` |
| `multi_profile`, `multi_profile_max_jobs` | a run with several parameter sets, never built |
| `locus_ownership_synteny_window` | read by nothing |
| `min_exon_query_cov` | a coverage floor for exon hits that no code applied |
| `pub_width_mm`, `pub_palette` | read by nothing (the print figure has a fixed size) |
| `seed_flanking_min_count`, `seed_flanking_min_qcov` | floors of `seed_on_flanking_support`, removed with it |
| `classify_context_min_identity`, `classify_context_min_qcov` | floors of the context promotion, removed with it |
| `classify_fallback_strong_min_identity_floor` | the identity floor of a route by which a hit chain reached MEDIUM; it went with the identity-based label |
| `strict_goi_family` | demoted hit-chain calls whose target annotation lacked a family name. The family-name check is annotation only: `--goi_family_tokens` writes `GoiFamilyConsistent` and `GoiFamilyReason` and changes no label |
| `distance_autotune_close_pct`, `distance_autotune_far_pct`, `distance_autotune_max_relax`, `distance_autotune_min_flanking`, `distance_autotune_high_floor`, `distance_autotune_medium_floor` | settings of the relaxation of the identity thresholds for distant genomes, removed with those thresholds |
| `disable_distance_autotune` | the switch of that relaxation. The judge computes the E-value of each call in the expected place of its own genome |
| `disable_ambiguous_tier` | took AMBIGUOUS out of the label. AMBIGUOUS ("needs curation", never counted) is part of the one label, with `Basis=not_specific` or `gene_order_only`, and cannot be switched off |
| `phylo_placement_promote` | raised an AMBIGUOUS call with a concordant phylogenetic placement to MEDIUM. It could no longer act, because the placement check reads only models labelled HIGH or MEDIUM. The verdicts stay advisory (`phylo_verdict` on each record) |
| `discovery_paralog_margin` (after 1.0.1) | the margin of the structural discovery's guard against a home relative. The guard acts only with a panel of home relatives, and no step gives the search one |
| `structural_max_length` (after 1.0.1) | a length above which ESMFold would skip a protein. The structural search folds every protein in overlapping pieces ([`enable_structural_search`](#enable_structural_search--structural_device)) |
| `disable_paralog_check`, `paralog_confusion_min_gap` (after 1.0.1) | the switch and the flag margin of the reciprocal-best paralog check, which went with them: it was to compare each call with several home paralogs from the query file, but the query is one sequence, so it could flag nothing. Locus ownership compares each call with the home genes ([§ 14](#disable_locus_ownership)) |

**Stops when on, warns when off**

| Parameter | Why it went |
|---|---|
| `goi_chain_query_aware` | a hit chain that follows the query (off by default). Switched on, it refused a continuation of the query in the same reading frame 24 bp on, and so lost true calls in a test. Its link rule lives on in the leftover join (`goi_join_leftover_hits`) |
| `allow_missing_smith_waterman` (after 1.0.1) | let the search start without parasail. Every search needs parasail: the bit score and the specificity terms of each call's label, the hull rescue and the locus-ownership check are parasail alignments, and such a run stopped at the ownership check after the whole search. Install parasail instead |

**Warns when on, stops when off**

| Parameter | Why it went |
|---|---|
| `goi_genome_wide_without_blocks` | removed on 2026-11-16 (after 1.0.1). With `goi_genome_wide_path` it ran the genome-wide path also in a genome without a significant block and marked retrocopy candidates. The path now runs in every searched genome, with or without a block (the genome place does not depend on the blocks), and its retrocopy test is always on. To search no genome-wide places, use `--goi_genome_wide_path false` |

**The LLM parameter advisor** (read only by the advisor, which was removed with them)

| Parameter | What it did |
|---|---|
| `auto_params` | ran the advisor, which suggested values for up to 10 search parameters from the kingdom of the home species and the distance to the targets, with a cloud LLM or built-in rules. On the Nextflow versions SynVoy runs on, none could be applied after launch, so it only printed flags to run again with. The settings that depend on the gene now come from the preset modules ([§ 1](#1-presets-and-preset-modules)), the distances from the run itself ([§ 7](#7-distances-in-gene-units)) |
| `force_disable_advanced_search` | kept the advisor from suggesting the language-model or structural search |
| `llm_provider`, `llm_api_key`, `llm_api_base_url`, `llm_model` | provider, key, endpoint and model of the LLM. The containers no longer receive `LLM_API_KEY`, `GOOGLE_API_KEY` or `OPENAI_API_KEY` |

**Options of the figure script.** `plot_synteny.py` no longer accepts `--target_names`,
`--max_legend_entries`, `--ribbon_alpha_dense`, `--pub_width`, `--pub_palette`, `--pub_svg`,
`--grid_goi_style` and `--grid_goi_max_models`. A run that passes one of them in
`--plot_extra_args` ends at the plot step with `unrecognized arguments`: take the option out and
run again with `-resume`. `--grid_goi_style` chose one of four ways to draw a GOI cell; the anchor
grid now has one GOI layout ([OUTPUT.md](OUTPUT.md#the-goi-column-of-the-anchor-grid)). The
nearest to `--grid_goi_max_models` is `--grid_max_models`. `plot_synteny_matrix.py` no longer
exists.

**Figures.** The anchor-positions figure (`*_anchor_positions.html` and `.svg`; it repeated the
anchor grid and the gene-position map), the publication SVG of the ribbon plot
(`*_synteny_plot.svg`) and the matrix (`*_synteny_matrix.svg`) are no longer drawn. A run draws the
anchor grid with its print version, the ribbon plot, the gene-position map and the tree.
