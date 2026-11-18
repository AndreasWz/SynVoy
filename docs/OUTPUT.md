# SynVoy output

What a run writes to the folder given with `--outdir` (default `results`), which file answers which question, and what each file and field means. The first sections are for reading a result; from [The report](#the-report-synvoy_reportjson) on, this document is a reference.

Some outputs depend on the profile. The launcher's default profile, `auto,low_mem`, builds no gene tree: the tree file is then a placeholder, and the tree page and the clade-coloured grid are not written. With `-profile standard` they are ([USAGE.md](USAGE.md#2-execution-profiles)).

Contents: [The output folder](#the-output-folder) · [Confidence levels](#confidence-levels) · [Reading the result](#reading-the-result) · [The figures](#the-figures) · [The report](#the-report-synvoy_reportjson) · [Gene models](#gene-models) · [Regions and blocks](#regions-and-blocks) · [Decoy calls](#decoy-calls---decoy_mode) · [Limits](#what-the-searchs-limits-dropped-limits) · [Intermediate files and logs](#intermediate-files-and-logs) · [Quick reference](#quick-reference) · [Reproducibility](#reproducibility)

---

## The output folder

```
<outdir>/
├── synvoy_report.json                                     the result: start here
├── synteny_block_locus_<N>_anchor_grid.html, .svg         the main figure
├── synteny_block_locus_<N>_anchor_grid_print.svg          the same at print size, with its legend text
│       (+ _print.legend.txt, _print.source.tsv)           and a table of every cell
├── synteny_block_locus_<N>_anchor_grid_clades_print.svg   copies coloured by gene-tree clade (needs the tree)
├── synteny_block_locus_<N>_synteny_plot.html              ribbon plot (+ a picture, _synteny_plot_view.svg)
├── synteny_block_locus_<N>_gene_positions.html, .svg      gene-position map
├── synteny_block_locus_<N>_*_with_fragments.*             ribbon plot and grid with fragment models drawn
├── synteny_block_locus_<N>_gene_names.tsv                 the gene labels of the figures
├── synteny_block_locus_<N>_species_names.tsv              the common names in the species labels
├── synteny_block_locus_<N>_tree.html                      the tree page (only with a tree)
├── locus_<N>_tree.nwk                                     the gene tree (a placeholder without a tree)
├── locus_<N>_tree.info.json                               how the tree was made (only with a tree)
├── plot_inputs_synteny_block_locus_<N>/                   gene models and tables per genome, and all a re-draw needs
├── regions/locus_<N>/                                     candidate regions per genome, with their scores
├── blocks/locus_<N>/                                      the blocks the search searched, per genome
├── rescue/locus_<N>/                                      gene models of the rescue passes (only when they found one)
├── decoys/locus_<N>/                                      decoy calls (only with --decoy_mode)
├── qc/genome_qc_summary.json                              assembly statistics of each target genome
├── intermediate/                                          outputs of the steps before the search
└── logs/                                                  the log of each task
```

`<N>` numbers the home loci: the places in the home genome where your gene aligns. Most runs have one, `locus_1`, which is the annotated home gene when there is one.

A file per target genome is named after the genome's FASTA file, extension included. In Easy Mode the genomes are `<accession>.fna`, so the gene models of one genome are in `plot_inputs_synteny_block_locus_1/GCF_003254395.2.fna.gff`. In the report, a genome is named by its assembly accession when its file name holds one (`GCF_003254395.2`), otherwise by its file name without the extensions.

Easy Mode also writes `query/` (the protein as fetched), `home_genome/` (the home genome and its annotation) and `downloaded_genomes/easy_mode_genomes/` (the target genomes, with `species_mapping.tsv`, which names the species of each file, and `download_record.tsv`).

A run into a folder that holds an earlier run's results overwrites the files it writes and leaves the others. It warns about this when it starts, and at the end of a successful run that was not resumed, the summary says how many files are older than the run and names up to three of them. Use an empty folder, or `-resume`.

**The home genome in Easy Mode.** When the downloaded home assembly holds chromosomes (an `NC_` or `CM` record, or a record of at least 1 Mb that is not `NW_` or `NZ_`; an organelle genome does not count), only the chromosomes and the organelle genomes are searched. `home_genome/home_genome/home_genome.removed.tsv` lists the other records (id and length; only its header line when nothing was removed), and `home_genome.removed.fna` holds them. When the home gene of your query lies on a removed record and no kept record holds as good a match, the run stops: the message names the record and says to run in Pro Mode with the whole assembly. The report's `home_genome_filter` gives the count and the bp removed. Target genomes are never filtered.

---

## Confidence levels

Every GOI call, that is every gene model of your gene in a target genome, carries one of four labels, and `Basis` says why. Only HIGH and MEDIUM are counted as orthologs.

| Label | `Basis` | Rule | In words |
|---|---|---|---|
| `HIGH` | `significant_specific` | E ≤ 1e-6 in the expected place, specific with a margin of ≥ 10 bits (or no real competitor), and the model covers ≥ 70 % of the query | ortholog, high confidence |
| `MEDIUM` | `significant_specific` | E ≤ 1e-3 in the expected place, and specific | ortholog, medium confidence |
| `AMBIGUOUS` | `not_specific` | E ≤ 1e-3, but another home gene explains the model as well | needs curation: another home gene fits as well |
| `AMBIGUOUS` | `gene_order_only` | 1e-3 < E ≤ 1 in the slot between the gene's neighbours, in a block whose gene order is significant, and specific. Offered only in a genome without a counted call; elsewhere it is left out of the report, the per-genome GFFs and the figures | needs curation: weak match where the gene should be |
| `LOW` | `not_significant` | everything else: a fragment, a weak hit, a chance alignment | no significant match |

The report's headline uses these words, with `gene order only` for a weak match where the gene should be ([The headline](#the-headline)). The figure legends say high confidence, medium confidence, low confidence and needs curation. The web pages call a HIGH call Found, a MEDIUM call Probable and an AMBIGUOUS call Candidate, and do not list LOW models as calls.

**The E-value** (`Eeff`) is that of the home protein against the model, computed in the smallest place SynVoy predicted for the gene: the slot between the nearest found home neighbours (± one gene pitch), else the core of the block, else the searched window (each only in a block with `BlockEsyn` ≤ 0.1, that is a block whose gene order is unlikely by chance), else the whole genome. It is multiplied by the number of such blocks in the genome.

**Specific** means that no other protein of the home proteome matches the model as well: there is no real competitor (`RealCompetitor=0`), or the home protein wins by at least 3.3 bits (`SpecMargin`). Within ±3.3 bits, a model in the slot or core of a significant block is decided for your gene, unless the competitor lies in the gene's own home neighbourhood (`BestOtherInNeighbourhood=1`, as in a tandem array).

With the `signal_anchor` module (a query with a signal peptide or a GPI anchor), the 70 % for HIGH may also be reached on the query's core, without those signals (`CoreCoverage`).

**Not counted either:** a HIGH or MEDIUM record that the report's paralog check gives to another home gene (`goi_class` `paralog_not_goi`), and one with a high identity over a short part of the query (`identity_coverage_decoupled`: identity ≥ 50 % and coverage below 0.35). The rules in full: [ALGORITHM.md § 7.8](ALGORITHM.md#78-confidence).

---

## Reading the result

Open **`synvoy_report.json`**. It holds the adjudicated result: the calls of every home locus merged into one record per target gene, checked against the home genome's paralogs, with the models of the rescue passes included.

### The headline

```bash
jq -r '.summary.headline' synvoy_report.json
```

```
3 ortholog(s) (2 high, 1 medium confidence); 1 need curation (1 gene order only); 4 no significant match; across 5 genome(s).
```

- **ortholog(s)**: the HIGH and MEDIUM records that are counted (`summary.high_confidence_goi`, `summary.medium_confidence_goi`). These are the numbers to quote. Without one, the headline starts with `NO ortholog`.
- **need curation**: the AMBIGUOUS records, and why: `gene order only` (a weak match where the gene should be) or `another home gene fits as well` (`needs_curation_goi`, `needs_curation_by_basis`).
- **of another home gene**: AMBIGUOUS records that the paralog check gave to another home gene. They need no curation for your gene (`ambiguous_other_home_gene_goi`).
- **no significant match**: LOW models (`low_confidence_goi`). They are counted per home locus, and a LOW model of either rescue pass counts twice (both passes write a model as a gene row and an mRNA row), so read this as a rough number of weak leads, not of genes.
- **across N genome(s)**: the genomes that hold at least one model of your gene, LOW included. It is not the number of genomes searched.

A run that built no model of your gene at all says `No GOI ortholog annotations were produced.`

### One row per gene

```bash
jq -r '.goi_dedup.records[] | [.genome, .chrom, .start, .end, .confidence, .goi_class, .identity, .query_cov] | @tsv' \
    synvoy_report.json
```

`goi_dedup.records` lists the HIGH, MEDIUM and AMBIGUOUS calls, one record per target gene, best first. A gene-order-only call in a genome with a counted call is left out (`gene_order_only_dropped_records`), and LOW models are not listed. A record is an ortholog when its `confidence` is HIGH or MEDIUM and its `goi_class` is neither `paralog_not_goi` nor `identity_coverage_decoupled`. Every field is described under [The records](#the-records-goi_deduprecords).

The exons of a call are in the per-genome GFF, `plot_inputs_synteny_block_locus_<N>/<genome>.gff`, under the record's `mrna_id`; those of a `synteny_hull_rescue` record are in `rescue/locus_<N>/` ([Gene models](#gene-models)).

### Genomes without a call

A genome without a call is not proof that the gene is absent there. The summary sorts these genomes by what the search saw:

| `summary` field | Genomes |
|---|---|
| `absent_context_complete` | without a counted call; the best block was searched with its whole window inside the scaffold, so the neighbourhood is there and the gene was not found in it |
| `unknown_context_truncated` | without a counted call; the best block's window ran past a scaffold end, so the place of the gene is partly missing from the assembly: unknown, not absent |
| `absent_no_block` | without a counted call; the search found no block of the gene's neighbours |
| `nearest_paralog_only` | without a counted call, whose best call is a model that another home gene explains as well: only the gene's nearest relative was found. Such a genome is also in one of the three lists above |
| `goi_found_but_not_syntenic` | with a strong hit of your gene that no synteny rule accepted, at a place that no counted call of the run (from any home locus or rescue pass) overlaps. The hits are in `rejected_candidates` |
| `goi_absent_genomes` | whose search built gene models of the neighbours but none of your gene, and found no strong hit elsewhere |
| `genomes_not_searched` | that the search could not search, each with the reason (for example MMseqs2 ran out of memory after its retries). They have no result at all, and `warnings` names them |

```bash
python3 -c "
import json; s = json.load(open('synvoy_report.json'))['summary']
print(s['headline'])
for key in ('absent_context_complete', 'unknown_context_truncated', 'absent_no_block',
            'nearest_paralog_only', 'goi_found_but_not_syntenic', 'genomes_not_searched'):
    print(key, s.get(key))
"
```

Before any number, read `warnings` in the report: it names everything that limits the result, such as a target genome that failed the assembly check and was searched anyway, a step that failed, a genome that could not be downloaded or searched, a home annotation that is only partial, or a UniProt lookup that did not answer. It is empty on a clean run.

---

## The figures

| File | What it shows |
|---|---|
| `*_anchor_grid.html`, `.svg` | **The main figure.** One row per genome (ordered by the species tree from the NCBI taxonomy; without the taxonomy, by the gene tree), one column per home gene, your gene's column in red. An arrow is the gene placed in that genome's neighbourhood: its colour is the home gene, the number its % identity, a number in brackets under it its query coverage when that is below 80 % (`100` over `(7)` = 100 % identity over 7 % of the home protein), its outline the confidence (solid HIGH, dashed MEDIUM, dotted LOW, palest and most broken AMBIGUOUS). A second arrow behind it: more than one copy at the locus. The axis above the grid shows where the home genes lie, and the column on the right gives each row's scaffold and position. Hover for details. An empty cell means that nothing was placed in that neighbourhood, not that the gene is absent. See [The GOI column](#the-goi-column-of-the-anchor-grid). |
| `*_anchor_grid_print.svg` | **The same grid at print size:** 183 mm wide, at most 170 mm tall, text of 7, 6 and 5 pt in Arial, as a journal asks for a two-column figure. It has no title, coordinates or names table: `_print.legend.txt` holds the sentences of a figure legend, and `_print.source.tsv` one row per cell and GOI call (identity, coverage, confidence, scaffold, coordinates, exons). When the grid is too wide, the spacing shrinks first; then the columns with the fewest orthologs are left out, and the log and the legend text say which. Options (through `--plot_extra_args`): `--print_grid_numbers` (the identity in every arrow), `--print_width_mm 89` (one column), `--print_height_mm`, `--no_print_grid`. |
| `*_anchor_grid_clades_print.svg` | **The print grid coloured by the GOI gene tree.** Flanking genes in grey; every copy of your gene is filled with the colour of its clade in the gene tree (`--clade_count`, default 4), and its exon structure is drawn in the panel on the right in the same colour. It shows which copies in different genomes group together, and which are duplications of one lineage or relatives. White marks a copy that is not a leaf of the tree (the tree is built from the search's HIGH and MEDIUM models, before the rescue passes and the report's checks). It needs the gene tree, so the default `low_mem` profile does not write it, and it is written only when it shows something the species tree does not: a genome with two or more copies in the tree, or copies that fall outside the clade of their neighbours in the species tree. With `_clades_print.legend.txt` and `_clades_print.source.tsv`. |
| `*_anchor_grid_compact_print.svg`, `*_compact_panel_print.svg`, `*_clades_muted_print.svg` | **Other layouts of the print grid**, written only on request (`--plot_extra_args '--grid_variants all'`, or a comma-separated choice). `compact`: the flanking genes of each side are one grey stack with how many of them were placed, and the width goes to the copies of your gene. `compact_panel`: the same stacks, with the gene models in a panel on the right. `clades_muted`: the print grid with the flanking genes in muted colours and the copies in gene-tree clade colours; it is written only when the gene tree has three or more leaves. `compact` and `compact_panel` colour the copies by clade when there is such a tree, and are drawn without clade colours otherwise (as under the default `low_mem` profile, which builds no tree). |
| `*_synteny_plot.html` | **Ribbon plot** of one locus, for exploring on screen: one track per genome, every gene at its real size with its exons, ribbons between orthologous genes, your gene's ribbons in red. Scaffolds that read in the opposite direction to the home genome are flipped. Your gene is named on the home track; another track carries its copy count (`×2`, `+1 weak`). `*_synteny_plot_view.svg` is a picture of it, not laid out for print. |
| `*_with_fragments.html`, `.svg` | The ribbon plot and the anchor grid again, **with** the fragment models of your gene (`ModelStatus=fragment`), which the main figures hide. Written only when a run has fragments. |
| `*_gene_positions.html`, `.svg` | **Gene-position map.** Each genome is a line with one dot per gene at its true position: where the gene sits, and whether the neighbourhood is rearranged. |
| `*_tree.html`, `locus_<N>_tree.nwk` | **The gene tree** of your gene's models (MAFFT and IQ-TREE). With the default `low_mem` profile no tree is built: the `.nwk` is a placeholder and there is no tree page. Use `-profile standard`. |

Every figure has its legend under the plot.

### The GOI column of the anchor grid

SynVoy reports GOI *calls*: everything in a genome's neighbourhood that looks like your gene. The figures count a call as a **copy** when it is the best call of its genome or has at least MEDIUM confidence (fragments excepted); every other call is a **weak call**.

- **One small arrow per copy** in the GOI column, in genomic order and pointing in its coding direction, numbered along the chromosome when a genome has several.
- **A panel right of the grid** draws the exon structure of each copy on a line of its own, 5′→3′, under the same number: coding exons to scale, `∧` for an intron, the identity to the query after it (`56 (44)` = 56 % identity over 44 % of the query; for a model found from a tandem copy, `SynVoy_Parent=GOI_copy_<n>`, the identity to that copy). A call built from search hits, not a gene model, is a thin bar.
- **A small pale mark per weak call** at its place (a box for a gene model, a bar for aligned hits); at most eight per cell, then `+N`.
- **`//`** stands before what is not at the locus: another scaffold, or farther away than the neighbourhood is long.
- **When no genome has two copies** (the home genome included), there is no panel: the gene model is drawn in the GOI column itself.
- A genome with more than ten copies shows the first ten models and `+N more` (`--grid_max_models`); every copy keeps its arrow and its row in the source table.
- **`–`** stands in the cell of a genome that was searched and has no GOI call. Its row is drawn all the same, with the flanking genes that were placed; a row in which nothing was placed reads "no flanking gene placed". The legend text names these genomes. `--hide_goi_absent_tracks true` leaves them out.

The same rule counts copies elsewhere: `×N` in the ribbon plot's GOI label and in the gene-position map is the number of copies, `+N weak` the weak calls. The gene model shows coding exons only, so it can have fewer exons than NCBI's exon count, which includes the non-coding exons of all transcripts.

### Gene and species names

Every figure labels a home gene with a gene name, not with its annotation ID. The static figures (anchor grid, gene positions, tree) list the names in a **Gene names** table under the plot: label, full name, gene ID and where the label comes from. In the ribbon plot, hover over a gene for the same. The table is also written to `synteny_block_locus_<N>_gene_names.tsv` and copied to `plot_inputs_synteny_block_locus_<N>/gene_names.tsv`.

A label comes from the first source that names the gene:

| Source (as printed in the table; code in the TSV) | What it is |
|---|---|
| `names table` (`user`) | A table you supply with `--gene_names_tsv` (see below). |
| `home GFF` (`gff_symbol`) | The gene symbol in the home GFF (`gene=` or `Name=`), e.g. `Adh`, `ACT1`. |
| `NCBI Gene` (`ncbi_gene`) | The gene's current symbol at NCBI Gene, looked up by the `GeneID` in the GFF, only for genes the GFF leaves without a symbol. An annotation release keeps the names of its date, while NCBI goes on naming genes, so this recovers many names where a GFF has mostly `LOC…` ids. |
| `product name` (`product_name`) | No symbol anywhere, but a product name: SynVoy abbreviates it. **These labels end in `*`**: they are not official gene symbols. |
| `no name annotated` (`id_only`) | "uncharacterized LOC…" and the like: the gene keeps its ID. |

How a product name is abbreviated:

1. A gene symbol at the end of the name is used as it stands (`E3 ubiquitin-protein ligase MYLIP` → `MYLIP*`, `segmentation protein Runt` → `Runt*`).
2. A short name with one real word is kept whole (`protein lozenge` → `lozenge*`, `formin-2` → `formin-2*`, `annexin B9` → `annexin B9*`): up to 16 characters for a single word, 12 with a designator after it.
3. Otherwise the name is split at spaces, hyphens and slashes. Filler words (protein, of, and, …) are dropped; numbers, single letters and existing abbreviations are kept whole; every other word gives its first letter, and "like" gives `L`. The result is upper case: `pancreatic triacylglycerol lipase` → `PTL*`, `coiled-coil domain-containing protein 170` → `CCDC170*`, `zinc finger protein 629-like` → `ZF629L*`. If that leaves fewer than four characters and the name has only one or two words, the first word gives its first three letters: `arginine kinase` → `ARGK*`.
4. Two genes with the same abbreviation are numbered in genomic order (`lozenge-1*`, `lozenge-2*`). No label is longer than 16 characters, not counting the `*`.

In the ribbon plot, a gene of a *target* genome keeps the name its own annotation gives it, when target GFFs were supplied: its symbol, else the abbreviation of its product name (also marked `*`). Otherwise it carries the label of the home gene it is the ortholog of, and its tooltip names that home gene.

To change a label, edit a copy of `synteny_block_locus_<N>_gene_names.tsv` (columns `gene_id`, `label`, `full_name`, `source`, `ncbi_gene_id`; set `source` to `user` on the rows you change) and draw the figures again with `--gene_names_tsv <file>`. A two-column table (`gene_id`, `label`) is enough; genes it does not list are named as usual. Without network access the NCBI lookup fails and the figures are drawn from the GFF alone; the log and the line under the table say so. `--no_network` skips the lookup on purpose, and `--no_gene_legend` leaves the table out.

A species label carries the species' common name (`--common_names both`, the default) when one is known: from `--common_names_tsv`, from the per-user cache of earlier lookups (`~/.cache/synvoy/common_names.tsv`), or from NCBI during the run. `synteny_block_locus_<N>_species_names.tsv` (also in `plot_inputs_*`) lists every common name with its source: `table`, `cache`, `ncbi`, `none` (NCBI has none), `unreachable` (NCBI did not answer) or `offline` (not looked up). A lookup that NCBI does not answer is asked once more after 2 s, at most twice per figure run; if it still fails, the label is the scientific name alone and the plot log has a warning. A label from the cache can differ between machines.

### Drawing the figures again

Figure options without a pipeline parameter of their own (`--print_grid_numbers`, `--print_width_mm`, `--grid_max_models`, `--grid_variants`, `--gene_names_tsv`, `--no_gene_legend`, `--no_orient_to_home`, `--no_anchor_grid`, `--no_fragment_variant`, …) go through `--plot_extra_args '…'` ([PARAMETERS.md](PARAMETERS.md#plot_extra_args)). To try them without running the search again, draw the figures of a finished run again; the originals are not touched:

```bash
conda activate synvoy_env
scripts/replot.sh results/my_run/plot_inputs_synteny_block_locus_1 results/my_run/replot \
    --home_species "Apis mellifera" --print_grid_numbers
```

Every option of `python3 bin/plot_synteny.py --help` can follow the two folders. A re-draw reuses the run's `gene_names.tsv`, so the labels stay the same; for a run without it, or one drawn without network access, it looks the names up at NCBI Gene. It also reuses the run's `effective_params.json`, so the GOI models that steer each row's view (`max_goi_per_genome`) are the run's; a folder without it is drawn with the default of 10.

To draw every locus of a run again at once, also of a run made with an older version: `scripts/replot_run.sh results/my_run [options]`. The figures go to `results/my_run/figures_redrawn_<date>/`, with a `REDRAWN.txt` that records the run's folder, the version and commit of the figure code, and the home species. The calls and labels are those of the run.

A re-draw follows the rule the run followed: a plot-input folder with block tables (`--plot_block_selection true`) is drawn by the blocks, one without by the candidate regions. To draw a run by its blocks whose tables are only in `blocks/`, pass them: `scripts/replot_run.sh results/my_run --block_tables results/my_run/blocks/*/*.blocks.tsv` (each locus takes its own tables).

---

## The report: synvoy_report.json

### The records (`goi_dedup.records`)

One record per target gene with a HIGH, MEDIUM or AMBIGUOUS call. Calls of one genome and scaffold are one gene when they overlap reciprocally by more than 80 %, when the shorter one lies at least 90 % inside the other, or when they are the two halves of one reading frame. The call that ranks first is kept (by label, then `eeff`, with values below 1e-10 equal, then `block_esyn`, `bits` and the shorter span), except where `superseded_wide_model` or `overruled_counted_piece` below says otherwise.

| Field | Meaning |
|---|---|
| `genome`, `chrom`, `start`, `end` | Where the gene is: the genome, the scaffold and the span (1-based, inclusive). |
| `confidence`, `basis` | The label and why ([Confidence levels](#confidence-levels)). |
| `eeff`, `place`, `block_esyn` | The E-value in the expected place, the place (`slot`, `core`, `window`, `genome`), and the E-value of the block's gene order. |
| `goi_class` | Follows the label: `confident_goi` (HIGH), `probable_goi` (MEDIUM), `syntenic_candidate_unconfirmed` (AMBIGUOUS). A copy from the tandem scan (`EvidenceType=tandem_copy`) is `tandem_goi_copy`, and a model of the hull rescue `synteny_hull_rescue`, whatever their label; the other copies of a tandem array are gene models (`l1m2`, …) and follow the label. The report's checks set `paralog_not_goi` (it belongs to another home gene; not counted), `identity_coverage_decoupled` (high identity over a short part of the query; not counted) and `cross_locus_duplicate` (found from two or more home loci, every merged model with the same identity; counted once). `original_goi_class` keeps the class before the checks. |
| `identity`, `query_cov` | % identity to the home protein, and the share of it that the model covers (`core_cov` too, with the `signal_anchor` module). |
| `coverage_demoted`, `query_coverage` | Only on a record that the coverage rule relabelled `identity_coverage_decoupled`: `true`, and the coverage it was judged on (`query_cov`, or with the `signal_anchor` module the larger of `query_cov` and `core_cov`). |
| `bits` | Bit score of the model's protein against the home protein (Smith-Waterman, BLOSUM62, gaps 11/1). |
| `model_status` | `complete` (≥ 70 % coverage with at least two exons, or a tandem copy), `partial`, `fragment` (below 40 % coverage), or `rescue` (a model of a rescue pass). |
| `refusal_lifted` | Only on a hit chain the fallback gate refused that the search wrote because its block is significant (`RefusalLifted`): the tests it failed. Such a call is labelled and counted like every other; it never seeds and is not in the tree. |
| `mrna_id` | The model's ID in the per-genome GFF, or in `rescue/locus_<N>/`. |
| `target_gene` | The target genome's own gene at that position, when target GFFs were given. For a model of a rescue pass it is the query id. |
| `provenance`, `n_source_loci`, `n_merged_hits` | The home loci that found the call, how many, and how many model rows were merged into the record. |
| `cross_locus_duplicate` | `true` when the record was found from two or more home loci and every merged model has the same identity; its `goi_class` is then `cross_locus_duplicate`. |
| `owning_gene`, `owning_locus`, `owning_is_goi_gene`, `ownership_bit`, `ownership_gap` | The paralog check (locus ownership): the protein of the home paralog panel that the call matches best, the home locus it belongs to, whether it counts as your gene, its raw Smith-Waterman score (BLOSUM62, gaps 11/1) and its lead over the next entry. When no home gene is your gene's own (for example a home genome without a GFF), your gene's home protein stands for it, named `GOI_<query id>`; `owning_locus` is then the primary home locus. A home copy near-identical to your gene counts as your gene. |
| `ownership_near_tie` | `true` (absent otherwise) when a paralog is the best panel entry but beats your gene's best entry by less than `locus_ownership_tiebreak_gap` (10): the call stays your gene's, and `owning_gene` and `owning_locus` name your gene's entry. |
| `ownership_not_checked` | `true` (absent otherwise) when the paralog check of this call's genome and home locus failed, or the home paralog panel did (`failed_steps`): the call may belong to another home gene. |
| `locus_reattributed`, `inferred_paralog`, `ownership_tiebroken` | The call belongs to another home locus than the one that found it; the paralog it was given to; whether the gene order of the neighbourhood decided a near tie between two loci. |
| `phylo_verdict`, `phylo_z` | Advisory: whether the call diverged more than its species should have (`phylo_discordant`, `phylo_concordant`, `insufficient_data`), from `intermediate/phylo_placement/locus_<N>.phylo_placement.tsv`. |
| `superseded_wide_model` | `true` when this model replaced a wider model of the same gene, at least twice as long but covering less than 0.15 more of the query. |
| `same_orf_merged` | `true` (absent otherwise) when two gene models that are the two halves of one reading frame were merged into this record (miniprot can model the two halves of a gene apart where the query aligns twice inside it). |
| `overruled_counted_piece` | `true` (absent otherwise) when this AMBIGUOUS wider model (`basis` `not_specific`) kept the record of a HIGH or MEDIUM fragment inside it: another home protein explains the gene as well, so its piece is no ortholog call. The piece keeps its own row in the per-genome GFF, and in a gene tree it is drawn hollow, "not counted". |
| `own_copy_set_aside`, `own_copy_bits` | The model's `OwnCopySetAside` and `OwnCopyBits` (absent when it has none): see [the GFF attributes](#the-per-genome-gff). |
| `low_quality_assembly` | `true` when the genome failed the assembly check and was searched anyway. |
| `copies_in_slot`, `slot_rank` | Only with `--report_copies_in_slot`: the model's `CopiesInSlot` and `SlotRank`. |
| `dispersed_rescue`, `rescue_rule`, `rescue_anchors`, `rescue_anchor_bp` | Only on a model found because the dispersed-GOI rescue made the search look there (`DispersedRescue=1`): `true`, the rule that seeded its block (`envelope`, `near_neighbour`, or `no_flanking_hits` in a genome without any hit of a flanking gene), the anchors on the hit's sequence and the distance to the nearest ([ALGORITHM.md § 7.11](ALGORITHM.md#711-goi-hits-far-from-their-anchors)). |

### `summary`

The fields of the headline and of [Genomes without a call](#genomes-without-a-call) are described above. The other fields:

| Field | Meaning |
|---|---|
| `own_copy_set_aside_goi` | Records whose model carries `OwnCopySetAside`: calls whose specificity left a near-identical home copy of your gene out of the competitors. 0 when the home locus holds no such copy. |
| `gene_order_only_dropped` | Gene-order-only calls left out because their genome has a counted call; `goi_dedup.gene_order_only_dropped_records` lists them. The per-genome GFFs in `plot_inputs_*` and the figures leave them out too. |
| `headline_metric` | The same number as `high_confidence_goi`. |
| `coverage_demoted_calls` | HIGH and MEDIUM calls relabelled `identity_coverage_decoupled`. |
| `coverage_not_recorded_calls` | HIGH and MEDIUM calls of at least 50 % identity that the coverage rule could not test, because no coverage was recorded; they stay counted. |
| `goi_hits_collapsed_by_dedup` | HIGH, MEDIUM and AMBIGUOUS model rows merged into another record. |
| `nonsyntenic_candidate_count`, `nearest_paralog_only_genomes` | How many `rejected_candidates` rows there are, and how many genomes are in `nearest_paralog_only`. |
| `nonsyntenic_place_count`, `nonsyntenic_unmodelled_place_count` | The places of `rejected_candidates`, and those that no counted call of the run overlaps (the places behind `goi_found_but_not_syntenic`). |
| `dispersed_rescue_goi`, `dispersed_rescue_by_rule` | Counted calls (HIGH or MEDIUM, not demoted) found only because the dispersed-GOI rescue made the search look there (`DispersedRescue=1`), and how many per rule. |
| `counted_fragment_calls`, `refusal_lifted_calls` | Counted calls (HIGH or MEDIUM, not demoted) whose model is a fragment (`model_status` `fragment`), and those that are hit chains written in a significant block although the fallback gate refused them (`refusal_lifted`). Both are in the headline counts. |
| `genomes_without_flanking_hits` | Genomes where a search found no hit of any flanking gene and still searched blocks: the gene's own hits seeded them, with no neighbour test (`RescueRule=no_flanking_hits`), and `warnings` names them. |
| `failed_qc_genomes_with_downstream_results`, `staging_empty` | Genomes that failed the assembly check and still have results; `true` when the report found none of its staged inputs (`staging_diagnostics`). |
| `paralog_misassignments`, `locus_reattributed_orthologs`, `cross_locus_duplicate_goi`, `goi_owner_search_fumbles` | What the paralog check and the merging of home loci changed. |
| `self_consistency_flag_count`, `phylo_discordant_calls`, `paralog_confusion_flags`, `strong_synteny_no_goi_flags` | Advisory flags (`self_consistency`). |
| `gene_family_query` | True when the query looks like a member of a large family (`gene_family_advisory`). |
| `low_quality_genomes_searched`, `ownership_not_checked_calls` | Genomes searched although they failed the assembly check; calls whose paralog check failed. |
| `*_pre_dedup`, `*_pre_ownership` | The HIGH and MEDIUM counts before the merging of records and before the paralog check. |
| `ambiguous_goi` | The same number as `needs_curation_goi`, under its old name. |
| `total_annotations`, `total_goi_annotations`, `confident_goi_annotations`, `resolved_goi_annotations`, `fallback_goi_annotations`, `genomes_with_annotations`, `genomes_with_goi_annotations`, `total_new_genes`, `low_confidence_regions` | The counts behind `annotations`, `synteny_results` and `regions`: GFF rows (as in `annotations`), model proteins and LOW regions, not genes. `genomes_with_goi_annotations` is the N of the headline's "across N genome(s)". |
| `ambiguous_goi_annotations`, `goi_ambiguous_only_genomes` | Despite their names, these count the class `ambiguous_goi_family_member`, which is the class of LOW models. |
| `total_raw_search_hits`, `genomes_with_raw_search_hits` | Always 0: the raw search hits do not reach the report. Not a result. |

### The other top-level keys

| Key | What it holds |
|---|---|
| `warnings` | Sentences to read before any number (see [Genomes without a call](#genomes-without-a-call)). |
| `goi_dedup` | `records` (above), the counts of the merging, and `gene_order_only_dropped_records`. |
| `self_consistency` | `flags`, things the run noticed about its own result: `identity_coverage_decoupled`, `paralog_misassignment`, `paralog_confusion` (a call matches another paralog better than most calls of its locus do), `locus_reattributed`, `goi_owner_search_fumble` (a call that belongs to your gene's own home locus was found only from another one), `cross_locus_duplicate`, `strong_synteny_no_goi` (a strong neighbourhood without a HIGH call of your gene) and `phylo_discordant`. `thresholds` gives the values the checks used. `locus_ownership_summary` counts what the paralog check did, and its `goi_entry` names the panel entry that stands for your gene and how it was found (`found_by`: `goi_info`, the home gene the home step took; `identity_rule`, a home protein ≥ 95 % identical over ≥ 90 % of your gene's home protein; `query`, your gene's home protein itself; null when no panel was built). |
| `home_goi` | Only when the home step took an annotated home gene (`gene_id` in `goi_info.json`): the home gene your query was matched to, with `query_id`, `gene_id`, `transcript_id`, `method` and `found_by`. |
| `rejected_candidates` | Strong GOI hits that no synteny rule accepted (`candidates`, one row per hit; from `hits/<locus>_<genome>.nonsyntenic.tsv` of the search), with the rule (`rejected_by`), the nearest flanking gene on the hit's sequence (`nearest_anchor_bp`; null, with `nearest_anchor_same_chrom` false, when that sequence carries none), the hit's `strand`, `qstart`–`qend` (the query residues), `bits`, `alnlen` (the search engine's alignment length: for MMseqs2, the hit's target span in nucleotides; the rescue's length test reads `qend` − `qstart` + 1), `locus` (the home locus whose search refused it), `envelope` (`inside` or `outside` the envelope of its sequence, `chromosome_without`, `none_in_genome`), `anchors_on_sequence` and `genes_on_sequence` (the flanking genes whose best hit lies on that sequence, which the rule counts, and those with any hit there), and `modelled_here` / `modelled_by` (the best label and id of a model of the same search over the hit). Coordinates are 0-based, end exclusive. `places` merges the rows of one home locus, genome and sequence that lie within 5 kb of each other, with `n_rows`, `queries`, the best row's strand, identity, E and bits, and `modelled_in_run`: the counted call of the run that overlaps the place (its `mrna_id`, `confidence`, `provenance`, and `other_locus`, true when it was found from another home locus), else null. `best_per_genome` and `goi_found_but_not_syntenic` read only the places without one. `dropped` lists the GOI hits outside every searched GOI zone that the rescue dropped, with the test each failed (`rejected_by` `dispersed_identity`, `dispersed_evalue`, `dispersed_alnlen`, or `dead_zone` for a hit whose block the cap dropped): the strongest 200 per genome search, `hits` counting them all. Not orthology calls: they are listed so that a genome with a strong hit that could not be placed is not read as a genome without a hit. |
| `gene_family_advisory` | Always present. `is_gene_family_query` (with `reasons` and `advice`) is true when at least 5 GOI rows, and at least 40 % of them, are LOW, or when the query has at least 3 home loci. Calls in a large family are candidates rather than orthologs. |
| `annotations` | Counts of the `gene` and `mRNA` rows in the GFFs of the search and of the rescue passes, by role, confidence, class and evidence type, in all and per genome. A model of a rescue pass counts twice. |
| `regions` | Counts of the candidate regions per genome and confidence. |
| `synteny_results` | `genes_discovered`: the model proteins per genome. |
| `target_policy` | Per target genome: the median identity of the best hit of each flanking gene (how far the genome is from home; it changes nothing), the blocks searched, how many are significant, the best block's E_syn and whether its window ran past a scaffold end, the genetic code used, and `rescue_dropped_outside_zones` (the GOI hits outside every searched GOI zone that the dispersed-GOI rescue dropped). With any `--modules` other than `legacy` also the genome-axis values (`anchor_max_places`, its source, `blocks_before_cap`, `blocks_rebuilt`); in a decoy run the decoy counts. With several home loci, of the row kept for the genome. |
| `limits` | What each fixed limit of the search dropped: see [below](#what-the-searchs-limits-dropped-limits). |
| `genome_qc`, `qc_summary` | The assembly quality of each target genome. |
| `failed_steps` | Only when a step that may fail without stopping the run did: one row per failure with `step` (`RESCUE_GOI_HULL`, `RESCUE_STRONG_SYNTENY`, `RECIPROCAL_BEST_PARALOG`, `ASSIGN_LOCUS_OWNERSHIP`, or `BUILD_HOME_PARALOG_PANEL` with `--require_paralog_panel false`), `genome`, `locus`, `exit` and `reason` (the last line it wrote to stderr), and a sentence in `warnings`. A failed paralog check marks the calls it would have checked (`ownership_not_checked`). |
| `trees` | Only when a gene tree was built: per home locus, how (from `locus_<N>_tree.info.json`: IQ-TREE's fixed `seed`, `threads`, `sequences`, and `bootstrap`, with `status` `ultrafast`, `none_fewer_than_4_sequences`, `dropped_after_failure` or `placeholder_tree`), its `leaves`, and `not_counted`: the leaves the report does not count, each with its `reason` (`merged` with another leaf of the same gene, `overruled_piece`, `paralog_not_goi`, `identity_coverage_decoupled`, `not_in_report`) and its `record`. The tree page draws those leaves hollow. |
| `genome_order` | How the target genomes were ordered into the search's waves: the `strategy` (`ete3_taxonomy` from the local taxonomy database, `accession_taxonomy` looked up at NCBI, `minhash` or `alphabetical`) and `network_failures`, the NCBI lookups that did not answer (with a sentence in `warnings`). The order decides which genomes seed the later waves. |
| `downloads` | Easy Mode only: every genome the run chose, with `accession`, `species`, `genome` (`ok` or `failed`) and `gff` (`ok` or `none`). A failed download and a genome without a GFF have a sentence in `warnings`. |
| `home_gff_unmatched` | Only when the home GFF has genes on sequences that the home FASTA does not hold: how many, and up to five examples. Those genes are neither flanking genes nor competitors of a call. A GFF that shares no sequence with the FASTA stops the run. |
| `home_genome_filter` | Easy Mode only: what the chromosome filter removed from the home assembly (`removed_sequences`, `removed_bp`). |
| `home_annotation` | Only in a run without a home annotation: where the home annotation and proteome came from (`source` `predicted` or `predicted+borrowed`, the target it was borrowed from, the counts) and `partial: true`, with a sentence in `warnings`. The home proteome then holds only the genes predicted near the query's hits and the genes borrowed from one target, so `BestOtherHome` and `RealCompetitor` compare a call with these alone. Borrowed genes have ids `borrowed_<target gene>`. |
| `run_settings` | What decided the run's settings ([below](#run_settings)). |
| `copies_in_slot` | Only with `--report_copies_in_slot`: per genome the copies of the gene family at home and in the slot of the genome's best call (`unit`, e.g. `2:3`). Advisory: equal rank is not positional orthology. |
| `genome_wide_path` | Only with `--goi_genome_wide_path`: the places of GOI hits outside every searched window, with the home protein's E-value over the genome, the specificity terms, the rank, the cap (home family size + 1) and what became of each (`modelled`, `cut_by_cap`, `not_specific`, `not_significant`, or `no_window` for a kept place where no window could be cut). Their calls carry `GenomeWide=1` and `Place=genome`; this place has no decoy calibration. |
| `decoy_rate` | Only with `--decoy_mode`: see [Decoy calls](#decoy-calls---decoy_mode). |
| `staging_diagnostics` | Plumbing diagnostics, useful when something looks wrong. |

### `run_settings`

What decided the run's settings, from `intermediate/locate_gene/effective_params.json`; null when that file is missing.

- `modules_mode` (the value of `--modules`), `modules_applied` and `modules_would_fire` (the preset modules applied, and every module whose trigger fires, each with its reason), `preset_applied` and `apply_reason`, and `legacy_preset`, the preset the old chooser would pick.
- `settings` and `source`: the value of every parameter the modules may set, and where it came from (`default`, `module:<id>`, `preset:<name>` or `user_set`). `user_set_kept` lists values of yours that were kept although a module or preset wanted another.
- `features` (what the triggers read: the query's features, the home loci, the home neighbourhood) and `home_copies` (the copies of the gene family at home).
- `profile` (the quantities the module rules read), `derivation` (for every value that a module rule, the preset or you set: the value, the default, and for a rule its formula and inputs, e.g. `max_goi_per_genome 11 = max(10, c + 1), c = 10 (A1)`) and `constraints` (the checks of how the settings relate, each with `applies`, `ok` and `detail`).
- `flags`: what the triggers found (`domain_shared`, `dispersed_family_loci`, `home_contig_end`, `home_subtelomeric`, `dist_to_scaffold_end`).
- `smith_waterman`: what ran the Smith-Waterman search (`backend`, `parasail_version`, and `simd`, the instruction set parasail used on that machine); `{"mixed": [...]}` when the search tasks ran on different machines.
- `run_resources`: the resource values a profile sets, such as `nextflow_profile`, `executor`, `mmseqs_sensitivity`, `iterative_search_cpus`, `iterative_search_memory`, `compute_tree_cpus`, `skip_tree` and `nextflow_version`. The MMseqs2 sensitivity changes which hits are found; the thread count and memory do not.
- `tool_versions`: the versions of the tools the run found (`python`, `parasail`, `mmseqs`, `tblastn`, `blastp`, `miniprot`, `samtools`, `mafft`, `iqtree2`, `iqtree`, `prodigal`, `augustus`, `datasets`; `not found` for a tool that is not installed). The tree step tries `iqtree2` first, then `iqtree`; the pinned IQ-TREE 3 installs `iqtree`, so `iqtree2` reads `not found` there. `environment.yml` pins each of the others, so an environment built from it runs these versions.
- With any `--modules` other than `legacy` and a home genome without a GFF, `settings`, `source` and `derivation` also hold `augustus_species`, the Augustus model of the home's gene prediction (`home_lineage`: the model nearest the home species in the NCBI taxonomy; `user_set`; or `default`, `fly`, without the taxonomy).

How the modules work: [PARAMETERS.md](PARAMETERS.md#modules) and [ALGORITHM.md § 5.3](ALGORITHM.md#53-parameter-preset-and-preset-modules).

---

## Gene models

### The per-genome GFF

**`plot_inputs_synteny_block_locus_<N>/<genome>.gff`** holds every model the search built in that genome, of your gene and of its neighbours: one `mRNA` row and its `CDS` rows per model (a tandem copy is a `gene` row). The CDS rows of a miniprot model carry its reading-frame phase; those of a model built from search hits carry phase 0. The CDS rows of a miniprot model of your gene (`EvidenceType=exon_annotation`) also carry `SpliceAcceptor` and `SpliceDonor`: the two bases before and after the segment, on the gene's strand. A target genome searched with a genetic code other than the standard one carries the header line `##synvoy-genetic-code N` (NCBI table N). The gene-order-only calls that the report leaves out are not in these copies.

The attributes of a model's `mRNA` row:

| Attribute | Meaning |
|---|---|
| `SynVoyRole` | `goi` or `flanking` |
| `Confidence`, `GOIClass`, `EvidenceType`, `ModelStatus` | As in the [per-genome table](#the-per-genome-table) |
| `Identity`, `QueryCoverage`, `Exons` | % identity, share of the query covered, exon count; both refer to the **home protein** (for a model whose `SynVoy_Parent` is a tandem copy, `GOI_copy_<n>`, to that copy) |
| `Bits` | GOI models: bit score of the model's protein against the home protein (Smith-Waterman, BLOSUM62, gaps 11/1). Among overlapping models it decides after the label, `Eeff` and `BlockEsyn` |
| `Place`, `BlockEsyn`, `EPlace`, `EPlaceTblastn`, `SigBlocks`, `Eeff`, `EeffTblastn` | GOI models: the smallest expected place that holds the model (`slot` between its nearest home neighbours, `core`, `window`, or `genome` when its block's E_syn is above 0.1); the block's E-value of gene order; the E-value of the home protein against the model in that place (blastp with composition statistics), and that of the window's best tblastn hit reaching into the model, rescaled to the place; the number of significant blocks searched in the genome; and the two E-values multiplied by it ([ALGORITHM.md § 7.13](ALGORITHM.md#713-e-values-in-the-expected-place)). The label rests on `Eeff` |
| `BestOtherHome`, `SpecMargin`, `RealCompetitor`, `BestOtherInNeighbourhood` | GOI models: the home protein other than your gene's own that matches the model best; the bit-score margin of your gene's home protein over it; whether it is a real competitor (`1`: it aligns to the model at E ≤ 0.01 over the home proteome, or is a home family member of your gene); and whether it lies in the gene's home neighbourhood (between the outermost flanking genes of the home locus, or within `n_flanking_genes` genes of your gene's own home gene). Your gene's own gene is the home gene in `goi_info.json`, else the home protein near-identical to the query; a near-identical copy of it at the home locus is no competitor unless the model favours it by 10 bits. A call seeds the next wave only with no real competitor or a margin of at least 10 bits |
| `OwnCopySetAside`, `OwnCopyBits` | GOI models, only when it happened: the near-identical home copy of your gene that matches the model best and was left out of `BestOtherHome` because the model favours it by less than 10 bits, and by how many bits (signed: `+6.2`, `-1.0`). A trace, never a gate: a call that rests on it is one to review in a home tandem array |
| `SearchQuery`, `SearchQueryIdentity` | Only on a model that the home protein did not find and the gene model of a previously searched genome did: that model's id, and the identity to it. `Identity` is still measured against the home protein |
| `StartCodon`, `StopCodon` | `yes` or `no`: whether the model begins with ATG and ends at a stop codon. SynVoy extends a miniprot model's two ends to the nearest in-frame start and stop when it can do so without crossing a splice site |
| `SyntenyContext`, `BlockFlankingSupport`, `BlockCollinearSupport` | The neighbourhood: how many distinct flanking genes its block holds, and the longest run of them in home order |
| `InferenceReason` | Why it got its confidence: for a GOI model `<Basis>_in_<Place>` (e.g. `significant_specific_in_slot`); `awaiting_label` on a GOI row that the judge could not label (no `Eeff`, so LOW) |
| `Basis` | GOI models: why the label (`significant_specific`, `not_specific`, `gene_order_only`, `not_significant`). `FrameFixed=1` or `FrameBroken=1`: the model's CDS was not in one clean frame; another frame over the same exons was used, or none was clean and the model was judged on its hit regions alone |
| `SynVoy_Parent` | The query or home gene the model derives from |
| `TargetGene`, `TargetProduct`, `TargetID` | The target genome's own annotation at that position (needs target GFFs) |
| `GoiFamilyConsistent`, `GoiFamilyReason` | Family-name check, only with `--goi_family_tokens`; it changes no label |
| `ContextTruncated` | GOI models: `left`, `right` or `both` when the search window of the model's block ran past a scaffold end on that side, so part of the neighbourhood is not in the assembly |
| `ExonsKept`, `CysKept` | GOI gene models: exons of the model / exons of the home gene, and cysteines of the home protein aligned to a cysteine of the model / cysteines of the home protein. Advisory |
| `CoreCoverage` | GOI gene models and hull-rescue models, only with `core_completeness` (set by the `signal_anchor` module): the share of the home protein's core (without signal peptide and GPI tail) that the model covers. Every coverage test (the label, `ModelStatus`, the hull rescue's filter, the report's coverage check) then reads the larger of `QueryCoverage` and this |
| `LowComplexity` | GOI models, only with `report_low_complexity` (set by the `low_complexity` module): the share of the model protein in low-complexity windows. Advisory |
| `CopiesInSlot`, `SlotRank` | GOI models, only with `report_copies_in_slot` (set by the `local_family` module): the GOI models of the block's slot labelled AMBIGUOUS or better, and this model's rank among them, counted from the home neighbour on the left. Advisory |
| `ChainQuery` | GOI hit chains (`fallback_hit_span`): the query interval of each hit of the chain (`1-40,41-88`) |
| `RefusalLifted` | GOI hit chains in a block with E_syn ≤ 0.1 that the fallback gate refused, written all the same: the tests they failed (`short_aln`, `short_bits`, `long_qcov`, `long_aln`, `long_bits`; [ALGORITHM.md § 7.7](ALGORITHM.md#77-from-hits-to-gene-models)). The judge labels them like every call; they never seed and are not in the tree |
| `MissingQueryStart`, `FirstExon`, `MissingQueryEnd`, `LastExon` | Only with `--goi_first_exon_guard`: the query residues before (after) the model where a short home first (last) exon is missing; the model keeps miniprot's start (stop) |
| `JoinedHits` | Only with `--goi_join_leftover_hits`: the query intervals of the leftover hits joined into this model |
| `GenomeWide` | Only with `--goi_genome_wide_path`: `1` on a call of the genome-wide path (`Place=genome`) |
| `DispersedRescue`, `RescueRule`, `RescueAnchors`, `RescueAnchorBp` | GOI models of a block that the dispersed-GOI rescue seeded, outside every GOI zone of the blocks built without it: `1`; the rule (`envelope`, `near_neighbour`, or `no_flanking_hits` in a genome without any hit of a flanking gene); the flanking genes whose best hit lies on that sequence; the distance in bp from the rescued hit to the nearest of them (absent when there is none) |
| `RetrocopyCandidate` | Only with `--goi_genome_wide_without_blocks`: `1` when a CDS or hit of a genome-wide call runs across an exon junction of the home gene, with 10 residues on each side, as a retrocopy would. A candidate only: an ortholog that lost the intron looks the same |
| `Rearranged_from` | Flanking models placed on another scaffold than their block |

### Model ids

`GOI_<name>|<genome>_b<block>_l<locus>_<kind>`: the block of flanking genes, the candidate locus inside its window (`l1` = the best), and how the model was built (`exon_ann` = miniprot model, `fallback` = hit chain). `l1m2`, `l1m3` … are further gene models of the same locus, such as the other copies of a tandem array. `r<n>` in place of `l<n>` marks a model found only with a relative's model as the query (see `SearchQuery`).

A single-segment tandem copy (`EvidenceType=tandem_copy`, a `gene` row) has `Name=GOI_copy_<n>` (the figures label it "GOI #n") and `SynVoy_Parent=<parent>`, the query it was found with. Its id is `GOI_copy_<n>|<genome>_b<block>_l<locus>_by_<parent>` with `--tandem_copy_parent seed` (the default), or `<parent>|<genome>_b<block>_l<locus>_tc<n>` with `--tandem_copy_parent home`.

### The distances of the search

The first lines of each per-genome GFF hold the distances used for that genome, in base pairs, and how they were obtained:

```
##synvoy-distances {"cluster_distance":21897,"dispersed_envelope_margin":28153,"expansion":1.301,"expansion_anchors":15,"expansion_source":"anchors","goi_margin":6256,"home_anchors":20,"home_pitch_bp":2404,"home_span_bp":45681,"region_padding":21897,"search_margin":21897,"source":{…},"synteny_bridge_max_gap":62562}
```

`home_pitch_bp` is the mean step between neighbouring flanking genes at home, and `expansion` how much larger the same neighbourhood is in this genome (`expansion_source`: `anchors`, measured over a run of flanking genes in home order, of `expansion_anchors` genes, at least 5; `genome size`, the ratio of the genome lengths; `none`, 1). `source` says for each distance whether it was `derived`, `set in bp` by you, or `fixed (no home coordinates)`. `region_padding` is the extent of the neighbourhood on each side of a block, and `search_margin` the DNA that was cut out and searched on each side; it is narrower where a wide neighbourhood holds DNA whose hits could not be used. How the distances are derived: [ALGORITHM.md](ALGORITHM.md#distances-in-gene-units).

**The genome axis.** With any `--modules` other than `legacy`, the next line holds the values the search set for this genome, each with its source:

```
##synvoy-genome-axis {"anchor_max_places":{"blocks_before_cap":112,"blocks_rebuilt":64,"source":"genome_axis","value":5}}
```

`source` is `default` (the value the genome started from), `genome_axis` (the search set it: here the block cap bound at 112 blocks, so the blocks were built again from the 5 best places of each flanking gene, which gave 64) or `user_set`. Under `legacy` the line is not written.

### The per-genome table

**`plot_inputs_synteny_block_locus_<N>/<genome>.homology.tsv`** lists every model the search built in that genome, of your gene and of its neighbours, before the report's checks. Use it to see the evidence behind a call or to look at LOW leads. It does not hold the models of the rescue passes or the report's verdicts.

| Column | Meaning |
|---|---|
| `target_id` | The model's id (e.g. `GOI_Melt\|Apis_florea_fna_b0_l1_exon_ann`) |
| `home_id` | The home gene it maps to (`GOI_<query>` for your gene, e.g. `gene-LOC<NNN>` for a flanking gene) |
| `role` | `goi` or `flanking` |
| `confidence` | `HIGH`, `MEDIUM`, `AMBIGUOUS` or `LOW` |
| `goi_class` | GOI rows only; it follows the label: `confident_goi` (HIGH), `probable_goi` (MEDIUM), `syntenic_candidate_unconfirmed` (AMBIGUOUS), `ambiguous_goi_family_member` (LOW, despite its name). A single-segment tandem copy is `tandem_goi_copy` whatever its label |
| `model_status` | `complete`, `partial` or `fragment`: how much of the gene the model recovers |
| `evidence_type` | How the model was built. GOI: `exon_annotation` (a miniprot gene model), `tandem_copy`, `fallback_hit_span` (a chain of search hits), `rescued_exon`, `raw_hit`; all kinds are labelled by the same rule. Flanking: `flanking_miniprot`, `flanking_hit_span`, `rearranged_flanking`, `rearranged_flanking_fallback` |
| `identity` | % amino-acid identity to the home protein |
| `n_exons` | Exon count of the model |
| `synteny_context` | `candidate_region_anchor`, `strong_flanking_support`, `weak_flanking_support`, `cross_chromosome_rearranged`, `no_flanking_support` |
| `block_flanking_support` | Distinct flanking genes found in the model's search block |
| `query_coverage` | Share of the query that the model covers (0 to 1) |
| `target_gene`, `target_product` | The target genome's own annotation at that position, if it had a GFF |
| `embedding_similarity`, `structural_similarity` | ProtT5 and Foldseek scores, when those optional steps are on |

```bash
# HIGH GOI models of all genomes:
awk -F'\t' 'NR==1 || ($3=="goi" && $4=="HIGH")' plot_inputs_*/*.homology.tsv

# HIGH and MEDIUM gene models only (no hit chains):
awk -F'\t' 'NR==1 || ($3=="goi" && ($4=="HIGH"||$4=="MEDIUM") && $7=="exon_annotation")' \
    plot_inputs_*/*.homology.tsv
```

### Rescue models

Two passes search again where the main search found no good model ([ALGORITHM.md § 9](ALGORITHM.md#9-rescue-passes)). Their models are in `rescue/locus_<N>/`, a file only when a pass found a model:

- **The hull rescue**: `<genome>.hull_rescue.gff` and `.hull_rescue.faa` (the protein), with the genome's file name without its FASTA extension (`GCF_003254395.2.hull_rescue.gff`). Its models are judged like the main search's and can be counted; a `synteny_hull_rescue` record in the report has its exons here, not in `plot_inputs_*`. The pass writes no model that repeats a model of the main search (same strand, each covering more than 80 % of the other) with an equal or better label.
- **The strong-synteny rescue**: `<genome>.rescue.gff` (`GCF_003254395.2.fna.rescue.gff`). Its models are always LOW (`GOIClass=relaxed_rescue`): leads, never records.

For a genome browser:

```
Gene models   →  plot_inputs_*/<genome>.gff  and  rescue/locus_*/*.gff
Region spans  →  regions/locus_*/<genome>.regions.bed
Home gene     →  the home GFF (your --home_gff, or home_genome/ in Easy Mode)
```

---

## Regions and blocks

### Region scores

**`regions/locus_<N>/<genome>.scores.tsv`** has one row per candidate region of that genome, and **`<genome>.regions.bed`** the same regions as BED rows. They are named like `Apis_florea|Reg1_G10_CHIGH_S0.70`: the species of the genome (from `species_mapping.tsv`, with `_` for spaces), then the rank, the flanking genes, the confidence and the score. When the species is not known, the name has no prefix. The regions are drawn in the figures and read by the strong-synteny rescue; they do not change the label of any call ([ALGORITHM.md § 8](ALGORITHM.md#8-region-scores)).

| Column | Meaning |
|---|---|
| `region_rank`, `region_name`, `species`, `chrom`, `start`, `end`, `strand` | The region |
| `score` | The ranking score: `synteny_score`, plus 0.15 when the region holds a GOI call, and at least 0.45 for a distant-macrosynteny rescue region |
| `synteny_score` | The neighbourhood evidence alone: `quality_score × coverage_score`. This is what the p-value tests |
| `quality_score` | 0.4 × `coverage_score` + 0.3 × `consistency` + 0.3 × `strand_consistency` (the weights are `synteny_weight_base`, `synteny_weight_consistency` and `synteny_weight_strand`) |
| `coverage_score` | The share of the expected flanking genes found: `unique_genes` / `total_genes_expected` |
| `unique_genes`, `total_genes_expected` | The flanking genes found in the region, or the flanking support of the GOI's search block (`goi_block_flanking`) when that is larger; and the flanking genes of the home locus, both sides |
| `consistency` | Gene-order agreement with the home genome (longest collinear run / genes) |
| `strand_consistency` | Strand agreement with the home genes; a uniformly inverted neighbourhood scores 1.0 |
| `p_value` | Permutation test of `synteny_score` (200 draws, so never below 1/201, about 0.005). It orders regions; it is not corrected for the number of regions tested |
| `goi_overlap`, `is_goi_anchor` | The region holds a GOI call; the region was made from a GOI call alone |
| `high_flanking_count`, `goi_block_flanking` | HIGH flanking models in the region; flanking support of the GOI's search block (the first sort key) |
| `region_class`, `goi_missing` | e.g. `goi_missing_but_strong_synteny`: the neighbourhood is there, but no GOI model (not LOW) or GOI hit lies near it |
| `confidence`, `selection_reason` | The region's own label and why the region was kept |

A row with `is_goi_anchor` true (`region_class` `goi_anchor`) is no synteny cluster: when no cluster overlaps a GOI call, the call itself is written as a region. Its `p_value` (0), `consistency`, `strand_consistency` and scores are placeholders, set from the call's label, not measured.

```bash
# The best region of each genome, by score:
awk -F'\t' 'NR>1 && $1==1 {print $3, $4, $8, $9, $16}' OFS='\t' regions/*/*.scores.tsv | sort -t$'\t' -k3,3nr
```

### Block tables

**`blocks/locus_<N>/locus_<N>_<genome>.blocks.tsv`** is the search's own view of the same question: one row per block it searched in that genome, also one without a call (a genome the search reached without a block has the header only). With `--plot_block_selection true` the figures choose what a target row shows from this table, and a copy is kept in `plot_inputs_*` for a re-draw.

| Column | Meaning |
|---|---|
| `genome` | The target genome (its file name) |
| `block` | The block number: the `N` of `_b<N>_` in the ids of the models found in it |
| `chrom`, `start`, `end` | The block's span, from its first to its last seed hit, 0-based, end exclusive |
| `hood_start`, `hood_end` | The neighbourhood its anchors were counted in (block ± padding); E_syn is computed for this length |
| `flanking_support`, `collinear_support` | Flanking genes with hits in that neighbourhood, and the longest run of them in home order |
| `esyn`, `significant` | The block's E_syn at full precision (`inf` when the window had no GOI query), and 1 when it is at most 0.1: a predicted place, counted in `SigBlocks` |
| `goi_seeded` | 1 when your gene seeded the block (your gene's own home gene among the anchors, or a rescued GOI hit): kept whatever its number of genes |
| `context` | `left`, `right` or `both`: the window ran past a scaffold end there |
| `best_call`, `best_call_id`, `n_high` … `n_low` | The block's best GOI model in the label's order, and its GOI models per label, as the search wrote them: decoy models are not counted, and a gene-order-only call that the report leaves out is still counted here |
| `dispersed_rescue` | The rule of the rescued GOI hit that seeded the block (`envelope`, `near_neighbour`, `no_flanking_hits`); empty for any other block |

Windows of the genome-wide path (`--goi_genome_wide_path`) are not blocks and have no row.

---

## Decoy calls (`--decoy_mode`)

With `--decoy_mode true`, the search also looks for the reversed home protein and two shuffles of it in every window, and for those of the relative's query where a window has one (filed under the home decoy of the same kind; `SearchQuery` names the decoy query, e.g. `GOI_x~rev~rel`). Their models are in `decoys/locus_<N>/locus_<N>_<genome>.decoy.gff`, and one row per call in `locus_<N>_<genome>.decoy.tsv`: the decoy kind, the label, the place, `Eeff`, `EeffTblastn`, bits, specificity, `SigBlocks`, then `PlaceBp` and the E-value in each nested place (`EPlace`, `EPlaceSlot`, `EPlaceCore`, `EPlaceWindow`, `EPlaceGenome`), then `SearchQuery`, then `DispersedRescue` and `RescueRule` (a decoy found only through the dispersed-GOI rescue), then `RefusalLifted` (a decoy hit chain written in a significant block although the fallback gate refused it), so the false-call rate of such calls can be read apart. Each decoy case (home gene × decoy kind) is chosen among its own candidates. No decoy is in the result GFFs or the figures.

Every decoy call is false by construction, so the share of decoy cases that reach a label is the false-call rate of that label on your data. The report gives it as `decoy_rate`: `cases`, the decoy cases searched (home gene × target genome with a search window × decoy kind, also one without a decoy model; the home genome is never searched), and for `high`, `medium_or_better` and `ambiguous_or_better` the `calls` (cases with such a decoy call) and `ci95_upper`, the upper end of their exact (Clopper–Pearson) 95 % interval. 0 calls means "at most that bound", not 0 %. The counts per genome are the `decoy_cases`, `decoy_high`, `decoy_medium_plus` and `decoy_ambiguous_plus` fields of `target_policy`.

---

## What the search's limits dropped (`limits`)

Several fixed limits of the search can drop a hit, a candidate or a block without a trace in the results. Each drop is counted per genome, in one INFO line at the end of each genome in the search log (`Limit counts (what each limit dropped; docs/OUTPUT.md): ...`, only the counts that are not 0, else `none`), and summed in the report: `limits.counts.<counter>.total`, `limits.counts.<counter>.genomes` (the genomes where it is not 0) and `limits.genome_searches` (home loci × genomes). No limit is changed by these counts.

A Smith-Waterman or miniprot run of the search that does not finish within its time guard (`sw_timeout_seconds`, `goi_miniprot_timeout`, `flanking_miniprot_timeout`) stops the task instead of giving fewer hits or no model, and so does a miniprot run killed for memory (exit 137, which Nextflow retries with more memory). An MMseqs2 search that fails for lack of memory is retried with one thread and a sensitivity of at most 7.0, then at most 6.0, so its hits can differ: `mmseqs_genome_lower_rung` and `mmseqs_window_lower_rung` count the searches that ran so, the genome's GFF names them (`##synvoy-reduced-search`), and `warnings` has a sentence. A window whose MMseqs2 search still fails stops the search. The other counts are decided by the inputs. Hits and candidates of decoy queries are left out, except in `regions_cap_decoy`. A limit that never shows a count on your data does not bind there.

| Counter | Counts |
|---|---|
| `miniprot_timeout_goi`, `miniprot_timeout_flanking`, `sw_timeout`, `goi_model_rounds_ended_by_timeout` | Always 0: a run past its time guard stops the task |
| `miniprot_failed_*`, `miniprot_no_alignment_*` | miniprot calls that failed, and those that found no alignment |
| `goi_model_loci`, `goi_models_found`, `goi_models_kept` | GOI loci modelled, the alignments miniprot returned there over all rounds, and the models kept (`goi_max_models_per_locus`, 10) |
| `goi_model_cap_stopped` | GOI loci where the model cap stopped the rounds while an alignment that would have been kept, or a window hit that no model explains, was left |
| `goi_models_below_outs` | Alignments of a masked round below 0.3 × the best score (miniprot's `--outs` bound, applied again) |
| `goi_models_near_outs` | Kept models scoring below 0.4 × the best: close to that bound, a hint at what `--outs` leaves out |
| `window_len_floor_mmseqs`, `_tblastn`, `_sw` | Window hits that pass the E-value and identity floors but align fewer query residues (qend − qstart + 1) than the length floor (15; Smith-Waterman 10) |
| `window_len_floor_lt15_long_query` | Of those, hits of fewer than 15 query residues on a query longer than 100 aa |
| `mmseqs_maxseqs_window`, `mmseqs_maxseqs_genome` | Queries with 300 result rows in a window or in the genome-wide search: MMseqs2 reports at most 300 targets per query (its `--max-seqs` default), so such a query may have lost more |
| `sw_long_query_skipped` | Query × window pairs that Smith-Waterman skipped as longer than 10,000 aa |
| `goi_core_filter` | Window GOI hits outside the block's GOI core (block ± GOI margin), which are not used |
| `regions_cap`, `regions_cap_rel`, `regions_cap_decoy` | GOI loci of one query in a window beyond `goi_candidate_regions_per_query` (2): of the home protein and models of earlier genomes, of the relative's query, of a decoy. Each cut is also logged with its loci (one INFO line per window and query, decoys left out) |
| `fallback_refused`, `fallback_short_aln`, `_short_bits`, `_long_qcov`, `_long_aln`, `_long_bits` | Hit chains refused by the fallback gate outside a block with E_syn ≤ 0.1, and the tests they failed (a query below `fallback_short_query_len`: aligned query residues, bit score; a longer one: query coverage 0.25, aligned query residues, bit score) |
| `fallback_lifted` | Not a drop: hit chains in a block with E_syn ≤ 0.1 that the fallback gate would have refused, written as candidates (`RefusalLifted=`) and judged |
| `fallback_span_cap` | Hit chains refused because their hits span more than 150 bp per query residue, at least 100 kb and at least 1.5 × your gene's span at home × the target's expansion |
| `flanking_loci_cap` | Loci of a flanking gene in a window beyond the two strongest |
| `blocks_capped`, `blocks_capped_goi_seeded`, `blocks_capped_without_home_order` | Blocks dropped by `max_blocks_per_genome`, of them GOI-seeded ones, and those dropped without the home order (the cap keeps the GOI-seeded blocks and those whose GOI zone holds a strong GOI hit first, with or without the home order). Each dropped block also has an INFO line with its position in the list the cap read |
| `dispersed_identity`, `dispersed_evalue`, `dispersed_alnlen`, `dispersed_in_goi_zone` | Genome-wide GOI hits that the rescue of hits far from their anchors (`enable_dispersed_goi_rescue`) does not take, by the first test they fail, in every genome with a flanking hit: identity < `dispersed_goi_min_identity`, E > `dispersed_goi_max_evalue`, fewer aligned query residues than `dispersed_goi_min_alnlen`, or lying in a block's GOI zone (left to that block). The hits among them that lie outside every searched GOI zone are listed in `rejected_candidates.dropped` |
| `dead_zone` | Of the hits left to a GOI zone, those that no searched block's GOI core holds, because the cap dropped their block: strong GOI hits that neither the blocks nor the rescue use |
| `dispersed_refused_per_piece_pass` | A measure only: hits the rescue refused that it would have seeded if each flanking gene counted on every sequence where it has a hit, not only at its best hit in the genome |

One limit outside the search is logged instead: `split_loci.py` warns whenever `--max_loci` drops a home locus (`WARNING: --max_loci N dropped K home locus/loci ...`). The rescue passes count nothing: a miniprot run there that does not finish within its time guard (`strong_synteny_rescue_miniprot_timeout`, and `goi_miniprot_timeout` for the hull rescue) ends the run, and so does a rescue task that the scheduler killed for memory or time a third time.

---

## Intermediate files and logs

`intermediate/` holds the outputs of the steps before the search and is always written.

| Folder | Contents |
|---|---|
| `query/` | `normalized_query.faa`: the cleaned query (translated if it was DNA; only the first record of a file with several) |
| `locate_gene/` | The query's hits on the home genome, `hit_profile.json`, `auto_preset.json`, `home_context.json`, `effective_params.json` |
| `annotate_goi/` | Your gene's exon structure on the home genome (`goi_exons.faa`, `goi_annotation.bed`, `goi_info.json`) |
| `split_loci/` | One BED per home locus (the annotated home gene first; `--home_locus_mode`). A row is a hit interval, one copy of the gene, snapped to the coding span of the annotated gene that covers at least half of it (the hit span then in columns 8 and 9); with several rows, the first, named `best_interval`, is where the walk to the flanking genes starts |
| `flanking/` | The flanking genes of each home locus (`synteny_block_locus_<N>.bed`, proteins) |
| `initial_db/` | The starting query set: the flanking proteins and your gene's protein |
| `phylo_sort/` | The target genomes, closest first, with distances (`sorted_genomes.txt`), and `sort_provenance.json`: the strategy that ordered them, the taxonomy database ETE3 opened (path, size, date, md5), the home TaxID and the order. Two machines with different taxonomy databases can order the search waves differently |
| `qc/` | `filtered_sorted_genomes.txt`: the target genomes the search takes, closest first. A genome that failed the assembly check is left out under `--qc_fail_policy drop`, and under `auto` (the default) when SynVoy chose it, unless every genome failed; one you chose stays in, flagged ([`qc_fail_policy`](PARAMETERS.md#qc_fail_policy)) |
| `phylo_placement/` | The phylogenetic-placement table of each home locus |

Files worth knowing on a healthy run:

| File | Why |
|---|---|
| `annotate_goi/goi_info.json` | How the home gene of the query was found. `method` is `gff_annotation` when an annotated gene was taken, else `hit_annotation` (a miniprot model built at the hits; `tandem_duplication` when that model holds tandem copies of the gene) or `full_protein_only`. `found_by` names the path: `name` (a name of the query in the home GFF), `position` (an annotated gene at a significant hit that agrees with the query), `miniprot` or `none`. For an annotated gene also `gene_id`, `transcript_id`, `identity`, `score_fraction` (the Smith-Waterman score as a share of the query's self-score), `overlaps_hit`, and the best other gene compared (`runner_up_gene_id`, `runner_up_score_fraction`). `home_gene_ambiguous` is true when that runner-up came within 0.05 of the gene taken by position: check which of the two your query is. `gene_id` is your gene's own gene in the specificity terms (`BestOtherHome` is never that gene); with `home_gene_ambiguous`, or without a `gene_id`, the search takes the home protein near-identical to the query instead. `uniprot_lookup`, only when the query id is a UniProt accession: `query_id`, `status` (`ok`, `not_found` or `failed`), and the `names` UniProt gave or the `error`. These names can decide which home gene is taken; when UniProt did not answer, `warnings` in the report says so. |
| `locate_gene/hit_profile.json` | The hits SynVoy read off the home genome (how many loci, how strong, how much of the query each covers) and the query's own features (`query_features`: length, low complexity, signal peptide, GPI tail, cysteines). `locate_retry`, only when it happened: the home genome was searched once more, for a short query without hits, with tblastn without its low-complexity filter (`short_query` module). |
| `locate_gene/home_context.json` | The home gene's neighbourhood: family members of the query among the genes next to it (`relatives_nearby`), and its distance to the scaffold end. |
| `locate_gene/effective_params.json` | The preset modules (or preset) applied and the values the run used ([`run_settings`](#run_settings) copies it). `source` says per parameter whether a value is the default, a module's (`module:A1`), a preset's or yours (`user_set`); `user_set_kept` lists your values that were kept although a module or preset wanted another. Look here first when a parameter seems to be ignored. |

`qc/genome_qc_summary.json` (in the output folder itself) has the assembly statistics of each target genome.

**Logs.** At the end of a run, also a failed one, the logs of the tasks are collected in `logs/`: one folder per task that wrote any output, named after the process, its tag and its first output (`<PROCESS>__<tag>___<output>`), with `stdout.log`, `stderr.log` and `script.sh` (empty files are left out). They are collected from the whole work folder, so a work folder shared with earlier runs adds their tasks too. They are the first place to look when a step failed.

---

## Quick reference

| Question | File |
|---|---|
| Which genomes have the ortholog? | `synvoy_report.json` → `summary.headline`, `goi_dedup.records` |
| Where in the target genome, with exons? | `plot_inputs_*/<genome>.gff`; for a `synteny_hull_rescue` record `rescue/locus_<N>/<genome>.hull_rescue.gff` |
| What is the evidence behind a call, and what else was found? | `plot_inputs_*/<genome>.homology.tsv` |
| Why is there no call in a genome? | `summary` → `absent_*`, `unknown_context_truncated`, `goi_found_but_not_syntenic`, `genomes_not_searched` |
| How good is the synteny in a genome? | `regions/locus_<N>/<genome>.scores.tsv`, `blocks/locus_<N>/` |
| An overview of all genomes | `synteny_block_locus_<N>_anchor_grid.html` |
| One locus in detail | `synteny_block_locus_<N>_synteny_plot.html` |
| The gene tree | `synteny_block_locus_<N>_tree.html`, `locus_<N>_tree.nwk` (needs `-profile standard`) |
| Was a parameter of mine applied? | `intermediate/locate_gene/effective_params.json` |
| Why did a step fail? | `logs/`, `warnings` and `failed_steps` in the report, `qc/genome_qc_summary.json` |

---

## Reproducibility

The same inputs, parameters and profile give the same gene models, labels and report: `PYTHONHASHSEED` is fixed in `nextflow.config`, and the search inside a window runs in one thread. Three things can still differ:

- **The gene tree.** IQ-TREE runs with a fixed seed, but with more than one thread (`compute_tree_cpus`, 2 by default) its branch lengths and bootstrap supports can differ slightly from run to run. With `--compute_tree_cpus 1` the tree is the same every time.
- **The profile.** `low_mem`, the launcher's default, searches with MMseqs2 at sensitivity 7.0 instead of 9.5, which finds slightly different hits.
- **Easy Mode** can differ between dates, because the genomes NCBI offers change.

Record the version you ran (the banner of `./run_synvoy.sh`, for example `SynVoy v1.0.1` and `main@1281e47`): results of different versions are not comparable in detail ([CHANGELOG.md](../CHANGELOG.md)). To continue a run, start it again with `-resume` ([USAGE.md](USAGE.md#7-resuming-and-caching)).
