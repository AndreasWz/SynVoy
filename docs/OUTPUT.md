# SynVoy Output Guide

After a successful run, `--outdir` contains a layered tree of files —
some are headline results you'll cite, others are intermediate artifacts
useful only for debugging. This guide tells you **which file to open for
which question**, in order of how often they're useful.

```
<outdir>/
├── synvoy_report.json                      the adjudicated result (start here)
├── synteny_block_locus_<N>_anchor_grid.html / .svg      main figure (full version)
├── synteny_block_locus_<N>_anchor_grid_print.svg        the same at print size (183 mm)
├── synteny_block_locus_<N>_anchor_grid_clades_print.svg GOI copies coloured by gene-tree clade (when informative)
├── synteny_block_locus_<N>_anchor_grid_<variant>_print.svg try-out versions: compact, compact_panel, clades_muted (opt-in)
│       + _print.legend.txt / _print.source.tsv          its legend text and cell values
├── synteny_block_locus_<N>_synteny_plot.html            ribbon plot (exploration)
├── synteny_block_locus_<N>_*_with_fragments.*           same figures, fragments drawn
├── synteny_block_locus_<N>_gene_positions.html / .svg   gene-position map
├── synteny_block_locus_<N>_gene_names.tsv               gene labels used in the figures
├── locus_<N>_tree.nwk  (+ synteny_block_locus_<N>_tree.html with a real tree)
├── plot_inputs_synteny_block_locus_<N>/    per-genome GFF + homology table
├── regions/<locus>/                        per-genome region BED + scores, one folder per home locus
├── rescue/locus_<N>/                       gene models found by the rescue passes
├── qc/genome_qc_summary.json
├── intermediate/                           per-stage artifacts
└── logs/                                   stdout / stderr / script of every task
```

`<N>` is the home locus (most runs have one). Easy Mode adds `query/`,
`home_genome/` and `downloaded_genomes/`.

---

## Confidence levels

Every GOI call carries one of four labels, and `Basis` says why (since 2026-10-14). Only
the first two are counted as orthologs.

| Label | `Basis` | Rule | Report, figures and web say |
|---|---|---|---|
| `HIGH` | `significant_specific` | E ≤ 1e-6 in the expected place, specific with a margin of ≥ 10 bits (or no real competitor), and the model covers ≥ 70 % of the query | ortholog, high confidence |
| `MEDIUM` | `significant_specific` | E ≤ 1e-3 in the expected place, and specific | ortholog, medium confidence |
| `AMBIGUOUS` | `not_specific` | E ≤ 1e-3, but another home gene explains the model as well | needs curation: another home gene fits as well |
| `AMBIGUOUS` | `gene_order_only` | 1e-3 < E ≤ 1 in the slot between the gene's neighbours in a block whose gene order is significant, and specific. Offered only in a genome without a counted call | needs curation: weak match where the gene should be |
| `LOW` | `not_significant` | everything else: a fragment, a weak hit, a chance alignment | no significant match |

**The E-value** (`Eeff`) is that of the home protein against the model, computed in the
smallest place SynVoy predicted for the gene: the slot between the nearest found home
neighbours (± one gene pitch), else the core of the block, else the searched window (each
only in a block with `BlockEsyn` ≤ 0.1, the chance of that much gene order), else the
whole genome; times the number of such blocks in the genome. **Specific** means that no
other protein of the home proteome matches the model as well: no real competitor
(`RealCompetitor=0`), or a margin of at least 3.3 bits (`SpecMargin`); within ±3.3 bits a
model in the slot or core of a significant block is decided for the gene, unless the
competitor lies in the gene's own home neighbourhood (`BestOtherInNeighbourhood=1`, a
tandem array). A model of the report's paralog check (`paralog_not_goi`) or of
high identity over a short slice of the query (`identity_coverage_decoupled`) is not
counted either. Docs: ALGORITHM.md §7.8.

---

## "I just want the orthologs"

Open **`synvoy_report.json`**. It holds the *adjudicated* result: calls from every
home locus merged into one record per gene, checked against the home genome's
paralogs, with the rescue-pass models included.

```bash
# The one-line answer
jq -r '.summary.headline' synvoy_report.json

# One row per adjudicated gene
jq -r '.goi_dedup.records[] | [.genome, .chrom, .start, .end, .confidence, .goi_class, .identity, .query_cov] | @tsv' \
    synvoy_report.json
```

Fields of a record in `goi_dedup.records`:

| Field | Meaning |
|---|---|
| `genome`, `chrom`, `start`, `end` | Where the gene is (1-based, inclusive). |
| `confidence` | `HIGH` / `MEDIUM` / `AMBIGUOUS` (LOW calls are not listed here; they are counted in `summary`). |
| `goi_class` | `confident_goi`, `probable_goi`, `synteny_hull_rescue` (found by the hull rescue), `syntenic_candidate_unconfirmed`, `paralog_not_goi`, `cross_locus_duplicate` (the same gene reached from two home loci), `identity_coverage_decoupled` (high identity over a short slice of the query). `original_goi_class` keeps the class before adjudication. |
| `identity`, `query_cov` | % identity to the query and the fraction of the query the model covers. |
| `bits` | Bit score of the model's protein against the home protein (the `Bits` attribute; empty for runs before 2026-10-14). Of several records of one gene the one with the best confidence, then the best bit score, is kept. |
| `model_status` | `complete` / `partial` / `fragment`. |
| `mrna_id` | The model's ID in the per-genome GFF (or in `rescue/`). |
| `target_gene` | The target genome's own gene name at that position, when target annotations were available. |
| `provenance`, `n_source_loci`, `n_merged_hits` | Which home loci produced the call and how many hits were merged into it. |
| `owning_gene`, `owning_locus`, `owning_is_goi_gene`, `ownership_bit`, `ownership_gap` | Result of the paralog check: the home gene the call matches best, and by what margin. |
| `phylo_verdict`, `phylo_z` | Advisory phylogenetic-placement verdict. |
| `superseded_wide_model` | `true` when a compact model replaced a sprawling rescue model of the same gene. |
| `basis`, `eeff`, `place`, `block_esyn` | Why the call has its label (`significant_specific`, `not_specific`, `gene_order_only`), its E-value in its place, the place (`slot`, `core`, `window`, `genome`) and the E-value of its block's gene order. Empty for runs before 2026-10-14. |

### The per-genome table

**`plot_inputs_synteny_block_locus_<N>/<species>.homology.tsv`** lists *every*
model the search built in that genome — GOI and flanking — before adjudication.
Use it to see the evidence behind a call or to look at LOW-confidence leads. It
does **not** contain the rescue-pass models or the paralog / dedup verdicts; those
are only in the report.

| Column | Meaning |
|---|---|
| `target_id` | Identifier of the model (e.g. `GOI_Melt\|Apis_florea_fna_b0_l1_exon_ann`) |
| `home_id` | Home gene it maps to (`GOI_<query>` for the GOI, `gene-LOC<NNN>` for a flanking gene) |
| `role` | `goi` or `flanking` |
| `confidence` | `HIGH` / `MEDIUM` / `AMBIGUOUS` / `LOW` |
| `goi_class` | GOI rows only: `confident_goi`, `probable_goi`, `syntenic_candidate_unconfirmed`, `tandem_goi_copy`, `ambiguous_goi_family_member` |
| `model_status` | `complete` / `partial` / `fragment` — how much of the gene is recovered |
| `evidence_type` | How the model was built. GOI: `exon_annotation` (a miniprot gene model — the only kind that reaches HIGH), `tandem_copy`, `fallback_hit_span` (a chain of search hits), `rescued_exon`, `raw_hit`. Flanking: `flanking_miniprot`, `flanking_hit_span`, `rearranged_flanking`, `rearranged_flanking_fallback`. |
| `identity` | % amino-acid identity to the home protein |
| `n_exons` | Exon count of the model |
| `synteny_context` | `candidate_region_anchor`, `strong_flanking_support`, `weak_flanking_support`, `cross_chromosome_rearranged`, `no_flanking_support` |
| `block_flanking_support` | Number of distinct flanking genes found in the model's search block |
| `query_coverage` | Fraction of the query covered by the model (0–1) |
| `target_gene` / `target_product` | The target genome's own annotation at that position, if it had a GFF |
| `embedding_similarity` / `structural_similarity` | Optional ProtT5 / Foldseek scores when those layers are enabled |

```bash
# HIGH-confidence GOI models across all species:
awk -F'\t' 'NR==1 || ($3=="goi" && $4=="HIGH")' plot_inputs_*/*.homology.tsv

# HIGH+MEDIUM, gene models only (no hit chains):
awk -F'\t' 'NR==1 || ($3=="goi" && ($4=="HIGH"||$4=="MEDIUM") && $7=="exon_annotation")' \
    plot_inputs_*/*.homology.tsv
```

---

## "I want coordinates of the ortholog (for IGV / liftover / curation)"

- **Gene models with exons**: `plot_inputs_synteny_block_locus_<N>/<species>.gff`.
  GFF3 with one `mRNA` row and its `CDS` rows per model (tandem copies are `gene`
  rows). CDS rows carry the real reading-frame phase, and flanking models carry
  their real exons.
- **Models found by a rescue pass**: `rescue/locus_<N>/<species>.hull_rescue.gff`
  (+ `.faa`, the protein) and `<species>.rescue.gff`. A file is written only when a
  rescue produced a model. A `synteny_hull_rescue` record in the report has its
  exons here, not in `plot_inputs_*`.
- **The syntenic region** (not the gene): `regions/<locus>/<species>.regions.bed`, one row
  per candidate region, named like `<species>|Reg1_G10_CHIGH_S0.70` (rank, flanking
  genes, confidence, score).

Attributes on a model's `mRNA` row:

| Attribute | Meaning |
|---|---|
| `SynVoyRole` | `goi` or `flanking` |
| `Confidence`, `GOIClass`, `EvidenceType`, `ModelStatus` | As in the table above |
| `Identity`, `QueryCoverage`, `Exons` | % identity, fraction of the query covered, exon count. Both always refer to the **home protein** |
| `Bits` | GOI models: bit score of the model's protein against the home protein (Smith-Waterman, BLOSUM62, gaps 11/1). The search, the report and the figures choose among overlapping models by it |
| `Place`, `BlockEsyn`, `EPlace`, `EPlaceTblastn`, `SigBlocks`, `Eeff`, `EeffTblastn` | GOI models: the smallest expected place that holds the model (`slot` between its nearest home neighbours, `core`, `window`, or `genome` when its block's E_syn is above 0.1), the block's E-value of gene order, the E-value of the home protein against the model in that place (blastp with composition statistics; and the window's best tblastn hit rescaled to the place), the number of significant blocks searched in the genome, and the two E-values multiplied by it (docs/ALGORITHM.md §7.13). The label rests on `Eeff` since 2026-10-14 (step 4 of the label rework) |
| `BestOtherHome`, `SpecMargin`, `RealCompetitor`, `BestOtherInNeighbourhood` | GOI models: the home protein other than the GOI's own gene that matches the model best, the bit-score margin of the home GOI protein over it, whether it is a real competitor (`1`: it aligns to the model at E ≤ 0.01 over the home proteome, or is a home family member of the GOI), and whether it lies between the outermost flanking genes of the home locus. A call seeds the next wave only with no real competitor or a margin of at least 10 bits |
| `SearchQuery`, `SearchQueryIdentity` | Only on a model that the home protein did not find and the gene model of a previously searched genome did: that model's id, and the identity to it. `Identity` is still measured against the home protein |
| `StartCodon`, `StopCodon` | `yes` / `no`: whether the model begins with ATG and ends at a stop codon. SynVoy extends a miniprot model's two ends to the nearest in-frame start / stop when it can do so without crossing a splice site. |
| `SyntenyContext`, `BlockFlankingSupport`, `BlockCollinearSupport` | The neighbourhood: how many distinct flanking genes its block holds, and the longest run of them in home order |
| `InferenceReason` | Why it got its confidence: for a GOI model `<Basis>_in_<Place>` (e.g. `significant_specific_in_slot`); `awaiting_label` on a GOI row the judge could not label (no `Eeff`, so LOW) |
| `Basis` | GOI models: why the label (`significant_specific`, `not_specific`, `gene_order_only`, `not_significant`; table above). `FrameFixed=1` / `FrameBroken=1`: the model's CDS was not in one clean frame; another frame over the same exons was used, or none was clean and the model was judged on its hit regions only |
| `SynVoy_Parent` | The query or home gene the model derives from |
| `TargetGene`, `TargetProduct`, `TargetID` | The target genome's own annotation at that position (needs target GFFs) |
| `GoiFamilyConsistent`, `GoiFamilyReason` | Family-name check, only when family tokens are set (`--goi_family_tokens`); annotation only, it changes no label |
| `ContextTruncated` | GOI models: `left` / `right` / `both` when the search window of the model's block ran past a scaffold end on that side, so part of the neighbourhood that was meant to be searched is not in the assembly. Absent when the window is complete |
| `ExonsKept`, `CysKept` | GOI gene models (miniprot): exons of the model / exons of the home gene (tandem copies of the home model not counted; absent when the home gene model is not known), and cysteines of the home protein aligned to a cysteine of the model protein / cysteines of the home protein. Advisory: no label reads them (an intron-poor lineage has fewer exons by nature) |
| `CopiesInSlot`, `SlotRank` | GOI models, only with `--report_copies_in_slot` (on with the `local_family` module): the GOI models of the block's slot (`Place=slot`) labelled AMBIGUOUS or better, and this model's rank among them counted from the home neighbour on the left. Absent in a block without a slot (a neighbour missing on one side). Advisory: equal rank is not positional orthology |
| `Rearranged_from` | Flanking models placed on another scaffold than their block |

`CDS` rows add `SpliceDonor` / `SpliceAcceptor`. A target genome searched with a genetic code other than the standard one carries the header line `##synvoy-genetic-code N` (NCBI table N).

**Model ids.** `GOI_<name>|<genome>_b<block>_l<locus>_<kind>`: the block of flanking
genes, the candidate locus inside its window (`l1` = the best), and how the model was
built (`exon_ann` = miniprot model, `fallback` = hit chain). Two additions since
2026-10-09: `l1m2`, `l1m3` … are further gene models of the same locus (the other
copies of a tandem array; before, only the first was reported), and `r<n>` in place of
`l<n>` marks a model found only with a relative's model as query (see `SearchQuery`).

**The distances of the search.** The first lines of each `<species>.gff` hold a
comment with the distances used for that genome, in base pairs, and how they were
obtained:

```
##synvoy-distances {"cluster_distance":21897,"dispersed_envelope_margin":28153,"expansion":1.301,"expansion_anchors":15,"expansion_source":"anchors","goi_margin":6256,"home_anchors":20,"home_pitch_bp":2404,"home_span_bp":45681,"region_padding":21897,"search_margin":21897,"source":{…},"synteny_bridge_max_gap":62562}
```

`home_pitch_bp` is the mean step between neighbouring flanking genes at home,
`expansion` how much larger the same neighbourhood is in this genome
(`expansion_source`: `anchors` = measured over a run of flanking genes in home order,
`genome size` = ratio of the genome lengths, `none` = 1; `expansion_anchors` is the
number of flanking genes in the longest such run, at least 5 for a measurement), and
`source` says for each
distance whether it was `derived`, `set in bp` by the user, or `fixed (no home
coordinates)`. `region_padding` is the extent of the neighbourhood on each side of a
block and `search_margin` the DNA that was cut out and searched on each side (the same
here; narrower where a wide neighbourhood holds DNA whose hits could not be used). See
[PARAMETERS.md § 5](PARAMETERS.md#distances-in-gene-units).

For a genome browser:
```
Gene models   →  plot_inputs_*/<species>.gff  and  rescue/locus_*/<species>.*.gff
Region spans  →  regions/<locus>/<species>.regions.bed
GOI on home   →  the home-genome GFF you supplied with --home_gff
```

---

## "How good is the synteny in each species?"

**`regions/<locus>/<species>.scores.tsv`** — one row per candidate region (until 2026-10-14 the files of all home loci shared one `regions/` folder and overwrote each other).

| Column | Meaning |
|---|---|
| `region_rank`, `region_name`, `species`, `chrom`, `start`, `end`, `strand` | The region |
| `score` | The ranking score: `synteny_score` plus a bonus when the region holds a GOI call |
| `synteny_score` | The neighbourhood evidence alone: `quality_score × coverage_score`. This is what the p-value tests. |
| `coverage_score` | Fraction of the expected flanking genes found (`unique_genes / total_genes_expected`) |
| `consistency` | Gene-order agreement with the home genome (longest collinear run / genes) |
| `strand_consistency` | Strand agreement with the home genes. A uniformly inverted neighbourhood scores 1.0. |
| `p_value` | Permutation test of `synteny_score` (200 draws, so never below 0.005). It orders regions; it is not corrected for the number of regions tested. |
| `goi_overlap`, `is_goi_anchor` | The region holds a GOI call / was anchored by one |
| `high_flanking_count`, `goi_block_flanking` | HIGH-confidence flanking genes in the region; flanking support of the GOI's search block (the primary sort key) |
| `region_class`, `goi_missing` | e.g. `goi_missing_but_strong_synteny`: the neighbourhood is there, no GOI model was built |
| `confidence`, `selection_reason` | Region-level label and why the region was kept |

```bash
# Best region per species, by score:
awk -F'\t' 'NR>1 && $1==1 {print $3, $4, $8, $9, $16}' OFS='\t' regions/*/*.scores.tsv | sort -k3,3nr
```

Scores from runs before 2026-07-26 are not comparable: the strand and order
terms were redefined then.

---

## "I want the big-picture summary"

**`synvoy_report.json`**, top-level keys:

| Key | What's inside |
|---|---|
| `summary` | The headline and the counts (see below) |
| `goi_dedup` | `records` — the adjudicated calls described above — plus dedup statistics |
| `self_consistency` | `flags`: things the run noticed about its own result — `identity_coverage_decoupled`, `paralog_misassignment`, `locus_reattributed`, `cross_locus_duplicate`, `strong_synteny_no_goi`, `phylo_discordant` |
| `rejected_candidates` | GOI hits that cleared the quality bar but were refused by a synteny gate, with the gate and the distance to the nearest flanking gene. Not orthology calls — listed so a genome with a strong unplaceable hit is not read as a genome with no hit. |
| `gene_family_advisory` | Set when the query looks like a member of a large family, where calls are candidates rather than orthologs |
| `annotations` | Counts by role, confidence, class and evidence type; `per_genome` breakdown |
| `regions`, `synteny_results` | Region statistics per genome |
| `warnings` | Sentences to read before any number: a target genome that failed the assembly quality check and was searched anyway (a genome you chose is kept, see `qc_fail_policy`). Empty on a clean run. Such a genome is also in `summary.low_quality_genomes_searched`, its calls carry `low_quality_assembly: true`, and its row has a dagger (†) in the figures |
| `genome_qc`, `qc_summary` | Assembly quality of each target |
| `run_settings` | What decided the run's settings (from `effective_params.json`): `modules_mode` (`--modules`), `modules_applied` and `modules_would_fire` (preset modules with the reason each fired), `preset_applied` and `apply_reason`, the chooser's `legacy_preset`, `settings` and `source` (per parameter: default, `module:<id>`, `preset:<name>`, `user_set`), `user_set_kept`, the `features` the triggers read (query features, home loci, home neighbourhood), and `home_copies` (copies of the gene family at home: the query's coverage at its home locus, rounded half up, or 1 + the family members next to the home gene, whichever is more). `flags` says what the triggers found: `domain_shared` (relatives share only a domain), `dispersed_family_loci` (the home loci of a family on several places, with the query coverage of each), `home_contig_end` (the home gene lies near a scaffold end), `home_subtelomeric` (and a telomeric repeat sits at that end), `dist_to_scaffold_end`. `null` when the file is missing |
| `copies_in_slot` | Only with `--report_copies_in_slot`: per genome `home_copies`, `copies_in_slot` of the slot that holds the genome's best call with a slot, `unit` (`2:3`) and that call (`slot_of`) |
| `target_policy` | Per target genome (`hits/<genome>.target_policy.tsv` of the search): the median identity of the best hit of each flanking gene (how far the genome is from home; it changes nothing), the blocks searched, how many are significant, the best block's E_syn and whether its window ran past a scaffold end, and the genetic code its translations used (`genetic_code`; see `target_genetic_codes`) |
| `staging_diagnostics` | Plumbing diagnostics — useful when something looks wrong |

The fields of `summary` worth reading:

| Field | Meaning |
|---|---|
| `headline` | One sentence, e.g. `12 ortholog(s) (9 high, 3 medium confidence); 4 need curation (3 gene order only, 1 another home gene fits as well); 40 no significant match; across 33 genome(s).` |
| `high_confidence_goi`, `medium_confidence_goi` | Counts after dedup and the paralog check — the numbers to quote |
| `needs_curation_goi`, `needs_curation_by_basis` | The AMBIGUOUS calls ("needs curation"), never part of the ortholog count, and how many rest on `gene_order_only` and how many are `not_specific`. `ambiguous_goi` is the old name of `needs_curation_goi` (kept for one release) |
| `gene_order_only_dropped` | Gene-order-only calls not offered because their genome has a counted call of the gene |
| `low_confidence_goi` | LOW calls, named separately, never part of the ortholog count |
| `goi_absent_genomes` | Genomes in which no GOI call was **placed**. Not proof of absence. |
| `goi_found_but_not_syntenic` | Genomes with a strong GOI hit that no synteny gate accepted (see `rejected_candidates`) |
| `absent_context_complete`, `unknown_context_truncated`, `absent_no_block` | Genomes without a counted call, by what the search saw there: the best block's window lay inside its scaffold; it ran past a scaffold end (the place of the gene is partly missing from the assembly: unknown, not absent); the search found no block at all |
| `nearest_paralog_only`, `nearest_paralog_only_genomes` | Genomes without a counted call whose best call is a model that another home gene explains as well (`Basis=not_specific`): only the gene's nearest relative was found |
| `paralog_misassignments`, `locus_reattributed_orthologs`, `cross_locus_duplicate_goi` | What adjudication changed |
| `self_consistency_flag_count`, `phylo_discordant_calls` | Advisory flags |
| `*_pre_dedup`, `*_pre_ownership` | The same counts before each adjudication step |
| `total_raw_search_hits` | A plumbing diagnostic, often `0` on a perfectly good run. Not a result. |

```bash
python3 -c "
import json; d=json.load(open('synvoy_report.json')); s=d['summary']
print(s['headline'])
print('no GOI call placed in:', s['goi_absent_genomes'])
print('strong hit, not placed:', s.get('goi_found_but_not_syntenic', []))
print('self-consistency flags:', s.get('self_consistency_flag_count', 0))
"
```

---

## Visualisations

| File | When to use |
|---|---|
| `<locus>_anchor_grid.html` / `.svg` | **The headline figure, full version.** Species × gene grid — rows are species (species tree at left), columns are home genes, the GOI column in red. An arrow is the ortholog placed in that neighbourhood: its colour is the home gene, the number its % identity, a number in brackets under the arrow its query coverage when that is below 80 % (`100` over `(7)` = 100 % identity over 7 % of the home protein), its outline the confidence (solid high, dashed medium, dotted low, most broken ambiguous), a second arrow behind it more than one copy at the locus. An empty cell means **no ortholog was placed in this neighbourhood**, not that the gene is absent. The axis above the grid shows where the home genes really lie; the right-hand column gives each row's scaffold and position. See [The GOI column](#the-goi-column-of-the-anchor-grid) below. |
| `<locus>_anchor_grid_print.svg` | **The same grid at print size**: 183 mm wide, at most 170 mm tall, text of 7, 6 and 5 pt in Arial, as a journal asks for a two-column figure. No title, coordinates or names table: `_print.legend.txt` holds the sentences of a figure legend, `_print.source.tsv` one row per cell and GOI call (identity, coverage, confidence, scaffold, coordinates, exons). When a grid is too wide the spacing shrinks first; then the columns with the fewest orthologs are left out, and the log, the legend text and nothing else say which. Options: `--print_grid_numbers` (identity in every arrow, and low query coverage in brackets under it; where that extra line would break the height limit, the identity gets a `*` instead and the value stays in the source table), `--print_width_mm 89` (one column), `--print_height_mm`, `--no_print_grid`. |
| `<locus>_anchor_grid_clades_print.svg` | **The print grid coloured by the GOI gene tree.** Flanking genes in one grey and in narrower columns; every GOI copy is a small arrow filled with the colour of its clade in the GOI tree, cut and coloured as on the tree page (`--clade_count`, default 4), and its exon structure is drawn in the panel on the right, under the same number and in the same colour, as in the print grid. It is written next to the print grid, not instead of it: the print grid shades by identity, this one groups the copies; with the panel it is about as tall as the print grid. It shows which copies in different genomes group together and which are lineage-specific duplications or relatives. White = a copy that is not a leaf of the tree (a call without a protein model, or one the tree was not built from: the tree is built from the search's HIGH and MEDIUM models, before the rescue passes and the report's checks). Weak calls are pale boxes in their clade colour; the outline is the confidence. A small tree in the legend shows how the drawn clades are related in the GOI tree (topology only). Written only when it shows something the species tree beside the rows does not: a genome with two or more copies in the tree, or a genome whose copies fall in no clade of the placed genomes above and below it while those two share one. With `_clades_print.legend.txt` and `_clades_print.source.tsv` (the print-grid table plus a `clade` column). `--no_print_grid` switches it off too. |
| `<locus>_anchor_grid_compact_print.svg`, `_compact_panel_print.svg`, `_clades_muted_print.svg` | **Try-out versions of the print grid** (opt-in: `--plot_extra_args '--grid_variants all'`, or a comma-separated choice). `compact`: the flanking genes of each side of the GOI are one grey stack with how many of them were placed in that genome's neighbourhood (×N; in the home row all of that side), and the width goes to the GOI copies, whose exon models stand in the GOI column in genomic order, three to a line, so a genome with six copies takes two lines instead of six. `compact_panel`: the same stacks, small numbered arrows in the column and the models in the panel on the right. `clades_muted`: the print grid itself with the flanking genes in muted home-gene colours. All three colour the GOI copies by GOI-tree clade (the compact ones fall back to GOI red when the run has no GOI tree; `clades_muted` is then not written). What the stacks do not show: which neighbours were placed and in what order. With `.legend.txt` and `.source.tsv` (in the compact ones every flanking gene keeps its row, `drawn_as = stack`). |
| `<locus>_synteny_plot.html` | Interactive ribbon plot: one track per genome (at most three lanes), every gene at its real size with its exons, one scale for all genomes, ribbons between orthologous flanking genes, GOI ribbons in red. Scaffolds that read in the opposite direction to the home genome are flipped so ribbons run straight. The GOI is named once, on the home track; another track carries only its copy count (`×2`, `+1 weak`) when it has one. For exploring one locus on screen. It has no print version: the figure for a paper is the anchor grid. |
| `<locus>_synteny_plot_view.svg` | The ribbon plot as a picture (always written), for a preview or a README. Not laid out for print. |
| `*_with_fragments.html` / `.svg` | The ribbon plot and the anchor grid again, **with** `ModelStatus=fragment` GOI models drawn. The main figures hide them: they are single-exon hits that otherwise fill the plot. Written only when a run has fragments. |
| `<locus>_gene_positions.html` / `.svg` | **Gene-position map.** Each genome is a line with one dot per gene at its true position. The grid answers "is the ortholog there?"; this answers "where does it sit, and is the neighbourhood rearranged?" |
| `<locus>_tree.html`, `locus_<N>_tree.nwk` | Tree of the GOI sequences (MAFFT + IQ-TREE). With the launcher's default `low_mem` profile the tree step is skipped and the `.nwk` is a placeholder; use `-profile standard` for the real tree. |

Flanking genes are drawn with their real exons. Runs made before 2026-09-17 stored
a single block per flanking gene, and the plot then drew evenly spaced exons that do
not exist — re-run the search to get real ones.

Every figure has its legend under the plot, in titled groups: **Layout** (what a
row, a column, a line is), **Flanking genes**, **GOI** and **Numbers**.

### The GOI column of the anchor grid

SynVoy reports GOI *calls*: everything in a genome's neighbourhood that looks like
the GOI. Many are weak hits on another scaffold. The figures count a call as a
**copy** when it is the best call of its genome or has at least medium confidence;
every other call is a **weak call**.

- **One small arrow per copy** in the GOI column, in genomic order and pointing in
  its coding direction, numbered along the chromosome when a genome has several.
- **A panel right of the grid** draws the exon structure of each copy on a line of
  its own, 5′→3′, under the same number: coding exons to scale, `∧` for an intron,
  the identity to the query after it (`56 (44)` = 56 % identity over 44 % of the
  query). A call built from search hits, not a gene model, is a thin bar.
- **A small pale mark per weak call** at its place (a box for a gene model, a bar
  for aligned hits); at most eight per cell, then `+N`.
- **`//`** stands before what is not at the locus: another scaffold, or farther
  away than the neighbourhood is long.
- **When no genome has two copies** (the home genome included) there is no panel:
  the gene model is drawn in the GOI column itself, in its genomic orientation.
- A genome with more than ten copies shows the first ten models and `+N more`
  (`--grid_max_models`); every copy keeps its small arrow and its row in the
  source table.
- **`–`** stands in the cell of a genome that was searched and has no GOI call.
  Its row is drawn all the same, with the flanking genes that were placed; a row
  in which nothing was placed reads "no flanking gene placed". The legend text
  (`*_print.legend.txt`) names these genomes and how many flanking genes each
  holds. A missing model is not evidence that the gene is absent. Until
  2026-10-09 such genomes were left out of every figure; `--hide_goi_absent_tracks
  true` brings that back.

The same rule counts copies elsewhere: `×N` in the ribbon plot's GOI label and in
the gene-position map is the number of copies, `+N weak` the weak calls. Earlier
versions drew one model per cell and `×N` for all calls, which read as N copies;
their options `--grid_goi_style` and `--grid_goi_max_models` no longer exist.
The GOI model shows coding exons only, so it can have fewer exons than NCBI's exon
count, which includes non-coding UTR exons of all transcripts.

### Gene names in the figures

Every figure labels a home gene with a gene name, not with its annotation ID. The
static figures (anchor grid, gene positions, tree) list the names in a **Gene names**
table under the plot: label,
full name, gene ID and where the label comes from. The interactive ribbon plot
(`*_synteny_plot.html` and its `_view.svg` copy) has no table: hover over a gene
for its label, full name and gene ID. The table is also written to
`<locus>_gene_names.tsv` (and copied to `plot_inputs_<locus>/gene_names.tsv`).

A label comes from the first source that names the gene:

| Source (as printed in the table) | What it is |
|---|---|
| `names table` | A table you supply with `--gene_names_tsv` (see below). |
| `home GFF` | The gene symbol in the home GFF (`gene=` / `Name=`), e.g. `Adh`, `ACT1`. |
| `NCBI Gene` | The gene's current symbol at NCBI Gene, looked up by the `GeneID` in the GFF. Asked only for genes the GFF leaves without a symbol (one request per figure run). Annotation releases are frozen while NCBI keeps naming genes, so this recovers many names for genomes whose GFF is mostly `LOC…` ids. |
| `product name` | No symbol anywhere, but the annotation has a product name: SynVoy abbreviates it. **These labels end in `*`**: they are not official gene symbols. |
| `no name annotated` | "uncharacterized LOC…" and the like: the gene keeps its ID. |

How a product name is abbreviated:

1. A gene symbol at the end of the name is used as it stands
   (`E3 ubiquitin-protein ligase MYLIP` → `MYLIP*`, `segmentation protein Runt` → `Runt*`).
2. A short name with one real word is kept whole (`protein lozenge` → `lozenge*`,
   `formin-2` → `formin-2*`, `annexin B9` → `annexin B9*`): up to 16 characters
   for a single word, 12 with a designator after it.
3. Otherwise the name is split at spaces, hyphens and slashes. Filler words
   (protein, of, and, …) are dropped; numbers, single letters and existing
   abbreviations are kept whole; every other word gives its first letter, and
   "like" gives `L`. The result is upper case:
   `pancreatic triacylglycerol lipase` → `PTL*`,
   `coiled-coil domain-containing protein 170` → `CCDC170*`,
   `zinc finger protein 629-like` → `ZF629L*`.
   If that leaves fewer than four characters and the name has only one or two
   words, the first word gives its first three letters: `arginine kinase` → `ARGK*`.
4. Two genes with the same abbreviation are numbered in genomic order
   (`lozenge-1*`, `lozenge-2*`). No abbreviation is longer than 12 characters.

In the ribbon plot a gene of a *target* genome keeps the name its own annotation
gives it, when target GFFs were supplied: its symbol, else the abbreviation of its
product name (also marked `*`); otherwise it carries the label of the home gene it
is the ortholog of. Its tooltip names that home gene.

To change a label, edit a copy of `<locus>_gene_names.tsv` (columns `gene_id`,
`label`, `full_name`, `source`; set `source` to `user` on the rows you change) and
re-draw with `--gene_names_tsv <file>`. A two-column table (`gene_id`, `label`) is
enough. Genes the table does not list are named as usual.

Without network access the NCBI lookup is skipped or fails; the figures are still
drawn, from the GFF alone, and both the log and the line under the table say so.
`--no_network` skips the lookup on purpose, `--no_gene_legend` leaves the table out.

Figure options that have no pipeline parameter (`--print_grid_numbers`,
`--print_width_mm`, `--grid_max_models`, `--no_fragment_variant`, `--no_orient_to_home`,
`--no_anchor_grid`, `--gene_names_tsv`, `--no_gene_legend`, …) go through
`--plot_extra_args '…'`; see
[USAGE.md § Visualization](USAGE.md#visualization). To try figure options without
re-running the search, re-draw the figures of a finished run (the originals are
not touched):

```bash
conda activate synvoy_env
scripts/replot.sh results/my_run/plot_inputs_synteny_block_locus_1 results/my_run/replot \
    --home_species "Apis mellifera" --print_grid_numbers
```

Every option of `python3 bin/plot_synteny.py --help` can follow the two folders.
A re-draw reuses the run's own `gene_names.tsv`, so the labels do not change. For a
run without that table, or one whose figures were drawn without network access (a
cluster node, for instance), the re-draw looks the names up at NCBI Gene.

---

## Decoy calls (`--decoy_mode`)

With `--decoy_mode true` the search also looks for the reversed home protein and two
shuffles of it in every window. Their models are in `decoys/<locus>/<genome>.decoy.gff`,
one row per call in `<genome>.decoy.tsv` (decoy kind, label, place, `Eeff`, `EeffTblastn`,
bits, specificity). None of them is in the result GFFs, the report or the figures. Each
is false by construction: the share of decoy cases that reach a label (or an E-value
threshold) is the false-call rate of that rule on this data.

## Intermediate / debugging artifacts (most users skip)

`intermediate/` holds per-stage outputs and is always written.

| Directory | Contents |
|---|---|
| `query/` | `normalized_query.faa` — the cleaned query (translated if the input was DNA; only the first record of a multi-sequence file) |
| `locate_gene/` | The query's hits on the home genome, `hit_profile.json`, `auto_preset.json`, `home_context.json`, `effective_params.json` |
| `split_loci/` | One BED per home locus (the annotated home gene is the first one; `--home_locus_mode`) |
| `flanking/` | The flanking genes of each home locus (`synteny_block_locus_<N>.bed`, proteins) |
| `annotate_goi/` | The GOI's exon structure on the home genome (`goi_exons.faa`, `goi_annotation.bed`, `goi_info.json`) |
| `phylo_sort/` | Target genomes ordered closest first, with distances |
| `initial_db/` | The starting query set (flanking proteins + GOI) |
| `qc/` | The targets that passed assembly QC |
| `phylo_placement/` | Per-locus phylogenetic-placement table |
| `estimate_params/` | Only with `--auto_params true`: `estimated_params.json` |

Two files are worth knowing on a healthy run:

| File | Why |
|---|---|
| `locate_gene/hit_profile.json` | The hit distribution SynVoy read off the home genome — how many loci, how strong, how much of the query each locus covers — and the query's own features (`query_features`: length, low complexity, signal peptide, GPI tail, cysteines). |
| `locate_gene/home_context.json` | The home gene's neighbourhood: family members of the query among the genes next to it (`relatives_nearby`, tested with the label's family test), and its distance to the scaffold end. |
| `locate_gene/effective_params.json` | Which preset or preset modules were applied and the parameter values the search actually ran with. `source` says, per parameter, whether a value is the default, a module's (`module:A1`), a preset's, or yours (`user_set`); `user_set_kept` lists values of yours that were kept although a preset or module wanted something else; `modules_would_fire` lists every module whose trigger fires, also under `--modules legacy`. Look here first when a parameter seems ignored. |

`qc/genome_qc_summary.json` (top level) has the assembly statistics of each target.
`logs/<PROCESS>__<tag>___<output>/` holds `stdout.log`, `stderr.log` and `script.sh`
for every task — the first place to look when a step failed.

---

## Quick reference — "I want X, open Y"

| Question | File |
|---|---|
| "Which species have the ortholog?" | `synvoy_report.json` → `summary.headline`, `goi_dedup.records` |
| "Where in the target genome, with exons?" | `plot_inputs_*/X.gff`; for a `synteny_hull_rescue` record, `rescue/locus_*/X.hull_rescue.gff` |
| "What is the evidence behind a call / what else was found?" | `plot_inputs_*/X.homology.tsv` |
| "How good is the synteny in species Y?" | `regions/<locus>/Y.scores.tsv` |
| "Overview of all species" | `<locus>_anchor_grid.html` |
| "One locus in detail" | `<locus>_synteny_plot.html` |
| "Ortholog phylogeny" | `<locus>_tree.html` / `locus_<N>_tree.nwk` (needs `-profile standard`) |
| "Was a parameter of mine applied?" | `intermediate/locate_gene/effective_params.json` |
| "Why did a step fail on species Z?" | `logs/`, `qc/genome_qc_summary.json` |

---

## Reproducibility

The same inputs and parameters give the same result: two cluster runs of the
19-genome melittin benchmark under different Python hash seeds produced 20 of 20
identical GFFs (2026-09-10). `PYTHONHASHSEED` is pinned in `nextflow.config`.
Easy Mode can still differ between dates, because the genomes NCBI offers change.

Re-run with the same `--outdir` and `-resume` to pick up where you left off.
Record the version you ran with (`git rev-parse --short HEAD`, or the banner
`./run_synvoy.sh` prints): results from different versions are not comparable in
detail.
