# SynVoy Output Guide

After a successful run, `--outdir` contains a layered tree of files —
some are headline results you'll cite, others are intermediate artifacts
useful only for debugging. This guide tells you **which file to open for
which question**, in order of how often they're useful.

```
<outdir>/
├── synvoy_report.json                      the adjudicated result (start here)
├── synteny_block_locus_<N>_anchor_grid.html / .svg      main figure
├── synteny_block_locus_<N>_synteny_plot.html            ribbon plot (exploration)
├── synteny_block_locus_<N>_*_with_fragments.*           same figures, fragments drawn
├── synteny_block_locus_<N>_gene_positions.* / _anchor_positions.*
├── locus_<N>_tree.nwk  (+ synteny_block_locus_<N>_tree.html with a real tree)
├── plot_inputs_synteny_block_locus_<N>/    per-genome GFF + homology table
├── regions/                                per-genome region BED + scores
├── rescue/locus_<N>/                       gene models found by the rescue passes
├── qc/genome_qc_summary.json
├── intermediate/                           per-stage artifacts
└── logs/                                   stdout / stderr / script of every task
```

`<N>` is the home locus (most runs have one). Easy Mode adds `query/`,
`home_genome/` and `downloaded_genomes/`.

---

## Confidence levels

Every GOI call carries one of four labels. They mean different things, and only
the first two are counted as orthologs.

| Label | Meaning |
|---|---|
| `HIGH` | A gene model (normally multi-exon) at ≥ 50 % identity in a neighbourhood with conserved, collinear flanking genes. |
| `MEDIUM` | A real but weaker case: lower identity, a single exon, or a hit chain with strong support. |
| `AMBIGUOUS` | A candidate in a conserved neighbourhood whose **own** sequence evidence does not show it is the gene (`syntenic_candidate_unconfirmed`), or a call that matches a home *paralog* better than the query (`paralog_not_goi`). Named in the report, drawn in the plots, **not counted as an ortholog**. A gene lost from a conserved locus leaves the neighbourhood intact, so flanking support alone cannot tell "present" from "lost". |
| `LOW` | A lead: a fragment, a weak hit, a chance alignment. Never cite without manual review. |

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
| `model_status` | `complete` / `partial` / `fragment`. |
| `mrna_id` | The model's ID in the per-genome GFF (or in `rescue/`). |
| `target_gene` | The target genome's own gene name at that position, when target annotations were available. |
| `provenance`, `n_source_loci`, `n_merged_hits` | Which home loci produced the call and how many hits were merged into it. |
| `owning_gene`, `owning_locus`, `owning_is_goi_gene`, `ownership_bit`, `ownership_gap` | Result of the paralog check: the home gene the call matches best, and by what margin. |
| `phylo_verdict`, `phylo_z` | Advisory phylogenetic-placement verdict. |
| `superseded_wide_model` | `true` when a compact model replaced a sprawling rescue model of the same gene. |

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
- **The syntenic region** (not the gene): `regions/<species>.regions.bed`, one row
  per candidate region, named like `<species>|Reg1_G10_CHIGH_S0.70` (rank, flanking
  genes, confidence, score).

Attributes on a model's `mRNA` row:

| Attribute | Meaning |
|---|---|
| `SynVoyRole` | `goi` or `flanking` |
| `Confidence`, `GOIClass`, `EvidenceType`, `ModelStatus` | As in the table above |
| `Identity`, `QueryCoverage`, `Exons` | % identity, fraction of the query covered, exon count |
| `StartCodon`, `StopCodon` | `yes` / `no`: whether the model begins with ATG and ends at a stop codon. SynVoy extends a miniprot model's two ends to the nearest in-frame start / stop when it can do so without crossing a splice site. |
| `SyntenyContext`, `BlockFlankingSupport`, `BlockCollinearSupport` | The neighbourhood: how many distinct flanking genes its block holds, and the longest run of them in home order |
| `InferenceReason` | Why it got its confidence (e.g. `goi_tandem_copy_weak_context`, `fragment_model_with_flanking_only`) |
| `SynVoy_Parent` | The query or home gene the model derives from |
| `TargetGene`, `TargetProduct`, `TargetID` | The target genome's own annotation at that position (needs target GFFs) |
| `GoiFamilyConsistent`, `GoiFamilyReason` | Family-name check, only when family tokens are set |
| `Rearranged_from` | Flanking models placed on another scaffold than their block |

`CDS` rows add `SpliceDonor` / `SpliceAcceptor`.

For a genome browser:
```
Gene models   →  plot_inputs_*/<species>.gff  and  rescue/locus_*/<species>.*.gff
Region spans  →  regions/<species>.regions.bed
GOI on home   →  the home-genome GFF you supplied with --home_gff
```

---

## "How good is the synteny in each species?"

**`regions/<species>.scores.tsv`** — one row per candidate region.

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
awk -F'\t' 'NR>1 && $1==1 {print $3, $4, $8, $9, $16}' OFS='\t' regions/*.scores.tsv | sort -k3,3nr
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
| `genome_qc`, `qc_summary` | Assembly quality of each target |
| `staging_diagnostics` | Plumbing diagnostics — useful when something looks wrong |

The fields of `summary` worth reading:

| Field | Meaning |
|---|---|
| `headline` | One sentence, e.g. `3 high-confidence GOI ortholog annotation(s) (+13 low-confidence) across 3 genome(s).` |
| `high_confidence_goi`, `medium_confidence_goi` | Counts after dedup and the paralog check — the numbers to quote |
| `ambiguous_goi`, `low_confidence_goi` | Named separately, never part of the ortholog count |
| `goi_absent_genomes` | Genomes in which no GOI call was **placed**. Not proof of absence. |
| `goi_found_but_not_syntenic` | Genomes with a strong GOI hit that no synteny gate accepted (see `rejected_candidates`) |
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
| `<locus>_anchor_grid.html` / `.svg` | **The headline figure.** Species × gene grid — rows are species (species tree at left), columns are home genes, the GOI column in red. Arrows are shaded by % identity and styled by confidence (HIGH solid, MEDIUM dashed, AMBIGUOUS palest and most broken, LOW striped). An empty cell means **no ortholog was placed in this neighbourhood**, not that the gene is absent. The GOI column draws each gene model's exons to scale; a call built from search hits, not a gene model, is a thin bar; a cell with many models shows the best one plus `×N`. `×N` above a flanking arrow means N copies of that gene at this locus (hits elsewhere in the genome are not counted; the synteny plot shows them). The GOI model shows coding exons only, so it can have fewer exons than NCBI's exon count, which includes non-coding UTR exons of all transcripts. |
| `<locus>_synteny_plot.html` | Interactive ribbon plot: one track per genome (at most three lanes), ribbons between orthologous flanking genes. Scaffolds that read in the opposite direction to the home genome are flipped so ribbons run straight. For exploring one locus. |
| `<locus>_synteny_plot_view.svg` | Static mirror of the ribbon plot (always written). |
| `*_with_fragments.html` / `.svg` | The ribbon plot and the anchor grid again, **with** `ModelStatus=fragment` GOI models drawn. The main figures hide them: they are single-exon hits that otherwise fill the plot. Written only when a run has fragments. |
| `<locus>_gene_positions.html` / `.svg` | **Gene-position map.** Each genome is a line with one dot per gene at its true position. The grid answers "is the ortholog there?"; this answers "where does it sit, and is the neighbourhood rearranged?" |
| `<locus>_anchor_positions.html` / `.svg` | The two views combined: the aligned grid plus each row's real-position line. |
| `<locus>_tree.html`, `locus_<N>_tree.nwk` | Tree of the GOI sequences (MAFFT + IQ-TREE). With the launcher's default `low_mem` profile the tree step is skipped and the `.nwk` is a placeholder; use `-profile standard` for the real tree. |
| `<locus>_synteny_plot.svg` | Narrow publication SVG. Opt-in: `--pub_svg true`. |
| `<locus>_synteny_matrix.svg` | Presence/absence matrix. Opt-in: `--enable_matrix_plot true`. |

Flanking genes are drawn with their real exons. Runs made before 2026-09-17 stored
a single block per flanking gene, and the plot then drew evenly spaced exons that do
not exist — re-run the search to get real ones.

Figure options that have no pipeline parameter (`--grid_goi_style`,
`--grid_goi_max_models`, `--no_fragment_variant`, `--no_orient_to_home`,
`--no_anchor_grid`, …) go through `--plot_extra_args '…'`; see
[USAGE.md § Visualization](USAGE.md#visualization). To try figure options without
re-running the search, re-draw the figures of a finished run (the originals are
not touched):

```bash
conda activate synvoy_env
scripts/replot.sh results/my_run/plot_inputs_synteny_block_locus_1 results/my_run/replot \
    --home_species "Apis mellifera" --grid_goi_style genomic
```

Every option of `python3 bin/plot_synteny.py --help` can follow the two folders.

---

## Intermediate / debugging artifacts (most users skip)

`intermediate/` holds per-stage outputs and is always written.

| Directory | Contents |
|---|---|
| `query/` | `normalized_query.faa` — the cleaned query (translated if the input was DNA; only the first record of a multi-sequence file) |
| `locate_gene/` | The query's hits on the home genome, `hit_profile.json`, `auto_preset.json`, `effective_params.json` |
| `split_loci/` | One BED per home locus |
| `flanking/` | The flanking genes of each home locus (`synteny_block_locus_<N>.bed`, proteins) |
| `annotate_goi/` | The GOI's exon structure on the home genome (`goi_exons.faa`, `goi_annotation.bed`, `goi_info.json`) |
| `phylo_sort/` | Target genomes ordered closest first, with distances |
| `initial_db/` | The starting query set (flanking proteins + GOI) |
| `qc/` | The targets that passed assembly QC |
| `phylo_placement/` | Per-locus phylogenetic-placement table |
| `resolve_home_locus/` | Only with `--enable_name_locus` / `--home_goi_gene` |
| `estimate_params/` | Only with `--auto_params true`: `estimated_params.json` |

Two files are worth knowing on a healthy run:

| File | Why |
|---|---|
| `locate_gene/hit_profile.json` | The hit distribution SynVoy read off the home genome — how many loci, how strong. |
| `locate_gene/effective_params.json` | Which preset was applied and the parameter values the search actually ran with. `source` says, per parameter, whether a value is the default, a preset's, or yours (`user_set`); `user_set_kept` lists values of yours that were kept although the preset wanted something else. Look here first when a parameter seems ignored. |

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
| "How good is the synteny in species Y?" | `regions/Y.scores.tsv` |
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
