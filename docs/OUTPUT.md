# SynVoy Output Guide

After a successful run, `--outdir` contains a layered tree of files —
some are headline results you'll cite, others are intermediate artifacts
useful only for debugging. This guide tells you **which file to open for
which question**, in order of how often they're useful.

---

## "I just want the orthologs"

Open **`plot_inputs_synteny_block_locus_<N>/<species>.homology.tsv`** for each species.
That's the canonical per-target table: every gene SynVoy associated with a
home-genome ortholog (GOI or flanking), with confidence, identity, and evidence.

Columns:

| Column | Meaning |
|---|---|
| `target_id` | Internal identifier of the call (e.g. `GOI_Melt\|Apis_florea_fna_b0_l1_exon_ann`) |
| `home_id` | Home-genome gene this ortholog maps to (`GOI_<query>` for the GOI itself, `gene-LOC<NNN>` for flanking genes) |
| `role` | `goi` or `flanking` |
| `confidence` | **`HIGH` / `MEDIUM` / `LOW`** — primary filter for downstream analysis |
| `goi_class` | Sub-class for GOI rows: `confident_goi`, `probable_goi`, `tandem_goi_copy`, `ambiguous_goi_family_member` |
| `model_status` | `complete` / `partial` / `fragment` — how much of the gene is recovered |
| `evidence_type` | `exon_annotation` (best), `flanking_miniprot`, `tandem_copy`, `fallback_hit_span`, `rescued_exon` |
| `identity` | % amino-acid identity to the home query |
| `n_exons` | Predicted exon count |
| `synteny_context` | `candidate_region_anchor` (best), `strong_flanking_support`, `weak_flanking_support` |
| `block_flanking_support` | Number of conserved flanking genes supporting this call |
| `query_coverage` | Fraction of the home query covered by the recovered model (0–1) |
| `target_gene` / `target_product` | Pre-existing annotation if the target genome had a GFF |
| `embedding_similarity` / `structural_similarity` | Optional ProtT5 / Foldseek scores when those features are enabled |

### Filtering for paper-quality calls

```bash
# All HIGH-confidence GOI orthologs across all species:
awk -F'\t' 'NR==1 || ($3=="goi" && $4=="HIGH")' plot_inputs_*/*.homology.tsv

# HIGH+MEDIUM, exon-annotation only (excludes rescue paths):
awk -F'\t' 'NR==1 || ($3=="goi" && ($4=="HIGH"||$4=="MEDIUM") && $7=="exon_annotation")' \
    plot_inputs_*/*.homology.tsv
```

LOW-confidence rows are useful for *exploratory* analysis (e.g. divergent
toxin families like the U11-myrmicitoxin case where 26%-identity calls were
genuine orthologs that BLAST missed) but should never be cited without
manual review.

---

## "I want coordinates of the ortholog (for IGV / liftover / curation)"

Two options depending on granularity:

- **Per-gene exon coordinates**: open the per-target GFF
  `plot_inputs_synteny_block_locus_<N>/<species>.gff`. Standard GFF3 with
  `gene` / `mRNA` / `CDS` features and SynVoy attribute fields (`Confidence`,
  `Identity`, `SynVoy_Parent`, etc.) on every row.
- **Per-region BED** (the syntenic block, not the gene): `regions/<species>.regions.bed`
  has one row per candidate region with start/end/strand and a name like
  `<species>|Reg1_G10_CHIGH_S0.70` (region rank, gene count, confidence, score).

For loading into a genome browser:
```
Tracks/genes  →  plot_inputs_*/<species>.gff
Region spans  →  regions/<species>.regions.bed
GOI on home   →  the home-genome GFF you supplied with --home_gff
```

---

## "How well did SynVoy do per species?"

Open **`regions/<species>.scores.tsv`** — one row per candidate region, with
all the synteny-scoring components (coverage, uniqueness, consistency,
strand consistency, p-value).

Useful columns:

| Column | Meaning |
|---|---|
| `score` | Composite synteny score (0 – 1) |
| `quality_score` / `coverage_score` | Component scores |
| `unique_genes` / `total_genes_expected` | Flanking-gene recovery (e.g. 10/11) |
| `consistency` | Order-preservation of flanking genes |
| `strand_consistency` | Fraction of flanking genes on the expected strand |
| `p_value` | Permutation-test p-value for the synteny score |
| `goi_overlap` | `True` if the GOI was located inside this region |
| `is_goi_anchor` | `True` if this region's score was anchored by a GOI hit |
| `confidence` | `HIGH` / `MEDIUM` / `LOW` (region-level) |

For ranking species by quality:
```bash
# Top region per species, sorted by score descending:
awk -F'\t' 'NR>1 && $1=="1" {print $3, $7, $8, $14}' OFS='\t' regions/*.scores.tsv | sort -k2,2nr
```

---

## "I want the big-picture summary"

**`synvoy_report.json`** — one machine-readable JSON with the whole run:

| Top-level key | What's inside |
|---|---|
| `summary` | Counts of GOI annotations by confidence/class, list of `goi_absent_genomes`, etc. |
| `annotations` | `role_counts`, `goi_confidence_counts`, `goi_class_counts`, `goi_evidence_counts`, per-genome breakdown |
| `regions` | Per-genome region statistics (number of candidate regions, top score, etc.) |
| `synteny_results` | Per-genome synteny-scoring outcomes |
| `genome_qc` / `qc_summary` | Genome-quality-control results from `ASSESS_GENOME_QUALITY` |
| `staging_diagnostics` | Plumbing diagnostics — useful when something looks wrong |

For a quick health check:
```bash
python3 -c "
import json; d=json.load(open('synvoy_report.json'))
s=d['summary']
print(f'Genomes searched: 19, GOI absent in: {len(s[\"goi_absent_genomes\"])}')
print(f'Confidence: {d[\"annotations\"][\"goi_confidence_counts\"]}')
"
```

---

## Visualisations

| File | When to use |
|---|---|
| `<locus>_synteny_plot.html` | Interactive track-style plot (legacy). Best for an exploratory single-locus view with ribbons between flanking orthologs. Resolved-vs-ambiguous GOIs distinguished by hatched fill. |
| `<locus>_synteny_plot_view.svg` | **Static SVG mirror of the HTML** (auto-generated). Same layout, same colours, drops cleanly into READMEs / Word / Inkscape. CSS is CDATA-embedded so the file is fully standalone. |
| `<locus>_synteny_plot.svg` | Narrow publication-format SVG. Render with `--pub_svg`. Different layout (vertical, condensed) optimised for two-column journal figures. |
| `<locus>_tree.html` | Standalone phylogenetic tree of the GOI sequences. Midpoint-rooted, with subclade colouring matching the matrix plot. |
| **(matrix plot)** | Run `bin/plot_synteny_matrix.py` separately on the same `plot_inputs_*` directory — gives a phylogeny-anchored matrix view of all species at once, with HIGH/MEDIUM/LOW confidence visually distinct. The fastest "is the gene present in this clade?" view. |

For paper figures, prefer the SVG export (`--pub_svg` flag) and the matrix
plot — both render cleanly at print size.

---

## Intermediate / debugging artifacts (most users skip)

`intermediate/` contains per-stage outputs of the pipeline. Useful if a run
behaved unexpectedly:

| Directory | Contents |
|---|---|
| `query/` | Resolved query (after auto-translation if input was DNA) |
| `locate_gene/` | The GOI's hit on the home genome (`tblastn` + `MMseqs2` results) |
| `split_loci/` | Multi-locus splitting evidence (when the query maps to multiple positions) |
| `flanking/` | The `n` flanking genes either side of each GOI hit |
| `annotate_goi/` | GOI exon structure (`goi_exons.faa`, `goi_annotation.bed`, `goi_info.json`) |
| `phylo_sort/` | Phylogenetic ordering of target genomes (closest first) |
| `initial_db/` | The starting MMseqs2 query DB (flanking genes + GOI) |
| `qc/` | Per-genome quality-control flags |

---

## Quick reference — "I want X, open Y"

| Question | File |
|---|---|
| "Which species have the GOI ortholog?" | `plot_inputs_*/X.homology.tsv` filtered by `role=goi` and `confidence=HIGH` |
| "Where in the target genome?" | `plot_inputs_*/X.gff` (gene/mRNA/CDS coordinates) |
| "How good is the synteny in species Y?" | `regions/Y.scores.tsv` |
| "Visual overview of all species at once" | matrix plot (run `plot_synteny_matrix.py`) |
| "Detailed per-locus visual" | `<locus>_synteny_plot.html` |
| "Ortholog phylogeny" | `<locus>_tree.html` (or `*_tree.nwk` for tools) |
| "One-line health check" | `synvoy_report.json` → `summary` |
| "Why did SynVoy crash on species Z?" | `intermediate/qc/genome_qc_summary.json` and `logs/` |

---

## Reproducibility

Every run captures the exact parameters used in the Nextflow `params` block.
Re-run with the same `--outdir` and `-resume` to pick up where you left off
without redoing finished steps. The pipeline-level cache key includes input
file hashes + parameter values — changing a parameter invalidates downstream
tasks correctly.
