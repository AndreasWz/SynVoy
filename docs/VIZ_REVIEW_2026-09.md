# Visualisation review for Nature submission (2026-09-16)

> **Progress, 2026-09-17.** Frank's decisions: hide fragments, with a `_with_fragments` sibling;
> keep per-gene colours; draw exons inside the grid from the final GOI model; the grid is the
> main figure and the ribbon plot the exploratory companion. Done in `bin/plot_synteny.py`
> (not committed):
>
> - **0.1** home GOI resolved to annotated gene(s) (`_resolve_home_goi_models`);
> - **0.2** fragments hidden, with `*_with_fragments` ribbon plot and grid
>   (`--no_fragment_variant` turns it off; module emit `with_fragments`);
> - **1.1–1.4** at most 3 lanes with priority, lanes assigned before widening,
>   `--orient_to_home` on by default (`--no_orient_to_home`);
> - **2.1/2.4** GOI column drawn from CDS blocks on one shared scale
>   (`--grid_goi_style cds|cds_aligned|genomic|notched`; option removed 2026-10-07, when
>   the grid got one GOI layout). Hit-based calls get a thin glyph;
>   row labels are italic binomials sized from Arial metrics;
> - AMBIGUOUS GOI calls were styled as HIGH in the grid; fixed;
> - the location column of all three grids printed gap-compressed (and, with orientation,
>   negative) plot coordinates; it now prints genomic ones.
>
> Tests: `tests/test_plot_goi_models_and_lanes.py` (17). Full suite: 773 passed, 4 skipped.
> Drafts for review: `tmp/plot_review/2026-09-17_goi_drafts/README.md`.
>
> **Round 2 (same day).** Decisions: `cds` is the GOI style (now the default); paralogs stay in
> the home row; everything on a reversed scaffold is flipped.
>
> - **0.3 root cause found and fixed.** `iterative_search_runner.py:collapse_flanking_cds_to_gene_span`
>   discarded every in-block flanking model's per-exon CDS rows (while `Exons=N` stayed). The plot
>   then invented N evenly spaced exons. Rearranged flanking models, written after that step,
>   kept theirs. The collapse is removed; no consumer needed it (checked `bin/`, the modules and
>   `scripts/benchmark`). The plot no longer synthesises exons, and hit chains are drawn thinner.
>   `scripts/reproduce_annotation.py` also imported the function; that script was already broken
>   (`flanking_query_utils` no longer exists).
> - **Orientation is decided per scaffold** (`_orient_tracks_to_home`): strand agreement plus
>   order concordance per scaffold. The old per-genome order test left multi-scaffold rows
>   mirrored (DPP4: 3 rows with 8/8–14/14 flanking genes on the opposite strand).
>
> Tests: 21 in the new file; full suite 777 passed. Round-2 renders:
> `tmp/plot_review/2026-09-17_round2/`.
>
> **Round 3: many GOI models.** A GOI cell with more than 3 models draws the best one plus `×N`
> (`--grid_goi_max_models`; option removed 2026-10-07, replaced by the copy rule and
> `--grid_max_models`); the synteny plot draws them all. The 10-per-genome GOI cap used to
> truncate every figure; it now only steers the neighbourhood choice, and all models inside the
> chosen view are restored. Renders: `tmp/plot_review/2026-09-17_round3_many_goi/` (SPIN, SP
> locus 3). Tests 26; full suite 782 passed. Still open: 1.5, Phase 3 (print styling), Phase 4.

All renders use headless Chrome (cairosvg ignores CSS, so it cannot be used to check figures).
The screenshots live in `tmp/plot_review/2026-09-16/` (git-ignored, numbered `01`–`10`).
`regen.sh` there rebuilds every figure from a run's `plot_inputs_*` directory with the local code.

Data reviewed (newest on disk, `local_runs/family_bench/lrz_fix5/`):

| Run | Genomes drawn | Ribbon canvas | Grid canvas |
|---|---|---|---|
| melittin_fix5 | 10 (10 GOI-absent hidden) | 1500 × 1699 px | 1258 × 978 px |
| APYR (apyrase, fix5_ladder) | 33 | 1500 × **8984** px | 1963 × 2174 px |
| DPP4, SP | 34 / 34 | – | – |

The local code gives the same output as the LRZ renders. The plot test files pass (31), and so
does the full suite (756 passed, 4 skipped).

---

## 1. Synteny (ribbon) plot: why genomes stack to 5–7 rows

A "row" here is one greedy interval-scheduling lane per genome track (`_assign_sub_tracks`,
`bin/plot_synteny.py:1465`). Rows per track, measured by instrumenting that function:

| Run | as rendered | GOI fragment models collapsed | + no gene widening |
|---|---|---|---|
| melittin | 1×1, 2×8, **4×1** (home) | same | 1×1, 2×8, 3×1 |
| APYR | 2×1, 3×1, **4×9, 5×15, 6×5, 7×2** | **2×30, 3×3** | 1×7, 2×26 |
| DPP4 | 2×4, 3×27, **4×3** | 2×5, 3×28, **4×1** | ≤ 3 |
| SP | 2×11, 3×23 | 2×14, 3×20 | ≤ 2 |

**Three causes, in order of impact:**

1. **Fragment GOI models are drawn as separate genes.** The APYR target GFFs hold
   **297 `EvidenceType=rescued_exon` / `ModelStatus=fragment` LOW models against 32 HIGH GOI
   models**. They are single-exon hits, 53–600 bp each, 3–35 % query coverage, and many overlap
   one another. In *Cataglyphis* there are 12 of them on the + strand, 6–13 kb from the real
   (− strand, 11-exon, 95 % coverage) apyrase model. Every lane beyond the 3rd in APYR is one of
   these. They are exon-level hits on the ortholog of the neighbouring paralog, which the
   home query also hits (see §2.1). The plot caps GOI entries at `_MAX_GOI_PER_GENOME = 10`, so
   the cap is filled by fragments.
2. **`_widen_sparse_plot` runs before lane assignment.** It widens every gene by up to 4×
   (1.76× melittin, 1.9× APYR) around its centre, then `_assign_sub_tracks` bumps the genes that
   now collide. This is the whole reason the melittin home track needs 4 lanes.
3. **There is no lane cap and no priority.** An unmatched home gene without a ribbon, or a LOW
   fragment, can take a lane as easily as the GOI.

**Other defects in this figure:**

- **Mirrored tracks.** `--orient_to_home` exists but `modules/plot_synteny.nf` never passes it.
  In APYR, 15 of 33 tracks are drawn reverse-complemented relative to home, so the ribbons swing
  in S-curves (screenshot 07). With the flag on, the ribbons run straight down (screenshot 08).
- **Canvas does not fit a page.** At Nature's 183 mm double-column width, APYR is 1,096 mm tall,
  4.4 pages. 33 genomes cannot be shown as a ribbon plot in a main figure.
- **Wasted width.** About a third of the width is the right-aligned label gutter. The x range
  is symmetric around the GOI, so when the GOI sits at the edge of its block, half the plot is
  empty (APYR).
- **Invented exon positions.** Flanking target genes carry `Exons=6` but only one CDS row
  (35 of 36 in *Cataglyphis*). The plot then draws evenly spaced exon blocks
  (`has_synth_exons`, `:2231`), which look like a barcode. This question was left open in June.
  For a paper, drawing exons at made-up positions is not acceptable.
- A rotated label (`★ GOI <gene> ×10`) is repeated on every track. The GOI ribbon is lime
  (a clade colour) while the GOI itself is red. There is a drop shadow, and 52 genes use the
  `geneGloss` gradient. The sub-label repeats the species (`Tetramorium_bicarinatum • OV788322.1`).

## 2. Anchor grid: the "multiple exon" GOI column

Screenshots 03 and 05.

### 2.1 The home GOI is split into query hits, on the wrong strand

The home BED (`synteny_block_locus_1.bed`) carries **one `GOI_` row per tblastn hit** (4 in
APYR), all with strand `+`. The home GFF shows the real picture:

- hits 1 and 2 are exons 1–10 and 11 of the **apyrase gene** (**− strand**);
- hits 3 and 4 (e-49 and e-17) lie in **the next gene** (+ strand), almost certainly a tandem
  paralog.

So the home row draws 4 red arrows pointing right for what are 2 genes, one of them pointing
left. The targets' HIGH model points left, so the grid tells readers the GOI is inverted in
every ant.

### 2.2 The GOI column packs every model into one cell

`draw_goi_array_grid` (`:4407`) lays out every GOI model left to right inside the column as a
"notched arrow" (`_notched_arrow`, `:4067`):

- **The notches read as a bow-tie or saw blade.** Exon junctions are white V-cuts at *evenly
  spaced* positions, not real ones. An 11-exon apyrase gets 10 cuts (the red saw blade in
  APYR); a 2-exon melittin looks like an hourglass.
- **Copies overflow the column.** APYR draws the HIGH model plus up to 12 fragment cells, so
  the GOI "column" spans about 4 columns of width and pushes into the neighbours.
- **The numbers collide.** The font floor is 8.5 px, and a label like `56 (44)` is wider than
  its arrow. Neighbouring labels overprint (`56 (3(95 (56)` in *Apis cerana*,
  `35 (7(56 (47)` in *Cardiocondyla*). At print size 8.5 px is **3.5 pt** (melittin) and
  **2.2 pt** (APYR).

### 2.3 Rest of the grid

- **Row labels run into the cells.** This happens in 7 of 10 rows (melittin) and 29 of 33
  (APYR) with the font that actually renders. Even measured with Arial it is still 3 and 23,
  so it is not only a font problem. The label repeats the genome id whenever the common-name
  lookup is offline (`Bombus impatiens (Bombus_impatiens)`), and the home row embeds the
  accession (`(NC_037641.1)`), which already sits in the right-hand column.
- **The requested font is missing.** The SVG asks for `"IBM Plex Sans"`, which is not
  installed. It falls back to DejaVu Sans, a wide face, and label widths are estimated from
  character counts.
- **Legend clipped.** The last line is cut at the right edge in the melittin grid and in the
  anchor-positions figure (`…listed in synvo`, `…arrow points in cod`).
- **Wasted space.** The column-header box is about 200 px tall of mostly empty panel. The
  "Species tree" caption floats, and species names are not italic.

## 3. Nature readiness (all figures)

Nature's figure rules: 89 mm (single) or 183 mm (double) wide, at most ~247 mm deep;
Arial/Helvetica at **5–7 pt**; editable text; no titles inside the figure (they belong in the
legend); no shadows, gloss or 3-D effects; colour-blind-safe colours, avoiding red/green as
the only contrast. Check the current guide once more at export time.

| Figure | Text at 183 mm width | Main problems |
|---|---|---|
| anchor grid, melittin | labels 7.4 pt ✓, cells 6.8 pt ✓, GOI cells 3.5 pt ✗, legend 4.1–4.9 pt ✗ | label overflow, GOI column, in-figure title and prose legend |
| anchor grid, APYR | labels 4.8 pt ✗, cells 4.4 pt ✗, GOI cells 2.2 pt ✗, legend 2.6 pt ✗ | same, and too much to fit |
| ribbon, melittin | 3.5–4.1 pt ✗ | fits the page height (207 mm), but shadows, gloss, lanes |
| ribbon, APYR | 3.5 pt ✗ | 1,096 mm tall: not a main-figure format |
| gene positions / anchor positions | ~4 pt ✗ | 22–30 px titles; duplicates what the grid shows |

**Figures that overlap.** Every run writes a ribbon plot, anchor grid, anchor positions, gene
positions and a tree page (the matrix is opt-in). Anchor positions and gene positions answer
the same question as the grid plus its right-hand span column. A paper needs one main panel
and at most one Extended Data companion.

---

## 4. Fix plan

The order is deliberate: first make the figures show correct content, then control the
layout, then style for print.

### Phase 0: correct content (prerequisite for every figure)

- **0.1 Home GOI = annotated gene(s).** Map `GOI_*` hit rows to the overlapping home GFF
  gene(s): real strand, real CDS. Hits that fall in a *different* gene (the paralog) become that
  gene's own column, marked as GOI-family if needed, instead of GOI copies. Test: the APYR home
  row has 1 GOI of strand − with 11 exons.
- **0.2 Collapse fragment GOI models for display.** Per genome, a `ModelStatus=fragment` GOI
  model does not become a glyph or count as a copy. It is folded into the nearest non-fragment
  GOI locus (same strand, within the block), or into one "N fragment hits" marker if no such
  locus exists. Details stay in the tooltip. Apply the rule before the
  `_MAX_GOI_PER_GENOME` cap. Ideally use the same rule as
  `generate_report.dedupe_goi_annotations`, so the figure counts what the headline counts.
- **0.3 No invented exons.** A gene without ≥ 2 real CDS rows is drawn as a solid arrow. The
  upstream alternative is to emit per-exon CDS for flanking miniprot models in
  `iterative_search_runner.py`. That changes the GFF output and cache hashes, so it is a
  separate decision.

### Phase 1: ribbon plot, at most 3 lanes

- **1.1** Phase 0.2 alone brings APYR to at most 3 lanes (measured: 30 tracks with 2, 3 with 3).
- **1.2** Assign lanes on **unwidened** coordinates, then widen each gene only as far as its
  lane neighbours allow. Measured: melittin home goes from 4 to 3, DPP4 from 4 to 3.
- **1.3** Hard cap of `MAX_LANES = 3` with priority: GOI and ribbon-bearing orthologs first,
  then unmatched home genes. A gene that still doesn't fit goes into the lane where it overlaps
  least and is drawn translucent; the GOI or an ortholog is never dropped. Regression test:
  every track uses ≤ 3 lanes on the 4 fixture runs.
- **1.4** Pass `--orient_to_home` from `modules/plot_synteny.nf` and make it the default.
- **1.5** One GOI label, placed on the home track only. The GOI ribbon becomes GOI-red. Trim
  the x range to the actual gene extent, not a symmetric range around the GOI.

### Phase 2: anchor-grid GOI column

- **2.1** One glyph per locus, the same arrow size as every other column, with no notches.
  The best model gets the identity number, and real extra copies (non-fragment,
  non-overlapping) get the existing `×N` badge. The tooltip lists all models.
- **2.2** Home GOI strand and exon count come from Phase 0.1.
- **2.3** If exon structure matters for the paper (melittin 2- vs 3-exon, Tb1a), it belongs
  in a **separate small gene-model panel**: exons drawn to scale from real CDS, aligned on the
  home model. It does not belong inside a grid cell.
- **2.4** Row label is the italic binomial only. No `(Genome_id)`, and the accession stays in
  the right-hand column. The label column is sized from real text width (see 3.1), so it can
  never overlap the cells.

### Phase 3: print styling (one shared "pub" theme for every static SVG)

- **3.1** Fonts: `Arial, Helvetica, "Liberation Sans", sans-serif`. Measure text with
  Liberation Sans metrics (Arial-compatible, installed; PIL `getlength`) instead of `len × 8`.
- **3.2** Design at print size: the SVG `viewBox` is set in pt for 183 or 89 mm. Text is 7 pt
  for labels and 6 pt for cell numbers, never below 5 pt. If a number doesn't fit at 5 pt, drop
  it (the shade and tooltip still carry it) rather than shrinking it to 2 pt.
- **3.3** Remove the title, subtitle, prose legend lines, drop shadows, gloss gradients,
  rounded cards and box shadows from static exports. Write the legend text to a
  `*_legend.txt` next to the SVG, for the figure legend.
- **3.4** Colour: GOI in Okabe-Ito vermillion; flanking genes in a muted palette. In the grid,
  the column already identifies the gene, so per-gene colour is redundant (decision below).
- **3.5** Static SVGs use presentation attributes, not CSS classes, so Illustrator, Inkscape
  and PDF converters render them the same as Chrome. Also export PDF via Chrome
  `--print-to-pdf`.

### Phase 4: figure set

- **Main panel:** species tree + anchor grid + location column.
- **Extended Data:** ribbon plot for ≤ 12 representative genomes (`--pub_genomes`, or one per
  genus from the species tree).
- **Interactive HTML only:** anchor positions and gene positions; they are no longer written
  as static SVG.

### Phase 5: verification

Tests for: ≤ 3 lanes; home GOI strand and exon count; fragments not drawn as copies; no text
box overlapping another text box or a cell (metric widths); every text ≥ 5 pt at the declared
print width. Visual pass with headless Chrome at 183 mm on melittin, APYR, DPP4, SP and DCN.

---

## 5. Decisions needed before implementing

1. **Fragment GOI models in figures:** hide them entirely (count only in the tooltip)
   *(recommended)*, or show one faint "N fragment hits" marker.
2. **Anchor-grid colour:** keep per-gene colours, or use neutral grey shaded by identity for
   flanking genes and red for the GOI *(recommended: cleaner and colour-blind-safe; the column
   already names the gene)*.
3. **Exon structure:** drop it from the grid and add a separate gene-model panel
   *(recommended)*, or keep a to-scale mini-model inside the GOI cell.
4. **Invented exons (0.3):** solid arrows in the plot now *(recommended)*, and/or real per-exon
   CDS from the pipeline later.
5. **Main figure:** grid as the main panel with the ribbon as Extended Data *(recommended)*,
   or the other way round.
