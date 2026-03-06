# Plot Review: Recent SynTerra Runs (2026-03-06)

## Scope
Reviewed the most recent synteny/tree plot outputs:

- `local_runs/results_mrjp3/synteny_block_locus_1_synteny_plot.html`
- `local_runs/results_mrjp3/synteny_block_locus_1_tree.html`
- `local_runs/results_melittin_koludarov_no_gff/synteny_block_locus_1_synteny_plot_v2.html`
- `local_runs/results_melittin_koludarov_no_gff/synteny_block_locus_1_synteny_plot_v2_tree.html`
- `local_runs/results_human_ly6e_tuned/synteny_plot_v3.html`
- `local_runs/results_human_ly6e_tuned/synteny_plot_v3_tree.html`

## Evidence Snapshot

| Run/file | Tracks (synteny) | Legend entries | Tree leaf labels | Canvas width |
|---|---:|---:|---:|---:|
| `results_mrjp3` | 5 | 39 | 5 | Synteny: 1500, Tree: 900 |
| `results_melittin_koludarov_no_gff` | 13 | 26 | 13 | Synteny: 1500, Tree: 900 |
| `results_human_ly6e_tuned` | 16 | 42 | 46 | Synteny: 1500, Tree: 900 |

Additional observations from parsed plot payloads:

- In dense synteny plots, many labels are rotated (`~121-122` rotated labels).
- In dense synteny plots, most genes are visually tiny after compression:
  - `results_melittin_koludarov_no_gff`: `51.3%` of genes have visual width `<1.5 kb` equivalent.
  - `results_human_ly6e_tuned`: `46.5%` of genes have visual width `<1.5 kb` equivalent.
- `results_human_ly6e_tuned` includes at least one track explicitly marked `GOI absent`, but still contributes to visual clutter.

## Main Issues

1. **Readability is too low in dense plots**
- Gene labels are very small and highly rotated in crowded regions.
- Tree labels are long and become hard to scan when many leaves are present (notably 46 leaves in `synteny_plot_v3_tree.html`).

2. **Signal-to-noise ratio is low**
- Ribbon crossings dominate in dense blocks, making ortholog paths hard to follow.
- Legends become oversized relative to track count (example: 39 legend entries for 5 tracks).

3. **Layout does not scale with data size**
- Synteny width is fixed to `args.plot_width` (default `1500`).
- Tree width is hard-coded to `900`, which leaves little horizontal room for long labels.

## Suggested Fixes (Prioritized)

### P0 (high impact, low risk)
1. **Adaptive font/label policy**
- Increase minimum label font for key labels (GOI, track labels, tree labels).
- For dense tracks, show text only for GOI + top-scoring homologs; keep others hover-only.
- Keep rotated labels only above a stricter width threshold, otherwise hide text.

2. **Adaptive canvas sizing**
- Make tree width depend on max leaf-label length and leaf count (instead of fixed `900`).
- Make synteny width scale with number of tracks/genes (or provide wider default for desktop).

3. **Legend reduction**
- Limit legend to most frequent/important groups plus GOI.
- Optionally add `--max_legend_entries` and hide low-frequency entries by default.

### P1 (quality improvements)
1. **Ribbon clutter controls**
- Add thresholding (for example by identity score or top-N links per target gene).
- Lower default ribbon alpha in dense plots and optionally draw only GOI-connected ribbons first.

2. **Tree label cleanup**
- Replace long fallback IDs in display labels with compact species/accession labels.
- Collapse duplicate GOI copies per genome in default tree mode (or provide toggle).

3. **Track filtering**
- Optional flag to hide tracks with GOI absent or very low informative content.
- Keep those tracks in a separate summary table instead of main panel.

### P2 (usability polish)
1. **Separate “overview” and “detail” outputs**
- Overview plot: fewer labels/ribbons, cleaner comparative view.
- Detail plot: full labels for manual inspection.

2. **Add explicit density diagnostics into subtitle**
- Example: tracks, genes, ribbons, hidden labels, and filtering thresholds used.

## Code Touchpoints

Main places to implement these changes in `bin/plot_synteny.py`:

- Tree layout sizing: `_render_tree_html()` (`width=900`, margins, tree text font).
- Synteny width default arg: `--plot_width` (default `1500`).
- Synteny layout and legend fonts: `fig.update_layout(...)` for legend/font/margins.
- Label density logic: section `# -- 5d. Gene labels`.
- Ribbon rendering and alpha: section `# -- 5b. Ribbons`.

## Recommended First Patch Set

1. Add CLI knobs:
- `--max_legend_entries` (default `25`)
- `--min_label_gene_width` (default stricter than current)
- `--ribbon_alpha_dense` and `--max_ribbons_per_track_pair`
- `--tree_width_auto` (on by default)

2. Update defaults:
- Larger tree label font and dynamic tree width.
- Smaller default ribbon alpha in dense plots.
- Hide non-priority labels when `n_tracks` exceeds a threshold.

3. Keep backward compatibility:
- Preserve old behavior behind explicit flags where needed.
