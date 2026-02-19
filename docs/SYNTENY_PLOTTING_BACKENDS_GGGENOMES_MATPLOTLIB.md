# SynTerra Synteny Plotting: gggenomes vs Matplotlib vs Current Plotly

## 1. Purpose

This document evaluates how SynTerra could support publication-ready synteny figures using:
- `gggenomes` (R)
- `matplotlib` (Python)

and compares both against the current Plotly implementation.

The goal is not to replace interactive plotting blindly. The goal is a practical plotting stack where:
- interactive plots remain useful for QC and debugging,
- static figures are high quality and reproducible for manuscripts,
- implementation effort stays realistic.

## 2. Current SynTerra Plotting State

### 2.1 Current pipeline wiring

`main.nf` groups locus-specific plotting inputs and calls `PLOT_SYNTENY` (`main.nf:468`, `main.nf:495`).

`modules/plot_synteny.nf` currently invokes `plot_synteny.py` and emits:
- `*_synteny_plot.html`
- `*_tree.html` (optional)
- staged input bundle `plot_inputs_*`

(`modules/plot_synteny.nf:16-18`, `modules/plot_synteny.nf:58-68`).

`bin/plot_synteny.py` currently:
- parses home BED, target GFF, candidate BED, homology TSV, and tree,
- filters target genes to candidate regions,
- draws arrows/ribbons with Plotly,
- writes an HTML plot.

### 2.2 Practical limitation observed in recent runs

In your recent run artifacts, the plotting process executed with all target files present, but output showed only home genome tracks (`Tracks: 1`).

Evidence:
- `work/7f/86e26b454ec07c30f6cd7190f57d67/.command.out`
- `work/eb/8a710b6da0534adcabf958015a9fb6/.command.out`

Root behavior in code:
- genes are filtered by candidate BED overlap (`bin/plot_synteny.py:861`),
- empty target tracks are dropped (`bin/plot_synteny.py:865-866`).

Important implication: changing plotting library alone will not fix missing tracks. Some data-contract and filtering diagnostics need to be improved in parallel.

## 3. Decision Criteria

For SynTerra, backend choice should be based on:
- publication quality (vector export, typographic control),
- compatibility with current inputs (BED/GFF/links/tree),
- effort to integrate in Nextflow profiles (conda/docker/singularity),
- scalability for many genomes/tracks,
- maintainability by your team.

## 4. gggenomes Path (R)

## 4.1 Why gggenomes is attractive for SynTerra

`gggenomes` is purpose-built for comparative genomics and synteny-style maps.

Key strengths from docs:
- Multi-track genomic plotting model (`seqs`, `genes/feats`, `links`) rather than forcing one table.
- Native geoms for sequence tracks, gene arrows, links, labels, and sequence break decorators.
- Built-in manipulation verbs to align, focus, flip/sync orientation, and pick subsets.

Relevant references:
- package site and installation guidance
- track-oriented API (`gggenomes()`, `geom_seq()`, `geom_gene()`, `geom_link()`)
- manipulation verbs (`focus()`, `sync()`, `align()`)

## 4.2 Data model fit with SynTerra outputs

SynTerra already emits almost everything gggenomes needs.

### Mapping proposal

| SynTerra output | gggenomes track | Notes |
|---|---|---|
| `home_bed`, target-derived seq spans | `seqs` | one row per sequence context with `seq_id,length` (+ optional `bin_id,start,end,strand`) |
| target/home gene models (GFF/BED-derived) | `genes` (or `feats`) | include `seq_id,start,end,strand,gene_id,label,home_gene_id,is_goi` |
| homology mapping + overlaps | `links` | include `seq_id,start,end,seq_id2,start2,end2,link_type,score` |
| candidate regions | either `feats` track or pre-filter logic | can be shown as context ribbons/boxes |
| tree | separate panel (`ggtree`) or ordering input | `ggtree` is optional per gggenomes install docs |

### Coordinate caveat

- `read_bed()` in gggenomes converts BED 0-based starts to 1-based starts.
- SynTerra must keep coordinate conversion explicit and tested to avoid off-by-one errors.

## 4.3 What you can do in gggenomes that maps directly to current needs

- `geom_seq()` for tracks,
- `geom_gene()` for gene arrows,
- `geom_link()` for synteny/homology links,
- `focus()` to zoom loci around GOI-like features,
- `geom_seq_break()` for truncated/zoomed sequence decorations,
- `sync()` to auto-orient bins by link support (helps reduce crossing links).

These are close to the semantics you currently hand-code in Plotly (track layout, arrows, ribbons, gap markers), but implemented in a genomics-specific plotting grammar.

## 4.4 Integration design for SynTerra

### Minimal architecture

1. Add a preprocessing step (Python or R) that writes normalized TSV tracks per locus:
- `plot_tracks/locus_X.seqs.tsv`
- `plot_tracks/locus_X.genes.tsv`
- `plot_tracks/locus_X.links.tsv`
- optional `plot_tracks/locus_X.candidates.tsv`

2. Add new plotting module:
- `modules/plot_synteny_gggenomes.nf`
- calls `scripts/plot_synteny_gggenomes.R`

3. Render static outputs per locus:
- `*_synteny_plot.gggenomes.pdf`
- `*_synteny_plot.gggenomes.svg`
- optional PNG

### Suggested R rendering skeleton

```r
library(gggenomes)
library(ggplot2)
library(readr)

seqs  <- read_tsv("locus_1.seqs.tsv", show_col_types = FALSE)
genes <- read_tsv("locus_1.genes.tsv", show_col_types = FALSE)
links <- read_tsv("locus_1.links.tsv", show_col_types = FALSE)

p <- gggenomes(seqs = seqs, genes = genes, links = links, adjacent_only = TRUE) +
  geom_seq() +
  geom_gene(aes(fill = home_gene_id), position = "strand") +
  geom_link(aes(fill = home_gene_id), alpha = 0.25) +
  geom_seq_label() +
  theme_gggenomes_clean()

# optional locus focusing if loci/GOI metadata exists
# p <- p %>% focus(is_goi, .track_id = "genes", .expand = 5000)

ggsave("synteny_block_locus_1.gggenomes.pdf", p, width = 12, height = 6, units = "in")
ggsave("synteny_block_locus_1.gggenomes.svg", p, width = 12, height = 6, units = "in")
```

## 4.5 Pros and cons of gggenomes for SynTerra

### Pros

- Best domain fit for comparative genomic maps.
- Better static figure defaults than hand-rolled Matplotlib for this problem class.
- Strong composability with ggplot ecosystem.
- Natural path to journal-friendly PDF/SVG figures.

### Cons

- Introduces an R runtime and package management burden.
- Tree side-by-side plotting usually needs `ggtree` integration.
- Another language in the pipeline increases maintenance surface.

## 4.6 Dependency and reproducibility strategy

Recommended in pipeline context:
- Do not install R packages at runtime from the internet.
- Use either:
  - a dedicated conda environment file for gg backend, or
  - a pinned container (Docker/Singularity) with R + gggenomes preinstalled.

This is critical for reproducible Nextflow runs, especially on HPC.

## 5. Matplotlib Path (Python)

## 5.1 Why Matplotlib is still useful

Matplotlib is already in your Python ecosystem and gives full low-level control over static export.

From official docs:
- `savefig()` supports explicit format, dpi, metadata, and backend control.
- Static backends produce vector outputs (`pdf`, `svg`, `ps`, `pgf`) and raster (`png`).
- style sheets and `rcParams` allow publication presets.

## 5.2 Fit with current SynTerra code

You already encode layout logic in Python (`plot_synteny.py`):
- track spacing,
- gene polygon coordinates,
- ribbons,
- labels,
- gap compression.

This can be ported incrementally to Matplotlib without data-model redesign.

### Two Matplotlib implementation levels

1. Quick migration
- Keep current data parsing and layout.
- Replace Plotly trace construction with Matplotlib patches.
- Output static `pdf/svg/png`.

2. Performance-optimized migration
- Batch gene polygons with `PolyCollection`.
- Batch ribbons/links with collections.
- Reduce artist count for large genome panels.

## 5.3 Publication-quality control in Matplotlib

Useful controls:
- style sheets / `rcParams` for consistent typography and line weights,
- explicit vector export via `savefig(..., format='pdf'|'svg')`,
- optional raster export for slides.

You can maintain two presets:
- `synterra_interactive_like.mplstyle` (fast QC)
- `synterra_publication.mplstyle` (journal style)

## 5.4 Pros and cons of Matplotlib for SynTerra

### Pros

- No new language stack.
- Full control over rendering and export.
- Easy integration with existing Python scripts/tests.

### Cons

- More custom code to maintain for genomics-specific operations that gggenomes already provides.
- Higher effort to match ggplot-level aesthetics out of the box.
- Harder to replicate gggenomes conveniences like focus/sync grammar without extra code.

## 6. Comparison Against Current Plotly Method

| Criterion | Current Plotly | gggenomes | Matplotlib |
|---|---|---|---|
| Interactive exploration | Strong | Weak/native static | Weak/native static |
| Publication static output | Moderate (HTML-first, can export but not ideal workflow) | Strong (ggplot workflow to PDF/SVG) | Strong (PDF/SVG via savefig) |
| Domain-specific genomics grammar | Custom/manual | Strong built-in | Custom/manual |
| Integration effort from current code | Already done | Medium-high | Medium |
| Long-term maintainability | Medium (custom logic) | Medium (R dependency) | Medium-high (more custom plotting code) |
| Tree + synteny layout options | Custom | Good with optional ggtree integration | Custom |

## 7. Recommended Strategy

## 7.1 Short answer

Use a hybrid strategy:
- keep Plotly for interactive QC,
- add static publication backend,
- prioritize gggenomes first for static comparative maps,
- keep Matplotlib as a fallback or for custom figure variants.

## 7.2 Why this is the best fit

- Your users need publication-ready figures now.
- gggenomes gives faster path to high-quality static comparative plots.
- Plotly remains valuable during method tuning and troubleshooting.
- Matplotlib remains useful where you need strict custom rendering control or avoid R.

## 8. Concrete SynTerra Implementation Plan

## Phase 0 (must do first, independent of backend)

Improve plotting data diagnostics:
- print per-target counts before/after candidate filtering,
- warn when all targets are removed,
- optionally fall back to unfiltered tracks if everything is filtered away.

This addresses current "home-only" silent failure mode.

## Phase 1 (backend-agnostic plotting contract)

Introduce canonical plotting tables per locus:
- `seqs.tsv`, `genes.tsv`, `links.tsv`, `meta.tsv`

and keep current Plotly renderer consuming these tables.

Benefits:
- one source of truth,
- backend swapping becomes straightforward,
- easier regression tests.

## Phase 2 (gggenomes backend)

Add:
- `scripts/plot_synteny_gggenomes.R`
- `modules/plot_synteny_gggenomes.nf`
- parameters:
  - `--plot_backend plotly|gggenomes|matplotlib|all` (default `plotly`)
  - `--plot_formats html,pdf,svg,png` (backend-dependent)

## Phase 3 (matplotlib backend)

Add:
- `scripts/plot_synteny_mpl.py`
- publication style presets (`.mplstyle`)
- optional high-performance collections mode.

## 9. Suggested CLI and Output Contract

## New user-facing parameters

- `--plot_backend` default `plotly`
- `--plot_static` default `false`
- `--plot_formats` default `pdf,svg`
- `--plot_style` default `synterra`

## Output naming

- `*_synteny_plot.html` (Plotly)
- `*_synteny_plot.gggenomes.pdf`
- `*_synteny_plot.gggenomes.svg`
- `*_synteny_plot.mpl.pdf`
- `*_synteny_plot.mpl.svg`

## 10. Effort Estimate

- Phase 0 diagnostics/fallbacks: 0.5-1 day
- Phase 1 canonical plotting tables: 1-2 days
- Phase 2 gggenomes backend: 2-4 days
- Phase 3 matplotlib backend: 2-4 days
- integration tests + docs: 1-2 days

Total practical path:
- gggenomes-first static support in about 1 week,
- full dual static backends in about 2 weeks.

## 11. Validation Checklist

- Compare one known-good locus between Plotly and gggenomes:
  - same genomes present,
  - same GOI marking,
  - same homolog links.
- Compare one difficult locus with sparse targets.
- Ensure PDF/SVG pass visual QC in Illustrator/Inkscape.
- Ensure no internet required at run-time for plotting.
- Ensure `-resume` works with backend switch.

## 12. Recommendation for your project now

1. Keep Plotly as default interactive backend.
2. Implement backend-agnostic plot tables.
3. Add gggenomes static backend first.
4. Add Matplotlib backend only if you need fully Python-native publication control or special custom layouts not easy in gggenomes.

This gives the best quality-to-effort ratio for SynTerra.

## References

### SynTerra local implementation
- `main.nf`
- `modules/plot_synteny.nf`
- `bin/plot_synteny.py`
- run artifacts in `work/*/.command.out` (track-count observations)

### gggenomes documentation
- https://thackl.github.io/gggenomes/
- https://thackl.github.io/gggenomes/reference/gggenomes.html
- https://thackl.github.io/gggenomes/reference/index.html
- https://thackl.github.io/gggenomes/reference/read_gff3.html
- https://thackl.github.io/gggenomes/reference/read_bed.html
- https://thackl.github.io/gggenomes/reference/focus.html
- https://thackl.github.io/gggenomes/reference/flip.html
- https://thackl.github.io/gggenomes/reference/geom_seq_break.html
- https://rdrr.io/cran/gggenomes/

### Matplotlib documentation
- https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html
- https://matplotlib.org/stable/users/explain/figure/backends.html
- https://matplotlib.org/stable/users/explain/customizing.html
- https://matplotlib.org/stable/users/explain/artists/performance.html
- https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.PathPatch.html
- https://matplotlib.org/stable/api/collections_api.html
