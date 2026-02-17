Findings and Patches (2026-02-14)

Scope
- Reviewed local changes vs current repo state and `changes_done_maybe_bullshit.txt`.
- Focused on: fallback explosion, coordinate mismatch vs ground truth, split GOI/flanking inputs, and pipeline wiring issues.

Patched

1) Smith-Waterman parser/output delimiter bug
- File: `bin/smith_waterman_search.py`
- Fixed literal escaped delimiters in the ssearch path (`\\n` / `\\t`) to real newline/tab handling.
- Result: ssearch output now parses and writes valid m8-style rows.

2) Easy mode local FASTA compatibility
- File: `bin/resolve_gene_input.py`
- Local file input now gets copied to `resolved_query/input_query.fasta`.
- This guarantees `modules/resolve_query.nf` output contract (`resolved_query/*.fasta`) works for file-based queries.

3) Query normalization restored in pipeline
- File: `main.nf`
- Query is now routed through `NORMALIZE_QUERY` before `LOCATE_GENE` and `ANNOTATE_GOI`.
- This restores protein-space behavior for DNA-like query inputs.

4) `max_genomes=0` auto mode preserved
- File: `main.nf`
- Replaced `?:` coercion with explicit null check so value `0` is passed through instead of silently becoming `10`.

5) Iterative search fallback/synteny tightening
- File: `bin/iterative_search_runner.py`
- Added GOI detection helper (`GOI_` / `GOI_copy_`) and seeded synteny blocks from GOI hits when available.
- Split per-parent hits into local loci before annotation to avoid collapsing distant paralog regions.
- Limited fallback output to GOI parents only; non-GOI flanking queries no longer flood raw fallback calls.
- Added compatibility wrapper `identify_best_synteny_block()` to prevent import/test breakage.
- Kept `bits` from parsed hits for downstream candidate ranking.

6) RBH helper robustness and API cleanup
- File: `bin/iterative_search_runner.py`
- `batch_rbh_check()` no longer references undefined global `args`.
- Supports both dict candidates and SeqRecord-like inputs used in tests.
- Normalizes qcov/tcov if reported as fraction vs percent.

7) Species mapping format alignment
- File: `bin/fetch_related_genomes.py`
- Target-species mode now writes 3-column mapping (`accession<TAB>species<TAB>tax_level`) like automatic mode.
- File: `bin/plot_synteny.py`
- Species-map loader already accepts both 2-column and 3-column lines and uses column 2 as display name.

8) Ground-truth harness with split query inputs
- File: `scripts/reproduce_annotation.py`
- Supports either:
  - `--queries` (combined FASTA), or
  - `--goi_queries` + `--flanking_queries` (auto-combined).
- GOI records are prefixed with `GOI_` in the combined FASTA so GOI-specific logic is preserved.
- BED writing now uses BED coordinate convention (0-based start), applies `--offset`, and can override chrom with ground-truth chrom.
- Synteny block seeding in harness now prioritizes GOI-prefixed hits.
- Chrom matching in comparison now normalizes common naming variants (`OV788327` vs `OV788327.1`).

Validation
- Command:
  - `pytest -q tests/test_core_functions.py tests/test_coordinates.py tests/test_rbh_batch.py tests/test_strand.py`
- Result:
  - `31 passed`

Notes
- Ground-truth BED coordinates will still look local/near-zero if `--offset` is not provided for extracted-region runs.
- The broader refactor in `bin/iterative_search_runner.py` remains; this patch pass focused on correctness and output control rather than redesigning the full search strategy.  


UPDATES::

Implemented and documented. I updated `CODEX_FINDINGS.md` with the findings + patch status.

Main fixes applied:

1. Smith-Waterman parser/writer delimiter bug fixed  
- `bin/smith_waterman_search.py`

2. Easy-mode local FASTA now correctly staged into `resolved_query/*.fasta`  
- `bin/resolve_gene_input.py`

3. Query normalization restored into active flow (DNA query -> protein ORF), and `max_genomes=0` no longer gets overwritten  
- `main.nf`

4. Fallback explosion reduction + GOI-focused block seeding + per-locus splitting + compatibility wrapper restored  
- `bin/iterative_search_runner.py`  
- added/used: GOI detection, `split_hits_into_loci`, `identify_best_synteny_block` wrapper, GOI-only fallback behavior

5. Species map format alignment (3-column writer; parser tolerant)  
- `bin/fetch_related_genomes.py`  
- `bin/plot_synteny.py`

6. Ground-truth harness now supports split GOI/flanking input and emits offset-aware BED coordinates  
- `scripts/reproduce_annotation.py`

Validation run:

- `pytest -q tests/test_core_functions.py tests/test_coordinates.py tests/test_rbh_batch.py tests/test_strand.py`
- Result: `31 passed`

Important note for your coordinate issue:
- If benchmarking on extracted region FASTA, you still need `--offset` in `scripts/reproduce_annotation.py` to place predictions into genomic coordinates (otherwise they stay near local region coordinates).


Additional updates (2026-02-14, flanking + pipeline audit)

1) Flanking query normalization before iterative search
- New file: `bin/flanking_query_utils.py`
- Added `collapse_flanking_query_records()`:
  - Reconstructs explicit exon split inputs (`parent|exon_N`) to one parent protein.
  - Collapses repeated-ID fragment inputs to one representative protein per parent.
  - Removes stop symbols (`*`) and non-AA artifacts from query records.

- Wired into:
  - `modules/prepare_initial_db.nf` (pipeline DB build now sees normalized full flanking genes)
  - `scripts/reproduce_annotation.py` (little-test query assembly now matches pipeline behavior)

2) Ground-truth harness output shape for annotation inspection
- `scripts/reproduce_annotation.py` now writes:
  - `predicted_annotations.bed` (GOI CDS predictions used for metrics; unchanged scope)
  - `predicted_annotations_all.bed` (gene-span/mRNA annotations for GOI + flanking)
  - `predicted_annotations_all_cds.bed` (all CDS annotations)

3) Little test rerun (syntenyfinder env)
- Command:
  - `conda run --no-capture-output -n syntenyfinder python scripts/reproduce_annotation.py --outdir tests/ground_truth_test/output_recheck`
- Result metrics (`tests/ground_truth_test/output_recheck/metrics.json`):
  - precision: `0.8888888888888888`
  - recall: `0.42105263157894735`
  - f1: `0.5714285714285714`
- Output sizes:
  - `predicted_annotations.bed`: 18
  - `predicted_annotations_all.bed`: 53
  - `predicted_annotations_all_cds.bed`: 77

4) Pipeline logic audit patch
- File: `modules/compute_tree.nf`
- Fixed redirection bug that could contaminate `goi_only.faa` with stderr text from dedup logging.
  - Old: `> goi_only.faa 2>&1 | head -5`
  - New: `> goi_only.faa 2> goi_dedup.log` + `head -5 goi_dedup.log`

5) Pipeline completion message accuracy
- File: `main.nf`
- Updated final “Key outputs” log lines to match actual generated artifacts (`*_synteny_plot.html`, `*_tree.nwk`, `regions/*.regions.bed`, `intermediate/`).

6) Documentation expansion
- Updated `README.md` and `USAGE.md`:
  - full phase-by-phase workflow
  - flanking normalization behavior
  - rerun commands (little test + big easy-mode run)
  - refreshed wired/reserved parameter documentation

Additional updates (2026-02-14, per-locus flanking + docs overhaul)

1) Flanking dedup switched to per-locus (not global)
- File: `bin/iterative_search_runner.py`
- Added:
  - `_split_models_into_loci(...)`
  - per-locus selection in `deduplicate_flanking_models(...)`
- Keeps one best flanking model per parent ID per genomic locus, using ranking by:
  1. exon-chain consistency
  2. source preference (`flanking_annotation` > `flanking_hits`)
  3. exon count
  4. total CDS length
  5. transcript span
  6. identity

2) Parent-ID handling made format-agnostic
- File: `bin/iterative_search_runner.py`
- Added `_select_parent_id(...)` fallback order:
  - `SynTerra_Parent`, `ParentProtein`, `protein_id`, `gene_id`, `gene`, `locus_tag`, `Name`, then model-ID fallback.
- No dependency on `XP_` naming.

3) Pro-mode query resolution logic fix
- File: `main.nf`
- In pro mode, `--query_id` now uses `RESOLVE_GENE_INPUT` (same robust resolver as easy mode), instead of UniProt-only fetch module.
- This improves support for NCBI IDs and mixed input styles.

4) Little test rerun after per-locus change
- Command:
  - `conda run --no-capture-output -n syntenyfinder python scripts/reproduce_annotation.py --outdir tests/ground_truth_test/output_recheck`
- Log summary:
  - `Flanking model dedup (per-locus): kept 34/35 models across 10 parent IDs (dropped 1 duplicates).`
- Metrics unchanged for GOI benchmark:
  - precision: `0.8888888888888888`
  - recall: `0.42105263157894735`
  - f1: `0.5714285714285714`

5) Documentation rewritten
- Rewritten: `README.md`
- Rewritten: `USAGE.md`
- Added detailed method doc: `PIPELINE_DETAILED.md`

