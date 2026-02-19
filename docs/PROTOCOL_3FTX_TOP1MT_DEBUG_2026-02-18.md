# SynTerra Debug Protocol: 3FTx/TOP1MT Context (2026-02-18)

## Goal
Audit why the snake 3FTx synteny plots did not recover expected `TOP1MT`/Ly6 context, verify annotation usage, and implement only fixes that improve biological interpretability.

## Inputs Audited
- Run directory: `results_3snake_3ftx_v2`
- New plots:
  - `results_3snake_3ftx_v2/synteny_block_locus_1_synteny_plot.html`
  - `results_3snake_3ftx_v2/synteny_block_locus_2_synteny_plot.html`
- Ground-truth reference: `supp_mat_1_nature_paper.pdf`
- Pipeline code paths:
  - `bin/cluster_grs.py`
  - `bin/plot_synteny.py`
  - `bin/fetch_related_genomes.py`
  - `modules/cluster_regions.nf`
  - `main.nf`

## Diagnostic Steps
1. Verified run completion and extracted cluster/plot inputs.
2. Checked whether target native annotations were propagated into plot-input GFFs.
3. Compared candidate BED regions against GOI coordinates from iterative target GFFs.
4. Cross-checked target/home GFF gene names for `TOP1MT`, `TOP1`, `LY6*`, `3FTx`.
5. Mapped root causes to pipeline stage (download, clustering, plotting).

## Findings
1. Annotation usage is asymmetric:
- Ophiophagus (`GCA_000516915.1`) had native annotation labels (`TargetGene`, `TargetProduct`) in plot-input GFF.
- Bungarus (`GCA_023653725.1`) had no native target labels because no target GFF was downloaded for that accession in easy mode.

2. Candidate-region selection was dropping true GOI loci:
- Region clustering was based on raw MMseqs hit density, not final iterative GOI models.
- Result: top candidate regions often did **not** overlap GOI annotations in final target GFF.
- Consequence: plot filtering removed GOI context and produced unusable multi-contig clutter.

3. `TOP1MT` naming mismatch in available annotations:
- In audited Ophiophagus annotation, nearby topoisomerase signal is labeled `TOP1`, not `TOP1MT`.
- Home Naja annotation is largely locus-tag based; flanking IDs are mostly `gene-E2320_*`, limiting direct biological label visibility.

## Fixes Implemented
1. GOI-aware region prioritization in clustering:
- File: `bin/cluster_grs.py`
- Added `--target_gff` input and GOI interval parsing from iterative GFF.
- Clusters overlapping GOI are prioritized.
- If no scored cluster overlaps GOI, GOI-anchor regions are injected so GOI loci cannot be dropped.

2. Wired iterative target GFF into region clustering:
- Files:
  - `main.nf`
  - `modules/cluster_regions.nf`
- Added locus+genome join to pass the correct target GFF to each `CLUSTER_REGIONS` task.

3. Plot made robust against bad candidate BEDs:
- File: `bin/plot_synteny.py`
- If candidate BED has GOI-overlapping regions, only those are used.
- If candidate BED misses GOI, fallback to GOI-centered context window from target GFF.
- If GOI is present, target track is restricted to GOI chromosome to reduce clutter.

4. Improved annotation retrieval fallback for GenBank targets:
- File: `bin/fetch_related_genomes.py`
- If `GCA_*` download lacks usable GFF, tool now tries `GCF_*` counterpart for annotation fallback.

## Validation Performed
1. Syntax checks:
- `python3 -m py_compile bin/cluster_grs.py bin/plot_synteny.py bin/fetch_related_genomes.py`

2. Functional sanity (manual reruns of clustering on existing run artifacts):
- New candidate beds now include GOI-overlapping regions for previously failing Bungarus and Ophiophagus cases.
- Example outputs created in `/tmp/l1_bung_new.bed`, `/tmp/l2_bung_new.bed`, `/tmp/l2_ophi_new.bed` show GOI-first regions.

## Residual Risks
1. If upstream assemblies truly lack informative annotations, `TOP1MT` may still appear as `TOP1`/locus-tag despite correct locus recovery.
2. Easy-mode species selection can still pick poor GenBank assemblies; new GCF fallback improves annotation retrieval but does not guarantee ideal assembly choice.
3. Full end-to-end rerun in your Nextflow runtime is still needed to regenerate final HTML with these code-level fixes.
