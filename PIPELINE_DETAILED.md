# SynTerra Detailed Pipeline Description

This document describes the implemented SynTerra pipeline logic in technical detail.
It is intended as the reference for method review, debugging, and scientific validation.

## 1. Purpose and Scope

SynTerra finds orthologous GOI neighborhoods across related genomes by combining:
- protein-to-genome homology search,
- synteny (flanking context),
- iterative database expansion,
- exon/intron-aware model reconstruction.

The workflow is designed to remain functional when genomes have:
- complete annotation,
- partial/inconsistent annotation,
- no annotation.

## 2. Terminology

- GOI: gene of interest (query protein concept).
- Home genome: reference genome where GOI is first located.
- Target genomes: genomes searched iteratively.
- Flanking genes: genes around the GOI locus used as synteny anchors.
- Synteny block: cluster of flanking anchor hits in a target genome.
- Locus: a genomic neighborhood for a parent query where nearby hits belong together.

## 3. High-Level Workflow Graph (main.nf)

1. Resolve query input (`RESOLVE_GENE_INPUT` / `FETCH_QUERY_FROM_ID` / local FASTA).
2. Normalize query to protein (`NORMALIZE_QUERY`).
3. Locate GOI in home genome (`LOCATE_GENE`).
4. Annotate GOI exons in home genome (`ANNOTATE_GOI`).
5. Split distinct home loci (`SPLIT_LOCI`).
6. Build effective home annotation source:
   - provided GFF OR
   - Prodigal home predictions + borrowed annotations.
7. Extract flanking genes (`EXTRACT_FLANKING`).
8. Build initial search DB = flanking + GOI exon/full entries (`PREPARE_INITIAL_DB`).
9. Sort targets by phylogenetic distance (`PHYLO_SORT`).
10. Iterative search on targets (`ITERATIVE_SEARCH`).
11. Cluster candidate regions (`CLUSTER_REGIONS`).
12. Build GOI trees (`COMPUTE_TREE`).
13. Plot integrated synteny (`PLOT_SYNTENY`).
14. Summarize report (`GENERATE_REPORT`).

## 4. Query Handling

### 4.1 Accepted query forms

- UniProt accession
- NCBI protein accession
- local FASTA

`resolve_gene_input.py` resolves ID-based inputs to FASTA and species metadata.

### 4.2 Protein-space normalization

`normalize_query.py` ensures downstream modules use a protein query representation.
This avoids mixing protein-vs-DNA search modes later.

## 5. Home GOI Localization and Exon Annotation

## 5.1 Home localization

`LOCATE_GENE` runs MMseqs and BLAST against home genome and merges hits to BED loci.

## 5.2 GOI exon annotation

`annotate_goi_exons.py` has two major branches:

- GFF-driven mode:
  - parse GFF/GTF features,
  - name/ID/region-based matching to GOI,
  - extract CDS model.

- Hit-driven mode (no usable GFF):
  - infer exon structure from protein-to-genome hits,
  - splice/start/stop-aware refinement,
  - output GOI full sequence plus exon entries.

Output:
- `goi_exons.faa`
- `goi_annotation.bed`
- `goi_info.json`

## 6. Flanking Gene Construction

`extract_flanking_genes.py` extracts flanking genes around each locus from effective home annotation.

When home GFF is unavailable/unusable:
- Prodigal-based fallback predicts proteins in GOI-centered windows.
- Borrowed annotation transfer can supplement home models.

Flanking query normalization (`flanking_query_utils.py`) collapses messy inputs into one parent-level protein per flanking gene concept:
- reconstruct from exon parts if present,
- collapse duplicate fragment records,
- keep strongest representative sequence.

Important:
- Flanking parent IDs are not required to be RefSeq (`XP_...`).
- Any stable ID string is accepted (NCBI, Ensembl, custom labels).

## 7. Initial DB Composition

`PREPARE_INITIAL_DB` combines:
- normalized flanking protein queries,
- GOI full sequence,
- GOI exon sequences,
- fallback fragments only if real exon/tandem evidence is absent.

This DB seeds iterative target-genome search.

## 8. Iterative Target Search (iterative_search_runner.py)

For each target genome (in phylogenetic order):

1. Search current DB vs target genome (MMseqs).
2. Parse hits and identify candidate synteny blocks.
3. Seed block detection with flanking anchors first when available.
4. For each block:
   - run augmented search in padded region (MMseqs + tblastn + optional SW),
   - build GOI candidate models,
   - build flanking candidate models,
   - output region proteins + GFF + homology table.

### 8.1 GOI vs flanking separation

The runner explicitly distinguishes GOI-derived query IDs and flanking-derived IDs.
This prevents flanking post-processing from altering GOI behavior.

### 8.2 Flanking model generation

Flanking candidates are built from:
- miniprot-guided exon chain annotation (`flanking_annotation`), or
- conservative hit-span fallback (`flanking_hits`) when no chain is resolved.

### 8.3 Per-locus flanking dedup (current behavior)

After all blocks for a target genome are processed, SynTerra applies a flanking-only dedup pass:

Function: `deduplicate_flanking_models(...)` in `bin/iterative_search_runner.py`.

Rules:
- Scope: sources `flanking_annotation` and `flanking_hits` only.
- Parent ID selection is attribute-agnostic:
  - tries `SynTerra_Parent`, `ParentProtein`, `protein_id`, `gene_id`, `gene`, `locus_tag`, `Name`, then model ID fallback.
- Models are split into loci by:
  - same chromosome,
  - same strand,
  - coordinate proximity (gap threshold).
- In each locus, keep one best model scored by:
  1. exon-chain strand consistency,
  2. source preference (`flanking_annotation` over `flanking_hits`),
  3. exon count,
  4. total CDS length,
  5. transcript span,
  6. identity score.

Result:
- one best flanking model per parent ID per locus,
- multiple loci for the same parent are retained when biologically distinct.

DB expansion behavior:
- GOI-derived models are used to expand the iterative query DB.
- Flanking models are kept in GFF/region outputs for contextual annotation and plotting, not as next-wave query seeds.

## 9. Region Clustering, Trees, and Plotting

- `CLUSTER_REGIONS` clusters hit regions using synteny score and distance constraints.
- `COMPUTE_TREE` filters to GOI entries and computes a GOI-centric phylogeny.
- `PLOT_SYNTENY` combines home/target beds, GFFs, homology mapping, species map, and tree into HTML plots.
- Plot step filters target genes to candidate BED regions before drawing, reducing off-locus clutter from non-selected regions.

## 10. Annotation Format Robustness

## 10.1 GFF/GTF parsing

`sequence_utils.py` parses:
- GFF3 key=value style,
- GTF key "value" style,
- partial/malformed lines with defensive fallback.

Feature ID retrieval falls back across multiple tags (`ID`, `gene_id`, `transcript_id`, `protein_id`, `Name`, `locus_tag`, etc.) and then synthetic IDs.

## 10.2 No-GFF mode

If GFF is absent, SynTerra still runs through:
- home prediction,
- hit-driven GOI exon inference,
- flanking extraction from fallback annotations,
- target iterative search and model reconstruction.

## 11. Logic Review Summary

The current full-pipeline logic is scientifically coherent for its target use case:

1. Anchor by local synteny, not GOI-only homology.
2. Annotate structure after finding candidate regions.
3. Keep GOI and flanking processing independent.
4. Deduplicate flanking models conservatively per locus.
5. Maintain fallback paths when annotations are missing.

Known practical caveats:

- Any strong homology pipeline can still produce paralog-rich false positives in large gene families.
- Tree quality depends on GOI sequence quality and alignment depth; sparse GOI recovery yields placeholder trees.
- Easy mode quality depends on external metadata/services (NCBI/UniProt availability).

These are method constraints, not single-code-path bugs.

## 12. Recommended Validation Workflow

1. Run a small ground-truth benchmark:

```bash
conda run --no-capture-output -n syntenyfinder \
  python scripts/reproduce_annotation.py \
  --outdir tests/ground_truth_test/output_recheck
```

2. Inspect:
- `predicted_annotations.bed` (GOI)
- `predicted_annotations_all.bed` (gene spans: GOI + flanking)
- `predicted_annotations_all_cds.bed` (CDS entries)

3. Run full Nextflow workflow:

```bash
nextflow run main.nf --gene P01501 --mode easy --outdir results -resume
```

4. Review:
- synteny HTML plot(s),
- GFF outputs in `iterative_results/regions/`,
- report JSON and QC summary.
