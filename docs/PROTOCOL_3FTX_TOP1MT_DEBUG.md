# SynTerra Debug Protocol: 3FTx-Ly6-TOP1 Context

## Goal
Ensure SynTerra uses real target annotations (not only synthetic IDs) so Ly6/TOP1-family context can be evaluated for 3FTx loci.

## Scope
- Inputs: `results_3snake_3ftx` run outputs (Naja home, Bungarus + Ophiophagus targets)
- Code touchpoints:
  - `bin/iterative_search_runner.py`
  - `bin/plot_synteny.py`
  - `bin/fetch_related_genomes.py`

## Protocol
1. Verify what annotations are actually available per target genome.
- Command:
```bash
find results_3snake_3ftx/downloaded_genomes/easy_mode_genomes -maxdepth 1 -type f | rg '\\.(fna|gff|gff3)(\\.gz)?$'
```
- Expected: each `.fna` has a matching `.gff/.gff3` when NCBI provides one.

2. Verify the old failure mode (synthetic labels dominating).
- Command:
```bash
rg -n 'SynTerra_Parent=|TargetGene=|TargetProduct=' results_3snake_3ftx/plot_inputs_synteny_block_locus_1/*.gff
```
- Old behavior: `SynTerra_Parent=gene-E2320_*` present, `TargetGene/TargetProduct` absent.

3. Verify native annotation signal exists in target reference GFF.
- Command:
```bash
rg -n 'LY6E|Lypd2|LYPD|TOP1|TOP1MT' results_3snake_3ftx/downloaded_genomes/easy_mode_genomes/GCA_000516915.1.gff
```
- Evidence already observed:
  - `AZIM01004932.1` has `LY6E` and `Lypd2`
  - `TOP1` exists on another contig (`AZIM01000244.1`)

4. Patch annotation propagation at model emission time.
- `iterative_search_runner.py` now:
  - finds adjacent native annotation files (`.gff/.gff3` and gz variants)
  - indexes native gene/transcript features
  - attaches `TargetGene`, `TargetProduct`, `TargetID` to emitted mRNA features when overlap/nearby matches exist

5. Patch plot label priority.
- `plot_synteny.py` now:
  - parses `TargetGene` and `TargetProduct`
  - uses informative target labels for non-home tracks
  - keeps `SynTerra_Parent` for homology grouping/ribbons

6. Patch downloader annotation discovery.
- `fetch_related_genomes.py` now scans:
  - `*.gff`, `*.gff3`, `*.gff.gz`, `*.gff3.gz`
- Compressed GFFs are decompressed to `<accession>.gff` for downstream use.

7. Re-run pipeline and validate improved biology-facing output.
- After rerun, check:
```bash
rg -n 'TargetGene=|TargetProduct=' results_3snake_3ftx/regions/*.gff
rg -n 'LY6E|Lypd2|TOP1|TOP1MT' results_3snake_3ftx/synteny_block_locus_*_synteny_plot.html
```
- Pass condition: target plots/hover and GFF carry real target labels where native annotation exists.

## Known Biological Limitation
- In this dataset, explicit `TOP1MT` annotation was not observed in the downloaded Ophiophagus GFF; `TOP1` was observed.
- Therefore, absence of `TOP1MT` label in a specific locus can be annotation-source limitation, not only pipeline logic.
