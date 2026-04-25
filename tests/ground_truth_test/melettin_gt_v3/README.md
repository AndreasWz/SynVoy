# Melittin GT v3 — regression fixture

This directory pins the known-good output of the melittin ground-truth run
(`results/melettin_gt_v3/`, 2026-03-28) as reference artifacts. Any change
to search logic (`iterative_search_runner.py`, `cluster_grs.py`,
`compute_tree.py`, flanking-gene extraction, scoring) must still reproduce
the coordinates and tree shape recorded here.

## Contents

| File | Source | What it pins |
|---|---|---|
| `regions/*.regions.bed` | `results/melettin_gt_v3/regions/` | Scaffold + coords + score + strand per target species |
| `regions/*.scores.tsv` | `results/melettin_gt_v3/regions/` | Full region-scoring breakdown per target species |
| `locus_1_tree.nwk` | `results/melettin_gt_v3/locus_1_tree.nwk` | Expected tree leaf set (query + 4 species) |
| `flanking_parents.tsv` | `local_data/ground_truth/melettin/melittin_flanking_genes.tsv` | 5 upstream + 5 downstream anchor genes used as parents |

The raw bee genomes (~1-2 GB each) are not committed — only the small
reference artifacts. The test suite ships only these files.

## Expected state (2026-03-28 baseline)

- 5 target species: `Colletes_gigas`, `Euglossa_dilemma`,
  `Melipona_beecheii`, `Tetragonula_carbonaria`, `Xylocopa_violacea`.
- Each species emits exactly one BED row at locus 1.
- Tree has 4 target leaves (Melipona's fallback hit is too weak to enter
  the tree — this is expected, documented here so regressions that
  accidentally *add* a 5th leaf are flagged).
- 10 flanking parents (5 up, 5 down).

## How to use

1. **Canary test** (runs in CI, no genomes needed):
   `pytest tests/test_melettin_gt_v3_fixture.py`
   Asserts the fixture is well-formed and has not been accidentally
   deleted or corrupted.

2. **Regression diff** (after any change to search logic):
   Re-run the melittin GT pipeline end-to-end, then diff against this
   fixture:
   ```
   python3 scripts/validate_melettin_gt_v3.py \
       --outdir results/melettin_gt_v3_rerun
   ```
   The script compares species set, per-species scaffold IDs, coordinate
   deltas, tree leaf set, and best-region confidence strings.

## How to regenerate

If a deliberate change to SynVoy's output format or scoring invalidates
this fixture (e.g. we decide to report a bigger region window, or the
tree-building switches to a method that produces different branch
lengths), regenerate with:

```bash
nextflow run main.nf \
    --mode pro \
    --query local_data/queries/mellettin/GOI_melittin.fasta \
    --home_genome local_data/ground_truth/melettin/home/Apis_mellifera.fa \
    --home_gff local_data/ground_truth/melettin/home/Apis_mellifera.gff \
    --target_genomes "local_data/ground_truth/melettin/targets/*.fa" \
    --outdir results/melettin_gt_v4 \
    --n_flanking_genes 5 \
    --auto_params false --multi_profile false
```

Then:
```bash
cp results/melettin_gt_v4/regions/*.regions.bed tests/ground_truth_test/melettin_gt_v3/regions/
cp results/melettin_gt_v4/regions/*.scores.tsv  tests/ground_truth_test/melettin_gt_v3/regions/
cp results/melettin_gt_v4/locus_1_tree.nwk      tests/ground_truth_test/melettin_gt_v3/
```

Note: only regenerate when the behavior change is intentional. Otherwise,
the pipeline is broken and the fix belongs in the code, not in the
fixture.

## Format notes

- Pre-2026-04-23 fixture captures use bare region names
  (`Reg1_G7_CMEDIUM_S0.47`). Regenerated fixtures will include a
  species prefix (`Colletes_gigas|Reg1_G7_CMEDIUM_S0.47`) because
  `cluster_grs.py` now consumes `species_mapping.tsv`.
- `scripts/validate_melettin_gt_v3.py` does not compare the `name`
  column directly, so both formats are acceptable for regression
  diffing.

## Fixture lineage and known divergences

This fixture was captured from `results/melettin_gt_v3/` on **2026-03-28**.
Pipeline behavior has shifted since then due to (notably) commit
`72bd7d2` "tighten GOI orthology calls and adaptive synteny filtering",
which:

1. **Adds GOI-anchor region injection in `cluster_grs.py`** when the
   iterative search finds GOI hits but no synteny cluster overlaps them.
   Result: species can now emit multiple BED rows per locus, with high
   confidence GOI-anchor regions sorting above the synteny-only region
   the v3 fixture captured. Example: Melipona_beecheii now emits 3
   regions (2 GOI-anchor on `contig_140` + 1 synteny on `contig_94`);
   v3 had only the `contig_94` row. The validator handles this via
   multi-row matching — the v3 region is still found, just no longer
   top-ranked.
2. **Tightens classification thresholds**, which can drop weak
   sequences from the tree-building step. Example: in v5 smoke runs,
   Euglossa's GOI sequence falls below the new tree-inclusion bar
   even though the BED region is still emitted.

Both are intentional changes, not regressions — but the validator
will report (1) as informational notes and (2) as a hard FAIL because
the leaf set genuinely changed. To clear the tree FAIL:

- If you accept the new pipeline behavior, regenerate `locus_1_tree.nwk`
  in this fixture.
- If you suspect over-tightening, loosen the relevant thresholds (see
  `bin/cluster_grs.py` `--classify_*_min_identity` flags or
  `bin/iterative_search_runner.py` GOI classification logic).
