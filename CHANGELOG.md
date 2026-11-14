# Changelog

What changed between SynVoy versions. Results of different versions are not comparable in detail, so keep the version and the commit that `./run_synvoy.sh` prints with your results (for example `SynVoy v1.0.1` and `main@1281e47`).

## 1.0.1 (11 November 2026)

Documentation only: working notes that were never meant to be public left the repository. A run computes exactly what it computes with 1.0.0.

## 1.0.0 (11 November 2026)

The first release, commit `b3e3202`.

The version file has said 1.0.0 since 13 October 2026, so a checkout from between 13 October and 11 November also prints `SynVoy v1.0.0`; the commit in the banner tells it apart from the release. These changes of that period matter when you compare with such a run:

- **Preset modules are the default** (10 November). `--modules auto` applies every preset module whose trigger fires. The chooser that picks one of the four presets is used only with `--modules legacy`; a preset you name yourself, with `--preset_override` or a `preset_*` profile, applies in any mode.
- **The tool versions are pinned** (10 November). `environment.yml` pins Nextflow, Python, parasail and the bioinformatics tools (MMseqs2, BLAST, miniprot, Augustus, Prodigal, MAFFT, IQ-TREE, samtools, NCBI datasets) to the tested versions; Java keeps a range, and Entrez Direct, ete3 and numpy have no version. Before, these tools had only lower bounds or none, so two machines could run different tool versions.
- **No silent failures** (30 October to 4 November). A Smith-Waterman, miniprot, Augustus or Prodigal run that does not finish within its time guard stops the task, instead of giving fewer hits or no model. A failed rescue, paralog check or ownership check is named in the report (`warnings`, `failed_steps`) instead of being skipped. The report also records the tool versions and resources of the run, and what each limit of the search dropped (`limits`).
- **One rule for every label** (15 October). HIGH, MEDIUM, AMBIGUOUS and LOW come from one judge: the E-value of the home protein in the place where the gene is expected, whether another home gene explains the model as well, and how much of the query the model covers. AMBIGUOUS means "needs curation" and is never counted. The labels of earlier runs rest on other rules.
- **One folder per home locus** (14 October) for the region files, `regions/<locus>/`. Before, the home loci of a run overwrote each other's region files.

## Before 13 October 2026

Until 13 October 2026 the version file said 2.0.0, so the launcher of that time printed `SynVoy v2.0.0`: such a run is older than 1.0.0, not newer. The checkouts of that time changed often. Two changes of that time matter most when you compare with such a run:

- **Distances in gene units** (11 October). The distances of the search follow the spacing of the flanking genes in each genome, instead of fixed numbers of base pairs. `-profile <execution>,legacy_distances` restores the fixed distances.
- **Figures** (5 to 11 October). Flanking genes are drawn with their real exons (5 October) and the anchor grid has a print version (7 October). From 11 October every searched genome gets a row, also one without a call.
