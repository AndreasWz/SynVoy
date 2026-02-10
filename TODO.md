# SynTerra TODO

## High Priority — Bugs / Broken
- [x] Fix GOIs missing in home genome track for Tetramorium run (multiple hits across loci — GOI re-injection works for single-hit locus but needs validation for multi-hit cases)
- [x] Fix borrowed flanking genes integration — ensure annotations from annotated target genomes are correctly propagated to the synteny plot
- [x] Account for evolutionary events: detect tandem duplications / tandem arrays and represent them appropriately in the analysis and visualization
- [ ] Prodigal on eukaryotic genomes: using metagenomic mode (`-p meta`) on eukaryotes is suboptimal — gene predictions miss small genes and multi-exon genes. Consider miniprot or GeMoMa as alternatives for home genome annotation when no GFF is available
- [x] Remove redundant `.first()` on `best_locus_id_ch` in main.nf — it's already a value channel from `.toSortedList().map()`, causes harmless WARN on every run
- [ ] Pipeline phase banners print multiple times in terminal output due to Nextflow's ANSI refresh — cosmetic but confusing

## Medium Priority — Improvements
- [ ] Show multiple loci in the plot — currently only the best locus (by e-value) is visualized; user might want to see all loci or pick one
- [ ] Better gene names for unannotated genomes — Prodigal names like `pred_OV788327.1_22401957` are not informative; could add functional annotation via InterProScan, eggNOG-mapper, or at least BLAST best-hit descriptions
- [ ] Gene inversions in synteny plot — currently genes on opposite strands are shown but inversions of syntenic blocks are not explicitly highlighted
- [ ] Support user-provided local target genomes — allow `--target_genomes /path/to/genome.fna` in addition to species names
- [ ] Improve tree integration — tree currently only includes GOIs from best locus; should ideally represent all loci or let user choose
- [ ] Add summary statistics panel to HTML plot — number of GOI copies per species, synteny scores, confidence metrics
- [ ] Graceful handling when no synteny is found — currently produces empty/confusing plots; should produce a clear "no syntenic regions found" message
- [ ] Make color scheme colorblind-friendly — current palette may not be distinguishable for all users

## Low Priority — Nice to Have
- [ ] Support multiple query genes simultaneously — run pipeline for a set of genes and produce combined report
- [ ] PDF/SVG export quality — Plotly's static export is mediocre; consider matplotlib/cairo backend for publication-quality figures
- [ ] Add CI/CD with automated test runs (melittin + tetramorium test cases)
- [ ] Clean up archive/ directory — contains old debugging scripts and test files
- [ ] Docker/Singularity container — ensure all dependencies (mmseqs2, prodigal, mafft, fasttree, miniprot) are in the container image
- [ ] Support non-insect taxa — verify fetch_related_genomes taxonomy walking works beyond Hymenoptera
- [ ] Configurable Prodigal top-N% filter — currently hardcoded to top 10% longest; should be a parameter or adaptive
- [ ] Cache downloaded genomes across runs — avoid re-downloading the same NCBI genomes for different queries on the same species




fix that if taxid valid but no genome found, no error happens, but instead it should just skip that species and continue with the next one, even if home species
