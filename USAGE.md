# SynTerra Usage Guide

This guide documents the current behavior of the pipeline **as implemented** in the codebase.

**Quickstart**

Easy mode (auto-fetch genomes):
```bash
nextflow run main.nf \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --outdir results/melittin
```

Pro mode (use your own data):
```bash
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome home.fna \
  --target_genomes "targets/*.fna" \
  --mode pro \
  --outdir results/my_run
```

**Inputs**

Query gene
- `--query_id` is a UniProt accession. The pipeline will fetch the protein from UniProt.
- `--gene` is a FASTA file. If it is nucleotide-only, SynTerra translates it to protein (longest ORF across six frames).

Home genome
- `--home_genome` is required in pro mode.
- `--home_species` is required in easy mode. The pipeline fetches the reference genome and GFF (if available) from NCBI.
- `--home_gff` is optional in pro mode. If absent or missing, the pipeline will use Prodigal and borrowed annotations.

Target genomes
- Easy mode: related genomes are downloaded from NCBI.
- Pro mode: use `--target_genomes "path/*.fna"`.

**Phylogenetic sorting**
- Sorting is performed using NCBI taxonomy via `ete3.NCBITaxa`.
- Set `TAXDB` to a local taxdump directory to avoid repeated downloads.
- If taxonomy lookup fails, the pipeline falls back to alphabetical ordering.

**Active parameters**
These are wired to the current pipeline.

| Parameter | Default | Meaning |
|---|---|---|
| `--mode` | `easy` | `easy` or `pro` |
| `--query_id` | null | UniProt accession to fetch |
| `--gene` | null | FASTA query (DNA or protein) |
| `--home_species` | null | Easy mode home species |
| `--home_genome` | null | Pro mode home genome FASTA |
| `--home_gff` | null | Optional home GFF |
| `--target_genomes` | null | Glob of target FASTAs (pro mode) |
| `--target_species` | null | Comma-separated list (easy mode override) |
| `--max_genomes` | 10 | Max related genomes to fetch |
| `--n_flanking_genes` | 10 | Flanking genes per side |
| `--min_flanking_size` | 500 | Min size for flanking genes |
| `--prefer_large_genes` | true | Prefer larger flanking genes |
| `--mmseqs_sensitivity` | 8.5 | MMseqs2 sensitivity |
| `--min_synteny_score` | 0.6 | Synteny score threshold |
| `--outdir` | `results` | Output directory |

**Parameters currently defined but not wired**
These appear in `nextflow.config` but are not used by the current pipeline.
- `cluster_distance`, `min_hit_identity`, `min_hit_length`
- `enable_smith_waterman`, `sw_method`, `sw_min_score`, `sw_min_identity`
- `region_padding`, `enable_splice_variants`, `enable_frameshifts`, `mutation_rate`, `num_mutant_variants`
- `expand_db_threshold`, `diamond_sensitivity`, `min_gene_identity`, `augustus_species`

If you want these activated, ask for a wiring pass.

**Outputs**
Primary outputs in `--outdir`:
- `*_synteny_plot.html` interactive visualization
- `*_tree.html` (optional) and `*.nwk` GOI phylogeny
- `synterra_report.json` summary report
- `qc/genome_qc_summary.json` assembly QC summary
- `regions/*.regions.bed` synteny region calls

Intermediate outputs in `--outdir/intermediate`:
- `locate_gene/` merged GOI hits and raw hit files
- `annotate_goi/` GOI exon FASTA, BED, JSON
- `split_loci/` `locus_*.bed`
- `flanking/` `synteny_block_*.bed`, `flanking_proteins_*.faa`
- `initial_db/` `initial_db_*.faa`
- `phylo_sort/` `sorted_genomes.txt`
- `query/` normalized query FASTA

The Nextflow `work/` directory is still used for execution and caching.

**Troubleshooting**
- No hits found: increase `--mmseqs_sensitivity` or reduce `--min_synteny_score`.
- Easy mode fails early: ensure `ncbi-datasets-cli` and Entrez Direct are installed.
- Phylo sort falls back to alphabetical: set `TAXDB` to a valid taxdump directory.

**Performance tips**
- Use `-profile docker` or `-profile singularity` for reproducible toolchains.
- For HPC, use the provided `slurm_submit.sh` template or the `hpc_*` profiles.

