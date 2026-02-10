# SynTerra

**Phylogenetically informed syntenic ortholog discovery**

SynTerra is a Nextflow pipeline that finds orthologs of a Gene of Interest (GOI) across related genomes. It is designed to work even when target genomes are unannotated by combining synteny (gene order conservation), iterative search, and exon-aware annotation.

**What it does**
- Locates the GOI in a home genome using MMseqs2 and BLAST.
- Annotates GOI exons (from GFF if available, otherwise from hits and splice-site logic).
- Extracts flanking genes around the GOI to define a synteny block.
- Orders target genomes by phylogenetic distance and searches them iteratively.
- Clusters syntenic regions, builds a GOI phylogeny, and generates an interactive synteny plot.

**How it works (high level)**
1. Validate inputs and normalize the query. Nucleotide queries are translated to protein.
2. Locate the GOI in the home genome and annotate its exon structure.
3. Extract flanking genes around each candidate locus.
4. Sort target genomes by taxonomy and run iterative search (MMseqs2 + optional Smith-Waterman).
5. Cluster regions by synteny, infer a GOI tree, and render the plot and report.

**Requirements**
- Nextflow >= 22.10.1
- Conda, Docker, or Singularity
- External tools (via `environment.yml`): mmseqs2, BLAST, prodigal, mafft, fasttree, ete3
- Easy mode only: `ncbi-datasets-cli` and Entrez Direct (`esearch`, `efetch`, `xtract`)

**Quickstart**

Easy mode (fetches genomes from NCBI):
```bash
nextflow run main.nf \
  --query_id P01501 \
  --home_species "Apis mellifera" \
  --outdir results/melittin
```

Pro mode (use your own files):
```bash
nextflow run main.nf \
  --gene my_gene.fasta \
  --home_genome home.fna \
  --target_genomes "targets/*.fna" \
  --mode pro \
  --outdir results/my_run
```

**Key outputs (in `--outdir`)**
- `*_synteny_plot.html` interactive synteny visualization
- `*_tree.html` (optional) and `*.nwk` GOI phylogeny
- `synterra_report.json` summary report
- `qc/genome_qc_summary.json` assembly quality summary
- `regions/*.regions.bed` synteny region calls
- `intermediate/` exported intermediate artifacts (loci, flanking blocks, etc.)

**Notes**
- If you provide a nucleotide query, SynTerra will translate it to protein and use the longest ORF across six frames.
- If multiple loci are detected, SynTerra will generate one synteny plot per locus.
- Phylogenetic sorting uses NCBI taxonomy. Set `TAXDB` to a local taxdump directory for best results.

**Docs**
- Full usage and parameter reference: `USAGE.md`
- Example configs: `conf/`
- Tests: `tests/`

If you want help or changes, open an issue or ask for a targeted tweak.
