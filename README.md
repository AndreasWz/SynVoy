# SynTerra

SynTerra is a Nextflow-based pipeline for the identification of orthologous genes across evolutionary distances using genomic synteny. 

Standard sequence similarity searches can fail to identify orthologs when sequences are highly divergent or consist of short, complex micro-exons. SynTerra attempts to address this by leveraging the conservation of gene order (macro-synteny). By identifying the conserved flanking genes surrounding a Gene of Interest (GOI) in a reference species, the pipeline locates the homologous genomic neighborhood in target species, followed by a localized sequence search to identify the GOI candidate.

The tool is currently in early development.

## License

SynTerra is distributed under the GNU AGPLv3 License.

## Pipeline Overview

The pipeline executes the following primary steps:
1. **Input Resolution:** Processes the input GOI from a UniProt/NCBI identifier or user-provided FASTA.
2. **Genome Staging:** Retrieves the home (reference) genome and target genomes.
3. **Home Locus Identification:** Maps the GOI to the home genome to establish its coordinates.
4. **Context Extraction:** Extracts the sequences of *n* flanking genes upstream and downstream of the GOI.
5. **Macro-Synteny Mapping:** Uses MMseqs2 to map the flanking genes against the target genomes, clustering hits to define candidate syntenic blocks.
6. **Iterative Local Search:** Performs a localized sequence search (using tblastn and miniprot) within the boundaries of the candidate syntenic blocks to identify the GOI.
7. **Phylogenetic Sorting & Visualization:** Sorts the resulting loci by phylogenetic distance and generates an HTML plot representing the syntenic alignments.

## Installation & Requirements

*   **OS:** Linux or macOS
*   **Workflow Manager:** Nextflow (>=22.10.1)
*   **Environment Integration:** Conda/Mamba, Docker, or Singularity

### Conda (Recommended)
```bash
git clone https://github.com/AndreasWz/SynTerra.git
cd SynTerra
conda env create -f environment.yml
conda activate syntenyfinder
```

### Docker
```bash
docker build -t synterra-local:latest .
```

## Execution Modes

SynTerra operates in two modes:

### 1. Easy Mode
Automates the retrieval of genomes and taxonomy metadata. The user provides a GOI identifier, and the pipeline fetches the reference genome and related target assemblies from NCBI.

```bash
nextflow run main.nf \
  --mode easy \
  --query_id Q16553 \
  --max_genomes 5 \
  --outdir results \
  -profile standard
```

### 2. Pro Mode
Designed for offline execution or custom datasets. The user must provide local FASTA files for the query, reference genome, and target genomes.

```bash
nextflow run main.nf \
  --mode pro \
  --query input/query.fasta \
  --home_genome input/reference.fna \
  --target_genomes "input/targets/*.fna" \
  --outdir results \
  -profile standard
```

## Output Structure

The `--outdir` will contain the following primary artifacts:
*   `*_synteny_plot.html`: Visual representation of mapped orthologs.
*   `regions/*.regions.bed`: Genomic coordinates of identified syntenic blocks.
*   `synterra_report.json`: Metrics and parameter summary of the run.
*   `*_tree.nwk`: Newick representation of the phylogenetic tree used for sorting.

## Documentation

For a comprehensive description of the algorithm, execution profiles, and configuration parameters, see the [USAGE.md](USAGE.md) manual.
