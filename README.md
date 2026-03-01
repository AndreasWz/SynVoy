<div align="center">
  <img src="assets/synterra_logo.png" alt="SynTerra Logo" width="300" onerror="this.src='https://via.placeholder.com/300x100?text=SynTerra'">
  <h1>SynTerra</h1>
  <p><strong>Synteny-Guided Evolutionary Ortholog Discovery & Visualization</strong></p>
  
  [![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://www.nextflow.io/)
  [![Conda](https://img.shields.io/badge/conda-supported-blue.svg)](https://docs.conda.io/en/latest/)
  [![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
</div>

---

**SynTerra** is an incredibly robust, massively scalable Nextflow pipeline designed to relentlessly hunt down heavily diverged orthologous genes across evolutionary vast distances. 

When standard sequence similarity searches (like BLAST) fail due to extreme sequence divergence, rapid mutation rates, or complex micro-exons, **SynTerra leverages the power of Genomic Synteny (gene order conservation)**. By identifying the genomic "neighborhood" of a Gene of Interest (GOI) in a home species, SynTerra can anchor onto highly conserved flanking genes to precisely triangulate the location of the missing GOI in target species—even when the GOI sequence itself has mutated beyond standard recognition.

## 🌟 Key Features

* **Synteny-Driven Targeting:** Uses conserved macro-synteny to identify regions likely to harbor highly divergent alleles or novel paralogs.
* **Iterative Deep Search:** Employs a multi-wave, dynamic-scoring engine powered by `MMseqs2`, `tblastn`, and `miniprot` to perform aggressive, block-by-block evaluations of candidate regions.
* **Zero-Annotation Compatibility:** Works seamlessly with fully annotated, partially annotated, or raw unannotated `.fna` assemblies securely falling back to `Prodigal` for dynamic locus annotations.
* **Fully Auto-Pilot "Easy Mode":** Simply provide a UniProt ID. SynTerra will automatically resolve the sequence, download the optimal reference genome from NCBI, auto-fetch phylogenetically related target genomes, and execute the entire pipeline.
* **Publication-Ready Visualizations:** Automatically generates stunning, interactivity-rich, visually compressed, and phylogenetically sorted HTML synteny block maps.

---

## 🚀 Quick Start

SynTerra requires **Nextflow (`>=22.10.1`)** and a package manager like Conda/Mamba or Docker.

### 1. Installation

The easiest way to install SynTerra dependencies is via Conda:

```bash
git clone https://github.com/AndreasWz/SynTerra.git
cd SynTerra
conda env create -f environment.yml
conda activate syntenyfinder
```

### 2. Auto-Pilot (Easy Mode)
Want to find orthologs for the Human LY6E protein (`Q16553`) across related mammals? *SynTerra will do everything for you.*

```bash
nextflow run main.nf \
  --mode easy \
  --query_id Q16553 \
  --outdir results_ly6e \
  -profile standard
```

### 3. Bring-Your-Own-Data (Pro Mode)
If you have your own local `.fasta` assemblies and you want to lock the search to your specific datasets:

```bash
nextflow run main.nf \
  --mode pro \
  --query input/my_gene.fasta \
  --home_genome input/reference_genome.fna.gz \
  --home_gff input/reference_annotation.gff \
  --target_genomes "input/target_clades/*.fna" \
  --outdir results_custom \
  -profile standard
```

---

## 📊 Pipeline Outputs

SynTerra builds a cleanly structured output directory (`--outdir`) containing everything you need for publication and downstream analysis:

| Directory / File | Description |
| :--- | :--- |
| `*_synteny_plot.html` | 🏆 **The core output.** Interactive, beautifully rendered synteny ribbons connecting homologous genes across species layers, ordered phylogenetically. |
| `*_tree.nwk` | Newick tree file representing the evolutionary distance calculated between the discovered loci. |
| `regions/*.regions.bed` | Condensed BED tracks defining the exact genomic coordinate boundaries of the successfully identified syntenic arrays. |
| `downloaded_genomes/` | *(Easy Mode only)* Contains all automatically fetched `.fna` assemblies, NCBI metadata manifests, and assembly quality `.tsv` reports. |
| `synterra_report.json` | Deep technical JSON file detailing execution states, metrics, parameters used, and search statistics. |

---

## 📖 Deep Documentation

This `README.md` is just a quick-start guide. 

For a complete breakdown of every conceivable parameter, execution mechanism, profile tuning (HPC vs Laptop), and advanced edge-case handling, you **must** consult the exhaustive `USAGE.md` manual.

👉 **[Read the Full SynTerra USAGE.md Manual Here](USAGE.md)**

---

## 🛠️ Architecture

SynTerra performs a directed sequence of operations:
1. **Resolution:** Extrapolates sequence data and taxonomies from single IDs.
2. **Staging:** Downloads or links heavy `fasta` libraries.
3. **Home Anchoring:** Identifies the precise boundaries of the GOI on the reference species.
4. **Context Building:** Extracts the `N` adjacent upstream and downstream flanking genes to build a syntenic "fingerprint".
5. **Macro-Synteny Block Hunting:** Sweeps target genomes using `MMseqs2` to locate clustering fields of orthologous flanking genes.
6. **Iterative Local Search:** Aggressively combs the gaps between clustered flanking genes in the target species utilizing relaxed-stringency `tblastn` and `miniprot` sweeps to unearth the hidden GOI.
7. **Phylo-Sorting:** Aligns successful targets on an evolutionary tree matrix.
8. **Rendering:** Compresses vast genomic empty spaces and plots color-coded orthology ribbons directly to an interactive dashboard.

---

## 💡 Citations & Support
If you use SynTerra in your research, please link back to this repository. For bugs, features, or questions, please open an Issue!
