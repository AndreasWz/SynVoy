# Ground Truth Validation Strategy using Ly6 Curated Annotations

## 1. Overview of the Ground Truth Data

### Ly6 / 3FTx dataset (`ground_truth/ly6_3ftx`)
- Paper: `ground_truth/ly6_3ftx/s41467-023-40550-0.pdf`.
- Curated annotations: `ground_truth/ly6_3ftx/Ly6/` contains 21 GFF files (GFF3 headers report `##source-version geneious 2025.1.2`).
- Feature + evidence mix observed across these GFFs: `gene`, `mRNA`, `CDS`, `exon`, `ncRNA`, `misc_RNA`, plus `source` / `extracted region` records. The `source` column includes `Geneious`, `maker`, `tblastx`, `venom_blast`, and `Long_blast`, indicating both curated models and evidence tracks.
- File list:
```text
CM019148_extraction.gff
CM019150_extraction.gff
CM019151_extraction.gff
Gaga_NC_006089.gff
Hosa_NC_000008.gff
NC_041731_extraction.gff
NW_001794279.gff
NW_003338769.gff
NW_005842480.gff
NW_006534208.gff
NW_006678715.gff
NW_015387915.gff
NW_015390125.gff
NW_018151385.gff
NW_020769337_extraction.gff
NW_020769389.gff
NW_020776291.gff
NW_025336455.gff
Phci_NW_018344010.gff
ScVE01q_1072;HRSCAF=1231.gff
Xetr_NC_030682.gff
```

### Melettin dataset (`ground_truth/melettin`)
- Paper: `ground_truth/melettin/s12915-023-01656-5.pdf`.
- Supplement: `ground_truth/melettin/8052397.zip` (46 files total; many are sub-archives with data tables and annotations).
- Genome annotation bundles for multiple bee species (contain `*.gff3`, `*.proteins.fa.gz`, `*.cds-transcripts.fa.gz`, `*.mrna-transcripts.fa.gz`, masked genomes, and repeat mask tracks; e.g., the Xylocopa bundle is a compact GFF/FASTA set, while Melipona is a full funannotate-style directory tree): `AdditionalFile41_GenomeAnnotation_Xylocopa_violacea.zip`, `AdditionalFile40_GenomAnnotation_Melipona_beecheii.zip`, `AdditionalFile39_GenomAnnotation_Tetragonula_carbonaria.zip`, `AdditionalFile38_GenomAnnotation_Euglossa_dilemma.zip`, `AdditionalFile37_GenomAnnotation_Colletes_gigas.zip`.
- Toxin gene annotation bundle: `AdditionalFile42_ToxinGeneAnnotations.zip` with many GFF3/GFF files (filenames include tags like `_SP`, `_VA`, `_APH`, `_PLA2`, `_DPP4`, `_SECA`, `_HYAL`, `_ICA`, plus `gff_modified_Koludarov_et_al_2021.gff3`).
- Melittin-like sequence resources: `AdditionalFile33_ML_Sequences_table.csv`, `AdditionalFile43_ML_sequences_full.fasta`, `AdditionalFile44_ML_sequences_mature.fasta`, plus peptide/FASTA text sets like `AdditionalFile35.txt`, `AdditionalFile34.txt`, `AdditionalFile30_Vollenhovia_melittin_like.txt`, and an alignment-style text file `AdditionalFile32.txt`.
- Comparative/analysis artifacts: `AdditionalFile45_SimilarityMatrices_for_PBVPs.zip` (SVG similarity matrices), `AdditionalFile46_PhylogeneticTrees_for_PBVPs_in ML_analysis.zip` (alignments + trees in `.phy`, `.nexus`, `.svg`), and `AdditionalFile11_Protein_3Dspaces_Fig3_Fig7.zip` (interactive HTML 3D PCA/UMAP spaces).
- Transcriptome/peptide assemblies: `AdditionalFile4_Xylocopa_VG_Assembly_Transdecoder.pep.zip`, `AdditionalFile5_Halictus_VG_Assembly_Transdecoder.pep.zip`, `AdditionalFile6_Apis_VG_Assembly_Transdecoder.pep.zip`.
- Additional docs/tables/scripts: multiple PDFs, spreadsheets, and a small script (`AdditionalFile36InHousePerlScript.txt`) that may be useful for provenance or reference tables.

---

## 2. Validation Objectives
Creating a test harness using these files allows us to quantify exactly how well SynVoy addresses the "magic parameter" problem. 

Our validation goals are twofold:
1. **Locus Discovery (Macro-Synteny):** Does `cluster_grs.py` correctly group the orthologous anchors to define the true locus spanning the Ly6 genes?
2. **GOI Gap-Filling (Micro-Synteny):** Does `iterative_search_runner.py` accurately find and stitch together the highly divergent/fragmented Ly6 query sequences, correctly aligning them precisely over the annotated ground truth `CDS` boundaries?

---

## 3. Systematic Testing Methodology

To leverage these GFFs as ground truth, we should implement the following test protocol:

### Step A: Test Data Preparation
- **Queries & Profiles:** Use a standard Human or Chicken Ly6 amino acid sequence (e.g., PSCA or LYPD2) as the query. Create a profile of the known flanking anchor genes from the GFF extractions.
- **Target Assemblies:** Fetch the full target genomes matching these extractions (e.g., the full *Naja naja* and *O. anatinus* references).

### Step B: Automated Intersection Script (`validate_ly6_results.py`)
We will develop a Python scoring script (using `pandas` and interval intersection logic like `pyranges` or `bedtools`) that accepts:
1. The SynVoy pipeline results (Locus bounds, GOI hit coordinates in TSV/BED).
2. The `Ly6/*.gff` ground truth files.

The script will automatically calculate the following metrics:

#### Metric 1: Locus Recall (Macro-Synteny)
* **Goal:** Determine if SynVoy found the right biological neighborhood.
* **Test:** Compare the start/end coordinates of the bounding anchors in the ground truth GFF to the final locus boundaries reported by SynVoy.
* **Score:** `Distance Offset (bp)` -> How tight is SynVoy's bounding box compared to the true anchor boundaries?

#### Metric 2: GOI Sensitivity (Exon-Level Recall)
* **Goal:** Determine if MMseqs2 iterative search missed any valid exons of the Ly6 gene.
* **Test:** Intersect SynVoy's `GOI_HITS` coordinates with the `CDS` coordinates inside the ground truth GFF for that specific target name.
* **Score:** `CDS Tracking (%)` -> (Total ground truth `CDS` bases covered by a SynVoy hit) / (Total ground truth `CDS` length). *Expectation: High sensitivity (e.g., >85%).*

#### Metric 3: GOI Specificity (Noise Reduction)
* **Goal:** Determine if SynVoy is producing spurious, false-positive hits from repetitive elements or intergenic noise.
* **Test:** Intersect SynVoy's `GOI_HITS` with regions outside the bounds of the ground truth `gene` features.
* **Score:** `False Positive Rate` -> Proportion of SynVoy candidate hits that do NOT overlap with any annotated Ly6 `CDS` or `mRNA`. *Expectation: Low FPR (e.g., <5%).*

---

## 4. Why This is Crucial for Future Development
Currently, when we alter SynVoy's internal logic—like modifying the `bad_quality_policy` timeouts or tuning the synteny scoring heuristic in `cluster_grs.py`—we rely on manual inspection (viewing the HTML plots) to verify if the pipeline still works.

By implementing this Ly6 validation suite:
1. **Regression Testing:** We can trigger this test automatically via a shell script. If a change drops the `CDS Tracking %` or completely fails to discover the locus in the platypus dataset, we immediately know the change broke something fundamental.
2. **Dynamic Parameter Tuning:** Since our next goal is *Dynamic Locus Parameterization* (adjusting search thresholds based on the target), we can use these ground truth scores as a fitness metric to programmatically optimize SynVoy's logic. If a dynamic parameter setting yields a perfectly bounded Ly6 locus with 95% CDS tracking, we know the heuristic succeeds.
