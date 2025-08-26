# Introduction
To comprehensively investigate the spatial organization and transcriptional heterogeneity of medulloblastoma (**MB**), we implemented a multi-scale analytical framework encompassing three spatial resolutions: macroscale, mesoscale, and microscale. At the macroscale, we compared spatial transcriptomic profiles across all medulloblastoma subtypes—Group 3, Group 4, and SHH—to identify subtype-specific spatial gene expression patterns and architectural signatures. The mesoscale focused on intra-patient comparisons by analyzing paired tumor samples from the same individual, enabling us to assess regional variation and transcriptional divergence within a shared genetic background. Finally, the microscale level dissected individual tumor samples at high spatial resolution, revealing localized patterns of gene expression, cellular neighborhoods, and microenvironmental niches that contribute to the spatial heterogeneity of the tumor. This hierarchical approach enabled the contextualization of  tumor biology across nested spatial levels, from subtype-level distinctions to fine-grained intra-tumoral architecture.

This document is provided to ensure the replications of the results presented; moreover, the detailed experimental procedures and data interpretation are available in extense at the [publication]().

## Sample description
Five patients with MB diagnoses were included in the study, the informed consent and assent from the patients and their families were signed prior to sample collection. The protocol was approved by the Institutional Ethics Committee (protocol INP 2022/003) CONBIOÉTICA-09-CEI-025-20161215, and the sample was handled in accordance with the Declaration of Helsinki. 

Our data included the spatial transcriptomics of seven paraffin-embedded samples from five pediatric medulloblastomas.

## Experimental protocol
The tissue tumor was fixed in formol and embedded in paraffin for histopathological diagnosis. The Visium Spatial Gene Expression workflow was performed according to the manufacturer's instructions (10x Genomics, Pleasanton, CA, USA) with slight modifications. The Visium Spatial Gene Expression libraries were sequenced on an Illumina sequencing platform NextSeq 2000 using paired-end sequencing with read lengths suitable for obtaining high-quality spatial transcriptomics data. For more detail, check the published document.

## Molecular Genotyping
The molecular classification of the MB samples was determined through gene expression profiling using microarray analysis (Affymetrix GeneChip U133 Plus 2.0). An unsupervised least squares model was applied to the expression profiles, which were compared against previously curated molecular classification data from  patients with medulloblastoma obtained through data mining.

 The samples were obtained from patients with the three most common and aggressive medulloblastoma subtypes: SHH, Group 4, and Group 3, with a patient distribution of 3, 1, and 1, respectively. The Group 4 samples showed transcriptional programs overlapping with SHH, but were classified as Group 4 due to their predominantly Group 4-like characteristics. For both the SHH and Group 4 subtypes, two samples were taken from the same patient. 

## Data
All the generated data are available in the GEO under the accession number **[GSE293081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293081)**.

Seven FFPE samples from five patients were secuenced under spatial transcriptmics protocols with the following name and molecular classficiation relationship:

| Sample Name    | Molecular Classification | Patient |
|----------------|--------------------------|---------|
| 09062023_G3-1  |        Group 3           |    1    |
| 09062023_SHH-1 |          SHH             |    2    |
| 09062023_SHH-2 |          SHH             |    3    |
|     SHH-3A     |          SHH             |    4    |
|     SHH-3B     |          SHH             |    4    |
|     G4-1A      |        Group 4           |    5    |
| 09062023_G4-1B |        Group 4           |    5    |

There were 1:3:1 patients for the MB subtypes Group 3:SHH: Group 4 MB, respectively. For a patient of the SHH and group 4 two tumor section were sequenced .

# Data Analysis
As mentioned, three approaches were used. The data analysis were performed using **R** v.4.4.3 and **Python** v. 3.12.

## Macro-Scale Analysis
To improve classification accuracy among patients and to identify the key differences between molecular subgroups, the samples were designated according to their molecular subtypes. We followed a traditional approach by identifying their differentially expressed genes; the count matrices for all samples were integrated under the Seurat functions in R.

### Data loading

You can download either the raw files or the processed files, if you choose for the raw files, we used the Spaceranger/2.0 implementation. The alignment and counting were generated via the count function with the STAR/2.7.2a aligner. The sequences were compared against the reference human genome (GRCh38). The sequenced data were correlated with the spatial coordinates in the H&E image on the basis of the information from the spatial identifiers. All reads and microsites with no associated tissue were removed; the resulting filtered matrix from SpaceRanger was used for subsequent analysis. So, in your working directory should look like the following.

```text
./
├── 09062023_G3-1
│   └── outs
│       ├── filtered_feature_bc_matrix
│       │    ├──barcodes.tsv.gz
│       │    ├──features.tsv.gz
│       │    └──matrix.mtx.gz
│       └── spatial
│           ├──aligned_fiducials.jpg
│           ├──detected_tissue_image.jpg
│           ├──scalefactors_json.json
│           ├──spatial_enrichment.csv
│           ├──tissue_hires_image.png
│           ├──tissue_lowres_image.png
│           └──tissue_positions.csv
├── 09062023_SHH-1
├── 09062023_SHH-2
├── SHH-3A
├── SHH-3B
├── G4-1A
└── 09062023_G4-1
```
If you choose to download the processed files, you must arrange your working folder as the previous tree to ensure the correct loading by the Seurat functions. So, header in file name must erase and plece it the corresponding folder.

### Data Objects
To improve manageability, we convert every data set in a Seurat Object and save it in the same folder.

```text
./Scripts/1_Data_Objects.R
```

### Data Integration
We loaded every MB object, filtered their spots by the percentage of mithocondrial genes, number of genes and total of transcripts, we aplied the SCT method according to seurat to integrate all samples and diminish the bash effect. An integrated object is saved. We used the 3000 most variable genes to perform the correction and anchoring.

```text
./Scripts/2_integration_analysis.R
```
### Dimensionality Reduction and DGE analysis
To perform dimension reduction we computed the firts 50 PCA components and a resolutions of 0.8 for the uMAP. A folder named Results is created with the following plots and results files. In this sections, samples were analyzed by MB molecular subtype. DGE genes were pre-filtered by a $log_2(FC) <= 0.25$.

```text
./Scripts/3_DR_DEG.R
```

```text
./Results
├── plots
│   └── uMAP plot (Figure S1)
└── data
    └── DGE analysis (Supplementary Table 1)
```
### Volcano Plots, Venn diagram and enrichment results
**Figure 2A** shows the volcano plots given de DEG analysis. It considered as differentially expresed gene with a $log_2|FC| \ge 1$ and and adjusted $p_{value}\le 0.01$.

**Figure 2B** shows a venn diagram for the postive differentially expresed genes wiht the criteria of $log_2(FC) \ge 1$ and and adjusted $p_{value}\le 0.01$.

Using the DEG we performed Gene Set Enrichment Analysis (GSEA), but we do not found a particular biological process of pathway representative for every molecular subgruop. 

```text
./Scripts/4_DEG_Macro.R
```

The generated files are depicted as follows:

```text
./Results
├── plots
│   ├── Volcano plots
│   └── Venn Diagram
└── data
    ├── GSEA results for every MB subtype
    └── Signature  macro file (DEGs)
 ```
### Signarutes plots

To assess how signature genes of each MB subtype are distributed across the tissue slides, we calculated an enrichment score for each spot using the differentially expressed genes. This analysis takes advantage of spatial transcriptomics, which integrates spatial location with subtype-specific gene expression patterns. For the definition of gene signatures, we focused exclusively on uniquely overexpressed genes, as these are more likely to serve as reliable markers.

```text
./Script/5_Signaures_Macro.R
```

As initial evaluation, we selected the top 20 genes for each subtype with the highest $log_2(FC) \ge 1$ values. Ideally, spots within tissue sections corresponding to the same MB subtype as the signature should exhibit higher enrichment scores.

Seven figures would be generated with three panels, one per sample, **Figure S2 - Figure S8**.  A. The association to the given signature with molecular classification. B. Heatmap representing the gene signature for the three considered MB subtypes. C. Correlation analysis for the signature scores.

``` text
./Results
└── plots
    ├── Figure S2
    ├── Figure S3
    ├── Figure S4
    ├── Figure S5
    ├── Figure S6
    ├── Figure S7
    └── Figure S8
```

### Percentage of classified spots

We calculated the percentage of spots classified according to the gene signatures in each sample. We verified that this results was not biased by the criteria used to select the DE genes by computing the percentage of classified spots under different $log_2(FC)$ thresholds. 
```text
./Script/6_Signatures_Macro_log2.R
```
The script outputs are a scores_signature.csv file for each sample, were a matrix contain the percentage of claasified spots for the Group 3, SHH and Group 4 signatures for every $log_2(FC)$ thresholds

```text
./Results
└── data
    ├── G3-1_scores_signature.csv
    ├── SHH-1_scores_signature.csv
    ├── SHH-2_scores_signature.csv
    ├── SHH-3A_scores_signature.csv
    ├── SHH-3B_scores_signature.csv
    ├── G4-1A_scores_signature.csv
    └── G4-1B_scores_signature.csv
 ```

### Bar plot of the percentage of classified spots

From the files ***scores_signature.csv** a plot is generated, one per sample showing the percentage of calssified spots for every $log_2(FC)$ thresholds, the results correspond to the **Figure S9**.

```text
./Script/7_barplot_giotto.R
```
As well, results showed that the best threshold was the top20 for the DEGs with $log_2(FC)\ge 1$, figure 2C.

## Machine Learning Models
We wondered if the classification could be improved, so we migrate to a Machine Learning (**ML**) approach. We used machine learning and explainable artificial intelligence algorithms to develop a classifier for clinical MB subtypes based on stRNA-seq data. Our approach involved four key steps: selecting data for model training and testing, constructing an eXtreme Gradient Boosting (XGBoost) model for classifying the samples, evaluating model performance, and pinpointing key variables influencing the classification.

### Data adaptation
Based on the integrated data object, we extract the count matrix to be exported into Python for the ML analysis.

```text
./Script/8_integration_ML_2.R
```
Five files should be saved, the data non integrated (Data_counts), data integrated (Data_data), the name of the genes, the name of every spot according of the molecular assingation for the samples, and classes file where the labels were transformed to 0, 1 or 2, to represent the three molecular MB subtypes in the data.

```text
./Results
└── ML
    ├── Data_counts.txt
    ├── genes.tsv
    ├── labels.tsv
    ├── Data_data.txt
    └── classes.csv
 ```

```text
Classes asignation
# 0 -> G4
# 1 -> SHH
# 2 -> G3
```

### XGBoost Model

```text
./Script/9_model.py
```
