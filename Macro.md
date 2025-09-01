# Macro-Scale Analysis
To improve classification accuracy among patients and to identify the key differences between molecular subgroups, the samples were designated according to their molecular subtypes. We followed a traditional approach by identifying their differentially expressed genes; the count matrices for all samples were integrated under the Seurat functions in R.

## Data loading

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

## Data Objects
To improve manageability, we convert every data set in a Seurat Object and save it in the same folder.

```text
./Scripts/1_Data_Objects.R
```

## Data Integration
We loaded every MB object, filtered their spots by the percentage of mitochondrial genes, number of genes and total of transcripts, we applied the SCT method according to seurat to integrate all samples and diminish the batch effect. An integrated object is saved. We used the 3000 most variable genes to perform the correction and anchoring.

```text
./Scripts/2_integration_analysis.R
```
## Dimensionality Reduction and DGE analysis
To perform dimension reduction we computed the first 50 PCA components and a resolution of 0.8 for the uMAP. A folder named Results is created with the following plots and results files. In this section, samples were analyzed by MB molecular subtype. DEGs were pre-filtered by a $log_2(FC) <= 0.25$.

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
## Volcano Plots, Venn diagram and enrichment results
**Figure 2A** shows the volcano plots given de DEG analysis. It considered as differentially expressed gene with a $log_2|FC| \ge 1$ and and adjusted $p_{value}\le 0.01$.

**Figure 2B** shows a Venn diagram for the positive differentially expressed genes with the criteria of $log_2(FC) \ge 1$ and and adjusted $p_{value}\le 0.01$.

Using the DEG we performed Gene Set Enrichment Analysis (GSEA), but we do not found a particular biological process of pathway representative for every molecular subgroup. 

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
## Signatures plots

To assess how signature genes of each MB subtype are distributed across the tissue slides, we calculated an enrichment score for each spot using the differentially expressed genes. This analysis takes advantage of spatial transcriptomics, which integrates spatial location with subtype-specific gene expression patterns. For the definition of gene signatures, we focused exclusively on uniquely overexpressed genes, as these are more likely to serve as reliable markers.

```text
./Script/5_Signaures_Macro.R
```

As an initial evaluation, we selected the top 20 genes for each subtype with the highest $log_2(FC) \ge 1$ values. Ideally, spots within tissue sections corresponding to the same MB subtype as the signature should exhibit higher enrichment scores.

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

## Percentage of classified spots

We calculated the percentage of spots classified according to the gene signatures in each sample. We verified that these results were not biased by the criteria used to select the DE genes by computing the percentage of classified spots under different $log_2(FC)$ thresholds. 
```text
./Script/6_Signatures_Macro_log2.R
```
The script outputs are a scores_signature.csv file for each sample, were a matrix contain the percentage of classified spots for the Group 3, SHH and Group 4 signatures for every $log_2(FC)$ thresholds

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

## Bar plot of the percentage of classified spots

From the files ***scores_signature.csv** a plot is generated, one per sample showing the percentage of classified spots for every $log_2(FC)$ thresholds, the results correspond to the **Figure S9**.

```text
./Script/7_barplot_giotto.R
```
As well, results showed that the best threshold was the top20 for the DEGs with $log_2(FC)\ge 1$, figure 2C.

# Machine Learning Models
We wondered if the classification could be improved, so we migrate to a Machine Learning (**ML**) approach. We used machine learning and explainable artificial intelligence algorithms to develop a classifier for clinical MB subtypes based on stRNA-seq data. Our approach involved four key steps: selecting data for model training and testing, constructing an eXtreme Gradient Boosting (XGBoost) model for classifying the samples, evaluating model performance, and pinpointing key variables influencing the classification.

## Data adaptation
Based on the integrated data object, we extract the count matrix to be exported into **Python** for the ML analysis.

```text
./Script/8_integration_ML_2.R
```
Five files should be saved, the data non integrated (Data_counts), data integrated (Data_data), the name of the genes, the name of every spot according of the molecular assignment for the samples, and classes file where the labels were transformed to 0, 1 or 2, to represent the three molecular MB subtypes in the data.

```text
./Results
└── ML
    ├── Data_counts.csv
    ├── genes.tsv
    ├── labels.tsv
    ├── Data_data.txt
    └── classes.csv
 ```

```text
Class assignment
# 0 -> G4
# 1 -> SHH
# 2 -> G3
```

## XGBoost Model(s)

To identify biomarkers across the three molecular subtypes of medulloblastoma, we performed a machine-learning analysis of spatial transcriptomic data. We used XGBoost, a widely applied tree-based classification method, to distinguish among SHH, Group 3 (G3), and Group 4 (G4) subtypes.

To minimize batch effects, all samples were integrated into a single dataset using Seurat’s IntegrateData function, and the resulting counts were exported for analysis in Python. Spots were labeled according to subtype (SHH = 0, G3 = 1, G4 = 2). We randomly split the dataset, assigning 75% for training and 25% for testing. The XGBoost model was trained on the training set and evaluated on the test set.

To ensure robust and unbiased results, we applied five-fold cross-validation (k = 5). The dataset was divided into five equal parts; in each iteration, four parts were used for training while the remaining part was used for testing. This process was repeated five times so that every subset served as a test set once. Model performance was evaluated with a confusion matrix.

To identify key genes driving classification, we calculated SHapley Additive exPlanations (SHAP) values, which provide interpretable, game-theoretical insights into model predictions. We then generated aggregate plots showing the top-ranked genes with their SHAP values and average expression levels. Biomarkers were defined as genes consistently identified across all cross-validation iterations

```text
./Script/9_model.py
```
The python script saves the following results, one per iteration.
```text
./Results
└── ML
    ├── Figures
    │   ├──Confusion Matrix (Figure S10)
    │   └──Shap plot 
    ├── Model 
    │   └── XGBoost models
    └── SHAP_values
        ├──SHAP_G3.csv
        ├──SHAP_G4.csv
        └──SHAP_SHH.csv
 ```

## Shap Plots
As we use a five-fold cross-validation, there are 5 lists for the most important genes. We integrate the list for every MB subtype taking the mean and standard deviation of the SHAPS values. The results are represented as summary shap plots, **Figure 2E**. A total of 19, 15, and 18 signature genes were identified for the SHH, Group 3, and Group 4 MB molecular subtypes, respectively.

```text
./Script/10_Shap_plots.R
```

Moreover, their expression can be positive or negative related to the classification. To correlate the directionality of relation expression-classification we used the dependency plot which relates the SHAP value with the RNA expression, **Figure S11-S13**.

```text
./Script/11_SHAP_dependency.R
```

As a result, two plots would be generated for the shap values.

```text
./Results
└── ML
    └── Figures
        ├──Spatial feature plot
        └──heatmap
```
