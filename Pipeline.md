# Introduction
To comprehensively investigate the spatial organization and transcriptional heterogeneity of medulloblastoma (**MB**), we implemented a multi-scale analytical framework encompassing three spatial resolutions: macroscale, mesoscale, and microscale. At the macroscale, we compared spatial transcriptomic profiles across all medulloblastoma subtypes—Group 3, Group 4, and SHH—to identify subtype-specific spatial gene expression patterns and architectural signatures. The mesoscale focused on intra-patient comparisons by analyzing paired tumor samples from the same individual, enabling us to assess regional variation and transcriptional divergence within a shared genetic background. Finally, the microscale level dissected individual tumor samples at high spatial resolution, revealing localized patterns of gene expression, cellular neighborhoods, and microenvironmental niches that contribute to the spatial heterogeneity of the tumor. This hierarchical approach enabled the contextualization of  tumor biology across nested spatial levels, from subtype-level distinctions to fine-grained intra-tumoral architecture.

This document is provided to ensure replication of the results presented; moreover, the detailed experimental procedures and data interpretation are available in in detail at the [publication]().

## Sample description
Five patients with MB diagnoses were included in the study, informed consent and assent were obtained from patients and their families prior to sample collection. The protocol was approved by the Institutional Ethics Committee (protocol INP 2022/003) CONBIOÉTICA-09-CEI-025-20161215, and the sample was handled in accordance with the Declaration of Helsinki. 

Our data included the spatial transcriptomics of seven paraffin-embedded samples from five pediatric medulloblastomas.

## Experimental protocol
The tumor tissue was fixed in formalin and embedded in paraffin for histopathological diagnosis. The Visium Spatial Gene Expression workflow was performed according to the manufacturer's instructions (10x Genomics, Pleasanton, CA, USA) with slight modifications. The Visium Spatial Gene Expression libraries were sequenced on an Illumina sequencing platform NextSeq 2000 using paired-end sequencing with read lengths suitable for obtaining high-quality spatial transcriptomics data. For more details, check the published document.

## Molecular Genotyping
The molecular classification of the MB samples was determined through gene expression profiling using microarray analysis (Affymetrix GeneChip U133 Plus 2.0). An unsupervised least squares model was applied to the expression profiles, which were compared against previously curated molecular classification data from  patients with medulloblastoma obtained through data mining.

 The samples were obtained from patients with the three most common and aggressive medulloblastoma subtypes: SHH, Group 4, and Group 3, with a patient distribution of 3, 1, and 1, respectively. The Group 4 samples showed transcriptional programs overlapping with SHH, but were classified as Group 4 due to their predominantly Group 4-like characteristics. For both the SHH and Group 4 subtypes, two samples were taken from the same patient. 

## Data
All the generated data are available in the GEO under the accession number **[GSE293081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293081)**.

Seven FFPE samples from five patients were sequenced under spatial transcriptomics protocols with the following name and molecular classification relationship:

| Sample Name    | Molecular Classification | Patient |
|----------------|--------------------------|---------|
| 09062023_G3-1  |        Group 3           |    1    |
| 09062023_SHH-1 |          SHH             |    2    |
| 09062023_SHH-2 |          SHH             |    3    |
|     SHH-3A     |          SHH             |    4    |
|     SHH-3B     |          SHH             |    4    |
|     G4-1A      |        Group 4           |    5    |
| 09062023_G4-1B |        Group 4           |    5    |

There were 1:3:1 patients for the MB subtypes Group 3:SHH: Group 4 MB, respectively. For one SHH and one Group 4 patient, two tumor sections were sequenced .

# Data Analysis
As mentioned, three approaches were used. The data analysis was performed using **R** v.4.4.3 and **Python** v. 3.12 depicted as follows:

* [Macro-Scale Analysis](Macro.md#macro-scale-analysis)
  - [DGE approach](Macro.md#dimensionality-reduction-and-dge-analysis)
  - [Machine Learning model](Macro.md#machine-learning-models)
* [Meso-Scale Analysis]()
* [Micro-Scale Analysis]()
