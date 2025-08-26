
## Atribution
This pipeline was developed by [Aarón Vázquez-Jiménez](https://scholar.google.com/citations?hl=en&user=p_LmnJcAAAAJ) in the [Human Systems Biology Group](https://resendislab.github.io/) at [INMEGEN](https://www.inmegen.gob.mx/). This pipeline was used to perfom a multi-scale functional analysis in spatial trancristomics data of medulloblastoma samples; [cite].

If you have any further question or comment, please do not hesitate to make contact:
  avazquez@inmegen.gob.mx

## Getting Started

In this study, we integrate spatial transcriptomics, machine learning, and explainable artificial intelligence to comprehensively analyze the spatial and transcriptional heterogeneity of the most prevalent and aggressive medulloblastoma subtypes (SHH, Group 3, and Group 4) in Mexican pediatric patients. By implementing a multiscale framework, we uncover extensive intra- and intertumoral variability that challenges traditional molecular classification.



The [pipeline](Pipeline.md) is structured as follows:

* [Introduction](Pipeline.md#Introduction)
  - [Sample Description](Pipeline.md#Sample-description)
  - [Experimental protocol](Pipeline.md#experimental-protocol)
  - [Molecular Genotyping](Pipeline.md#Molecular-genotyping)
  - [Data](Pipeline.md#Data)
* [Knee filter](Pipeline/pipeline.md#knee-filter)
* [Gene Over-dispersion](Pipeline/pipeline.md#gene-over-dispersion)
* [Gene number selection](Pipeline/pipeline.md#gene-number-selection)
* [Dimensionality Reduction](Pipeline/pipeline.md#Dimensionality-Reduction)
* [Clustering](Pipeline/pipeline.md#Clustering)
  - [Kmeans](Pipeline/pipeline.md#kmeans)
  - [EM](pipeline.md#expectation-maximization-algorithm)
