# Lipidomics Analysis Pipeline

## Description
This pipeline performs lipidomics data analysis for brain and serum samples, including data loading, filtering, imputation, normalization, clustering analysis, and bulk RNA-seq gene expression analysis. This analysis pipeline was developed for the article entitled "Brain and serum lipidomic profiles implicate Lands cycle acyl chain remodeling association with APOEε4 and mild cognitive impairment", published in Frontiers in Aging Neuroscience.
Please cite the article as follows if you use this pipeline:

## Citation:
Mares J, Costa AP, Dartora WJ, Wartchow KM, Lazarian A, Bennett DA, Nuriel T, Menon V and McIntire LBJ (2024) Brain and serum lipidomic profiles implicate Lands cycle acyl chain remodeling association with APOEε4 and mild cognitive impairment. Front. Aging Neurosci. 16:1419253. doi: 10.3389/fnagi.2024.1419253

## File Structure
- `01_Lipid_Inputs.R`: Loading and preparation of lipidomics data.
- `02_lipid_workflow_analysis.R`: Differential analyses and heatmap generation.
- `03_knn_WGCNA_pipeline.R`: KNN clustering and WGCNA analysis.
- `04_bulkRNA_seq_analysis.R`: Bulk RNA-seq gene expression analysis.

## Reference
This pipeline was used in the article published in Frontiers in Aging Neuroscience:
[Frontiers in Aging Neuroscience, 2024, DOI:10.3389/fnagi.2024.1419253](https://www.frontiersin.org/articles/10.3389/fnagi.2024.1419253/full).

## Requirements
- R
- tidyverse
- struct
- SummarizedExperiment
- pmp
- structToolbox
- limma
- lme4
- nlme
- expss
- RColorBrewer
- readr
- umap
- ggrepel
- gridExtra
- WGCNA
- ggpubr
- rstatix
- patchwork
- circlize

## Installation
To install the necessary dependencies, use the following command in R:

```r
install.packages(c("tidyverse", "struct", "SummarizedExperiment", "pmp", "structToolbox", "limma", "lme4", "nlme", "expss", "RColorBrewer", "readr", "umap", "ggrepel", "gridExtra", "WGCNA", "ggpubr", "rstatix", "patchwork", "circlize"))
