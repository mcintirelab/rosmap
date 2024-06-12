# Lipidomics Analysis Pipeline

## Description
This pipeline performs lipidomics data analysis for brain and serum samples, including data loading, filtering, imputation, normalization, clustering analysis, and bulk RNA-seq gene expression analysis.

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
