# TREZ and Oligodendrocyte-Lineage Analysis
## Overview

This repository contains the R code used for the paper's transcriptomic
analysis of the transition region (TREZ) and oligodendrocyte precursor
cell (OPC) / oligodendrocyte (OL) lineage in the HC and TN groups.

The workflow includes clustering and cluster annotation, TREZ-focused
analysis, marker-score analysis after lineage gating, comparison of
differential-expression methods, and oligodendrocyte-lineage stage
scoring.
## Repository contents
A suggested repository layout is:
``` text
.
├── README.md
└── code.R
```

`code.R` refers to the supplied analysis script. 
The `Matrix/` directory is expected to contain one filtered count matrix per sample. 
## Analysis workflow
The script is organized into sequential analysis sections:
### 1 Quality control and clustering

- Load the HC, TN1, TN2, and TN3 count matrices.
- Create Seurat objects, merge samples, and filter spots using
- `nFeature_RNA > 200` and `percent.mt < 20`.
- Run SCTransform, PCA, UMAP, and neighbor detection.
- Evaluate a range of clustering resolutions and assign 10
- clusters using k-means on the first 30 principal components.
### 2 Cluster mapping and TREZ selection
- Examine the predefined myelin- and root-zone marker sets.
- Export marker DotPlots and cluster composition tables to support cluster annotation.
- Extract these clusters for downstream TREZ analyses.

### 3 TREZ marker scoring and differential expression
- Calculate UCell scores for OPC- and OL-associated gene signatures.
- Generate group and per-sample visualizations.
- Perform TN-versus-HC differential-expression analysis using MAST on the SCT assay.
- Export differential-expression tables, volcano plots, and GO/KEGG enrichment results when available.
### 4 Lineage-gated OPC/OL analysis
- Gate spots positive for either `Olig1` or `Olig2`.
- Calculate OPC and OL signature scores and their ratio.
- Export per-sample summaries, visualizations, and MAST
- differential-expression results.
### 5 Additional visualization and robustness analyses
- Compare Wilcoxon-based differential-expression results with the MAST results for selected genes. 
- Generate volcano plots, gene-level effect-size plots, and spot-level bootstrap intervals.
### 6 Oligodendrocyte-lineage stage scoring
- Calculate UCell scores for OPC, COP, NFOL, MFOL, and MOL gene sets.
- Derive maturation and early-to-late indices.
- Assign a dominant stage per spot and export stage scores, sample summaries, and plots.
## Requirements
The script uses R and the following R/Bioconductor packages in different
sections:
- Seurat
- MAST
- UCell
- clusterProfiler
- org.Rn.eg.db
- enrichplot
- EnhancedVolcano
- dplyr
- tidyr
- ggplot2
- patchwork
- ggpubr
- ggrepel
- svglite
- scales
- BiocParallel (used when available)
Install packages in an R session before running the script. For example:
``` r
install.packages(c(
  "Seurat", "dplyr", "tidyr", "ggplot2", "patchwork",
  "ggpubr", "ggrepel", "svglite", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")}
BiocManager::install(c(
  "MAST", "UCell", "clusterProfiler", "org.Rn.eg.db",
  "enrichplot", "EnhancedVolcano", "BiocParallel"))
```
