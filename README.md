# 🧬 Complete RNA-seq Analysis Pipeline

A comprehensive R-based pipeline for bulk RNA-seq quality control, differential expression analysis, and enrichment analysis.

## 📋 Overview

This repository contains R scripts for end-to-end RNA-seq analysis, from raw count matrices and metadata file. The pipeline uses DESeq2 null model for quality control followed by Principal component analysis for visualizing the sample distribution.

## ✨ Key Features

### 🔍 Quality Control & Preprocessing

- **Multi-level QC**: Expression density plots, library size distributions, batch effect visualization
- **Outlier detection**: Identification of outlier samples (using ±2 SD threshold) with visualization
- **Batch assessment**: RNA.Batch and Library.Batch tracking with color-coded visualizations

### 📊 Dimensionality Reduction

- **PCA analysis**: Variance stabilizing transformation (VST) with top variable genes
- **Multi-variable encoding**: Integrated visualization of batch (color), treatment (shape), sex (shape), and strain (ellipses)
- **PC association testing**: ANOVA-based testing for metadata-PC relationships
- **Correlation analysis**: Spearman (categorical) and Pearson (continuous) correlations with PC2

### 🧪 Differential Expression

- **DESeq2 framework**: Standard Wald test with multiple model designs
- **Model-based contrasts**: Support for treatment, sex, strain, and interaction effects
- **Batch correction**: Employing DESeq2 batch correction covariate approach in its robust negative binomial generalized linear model
- **Deterministic tie-breaking**: Wald statistic ranking for gene set enrichment analysis

### 🔬 Enrichment Analysis

- **fGSEA pipeline**: Fast gene set enrichment analysis with 10,000 permutations
- **Multiple databases**: MSigDB Hallmark, GO-BP, KEGG, Reactome support, etc.
- **Cross-model comparison**: Heatmap visualization of pathway enrichment across models based on leading edge genes
- **Top pathway selection**: Mean absolute normalized enrichment score (NES) ranking

## 📁 Repository Structure
```
Complete_RNAseq_pipeline/
├── 📁 RNA_SEQ_Quality_Control/          
│   └── 📜 Time_point_60_full_QC_code.R       # Quality control scripts           
│    
├── 📁 RNA_SEQ_analysis_Deseq2/                          
│   ├── 📁 withoutPC2_as_covariate/
│   │   ├── 📜 biplot_v2.R
│   │   ├── 📜 cross_path_comparison_detailedv1.R
│   │   ├── 📜 no_shrinkage_fgsea_wald_ranking_withoutPC2.R
│   │   └── 📜 THC_wald_test_60_no_PC2v2_noOutlier_noShrinkage.R
│   │
│   └── 📁 withPC2_as_covariate/                
│       ├── 📜 biplot_v2.R
│       ├── 📜 cross_path_comparison_detailedv1.R
│       ├── 📜 no_shrinkage_fgsea_wald_ranking_withPC2.R
│       └── 📜 THC_wald_test_60_PC2_noOutlier_noShrinkage.R
│
└── 📄 README.md                                # This file
```

## 📦 Requirements

### R Packages
```r
# Core analysis
DESeq2

# Enrichment
fgsea, msigdbr

# Visualization
ggplot2, pheatmap, ggrepel, patchwork, RColorBrewer

# Data manipulation
dplyr, tidyr, matrixStats
```

## 🚀 Quick Start
```r
# 1. Load count and metadata
counts <- read.csv("count.csv", row.names = 1)
metadata <- read.csv("metv2.csv", row.names = 1)

# 2. Run QC pipeline
source("RNA_SEQ_Quality_Control/Time_point_60_full_QC_code.R")

# 3. Differential expression analysis
source("RNA_SEQ_analysis_Deseq2/withPC2_as_covariate/THC_wald_test_60_PC2_noOutlier_noShrinkage.R")

# 4. Functional enrichment
source("RNA_SEQ_analysis_Deseq2/withPC2_as_covariate/no_shrinkage_fgsea_wald_ranking_withPC2.R")
```

## 📥 Input Data Format

**Count matrix** (`count.csv`): Genes (rows) × Samples (columns)
- Integer counts from alignment tools (STAR, RSEM, Salmon)
- Gene IDs as row names

**Metadata** (`metv2.csv`): Samples (rows) × Variables (columns)
- Required: `Sample.ID`, `Treatment`, `Sex`, `Strain`, `Time.Point`
- Optional: `RNA.Batch`, `Library.Batch`, `RIN`, `Conc`

## 📤 Output Files

### Quality Control
- Expression density plots (with/without outliers)
- PCA plots (individual variables + combined multi-factor)
- Sample correlation heatmaps
- PC-variable association tables
- Outlier detection reports

### Differential Expression
- Normalized count matrices
- DESeq2 results tables with log2FC, p-values, FDR
- Volcano plots and MA plots
- Top differentially expressed gene lists

### Enrichment Analysis
- fGSEA results for multiple gene sets
- Enrichment heatmaps across conditions
- Leading edge gene lists
- Pathway-level statistics

## 🎯 Analysis Highlights

- **Filtering**: Genes with ≥10 counts in ≥3 samples
- **Normalization**: DESeq2 median-of-ratios
- **Transformation**: Variance stabilizing transformation (VST) for visualization
- **Multiple testing correction**: Benjamini-Hochberg FDR
- **Significance thresholds**: FDR < 0.05, |log2FC| > 1 (customizable)

## 📚 Citation

If you use this pipeline, please cite:
- Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology*
- Korotkevich G et al. (2021). "Fast gene set enrichment analysis." *bioRxiv*

