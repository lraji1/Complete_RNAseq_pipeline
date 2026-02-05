Complete RNA-seq Analysis Pipeline
A comprehensive R-based pipeline for bulk RNA-seq quality control, differential expression analysis, and enrichment analysis

Overview
This repository contains R scripts for end-to-end RNA-seq analysis, from raw count matrices and metadata file. 
The pipeline uses DESeq2 null model for quality control followed by Principal component analysis for visualing the sample distribution.

Key Features
Quality Control & Preprocessing

Multi-level QC: Expression density plots, library size distributions, batch effect visualization
Outlier detection: Identification of outlier samples( using±2 SD threshold) with visualization
Batch assessment: RNA.Batch and Library.Batch tracking with color-coded visualizations


Dimensionality Reduction

PCA analysis: Variance stabilizing transformation (VST) with top variable genes
Multi-variable encoding: Integrated visualization of batch (color), treatment (shape), sex (shape), and strain (ellipses)
PC association testing: ANOVA-based testing for metadata-PC relationships
Correlation analysis: Spearman (categorical) and Pearson (continuous) correlations with PC2

Differential Expression

DESeq2 framework: Standard Wald test with multiple model designs
Model based contrasts: Support for treatment, sex, strain, and interaction effects
Batch correction: Employing DESeq2 batch correction covariate approach in its robust negative binomial generalized linear model
Deterministic tie-breaking: Wald statistic ranking for gene set enrichment analysis

Enrichment Analysis

fGSEA pipeline: Fast gene set enrichment analysis with 10,000 permutations
Multiple databases: MSigDB Hallmark, GO-BP, KEGG, Reactome support . etc
Cross-model comparison: Heatmap visualization of pathway enrichment across models based on leading edge genes
Top pathway selection: Mean absolute normalized enrichment score (NES) ranking