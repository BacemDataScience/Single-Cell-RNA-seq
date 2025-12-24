# Single-Cell RNA-seq Analysis of In Vitro Expanded Immune Cells

This repository contains the R scripts and analysis workflow used in the study  
"In Vitro Expansion Reprograms Cytotoxic Immune States in Single Cell RNA Sequencing Data".

## Overview
The pipeline performs a state-focused single-cell RNA-seq analysis comparing freshly profiled and in vitro expanded immune cells using publicly available data (GSE211644). Analyses include quality control, normalization, dataset integration, clustering, differential abundance testing, differential expression analysis, functional module scoring, and gene-level visualization.

## Dataset
- Source: Gene Expression Omnibus (GEO)
- Accession: GSE211644

## Main analyses
- Quality control and filtering
- SCTransform normalization and dataset integration
- Dimensionality reduction (PCA, UMAP)
- Graph-based clustering
- Canonical marker visualization
- Functional module scoring (cytotoxicity, exhaustion, proliferation)
- Differential abundance analysis
- Differential expression analysis
- Gene-level spotlight analyses (e.g. CCL4, GZMK, GNLY, FGFBP2)

## Software requirements
- R (≥ 4.2)
- Seurat (v5 or compatible)
- tidyverse
- data.table
- ggplot2
- patchwork

## Reproducibility
All figures and tables reported in the manuscript were generated directly from the scripts provided in this repository. Intermediate outputs (CSV files and plots) are saved to disk to facilitate inspection and reuse.
