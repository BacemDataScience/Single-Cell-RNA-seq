# State-Focused Single-Cell RNA Sequencing Reveals Reproducible Cytotoxic and Transitional Immune Programs Following In Vitro Expansion

This repository contains the R scripts, validation workflows, and analysis framework used in the study:

"State-Focused Single-Cell RNA Sequencing Reveals Reproducible Cytotoxic and Transitional Immune Programs Following In Vitro Expansion"

## Overview
This repository provides a reproducible state-focused single-cell RNA-seq analysis framework for studying transcriptional remodeling of immune cells following in vitro expansion. The workflow integrates discovery and independent cross-dataset validation analyses using publicly available human immune-cell and melanoma tumor infiltrating lymphocyte (TIL) datasets.

The analytical framework emphasizes:

- functional-state interpretation,
- reproducibility across independent datasets,
- overlapping cytotoxic/exhaustion/transitional programs,
- and systems-level immune-state architecture rather than rigid cluster identity assignment.

## Dataset
- Source: Gene Expression Omnibus (GEO)
- GSE211644:	Discovery dataset comparing freshly isolated and in vitro expanded immune cells
- GSE120575:	Independent melanoma tumor infiltrating lymphocyte (TIL) validation dataset
- GSE123139:	Independent melanoma T-cell validation dataset

## Main analyses
Quality control and filtering
SCTransform normalization
Seurat v5 integration workflow
Dimensionality reduction (PCA, UMAP)
Graph-based clustering
Functional-state annotation
Canonical marker visualization
Functional module scoring
Cytotoxicity program analyses
Exhaustion-associated program analyses
Proliferation program analyses
Interferon signaling analyses
Differential abundance analysis
Differential expression analysis
Cross-dataset validation
Program overlap analyses
Correlation heatmaps
Cluster-level functional-state analyses
Patient-level heterogeneity analyses
Recurrent gene validation analyses
Publication-quality figure generation

## Software requirements
- R (≥ 4.2)
- Seurat (v5 or compatible)
- tidyverse
- data.table
- ggplot2
- patchwork
- pheatmap
- future
- future.apply
- ggrastr

## Study Focus
This project demonstrates that:
- in vitro expansion reshapes immune landscapes primarily through functional-state remodeling,
- recurrent cytotoxic and transitional transcriptional programs are reproducible across datasets,
- and immune single-cell datasets are better interpreted through overlapping functional-state architectures rather than rigid cluster-centric frameworks.
