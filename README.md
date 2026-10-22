# balf-mngs-network-pipeline

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Ubuntu%20Linux-lightgrey.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Reproducibility-Verified-success.svg)]()

A reproducible, high-throughput R analytical pipeline for **Bronchoalveolar Lavage Fluid (BALF) metagenomic next-generation sequencing (mNGS)** cohorts. 

This repository implements batch-effect correction, ecological diversity profiling, a **22-method multi-paradigm co-occurrence network framework**, 1,000-iteration bootstrap stability estimation, dual consensus/robust keystone validation, and factor-driven pathogen rewiring analysis across clinical and host strata.

---

## Table of Contents
- [Workflow Overview](#workflow-overview)
- [Repository Structure](#repository-structure)
- [Prerequisites & Dependencies](#prerequisites--dependencies)
- [Input Data Requirements](#input-data-requirements)
- [Quick Start](#quick-start)
- [Modular Figure Regeneration](#modular-figure-regeneration)
- [Figure & Deliverables Index](#figure--deliverables-index)
- [Performance & Hardware Considerations](#performance--hardware-considerations)
- [Citation & License](#citation--license)

---

## Workflow Overview

```text
[Raw Kraken2 Species Matrix + Clinical Metadata]
                       │
                       ▼ (S0: Preprocessing & Sensitivity Calibration)
  [Cohort Prevalence >= 5% Filter] ───► Figure 01: Prevalence Calibration & Sparsity
                       │
                       ▼ (S0b: ConQuR Conditional Quantile Regression)
  [Batch-Corrected Matrix + CLR] ──────► Figure 02: PCA Batch Correction Diagnostics
                       │
         ┌─────────────┴─────────────┐
         ▼                           ▼
(S1: Alpha Diversity)       (S2: Beta Diversity & PERMANOVA)
  Observed, Shannon, Simpson   Bray-Curtis Distance + adonis2 (999 perm)
  └──► Figure 03               └──► Figure 04
         │                           │
         └─────────────┬─────────────┘
                       ▼ (S3: Taxonomic Profiling)
       Top-20 Species Composition Heatmap Across Host Factors
       └──► Figure 05
                       │
                       ▼ (S4: 22-Method Network Construction with BH-FDR)
  ├── Correlation Family: Spearman, iNAP, PBCDM, TRECOR, kLDM, MDiNE
  ├── Precision / Inverse Covariance Family: GLASSO, SPIEC-EASI, SparCC
  ├── Covariate / Latent Family: Covariate-Adjusted, Earth Geo, LUPINE, GNN, MIHN
  └── Stratified Subgroup Networks: City, Sex, Age, Infection, Underlying Disease
                       │
                       ▼ (S5: Paradigm Clustering & TOP-6 Selection)
  Cross-Method Jaccard Dissimilarity Dendrogram
  └──► Figure 06
                       │
                       ▼ (S6: 1,000-Iteration Bootstrap Stability & Consensus)
  ├── Subsampling (80% Cohort, 1,000 Reps) ──► Figure 07: Stability Ranking
  └── Weighted Multi-Paradigm Consensus   ──► Figure 08: Top-50 Consensus Network
                       │
         ┌─────────────┴─────────────┐
         ▼ (Story 1: Global Core)    ▼ (Story 2: Factor-Driven Rewiring)
(S7–S8: Keystone Concordance) (S10: Pathogen Plasticity Dashboard)
  Dual-Validation (Consensus vs. Robust)  Top-30 Volatile Pathogen Trajectories
  └──► Figure 09                         └──► Figure 11
(S9: Modularity & Hotspots)   (S10: Differential Edge Rewiring)
  Louvain + Pathogen Fisher's Exact Test  Pairwise Jaccard Edge Conservation
  └──► Figure 10 & Supp Fig S1           └──► Figure 12






