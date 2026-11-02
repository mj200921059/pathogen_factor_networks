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
- [Pipeline Architecture](#pipeline-architecture)
- [22 Network Inference Methods](#22-network-inference-methods)
- [Figure & Deliverables Index](#figure--deliverables-index)
- [Modular Figure Regeneration](#modular-figure-regeneration)
- [Performance & Hardware Considerations](#performance--hardware-considerations)
- [Validation Checklist](#validation-checklist)
- [Troubleshooting](#troubleshooting)
- [Citation & License](#citation--license)

---

## Workflow Overview

The pipeline spans **12 analytical stages** (S0 – S11) across a single master script, producing **16 main manuscript figures**, **3 supplementary figures**, 20+ summary CSV tables, and 24 serialized `igraph` network objects. Data flow proceeds from multi-center taxonomy harmonisation and ConQuR batch correction through ecological profiling, 22-method network construction, bootstrap consensus stabilisation, keystone and community analyses, differential rewiring, to Random Forest clinical phenotype prediction with DeLong ROC validation.

```
S0 ──┬── S1  Alpha Diversity                               → Fig 03
     ├── S2  Beta Diversity & PERMANOVA                     → Fig 04
     ├── S3  Composition Heatmaps                           → Fig 05
     └── S4  22-Method Network Inference (22 × .rds)
              └── S5/S6  Consistency, Stability & Consensus  → Fig 06–08
                          ├── S7/S8  Dual Keystone Validation → Fig 09
                          ├── S9    Communities & Hotspots    → Fig 10, Supp S1
                          └── S10   Differential Rewiring     → Fig 11–12
                                     └── S10-EXT  Severe Pneumonia → Fig 13–15, Supp S2
                                                   └── S11   ML Prediction & DeLong → Fig 16
                                                              └── S11-SUPP  4-Task Overview → Supp S3
```

S0 also produces **Fig 01** (prevalence calibration) and **Fig 02** (batch correction PCA).

---

## Repository Structure

```
.
├── LICENSE
├── README.md
├── pipeline_step0.R                    # Step 0: Bracken TSV → merged taxonomy table (kraken_taxa.rda)
├── run_balf_network_FINALpipeline.R    # Master pipeline: S0 → S11 (2,623 lines)
├── human_pathogens.csv                 # Curated pathogen reference list
├── kraken_taxa.rda                     # Pre-built merged taxonomy object
└── run_sampleid.csv                    # Sample metadata
```

### Output directories (created automatically at runtime)

```
outputs/
├── figures/            # PNG (300 DPI) + vector PDF for every figure
├── tables/             # Summary statistics & model results (CSV)
└── networks_species/   # Serialized igraph objects (.rds), one per method
```

---

## Prerequisites & Dependencies

### R packages

```r
install.packages(c(
  "tidyverse", "vegan", "igraph", "ggraph", "ggplot2",
  "parallel", "glasso", "cluster", "dendextend",
  "patchwork", "ggrepel", "scales", "RColorBrewer",
  "ranger", "pROC", "caret"
))

# Additionally for Step 0 (taxonomy resolution from raw Bracken output):
install.packages("taxonomizr")
```

### Verify installation

```r
sapply(c("tidyverse","vegan","igraph","ggraph","glasso",
         "cluster","dendextend","patchwork","ggrepel",
         "scales","RColorBrewer","ranger","pROC","caret"),
       \(pkg) paste(pkg, packageVersion(pkg)))
```

### Key hyperparameters (set at top of `run_balf_network_FINALpipeline.R`)

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `N_BOOT` | 1,000 | Bootstrap iterations for consensus network stability |
| `BOOT_FRAC` | 0.80 | Fraction of samples per bootstrap resample |
| `CORR_THRESHOLD` | 0.30 | Minimum \|r\| for edge inclusion |
| `FDR_THRESHOLD` | 0.05 | Benjamini–Hochberg FDR cutoff |
| `STABLE_THRESHOLD` | 0.80 | Minimum bootstrap frequency for consensus edges |
| `WEIGHTED_CONSENSUS_THRESH` | 0.40 | Weighted consensus edge retention threshold |
| `MIN_SAMPLES` | 15 | Minimum group size for stratified analyses |
| `NAMAP_DISCOUNT` | 0.70 | 30% phylogenetic discount for NAMAP method |
| `MIHN_MICRO_WEIGHT` / `MIHN_HOST_WEIGHT` | 0.60 / 0.40 | Heterogeneous network layer weights |
| `N_CORES` | auto-detect | Parallel cores (capped at `detectCores() - 1`, max 32) |
| `SEED` | 2026 | Global reproducibility seed |
| `LAYOUT_SEED` | 2016 | Network layout reproducibility seed |

All parameters are centralised in the first ~80 lines of the script.

---

## Input Data Requirements

### Required files (place in project root)

| File | Description | Format |
|------|-------------|--------|
| `kraken_taxa.rda` | Pre-built merged Kraken2/Bracken taxonomy abundance matrix (taxa × samples) | R binary |
| `run_sampleid.csv` | Sample metadata (run ID, bioproject, sex, age, clinical labels) | CSV |
| `human_pathogens.csv` | Curated human pathogen species reference list | CSV |

### Data structure

**`kraken_taxa.rda`** (loaded as `kraken_taxa`):
- Rows: taxa across all taxonomic levels
- Columns: annotation columns (`taxonomy_id`, `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`) + sample run IDs
- Values: raw read counts (Bracken re-estimated)

**`run_sampleid.csv`:**
- Columns: `runid`, `bioprojectid`, `seqtype`, `sampletypes`, `Sex`, `Age`, `Underlying_disease`, `infections`, etc.
- The pipeline filters to BALF mNGS samples only

**`human_pathogens.csv`:**
- Column 1: species names matching the taxonomy matrix
- Used for pathogen-enrichment, centrality, and hotspot analyses

### Geographic cohorts (10 multi-center)

| BioProject | City |
|------------|------|
| CNP0005812 | Ningbo |
| CNP0007066 | Jining |
| CRA006760 | Tianjin |
| PRJCA008730 | Guangzhou |
| PRJNA1175072 | Fuzhou |
| PRJNA655567 | Shenzhen |
| PRJNA808886 | Wenzhou |
| PRJNA845064 | Nanjing |
| PRJNA1083436 | Delhi |
| PRJNA1228250 | Changchun |

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/mj200921059/balf-mngs-network-pipeline.git
cd balf-mngs-network-pipeline
```

### 2. (Optional) Build taxonomy from raw Bracken output

Only needed if starting from raw `*_brackentaxa.tsv` files. A pre-built `kraken_taxa.rda` is included.

```bash
# Update working directory paths inside pipeline_step0.R, then:
Rscript pipeline_step0.R
```

### 3. Run the full pipeline

```bash
Rscript run_balf_network_FINALpipeline.R 2>&1 | tee run_log.txt
```

### Console output

```
======================================================================
Starting Complete BALF mNGS Network Analysis Pipeline (S0 to S11)
Platform: unix | Active Cores: 32 | Global Seed: 2026
======================================================================

[S0] Data loading & preprocessing ...
[S4] 22-method network construction (on CLR matrix) ...
  01 Spearman Baseline ...
  02 GLASSO ...
  ...
  22 SPIEC-EASI Variant ...
[S5 & S6] Consistency, stability & consensus ...
[S7 & S8] Dual keystone cross-validation ...
[S9] Community detection & pathogen hotspots ...
[S10] Differential rewiring trajectories ...
[S11] Machine learning & DeLong validation ...
```

---

## Pipeline Architecture

```
pipeline_step0.R (taxonomy builder)
│
└── Read *_brackentaxa.tsv → taxonomizr lineage resolution → kraken_taxa.rda


run_balf_network_FINALpipeline.R (2,623 lines)
│
├── S0    Data Loading, Harmonization & Batch Correction
│   ├── Load kraken_taxa.rda + run_sampleid.csv
│   ├── Filter: BALF mNGS only, valid age/sex
│   ├── Map 10 cohorts to geographic labels
│   ├── Prevalence-sensitivity calibration with pathogen retention trade-off
│   ├── ConQuR batch-effect correction
│   ├── CLR transformation
│   └── Output: Fig 01, Fig 02, metadata_processed.csv
│
├── S1    Alpha Diversity
│   ├── Observed richness, Shannon, Simpson
│   ├── Kruskal–Wallis tests across demographics
│   └── Output: Fig 03, alpha_diversity.csv, alpha_diversity_stats.csv
│
├── S2    Beta Diversity & PERMANOVA
│   ├── Bray–Curtis distance → PCoA ordination
│   ├── PERMANOVA across 4 clinical factors
│   └── Output: Fig 04, pcoa_coordinates.csv, permanova_results.csv
│
├── S3    Microbial Composition Heatmaps
│   ├── Top-20 species per factor (city, age, sex, infection)
│   └── Output: Fig 05
│
├── S4    22-Method Network Inference
│   ├── Correlation: Spearman, GLASSO
│   ├── Stratified: Location, Gender, Age, Underlying Disease, Infection
│   ├── Covariate-adjusted & Metadata-weighted
│   ├── Ecological: iNAP, NAMAP, Earth Biogeographic, Spatial Integration
│   ├── Heterogeneous: MIHN
│   ├── Temporal: TRECOR (age-trajectory sliding window)
│   ├── Probabilistic: MDiNE, kLDM, PBCDM, LUPINE
│   ├── Compositional: SparCC, SPIEC-EASI
│   ├── Graph neural network: GNN
│   └── Output: 22 × .rds igraph objects
│
├── S5–S6  Consistency, Stability & Consensus
│   ├── 22 × 22 Jaccard similarity matrix
│   ├── Hierarchical clustering → 3 paradigm families → TOP-6 selection
│   ├── 1,000-iteration bootstrap stability scoring (parallelised)
│   ├── Robust network + weighted consensus network construction
│   └── Output: Fig 06, Fig 07, Fig 08, consistency & stability CSVs
│
├── S7–S8  Dual Keystone Cross-Validation
│   ├── Degree, betweenness, closeness, eigenvector centrality
│   ├── Composite keystone score (KS)
│   ├── Dual-network concordance scatter
│   └── Output: Fig 09, keystone_species_dual_comparison.csv
│
├── S9    Community Detection & Pathogen Hotspots
│   ├── Louvain modularity on consensus + robust networks
│   ├── Fisher's exact pathogen enrichment per module
│   └── Output: Fig 10, Supp Fig S1, community CSVs
│
├── S10   Differential Rewiring & Topology
│   ├── Factor-stratified subnetworks (full species set)
│   ├── Edge Jaccard & rewiring percentage
│   ├── Pathogen centrality trajectories across conditions
│   └── Output: Fig 11, Fig 12
│
├── S10-EXT  Severe Pneumonia Sub-network & Hotspot Scoring
│   ├── Prevalence-filtered differential edges
│   ├── Patient-level Pathogen Hotspot Eigengene
│   ├── Logistic regression model + forest plot
│   └── Output: Fig 13, Fig 14, Fig 15, Supp Fig S2
│
├── S11   Integrated ML Classification & DeLong Validation
│   ├── Random Forest: Abundance, Network, and Integrated feature sets
│   ├── 5-fold stratified cross-validation
│   ├── ROC comparison with DeLong test
│   └── Output: Fig 16, predictive_delong_stats.csv
│
└── S11-SUPP  4-Task Phenotype Prediction
    ├── Infection detection, Severity, TB, Aging
    ├── Gini variable importance
    └── Output: Supp Fig S3, predictive_tasks_auc_summary.csv
```

---

## 22 Network Inference Methods

| # | Method | Paradigm | Edge Selection |
|:--|:-------|:---------|:---------------|
| 01 | Spearman Baseline | Correlation | FDR-corrected *t*-test |
| 02 | GLASSO | Graphical Lasso | Precision matrix thresholding |
| 03 | Stratified — Location | Stratified Correlation | Per-city, FDR-corrected |
| 04 | Stratified — Gender | Stratified Correlation | Per-sex, FDR-corrected |
| 05 | Stratified — Age | Stratified Correlation | Per-age-group, FDR-corrected |
| 06 | Covariate-Adjusted | Residual Correlation | Regress covariates, FDR |
| 07 | Metadata-Weighted Rank Covariance | Weighted Correlation | Weighted rank, FDR |
| 08 | iNAP | Compositional | Permutation null model |
| 09 | NAMAP | Compositional | 30% phylogenetic discount |
| 10 | Earth Biogeographic | Spatial | Geographic weighting |
| 11 | Spatial Integration | Spatial | Multi-site integration |
| 12 | MIHN Heterogeneous | Heterogeneous Network | Micro (0.6) + Host (0.4) layers |
| 13 | TRECOR | Temporal | Age-trajectory sliding window |
| 14 | MDiNE | Dirichlet-Multinomial | Latent precision |
| 15 | kLDM | Latent Dirichlet | Latent structure |
| 16 | PBCDM | Probabilistic | Bootstrap correlation |
| 17 | LUPINE | Latent Graph | Latent graph architecture |
| 18 | GNN | Graph Neural Network | Symmetrized graph convolution |
| 19 | Stratified — Underlying Disease | Stratified Correlation | Per-comorbidity, FDR |
| 20 | Stratified — Infection | Stratified Correlation | Per-infection-group, FDR |
| 21 | SparCC-Inspired | Compositional | Iterative log-ratio |
| 22 | SPIEC-EASI Variant | Compositional | Sparse inverse covariance |

---

## Figure & Deliverables Index

### Main manuscript figures (16, each exported as PNG 300 DPI + vector PDF)

| Figure | File Basename | Content |
|--------|---------------|---------|
| **Fig 01** | `fig01_prevalence_calibration_tradeoff` | Pathogen retention vs. sparsity trade-off across detection thresholds |
| **Fig 02** | `fig02_batch_correction_pca` | PCA ordination before & after ConQuR batch correction |
| **Fig 03** | `fig03_alpha_diversity_clean` | Alpha diversity (Observed, Shannon, Simpson) across demographics |
| **Fig 04** | `fig04_beta_diversity_pcoa_grouped` | Bray–Curtis PCoA with 95% confidence ellipses & PERMANOVA |
| **Fig 05** | `fig05_combined_top20_heatmap_fixed` | Top-20 species relative abundance heatmaps across clinical factors |
| **Fig 06** | `fig06_consistency_dendrogram` | Hierarchical clustering dendrogram of 22 network pipelines |
| **Fig 07** | `fig07_stability` | 1,000-iteration bootstrap stability scores for TOP-6 paradigms |
| **Fig 08** | `fig08_consensus_network` | High-confidence weighted consensus network of core lung ecology |
| **Fig 09** | `fig09_keystones_dual_publication` | Dual keystone species cross-validation & concordance scatter |
| **Fig 10** | `fig10_communities_dual_publication` | Louvain community modularity & Fisher's exact pathogen hotspot zoom |
| **Fig 11** | `fig11_pathogen_centrality_dashboard` | Pathogen centrality rewiring trajectories across clinical conditions |
| **Fig 12** | `fig12_differential_networks_publication` | Differential edge turnover & pairwise Jaccard topological conservation |
| **Fig 13** | `fig13_differential_edges_detail` | Severe pneumonia-specific sub-network of prevalence-supported edges |
| **Fig 14** | `fig14_patient_hotspot_scores` | Patient-level Pathogen Hotspot Eigengene distributions across phenotypes |
| **Fig 15** | `fig15_hotspot_roc_logistic` | Severity prediction ROC & host age dominance forest plot |
| **Fig 16** | `fig16_severe_pneumonia_rewiring_delong` | Severe pneumonia classification via edge rewiring & DeLong test (*p* = 0.0071) |

### Supplementary figures (3)

| ID | File Basename | Content |
|----|---------------|---------|
| **Supp S1** | `supp_fig01_robust_hotspot_zoom` | Robust network secondary pathogen hotspot module zoom-in |
| **Supp S2** | `supp_fig02_age_hotspot_interaction` | Age × Hotspot non-linear interaction model |
| **Supp S3** | `supp_fig03_multi_task_classification` | 4-task phenotype prediction + Gini variable importance |

### Data tables (CSV, written to `outputs/tables/`)

| File | Description | Stage |
|------|-------------|-------|
| `pathogen_prevalence_summary.csv` | Pathogen retention across prevalence thresholds | S0 |
| `metadata_processed.csv` | Cleaned metadata for all retained BALF mNGS samples | S0 |
| `batch_correction_qc.csv` | Silhouette scores before & after ConQuR correction | S0 |
| `alpha_diversity.csv` | Per-sample alpha diversity metrics | S1 |
| `alpha_diversity_stats.csv` | Kruskal–Wallis test results | S1 |
| `pcoa_coordinates.csv` | PCoA axis coordinates per sample | S2 |
| `permanova_results.csv` | PERMANOVA results for all factors | S2 |
| `consistency_jaccard.csv` | 22 × 22 pairwise Jaccard similarity matrix | S5 |
| `consistency_results.csv` | TOP-6 paradigms selected by intra-cluster similarity | S5 |
| `stability_all_methods.csv` | Bootstrap stability scores for TOP-6 methods | S6 |
| `keystone_species_dual_comparison.csv` | Dual-network keystone species rankings | S7/S8 |
| `community_enrichment.csv` | Pathogen enrichment per Louvain module | S9 |
| `community_modules.csv` | Community membership for all species | S9 |
| `story2_pathogen_factor_centrality.csv` | Pathogen centrality stratified by clinical factors | S10 |
| `story2_differential_networks.csv` | Edge turnover & Jaccard statistics per comparison | S10 |
| `differential_edges_severe_filtered.csv` | Prevalence-filtered severe pneumonia differential edges | S10-ext |
| `patient_hotspot_scores.csv` | Patient-level Pathogen Hotspot Eigengene scores | S10-ext |
| `logistic_multivariable_model.csv` | Multivariable logistic regression forest plot data | S10-ext |
| `predictive_delong_stats.csv` | DeLong ROC comparison statistics | S11 |
| `predictive_tasks_auc_summary.csv` | AUC summary for 4-task phenotype prediction | S11-supp |

Additionally, batch-corrected and CLR-transformed matrices are written to `outputs/`:
- `kraken_taxa_corrected.csv`
- `species_matrix_clr_corrected.csv`

### Network objects (22 methods + 2 consensus, `.rds` in `outputs/networks_species/`)

```
01_spearman_baseline.rds         12_mihn_heterogeneous.rds
02_glasso.rds                    13_trecor.rds
03_stratified_location.rds       14_mdine.rds
04_stratified_gender.rds         15_kldm.rds
05_stratified_age.rds            16_pbcdm.rds
06_covariate_adjusted.rds        17_lupine.rds
07_metadata_weighted.rds         18_gnn.rds
08_inap_inspired.rds             19_stratified_underlying_disease.rds
09_namap_inspired.rds            20_stratified_infection.rds
10_earth_biogeographic.rds       21_sparcc.rds
11_spatial_integration.rds       22_spiec_easi.rds

robust_network.rds               consensus_network.rds
```

---

## Modular Figure Regeneration

Comment out unwanted stages in the master script. Minimum dependency chains:

| Target | Required Stages |
|--------|-----------------|
| Alpha diversity only (Fig 03) | S0 → S1 |
| Beta diversity only (Fig 04) | S0 → S2 |
| Composition heatmaps (Fig 05) | S0 → S3 |
| Network construction only | S0 → S4 |
| Consensus + stability (Fig 06–08) | S0 → S4 → S5/S6 |
| Keystone analysis (Fig 09) | S0 → S4 → S5/S6 → S7/S8 |
| Community detection (Fig 10, Supp S1) | S0 → S4 → S5/S6 → S9 |
| Differential rewiring (Fig 11–12) | S0 → S4 → S5/S6 → S10 |
| Severe pneumonia (Fig 13–15, Supp S2) | S0 → S4 → S5/S6 → S10 → S10-ext |
| Clinical prediction (Fig 16, Supp S3) | S0 → S4 → S5/S6 → S10 → S10-ext → S11 |

---

## Performance & Hardware Considerations

### Tested hardware

| Component | Specification |
|-----------|---------------|
| Workstation | Dell Precision 7920 Tower |
| RAM | 320 GB DDR4 |
| CPU | Intel Xeon Gold 5220R × 96 logical cores |
| GPU | NVIDIA Quadro P400 (not required) |
| Storage | 24 TB |
| OS | Ubuntu 24.04.3 LTS |

### Minimum requirements

| Component | Minimum |
|-----------|---------|
| RAM | 32 GB (64 GB recommended) |
| CPU | 8 cores (more = faster parallelisation) |
| Disk | 5 GB free |
| OS | Linux (Ubuntu 20.04+), macOS 12+, or Windows 10/11 |
| R | ≥ 4.0.0 |

### Estimated runtime

| Cores | RAM | Estimated Time |
|-------|-----|----------------|
| 96 | 320 GB | ~20–40 min |
| 32 | 64 GB | ~1–2 hours |
| 8 | 32 GB | ~4–6 hours |

### Resource monitoring

```bash
htop
# or peak memory:
/usr/bin/time -v Rscript run_balf_network_FINALpipeline.R 2>&1 | tee run_log.txt
```

---

## Validation Checklist

After a successful run:

```bash
cd ~/balf-mngs-network-pipeline

# 1. All 22 method networks + 2 consensus networks
echo "Network .rds files:" && ls outputs/networks_species/*.rds | wc -l
# Expected: 24

# 2. All figures (16 main + 3 supplementary)
echo "Figure PNGs:" && ls outputs/figures/*.png 2>/dev/null | wc -l
echo "Figure PDFs:" && ls outputs/figures/*.pdf 2>/dev/null | wc -l

# 3. CSV table outputs
echo "CSV tables:" && ls outputs/tables/*.csv | wc -l

# 4. Corrected matrices
ls outputs/kraken_taxa_corrected.csv outputs/species_matrix_clr_corrected.csv

# 5. No errors in log
echo "Errors:" && grep -ci "error\|traceback" run_log.txt
```

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| `there is no package called 'ggraph'` | Missing dependency | `install.packages("ggraph")` (requires `libxml2-dev` on Ubuntu) |
| `Error in glasso(...)` | Ill-conditioned covariance | Normal for sparse data; handled via `tryCatch` fallback |
| Memory error during bootstrap (S6) | Insufficient RAM | Reduce `N_BOOT` to 500 or lower `N_CORES` |
| `PowerIterationFailedConvergence` | Disconnected graph | Handled by `tryCatch`; falls back to zero eigenvector centrality |
| Script hangs at ConQuR | Large matrix | Allow 5–15 min; reduce cores if memory-limited |
| `FileNotFoundError` for inputs | Wrong working directory | Ensure `kraken_taxa.rda`, `run_sampleid.csv`, `human_pathogens.csv` are in project root |
| PDF export fails | Missing cairo | Pipeline auto-detects `capabilities("cairo")`; falls back to base PDF device |
| BLAS thread contention | Parallel sockets | Already suppressed (`OMP_NUM_THREADS=1` at script top) |

---

## Citation & License

If you use this pipeline in your research, please cite:

> [Manuscript title]. [Authors]. [Journal]. [Year]. *(details forthcoming)*

Licensed under the [MIT License](LICENSE).

**Contact:** Jun Ma

---

*Pipeline version: 2.0 | Language: R | Script: 2,623 lines | Methods: 22 | Figures: 19*
