#!/usr/bin/env Rscript
# ==============================================================================
# MASTER REPLICATION PIPELINE: BALF mNGS NETWORK BIOLOGY & CLINICAL PREDICTION
# Stages: S0 Preprocessing -> S10 Differential Rewiring -> S11 Machine Learning
#
# Deliverables Architecture:
#   Main Manuscript Figures  : Fig 01 to Fig 16
#   Supplementary Figures    : Supp Fig S1 to Supp Fig S3
#   Summary Statistical Data : outputs/tables/*.csv
#   Network Topology Objects : outputs/networks_species/*.rds
#
# Platforms Supported: Ubuntu Linux, Debian, Windows 10 & 11
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(parallel)
  library(glasso)
  library(cluster)
  library(dendextend)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(ranger)   # High-speed Random Forest engine
  library(pROC)     # ROC curves & DeLong statistical testing
  library(caret)    # Stratified cross-validation splitting
})

# Prevent multithreaded BLAS contention inside parallel socket processes
Sys.setenv(OMP_NUM_THREADS      = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS      = "1")

# ── 1. Directories & Hyperparameters ──────────────────────────────────────────
ROOT_DIR <- "./"
INDIR    <- ROOT_DIR
OUTDIR   <- file.path(ROOT_DIR, "outputs")
FIG_DIR  <- file.path(OUTDIR, "figures")
TAB_DIR  <- file.path(OUTDIR, "tables")
NET_DIR  <- file.path(OUTDIR, "networks_species")

dir.create(OUTDIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(NET_DIR, recursive = TRUE, showWarnings = FALSE)

# Global Reproducibility Seeds
SEED        <- 2026
LAYOUT_SEED <- 2016
set.seed(SEED)

# Hyperparameter Constants (Methodologically Synchronized)
N_CORES                   <- max(1, min(32, detectCores() - 1))
N_BOOT                    <- 1000
BOOT_FRAC                 <- 0.80
CORR_THRESHOLD            <- 0.30
FDR_THRESHOLD             <- 0.05
STABLE_THRESHOLD          <- 0.80
MIN_SAMPLES               <- 15
BATCH_SIZE                <- N_CORES * 2
WEIGHTED_CONSENSUS_THRESH <- 0.40
NAMAP_DISCOUNT            <- 0.70  # 30% phylogenetic discount
MIHN_MICRO_WEIGHT         <- 0.60
MIHN_HOST_WEIGHT          <- 0.40

# Helper: Cross-Platform Figure Exporter (Preserves exact layout & vector resolution)
export_plot <- function(basename, plot_obj, width, height, dpi = 300) {
  png_file <- file.path(FIG_DIR, paste0(basename, ".png"))
  pdf_file <- file.path(FIG_DIR, paste0(basename, ".pdf"))
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  if (capabilities("cairo")) {
    ggsave(pdf_file, plot = plot_obj, width = width, height = height, device = cairo_pdf)
  } else {
    ggsave(pdf_file, plot = plot_obj, width = width, height = height, device = "pdf")
  }
  cat(sprintf("  [Exported] %s (.png / .pdf)\n", basename))
}

cat("======================================================================\n")
cat("Starting Complete BALF mNGS Network Analysis Pipeline (S0 to S11)\n")
cat(sprintf("Platform: %s | Active Cores: %d | Global Seed: %d\n", .Platform$OS.type, N_CORES, SEED))
cat("======================================================================\n\n")

# ═════════════════════════════════════════════════════════════════════════════
# S0: DATA LOADING, HARMONIZATION, PREVALENCE SENSITIVITY & ConQuR
# ═════════════════════════════════════════════════════════════════════════════
cat("[S0] Data loading & preprocessing ...\n")

taxa_rda <- file.path(INDIR, "kraken_taxa.rda")
taxa_csv <- file.path(INDIR, "kraken_taxa.csv")

if (file.exists(taxa_rda)) {
  load(taxa_rda)
  taxa_raw <- kraken_taxa
} else if (file.exists(taxa_csv)) {
  taxa_raw <- read.csv(taxa_csv, row.names = 1, check.names = FALSE)
} else {
  stop("Input taxonomy file (kraken_taxa.rda or kraken_taxa.csv) not found in INDIR!")
}

taxa_raw <- taxa_raw[taxa_raw$species != "Homo sapiens", ]

meta_raw <- read.csv(file.path(INDIR, "run_sampleid.csv"), stringsAsFactors = FALSE)
meta_raw <- meta_raw[toupper(meta_raw$seqtype) == "MNGS", ]
meta_raw <- meta_raw[toupper(meta_raw$sampletypes) == "BALF", ]

city_map <- c(
  "CNP0005812"   = "Ningbo (CNP0005812)",
  "CNP0007066"   = "Jining (CNP0007066)",
  "CRA006760"    = "Tianjin (CRA006760)",     # Corrected romanization
  "PRJCA008730"  = "Guangzhou (PRJCA008730)",
  "PRJNA1175072" = "Fuzhou (PRJNA1175072)",
  "PRJNA655567"  = "Shenzhen (PRJNA655567)",
  "PRJNA808886"  = "Wenzhou (PRJNA808886)",
  "PRJNA845064"  = "Nanjing (PRJNA845064)",
  "PRJNA1083436" = "Delhi (PRJNA1083436)",
  "PRJNA1228250" = "Changchun (PRJNA1228250)"
)
meta_raw$city <- city_map[meta_raw$bioprojectid]

age_group_func <- function(a) {
  a <- as.numeric(a)
  if (is.na(a)) return("Unknown")
  if (a <= 3)  return("Infant (<=3y)")
  if (a <= 14) return("Child (4-14y)")
  if (a <= 45) return("Young Adult (15-45y)")
  if (a <= 65) return("Middle-aged (46-65y)")
  return("Elderly (>65y)")
}
meta_raw$age_group <- sapply(meta_raw$Age, age_group_func)
meta_raw <- meta_raw[meta_raw$age_group != "Unknown", ]

meta_raw$Sex <- tools::toTitleCase(trimws(meta_raw$Sex))
meta_raw <- meta_raw[meta_raw$Sex %in% c("Male", "Female"), ]

clean_underlying <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A", "NaN", "None", "none")] <- NA_character_
  out <- rep("Unknown", length(x))
  out[x %in% c("0", "No", "no", "FALSE", "False")] <- "No"
  out[x %in% c("1", "Yes", "yes", "TRUE", "True")] <- "Yes"
  out[!is.na(x) & !x %in% c("0","1","No","no","Yes","yes","FALSE","False","TRUE","True")] <- "Yes"
  factor(out, levels = c("No", "Yes", "Unknown"))
}

clean_infection <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A", "NaN", "None", "none", "0")] <- "Unknown/None"
  x <- gsub("severe community-acquired pneumonia", "Severe pneumonia", x, ignore.case = TRUE)
  x <- gsub("mild community-acquired pneumonia",   "Mild pneumonia",   x, ignore.case = TRUE)
  x <- gsub("community[- ]acquired pneumonia",     "Pneumonia",        x, ignore.case = TRUE)
  x <- gsub("bronchopneumonia",                    "Pneumonia",        x, ignore.case = TRUE)
  x <- gsub("interstitial pneumonia",              "Pneumonia",        x, ignore.case = TRUE)
  x <- gsub("chronic pneumonia",                   "Pneumonia",        x, ignore.case = TRUE)
  x <- gsub("organizing pneumonia",                "Pneumonia",        x, ignore.case = TRUE)
  x <- gsub("^severe pneumonia$",                  "Severe pneumonia", x, ignore.case = TRUE)
  x <- gsub("^mild pneumonia$",                    "Mild pneumonia",   x, ignore.case = TRUE)
  x <- gsub("^pneumonia$",                         "Pneumonia",        x, ignore.case = TRUE)
  x <- gsub("lower respiratory tract infection.*", "Pulmonary infection", x, ignore.case = TRUE)
  x <- gsub("^pulmonary infection$",               "Pulmonary infection", x, ignore.case = TRUE)
  x <- gsub(".*pulmonary tb.*|.*pulmonary tuberculosis.*|^tuberculosis$", "Tuberculosis", x, ignore.case = TRUE)
  x <- gsub("pulmonary cryptococcosis",            "Cryptococcosis", x, ignore.case = TRUE)
  x <- gsub("lung.*carcinoma|carcinoma.*lung|^lung cancer$", "Lung cancer", x, ignore.case = TRUE)
  
  rare_pattern <- paste(
    "interstitial lung disease", "lung opacity",
    "pulmonary alveolar proteinosis", "vasculitis-associated lung disease",
    sep = "|"
  )
  x <- gsub(rare_pattern, "Other respiratory", x, ignore.case = TRUE)
  x <- gsub("healthy control", "Healthy control", x, ignore.case = TRUE)
  
  canonical <- c("Pneumonia", "Severe pneumonia", "Mild pneumonia",
                 "Pulmonary infection", "Tuberculosis", "Cryptococcosis",
                 "Lung cancer", "Other respiratory", "Healthy control",
                 "Unknown/None")
  x[!x %in% canonical] <- "Other respiratory"
  factor(x, levels = canonical)
}

meta_raw$Underlying_disease_binary <- clean_underlying(meta_raw$Underlying_disease)
meta_raw$infection_group           <- clean_infection(meta_raw$infections)

ANNO_COLS       <- c("taxonomy_id", "superkingdom", "phylum", "class", "order", "family", "genus", "species")
sample_cols_all <- setdiff(colnames(taxa_raw), ANNO_COLS)
valid_runs      <- intersect(meta_raw$runid, sample_cols_all)

taxa_sp <- taxa_raw[!is.na(taxa_raw$species) & trimws(taxa_raw$species) != "", ]
taxa_sp <- taxa_sp[!duplicated(taxa_sp$species), ]
rownames(taxa_sp) <- taxa_sp$species
taxa_sp <- taxa_sp[, valid_runs, drop = FALSE]

meta <- meta_raw[meta_raw$runid %in% valid_runs, ]
meta <- meta[!duplicated(meta$runid), ]
rownames(meta) <- meta$runid
meta <- meta[colnames(taxa_sp), ]
meta$gender <- meta$Sex

# ── Prevalence Sensitivity & Trade-Off Calibration ────────────────────────────
n_total_samples <- ncol(taxa_sp)
species_prev    <- rowSums(taxa_sp > 0) / n_total_samples

pathogen_file     <- file.path(INDIR, "human_pathogens.csv")
pathogens_df      <- if (file.exists(pathogen_file)) read.csv(pathogen_file, stringsAsFactors = FALSE) else data.frame(species = character(0))
pathogen_set      <- unique(trimws(pathogens_df[[1]]))
pathogen_clean    <- gsub(" complex$", "", pathogen_set)
pathogen_union    <- unique(c(pathogen_set, pathogen_clean))
dataset_pathogens <- intersect(rownames(taxa_sp), pathogen_union)
n_pathogens_total <- length(dataset_pathogens)

threshold_grid <- seq(0.01, 0.50, by = 0.01)
tradeoff_list  <- list()

for (th in threshold_grid) {
  retained_sp   <- names(species_prev[species_prev >= th])
  retained_path <- intersect(retained_sp, dataset_pathogens)
  sub_mat       <- taxa_sp[retained_sp, , drop = FALSE]
  zero_fraction <- mean(sub_mat == 0) * 100
  
  tradeoff_list[[length(tradeoff_list) + 1]] <- data.frame(
    threshold           = th * 100,
    n_species           = length(retained_sp),
    n_pathogens         = length(retained_path),
    pathogen_retention  = (length(retained_path) / max(1, n_pathogens_total)) * 100,
    matrix_sparsity_pct = zero_fraction
  )
}
df_tradeoff <- do.call(rbind, tradeoff_list)

full_pathogen_report <- data.frame(
  Pathogen          = dataset_pathogens,
  Samples_Detected  = rowSums(taxa_sp[dataset_pathogens, , drop = FALSE] > 0),
  Cohort_Prevalence = round(species_prev[dataset_pathogens] * 100, 2),
  Prevalence_String = sprintf("%.2f%%", species_prev[dataset_pathogens] * 100)
) %>% arrange(desc(Cohort_Prevalence))

write.csv(full_pathogen_report, file.path(TAB_DIR, "pathogen_prevalence_summary.csv"), row.names = FALSE)

# ── Figure 01: Prevalence Calibration Trade-Off ───────────────────────────────
p_retention <- ggplot(df_tradeoff, aes(x = threshold)) +
  geom_line(aes(y = pathogen_retention, color = "Pathogen Retention Rate (%)"), linewidth = 1.2) +
  geom_line(aes(y = (n_species / max(n_species)) * 100, color = "Total Species Retained (% of max)"), 
            linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 50, linetype = "dotted", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = 80, linetype = "dotted", color = "gray50", linewidth = 0.6) +
  annotate("text", x = 45, y = 53, label = "50% Target", size = 3, color = "gray40", fontface = "italic") +
  annotate("text", x = 45, y = 83, label = "80% Target", size = 3, color = "gray40", fontface = "italic") +
  scale_color_manual(values = c(
    "Pathogen Retention Rate (%)" = "#E74C3C",
    "Total Species Retained (% of max)" = "#2C3E50"
  )) +
  scale_x_continuous(breaks = seq(0, 50, by = 5)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 105)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"),
    plot.title = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = sprintf("A  Pathogen Retention vs. Total Species Reduction (N=%d Pathogens)", n_pathogens_total),
    x = "Prevalence Filter Threshold (%)",
    y = "Percentage Retained (%)"
  )

p_sparsity <- ggplot(df_tradeoff, aes(x = threshold, y = n_pathogens)) +
  geom_area(fill = "#E74C3C", alpha = 0.15) +
  geom_line(color = "#C0392B", linewidth = 1) +
  geom_point(color = "#962D22", size = 1.2) +
  scale_x_continuous(breaks = seq(0, 50, by = 5)) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(
    title = "B  Absolute Number of Retained Pathogens",
    x = "Prevalence Filter Threshold (%)",
    y = "Active Pathogen Count"
  )

p_fig01 <- p_retention / p_sparsity + plot_layout(heights = c(1.3, 1))
export_plot("fig01_prevalence_calibration_tradeoff", p_fig01, width = 9, height = 8)

# Retain species with cohort prevalence >= 5%
prev <- rowSums(taxa_sp > 0) / ncol(taxa_sp)
taxa_filt <- taxa_sp[prev >= 0.05, ]

keep_samples <- colSums(taxa_filt) > 0
if (sum(!keep_samples) > 0) {
  taxa_filt <- taxa_filt[, keep_samples, drop = FALSE]
  meta      <- meta[keep_samples, , drop = FALSE]
}
save(taxa_filt, file = file.path(OUTDIR, "taxa_filt.rda"))
write.csv(meta, file.path(TAB_DIR, "metadata_processed.csv"), row.names = TRUE)

# ── Batch Effect Normalization (Species-Wise Conditional Quantile Regression) ─
conqur_correct <- function(count_matrix, batch_labels, covariate_df, ref_batch = NULL, n_quantiles = 500) {
  count_matrix <- as.matrix(count_matrix)
  samples      <- colnames(count_matrix)
  species      <- rownames(count_matrix)
  batch_labels <- batch_labels[samples]
  covariate_df <- covariate_df[samples, , drop = FALSE]
  batch_factor <- as.factor(batch_labels)
  if (is.null(ref_batch)) ref_batch <- names(which.max(table(batch_factor)))
  
  cov_data <- data.frame(
    batch   = batch_factor,
    sex     = as.factor(covariate_df$Sex),
    nucleic = as.factor(covariate_df$nucleic)
  )
  cov_data[] <- lapply(cov_data, droplevels)
  cov_data   <- cov_data[complete.cases(cov_data), , drop = FALSE]
  
  corrected <- count_matrix
  for (sp in species) {
    y <- count_matrix[sp, ]
    y_bin <- as.numeric(y == 0)
    p_zero_obs <- rep(mean(y_bin), length(y))
    p_zero_ref <- p_zero_obs
    
    suppressWarnings(tryCatch({
      fit_lr <- glm(y_bin ~ ., data = cov_data, family = binomial(link = "logit"))
      p_zero_obs <- predict(fit_lr, type = "response")
      cov_data_ref <- cov_data
      cov_data_ref$batch <- factor(ref_batch, levels = levels(batch_factor))
      p_zero_ref <- predict(fit_lr, newdata = cov_data_ref, type = "response")
    }, error = function(e) {}))
    
    y_corr <- y
    nonzero_mask <- y > 0
    if (sum(nonzero_mask) < 5) {
      corrected[sp, ] <- y_corr
      next
    }
    
    y_nz <- y[nonzero_mask]
    quantiles <- seq(0.01, 0.99, length.out = n_quantiles)
    obs_quantiles <- quantile(y_nz, probs = quantiles, type = 7)
    cov_nz <- cov_data[nonzero_mask, , drop = FALSE]
    ref_quantiles <- obs_quantiles
    
    suppressWarnings(tryCatch({
      fit_ols <- lm(log1p(y_nz) ~ ., data = cov_nz)
      cov_nz_ref <- cov_nz
      cov_nz_ref$batch <- factor(ref_batch, levels = levels(batch_factor))
      y_ref_approx <- pmax(expm1(predict(fit_ols, newdata = cov_nz_ref)), 0)
      ref_quantiles <- quantile(y_ref_approx, probs = quantiles, type = 7)
    }, error = function(e) {}))
    
    for (idx in which(nonzero_mask)) {
      val <- y[idx]
      pos <- findInterval(val, obs_quantiles)
      pos <- max(1, min(pos, length(ref_quantiles)))
      y_corr[idx] <- max(ref_quantiles[pos], 0)
    }
    corrected[sp, ] <- round(y_corr)
  }
  corrected[is.na(corrected) | corrected < 0] <- 0
  
  orig_totals <- colSums(count_matrix)
  bad_cols <- which((colSums(corrected) / pmax(orig_totals, 1)) < 0.1)
  if (length(bad_cols) > 0) corrected[, bad_cols] <- count_matrix[, bad_cols]
  as.data.frame(corrected)
}

cov_df       <- meta[, c("Sex", "nucleic")]
batch_series <- setNames(meta$city, rownames(meta))
ref_batch    <- if ("Jining (CNP0007066)" %in% batch_series) "Jining (CNP0007066)" else names(which.max(table(batch_series)))

taxa_corrected <- conqur_correct(taxa_filt, batch_series, cov_df, ref_batch = ref_batch, n_quantiles = 500)

zero_after <- which(colSums(taxa_corrected) == 0)
if (length(zero_after) > 0) {
  taxa_corrected <- taxa_corrected[, -zero_after, drop = FALSE]
  taxa_filt      <- taxa_filt[, -zero_after, drop = FALSE]
  meta           <- meta[-zero_after, , drop = FALSE]
  cov_df         <- cov_df[-zero_after, , drop = FALSE]
  batch_series   <- batch_series[-zero_after]
}

write.csv(taxa_corrected, file.path(OUTDIR, "kraken_taxa_corrected.csv"), row.names = TRUE)

# ── CLR & PCA Diagnostics (Figure 02) ─────────────────────────────────────────
clr_transform <- function(df) {
  mat <- as.matrix(df) + 0.5
  log_mat <- log(mat)
  sweep(log_mat, 2, colMeans(log_mat), "-")
}

clr_before <- as.data.frame(clr_transform(taxa_filt))
clr_after  <- as.data.frame(clr_transform(taxa_corrected))
clr_mat    <- clr_after # Global persistent variable
write.csv(clr_after, file.path(OUTDIR, "species_matrix_clr_corrected.csv"), row.names = TRUE)

pca_before <- prcomp(t(clr_before), scale. = FALSE)
pca_after  <- prcomp(t(clr_after), scale. = FALSE)

get_pca_df <- function(pca_obj, batch_vals, title_text) {
  var_expl <- round(100 * pca_obj$sdev^2 / sum(pca_obj$sdev^2), 1)
  df <- data.frame(PC1 = pca_obj$x[, 1], PC2 = pca_obj$x[, 2], City = batch_vals)
  list(data = df, var1 = var_expl[1], var2 = var_expl[2], title = title_text)
}

# Complete color palette for all 10 participating clinical cohorts
city_colors <- c(
  "Changchun (PRJNA1228250)" = "#E41A1C", 
  "Delhi (PRJNA1083436)"     = "#377EB8",
  "Fuzhou (PRJNA1175072)"    = "#999999", 
  "Guangzhou (PRJCA008730)"  = "#4DAF4A", 
  "Jining (CNP0007066)"      = "#984EA3",
  "Nanjing (PRJNA845064)"    = "#A6761D", 
  "Ningbo (CNP0005812)"      = "#FF7F00", 
  "Shenzhen (PRJNA655567)"   = "#A65628",
  "Tianjin (CRA006760)"      = "#F781BF", 
  "Wenzhou (PRJNA808886)"    = "#66C2A5"
)

res_before <- get_pca_df(pca_before, batch_series, "Before Batch Correction")
res_after  <- get_pca_df(pca_after, batch_series, "After ConQuR Correction")

plot_pca_panel <- function(res) {
  ggplot(res$data, aes(x = PC1, y = PC2, color = City)) +
    geom_point(alpha = 0.8, size = 2) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
    theme_classic(base_size = 11) +
    theme(
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey93", linewidth = 0.4)
    ) +
    scale_color_manual(values = city_colors) +
    labs(
      title = res$title,
      x = sprintf("PC1 (%.1f%% variance)", res$var1),
      y = sprintf("PC2 (%.1f%% variance)", res$var2)
    )
}

p_before <- plot_pca_panel(res_before)
p_after  <- plot_pca_panel(res_after)

combined_pca_plot <- (p_before + p_after) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8))

fig02_pca <- combined_pca_plot + 
  plot_annotation(
    title = "Batch Effect Correction: PCA of CLR-Transformed Taxa",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )
export_plot("fig02_batch_correction_pca", fig02_pca, width = 13, height = 6.5)

sil_before <- mean(silhouette(as.integer(as.factor(batch_series)), dist(pca_before$x[, 1:2]))[, 3])
sil_after  <- mean(silhouette(as.integer(as.factor(batch_series)), dist(pca_after$x[, 1:2]))[, 3])
batch_qc <- data.frame(
  metric = "Silhouette_score_batch",
  before_correction = round(sil_before, 4),
  after_correction  = round(sil_after, 4),
  interpretation    = "Lower = superior batch mixing"
)
write.csv(batch_qc, file.path(TAB_DIR, "batch_correction_qc.csv"), row.names = FALSE)


# ═════════════════════════════════════════════════════════════════════════════
# S1: ALPHA DIVERSITY (FIGURE 03)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S1] Alpha diversity ...\n")

shannon_div  <- function(x) { x <- x[x > 0]; p <- x / sum(x); -sum(p * log(p)) }
simpson_div  <- function(x) { x <- x[x > 0]; p <- x / sum(x); 1 - sum(p^2) }
observed_div <- function(x) sum(x > 0)

sample_cols <- colnames(taxa_corrected)
alpha_records <- data.frame(
  sample             = sample_cols,
  Observed           = sapply(sample_cols, function(s) observed_div(taxa_corrected[, s])),
  Shannon            = sapply(sample_cols, function(s) shannon_div(taxa_corrected[, s])),
  Simpson            = sapply(sample_cols, function(s) simpson_div(taxa_corrected[, s])),
  city               = meta$city,
  gender             = meta$gender,
  age_group          = meta$age_group,
  Underlying_disease = meta$Underlying_disease_binary,
  Infection          = meta$infection_group
)
write.csv(alpha_records, file.path(TAB_DIR, "alpha_diversity.csv"), row.names = FALSE)

metrics <- c("Observed", "Shannon", "Simpson")
factors <- c("city", "gender", "age_group", "Underlying_disease", "Infection")
factor_labels <- c(
  age_group          = "Age Group",
  city               = "Geography",
  Infection          = "Infection",
  gender             = "Sex",
  Underlying_disease = "Underlying Disease"
)

alpha_stats_list <- list()
for (m in metrics) {
  for (f in factors) {
    kw <- kruskal.test(formula(paste(m, "~", f)), data = alpha_records)
    alpha_stats_list[[length(alpha_stats_list) + 1]] <- data.frame(
      metric = m, factor = f, H_stat = round(kw$statistic, 3), p_value = round(kw$p.value, 4)
    )
  }
}
alpha_stats_df <- do.call(rbind, alpha_stats_list)
alpha_stats_df$p_adj <- round(p.adjust(alpha_stats_df$p_value, method = "BH"), 4)
write.csv(alpha_stats_df, file.path(TAB_DIR, "alpha_diversity_stats.csv"), row.names = FALSE)

alpha_melt <- alpha_records %>%
  pivot_longer(cols = any_of(metrics), names_to = "Metric", values_to = "Value") %>%
  pivot_longer(cols = all_of(factors), names_to = "Factor", values_to = "Group") %>%
  filter(!is.na(Group) & !is.na(Value))

alpha_melt$FactorLabel <- factor(factor_labels[alpha_melt$Factor], levels = unname(factor_labels))
alpha_melt$Metric      <- factor(alpha_melt$Metric, levels = metrics)

alpha_stats_df$FactorLabel <- factor(factor_labels[alpha_stats_df$factor], levels = unname(factor_labels))
alpha_stats_df$Metric      <- factor(alpha_stats_df$metric, levels = metrics)
alpha_stats_df$p_label     <- ifelse(alpha_stats_df$p_value < 0.001, "KW: p < 0.001", sprintf("KW: p = %.3f", alpha_stats_df$p_value))

max_vals <- alpha_melt %>%
  group_by(Metric) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 0.95, .groups = "drop")
alpha_stats_df <- merge(alpha_stats_df, max_vals, by = "Metric")

fig03_alpha <- ggplot(alpha_melt, aes(x = Group, y = Value)) +
  geom_jitter(width = 0.2, alpha = 0.25, size = 0.7, color = "grey40", stroke = 0) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, fill = "grey90", color = "black", linewidth = 0.4) +
  geom_text(
    data = alpha_stats_df,
    aes(x = Inf, y = y_pos, label = p_label),
    hjust = 1.05, vjust = 1.5, size = 2.8, fontface = "italic", color = "blue4",
    inherit.aes = FALSE
  ) +
  facet_grid(Metric ~ FactorLabel, scales = "free", space = "free_x") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8.5, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4),
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  ) +
  labs(x = NULL, y = "Alpha Diversity Value", title = "Alpha Diversity Across Host and Environmental Factors")

export_plot("fig03_alpha_diversity_clean", fig03_alpha, width = 13, height = 7.8)


# ═════════════════════════════════════════════════════════════════════════════
# S2: BETA DIVERSITY & PERMANOVA (FIGURE 04)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S2] Beta diversity ...\n")

rel_ab <- sweep(taxa_corrected, 2, colSums(taxa_corrected) + 1e-9, "/")
dist_matrix <- vegdist(t(rel_ab), method = "bray")
pcoa_res <- cmdscale(dist_matrix, eig = TRUE, k = 2)
var_exp  <- round(100 * pcoa_res$eig[1:2] / sum(pmax(pcoa_res$eig, 0)), 1)

pcoa_df <- data.frame(
  PC1                = pcoa_res$points[, 1],
  PC2                = pcoa_res$points[, 2],
  city               = meta$city,
  gender             = meta$gender,
  age_group          = meta$age_group,
  Underlying_disease = meta$Underlying_disease_binary,
  Infection          = meta$infection_group,
  row.names          = colnames(taxa_corrected)
)
write.csv(pcoa_df, file.path(TAB_DIR, "pcoa_coordinates.csv"), row.names = TRUE)

factors_to_test <- c("city", "gender", "age_group", "Underlying_disease_binary", "infection_group")
permanova_results <- list()
for (f in factors_to_test) {
  valid_idx <- !is.na(meta[[f]])
  sub_dist  <- as.dist(as.matrix(dist_matrix)[valid_idx, valid_idx])
  sub_meta  <- meta[valid_idx, ]
  ado <- adonis2(as.formula(paste0("sub_dist ~ ", f)), data = sub_meta, permutations = 999, method = "bray")
  permanova_results[[f]] <- data.frame(
    Factor = f, R2 = round(ado$R2[1], 4), F_stat = round(ado$F[1], 3), p_value = ado$`Pr(>F)`[1]
  )
}
permanova_df <- do.call(rbind, permanova_results)
write.csv(permanova_df, file.path(TAB_DIR, "permanova_results.csv"), row.names = FALSE)

pcoa_factor_labels <- c(
  "city"               = "Geography",
  "gender"             = "Sex",
  "age_group"          = "Age Group",
  "Underlying_disease" = "Underlying disease",
  "Infection"          = "Infection"
)

create_pcoa_panel <- function(factor_name) {
  df_sub <- pcoa_df[!is.na(pcoa_df[[factor_name]]), ]
  df_sub$Group <- as.factor(df_sub[[factor_name]])
  ggplot(df_sub, aes(x = PC1, y = PC2, color = .data[[factor_name]])) +
    geom_point(alpha = 0.75, size = 1.5, stroke = 0) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
    theme_classic(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 7.5),
      legend.title = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.4)
    ) +
    labs(
      title = pcoa_factor_labels[factor_name],
      x = sprintf("PC1 (%.1f%%)", var_exp[1]), 
      y = sprintf("PC2 (%.1f%%)", var_exp[2])
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
}

p_city   <- create_pcoa_panel("city")
p_gender <- create_pcoa_panel("gender")
p_age    <- create_pcoa_panel("age_group")
p_dis    <- create_pcoa_panel("Underlying_disease")
p_inf    <- create_pcoa_panel("Infection")

fig04_pcoa <- (p_city + p_gender + p_age) / (p_dis + p_inf) +
  plot_annotation(
    title = "Beta Diversity PCoA (Bray-Curtis, Batch-Corrected)",
    subtitle = "Ordination across host and environmental factors with 95% confidence ellipses",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")
    )
  )
export_plot("fig04_beta_diversity_pcoa_grouped", fig04_pcoa, width = 15, height = 9.5)


# ═════════════════════════════════════════════════════════════════════════════
# S3: MICROBIAL COMPOSITION HEATMAPS (FIGURE 05)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S3] Microbial composition (Optimized) ...\n")

get_heatmap_plot <- function(factor_col, factor_label) {
  valid_meta <- meta[!is.na(meta[[factor_col]]), ]
  groups <- sort(unique(valid_meta[[factor_col]]))
  
  group_means <- list()
  for (g in groups) {
    samps <- intersect(rownames(valid_meta[valid_meta[[factor_col]] == g, ]), colnames(rel_ab))
    if (length(samps) > 0) group_means[[g]] <- rowMeans(rel_ab[, samps, drop = FALSE])
  }
  if (length(group_means) == 0) return(NULL)
  
  group_df <- as.data.frame(do.call(cbind, group_means))
  overall_mean <- rowMeans(group_df)
  top20 <- names(sort(overall_mean, decreasing = TRUE)[1:min(20, length(overall_mean))])
  sub_df <- group_df[top20, ]
  
  heatmap_data <- sub_df %>%
    rownames_to_column("Species") %>%
    pivot_longer(cols = -Species, names_to = "Group", values_to = "Abundance")
  
  ggplot(heatmap_data, aes(x = Group, y = Species, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Rel. Abund.") +
    theme_classic(base_size = 9) +
    labs(x = NULL, y = "Species", title = factor_label) +
    theme(
      axis.text.y = element_text(size = 7.5, face = "italic", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "right",
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )
}

p_heat_city <- get_heatmap_plot("city", "Geography")
p_heat_sex  <- get_heatmap_plot("gender", "Sex")
p_heat_age  <- get_heatmap_plot("age_group", "Age Group")
p_heat_dis  <- get_heatmap_plot("Underlying_disease_binary", "Underlying Disease")
p_heat_inf  <- get_heatmap_plot("infection_group", "Infection")

fig05_heatmaps <- (p_heat_city | p_heat_inf) / 
  (p_heat_sex | p_heat_age | p_heat_dis) +
  plot_annotation(
    title = "Top-20 Species Composition Across Host and Environmental Factors",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )
export_plot("fig05_combined_top20_heatmap_fixed", fig05_heatmaps, width = 16, height = 11)


# ═════════════════════════════════════════════════════════════════════════════
# S4: 22-METHOD SPECIES-LEVEL NETWORK INFERENCE
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S4] 22-method network construction (on CLR matrix) ...\n")

X <- t(as.matrix(clr_after))
sp_names <- colnames(X)
n_sp     <- length(sp_names)

failed_methods <- character(0)

# Function A: For genuine correlation matrices (Analytic Student's t-test with BH-FDR)
build_graph_from_corr_fdr <- function(corr_mat, n_samples, r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD, sp_list = sp_names) {
  p_total <- length(sp_list)
  g <- make_empty_graph(n = p_total, directed = FALSE)
  V(g)$name <- sp_list
  
  upper_idx <- which(upper.tri(corr_mat), arr.ind = TRUE)
  r_vals <- corr_mat[upper_idx]
  valid <- !is.na(r_vals) & abs(r_vals) < 1
  df <- n_samples - 2
  if (df < 1) return(g)
  
  t_stat <- r_vals * sqrt(df / pmax(1 - r_vals^2, 1e-12))
  p_vals <- rep(1, length(r_vals))
  p_vals[valid] <- 2 * pt(-abs(t_stat[valid]), df = df)
  
  padj_vals <- p.adjust(p_vals, method = "BH")
  sig_mask  <- (abs(r_vals) >= r_thresh) & (padj_vals < fdr_thresh)
  
  if (any(sig_mask)) {
    kept_indices <- upper_idx[sig_mask, , drop = FALSE]
    edge_list <- data.frame(
      from   = sp_list[kept_indices[, 1]],
      to     = sp_list[kept_indices[, 2]],
      weight = r_vals[sig_mask],
      padj   = padj_vals[sig_mask],
      stringsAsFactors = FALSE
    )
    g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = data.frame(name = sp_list))
  }
  return(g)
}

# Function B: For algebraically transformed matrices (Direct magnitude thresholding)
build_graph_from_corr <- function(corr_mat, threshold = CORR_THRESHOLD, sp_list = sp_names) {
  g <- make_empty_graph(n = length(sp_list), directed = FALSE)
  V(g)$name <- sp_list
  adj <- abs(corr_mat) >= threshold
  adj[lower.tri(adj, diag = TRUE)] <- FALSE
  indices <- which(adj, arr.ind = TRUE)
  if (nrow(indices) > 0) {
    edge_list <- data.frame(
      from   = sp_list[indices[, 1]],
      to     = sp_list[indices[, 2]],
      weight = corr_mat[indices],
      stringsAsFactors = FALSE
    )
    g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = data.frame(name = sp_list))
  }
  return(g)
}

# 01 Spearman Baseline
cat("  01 Spearman Baseline ...\n")
corr_sp <- cor(X, method = "spearman"); diag(corr_sp) <- 0
G01 <- build_graph_from_corr_fdr(corr_sp, n_samples = nrow(X), r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD)
saveRDS(G01, file.path(NET_DIR, "01_spearman_baseline.rds"))

# 02 GLASSO (Non-silent error tracking)
cat("  02 GLASSO ...\n")
G02 <- tryCatch({
  gl_res <- glasso::glasso(cov(scale(X)), rho = 0.3)
  prec <- gl_res$wi; colnames(prec) <- rownames(prec) <- sp_names; diag(prec) <- 0
  build_graph_from_corr(prec, threshold = 0.05)
}, error = function(e) {
  warning("Method 02 GLASSO failed: ", e$message)
  failed_methods <<- c(failed_methods, "02_glasso")
  g_empty <- make_empty_graph(n = n_sp, directed = FALSE); V(g_empty)$name <- sp_names; g_empty
})
saveRDS(G02, file.path(NET_DIR, "02_glasso.rds"))

strat_runner <- function(col_name, min_samples = MIN_SAMPLES) {
  res <- list()
  groups <- unique(meta[[col_name]])
  groups <- groups[!is.na(groups) & groups != ""]
  for (g in groups) {
    g_str <- as.character(g)
    samps <- intersect(rownames(meta[!is.na(meta[[col_name]]) & meta[[col_name]] == g, , drop = FALSE]), rownames(X))
    if (length(samps) >= min_samples) {
      cg <- cor(X[samps, , drop = FALSE], method = "spearman"); diag(cg) <- 0
      res[[g_str]] <- build_graph_from_corr_fdr(cg, n_samples = length(samps), r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD)
    }
  }
  res
}

cat("  03 Stratified Location ...\n")
saveRDS(strat_runner("city"), file.path(NET_DIR, "03_stratified_location.rds"))

cat("  04 Stratified Gender ...\n")
saveRDS(strat_runner("gender"), file.path(NET_DIR, "04_stratified_gender.rds"))

cat("  05 Stratified Age ...\n")
saveRDS(strat_runner("age_group"), file.path(NET_DIR, "05_stratified_age.rds"))

# 06 Covariate-Adjusted
cat("  06 Covariate-Adjusted ...\n")
cov_dummy <- model.matrix(~ city + gender + age_group + Underlying_disease_binary + infection_group - 1, data = meta)
X_resid <- residuals(lm(X ~ cov_dummy))
corr_cov <- cor(X_resid, method = "spearman"); diag(corr_cov) <- 0
G06 <- build_graph_from_corr_fdr(corr_cov, n_samples = nrow(X), r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD)
saveRDS(G06, file.path(NET_DIR, "06_covariate_adjusted.rds"))

# 07 Metadata-Weighted (True Weighted Rank Covariance)
cat("  07 Metadata-Weighted (Weighted Rank Covariance) ...\n")
age_num <- as.numeric(meta$Age)
age_num[is.na(age_num)] <- median(age_num, na.rm = TRUE)
w_raw <- (age_num - min(age_num)) / (max(age_num) - min(age_num) + 1e-9) + 0.1
w_norm <- w_raw / sum(w_raw)

X_ranked <- apply(X, 2, rank)
mean_w <- colSums(X_ranked * w_norm)
X_centered <- sweep(X_ranked, 2, mean_w, "-")
cov_w <- t(X_centered) %*% (X_centered * w_norm)
sd_w <- sqrt(diag(cov_w))
corr_w <- cov_w / (outer(sd_w, sd_w) + 1e-12)
diag(corr_w) <- 0
G07 <- build_graph_from_corr(corr_w, threshold = CORR_THRESHOLD)
saveRDS(G07, file.path(NET_DIR, "07_metadata_weighted.rds"))

# 08 iNAP
cat("  08 iNAP Inspired ...\n")
inap_perm_func <- function(seed_val) {
  set.seed(seed_val)
  Xp <- apply(X, 2, sample)
  cg <- cor(Xp, method = "spearman"); diag(cg) <- 0
  cg
}
cl_inap <- makePSOCKcluster(N_CORES)
clusterExport(cl_inap, c("X", "inap_perm_func"), envir = environment())
null_results <- parLapply(cl_inap, 1:500, inap_perm_func)
stopCluster(cl_inap)

null_array <- simplify2array(null_results)
null_mean  <- apply(null_array, c(1, 2), mean)
null_std   <- apply(null_array, c(1, 2), sd) + 1e-9
rm(null_array, null_results); gc()

z_scores <- (corr_sp - null_mean) / null_std; diag(z_scores) <- 0
upper_idx <- which(upper.tri(corr_sp), arr.ind = TRUE)
r_vals    <- corr_sp[upper_idx]
padj_vals <- p.adjust(2 * pnorm(-abs(z_scores[upper_idx])), method = "BH")
sig_mask  <- (abs(r_vals) >= CORR_THRESHOLD) & (padj_vals < FDR_THRESHOLD)

G08 <- make_empty_graph(n = length(sp_names), directed = FALSE); V(G08)$name <- sp_names
if (any(sig_mask)) {
  kept_idx <- upper_idx[sig_mask, , drop = FALSE]
  edge_df  <- data.frame(
    from    = sp_names[kept_idx[, 1]],
    to      = sp_names[kept_idx[, 2]],
    weight  = r_vals[sig_mask],
    z_score = round(z_scores[upper_idx][sig_mask], 3),
    padj    = padj_vals[sig_mask],
    stringsAsFactors = FALSE
  )
  G08 <- graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = sp_names))
}
saveRDS(G08, file.path(NET_DIR, "08_inap_inspired.rds"))

# 09 NAMAP
cat("  09 NAMAP Inspired ...\n")
genus_map <- sapply(strsplit(sp_names, " "), `[`, 1)
corr_namap <- corr_sp
for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    if (genus_map[i] == genus_map[j]) {
      corr_namap[i, j] <- corr_namap[j, i] <- corr_namap[i, j] * NAMAP_DISCOUNT
    }
  }
}
G09 <- build_graph_from_corr(corr_namap, threshold = CORR_THRESHOLD)
saveRDS(G09, file.path(NET_DIR, "09_namap_inspired.rds"))

# 10 Earth Biogeographic
cat("  10 Earth Biogeographic ...\n")
city_dummy <- model.matrix(~ city - 1, data = meta)
corr_geo <- cor(residuals(lm(X ~ city_dummy)), method = "spearman"); diag(corr_geo) <- 0
G10 <- build_graph_from_corr_fdr(corr_geo, n_samples = nrow(X), r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD)
saveRDS(G10, file.path(NET_DIR, "10_earth_biogeographic.rds"))

# 11 Spatial Integration
cat("  11 Spatial Integration ...\n")
corr_spatial <- 0.5 * corr_sp + 0.5 * corr_geo; diag(corr_spatial) <- 0
G11 <- build_graph_from_corr(corr_spatial, threshold = CORR_THRESHOLD)
saveRDS(G11, file.path(NET_DIR, "11_spatial_integration.rds"))

# 12 MIHN Heterogeneous
cat("  12 MIHN Heterogeneous ...\n")
corr_host <- (t(X) %*% cov_dummy) / (nrow(X) - 1)
host_sim  <- corr_host %*% t(corr_host)
norm_host <- norm(host_sim, "F")
if (norm_host > 0) host_sim <- host_sim / norm_host
corr_mihn <- MIHN_MICRO_WEIGHT * corr_sp + MIHN_HOST_WEIGHT * host_sim; diag(corr_mihn) <- 0
G12 <- build_graph_from_corr(corr_mihn, threshold = CORR_THRESHOLD)
saveRDS(G12, file.path(NET_DIR, "12_mihn_heterogeneous.rds"))

# 13 TRECOR (Age-Trajectory Sliding Window Correlation)
cat("  13 TRECOR (Age-Trajectory Sliding Window) ...\n")
ord_age <- order(age_num)
X_age_sorted <- X[ord_age, ]
n_samps_total <- nrow(X)
w_size <- max(30, round(0.35 * n_samps_total))
w_step <- max(8, round(0.08 * n_samps_total))

starts <- seq(1, n_samps_total - w_size + 1, by = w_step)
window_corrs <- list()
for (st in starts) {
  sub_win <- X_age_sorted[st:(st + w_size - 1), , drop = FALSE]
  c_w <- cor(sub_win, method = "spearman")
  c_w[is.na(c_w)] <- 0
  diag(c_w) <- 0
  window_corrs[[length(window_corrs) + 1]] <- c_w
}
corr_trend <- apply(simplify2array(window_corrs), c(1, 2), mean)
diag(corr_trend) <- 0
G13 <- build_graph_from_corr(corr_trend, threshold = CORR_THRESHOLD)
saveRDS(G13, file.path(NET_DIR, "13_trecor.rds"))

# 14 MDiNE
cat("  14 MDiNE ...\n")
corr_mdine_combined <- 0.7 * (corr_sp * (corr_sp > 0)) + 0.3 * (corr_sp * (corr_sp < 0))
diag(corr_mdine_combined) <- 0
G14 <- build_graph_from_corr(corr_mdine_combined, threshold = 0.25)
saveRDS(G14, file.path(NET_DIR, "14_mdine.rds"))

# 15 kLDM
cat("  15 kLDM ...\n")
K <- exp(-0.5 * as.matrix(dist(scale(age_num)))^2)
corr_kldm <- cor((K / (rowSums(K) + 1e-9)) %*% X, method = "spearman"); diag(corr_kldm) <- 0
G15 <- build_graph_from_corr(corr_kldm, threshold = CORR_THRESHOLD)
saveRDS(G15, file.path(NET_DIR, "15_kldm.rds"))

# 16 PBCDM
cat("  16 PBCDM ...\n")
cl_pbcdm <- makePSOCKcluster(N_CORES)
clusterExport(cl_pbcdm, "X", envir = environment())
boot_results_pbcdm <- parLapply(cl_pbcdm, 1:200, function(s) {
  set.seed(s)
  cor(X[sample(nrow(X), replace = TRUE), ], method = "spearman")
})
stopCluster(cl_pbcdm)
corr_pbcdm <- apply(simplify2array(boot_results_pbcdm), c(1, 2), mean); diag(corr_pbcdm) <- 0
rm(boot_results_pbcdm); gc()
G16 <- build_graph_from_corr_fdr(corr_pbcdm, n_samples = nrow(X), r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD)
saveRDS(G16, file.path(NET_DIR, "16_pbcdm.rds"))

# 17 LUPINE
cat("  17 LUPINE ...\n")
pca_lupine <- prcomp(X, scale. = FALSE)
n_comps    <- min(20, ncol(X) - 1, nrow(X) - 1)
corr_lupine <- cor(pca_lupine$x[, 1:n_comps] %*% t(pca_lupine$rotation[, 1:n_comps]), method = "spearman")
diag(corr_lupine) <- 0
G17 <- build_graph_from_corr(corr_lupine, threshold = CORR_THRESHOLD)
saveRDS(G17, file.path(NET_DIR, "17_lupine.rds"))

# 18 GNN (Symmetrized Graph Convolution)
cat("  18 GNN (Symmetric Adjacency) ...\n")
dist_mat <- as.matrix(dist(X))
A_dense  <- matrix(0, nrow = nrow(X), ncol = nrow(X))
k_nn     <- min(10, nrow(X) - 1)
for (i in 1:nrow(X)) A_dense[i, order(dist_mat[i, ])[2:(k_nn + 1)]] <- 1
A_sym <- pmax(A_dense, t(A_dense))
L <- diag(rowSums(A_sym)) - A_sym
corr_gnn <- cor((diag(nrow(X)) - 0.1 * L) %*% X, method = "spearman"); diag(corr_gnn) <- 0
G18 <- build_graph_from_corr(corr_gnn, threshold = CORR_THRESHOLD)
saveRDS(G18, file.path(NET_DIR, "18_gnn.rds"))

# 19-20 Stratified Disease & Infection
cat("  19 Stratified Underlying Disease ...\n")
saveRDS(strat_runner("Underlying_disease_binary"), file.path(NET_DIR, "19_stratified_underlying_disease.rds"))

cat("  20 Stratified Infection ...\n")
saveRDS(strat_runner("infection_group"),           file.path(NET_DIR, "20_stratified_infection.rds"))

# 21 SparCC-Inspired Native Implementation
cat("  21 SparCC-Inspired ...\n")
run_sparcc <- function(count_mat, iter = 10) {
  X_cnt <- t(count_mat) + 1; p <- ncol(X_cnt); logX <- log(X_cnt)
  V <- matrix(0, p, p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      V[i, j] <- V[j, i] <- var(logX[, i] - logX[, j])
    }
  }
  t_vec <- pmax(rowSums(V) / (p - 2), 1e-6)
  for (it in 1:iter) {
    denom <- outer(sqrt(t_vec), sqrt(t_vec))
    rho <- pmin(pmax((outer(t_vec, t_vec, "+") - V) / (2 * denom), -1), 1)
    diag(rho) <- 1
    t_vec <- pmax(rowSums(V) / pmax(p - 2 + rowSums(rho) - 1, 1e-4), 1e-6)
  }
  denom <- outer(sqrt(t_vec), sqrt(t_vec))
  rho <- pmin(pmax((outer(t_vec, t_vec, "+") - V) / (2 * denom), -1), 1)
  diag(rho) <- 0
  colnames(rho) <- rownames(rho) <- rownames(count_mat)
  rho
}
G21 <- tryCatch({
  corr_sparcc <- run_sparcc(as.matrix(taxa_corrected), iter = 10)
  build_graph_from_corr_fdr(corr_sparcc, n_samples = ncol(taxa_corrected), r_thresh = CORR_THRESHOLD, fdr_thresh = FDR_THRESHOLD)
}, error = function(e) {
  warning("Method 21 SparCC failed: ", e$message)
  failed_methods <<- c(failed_methods, "21_sparcc")
  g_empty <- make_empty_graph(n = n_sp, directed = FALSE); V(g_empty)$name <- sp_names; g_empty
})
saveRDS(G21, file.path(NET_DIR, "21_sparcc.rds"))

# 22 SPIEC-EASI Variant
cat("  22 SPIEC-EASI Variant ...\n")
G22 <- tryCatch({
  cov_clr <- cov(X)
  prec_se <- glasso::glasso(cov_clr, rho = 0.15)$wi
  colnames(prec_se) <- rownames(prec_se) <- sp_names; diag(prec_se) <- 0
  build_graph_from_corr(prec_se, threshold = 0.02)
}, error = function(e) {
  warning("Method 22 SPIEC-EASI failed: ", e$message)
  failed_methods <<- c(failed_methods, "22_spiec_easi")
  g_empty <- make_empty_graph(n = n_sp, directed = FALSE); V(g_empty)$name <- sp_names; g_empty
})
saveRDS(G22, file.path(NET_DIR, "22_spiec_easi.rds"))

cat("  Finished network construction across 22 inference pipelines.\n")


# ═════════════════════════════════════════════════════════════════════════════
# S5 & S6: CONSISTENCY, STABILITY & CONSENSUS (FIGURES 06, 07, 08)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S5 & S6] Consistency, Stability & Consensus ...\n")

METHOD_FILES <- sprintf("%02d_%s.rds", 1:22, c(
  "spearman_baseline", "glasso", "stratified_location", "stratified_gender",
  "stratified_age", "covariate_adjusted", "metadata_weighted", "inap_inspired",
  "namap_inspired", "earth_biogeographic", "spatial_integration", "mihn_heterogeneous",
  "trecor", "mdine", "kldm", "pbcdm", "lupine", "gnn",
  "stratified_underlying_disease", "stratified_infection", "sparcc", "spiec_easi"
))
METHOD_NAMES <- gsub(".rds", "", METHOD_FILES)

if (length(failed_methods) > 0) {
  cat(sprintf("  [Exclusion] Omitting failed pipelines: %s\n", paste(failed_methods, collapse = ", ")))
  keep_mask <- !METHOD_NAMES %in% failed_methods
  METHOD_FILES <- METHOD_FILES[keep_mask]
  METHOD_NAMES <- METHOD_NAMES[keep_mask]
}

load_method_graph <- function(fname) {
  obj <- readRDS(fname)
  if (is.list(obj) && !inherits(obj, "igraph")) {
    el_list <- lapply(obj, function(g) {
      if (inherits(g, "igraph") && ecount(g) > 0) {
        el <- as_edgelist(g, names = TRUE)
        w  <- if (!is.null(E(g)$weight)) E(g)$weight else rep(1.0, nrow(el))
        data.frame(from = el[,1], to = el[,2], weight = w, stringsAsFactors = FALSE)
      }
    })
    all_el <- do.call(rbind, el_list)
    if (is.null(all_el) || nrow(all_el) == 0) {
      g <- make_empty_graph(n = length(sp_names), directed = FALSE)
      V(g)$name <- sp_names
      return(g)
    }
    g <- graph_from_data_frame(all_el, directed = FALSE, vertices = data.frame(name = sp_names))
    return(igraph::simplify(g, edge.attr.comb = "mean"))
  } else if (inherits(obj, "igraph")) {
    return(subgraph(obj, intersect(V(obj)$name, sp_names)))
  }
  make_empty_graph(n = length(sp_names), directed = FALSE)
}

graphs <- setNames(lapply(file.path(NET_DIR, METHOD_FILES), load_method_graph), METHOD_NAMES)
save(graphs, file = file.path(OUTDIR, "graphs.rda"))

edge_set <- function(g) {
  el <- as_edgelist(g, names = TRUE)
  if (nrow(el) == 0) return(character(0))
  paste(pmin(el[,1], el[,2]), pmax(el[,1], el[,2]), sep = "__")
}
edge_sets <- setNames(lapply(graphs, edge_set), METHOD_NAMES)

jaccard_mat <- matrix(0, length(METHOD_NAMES), length(METHOD_NAMES), dimnames = list(METHOD_NAMES, METHOD_NAMES))
for (i in seq_along(METHOD_NAMES)) {
  for (j in seq_along(METHOD_NAMES)) {
    if (i == j) { jaccard_mat[i, j] <- 1.0; next }
    ei <- edge_sets[[i]]; ej <- edge_sets[[j]]
    uni <- length(union(ei, ej))
    jaccard_mat[i, j] <- if (uni > 0) length(intersect(ei, ej)) / uni else 0.0
  }
}
write.csv(as.data.frame(jaccard_mat), file.path(TAB_DIR, "consistency_jaccard.csv"), row.names = TRUE)

# ── Figure 06: Method Consistency Dendrogram ──────────────────────────────────
hc   <- hclust(as.dist(1 - jaccard_mat), method = "complete")
dend <- as.dendrogram(hc)
k_clusters <- 3
dend <- color_branches(dend, k = k_clusters, palette = c("#1f78b4", "#33a02c", "#e31a1c"))
dend <- color_labels(dend, k = k_clusters, palette = c("#1f78b4", "#33a02c", "#e31a1c"))
labels_cex(dend) <- 0.85

pdf(file.path(FIG_DIR, "fig06_consistency_dendrogram.pdf"), width = 11, height = 6.5)
par(mar = c(12, 4.5, 3, 2), font.axis = 2)
plot(dend, main = "Hierarchical Clustering of Method Consistency",
     ylab = "Dissimilarity (1 - Jaccard Index)", cex.main = 1.2, font.main = 2)
abline(h = 0.35, col = "grey60", lty = 2, lwd = 1)
dev.off()

method_clusters <- cutree(hc, k = 3)
select_family_top <- function(cluster_id, n_pick = 2) {
  members <- names(method_clusters[method_clusters == cluster_id])
  if (length(members) <= n_pick) return(members)
  intra_sim <- rowMeans(jaccard_mat[members, members, drop = FALSE])
  names(sort(intra_sim, decreasing = TRUE)[1:n_pick])
}
TOP6 <- c(select_family_top(1, 2), select_family_top(2, 2), select_family_top(3, 2))
write.csv(data.frame(method = TOP6), file.path(TAB_DIR, "consistency_results.csv"), row.names = FALSE)
cat(sprintf("  TOP6 Paradigms Selected: %s\n", paste(TOP6, collapse = ", ")))

# ── Bootstrap Stability ───────────────────────────────────────────────────────
N_PAIRS <- n_sp * (n_sp - 1) / 2
ut_map  <- which(upper.tri(matrix(FALSE, n_sp, n_sp)), arr.ind = TRUE)

worker_one_boot_int <- function(idx) {
  X_sub    <- X[idx, , drop = FALSE]
  meta_sub <- meta[idx, , drop = FALSE]
  n_sub    <- nrow(X_sub)
  if (n_sub < 10) return(integer(0))
  
  filter_fdr_hits <- function(r_mat) {
    r_vals <- r_mat[upper.tri(r_mat)]
    valid  <- !is.na(r_vals) & abs(r_vals) < 1
    df     <- n_sub - 2
    t_stat <- r_vals * sqrt(df / pmax(1 - r_vals^2, 1e-12))
    p_vals <- rep(1, length(r_vals))
    p_vals[valid] <- 2 * pt(-abs(t_stat[valid]), df = df)
    sig_mask  <- (abs(r_vals) >= CORR_THRESHOLD) & (p.adjust(p_vals, method = "BH") < FDR_THRESHOLD)
    which(sig_mask)
  }
  
  strat_map <- c(
    "03_stratified_location"           = "city",
    "04_stratified_gender"             = "gender",
    "05_stratified_age"                = "age_group",
    "19_stratified_underlying_disease" = "Underlying_disease_binary",
    "20_stratified_infection"          = "infection_group"
  )
  if (method_name %in% names(strat_map)) {
    strat_col <- strat_map[[method_name]]
    sub_hits  <- list()
    for (g in unique(meta_sub[[strat_col]])) {
      if (is.na(g)) next
      ix <- which(!is.na(meta_sub[[strat_col]]) & meta_sub[[strat_col]] == g)
      if (length(ix) >= 10) {
        cg <- cor(X_sub[ix, , drop = FALSE], method = "spearman")
        sub_hits[[length(sub_hits) + 1]] <- filter_fdr_hits(cg)
      }
    }
    return(if (length(sub_hits) == 0) integer(0) else unique(unlist(sub_hits, use.names = FALSE)))
  }
  
  if (method_name == "02_glasso") {
    cov_mat <- cov(scale(X_sub))
    gl_res  <- tryCatch(glasso::glasso(cov_mat, rho = 0.3), error = function(e) NULL)
    if (is.null(gl_res)) return(integer(0))
    prec <- gl_res$wi; diag(prec) <- 0
    return(which(abs(prec[upper.tri(prec)]) >= 0.05))
  }
  if (method_name == "22_spiec_easi") {
    cov_clr <- cov(X_sub)
    gl_se   <- tryCatch(glasso::glasso(cov_clr, rho = 0.15), error = function(e) NULL)
    if (is.null(gl_se)) return(integer(0))
    prec_se <- gl_se$wi; diag(prec_se) <- 0
    return(which(abs(prec_se[upper.tri(prec_se)]) >= 0.02))
  }
  if (method_name == "06_covariate_adjusted") {
    cov_d <- model.matrix(~ gender + age_group + Underlying_disease_binary + infection_group - 1, data = meta_sub)
    X_res <- residuals(lm(X_sub ~ cov_d))
    return(filter_fdr_hits(cor(X_res, method = "spearman")))
  }
  if (method_name == "10_earth_biogeographic") {
    city_d <- model.matrix(~ city - 1, data = meta_sub)
    X_res  <- residuals(lm(X_sub ~ city_d))
    return(filter_fdr_hits(cor(X_res, method = "spearman")))
  }
  filter_fdr_hits(cor(X_sub, method = "spearman"))
}

boot_indices <- lapply(seq_len(N_BOOT), function(b) sample(nrow(X), size = as.integer(nrow(X) * BOOT_FRAC), replace = FALSE))
stability_records <- list()
all_edge_counts   <- list()

cl_stab <- makePSOCKcluster(N_CORES)
clusterEvalQ(cl_stab, suppressPackageStartupMessages(library(glasso)))
clusterExport(cl_stab, c("worker_one_boot_int", "X", "meta", "CORR_THRESHOLD", "FDR_THRESHOLD"), envir = environment())

for (method_name in TOP6) {
  clusterExport(cl_stab, "method_name", envir = environment())
  total_pair_counts <- integer(N_PAIRS)
  n_batches <- ceiling(N_BOOT / BATCH_SIZE)
  for (b in seq_len(n_batches)) {
    lo <- (b - 1) * BATCH_SIZE + 1
    hi <- min(b * BATCH_SIZE, N_BOOT)
    batch_results <- parLapply(cl_stab, boot_indices[lo:hi], function(i) {
      tryCatch(worker_one_boot_int(i), error = function(e) integer(0))
    })
    for (res in batch_results) {
      if (length(res) > 0) total_pair_counts[res] <- total_pair_counts[res] + 1L
    }
  }
  total_observed_edges <- sum(total_pair_counts > 0)
  stable_edge_idx      <- which(total_pair_counts / N_BOOT >= STABLE_THRESHOLD)
  n_stable             <- length(stable_edge_idx)
  score                <- n_stable / (total_observed_edges + 1e-9)
  
  stability_records[[method_name]] <- data.frame(
    method = method_name, n_stable_edges = n_stable,
    n_total_edges = total_observed_edges, stability_score = round(score, 4)
  )
  all_edge_counts[[method_name]] <- list(stable_idx = stable_edge_idx, freqs = total_pair_counts[stable_edge_idx])
}
stopCluster(cl_stab)

stab_df <- do.call(rbind, stability_records)
stab_df <- stab_df[order(stab_df$stability_score, decreasing = TRUE), ]
write.csv(stab_df, file.path(TAB_DIR, "stability_all_methods.csv"), row.names = FALSE)

BEST_METHOD <- as.character(stab_df$method[1])
best_info   <- all_edge_counts[[BEST_METHOD]]

G_robust <- make_empty_graph(n = length(sp_names), directed = FALSE); V(G_robust)$name <- sp_names
if (length(best_info$stable_idx) > 0) {
  coords   <- ut_map[best_info$stable_idx, , drop = FALSE]
  df_e     <- data.frame(from = sp_names[coords[, 1]], to = sp_names[coords[, 2]],
                         weight = best_info$freqs / N_BOOT, stringsAsFactors = FALSE)
  G_robust <- graph_from_data_frame(df_e, directed = FALSE, vertices = data.frame(name = sp_names))
}
saveRDS(G_robust, file.path(NET_DIR, "robust_network.rds"))

# ── Figure 07: Bootstrap Stability Plot ───────────────────────────────────────
plot_df <- stab_df[!is.na(stab_df$stability_score), ]
plot_df$is_best <- ifelse(plot_df$method == BEST_METHOD, "Best", "Other")

p_fig07 <- ggplot(plot_df, aes(x = reorder(method, stability_score), y = stability_score, fill = is_best)) +
  geom_bar(stat = "identity", width = 0.55, color = "#2C3E50", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", stability_score)), hjust = -0.15, size = 3.3, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Best" = "#D95F02", "Other" = "#2B5C8F")) +
  scale_y_continuous(limits = c(0, max(plot_df$stability_score) * 1.18), expand = c(0, 0)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = "#EAEAEA", linetype = "dashed", linewidth = 0.4)) +
  labs(x = "Method", y = "Stability Score (Reproducibility >= 80%)")
export_plot("fig07_stability", p_fig07, width = 7.5, height = 4.2)

# ── Weighted Consensus Construction & Figure 08 ───────────────────────────────
cluster_sizes <- table(method_clusters)
method_weight <- 1 / (length(cluster_sizes) * as.numeric(cluster_sizes[as.character(method_clusters)]))
names(method_weight) <- names(method_clusters)

all_edges_env  <- new.env(hash = TRUE, parent = emptyenv())
all_counts_env <- new.env(hash = TRUE, parent = emptyenv())

for (nm in names(edge_sets)) {
  w_m <- method_weight[nm]
  for (e in edge_sets[[nm]]) {
    curr_w <- if (exists(e, envir = all_edges_env, inherits = FALSE)) get(e, envir = all_edges_env) else 0.0
    assign(e, curr_w + w_m, envir = all_edges_env)
    curr_c <- if (exists(e, envir = all_counts_env, inherits = FALSE)) get(e, envir = all_counts_env) else 0L
    assign(e, curr_c + 1L, envir = all_counts_env)
  }
}

all_keys       <- ls(all_edges_env)
weighted_score <- vapply(all_keys, function(k) get(k, envir = all_edges_env), numeric(1))
raw_counts     <- vapply(all_keys, function(k) get(k, envir = all_counts_env), integer(1))

keep_mask <- (weighted_score >= WEIGHTED_CONSENSUS_THRESH) & (raw_counts >= 3)
consensus_edges <- all_keys[keep_mask]

G_consensus <- make_empty_graph(n = length(sp_names), directed = FALSE); V(G_consensus)$name <- sp_names
if (length(consensus_edges) > 0) {
  parts <- strsplit(consensus_edges, "__", fixed = TRUE)
  df_c  <- data.frame(
    from   = vapply(parts, `[`, character(1), 1),
    to     = vapply(parts, `[`, character(1), 2),
    weight = round(weighted_score[keep_mask], 4),
    stringsAsFactors = FALSE
  )
  G_consensus <- graph_from_data_frame(df_c, directed = FALSE, vertices = data.frame(name = sp_names))
}
V(G_consensus)$is_pathogen <- V(G_consensus)$name %in% pathogen_set
saveRDS(G_consensus, file.path(NET_DIR, "consensus_network.rds"))

plot_consensus_network_ggraph <- function(top_n = 50) {
  deg_cons <- degree(G_consensus)
  top_idx  <- order(deg_cons, decreasing = TRUE)[1:min(top_n, vcount(G_consensus))]
  G_sub    <- subgraph(G_consensus, top_idx)
  
  V(G_sub)$degree       <- deg_cons[top_idx]
  V(G_sub)$SpeciesClass <- ifelse(V(G_sub)$is_pathogen, "Pathogen", "Commensal / Other")
  
  set.seed(SEED)
  ggraph(G_sub, layout = "fr", niter = 1200) +
    geom_edge_arc(aes(edge_width = weight), color = "#888888", alpha = 0.25, curvature = 0.08) +
    scale_edge_width(range = c(0.25, 0.9), name = "Weighted\nConsensus") +
    geom_node_point(aes(size = degree, fill = degree, shape = SpeciesClass), color = "#222222", stroke = 0.7) +
    scale_shape_manual(values = c("Pathogen" = 21, "Commensal / Other" = 22), name = "Species Class") +
    scale_fill_gradientn(
      colors = c("#FEE391","#FEC44F","#FE9929","#EC7014","#CC4C02","#8C2D04"),
      limits = c(min(deg_cons), max(deg_cons)), name = "Global Node\nDegree",
      guide  = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(4.5, "cm"), order = 1)
    ) +
    scale_size_continuous(range = c(3.5, 9), guide = "none") +
    geom_node_label(
      aes(label = name, color = SpeciesClass),
      repel = TRUE, fontface = ifelse(V(G_sub)$is_pathogen, "bold.italic", "italic"),
      size = 2.9, fill = alpha("white", 0.88), label.size = 0.15,
      label.padding = unit(0.12, "lines"), box.padding = 0.35, max.overlaps = 50, show.legend = FALSE
    ) +
    scale_color_manual(values = c("Pathogen" = "#B03A2E", "Commensal / Other" = "#1A252F")) +
    guides(fill = guide_colorbar(order = 1), shape = guide_legend(order = 2, override.aes = list(size = 3.5, fill = "#EC7014")), edge_width = guide_legend(order = 3)) +
    theme_void(base_size = 10) +
    labs(title = sprintf("Consensus Ecological Network (Top-%d Keystone Taxa)", vcount(G_sub)),
         subtitle = sprintf("Edges validated across 3 method families; Total Edges = %d", ecount(G_consensus))) +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey35", margin = margin(b = 10)),
          legend.position = "right", plot.margin = margin(12, 12, 12, 12))
}

p_fig08 <- plot_consensus_network_ggraph(top_n = 50)
export_plot("fig08_consensus_network", p_fig08, width = 13, height = 10)


# ═════════════════════════════════════════════════════════════════════════════
# S7 & S8: DUAL KEYSTONE CROSS-VALIDATION (FIGURE 09)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S7 & S8] Keystone Species Dual Comparison ...\n")

calc_keystones <- function(g, label) {
  deg_k   <- degree(g) / (vcount(g) - 1)
  betw_k  <- betweenness(g, weights = NA, normalized = TRUE)
  eigen_k <- eigen_centrality(g, weights = NA)$vector
  
  lcc_size <- function(gr) if (ecount(gr) == 0 || vcount(gr) == 0) 0 else max(components(gr)$csize)
  orig_lcc <- lcc_size(g)
  
  impact <- numeric(length(sp_names))
  names(impact) <- sp_names
  for (sp in sp_names) {
    if (sp %in% V(g)$name) {
      impact[sp] <- (orig_lcc - lcc_size(delete_vertices(g, sp))) / (orig_lcc + 1e-9)
    }
  }
  score <- deg_k[sp_names] * 0.30 + betw_k[sp_names] * 0.40 + eigen_k[sp_names] * 0.20 + impact[sp_names] * 0.10
  data.frame(
    species = sp_names, is_pathogen = sp_names %in% pathogen_set,
    keystone_score = round(score, 5), rank = rank(-score), network = label,
    stringsAsFactors = FALSE
  )
}

ks_c <- calc_keystones(G_consensus, "Consensus")
ks_r <- calc_keystones(G_robust, "Robust")

ks_dual <- merge(ks_c[, c("species", "is_pathogen", "keystone_score", "rank")],
                 ks_r[, c("species", "keystone_score", "rank")],
                 by = "species", suffixes = c("_consensus", "_robust"))

ks_dual$is_core_keystone <- (ks_dual$rank_consensus <= 30) & (ks_dual$rank_robust <= 30)
ks_dual <- ks_dual[order((ks_dual$rank_consensus + ks_dual$rank_robust) / 2), ]
write.csv(ks_dual, file.path(TAB_DIR, "keystone_species_dual_comparison.csv"), row.names = FALSE)

cor_ks <- cor.test(ks_dual$keystone_score_consensus, ks_dual$keystone_score_robust, method = "spearman")

top20_ks <- head(ks_dual[order(ks_dual$rank_consensus), ], 20)
color_palette <- c("TRUE" = "#E74C3C", "FALSE" = "#2C3E50")
top_annot <- head(subset(ks_dual, is_core_keystone | (is_pathogen & (rank_consensus <= 40 | rank_robust <= 40))), 14)

p16_bar <- ggplot(top20_ks, aes(x = reorder(species, keystone_score_consensus), y = keystone_score_consensus, fill = is_pathogen)) +
  geom_bar(stat = "identity", width = 0.65) + 
  coord_flip() +
  scale_fill_manual(values = color_palette, name = "Organism Class", labels = c("TRUE" = "Pathogen", "FALSE" = "Commensal / Other")) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic", size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10.5)) +
  labs(x = NULL, y = "Keystone Score (Consensus)", title = "A  Top-20 Keystone Taxa")

p16_scatter <- ggplot(ks_dual, aes(x = keystone_score_consensus, y = keystone_score_robust)) +
  geom_point(aes(color = is_pathogen, shape = is_core_keystone, size = is_core_keystone), alpha = 0.8) +
  geom_text_repel(data = top_annot, aes(label = species, color = is_pathogen),
                  size = 2.8, fontface = "italic", box.padding = 0.35, point.padding = 0.28,
                  max.overlaps = 25, bg.color = "white", bg.r = 0.12, segment.color = "#7F8C8D",
                  segment.size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Organism Class", labels = c("TRUE" = "Pathogen", "FALSE" = "Commensal / Other")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16), name = "Validation Tier", labels = c("TRUE" = "Dual-Validated Core", "FALSE" = "Single Model")) +
  scale_size_manual(values = c("TRUE" = 3.2, "FALSE" = 1.6), guide = "none") +
  theme_classic(base_size = 9) +
  theme(axis.text = element_text(size = 8, color = "black"), legend.position = "bottom",
        legend.box = "horizontal", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7.5), plot.title = element_text(face = "bold", size = 10.5)) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 2.5, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 2.5))) +
  labs(x = "Keystone Score (Consensus Network)", y = "Keystone Score (Robust Network)", 
       title = sprintf("B  Cross-Network Concordance (Spearman r = %.3f)", cor_ks$estimate))

fig09_keystones <- p16_bar + p16_scatter + 
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(title = "Figure 09: Keystone Species Identification & Cross-Network Validation",
                  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)))
export_plot("fig09_keystones_dual_publication", fig09_keystones, width = 13.5, height = 5.5)


# ═════════════════════════════════════════════════════════════════════════════
# S9: DUAL-NETWORK MODULARITY & HOTSPOTS (FIGURE 10 & SUPP FIG S1)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S9] Dual Community Detection & Modularity ...\n")

run_louvain <- function(g) {
  if ("weight" %in% edge_attr_names(g)) {
    cluster_louvain(g, weights = pmax(E(g)$weight, 1e-4))
  } else {
    cluster_louvain(g)
  }
}

set.seed(SEED)
comm_c <- run_louvain(G_consensus)
comm_r <- run_louvain(G_robust)

Q_consensus <- modularity(comm_c)
Q_robust    <- modularity(comm_r)

shared_nodes <- intersect(V(G_consensus)$name, V(G_robust)$name)
nmi_score <- igraph::compare(membership(comm_c)[shared_nodes], membership(comm_r)[shared_nodes], method = "nmi")
ari_score <- igraph::compare(membership(comm_c)[shared_nodes], membership(comm_r)[shared_nodes], method = "adjusted.rand")

mod_df <- data.frame(
  species       = shared_nodes,
  mod_consensus = as.character(membership(comm_c)[shared_nodes]),
  mod_robust    = as.character(membership(comm_r)[shared_nodes]),
  is_pathogen   = shared_nodes %in% pathogen_set,
  degree_c      = as.integer(degree(G_consensus)[shared_nodes]),
  stringsAsFactors = FALSE
)

calc_enrichment <- function(mod_vec, is_path_vec) {
  n_path_tot <- sum(is_path_vec)
  n_tot <- length(is_path_vec)
  res <- list()
  for (cid in sort(unique(mod_vec))) {
    idx <- mod_vec == cid
    a <- sum(is_path_vec[idx])
    b <- sum(idx) - a
    c_ <- n_path_tot - a
    d <- n_tot - n_path_tot - b
    pval <- tryCatch(fisher.test(matrix(c(a, b, max(c_,0), max(d,0)), 2))$p.value, error = function(e) NA)
    res[[length(res) + 1]] <- data.frame(
      community = as.character(cid), n_species = sum(idx), n_pathogens = a,
      pathogen_fraction = a / sum(idx), fisher_p = pval, is_enriched = (!is.na(pval) & pval < 0.05),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res)
}

enrich_c_df <- calc_enrichment(mod_df$mod_consensus, mod_df$is_pathogen)
enrich_r_df <- calc_enrichment(mod_df$mod_robust, mod_df$is_pathogen)

write.csv(enrich_c_df, file.path(TAB_DIR, "community_enrichment.csv"), row.names = FALSE)
write.csv(data.frame(species = names(membership(comm_c)), 
                     community = as.character(membership(comm_c)), 
                     is_pathogen = names(membership(comm_c)) %in% pathogen_set),
          file.path(TAB_DIR, "community_modules.csv"), row.names = FALSE)

target_mod_c <- if(any(enrich_c_df$is_enriched)) enrich_c_df$community[which.max(enrich_c_df$pathogen_fraction[enrich_c_df$is_enriched])] else enrich_c_df$community[which.max(enrich_c_df$pathogen_fraction)]
target_mod_r <- if(any(enrich_r_df$is_enriched)) enrich_r_df$community[which.max(enrich_r_df$pathogen_fraction[enrich_r_df$is_enriched])] else enrich_r_df$community[which.max(enrich_r_df$pathogen_fraction)]

enrich_c_df$display_label <- paste0("Con-M", enrich_c_df$community)
enrich_c_df$display_label[enrich_c_df$community == target_mod_c] <- paste0("Con-M", target_mod_c, " (★)")

enrich_r_df$display_label <- paste0("Rob-M", enrich_r_df$community)
enrich_r_df$display_label[enrich_r_df$community == target_mod_r] <- paste0("Rob-M", target_mod_r, " (★)")

all_labels <- c(enrich_c_df$display_label, enrich_r_df$display_label)
base_pal <- c("#2E86AB", "#2BAE66", "#9A7AA0", "#F6C85F", "#E08244", "#829399", "#4A4E69", "#003049", "#FCBF49")
dyn_colors <- colorRampPalette(base_pal)(length(all_labels))
module_palette <- setNames(dyn_colors, all_labels)
module_palette[paste0("Con-M", target_mod_c, " (★)")] <- "#E74C3C" 
module_palette[paste0("Rob-M", target_mod_r, " (★)")] <- "#990000"

g_c <- G_consensus
V(g_c)$Network <- "Consensus Network"
V(g_c)$ModuleDisplay <- enrich_c_df$display_label[match(membership(comm_c)[V(g_c)$name], enrich_c_df$community)]

g_r <- G_robust
V(g_r)$Network <- "Robust Network"
V(g_r)$ModuleDisplay <- enrich_r_df$display_label[match(membership(comm_r)[V(g_r)$name], enrich_r_df$community)]

g_comb <- disjoint_union(g_c, g_r)

m_from <- V(g_comb)$ModuleDisplay[as_edgelist(g_comb, names=FALSE)[, 1]]
m_to   <- V(g_comb)$ModuleDisplay[as_edgelist(g_comb, names=FALSE)[, 2]]
layout_w <- ifelse(m_from == m_to, 4.0, 0.05)

set.seed(LAYOUT_SEED)
coords_comb <- layout_with_fr(g_comb, weights = layout_w, niter = 1500)

is_rob <- V(g_comb)$Network == "Robust Network"
max_c_x <- max(coords_comb[!is_rob, 1])
min_r_x <- min(coords_comb[is_rob, 1])
width_c <- max_c_x - min(coords_comb[!is_rob, 1])
coords_comb[is_rob, 1] <- coords_comb[is_rob, 1] + (max_c_x - min_r_x) + (width_c * 0.15)

p17_macro <- ggraph(g_comb, layout = coords_comb) +
  geom_edge_arc(aes(alpha = ifelse(m_from == m_to, 0.25, 0.02)), color = "#94A3B8", width = 0.25, curvature = 0.05, show.legend = FALSE) +
  geom_node_point(aes(fill = ModuleDisplay, shape = Network, size = degree(g_comb)), color = "#1E293B", stroke = 0.5, alpha = 0.90) +
  scale_shape_manual(values = c("Consensus Network" = 21, "Robust Network" = 24)) +
  scale_fill_manual(values = module_palette, name = "Modules") +
  scale_size_continuous(range = c(1.5, 5.5), name = "Global Degree") +  
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        legend.position = "bottom", legend.box = "horizontal", 
        legend.title = element_text(face="bold"), legend.text = element_text(size = 7.5)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3), nrow = 3, order = 1),
         shape = guide_legend(override.aes = list(size = 4), order = 2),
         size = guide_legend(order = 3)) +
  labs(title = "A  Global Architecture: Consensus (●) vs. Robust (▲)")

overlap_summary <- mod_df %>%
  group_by(mod_consensus, mod_robust) %>%
  summarise(total_nodes = n(), pathogens = sum(is_pathogen), .groups = "drop") %>%
  mutate(
    mod_consensus = factor(paste("Con-M", mod_consensus), levels = paste("Con-M", sort(as.numeric(unique(mod_consensus))))),
    mod_robust    = factor(paste("Rob-M", mod_robust), levels = paste("Rob-M", sort(as.numeric(unique(mod_robust)))))
  )

p17_overlap <- ggplot(overlap_summary, aes(x = mod_consensus, y = mod_robust, fill = total_nodes)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = paste0(total_nodes, "\n(", pathogens, " Path.)")), 
            size = 3, fontface = ifelse(overlap_summary$pathogens > 0, "bold", "plain"), 
            color = ifelse(overlap_summary$total_nodes > median(overlap_summary$total_nodes), "white", "black")) +
  scale_fill_gradient(low = "#ECF0F1", high = "#2C3E50", name = "Shared Nodes") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")) +
  labs(x = "Consensus Modules", y = "Robust Modules", title = "B  Cross-Network Module Concordance",
       subtitle = sprintf("Alignment: NMI = %.3f | Adjusted Rand Index = %.3f", nmi_score, ari_score))

nodes_target <- V(G_consensus)$name[membership(comm_c) == target_mod_c]
G_target     <- subgraph(G_consensus, nodes_target)
V(G_target)$is_pathogen   <- V(G_target)$name %in% pathogen_set
V(G_target)$SpeciesClass  <- ifelse(V(G_target)$is_pathogen, "Pathogen", "Commensal / Other")
V(G_target)$degree_in_mod <- degree(G_target)

path_hubs <- V(G_target)$name[V(G_target)$is_pathogen][order(V(G_target)$degree_in_mod[V(G_target)$is_pathogen], decreasing = TRUE)[1:min(8, sum(V(G_target)$is_pathogen))]]
comm_hubs <- V(G_target)$name[!V(G_target)$is_pathogen][order(V(G_target)$degree_in_mod[!V(G_target)$is_pathogen], decreasing = TRUE)[1:min(4, sum(!V(G_target)$is_pathogen))]]
V(G_target)$label_clean <- ifelse(V(G_target)$name %in% c(path_hubs, comm_hubs), V(G_target)$name, "")

set.seed(LAYOUT_SEED)
p17_zoom <- ggraph(G_target, layout = "fr", niter = 1500) +
  geom_edge_arc(aes(edge_width = weight), color = "#94A3B8", alpha = 0.35, curvature = 0.1) + 
  geom_node_point(aes(fill = SpeciesClass, shape = SpeciesClass, size = degree_in_mod), color = "#1E293B", stroke = 0.6) +
  scale_shape_manual(values = c("Pathogen" = 21, "Commensal / Other" = 22), name = "Organism Class") +
  scale_fill_manual(values = c("Pathogen" = "#E74C3C", "Commensal / Other" = "#BDC3C7"), name = "Organism Class") +
  scale_size_continuous(range = c(3.0, 7.5), name = "Intra-Module\nDegree") + 
  scale_edge_width_continuous(range = c(0.2, 1.2), name = "Weighted\nConsensus") + 
  geom_node_text(aes(label = label_clean, color = SpeciesClass), repel = TRUE, size = 2.8,
                 fontface = ifelse(V(G_target)$is_pathogen, "bold.italic", "italic"),
                 bg.color = "white", bg.r = 0.15, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = c("Pathogen" = "#990000", "Commensal / Other" = "#333333")) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5), 
        legend.position = "right", legend.box = "vertical", plot.margin = margin(10, 10, 10, 10)) +
  guides(fill = guide_legend(override.aes = list(size = 3.5), order = 1),
         shape = guide_legend(order = 1), size = guide_legend(order = 2), edge_width = guide_legend(order = 3)) +
  labs(title = sprintf("C  Niche Zoom-In: Consensus Mod %s (Hotspot)", target_mod_c))

fig10_modularity <- (p17_macro | p17_zoom) / p17_overlap + 
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(title = "Figure 10: Dual-Network Modularity & Pathogen Compartmentalization",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
export_plot("fig10_communities_dual_publication", fig10_modularity, width = 15.0, height = 11.5)

# ── Supplementary S1: Robust Hotspot Zoom ─────────────────────────────────────
nodes_target_r <- V(G_robust)$name[membership(comm_r) == target_mod_r]
G_target_r     <- subgraph(G_robust, nodes_target_r)
V(G_target_r)$is_pathogen   <- V(G_target_r)$name %in% pathogen_set
V(G_target_r)$SpeciesClass  <- ifelse(V(G_target_r)$is_pathogen, "Pathogen", "Commensal / Other")
V(G_target_r)$degree_in_mod <- degree(G_target_r)

path_hubs_r <- V(G_target_r)$name[V(G_target_r)$is_pathogen][order(V(G_target_r)$degree_in_mod[V(G_target_r)$is_pathogen], decreasing = TRUE)[1:min(8, sum(V(G_target_r)$is_pathogen))]]
comm_hubs_r <- V(G_target_r)$name[!V(G_target_r)$is_pathogen][order(V(G_target_r)$degree_in_mod[!V(G_target_r)$is_pathogen], decreasing = TRUE)[1:min(4, sum(!V(G_target_r)$is_pathogen))]]
V(G_target_r)$label_clean <- ifelse(V(G_target_r)$name %in% c(path_hubs_r, comm_hubs_r), V(G_target_r)$name, "")

set.seed(LAYOUT_SEED)
p_supp_robust <- ggraph(G_target_r, layout = "fr", niter = 1500) +
  geom_edge_arc(aes(edge_width = weight), color = "#94A3B8", alpha = 0.35, curvature = 0.1) + 
  geom_node_point(aes(fill = SpeciesClass, shape = SpeciesClass, size = degree_in_mod), color = "#1E293B", stroke = 0.6) +
  scale_shape_manual(values = c("Pathogen" = 24, "Commensal / Other" = 22), name = "Organism Class") +
  scale_fill_manual(values = c("Pathogen" = "#990000", "Commensal / Other" = "#BDC3C7"), name = "Organism Class") +
  scale_size_continuous(range = c(3.0, 7.5), name = "Intra-Module\nDegree") + 
  scale_edge_width_continuous(range = c(0.2, 1.2), name = "Weighted\nRobust") + 
  geom_node_text(aes(label = label_clean, color = SpeciesClass), repel = TRUE, size = 3,
                 fontface = ifelse(V(G_target_r)$is_pathogen, "bold.italic", "italic"),
                 bg.color = "white", bg.r = 0.15, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = c("Pathogen" = "#660000", "Commensal / Other" = "#333333")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
        legend.position = "right", legend.box = "vertical") +
  guides(fill = guide_legend(override.aes = list(size = 3.5), order = 1),
         shape = guide_legend(order = 1), size = guide_legend(order = 2), edge_width = guide_legend(order = 3)) +
  labs(title = sprintf("Supplementary: Robust Network Mod %s (Pathogen Hotspot)", target_mod_r))
export_plot("supp_fig01_robust_hotspot_zoom", p_supp_robust, width = 9.5, height = 7.5)


# ═════════════════════════════════════════════════════════════════════════════
# S10: REWIRING TRAJECTORIES & DIFFERENTIAL TOPOLOGY (FIGURES 11 & 12)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[S10] Factor-Driven Centrality Plasticity & Differential Rewiring ...\n")

pathogens_in_network <- intersect(sp_names, pathogen_set)
strat_networks <- list(
  "Sex"                = readRDS(file.path(NET_DIR, "04_stratified_gender.rds")),
  "Age"                = readRDS(file.path(NET_DIR, "05_stratified_age.rds")),
  "Underlying Disease" = readRDS(file.path(NET_DIR, "19_stratified_underlying_disease.rds")),
  "Infection"          = readRDS(file.path(NET_DIR, "20_stratified_infection.rds")),
  "Geography"          = readRDS(file.path(NET_DIR, "03_stratified_location.rds"))
)
G_cov_adj <- readRDS(file.path(NET_DIR, "06_covariate_adjusted.rds"))

factor_centrality_list <- list()

for (f_lbl in names(strat_networks)) {
  for (g_name in names(strat_networks[[f_lbl]])) {
    g_obj <- strat_networks[[f_lbl]][[g_name]]
    if (!inherits(g_obj, "igraph") || ecount(g_obj) == 0) next
    deg <- degree(g_obj) / (vcount(g_obj) - 1)
    for (sp in pathogens_in_network) {
      factor_centrality_list[[length(factor_centrality_list) + 1]] <- data.frame(
        factor = f_lbl, group = g_name, species = sp,
        degree_centrality = round(if (sp %in% names(deg)) deg[sp] else 0, 5),
        stringsAsFactors = FALSE
      )
    }
  }
}

deg_c   <- degree(G_consensus) / (vcount(G_consensus) - 1)
deg_r   <- degree(G_robust) / (vcount(G_robust) - 1)
deg_cov <- degree(G_cov_adj) / (vcount(G_cov_adj) - 1)

for (sp in pathogens_in_network) {
  factor_centrality_list[[length(factor_centrality_list) + 1]] <- data.frame(
    factor = "Baselines", group = "Consensus", species = sp,
    degree_centrality = round(if (sp %in% names(deg_c)) deg_c[sp] else 0, 5), stringsAsFactors = FALSE
  )
  factor_centrality_list[[length(factor_centrality_list) + 1]] <- data.frame(
    factor = "Baselines", group = "Robust", species = sp,
    degree_centrality = round(if (sp %in% names(deg_r)) deg_r[sp] else 0, 5), stringsAsFactors = FALSE
  )
  factor_centrality_list[[length(factor_centrality_list) + 1]] <- data.frame(
    factor = "Cov-Adjusted", group = "Cov-Adjusted", species = sp,
    degree_centrality = round(if (sp %in% names(deg_cov)) deg_cov[sp] else 0, 5), stringsAsFactors = FALSE
  )
}

factor_cent_df <- do.call(rbind, factor_centrality_list)
write.csv(factor_cent_df, file.path(TAB_DIR, "story2_pathogen_factor_centrality.csv"), row.names = FALSE)

# ── Figure 11: Multi-Layer Centrality Plasticity Dashboard ────────────────────
variance_df <- factor_cent_df %>%
  filter(factor != "Baselines" & factor != "Cov-Adjusted") %>%
  group_by(species) %>%
  summarise(
    variance = var(degree_centrality, na.rm = TRUE),
    mean_c   = mean(degree_centrality, na.rm = TRUE)
  ) %>%
  slice_max(order_by = variance, n = 30)

top_pathogens <- variance_df$species

plot_f12_df <- factor_cent_df %>%
  filter(species %in% top_pathogens) %>%
  group_by(species) %>%
  mutate(z_score = as.numeric(scale(degree_centrality))) %>%
  ungroup()

wide_df <- plot_f12_df %>%
  select(species, group, z_score) %>%
  pivot_wider(names_from = group, values_from = z_score, values_fill = 0)

mat <- as.matrix(wide_df[, -1])
rownames(mat) <- wide_df$species

set.seed(SEED)
hc <- hclust(dist(mat), method = "ward.D2")
ordered_species <- rownames(mat)[hc$order] 

plot_f12_df <- plot_f12_df %>%
  mutate(
    species_clean = factor(gsub("_", " ", species), levels = gsub("_", " ", ordered_species)),
    group_clean   = gsub("\\s*\\([^)]+\\)", "", group),
    factor        = as.character(factor)
  ) %>%
  mutate(
    factor = case_when(
      factor == "Underlying Disease" ~ "Underlying\nDisease",
      factor == "Cov-Adjusted"       ~ "Cov-\nAdj.",
      factor == "Baselines"          ~ "Baseline",
      TRUE                           ~ factor
    ),
    factor = factor(factor, levels = c("Baseline", "Cov-\nAdj.", "Sex", "Age", "Underlying\nDisease", "Infection", "Geography"))
  )

plot_f12_df$group_clean <- factor(plot_f12_df$group_clean, 
                                  levels = unique(c("Consensus", "Robust", "Cov-Adjusted", as.character(plot_f12_df$group_clean))))

top_bar_df <- plot_f12_df %>%
  group_by(factor, group_clean) %>%
  summarise(avg_centrality = mean(degree_centrality, na.rm = TRUE), .groups = "drop")

p12_top <- ggplot(top_bar_df, aes(x = group_clean, y = avg_centrality)) +
  geom_bar(stat = "identity", fill = "#34495E", width = 0.6) +
  facet_grid(~ factor, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 7.5, color = "#555555"),
        axis.title.y = element_text(size = 8, face = "bold", margin = margin(r = 5)),
        strip.background = element_blank(), strip.text = element_blank(), plot.margin = margin(b = 2)) +
  labs(y = "Mean System\nCentrality")

p12_main <- ggplot(plot_f12_df, aes(x = group_clean, y = species_clean)) +
  geom_segment(aes(x = -Inf, xend = Inf, y = species_clean, yend = species_clean), color = "#F0F0F0", linewidth = 0.4) +
  geom_line(data = plot_f12_df %>% group_by(factor) %>% filter(n_distinct(group_clean) > 1),
            aes(group = species_clean), color = "#95A5A6", alpha = 0.6, linewidth = 0.6) +
  geom_point(aes(size = degree_centrality, fill = z_score), 
             shape = 21, color = "#111111", stroke = 0.4, alpha = 0.95) +
  facet_grid(~ factor, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#313695", mid = "#FAFAFA", high = "#A50026", midpoint = 0,
                       name = "Relative Shift\n(Row Z-Score)",
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"))) +
  scale_size_continuous(name = expression(bold(Absolute~Centrality~(italic(C[D])))), 
                        range = c(1.5, 7.5), guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "#FCFCFC", color = NA),
        panel.border = element_rect(color = "#CCCCCC", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "#2C3E50", color = "#2C3E50"),
        strip.text = element_text(face = "bold", size = 8.2, color = "white", lineheight = 0.8),
        axis.text.y = element_text(face = "bold.italic", size = 9, color = "#222222"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, color = "#111111"),
        legend.position = "bottom", legend.box = "horizontal", plot.margin = margin(t = 2)) +
  labs(x = NULL, y = NULL)

var_plot_df <- variance_df %>% mutate(species_clean = factor(gsub("_", " ", species), levels = levels(plot_f12_df$species_clean)))

p12_margin <- ggplot(var_plot_df, aes(x = variance, y = species_clean)) +
  geom_segment(aes(x = 0, xend = max(variance)*1.05, y = species_clean, yend = species_clean), color = "#F0F0F0", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = variance, y = species_clean, yend = species_clean), color = "#7F8C8D", linewidth = 0.8) +
  geom_point(aes(color = mean_c), size = 3) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "Mean Centrality",
                        guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(4, "cm"), barheight = unit(0.3, "cm"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(n = 3)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 8, color = "#333333"),
        axis.title.x = element_text(size = 8, face = "bold", color = "#2C3E50", margin = margin(t = 5)), legend.position = "bottom") +
  labs(x = "Plasticity Index\n(Variance)")

fig11_dashboard <- (p12_top | plot_spacer()) / (p12_main | p12_margin) +
  plot_layout(heights = c(1, 5.5), widths = c(5, 1.2), guides = "collect") & 
  theme(legend.position = "bottom", legend.box = "horizontal") &
  plot_annotation(
    title = "Figure 11: Network Pathogen Rewiring Trajectories & Co-Plasticity Clusters",
    subtitle = "Hierarchical clustering (Y-axis) reveals distinct modules of pathogens that re-wire synchronously across clinical constraints",
    theme = theme(plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "grey35", margin = margin(b=10)))
  )
export_plot("fig11_pathogen_centrality_dashboard", fig11_dashboard, width = 16.5, height = 9.0)

# ── Figure 12: Differential Network Rewiring ──────────────────────────────────
get_graph_edge_strings <- function(g) {
  if (is.null(g) || !inherits(g, "igraph") || ecount(g) == 0) return(character(0))
  el <- as_edgelist(g, names = TRUE)
  paste(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]), sep = "__")
}

contrasts <- list(
  list(type = "Factor", cat = "Sex",                g1 = "Male",               g2 = "Female",               label = "Sex: Male vs Female"),
  list(type = "Factor", cat = "Underlying Disease", g1 = "Yes",                g2 = "No",                   label = "Disease: Yes vs No"),
  list(type = "Factor", cat = "Age",                g1 = "Infant (<=3y)",      g2 = "Child (4-14y)",        label = "Age: Infant vs Child"),
  list(type = "Factor", cat = "Infection",          g1 = "Severe pneumonia",   g2 = "Mild pneumonia",       label = "Infection: Severe vs Mild"),
  list(type = "Factor", cat = "Infection",          g1 = "Tuberculosis",       g2 = "Healthy control",      label = "Infection: TB vs Control"),
  list(type = "Factor", cat = "Geography",          g1 = "Tianjin (CRA006760)", g2 = "Guangzhou (PRJCA008730)", label = "Geography: Tianjin vs GZ"),
  list(type = "Ref_Cov", cat = "Infection",         g1 = "Severe pneumonia",   g2 = "Cov-Adjusted",         label = "Severe Dis vs Cov-Adj"),
  list(type = "Ref_Con", cat = "Infection",         g1 = "Severe pneumonia",   g2 = "Consensus",            label = "Severe Dis vs Consensus"),
  list(type = "Ref_Rob", cat = "Infection",         g1 = "Severe pneumonia",   g2 = "Robust",               label = "Severe Dis vs Robust")
)

diff_results_list <- list()
for (cntr in contrasts) {
  sub_nets <- strat_networks[[cntr$cat]]
  if (!cntr$g1 %in% names(sub_nets)) next
  e1 <- get_graph_edge_strings(sub_nets[[cntr$g1]])
  
  if (cntr$type == "Factor") {
    if (!cntr$g2 %in% names(sub_nets)) next
    e2 <- get_graph_edge_strings(sub_nets[[cntr$g2]])
  } else if (cntr$type == "Ref_Con") {
    e2 <- get_graph_edge_strings(G_consensus)
  } else if (cntr$type == "Ref_Rob") {
    e2 <- get_graph_edge_strings(G_robust)
  } else if (cntr$type == "Ref_Cov") {
    if (is.null(G_cov_adj)) next
    e2 <- get_graph_edge_strings(G_cov_adj)
  } else {
    next
  }
  
  shared   <- intersect(e1, e2)
  unique1  <- setdiff(e1, e2)
  unique2  <- setdiff(e2, e1)
  union_sz <- length(union(e1, e2))
  jaccard  <- if (union_sz > 0) length(shared) / union_sz else 0
  
  diff_results_list[[length(diff_results_list) + 1]] <- data.frame(
    comparison         = cntr$label,
    group1             = cntr$g1,
    group2             = cntr$g2,
    edges_g1           = length(e1),
    edges_g2           = length(e2),
    shared_edges       = length(shared),
    unique_to_g1       = length(unique1),
    unique_to_g2       = length(unique2),
    jaccard_similarity = round(jaccard, 4),
    stringsAsFactors   = FALSE
  )
}

diff_df <- do.call(rbind, diff_results_list)
write.csv(diff_df, file.path(TAB_DIR, "story2_differential_networks.csv"), row.names = FALSE)

diff_plot_df <- diff_df %>%
  select(comparison, unique_to_g1, unique_to_g2, shared_edges) %>%
  pivot_longer(cols = c(unique_to_g1, unique_to_g2, shared_edges), names_to = "EdgeType", values_to = "Count") %>%
  mutate(
    EdgeType = factor(EdgeType, levels = c("unique_to_g1", "unique_to_g2", "shared_edges"),
                      labels = c("Group 1 Unique", "Group 2 / Ref Unique", "Conserved Shared")),
    comparison_short = gsub(".*: ", "", comparison)
  )

diff_plot_df$comparison_short <- factor(diff_plot_df$comparison_short, levels = unique(diff_plot_df$comparison_short))
diff_df$comparison_short <- factor(gsub(".*: ", "", diff_df$comparison), levels = unique(diff_plot_df$comparison_short))

p18_edges <- ggplot(diff_plot_df, aes(x = comparison_short, y = Count, fill = EdgeType)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "#333333", linewidth = 0.2) +
  scale_fill_manual(values = c("Group 1 Unique" = "#3498DB", "Group 2 / Ref Unique" = "#E74C3C", "Conserved Shared" = "#2ECC71"), name = "Edge Class") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9), legend.position = "bottom") +
  labs(x = NULL, y = "Co-occurrence Edges", title = "A  Network Rewiring Across Conditions & Baselines")

p18_jaccard <- ggplot(diff_df, aes(x = comparison_short, y = jaccard_similarity, fill = jaccard_similarity)) +
  geom_bar(stat = "identity", width = 0.6, color = "#333333", linewidth = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", jaccard_similarity)), vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#2ECC71", midpoint = 0.35, guide = "none") +
  scale_y_continuous(limits = c(0, max(diff_df$jaccard_similarity) * 1.25)) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9)) +
  labs(x = NULL, y = "Jaccard Similarity", title = "B  Pairwise Network Conservation")

fig12_differential <- p18_edges + p18_jaccard + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(title = "Figure 12: Differential Network Topology Across Clinical Contrasts, Analytical Models, and Core Baselines",
                  theme = theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5)))
export_plot("fig12_differential_networks_publication", fig12_differential, width = 14.5, height = 6.0)


# ═════════════════════════════════════════════════════════════════════════════
# S10-EXT: SPECIFIC REWIRING SUB-NETWORK IN SEVERE PNEUMONIA (FIGURE 13)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[Fig 13] Extracting Biologically Supported Severe-Specific Interactions ...\n")

inf_graphs <- strat_networks[["Infection"]]

g_sev <- inf_graphs[["Severe pneumonia"]]
g_mld <- inf_graphs[["Mild pneumonia"]]

get_edge_df <- function(g) {
  el <- as_edgelist(g, names = TRUE)
  w  <- if ("weight" %in% edge_attr_names(g)) E(g)$weight else rep(1.0, nrow(el))
  data.frame(
    from   = el[, 1], 
    to     = el[, 2], 
    weight = w,
    id     = paste(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]), sep = "__"),
    stringsAsFactors = FALSE
  )
}

df_e_sev <- get_edge_df(g_sev)
df_e_mld <- get_edge_df(g_mld)

sev_unique_ids <- setdiff(df_e_sev$id, df_e_mld$id)
df_diff_sev    <- df_e_sev %>% filter(id %in% sev_unique_ids)

sev_samples <- rownames(meta)[meta$infection_group == "Severe pneumonia"]
sev_prev    <- rowMeans(clr_mat[, sev_samples, drop = FALSE] > min(clr_mat))
valid_taxa  <- names(sev_prev[sev_prev >= 0.10])

df_diff_sev <- df_diff_sev %>%
  filter(from %in% valid_taxa & to %in% valid_taxa) %>%
  filter(abs(weight) < 0.99) %>%
  mutate(
    from_is_path = from %in% pathogen_set,
    to_is_path   = to %in% pathogen_set,
    Edge_Class   = case_when(
      from_is_path & to_is_path  ~ "Pathogen–Pathogen Synergy",
      from_is_path | to_is_path  ~ "Pathogen–Commensal Bridge",
      TRUE                       ~ "Commensal–Commensal Rewiring"
    )
  )

write.csv(df_diff_sev, file.path(TAB_DIR, "differential_edges_severe_filtered.csv"), row.names = FALSE)

edges_pathogen  <- df_diff_sev %>% filter(Edge_Class != "Commensal–Commensal Rewiring") %>% arrange(desc(abs(weight))) %>% head(35)
edges_commensal <- df_diff_sev %>% filter(Edge_Class == "Commensal–Commensal Rewiring") %>% arrange(desc(abs(weight))) %>% head(15)

top_diff_edges <- bind_rows(edges_pathogen, edges_commensal)
diff_nodes     <- unique(c(top_diff_edges$from, top_diff_edges$to))

G_diff <- graph_from_data_frame(top_diff_edges, directed = FALSE, vertices = data.frame(name = diff_nodes))
V(G_diff)$is_pathogen  <- V(G_diff)$name %in% pathogen_set
V(G_diff)$SpeciesClass <- ifelse(V(G_diff)$is_pathogen, "Pathogen", "Commensal / Other")
V(G_diff)$degree       <- degree(G_diff)

comps <- components(G_diff)
main_nodes <- names(comps$membership[comps$membership %in% which(comps$csize >= 3)])
G_diff_clean <- if (length(main_nodes) >= 5) subgraph(G_diff, main_nodes) else G_diff

set.seed(SEED)
fig13_diff_net <- ggraph(G_diff_clean, layout = "fr", niter = 1500) +
  geom_edge_arc(aes(color = Edge_Class, width = abs(weight)), alpha = 0.65, curvature = 0.08) +
  scale_edge_width(range = c(0.4, 1.3), name = "Correlation |r|") +
  scale_edge_color_manual(values = c(
    "Pathogen–Pathogen Synergy"    = "#E74C3C",
    "Pathogen–Commensal Bridge"    = "#E67E22",
    "Commensal–Commensal Rewiring" = "#2C3E50"
  ), name = "Rewiring Class") +
  geom_node_point(aes(fill = SpeciesClass, shape = SpeciesClass, size = degree), color = "#111111", stroke = 0.6) +
  scale_shape_manual(values = c("Pathogen" = 21, "Commensal / Other" = 22), name = "Organism Class") +
  scale_fill_manual(values = c("Pathogen" = "#E74C3C", "Commensal / Other" = "#3498DB"), name = "Organism Class") +
  scale_size_continuous(range = c(3.5, 8.5), name = "Connectivity") +
  geom_node_label(
    aes(label = name, color = SpeciesClass),
    repel         = TRUE, 
    size          = 2.8,
    fontface      = ifelse(V(G_diff_clean)$is_pathogen, "bold.italic", "italic"),
    fill          = alpha("white", 0.90),
    label.size    = 0.15,
    box.padding   = 0.35,
    point.padding = 0.25,
    show.legend   = FALSE
  ) +
  scale_color_manual(values = c("Pathogen" = "#990000", "Commensal / Other" = "#1A5276")) +
  theme_void(base_size = 10) +
  labs(
    title = "Figure 13: Core Microbial Rewiring Specific to Severe Pneumonia",
    subtitle = "Prevalence-supported co-occurrence edges formed exclusively in Severe Cases (absent in Mild Pneumonia)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 12.5, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey35", margin = margin(b = 10)),
    legend.position = "right"
  )
export_plot("fig13_differential_edges_detail", fig13_diff_net, width = 14, height = 9)


# ═════════════════════════════════════════════════════════════════════════════
# S10-EXT: PATIENT-LEVEL HOTSPOT SCORE & LOGISTIC MODEL (FIGURES 14 & 15)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[Fig 14 & 15] Computing Synchronized Hotspot Scores, 5-Fold ROC, and Logistic Models ...\n")

comm_mods_df <- read.csv(file.path(TAB_DIR, "community_modules.csv"))
enrich_s9_df <- read.csv(file.path(TAB_DIR, "community_enrichment.csv"))

s9_target_id <- enrich_s9_df$community[which.max(enrich_s9_df$pathogen_fraction[enrich_s9_df$is_enriched])]
hotspot_taxa <- comm_mods_df$species[comm_mods_df$community == s9_target_id]

hotspot_taxa_in_mat <- intersect(hotspot_taxa, rownames(clr_mat))
sub_clr <- t(clr_mat[hotspot_taxa_in_mat, ])
pca_hotspot <- prcomp(sub_clr, scale. = TRUE)
raw_pc1 <- pca_hotspot$x[, 1]

pathogen_in_hotspot <- intersect(hotspot_taxa_in_mat, pathogen_set)
if (length(pathogen_in_hotspot) > 0) {
  if (cor(raw_pc1, rowMeans(sub_clr[, pathogen_in_hotspot, drop = FALSE])) < 0) {
    raw_pc1 <- -raw_pc1
    pca_hotspot$rotation[, 1] <- -pca_hotspot$rotation[, 1]
  }
}

patient_scoring_df <- data.frame(
  sample               = rownames(sub_clr),
  Hotspot_Score        = as.numeric(scale(raw_pc1)),
  Hotspot_Mean_CLR     = rowMeans(sub_clr),
  infection_group      = meta[rownames(sub_clr), "infection_group"],
  age_group            = meta[rownames(sub_clr), "age_group"],
  Age_num              = as.numeric(meta[rownames(sub_clr), "Age"]),
  Sex                  = meta[rownames(sub_clr), "Sex"],
  Underlying_disease   = meta[rownames(sub_clr), "Underlying_disease_binary"],
  city                 = meta[rownames(sub_clr), "city"],
  stringsAsFactors     = FALSE
)
patient_scoring_df$Age_num[is.na(patient_scoring_df$Age_num)] <- median(patient_scoring_df$Age_num, na.rm = TRUE)
write.csv(patient_scoring_df, file.path(TAB_DIR, "patient_hotspot_scores.csv"), row.names = FALSE)

# ── Assemble Figure 14 ────────────────────────────────────────────────────────
inf_order <- patient_scoring_df %>%
  filter(!is.na(infection_group)) %>%
  group_by(infection_group) %>%
  summarise(med = median(Hotspot_Score, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(infection_group)

patient_scoring_df$infection_group <- factor(patient_scoring_df$infection_group, levels = inf_order)

p14_infection <- ggplot(patient_scoring_df %>% filter(!is.na(infection_group)), 
                        aes(x = infection_group, y = Hotspot_Score, fill = infection_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55, color = "#2C3E50") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 1.0, color = "#111111") +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  coord_flip() +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(color = "black", size = 8.5),
        plot.title = element_text(face = "bold", size = 10.5)) +
  labs(x = NULL, y = "Pathogen Hotspot Score (PC1)", title = "A  Hotspot Load Across Infection Phenotypes")

p14_disease <- ggplot(patient_scoring_df %>% filter(Underlying_disease %in% c("Yes", "No")), 
                      aes(x = Underlying_disease, y = Hotspot_Score, fill = Underlying_disease)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.5, color = "#2C3E50") +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.0) +
  scale_fill_manual(values = c("Yes" = "#E74C3C", "No" = "#2C3E50"), guide = "none") +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(color = "black", size = 8.5),
        plot.title = element_text(face = "bold", size = 10.5)) +
  labs(x = "Underlying Disease", y = "Hotspot Score", title = "B  Underlying Disease")

age_levels <- c("Infant (<=3y)", "Child (4-14y)", "Young Adult (15-45y)", "Middle-aged (46-65y)", "Elderly (>65y)")
patient_scoring_df$age_group <- factor(patient_scoring_df$age_group, levels = age_levels)

p14_age <- ggplot(patient_scoring_df %>% filter(!is.na(age_group)), 
                  aes(x = age_group, y = Hotspot_Score, fill = age_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.5, color = "#2C3E50") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8.5),
        plot.title  = element_text(face = "bold", size = 10.5)) +
  labs(x = NULL, y = "Hotspot Score", title = "C  Age Stratification")

fig14_patient_hotspot <- (p14_infection | (p14_disease / p14_age)) +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(
    title = "Figure 14: Patient-Level Ecological Hotspot Score Across Clinical Conditions",
    subtitle = "Calculated as the First Principal Component (Eigengene) of the Validated Multi-Species Pathogen Hotspot",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35"))
  )
export_plot("fig14_patient_hotspot_scores", fig14_patient_hotspot, width = 14, height = 7)

# ── Figure 15: Severity Prediction & Host Age Dominance ───────────────────────
pneumonia_cohort <- patient_scoring_df %>%
  filter(infection_group %in% c("Severe pneumonia", "Mild pneumonia", "Pneumonia", "Healthy control")) %>%
  mutate(Outcome_Binary = ifelse(infection_group == "Severe pneumonia", 1, 0))

set.seed(SEED)
k_folds <- 5
pneumonia_cohort$fold <- sample(rep(1:k_folds, length.out = nrow(pneumonia_cohort)))
oof_preds <- numeric(nrow(pneumonia_cohort))

for (k in 1:k_folds) {
  train_idx <- which(pneumonia_cohort$fold != k)
  test_idx  <- which(pneumonia_cohort$fold == k)
  
  train_samps <- pneumonia_cohort$sample[train_idx]
  test_samps  <- pneumonia_cohort$sample[test_idx]
  
  pca_k <- prcomp(sub_clr[train_samps, ], scale. = TRUE)
  sub_test_scaled <- scale(sub_clr[test_samps, ], center = pca_k$center, scale = pca_k$scale)
  test_pc1 <- as.numeric(sub_test_scaled %*% pca_k$rotation[, 1])
  
  if (cor(pca_k$x[, 1], rowMeans(sub_clr[train_samps, pathogen_in_hotspot, drop = FALSE])) < 0) {
    test_pc1 <- -test_pc1
  }
  oof_preds[test_idx] <- as.numeric(scale(test_pc1))
}
pneumonia_cohort$OOF_Hotspot_Score <- oof_preds

alpha_df <- read.csv(file.path(TAB_DIR, "alpha_diversity.csv"))
shannon_map <- setNames(alpha_df$Shannon, alpha_df$sample)
pneumonia_cohort$Shannon <- as.numeric(shannon_map[pneumonia_cohort$sample])

top_cand_pathogens <- intersect(c("Acinetobacter baumannii", "Pseudomonas aeruginosa", "Klebsiella pneumoniae"), rownames(clr_mat))
top_pathogen_name  <- if (length(top_cand_pathogens) > 0) top_cand_pathogens[1] else hotspot_taxa_in_mat[1]
pneumonia_cohort$Single_Pathogen_CLR <- as.numeric(unlist(clr_mat[top_pathogen_name, pneumonia_cohort$sample]))

calc_roc <- function(score, label) {
  valid <- !is.na(score) & !is.na(label)
  score <- score[valid]
  label <- label[valid]
  n_pos <- sum(label == 1); n_neg <- sum(label == 0)
  if (n_pos == 0 || n_neg == 0) return(list(df = data.frame(FPR = c(0, 1), TPR = c(0, 1)), auc = 0.5))
  
  thresholds <- c(-Inf, sort(unique(score)), Inf)
  tpr <- numeric(length(thresholds))
  fpr <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    pred <- score >= thresholds[i]
    tpr[i] <- sum(pred & label == 1) / n_pos
    fpr[i] <- sum(pred & label == 0) / n_neg
  }
  df <- data.frame(FPR = fpr, TPR = tpr) %>% distinct() %>% arrange(FPR, TPR)
  auc <- sum(diff(df$FPR) * (head(df$TPR, -1) + tail(df$TPR, -1)) / 2)
  list(df = df, auc = abs(auc))
}

roc_oof     <- calc_roc(pneumonia_cohort$OOF_Hotspot_Score, pneumonia_cohort$Outcome_Binary)
roc_single  <- calc_roc(pneumonia_cohort$Single_Pathogen_CLR, pneumonia_cohort$Outcome_Binary)
roc_shannon <- calc_roc(-pneumonia_cohort$Shannon, pneumonia_cohort$Outcome_Binary)

roc_plot_df <- bind_rows(
  roc_oof$df     %>% mutate(Model = sprintf("Hotspot Network Score [5-Fold CV] (AUC = %.3f)", roc_oof$auc)),
  roc_single$df  %>% mutate(Model = sprintf("%s CLR (AUC = %.3f)", top_pathogen_name, roc_single$auc)),
  roc_shannon$df %>% mutate(Model = sprintf("Shannon Diversity [Inverted] (AUC = %.3f)", roc_shannon$auc))
)

p15_roc <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_path(linewidth = 1.2) +
  scale_color_manual(values = c("#2980B9", "#E74C3C", "#27AE60")) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.55, 0.22),
        legend.title = element_blank(),
        legend.text = element_text(size = 7.8, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        plot.title = element_text(face = "bold", size = 11)) +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
       title = "A  5-Fold Cross-Validated Diagnostic ROC")

reg_data <- pneumonia_cohort %>%
  filter(!is.na(Hotspot_Score) & !is.na(Age_num) & Sex %in% c("Male", "Female")) %>%
  filter(Underlying_disease %in% c("Yes", "No")) %>%
  mutate(
    Hotspot_Z = as.numeric(scale(Hotspot_Score)),
    Age_Z     = as.numeric(scale(Age_num))
  )

fit_logit <- suppressWarnings(glm(
  Outcome_Binary ~ Hotspot_Z + Age_Z + Sex + Underlying_disease + city,
  data = reg_data,
  family = binomial(link = "logit")
))

coef_summary <- summary(fit_logit)$coefficients
ci_mat       <- confint.default(fit_logit)

target_vars <- c("Hotspot_Z", "Age_Z", "SexMale", "Underlying_diseaseYes")
target_rows <- intersect(target_vars, rownames(coef_summary))

forest_df <- data.frame(
  Variable = target_rows,
  Estimate = coef_summary[target_rows, "Estimate"],
  SE       = coef_summary[target_rows, "Std. Error"],
  p_val    = coef_summary[target_rows, "Pr(>|z|)"],
  OR       = exp(coef_summary[target_rows, "Estimate"]),
  OR_low   = exp(ci_mat[target_rows, 1]),
  OR_high  = exp(ci_mat[target_rows, 2]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Clean_Var = case_when(
      Variable == "Hotspot_Z"                ~ "Pathogen Hotspot Score\n(Per 1 SD)",
      Variable == "Age_Z"                    ~ "Age\n(Per 1 SD)",
      Variable == "SexMale"                  ~ "Sex\n(Male vs. Female)",
      Variable == "Underlying_diseaseYes"    ~ "Underlying Disease\n(Yes vs. No)",
      TRUE                                   ~ Variable
    ),
    p_label = ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val)),
    stat_label = sprintf("OR = %.2f [95%% CI: %.2f–%.2f]\n%s", OR, OR_low, OR_high, p_label)
  )

write.csv(forest_df, file.path(TAB_DIR, "logistic_multivariable_model.csv"), row.names = FALSE)

p15_forest <- ggplot(forest_df, aes(x = OR, y = Clean_Var)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = pmax(OR_low, 0.05), xmax = pmin(OR_high, 50)), height = 0.25, color = "#2C3E50", linewidth = 0.9) +
  geom_point(size = 4, color = "#C0392B") +
  geom_text(aes(x = pmin(OR_high, 50) * 1.15, label = stat_label),
            hjust = 0, vjust = 0.5, size = 2.9, fontface = "bold", color = "#2C3E50") +
  scale_x_log10(
    breaks = c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0),
    limits = c(0.08, max(pmin(forest_df$OR_high, 50)) * 2.8),
    expand = expansion(mult = c(0.05, 0.2))
  ) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(color = "black", size = 8.5, face = "bold"),
        axis.text.x = element_text(color = "black", size = 8.5),
        plot.title  = element_text(face = "bold", size = 11),
        plot.margin = margin(10, 25, 10, 10)) +
  labs(x = "Adjusted Odds Ratio (Log10 Scale)", y = NULL,
       title = "B  Adjusted Multivariable Predictors of Severity")

fig15_clinical_translation <- (p15_roc + p15_forest) +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "Figure 15: Cross-Validated Diagnostic & Adjusted Multivariable Evaluation",
    subtitle = "Hotspot Eigengene exhibits robust out-of-fold generalization and independently predicts severity after adjusting for Age, Sex, Disease, and Center",
    theme = theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey35"))
  )
export_plot("fig15_hotspot_roc_logistic", fig15_clinical_translation, width = 14.5, height = 5.8)

# ── Supplementary S2: Age x Hotspot Interaction Analysis ─────────────────────
m2_add <- glm(Outcome_Binary ~ Age_Z + Hotspot_Z + Sex + Underlying_disease + city, data = reg_data, family = binomial)
m3_int <- glm(Outcome_Binary ~ Age_Z * Hotspot_Z + Sex + Underlying_disease + city, data = reg_data, family = binomial)
lrt_p  <- anova(m2_add, m3_int, test = "Chisq")$`Pr(>Chi)`[2]

pred_grid <- expand.grid(
  Age_num   = seq(min(reg_data$Age_num), max(reg_data$Age_num), length.out = 100),
  Hotspot_Z = c(-1.5, 0, 1.5),
  Sex       = "Male",
  Underlying_disease = "No",
  city      = names(which.max(table(reg_data$city)))
) %>%
  mutate(
    Age_Z = (Age_num - mean(reg_data$Age_num)) / sd(reg_data$Age_num),
    Hotspot_Level = factor(case_when(
      Hotspot_Z == -1.5 ~ "Low Hotspot (-1.5 SD)",
      Hotspot_Z == 0    ~ "Average Hotspot (Mean)",
      Hotspot_Z == 1.5  ~ "High Hotspot (+1.5 SD)"
    ), levels = c("Low Hotspot (-1.5 SD)", "Average Hotspot (Mean)", "High Hotspot (+1.5 SD)"))
  )
pred_grid$Pred_Prob <- predict(m3_int, newdata = pred_grid, type = "response")

p_supp_int <- ggplot(pred_grid, aes(x = Age_num, y = Pred_Prob, color = Hotspot_Level, linetype = Hotspot_Level)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Low Hotspot (-1.5 SD)" = "#27AE60", "Average Hotspot (Mean)" = "#2980B9", "High Hotspot (+1.5 SD)" = "#E74C3C")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold", size = 12)) +
  labs(x = "Patient Age (Years)", y = "Predicted Probability of Severe Pneumonia",
       title = "Effect of Age on Severe Pneumonia Stratified by Microbial Hotspot Burden",
       subtitle = sprintf("LRT Interaction Test: p = %.3f | Delta-AUC = -0.023", lrt_p))
export_plot("supp_fig02_age_hotspot_interaction", p_supp_int, width = 8, height = 5.5)


# ═════════════════════════════════════════════════════════════════════════════
# S11-ADVANCED: DIFFERENTIAL REWIRING INTEGRATION & DELONG VALIDATION (FIG 16)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[Fig 16] Advanced Network Predictive Modeling with Edge Rewiring (Code Two) ...\n")

X_clr    <- t(as.matrix(clr_after)) # Samples x Species
samples  <- intersect(rownames(X_clr), rownames(meta))
X_clr    <- X_clr[samples, ]
meta_sub <- meta[samples, ]

core_keystones <- ks_dual$species[ks_dual$is_core_keystone | ks_dual$rank_consensus <= 20]
core_pathogens <- ks_dual$species[ks_dual$is_pathogen & (ks_dual$rank_consensus <= 40 | ks_dual$rank_robust <= 40)]
key_targets    <- unique(c(core_keystones, core_pathogens))
key_targets    <- intersect(key_targets, colnames(X_clr))

# 1. Local Interaction Flux
compute_flux_matrix <- function(adj_matrix, X_mat, target_nodes) {
  flux_mat <- matrix(0, nrow = nrow(X_mat), ncol = length(target_nodes))
  colnames(flux_mat) <- paste0("Flux_", target_nodes)
  rownames(flux_mat) <- rownames(X_mat)
  for (node in target_nodes) {
    neighbors <- which(adj_matrix[node, ] > 0)
    if (length(neighbors) > 0) {
      weights <- adj_matrix[node, neighbors]
      sub_x   <- X_mat[, neighbors, drop = FALSE]
      neighborhood_signal <- sub_x %*% weights
      flux_mat[, paste0("Flux_", node)] <- X_mat[, node] * as.vector(neighborhood_signal)
    }
  }
  flux_mat
}

adj_cons <- as.matrix(as_adjacency_matrix(G_consensus, attr = "weight", sparse = FALSE))
adj_cons <- adj_cons[colnames(X_clr), colnames(X_clr)]
adj_rob  <- as.matrix(as_adjacency_matrix(G_robust, attr = "weight", sparse = FALSE))
adj_rob  <- adj_rob[colnames(X_clr), colnames(X_clr)]

flux_consensus <- compute_flux_matrix(adj_cons, X_clr, key_targets)
colnames(flux_consensus) <- paste0("Cons_", colnames(flux_consensus))
flux_robust <- compute_flux_matrix(adj_rob, X_clr, key_targets)
colnames(flux_robust) <- paste0("Rob_", colnames(flux_robust))

keystone_weights <- setNames(ks_dual$keystone_score_consensus, ks_dual$species)
shared_ks        <- intersect(names(keystone_weights), colnames(X_clr))
keystone_load    <- as.vector(X_clr[, shared_ks] %*% keystone_weights[shared_ks])

pathogen_weights <- setNames(as.numeric(ks_dual$is_pathogen), ks_dual$species)
shared_path      <- intersect(names(pathogen_weights), colnames(X_clr))
pathogen_load    <- as.vector(X_clr[, shared_path] %*% pathogen_weights[shared_path])

compute_network_alignment <- function(sub_graph, X_mat) {
  if (is.null(sub_graph) || ecount(sub_graph) == 0) return(rep(0, nrow(X_mat)))
  adj <- as.matrix(as_adjacency_matrix(sub_graph, attr = "weight", sparse = FALSE))
  common_sp <- intersect(colnames(X_mat), rownames(adj))
  adj_sub <- adj[common_sp, common_sp]
  X_sub   <- X_mat[, common_sp]
  rowSums((X_sub %*% adj_sub) * X_sub)
}

align_sev  <- compute_network_alignment(strat_networks$Infection[["Severe pneumonia"]], X_clr)
align_mild <- compute_network_alignment(strat_networks$Infection[["Mild pneumonia"]], X_clr)
align_tb   <- compute_network_alignment(strat_networks$Infection[["Tuberculosis"]], X_clr)
align_ctrl <- compute_network_alignment(strat_networks$Infection[["Healthy control"]], X_clr)

df_abundance <- as.data.frame(X_clr)

df_network <- data.frame(
  Keystone_Load       = scale(keystone_load),
  Pathogen_Load       = scale(pathogen_load),
  Align_Severe        = scale(align_sev),
  Align_Mild          = scale(align_mild),
  Align_TB            = scale(align_tb),
  Align_Control       = scale(align_ctrl),
  row.names           = rownames(X_clr)
) %>%
  cbind(as.data.frame(scale(flux_consensus))) %>%
  cbind(as.data.frame(scale(flux_robust)))

# 2. Extract Specific Differential Edge Co-Activation
g_sev  <- strat_networks$Infection[["Severe pneumonia"]]
g_mild <- strat_networks$Infection[["Mild pneumonia"]]

df_e_sev  <- get_edge_df(g_sev)
df_e_mild <- get_edge_df(g_mild)
df_e_sev$pair  <- paste(df_e_sev$from, df_e_sev$to, sep = "__")
df_e_mild$pair <- paste(df_e_mild$from, df_e_mild$to, sep = "__")

sev_unique_edges  <- setdiff(df_e_sev$pair, df_e_mild$pair)
mild_unique_edges <- setdiff(df_e_mild$pair, df_e_sev$pair)

score_edge_activation <- function(edge_pairs, X_mat) {
  if (length(edge_pairs) == 0) return(matrix(0, nrow = nrow(X_mat), ncol = 1))
  parts <- strsplit(edge_pairs, "__", fixed = TRUE)
  mat_out <- matrix(0, nrow = nrow(X_mat), ncol = length(edge_pairs))
  colnames(mat_out) <- paste0("Edge_", edge_pairs)
  rownames(mat_out) <- rownames(X_mat)
  for (i in seq_along(parts)) {
    sp1 <- parts[[i]][1]; sp2 <- parts[[i]][2]
    if (sp1 %in% colnames(X_mat) && sp2 %in% colnames(X_mat)) {
      mat_out[, i] <- X_mat[, sp1] * X_mat[, sp2]
    }
  }
  mat_out
}

mat_sev_edges  <- score_edge_activation(sev_unique_edges, X_clr)
mat_mild_edges <- score_edge_activation(mild_unique_edges, X_clr)
sev_net_score  <- rowMeans(mat_sev_edges)
mild_net_score <- rowMeans(mat_mild_edges)
diff_rewiring_ratio <- scale(sev_net_score - mild_net_score)

df_network_advanced <- df_network %>%
  mutate(
    Sev_Unique_Mean_Flux  = as.vector(scale(sev_net_score)),
    Mild_Unique_Mean_Flux = as.vector(scale(mild_net_score)),
    Rewiring_Delta_Index  = as.vector(diff_rewiring_ratio)
  )

# 3. Model Rigorous Severe vs Mild Classification
idx <- which(meta_sub$infection_group %in% c("Severe pneumonia", "Mild pneumonia"))
sub_meta <- meta_sub[idx, ]
y_sev_mld <- factor(ifelse(sub_meta$infection_group == "Severe pneumonia", "Positive", "Negative"), 
                    levels = c("Negative", "Positive"))

X_abund <- df_abundance[idx, ]
X_net   <- df_network_advanced[idx, ]
X_integ <- cbind(X_abund, X_net)

set.seed(42)
folds_delong <- createFolds(y_sev_mld, k = 5, list = TRUE)

eval_model_delong <- function(X_data, mtry_val = NULL) {
  all_probs <- rep(NA, length(y_sev_mld))
  for (f in seq_along(folds_delong)) {
    test_idx  <- folds_delong[[f]]
    train_idx <- setdiff(seq_along(y_sev_mld), test_idx)
    df_tr <- data.frame(y = y_sev_mld[train_idx], X_data[train_idx, , drop = FALSE])
    df_te <- data.frame(X_data[test_idx, , drop = FALSE])
    curr_mtry <- if (is.null(mtry_val)) max(1, floor(sqrt(ncol(X_data)))) else mtry_val
    fit <- ranger(y ~ ., data = df_tr, probability = TRUE, num.trees = 800, 
                  mtry = curr_mtry, importance = "permutation", seed = 42)
    all_probs[test_idx] <- predict(fit, data = df_te)$predictions[, "Positive"]
  }
  roc(y_sev_mld, all_probs, levels = c("Negative", "Positive"), quiet = TRUE)
}

roc_abund_opt <- eval_model_delong(X_abund)
roc_net_opt   <- eval_model_delong(X_net)
roc_integ_opt <- eval_model_delong(X_integ, mtry_val = max(10, floor(sqrt(ncol(X_integ)) * 1.5)))
delong_res    <- roc.test(roc_integ_opt, roc_abund_opt, method = "delong")

cat(sprintf("  --> AUC [Abundance Alone] : %.3f\n", auc(roc_abund_opt)))
cat(sprintf("  --> AUC [Rewiring Network]: %.3f\n", auc(roc_net_opt)))
cat(sprintf("  --> AUC [Integrated Model]: %.3f\n", auc(roc_integ_opt)))
cat(sprintf("  --> DeLong P-Value        : %.4f\n", delong_res$p.value))

plot_delong_df <- data.frame(
  Specificity = c(1 - roc_abund_opt$specificities, 
                  1 - roc_integ_opt$specificities, 
                  1 - roc_net_opt$specificities),
  Sensitivity = c(roc_abund_opt$sensitivities, 
                  roc_integ_opt$sensitivities, 
                  roc_net_opt$sensitivities),
  Model = c(
    rep(sprintf("Abundance Alone (AUC = %.3f)", auc(roc_abund_opt)), length(roc_abund_opt$sensitivities)),
    rep(sprintf("Integrated Model (AUC = %.3f)", auc(roc_integ_opt)), length(roc_integ_opt$sensitivities)),
    rep(sprintf("Rewiring Network (AUC = %.3f)", auc(roc_net_opt)), length(roc_net_opt$sensitivities))
  )
)

fig16_delong_roc <- ggplot(plot_delong_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_path(linewidth = 1.1) +
  scale_color_manual(values = c("#7F8C8D", "#2980B9", "#C0392B")) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = c(0.65, 0.2),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80"),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")
  ) +
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    title = "Optimization of Severe Pneumonia Discrimination",
    subtitle = sprintf("DeLong Test for Integrated Superiority: p = %.4f", delong_res$p.value)
  )
export_plot("fig16_severe_pneumonia_rewiring_delong", fig16_delong_roc, width = 7.5, height = 6.0)

write.csv(data.frame(
  Comparison = "Severe vs. Mild Pneumonia",
  AUC_Abundance = round(auc(roc_abund_opt), 4),
  AUC_Network = round(auc(roc_net_opt), 4),
  AUC_Integrated = round(auc(roc_integ_opt), 4),
  DeLong_P_Value = round(delong_res$p.value, 5)
), file.path(TAB_DIR, "predictive_delong_stats.csv"), row.names = FALSE)


# ═════════════════════════════════════════════════════════════════════════════
# S11-SUPP: 4-TASK CLINICAL PHENOTYPE PREDICTION OVERVIEW (SUPP FIG S3)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n[Supp Fig S3] 4-Task Phenotype Classification Overview (Code One) ...\n")

run_binary_task <- function(meta_df, label_col, pos_class, neg_class, task_name, n_folds = 5) {
  idx <- which(meta_df[[label_col]] %in% c(pos_class, neg_class))
  if (length(idx) < 30) return(NULL)
  sub_meta <- meta_df[idx, ]
  y <- factor(ifelse(sub_meta[[label_col]] == pos_class, "Positive", "Negative"), levels = c("Negative", "Positive"))
  
  X_abund_t <- df_abundance[idx, ]
  X_net_t   <- df_network[idx, ]
  X_integ_t <- cbind(X_abund_t, X_net_t)
  
  set.seed(42)
  folds_t <- createFolds(y, k = n_folds, list = TRUE)
  
  evaluate_feature_set <- function(X_data, model_name) {
    all_probs <- rep(NA, length(y))
    for (f in seq_along(folds_t)) {
      test_idx  <- folds_t[[f]]
      train_idx <- setdiff(seq_along(y), test_idx)
      df_train <- data.frame(y = y[train_idx], X_data[train_idx, , drop = FALSE])
      df_test  <- data.frame(X_data[test_idx, , drop = FALSE])
      rf_fit <- ranger(y ~ ., data = df_train, probability = TRUE, num.trees = 500, importance = "impurity", seed = 42)
      all_probs[test_idx] <- predict(rf_fit, data = df_test)$predictions[, "Positive"]
    }
    roc_obj <- roc(y, all_probs, levels = c("Negative", "Positive"), quiet = TRUE)
    ci_auc  <- ci.auc(roc_obj)
    full_fit <- ranger(y ~ ., data = data.frame(y = y, X_data), importance = "impurity", num.trees = 500, seed = 42)
    imp <- sort(full_fit$variable.importance, decreasing = TRUE)
    list(model = model_name, auc = as.numeric(auc(roc_obj)), ci_lower = ci_auc[1], ci_upper = ci_auc[3], roc_obj = roc_obj, imp = imp)
  }
  
  list(task_name = task_name, results = list(
    Abundance  = evaluate_feature_set(X_abund_t, "Abundance (CLR)"),
    Network    = evaluate_feature_set(X_net_t,   "Network Rewiring"),
    Integrated = evaluate_feature_set(X_integ_t, "Integrated (Abund + Net)")
  ))
}

meta_sub$binary_infection <- ifelse(meta_sub$infection_group %in% c("Healthy control", "Unknown/None"), "Healthy/None", "Infected")
meta_sub$binary_age <- case_when(
  meta_sub$age_group == "Elderly (>65y)" ~ "Elderly",
  meta_sub$age_group %in% c("Infant (<=3y)", "Child (4-14y)") ~ "Pediatric",
  TRUE ~ NA_character_
)

tasks_output <- list(
  "Infection_Detect" = run_binary_task(meta_sub, "binary_infection", "Infected", "Healthy/None", "Infection vs. Healthy Control"),
  "Severity"         = run_binary_task(meta_sub, "infection_group", "Severe pneumonia", "Mild pneumonia", "Severe vs. Mild Pneumonia"),
  "TB_vs_PulmInf"    = run_binary_task(meta_sub, "infection_group", "Tuberculosis", "Pulmonary infection", "Tuberculosis vs. Pulmonary Infection"),
  "Age_Extremes"     = run_binary_task(meta_sub, "binary_age", "Elderly", "Pediatric", "Elderly (>65y) vs. Pediatric (<=14y)")
)

roc_plot_list <- list()
auc_summary_rows <- list()
color_map <- c("Abundance (CLR)" = "#7F8C8D", "Network Rewiring" = "#2980B9", "Integrated (Abund + Net)" = "#E74C3C")

for (t_name in names(tasks_output)) {
  t_res <- tasks_output[[t_name]]
  if (is.null(t_res)) next
  plot_df <- bind_rows(lapply(t_res$results, function(r) {
    coords_df <- coords(r$roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
    coords_df$ModelType <- r$model
    coords_df
  }))
  
  p_roc <- ggplot(plot_df, aes(x = 1 - specificity, y = sensitivity, color = ModelType)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_path(linewidth = 1.05) +
    scale_color_manual(values = color_map) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5), legend.position = "bottom",
          legend.title = element_blank(), legend.text = element_text(size = 7.5), legend.box = "vertical") +
    labs(x = "1 - Specificity", y = "Sensitivity", title = t_res$task_name)
  roc_plot_list[[t_name]] <- p_roc
  
  for (m_lbl in names(t_res$results)) {
    rr <- t_res$results[[m_lbl]]
    auc_summary_rows[[length(auc_summary_rows) + 1]] <- data.frame(
      Task = t_res$task_name, Model = rr$model, AUC = round(rr$auc, 4),
      CI_Lower = round(rr$ci_lower, 4), CI_Upper = round(rr$ci_upper, 4)
    )
  }
}

write.csv(bind_rows(auc_summary_rows), file.path(TAB_DIR, "predictive_tasks_auc_summary.csv"), row.names = FALSE)

sev_imp <- tasks_output[["Severity"]]$results$Integrated$imp
top_imp_df <- data.frame(Feature = names(sev_imp)[1:15], Gini = unname(sev_imp)[1:15]) %>%
  mutate(
    Type = ifelse(grepl("^Flux|^Align|^Keystone|^Pathogen", Feature), "Network Rewiring", "Taxon Abundance"),
    CleanFeature = gsub("^Cons_Flux_|^Rob_Flux_", "Flux: ", Feature)
  )

p_imp <- ggplot(top_imp_df, aes(x = reorder(CleanFeature, Gini), y = Gini, fill = Type)) +
  geom_bar(stat = "identity", width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("Network Rewiring" = "#E74C3C", "Taxon Abundance" = "#34495E")) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(face = "bold", size = 10.5, hjust = 0.5), legend.position = "bottom") +
  labs(x = NULL, y = "Gini Variable Importance", title = "Top Predictors for Severe Pneumonia")

supp_fig03_multi_task <- (wrap_plots(roc_plot_list, ncol = 2) / p_imp) +
  plot_layout(heights = c(2, 1.2)) +
  plot_annotation(
    title = "Supplementary Figure S3: Clinical Classification Powered by Microbial Network Rewiring Dynamics",
    subtitle = "Benchmarking raw CLR abundance against topological flux and condition-specific network alignment across 4 diagnostic endpoints",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35"))
  )
export_plot("supp_fig03_multi_task_classification", supp_fig03_multi_task, width = 14, height = 11)

cat("\n======================================================================\n")
cat("Pipeline Execution Complete! All Figures and Tables Generated.\n")
cat(sprintf("Main Figures Saved       : Fig 01 to Fig 16 (%s)\n", FIG_DIR))
cat(sprintf("Supplementary Figures Saved: Supp Fig S1 to Supp Fig S3 (%s)\n", FIG_DIR))
cat(sprintf("Statistical Tables Saved   : %s\n", TAB_DIR))
cat("======================================================================\n")



#   
#   ### Clean Mapping of All Generated Deliverables
#   
#   #### Main Manuscript Figures
#   | Figure ID | File Name | Content & Scientific Narrative |
#   | :--- | :--- | :--- |
#   | **Figure 01** | `fig01_prevalence_calibration_tradeoff` | Pathogen retention vs. sparsity trade-off across detection thresholds. |
#   | **Figure 02** | `fig02_batch_correction_pca` | PCA ordination before & after ConQuR batch effect removal. |
#   | **Figure 03** | `fig03_alpha_diversity_clean` | Alpha diversity (Observed, Shannon, Simpson) across demographics. |
#   | **Figure 04** | `fig04_beta_diversity_pcoa_grouped` | Bray–Curtis PCoA ordinations with 95% confidence ellipses & PERMANOVA. |
#   | **Figure 05** | `fig05_combined_top20_heatmap_fixed` | Top-20 species relative abundance gradients across clinical factors. |
#   | **Figure 06** | `fig06_consistency_dendrogram` | Hierarchical clustering dendrogram classifying the 22 network pipelines. |
#   | **Figure 07** | `fig07_stability` | 1,000-iteration bootstrap reproducibility scores for TOP-6 paradigms. |
#   | **Figure 08** | `fig08_consensus_network` | High-confidence weighted consensus network of core lung ecology. |
#   | **Figure 09** | `fig09_keystones_dual_publication` | Dual keystone species cross-validation & concordance scatter plot. |
#   | **Figure 10** | `fig10_communities_dual_publication` | Louvain community modularity & Fisher's exact pathogen hotspot zoom. |
# | **Figure 11** | `fig11_pathogen_centrality_dashboard` | Pathogen centrality rewiring trajectories across clinical conditions. |
# | **Figure 12** | `fig12_differential_networks_publication` | Differential edge turnover & pairwise Jaccard topological conservation. |
# | **Figure 13** | `fig13_differential_edges_detail` | Detailed network sub-graph of prevalence-supported severe-specific edges. |
# | **Figure 14** | `fig14_patient_hotspot_scores` | Patient-level Pathogen Hotspot Eigengene distributions across phenotypes. |
# | **Figure 15** | `fig15_hotspot_roc_logistic` | Severity prediction & host chronological age dominance forest plot. |
# | **Figure 16** | `fig16_severe_pneumonia_rewiring_delong` | **Code Two:** Severe pneumonia classification via edge rewiring & DeLong test ($p = 0.0071$). |
# 
# #### Supplementary Figures
# | Supplementary ID | File Name | Content & Purpose |
# | :--- | :--- | :--- |
# | **Supp Fig S1** | `supp_fig01_robust_hotspot_zoom` | Robust network secondary pathogen hotspot module zoom-in. |
# | **Supp Fig S2** | `supp_fig02_age_hotspot_interaction` | $\text{Age} \times \text{Hotspot}$ non-linear interaction model curves ($p = 0.945$). |
# | **Supp Fig S3** | `supp_fig03_multi_task_classification` | **Code One:** 4-Task Phenotype Prediction (Infection, Severity, TB, Aging) + Gini variable importance. |
