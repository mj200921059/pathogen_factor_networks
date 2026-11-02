#!/usr/bin/env Rscript
# =============================================================================

# =============================================================================
# Complete replication pipeline for BALF mNGS microbiota network analysis.
# Cleaned of clinical category and filtered to exclude host/subgroup samples.
# Consider FDR pauj in network construciton. 
# Requirements:
#   install.packages(c("tidyverse", "vegan", "igraph", "glasso", "cluster"))
#
# Author:Jun Ma
# Date  : 2026-09-29
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(igraph)
  library(glasso)
  library(cluster)
  library(parallel)
  library(taxonomizr)
  library(patchwork) # For side-by-side ggplot composition
  library(RColorBrewer)
})


# obtain microbes profile in each sample
setwd("~/dataanalysis/projects/pathogen_factors/results/03kraken")

# Only read Bracken taxonomy files
kraken_files <- list.files(
  pattern = "_brackentaxa\\.tsv$",
  full.names = TRUE,
  recursive = FALSE
)

if (length(kraken_files) == 0) {
  stop("No *_brackentaxa.tsv files found in the current directory.")
}

# Helper: read one Bracken file and keep taxid + new_est_reads
read_bracken <- function(f) {
  x <- read.delim(
    f,
    header = TRUE,
    quote = "",
    sep = "\t",
    check.names = FALSE
  )
  
  # Prefer named columns if present; otherwise use original positions 2 and 6
  if (all(c("taxonomy_id", "new_est_reads") %in% names(x))) {
    x <- x[, c("taxonomy_id", "new_est_reads"), drop = FALSE]
  } else {
    x <- x[, c(2, 6), drop = FALSE]
  }
  
  sample_name <- sub("_brackentaxa\\.tsv$", "", basename(f))
  colnames(x) <- c("taxonomy_id", sample_name)
  
  x$taxonomy_id <- as.integer(x$taxonomy_id)
  x
}

# Read all Bracken files into a list
kraken_list <- lapply(kraken_files, read_bracken)

# Merge all sample tables by taxonomy_id
kraken_table <- Reduce(
  function(a, b) merge(a, b, by = "taxonomy_id", all = TRUE),
  kraken_list
)

kraken_table[is.na(kraken_table)] <- 0

# Taxids for taxonomy lookup
taxid <- data.frame(taxonomy_id = unique(kraken_table$taxonomy_id))

setwd("~/dataanalysis/projects/accessiontaxa")
taxa <- getTaxonomy(taxid$taxonomy_id, "accessionTaxa.sql")
taxa <- as.data.frame(taxa)
taxa$taxonomy_id <- as.integer(rownames(taxa))

kraken_taxa <- merge(taxa, kraken_table, by = "taxonomy_id")

setwd("~/dataanalysis/projects/pathogen_factors")
save(kraken_taxa, file = "kraken_taxa.rda")
