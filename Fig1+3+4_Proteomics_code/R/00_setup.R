## =============================================================================
## 00_setup.R
## -----------------------------------------------------------------------------
## Loads every package the pipeline needs and defines the constants that are
## shared across all the other R/ scripts:
##   - file paths (raw data location + where figures/tables get saved)
##   - which sample(s) to drop as batch/QC outliers
##   - significance cutoffs used throughout
##   - the shared color palette
##
## Nothing in this file does any analysis -- it only sets things up.
##
## Usage:
##     source("R/00_setup.R")   # from the project root, first, before any
##                               # other R/ script
## =============================================================================

## ===================== PACKAGES =============================================
# Data wrangling
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)

  # Plotting
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
  library(gridExtra)
  library(gtable)

  # Proteomics / stats
  library(MSstats)
  library(MSstatsConvert)

  # Heatmaps
  library(ComplexHeatmap)
  library(circlize)

  # GO enrichment
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## ===================== PATHS ================================================
## EDIT THESE TWO PATHS if you move the raw data or this project folder.
## Everything else in the pipeline is relative to project_dir, so the whole
## project folder can be moved/renamed without breaking anything else.

# Where the four MSstats output CSVs live.
# EDIT ME: point this at your own copy of the KO07_53min_data folder (e.g.
# the output of KO07_CD34_Paper_Proteomics's <data.name>_data/ folder).
data_dir <- "path/to/KO07_53min_data"

# Root of this project (figures/ and output_tables/ live under here).
# `here::here()` isn't used on purpose, so this works whether you run the
# .Rmd via knit, source(), or an interactive session -- just make sure your
# R working directory is the project folder before sourcing this file.
project_dir <- getwd()

fig_dir    <- file.path(project_dir, "figures")
table_dir  <- file.path(project_dir, "output_tables")

dir.create(file.path(fig_dir, "PCA"),                     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "Correlation_Heatmaps"),     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "Temporal_Clusters"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "GO_Enrichment"),            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "Volcano_Plots"),            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "Selected_Gene_Sets"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(table_dir, "cluster_gene_lists"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(table_dir, "go_enrichment_tables"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(table_dir, "zscore_matrices"),         recursive = TRUE, showWarnings = FALSE)

## ===================== ANALYSIS CONSTANTS ===================================

# Sample(s) excluded from analysis. D9.2 is dropped for a documented lab/QC
# issue (culture contamination), not a statistical outlier call -- raw MS
# data is untouched in data_dir, this only filters it at the analysis step
# (protein_qc <- protein[!(SUBJECT %in% outlier_samples)]), applied
# consistently everywhere. Set to character(0) to keep every sample.
outlier_samples <- c("D9.2")

# Pairwise comparisons used for the temporal-trajectory figures
comparisons_of_interest <- c("D9-D5", "D14-D9", "D14-D5")

# Significance thresholds (shared across volcano plots, GO input gene lists,
# etc). ALL CAPS on purpose -- avoids R's "promise already under evaluation"
# error from a function parameter defaulting to a global of the same name.
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1     # used where a fold-change cutoff is wanted in addition to padj
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.2

# Number of temporal clusters (k) for hierarchical clustering of z-scored
# abundance profiles
n_temporal_clusters <- 3

## ===================== SHARED COLOR PALETTE =================================
## Every "pinkish" plot in the pipeline (Pearson correlation heatmaps, the
## z-scored cluster heatmap, cluster trend lines, and GO dot plots) pulls its
## colors from this single palette so they read as one consistent figure set.
## Change the three hex codes here and the whole figure set updates together.

pink_palette <- c(low = "#F3E79B", mid = "#D8B4D8", high = "#C03A87")

# A ready-to-use diverging color function for ComplexHeatmap col= arguments,
# e.g. z-scored intensity data that ranges roughly from -2 to +2.
pink_ramp_diverging <- function(range = c(-2, 0, 2)) {
  circlize::colorRamp2(range, unname(pink_palette))
}

# A ready-to-use sequential color function for correlation-style data that
# ranges from 0 to 1 (e.g. Pearson r).
pink_ramp_sequential <- function(range = c(0, 0.5, 1)) {
  circlize::colorRamp2(range, unname(pink_palette))
}

# ggplot2 equivalent (scale_color_gradientn / scale_fill_gradientn colours=)
pink_gradient_colors <- unname(pink_palette)

cat("00_setup.R loaded. Project directory:", project_dir, "\n")
