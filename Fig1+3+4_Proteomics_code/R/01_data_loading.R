## =============================================================================
## 01_data_loading.R
## -----------------------------------------------------------------------------
## Loads the four MSstats output CSVs and builds the two data structures every
## downstream analysis is built from:
##   - protein     : long-format per-sample protein abundance (ProteinLevelData)
##   - groupcomp   : long-format pairwise group comparisons (GroupComparisonsData)
##   - intensity_matrix : wide Protein x Subject matrix of log intensities
##
## Reads the four CSVs from disk (no MSstats re-run needed) rather than
## depending on in-memory dataProcess()/groupComparison() objects.
##
## Usage:
##     source("R/00_setup.R")
##     source("R/01_data_loading.R")
## =============================================================================

## ===================== LOAD RAW TABLES ======================================

load_proteomics_data <- function(data_dir) {
  list(
    protein   = fread(file.path(data_dir, "20260325_ProteinLevelData.csv")),
    groupcomp = fread(file.path(data_dir, "20260325_GroupComparisonsData.csv")),
    feature   = fread(file.path(data_dir, "20260325_FeatureLevelData.csv")),
    cleaned   = fread(file.path(data_dir, "20260325_CleanedPreprocessedData.csv"))
  )
}

## ===================== GENE-NAME / COLUMN STANDARDIZATION ===================

# Upper-cases and trims gene.name, and makes sure GROUP/SUBJECT/Label are
# plain character columns (data.table sometimes reads these as factors,
# which breaks %in% comparisons downstream). Modifies `dt` in place and
# returns it invisibly, matching data.table convention.
standardize_protein_dt <- function(dt) {
  dt <- as.data.table(dt)
  if ("gene.name" %in% names(dt)) dt[, gene.name := toupper(trimws(gene.name))]
  if ("GROUP"     %in% names(dt)) dt[, GROUP     := as.character(GROUP)]
  if ("SUBJECT"   %in% names(dt)) dt[, SUBJECT   := as.character(SUBJECT)]
  dt[]
}

standardize_groupcomp_dt <- function(dt) {
  dt <- as.data.table(dt)
  if ("gene.name" %in% names(dt)) dt[, gene.name := toupper(trimws(gene.name))]
  if ("Label"     %in% names(dt)) dt[, Label     := as.character(Label)]
  dt[]
}

## ===================== WIDE INTENSITY MATRIX ================================

# Builds a Protein x SUBJECT matrix of LogIntensities from the long-format
# protein table (rows = proteins, columns = individual runs/subjects).
build_intensity_matrix <- function(protein_dt) {
  wide <- dcast(
    as.data.table(protein_dt),
    Protein ~ SUBJECT,
    value.var = "LogIntensities"
  )
  as.matrix(wide, rownames = "Protein")
}

# Drops named columns (samples) from a Protein x Subject matrix. Use this for
# the centralized outlier_samples list defined in 00_setup.R so every figure
# in the pipeline treats outliers the same way.
drop_samples <- function(mat, samples_to_drop) {
  if (length(samples_to_drop) == 0) return(mat)
  keep <- !(colnames(mat) %in% samples_to_drop)
  mat[, keep, drop = FALSE]
}

# Keeps only rows with no missing values -- required before PCA/correlation/
# clustering, all of which need complete cases.
complete_rows_only <- function(mat) {
  mat[complete.cases(mat), , drop = FALSE]
}
