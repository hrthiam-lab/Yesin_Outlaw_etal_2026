## =============================================================================
## 02_pca_functions.R
## -----------------------------------------------------------------------------
## PCA on the sample x protein intensity matrix. Replaces the two nearly
## identical PCA chunks in the original script (all samples, then again after
## dropping the D9.2 outlier) with one function called twice.
##
## Usage:
##     source("R/00_setup.R"); source("R/01_data_loading.R")
##     source("R/02_pca_functions.R")
## =============================================================================

# Runs PCA on a complete-case Protein x Subject matrix and returns everything
# the plotting function needs. Sample names are expected in "Condition.Batch"
# form (e.g. "D9.2"), matching how MSstats names SUBJECT columns.
run_pca <- function(intensity_matrix) {
  mat <- complete_rows_only(intensity_matrix)
  pca_out <- prcomp(t(mat))

  pca_dt <- as.data.table(pca_out$x, keep.rownames = TRUE)
  pca_dt[, Condition := tstrsplit(rn, "\\.")[[1]]]
  pca_dt[, Batch := sapply(strsplit(rn, "\\."), function(x) tail(x, 1))]

  percent_var <- round(100 * (pca_out$sdev^2) / sum(pca_out$sdev^2), 1)

  list(
    pca_out     = pca_out,
    pca_dt      = pca_dt,
    percent_var = percent_var,
    n_proteins  = nrow(mat)
  )
}

# Builds the PC1 vs PC2 scatter plot with sample labels, colored by Condition.
plot_pca <- function(pca_result, title = NULL) {
  if (is.null(title)) {
    title <- sprintf("PCA using %d Proteins (log intensity)", pca_result$n_proteins)
  }

  ggplot(pca_result$pca_dt, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(alpha = 1.0, size = 4) +
    ggrepel::geom_text_repel(aes(label = rn), show.legend = FALSE, size = 3) +
    theme_bw() +
    xlab(sprintf("PC1, %.1f%%", pca_result$percent_var[1])) +
    ylab(sprintf("PC2, %.1f%%", pca_result$percent_var[2])) +
    ggtitle(title)
}
