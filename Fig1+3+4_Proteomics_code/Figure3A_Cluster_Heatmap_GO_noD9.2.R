## =============================================================================
## Figure3A_Cluster_Heatmap_GO_noD9.2.R
## -----------------------------------------------------------------------------
## Figure 3A: temporal cluster z-score heatmap (blue/white/red) + mean trend
## lines + per-cluster top-5 GO Biological Process terms (simplify()-collapsed,
## Wang, cutoff 0.6; pink/yellow dot color), with the D9.2 outlier sample
## dropped from the abundance data before z-scoring/clustering. Significance
## calls (`sig_genes`) still come from `groupcomp` as MSstats computed them --
## dropping D9.2 here does not re-run groupComparison() without that sample.
## Uses k = 4 clusters (N_CLUSTERS_NO_D92), not the shared default of 3.
## Heatmap blocks are renumbered largest-to-smallest (1 = largest) for
## readability; the tree's own physical order is untouched. A console
## "heatmap cluster N = table cluster M" translation keeps the displayed
## numbers traceable to cluster_result$row_groups, since ComplexHeatmap cuts
## its own dendrogram independently once row_split uses a real tree.
##
## Usage:
##     Rscript Figure3A_Cluster_Heatmap_GO_noD9.2.R      # standalone
##     source("Figure3A_Cluster_Heatmap_GO_noD9.2.R")    # or sourced from
##                                                        # KO07CD34_Paper_Script.Rmd
## =============================================================================

## ===================== 0. Make sure R is running from the project folder ====
if (!file.exists("R/00_setup.R")) {
  # Find this script's own path (source() call, Rscript --file=, or the active
  # RStudio tab, in that order of reliability) and setwd() to its folder.
  this_file <- NULL
  for (i in rev(seq_along(sys.frames()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile) && nzchar(ofile)) { this_file <- ofile; break }
  }
  if (is.null(this_file)) {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) > 0) this_file <- sub("^--file=", "", file_arg[1])
  }
  if (is.null(this_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
  }
  if (!is.null(this_file) && nzchar(this_file)) setwd(dirname(this_file))
}
if (!file.exists("R/00_setup.R")) {
  stop(
    "\n\nCan't find R/00_setup.R from the current working directory:\n  ", getwd(),
    "\n\nOpen this script from inside AY_KO_Figure1_Proteomics_Analysis/ (the ",
    "folder this file and the R/ folder both live in), or setwd() there first.\n",
    call. = FALSE
  )
}

cat("Working directory:", getwd(), "\n")

suppressPackageStartupMessages(library(GOSemSim))

source("R/00_setup.R")
source("R/01_data_loading.R")
source("R/04_temporal_cluster_functions.R")
source("R/05_go_enrichment_functions.R")

out_pdf_dir <- file.path(fig_dir, "Temporal_Clusters")
dir.create(out_pdf_dir, recursive = TRUE, showWarnings = FALSE)

## ===================== 1. Load + standardize data, then drop the D9.2 outlier ====

raw       <- load_proteomics_data(data_dir)
protein   <- standardize_protein_dt(raw$protein)
groupcomp <- standardize_groupcomp_dt(raw$groupcomp)

# `outlier_samples` (currently c("D9.2")) is defined once in R/00_setup.R.
protein_qc <- protein[!(SUBJECT %in% outlier_samples)]
cat("Dropped as outliers:", paste(outlier_samples, collapse = ", "), "\n")
cat("Remaining D9 subjects:", paste(sort(unique(protein_qc[GROUP == "D9", SUBJECT])), collapse = ", "), "\n")

N_CLUSTERS_NO_D92 <- 4

## ===================== 2. Colors for THIS figure only =======================

HEAT_COLORS <- c(low = "#2166AC", mid = "white", high = "#B2182B")
heat_ramp   <- circlize::colorRamp2(c(-2, 0, 2), unname(HEAT_COLORS))

# Blue sequential strip for the Day 5/9/14 top annotation blocks, ending on
# the heatmap's own low (negative z-score) color.
day_groups <- c("D5", "D9", "D14")
day_block_colors <- c("D5" = "#DEEBF7", "D9" = "#6BAED6", "D14" = "#2166AC")
day_titles <- c("D5" = "Day 5", "D9" = "Day 9", "D14" = "Day 14")

TREND_LINE_COLOR <- unname(HEAT_COLORS["high"])   # matches the heatmap's high (positive) red

# GO dot color uses pink_gradient_colors from R/00_setup.R as-is.

## ===================== 3. Significant genes -> z-scored abundance matrix (D9.2 excluded) ====

sig_genes <- significant_genes(groupcomp, comparisons_of_interest,
                                padj_cutoff = PADJ_CUTOFF, lfc_cutoff = 0)
cat("Number of significant genes retained:", length(sig_genes), "\n")

abundance_mat   <- build_group_wide_matrix(protein_qc, day_groups, genes = sig_genes)
abundance_mat_z <- zscore_rows(abundance_mat)
cat("Proteins retained in abundance heatmap:", nrow(abundance_mat_z), "\n")

cluster_result <- cluster_temporal_profiles(abundance_mat_z, k = N_CLUSTERS_NO_D92)
print(cluster_result$cluster_sizes)

day_vec <- sub("_.*$", "", colnames(abundance_mat_z))

## ===================== 4. Heatmap: blue/white/red, day blocks instead of per-sample labels ====

col_split <- factor(day_vec, levels = day_groups)

top_ha <- HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = unname(day_block_colors[day_groups]), col = NA)),
  show_annotation_name = FALSE
)

ht_clusters <- Heatmap(
  abundance_mat_z,
  name = "z-score",
  col = heat_ramp,
  top_annotation = top_ha,
  column_split = col_split,
  column_title = unname(day_titles[day_groups]),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  cluster_rows = cluster_result$hclust,
  show_row_dend = TRUE,
  row_dend_side = "right",     # keeps the tree clear of the "Cluster" label on the left
  cluster_row_slices = TRUE,   # also draws the connecting tree BETWEEN clusters
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_split = N_CLUSTERS_NO_D92,   # ComplexHeatmap only accepts a single number once cluster_rows is a dendrogram
  row_gap = unit(2, "mm"),
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  show_row_names = FALSE,
  show_heatmap_legend = FALSE   # drawn separately below, stacked with the GO legend
)

# Probe-draw into a throwaway device to read back which rows landed in which
# visual split block, then map that back to cluster_result$row_groups --
# membership is guaranteed identical (same hclust, same k), the block number
# ComplexHeatmap assigns is not.
pdf(NULL)
ht_probe <- draw(ht_clusters, heatmap_legend_side = "left")
row_blocks <- row_order(ht_probe)
dev.off()

# Renumber displayed blocks largest -> smallest (1 = largest); the tree's own
# physical top-to-bottom order is untouched.
block_sizes <- vapply(row_blocks, length, integer(1))
heatmap_cluster_labels <- as.character(rank(-block_sizes, ties.method = "first"))

ht_clusters <- Heatmap(
  abundance_mat_z,
  name = "z-score",
  col = heat_ramp,
  top_annotation = top_ha,
  column_split = col_split,
  column_title = unname(day_titles[day_groups]),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  cluster_rows = cluster_result$hclust,
  show_row_dend = TRUE,
  row_dend_side = "right",
  cluster_row_slices = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_split = N_CLUSTERS_NO_D92,
  row_title = heatmap_cluster_labels,   # size-ranked labels, same physical block order as the first build
  row_gap = unit(2, "mm"),
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)

pdf(NULL)
ht_probe2 <- draw(ht_clusters, heatmap_legend_side = "left")
row_blocks_check <- row_order(ht_probe2)
dev.off()

if (!identical(vapply(row_blocks_check, length, integer(1)), block_sizes)) {
  warning(
    "\n\nThe rebuilt heatmap's block sizes/order don't match the first probe draw -- ",
    "the size-ranked numbers in row_title may not be aligned with the actual blocks. ",
    "Compare the two 'Cluster sizes' printouts below carefully before trusting this figure.\n",
    call. = FALSE
  )
  cat("First probe block sizes: ");  print(block_sizes)
  cat("Rebuilt probe block sizes: "); print(vapply(row_blocks_check, length, integer(1)))
} else {
  cat("Sanity check passed: rebuilt heatmap's block order/sizes match the first probe.\n")
}

heatmap_to_table_cluster <- vapply(seq_along(row_blocks), function(i) {
  genes_in_block <- rownames(abundance_mat_z)[row_blocks[[i]]]
  table_cls <- cluster_result$row_groups[genes_in_block]
  as.integer(names(sort(table(table_cls), decreasing = TRUE))[1])
}, integer(1))

cat("\nHeatmap cluster number (size-ranked, 1 = largest) -> table (CSV) cluster number:\n")
for (i in seq_along(row_blocks)) {
  cat("  Heatmap cluster", heatmap_cluster_labels[i], "(n =", block_sizes[i],
      ") = table cluster", heatmap_to_table_cluster[i], "\n")
}
cat("(Cross-referencing cluster_gene_lists/ or",
    "Cluster_GO_BP_full_results_matched_noD9.2.csv against a specific block",
    "on the heatmap? Use this mapping -- the two numbering schemes are not",
    "the same.)\n\n")

# cluster_row_slices = TRUE lets ComplexHeatmap reorder blocks top-to-bottom,
# so the trend-line and GO panels below use heatmap_cluster_labels /
# heatmap_to_table_cluster (not cluster_result$row_groups's own 1..k order)
# to stay aligned with whatever ComplexHeatmap actually drew.
table_to_heatmap_label <- setNames(heatmap_cluster_labels, heatmap_to_table_cluster)
cluster_levels <- paste0("Cluster ", heatmap_cluster_labels)   # heatmap's own top-to-bottom order

## ===================== 5. Trend lines (n = ... boxes), one per cluster ======

cluster_day_mean <- compute_cluster_day_means(cluster_result, day_vec, day_levels = day_groups)
cluster_day_mean$cluster <- factor(
  paste0("Cluster ", table_to_heatmap_label[sub("^Cluster ", "", as.character(cluster_day_mean$cluster))]),
  levels = cluster_levels
)
cluster_sizes_by_heatmap_label <- setNames(
  as.integer(cluster_result$cluster_sizes[as.character(heatmap_to_table_cluster)]),
  heatmap_cluster_labels
)
trend_grobs <- make_trend_grobs(cluster_day_mean, cluster_sizes_by_heatmap_label,
                                 line_color = TREND_LINE_COLOR)

## ===================== 6. GO Biological Process, top 5 terms per cluster, simplify()-collapsed ====
## Universe = every protein in the clustered heatmap (D9.2-excluded matrix).

SIMPLIFY_MEASURE <- "Wang"
SIMPLIFY_CUTOFF  <- 0.6

go.semdata <- GOSemSim::godata(
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  computeIC = FALSE   # Wang is graph-topology-based, doesn't need Information Content
)

run_go_by_cluster_simplified <- function(cluster_genes, universe) {
  go_list <- lapply(seq_along(cluster_genes), function(i) {
    ego <- run_go_enrichment(cluster_genes[[i]], universe = universe)
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)

    ego <- clusterProfiler::simplify(
      ego,
      cutoff = SIMPLIFY_CUTOFF,
      by = "p.adjust",
      select_fun = min,
      measure = SIMPLIFY_MEASURE,
      semData = go.semdata
    )

    df <- as.data.frame(ego)
    if (nrow(df) == 0) return(NULL)
    df$Cluster <- paste0("Cluster ", i)
    df
  })
  rbindlist(go_list, fill = TRUE)
}

go_df <- run_go_by_cluster_simplified(cluster_result$cluster_genes, universe = rownames(cluster_result$mat_z))

if (nrow(go_df) > 0) {
  go_df$Cluster <- paste0("Cluster ", table_to_heatmap_label[sub("^Cluster ", "", go_df$Cluster)])
  go_top   <- top_go_terms(go_df, n = 5, cluster_levels = cluster_levels)
  go_grobs <- make_go_grobs(go_top, cluster_levels)
  go_legend_gtable <- build_go_legend(go_top)   # Count (size) + -Log10(P.adjust) (pink/yellow color)
} else {
  message("No GO terms enriched in any cluster -- GO column and legend will be blank.")
  go_grobs <- lapply(cluster_levels, function(x) ggplotGrob(ggplot() + theme_void()))
  go_legend_gtable <- NULL
}

## ===================== 7. Standalone z-score legend (stacked with the GO legend on the right) ====

z_legend <- ComplexHeatmap::Legend(
  col_fun = heat_ramp,
  title = "z-score",
  at = c(-2, -1, 0, 1, 2)
)

## ===================== 8. Assemble the combined page ========================
## Four columns: heatmap | trend-line boxes | GO dot terms | legends.

assemble_figure <- function(outfile, width = 14, height = 9,
                             cluster_label_x = unit(3, "mm"),
                             legend_heights = c(0.62, 0.38)) {
  widths <- unit.c(unit(0.36, "npc"), unit(0.18, "npc"), unit(0.34, "npc"), unit(0.12, "npc"))

  pdf(outfile, width = width, height = height)
  on.exit(dev.off(), add = TRUE)

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 4, widths = widths)))

  # -- Column 1: heatmap --
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_clusters, newpage = FALSE)
  grid.text("Cluster", x = cluster_label_x, y = unit(0.5, "npc"), rot = 90,
             just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
  upViewport()

  # -- Column 2: trend-line boxes, one per cluster --
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                         layout = grid.layout(nrow = length(trend_grobs), ncol = 1)))
  for (i in seq_along(trend_grobs)) {
    grid.draw(editGrob(trend_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()

  # -- Column 3: GO dot terms, one panel per cluster --
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3,
                         layout = grid.layout(nrow = length(go_grobs), ncol = 1)))
  for (i in seq_along(go_grobs)) {
    grid.draw(editGrob(go_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()

  # -- Column 4: legends -- GO (Count above -Log10(P.adjust)) stacked above z-score --
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4,
                         layout = grid.layout(nrow = 2, ncol = 1, heights = legend_heights)))
  if (!is.null(go_legend_gtable)) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.draw(go_legend_gtable)
    upViewport()
  }
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(z_legend, x = unit(0.5, "npc"), y = unit(0.9, "npc"), just = "top")
  upViewport()
  upViewport()

  invisible(outfile)
}

outfile_pdf <- file.path(out_pdf_dir, sprintf("Figure3A_Cluster_Heatmap_GO_matched_noD9.2_k%d.pdf", N_CLUSTERS_NO_D92))
assemble_figure(outfile_pdf)
cat("Saved:", outfile_pdf, "\n")

# PNG version at the same proportions, for quick viewing/sharing outside R.
outfile_png <- sub("\\.pdf$", ".png", outfile_pdf)
png(outfile_png, width = 14, height = 9, units = "in", res = 300)
tryCatch({
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 4,
                                              widths = unit.c(unit(0.36, "npc"), unit(0.18, "npc"),
                                                               unit(0.34, "npc"), unit(0.12, "npc")))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_clusters, newpage = FALSE)
  grid.text("Cluster", x = unit(3, "mm"), y = unit(0.5, "npc"), rot = 90,
             just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                         layout = grid.layout(nrow = length(trend_grobs), ncol = 1)))
  for (i in seq_along(trend_grobs)) {
    grid.draw(editGrob(trend_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3,
                         layout = grid.layout(nrow = length(go_grobs), ncol = 1)))
  for (i in seq_along(go_grobs)) {
    grid.draw(editGrob(go_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4,
                         layout = grid.layout(nrow = 2, ncol = 1, heights = c(0.62, 0.38))))
  if (!is.null(go_legend_gtable)) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.draw(go_legend_gtable)
    upViewport()
  }
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(z_legend, x = unit(0.5, "npc"), y = unit(0.9, "npc"), just = "top")
  upViewport()
  upViewport()
}, finally = dev.off())
cat("Saved:", outfile_png, "\n")

## ===================== 9. Save the underlying tables too (gene lists + full GO results) ====

for (i in seq_along(cluster_result$cluster_genes)) {
  write.csv(
    data.frame(gene = cluster_result$cluster_genes[[i]]),
    file.path(table_dir, "cluster_gene_lists", sprintf("Cluster_%d_genes_noD9.2.csv", i)),
    row.names = FALSE
  )
}
if (nrow(go_df) > 0) {
  setorder(go_df, Cluster, p.adjust)
  go_out_path <- file.path(table_dir, "go_enrichment_tables", "Cluster_GO_BP_full_results_matched_noD9.2.csv")
  fwrite(go_df, go_out_path)
  cat("Full GO term table (all terms, ranked by p.adjust within cluster) saved to:\n ", go_out_path, "\n")
}

## ===================== 10. Publication-ready z-score table ==================
## One row per protein, one column per sample, plus both cluster numbering
## schemes (Table_Cluster, Heatmap_Cluster).

gene_to_heatmap_cluster <- setNames(
  rep(heatmap_cluster_labels, times = vapply(row_blocks, length, integer(1))),
  rownames(abundance_mat_z)[unlist(row_blocks)]
)

zscore_table <- data.frame(
  Protein = rownames(abundance_mat_z),
  abundance_mat_z,
  check.names = FALSE
)
zscore_table$Table_Cluster   <- cluster_result$row_groups[zscore_table$Protein]
zscore_table$Heatmap_Cluster <- gene_to_heatmap_cluster[zscore_table$Protein]
zscore_table <- zscore_table[order(as.integer(zscore_table$Heatmap_Cluster), zscore_table$Protein), ]

zscore_out_path <- file.path(table_dir, "zscore_matrices", "Figure3A_zscore_matrix_noD9.2.csv")
write.csv(zscore_table, zscore_out_path, row.names = FALSE)
cat("Per-protein z-score matrix (one column per sample, plus cluster assignment) saved to:\n ", zscore_out_path, "\n")
