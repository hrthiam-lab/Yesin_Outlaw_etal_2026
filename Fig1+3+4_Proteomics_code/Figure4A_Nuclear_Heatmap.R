## =============================================================================
## Figure4A_Nuclear_Heatmap.R
## -----------------------------------------------------------------------------
## Figure 4A: HPA nuclear gene-set heatmap (gene_lists/nuclear2.tsv), styled
## like Figure 3A (blue/white/red heat, Day 5/9/14 blocks). Sparse row labels:
## top 5 genes per cluster by score = -log10(best adj.pvalue) * best |log2FC|
## (best across the 3 pairwise comparisons); LBR/LMNA/LMNB1/LMNB2 are always
## included and always labeled (bold, red), as extras if they miss the top 5.
## Each of those 4 also gets its own trajectory panel (real per-gene z-score +
## jittered replicate points, not a cluster average). D9.2 excluded (see
## R/00_setup.R); k = 4 clusters (N_CLUSTERS_FIG4A) instead of the shared
## default of 3, so LMNA can split from LMNB1/LMNB2.
##
## Usage:
##     Rscript Figure4A_Nuclear_Heatmap.R      # standalone
##     source("Figure4A_Nuclear_Heatmap.R")    # or sourced from
##                                              # KO07CD34_Paper_Script.Rmd
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
    "folder this file, the R/ folder, and gene_lists/ all live in), or setwd() there first.\n",
    call. = FALSE
  )
}
if (!file.exists("gene_lists/nuclear2.tsv")) {
  stop("\n\nCan't find gene_lists/nuclear2.tsv (the Human Protein Atlas nuclear list) from:\n  ",
       getwd(), "\n", call. = FALSE)
}
cat("Working directory:", getwd(), "\n")

source("R/00_setup.R")
source("R/01_data_loading.R")
source("R/04_temporal_cluster_functions.R")

out_fig_dir <- file.path(fig_dir, "Selected_Gene_Sets")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
out_tbl_dir <- file.path(table_dir, "nuclear_envelope_gene_sets")
dir.create(out_tbl_dir, recursive = TRUE, showWarnings = FALSE)

## ===================== 1. Load + standardize data, drop D9.2 ================

raw       <- load_proteomics_data(data_dir)
protein   <- standardize_protein_dt(raw$protein)
groupcomp <- standardize_groupcomp_dt(raw$groupcomp)
protein_qc <- protein[!(SUBJECT %in% outlier_samples)]
cat("Dropped as outliers:", paste(outlier_samples, collapse = ", "), "\n")

## ===================== 2. Colors (match Figure3A_Cluster_Heatmap_GO_noD9.2.R exactly) ====

HEAT_COLORS <- c(low = "#2166AC", mid = "white", high = "#B2182B")
heat_ramp   <- circlize::colorRamp2(c(-2, 0, 2), unname(HEAT_COLORS))

day_groups <- c("D5", "D9", "D14")
day_block_colors <- c("D5" = "#DEEBF7", "D9" = "#6BAED6", "D14" = "#2166AC")
day_titles <- c("D5" = "Day 5", "D9" = "Day 9", "D14" = "Day 14")

## ===================== MANUAL LABEL OVERRIDES (OPTIONAL) ====================
## Hand-pick a cluster's 5 labels instead of the automatic top-5-by-score,
## e.g. MANUAL_LABEL_OVERRIDES <- list(`1` = c("RAB7A", "VPS35")). Unlisted
## clusters still use the automatic ranking; genes not actually in that
## cluster are skipped with a message.
MANUAL_LABEL_OVERRIDES <- list()

TRAJ_MEAN_COLOR  <- unname(HEAT_COLORS["high"])   # "#B2182B" -- day-mean line/point
TRAJ_POINT_COLOR <- unname(HEAT_COLORS["low"])    # "#2166AC" -- individual replicate points

LAMINA_GENES <- c("LBR", "LMNA", "LMNB1", "LMNB2")

# k for this figure only -- local override of the shared n_temporal_clusters
# (= 3, R/00_setup.R); doesn't affect Figure 3A or anything else that sources
# R/00_setup.R directly.
N_CLUSTERS_FIG4A <- 4

## ===================== 3. Gene list + significance -> candidate genes for the heatmap ====

hpa_dt    <- fread("gene_lists/nuclear2.tsv")
genes_hpa <- unique(toupper(trimws(as.character(hpa_dt$Gene))))
genes_hpa <- genes_hpa[!is.na(genes_hpa) & genes_hpa != ""]
cat("Unique gene symbols from nuclear2.tsv (Human Protein Atlas):", length(genes_hpa), "\n")

sig_genes <- significant_genes(groupcomp, comparisons_of_interest,
                                padj_cutoff = PADJ_CUTOFF, lfc_cutoff = 0)
cat("Genes significant in >=1 pairwise comparison, whole dataset:", length(sig_genes), "\n")

hpa_significant <- intersect(genes_hpa, sig_genes)
cat("HPA nuclear genes significant in >=1 comparison:", length(hpa_significant), "\n")

candidate_genes <- unique(c(hpa_significant, LAMINA_GENES))
cat("Candidate genes for heatmap (HPA-significant + LBR/LMNA/LMNB1/LMNB2 forced in):",
    length(candidate_genes), "\n")

## ===================== 4. Build matrix (D9.2 excluded), z-score, complete-case filter ====

abundance_mat   <- build_group_wide_matrix(protein_qc, day_groups, genes = candidate_genes)
abundance_mat_z <- zscore_rows(abundance_mat)
dropped_missing <- setdiff(rownames(abundance_mat), rownames(abundance_mat_z))
cat("Retained after requiring complete data in every replicate:", nrow(abundance_mat_z), "\n")
if (length(dropped_missing) > 0) {
  cat("Dropped for missing values in >=1 replicate:", paste(dropped_missing, collapse = ", "), "\n")
}

lamina_present <- intersect(LAMINA_GENES, rownames(abundance_mat_z))
lamina_missing <- setdiff(LAMINA_GENES, rownames(abundance_mat_z))
cat("LBR/LMNA/LMNB1/LMNB2 on the heatmap:", paste(lamina_present, collapse = ", "), "\n")
if (length(lamina_missing) > 0) {
  message("LBR/LMNA/LMNB1/LMNB2 NOT on the heatmap (not detected, or missing a replicate value): ",
          paste(lamina_missing, collapse = ", "))
}

if (nrow(abundance_mat_z) < 2) {
  stop("Fewer than 2 genes survived the significance + complete-case filters -- nothing to cluster.", call. = FALSE)
}

## ===================== 5. Cluster ===========================================

k_use <- min(N_CLUSTERS_FIG4A, nrow(abundance_mat_z))
cat("Using k =", k_use, "clusters for Figure 4A (N_CLUSTERS_FIG4A =", N_CLUSTERS_FIG4A, ")\n")
cluster_result <- cluster_temporal_profiles(abundance_mat_z, k = k_use)
print(cluster_result$cluster_sizes)

day_vec <- sub("_.*$", "", colnames(abundance_mat_z))

## ComplexHeatmap only accepts a single integer for row_split once cluster_rows
## is a dendrogram, so it re-cuts the same hclust object itself below (row_split
## = k_use) -- membership matches cluster_result$row_groups exactly, but the
## printed block number may not; section 7 prints a translation table.

## ===================== 6. Rank genes per cluster: ONE combined score, top 5 per cluster ====
## score = -log10(best_adj_pvalue) * best_abs_log2fc, ranked among every gene
## on the heatmap (LAMINA_GENES included, competing on equal footing).

rank_pool <- rownames(cluster_result$mat_z)

rank_dt <- groupcomp[
  Label %in% comparisons_of_interest & gene.name %in% rank_pool &
    !is.na(adj.pvalue) & !is.na(log2FC),
  .(best_adj_pvalue = min(adj.pvalue, na.rm = TRUE),
    best_abs_log2fc = max(abs(log2FC), na.rm = TRUE)),
  by = gene.name
]
rank_dt[, Cluster := cluster_result$row_groups[gene.name]]
rank_dt[, score := -log10(pmax(best_adj_pvalue, .Machine$double.xmin)) * best_abs_log2fc]

cluster_ids <- sort(unique(cluster_result$row_groups))

# Per cluster: MANUAL_LABEL_OVERRIDES if set, else automatic top 5 by score --
# plus any LAMINA_GENES in that cluster not already picked, added as extras.
label_by_cluster <- setNames(lapply(cluster_ids, function(cl) {
  cluster_genes_this <- names(cluster_result$row_groups)[cluster_result$row_groups == cl]
  override <- MANUAL_LABEL_OVERRIDES[[as.character(cl)]]

  if (!is.null(override)) {
    picked  <- intersect(override, cluster_genes_this)
    skipped <- setdiff(override, cluster_genes_this)
    if (length(skipped) > 0) {
      message("Cluster ", cl, ": manual label(s) not in this cluster, skipped: ", paste(skipped, collapse = ", "))
    }
    top5   <- picked
    source <- "manual"
  } else {
    top5   <- head(rank_dt[Cluster == cl][order(-score), gene.name], 5)
    source <- "auto_score"
  }

  lamina_here  <- intersect(LAMINA_GENES, cluster_genes_this)
  lamina_extra <- setdiff(lamina_here, top5)

  list(top5 = top5, lamina_extra = lamina_extra, all = unique(c(top5, lamina_extra)), source = source)
}), cluster_ids)

for (cl in cluster_ids) {
  info <- label_by_cluster[[as.character(cl)]]
  cat("Cluster", cl, "- ", if (info$source == "manual") "manual picks" else "top 5", ":",
      paste(info$top5, collapse = ", "))
  if (length(info$lamina_extra) > 0) {
    cat(" | + extra (lamina gene, not already in that list):", paste(info$lamina_extra, collapse = ", "))
  }
  cat("\n")
}

label_genes_all <- unique(unlist(lapply(label_by_cluster, `[[`, "all")))

# Why each gene was labeled, for the saved table (top5 / manual / lamina_extra).
label_lookup <- rbindlist(lapply(cluster_ids, function(cl) {
  info <- label_by_cluster[[as.character(cl)]]
  top5_reason <- if (info$source == "manual") "manual" else "top5"
  rbind(
    if (length(info$top5) > 0) data.table(gene.name = info$top5, label_reason = top5_reason) else NULL,
    if (length(info$lamina_extra) > 0) data.table(gene.name = info$lamina_extra, label_reason = "lamina_extra") else NULL
  )
}), fill = TRUE)
rank_dt <- merge(rank_dt, label_lookup, by = "gene.name", all.x = TRUE)
rank_dt[is.na(label_reason), label_reason := "not_labeled"]

ranking_csv_path <- file.path(out_tbl_dir, "Figure4A_gene_ranking_by_cluster.csv")
fwrite(rank_dt[order(Cluster, -score)], ranking_csv_path)
cat("Saved gene ranking table to:", ranking_csv_path, "\n")

## ===================== 7. Heatmap: blue/white/red, day blocks, sparse labels ====

col_split <- factor(day_vec, levels = day_groups)
top_ha <- HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = unname(day_block_colors[day_groups]), col = NA)),
  show_annotation_name = FALSE
)

# anno_mark()'s `at` indexes the ORIGINAL row order (abundance_mat_z);
# ComplexHeatmap remaps it internally after clustering.
label_idx <- match(label_genes_all, rownames(abundance_mat_z))
is_lamina <- label_genes_all %in% LAMINA_GENES

row_anno <- rowAnnotation(
  link = anno_mark(
    at = label_idx,
    labels = label_genes_all,
    labels_gp = gpar(
      fontsize = 8,
      fontface = ifelse(is_lamina, "bold", "plain"),
      col = ifelse(is_lamina, TRAJ_MEAN_COLOR, "black")
    )
  )
)

# cluster_rows = the same hclust object cluster_temporal_profiles() computed
# (no re-clustering). Passing the original abundance_mat_z (not the
# pre-sorted cluster_result$mat_z) lets ComplexHeatmap reorder rows from the
# dendrogram's own leaf order; row_split = k_use is the only form it accepts
# once cluster_rows is a dendrogram (see the numbering note above section 6).
ht_nuclear <- Heatmap(
  abundance_mat_z,
  name = "z-score",
  col = heat_ramp,
  top_annotation = top_ha,
  right_annotation = row_anno,
  column_split = col_split,
  column_title = unname(day_titles[day_groups]),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  cluster_rows = cluster_result$hclust,
  show_row_dend = TRUE,
  cluster_row_slices = TRUE,   # also draws the connecting tree BETWEEN clusters
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_split = k_use,
  row_gap = unit(2, "mm"),
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  show_row_names = FALSE
)

# Draw once to a throwaway device just to read back row_order() -- maps
# "heatmap cluster N" to "table (CSV) cluster M" since the two numbering
# schemes aren't guaranteed to agree (same hclust object, same k, so gene
# MEMBERSHIP is identical either way).
pdf(NULL)
ht_probe <- draw(ht_nuclear, heatmap_legend_side = "left")
row_blocks <- row_order(ht_probe)
dev.off()

# Use ComplexHeatmap's own slice names (the printed row title) where
# available, rather than assuming list position lines up with the display.
heatmap_cluster_labels <- names(row_blocks)
if (is.null(heatmap_cluster_labels)) heatmap_cluster_labels <- as.character(seq_along(row_blocks))

heatmap_to_table_cluster <- vapply(seq_along(row_blocks), function(i) {
  genes_in_block <- rownames(abundance_mat_z)[row_blocks[[i]]]
  table_cls <- cluster_result$row_groups[genes_in_block]
  as.integer(names(sort(table(table_cls), decreasing = TRUE))[1])
}, integer(1))

cat("\nHeatmap cluster number -> table (CSV) cluster number:\n")
for (i in seq_along(row_blocks)) {
  cat("  Heatmap cluster", heatmap_cluster_labels[i], "= table cluster", heatmap_to_table_cluster[i], "\n")
}
cat("(If you're cross-referencing Figure4A_gene_clusters.csv or",
    "Figure4A_gene_ranking_by_cluster.csv against a specific block on the",
    "heatmap, use this mapping -- the two numbering schemes are not the same.)\n\n")

## ===================== 8. Trajectory panels for LBR/LMNA/LMNB1/LMNB2 ========
## Real per-gene z-score (same matrix the heatmap uses), individual replicate
## points jittered so all 3 are visible, connected through each day's mean --
## NOT a cluster-average trend line.

plot_gene_trajectory <- function(gene, mat_z, day_levels) {
  row_vals <- mat_z[gene, ]
  df <- data.frame(sample = names(row_vals), z = as.numeric(row_vals))
  df$Day <- factor(sub("_.*$", "", df$sample), levels = day_levels)

  day_mean <- aggregate(z ~ Day, df, mean)

  ggplot() +
    geom_jitter(data = df, aes(x = Day, y = z), width = 0.08, height = 0,
                color = TRAJ_POINT_COLOR, size = 1.6, alpha = 0.75) +
    geom_line(data = day_mean, aes(x = Day, y = z, group = 1),
              color = TRAJ_MEAN_COLOR, linewidth = 1) +
    geom_point(data = day_mean, aes(x = Day, y = z),
               color = TRAJ_MEAN_COLOR, size = 2.2) +
    theme_bw() +
    labs(title = gene, x = NULL, y = "z-score") +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

blank_gene_panel <- function(gene) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste0(gene, "\nnot on heatmap\n(missing data)"),
              size = 3, hjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

gene_trajectory_grobs <- lapply(LAMINA_GENES, function(g) {
  p <- if (g %in% lamina_present) {
    plot_gene_trajectory(g, cluster_result$mat_z, day_groups)
  } else {
    blank_gene_panel(g)
  }
  ggplotGrob(p)
})

## ===================== 9. Assemble the combined page: heatmap | 4 gene trajectory panels ====

assemble_figure <- function(outfile, width = 12, height = 9) {
  widths <- unit.c(unit(0.72, "npc"), unit(0.28, "npc"))

  pdf(outfile, width = width, height = height)
  on.exit(dev.off(), add = TRUE)

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = widths)))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_nuclear, newpage = FALSE, heatmap_legend_side = "left")
  upViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                         layout = grid.layout(nrow = length(gene_trajectory_grobs), ncol = 1)))
  for (i in seq_along(gene_trajectory_grobs)) {
    grid.draw(editGrob(gene_trajectory_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()

  upViewport()
  invisible(outfile)
}

outfile_pdf <- file.path(out_fig_dir, "Figure4A_Nuclear_Heatmap_Labeled.pdf")
assemble_figure(outfile_pdf)
cat("Saved:", outfile_pdf, "\n")

outfile_png <- sub("\\.pdf$", ".png", outfile_pdf)
png(outfile_png, width = 12, height = 9, units = "in", res = 300)
tryCatch({
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2,
                                              widths = unit.c(unit(0.72, "npc"), unit(0.28, "npc")))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_nuclear, newpage = FALSE, heatmap_legend_side = "left")
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                         layout = grid.layout(nrow = length(gene_trajectory_grobs), ncol = 1)))
  for (i in seq_along(gene_trajectory_grobs)) {
    grid.draw(editGrob(gene_trajectory_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()
  upViewport()
}, finally = dev.off())
cat("Saved:", outfile_png, "\n")

## ===================== 10. Save tables ======================================

cluster_df <- data.frame(
  gene = names(cluster_result$row_groups),
  cluster = unname(cluster_result$row_groups),
  labeled = names(cluster_result$row_groups) %in% label_genes_all,
  is_lamina_gene = names(cluster_result$row_groups) %in% LAMINA_GENES
)
write.csv(cluster_df, file.path(out_tbl_dir, "Figure4A_gene_clusters.csv"), row.names = FALSE)
cat("Saved gene cluster table to:", file.path(out_tbl_dir, "Figure4A_gene_clusters.csv"), "\n")

## ===================== PUBLICATION-READY Z-SCORE TABLE ======================
## One row per protein, one column per sample, plus both cluster numbering
## schemes (Table_Cluster, Heatmap_Cluster) and Is_Lamina_Gene.

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
zscore_table$Is_Lamina_Gene  <- zscore_table$Protein %in% LAMINA_GENES
zscore_table <- zscore_table[order(zscore_table$Heatmap_Cluster, zscore_table$Protein), ]

zscore_out_path <- file.path(table_dir, "zscore_matrices", "Figure4A_zscore_matrix.csv")
write.csv(zscore_table, zscore_out_path, row.names = FALSE)
cat("Per-protein z-score matrix (one column per sample, plus cluster assignment) saved to:\n ", zscore_out_path, "\n")
