## =============================================================================
## 04_temporal_cluster_functions.R
## -----------------------------------------------------------------------------
## Builds the "Figure 3A"-style temporal cluster figure: a z-scored abundance
## heatmap split into k hierarchical clusters, a small line plot of each
## cluster's mean trend across days, and (optionally, see 05_) a GO panel next
## to it. Consolidates the big single chunk from the original script plus the
## smaller repeat chunks that reused pieces of it, and fixes a bug where a
## later chunk referenced `row_clusters`, an object that was never created
## (only `row_groups` / `row_groups_renamed` existed) -- every function below
## consistently uses the renamed cluster labels.
##
## The z-scored heatmap and the trend-line color now both pull from
## `pink_palette` (00_setup.R) instead of the original blue/red diverging
## scale and plain red trend lines, per your request to match the Pearson
## figure's pink palette.
##
## Usage:
##     source("R/00_setup.R"); source("R/01_data_loading.R")
##     source("R/04_temporal_cluster_functions.R")
## =============================================================================

## ===================== 1. Build a GROUP x SUBJECT abundance matrix for chosen genes ====

# `groups`: character vector giving both which GROUPs to keep AND their
# left-to-right column order (e.g. c("D5", "D9", "D14")).
# `genes`: optional character vector to restrict rows to (e.g. significant
# genes only). NULL keeps every gene.
build_group_wide_matrix <- function(protein_dt, groups, genes = NULL) {
  dt <- protein_dt[GROUP %in% groups & !is.na(LogIntensities) & !is.na(gene.name)]
  if (!is.null(genes)) dt <- dt[gene.name %in% genes]

  dt <- copy(dt)
  dt[, sample_name := paste(GROUP, SUBJECT, sep = "_")]

  wide <- dcast(dt, gene.name ~ sample_name, value.var = "LogIntensities", fun.aggregate = mean)
  mat <- as.data.frame(wide)
  rownames(mat) <- mat$gene.name
  mat$gene.name <- NULL
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  # order columns to match `groups`, samples sorted within each group
  ordered_cols <- unlist(lapply(groups, function(g) {
    cols <- colnames(mat)[grepl(paste0("^", g, "_"), colnames(mat))]
    sort(cols)
  }))
  mat[, ordered_cols, drop = FALSE]
}

# Genes significant in any of `labels` (adj.pvalue < padj_cutoff, |log2FC| > lfc_cutoff)
significant_genes <- function(groupcomp_dt, labels, padj_cutoff = 0.05, lfc_cutoff = 0) {
  unique(groupcomp_dt[
    Label %in% labels & !is.na(gene.name) & !is.na(adj.pvalue) & !is.na(log2FC) &
      adj.pvalue < padj_cutoff & abs(log2FC) > lfc_cutoff,
    gene.name
  ])
}

## ===================== 2. Row z-score + hierarchical clustering =============

zscore_rows <- function(mat) {
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  z <- t(scale(t(mat)))
  z[apply(z, 1, function(x) all(is.finite(x))), , drop = FALSE]
}

# Clusters rows of a z-scored matrix into k groups, reorders rows by cluster,
# and relabels clusters 1..k in the order they appear top-to-bottom (matching
# the original script's "cluster_map" relabeling step).
cluster_temporal_profiles <- function(mat_z, k = 3, method = "complete") {
  row_hc <- hclust(dist(mat_z), method = method)
  row_groups <- cutree(row_hc, k = k)

  ord <- order(row_groups)
  mat_z_ordered <- mat_z[ord, , drop = FALSE]
  row_groups <- row_groups[rownames(mat_z_ordered)]

  cluster_levels <- sort(unique(row_groups))
  cluster_map <- setNames(seq_along(cluster_levels), cluster_levels)
  row_groups_renamed <- unname(cluster_map[as.character(row_groups)])
  names(row_groups_renamed) <- names(row_groups)

  list(
    mat_z = mat_z_ordered,
    row_groups = row_groups_renamed,
    cluster_sizes = table(row_groups_renamed),
    cluster_genes = split(names(row_groups_renamed), row_groups_renamed),
    hclust = row_hc
  )
}

## ===================== 3. Cluster heatmap (pink diverging palette) ==========

# `day_vec` should be a character vector, one entry per column of mat_z,
# giving that column's day/group label (e.g. "D5"). `day_colors` is a named
# vector mapping each day to a color for the top annotation strip.
plot_cluster_heatmap <- function(cluster_result, day_vec, day_colors,
                                  title = "Proteomics temporal clusters (significant proteins only)") {
  top_ha <- HeatmapAnnotation(
    Day = day_vec,
    col = list(Day = day_colors),
    annotation_name_side = "left"
  )

  Heatmap(
    cluster_result$mat_z,
    name = "z-score",
    top_annotation = top_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(cluster_result$row_groups, levels = sort(unique(cluster_result$row_groups))),
    row_gap = unit(2, "mm"),
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_title = title,
    col = pink_ramp_diverging()
  )
}

## ===================== 4. Cluster mean trend lines ==========================

# Returns a tidy data.frame of mean z-score per cluster per day, ready for
# ggplot faceting.
compute_cluster_day_means <- function(cluster_result, day_vec, day_levels) {
  mat_z <- cluster_result$mat_z
  row_groups <- cluster_result$row_groups

  cluster_mean_df <- lapply(sort(unique(row_groups)), function(cl) {
    genes <- names(row_groups)[row_groups == cl]
    submat <- mat_z[genes, , drop = FALSE]
    data.frame(
      sample = colnames(submat),
      mean_z = colMeans(submat, na.rm = TRUE),
      cluster = paste0("Cluster ", cl),
      Day = day_vec
    )
  }) %>% bind_rows()

  out <- cluster_mean_df %>%
    group_by(cluster, Day) %>%
    summarise(mean_z = mean(mean_z), .groups = "drop")

  out$Day <- factor(out$Day, levels = day_levels)
  out$cluster <- factor(out$cluster, levels = paste0("Cluster ", sort(unique(row_groups))))
  out
}

# One combined facet plot, all clusters stacked vertically.
plot_cluster_trends <- function(cluster_day_mean, line_color = pink_palette[["high"]]) {
  ggplot(cluster_day_mean, aes(x = Day, y = mean_z, group = 1)) +
    geom_line(color = line_color, linewidth = 1) +
    geom_point(color = line_color, size = 2) +
    facet_wrap(~cluster, ncol = 1, scales = "fixed") +
    theme_bw() +
    labs(title = "Cluster mean patterns", x = NULL, y = "Mean z-score") +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
}

# Small per-cluster grobs, used when assembling the combined multi-panel PDF
# page (heatmap + trends + GO side by side) in make_combined_cluster_figure().
make_trend_grobs <- function(cluster_day_mean, cluster_sizes, line_color = pink_palette[["high"]]) {
  clusters <- levels(cluster_day_mean$cluster)
  lapply(clusters, function(cl_name) {
    df <- cluster_day_mean %>% filter(cluster == cl_name)
    cl_num <- sub("Cluster ", "", cl_name)
    n_cl <- cluster_sizes[[cl_num]]
    p <- ggplot(df, aes(x = Day, y = mean_z, group = 1)) +
      geom_line(color = line_color, linewidth = 0.9) +
      geom_point(color = line_color, size = 1.8) +
      theme_bw() +
      labs(title = cl_name, subtitle = paste0("n = ", n_cl), x = NULL, y = NULL) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7, hjust = 0.5),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2)
      )
    ggplotGrob(p)
  })
}

## ===================== 5. Combined heatmap + trends + GO page (matches original Figure A) ====

# `go_grobs` and `go_legend` come from make_go_grobs()/build_go_legend() in
# 05_go_enrichment_functions.R. Pass go_grobs = NULL to skip the GO column.
make_combined_cluster_figure <- function(outfile, ht, trend_grobs, go_grobs = NULL, go_legend = NULL,
                                          title = "Proteomics temporal clusters (significant proteins only)",
                                          width = 14, height = 10) {
  ncol_layout <- if (is.null(go_grobs)) 2 else if (is.null(go_legend)) 3 else 4
  widths <- switch(
    as.character(ncol_layout),
    "2" = unit.c(unit(0.6, "npc"), unit(0.4, "npc")),
    "3" = unit.c(unit(0.4, "npc"), unit(0.2, "npc"), unit(0.4, "npc")),
    "4" = unit.c(unit(0.34, "npc"), unit(0.18, "npc"), unit(0.34, "npc"), unit(0.14, "npc"))
  )

  pdf(outfile, width = width, height = height)
  on.exit(dev.off(), add = TRUE)

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2, ncol = ncol_layout,
    heights = unit.c(unit(0.6, "in"), unit(1, "npc") - unit(0.6, "in")),
    widths = widths
  )))

  grid.text(
    title, x = 0.02, y = 0.5, just = c("left", "center"),
    gp = gpar(fontsize = 18, fontface = "bold"),
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol_layout)
  )

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(ht, newpage = FALSE, heatmap_legend_side = "left", annotation_legend_side = "left")
  upViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2,
                         layout = grid.layout(nrow = length(trend_grobs), ncol = 1)))
  for (i in seq_along(trend_grobs)) {
    grid.draw(editGrob(trend_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
  }
  upViewport()

  if (!is.null(go_grobs)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3,
                           layout = grid.layout(nrow = length(go_grobs), ncol = 1)))
    for (i in seq_along(go_grobs)) {
      grid.draw(editGrob(go_grobs[[i]], vp = viewport(layout.pos.row = i, layout.pos.col = 1)))
    }
    upViewport()
  }

  if (!is.null(go_legend)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
    grid.draw(go_legend)
    upViewport()
  }

  invisible(outfile)
}
