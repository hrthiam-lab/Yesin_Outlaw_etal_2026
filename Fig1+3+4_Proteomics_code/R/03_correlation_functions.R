## =============================================================================
## 03_correlation_functions.R
## -----------------------------------------------------------------------------
## Pearson correlation heatmaps. Consolidates FOUR near-duplicate blocks from
## the original script into reusable functions:
##   1. sample-level correlation heatmap (run x run)
##   2. condition-level correlation, manually ordered ggplot tile heatmap
##   3. condition-level correlation, hierarchically-ordered ggplot tile heatmap
##   4. condition-level correlation, ComplexHeatmap with dendrograms + cell labels
## All four now share the same pink palette (pink_ramp_sequential /
## pink_gradient_colors from 00_setup.R) -- this is the "Pearson" color scheme
## the cluster heatmap, trend lines, and GO plots are matched to.
##
## Usage:
##     source("R/00_setup.R"); source("R/01_data_loading.R")
##     source("R/03_correlation_functions.R")
## =============================================================================

## ===================== CORRELATION MATRICES =================================

# Run-level (sample-level) Pearson correlation on a Protein x Subject matrix.
compute_sample_correlation <- function(intensity_matrix, method = "pearson") {
  cor(intensity_matrix, use = "pairwise.complete.obs", method = method)
}

# Condition-level Pearson correlation: averages replicates within each GROUP
# first, then correlates the condition means to each other. Only proteins
# with a finite value in every requested condition are used.
compute_condition_correlation <- function(protein_dt, conditions, method = "pearson") {
  cond_mean_dt <- protein_dt[
    GROUP %in% conditions & !is.na(LogIntensities) & !is.na(gene.name),
    .(mean_logint = mean(LogIntensities, na.rm = TRUE)),
    by = .(gene.name, GROUP)
  ]

  cond_wide <- dcast(cond_mean_dt, gene.name ~ GROUP, value.var = "mean_logint")
  cond_mat  <- as.data.frame(cond_wide)
  rownames(cond_mat) <- cond_mat$gene.name
  cond_mat$gene.name <- NULL
  cond_mat <- as.matrix(cond_mat)
  mode(cond_mat) <- "numeric"
  cond_mat <- complete_rows_only(cond_mat)

  cor(cond_mat, method = method)
}

## ===================== SAVING A COMPLEXHEATMAP OBJECT TO DISK ===============

save_heatmap <- function(ht, path, width, height, ...) {
  ext <- tools::file_ext(path)
  if (ext == "pdf") {
    pdf(path, width = width, height = height)
  } else if (ext == "png") {
    png(path, width = width, height = height, units = "in", res = 300)
  } else {
    stop("save_heatmap: path must end in .pdf or .png")
  }
  ComplexHeatmap::draw(ht, ...)
  dev.off()
  invisible(path)
}

## ===================== PLOT 1: PLAIN SAMPLE-LEVEL CORRELATION HEATMAP =======

plot_sample_correlation_heatmap <- function(cor_mat, name = "Pearson r") {
  Heatmap(
    cor_mat,
    name = name,
    col = pink_ramp_sequential(),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_side = "left",
    column_names_side = "bottom"
  )
}

## ===================== PLOT 2/3: CONDITION-LEVEL CORRELATION AS A GGPLOT TILE HEATMAP ====

# order_type: "manual" uses `order` as given; "cluster" derives the order
# from hierarchical clustering of (1 - correlation) distance.
plot_condition_correlation_tile <- function(cor_mat, order_type = c("manual", "cluster"),
                                             order = NULL) {
  order_type <- match.arg(order_type)

  if (order_type == "cluster") {
    hc <- hclust(as.dist(1 - cor_mat), method = "complete")
    order <- colnames(cor_mat)[hc$order]
  } else if (is.null(order)) {
    order <- colnames(cor_mat)
  }
  order <- order[order %in% colnames(cor_mat)]

  mat_ord <- cor_mat[order, order, drop = FALSE]
  df <- as.data.frame(as.table(mat_ord))
  colnames(df) <- c("Var1", "Var2", "value")
  df$Var1 <- factor(df$Var1, levels = rev(order))
  df$Var2 <- factor(df$Var2, levels = order)

  ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "grey70", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 4) +
    scale_fill_gradientn(colours = pink_gradient_colors, limits = c(0, 1), name = NULL) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_fixed()
}

## ===================== PLOT 4: CONDITION-LEVEL CORRELATION WITH DENDROGRAMS (COMPLEXHEATMAP) ====

## ===================== PAIRWISE CONDITION SCATTER PLOTS (E.G. PRIMARY VS D14) ====

# Per-protein mean log-intensity in two conditions, finite values only.
compute_condition_scatter_data <- function(protein_dt, cond_x, cond_y) {
  mean_dt <- protein_dt[
    GROUP %in% c(cond_x, cond_y) & !is.na(LogIntensities) & !is.na(gene.name),
    .(mean_logint = mean(LogIntensities, na.rm = TRUE)),
    by = .(gene.name, GROUP)
  ]
  wide <- dcast(mean_dt, gene.name ~ GROUP, value.var = "mean_logint")
  setnames(wide, c(cond_x, cond_y), c("x", "y"))
  wide[is.finite(x) & is.finite(y)]
}

# Scatter + linear fit + Pearson r annotation, with shared axis limits so
# multiple comparisons can be placed side by side (as in the original
# script's Primary-vs-D14 / Primary-vs-dHL60 panel pair).
plot_condition_scatter <- function(scatter_dt, cond_x, cond_y, x_lim, y_lim,
                                    line_color = "#2C8ED6", title = NULL) {
  r <- cor(scatter_dt$x, scatter_dt$y, method = "pearson")
  if (is.null(title)) title <- sprintf("%s vs %s (protein abundances)", cond_x, cond_y)

  ggplot(scatter_dt, aes(x = x, y = y)) +
    geom_point(color = "grey40", alpha = 0.45, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = line_color, linewidth = 1.4) +
    coord_fixed(xlim = x_lim, ylim = y_lim) +
    labs(
      title = title,
      subtitle = paste0("Pearson = ", sprintf("%.3f", r)),
      x = paste(cond_x, "mean log-intensity"),
      y = paste(cond_y, "mean log-intensity")
    ) +
    theme_classic(base_size = 12) +
    theme(
      aspect.ratio = 1, panel.grid = element_blank(),
      plot.title = element_text(size = 11), plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}

# Adds a bold panel letter (e.g. "G") to the top-left corner of a plot, for
# assembling lettered multi-panel figures with patchwork.
add_panel_label <- function(p, label) {
  p + annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.2, vjust = 1.5,
               size = 6, fontface = "bold")
}

plot_condition_correlation_dendro <- function(cor_mat, condition_colors = NULL) {
  row_hc <- hclust(as.dist(1 - cor_mat), method = "complete")
  col_hc <- row_hc  # symmetric matrix -> same clustering both ways

  if (is.null(condition_colors)) {
    # Fall back to a pink-family qualitative palette sized to the conditions
    pal <- grDevices::colorRampPalette(c(pink_palette["low"], pink_palette["high"]))
    condition_colors <- setNames(pal(ncol(cor_mat)), colnames(cor_mat))
  }

  ha_top <- HeatmapAnnotation(
    Condition = colnames(cor_mat),
    col = list(Condition = condition_colors),
    annotation_name_side = "left"
  )
  ha_left <- rowAnnotation(
    Condition = rownames(cor_mat),
    col = list(Condition = condition_colors),
    annotation_name_side = "top"
  )

  Heatmap(
    cor_mat,
    name = "Pearson",
    top_annotation = ha_top,
    left_annotation = ha_left,
    cluster_rows = as.dendrogram(row_hc),
    cluster_columns = as.dendrogram(col_hc),
    row_names_side = "left",
    column_names_side = "bottom",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    col = pink_ramp_sequential(),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 10))
    },
    rect_gp = gpar(col = "grey70")
  )
}
