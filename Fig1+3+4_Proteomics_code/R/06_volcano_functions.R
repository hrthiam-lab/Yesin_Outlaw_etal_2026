## =============================================================================
## 06_volcano_functions.R
## -----------------------------------------------------------------------------
## Volcano plots for pairwise comparisons. The original script had two
## versions of essentially the same plot: an earlier draft that auto-labeled
## the top N significant genes, and a later "Figure 3" version with hand-picked
## gene labels and extra styling (axis hooks, wider x-range). Both behaviors
## are kept here as options on one function instead of duplicated code.
##
## Usage:
##     source("R/00_setup.R"); source("R/01_data_loading.R")
##     source("R/06_volcano_functions.R")
## =============================================================================

# Adds a Significance column (Up/Down/Not) and -log10(chosen p-value metric)
# to a groupcomp subset for one comparison Label. y_metric picks which column
# drives the y-axis: "pvalue" (raw, the historical default) or "adj.pvalue"
# (BH-adjusted -- matches the significance calling, which always uses
# adj.pvalue regardless of y_metric).
add_volcano_columns <- function(dt, padj_cutoff = PADJ_CUTOFF, lfc_cutoff = LFC_CUTOFF,
                                 y_metric = c("pvalue", "adj.pvalue")) {
  y_metric <- match.arg(y_metric)
  copy(dt)[
    ,
    Significance := fifelse(
      adj.pvalue < padj_cutoff & log2FC > lfc_cutoff, "Up",
      fifelse(adj.pvalue < padj_cutoff & log2FC < -lfc_cutoff, "Down", "Not")
    )
  ][, neglog10p := -log10(get(y_metric))]
}

# label_mode = "top_n"  -> auto-labels the top `top_n` Up and Down genes by adj.pvalue
# label_mode = "custom" -> labels exactly the genes in `label_genes`
# y_metric   = "pvalue" (default, unchanged historical behavior) or
#              "adj.pvalue" (plots -log10(adjusted p-value) instead, so the
#              axis matches the Up/Down significance calling, which is always
#              adjusted-p-value-based regardless of this setting)
plot_volcano <- function(comp_label, groupcomp_dt,
                          padj_cutoff = PADJ_CUTOFF, lfc_cutoff = LFC_CUTOFF,
                          label_mode = c("top_n", "custom"),
                          top_n = 20, label_genes = NULL,
                          x_lim = NULL, y_lim = NULL,
                          title = comp_label, x_label = NULL,
                          y_metric = c("pvalue", "adj.pvalue"),
                          sig_colors = c(Down = "#7DB2D9", Not = "#BFBFBF", Up = "#D9857B")) {
  label_mode <- match.arg(label_mode)
  y_metric <- match.arg(y_metric)

  volc <- groupcomp_dt[Label == comp_label & !is.na(log2FC) & !is.na(pvalue) & !is.na(adj.pvalue)]
  volc <- add_volcano_columns(volc, padj_cutoff, lfc_cutoff, y_metric = y_metric)

  if (is.null(x_lim)) x_lim <- c(floor(min(volc$log2FC, na.rm = TRUE)), ceiling(max(volc$log2FC, na.rm = TRUE)))
  if (is.null(y_lim)) y_lim <- c(0, ceiling(max(volc$neglog10p, na.rm = TRUE)))
  if (is.null(x_label)) x_label <- sprintf("Log2FC (%s)", comp_label)

  y_axis_label <- if (y_metric == "adj.pvalue") {
    expression(-Log[10]("adjusted p-value"))
  } else {
    expression(-Log[10](p - value))
  }
  ref_line_label <- if (y_metric == "adj.pvalue") {
    paste0("adj. pval = ", padj_cutoff)
  } else {
    paste0("pval = ", padj_cutoff)
  }

  lab_df <- if (label_mode == "custom") {
    volc[gene.name %in% label_genes]
  } else {
    rbind(
      volc[Significance == "Up"][order(adj.pvalue, -abs(log2FC))][1:min(.N, top_n)],
      volc[Significance == "Down"][order(adj.pvalue, -abs(log2FC))][1:min(.N, top_n)]
    )
  }

  base_theme <- theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.line = element_line(linewidth = 0.8, lineend = "square"),
      plot.margin = margin(t = 10, r = 25, b = 10, l = 10)
    )

  ggplot(volc, aes(x = log2FC, y = neglog10p, color = Significance)) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_manual(values = sig_colors) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60") +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = gene.name), size = 3, max.overlaps = Inf,
      show.legend = FALSE, box.padding = 0.4, point.padding = 0.2,
      min.segment.length = 0, seed = 123
    ) +
    coord_cartesian(xlim = x_lim, ylim = y_lim, clip = "off") +
    labs(title = title, x = x_label, y = y_axis_label, color = "Significance") +
    annotate("text", x = -lfc_cutoff, y = 0.2, label = paste0("log2FC = -", lfc_cutoff), size = 4) +
    annotate("text", x = lfc_cutoff, y = 0.2, label = paste0("log2FC = ", lfc_cutoff), size = 4) +
    annotate("text", x = x_lim[2], y = 0.2, label = ref_line_label, hjust = 1, size = 4) +
    base_theme
}

## ===================== TOP HITS TABLE =======================================

# Top `n` up- and down-regulated genes for a comparison, ranked by log2FC
# among adj.pvalue-significant hits. Direct replacement for the original
# script's loop over `global.pairwiseComparison` (an object that was never
# defined in the pasted code) -- uses `groupcomp` instead, which holds the
# same data (it's GroupComparisonsData.csv, loaded in 01_data_loading.R).
top_hits_table <- function(groupcomp_dt, comp_label, padj_cutoff = PADJ_CUTOFF, n = 50) {
  sig <- groupcomp_dt[
    Label == comp_label & !is.na(adj.pvalue) & adj.pvalue < padj_cutoff & is.finite(log2FC)
  ]
  list(
    up   = head(sig[log2FC > 0][order(-log2FC)], n)[, .(gene.name, log2FC, adj.pvalue)],
    down = head(sig[log2FC < 0][order(log2FC)], n)[, .(gene.name, log2FC, adj.pvalue)]
  )
}
