## =============================================================================
## 05_go_enrichment_functions.R
## -----------------------------------------------------------------------------
## GO Biological Process enrichment. Consolidates the ~5 overlapping GO chunks
## in the original script (one embedded inside the Figure 3A chunk, one
## standalone "per-cluster enrichGO + dotplot" chunk, one final polished
## pink-colored dotplot chunk, plus the up/down bar-plot version for pairwise
## comparisons) into a small set of reusable functions.
##
## The per-cluster / per-comparison DOT PLOTS (size = gene count, color =
## -log10 adjusted p-value) now use the shared pink_gradient_colors palette
## from 00_setup.R, per your request to match the Pearson heatmap.
##
## The up/down BAR plots (make_go_barplot equivalent below) intentionally
## keep their two-color Up/Down encoding rather than switching to the pink
## gradient -- that color there encodes direction of change, not p-value, so
## collapsing it to one palette would make the plot unreadable. Let me know if
## you'd like those recolored too.
##
## Usage:
##     source("R/00_setup.R"); source("R/01_data_loading.R")
##     source("R/05_go_enrichment_functions.R")
## =============================================================================

## ===================== CORE ENRICHMENT CALL =================================

run_go_enrichment <- function(genes, universe = NULL, ont = "BP",
                               pvalueCutoff = GO_PVALUE_CUTOFF,
                               qvalueCutoff = GO_QVALUE_CUTOFF) {
  if (length(genes) == 0) return(NULL)
  enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
}

# Runs GO enrichment separately for each cluster in `cluster_genes` (a named
# list of gene vectors, e.g. cluster_result$cluster_genes from
# cluster_temporal_profiles()) and returns one combined data.table with a
# Cluster column, ready for top_go_terms()/plot_go_dotplot().
run_go_by_cluster <- function(cluster_genes, universe) {
  go_list <- lapply(seq_along(cluster_genes), function(i) {
    ego <- run_go_enrichment(cluster_genes[[i]], universe = universe)
    if (is.null(ego)) return(NULL)
    df <- as.data.frame(ego)
    if (nrow(df) == 0) return(NULL)
    df$Cluster <- paste0("Cluster ", i)
    df
  })
  rbindlist(go_list, fill = TRUE)
}

## ===================== RESHAPING FOR PLOTS ==================================

# Keeps the top `n` GO terms per Cluster (by p.adjust) and adds the columns
# the dot plot needs (gene ratio as a number, -log10 padj, wrapped labels).
top_go_terms <- function(go_df, n = 5, group_col = "Cluster", cluster_levels = NULL) {
  if (nrow(go_df) == 0) return(go_df)

  go_top <- as.data.frame(go_df) %>%
    group_by(.data[[group_col]]) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    mutate(
      GeneRatioNum = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      negLog10Padj = -log10(p.adjust),
      Description = str_wrap(Description, width = 35)
    )

  if (!is.null(cluster_levels)) {
    go_top[[group_col]] <- factor(go_top[[group_col]], levels = cluster_levels)
  }
  go_top
}

## ===================== DOT PLOT (PINK PALETTE) ==============================

plot_go_dotplot <- function(go_top, facet_col = "Cluster",
                             title = "GO Biological Process by cluster") {
  ggplot(go_top, aes(x = GeneRatioNum, y = Description, size = Count, color = negLog10Padj)) +
    geom_point() +
    facet_wrap(stats::as.formula(paste("~", facet_col)), ncol = 1, scales = "free_y") +
    scale_color_gradientn(colours = pink_gradient_colors, name = expression(-Log[10]("(P.adjust)"))) +
    labs(title = title, x = "Gene ratio", y = NULL, size = "Count") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
}

# Small per-cluster grobs (no legend, no title) for the combined multi-panel
# figure built by make_combined_cluster_figure() in 04_temporal_cluster_functions.R.
make_go_grobs <- function(go_top, clusters) {
  lapply(clusters, function(cl_name) {
    df <- go_top %>% filter(Cluster == cl_name)
    if (nrow(df) == 0) {
      p <- ggplot() + theme_void()
    } else {
      p <- ggplot(df, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum),
                           size = Count, color = negLog10Padj)) +
        geom_point() +
        scale_color_gradientn(colours = pink_gradient_colors) +
        theme_bw() +
        labs(x = NULL, y = NULL) +
        theme(
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.grid.minor = element_blank(),
          plot.margin = margin(2, 2, 2, 2),
          legend.position = "none"
        )
    }
    ggplotGrob(p)
  })
}

# Extracts just the color/size legend as a standalone grob, for placement in
# the combined figure's 4th column.
build_go_legend <- function(go_top) {
  p <- ggplot(go_top, aes(x = GeneRatioNum, y = Description, size = Count, color = negLog10Padj)) +
    geom_point() +
    scale_color_gradientn(colours = pink_gradient_colors) +
    theme_bw() +
    labs(color = expression(-Log[10]("(P.adjust)")), size = "Count")
  gtable::gtable_filter(ggplotGrob(p), "guide-box")
}

## ===================== UP/DOWN GO BAR+DOT PLOT FOR A SINGLE PAIRWISE COMPARISON ====

# Direct equivalent of `make_go_barplot()` in the original script. Splits a
# comparison's significant genes into Up/Down, runs GO BP on each half, and
# returns both the combined term table and the bar+dot plot.
plot_go_updown <- function(comp_label, groupcomp_dt, universe_genes,
                            padj_cutoff = PADJ_CUTOFF, top_n = 5) {
  sub_dt <- groupcomp_dt[
    Label == comp_label & !is.na(gene.name) & !is.na(log2FC) & !is.na(adj.pvalue)
  ]
  up_genes   <- unique(sub_dt[adj.pvalue < padj_cutoff & log2FC > 0, gene.name])
  down_genes <- unique(sub_dt[adj.pvalue < padj_cutoff & log2FC < 0, gene.name])

  up_df <- data.frame()
  down_df <- data.frame()

  ego_up <- run_go_enrichment(up_genes, universe = universe_genes)
  if (!is.null(ego_up)) {
    up_df <- as.data.frame(ego_up)
    if (nrow(up_df) > 0) up_df$Regulation <- "Up"
  }
  ego_down <- run_go_enrichment(down_genes, universe = universe_genes)
  if (!is.null(ego_down)) {
    down_df <- as.data.frame(ego_down)
    if (nrow(down_df) > 0) down_df$Regulation <- "Down"
  }

  res_df <- rbindlist(list(as.data.table(up_df), as.data.table(down_df)), fill = TRUE)
  if (nrow(res_df) == 0) return(list(table = res_df, plot = NULL))

  plot_df <- res_df %>%
    as.data.frame() %>%
    group_by(Regulation) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      negLog10Padj = -log10(p.adjust),
      Description = str_wrap(Description, width = 34),
      x_bar = ifelse(Regulation == "Down", -negLog10Padj, negLog10Padj)
    ) %>%
    arrange(Regulation, abs(x_bar))
  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)

  p <- ggplot(plot_df, aes(y = Description)) +
    geom_col(aes(x = x_bar, fill = Regulation), width = 0.8) +
    geom_point(aes(x = x_bar, size = Count), color = "black") +
    scale_fill_manual(values = c("Down" = "#80B1D3", "Up" = "#D97A73")) +
    theme_bw() +
    labs(
      title = paste0(comp_label, " GO Biological Process"),
      x = expression(-log[10]("adjusted p-value")), y = NULL,
      fill = "Regulation", size = "Count"
    ) +
    theme(panel.grid.minor = element_blank())

  list(table = plot_df, plot = p)
}

## ===================== UTILITY: FLATTEN A NESTED LIST OF ENRICHRESULT OBJECTS ====

# General-purpose helper for turning a nested list of enrichGO() results
# (however deeply nested) into one flat data.frame with a `source_path`
# column recording where each row came from. Not wired into the main
# pipeline automatically -- point it at whatever named list of enrichResult
# objects you build (e.g. from run_go_by_cluster() variants), for example:
#   flatten_go_results(list(D9_D5 = ego_up, D14_D9 = ego_down))
flatten_go_results <- function(x, path = "root") {
  out <- list()
  df <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (!is.null(df) && nrow(df) > 0 && all(c("Description", "geneID") %in% colnames(df))) {
    df$source_path <- path
    return(list(df))
  }
  if (is.list(x)) {
    element_names <- names(x)
    if (is.null(element_names)) element_names <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      nm <- element_names[i]
      if (is.na(nm) || nm == "") nm <- as.character(i)
      out <- c(out, flatten_go_results(x[[i]], path = paste(path, nm, sep = " / ")))
    }
  }
  out
}
