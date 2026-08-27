## =============================================================================
## GO_Simplify_UpDown.R
## -----------------------------------------------------------------------------
## Up/Down GO Biological Process enrichment with clusterProfiler::simplify()
## collapsing semantically-redundant terms, then a top-5-per-direction bar
## plot (coral = Up, blue = Down, dot size = protein Count). Runs separately
## for each comparison (D9-D5, D14-D9, D14-D5): gene selection is
## adj.pvalue < PADJ_CUTOFF with direction from the sign of log2FC only (no
## fold-change cutoff), tested against a comparison-specific universe
## (D9-D5/D14-D9 pooled, D14-D5 separate) rather than the whole proteome.
##
## simplify() measure/cutoff: SIMPLIFY_MEASURE/SIMPLIFY_CUTOFF below. Wang,
## Lin, Jiang, Rel land roughly in [0, 1]; Resnik is raw Information Content
## and is NOT bounded to [0, 1], so a cutoff tuned for Wang won't transfer.
## The per-comparison pairwise-similarity diagnostic printed below shows the
## real numbers your data lands on for whichever measure is set.
##
## Usage:
##     Rscript GO_Simplify_UpDown.R      # standalone
##     source("GO_Simplify_UpDown.R")    # or sourced from
##                                        # KO07CD34_Paper_Script.Rmd
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
})

## ===================== 0. Make sure R is running from THIS script's folder, and that an output/ subfolder exists to save the plots into ====
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
cat("Working directory:", getwd(), "\n")

out_dir <- file.path(getwd(), "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ===================== 1. Load data -- same data_dir/files as R/00_setup.R ====

# EDIT ME: point this at your own copy of the KO07_53min_data folder (keep in
# sync with data_dir in R/00_setup.R).
data_dir <- "path/to/KO07_53min_data"

groupcomp_path <- file.path(data_dir, "20260325_GroupComparisonsData.csv")
protein_path   <- file.path(data_dir, "20260325_ProteinLevelData.csv")

if (!file.exists(groupcomp_path) || !file.exists(protein_path)) {
  stop(
    "\n\nCan't find the data files in:\n  ", data_dir,
    "\n\nExpected to find:\n  ", basename(groupcomp_path), "\n  ", basename(protein_path),
    "\n\nIf data_dir has moved, update it near the top of this script (same ",
    "value as data_dir in R/00_setup.R of the main project).\n",
    call. = FALSE
  )
}

groupcomp <- fread(groupcomp_path)
groupcomp[, gene.name := toupper(trimws(gene.name))]
groupcomp[, Label := as.character(Label)]

protein <- fread(protein_path)
protein[, gene.name := toupper(trimws(gene.name))]

## ===================== 2. Comparisons to run, and the comparison-specific universes ====
## D9-D5 and D14-D9 share a pooled universe (every gene tested in either);
## D14-D5 gets its own. This is NOT the whole detected proteome -- using that
## as a background is a different, larger null than "every gene MSstats
## actually tested in this comparison."
PADJ_CUTOFF <- 0.05

comparisons <- c("D9-D5", "D14-D9", "D14-D5")

comparisons_updown   <- c("D9-D5", "D14-D9")
gene_universe_updown <- unique(groupcomp[Label %in% comparisons_updown & !is.na(gene.name), gene.name])
gene_universe_D14D5  <- unique(groupcomp[Label == "D14-D5" & !is.na(gene.name), gene.name])

universe_for <- function(comp_label) {
  if (comp_label == "D14-D5") gene_universe_D14D5 else gene_universe_updown
}

## ===================== 3. Semantic similarity data for simplify() ===========
## Built once, reused across comparisons/directions (same org.Hs.eg.db BP
## annotation regardless of which genes are tested).

SIMPLIFY_MEASURE <- "Wang"   # "Wang", "Resnik", "Lin", "Jiang", or "Rel"
SIMPLIFY_CUTOFF  <- 0.7

go.semdata <- GOSemSim::godata(
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  computeIC = (SIMPLIFY_MEASURE != "Wang")  # only the IC-based measures need this
)

## ===================== 4. Per-comparison: build Up/Down gene sets, run enrichGO() + simplify(), and make/save the top-5 bar plot ====

top_n <- 5

run_comparison <- function(comp_label) {
  cat("\n==== ", comp_label, " ====\n", sep = "")

  comp_dt <- groupcomp[
    Label == comp_label & !is.na(gene.name) & !is.na(log2FC) & !is.na(adj.pvalue)
  ]
  universe_genes <- universe_for(comp_label)

  target.genes <- list(
    UP = unique(comp_dt[adj.pvalue < PADJ_CUTOFF & log2FC > 0, gene.name]),
    DOWN = unique(comp_dt[adj.pvalue < PADJ_CUTOFF & log2FC < 0, gene.name])
  )
  cat("UP genes:", length(target.genes$UP), " | DOWN genes:", length(target.genes$DOWN), "\n")

  enrichment.list <- lapply(names(target.genes), function(direction) {
    genes <- target.genes[[direction]]
    if (length(genes) == 0) return(data.table())

    enrichment <- clusterProfiler::enrichGO(
      gene = genes,
      universe = universe_genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = GO_PVALUE_CUTOFF,
      qvalueCutoff = GO_QVALUE_CUTOFF,
      readable = TRUE
    )
    if (is.null(enrichment) || nrow(as.data.frame(enrichment)) == 0) return(data.table())

    # Diagnostic: raw (pre-simplify) top 5 terms, for comparing against
    # plot_go_updown()'s output for the same comparison/direction.
    raw_df <- as.data.frame(enrichment)
    raw_top <- head(raw_df[order(raw_df$p.adjust), c("Description", "p.adjust", "Count")], 5)
    cat("  [", comp_label, "/", direction, "] RAW top 5 (pre-simplify):\n", sep = "")
    print(raw_top, row.names = FALSE)

    # Diagnostic: actual pairwise similarity range for SIMPLIFY_MEASURE
    # (computed the way simplify() computes it internally), so
    # SIMPLIFY_CUTOFF can be chosen from real numbers rather than assumption.
    if (nrow(raw_df) >= 2) {
      sim_mat <- tryCatch(
        GOSemSim::mgoSim(raw_df$ID, raw_df$ID, semData = go.semdata,
                          measure = SIMPLIFY_MEASURE, combine = NULL),
        error = function(e) NULL
      )
      if (!is.null(sim_mat)) {
        off_diag <- sim_mat[upper.tri(sim_mat)]
        off_diag <- off_diag[!is.na(off_diag)]
        if (length(off_diag) > 0) {
          cat("  [", comp_label, "/", direction, "] ", SIMPLIFY_MEASURE,
              " pairwise similarity over ", length(off_diag), " term pairs -- ",
              "min=", round(min(off_diag), 3),
              " median=", round(median(off_diag), 3),
              " max=", round(max(off_diag), 3),
              " (current SIMPLIFY_CUTOFF=", SIMPLIFY_CUTOFF, ")\n", sep = "")
        } else {
          cat("  [", comp_label, "/", direction, "] ", SIMPLIFY_MEASURE,
              ": no valid pairwise similarities (all NA) -- simplify() will collapse nothing.\n", sep = "")
        }
      }
    }

    # Reduce redundancy among semantically similar GO terms
    enrichment <- clusterProfiler::simplify(
      enrichment,
      cutoff = SIMPLIFY_CUTOFF,
      by = "p.adjust",
      select_fun = min,
      measure = SIMPLIFY_MEASURE,
      semData = go.semdata
    )
    setDT(as.data.table(enrichment))
  })
  names(enrichment.list) <- names(target.genes)

  cat("Terms after simplify() -- UP:", nrow(enrichment.list$UP), " | DOWN:", nrow(enrichment.list$DOWN), "\n")

  plot_df <- rbindlist(lapply(names(enrichment.list), function(direction) {
    dt <- enrichment.list[[direction]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt <- copy(dt)
    dt[, Regulation := ifelse(direction == "UP", "Up", "Down")]
    dt
  }), fill = TRUE)

  slug <- gsub("[^A-Za-z0-9]+", "", comp_label)

  if (nrow(plot_df) == 0) {
    message("No GO terms survived simplify() in either direction for ", comp_label, " -- nothing to plot.")
    return(invisible(NULL))
  }

  plot_df <- plot_df %>%
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

  p_go_top5 <- ggplot(plot_df, aes(y = Description)) +
    geom_col(aes(x = x_bar, fill = Regulation), width = 0.8) +
    geom_point(aes(x = x_bar, size = Count), color = "black") +
    scale_fill_manual(values = c("Down" = "#80B1D3", "Up" = "#D97A73")) +
    theme_bw() +
    labs(
      title = paste0(comp_label, ": top 5 enriched GO Biological Process terms (simplify()-collapsed)"),
      x = expression(-log[10]("adjusted p-value")), y = NULL,
      fill = "Regulation", size = "Count"
    ) +
    theme(panel.grid.minor = element_blank())

  print(p_go_top5)

  ggsave(file.path(out_dir, paste0("Top5_GO_BP_simplified_updown_barplot_", slug, ".pdf")), p_go_top5, width = 9, height = 6)
  ggsave(file.path(out_dir, paste0("Top5_GO_BP_simplified_updown_barplot_", slug, ".png")), p_go_top5, width = 9, height = 6, dpi = 300)
  cat("Saved plots to:", out_dir, "\n")

  invisible(list(enrichment = enrichment.list, plot = p_go_top5, plot_df = plot_df))
}

results_by_comparison <- setNames(lapply(comparisons, run_comparison), comparisons)
