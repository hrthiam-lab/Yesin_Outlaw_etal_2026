# Proteomics_code

R pipeline for the quantitative proteomics analysis of HSC-to-neutrophil
differentiation: raw Spectronaut/MSstats processing through PCA, sample
correlation, temporal clustering, GO enrichment, and the nuclear-envelope
gene-set heatmap. This code is interwoven across several figures/panels in
the paper rather than a single one.

## Requirements

R 4.x, with:

```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "purrr",
                    "ggplot2", "ggrepel", "gridExtra", "gtable", "patchwork",
                    "magrittr", "circlize", "devtools", "R.utils"))

BiocManager::install(c("MSstats", "MSstatsConvert", "ComplexHeatmap",
                        "clusterProfiler", "org.Hs.eg.db", "GOSemSim",
                        "UniProt.ws"))
```

Needs network access at knit/run time: the pipeline `source()`s the
Huttenhain Lab's shared `bp_utils`/`drb_utils` helper scripts live from
GitHub, and queries UniProt / `org.Hs.eg.db` for gene-name mapping.

## Input

The four MSstats output tables (`ProteinLevelData`, `GroupComparisonsData`,
`FeatureLevelData`, `CleanedPreprocessedData`) — or, optionally, the raw
Spectronaut `.tsv` report itself, which the pipeline can process from
scratch (see `RUN_RAW_IMPORT` below).

Before running anything, set these paths for your own machine:

- `data_dir` in `R/00_setup.R` and (same value) in `GO_Simplify_UpDown.R`
  — folder containing the four MSstats CSVs.
- `spectronaut.lfq.report` in `KO07CD34_Paper_Script.Rmd` — only needed if
  `RUN_RAW_IMPORT <- TRUE`.

## Running

**Entry point — `KO07CD34_Paper_Script.Rmd`.** One file, start to finish:

```
raw Spectronaut report (optional, RUN_RAW_IMPORT toggle)
        │
        ▼
MSstats dataProcess() / groupComparison()  ──or──  pre-computed CSVs
        │
        ▼
D9.2 exclusion  ──►  PCA  ──►  Pearson & Spearman correlation
        │
        ▼
Temporal cluster heatmap  ──►  volcano plots (raw + adjusted p)
        │
        ▼
GO enrichment + simplify()  ──►  nuclear-envelope gene-set heatmap
        │
        ▼
output tables: group comparisons, cleaned/processed data,
protein-level & feature-level tables
```

```r
rmarkdown::render("KO07CD34_Paper_Script.Rmd")   # from this folder
# or open in RStudio and knit top to bottom
```

`RUN_RAW_IMPORT` (near the top, section 1) chooses where it starts: `TRUE`
runs the raw Spectronaut report all the way through
`SpectronauttoMSstatsFormat()` → `dataProcess()` → `groupComparison()`
yourself; `FALSE` (default) skips straight to the four pre-computed CSVs.

**Required dependencies**, called via `source()`/`library()` rather than
duplicated inline — not optional extras:

| File | Role |
| --- | --- |
| `R/00_setup.R` – `R/06_volcano_functions.R` | Shared packages, paths, constants, and plotting/analysis functions, sourced at the top. |
| `Figure3A_Cluster_Heatmap_GO_noD9.2.R` | Temporal cluster heatmap + trend lines + GO terms (D9.2 excluded) — `source()`d for the Figure 3A section. |
| `GO_Simplify_UpDown.R` | `simplify()`-collapsed top-5 Up/Down GO terms per pairwise comparison — `source()`d for the GO section. |
| `Figure4A_Nuclear_Heatmap.R` | HPA nuclear gene-set heatmap with labeled lamina genes (LBR/LMNA/LMNB1/LMNB2) — `source()`d for the Figure 4A section; needs `gene_lists/nuclear2.tsv`. |

Each of the four scripts above also runs standalone (`Rscript <file>.R` or
open in RStudio), against the same four MSstats CSVs, if you want to
regenerate one figure on its own.

## Output

- **Tables** (`output_tables/`): group-comparison results, cleaned/processed
  protein- and feature-level data, per-cluster gene lists, GO enrichment
  results, and a per-protein z-score matrix for each heatmap figure.
- **Figures** (`figures/`): one subfolder per analysis (PCA, correlation
  heatmaps, temporal clusters, GO enrichment, volcano plots, gene-set
  heatmaps).

Both folders are created automatically on first run.
