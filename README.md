# Yesin_Outlaw_etal_2026

Codes for analyzing the data used in Yesin et al., paper.

> Yesin A*, Outlaw K*, Sanchez-Lopez R, Doucoure A, Hüttenhain R, Thiam HR (2026). Nuclear remodeling and high deformability emerge early while migration is progressively optimized during HSCs to neutrophils differentiation. 

*Shared first authorship

Corresponding author: Hawa Racine Thiam (hrthiam@stanford.edu), Stanford University

Keywords: Neutrophils, Neutrophils differentiation, Quantitative proteomics, Nucleus, Confined migration


## Contents

- **Proteomics_code** — R pipeline for the quantitative proteomics analysis
  of HSC-to-neutrophil differentiation (raw Spectronaut/MSstats processing
  through PCA, correlation, temporal clustering, GO enrichment, and the
  nuclear-envelope gene-set heatmap); interwoven across several figures
- **Fig4_NuclearEnvelopeIF_code** — Python pipeline for quantifying nuclear
  envelope protein intensity and nuclear shape in fixed, immunostained cells
- **Fig5_MigrationTracking_code** — Python pipeline for segmenting, tracking
  and measuring the speed of nuclei migrating in PDMS microfluidic channels

## Repository structure

```
Yesin_Outlaw_etal_2026/
├── Proteomics_code/
│   ├── KO07CD34_Paper_Script.Rmd   <- entry point, full pipeline
│   ├── Figure3A_Cluster_Heatmap_GO_noD9.2.R
│   ├── Figure4A_Nuclear_Heatmap.R
│   ├── GO_Simplify_UpDown.R
│   ├── R/
│   ├── gene_lists/
│   └── README.md
├── Fig4_NuclearEnvelopeIF_code/
│   ├── nuclear_IF_analysis.py
│   └── README.md
├── Fig5_MigrationTracking_code/
│   ├── 01_segment_and_track.py
│   ├── 02_speed_whole_channel.py
│   ├── 03_speed_by_region.py
│   └── README.md
├── requirements.txt
├── LICENSE
└── README.md
```

## Requirements

Proteomics_code: R 4.x — see that folder's README for the full package list
(CRAN + Bioconductor) and setup notes.

Fig4_NuclearEnvelopeIF_code: Python 3.12 or later, `pip install numpy
scikit-image scipy tifffile openpyxl`.

Fig5_MigrationTracking_code: Python 3.12 or later, `pip install -r requirements.txt`.


## License

[MIT](LICENSE)
