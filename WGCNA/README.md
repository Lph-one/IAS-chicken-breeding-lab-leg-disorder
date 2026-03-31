# WGCNA Analysis and Visualization for CA and CT Groups

A reproducible workflow for Weighted Gene Co-expression Network Analysis (WGCNA) and cross-group module comparison using metabolomics data.

---

## Overview

This repository provides a reproducible workflow for performing Weighted Gene Co-expression Network Analysis (WGCNA) on metabolomics data from two groups (CA and CT), together with downstream visualization including clustered heatmaps and Sankey plots for cross-group module comparison.

The analysis starts from two input datasets:

- `CAmeta.csv`
- `CTmeta.csv`

In these files:
- Rows represent metabolites  
- Columns represent samples  

---

## Workflow

The pipeline consists of the following steps:

1. Perform WGCNA separately for CA and CT  
2. Generate TOM-based heatmaps  
3. Compare modules across groups using Sankey plots  

---

## Project Structure
run_CA_WGCNA.R # WGCNA analysis for CA
run_CT_WGCNA.R # WGCNA analysis for CT
heatmap_by_group.R # TOM heatmap visualization
sankey_plot.R # Cross-group module comparison (Sankey)

---

## WGCNA Analysis
run_CA_WGCNA.R
run_CT_WGCNA.R

---

## Heatmap Visualization
heatmap_by_group.R

---

## Sankey Plot (Module Comparison)
sankey_plot.R

---

## Requirements
The following R packages are required:

WGCNA
pheatmap
dplyr
grid
tidyr
ggplot2
ggalluvial

---

## Usage
source("run_CA_WGCNA.R")

source("run_CT_WGCNA.R")

source("heatmap_by_group.R")

source("sankey_plot.R")





## License
For academic and research use only.

