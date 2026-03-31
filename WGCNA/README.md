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

## Workflow
### 1. Perform WGCNA separately for CA and CT 
run_CA_WGCNA.R
run_CT_WGCNA.R

### 2. Generate TOM-based heatmaps  
heatmap_by_group.R

### 3. Compare modules across groups using Sankey plots  
sankey_plot.R

---

## Output
Results are saved in:
output/CA_WGCNA/
output/CT_WGCNA/
Including:
TOM matrices
Module color assignments
Module eigengenes
Cytoscape network files







## License
For academic and research use only.

