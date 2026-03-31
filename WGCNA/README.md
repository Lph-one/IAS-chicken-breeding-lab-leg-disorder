# WGCNA Analysis and Visualization for CA and CT Groups

This repository provides a reproducible workflow for performing Weighted Gene Co-expression Network Analysis (WGCNA) on metabolomics data from two groups (CA and CT), together with downstream visualization including clustered heatmaps and Sankey plots for cross-group module comparison.

The analysis starts from two input datasets (`CAmeta.csv` and `CTmeta.csv`), where rows represent metabolites and columns represent samples. The scripts `run_CA_WGCNA.R` and `run_CT_WGCNA.R` are used to perform WGCNA separately for each group, including data quality control, soft-threshold selection, network construction, module detection, and module merging. The resulting outputs (TOM matrices, module color annotations, eigengenes, and Cytoscape files) are stored under the `output/CA_WGCNA/` and `output/CT_WGCNA/` directories.

To visualize network structure, the script `heatmap_by_group.R` generates TOM-based heatmaps for both groups. Metabolites are ordered according to their module assignments, and module boundaries are explicitly highlighted on the heatmap. For efficiency, partial TOM matrices (`*_part.csv`) are used for visualization. This step also produces row annotation files (`*_Heatmap_row_info.csv`), which record the ordering and module membership of each metabolite and serve as input for downstream comparison.

The script `sankey_plot.R` integrates the results from both groups and generates a Sankey diagram to visualize how metabolite modules correspond between CA and CT. Only shared metabolites are considered, and module transitions are displayed using consistent color mapping based on WGCNA module colors. This allows intuitive assessment of module preservation, splitting, or reorganization across conditions.

All scripts are written using relative paths and can be executed in any environment without modification. The required R packages include `WGCNA`, `pheatmap`, `dplyr`, `grid`, `tidyr`, `ggplot2`, and `ggalluvial`.

A typical workflow is to first run the WGCNA scripts for both groups, followed by the heatmap generation, and finally the Sankey plot:

```r
source("run_CA_WGCNA.R")
source("run_CT_WGCNA.R")
source("heatmap_by_group.R")
source("sankey_plot.R")