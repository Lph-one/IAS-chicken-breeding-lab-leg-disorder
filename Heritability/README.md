# Genotype-Based PCA and Heritability Analysis Pipeline (CA vs CT)

A reproducible pipeline for genotype-based PCA and large-scale heritability estimation using GCTA and HIBLUP.

This repository provides a reproducible workflow for genotype-based analysis, including group-wise data splitting, principal component analysis (PCA), and heritability estimation using both GCTA and HIBLUP. The pipeline is designed for comparative analysis between two groups (CA and CT) and supports large-scale multi-trait analysis in a parallelized manner.

---

## Overview

The pipeline consists of four main steps:

1. Split genotype data into CA and CT groups  
2. Perform PCA for each group  
3. Estimate heritability using GCTA  
4. Estimate heritability using HIBLUP  

All steps are modular and can be executed independently or as a full workflow.

---

## Project Structure

.
├── split_genotype.sh # Split genotype data into CA and CT
├── run_pca.sh # Perform PCA for each group
├── run_CA_h2_GCTA.sh # GCTA heritability for CA
├── run_CT_h2_GCTA.sh # GCTA heritability for CT
├── run_CA_h2_Hiblup.sh # HIBLUP heritability for CA
├── run_CT_h2_Hiblup.sh # HIBLUP heritability for CT
├── run_h2_by_group_GCTA.py # Core GCTA pipeline (parameterized)
├── run_h2_by_group_Hiblup.sh # Core HIBLUP pipeline (parameterized)


---

## Requirements

The following software is required:

- PLINK (v1.9 or later)
- GCTA (v1.93 or later)
- HIBLUP
- GNU Parallel
- Python ≥ 3.7

Make sure all tools are available in your environment (`$PATH`).  
If not, you can specify their paths in the scripts.

---

## Workflow
### 1. Split genotype data

bash split_genotype.sh

### 2. Run PCA
bash run_pca.sh

### 3. GCTA heritability analysis
bash run_CA_h2_GCTA.sh
bash run_CT_h2_GCTA.sh

### 4. HIBLUP heritability analysis
bash run_CA_h2_Hiblup.sh
bash run_CT_h2_Hiblup.sh


---

## Output
For each group (CA and CT), the pipeline generates:

PCA results
GRM files (GCTA)
Trait-level analysis outputs
Summary heritability files





## License
For academic and research use only.