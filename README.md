

## 📋 Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Data Requirements](#data-requirements)
- [Analysis Workflow](#analysis-workflow)
- [File Organization](#file-organization)
- [Helper Functions](#helper-functions)
- [External Datasets](#external-datasets)
- [License](#license)
- [Support](#support)

## 🔬 Overview

This project analyzes immune system changes during aging and vaccination responses using:

- **Multi-omics data**: scRNA-seq, flow cytometry, Olink proteomics, total IgG, HAI assay
- **Multiple cohorts**: BRI s with longitudinal sampling and SF4 cohort
- **Comprehensive cell typing**: 71 distinct immune cell types
- **External validation**: OneK1K, Terekhova, AIFI RA, and Tissue Aging datasets
- **Time series analysis**: Multiple timepoints (Y1D0, Y1D7, Y2D0, etc.)

## 📁 Project Structure

```
├── Celltype_Mapping/           # Cell type classification hierarchy
├── Color_hex_codes/            # Standardized color schemes
├── Conda_Environment/          # Python and R environment files
├── Dataset/                    # Data storage and download scripts
├── Demographics_table/         # Cohort demographic analysis
├── Figure1-5/                  # Main manuscript figures
├── Extended_Figure1-9/         # Extended analysis figures
├── Supplementary_Figure1-2/    # Supplementary figures
├── Supplementary_Table1/       # Supplementary tables
├── helper_function/            # Reusable analysis functions
├── README.md                   # This file
└── LICENSE.txt                 # Allen Institute Software License
```

## 📊 Data Requirements

### Primary Datasets
- **BRI Cohort**: Single-cell RNA-seq, flow cytometry, Olink proteomics, HAI, total IgG assay
- **SF4 Cohort**: Single-cell RNA-seq, flow cytometry, Olink proteomics

### External Validation Datasets
- **OneK1K**: Population-scale scRNA-seq data
- **Terekhova Dataset**: Independent aging study
- **RA Dataset**: Rheumatoid arthritis cohort in AIFI
- **Tissue Aging Dataset**: Tissue aging dataset
  
### Data Formats
- **scRNA-seq**: H5AD files with AnnData format
- **Flow cytometry**: FCS files processed through standard gating
- **Proteomics**: Olink NPX values
- **Metadata**: CSV files with sample annotations
- **Total IgG/HAI Assay**: CSV files




### Key Dependencies

**Python:**
- scanpy, pandas, numpy, matplotlib, seaborn
- scikit-learn, scipy, statsmodels
- celltypist, milopy, harmonypy

**R:**
- DESeq2, limma, edgeR
- ggplot2, dplyr, tidyr


## 📂 File Organization

### Main Figures (Figure1-5)


### Extended Figures (Extended_Figure1-9)


### Dataset Structure
```
Dataset/
├── scRNA/
│   ├── BRI/                    # BRI cohort scRNA-seq data
│   └── SF4/                    # SF4 cohort scRNA-seq data
├── FlowCyto/                   # Flow cytometry data
├── Olink/                      # Proteomics data
├── HAI/                        # Hemagglutination inhibition data
└── MSD/                        # Flu specific total IgG
```

## 🛠️ Helper Functions

### Python Functions (`helper_function_IHA.py`)
- `hex_to_rgb()`: Convert hex colors to RGB
- `create_cmap_from_hex()`: Create custom colormaps
- `plot_nhood_graph()`: Plot neighborhood graphs for spatial analysis
- `grouped_obs_sum_raw()`: Aggregate raw counts by group
- `grouped_obs_mean()`: Calculate mean expression by group

### R Functions (`helper_function_IHA.r`)
- `read_pseudobulk_expression()`: Parallel data reading
- `filter_genes_and_celltype()`: Data filtering utilities
- `deseq2_analysis()`: DESeq2 differential expression wrapper
- `clr_transform()`: Centered log-ratio transformation


## 🎨 Visualization

### Color Schemes
Standardized color palettes are defined in `Color_hex_codes/`:
- **Cell types**: 71 distinct colors for immune cell types
- **Demographics**: Age groups, sex, CMV status
- **Consistent theming**: Across all figures and analyses

### Cell Type Hierarchy
The `AIFI_Reference.json` file defines a hierarchical cell type classification:
- **Major categories**: B cell, T cell, NK cell, Monocyte, DC, etc.
- **Subcategories**: Memory, Naive, Effector subtypes
- **Detailed annotations**: 71 specific cell type labels

## 📈 Key Analyses

### 1. Aging Analysis
- Composite aging scores
- Cell type frequency changes
- Gene expression alterations
- Cross-cohort validation

### 2. Vaccination Response
- Pre/post vaccination comparisons
- Vaccine-specific effects
- Time course analysis

### 3. Multi-omics Integration
- scRNA-seq + flow cytometry
- scRNA-seq + proteomics
- Cross-platform validation
- Pathway analysis


## 📄 License

This project is licensed under the Allen Institute Software License - a 2-clause BSD license with additional restrictions on commercial use. See `LICENSE.txt` for details.

**Key restrictions:**
- Commercial use requires written permission from Allen Institute
- Redistribution must include copyright notice
- No warranty provided

## 🤝 Support

This code is released AS IS without active support. While we welcome issues and questions, please understand that active responses are not guaranteed.


### Citation
If you use this code/dataset in your research, please cite our manuscript.

Gong, Q., Sharma, M., Glass, M.C. et al. Multi-omic profiling reveals age-related immune dynamics in healthy adults. Nature (2025). https://doi.org/10.1038/s41586-025-09686-5

---

**Last updated**: Oct 2025  
**Maintainer**: Qiuyu Gong
