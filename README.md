# Cd70 Deficiency Accelerates Bcl6-driven Diffuse Large B Cell Lymphoma and Impairs CD4⁺ Cytotoxic T-cell Immune Surveillance

## I. Installation and Environment Setup

To ensure full reproducibility of the analyses described in this repository, please follow the setup instructions below.

### 1. System Requirements

- Operating System: Linux (tested on CentOS 7)

- R version: 4.1.1

- Platform: x86_64-conda-linux-gnu (64-bit)


### 2. Create and Activate the Conda Environment
You can create a compatible environment with R 4.1.1 as follows:

```{bash}
conda create -n cd70_dlbcl_env r-base=4.1.1 -y
conda activate cd70_dlbcl_env
conda install -c conda-forge r-essentials r-devtools r-tidyverse -y
```

Some analyses in this repository use Python-based Jupyter notebooks for visualization and statistical analyses.

```{bash}
conda install -c conda-forge python=3.9 jupyterlab numpy pandas scipy seaborn matplotlib scikit-learn -y

pip install pvalannot heatmapannot
```

### 3. Install Required R Packages

```{R}
# Base packages
install.packages(c(
  "tidyverse", "data.table", "ggplot2", "ggpubr", "ComplexHeatmap", 
  "RColorBrewer", "viridis", "viridisLite", "Polychrome", 
  "Rcpp", "here", "sp", "forcats", "stringr", "dplyr", 
  "purrr", "readr", "tidyr", "tibble"
))

# Additional dependencies
install.packages(c(
  "Seurat", "SeuratObject", "harmony", "presto", "devtools", "usethis", 
  "infercnv", "edgeR", "limma", "circlize", "patchwork", "rstatix", 
  "cowplot", "ComplexHeatmap"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges", "SummarizedExperiment", "SingleCellExperiment",
  "Biobase", "MatrixGenerics", "GenomeInfoDb", "S4Vectors"
))
```
### 4. Optional: Reproduce the Original Environment
```{bash}
conda env create -f cd70_dlbcl_env.yml
```

## II. Explanation of the scripts
- **0_data_preprocessing.ipynb**

    Reads the count matrix and removes low quality cells

- **1_split_merged_data_object_and_find_marker.ipynb**

    Splits Cd3⁺ Cd4⁻, Cd3⁺ Cd8⁻, and B-cell subsets based on marker expression, and identifies marker genes for each cluster

- **2_remove_clusters_and_subclustering.ipynb**

    Removes underrepresented clusters and performs additional clustering.

- **3_add_celltype_HEX_color_and_annotation.ipynb**
    
    Annotates each cluster manually based on marker expression and stores cell type information in the seurat.

- **4_plot_DE_genes_between_clusters_on_final_seurat**
    
    Plots top marker genes for the final Seurat object.

- **5_plot_clustered_doplot_{space}.ipynb**
    
    Plots a subset of cell-type markers. This includes final version of clustered dot plots shown in the manuscript figures.

- **6_trust4_pl/**
    
    Contains the original Perl script used for submitting TRUST4 SLURM array jobs.

- **6_trust4_diversity/**
    
    Calculates clonality, entropy, and richness for TCR/BCR inferred by TRUST4.

- **7_trus4_analysis_{space}.ipynb**
    
    Overlays clonality, entropy, and richness on UMAP visualizations. This includes the the UMAP shown in the manuscript figures.

- **8_calculate_celltype_enrichment_per_mouse.ipynb**

    Plots the fraction of each cell type per sample.

- **9_cellphoneDB_specific_celltype_analysis_specific_B_cell_cluster**

    Prepares the input data for cellphonedb and visualizes cellphonedb output.

- **9_cellphoneDB/**

    Includes SLURM scripts for running CellPhoneDB analyses for each genotype and timepoint combination.

- **10_infercnv/**

    Contains R and SLURM scripts for running inferCNV on individual samples.

- **10_define_inferCNV_pos_cells_in_each_samples.ipynb**

    Identifies inferCNV-positive and inferCNV-negative cells from inferCNV output files.

- **11_R_package_versions.ipynb**

    Lists the R package versions used in the analyses.

- **12_percent_of_infercnv_positive_cells_per_sample.ipynb**

    Calculates the percentage of inferCNV-positive cells per sample and compares results between genotypes and timepoints.

- **12_inferCNV_with_B_cell_clonal_expansion.ipynb**

    Determines the percentage of inferCNV-positive cells within clonally expanded cell populations.


## III. Reproducing Data and Figures from the Manuscript

### 1. Access full dataset from EGA and build Seurat object
The single-cell RNA sequencing count matrices are available in the European Genome-phenome Archive (EGA) under accession number EGAD00001012104. You can apply for data access through the EGA portal.

### 2. Clone github repository
```{bash}
git clone https://github.com/baigal628/cd70_dlbcl_manuscript.git

cd cd70_dlbcl_manuscript/

cd data/
```

### 3. Demo dataset
We have provided one dataset to test the code. But the result will not reproduce what is shown in the manuscript.

```{bash}
ls data/all_analysis/Pool105_17/

barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

```

### 4. Analyze downloaded data set

Download the file all_analysis.tar.gz and place it inside the data/ folder:

```{bash}
cd data/

tar xvzf all_analysis.tar.gz

ls all_analysis/ 
```
You should see multiple data folders named according to sample pool IDs (e.g., Pool105_1/).


Next, navigate to the scripts/ folder and open the corresponding Jupyter notebooks to construct Seurat objects and perform downstream analyses.
In each notebook, set the working directory to the path where the repository is cloned. For example:
```{bash}
setwd("/path/to/cd70_dlbcl_manuscript/")
```