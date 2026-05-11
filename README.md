# Tumor Educated Platelets (TEPs) Differential Centrality Measures

This repository contains the code for the analysis of differential centrality measures in tumor educated platelets

The file `main.R` contains the code for conducting the analysis and creating train and test data.

The file `randomization.R` contains the code for randomization of AUC values. 

## System Requirements

### Operating System

The model was developed and tested on a machine running **Windows 11**.

Training and test were performed on a system equipped with a 12th Gen Intel(R) Core(TM) i7-12700H @ 2.30 GHz CPU and 16 GB of RAM.

### R packages

The analysis was conducted using **R version 4.4.0** with the following package versions:

- `clusterProfiler`: 4.14.6  
- `org.Hs.eg.db`: 3.20.0  
- `stringr`: 1.6.0  
- `dplyr`: 1.1.4  
- `pathview`: 1.46.0  
- `igraph`: 2.2.1
- `network`: 1.19.0  
- `DESeq2`: 1.46.0  
- `DT`: 0.34.0  
- `ggplot2`: 4.0.1  
- `NetworkToolbox`: 1.4.4  
- `psych`: 2.5.6  
- `GGally`: 2.4.0  
- `sna`: 2.8  
- `BiocGenerics`: 0.52.0  
- `SummarizedExperiment`: 1.36.0  
- `reshape2`: 1.4.5  
- `vegan`: 2.7.2  
- `ggrepel`: 0.9.6  
- `ReactomePA`: 1.50.0  
- `readxl`: 1.4.5  
- `tibble`: 3.3.0  

### Python Packages

The model was trained and tested using **Python 3.12.3** with the following package versions:

- `python`: 3.12.3  
- `pandas`: 2.2.2  
- `numpy`: 1.26.4  
- `scikit-learn`: 1.5.1  
- `xgboost`: 3.0.0  
- `catboost`: 1.2.8  
- `matplotlib`: 3.9.2
