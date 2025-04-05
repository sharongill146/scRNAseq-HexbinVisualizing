# scRNAseq-HexbinVisualizing
Data set on immune cells in a blood sample undergone chemotherapy with improved visualizations using Hexbin Plots and Hexbin.R Functions

# HexnicPlot: Enhanced scRNA-seq Visualization with Hexagonal Binning  

HexnicPlot is an R-based tool for improving UMAP visualizations in single-cell RNA sequencing (scRNA-seq) data. Instead of plotting individual cells, it aggregates them into hexagonal bins, allowing users to visualize majority cell types, Shannon entropy, and gene expression patterns while reducing visual clutter.

## Features  
✔️ Aggregates scRNA-seq data into hexagonal bins for better visualization  
✔️ Supports Shannon entropy calculation to highlight cellular heterogeneity  
✔️ Integrates seamlessly with Seurat for easy workflow compatibility  
✔️ Customizable parameters (bin density, color palettes, visualization modes)  

## Dataset Information

### Biological Dataset  

This tool is demonstrated using an scRNA-seq dataset on immune cells exposed to low-dose radiation (LDR), including:  
- Natural Killer Cells  
- CD4+ Effector Memory T Cells  
- CD8+ Effector Memory T Cells  
- Classical Monocytes  
- Memory B Cells  

### Dataset Source 

The dataset used in this project is publicly available from CellxGene.
https://cellxgene.cziscience.com/e/9c1b5626-58df-4401-ae7b-f66d068c1551.cxg/

### Future Directions  

- Extend HexnicPlot to multi-omics datasets (epigenetics, proteomics, metabolomics)  
- Improve entropy-based clustering to refine cellular heterogeneity analysis  
- Add interactive visualizations for real-time exploration  

## Contributors  

Developed by Sharon Gill
With Motivation for HexnicPlot Functions by Nicolas Ho (Faculty of Medicine, University of Ottawa)
