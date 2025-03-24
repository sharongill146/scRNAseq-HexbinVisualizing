# Install and Load Necessary Packages
#=====================================
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggdensity", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("hexbin", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("entropy", quietly = TRUE)) install.packages("dplyr")

library(Seurat)   # For scRNA-seq analysis
library(ggplot2)  # For visualizations
library(ggdensity) #For density plots
library(dplyr)    # For data manipulation
library(patchwork) # For combining plots
library(hexbin) # For hexagonal binning 
library(ggrepel) #For labelling cell types
library(viridisLite)
library(reshape2)
library(entropy)

# Generate a distinct color palette using RColorBrewer
#=====================================================
library(RColorBrewer)
num_classes <- length(unique(df_clean$General_type))

# Use Set3 or other palettes that support more colors
colors <- brewer.pal(max(min(num_classes, 12), 3), "Set3")

# Load Single-Cell Dataset (Adjust Path Accordingly)
#====================================================
dataset <- "/Users/sharongill/Desktop/cac9d341-6b8d-4ed3-8eaa-054688787840.rds"
seurat_obj <- readRDS(dataset)

# Inspect Metadata (Sample Information)
#======================================
# Metadata contains annotations for each cell (e.g., cell type, sample origin)
meta_data <- seurat_obj@meta.data
head(meta_data)  # View first few rows

# Summary statistics of metadata
summary(meta_data)

# Extract UMAP Coordinates (for visualization)
if ("umap" %in% names(seurat_obj@reductions)) {
  umap_coords <- seurat_obj@reductions$umap@cell.embeddings
  head(umap_coords)
} else {
  message("UMAP coordinates not found in dataset.")
}

# Visualize Cell Type Distribution (if available) 
# Output of bar graph with cell count over cell type.
#====================================================

if ("cell_type" %in% colnames(meta_data)) {
  ggplot(meta_data, aes(x = cell_type)) +
    geom_bar(fill = "steelblue") +
    theme_minimal() +
    labs(title = "Cell Type Distribution", x = "Cell Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  message("No 'cell_type' column found in metadata.")
}

# Save Processed Metadata for Reference
#======================================

write.csv(meta_data, "metadata_output.csv", row.names = TRUE)

