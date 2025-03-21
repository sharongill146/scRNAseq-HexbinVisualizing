
# Install and Load Necessary Packages
#=====================================
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggdensity", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("hexbin", quietly = TRUE)) install.packages("dplyr")

library(Seurat)   # For scRNA-seq analysis
library(ggplot2)  # For visualizations
library(ggdensity) #For density plots
library(dplyr)    # For data manipulation
library(patchwork) # For combining plots
library(hexbin) # For hexagonal binning 

# Generate a distinct color palette using RColorBrewer
#=====================================================
library(RColorBrewer)
num_classes <- length(unique(df_clean$General_type))

# Use Set3 or other palettes that support more colors
colors <- brewer.pal(min(num_classes, 12), "Set3")

ggplot(df_clean, aes(x = UMAP_1, y = UMAP_2, fill = General_type)) +
  geom_hdr() +
  geom_point(shape = 21, alpha = 0.6) +
  scale_fill_manual(values = colors) +
  labs(title = "UMAP Subpopulation Visualization",
       x = "UMAP Dimension 1", y = "UMAP Dimension 2", fill = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "right", axis.text = element_text(size = 12))

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


