# Majority plot: Identifies the dominant cell type in each hexbin
# ================================================================
HexnicPlot(seurat_obj, label_col = "cell_type", plot_type = "majority")

# Shannon entropy plot: Measures the diversity of cell types in each hexbin
# =========================================================================
HexnicPlot(seurat_obj, label_col = "cell_type", plot_type = "entropy")

# Gene expression plot: Visualizes the expression of a specific gene 
# ===================================================================
gene_name <- "ENSG00000238009"
HexnicPlot(seurat_obj, label_col = "cell_type", plot_type = "gene", gene_name = gene_name, assay = "RNA", layer = "counts")

