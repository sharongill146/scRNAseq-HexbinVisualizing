#Working with Hexbins - Centroids, Shannon Entropy, Gene Expression
#==========================================================================================
# Ensure df_clean contains the correct data
if (!exists("df_clean")) {
    stop("Error: df_clean does not exist. Ensure you have loaded the correct dataset.")
}

# Extract UMAP coordinates and cell types from df_clean
x <- df_clean$UMAP_1
y <- df_clean$UMAP_2
cell_types <- seurat_obj@meta.data$cell_type
cell_type <- df_clean$cell_type  # Ensure 'cell_type' exists

# Convert to data frame
df <- data.frame(UMAP_1 = x, UMAP_2 = y, cell_type = cell_types)

# Step 1: Calculate optimal bins based on the number of cells
num_points <- nrow(df)
optimal_bins <- round(sqrt(num_points / 2))  # Square root heuristic for bin count

#Compute centroid of each cell type
centroids <- df %>%
    group_by(cell_type) %>%
    summarise(
        centroid_x = mean(UMAP_1, na.rm = TRUE),
        centroid_y = mean(UMAP_2, na.rm = TRUE)
    )

#Create Hexbin plot with centroids and bin count annotation
hex_plot <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_hex(bins = optimal_bins) +  # Hexbin plot with optimal bin count
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Cell Count") +  # Color scale
    geom_point(data = centroids, aes(x = centroid_x, y = centroid_y), color = "red", size = 2) +  # Centroid points
    geom_text(data = centroids, aes(x = centroid_x, y = centroid_y, label = cell_type), 
              color = "black", vjust = -1, size = 4) +  # Label cell type names at centroids
    annotate("text", x = max(df$UMAP_1) * 0.9, y = max(df$UMAP_2) * 0.1, 
             label = paste("Hexbins Used:", optimal_bins), 
             color = "black", size = 5, fontface = "bold", hjust = 1) +  # Display hexbin count
    labs(title = "Hexbin Plot of UMAP Coordinates",
         x = "UMAP Dimension 1", 
         y = "UMAP Dimension 2") +
    theme_minimal()

# Print the plot
print(hex_plot)

# Step 4: Create Hexbin plot, coloring by cell type

#Generate distinct colors dynamically based on the number of cell types
num_cell_types <- length(unique(df$cell_type))
palette_colors <- viridisLite::viridis(num_cell_types + 1, option = "C")  # "C" gives a good color range


hex_plot <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
    stat_binhex(aes(fill = cell_type), bins = optimal_bins, alpha = 0.7) +  # Color by cell type
    scale_fill_manual(values = palette_colors[1:num_cell_types], name = "Cell Type") +  # Use dynamic color palette
    geom_point(data = centroids, aes(x = centroid_x, y = centroid_y, fill = cell_type), 
        shape = 21, size = 4, stroke = 1.2, color = "black") +  # Centroid points (correct color)
    scale_color_manual(values = palette_colors[1:num_cell_types], guide = "none") +  # Match colors for centroids
        annotate("text", x = max(df$UMAP_1) * 0.9, y = min(df$UMAP_2) * 1.1, 
            label = paste("Hexbins Used:", optimal_bins), 
            color = "black", size = 5, fontface = "bold", hjust = 1) +  # Display hexbin count
            labs(title = "Hexbin Plot of UMAP Coordinates by Cell Type",
            x = "UMAP Dimension 1", 
            y = "UMAP Dimension 2") +
        theme_minimal() +
        theme(legend.position = "right")
 
# Print the plot
 print(hex_plot)

#Step 5: Calculate Shannon Entropy 

# Ensure df_clean exists and contains UMAP data
if (!exists("df_clean")) stop("df_clean does not exist!")

#Extract UMAP coordinates
x <- df_clean$UMAP_1
y <- df_clean$UMAP_2

#Compute hexbins
num_points <- nrow(df_clean)
optimal_bins <- round(sqrt(num_points / 2))  # Optimal number of bins

hbin <- hexbin(x, y, xbins = optimal_bins)  # Create hexbins

#Compute Shannon Entropy
counts <- hbin@count  # Get the number of cells in each hexbin
total_cells <- sum(counts)  # Total number of cells
probabilities <- counts / total_cells  # Convert counts to probabilities

#Shannon entropy formula: H = -sum(p * log2(p)), ignoring zero probabilities
shannon_entropy <- -sum(probabilities[probabilities > 0] * log2(probabilities[probabilities > 0]))

#Print the result
cat("Shannon Entropy of Hexbin Distribution:", shannon_entropy, "\n")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Function to calculate Shannon entropy
calculate_shannon_entropy <- function(probabilities) {
    -sum(probabilities * log2(probabilities + .Machine$double.xmin))
}

# Extract gene expression data for your first 20 genes (or as needed)
genes_of_interest <- rownames(seurat_obj)[1:20]
expression_data <- FetchData(seurat_obj, vars = genes_of_interest)

# Add cell type information from Seurat object
expression_data$cell_type <- seurat_obj$cell_type

# Reshape the data for entropy calculation
expression_data_long <- melt(expression_data, id.vars = "cell_type", variable.name = "gene", value.name = "expression")

# Ensure expression values are numeric
expression_data_long$expression <- as.numeric(expression_data_long$expression)

# Calculate Shannon entropy for each cell type
entropy_data <- expression_data_long %>%
    group_by(cell_type) %>%
    summarize(
        entropy = calculate_shannon_entropy(table(expression) / sum(table(expression))),
        .groups = "drop"
    )

# Plot Shannon entropy
ggplot(entropy_data, aes(x = cell_type, y = entropy, fill = cell_type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Shannon Entropy of Gene Expression by Cell Type",
         x = "Cell Type",
         y = "Shannon Entropy") +
    scale_fill_viridis_d() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          legend.position = "none")

# Step 6: Plotting Gene Expression

# Check for Expression in the Data Set
head(seurat_obj@assays$RNA@data)

# Define the genes of interest (first 20 genes)
genes_of_interest <- rownames(seurat_obj)[1:20]

# Extract expression data for the genes of interest
expression_data <- FetchData(seurat_obj, vars = genes_of_interest)

# Add cell type information
expression_data$cell_type <- seurat_obj$cell_type

# Reshape the data for ggplot
expression_data_long <- melt(expression_data, id.vars = "cell_type", variable.name = "gene", value.name = "expression")

# Create the violin plot
ggplot(expression_data_long, aes(x = cell_type, y = expression, fill = gene)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(size = 0.1, aes(color = gene), alpha = 0.5) +
  scale_fill_manual(values = viridis(length(genes_of_interest))) +
  scale_color_manual(values = viridis(length(genes_of_interest))) +
  theme_minimal() +
  labs(x = "Cell Type", y = "Expression Levels", title = "Violin Plot of Gene Expression by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_text(text = "Cell Type"))