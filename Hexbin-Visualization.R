#Automatically Optimize Number of Bins and Organize Cell Density into Hexbin with Centroids
#==========================================================================================
# Ensure df_clean contains the correct data
if (!exists("df_clean")) {
    stop("Error: df_clean does not exist. Ensure you have loaded the correct dataset.")
}

# Extract UMAP coordinates and cell types from df_clean
x <- df_clean$UMAP_1
y <- df_clean$UMAP_2
cell_types <- df_clean$cell_type  # Ensure 'cell_type' exists

# Convert to data frame
df <- data.frame(UMAP_1 = x, UMAP_2 = y, cell_type = cell_types)

# Step 1: Calculate optimal bins based on the number of cells
num_points <- nrow(df)
optimal_bins <- round(sqrt(num_points / 2))  # Square root heuristic for bin count

# Step 2: Compute centroid of each cell type
centroids <- df %>%
    group_by(cell_type) %>%
    summarise(
        centroid_x = mean(UMAP_1, na.rm = TRUE),
        centroid_y = mean(UMAP_2, na.rm = TRUE)
    )

# Step 3: Create Hexbin plot with centroids and bin count annotation
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
