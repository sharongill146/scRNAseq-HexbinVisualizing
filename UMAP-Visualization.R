#Visualize the joint distribution - Generating a 2D density plot
#===============================================================
# Ensure UMAP coordinates are extracted correctly
if (exists("umap_coords") && !is.null(umap_coords)) {
  
  # Convert matrix to data frame
  umap_df <- as.data.frame(umap_coords)
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")  # Ensure proper column names
  
  # Generate a 2D density plot
  p_1 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    coord_equal() +  # For equal scaling on both axes
    geom_density_2d_filled(alpha = 0.6) +  # Filled 2D density plot with transparency
    geom_density_2d(color = "black", linewidth = 0.2) +  # Outline for better visibility
    labs(title = "UMAP 2D Density", x = "UMAP 1", y = "UMAP 2")
  
  # Print the plot
  print(p_1)
  
} else {
  message("UMAP coordinates not found or not properly extracted.")
}

#Highlight High-Density Regions for immediately interpresetable regions
p_1 + geom_hdr()

# Visualizing Subpopulations
#============================
# Ensure UMAP coordinates are extracted correctly
if (exists("umap_coords") && !is.null(umap_coords)) {
  
  # Convert UMAP coordinates matrix to a data frame
  df_clean <- as.data.frame(umap_coords)
  colnames(df_clean) <- c("UMAP_1", "UMAP_2")  # Assign correct column names

  # Generate a distinct color palette using RColorBrewer
  library(RColorBrewer)
  num_classes <- length(unique(df_clean$General_type))

  # Use Set3 or other palettes that support more colors
colors <- brewer.pal(max(min(num_classes, 12), 3), "Set3")

  # Add metadata information (e.g., cell types) from Seurat object
  df_clean$cell_type <- seurat_obj@meta.data$cell_type  # Adjust column name if needed

  # Remove NA values in UMAP coordinates
  df_clean <- df_clean[!is.na(df_clean$UMAP_1) & !is.na(df_clean$UMAP_2), ]
  
  # Replace remaining NA values (if any) with column means
  df_clean$UMAP_1[is.na(df_clean$UMAP_1)] <- mean(df_clean$UMAP_1, na.rm = TRUE)
  df_clean$UMAP_2[is.na(df_clean$UMAP_2)] <- mean(df_clean$UMAP_2, na.rm = TRUE)

  # Handle missing values in 'cell_type'
  if ("cell_type" %in% colnames(df_clean)) {
    if (any(is.na(df_clean$cell_type))) {
      df_clean$cell_type[is.na(df_clean$cell_type)] <- "Unknown"
    }
  } else {
    message("Column 'cell_type' not found in metadata.")
  }

  # Adjust xlim and ylim based on data distribution
  x_range <- range(df_clean$UMAP_1, na.rm = TRUE)
  y_range <- range(df_clean$UMAP_2, na.rm = TRUE)
  
  # Example adjusted limits
  x_lim <- c(x_range[1] + 0.1 * diff(x_range), x_range[2] - 0.1 * diff(x_range))
  y_lim <- c(y_range[1] + 0.1 * diff(y_range), y_range[2] - 0.1 * diff(y_range))
  
  # Generate the visualization
  umap_plot <- ggplot(df_clean, aes(x = UMAP_1, y = UMAP_2, fill = cell_type)) +
    geom_hdr(xlim = x_range, ylim = y_range) +  # High-density region visualization
    geom_point(shape = 21, alpha = 0.6) +  # Semi-transparent points for better visualization
    scale_fill_brewer(palette = "Set3") +  # Better color mapping
    labs(title = "UMAP Subpopulation Visualization",
         x = "UMAP Dimension 1", y = "UMAP Dimension 2", fill = "Cell Type") +
    theme_minimal() +
    theme(legend.position = "right", axis.text = element_text(size = 12))
  
  # Print the plot
  print(umap_plot)

} else {
  message("UMAP coordinates not found or not properly extracted.")
}

#UMAP Embedding and Density Visualization with Cell-Type Distribution
#====================================================================
  ggplot(df_clean, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_hdr_lines(xlim = x_lim, ylim = y_lim) +  # High-density region contours
  geom_hdr(alpha = 0.3, xlim = x_lim, ylim = y_lim) +  # Highlight density, reduce alpha for clarity
  geom_point(size = 1, alpha = 0.8, shape = 21, fill = "white", color = "black", stroke = 0.3) +  # Plot individual points, outlined points for contrast
  geom_jitter(size = 0.8, width = 0.2, height = 0.2, alpha = 0.7) +  # More jitter to separate overlapping points
  facet_wrap(~cell_type) +  # Facet by cell type
  scale_color_viridis_d(option = "plasma") +  
  theme_minimal() +
  theme(
    strip.text = element_blank(),  
    panel.spacing = unit(2, "lines"),  
    plot.margin = margin(10, 10, 10, 10), 
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.ticks.x = element_line()  
  )