#Create Hexagonal Map Visualization for Single-cell Data
#'
#' @param seurat_obj Seurat object
#' @param reduction Dimension reduction to use (default: "umap")
#' @param label_col Column name from meta.data for cell labels
#' @param hex_density Controls hex bin size (default: 50)
#' @param plot_type Type of plot: "majority", "entropy", or "gene"
#' @param gene_name Gene name for gene expression plot
#' @param assay Assay to use for gene expression (default: "SCT")
#' @param layer Data layer for gene expression (default: "scale.data")
#' @param min_cells Minimum cells per hex (default: 3)
#' @param color_palette ggsci palette: "npg","aaas","nejm","jama","jco","lancet","d3" (default: "npg")
#' @param alpha Transparency level (default: 0.3)
#' @param viridis Viridis palette for entropy plot (default: "magma")
#'
#' @return A ggplot object
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import ggsci
#'
#' @export

HexnicPlot <- function(
    seurat_obj,
    reduction = "umap",
    label_col,
    hex_density = 62,
    plot_type = c("majority", "entropy", "gene"),
    gene_name = NULL,
    assay = "SCT",
    layer = "scale.data",
    min_cells = 3,
    color_palette = "npg",
    alpha = 0.3,
    viridis = "magma"
) {
  plot_type <- match.arg(plot_type)
  if (plot_type == "gene" && is.null(gene_name)) {
    stop("gene_name must be provided when plot_type is 'gene'")
  }

  get_color_scale <- function(palette_name, continuous = FALSE) {
    adaptive_pal_inner <- function(values) {
      force(values)
      n_colors <- length(values)
      function(n) {
        if (n <= n_colors) {
          values[seq_len(n)]
        } else {
          colorRampPalette(values, alpha = TRUE)(n)
        }
      }
    }

    raw_cols <- switch(palette_name,
                       "npg" = ggsci:::ggsci_db$"npg"[["nrc"]],
                       "aaas" = ggsci:::ggsci_db$"aaas"[["default"]],
                       "nejm" = ggsci:::ggsci_db$"nejm"[["default"]],
                       "jama" = ggsci:::ggsci_db$"jama"[["default"]],
                       "jco" = ggsci:::ggsci_db$"jco"[["default"]],
                       "lancet" = ggsci:::ggsci_db$"lancet"[["lanonc"]],
                       "d3" = ggsci:::ggsci_db$"d3"[["category20"]]
    )

    raw_cols_rgb <- col2rgb(raw_cols)
    alpha_cols <- rgb(
      raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
      alpha = 255L, names = names(raw_cols),
      maxColorValue = 255L
    )

    if (continuous) {
      scale_fill_gradientn(colors = raw_cols)
    } else {
      discrete_scale("fill", palette_name, adaptive_pal_inner(unname(alpha_cols)))
    }
  }

  umap_coords <- seurat_obj@reductions[[reduction]]@cell.embeddings

  hb_umap_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    class = seurat_obj@meta.data[[label_col]]
  )

  x_range <- diff(range(hb_umap_df$UMAP1))
  y_range <- diff(range(hb_umap_df$UMAP2))
  hex_size <- min(x_range, y_range) / hex_density

  hex_coords <- function(x, y, size) {
    sqrt3 <- sqrt(3)
    col <- round(x / (size * 1.5))
    row <- round((y - (col %% 2) * size * sqrt3/2) / (size * sqrt3))
    list(
      x = col * size * 1.5,
      y = row * size * sqrt3 + (col %% 2) * size * sqrt3/2
    )
  }

  coords <- hex_coords(hb_umap_df$UMAP1, hb_umap_df$UMAP2, hex_size)
  hb_umap_df$hex_x <- coords$x
  hb_umap_df$hex_y <- coords$y

  if (plot_type == "gene") {
    gene_expr <- GetAssayData(seurat_obj, assay = assay, layer = layer)[gene_name,]
    hb_umap_df$gene_expr <- gene_expr[rownames(hb_umap_df)]
  }

  hex_summary <- switch(
    plot_type,
    "majority" = {
      hb_umap_df %>%
        group_by(hex_x, hex_y) %>%
        summarize(
          cell_count = n(),
          majority = names(sort(table(class), decreasing = TRUE))[1],
          prop_majority = max(table(class)) / n(),
          .groups = "drop"
        )
    },
    "entropy" = {
      hb_umap_df %>%
        group_by(hex_x, hex_y) %>%
        summarize(
          cell_count = n(),
          majority = names(sort(table(class), decreasing = TRUE))[1],
          entropy = {
            props <- table(class)/n()
            -sum(props * log2(props + .Machine$double.xmin))
          },
          .groups = "drop"
        )
    },
    "gene" = {
      hb_umap_df %>%
        group_by(hex_x, hex_y) %>%
        summarize(
          cell_count = n(),
          majority = names(sort(table(class), decreasing = TRUE))[1],
          mean_expr = mean(gene_expr, na.rm = TRUE),
          .groups = "drop"
        )
    }
  ) %>%
    filter(cell_count >= min_cells)

  if (nrow(hex_summary) < 2) {
    stop("Not enough rows in hex_summary to perform interpolation.")
  }

  p <- switch(
    plot_type,
    "majority" = {
      ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = majority, alpha = prop_majority), stat = "identity") +
        get_color_scale(color_palette) +
        scale_alpha_continuous(range = c(alpha, 1)) +
        labs(alpha = "Proportion\nMajority", fill = "Cell Type")
    },
    "entropy" = {
      ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = entropy), stat = "identity") +
        scale_fill_viridis_c(option = viridis) +
        labs(fill = "Shannon\nEntropy")
    },
    "gene" = {
      ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = mean_expr), stat = "identity") +
        scale_fill_viridis_c(option = viridis) +  # Use magma from viridis for color scale
        scale_alpha_continuous(range = c(alpha, 1)) +
        labs(fill = paste0(gene_name, "\nExpression"))
    }
  ) +
    coord_equal() +
    theme_minimal() +
    theme(panel.grid = element_blank())

  return(p)
}

# Majority plot: Identifies the dominant cell type in each hexbin
# ================================================================
HexnicPlot(seurat_obj, label_col = "cell_type", plot_type = "majority")

# Shannon entropy plot: Measures the diversity of cell types in each hexbin
# =========================================================================
HexnicPlot(seurat_obj, label_col = "cell_type", plot_type = "entropy")

# Gene expression plot: Visualizes the expression of a specific gene 
# ===================================================================
gene_name <- "ENSG00000129757"

# Create the gene expression plot
p <- HexnicPlot(seurat_obj, label_col = "cell_type", plot_type = "gene", gene_name = gene_name, assay = "RNA", layer = "counts")

# Use the magma color palette for the gene expression plot
p + scale_fill_viridis_c(option = "magma")


#A Second Hexnic Plot Function Specific to Calculating Methylation
#==================================================================
HexnicPlotMethylation <- function(
    seurat_obj,
    reduction = "umap",
    label_col,
    hex_density = 62,
    assay = "RNA",
    layer = "data",
    min_cells = 3,
    alpha = 0.8,
    viridis = "magma"
) {
    umap_coords <- seurat_obj@reductions[[reduction]]@cell.embeddings
    
    hb_umap_df <- data.frame(
        UMAP1 = umap_coords[, 1],
        UMAP2 = umap_coords[, 2],
        class = seurat_obj@meta.data[[label_col]]
    )
    
    # Calculate hex bin size
    x_range <- diff(range(hb_umap_df$UMAP1))
    y_range <- diff(range(hb_umap_df$UMAP2))
    hex_size <- min(x_range, y_range) / hex_density
    
    # Assign hexagon bins
    hb_umap_df <- hb_umap_df %>%
        mutate(hex_x = round(UMAP1 / hex_size) * hex_size,
               hex_y = round(UMAP2 / hex_size) * hex_size)
    
    # Compute methylation score (mean expression of first 20 genes)
    gene_list <- rownames(GetAssayData(seurat_obj, assay = assay, layer = layer))[1:20]
    methylation_score <- colMeans(GetAssayData(seurat_obj, assay = assay, layer = layer)[gene_list, ], na.rm = TRUE)
    hb_umap_df$methylation_score <- methylation_score[rownames(hb_umap_df)]
    
    # Compute summary statistics per hex
    hex_summary <- hb_umap_df %>%
        group_by(hex_x, hex_y) %>%
        summarize(
            cell_count = n(),
            majority = names(sort(table(class), decreasing = TRUE))[1],  # Most common cell type
            mean_methylation = mean(methylation_score, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        filter(cell_count >= min_cells)
    
    # Compute cluster centroids for labels
    cluster_labels <- hb_umap_df %>%
        group_by(class) %>%
        summarize(
            centroid_x = median(UMAP1),
            centroid_y = median(UMAP2),
            .groups = "drop"
        )
    
    # Plot with light colors and dark grey background
    p <- ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = mean_methylation), stat = "identity") +
        scale_fill_viridis_c(option = viridis, begin = 0, end = 0.7) +  # Use light colors by controlling the color range
        geom_text(data = cluster_labels, aes(x = centroid_x, y = centroid_y, label = class), 
                  color = "white", size = 1.5, fontface = "bold") +  # Change label color to white
        labs(fill = "Methylation Score") +
        coord_equal() +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            plot.background = element_rect(fill = "grey60"),  # Dark grey background
            panel.background = element_rect(fill = "grey60"),  # Dark grey panel background
            axis.text = element_text(color = "white"),         # White axis text
            axis.title = element_text(color = "white"),        # White axis titles
            legend.background = element_rect(fill = "grey60"), # Dark grey legend background
            legend.title = element_text(color = "white"),      # White legend title
            legend.text = element_text(color = "white")        # White legend text
        )
    
    print(p)
}

HexnicPlotMethylation(seurat_obj, label_col = "cell_type", viridis = "magma")

#Solving for most methylated gene in a certain cell type
#========================================================
FindMostMethylatedGene <- function(seurat_obj, cell_type, assay = "RNA", layer = "data") {
    # Check if the cell_type column exists
    if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
        stop("Cell type information is missing from metadata.")
    }
    
    # Filter cells based on the specified cell type
    filtered_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$cell_type == cell_type, ])
    
    # Check if there are any filtered cells
    if (length(filtered_cells) == 0) {
        stop("No cells match the specified cell type.")
    }
    
    # Extract methylation data for the filtered cells
    methylation_data <- GetAssayData(seurat_obj, assay = assay, slot = layer)
    
    # Calculate the mean methylation for each gene across the filtered cells
    mean_methylation <- apply(methylation_data[, filtered_cells], 1, mean, na.rm = TRUE)
    
    # Find the gene with the highest mean methylation
    most_methylated_gene <- names(mean_methylation)[which.max(mean_methylation)]
    
    return(most_methylated_gene)
}

# Example of using the function for "CD14-low, CD16-positive monocyte" cell type
most_methylated_gene <- FindMostMethylatedGene(seurat_obj, cell_type = "CD14-low, CD16-positive monocyte")
print(most_methylated_gene)