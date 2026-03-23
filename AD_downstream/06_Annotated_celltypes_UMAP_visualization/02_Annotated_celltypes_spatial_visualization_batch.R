# 02_Annotated_celltypes_spatial_visualization_batch
rm(list = ls())
library(Seurat)
library(anndata)
library(reticulate)
library(dplyr)
library(matrixStats) 
library(ggplot2)
library(viridis)
setwd("~/AD_downstream/06_Annotated_celltypes_UMAP_visualization/")

states <- readRDS("states_with_plaque_info.rds")
sc_ad  <- read_h5ad("ADcombined_mousebrain.h5ad")

type_levels <- sort(unique(states$type))
print(type_levels)

celltype_colors_l2 <- c(
  "OPC"       = "#667872",     
  "DE/MEN"    = "#b274e8",  
  "CHO/PEP"   = "#97f4f7",  
  "OLG"       = "#e6db17",      
  "CHOR/EPEN" = "#F1B2E9", 
  "TEPN"      = "#6C9F35",     
  "MLG"       = "#8597c6",     
  "VAS"       = "#17BED0",      
  "INH"       = "#3182BD",     
  "AC"        = "#FDC06F",
  "Mix"       = "#F5F5F5"
)

celltype_colors_l3 <- c(
  "AC1" = "#ccba33",     
  "AC2" = "#ffbe85",   
  "AC3" = "#e3782b",  
  "AC4" = "#f0dea3",
  "CHOR" = "#7f52a9",     
  "EPEN" = "#c4b0d4",    
  "CHO/PEP" = "#97f4f7", 
  "INH_Sst" = "#96abeb",      
  "INH_Pvalb" = "#96cad4", 
  "INH_Cnr1_Vip" = "#a8e1eb",
  "MLG1" = "#A9B4D8",
  "MLG2" = "#8597C6",
  "MLG3" = "#6274A3",  
  "OPC" = "#667872", 
  "OLG1" = "#e4f768",      
  "OLG2" = "#e6db17",   
  "VLMC" = "#1f76b3",
  "VSMC" = "#00aeef",
  "Peri/VEC" = "#d3a59c",
  "DE/MEN" = "#b274e8", 
  "MSN" = "#D96DA1",
  "DGGRC" = "#a6e8a6",  
  "TEGLU CA1"  = "#77ed8f",      
  "TEGLU CA2"  = "#82ad2d",      
  "TEGLU CA3"  = "#28330b",      
  "TEGLU L2/3" = "#cbfc60",
  "TEGLU L4"   = "#96db00",
  "TEGLU L4/5" = "#04b361",
  "TEGLU L5"   = "#40d102",
  "TEGLU L6"   = "#32a630",
  "TEGLU Mix"  = "#c5fcc5",      
  "Mix"        = "#F5F5F5"
)

for (tp in type_levels) {
  message("Processing type: ", tp)
  
  states_sub <- subset(states, subset = type == tp)
  if (ncol(states_sub) == 0) {
    message("  No cells for type ", tp, ", skip.")
    next
  }
  
  common_cells <- intersect(rownames(states_sub@meta.data), rownames(sc_ad$obs))
  if (length(common_cells) < 10) {
    message("  Too few cells after intersect for type ", tp, ", skip.")
    next
  }
  states_sub <- subset(states_sub, cells = common_cells)

  states_sub@meta.data$spatial_x <- as.numeric(as.character(sc_ad$obs[common_cells, "column"]))
  states_sub@meta.data$spatial_y <- as.numeric(as.character(sc_ad$obs[common_cells, "row"]))
  
  plot_data <- states_sub@meta.data[
    !is.na(states_sub@meta.data$spatial_x) &
      !is.na(states_sub@meta.data$spatial_y),
  ]
  
  if (nrow(plot_data) == 0) {
    message("  No valid spatial coords for type ", tp, ", skip.")
    next
  }

  plot_data$states_nn_alg1_label2 <- factor(plot_data$states_nn_alg1_label2)
  plot_data$states_nn_alg1_label3 <- factor(plot_data$states_nn_alg1_label3)
  
  lev_l2  <- levels(plot_data$states_nn_alg1_label2)
  miss_l2 <- setdiff(lev_l2, names(celltype_colors_l2))
  
  if (length(miss_l2) > 0) {
    extra_cols_l2 <- viridis::viridis(length(miss_l2))
    names(extra_cols_l2) <- miss_l2
    celltype_colors_l2_all <- c(celltype_colors_l2, extra_cols_l2)
  } else {
    celltype_colors_l2_all <- celltype_colors_l2
  }
  
  p_l2 <- ggplot(
    plot_data,
    aes(x = spatial_x, y = spatial_y, color = states_nn_alg1_label2)
  ) +
    geom_point(size = 0.5) +
    scale_color_manual(values = celltype_colors_l2_all) +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_classic() +
    labs(
      x = "Column",
      y = "Row",
      title = paste0("Label2 spatial (type = ", tp, ")")
    ) +
    NoLegend()
  
  ggsave(
    filename = paste0("spatial_label2_", tp, ".pdf"),
    plot     = p_l2,
    width    = 8,
    height   = 7,
    device   = cairo_pdf
  )
  
  lev_l3  <- levels(plot_data$states_nn_alg1_label3)
  miss_l3 <- setdiff(lev_l3, names(celltype_colors_l3))
  
  if (length(miss_l3) > 0) {
    extra_cols_l3 <- viridis::viridis(length(miss_l3))
    names(extra_cols_l3) <- miss_l3
    celltype_colors_l3_all <- c(celltype_colors_l3, extra_cols_l3)
  } else {
    celltype_colors_l3_all <- celltype_colors_l3
  }
  
  p_l3 <- ggplot(
    plot_data,
    aes(x = spatial_x, y = spatial_y, color = states_nn_alg1_label3)
  ) +
    geom_point(size = 0.5) +
    scale_color_manual(values = celltype_colors_l3_all) +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_classic() +
    labs(
      x = "Column",
      y = "Row",
      title = paste0("Label3 spatial (type = ", tp, ")")
    ) +
    NoLegend()
  
  ggsave(
    filename = paste0("spatial_label3_", tp, ".pdf"),
    plot     = p_l3,
    width    = 8,
    height   = 7,
    device   = cairo_pdf
  )
}

message("Done: saved whole-sample label2 and label3 spatial plots for each type.")