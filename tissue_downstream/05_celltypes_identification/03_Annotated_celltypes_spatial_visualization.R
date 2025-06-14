# 03_Annotated_celltypes_spatial_visualization
# Load libraries and set environment
library(Seurat)
library(anndata)
library(reticulate)
library(dplyr)
library(matrixStats) 
library(ggplot2)
rm(list = ls())
setwd("~/tissue_downstream/05_celltypes_identification/")
use_condaenv("scanpy_env", required = TRUE)
py_config()

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Select the brain slices that you wish to display, taking STATES-C1 as an example
states_sub <- subset(states, subset = `protocol.replicate` == "STATES-C1")

# Check if the rownames match
print(head(rownames(states_sub@meta.data)))
print(head(rownames(sc_ad$obs)))

# Intersect cell names to ensure consistent matching.
common_cells <- intersect(rownames(states_sub@meta.data), rownames(sc_ad$obs))
states_sub <- subset(states_sub, cells = common_cells)

# Add spatial coordinates
states_sub@meta.data$spatial_x <- as.numeric(as.character(sc_ad$obs[common_cells, "column"]))
states_sub@meta.data$spatial_y <- as.numeric(as.character(sc_ad$obs[common_cells, "row"]))
states_sub@meta.data$spatial_x <- max(states_sub@meta.data$spatial_x, na.rm = TRUE) - states_sub@meta.data$spatial_x

# Filter NA values in spatial coordinates to ensure ggplot compatibility.
plot_data <- states_sub@meta.data[!is.na(states_sub@meta.data$spatial_x) & !is.na(states_sub@meta.data$spatial_y), ]

# label2 spatial visualization
celltype_colors <- c(
  "OPC" = "#667872",     
  "DE/MEN" = "#b274e8",  
  "CHO/PEP" = "#97f4f7",  
  "OLG" = "#e6db17",      
  "CHOR/EPEN" = "#F1B2E9", 
  "TEPN" = "#6C9F35",     
  "MLG" = "#8597c6",     
  "VAS" = "#17BED0",      
  "INH" = "#3182BD",     
  "AC" = "#FDC06F",
  "Mix"="#F5F5F5"
)
p_spatial <- ggplot(plot_data, aes(x = spatial_x, y = spatial_y, color = states_nn_alg1_label2)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = celltype_colors) +
  theme_classic() +
  labs(x = "Column", y = "Row")+ NoLegend()
p_spatial
ggsave("spatial_STATESC1_label2.pdf", p_spatial, width = 8, height = 7, device = cairo_pdf)

# label3 spatial visualization
celltype_colors <- c(
  "AC1" = "#ccba33",      
  "AC2" = "#ffbe85",   
  "AC3" = "#e3782b",  
  "CHOR" = "#7f52a9",     
  "EPEN" = "#c4b0d4",    
  "CHO/PEP" = "#97f4f7", 
  "INH_Sst" = "#96abeb",      
  "INH_Pvalb" = "#96cad4", 
  "INH_Cnr1_Vip" = "#a8e1eb",
  "MLG" = "#8597c6",     
  "OPC" = "#667872", 
  "OLG1" = "#e4f768",      
  "OLG2" = "#e6db17",   
  "VLMC" = "#1f76b3" ,      
  "VSMC" = "#00aeef" ,
  "Peri/VEC" = "#d3a59c",
  "DE/MEN" = "#b274e8", 
  "MSN" = "#D96DA1" ,
  "DGGRC" = "#a6e8a6",  
  "TEGLU CA1" = "#77ed8f",      
  "TEGLU CA2" = "#82ad2d",      
  "TEGLU CA3" = "#28330b",      
  "TEGLU L2/3" = "#cbfc60",
  "TEGLU L4" = "#96db00",
  "TEGLU L5" = "#04b361",
  "TEGLU L5/6" = "#40d102", 
  "TEGLU L6" = "#32a630",
  "TEGLU L6b" = "#406e27",
  "TEGLU Mix" = "#c5fcc5",      
  "Mix"="#F5F5F5"
)
p_spatial <- ggplot(plot_data, aes(x = spatial_x, y = spatial_y, color = states_nn_alg1_label3)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = celltype_colors) +
  theme_classic() +
  labs(x = "Column", y = "Row")+ NoLegend()
p_spatial
ggsave("spatial_STATESC1_label3.pdf", p_spatial, width = 8, height = 7, device = cairo_pdf)

# Select specific cell types for visualization.
# Define color mapping by assigning colors to all cell types using a hue palette.
levels_labels <- levels(plot_data$states_nn_alg1_label3)

# Initialize all cell type colors to gray.
default_palette <- setNames(rep("lightgrey", length(levels_labels)), levels_labels)

# Specify cell types for visualization.
# TEPN
default_palette["TEGLU L2/3"] <- "#cbfc60"#
default_palette["TEGLU L4"] <- "#96db00"#                                         
default_palette["TEGLU L5"] <- "#04b361"
default_palette["TEGLU L5/6"] <- "#40d102"
default_palette["TEGLU L6"] <- "#32a630"
default_palette["TEGLU L6b"] <- "#406e27"
default_palette["TEGLU Mix"] <- "#c5fcc5"
default_palette["TEGLU CA1"] <- "#77ed8f"
default_palette["TEGLU CA2"] <- "#82ad2d"
default_palette["TEGLU CA3"] <- "#28330b"
default_palette["DGGRC"] <- "#a6e8a6"
default_palette["MSN"] <- "#D96DA1"

# AC
default_palette["AC1"] <- "#ccba33"
default_palette["AC2"] <- "#ffbe85"
default_palette["AC3"] <- "#e3782b"

# CHOR/EPEN
default_palette["CHOR"] <- "#7f52a9"
default_palette["EPEN"] <- "#c4b0d4"

# CHO/PEP
default_palette["CHO/PEP"] <- "#97f4f7"

# INH
default_palette["INH_Sst"] <- "#96abeb"
default_palette["INH_Pvalb"] <- "#96cad4"
default_palette["INH_Cnr1_Vip"] <- "#a8e1eb"

# MLG
default_palette["MLG"] <- "#8597c6"

# OPC
default_palette["OPC"] <- "#667872"

# OLG
default_palette["OLG1"] <- "#e4f768"
default_palette["OLG2"] <- "#e6db17"

# VAS
default_palette["VLMC"] <- "#1f76b3"
default_palette["VSMC"] <- "#00aeef"
default_palette["Peri/VEC"] <- "#d3a59c"

# DE/MEN
default_palette["DE/MEN"] <- "#b274e8"

# Visualize the selected cell types.
p_spatial <- ggplot(plot_data, aes(x = spatial_x, y = spatial_y, color = states_nn_alg1_label3)) +
  geom_point(size = 0.25) +
  scale_color_manual(values = default_palette) +
  theme_classic() +
  labs(x = "Column", y = "Row", title = "Selected Cell Types Spatial Distribution for STATES-C1")+ NoLegend()
p_spatial
ggsave("spatial_STATESC1_selected_celltypes.pdf", p_spatial, width = 5.5, height = 5, device = cairo_pdf)


