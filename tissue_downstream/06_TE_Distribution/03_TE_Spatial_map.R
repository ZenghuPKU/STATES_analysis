# 03_TE_Spatial_map
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(Seurat)
library(anndata)
library(reticulate)
library(ggplot2)
library(viridis)
library(circlize)
rm(list = ls())
setwd("~/tissue_downstream/06_TE_Distribution/")
use_condaenv("scanpy_env", required = TRUE)
py_config()

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Construct Seurat Objects and Calculate TE
totalRNA_matrix <- t(sc_ad$layers['totalRNA_raw'])
rbRNA_matrix <- t(sc_ad$layers['rbRNA'])
metadata <- sc_ad$obs
totalRNA <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata)
rbRNA <- CreateSeuratObject(counts = rbRNA_matrix, meta.data = metadata)
common_cells <- intersect(rownames(totalRNA@meta.data), rownames(rbRNA@meta.data))
totalRNA <- subset(totalRNA, cells = common_cells)
rbRNA <- subset(rbRNA, cells = common_cells)
totalRNA_expr <- as.matrix(GetAssayData(totalRNA, slot = "counts"))
rbRNA_expr <- as.matrix(GetAssayData(rbRNA, slot = "counts"))
common_genes <- intersect(rownames(totalRNA_expr), rownames(rbRNA_expr))
totalRNA_expr <- totalRNA_expr[common_genes, ]
rbRNA_expr <- rbRNA_expr[common_genes, ]
totalRNA_total <- colSums(totalRNA_expr)
rbRNA_total <- colSums(rbRNA_expr)
translation_efficiency <- rbRNA_total / totalRNA_total
translation_efficiency[is.na(translation_efficiency) | is.infinite(translation_efficiency)] <- 0

# STATES-C1 Subset, Coordinates & Add TE
states_sub <- subset(states, subset = `protocol.replicate` == "STATES-C1")
common_cells_c1 <- intersect(rownames(states_sub@meta.data), rownames(sc_ad$obs))
states_sub <- subset(states_sub, cells = common_cells_c1)

states_sub@meta.data$spatial_x <- as.numeric(as.character(sc_ad$obs[common_cells_c1, "column"]))
states_sub@meta.data$spatial_y <- as.numeric(as.character(sc_ad$obs[common_cells_c1, "row"]))
states_sub@meta.data$spatial_x <- max(states_sub@meta.data$spatial_x, na.rm = TRUE) - states_sub@meta.data$spatial_x
states_sub@meta.data$Translation_Efficiency <- translation_efficiency[rownames(states_sub@meta.data)]

# Prepare Plotting Data
plot_data_te <- states_sub@meta.data[!is.na(states_sub@meta.data$spatial_x) &
                                     !is.na(states_sub@meta.data$spatial_y) &
                                     !is.na(states_sub@meta.data$Translation_Efficiency), ]

# Original TE Spatial Map (viridis)
# Obtain inferno color palette
inferno_colors <- viridis::inferno(256)

# Calculate the 10% and 90% quantiles of TE
te_range <- quantile(plot_data_te$Translation_Efficiency, probs = c(0.1, 0.9), na.rm = TRUE)

# Use inferno color palette with reversed direction (high values dark)
p1 <- ggplot(plot_data_te, aes(x = spatial_x, y = spatial_y, color = Translation_Efficiency)) +
  geom_point(size = 0.25) +
  scale_color_viridis(
    option = "C",
    direction = -1,
    limits = te_range,
    oob = scales::squish,
    name = "TE"
  ) +
  coord_fixed() +
  theme_classic() +
  labs(title = "Original TE", x = "Column", y = "Row") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p1

# Save result
ggsave("TE_Spatial_map.pdf", p1, width = 6, height = 5)




