# 02_TE_Distribution_label3
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(Seurat)
library(anndata)
library(reticulate)
library(ggplot2)
library(ggridges)
library(viridis)
rm(list = ls())
setwd("~/tissue_downstream/06_TE_Distribution/")
use_condaenv("scanpy_env", required = TRUE)
py_config()

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Define colors for cell types
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

# Process expression matrices
totalRNA_matrix <- t(sc_ad$layers['totalRNA_raw'])
rbRNA_matrix <- t(sc_ad$layers['rbRNA'])
metadata <- sc_ad$obs
totalRNA <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata)
rbRNA <- CreateSeuratObject(counts = rbRNA_matrix, meta.data = metadata)

# Obtain intersection of cells
common_cells <- intersect(rownames(totalRNA@meta.data), rownames(rbRNA@meta.data))
totalRNA <- subset(totalRNA, cells = common_cells)
rbRNA <- subset(rbRNA, cells = common_cells)

# Extract label3 annotations from states
if (!exists("states")) stop("Error: states object not properly loaded.")
if (!"states_nn_alg1_label3" %in% colnames(states@meta.data)) stop("Error: label3 annotations missing.")
states_metadata <- states@meta.data[match(common_cells, rownames(states@meta.data)), ]
cell_labels_label3 <- states_metadata$states_nn_alg1_label3

# Retrieve expression matrices & common genes
totalRNA_expr <- as.matrix(GetAssayData(totalRNA, slot = "counts"))
rbRNA_expr <- as.matrix(GetAssayData(rbRNA, slot = "counts"))
common_genes <- intersect(rownames(totalRNA_expr), rownames(rbRNA_expr))
totalRNA_expr <- totalRNA_expr[common_genes, ]
rbRNA_expr <- rbRNA_expr[common_genes, ]

# Calculate TE (cell level)
totalRNA_total_counts <- colSums(totalRNA_expr)
rbRNA_total_counts <- colSums(rbRNA_expr)
translation_efficiency <- rbRNA_total_counts / totalRNA_total_counts
translation_efficiency[is.na(translation_efficiency) | is.infinite(translation_efficiency)] <- 0

# Create TE data frame (cell level)
cor_df_label3 <- data.frame(
  Cell = common_cells,
  CellType = cell_labels_label3,
  Translation_Efficiency = translation_efficiency
)
cor_df_label3 <- subset(cor_df_label3, CellType != "Mix")
cor_df_label3 <- subset(cor_df_label3, CellType != "TEGLU Mix")

# Sort
celltype_means_label3 <- cor_df_label3 %>%
  group_by(CellType) %>%
  summarise(mean_TE = mean(Translation_Efficiency)) %>%
  arrange(mean_TE)
celltype_levels_label3 <- celltype_means_label3$CellType

# Plot: TE density plot (cell level)
pdf("translation_efficiency_density_by_celltype_label3.pdf", width = 9, height = 14)
ggplot(cor_df_label3, aes(x = Translation_Efficiency, y = factor(CellType, levels = celltype_levels_label3), fill = CellType)) +
  geom_density_ridges(scale = 1.8, alpha = 0.8) + 
  scale_y_discrete(expand = expansion(mult = c(0, 0.4))) +  
  scale_x_continuous(
    limits = c(min(cor_df_label3$Translation_Efficiency), max(cor_df_label3$Translation_Efficiency)),
    breaks = seq(0, max(cor_df_label3$Translation_Efficiency), by = 0.05)
  ) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  labs(title = "TE Distribution Each Cell(Label3)", 
       x = "TE_cell_mean", 
       y = "Cell Type") +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.2, "cm")
  )
dev.off()
