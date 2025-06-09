# 01_TE_Distribution_label2
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(Seurat)
library(anndata)
library(ggridges)
library(viridis)
library(ggplot2)
rm(list = ls())
setwd("~/tissue_downstream/06_TE_Distribution/")

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Define colors for cell types
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

# Retrieve matrix data
totalRNA_matrix <- t(sc_ad$layers['totalRNA'])
rbRNA_matrix <- t(sc_ad$layers['rbRNA'])

# Load metadata
metadata <- sc_ad$obs

# Create Seurat objects
totalRNA <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata)
rbRNA <- CreateSeuratObject(counts = rbRNA_matrix, meta.data = metadata)

# Extract cell type information from states object
if (!exists("states")) stop("Error: states object not properly loaded. Please check the RData file.")
if (!"states_nn_alg1_label2" %in% colnames(states@meta.data)) stop("Error: 'states_nn_alg1_label2' column not found in states@meta.data.")
cell_labels <- states@meta.data$states_nn_alg1_label2

# Ensure consistent cell names (take intersection)
common_cells <- intersect(rownames(totalRNA@meta.data), rownames(rbRNA@meta.data))
totalRNA <- subset(totalRNA, cells = common_cells)
rbRNA <- subset(rbRNA, cells = common_cells)

# Extract cell types from states metadata
states_metadata <- states@meta.data[match(common_cells, rownames(states@meta.data)), ]
cell_labels <- states_metadata$states_nn_alg1_label2

# Obtain expression matrices
totalRNA_expr <- as.matrix(GetAssayData(totalRNA, slot = "counts"))
rbRNA_expr <- as.matrix(GetAssayData(rbRNA, slot = "counts"))

# Ensure matching gene names
common_genes <- intersect(rownames(totalRNA_expr), rownames(rbRNA_expr))
totalRNA_expr <- totalRNA_expr[common_genes, ]
rbRNA_expr <- rbRNA_expr[common_genes, ]

# Calculate total gene expression sums for totalRNA and rbRNA per cell
totalRNA_total_counts <- colSums(totalRNA_expr)
rbRNA_total_counts <- colSums(rbRNA_expr)

# Calculate translation efficiency (TE = rbRNA / totalRNA)
translation_efficiency <- rbRNA_total_counts / totalRNA_total_counts
translation_efficiency[is.na(translation_efficiency) | is.infinite(translation_efficiency)] <- 0

# Create a data frame containing translation efficiency and cell type for each cell
cor_df <- data.frame(
  Cell = common_cells,
  CellType = cell_labels,
  Translation_Efficiency = translation_efficiency
)

# Remove cells with "Mix" cell type
cor_df <- subset(cor_df, CellType != "Mix")

# Calculate mean TE for each cell type for sorting
celltype_means <- cor_df %>%
  group_by(CellType) %>%
  summarise(mean_TE = mean(Translation_Efficiency)) %>%
  arrange(mean_TE)

# Set cell type order (sorted by mean TE in descending order)
celltype_levels <- celltype_means$CellType

# Generate density plot (TE at the cell level)
pdf("translation_efficiency_density_by_celltype_label2.pdf", width = 8, height = 10)
ggplot(cor_df, aes(x = Translation_Efficiency, y = factor(CellType, levels = celltype_levels), fill = CellType)) +
  geom_density_ridges(scale = 1.8, alpha = 0.8) + 
  scale_y_discrete(expand = expansion(mult = c(0, 0.4))) +  
  scale_x_continuous(
    limits = c(min(cor_df$Translation_Efficiency), max(cor_df$Translation_Efficiency)),
    breaks = seq(0, max(cor_df$Translation_Efficiency), by = 0.05)
  ) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  labs(title = "TE Distribution Each Cell", 
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
