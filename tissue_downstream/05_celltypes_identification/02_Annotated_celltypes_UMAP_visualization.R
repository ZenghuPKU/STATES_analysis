# 02_Annotated_celltypes_UMAP_visualization
# Load libraries and set environment
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
rm(list = ls())
setwd("~/tissue_downstream/05_celltypes_identification/")

# Load the annotated object (states_celltypes_identification.RData)
load("states_celltypes_identification.RData")

# Extract UMAP coordinates
umap_coords <- Embeddings(states, reduction = "states.umap")

# Create a data frame that includes coordinates and labels
umap_df <- data.frame(
  cell_id = rownames(umap_coords),
  umap_coords,
  label1 = states@meta.data$states_nn_alg1_label1,
  label2 = states@meta.data$states_nn_alg1_label2,
  label3 = states@meta.data$states_nn_alg1_label3
)

# Save as CSV file
write.csv(umap_df, file = "mousebrain_umap_embeddings_with_labels.csv", row.names = FALSE)

# label1 UMAP visualization
Idents(states) <- "states_nn_alg1_label1"
table(states@active.ident)
celltype_colors <- c(
  "Non_Neuron" = "#66C2A5",     
  "Neuron" = "#F781BF",
  "Mix"="#F5F5F5"
)
P1=DimPlot(states, reduction = "states.umap", label = T, group.by = "states_nn_alg1_label1") +
  scale_color_manual(values = celltype_colors) +
  NoLegend()+
  ggtitle(NULL)
P1
ggsave("states_Umap_label1.png", P1, width = 8, height = 8)


# label2 UMAP visualization
Idents(states) <- "states_nn_alg1_label2"
table(states@active.ident)
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
P2=DimPlot(states, reduction = "states.umap", label = T, group.by = "states_nn_alg1_label2") +
  scale_color_manual(values = celltype_colors) +
  NoLegend() +
  ggtitle(NULL)
P2
ggsave("states_Umap_label2.pdf", P2, width = 8, height = 8, device = cairo_pdf)

# label3 UMAP visualization
Idents(states) <- "states_nn_alg1_label3"
table(states@active.ident)
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
P3=DimPlot(states, reduction = "states.umap", label = T, group.by = "states_nn_alg1_label3") +
  scale_color_manual(values = celltype_colors)+
  NoLegend() +
  ggtitle(NULL)
P3
ggsave("states_Umap_label3.pdf", P3, width = 8, height = 8, device = cairo_pdf)


