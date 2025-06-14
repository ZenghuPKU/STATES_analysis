# 03_Normalization_Visualization
# Load libraries and set environment
library(devtools)
library(Seurat)
library(anndata)
library(reticulate)
library(dplyr)
library(matrixStats)
library(ggplot2)
rm(list = ls())
setwd("~/tissue_downstream/")
use_condaenv("scanpy_env", required = TRUE)
py_config()

# Read h5ad data, create Seurat object, and add 'rbRNA' assay
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")
hvg_genes <- rownames(sc_ad$var)[sc_ad$var$highly_variable == TRUE]
totalRNA_matrix <- t(as.matrix(sc_ad$layers[['totalRNA_raw']]))[hvg_genes, , drop = FALSE]
rbRNA_matrix <- t(as.matrix(sc_ad$layers[['rbRNA']]))[hvg_genes, , drop = FALSE]
metadata <- sc_ad$obs
states <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata, assay = "totalRNA")
states[["rbRNA"]] <- CreateAssayObject(counts = rbRNA_matrix)

# Normalize 'totalRNA' assay data
DefaultAssay(states) <- "totalRNA"
raw_counts <- GetAssayData(states, assay = "totalRNA", slot = "counts")
scale_factor <- median(Matrix::colSums(raw_counts))
states <- NormalizeData(states, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
norm_counts <- GetAssayData(states, assay = "totalRNA", slot = "data")
states[["totalRNA"]]@data <- log1p(states[["totalRNA"]]@data)

# Normalize rbRNA data using TE matrix and log1p transform
DefaultAssay(states) <- "rbRNA"
te_matrix <- t(as.matrix(sc_ad$layers[['TE']]))[hvg_genes, , drop = FALSE]
norm_rbRNA <- norm_counts * te_matrix
states <- SetAssayData(states, assay = "rbRNA", slot = "data", new.data = norm_rbRNA)
states[["rbRNA"]]@data <- log1p(states[["rbRNA"]]@data)

# Visualization
DefaultAssay(states) <- "totalRNA"
states <- FindVariableFeatures(states, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
states <- ScaleData(states, features = rownames(states), assay = "totalRNA")
states <- RunPCA(states, reduction.name = "pca", npcs = 100)
states <- FindNeighbors(states, reduction = "pca", dims = 1:25, graph.name = "totalRNA_nn")
states <- FindClusters(states, graph.name = "totalRNA_nn", resolution = 2, algorithm = 1, method = "igraph")
states[["states_nn_alg1"]] <- Idents(states)

# Two approaches to generate UMAP embeddings for visualization:
# Method 1 loads UMAP embedding from the "mousebrain_umap_embeddings_with_labels.csv" file.
# Method 2 demonstrates how these embeddings were generated using dimensionality reduction with RunUMAP based on PCA results.

# Method 1: Load UMAP embedding from the "mousebrain_umap_embeddings_with_labels.csv" file
umap_df <- read.csv("mousebrain_umap_embeddings_with_labels.csv", row.names = 1)
head(colnames(states))
head(rownames(umap_df))   
custom_umap <- CreateDimReducObject(
  embeddings = as.matrix(umap_df[, c("UMAP_1", "UMAP_2")]),
  key = "UMAP_",
  assay = DefaultAssay(states)
)
states[["states.umap"]] <- custom_umap
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

#Output the normalized and visualized results
save(states, file = "states_Normalized_Visualized_output.RData")

# Method 2: Demonstrate how the UMAP embedding was generated using RunUMAP and PCA
states <- RunUMAP(states, reduction = "pca", dims = 1:25,
                       reduction.name = "states.umap", reduction.key = "UMAP_", n.neighbors = 50,
                       min.dist = 0.001, spread = 3)
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

#Output the normalized and visualized results
save(states, file = "states_Normalized_Visualized_output.RData")
