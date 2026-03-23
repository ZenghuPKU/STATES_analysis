# 02_celltypes_identification_Subcelltype
rm(list = ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyr)
setwd("~/AD_downstream/04_celltypes_identification/")

states=readRDS("~/AD_downstream/04_celltypes_identification/states_Major_celltype.rds")

Idents(states) <- "states_nn_alg1_label2"

DefaultAssay(states) <- "totalRNA"
raw_counts <- GetAssayData(states, assay = "totalRNA", slot = "counts")
total_counts <- Matrix::colSums(raw_counts)
scale_factor <- median(total_counts)
cat("RC：", scale_factor, "\n")

# MLG Normalize
MLG_cells <- WhichCells(states, idents = "MLG")
seurat_sub <- subset(states, cells = MLG_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("MLG_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[MLG_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$MLG_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "MLG_subcluster"

# VAS Normalize
VAS_cells <- WhichCells(states, idents = "VAS")
seurat_sub <- subset(states, cells = VAS_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("VAS_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[VAS_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$VAS_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "VAS_subcluster"

# AC Normalize
AC_cells <- WhichCells(states, idents = "AC")
seurat_sub <- subset(states, cells = AC_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("AC_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[AC_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$AC_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "AC_subcluster"

# AC_0 Normalize
AC_0_cells <- WhichCells(states, idents = "AC_0")
seurat_sub <- subset(states, cells = AC_0_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("AC_0_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[AC_0_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$AC_0_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "AC_0_subcluster"

# AC_2 Normalize
AC_2_cells <- WhichCells(states, idents = "AC_2")
seurat_sub <- subset(states, cells = AC_2_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("AC_2_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[AC_2_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$AC_2_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "AC_2_subcluster"


# CHOR/EPEN Normalize
CHOR_EPEN_cells <- WhichCells(states, idents = "CHOR/EPEN")
seurat_sub <- subset(states, cells = CHOR_EPEN_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("CHOR/EPEN_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[CHOR_EPEN_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$CHOR_EPEN_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "CHOR_EPEN_subcluster"

# INH Normalize
INH_cells <- WhichCells(states, idents = "INH")
seurat_sub <- subset(states, cells = INH_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("INH_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[INH_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$INH_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "INH_subcluster"

# OLG Normalize
OLG_cells <- WhichCells(states, idents = "OLG")
seurat_sub <- subset(states, cells = OLG_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("OLG_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[OLG_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$OLG_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "OLG_subcluster"

# TEPN Normalize
TEPN_cells <- WhichCells(states, idents = "TEPN")
seurat_sub <- subset(states, cells = TEPN_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_subcluster"

# TEPN_0 Normalize
TEPN_0_cells <- WhichCells(states, idents = "TEPN_0")
seurat_sub <- subset(states, cells = TEPN_0_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_0_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_0_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_0_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_0_subcluster"

# TEPN_1 Normalize
TEPN_1_cells <- WhichCells(states, idents = "TEPN_1")
seurat_sub <- subset(states, cells = TEPN_1_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_1_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_1_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_1_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_1_subcluster"

# TEPN_5 Normalize
TEPN_5_cells <- WhichCells(states, idents = "TEPN_5")
seurat_sub <- subset(states, cells = TEPN_5_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_5_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_5_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_5_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_5_subcluster"

# TEPN_7 Normalize
TEPN_7_cells <- WhichCells(states, idents = "TEPN_7")
seurat_sub <- subset(states, cells = TEPN_7_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_7_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_7_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_7_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_7_subcluster"

# TEPN_0_0 Normalize
TEPN_0_0_cells <- WhichCells(states, idents = "TEPN_0_0")
seurat_sub <- subset(states, cells = TEPN_0_0_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_0_0_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_0_0_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_0_0_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_0_0_subcluster"

# TEPN_1_0 Normalize
TEPN_1_0_cells <- WhichCells(states, idents = "TEPN_1_0")
seurat_sub <- subset(states, cells = TEPN_1_0_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_1_0_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_1_0_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_1_0_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_1_0_subcluster"

# TEPN_5_1 Normalize
TEPN_5_1_cells <- WhichCells(states, idents = "TEPN_5_1")
seurat_sub <- subset(states, cells = TEPN_5_1_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_5_1_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_5_1_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_5_1_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_5_1_subcluster"

# TEPN_5_3 Normalize
TEPN_5_3_cells <- WhichCells(states, idents = "TEPN_5_3")
seurat_sub <- subset(states, cells = TEPN_5_3_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:25,graph.name = "sub_nn")
seurat_sub <- FindClusters(seurat_sub,graph.name = "sub_nn")
table(seurat_sub$group,Idents(seurat_sub))
table(seurat_sub$type,Idents(seurat_sub))
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_5_3_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_5_3_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)
seurat_sub$TEPN_5_3_subcluster <- Idents(seurat_sub)
Idents(seurat_sub) <- "TEPN_5_3_subcluster"

# Identification of neurons and non-neurons (label1)
Idents(states) <- "states_nn_alg1_label2_new"
new.cluster.ids <- c(  "MLG_0" = "Non_Neuron",
                       "MLG_1" = "Non_Neuron",
                       "MLG_2" = "Non_Neuron",
                       "MLG_3" = "Non_Neuron",
                       "MLG_4" = "Non_Neuron",
                       "MLG_5" = "Non_Neuron",
                       "MLG_6" = "Non_Neuron",
                       "MLG_7" = "Non_Neuron",
                       "VAS_0" = "Non_Neuron",
                       "VAS_1" = "Non_Neuron",
                       "VAS_2" = "Non_Neuron",
                       "VAS_3" = "Non_Neuron",
                       "VAS_4" = "Non_Neuron",
                       "VAS_5" = "Non_Neuron",
                       "VAS_6" = "Non_Neuron",
                       "VAS_7" = "Non_Neuron",
                       "VAS_8" = "Non_Neuron",
                       "CHOR/EPEN_0" = "Non_Neuron",
                       "CHOR/EPEN_1" = "Non_Neuron",
                       "CHOR/EPEN_2" = "Non_Neuron",
                       "CHOR/EPEN_3" = "Non_Neuron",
                       "CHOR/EPEN_4" = "Non_Neuron",
                       "INH_0" = "Neuron",
                       "INH_1" = "Neuron",
                       "INH_2" = "Neuron",
                       "INH_3" = "Neuron",
                       "INH_4" = "Neuron",
                       "INH_5" = "Neuron",
                       "INH_6" = "Neuron",
                       "OLG_0" = "Non_Neuron",
                       "OLG_1" = "Non_Neuron",
                       "OLG_2" = "Non_Neuron",
                       "OLG_3" = "Non_Neuron",
                       "OLG_4" = "Non_Neuron",
                       "OLG_5" = "Non_Neuron",
                       "OLG_6" = "Non_Neuron",
                       "OLG_7" = "Non_Neuron",
                       "OLG_8" = "Non_Neuron",
                       "OLG_9" = "Non_Neuron",
                       "OLG_10" = "Non_Neuron",
                       "TEPN_0_0_0" = "Neuron",
                       "TEPN_0_0_1" = "Neuron",
                       "TEPN_0_0_2" = "Neuron",
                       "TEPN_0_0_3" = "Neuron",
                       "TEPN_0_1" = "Neuron",
                       "TEPN_0_2" = "Neuron",
                       "TEPN_0_3" = "Neuron",
                       "TEPN_1_0_0" = "Neuron",
                       "TEPN_1_0_1" = "Neuron",
                       "TEPN_1_0_2" = "Neuron",
                       "TEPN_1_0_3" = "Neuron",
                       "TEPN_1_1" = "Neuron",
                       "TEPN_1_2" = "Neuron",
                       "TEPN_1_3" = "Neuron",
                       "TEPN_1_4" = "Neuron",
                       "TEPN_1_5" = "Neuron",
                       "TEPN_2" = "Neuron",
                       "TEPN_3" = "Neuron",
                       "TEPN_4" = "Neuron",
                       "TEPN_5_0" = "Neuron",
                       "TEPN_5_1_0" = "Neuron",
                       "TEPN_5_1_1" = "Neuron",
                       "TEPN_5_1_2" = "Neuron",
                       "TEPN_5_1_3" = "Neuron",
                       "TEPN_5_1_4" = "Neuron",
                       "TEPN_5_2" = "Neuron",
                       "TEPN_5_3_0" = "Neuron",
                       "TEPN_5_3_1" = "Neuron",
                       "TEPN_5_3_2" = "Neuron",
                       "TEPN_5_3_3" = "Neuron",
                       "TEPN_5_4" = "Neuron",
                       "TEPN_5_5" = "Neuron",
                       "TEPN_6" = "Neuron",
                       "TEPN_7_0" = "Neuron",
                       "TEPN_7_1" = "Neuron",
                       "TEPN_7_2" = "Neuron",
                       "TEPN_7_3" = "Neuron",
                       "TEPN_7_4" = "Neuron",
                       "TEPN_8" = "Neuron",
                       "TEPN_9" = "Neuron",
                       "TEPN_10" = "Neuron",
                       "TEPN_11" = "Neuron",
                       "TEPN_12" = "Neuron",
                       "TEPN_13" = "Neuron",
                       "TEPN_14" = "Neuron",
                       "TEPN_15" = "Neuron",
                       "AC_0_0" = "Mix",
                       "AC_0_1" = "Mix",
                       "AC_0_2" = "Non_Neuron",
                       "AC_0_3" = "Mix",
                       "AC_0_4" = "Non_Neuron",
                       "AC_0_5" = "Mix",
                       "AC_1" = "Non_Neuron",
                       "AC_2_0" = "Non_Neuron",
                       "AC_2_1" = "Mix",
                       "AC_2_2" = "Non_Neuron",
                       "AC_2_3" = "Non_Neuron",
                       "AC_2_4" = "Mix",
                       "AC_3" = "Non_Neuron",
                       "AC_4" = "Non_Neuron",
                       "AC_5" = "Non_Neuron",
                       "AC_6" = "Non_Neuron",
                       "OPC" = "Non_Neuron",
                       "CHO/PEP" = "Neuron",
                       "DE/MEN" = "Neuron"
)
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label1<- states@active.ident
Idents(states) <- "states_nn_alg1_label1"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Identification of major celltypes (label2)
Idents(states) <- "states_nn_alg1_label2_new"
new.cluster.ids <- c(  "MLG_0" = "MLG",
                       "MLG_1" = "MLG",
                       "MLG_2" = "MLG",
                       "MLG_3" = "MLG",
                       "MLG_4" = "MLG",
                       "MLG_5" = "MLG",
                       "MLG_6" = "MLG",
                       "MLG_7" = "MLG",
                       "VAS_0" = "VAS",
                       "VAS_1" = "VAS",
                       "VAS_2" = "VAS",
                       "VAS_3" = "VAS",
                       "VAS_4" = "VAS",
                       "VAS_5" = "VAS",
                       "VAS_6" = "VAS",
                       "VAS_7" = "VAS",
                       "VAS_8" = "VAS",
                       "CHOR/EPEN_0" = "CHOR/EPEN",
                       "CHOR/EPEN_1" = "CHOR/EPEN",
                       "CHOR/EPEN_2" = "CHOR/EPEN",
                       "CHOR/EPEN_3" = "CHOR/EPEN",
                       "CHOR/EPEN_4" = "CHOR/EPEN",
                       "INH_0" = "INH",
                       "INH_1" = "INH",
                       "INH_2" = "INH",
                       "INH_3" = "INH",
                       "INH_4" = "INH",
                       "INH_5" = "INH",
                       "INH_6" = "INH",
                       "OLG_0" = "OLG",
                       "OLG_1" = "OLG",
                       "OLG_2" = "OLG",
                       "OLG_3" = "OLG",
                       "OLG_4" = "OLG",
                       "OLG_5" = "OLG",
                       "OLG_6" = "OLG",
                       "OLG_7" = "OLG",
                       "OLG_8" = "OLG",
                       "OLG_9" = "OLG",
                       "OLG_10" = "OLG",
                       "TEPN_0_0_0" = "TEPN",
                       "TEPN_0_0_1" = "TEPN",
                       "TEPN_0_0_2" = "TEPN",
                       "TEPN_0_0_3" = "TEPN",
                       "TEPN_0_1" = "TEPN",
                       "TEPN_0_2" = "TEPN",
                       "TEPN_0_3" = "TEPN",
                       "TEPN_1_0_0" = "TEPN",
                       "TEPN_1_0_1" = "TEPN",
                       "TEPN_1_0_2" = "TEPN",
                       "TEPN_1_0_3" = "TEPN",
                       "TEPN_1_1" = "TEPN",
                       "TEPN_1_2" = "TEPN",
                       "TEPN_1_3" = "TEPN",
                       "TEPN_1_4" = "TEPN",
                       "TEPN_1_5" = "TEPN",
                       "TEPN_2" = "TEPN",
                       "TEPN_3" = "TEPN",
                       "TEPN_4" = "TEPN",
                       "TEPN_5_0" = "TEPN",
                       "TEPN_5_1_0" = "TEPN",
                       "TEPN_5_1_1" = "TEPN",
                       "TEPN_5_1_2" = "TEPN",
                       "TEPN_5_1_3" = "TEPN",
                       "TEPN_5_1_4" = "TEPN",
                       "TEPN_5_2" = "TEPN",
                       "TEPN_5_3_0" = "TEPN",
                       "TEPN_5_3_1" = "TEPN",
                       "TEPN_5_3_2" = "TEPN",
                       "TEPN_5_3_3" = "TEPN",
                       "TEPN_5_4" = "TEPN",
                       "TEPN_5_5" = "TEPN",
                       "TEPN_6" = "TEPN",
                       "TEPN_7_0" = "TEPN",
                       "TEPN_7_1" = "TEPN",
                       "TEPN_7_2" = "TEPN",
                       "TEPN_7_3" = "TEPN",
                       "TEPN_7_4" = "TEPN",
                       "TEPN_8" = "TEPN",
                       "TEPN_9" = "TEPN",
                       "TEPN_10" = "TEPN",
                       "TEPN_11" = "TEPN",
                       "TEPN_12" = "TEPN",
                       "TEPN_13" = "TEPN",
                       "TEPN_14" = "TEPN",
                       "TEPN_15" = "TEPN",
                       "AC_0_0" = "Mix",
                       "AC_0_1" = "Mix",
                       "AC_0_2" = "AC",
                       "AC_0_3" = "Mix",
                       "AC_0_4" = "AC",
                       "AC_0_5" = "Mix",
                       "AC_1" = "AC",
                       "AC_2_0" = "AC",
                       "AC_2_1" = "Mix",
                       "AC_2_2" = "AC",
                       "AC_2_3" = "AC",
                       "AC_2_4" = "Mix",
                       "AC_3" = "AC",
                       "AC_4" = "AC",
                       "AC_5" = "AC",
                       "AC_6" = "AC"
)
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label2<- states@active.ident
Idents(states) <- "states_nn_alg1_label2"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()


# Identification of subcelltypes (label3)
Idents(states) <- "states_nn_alg1_label2_new"
table(states$states_nn_alg1_label2_new)
new.cluster.ids <- c(  "MLG_0" = "MLG3",
                       "MLG_1" = "MLG1",
                       "MLG_2" = "MLG1",
                       "MLG_3" = "MLG2",
                       "MLG_4" = "MLG2",
                       "MLG_5" = "MLG1",
                       "MLG_6" = "MLG1",
                       "MLG_7" = "MLG3",
                       "VAS_0" = "Peri/VEC",
                       "VAS_1" = "VLMC",
                       "VAS_2" = "Peri/VEC",
                       "VAS_3" = "VLMC",
                       "VAS_4" = "VSMC",
                       "VAS_5" = "Peri/VEC",
                       "VAS_6" = "Peri/VEC",
                       "VAS_7" = "Peri/VEC",
                       "VAS_8" = "VLMC",
                       "CHOR/EPEN_0" = "EPEN",
                       "CHOR/EPEN_1" = "CHOR",
                       "CHOR/EPEN_2" = "CHOR",
                       "CHOR/EPEN_3" = "CHOR",
                       "CHOR/EPEN_4" = "EPEN",
                       "INH_0" = "INH_Sst",
                       "INH_1" = "INH_Pvalb",
                       "INH_2" = "INH_Sst",
                       "INH_3" = "INH_Pvalb",
                       "INH_4" = "INH_Cnr1_Vip",
                       "INH_5" = "INH_Pvalb",
                       "INH_6" = "INH_Sst",
                       "OLG_0" = "OLG1",
                       "OLG_1" = "OLG2",
                       "OLG_2" = "OLG2",
                       "OLG_3" = "OLG2",
                       "OLG_4" = "OLG2",
                       "OLG_5" = "OLG2",
                       "OLG_6" = "OLG1",
                       "OLG_7" = "OLG2",
                       "OLG_8" = "OLG2",
                       "OLG_9" = "OLG2",
                       "OLG_10" = "OLG2",
                       "TEPN_0_0_0" = "TEGLU L4",
                       "TEPN_0_0_1" = "TEGLU L2/3",
                       "TEPN_0_0_2" = "TEGLU L4",
                       "TEPN_0_0_3" = "TEGLU L2/3",
                       "TEPN_0_1" = "TEGLU L6",
                       "TEPN_0_2" = "TEGLU L4",
                       "TEPN_0_3" = "TEGLU L2/3",
                       "TEPN_1_0_0" = "TEGLU L4/5",
                       "TEPN_1_0_1" = "TEGLU L4/5",
                       "TEPN_1_0_2" = "TEGLU L6",
                       "TEPN_1_0_3" = "TEGLU Mix",
                       "TEPN_1_1" = "TEGLU L2/3",
                       "TEPN_1_2" = "TEGLU L4/5",
                       "TEPN_1_3" = "TEGLU Mix",
                       "TEPN_1_4" = "TEGLU Mix",
                       "TEPN_1_5" = "TEGLU Mix",
                       "TEPN_2" = "TEGLU L6",
                       "TEPN_3" = "DGGRC",
                       "TEPN_4" = "TEGLU L2/3",
                       "TEPN_5_0" = "TEGLU CA3",
                       "TEPN_5_1_0" = "TEGLU CA1",
                       "TEPN_5_1_1" = "TEGLU CA2",
                       "TEPN_5_1_2" = "TEGLU CA2",
                       "TEPN_5_1_3" = "TEGLU Mix",
                       "TEPN_5_1_4" = "TEGLU Mix",
                       "TEPN_5_2" = "TEGLU CA3",
                       "TEPN_5_3_0" = "TEGLU CA2",
                       "TEPN_5_3_1" = "TEGLU CA1",
                       "TEPN_5_3_2" = "TEGLU CA1",
                       "TEPN_5_3_3" = "TEGLU CA2",
                       "TEPN_5_4" = "TEGLU CA3",
                       "TEPN_5_5" = "TEGLU Mix",
                       "TEPN_6" = "TEGLU L5",
                       "TEPN_7_0" = "TEGLU L4/5",
                       "TEPN_7_1" = "TEGLU L4",
                       "TEPN_7_2" = "TEGLU L6",
                       "TEPN_7_3" = "TEGLU Mix",
                       "TEPN_7_4" = "TEGLU Mix",
                       "TEPN_8" = "MSN",
                       "TEPN_9" = "TEGLU CA1",
                       "TEPN_10" = "TEGLU Mix",
                       "TEPN_11" = "MSN",
                       "TEPN_12" = "DGGRC",
                       "TEPN_13" = "TEGLU Mix",
                       "TEPN_14" = "TEGLU Mix",
                       "TEPN_15" = "MSN",
                       "AC_0_0" = "Mix",
                       "AC_0_1" = "Mix",
                       "AC_0_2" = "AC3",
                       "AC_0_3" = "Mix",
                       "AC_0_4" = "AC3",
                       "AC_0_5" = "Mix",
                       "AC_1" = "AC2",
                       "AC_2_0" = "AC4",
                       "AC_2_1" = "Mix",
                       "AC_2_2" = "AC4",
                       "AC_2_3" = "AC4",
                       "AC_2_4" = "Mix",
                       "AC_3" = "AC1",
                       "AC_4" = "AC4",
                       "AC_5" = "AC1",
                       "AC_6" = "AC2"
)

states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label3<- states@active.ident
Idents(states) <- "states_nn_alg1_label3"
table(states@active.ident)
table(states$type, Idents(states))
table(states$group, Idents(states))
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Output the cell types identification result
saveRDS(states, file = "states_celltypes_identification.rds")
