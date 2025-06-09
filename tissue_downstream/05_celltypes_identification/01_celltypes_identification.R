# 01_celltypes_identification
# Load libraries and set environment
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("~/tissue_downstream/05_celltypes_identification/")
rm(list = ls())

# Load the Mixfind object (states_mixsfind.Rdata)
load("states_mixsfind.Rdata")

# Find highly variable genes
Idents(states) <- "states_nn_alg1_new"
markers = FindAllMarkers(object = states, graph.name = "totalRNA_nn",test.use="wilcox" ,only.pos = TRUE,logfc.threshold = 0.1,min.pct = 0.25)   
all.markers = markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox_alg1_sub.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox_alg1_sub.csv", row.names = F)

# Draw dot plot showing the expression status of certain marker genes
C=c("Pdgfra","Olig2")
DotPlot(states,features = C)+RotatedAxis()

# Initial Identification of neurons and non-neurons (label1)
new.cluster.ids <- c( "0" = "Neuron",
                      "1" = "Non_Neuron",
                      "2" = "Neuron",
                      "3" = "Neuron",
                      "4" = "Non_Neuron",
                      "5" = "Non_Neuron",
                      "6" = "Non_Neuron",
                      "7" = "Neuron",
                      "8" = "Neuron",
                      "9" = "Non_Neuron",
                      "10" = "Neuron",
                      "11" = "Neuron",
                      "12" = "Neuron",
                      "13" = "Neuron",
                      "14" = "Non_Neuron",
                      "15" = "Neuron",
                      "16" = "Neuron",
                      "17" = "Non_Neuron",
                      "18" = "Non_Neuron",
                      "19" = "Neuron",
                      "21" = "Non_Neuron",
                      "22" = "Neuron",
                      "23" = "Non_Neuron",
                      "24" = "Non_Neuron")
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label1<- states@active.ident
Idents(states) <- "states_nn_alg1_label1"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Initial identification of the main cell types
Idents(states) <- "states_nn_alg1_new"
new.cluster.ids <- c( "0" = "TEPN",
                      "1" = "AC",
                      "2" = "TEPN",
                      "3" = "TEPN",
                      "4" = "VAS",
                      "5" = "OLG",
                      "6" = "OLG",
                      "7" = "TEPN",
                      "8" = "TEPN",
                      "9" = "OLG",
                      "10" = "INH",
                      "11" = "TEPN",
                      "12" = "TEPN",
                      "13" = "TEPN",
                      "14" = "VAS",
                      "15" = "TEPN",
                      "16" = "TEPN",
                      "17" = "MLG",
                      "18" = "CHOR/EPEN",
                      "19" = "DE/MEN",
                      "21" = "CHOR/EPEN",
                      "22" = "CHO/PEP",
                      "23" = "OPC",
                      "24" = "VAS")
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label2<- states@active.ident
Idents(states) <- "states_nn_alg1_label2"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Identification of major cell types and cell subtypes (Major celltypes: label2; Subcelltypes: label3)
# Normalize data (method: RC)
DefaultAssay(states) <- "totalRNA"
raw_counts <- GetAssayData(states, assay = "totalRNA", slot = "counts")
total_counts <- Matrix::colSums(raw_counts)
scale_factor <- median(total_counts)
cat("中位数归一化因子：", scale_factor, "\n")

# AC Normalize
AC_cells <- WhichCells(states, idents = "AC")
seurat_sub <- subset(states, cells = AC_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("AC_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[AC_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident

# VAS Normalize
VAS_cells <- WhichCells(states, idents = "VAS")
seurat_sub <- subset(states, cells = VAS_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("VAS_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[VAS_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident

# OLG Normalize
OLG_cells <- WhichCells(states, idents = "OLG")
seurat_sub <- subset(states, cells = OLG_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("OLG_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[OLG_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident

# INH Normalize
INH_cells <- WhichCells(states, idents = "INH")
seurat_sub <- subset(states, cells = INH_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("INH_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[INH_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident

# CHOR/EPEN Normalize
CHOR_EPEN_cells <- WhichCells(states, idents = "CHOR/EPEN")
seurat_sub <- subset(states, cells = CHOR_EPEN_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("CHOR/EPEN_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[CHOR_EPEN_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident

# TEPN Normalize
TEPN_cells <- WhichCells(states, idents = "TEPN")
seurat_sub <- subset(states, cells = TEPN_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# TEPN_0 Normalize
TEPN_0_cells <- WhichCells(states, idents = "TEPN_0")
seurat_sub <- subset(states, cells = TEPN_0_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub,resolution = 1)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_0_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_0_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# TEPN_1 Normalize
TEPN_1_cells <- WhichCells(states, idents = "TEPN_1")
seurat_sub <- subset(states, cells = TEPN_1_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub,resolution = 1)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_1_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_1_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# TEPN_2 Normalize
TEPN_2_cells <- WhichCells(states, idents = "TEPN_2")
seurat_sub <- subset(states, cells = TEPN_2_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub,resolution = 1)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_2_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_2_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# TEPN_3 Normalize
TEPN_3_cells <- WhichCells(states, idents = "TEPN_3")
seurat_sub <- subset(states, cells = TEPN_3_cells)
seurat_sub <- NormalizeData(seurat_sub, assay = "totalRNA", normalization.method = "RC", scale.factor = scale_factor)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 1500, assay = "totalRNA")
seurat_sub <- ScaleData(seurat_sub)
seurat_sub <- RunPCA(seurat_sub)
ElbowPlot(seurat_sub, ndims = 50, reduction = "pca")
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:30)
seurat_sub <- FindClusters(seurat_sub,resolution = 1)
new_clusters <- Idents(seurat_sub)
levels(new_clusters) <- paste0("TEPN_3_", seq(0, length(levels(new_clusters))-1))
all_idents <- Idents(states)
levels(all_idents) <- union(levels(all_idents), levels(new_clusters))
all_idents[TEPN_3_cells] <- new_clusters
Idents(states) <- all_idents
table(Idents(states))
states$states_nn_alg1_label2_new<- states@active.ident
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Identification of neurons and non-neurons (label1)
Idents(states) <- "states_nn_alg1_label2_new"
new.cluster.ids <- c( "TEPN_0_0" = "Neuron",
                      "TEPN_0_1" = "Neuron",
                      "TEPN_0_2" = "Neuron",
                      "TEPN_0_3" = "Neuron",
                      "TEPN_0_4" = "Neuron",
                      "TEPN_1_0" = "Neuron",
                      "TEPN_1_1" = "Neuron",
                      "TEPN_1_2" = "Neuron",
                      "TEPN_1_3" = "Neuron",
                      "TEPN_1_4" = "Neuron",
                      "TEPN_1_5" = "Neuron",
                      "TEPN_2_0" = "Neuron",
                      "TEPN_2_1" = "Neuron",
                      "TEPN_2_2" = "Neuron",
                      "TEPN_2_3" = "Neuron",
                      "TEPN_2_4" = "Neuron",
                      "TEPN_2_5" = "Neuron",
                      "TEPN_3_0" = "Neuron",
                      "TEPN_3_1" = "Neuron",
                      "TEPN_3_2" = "Neuron",
                      "TEPN_3_3" = "Neuron",
                      "TEPN_4" = "Neuron",
                      "TEPN_5" = "Neuron",
                      "TEPN_6" = "Neuron",
                      "TEPN_7" = "Neuron",
                      "TEPN_8" = "Non_Neuron",
                      "TEPN_9" = "Neuron",
                      "TEPN_10" = "Neuron",
                      "TEPN_11" = "Neuron",
                      "AC_0" = "Non_Neuron",
                      "AC_1" = "Non_Neuron",
                      "AC_2" = "Non_Neuron",
                      "AC_3" = "Non_Neuron",
                      "VAS_0" = "Non_Neuron",
                      "VAS_1" = "Non_Neuron",
                      "VAS_2" = "Non_Neuron",
                      "VAS_3" = "Non_Neuron",
                      "VAS_4" = "Non_Neuron",
                      "VAS_5" = "Non_Neuron",
                      "VAS_6" = "Non_Neuron",
                      "VAS_7" = "Non_Neuron",
                      "INH_0" = "Neuron",
                      "INH_1" = "Neuron",
                      "INH_2" = "Neuron",
                      "INH_3" = "Neuron",
                      "DE/MEN" = "Neuron",
                      "CHO/PEP" = "Neuron",
                      "OPC" = "Non_Neuron",
                      "MLG" = "Non_Neuron",
                      "CHOR/EPEN_0" = "Non_Neuron",
                      "CHOR/EPEN_1" = "Non_Neuron",
                      "CHOR/EPEN_2" = "Non_Neuron",
                      "CHOR/EPEN_3" = "Non_Neuron",
                      "OLG_0" = "Non_Neuron",
                      "OLG_1" = "Non_Neuron",
                      "OLG_2" = "Non_Neuron",
                      "OLG_3" = "Non_Neuron",
                      "OLG_4" = "Non_Neuron"
)
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label1<- states@active.ident
Idents(states) <- "states_nn_alg1_label1"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Identification of major celltypes (label2)
Idents(states) <- "states_nn_alg1_label2_new"
new.cluster.ids <- c("TEPN_0_0" = "TEPN",
                     "TEPN_0_1" = "TEPN",
                     "TEPN_0_2" = "TEPN",
                     "TEPN_0_3" = "TEPN",
                     "TEPN_0_4" = "TEPN",
                     "TEPN_1_0" = "TEPN",
                     "TEPN_1_1" = "TEPN",
                     "TEPN_1_2" = "TEPN",
                     "TEPN_1_3" = "TEPN",
                     "TEPN_1_4" = "TEPN",
                     "TEPN_1_5" = "TEPN",
                     "TEPN_2_0" = "TEPN",
                     "TEPN_2_1" = "TEPN",
                     "TEPN_2_2" = "TEPN",
                     "TEPN_2_3" = "TEPN",
                     "TEPN_2_4" = "TEPN",
                     "TEPN_2_5" = "TEPN",
                     "TEPN_3_0" = "TEPN",
                     "TEPN_3_1" = "TEPN",
                     "TEPN_3_2" = "TEPN",
                     "TEPN_3_3" = "TEPN",
                     "TEPN_4" = "TEPN",
                     "TEPN_5" = "TEPN",
                     "TEPN_6" = "TEPN",
                     "TEPN_7" = "TEPN",
                     "TEPN_8" = "AC",
                     "TEPN_9" = "TEPN",
                     "TEPN_10" = "TEPN",
                     "TEPN_11" = "TEPN",
                     "AC_0" = "AC",
                     "AC_1" = "AC",
                     "AC_2" = "AC",
                     "AC_3" = "AC",
                     "VAS_0" = "VAS",
                     "VAS_1" = "VAS",
                     "VAS_2" = "VAS",
                     "VAS_3" = "VAS",
                     "VAS_4" = "VAS",
                     "VAS_5" = "VAS",
                     "VAS_6" = "VAS",
                     "VAS_7" = "VAS",
                     "INH_0" = "INH",
                     "INH_1" = "INH",
                     "INH_2" = "INH",
                     "INH_3" = "INH",
                     "CHOR/EPEN_0" = "CHOR/EPEN",
                     "CHOR/EPEN_1" = "CHOR/EPEN",
                     "CHOR/EPEN_2" = "CHOR/EPEN",
                     "CHOR/EPEN_3" = "CHOR/EPEN",
                     "OLG_0" = "OLG",
                     "OLG_1" = "OLG",
                     "OLG_2" = "OLG",
                     "OLG_3" = "OLG",
                     "OLG_4" = "OLG"
)
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label2<- states@active.ident
Idents(states) <- "states_nn_alg1_label2"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Identification of subcelltypes (label3)
Idents(states) <- "states_nn_alg1_label2_new"
new.cluster.ids <- c( "TEPN_0_0" = "TEGLU L2/3/4",
                      "TEPN_0_1" = "TEGLU L2/3/4",
                      "TEPN_0_2" = "TEGLU L5/6",
                      "TEPN_0_3" = "TEGLU L2/3",
                      "TEPN_0_4" = "TEGLU Mix",
                      "TEPN_1_0" = "TEGLU Mix",
                      "TEPN_1_1" = "TEGLU L2/3",
                      "TEPN_1_2" = "TEGLU L2/3",
                      "TEPN_1_3" = "TEGLU L2/3",
                      "TEPN_1_4" = "TEGLU L2/3",
                      "TEPN_1_5" = "TEGLU Mix",
                      "TEPN_2_0" = "TEGLU L6",
                      "TEPN_2_1" = "TEGLU L6b",
                      "TEPN_2_2" = "TEGLU L6",
                      "TEPN_2_3" = "TEGLU L5/6",
                      "TEPN_2_4" = "TEGLU L5/6",
                      "TEPN_2_5" = "TEGLU Mix",
                      "TEPN_3_0" = "TEGLU L2/3",
                      "TEPN_3_1" = "TEGLU L5/6",
                      "TEPN_3_2" = "TEGLU L5",
                      "TEPN_3_3" = "TEGLU L5",
                      "TEPN_4" = "DGGRC",
                      "TEPN_5" = "TEGLU CA1",
                      "TEPN_6" = "TEGLU CA3",
                      "TEPN_7" = "MSN",
                      "TEPN_8" = "AC3",
                      "TEPN_9" = "TEGLU CA2",
                      "TEPN_10" = "TEGLU Mix",
                      "TEPN_11" = "TEGLU Mix",
                      "AC_0" = "AC1",
                      "AC_1" = "AC2",
                      "AC_2" = "AC3",
                      "AC_3" = "AC1",
                      "VAS_0" = "Peri/VEC",
                      "VAS_1" = "Peri/VEC",
                      "VAS_2" = "VLMC",
                      "VAS_3" = "VLMC",
                      "VAS_4" = "VSMC",
                      "VAS_5" = "Peri/VEC",
                      "VAS_6" = "VLMC",
                      "VAS_7" = "VLMC",
                      "INH_0" = "INH_Sst",
                      "INH_1" = "INH_Cnr1_Vip",
                      "INH_2" = "INH_Pvalb",
                      "INH_3" = "INH_Pvalb",
                      "CHOR/EPEN_0" = "CHOR",
                      "CHOR/EPEN_1" = "EPEN",
                      "CHOR/EPEN_2" = "CHOR",
                      "CHOR/EPEN_3" = "CHOR",
                      "OLG_0" = "OLG1",
                      "OLG_1" = "OLG1",
                      "OLG_2" = "OLG1",
                      "OLG_3" = "OLG2",
                      "OLG_4" = "OLG2"
)
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label3<- states@active.ident
Idents(states) <- "states_nn_alg1_label3"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

# Output the cell types identification result
save(states,file = "states_celltypes_identification.RData")
