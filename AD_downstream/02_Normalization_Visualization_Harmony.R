# 02_Normalization_Visualization_Harmony.R
rm(list = ls())
library(Seurat)
library(anndata)
library(dplyr)
library(ggplot2)
library(Matrix)
library(harmony)
setwd("~/AD_downstream/")

h5ad_file <- "ADcombined_mousebrain.h5ad"
pdf_file  <- "ADcombined_harmony_UMAP.pdf"
rds_file  <- "states_ADcombined_harmony.rds"

sc_ad      <- read_h5ad(h5ad_file)
meta       <- sc_ad$obs
cell_names <- rownames(meta)
gene_names <- rownames(sc_ad$var)

totalRNA_raw <- t(as.matrix(sc_ad$layers[["totalRNA_raw"]]))      
rownames(totalRNA_raw) <- gene_names
colnames(totalRNA_raw) <- cell_names

states <- CreateSeuratObject(
  counts    = totalRNA_raw,
  meta.data = meta,
  assay     = "totalRNA"
)
DefaultAssay(states) <- "totalRNA"

totalRNA_log <- t(as.matrix(sc_ad$layers[["totalRNA_log"]]))
rownames(totalRNA_log) <- gene_names
colnames(totalRNA_log) <- cell_names
states[["totalRNA"]]@data <- totalRNA_log

rbRNA_raw <- t(as.matrix(sc_ad$layers[["rbRNA_raw"]]))
rbRNA_log <- t(as.matrix(sc_ad$layers[["rbRNA_log"]]))

rownames(rbRNA_raw) <- gene_names
colnames(rbRNA_raw) <- cell_names
rownames(rbRNA_log) <- gene_names
colnames(rbRNA_log) <- cell_names

rbRNA_assay        <- CreateAssayObject(counts = rbRNA_raw)
rbRNA_assay@data   <- rbRNA_log
states[["rbRNA"]]  <- rbRNA_assay

TE_mat <- t(as.matrix(sc_ad$layers[["TE"]]))
rownames(TE_mat) <- gene_names
colnames(TE_mat) <- cell_names

TE_assay       <- CreateAssayObject(counts = TE_mat)
TE_assay@data  <- TE_mat
states[["TE"]] <- TE_assay

DefaultAssay(states) <- "totalRNA"

states <- FindVariableFeatures(
  states,
  selection.method = "vst",
  nfeatures        = 1500,
  assay            = "totalRNA"
)

states <- ScaleData(states, features = rownames(states), assay = "totalRNA")

states <- RunPCA(
  states,
  features = VariableFeatures(states),
  npcs     = 50,
  verbose  = FALSE
)

ElbowPlot(states, reduction = "pca", ndims = 50)

states <- RunHarmony(
  object        = states,
  group.by.vars = "type",
  verbose       = TRUE
)

states <- FindNeighbors(
  states,
  reduction = "harmony",
  dims      = 1:25,
  graph.name = "totalRNA_nn"
)

states <- FindClusters(
  states,
  graph.name = "totalRNA_nn",
  resolution = 2
)
states[["states_nn_alg1"]] <- Idents(states)

states <- RunUMAP(
  states,
  reduction      = "harmony",
  dims           = 1:25,
  reduction.name = "states.umap",
  reduction.key  = "UMAP_",
  n.neighbors    = 50,
  min.dist       = 0.001,
  spread         = 3,
  seed.use       = 42
)

saveRDS(states, file = rds_file)
