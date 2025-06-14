library(Seurat)
library(anndata)
library(reticulate)
library(dplyr)
library(matrixStats) 
library(ggplot2)

rm(list = ls())
setwd("/storage/lingyuan2/STATES_data")
use_python("/home/lingyuan2/mambaforge/envs/bigfish_env/bin/python", required = TRUE)
py_config()
sc_ad <- read_h5ad("cellcyclescore.h5ad")

totalRNA_matrix <- t(sc_ad$X)
rbRNA_matrix <- t(sc_ad$layers['rbRNA'])
metadata <- sc_ad$obs
metadata['phase']
metadata['sample']
totalRNA <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata)

totalRNA[["rbRNA"]] <- CreateAssayObject(counts = rbRNA_matrix)

DefaultAssay(totalRNA) <- "RNA"

# Normalize RNA counts (no log)
totalRNA <- NormalizeData(totalRNA, assay = "RNA", normalization.method = "RC", scale.factor = 2203)
norm_rna <- GetAssayData(totalRNA, slot = "data", assay = "RNA")
te_matrix <- t(sc_ad$layers['TE'])

if(!all(dim(norm_rna) == dim(te_matrix))) {
  stop("Normalized RNA matrix and TE matrix dimensions do not match!")
}

norm_rbRNA <- norm_rna * te_matrix
norm_rbRNA[is.nan(norm_rbRNA)] <- 0
totalRNA[["rbRNA"]]@data <- as.matrix(norm_rbRNA)

# Log1p transform
totalRNA <- SetAssayData(totalRNA, assay = "RNA", slot = "data", new.data = log1p(norm_rna))
totalRNA <- SetAssayData(totalRNA, assay = "rbRNA", slot = "data", new.data = log1p(totalRNA[["rbRNA"]]@data))

# Select highly variable genes (top 700)
totalRNA <- FindVariableFeatures(totalRNA, selection.method = "vst", nfeatures = 700)

# Regress out cell cycle
totalRNA <- ScaleData(totalRNA, features = rownames(totalRNA), assay = "RNA", vars.to.regress = c("S_score", "G2M_score"))
totalRNA <- RunPCA(totalRNA, reduction.name = "pca")

DefaultAssay(totalRNA) <- "rbRNA"
totalRNA <- ScaleData(totalRNA, features = rownames(totalRNA), assay = "rbRNA", vars.to.regress = c("S_score", "G2M_score"))
totalRNA <- FindVariableFeatures(totalRNA, selection.method = "vst", nfeatures = 700)
totalRNA <- RunPCA(totalRNA, reduction.name = "rbRNA_pca")

totalRNA <- FindMultiModalNeighbors(
  totalRNA, reduction.list = list("pca", "rbRNA_pca"),
  dims.list = list(1:30, 1:20), modality.weight.name = "RNA.weight"
)

totalRNA <- RunUMAP(totalRNA, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                n.neighbors = 30, min.dist = 1, spread = 3)

totalRNA <- FindClusters(totalRNA, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

p1_sample <- DimPlot(totalRNA, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 2.5)
p1_sample

# Save weights before metadata update
RNA_weight <- totalRNA$RNA.weight
rbRNA_weight <- totalRNA$rbRNA.weight

print("Before updating metadata:")
print(head(totalRNA@meta.data))

metadata$sample <- factor(metadata$sample, levels = c("C3control", "B4Tg15min", "B5Tg30min", "B6Tg1h", "C4Tg2h", "C5Tg4h"))
metadata$phase <- factor(metadata$phase, levels = c("G1", "S", "G2M"))

totalRNA@meta.data <- metadata
totalRNA$RNA.weight <- RNA_weight
totalRNA$rbRNA.weight <- rbRNA_weight

print("After updating metadata:")
print(head(totalRNA@meta.data))

totalRNA$sample <- factor(totalRNA$sample, levels = c("C3control", "B4Tg15min", "B5Tg30min", "B6Tg1h", "C4Tg2h", "C5Tg4h"))

print("Sample column after factor conversion:")
print(table(totalRNA$sample))

print("Checking weights after metadata update:")
print(colnames(totalRNA@meta.data))

output_dir <- "/storage/lingyuan2/STATES_data"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

add_sample_legend <- function(plot) {
  plot + theme(legend.position = "right") +
    guides(color = guide_legend(title = "Sample", override.aes = list(size = 5)))
}

add_phase_legend <- function(plot) {
  plot + theme(legend.position = "right") +
    guides(color = guide_legend(title = "Phase", override.aes = list(size = 5)))
}

# Color schemes
sample_colors <- c(
  "C3control" = "#FDE725FF",
  "B4Tg15min" = "#7AD151FF",
  "B5Tg30min" = "#22A884FF",
  "B6Tg1h" = "#2A788EFF",
  "C4Tg2h" = "#414487FF",
  "C5Tg4h" = "#440154FF"
)

phase_colors <- c(
  "G1" = "#1f77b4",
  "S" = "#2ca02c",
  "G2M" = "#ff7f0e"
)

# WNN UMAP by sample
p1_sample <- DimPlot(totalRNA, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 2.5, 
                     group.by = "sample", cols = sample_colors) + labs(subtitle = "WNN by Sample") + coord_fixed() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
p1_sample <- add_sample_legend(p1_sample)
ggsave(filename = file.path(output_dir, "wnn_by_sample.pdf"), plot = p1_sample, width = 10, height = 5)
p1_sample

p1_phase <- DimPlot(totalRNA, reduction = 'wnn.umap', label = FALSE, repel = TRUE, label.size = 2.5, 
                    group.by = "phase", cols = phase_colors) + labs(subtitle = "WNN by Phase")
p1_phase <- add_phase_legend(p1_phase)
ggsave(filename = file.path(output_dir, "wnn_by_phase.pdf"), plot = p1_phase, width = 9, height = 12)
p1_phase

# totalRNA UMAP
totalRNA <- RunUMAP(totalRNA, reduction = 'pca', dims = 1:30, assay = 'RNA',
                 reduction.name = 'rna.umap', reduction.key = 'totalRNAUMAP_',
                n.neighbors = 30, min.dist = 2, spread = 3)
p3_sample <- DimPlot(
  totalRNA,
  reduction = 'rna.umap',
  label = FALSE,
  repel = TRUE,
  label.size = 2.5,
  group.by = "sample",
  cols = sample_colors
) +
  labs(subtitle = "totalRNA by Sample") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p3_sample <- DimPlot(totalRNA, reduction = 'rna.umap', label = FALSE, repel = TRUE, label.size = 2.5, 
                     group.by = "sample", cols = sample_colors) + labs(subtitle = "totalRNA by Sample") + coord_fixed() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
p3_sample <- add_sample_legend(p3_sample)
ggsave(filename = file.path(output_dir, "totalRNA_by_sample.pdf"), plot = p3_sample, width = 5.5, height = 5)
p3_sample

p3_phase <- DimPlot(totalRNA, reduction = 'rna.umap', label = FALSE, repel = TRUE, label.size = 2.5, 
                    group.by = "phase", cols = phase_colors) + labs(subtitle = "totalRNA by Phase")
p3_phase <- add_phase_legend(p3_phase)
ggsave(filename = file.path(output_dir, "totalRNA_by_phase.pdf"), plot = p3_phase, width = 9, height = 12)
p3_phase 

# rbRNA UMAP
totalRNA <- RunUMAP(totalRNA, reduction = 'rbRNA_pca', dims = 1:20, assay = 'rbRNA',
                 reduction.name = 'rbRNA.umap', reduction.key = 'rbRNAUMAP_',
                n.neighbors = 30, min.dist = 2, spread = 3)
p4_sample <- DimPlot(
  totalRNA,
  reduction = 'rbRNA.umap',
  label = FALSE,
  repel = TRUE,
  label.size = 2.5,
  group.by = "sample",
  cols = sample_colors
) +
  labs(subtitle = "rbRNA by Sample") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p4_sample <- DimPlot(totalRNA, reduction = 'rbRNA.umap', label = FALSE, repel = TRUE, label.size = 2.5, 
                     group.by = "sample", cols = sample_colors) + labs(subtitle = "rbRNA by Sample") + coord_fixed() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
p4_sample <- add_sample_legend(p4_sample)
ggsave(filename = file.path(output_dir, "rbRNA_by_sample.pdf"), plot = p4_sample, width = 5.5, height = 5)
p4_sample

p4_phase <- DimPlot(totalRNA, reduction = 'rbRNA.umap', label = FALSE, repel = TRUE, label.size = 2.5, 
                    group.by = "phase", cols = phase_colors) + labs(subtitle = "rbRNA by Phase")
p4_phase <- add_phase_legend(p4_phase)
ggsave(filename = file.path(output_dir, "rbRNA_by_phase.pdf"), plot = p4_phase, width = 9, height = 12)

# Save UMAP embeddings
wnn_umap_embeddings <- Embeddings(totalRNA, reduction = "wnn.umap")
write.csv(wnn_umap_embeddings, file = file.path(output_dir, "wnn_umap_embeddings.csv"), row.names = TRUE)

rna_umap_embeddings <- Embeddings(totalRNA, reduction = "rna.umap")
write.csv(rna_umap_embeddings, file = file.path(output_dir, "totalRNA_umap_embeddings.csv"), row.names = TRUE)

rbRNA_umap_embeddings <- Embeddings(totalRNA, reduction = "rbRNA.umap")
write.csv(rbRNA_umap_embeddings, file = file.path(output_dir, "rbRNA_umap_embeddings.csv"), row.names = TRUE)

library(monocle3)
ls("package:monocle3")

# Create Monocle3 object
cds <- new_cell_data_set(
  expression_data = GetAssayData(totalRNA, assay = "rbRNA", slot = "counts"),
  cell_metadata = totalRNA@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(totalRNA), row.names = rownames(totalRNA))
)
# Use Seurat UMAP embedding
reducedDims(cds)$UMAP <- Embeddings(totalRNA, reduction = "wnn.umap")
set.seed(2025)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

# Pseudotime plot
p_pseudotime <- plot_cells(cds, color_cells_by = "pseudotime")
library(ggplot2)

p_pseudotime <- p_pseudotime +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
ggsave(filename = file.path(output_dir, "pseudotime_initial_plot.pdf"), plot = p_pseudotime, width = 5.5, height = 5)

library(ggplot2)
library(RColorBrewer)

plot_cells(cds, show_trajectory_graph = TRUE, color_cells_by = "pseudotime", cell_size = 0.6) + 
  coord_fixed() +
  scale_color_gradientn(colors = brewer.pal(11, "Blues")[1:9]) + 
  theme_classic() + 
  theme(
    line = element_line(linewidth = 0.3),
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    aspect.ratio = 1 / 1
  ) +
  labs(color = "Pseudotime") -> p2
ggsave("pseudotime_plot_size2.pdf", plot = p2, width = 5.5, height = 5)

# Reverse sample order for ridge plot if needed
#plot_data$sample <- factor(plot_data$sample, levels = rev(unique(plot_data$sample)))
#plot_data$sample <- factor(plot_data$sample, levels = rev(levels(plot_data$sample)))

pseudotime_values <- pseudotime(cds)
pseudotime_data <- data.frame(
  Cell = names(pseudotime_values),
  Pseudotime = pseudotime_values
  )
write.csv(pseudotime_data, file = file.path(output_dir, "pseudotime_values2.csv"), row.names = FALSE)
library(ggridges)
plot_data <- as.data.frame(cds@colData)
plot_data["pseudotime"] = cds@principal_graph_aux@listData$UMAP$pseudotime

p <- ggplot(plot_data, aes(x = pseudotime, y = sample, color = sample, fill = sample)) +
  geom_density_ridges(alpha = 0.6, size = 0.6) +
  labs(x = "Pseudotime", y = NULL) +
  theme_classic() +
  theme(
    line = element_line(linewidth = 0.3),
    text = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    axis.title = element_blank()
  ) +
  scale_color_manual(values = sample_colors) +
  scale_fill_manual(values = sample_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_discrete(limits = rev(levels(plot_data$sample)))

print(p)
ggsave(filename = file.path(output_dir, "pseudotime_distribution_ridge_plot2.pdf"), p, width = 6, height = 6)
