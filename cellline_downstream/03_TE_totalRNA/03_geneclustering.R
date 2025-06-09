rm(list=ls())

library(Seurat)
library(anndata)
library(reticulate)
library(ComplexHeatmap)
library(circlize)

use_python("/home/lingyuan2/mambaforge/envs/bigfish_env/bin/python", required = TRUE)
py_config()
sc_ad <- read_h5ad("filtered_data.h5ad")

totalRNA_matrix <- t(sc_ad$X)
rbRNA_matrix <- t(sc_ad$layers['rbRNA'])
metadata <- sc_ad$obs

up_genes <- read.csv("early_upregulated_genes.csv", header = FALSE)

selected_genes_df <- rbind(up_genes)
colnames(selected_genes_df) <- c("gene")
selected_genes <- selected_genes_df$gene

totalRNA_matrix <- totalRNA_matrix[selected_genes, , drop = FALSE]
rbRNA_matrix <- rbRNA_matrix[selected_genes, , drop = FALSE]

totalRNA <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata)
rbRNA <- CreateSeuratObject(counts = rbRNA_matrix, meta.data = metadata)

samples <- unique(totalRNA@meta.data$sample)

expr_matrix <- data.frame(row.names = selected_genes)

for(sample_name in samples) {
  sample_cells <- rownames(totalRNA@meta.data)[totalRNA@meta.data$sample == sample_name]
  totalRNA_counts <- GetAssayData(totalRNA[, sample_cells], slot = "counts")
  rbRNA_counts <- GetAssayData(rbRNA[, sample_cells], slot = "counts")
  totalRNA_counts <- as(totalRNA_counts, "dgCMatrix")
  rbRNA_counts <- as(rbRNA_counts, "dgCMatrix")
  sumtotalRNA <- rowSums(totalRNA_counts)
  sumrbRNA <- rowSums(rbRNA_counts)
  originalTE <- sumrbRNA / sumtotalRNA
  temp_totalRNA <- CreateSeuratObject(counts = totalRNA_counts)
  temp_totalRNA <- NormalizeData(temp_totalRNA, normalization.method = "RC", scale.factor = 2203)
  normalizedtotalRNAExpr <- GetAssayData(temp_totalRNA, slot = "data")
  normalizedrbRNAExpr <- normalizedtotalRNAExpr * originalTE
  finaltotalRNA <- rowMeans(log1p(normalizedtotalRNAExpr))
  finalrbRNA <- rowMeans(log1p(normalizedrbRNAExpr))
  finalTE <- log1p(originalTE)
  expr_matrix[[paste0("totalRNA_", sample_name)]] <- finaltotalRNA
  expr_matrix[[paste0("rbRNA_", sample_name)]] <- finalrbRNA
  expr_matrix[[paste0("te_", sample_name)]] <- finalTE
}

totalRNA_cols <- grep("^totalRNA_", colnames(expr_matrix), value = TRUE)
rbRNA_cols <- grep("^rbRNA_", colnames(expr_matrix), value = TRUE)
te_cols <- grep("^te_", colnames(expr_matrix), value = TRUE)

all_totalRNA <- as.matrix(expr_matrix[, totalRNA_cols])
all_rbRNA <- as.matrix(expr_matrix[, rbRNA_cols])
all_te <- as.matrix(expr_matrix[, te_cols])

scaled_totalRNA <- t(scale(t(all_totalRNA)))
scaled_rbRNA <- t(scale(t(all_rbRNA)))
scaled_te <- t(scale(t(all_te)))

cluster_data <- cbind(scaled_totalRNA, scaled_te)
rownames(cluster_data) <- rownames(expr_matrix)

gene_dist <- dist(cluster_data)
gene_hclust <- hclust(gene_dist)
gene_order <- rownames(cluster_data)[gene_hclust$order]

ordered_totalRNA <- scaled_totalRNA[gene_order, , drop = FALSE]
ordered_rbRNA <- scaled_rbRNA[gene_order, , drop = FALSE]
ordered_te <- scaled_te[gene_order, , drop = FALSE]

col_fun <- colorRamp2(
  quantile(c(ordered_totalRNA, ordered_te, ordered_rbRNA), probs = c(0.1, 0.5, 0.9)),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

pdf("clustered_heatmaps_consistent_162gene_order.pdf", width = 15, height = 15)

ht_totalRNA <- Heatmap(
  matrix = ordered_totalRNA,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  col = col_fun,
  name = "totalRNA",
  column_title = "totalRNA Expression",
  column_title_gp = gpar(fontsize = 16)
)

ht_te <- Heatmap(
  matrix = ordered_te,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  col = col_fun,
  name = "TE",
  column_title = "TE",
  column_title_gp = gpar(fontsize = 16)
)

ht_rbRNA <- Heatmap(
  matrix = ordered_rbRNA,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  col = col_fun,
  name = "rbRNA",
  column_title = "rbRNA Expression",
  column_title_gp = gpar(fontsize = 16)
)

draw(ht_totalRNA + ht_te + ht_rbRNA, 
     column_title = "Gene Expression Patterns (Clustered by totalRNA & TE)", 
     column_title_gp = gpar(fontsize = 20),
     heatmap_legend_side = "bottom",
     gap = unit(5, "mm"))

dev.off()
