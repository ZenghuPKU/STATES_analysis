## ===============================
## environment
## ===============================
rm(list = ls())

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/12_SDEG")

states <- readRDS("../05_plaque/states_with_plaque_info.rds")

## ===============================
## 14mAD
## ===============================

obj <- subset(states, subset = type == "14mAD")
tp <- "14mAD"

message("Processing type: ", tp)

totalRNA_counts_mat <- GetAssayData(obj, assay = "totalRNA", slot = "counts")
totalRNA_data_mat   <- GetAssayData(obj, assay = "totalRNA", slot = "data")

totalRNA_counts_mat <- as(totalRNA_counts_mat, "dgCMatrix")
totalRNA_data_mat   <- as.matrix(totalRNA_data_mat)

meta <- obj@meta.data

# ====================== inner, outer ======================
inner_rings <- c("ring1_0_10um", "ring2_10_20um", "ring3_20_30um")
outer_rings <- c("ring4_30_40um", "ring5_40_50um", "outer_50um_plus")

inner_cells <- rownames(meta)[meta$ring %in% inner_rings]
outer_cells <- rownames(meta)[meta$ring %in% outer_rings]

message("Inner cells: ", length(inner_cells), ", Outer cells: ", length(outer_cells))

if (length(inner_cells) < 10 | length(outer_cells) < 10) {
  stop("Not enough cells in inner or outer group.")
}

## --------------------------------
## gene filter
## --------------------------------
total_inner_counts <- totalRNA_counts_mat[, inner_cells, drop = FALSE]
valid_gene <- Matrix::rowSums(total_inner_counts >= 1) >= 10
genes_use <- rownames(totalRNA_counts_mat)[valid_gene]

message("Genes tested after inner ring filter: ", length(genes_use))

## Wilcoxon test
res_list <- lapply(genes_use, function(g) {
  n_inner  <- length(inner_cells)
  n_outer  <- length(outer_cells)

  expr_inner  <- as.numeric(totalRNA_data_mat[g, inner_cells])
  expr_outer  <- as.numeric(totalRNA_data_mat[g, outer_cells])

  wt <- wilcox.test(expr_inner, expr_outer, exact = FALSE)

  mean_inner <- mean(expr_inner)
  mean_outer <- mean(expr_outer)

  pseudo <- 1e-6
  log2FC <- log2((mean_inner + pseudo) / (mean_outer + pseudo))

  data.frame(
    gene        = g,
    group1      = "inner",
    group2      = "outer",
    type        = tp,
    n_inner     = n_inner,
    n_outer     = n_outer,
    mean_inner  = mean_inner,
    mean_outer  = mean_outer,
    delta_mean  = mean_inner - mean_outer,
    log2FC      = log2FC,
    p_value     = wt$p.value,
    stringsAsFactors = FALSE
  )
})

res_df <- bind_rows(res_list)

if (nrow(res_df) == 0) {
  stop("No genes left after inner ring filter.")
}

res_df$padj <- p.adjust(res_df$p_value, method = "BH")

## save
out_file <- paste0(
  "totalRNA_Wilcox_dataSlot_countsFilter_",
  tp, "_inner_vs_outer_30ring_14mAD.csv"
)
write.csv(res_df, out_file, row.names = FALSE)
message("Results saved: ", out_file)
