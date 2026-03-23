rm(list = ls())

library(Seurat)
library(Matrix)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readr)

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/12_SDEG")


states <- readRDS("../05_plaque/states_with_plaque_info.rds")

target_type <- "14mAD"

ordered_rings <- c("ring1_0_10um", "ring2_10_20um", "ring3_20_30um", "ring4_30_40um", "ring5_40_50um", "outer_50um_plus")

get_ring_label <- function(ring_names) {
  sapply(ring_names, function(x) {
    if (grepl("^ring[0-9]+", x)) {
      sub("^((ring[0-9]+)).*", "\\1", x)
    } else if (grepl("^outer", x)) {
      "outer"
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}


gene_list_path <- "ring30_14mAD_SDEG_totalRNA.csv"
genes_to_keep <- unique(read.csv(gene_list_path, stringsAsFactors = FALSE)[[1]])  # first column as gene name


ad_obj <- subset(states, subset = type == target_type)

rings <- intersect(ordered_rings, unique(as.character(ad_obj@meta.data$ring)))


total_mat <- GetAssayData(ad_obj, assay = "totalRNA", slot = "data")  # 用 data slot 求均值
rb_mat    <- GetAssayData(ad_obj, assay = "rbRNA",    slot = "counts")

total_mat <- as.matrix(total_mat)
rb_mat    <- as(rb_mat, "dgCMatrix")

total_mat <- total_mat[rownames(total_mat) %in% genes_to_keep, , drop = FALSE]
rb_mat <- rb_mat[rownames(rb_mat) %in% rownames(total_mat), , drop = FALSE]  # 保证TE和total行顺序一致

n_gene <- nrow(total_mat)
n_ring <- length(rings)

total_mean_mat <- matrix(
  NA,
  nrow = n_gene,
  ncol = n_ring,
  dimnames = list(rownames(total_mat), rings)
)
TE_mat <- matrix(
  NA,
  nrow = n_gene,
  ncol = n_ring,
  dimnames = list(rownames(total_mat), rings)
)

total_mat_counts <- GetAssayData(ad_obj, assay = "totalRNA", slot = "counts")
total_mat_counts <- as(total_mat_counts, "dgCMatrix")
total_mat_counts <- total_mat_counts[rownames(total_mat), , drop = FALSE]  

rb_mat <- rb_mat[rownames(total_mat_counts), , drop = FALSE]  

for (r in rings) {
  cells_r <- rownames(ad_obj@meta.data)[as.character(ad_obj@meta.data$ring) == r]
  if (length(cells_r) == 0) next

  mean_total <- rowMeans(total_mat[, cells_r, drop = FALSE])
  sum_total_counts <- Matrix::rowSums(total_mat_counts[, cells_r, drop = FALSE])
  sum_rb    <- Matrix::rowSums(rb_mat[, cells_r, drop = FALSE])

  total_mean_mat[, r] <- mean_total
  TE_mat[, r] <- sum_rb / sum_total_counts
}


TE_mat[!is.finite(TE_mat)] <- NA
total_mean_mat[!is.finite(total_mean_mat)] <- NA

valid_genes <- (rowSums(is.na(total_mean_mat)) == 0) & (rowSums(is.na(TE_mat)) == 0)
total_mean_mat <- total_mean_mat[valid_genes, , drop = FALSE]
TE_mat <- TE_mat[valid_genes, , drop = FALSE]

message("Genes retained for heatmap: ", nrow(total_mean_mat))


mat_total <- total_mean_mat
mat_te <- TE_mat

mat_total_scaled <- t(scale(t(mat_total)))
mat_te_scaled <- t(scale(t(mat_te)))

show_rings <- ordered_rings[ordered_rings %in% colnames(mat_total_scaled)]
mat_total_scaled <- mat_total_scaled[, show_rings, drop = FALSE]
mat_te_scaled <- mat_te_scaled[, show_rings, drop = FALSE]

ring_labels <- get_ring_label(colnames(mat_total_scaled))
colnames(mat_total_scaled) <- ring_labels
colnames(mat_te_scaled) <- ring_labels


if (nrow(mat_total_scaled) < 2) {
  warning("Too few genes for clustering in 14mAD")
} else {
  cluster_data <- cbind(mat_total_scaled, mat_te_scaled)
  gene_dist <- dist(cluster_data)
  gene_hclust <- hclust(gene_dist, method = "ward.D2")
  clusters <- cutree(gene_hclust, k = min(4, nrow(mat_total_scaled)))
  row_split <- factor(clusters, levels = sort(unique(clusters)))

  gene_order <- rownames(cluster_data)[gene_hclust$order]
  mat_total_scaled <- mat_total_scaled[gene_order, , drop = FALSE]
  mat_te_scaled <- mat_te_scaled[gene_order, , drop = FALSE]
  row_split <- row_split[gene_order]


  zlim <- c(-1, 0, 1)

  col_fun_total <- colorRamp2(
    zlim,
    c("#2166ac", "#ffffff", "#b2182b")   
  )

  col_fun_te <- colorRamp2(
    zlim,
    c("#2166ac", "#ffffff", "#b2182b")  
  )


  ht_total <- Heatmap(
    matrix = mat_total_scaled,
    name = "totalRNA",
    col = col_fun_total,
    show_row_names = TRUE,    
    row_names_gp = gpar(fontsize = 6),  
    cluster_rows = FALSE,          
    cluster_columns = FALSE,
    column_names_rot = 45,
    column_title = "totalRNA (14mAD)",
    column_title_side = "top",
    row_split = row_split,
    row_gap = unit(1.5, "mm"),
    column_gap = unit(4, "mm"),
    border = TRUE
  )

  ht_te <- Heatmap(
    matrix = mat_te_scaled,
    name = "TE",
    col = col_fun_te,
    show_row_names = TRUE,    
    row_names_gp = gpar(fontsize = 6),  
    cluster_rows = FALSE,         
    cluster_columns = FALSE,
    column_names_rot = 45,
    column_title = "TE (14mAD)",
    column_title_side = "top",
    row_split = row_split,
    row_gap = unit(1.5, "mm"),
    column_gap = unit(4, "mm"),
    border = TRUE
  )

  ht_list <- ht_total + ht_te


  pdf("totalRNA_TE_ring_heatmap_14mAD_ring1-3SDEG_genelabel_split.pdf", width = 4, height = 8)
  draw(ht_list, heatmap_legend_side = "right")
  dev.off()


  gene_cluster_df <- data.frame(
    gene = rownames(mat_te_scaled),
    cluster = row_split
  )

  write.csv(
    gene_cluster_df,
    "totalRNA_TE_ring_gene_clusters_14mAD_ring1-3SDEG_split.csv",
    row.names = FALSE
  )

  message("Finished 14mAD")
}
