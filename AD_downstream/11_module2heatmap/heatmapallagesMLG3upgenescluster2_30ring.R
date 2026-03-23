rm(list=ls())

library(Seurat)
library(dplyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/11_module2heatmap")
states_rds <- "../05_plaque/states_with_plaque_info.rds"
states <- readRDS(states_rds)

selected_genes_file <- "MLGupgenescluster2.txt"
selected_genes <- readLines(selected_genes_file)
selected_genes <- unique(selected_genes[selected_genes != ""])

inner_rings <- c("ring1_0_10um", "ring2_10_20um", "ring3_20_30um")
outer_rings <- c("ring4_30_40um", "ring5_40_50um", "outer_50um_plus")

mlg_cells <- subset(states, subset = (states_nn_alg1_label2 == "MLG") & (group == "AD"))
meta <- mlg_cells@meta.data

required_cols <- c("ring_assignment", "states_nn_alg1_label3")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "None of the requested variables were found in meta.data: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Assign inner/outer
meta$ring_group <- dplyr::case_when(
  meta$ring_assignment %in% inner_rings ~ "inner",
  meta$ring_assignment %in% outer_rings ~ "outer",
  TRUE ~ NA_character_
)

meta_valid <- meta[!is.na(meta$ring_group), , drop=FALSE]

if (nrow(meta_valid) == 0) {
  stop("No valid cells found with defined inner/outer ring_group assignments.")
}

inner_cells <- rownames(meta_valid)[meta_valid$ring_group == "inner"]
if (length(inner_cells) == 0) {
  stop("No inner cells found for analysis.")
}

# Slot extraction
totalRNA_counts_all <- tryCatch({
  GetAssayData(mlg_cells, assay = "totalRNA", slot = "counts")
}, error = function(e) {
  stop("Cannot retrieve 'totalRNA' counts slot: ", e$message)
})

genes_can_use <- intersect(selected_genes, rownames(totalRNA_counts_all))
if (length(genes_can_use) == 0) stop("No selected genes found in totalRNA counts!")

ordered_genes <- selected_genes[selected_genes %in% genes_can_use]

if (length(ordered_genes) == 0) {
  stop("No ordered genes found for plotting.")
}

totalRNA_counts <- totalRNA_counts_all[ordered_genes, inner_cells, drop=FALSE]
gene_filter <- rowSums(totalRNA_counts >= 1) / ncol(totalRNA_counts) > 0.05
genes_keep <- names(gene_filter)[gene_filter]

if (length(genes_keep) == 0) {
  stop("No genes found passing the filter (>5% of inner cells with totalRNA count >= 1) from selected gene list.")
}

ordered_genes_final <- ordered_genes[ordered_genes %in% genes_keep]
if (length(ordered_genes_final) == 0) {
  stop("No genes remain after filtering by expression.")
}

# Build group variable (label3 + ring_group)
meta_valid$group2 <- paste(
  meta_valid$states_nn_alg1_label3,
  meta_valid$ring_group,
  sep = "_"
)

desired_group_order <- c(
  "MLG1_inner", "MLG2_inner", "MLG3_inner",
  "MLG1_outer", "MLG2_outer", "MLG3_outer"
)

# Try extract other assay slots for selected genes, catch errors if slots do not exist
totalRNA_data <- tryCatch({
  GetAssayData(mlg_cells, assay = "totalRNA", slot = "data")[ordered_genes_final, , drop=FALSE]
}, error = function(e) {
  stop("Cannot retrieve 'totalRNA' data slot: ", e$message)
})
totalRNA_counts_all <- totalRNA_counts_all[ordered_genes_final, , drop=FALSE]

rbRNA_counts_all <- tryCatch({
  GetAssayData(mlg_cells, assay = "rbRNA", slot = "counts")[ordered_genes_final, , drop=FALSE]
}, error = function(e) {
  stop("Cannot retrieve 'rbRNA' counts slot: ", e$message)
})

all_groups_in_data <- unique(meta_valid$group2)
all_groups_in_order <- desired_group_order[desired_group_order %in% all_groups_in_data]
if (length(all_groups_in_order) == 0) {
  stop("None of the desired column groups exist in your data!")
}

totalRNA_avg <- matrix(NA, nrow=length(ordered_genes_final), ncol=length(desired_group_order),
                      dimnames=list(ordered_genes_final, desired_group_order))
TE_avg <- matrix(NA, nrow=length(ordered_genes_final), ncol=length(desired_group_order),
                 dimnames=list(ordered_genes_final, desired_group_order))

for (g in desired_group_order) {
  if(!(g %in% all_groups_in_data)) next
  group_cells <- rownames(meta_valid)[meta_valid$group2 == g]
  if(length(group_cells) == 0) next

  totalRNA_avg[, g] <- if (length(group_cells) > 1) {
    rowMeans(totalRNA_data[, group_cells, drop=FALSE])
  } else {
    as.vector(totalRNA_data[, group_cells, drop=FALSE])
  }
  rbRNA_sums <- Matrix::rowSums(rbRNA_counts_all[, group_cells, drop=FALSE])
  totalRNA_counts_sums <- Matrix::rowSums(totalRNA_counts_all[, group_cells, drop=FALSE])
  TE_g <- ifelse(totalRNA_counts_sums > 0, rbRNA_sums / totalRNA_counts_sums, NA_real_)
  TE_g[!is.finite(TE_g)] <- NA
  TE_avg[, g] <- TE_g
}

# Remove genes with missing values everywhere
valid_genes <- (rowSums(is.na(totalRNA_avg[, all_groups_in_order, drop=FALSE])) == 0) & 
               (rowSums(is.na(TE_avg[, all_groups_in_order, drop=FALSE])) == 0)
totalRNA_avg <- totalRNA_avg[valid_genes, , drop=FALSE]
TE_avg <- TE_avg[valid_genes, , drop=FALSE]

if (nrow(totalRNA_avg) == 0 || nrow(TE_avg) == 0) {
  stop("No valid genes with sufficient data after filtering NAs from totalRNA and TE matrices.")
}
cols_to_plot <- desired_group_order[desired_group_order %in% colnames(totalRNA_avg)]
totalRNA_avg <- totalRNA_avg[, cols_to_plot, drop=FALSE]
TE_avg <- TE_avg[, cols_to_plot, drop=FALSE]

group_info <- data.frame(
  group = cols_to_plot,
  label3 = sapply(strsplit(cols_to_plot, "_"), `[`, 1),
  ring_group = sapply(strsplit(cols_to_plot, "_"), `[`, 2),
  stringsAsFactors = FALSE
)

ring_group_palette <- c(inner="#1F77B4", outer="#FF7F0E")
label3_palette <- c(MLG1="#A1D99B", MLG2="#FC9272", MLG3="#9ECAE1")
column_anno = HeatmapAnnotation(
  df = group_info[, c("ring_group", "label3")],
  col = list(
    ring_group = ring_group_palette,
    label3 = label3_palette
  )
)

zscore_clip_global <- function(mat, z=1) {
  mat_z <- t(scale(t(mat)))
  mat_z[mat_z >  z] <-  z
  mat_z[mat_z < -z] <- -z
  mat_z
}
mat_total_scaled <- zscore_clip_global(totalRNA_avg, 1)
mat_te_scaled <- zscore_clip_global(TE_avg, 1)

zlim <- c(-1,0,1)
col_fun_total <- colorRamp2(zlim, c("#2166ac", "#ffffff", "#b2182b"))
col_fun_te <- colorRamp2(zlim, c("#2166ac", "#ffffff", "#b2182b"))

ht_total <- Heatmap(
  matrix = mat_total_scaled,
  name = "totalRNA",
  col = col_fun_total,
  show_row_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_title = "totalRNA",
  column_title_side = "top",
  row_names_side = "left",
  row_gap = unit(1.5, "mm"),
  column_gap = unit(3, "mm"),
  border = TRUE,
  top_annotation = column_anno
)
ht_te <- Heatmap(
  matrix = mat_te_scaled,
  name = "TE",
  col = col_fun_te,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_title = "TE",
  column_title_side = "top",
  row_gap = unit(1.5, "mm"),
  column_gap = unit(3, "mm"),
  border = TRUE,
  top_annotation = column_anno
)
ht_list <- ht_total + ht_te

pdf("MLGupgenescluster2_totalRNA_TE_inner_outer_heatmap_allAD_30ring.pdf", width=7, height=6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()
