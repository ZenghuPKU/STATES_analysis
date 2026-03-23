# 08_MLG_Slingshot
rm(list = ls())
library(Seurat)
library(tidyverse)
library(destiny)        
library(SingleCellExperiment)
library(slingshot)
library(viridis)
library(ggplot2)
library(plyr)
library(fields)     
setwd("~/AD_downstream/")

states <- readRDS("states_with_plaque_info.rds")
DefaultAssay(states) <- "totalRNA"

mlg_lineage <- c("MLG1", "MLG2", "MLG3")

mlg_cells <- rownames(states@meta.data)[
  states$states_nn_alg1_label3 %in% mlg_lineage
]

mlg_subset <- subset(states, cells = mlg_cells)

mlg_subset$MLG_subtype <- mlg_subset$states_nn_alg1_label3
Idents(mlg_subset) <- factor(mlg_subset$MLG_subtype,
                             levels = mlg_lineage)

mlg_levels <- levels(Idents(mlg_subset))
print(mlg_levels)
table(Idents(mlg_subset))

mlg_subset <- NormalizeData(mlg_subset, verbose = FALSE)
mlg_subset <- FindVariableFeatures(
  mlg_subset,
  selection.method = "vst",
  nfeatures        = 1500,
  verbose          = FALSE
)
mlg_subset <- ScaleData(
  mlg_subset,
  features = VariableFeatures(mlg_subset),
  verbose  = FALSE
)
mlg_subset <- RunPCA(
  mlg_subset,
  features = VariableFeatures(mlg_subset),
  npcs     = 50,
  verbose  = FALSE
)

pdf("MLG_ElbowPlot_PCA.pdf", width = 6, height = 5)
ElbowPlot(mlg_subset, ndims = 50)
dev.off()
pca_embeddings <- Embeddings(mlg_subset, reduction = "pca")

pca_for_dm <- pca_embeddings[, 1:25, drop = FALSE]

dup_flag <- duplicated(pca_for_dm)
cat("Duplicated rows in PC space:", sum(dup_flag), "\n")

pca_for_dm_unique <- pca_for_dm[!dup_flag, , drop = FALSE]

dm <- DiffusionMap(
  pca_for_dm_unique,
  sigma    = "local",
  distance = "euclidean"
)

dm_embeddings <- as.data.frame(eigenvectors(dm))
colnames(dm_embeddings) <- paste0("DC", seq_len(ncol(dm_embeddings)))
rownames(dm_embeddings) <- rownames(pca_for_dm_unique)

sce <- as.SingleCellExperiment(mlg_subset)
sce$cluster <- Idents(mlg_subset)

if (!"type" %in% colnames(colData(sce))) {
  sce$type <- mlg_subset$type
}

dm_mat <- as.matrix(dm_embeddings[, 1:2])
common_cells <- intersect(rownames(dm_mat), colnames(sce))
sce    <- sce[, common_cells]
dm_mat <- dm_mat[common_cells, , drop = FALSE]

reducedDim(sce, "Diffusion") <- dm_mat

diffusion_coords <- reducedDim(sce, "Diffusion")

xlim_all <- range(diffusion_coords[, 1])
ylim_all <- range(diffusion_coords[, 2])

start_cluster <- "MLG1"
end_cluster   <- "MLG3"

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim    = "Diffusion",
  start.clus    = start_cluster,
  end.clus      = end_cluster
)

pt_mat <- slingPseudotime(sce)
head(pt_mat)

pseudotime_lineage1 <- pt_mat[, 1]

pt_min <- min(pseudotime_lineage1, na.rm = TRUE)
pt_max <- max(pseudotime_lineage1, na.rm = TRUE)
pseudotime_scaled <- (pseudotime_lineage1 - pt_min) / (pt_max - pt_min)

pal_pt <- viridis(100, option = "plasma")
cols_pt <- pal_pt[
  cut(pseudotime_scaled, breaks = 100, labels = FALSE)
]

pdf("Diffusion_MLG_pseudotime.pdf", width = 8.5, height = 8)
par(mar = c(5, 4, 4, 5) + 0.1)  

plot(
  diffusion_coords,
  col  = cols_pt,
  pch  = 16,
  cex  = 0.4,
  asp  = 1,
  xlab = "Diffusion Component 1",
  ylab = "Diffusion Component 2",
  main = "MLG lineage: Slingshot pseudotime",
  xlim = xlim_all,
  ylim = ylim_all
)
lines(SlingshotDataSet(sce), lwd = 2)

image.plot(
  legend.only = TRUE,
  zlim        = c(0, 1),
  col         = pal_pt,
  legend.lab  = "Pseudotime (0–1)",
  legend.line = 2
)

dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)  

cols_clu <- c(
  "MLG1" = "#A9B4D8",
  "MLG2" = "darkgreen",
  "MLG3" = "darkorange"
)

pdf("MLG_DiffusionMap_QC_byCluster.pdf", width = 8, height = 8)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(
  diffusion_coords,
  col  = cols_clu[as.character(sce$cluster)],
  pch  = 16,
  cex  = 0.4,
  asp  = 1,
  xlab = "Diffusion Component 1",
  ylab = "Diffusion Component 2",
  main = "MLG lineage: Diffusion Map (by cluster)",
  xlim = xlim_all,
  ylim = ylim_all
)
lines(SlingshotDataSet(sce), lwd = 2)
legend(
  "topright",
  legend = mlg_levels,
  col    = cols_clu[mlg_levels],
  pch    = 16,
  cex    = 0.6
)
dev.off()

genes_to_plot <- c("Trem2","Apoe")

cells_use <- colnames(sce)

totalRNA_norm_mlg <- GetAssayData(mlg_subset, assay = "totalRNA", slot = "data")
totalRNA_norm_use <- totalRNA_norm_mlg[, cells_use, drop = FALSE]

totalRNA_counts_all <- GetAssayData(states, assay = "totalRNA", slot = "counts")
rbRNA_counts_all    <- GetAssayData(states, assay = "rbRNA",    slot = "counts")

totalRNA_counts_use <- totalRNA_counts_all[, cells_use, drop = FALSE]
rbRNA_counts_use    <- rbRNA_counts_all[,    cells_use, drop = FALSE]

library(tidyr)

pt_all <- pseudotime_scaled          
pt_all <- pt_all[colnames(sce)]      

valid_pt    <- is.finite(pt_all)
pt_valid    <- pt_all[valid_pt]
cells_valid <- names(pt_valid)

length(pt_valid)

totalRNA_norm_valid   <- totalRNA_norm_use[,  cells_valid, drop = FALSE]
totalRNA_counts_valid <- totalRNA_counts_use[, cells_valid, drop = FALSE]
rbRNA_counts_valid    <- rbRNA_counts_use[,    cells_valid, drop = FALSE]

n_bins        <- 30   
min_cells_bin <- 20   

bin_id <- cut(
  pt_valid,
  breaks = n_bins,
  labels = FALSE,
  include.lowest = TRUE
)

bin_pt_center <- tapply(pt_valid, bin_id, mean, na.rm = TRUE)

traj_list <- list()

for (gn in genes_to_plot) {
  message("Computing pseudotime trajectory for gene: ", gn)
  
  if (!gn %in% rownames(totalRNA_norm_valid)) {
    warning("  Gene ", gn, " not in totalRNA_norm_valid, skip.")
    next
  }
  if (!gn %in% rownames(totalRNA_counts_valid) ||
      !gn %in% rownames(rbRNA_counts_valid)) {
    warning("  Gene ", gn, " not in raw counts, skip TE.")
    next
  }
  
  for (b in sort(unique(bin_id))) {
    cells_b <- cells_valid[bin_id == b]
    n_b     <- length(cells_b)
    if (n_b < min_cells_bin) next
    
    gene_norm <- as.numeric(totalRNA_norm_valid[gn, cells_b])
    totalRNA_mean_b <- mean(gene_norm, na.rm = TRUE)
    
    gene_total_raw <- as.numeric(totalRNA_counts_valid[gn, cells_b])
    gene_rb_raw    <- as.numeric(rbRNA_counts_valid[gn, cells_b])
    
    sum_totalRNA <- sum(gene_total_raw, na.rm = TRUE)
    sum_rbRNA    <- sum(gene_rb_raw,    na.rm = TRUE)
    
    if (!is.finite(sum_totalRNA) || sum_totalRNA <= 0) {
      original_te <- NA_real_
    } else {
      original_te <- sum_rbRNA / sum_totalRNA
    }
    if (!is.finite(original_te) || original_te < 0) {
      original_te <- NA_real_
    }
    
    log_te <- log1p(original_te)
    
    traj_list[[length(traj_list) + 1]] <- data.frame(
      Gene          = gn,
      Bin           = as.integer(b),
      PT_center     = bin_pt_center[as.character(b)],  # 这里就是 0–1 中心
      totalRNA_mean = totalRNA_mean_b,
      TE_raw        = original_te,
      TE_log1p      = log_te,
      stringsAsFactors = FALSE
    )
  }
}

traj_df <- dplyr::bind_rows(traj_list)
table(traj_df$Gene)

out_dir_traj <- "Pseudotime_Trajectories_dualAxis_smooth"
dir.create(out_dir_traj, showWarnings = FALSE)

for (gn in unique(traj_df$Gene)) {
  df_g <- traj_df[traj_df$Gene == gn, , drop = FALSE]
  if (nrow(df_g) == 0) next
  
  df_g <- df_g %>%
    dplyr::filter(is.finite(totalRNA_mean),
                  is.finite(TE_log1p),
                  is.finite(PT_center))
  if (nrow(df_g) < 3) next
  
  fit_total <- stats::smooth.spline(
    x    = df_g$PT_center,
    y    = df_g$totalRNA_mean,
    spar = 0.5
  )
  fit_te    <- stats::smooth.spline(
    x    = df_g$PT_center,
    y    = df_g$TE_log1p,
    spar = 0.5
  )
  
  x_new <- seq(
    min(df_g$PT_center),
    max(df_g$PT_center),
    length.out = 400
  )
  
  total_smooth    <- predict(fit_total, x_new)$y
  te_log1p_smooth <- predict(fit_te,    x_new)$y
  
  ok <- is.finite(total_smooth) & is.finite(te_log1p_smooth)
  x_new           <- x_new[ok]
  total_smooth    <- total_smooth[ok]
  te_log1p_smooth <- te_log1p_smooth[ok]
  if (length(x_new) < 3) next
  
  df_plot <- data.frame(
    PT_center     = x_new,            
    totalRNA_mean = total_smooth,
    TE_log1p      = te_log1p_smooth
  )
  
  total_min  <- min(df_plot$totalRNA_mean, na.rm = TRUE)
  total_max  <- max(df_plot$totalRNA_mean, na.rm = TRUE)
  te_min     <- min(df_plot$TE_log1p,      na.rm = TRUE)
  te_max     <- max(df_plot$TE_log1p,      na.rm = TRUE)
  total_span <- total_max - total_min
  te_span    <- te_max    - te_min
  if (total_span == 0 || te_span == 0) next
  
  df_plot$TE_scaled <- (df_plot$TE_log1p - te_min) *
    (total_span / te_span) + total_min
  
  inv_trans <- function(y) {
    (y - total_min) * (te_span / total_span) + te_min
  }
  
  p_traj <- ggplot(df_plot, aes(x = PT_center)) +
    geom_line(aes(y = totalRNA_mean, color = "totalRNA"), size = 0.7) +
    geom_line(aes(y = TE_scaled,     color = "TE"),       size = 0.7) +
    scale_y_continuous(
      name = "totalRNA (Seurat data slot, log-normalized, smoothed)",
      sec.axis = sec_axis(
        trans = inv_trans,
        name  = "TE (log1p(sum rbRNA / sum totalRNA), smoothed)"
      )
    ) +
    scale_color_manual(
      values = c(
        "totalRNA" = "#1f78b4",
        "TE"       = "#e31a1c"
      ),
      breaks = c("totalRNA", "TE"),
      labels = c("totalRNA (mean, smooth)", "TE (log1p, smooth)")
    ) +
    theme_classic() +
    theme(
      axis.title.y.right = element_text(color = "#e31a1c"),
      axis.title.y.left  = element_text(color = "#1f78b4"),
      legend.title       = element_blank()
    ) +
    labs(
      title = paste0("Pseudotime trajectory of ", gn, " (MLG lineage)"),
      x     = "Pseudotime (scaled 0–1, Slingshot lineage 1)"
    )
  
  out_pdf <- file.path(
    out_dir_traj,
    paste0("Pseudotime_trajectory_", gn, "_dualAxis_smooth.pdf")
  )
  ggsave(out_pdf, p_traj, width = 8, height = 4)
  
  message("  Saved SMOOTH dual-axis trajectory for ", gn, " : ", out_pdf)
}
