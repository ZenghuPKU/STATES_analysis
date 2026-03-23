## 04_MLG123_marker_dot_plot
rm(list = ls())
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(viridis)
library(Matrix)
setwd("~/AD_downstream/06_Annotated_celltypes_UMAP_visualization/")

states <- readRDS("states_with_plaque_info.rds")

cells_to_keep <- which(
  !(states$states_nn_alg1_label2 %in% c("Mix") |
      states$states_nn_alg1_label3 %in% c("Mix", "TEGLU Mix"))
)
states <- subset(states, cells = colnames(states)[cells_to_keep])

Idents(states) <- "states_nn_alg1_label3"
mlg_levels <- c("MLG1", "MLG2", "MLG3")
cells_mlg123 <- WhichCells(states, idents = mlg_levels)
states <- subset(states, cells = cells_mlg123)

Idents(states) <- "states_nn_alg1_label3"
states@active.ident <- factor(states@active.ident, levels = mlg_levels)
states$states_nn_alg1_label3 <- states@active.ident

markers_mlg <- c(
  "Lyz2","Ccl6","Cd9","Spp1","Ctsl","Cst7","Trem2","Apoe",
  "Sparc","Gpr34",
  "Tmem119","P2ry12"
)

make_dotdata <- function(obj, assay_use, markers){
  DefaultAssay(obj) <- assay_use
  dp <- DotPlot(obj, features = markers) + coord_flip()
  dp[["data"]]
}

plot_from_dot <- function(dd, max_pct, size_range = c(0,6),
                          color_limits = c(-1, 1),
                          title = "") {
  ggplot(dd, aes(x = id, y = features.plot,
                 color = avg.exp.scaled, size = pct.exp)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
          plot.title  = element_text(hjust = 0.5)) +
    scale_color_gradientn(
      colours = c("blue","white","red"),
      limits  = color_limits,
      oob     = scales::squish,
      name    = "expression\n(z-score)"
    ) +
    scale_size(range = size_range, limits = c(0, max_pct),
               name = "% cells") +
    labs(title = title, x = "Cell type", y = "Marker") +
    guides(size = guide_legend(order = 2), color = guide_colorbar(order = 1))
}

Idents(states) <- "states_nn_alg1_label3"
states@active.ident <- factor(states@active.ident, levels = mlg_levels)
states$states_nn_alg1_label3 <- states@active.ident

# totalRNA / rbRNA dot data
dot_totalRNA_mlg <- make_dotdata(states, "totalRNA", markers_mlg)
dot_rbRNA_mlg    <- make_dotdata(states, "rbRNA",    markers_mlg)

max_pct_mlg <- max(c(dot_totalRNA_mlg$pct.exp, dot_rbRNA_mlg$pct.exp))

P_mlg_totalRNA <- plot_from_dot(dot_totalRNA_mlg, max_pct_mlg,
                                color_limits = c(-1, 1),
                                title = "MLG1/2/3 totalRNA")
P_mlg_rbRNA    <- plot_from_dot(dot_rbRNA_mlg,    max_pct_mlg,
                                color_limits = c(-1, 1),
                                title = "MLG1/2/3 rbRNA")

ggsave("MLG123_label3_totalRNA_scaled.pdf", P_mlg_totalRNA, width = 6, height = 6)

totalRNA_counts <- GetAssayData(states, assay = "totalRNA", slot = "counts")
rbRNA_counts    <- GetAssayData(states, assay = "rbRNA",    slot = "counts")

genes_used <- markers_mlg[markers_mlg %in% rownames(totalRNA_counts)]

ct_levels <- mlg_levels
te_mat <- matrix(
  NA_real_,
  nrow = length(genes_used),
  ncol = length(ct_levels),
  dimnames = list(genes_used, ct_levels)
)

for (ct in ct_levels) {
  cells_ct <- WhichCells(states, idents = ct)
  if (length(cells_ct) == 0) next
  
  total_sub <- totalRNA_counts[genes_used, cells_ct, drop = FALSE]
  rb_sub    <- rbRNA_counts[genes_used,    cells_ct, drop = FALSE]
  
  total_sum <- Matrix::rowSums(total_sub)
  rb_sum    <- Matrix::rowSums(rb_sub)
  
  te_vec <- rb_sum / total_sum
  te_vec[!is.finite(te_vec)] <- NA_real_
  
  te_mat[, ct] <- te_vec
}

make_te_dotdata <- function(te_mat, ref_dotdata) {
  df <- as.data.frame(te_mat)
  df$features.plot <- rownames(df)
  df_long <- pivot_longer(df, -features.plot,
                          names_to = "id", values_to = "avg.exp")
  
  gene_order     <- unique(ref_dotdata$features.plot)
  celltype_order <- unique(ref_dotdata$id)
  
  df_long <- df_long %>%
    mutate(
      features.plot = factor(features.plot, levels = gene_order),
      id            = factor(id,            levels = celltype_order)
    ) %>%
    group_by(features.plot) %>%
    mutate(avg.exp.scaled = scale(avg.exp)[,1]) %>%
    ungroup()
  
  df_long
}

plot_te_dot <- function(df_te, totalRNA_dot, max_pct, title = "") {
  df_te <- dplyr::left_join(
    df_te,
    totalRNA_dot[, c("features.plot", "id", "pct.exp")],
    by = c("features.plot", "id")
  )
  
  ggplot(df_te, aes(
    x     = id,
    y     = features.plot,
    color = avg.exp.scaled,
    size  = pct.exp
  )) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x  = element_text(angle = 90, vjust = .5, hjust = 1),
      plot.title   = element_text(hjust = 0.5)
    ) +
    scale_color_gradientn(
      colours = c("blue", "white", "red"),
      limits  = c(-1, 1),   # TE z-score 截 [-1, 1]
      oob     = scales::squish,
      name    = "TE\n(z-score)"
    ) +
    scale_size(range = c(0,6), limits = c(0, max_pct), name = "% cells") +
    labs(title = title, x = "Cell type (label3 MLG)", y = "Marker") +
    guides(size = guide_legend(order = 2), color = guide_colorbar(order = 1))
}

te_dotdata_mlg <- make_te_dotdata(te_mat, dot_totalRNA_mlg)
P_mlg_te_scaled <- plot_te_dot(
  te_dotdata_mlg,
  dot_totalRNA_mlg,
  max_pct = max_pct_mlg,
  title = "MLG1/2/3 TE (sum rbRNA / sum totalRNA)"
)

ggsave("MLG123_label3_TE_scaled.pdf", P_mlg_te_scaled, width = 6, height = 6)
