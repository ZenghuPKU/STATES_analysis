# 04_Marker_gene_dot_plot
# Load libraries and set environment
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(viridis)
library(Matrix)
rm(list = ls())
setwd("~/tissue_downstream/05_celltypes_identification/")

# Load the annotated object (states_celltypes_identification.RData)
load("states_celltypes_identification.RData")

# Remove Mix cells and set groups
cells_to_keep <- which(!(states$states_nn_alg1_label2 %in% c("Mix") |
                           states$states_nn_alg1_label3 %in% c("Mix", "TEGLU Mix")))
states <- subset(states, cells = colnames(states)[cells_to_keep])

# label3 order
Idents(states) <- "states_nn_alg1_label3"
order_l3 <- c("TEGLU CA1","TEGLU CA2","TEGLU CA3","TEGLU L2/3","TEGLU L4","TEGLU L5",
              "TEGLU L5/6","TEGLU L6","TEGLU L6b","DGGRC","MSN","INH_Pvalb","INH_Sst",
              "INH_Cnr1_Vip","DE/MEN","CHO/PEP","AC1","AC2","AC3","OLG1","OLG2",
              "OPC","MLG","CHOR","EPEN","Peri/VEC","VLMC","VSMC")
states@active.ident <- factor(states@active.ident, levels = order_l3)
states$states_nn_alg1_label3 <- states@active.ident

# label2 order
order_l2 <- c("TEPN","INH","DE/MEN","CHO/PEP","AC","OLG","OPC","MLG","CHOR/EPEN","VAS")

# Standardize colors & dot size functions
make_dotdata <- function(obj, assay_use, markers){
  DefaultAssay(obj) <- assay_use
  dp <- DotPlot(obj, features = markers) + coord_flip()
  dp[["data"]]
}

plot_from_dot <- function(dd, max_pct, size_range = c(0,6)){
  ggplot(dd, aes(x = id, y = features.plot,
                 color = avg.exp.scaled, size = pct.exp)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    scale_color_gradientn(colours = c("blue","white","red"),
                          limits = c(-2,2), oob = scales::squish,
                          name = "expression\n(z‑score)") +
    scale_size(range = size_range, limits = c(0, max_pct),
               name = "% cells") +
    guides(size = guide_legend(order = 2), color = guide_colorbar(order = 1))
}

# label3: totalRNA & rbRNA
markers_l3 <- c("Myh11","Myl9","Ptgds","Gjb2","Rgs5","Ly6c1","Vim","Rarres2","Ttr","Enpp2",
                "Hexb","C1qa","Pdgfra","Cacng4","Mbp","Mobp","Mag","Mog","Gfap","Mt2",
                "Apoe","Clu","Slc7a10","Mfge8","Gng8","Cadps2","Prkcd","Ntng1","Cnr1","Vip",
                "Sst","Gad2","Gad1","Pvalb","Ppp1r1b","Adora2a","C1ql2","Cplx3","Rprm",
                "Sncb","Cplx1","Fezf2","Pde1a","Rgs4","Lamp5","Cux2","Chgb","Rgs14","Ppp3r1","Wfs1")

dot_totalRNA_l3  <- make_dotdata(states, "totalRNA",  markers_l3)
dot_rbRNA_l3  <- make_dotdata(states, "rbRNA", markers_l3)
max_pct_l3   <- max(c(dot_totalRNA_l3$pct.exp, dot_rbRNA_l3$pct.exp))

P_l3_totalRNA  <- plot_from_dot(dot_totalRNA_l3, max_pct_l3)
P_l3_rbRNA  <- plot_from_dot(dot_rbRNA_l3, max_pct_l3)

ggsave("label3_totalRNA.pdf",  P_l3_totalRNA, width = 10, height = 15)
ggsave("label3_rbRNA.pdf",  P_l3_rbRNA, width = 10, height = 15)

# label2: totalRNA & rbRNA
Idents(states) <- "states_nn_alg1_label2"
states@active.ident <- factor(states@active.ident, levels = order_l2)
states$states_nn_alg1_label2 <- states@active.ident

markers_l2 <- c("Flt1","Ly6c1","Ly6e","Ptgds","Rgs5","Rarres2","Clic6","Folr1","Enpp2","Ttr",
                "Ctss","C1qc","Csf1r","C1qa","Hexb","Cacng4","Olig2","Sox10","Pdgfra","Ptprz1",
                "Trf","Cnp","Fth1","Mobp","Plp1","Clu","Mt2","Atp1b2","Gja1","Aldoc",
                "Resp18","Hap1","Cadps2","Gng8","Tac2","Tubb5","Ntng1","Prkcd","Sparc","Sparcl1",
                "Npy","Sst","Pvalb","Gad2","Gad1","Hpca","Nell2","Mapk1","Ppp3r1","Nrgn")

dot_totalRNA_l2 <- make_dotdata(states, "totalRNA",  markers_l2)
dot_rbRNA_l2 <- make_dotdata(states, "rbRNA", markers_l2)
max_pct_l2  <- max(c(dot_totalRNA_l2$pct.exp, dot_rbRNA_l2$pct.exp))

P_l2_totalRNA <- plot_from_dot(dot_totalRNA_l2, max_pct_l2)
P_l2_rbRNA <- plot_from_dot(dot_rbRNA_l2, max_pct_l2)

ggsave("label2_totalRNA.pdf", P_l2_totalRNA, width = 10, height = 15)
ggsave("label2_rbRNA.pdf", P_l2_rbRNA, width = 10, height = 15)

# Calculate TE matrix (label2 and label3)
# Extract expression matrix
library(anndata)
library(reticulate)
use_condaenv("scanpy_env", required = TRUE)
py_config()
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")
totalRNA_matrix <- t(sc_ad$layers[["totalRNA_raw"]])
rbRNA_matrix <- t(sc_ad$layers[["rbRNA"]])
meta <- sc_ad$obs

# Remove mixed cells
meta$label2 <- states$states_nn_alg1_label2[rownames(meta)]
meta$label3 <- states$states_nn_alg1_label3[rownames(meta)]
valid_cells <- which(!(meta$label2 %in% c("Mix") | meta$label3 %in% c("Mix", "TEGLU Mix")))
meta <- meta[valid_cells, ]
totalRNA_matrix <- totalRNA_matrix[, valid_cells]
rbRNA_matrix <- rbRNA_matrix[, valid_cells]

# Define TE calculation function
calculate_TE_matrix <- function(celltype_label, marker_list, celltype_order) {
  meta$CellType <- meta[[celltype_label]]
  celltypes <- intersect(celltype_order, unique(meta$CellType))
  genes <- intersect(marker_list, rownames(totalRNA_matrix))
  te_mat <- matrix(NA, nrow = length(genes), ncol = length(celltypes),
                   dimnames = list(genes, celltypes))
  for (ct in celltypes) {
    ct_cells <- rownames(meta)[meta$CellType == ct]
    ct_cells <- intersect(ct_cells, colnames(totalRNA_matrix))
    if (length(ct_cells) == 0) next
    totalRNA_sub <- totalRNA_matrix[genes, ct_cells, drop = FALSE]
    rbRNA_sub <- rbRNA_matrix[genes, ct_cells, drop = FALSE]
    totalRNA_sum <- rowMeans(totalRNA_sub)
    rbRNA_sum <- rowMeans(rbRNA_sub)
    te <- rbRNA_sum / totalRNA_sum
    te[is.infinite(te) | is.na(te)] <- NA
    te_mat[genes, ct] <- te
  }
  return(te_mat)
}

te_label2 <- calculate_TE_matrix("label2", markers_l2, order_l2)
te_label3 <- calculate_TE_matrix("label3", markers_l3, order_l3)

# TE dot plot functions
make_te_dotdata <- function(te_mat, ref_dotdata) {
  df <- as.data.frame(te_mat)
  df$features.plot <- rownames(df)
  df_long <- pivot_longer(df, -features.plot, names_to = "id", values_to = "avg.exp")
  
  # Maintain consistent order with totalRNA plot
  gene_order <- unique(ref_dotdata$features.plot)
  celltype_order <- unique(ref_dotdata$id)
  
  df_long <- df_long %>%
    mutate(
      features.plot = factor(features.plot, levels = gene_order), 
      id = factor(id, levels = celltype_order)
    ) %>%
    group_by(features.plot) %>%
    mutate(avg.exp.scaled = scale(avg.exp)[,1]) %>%
    ungroup()
  
  df_long
}

plot_te_dot <- function(df_te, totalRNA_dot, max_pct, scale_color = TRUE, title = "") {
  df_te <- left_join(df_te, totalRNA_dot[, c("features.plot", "id", "pct.exp")], 
                     by = c("features.plot", "id"))
  ggplot(df_te, aes(x = id, y = features.plot,
                    color = if (scale_color) avg.exp.scaled else avg.exp,
                    size = pct.exp)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_gradientn(
      colours = if (scale_color) c("blue", "white", "red") else c("white", "#FDBE85", "#D94701"),
      limits = if (scale_color) c(-2, 2) else c(0, NA),
      oob = scales::squish,
      name = if (scale_color) "TE\n(z‑score)" else "TE") +
    scale_size(range = c(0,6), limits = c(0, max_pct), name = "% cells") +
    
    guides(size = guide_legend(order = 2), color = guide_colorbar(order = 1))
}

# label3 TE dot plot
te_dotdata_l3 <- make_te_dotdata(te_label3, dot_totalRNA_l3)
P_l3_te_scaled <- plot_te_dot(te_dotdata_l3, dot_totalRNA_l3, max_pct_l3, scale_color = TRUE)
ggsave("label3_TE_scaled.pdf", P_l3_te_scaled, width = 10, height = 15)

# label2 TE dot plot
te_dotdata_l2 <- make_te_dotdata(te_label2, dot_totalRNA_l2)
P_l2_te_scaled <- plot_te_dot(te_dotdata_l2, dot_totalRNA_l2, max_pct_l2, scale_color = TRUE)
ggsave("label2_TE_scaled.pdf", P_l2_te_scaled, width = 10, height = 15)


