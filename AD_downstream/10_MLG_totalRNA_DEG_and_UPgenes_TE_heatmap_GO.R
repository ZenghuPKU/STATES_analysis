# 10_MLG_totalRNA_DEG_and_UPgenes_TE_heatmap_GO
rm(list = ls())
library(Seurat)
library(Matrix)
library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("~/AD_downstream/")
out_dir <- getwd()
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

k_clust  <- 3
prop_min <- 0.05    
min_cells_grp <- 5  
min_genes_deg <- 10 
lfc_cut  <- 0.4
fdr_cut  <- 0.05
eps      <- 1e-6
max_logFC     <- 4
max_negLogFDR <- 50

states <- readRDS("states_with_plaque_info.rds")
meta <- states@meta.data

mlg_cts   <- c("MLG1", "MLG2", "MLG3")
cells_mlg <- rownames(meta)[meta$states_nn_alg1_label3 %in% mlg_cts]

length(cells_mlg)
table(meta$states_nn_alg1_label3[cells_mlg])

totalRNA_counts_all <- GetAssayData(states, assay = "totalRNA", slot = "counts")
rbRNA_counts_all    <- GetAssayData(states, assay = "rbRNA",    slot = "counts")
totalRNA_norm_all   <- GetAssayData(states, assay = "totalRNA", slot = "data")

totalRNA_counts_mlg <- totalRNA_counts_all[, cells_mlg, drop = FALSE]
rbRNA_counts_mlg    <- rbRNA_counts_all[,    cells_mlg, drop = FALSE]
totalRNA_norm_mlg   <- totalRNA_norm_all[,   cells_mlg, drop = FALSE]

common_genes <- Reduce(intersect, list(
  rownames(totalRNA_counts_mlg),
  rownames(rbRNA_counts_mlg),
  rownames(totalRNA_norm_mlg)
))
totalRNA_counts_mlg <- totalRNA_counts_mlg[common_genes, , drop = FALSE]
rbRNA_counts_mlg    <- rbRNA_counts_mlg[common_genes,    , drop = FALSE]
totalRNA_norm_mlg   <- totalRNA_norm_mlg[common_genes,   , drop = FALSE]

ct_vec <- meta[cells_mlg, "states_nn_alg1_label3"]
ct_vec <- droplevels(factor(ct_vec, levels = mlg_cts))

n_cells <- ncol(totalRNA_counts_mlg)
expr_ge1       <- totalRNA_counts_mlg >= 1
cell_count_ge1 <- Matrix::rowSums(expr_ge1)
prop_ge1       <- cell_count_ge1 / n_cells

genes_keep <- names(prop_ge1)[prop_ge1 >= prop_min]
message("Genes pass raw>=1 in ≥", prop_min*100, "% cells: ", length(genes_keep))

totalRNA_counts_mlg <- totalRNA_counts_mlg[genes_keep, , drop = FALSE]
rbRNA_counts_mlg    <- rbRNA_counts_mlg[genes_keep,    , drop = FALSE]
totalRNA_norm_mlg   <- totalRNA_norm_mlg[genes_keep,   , drop = FALSE]

cells_A <- cells_mlg[ct_vec == "MLG3"]
cells_B <- cells_mlg[ct_vec %in% c("MLG1","MLG2")]

cond_vec <- rep(NA_character_, length(cells_mlg))
names(cond_vec) <- cells_mlg
cond_vec[cells_A] <- "A"
cond_vec[cells_B] <- "B"
cond_vec <- cond_vec[!is.na(cond_vec)]
cond_vec <- factor(cond_vec, levels = c("B","A"))  

cells_use_deg <- names(cond_vec)

message("DEG cell table:")
print(table(cond_vec))

total_deg <- totalRNA_norm_mlg[, cells_use_deg, drop = FALSE]

logfc_rna_vec <- numeric(length(genes_keep))
pval_rna_vec  <- numeric(length(genes_keep))
names(logfc_rna_vec) <- genes_keep
names(pval_rna_vec)  <- genes_keep

for (gn in genes_keep) {
  expr_g <- total_deg[gn, ]
  expr_A <- as.numeric(expr_g[cond_vec == "A"])
  expr_B <- as.numeric(expr_g[cond_vec == "B"])
  
  expr_A <- expr_A[!is.na(expr_A)]
  expr_B <- expr_B[!is.na(expr_B)]
  
  if (length(expr_A) >= min_cells_grp && length(expr_B) >= min_cells_grp) {
    logfc_rna_vec[gn] <- log2(mean(expr_A) + eps) - log2(mean(expr_B) + eps)
    pval_rna_vec[gn]  <- tryCatch(
      wilcox.test(expr_A, expr_B)$p.value,
      error = function(e) NA_real_
    )
  } else {
    logfc_rna_vec[gn] <- NA_real_
    pval_rna_vec[gn]  <- NA_real_
  }
}

deg_rna_df <- data.frame(
  Gene   = genes_keep,
  logFC  = logfc_rna_vec,
  pvalue = pval_rna_vec
) %>%
  filter(!is.na(logFC) & !is.na(pvalue)) %>%
  mutate(
    FDR       = p.adjust(pvalue, method = "BH"),
    negLogFDR = -log10(FDR + 1e-300),
    Sig = case_when(
      logFC >  lfc_cut & FDR < fdr_cut ~ "Up",
      logFC < -lfc_cut & FDR < fdr_cut ~ "Down",
      TRUE                             ~ "NS"
    )
  )

out_deg_csv <- file.path(out_dir, "totalRNA_DEG_MLG3_vs_MLG1_2.csv")
write.csv(deg_rna_df, out_deg_csv, row.names = FALSE)

if (nrow(deg_rna_df) > 0) {
  df <- deg_rna_df %>%
    mutate(
      logFC_capped   = pmax(pmin(logFC,  max_logFC), -max_logFC),
      negLogFDR_cap  = pmin(negLogFDR, max_negLogFDR)
    )
  
  n_up   <- sum(df$Sig == "Up")
  n_down <- sum(df$Sig == "Down")
  
  top_genes <- df %>%
    filter(Sig != "NS") %>%
    arrange(FDR) %>%
    group_by(Sig) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  p_vol <- ggplot(df, aes(x = logFC_capped, y = negLogFDR_cap)) +
    geom_point(aes(color = Sig), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Up" = "#B2182B",
                                  "Down" = "#2166AC",
                                  "NS" = "gray")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(fdr_cut),
               linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(-max_logFC, max_logFC)) +
    scale_y_continuous(limits = c(0, max_negLogFDR)) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene),
      size = 3, max.overlaps = 100
    ) +
    annotate("text", x =  max_logFC * 0.9, y = max_negLogFDR,
             label = paste0("Up: ", n_up),
             color = "#B2182B", size = 4, hjust = 1) +
    annotate("text", x = -max_logFC * 0.9, y = max_negLogFDR,
             label = paste0("Down: ", n_down),
             color = "#2166AC", size = 4, hjust = 0) +
    theme_classic() +
    labs(
      title = paste0("totalRNA Volcano: MLG3 vs (MLG1+MLG2)\n",
                     "(", nrow(df), " genes; raw≥1 in ≥",
                     prop_min*100, "% cells)"),
      x = "log2 Fold Change (MLG3 vs MLG1+MLG2)",
      y = "-log10(FDR)",
      color = "Significance"
    )
  
  out_vol <- file.path(out_dir, "Volcano_totalRNA_MLG3_vs_MLG1_2.pdf")
  ggsave(out_vol, p_vol, width = 7, height = 6)
  message("Volcano saved: ", out_vol)
}

genes_up <- deg_rna_df %>%
  filter(Sig == "Up") %>%
  pull(Gene) %>%
  unique()

message("Up genes (logFC>", lfc_cut, " & FDR<", fdr_cut, "): ", length(genes_up))

genes_use <- intersect(genes_up, rownames(totalRNA_norm_mlg))
length(genes_use)

ct_levels <- mlg_cts
n_ct      <- length(ct_levels)

totalRNA_norm_byCT <- matrix(
  NA_real_, nrow = length(genes_use), ncol = n_ct,
  dimnames = list(genes_use, ct_levels)
)
totalRNA_sum_byCT  <- totalRNA_norm_byCT
rbRNA_sum_byCT     <- totalRNA_norm_byCT

for (ct in ct_levels) {
  cells_ct <- cells_mlg[ct_vec == ct]
  if (length(cells_ct) == 0) next
  
  totalRNA_norm_byCT[, ct] <- Matrix::rowMeans(
    totalRNA_norm_mlg[genes_use, cells_ct, drop = FALSE]
  )
  totalRNA_sum_byCT[, ct]  <- Matrix::rowSums(
    totalRNA_counts_mlg[genes_use, cells_ct, drop = FALSE]
  )
  rbRNA_sum_byCT[, ct]     <- Matrix::rowSums(
    rbRNA_counts_mlg[genes_use, cells_ct, drop = FALSE]
  )
}

TE_byCT <- rbRNA_sum_byCT / totalRNA_sum_byCT
TE_byCT[!is.finite(TE_byCT)] <- NA_real_

valid_genes <- (rowSums(is.na(totalRNA_norm_byCT)) == 0) &
  (rowSums(is.na(TE_byCT))            == 0)

totalRNA_norm_byCT <- totalRNA_norm_byCT[valid_genes, , drop = FALSE]
TE_byCT            <- TE_byCT[valid_genes,            , drop = FALSE]

message("Up genes kept after removing NA rows: ", nrow(totalRNA_norm_byCT))

scale_rows <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- matrixStats::rowSds(mat, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- 1
  (mat - m) / s
}

totalRNA_z <- scale_rows(totalRNA_norm_byCT)
TE_z       <- scale_rows(TE_byCT)

clust_mat <- cbind(totalRNA_z, TE_z)
colnames(clust_mat) <- c(
  paste0("RNA_", colnames(totalRNA_z)),
  paste0("TE_",  colnames(TE_z))
)

gene_dist   <- dist(clust_mat)
gene_hclust <- hclust(gene_dist, method = "ward.D2")
gene_clusters <- cutree(gene_hclust, k = k_clust)
row_split    <- factor(gene_clusters, levels = 1:k_clust)

gene_order       <- rownames(clust_mat)[gene_hclust$order]
totalRNA_z_ord   <- totalRNA_z[gene_order, , drop = FALSE]
TE_z_ord         <- TE_z[gene_order,       , drop = FALSE]
row_split_ord    <- row_split[gene_order]

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
  matrix           = totalRNA_z_ord,
  name             = "totalRNA_Z",
  col              = col_fun_total,
  show_row_names   = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_rot = 45,
  column_title     = paste0("Up genes (MLG3 vs MLG1+2)\nMLG totalRNA (Z, k=", k_clust, ")"),
  column_title_side = "top",
  row_split        = row_split_ord,
  row_gap          = unit(1.5, "mm"),
  column_gap       = unit(4, "mm"),
  border           = TRUE
)

ht_te <- Heatmap(
  matrix           = TE_z_ord,
  name             = "TE_Z",
  col              = col_fun_te,
  show_row_names   = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_rot = 45,
  column_title     = paste0("Up genes (MLG3 vs MLG1+2)\nMLG TE (Z, k=", k_clust, ")"),
  column_title_side = "top",
  row_split        = row_split_ord,
  row_gap          = unit(1.5, "mm"),
  column_gap       = unit(4, "mm"),
  border           = TRUE
)

ht_list <- ht_total + ht_te

out_ht_pdf <- file.path(out_dir,
                        paste0("MLG_UpGenes_totalRNA_TE_ComplexHeatmap_k", k_clust, ".pdf"))
pdf(out_ht_pdf, width = 6, height = 10)
draw(ht_list, heatmap_legend_side = "right")
dev.off()
message("Heatmap saved: ", out_ht_pdf)

cluster_df <- data.frame(
  gene    = gene_order,
  cluster = as.integer(row_split_ord)
)

totalRNA_z_export <- as.data.frame(totalRNA_z_ord)
colnames(totalRNA_z_export) <- paste0("totalRNA_Z_", colnames(totalRNA_z_export))

TE_z_export <- as.data.frame(TE_z_ord)
colnames(TE_z_export) <- paste0("TE_Z_", colnames(TE_z_export))

heatmap_matrix_df <- data.frame(
  Gene    = gene_order,
  Cluster = as.integer(row_split_ord),
  totalRNA_z_export,
  TE_z_export,
  check.names = FALSE
)

out_heatmap_matrix_csv <- file.path(
  out_dir,
  paste0("MLG_UpGenes_totalRNA_TE_heatmap_matrix_k", k_clust, ".csv")
)
write.csv(heatmap_matrix_df, out_heatmap_matrix_csv, row.names = FALSE)
message("Heatmap matrix table saved: ", out_heatmap_matrix_csv)

perform_go_analysis <- function(genes, ont) {
  enrichGO(
    gene          = genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
}

process_go_data <- function(go_data, category) {
  if (is.null(go_data)) return(NULL)
  go_df_raw <- as.data.frame(go_data)
  if (nrow(go_df_raw) == 0) return(NULL)
  
  go_df <- go_df_raw %>%
    tidyr::separate(GeneRatio,
                    into = c("GeneInTerm", "GeneInBackground"),
                    sep  = "/") %>%
    mutate(
      GeneRatio = as.numeric(GeneInTerm) / as.numeric(GeneInBackground),
      pval_log  = -log10(p.adjust),
      Category  = category
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5)
  go_df
}

create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean) %>% na.omit()
  if (nrow(go_combined) == 0) return(NULL)
  
  category_order <- c("Biological Process", "Molecular Function", "Cellular Component")
  go_combined$Category <- factor(go_combined$Category, levels = category_order)
  
  go_sorted <- bind_rows(lapply(category_order, function(cat) {
    go_combined %>%
      filter(Category == cat) %>%
      arrange(desc(pval_log))
  }))
  
  go_sorted$Description <- factor(go_sorted$Description, levels = rev(go_sorted$Description))
  scale_factor <- max(go_sorted$pval_log) / max(go_sorted$GeneRatio)
  
  category_colors_bar    <- c("Biological Process" = "#9FD4EC",
                              "Molecular Function" = "#FFDFAD",
                              "Cellular Component" = "#D0BCDF")
  category_colors_points <- c("Biological Process" = "#0F89CA",
                              "Molecular Function" = "#FCA828",
                              "Cellular Component" = "#74509C")
  
  p <- ggplot(go_sorted, aes(y = Description)) +
    geom_bar(aes(x = pval_log, fill = Category),
             stat = "identity", alpha = 0.6) +
    scale_fill_manual(values = category_colors_bar) +
    geom_path(aes(x = GeneRatio * scale_factor, group = Category),
              color = "black", size = 0.75) +
    geom_point(aes(x = GeneRatio * scale_factor, color = Category),
               size = 3) +
    scale_color_manual(values = category_colors_points) +
    scale_x_continuous(
      name = "-log10(adjusted p-value)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Gene Ratio")
    ) +
    labs(y = "", title = title) +
    theme_minimal() +
    theme(
      panel.border    = element_rect(colour = "black", fill = NA),
      axis.text       = element_text(color = "black", size = 12),
      axis.title      = element_text(size = 12),
      legend.position = "right",
      legend.text     = element_text(size = 12),
      legend.title    = element_text(size = 12)
    )
  p
}

go_dir <- file.path(out_dir, paste0("GO_k", k_clust))
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

clusters_vec <- cluster_df$cluster
names(clusters_vec) <- cluster_df$gene

for (cl in sort(unique(clusters_vec))) {
  
  gene_list <- names(clusters_vec)[clusters_vec == cl]
  
  if (length(gene_list) < 5) {
    cat("Cluster", cl, "genes <5, skip GO\n")
    next
  }
  
  cat(">>> Cluster", cl, " GO\n")
  
  go_bp <- perform_go_analysis(gene_list, "BP")
  go_mf <- perform_go_analysis(gene_list, "MF")
  go_cc <- perform_go_analysis(gene_list, "CC")
  
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  go_top <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean)
  
  p_go <- create_combined_plot(
    go_bp, go_mf, go_cc,
    title = paste0("GO - Cluster ", cl)
  )
  
  if (!is.null(p_go)) {
    ggsave(
      filename = file.path(go_dir,
                           paste0("GO_k", k_clust, "_Cluster", cl, ".pdf")),
      plot     = p_go,
      width    = 12,
      height   = 6
    )
  }
}

genes_cluster2 <- cluster_df$gene[cluster_df$cluster == 2]
genes_cluster2 <- gene_order[gene_order %in% genes_cluster2]

type_vec_mlg <- as.character(meta[cells_mlg, "type"])

type_levels_pref <- c("4mWT","8mWT","14mWT","4mAD","8mAD","14mAD")
type_levels_use  <- type_levels_pref[type_levels_pref %in% unique(type_vec_mlg)]

totalRNA_norm_byType <- matrix(
  NA_real_,
  nrow = length(genes_cluster2),
  ncol = length(type_levels_use),
  dimnames = list(genes_cluster2, type_levels_use)
)
totalRNA_sum_byType <- totalRNA_norm_byType
rbRNA_sum_byType    <- totalRNA_norm_byType

for (tp in type_levels_use) {
  cells_tp <- cells_mlg[type_vec_mlg == tp]
  if (length(cells_tp) == 0) next
  
  totalRNA_norm_byType[, tp] <- Matrix::rowMeans(
    totalRNA_norm_mlg[genes_cluster2, cells_tp, drop = FALSE]
  )
  totalRNA_sum_byType[, tp]  <- Matrix::rowSums(
    totalRNA_counts_mlg[genes_cluster2, cells_tp, drop = FALSE]
  )
  rbRNA_sum_byType[, tp]     <- Matrix::rowSums(
    rbRNA_counts_mlg[genes_cluster2, cells_tp, drop = FALSE]
  )
}

TE_byType <- rbRNA_sum_byType / totalRNA_sum_byType
TE_byType[!is.finite(TE_byType)] <- NA_real_

valid_genes_c2 <- (rowSums(is.na(totalRNA_norm_byType)) == 0) &
  (rowSums(is.na(TE_byType))            == 0)

totalRNA_norm_byType_c2 <- totalRNA_norm_byType[valid_genes_c2, , drop = FALSE]
TE_byType_c2            <- TE_byType[valid_genes_c2,            , drop = FALSE]

totalRNA_z_type_c2 <- scale_rows(totalRNA_norm_byType_c2)
TE_z_type_c2       <- scale_rows(TE_byType_c2)

ht_total_c2_type <- Heatmap(
  matrix           = totalRNA_z_type_c2,
  name             = "totalRNA_Z",
  col              = col_fun_total,
  show_row_names   = TRUE,                
  row_names_side   = "left",                 
  row_names_gp     = grid::gpar(fontsize = 6), 
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_rot = 45,
  column_title     = "Cluster2 genes\n(MLG1/2/3, totalRNA by type, Z)",
  column_title_side = "top",
  row_gap          = unit(1.5, "mm"),
  column_gap       = unit(4, "mm"),
  border           = TRUE
)

ht_te_c2_type <- Heatmap(
  matrix           = TE_z_type_c2,
  name             = "TE_Z",
  col              = col_fun_te,
  show_row_names   = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_rot = 45,
  column_title     = "Cluster2 genes\n(MLG1/2/3, TE by type, Z)",
  column_title_side = "top",
  row_gap          = unit(1.5, "mm"),
  column_gap       = unit(4, "mm"),
  border           = TRUE
)

ht_list_c2_type <- ht_total_c2_type + ht_te_c2_type

out_ht_c2_type_pdf <- file.path(
  out_dir,
  paste0("MLG_Cluster2_genes_totalRNA_TE_byType_ComplexHeatmap.pdf")
)
pdf(out_ht_c2_type_pdf, width = 7, height = 10)
draw(ht_list_c2_type, heatmap_legend_side = "right")
dev.off()