# 07_Gene_clustering_based_on_TE_and_GO
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(Seurat)
library(anndata)
library(matrixStats)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
rm(list = ls())
setwd("~/tissue_downstream/")

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Load expression matrix
totalRNA_matrix <- t(sc_ad$layers['totalRNA'])
rbRNA_matrix <- t(sc_ad$layers['rbRNA'])
metadata <- sc_ad$obs

# Filter highly variable genes
selected_genes <- rownames(sc_ad$var)[sc_ad$var$highly_variable == TRUE]
totalRNA_matrix <- totalRNA_matrix[selected_genes, , drop = FALSE]
rbRNA_matrix <- rbRNA_matrix[selected_genes, , drop = FALSE]

# Create Seurat objects and add labels
totalRNA <- CreateSeuratObject(counts = totalRNA_matrix, meta.data = metadata)
rbRNA <- CreateSeuratObject(counts = rbRNA_matrix, meta.data = metadata)
cell_labels <- states$states_nn_alg1_label2
totalRNA@meta.data$CellType <- cell_labels
rbRNA@meta.data$CellType <- cell_labels
cell_types <- setdiff(unique(cell_labels), "Mix")

# Calculate translation efficiency (TE)
expr_matrix <- data.frame(gene = rownames(totalRNA), stringsAsFactors = FALSE)

for (ct in cell_types) {
  ct_cells <- rownames(totalRNA@meta.data)[totalRNA@meta.data$CellType == ct]
  totalRNA_counts <- GetAssayData(totalRNA[, ct_cells], slot = "counts")
  rbRNA_counts <- GetAssayData(rbRNA[, ct_cells], slot = "counts")
  totalRNA_counts <- as(totalRNA_counts, "dgCMatrix")
  rbRNA_counts <- as(rbRNA_counts, "dgCMatrix")
  
  sumtotalRNA <- rowSums(totalRNA_counts)
  sumrbRNA <- rowSums(rbRNA_counts)
  te_vec <- sumrbRNA / sumtotalRNA
  te_vec[is.na(te_vec) | is.infinite(te_vec)] <- 0
  
  expr_matrix[[paste0("te_", ct)]] <- log1p(te_vec)
}

# Z-score normalize TE
te_cols <- grep("^te_", colnames(expr_matrix), value = TRUE)
scaled_te <- t(scale(t(as.matrix(expr_matrix[, te_cols]))))
expr_matrix[, te_cols] <- scaled_te

# Preliminary clustering (TE)
mat_te <- as.matrix(expr_matrix[, te_cols])
rownames(mat_te) <- expr_matrix$gene
hc <- hclust(dist(mat_te), method = "ward.D2")
clusters <- cutree(hc, k = 10)
expr_matrix$cluster <- factor(clusters, levels = 1:10)

# Sort clusters by TE mean and reassign labels
cluster_means <- expr_matrix %>%
  group_by(cluster) %>%
  summarise(across(all_of(te_cols), mean, na.rm = TRUE))
cluster_order_val <- apply(cluster_means[, te_cols], 1, max, na.rm = TRUE)
cluster_order <- cluster_means$cluster[order(cluster_order_val, decreasing = TRUE)]
new_cluster_labels <- setNames(1:10, cluster_order)
expr_matrix_sorted <- expr_matrix %>%
  mutate(cluster = factor(new_cluster_labels[as.character(cluster)], levels = 1:10))

# Prepare heatmap data
manual_order <- c("TEPN", "INH", "DE/MEN", "CHO/PEP", "AC", "OLG", "OPC", "MLG", "CHOR/EPEN", "VAS")
new_te_cols <- paste0("te_", manual_order)
mat_te <- as.matrix(expr_matrix_sorted[, new_te_cols])
rownames(mat_te) <- expr_matrix_sorted$gene
row_split <- expr_matrix_sorted$cluster

# Construct RdBu color mapping (quantiles)
col_fun_te <- colorRamp2(
  quantile(mat_te, probs = c(0.2, 0.5, 0.8), na.rm = TRUE),
  c("#2166ac", "#ffffff", "#b2182b")
)

# Generate heatmap
ht_te <- Heatmap(
  matrix = mat_te,
  name = "TE",
  col = col_fun_te,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_title = "TE",
  column_title_side = "top",
  column_gap = unit(5, "mm"),
  row_split = row_split,
  border = TRUE
)

pdf("s10_te_only_heatmap_clusterSorted.pdf", width = 5, height = 8)
draw(ht_te)
dev.off()

# Export to CSV
write.csv(expr_matrix_sorted, "s10_te_only_heatmap_clusterSorted.csv")
print("Gene counts per TE cluster (sorted):")
print(table(expr_matrix_sorted$cluster))

# GO enrichment analysis for TE clustering results
dir.create("s10_TE_only_GO", showWarnings = FALSE)

# GO enrichment function
perform_go_analysis <- function(genes, ont) {
  enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    maxGSSize = 200,
    pvalueCutoff = 0.05
  )
}

# GO results processing function
process_go_data <- function(go_data, category) {
  if (nrow(go_data) == 0) return(NULL)
  go_df <- as.data.frame(go_data) %>%
    tidyr::separate(GeneRatio, into = c("GeneInTerm", "GeneInBackground"), sep = "/") %>%
    dplyr::mutate(
      GeneRatio = as.numeric(GeneInTerm) / as.numeric(GeneInBackground),
      pval_log = -log10(p.adjust),
      Category = category
    ) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 5)
  return(go_df)
}

# Plotting function
create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  go_combined <- dplyr::bind_rows(go_bp_clean, go_mf_clean, go_cc_clean) %>% na.omit()
  if (nrow(go_combined) == 0) return(NULL)
  
  category_order <- c("Biological Process", "Molecular Function", "Cellular Component")
  go_combined$Category <- factor(go_combined$Category, levels = category_order)
  
  go_sorted <- dplyr::bind_rows(lapply(category_order, function(cat) {
    go_combined %>% dplyr::filter(Category == cat) %>% dplyr::arrange(desc(pval_log))
  }))
  
  go_sorted$Description <- factor(go_sorted$Description, levels = rev(go_sorted$Description))
  scale_factor <- max(go_sorted$pval_log) / max(go_sorted$GeneRatio)
  
  category_colors_bar <- c("Biological Process" = "#9FD4EC", "Molecular Function" = "#FFDFAD", "Cellular Component" = "#D0BCDF")
  category_colors_points <- c("Biological Process" = "#0F89CA", "Molecular Function" = "#FCA828", "Cellular Component" = "#74509C")
  
  p <- ggplot(go_sorted, aes(y = Description)) +
    geom_bar(aes(x = pval_log, fill = Category), stat = "identity", alpha = 0.6) +
    scale_fill_manual(values = category_colors_bar) +
    geom_path(aes(x = GeneRatio * scale_factor, group = Category), color = "black", size = 0.75) +
    geom_point(aes(x = GeneRatio * scale_factor, fill = Category),
               shape = 21, color = "black", size = 3) +
    scale_fill_manual(values = category_colors_points) +
    scale_x_continuous(
      name = "-log10(adjusted p-value)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Gene Ratio")
    ) +
    labs(y = "", fill = "Category", title = title) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
  return(p)
}

# Perform GO analysis for each cluster
for (cl in sort(unique(expr_matrix_sorted$cluster))) {
  cat("Processing GO for TE cluster:", cl, "\n")
  gene_list <- expr_matrix_sorted$gene[expr_matrix_sorted$cluster == cl]
  if (length(gene_list) < 5) next
  
  go_bp <- perform_go_analysis(gene_list, "BP")
  go_mf <- perform_go_analysis(gene_list, "MF")
  go_cc <- perform_go_analysis(gene_list, "CC")
  
  plot_title <- paste0("GO Enrichment - TE Cluster ", cl)
  p <- create_combined_plot(go_bp, go_mf, go_cc, plot_title)
  
  if (!is.null(p)) {
    ggsave(paste0("s10_TE_only_GO/TE_GO_Cluster", cl, ".pdf"), plot = p, width = 12, height = 6)
  }
}
