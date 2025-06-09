# 08_TE_logFC_major_celltypes
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(anndata)
rm(list = ls())
setwd("~/tissue_downstream/")

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Prepare expression matrix
totalRNA_matrix <- t(sc_ad$layers[["totalRNA"]])
rbRNA_matrix <- t(sc_ad$layers[["rbRNA"]])
metadata <- sc_ad$obs
selected_genes <- rownames(sc_ad$var)[sc_ad$var$highly_variable == TRUE]
totalRNA_matrix <- totalRNA_matrix[selected_genes, , drop = FALSE]
rbRNA_matrix <- rbRNA_matrix[selected_genes, , drop = FALSE]

# Calculate translation efficiency (TE)
te_matrix <- rbRNA_matrix / totalRNA_matrix
te_matrix[is.na(te_matrix) | is.infinite(te_matrix)] <- NA
te_matrix[totalRNA_matrix < 2] <- NA

# Set up loop & summary list
cell_labels <- states$states_nn_alg1_label2
cell_types <- unique(cell_labels)
summary_list <- list()

# Iterate through each cell type
for (target_celltype in cell_types) {
  cat("\n========= Processing:", target_celltype, "=========\n")
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", target_celltype)
  
  target_cells <- colnames(te_matrix)[cell_labels == target_celltype]
  other_cells <- colnames(te_matrix)[cell_labels != target_celltype]
  
  # Filter genes: ≥5% of target cells with totalRNA ≥ 1
  states_sub <- totalRNA_matrix[selected_genes, target_cells]
  min_cells_expr <- ceiling(length(target_cells) * 0.05)
  genes_pass_filter <- rownames(states_sub)[rowSums(states_sub >= 1, na.rm = TRUE) >= min_cells_expr]
  
  cat("Number of genes passing filter (totalRNA ≥ 1 in ≥5% of cells):", length(genes_pass_filter), "\n")
  if (length(genes_pass_filter) < 10) next
  
  # Calculate logFC and p-values
  logfc_vec <- numeric(length(genes_pass_filter))
  pval_vec <- numeric(length(genes_pass_filter))
  names(logfc_vec) <- genes_pass_filter
  names(pval_vec) <- genes_pass_filter
  
  for (gene in genes_pass_filter) {
    te_target <- te_matrix[gene, target_cells]
    te_other <- te_matrix[gene, other_cells]
    
    te_target <- te_target[!is.na(te_target)]
    te_other <- te_other[!is.na(te_other)]
    
    if (length(te_target) >= 5 && length(te_other) >= 5) {
      logfc_vec[gene] <- log2(mean(te_target) + 1e-6) - log2(mean(te_other) + 1e-6)
      pval_vec[gene] <- wilcox.test(te_target, te_other)$p.value
    } else {
      logfc_vec[gene] <- NA
      pval_vec[gene] <- NA
    }
  }
  
  # Organize results
  result_df <- data.frame(
    Gene = genes_pass_filter,
    logFC = logfc_vec,
    pvalue = pval_vec
  ) %>%
    mutate(
      negLogP = -log10(pvalue),
      Sig = case_when(
        logFC > 0.5 & pvalue < 0.05 ~ "Up",
        logFC < -0.5 & pvalue < 0.05 ~ "Down",
        TRUE ~ "NS"
      )
    )%>%
    filter(!is.na(logFC) & !is.na(pvalue))
  
  write.csv(result_df, paste0("TE_logFC_", safe_name, ".csv"), row.names = FALSE)
  n_genes_used <- nrow(result_df)
  
  # Constrain value ranges
  result_df <- result_df %>%
    mutate(
      logFC_capped = pmax(pmin(logFC, 4), -4),
      negLogP_capped = pmin(negLogP, 50)
    )
  
  # Annotate top genes
  top_genes <- result_df %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 5) %>%
    bind_rows(
      result_df %>%
        filter(pvalue < 0.05) %>%
        arrange(logFC) %>%
        slice_head(n = 5)
    )
  
  # Count up- and down-regulated genes
  n_up <- sum(result_df$Sig == "Up", na.rm = TRUE)
  n_down <- sum(result_df$Sig == "Down", na.rm = TRUE)
  
  # Generate volcano plot
  p <- ggplot(result_df, aes(x = logFC_capped, y = negLogP_capped)) +
    geom_point(aes(color = Sig), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(-4, 4)) +
    scale_y_continuous(limits = c(0, 55)) +
    geom_text_repel(data = top_genes,
                    aes(x = logFC_capped, y = negLogP_capped, label = Gene),
                    size = 3, max.overlaps = 100) +
    annotate("text", x = 3.5, y = 55, label = paste0("Up: ", n_up), color = "#B2182B", size = 4, hjust = 1) +
    annotate("text", x = -3.5, y = 55, label = paste0("Down: ", n_down), color = "#2166AC", size = 4, hjust = 0) +
    theme_classic() +
    labs(
      title = paste0("TE Volcano Plot: ", target_celltype,
                     "\n(", n_genes_used, " genes used, totalRNA≥1 in ≥5%)"),
      x = "log2 Fold Change (TE)",
      y = "-log10(P-value)",
      color = "Significance"
    )
  
  ggsave(paste0("Volcano_TE_", safe_name, ".pdf"), p, width = 7, height = 6)
  
  # Record up- and down-regulated genes
  up_genes <- result_df$Gene[result_df$Sig == "Up"]
  down_genes <- result_df$Gene[result_df$Sig == "Down"]
  
  if (length(up_genes) > 0) {
    summary_list[[paste0(target_celltype, "_Up")]] <- data.frame(
      CellType = target_celltype,
      Gene = up_genes,
      Direction = "Up"
    )
  }
  if (length(down_genes) > 0) {
    summary_list[[paste0(target_celltype, "_Down")]] <- data.frame(
      CellType = target_celltype,
      Gene = down_genes,
      Direction = "Down"
    )
  }
}

# Consolidate all up- and down-regulated genes into a single table
summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, "TElogFC_AllCellTypes_UpDownGenes_forGO.csv", row.names = FALSE)
