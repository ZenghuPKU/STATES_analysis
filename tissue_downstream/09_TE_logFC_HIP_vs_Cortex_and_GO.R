# 09_TE_logFC_HIP_vs_Cortex_and_GO
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(anndata)
rm(list = ls())
setwd("~/tissue_downstream/")

# Load the annotated object (states_celltypes_identification.RData), Read h5ad data
load("states_celltypes_identification.RData")
sc_ad <- read_h5ad("mousebrain_harmony.h5ad")

# Prepare expression matrix
totalRNA_matrix <- t(sc_ad$layers[["totalRNA"]])
rbRNA_matrix <- t(sc_ad$layers[["rbRNA"]])
selected_genes <- rownames(sc_ad$var)[sc_ad$var$highly_variable == TRUE]
totalRNA_matrix <- totalRNA_matrix[selected_genes, , drop = FALSE]
rbRNA_matrix <- rbRNA_matrix[selected_genes, , drop = FALSE]

# Calculate translation efficiency (TE)
te_matrix <- rbRNA_matrix / totalRNA_matrix
te_matrix[is.na(te_matrix) | is.infinite(te_matrix)] <- NA
te_matrix[totalRNA_matrix < 2] <- NA

# Set cell labels & groups
cell_labels <- states$states_nn_alg1_label3
table(states$states_nn_alg1_label3)
cell_groups <- list(
  HIP = c("TEGLU CA1", "TEGLU CA2", "TEGLU CA3", "DGGRC"),
  Cortex = c("TEGLU L2/3", "TEGLU L4", "TEGLU L5", "TEGLU L5/6","TEGLU L6","TEGLU L6b","TEGLU Mix")
)

# Initialize variables
group_names <- names(cell_groups)
summary_list <- list()

# Pairwise group comparisons
for (i in 1:(length(group_names) - 1)) {
  for (j in (i + 1):length(group_names)) {
    name1 <- group_names[i]
    name2 <- group_names[j]
    
    group1 <- cell_groups[[name1]]
    group2 <- cell_groups[[name2]]
    
    comp_name <- paste0(name1, "_vs_", name2)
    safe_name <- gsub("[^A-Za-z0-9_]+", "_", comp_name)
    
    cat("\n========= Comparing:", name1, "vs", name2, "=========\n")
    
    target_cells <- colnames(te_matrix)[cell_labels %in% group1]
    other_cells <- colnames(te_matrix)[cell_labels %in% group2]
    
    # Filter genes: ≥5% of group cells with totalRNA ≥ 1
    states_sub <- totalRNA_matrix[selected_genes, target_cells]
    min_cells_expr <- ceiling(length(target_cells) * 0.05)
    genes_pass_filter <- rownames(states_sub)[rowSums(states_sub >= 1, na.rm = TRUE) >= min_cells_expr]
    
    cat("Number of genes passing filter", length(genes_pass_filter), "\n")
    if (length(genes_pass_filter) < 10) next
    
    # Calculate logFC and p 
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
      )
    
    write.csv(result_df, paste0("TE_logFC_", safe_name, ".csv"), row.names = FALSE)
    n_genes_used <- nrow(result_df)
    
    # Cap values
    result_df <- result_df %>%
      mutate(
        logFC_capped = pmax(pmin(logFC, 4), -4),
        negLogP_capped = pmin(negLogP, 50)
      )
    
    # Top genes
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
    
    # Volcano plot
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
        title = paste0("TE Volcano Plot: ", comp_name,
                       "\n(", n_genes_used, " genes used, totalRNA≥1 in ≥5%)"),
        x = "log2 Fold Change (TE)",
        y = "-log10(P-value)",
        color = "Significance"
      )
    
    ggsave(paste0("Volcano_TE_", safe_name, ".pdf"), p, width = 7, height = 6)
    
    # Collect up- and down-regulated genes for GO analysis
    up_genes <- result_df$Gene[result_df$Sig == "Up"]
    down_genes <- result_df$Gene[result_df$Sig == "Down"]
    
    if (length(up_genes) > 0) {
      summary_list[[paste0(comp_name, "_Up")]] <- data.frame(
        CellType = comp_name,
        Gene = up_genes,
        Direction = "Up"
      )
    }
    if (length(down_genes) > 0) {
      summary_list[[paste0(comp_name, "_Down")]] <- data.frame(
        CellType = comp_name,
        Gene = down_genes,
        Direction = "Down"
      )
    }
  }
}
# Save consolidated up/down-regulated genes
summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, "TElogFC_MultiGroupComparisons_UpDownGenes_forGO.csv", row.names = FALSE)

# Set up GO output directory
go_dir <- "./GO/"
dir.create(go_dir, showWarnings = FALSE)

# Define GO analysis function
perform_go_analysis <- function(genes, ont) {
  enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    maxGSSize = 200
  )
}

# GO results processing function (limit to 10 terms)
process_go_data <- function(go_data, category) {
  if (nrow(go_data) == 0) return(NULL)
  as.data.frame(go_data) %>%
    separate(GeneRatio, into = c("GeneInTerm", "GeneInBackground"), sep = "/") %>%
    mutate(
      GeneRatio = as.numeric(GeneInTerm) / as.numeric(GeneInBackground),
      pval_log = -log10(p.adjust),
      Category = category
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5)
}

# Plotting function
create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean) %>% na.omit()
  if (nrow(go_combined) == 0) return(NULL)
  
  category_order <- c("Biological Process", "Molecular Function", "Cellular Component")
  go_combined$Category <- factor(go_combined$Category, levels = category_order)
  
  go_sorted <- bind_rows(lapply(category_order, function(cat) {
    go_combined %>% filter(Category == cat) %>% arrange(desc(pval_log))
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
# Batch process GO analysis for up/down-regulated groups
group_info <- unique(summary_df[, c("CellType", "Direction")])

for (i in seq_len(nrow(group_info))) {
  celltype <- group_info$CellType[i]
  direction <- group_info$Direction[i]
  cat("🔍 Processing GO for:", celltype, "-", direction, "\n")
  
  gene_list <- summary_df %>%
    filter(CellType == celltype, Direction == direction) %>%
    pull(Gene) %>% unique()
  
  if (length(gene_list) < 5) next
  
  go_bp <- perform_go_analysis(gene_list, "BP")
  go_mf <- perform_go_analysis(gene_list, "MF")
  go_cc <- perform_go_analysis(gene_list, "CC")
  
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean) %>% na.omit()
  
  plot_title <- paste0("GO Enrichment: ", celltype, " - ", direction)
  p <- create_combined_plot(go_bp, go_mf, go_cc, plot_title)
  
  if (!is.null(p)) {
    filename <- paste0("GO_", gsub("[^A-Za-z0-9_]+", "_", celltype), "_", direction, ".pdf")
    ggsave(file.path(go_dir, filename), plot = p, width = 12, height = 6)
  }
  
  if (nrow(go_combined) > 0) {
    filename_csv <- paste0("GO_", gsub("[^A-Za-z0-9_]+", "_", celltype), "_", direction, ".csv")
    write.csv(go_combined, file = filename_csv, row.names = FALSE)
  }
}

