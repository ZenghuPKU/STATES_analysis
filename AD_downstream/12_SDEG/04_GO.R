rm(list = ls())

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/12_SDEG")

genes_df <- read.csv("ring1-3_SDEG_totalRNA.csv", stringsAsFactors = FALSE)
genes_to_use <- unique(genes_df$gene) 

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)

perform_go_analysis <- function(genes, ont) {
  enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
}

process_go_data <- function(go_data, category) {
  if(is.null(go_data) || nrow(go_data) == 0) return(NULL)
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

create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  category_order <- c("Cellular Component", "Molecular Function", "Biological Process")
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean) %>% na.omit()
  
  if (nrow(go_combined) == 0) {
    warning("No enrichment terms for plot: ", title)
    return(NULL)
  }
  
  go_combined$Category <- factor(go_combined$Category, levels = category_order)
  
  go_sorted <- data.frame()
  for(cat in category_order) {
    cat_data <- go_combined %>%
      filter(Category == cat)
    if (!("pval_log" %in% colnames(cat_data))) next
    if(nrow(cat_data)==0) next
    cat_data <- cat_data %>% arrange(pval_log)
    go_sorted <- bind_rows(go_sorted, cat_data)
  }
  go_combined <- go_sorted
  if ("Description" %in% colnames(go_combined)) {
    go_combined$Description <- factor(go_combined$Description, levels = go_combined$Description)
  }
  
  category_colors_bar <- c(
    "Biological Process" = "#9FD4EC",
    "Molecular Function" = "#FFDFAD", 
    "Cellular Component" = "#D0BCDF"
  )
  
  category_colors_points <- c(
    "Biological Process" = "#0F89CA",
    "Molecular Function" = "#FCA828",
    "Cellular Component" = "#74509C"
  )
  
  if (!("pval_log" %in% colnames(go_combined)) || !("GeneRatio" %in% colnames(go_combined)) ||
      all(is.na(go_combined$pval_log)) || all(is.na(go_combined$GeneRatio)) || 
      max(go_combined$GeneRatio, na.rm = TRUE) == 0) {
    warning("Skip plot: insufficient data.")
    return(NULL)
  }
  
  scale_factor <- max(go_combined$pval_log, na.rm = TRUE) / max(go_combined$GeneRatio, na.rm = TRUE)
  
  p <- ggplot(go_combined, aes(y = Description)) +
    geom_bar(aes(x = pval_log, fill = Category), stat = "identity", alpha = 0.6) +
    scale_fill_manual(values = category_colors_bar) +
    geom_path(aes(x = GeneRatio * scale_factor, group = Category), color = "black", size = 0.75) +
    geom_point(aes(x = GeneRatio * scale_factor, fill = Category),
               shape = 21,
               color = "black", 
               size = 3) +
    scale_fill_manual(values = category_colors_points) +
    scale_x_continuous(
      name = "-log10(adjusted p-value)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Gene Ratio")
    ) +
    labs(y = "", fill = "Category", title = title) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 12),
      axis.ticks = element_line(colour = "black", linewidth = 1),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
  return(p)
}

cat("GO enrichment (no background universe)\n")
go_bp <- perform_go_analysis(genes_to_use, "BP")
go_mf <- perform_go_analysis(genes_to_use, "MF")
go_cc <- perform_go_analysis(genes_to_use, "CC")

if (!is.null(go_bp) && nrow(go_bp) > 0) write.csv(as.data.frame(go_bp), "ring1-3_SDEG_totalRNA_GO_BP.csv", row.names = FALSE)
if (!is.null(go_mf) && nrow(go_mf) > 0) write.csv(as.data.frame(go_mf), "ring1-3_SDEG_totalRNA_GO_MF.csv", row.names = FALSE)
if (!is.null(go_cc) && nrow(go_cc) > 0) write.csv(as.data.frame(go_cc), "ring1-3_SDEG_totalRNA_GO_CC.csv", row.names = FALSE)

plot_title <- "ring1-3 SDEG totalRNA GO Enrichment(14mAD)"
p <- create_combined_plot(go_bp, go_mf, go_cc, plot_title)
if (!is.null(p)) {
  print(p)
  ggsave("ring1-3_SDEG_totalRNA_GO_enrichment(14mAD).pdf", plot = p, width = 10, height = 4)
}
