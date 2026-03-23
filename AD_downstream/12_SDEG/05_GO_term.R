rm(list = ls())

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/12_SDEG")

genes_df <- read.csv("ring30_14mAD_SDEG_totalRNA.csv", stringsAsFactors = FALSE)
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


get_top5_terms <- function(go_obj) {
  if (is.null(go_obj) || nrow(go_obj) == 0) return(NULL)
  
  go_df <- as.data.frame(go_obj) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5)
  
  go_df$geneID <- gsub("/", ",", go_df$geneID)
  colnames(go_df)[colnames(go_df) == "geneID"] <- "Genes"
  
  return(go_df)
}


process_go_data <- function(go_data, category) {
  if(is.null(go_data) || nrow(go_data) == 0) return(NULL)
  
  as.data.frame(go_data) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    separate(GeneRatio, into = c("GeneInTerm", "GeneInBackground"), sep = "/") %>%
    mutate(
      GeneRatio = as.numeric(GeneInTerm) / as.numeric(GeneInBackground),
      pval_log = -log10(p.adjust),
      Category = category
    )
}

create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean)
  
  if (is.null(go_combined) || nrow(go_combined) == 0) {
    warning("No enrichment terms for plot.")
    return(NULL)
  }
  
  go_combined$Description <- factor(go_combined$Description,
                                    levels = go_combined$Description)
  
  scale_factor <- max(go_combined$pval_log) / max(go_combined$GeneRatio)
  
  p <- ggplot(go_combined, aes(y = Description)) +
    geom_bar(aes(x = pval_log, fill = Category),
             stat = "identity", alpha = 0.6) +
    geom_path(aes(x = GeneRatio * scale_factor, group = Category),
              color = "black", linewidth = 0.75) +
    geom_point(aes(x = GeneRatio * scale_factor, fill = Category),
               shape = 21,
               color = "black", size = 3) +
    scale_x_continuous(
      name = "-log10(adjusted p-value)",
      sec.axis = sec_axis(~ . / scale_factor, name = "Gene Ratio")
    ) +
    labs(y = "", fill = "Category", title = title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text = element_text(color = "black", size = 12)
    )
  
  return(p)
}


cat("GO enrichment (no background universe)\n")

go_bp <- perform_go_analysis(genes_to_use, "BP")
go_mf <- perform_go_analysis(genes_to_use, "MF")
go_cc <- perform_go_analysis(genes_to_use, "CC")


top5_bp <- get_top5_terms(go_bp)
top5_mf <- get_top5_terms(go_mf)
top5_cc <- get_top5_terms(go_cc)

if (!is.null(top5_bp))
  write.csv(top5_bp,
            "ring30_SDEG_totalRNA_GO_BP_top5.csv",
            row.names = FALSE)

if (!is.null(top5_mf))
  write.csv(top5_mf,
            "ring30_SDEG_totalRNA_GO_MF_top5.csv",
            row.names = FALSE)

if (!is.null(top5_cc))
  write.csv(top5_cc,
            "ring30_SDEG_totalRNA_GO_CC_top5.csv",
            row.names = FALSE)


plot_title <- "ring30 SDEG totalRNA GO Enrichment (14mAD)"

p <- create_combined_plot(go_bp, go_mf, go_cc, plot_title)

if (!is.null(p)) {
  print(p)
  ggsave("ring30_SDEG_totalRNA_GO_enrichment_14mAD.pdf",
         plot = p, width = 10, height = 4)
}