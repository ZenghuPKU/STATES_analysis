rm(list=ls())

bottom_10_percent <- read.csv("late_downregulated_genes.csv")
top_10_percent <- read.csv("early_upregulated_genes.csv")

background_genes <- read.csv("back.csv")
background_genes <- unique(background_genes$GeneSymbol)

library(clusterProfiler)
library(org.Hs.eg.db)

bottom_genes <- bottom_10_percent[,1]
top_genes <- top_10_percent[,1]
perform_go_analysis <- function(genes, ont) {
  enrichGO(gene = genes,
           OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL",
           ont = ont,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           universe = background_genes)
}

bottom_go_bp <- perform_go_analysis(bottom_genes, "BP")
bottom_go_mf <- perform_go_analysis(bottom_genes, "MF")
bottom_go_cc <- perform_go_analysis(bottom_genes, "CC")

top_go_bp <- perform_go_analysis(top_genes, "BP")
top_go_mf <- perform_go_analysis(top_genes, "MF")
top_go_cc <- perform_go_analysis(top_genes, "CC")


library(ggplot2)
library(dplyr) 
library(tidyr)

# Function to process GO results
process_go_data <- function(go_data, category) {
  if (nrow(go_data) == 0) {
    return(NULL)
  }
  
  go_select <- as.data.frame(go_data) %>%
    separate(GeneRatio, into = c("GeneInTerm", "GeneInBackground"), sep = "/") %>%
    mutate(
      GeneRatio = as.numeric(GeneInTerm) / as.numeric(GeneInBackground),
      pval_log = -log10(p.adjust),
      Category = category
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) 
  return(go_select)
}

# Function to create combined plot
create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  # Process data
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function") 
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  
  # Combine data
  category_order <- c("Cellular Component", "Molecular Function", "Biological Process")
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean)%>% na.omit()
  go_combined$Category <- factor(go_combined$Category, levels = category_order)
  
  
  # Sort within categories
  go_sorted <- data.frame()
  for(cat in category_order) {
    cat_data <- go_combined %>%
      filter(Category == cat) %>%
      arrange(pval_log)
    go_sorted <- bind_rows(go_sorted, cat_data)
  }
  go_combined <- go_sorted
  go_combined$Description <- factor(go_combined$Description, levels = go_combined$Description)
  
  
  # Define colors
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
  
  # Calculate scale factor
  scale_factor <- max(go_combined$pval_log, na.rm = TRUE) / max(go_combined$GeneRatio, na.rm = TRUE)
  
  # Create plot
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

# Create plots for each group
p_bottom <- create_combined_plot(bottom_go_bp, bottom_go_mf, bottom_go_cc, "Bottom 10%")
p_bottom
p_top <- create_combined_plot(top_go_bp, top_go_mf, top_go_cc, "Top 10%")  
p_top
# Save plots
ggsave("late_downregulated_go.pdf", plot = p_bottom, width = 10, height = 4) 
ggsave("early_upregulated_go.pdf", plot = p_top, width = 10, height = 3) 
