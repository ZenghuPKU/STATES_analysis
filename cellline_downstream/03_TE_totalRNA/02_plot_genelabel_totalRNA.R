library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("/storage/lingyuan2/STATES_analysis/cellline_downstream/03_TE_totalRNA/totalRNAresults")

data <- read.csv("wilcoxon_results_last_four_columns.csv")

genes_to_label_list <- c("DDIT3", "MANF", "HSPA5", "SEC61B", "IL27RA","PARP14","HERPUD1")

plot_volcano <- function(data, pval_col, lfc_col, title) {
  plot_data <- data.frame(
    feature_name = data$feature_name,
    p_value = pmin(-log10(data[[pval_col]]), 60),  
    log2_fold_change = data[[lfc_col]],
    significant = case_when(
      data[[lfc_col]] > 0.4 & data[[pval_col]] < 0.05 ~ "Up",
      data[[lfc_col]] < -0.4 & data[[pval_col]] < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )
  
  x_max <- max(abs(plot_data$log2_fold_change))
  
  genes_to_label <- plot_data %>% filter(feature_name %in% genes_to_label_list)
  
  my_colors <- c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray")
  
  p <- ggplot(plot_data, aes(x = log2_fold_change, y = p_value, color = significant)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = my_colors) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(Adjusted p-value)"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5)
    ) +
    geom_text_repel(
      data = genes_to_label,
      aes(label = feature_name),
      size = 4,
      box.padding = 0.6,
      max.overlaps = Inf,
      segment.color = "grey50",
      min.segment.length = 0.2
    ) +
    xlim(c(-x_max, x_max)) 
  
  return(p)
}

late_plot <- plot_volcano(
  data, 
  "late_adjusted_p_value", 
  "late_log2_fold_change", 
  "TotalRNA - Late vs Control"
)

print(late_plot)

ggsave("totalRNA_late_final.pdf", late_plot, width = 4, height = 4, bg = "white")