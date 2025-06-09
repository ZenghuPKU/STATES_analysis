rm(list=ls())
setwd('/media/zenglab/result/lingyuan/mousebraindownstream/soma_process')
fcdfraw <- read.csv("eachgenefinal_all0503.csv")
fcdf <- fcdfraw[fcdfraw$adj_pvalue < 0.05,]
fcdf_top <- fcdf[fcdf$log2fc < -1,]

write.csv(fcdf_top, "fcdf_1.csv", row.names = FALSE)

###########GO

background_genes <- read.csv("2323detected_genes.csv")
background_genes <- unique(background_genes$Gene)

library(clusterProfiler)
library(org.Mm.eg.db)

top_genes <- fcdf_top$gene

min_goterm_size <- 50
max_goterm_size <- 500

perform_go_analysis <- function(genes, ont, minSize = 10, maxSize = 500) {
  enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    universe = background_genes,
    minGSSize = minSize,
    maxGSSize = maxSize
  )
}

top_go_bp <- perform_go_analysis(top_genes, "BP", min_goterm_size, max_goterm_size)
top_go_mf <- perform_go_analysis(top_genes, "MF", min_goterm_size, max_goterm_size)
top_go_cc <- perform_go_analysis(top_genes, "CC", min_goterm_size, max_goterm_size)

library(ggplot2)
library(dplyr)
library(tidyr)

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
    slice_head(n = 3)
  return(go_select)
}

create_combined_plot <- function(go_bp, go_mf, go_cc, title) {
  go_bp_clean <- process_go_data(go_bp, "Biological Process")
  go_mf_clean <- process_go_data(go_mf, "Molecular Function")
  go_cc_clean <- process_go_data(go_cc, "Cellular Component")
  category_order <- c("Cellular Component", "Molecular Function", "Biological Process")
  go_combined <- bind_rows(go_bp_clean, go_mf_clean, go_cc_clean)%>% na.omit()
  go_combined$Category <- factor(go_combined$Category, levels = category_order)
  go_sorted <- data.frame()
  for(cat in category_order) {
    cat_data <- go_combined %>%
      filter(Category == cat) %>%
      arrange(pval_log)
    go_sorted <- bind_rows(go_sorted, cat_data)
  }
  go_combined <- go_sorted
  go_combined$Description <- factor(go_combined$Description, levels = go_combined$Description)
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

p_top <- create_combined_plot(top_go_bp, top_go_mf, top_go_cc, "Somata-biased translation")
p_top
ggsave("somata_biased_all_GO.pdf", plot = p_top, device = "pdf", width = 12, height = 4)

###########volcano plot

library(ggplot2)
library(dplyr)
library(ggrepel)

fcdfraw$Group <- "Other"
fcdfraw$Group[fcdfraw$gene %in% fcdf_top$gene] <- "Somata-biased"

genes_to_label <- c("olfm1", "kcna6", "cnr1", "mrpl43", "smoc1", "spink8", "tcf20")
label_genes <- fcdfraw %>% filter(tolower(gene) %in% genes_to_label)

library(ggplot2)
library(ggrepel)

fcdfraw$plot_y <- -log10(fcdfraw$adj_pvalue)
fcdfraw$plot_y[fcdfraw$plot_y > 200] <- 205

ymax <- max(fcdfraw$plot_y, na.rm = TRUE)

p <- ggplot(fcdfraw, aes(x = log2fc, y = plot_y, color = Group)) +
  geom_point(alpha = 0.7, size = 1.5) +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 2.5,
    max.overlaps = Inf,
    segment.size = 0.3,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#2166AC", linewidth = 0.6) +
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = seq(-2, 2, by = 0.5)
  ) +
  scale_color_manual(
    values = c(
      "Somata-biased" = "#B2182B",
      "Other" = "gray"
    ),
    breaks = c("Somata-biased top 10%", "Other")
  ) +
  labs(
    x = "log2(ProcessesTE / SomataTE)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 12),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 200))
p
ggsave("volcano_plot_all_label.pdf", p, width = 4, height = 4, device = "pdf")



###########violin plot
gene_summary <- read.csv('soma_process_eachgene_all_0503.csv')

row.names(gene_summary) <- gene_summary$gene
gene_summary$gene <- NULL
library(ggpubr)
te_data <- data.frame(
  TE = c(gene_summary$somata_TE, gene_summary$processes_TE),
  Region = rep(c("Somata", "Processes"), each = nrow(gene_summary))
)
library(ggplot2)
library(ggpubr)

te_data$Region <- factor(te_data$Region, levels = c("Somata", "Processes"))

p <- ggplot(te_data, aes(x = Region, y = TE, fill = Region)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("Somata" = "#B2182B", "Processes" = "#2166AC")) +
  labs(
    x = "Region",
    y = "TE Value"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("Somata", "Processes")),
    label = "p.format",
    label.y = max(te_data$TE) * 1.1,
    size = 5
  )

print(p)
ggsave(
  filename = "violin_plot_all.pdf",
  plot = p,
  device = "pdf",
  width = 4,
  height = 4,
  units = "in",
  dpi = 300
)
label_genes['gene']
