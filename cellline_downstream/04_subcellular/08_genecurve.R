rm(list = ls())
library(tidyverse)
library(ggplot2)
setwd('/storage/lingyuan2/STATES_data')

df <- read.csv("te_by_dr_bin_gene_control1021_3bin.csv")
result <- df %>%
  group_by(gene) %>%
  reframe(
    TE_diff = TE[DR_bin == "[0.66, 1.0)"] - mean(TE[DR_bin != "[0.66, 1.0)"], na.rm = TRUE)
  ) %>%
  filter(!is.na(TE_diff))
te_diff_df <- as.data.frame(result)
te_diff_df <- te_diff_df %>% arrange(desc(TE_diff))
te_diff_df <- te_diff_df %>%
  mutate(gene_order = row_number())
top_genes <- head(te_diff_df$gene, 3)
bottom_genes <- tail(te_diff_df$gene, 3)

top_plot_data <- df %>%
  filter(gene %in% top_genes)

bottom_plot_data <- df %>%
  filter(gene %in% bottom_genes)
p1 <- ggplot(top_plot_data, aes(x = DR_bin, y = TE, color = gene, group = gene)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text()
  ) +
  scale_y_continuous(limits = c(0.25, 0.65), breaks = seq(0.25, 0.65, 0.1)) +
  labs(
    x = "DR bin", 
    y = "Translation Efficiency (TE)",
    color = "Gene"
  )
p1
p2 <- ggplot(bottom_plot_data, aes(x = DR_bin, y = TE, color = gene, group = gene)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text()
  ) +
  scale_y_continuous(limits = c(0.25, 0.65), breaks = seq(0.25, 0.65, 0.1)) +
  labs(
    x = "DR bin",
    y = "Translation Efficiency (TE)",
    color = "Gene"
  )
p2

ggsave("te_by_dr_bin_extreme_genes_top3.pdf", plot = p1, width = 4, height = 3)
ggsave("te_by_dr_bin_extreme_genes_bottom3.pdf", plot = p2, width = 4, height = 3)
