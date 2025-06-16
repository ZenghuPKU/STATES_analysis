rm(list = ls())
setwd('/storage/lingyuan2/STATES_data')
library(dplyr)
df = read.csv('te_by_dr_bin_gene_control1021_3bin.csv')

str(df)
print(head(df))

result <- df %>%
  group_by(gene) %>%
  reframe(
    TE_diff = TE[DR_bin == "[0.66, 1.0)"] - mean(TE[DR_bin != "[0.66, 1.0)"], na.rm = TRUE)
  ) %>%
  filter(!is.na(TE_diff))

te_diff_df <- as.data.frame(result)

print(head(te_diff_df))
print(paste("Number of genes:", nrow(te_diff_df)))

te_diff_df <- te_diff_df %>% arrange(desc(TE_diff))

te_diff_df <- te_diff_df %>%
  mutate(gene_order = row_number())

genes_to_label <- unique(c(
  head(te_diff_df$gene, 3),
  tail(te_diff_df$gene, 3),
  "PXDN", "SLC16A1"
))
label_data <- te_diff_df %>% filter(gene %in% genes_to_label)

library(ggrepel)

high_te_genes <- te_diff_df[te_diff_df$TE_diff >= 0.1, "gene", drop = FALSE]
low_te_genes <- te_diff_df[te_diff_df$TE_diff < -0.1, "gene", drop = FALSE]
middle_te_genes <- te_diff_df[te_diff_df$TE_diff > -0.1 & te_diff_df$TE_diff < 0.1, "gene", drop = FALSE]
write.csv(
  low_te_genes,
  file = "/storage/lingyuan2/STATES_data/3bin_lower_genes0_1.csv",
  row.names = FALSE
)
p <- ggplot(te_diff_df, aes(x = gene_order, y = TE_diff)) +
  geom_point(aes(color = ifelse(TE_diff >= 0.1, "high", ifelse(TE_diff < -0.1, "low", "middle"))), alpha = 0.8) +
  scale_color_manual(values = c("low" = "red", "middle" = "gray", "high" = "blue")) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(data = label_data,
                  aes(label = gene),
                  size = 3,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.2,
                  max.overlaps = Inf) +
  scale_y_continuous(breaks = seq(floor(min(te_diff_df$TE_diff)),
                                ceiling(max(te_diff_df$TE_diff)), 
                                by = 0.1),
                    minor_breaks = seq(floor(min(te_diff_df$TE_diff)),
                                     ceiling(max(te_diff_df$TE_diff)),
                                     by = 0.05)) +
  scale_x_continuous(limits = c(0, 1022), 
                    breaks = seq(0, 1022, by = 100),
                    minor_breaks = seq(0, 1022, by = 50)) +
  labs(x = "Genes",
       y = "Cell Edge Bin TE - Other Bins' Average TE",
       title = "TE_diff across genes") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm")) +
  coord_cartesian(clip = "off")
p
ggsave("TE_diff_across_genes.pdf", p, width = 4.5, height = 4, units = "in", dpi = 300)
