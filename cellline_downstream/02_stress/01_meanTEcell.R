rm(list=ls())
setwd("/storage/lingyuan2/20250101rj/downstream/alldatafinalzarr")

library(Seurat)
library(anndata)
library(reticulate)

use_python("/home/lingyuan2/mambaforge/envs/bigfish_env/bin/python", required = TRUE)
py_config()
sc_ad <- read_h5ad("filtered_data.h5ad")

metadata <- sc_ad$obs

library(ggplot2)
library(ggridges)
library(forcats)

sample_order <- c("C5Tg4h", "C4Tg2h", "B6Tg1h", "B5Tg30min", "B4Tg15min", "C3control")

metadata <- metadata[metadata$sample %in% sample_order, ]
metadata$sample <- factor(metadata$sample, levels = sample_order)

label_map <- c("C5Tg4h" = "Tg4h", 
               "C4Tg2h" = "Tg2h", 
               "B6Tg1h" = "Tg1h", 
               "B5Tg30min" = "Tg30min", 
               "B4Tg15min" = "Tg15min", 
               "C3control" = "control")

metadata$sample_label <- label_map[as.character(metadata$sample)]

library(viridis)
colors <- viridis_pal(option = "D")(length(sample_order))
names(colors) <- sample_order

cat("Sample colors:\n")
for (sample in sample_order) {
  cat(paste0(sample, " (", label_map[sample], "): ", colors[sample], "\n"))
}

p <- ggplot(metadata, aes(x = TE_cell_mean, y = sample, fill = sample)) +
  geom_density_ridges(scale = 1.8, alpha = 0.8) +
  scale_y_discrete(labels = label_map, expand = expansion(mult = c(0, 0.4))) +
  scale_x_continuous(
    limits = c(min(metadata$TE_cell_mean), max(metadata$TE_cell_mean)),
    breaks = seq(0, max(metadata$TE_cell_mean), by = 0.05)
  ) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribution of cell TE by sample", x = "TE_cell_mean", y = "Sample") +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.2, "cm")
  )

p
ggsave("TE_distribution_by_sample.pdf", plot = p, width = 8, height = 6)