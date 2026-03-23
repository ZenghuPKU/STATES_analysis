setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/09_MLG_cell_density")

rm(list = ls())

library(dplyr)
library(ggplot2)
library(rstatix)

#========================
# 1. Read data
#========================
pooled <- read.csv("MLG_8mAD_ring_label3_density_pooled.csv", check.names = FALSE)
outer  <- read.csv("MLG_8mAD_outer_label3_density.csv", check.names = FALSE)

#========================
# 2. Outer 50 µm+ baseline
#========================
outer_ref <- outer %>%
  dplyr::select(states_nn_alg1_label3, density_per_um2) %>%
  dplyr::rename(outer_density = density_per_um2)

#========================
# 4. Pooled density for plotting
#========================
plot_df <- pooled %>%
  dplyr::rename(area_um2 = total_area_um2) %>%
  dplyr::bind_rows(
    outer %>% dplyr::mutate(ring = "outer_50um_plus")
  ) %>%
  dplyr::mutate(
    density_mm2 = density_per_um2 * 1e6,
    ring = factor(ring, levels = c(sort(unique(pooled$ring)), "outer_50um_plus"))
  )

#========================
# 6. Plot and save to PDF with title
#========================
pdf("MLG_8mAD_density_barplot.pdf", width = 5, height = 5)
p <- ggplot(plot_df, aes(x = ring, y = density_mm2, fill = states_nn_alg1_label3)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    y = expression("Cell density (per " * mm^2 * ")"),
    x = "Ring",
    fill = "Cell state (label3)",
    title = "8mAD"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
print(p)
dev.off()