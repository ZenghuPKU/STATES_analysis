rm(list=ls())

library(ggplot2)
library(ggpubr)

te_by_dr_bin <- read.csv("te_by_dr_bin_cell_1021genes.csv")

suffixes <- c('-0', '-1', '-2', '-3', '-4', '-5')
titles <- c('Control', 'Tg15min', 'Tg30min', 'Tg1h', 'Tg2h', 'Tg4h')

for (i in 1:length(suffixes)) {
  te_by_dr_bin_filtered <- te_by_dr_bin[grep(paste0(suffixes[i], "$"), te_by_dr_bin$cell_idx), ]
  
  bins <- sort(unique(te_by_dr_bin_filtered$DR_bin))
  
  comparisons <- list()
  for(j in 1:(length(bins)-1)) {
    for(k in (j+1):length(bins)) {
      comparisons[[length(comparisons) + 1]] <- c(bins[j], bins[k])
    }
  }
  
  p <- ggplot(te_by_dr_bin_filtered, aes(x = DR_bin, y = TE, fill = DR_bin)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    scale_fill_manual(values = c("#1f78b4", "#a6cee3", "#b2df8a")) +
    stat_compare_means(comparisons = comparisons,
                      method = "t.test", 
                      label = "p.signif",
                      p.adjust.method = "bonferroni",
                      hide.ns = FALSE) +
    coord_cartesian(ylim = c(0, 0.8), clip = "off") +
    scale_y_continuous(limits = c(0, 1.0)) +
    labs(title = paste("TE Distribution (", titles[i], ")"),
         x = "DR Bins",
         y = "TE") +
    theme_classic() +
    theme(axis.text.x = element_text(),
          legend.position = "none",
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  
  ggsave(paste0("violin_", titles[i], ".pdf"), p, width = 5, height = 5)
}
