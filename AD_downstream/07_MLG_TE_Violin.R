# 07_MLG_TE_Violin
rm(list = ls())
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggpubr)   
setwd("~/AD_downstream/")

states <- readRDS("states_with_plaque_info.rds")
meta <- states@meta.data

total_counts <- GetAssayData(states, assay = "totalRNA", slot = "counts")
rb_counts    <- GetAssayData(states, assay = "rbRNA",    slot = "counts")

total_sum <- Matrix::colSums(total_counts)
rb_sum    <- Matrix::colSums(rb_counts)

te_cell <- rb_sum / total_sum
te_cell[!is.finite(te_cell)] <- NA_real_

df_te <- data.frame(
  cell   = colnames(states),
  label3 = meta$states_nn_alg1_label3,
  TE     = as.numeric(te_cell),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(TE) & !is.na(label3))

mlg_levels <- c("MLG1", "MLG2", "MLG3")

df_mlg <- df_te %>%
  dplyr::filter(label3 %in% mlg_levels) %>%
  dplyr::mutate(label3 = factor(label3, levels = mlg_levels))
print(table(df_mlg$label3))

comp_mlg <- list(
  c("MLG1", "MLG2"),   
  c("MLG2", "MLG3"),
  c("MLG1", "MLG3")
)

mlg_tests <- lapply(comp_mlg, function(cc) {
  x <- df_mlg$TE[df_mlg$label3 == cc[1]]
  y <- df_mlg$TE[df_mlg$label3 == cc[2]]
  t_res <- t.test(x, y, var.equal = FALSE)
  data.frame(
    group1  = cc[1],
    group2  = cc[2],
    p.value = t_res$p.value,
    stringsAsFactors = FALSE
  )
})
mlg_tests <- do.call(rbind, mlg_tests)
mlg_tests$FDR <- p.adjust(mlg_tests$p.value, method = "BH")
print(mlg_tests)

p_mlg <- ggplot(df_mlg, aes(x = label3, y = TE, fill = label3)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    alpha = 0.8
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "TE per cell in MLG1–3",
    x = "MLG subtype (label3)",
    y = "TE = sum(rbRNA) / sum(totalRNA)",
    fill = "MLG"
  ) +
  ggpubr::stat_compare_means(
    comparisons       = comp_mlg,
    method            = "t.test",
    p.adjust.method   = "BH",
    label             = "p.signif"
  )

ggsave("TE_violin_MLG1-3_box.pdf", p_mlg, width = 5, height = 4)