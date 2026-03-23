# 01_Annotated_celltypes_UMAP_visualization
rm(list = ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("~/AD_downstream/06_Annotated_celltypes_UMAP_visualization")

states <- readRDS("states_with_plaque_info.rds")

umap_coords <- Embeddings(states, reduction = "states.umap")
umap_df <- data.frame(
  cell_id = rownames(umap_coords),
  umap_coords,
  label1 = states@meta.data$states_nn_alg1_label1,
  label2 = states@meta.data$states_nn_alg1_label2,
  label3 = states@meta.data$states_nn_alg1_label3
)
write.csv(umap_df, file = "5xFAD_mousebrain_umap_embeddings_with_labels.csv", row.names = FALSE)

Idents(states) <- "states_nn_alg1_label1"
table(states@active.ident)

celltype_colors_l1 <- c(
  "Non_Neuron" = "#66C2A5",
  "Neuron"     = "#F781BF",
  "Mix"        = "#FAFAFA"
)

P1_pdf <- DimPlot(
  states,
  reduction = "states.umap",
  label     = TRUE,
  group.by  = "states_nn_alg1_label1"
) +
  scale_color_manual(values = celltype_colors_l1) +
  NoLegend() +
  ggtitle(NULL)

ggsave(
  "states_Umap_label1.pdf",
  P1_pdf,
  width  = 8,
  height = 8,
  device = cairo_pdf
)

P1_png <- DimPlot(
  states,
  reduction = "states.umap",
  label     = FALSE,
  group.by  = "states_nn_alg1_label1"
) +
  scale_color_manual(values = celltype_colors_l1) +
  ggtitle(NULL)

P1_png_clean <- P1_png +
  theme_void() +
  theme(legend.position = "none")   # ★ 强制去 legend

ggsave(
  "states_Umap_label1.png",
  P1_png_clean,
  width  = 8,
  height = 8,
  units  = "in",
  dpi    = 300,
  bg     = "transparent"
)

Idents(states) <- "states_nn_alg1_label2"
table(states@active.ident)

celltype_colors_l2 <- c(
  "OPC"       = "#667872",
  "DE/MEN"    = "#b274e8",
  "CHO/PEP"   = "#97f4f7",
  "OLG"       = "#e6db17",
  "CHOR/EPEN" = "#F1B2E9",
  "TEPN"      = "#6C9F35",
  "MLG"       = "#8597c6",
  "VAS"       = "#17BED0",
  "INH"       = "#3182BD",
  "AC"        = "#FDC06F",
  "Mix"       = "#FAFAFA"
)

P2_pdf <- DimPlot(
  states,
  reduction = "states.umap",
  label     = TRUE,
  group.by  = "states_nn_alg1_label2"
) +
  scale_color_manual(values = celltype_colors_l2) +
  NoLegend() +
  ggtitle(NULL)

ggsave(
  "states_Umap_label2.pdf",
  P2_pdf,
  width  = 8,
  height = 8,
  device = cairo_pdf
)

P2_png <- DimPlot(
  states,
  reduction = "states.umap",
  label     = FALSE,
  group.by  = "states_nn_alg1_label2"
) +
  scale_color_manual(values = celltype_colors_l2) +
  ggtitle(NULL)

P2_png_clean <- P2_png +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "states_Umap_label2.png",
  P2_png_clean,
  width  = 8,
  height = 8,
  units  = "in",
  dpi    = 300,
  bg     = "transparent"
)

Idents(states) <- "states_nn_alg1_label3"
table(states@active.ident)

celltype_colors_l3 <- c(
  "AC1"          = "#ccba33",
  "AC2"          = "#ffbe85",
  "AC3"          = "#e3782b",
  "AC4"          = "#f0dea3",
  "CHOR"         = "#7f52a9",
  "EPEN"         = "#c4b0d4",
  "CHO/PEP"      = "#97f4f7",
  "INH_Sst"      = "#96abeb",
  "INH_Pvalb"    = "#96cad4",
  "INH_Cnr1_Vip" = "#a8e1eb",
  "MLG1"         = "#A9B4D8",
  "MLG2"         = "#8597C6",
  "MLG3"         = "#445172",
  "OPC"          = "#667872",
  "OLG1"         = "#e4f768",
  "OLG2"         = "#e6db17",
  "VLMC"         = "#1f76b3",
  "VSMC"         = "#00aeef",
  "Peri/VEC"     = "#d3a59c",
  "DE/MEN"       = "#b274e8",
  "MSN"          = "#D96DA1",
  "DGGRC"        = "#a6e8a6",
  "TEGLU CA1"    = "#77ed8f",
  "TEGLU CA2"    = "#82ad2d",
  "TEGLU CA3"    = "#28330b",
  "TEGLU L2/3"   = "#cbfc60",
  "TEGLU L4"     = "#96db00",
  "TEGLU L4/5"   = "#04b361",
  "TEGLU L5"     = "#40d102",
  "TEGLU L6"     = "#32a630",
  "TEGLU Mix"    = "#c5fcc5",
  "Mix"          = "#FAFAFA"
)

P3_pdf <- DimPlot(
  states,
  reduction = "states.umap",
  label     = FALSE,
  group.by  = "states_nn_alg1_label3"
) +
  scale_color_manual(values = celltype_colors_l3) +
  NoLegend() +
  ggtitle(NULL)

ggsave(
  "states_Umap_label3.pdf",
  P3_pdf,
  width  = 8,
  height = 8,
  device = cairo_pdf
)

P3_png <- DimPlot(
  states,
  reduction = "states.umap",
  label     = FALSE,
  group.by  = "states_nn_alg1_label3"
) +
  scale_color_manual(values = celltype_colors_l3) +
  ggtitle(NULL)

P3_png_clean <- P3_png +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "states_Umap_label3.png",
  P3_png_clean,
  width  = 8,
  height = 8,
  units  = "in",
  dpi    = 300,
  bg     = "transparent"
)

table(states$type)

states$group_ADWT <- NA_character_
states$group_ADWT[grepl("AD", states$type)] <- "AD"
states$group_ADWT[grepl("WT", states$type)] <- "WT"

states$group_color <- ifelse(
  states$states_nn_alg1_label1 == "Mix",
  "Mix",
  states$group_ADWT
)
states$group_color <- factor(states$group_color,
                             levels = c("WT", "AD", "Mix"))

group_colors_adwt <- c(
  "WT"  = "#00aeef",
  "AD"  = "red",
  "Mix" = "#FAFAFA"
)

P_ADWT_pdf <- DimPlot(
  states,
  reduction = "states.umap",
  group.by  = "group_color",
  label     = FALSE
) +
  scale_color_manual(values = group_colors_adwt) +
  NoLegend() +
  ggtitle(NULL)

ggsave(
  "states_Umap_ADWT.pdf",
  P_ADWT_pdf,
  width  = 8,
  height = 8,
  device = cairo_pdf
)

P_ADWT_png <- DimPlot(
  states,
  reduction = "states.umap",
  group.by  = "group_color",
  label     = FALSE
) +
  scale_color_manual(values = group_colors_adwt) +
  ggtitle(NULL)

P_ADWT_png_clean <- P_ADWT_png +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "states_Umap_ADWT.png",
  P_ADWT_png_clean,
  width  = 8,
  height = 8,
  units  = "in",
  dpi    = 300,
  bg     = "transparent"
)
