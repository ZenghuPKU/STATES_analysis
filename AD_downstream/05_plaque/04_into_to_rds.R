library(Seurat)
library(dplyr)

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/05_plaque")

states <- readRDS("../04_celltypes_identification/states_celltypes_identification.rds")

head(states@meta.data)
table(states$states_nn_alg1_label1)
table(states$states_nn_alg1_label2)
table(states$states_nn_alg1_label3)

plaque_4m <- read.csv("4mAD_mousebrain_with_plaque_info_obs.csv", header = TRUE)
plaque_8m <- read.csv("8mAD_mousebrain_with_plaque_info_obs.csv", header = TRUE)
plaque_14m <- read.csv("14mAD_mousebrain_with_plaque_info_obs.csv", header = TRUE)

plaque_4m$age <- "4m"
plaque_8m$age <- "8m"
plaque_14m$age <- "14m"

plaque_info <- rbind(plaque_4m, plaque_8m, plaque_14m)

plaque_info$cell_barcode_full <- paste0(plaque_info$age, "AD-", plaque_info$cell_barcode)

plaque_info2 <- plaque_info[plaque_info$cell_barcode_full %in% rownames(states@meta.data), ]

md <- states@meta.data
plaque_info2 <- plaque_info2[match(rownames(md), plaque_info2$cell_barcode_full), ]

md$assigned_plaque_id <- plaque_info2$assigned_plaque_id
md$min_border_dist_um <- plaque_info2$min_border_dist_um
md$ring_assignment <- plaque_info2$ring_assignment

states@meta.data <- md

saveRDS(states, file = "states_with_plaque_info.rds")
