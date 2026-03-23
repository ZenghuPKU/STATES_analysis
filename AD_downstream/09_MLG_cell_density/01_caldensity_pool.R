setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/09_MLG_cell_density")

library(Seurat)
library(dplyr)
library(tidyr)
library(readr)

states_rds <- "../05_plaque/states_with_plaque_info.rds"
area_csv <- "../05_plaque/14mAD_plaque_ring_area_and_cells.csv"

rings_to_use <- c(
  "ring1_0_10um",
  "ring2_10_20um",
  "ring3_20_30um",
  "ring4_30_40um",
  "ring5_40_50um"
)
outer_ring <- "outer_50um_plus"

states <- readRDS(states_rds)
meta <- states@meta.data

meta_mlg <- meta %>%
  filter(type == "14mAD", states_nn_alg1_label2 == "MLG")

area_df <- read_csv(area_csv, show_col_types = FALSE)
colnames(area_df)[colnames(area_df) == "plaque_id"] <- "assigned_plaque_id"


all_rings_cells <- meta_mlg %>%
  filter(ring_assignment %in% rings_to_use)

ring_label3_counts <- all_rings_cells %>%
  group_by(ring = ring_assignment, states_nn_alg1_label3) %>%
  summarise(cell_count = n(), .groups = "drop")

ring_areas <- area_df %>%
  filter(ring %in% rings_to_use) %>%
  group_by(ring) %>%
  summarise(total_area_um2 = sum(area_um2, na.rm = TRUE), .groups = "drop")

ring_label3_density <- left_join(ring_label3_counts, ring_areas, by = "ring") %>%
  mutate(density_per_um2 = cell_count / total_area_um2)

write.csv(ring_label3_density,
          file = "MLG_14mAD_ring_label3_density_pooled.csv",
          row.names = FALSE)

outer_cells <- meta_mlg %>%
  filter(ring_assignment == outer_ring)

outer_counts <- outer_cells %>%
  group_by(states_nn_alg1_label3) %>%
  summarise(cell_count = n(), .groups = "drop")

outer_area_df <- area_df %>%
  filter(ring == outer_ring)
total_outer_area <- sum(outer_area_df$area_um2, na.rm = TRUE)

outer_counts <- outer_counts %>%
  mutate(total_area_um2 = total_outer_area,
         density_per_um2 = cell_count / total_area_um2)

write.csv(outer_counts,
          file = "MLG_14mAD_outer_label3_density.csv",
          row.names = FALSE)

