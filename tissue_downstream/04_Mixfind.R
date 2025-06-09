# 04_Mixfind
# Load libraries and set environment
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reticulate)
library(Seurat)
library(grid)
setwd("~/tissue_downstream/")
rm(list = ls())

# Load the Normalized and Visualized object (states_Normalized_Visualized_output.RData)
load("states_Normalized_Visualized_output.RData")

# Initialize a column for storing distance from cluster centroid
states@meta.data$distance2centroid_leiden_teglu <- 9999

# UMAP Coordinates & Centroid Calculation
umap_coords <- Embeddings(states, "states.umap")
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Use the existing cluster assignment (e.g. "states_nn_alg1") for grouping
umap_df$states_nn_alg1 <- states@meta.data$states_nn_alg1

# Compute centroids for each cluster
centroids <- umap_df %>%
  group_by(states_nn_alg1) %>%
  summarise(centroid_x = mean(UMAP1),
            centroid_y = mean(UMAP2))

# Plot UMAP with cluster centroids
p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = states_nn_alg1)) +
  geom_point(size = 1) +
  geom_point(data = centroids, aes(x = centroid_x, y = centroid_y),
             color = "red", size = 3) +
  theme_minimal() +
  ggtitle("UMAP with Cluster Centroids")
print(p_umap)
ggsave("UMAP_with_Centroids.png", p_umap, width = 8, height = 6, dpi = 300)

# Compute Distance to Centroid for Each Cluster & Plot
unique_labels <- sort(unique(states@meta.data$states_nn_alg1))
for(i in seq_along(unique_labels)){
  current_label <- unique_labels[i]
  safe_label <- gsub("/", "_", current_label)  # replace "/" with "_" for filename
  current_centroid <- centroids %>% filter(states_nn_alg1 == current_label)
  
  # Get indices and UMAP coordinates for cells in current cluster
  idx <- which(states@meta.data$states_nn_alg1 == current_label)
  umap_current <- umap_coords[idx, , drop = FALSE]
  
  # Calculate Euclidean distances to the centroid
  dm <- sqrt((umap_current[, 1] - current_centroid$centroid_x)^2 +
               (umap_current[, 2] - current_centroid$centroid_y)^2)
  
  # Update meta.data with the computed distances
  states@meta.data$distance2centroid_leiden_teglu[idx] <- dm
  
  # Plot distance histogram for the current cluster
  p_hist <- ggplot(data.frame(distance = dm), aes(x = distance)) +
    geom_histogram(binwidth = diff(range(dm)) / 30, fill = "grey", color = "black") +
    ggtitle(paste("Distance Distribution for Cluster", current_label))
  ggsave(filename = paste0("Cluster_", safe_label, "_hist.png"), 
         plot = p_hist, width = 8, height = 6, dpi = 300)
  
  # Create UMAP plots:
  # Left: all cells with current cluster and centroid highlighted
  p_left <- ggplot() +
    geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2),
               color = "#111111", size = 0.5) +
    geom_point(data = umap_df[umap_df$states_nn_alg1 == current_label, ],
               aes(x = UMAP1, y = UMAP2), size = 0.5) +
    geom_point(data = current_centroid,
               aes(x = centroid_x, y = centroid_y),
               color = "red", size = 3) +
    ggtitle(paste("Cluster", current_label, "UMAP"))
  
  # Right: current cluster colored by distance
  umap_current_df <- data.frame(UMAP1 = umap_current[, 1],
                                UMAP2 = umap_current[, 2],
                                distance = dm)
  p_right <- ggplot() +
    geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2),
               color = "#111111", size = 0.5) +
    geom_point(data = umap_current_df,
               aes(x = UMAP1, y = UMAP2, color = distance),
               size = 0.5) +
    scale_color_viridis_c() +
    geom_point(data = current_centroid,
               aes(x = centroid_x, y = centroid_y),
               color = "red", size = 3) +
    ggtitle(paste("Cluster", current_label, "Colored by Distance"))
  
  combined_plot <- arrangeGrob(p_left, p_right, ncol = 2)
  ggsave(filename = paste0("Cluster_", safe_label, "_umap_combined.png"),
         plot = combined_plot, width = 12, height = 6, dpi = 300)
}

# Overall Distance Distribution & Mixed Cell Classification
# Plot overall distance distribution for all cells
p_overall <- ggplot(states@meta.data, aes(x = distance2centroid_leiden_teglu)) +
  geom_histogram(binwidth = diff(range(states@meta.data$distance2centroid_leiden_teglu)) / 30,
                 fill = "grey", color = "black") +
  ggtitle("Overall Distance Distribution")
print(p_overall)
ggsave(filename = "Overall_Distance_Distribution.png", plot = p_overall,
       width = 8, height = 6, dpi = 300)

# Define manual thresholds for each cluster (order must match sorted unique_labels)
manual_threshold <- c(5,5,7,6,5,4,
                      5,10,5,5,5,5,
                      5,5,5,4,5,4,
                      5,4,0,4,5,3,
                      5)

# Initialize mix cell flag as "False"
states@meta.data$is_mix_teglu <- "False"

# For each cluster, plot histogram with threshold and mark cells above threshold as "True"
for(i in seq_along(unique_labels)){
  current_label <- unique_labels[i]
  safe_label <- gsub("/", "_", current_label)
  
  idx <- which(states@meta.data$states_nn_alg1 == current_label)
  current_distances <- states@meta.data$distance2centroid_leiden_teglu[idx]
  
  p_cluster <- ggplot(data.frame(distance = current_distances), aes(x = distance)) +
    geom_histogram(binwidth = diff(range(current_distances)) / 30,
                   fill = "grey", color = "black") +
    geom_vline(xintercept = manual_threshold[i], color = "red") +
    ggtitle(paste("Cluster", current_label, "Distance Distribution"))
  ggsave(filename = paste0("Cluster_", safe_label, "_distance.png"),
         plot = p_cluster, width = 8, height = 6, dpi = 300)
  
  # Mark cells with distance greater than the manual threshold as "True"
  states@meta.data$is_mix_teglu[idx][ current_distances > manual_threshold[i] ] <- "True"
}
states@meta.data$is_mix_teglu <- factor(states@meta.data$is_mix_teglu)

# Plot combined UMAP colored by mix classification
umap_df$is_mix_teglu <- states@meta.data$is_mix_teglu
p_final <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = is_mix_teglu)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("False" = "blue", "True" = "red")) +
  ggtitle("Combined UMAP: Mix Cell Classification")
print(p_final)
ggsave(filename = "Combined_UMAP_Mix_Classification.png", plot = p_final,
       width = 8, height = 6, dpi = 300)

# Update Cluster Labels with Mix Classification
# Update the original cluster assignment: for cells flagged as "True", set label to "Mix"
states@meta.data$states_nn_alg1_new <- as.character(states@meta.data$states_nn_alg1)
states@meta.data$states_nn_alg1_new[states@meta.data$is_mix_teglu == "True"] <- "Mix"
states@meta.data$states_nn_alg1_new <- factor(states@meta.data$states_nn_alg1_new)
umap_df$states_nn_alg1_new <- states@meta.data$states_nn_alg1_new

# Define colors: non-"Mix" classes use Set1 palette; "Mix" is lightgrey
library(RColorBrewer)
all_levels <- levels(umap_df$states_nn_alg1_new)
non_mix_levels <- all_levels[all_levels != "Mix"]
non_mix_colors <- brewer.pal(n = max(3, length(non_mix_levels)), name = "Set1")[1:length(non_mix_levels)]
color_vector <- setNames(non_mix_colors, non_mix_levels)
color_vector["Mix"] <- "lightgrey"

# Plot updated UMAP with new classification
p_updated_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = states_nn_alg1_new)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = color_vector) +
  ggtitle("UMAP with Updated states_nn_alg1 (Mix as New Class)") +
  theme_minimal()
print(p_updated_umap)
ggsave(filename = "Updated_UMAP_states_nn_alg1_Lightgrey.png", 
       plot = p_updated_umap, width = 8, height = 6, dpi = 300)
table(states@meta.data$states_nn_alg1_new)

# Save the updated object
save(states, file = "states_mixsfind.Rdata")
