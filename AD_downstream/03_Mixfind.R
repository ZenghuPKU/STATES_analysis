# 03_Mixfind
rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reticulate)
library(Seurat)
library(grid)
setwd("~/AD_downstream/")

states <- readRDS("~/AD_downstream/states_ADcombined_harmony.rds")

states@meta.data$distance2centroid_leiden_teglu <- 9999

umap_coords <- Embeddings(states, "states.umap")
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df$states_nn_alg1 <- states@meta.data$states_nn_alg1

centroids <- umap_df %>%
  group_by(states_nn_alg1) %>%
  summarise(centroid_x = mean(UMAP1),
            centroid_y = mean(UMAP2))

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = states_nn_alg1)) +
  geom_point(size = 1) +
  geom_point(data = centroids, aes(x = centroid_x, y = centroid_y),
             color = "red", size = 3) +
  theme_minimal() +
  ggtitle("UMAP with Cluster Centroids")
print(p_umap)
ggsave("UMAP_with_Centroids.png", p_umap, width = 8, height = 6, dpi = 300)

unique_labels <- sort(unique(states@meta.data$states_nn_alg1))
for(i in seq_along(unique_labels)){
  current_label <- unique_labels[i]
  safe_label <- gsub("/", "_", current_label)  
  current_centroid <- centroids %>% filter(states_nn_alg1 == current_label)
  
  idx <- which(states@meta.data$states_nn_alg1 == current_label)
  umap_current <- umap_coords[idx, , drop = FALSE]
  
  dm <- sqrt((umap_current[, 1] - current_centroid$centroid_x)^2 +
               (umap_current[, 2] - current_centroid$centroid_y)^2)
  
  states@meta.data$distance2centroid_leiden_teglu[idx] <- dm
  
  p_hist <- ggplot(data.frame(distance = dm), aes(x = distance)) +
    geom_histogram(binwidth = diff(range(dm)) / 30, fill = "grey", color = "black") +
    ggtitle(paste("Distance Distribution for Cluster", current_label))
  ggsave(filename = paste0("Cluster_", safe_label, "_hist.png"), 
         plot = p_hist, width = 8, height = 6, dpi = 300)
  
  p_left <- ggplot() +
    geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2),
               color = "#111111", size = 0.5) +
    geom_point(data = umap_df[umap_df$states_nn_alg1 == current_label, ],
               aes(x = UMAP1, y = UMAP2), size = 0.5) +
    geom_point(data = current_centroid,
               aes(x = centroid_x, y = centroid_y),
               color = "red", size = 3) +
    ggtitle(paste("Cluster", current_label, "UMAP"))
  
  umap_current_df <- data.frame(UMAP1 = umap_current[, 1],
                                UMAP2 = umap_current[, 2],
                                distance = dm)
  p_right <- ggplot() +
    geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2),
               color = "#111111", size = 0.3) +
    geom_point(data = umap_current_df,
               aes(x = UMAP1, y = UMAP2, color = distance),
               size = 0.3) +
    scale_color_viridis_c() +
    geom_point(data = current_centroid,
               aes(x = centroid_x, y = centroid_y),
               color = "red", size = 3) +
    ggtitle(paste("Cluster", current_label, "Colored by Distance"))
  
  combined_plot <- arrangeGrob(p_left, p_right, ncol = 1)
  ggsave(filename = paste0("Cluster_", safe_label, "_umap_combined.png"),
         plot = combined_plot, width = 8, height = 12, dpi = 300)
}

p_overall <- ggplot(states@meta.data, aes(x = distance2centroid_leiden_teglu)) +
  geom_histogram(binwidth = diff(range(states@meta.data$distance2centroid_leiden_teglu)) / 30,
                 fill = "grey", color = "black") +
  ggtitle("Overall Distance Distribution")
print(p_overall)

manual_threshold <- c(5,5,5,4,4,
                      6,6,5,4,4,
                      5,3,3,5,3,
                      6,5,5,7.5,5,
                      3,4,2.5,3,2.5,4
)

states@meta.data$is_mix_teglu <- "False"

for(i in seq_along(unique_labels)){
  current_label <- unique_labels[i]
  safe_label <- gsub("/", "_", current_label)
  
  idx <- which(states@meta.data$states_nn_alg1 == current_label)
  current_distances <- states@meta.data$distance2centroid_leiden_teglu[idx]
  
  p_cluster <- ggplot(data.frame(distance = current_distances), aes(x = distance)) +
    geom_histogram(binwidth = diff(range(current_distances)) / 30,
                   fill = "grey", color = "black") +
    geom_vline(xintercept = manual_threshold[i], color = "red") +
    ggtitle(paste("Cluster", current_label, "Distance Distribution",
                  "\nThreshold =", manual_threshold[i]))
  ggsave(filename = paste0("Cluster_", safe_label, "_distance.png"),
         plot = p_cluster, width = 8, height = 6, dpi = 300)
  
  current_centroid <- centroids %>% dplyr::filter(states_nn_alg1 == current_label)
  
  theta <- seq(0, 2 * pi, length.out = 200)
  circle_df <- data.frame(
    x = current_centroid$centroid_x + manual_threshold[i] * cos(theta),
    y = current_centroid$centroid_y + manual_threshold[i] * sin(theta)
  )
  
  umap_current_df <- umap_df[umap_df$states_nn_alg1 == current_label, ]
  
  p_umap_circle <- ggplot() +
    
    geom_point(data = umap_df,
               aes(x = UMAP1, y = UMAP2),
               color = "#dddddd", size = 0.3) +
    
    geom_point(data = umap_current_df,
               aes(x = UMAP1, y = UMAP2),
               color = "blue", size = 0.4) +
    
    geom_path(data = circle_df,
              aes(x = x, y = y),
              color = "red", size = 0.6) +
    
    geom_point(data = current_centroid,
               aes(x = centroid_x, y = centroid_y),
               color = "red", size = 2) +
    ggtitle(paste("Cluster", current_label,
                  "- UMAP with Radius =", manual_threshold[i]))
  
  ggsave(filename = paste0("Cluster_", safe_label, "_umap_threshold_circle.png"),
         plot = p_umap_circle, width = 6, height = 6, dpi = 300)
  
  states@meta.data$is_mix_teglu[idx][ current_distances > manual_threshold[i] ] <- "True"
}

states@meta.data$is_mix_teglu <- factor(states@meta.data$is_mix_teglu)

umap_df$is_mix_teglu <- states@meta.data$is_mix_teglu
p_final <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = is_mix_teglu)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("False" = "blue", "True" = "red")) +
  ggtitle("Combined UMAP: Mix Cell Classification")
print(p_final)
ggsave(filename = "Combined_UMAP_Mix_Classification.png", plot = p_final,
       width = 8, height = 6, dpi = 300)

states@meta.data$states_nn_alg1_new <- as.character(states@meta.data$states_nn_alg1)
states@meta.data$states_nn_alg1_new[states@meta.data$is_mix_teglu == "True"] <- "Mix"
states@meta.data$states_nn_alg1_new <- factor(states@meta.data$states_nn_alg1_new)
umap_df$states_nn_alg1_new <- states@meta.data$states_nn_alg1_new

library(RColorBrewer)
all_levels <- levels(umap_df$states_nn_alg1_new)
non_mix_levels <- all_levels[all_levels != "Mix"]
non_mix_colors <- brewer.pal(n = max(3, length(non_mix_levels)), name = "Set1")[1:length(non_mix_levels)]
color_vector <- setNames(non_mix_colors, non_mix_levels)
color_vector["Mix"] <- "lightgrey"

p_updated_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = states_nn_alg1_new)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = color_vector) +
  ggtitle("UMAP with Updated states_nn_alg1 (Mix as New Class)") +
  theme_minimal()
print(p_updated_umap)
ggsave(filename = "Updated_UMAP_states_nn_alg1_Lightgrey.png", 
       plot = p_updated_umap, width = 8, height = 6, dpi = 300)
table(states@meta.data$states_nn_alg1_new)

saveRDS(states, file = "states_mixsfind.rds")

