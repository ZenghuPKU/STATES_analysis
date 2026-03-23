# Run Plaque Segmentation

# Load required libraries
library(EBImage)

# Get command line arguments or use defaults
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments: image_path, output_dir, sigma, threshold, min_area, save_images, work_dir
# Or use environment variables if arguments not provided
if (length(args) >= 1 && args[1] != "") {
  image_path <- args[1]
} else {
  image_path <- Sys.getenv("IMAGE_PATH", unset = "14mADabeta.tif")
}

if (length(args) >= 2 && args[2] != "") {
  output_dir <- args[2]
} else {
  output_dir <- Sys.getenv("OUTPUT_DIR", unset = "14mADabeta_plaque_results")
}

if (length(args) >= 3 && args[3] != "") {
  sigma <- as.numeric(args[3])
} else {
  sigma <- as.numeric(Sys.getenv("SIGMA", unset = "2"))
}

if (length(args) >= 4 && args[4] != "") {
  threshold <- as.numeric(args[4])
} else {
  threshold <- as.numeric(Sys.getenv("THRESHOLD", unset = "30"))
}

if (length(args) >= 5 && args[5] != "") {
  min_area <- as.numeric(args[5])
} else {
  min_area <- as.numeric(Sys.getenv("MIN_AREA", unset = "400"))
}

if (length(args) >= 6 && args[6] != "") {
  max_eccentricity <- as.numeric(args[6])
} else {
  max_eccentricity <- as.numeric(Sys.getenv("MAX_ECCENTRICITY", unset = "0.7"))
}

if (length(args) >= 7 && args[7] != "" && tolower(args[7]) != "null" && tolower(args[7]) != "na") {
  max_radius_cv <- as.numeric(args[7])
} else {
  max_radius_cv_env <- Sys.getenv("MAX_RADIUS_CV", unset = "")
  if (max_radius_cv_env != "" && tolower(max_radius_cv_env) != "null" && tolower(max_radius_cv_env) != "na") {
    max_radius_cv <- as.numeric(max_radius_cv_env)
  } else {
    max_radius_cv <- NULL
  }
}

if (length(args) >= 8 && args[8] != "" && tolower(args[8]) != "null" && tolower(args[8]) != "na") {
  max_radius_ratio <- as.numeric(args[8])
} else {
  max_radius_ratio_env <- Sys.getenv("MAX_RADIUS_RATIO", unset = "")
  if (max_radius_ratio_env != "" && tolower(max_radius_ratio_env) != "null" && tolower(max_radius_ratio_env) != "na") {
    max_radius_ratio <- as.numeric(max_radius_ratio_env)
  } else {
    max_radius_ratio <- NULL
  }
}

if (length(args) >= 9 && args[9] != "") {
  save_images <- as.logical(args[9])
} else {
  save_images <- as.logical(Sys.getenv("SAVE_IMAGES", unset = "FALSE"))
}

if (length(args) >= 10 && args[10] != "") {
  work_dir <- args[10]
} else {
  work_dir <- Sys.getenv("WORK_DIR", unset = "/storage/lingyuan2/STATES_analysis/AD_downstream/05_plaque")
}

# Set working directory
setwd(work_dir)

# Source the segmentation function
source("plaque_segmentation.R")

# Check if file exists (handle both absolute and relative paths)
if (!file.exists(image_path)) {
  # Try relative to work_dir
  image_path_full <- file.path(work_dir, image_path)
  if (file.exists(image_path_full)) {
    image_path <- image_path_full
  } else {
    stop("Image file not found: ", image_path)
  }
}

# Convert to absolute path for clarity
image_path <- normalizePath(image_path)

message("========================================")
message("Running plaque segmentation...")
message(paste("Working directory:", work_dir))
message(paste("Image file:", image_path))
message(paste("Output directory:", output_dir))
param_str <- paste("sigma =", sigma, ", threshold =", threshold, 
                   ", min_area =", min_area, ", max_eccentricity =", max_eccentricity)
if (!is.null(max_radius_cv)) {
  param_str <- paste(param_str, ", max_radius_cv =", max_radius_cv)
}
if (!is.null(max_radius_ratio)) {
  param_str <- paste(param_str, ", max_radius_ratio =", max_radius_ratio)
}
message(paste("Parameters:", param_str))
message(paste("Save images:", save_images))
message(paste("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
message("========================================")

# Run segmentation
results <- plaque_segmentation(
  image_path = image_path,
  sigma = sigma,
  threshold = threshold,
  min_area = min_area,
  max_eccentricity = max_eccentricity,
  max_radius_cv = max_radius_cv,
  max_radius_ratio = max_radius_ratio,
  show_plots = FALSE,  # Disable plots for background execution
  max_plot_size = 5000 # Maximum dimension for plotting
)

# Save results
save_plaque_results(results, output_dir, save_images = save_images)

# Print summary
message("\n========================================")
message("=== FINAL RESULTS ===")
message("========================================")
message(paste("Total plaques detected:", nrow(results$plaque_features)))
message(paste("Plaques after filtering (area >", min_area, "px, eccentricity <", max_eccentricity, "):", nrow(results$plaque_features_filtered)))

if (nrow(results$plaque_features_filtered) > 0) {
  message("\nPlaque statistics:")
  message(paste("  Mean area:", round(mean(results$plaque_features_filtered$s.area), 2), "pixels"))
  message(paste("  Median area:", round(median(results$plaque_features_filtered$s.area), 2), "pixels"))
  area_range <- range(results$plaque_features_filtered$s.area)
  message(paste("  Area range:", round(area_range[1], 2), "-", round(area_range[2], 2), "pixels"))
}

message("\n========================================")
message(paste("Job completed successfully at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
message("========================================")
