# Plaque Segmentation Function

library(EBImage)
library(ggplot2)
library(tiff)

plaque_segmentation <- function(image_path, 
                                sigma = 2, 
                                threshold = 30,  # Threshold in 0-255 range (auto-converted if image is 0-1)
                                min_area = 400,
                                max_eccentricity = 0.7,  # Maximum eccentricity (0-1, where 0 is circle, 1 is line)
                                max_radius_cv = NULL,  # Optional: max (s.radius.sd / s.radius.mean), NULL to disable
                                max_radius_ratio = NULL,  # Optional: max (s.radius.max / s.radius.min), NULL to disable
                                show_plots = TRUE,
                                max_plot_size = 5000,
                                auto_downsample = TRUE) {
  
  # Step 1: Load the image
  start_time <- Sys.time()
  message(paste("[", format(start_time, "%H:%M:%S"), "] Loading image..."))
  if (is.character(image_path)) {
    img <- readImage(image_path)
  } else {
    img <- image_path
  }
  
  # Convert to grayscale if needed
  if (colorMode(img) == Color) {
    img <- channel(img, "gray")
  }
  
  # Check image size and warn if too large
  img_dims <- dim(img)
  img_width <- img_dims[1]
  img_height <- img_dims[2]
  total_pixels <- img_width * img_height
  
  load_time <- difftime(Sys.time(), start_time, units = "secs")
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Image loaded in", round(load_time, 2), "seconds"))
  message(paste("Image dimensions:", paste(img_dims, collapse = " x ")))
  message(paste("Total pixels:", format(total_pixels, scientific = FALSE, big.mark = ",")))
  
  # Warn if image is very large
  if (total_pixels > 500000000) {  # > 500M pixels
    message("WARNING: Image is very large (>500M pixels). Processing may take a long time.")
    message("Consider using show_plots=FALSE to skip plotting for faster processing.")
  }
  
  # Determine if we need to downsample for plotting
  needs_downsample <- (img_width > max_plot_size || img_height > max_plot_size)
  if (needs_downsample && show_plots) {
    if (auto_downsample) {
      # Calculate downsampling factor to fit within max_plot_size
      downsample_factor <- max(img_width, img_height) / max_plot_size
      message(paste("Image is too large for plotting. Will downsample by factor", 
                    round(downsample_factor, 2), "for display only."))
    } else {
      message("WARNING: Image is very large. Plotting may fail or be very slow.")
      message("Consider setting auto_downsample=TRUE or show_plots=FALSE")
    }
  }
  
  # Step 2: Apply Gaussian filter with σ = 2
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Applying Gaussian filter (σ =", sigma, ")..."))
  message("  This may take a while for large images...")
  filter_start <- Sys.time()
  img_filtered <- gblur(img, sigma = sigma)
  filter_time <- difftime(Sys.time(), filter_start, units = "secs")
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Gaussian filter completed in", round(filter_time, 2), "seconds"))
  
  # Diagnostic: Check image value range
  img_min <- min(img_filtered, na.rm = TRUE)
  img_max <- max(img_filtered, na.rm = TRUE)
  img_mean <- mean(img_filtered, na.rm = TRUE)
  img_median <- median(img_filtered, na.rm = TRUE)
  message(paste("Filtered image statistics:"))
  message(paste("  Min:", round(img_min, 2), "Max:", round(img_max, 2), 
                "Mean:", round(img_mean, 2), "Median:", round(img_median, 2)))
  
  # Step 3: Binarize with intensity threshold
  # Note: threshold parameter is always specified in 0-255 range for user convenience.
  # If image is normalized to 0-1, we automatically convert threshold/255.
  # This way users can always think in terms of 0-255 values regardless of image format.
  message(paste("Binarizing with threshold =", threshold, "(0-255 range)..."))
  # Check if image is normalized (0-1 range) with tolerance for floating point precision
  # If max value is <= 1.1, we assume it's normalized (allowing some margin for rounding)
  is_normalized <- (img_max <= 1.1)
  if (is_normalized) {
    # Image is normalized to 0-1: convert threshold from 0-255 to 0-1
    threshold_normalized <- threshold / 255
    message(paste("  Image is normalized (0-1, max =", round(img_max, 4), 
                  "), auto-converting threshold:", threshold, "->", round(threshold_normalized, 4)))
    img_binary <- img_filtered > threshold_normalized
  } else {
    # Image is in 0-255 range: use threshold directly
    message(paste("  Image range is 0-", round(img_max, 2), ", using threshold =", threshold))
    img_binary <- img_filtered > threshold
  }
  
  # Diagnostic: Check binary image
  foreground_pixels <- sum(img_binary, na.rm = TRUE)
  total_pixels_binary <- length(img_binary)
  foreground_percent <- (foreground_pixels / total_pixels_binary) * 100
  message(paste("  Binary image: foreground pixels =", format(foreground_pixels, big.mark = ","), 
                "(", round(foreground_percent, 2), "%)"))
  
  # Step 4: Segment plaques using bwlabel
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Segmenting plaques with bwlabel..."))
  message("  This may take a while for large images...")
  label_start <- Sys.time()
  plaque_labels <- bwlabel(img_binary)
  label_time <- difftime(Sys.time(), label_start, units = "secs")
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] bwlabel completed in", round(label_time, 2), "seconds"))
  
  # Diagnostic: Count connected components
  num_labels <- max(plaque_labels, na.rm = TRUE)
  message(paste("  Found", num_labels, "connected components"))
  
  # Validate image dimensions to prevent plotting errors
  img_dims <- dim(img)
  if (length(img_dims) < 2 || any(img_dims <= 0) || any(is.infinite(img_dims)) || any(is.na(img_dims))) {
    stop("Invalid image dimensions detected. Cannot proceed with plotting.")
  }
  message(paste("Image dimensions:", paste(img_dims, collapse = " x ")))
  
  # Step 5: Calculate size and center of each plaque
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Computing plaque features..."))
  features_start <- Sys.time()
  plaque_moments <- computeFeatures.moment(plaque_labels)
  plaque_shapes <- computeFeatures.shape(plaque_labels)
  features_time <- difftime(Sys.time(), features_start, units = "secs")
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Feature computation completed in", round(features_time, 2), "seconds"))
  
  # Combine features and convert to data.frame for safe column access
  plaque_features <- as.data.frame(cbind(plaque_moments, plaque_shapes),
                                   stringsAsFactors = FALSE)
  # Add an explicit plaque_id column
  row_ids <- rownames(plaque_features)
  if (is.null(row_ids)) {
    row_ids <- seq_len(nrow(plaque_features))
  }
  plaque_features$plaque_id <- as.integer(row_ids)
  
  # Diagnostic: Show plaque statistics before filtering
  if (nrow(plaque_features) > 0) {
    message(paste("  Total plaques detected (before filtering):", nrow(plaque_features)))
    if ("s.area" %in% colnames(plaque_features)) {
      message(paste("  Area statistics (pixels):"))
      message(paste("    Min:", round(min(plaque_features$s.area), 2),
                    "Max:", round(max(plaque_features$s.area), 2),
                    "Mean:", round(mean(plaque_features$s.area), 2),
                    "Median:", round(median(plaque_features$s.area), 2)))
    }
  } else {
    message("  WARNING: No plaques detected at all! Check threshold settings.")
  }
  
  # Step 6: Filter plaques with multiple criteria
  filter_criteria <- c(paste("area >", min_area))
  if (!is.null(max_eccentricity)) {
    filter_criteria <- c(filter_criteria, paste("eccentricity <", max_eccentricity))
  }
  if (!is.null(max_radius_cv)) {
    filter_criteria <- c(filter_criteria, paste("radius_cv <", max_radius_cv))
  }
  if (!is.null(max_radius_ratio)) {
    filter_criteria <- c(filter_criteria, paste("radius_ratio <", max_radius_ratio))
  }
  message(paste("Filtering plaques with:", paste(filter_criteria, collapse = ", ")))
  
  # Start with all plaques
  plaque_features_filtered <- plaque_features
  
  # Filter 1: Area
  plaque_features_filtered <- plaque_features_filtered[
    plaque_features_filtered$s.area > min_area, , drop = FALSE]
  message(paste("  After area filtering (>", min_area, "px):", nrow(plaque_features_filtered), "plaques"))
  
  # Filter 2: Eccentricity
  if (!is.null(max_eccentricity) && "m.eccentricity" %in% colnames(plaque_features_filtered)) {
    valid_ecc <- !is.na(plaque_features_filtered$m.eccentricity)
    before_count <- nrow(plaque_features_filtered)
    plaque_features_filtered <- plaque_features_filtered[
      valid_ecc & plaque_features_filtered$m.eccentricity < max_eccentricity, , drop = FALSE]
    
    # Diagnostic: Show eccentricity statistics
    if (before_count > 0) {
      num_filtered_by_ecc <- before_count - nrow(plaque_features_filtered)
      if (nrow(plaque_features_filtered) > 0) {
        ecc_stats <- plaque_features_filtered$m.eccentricity[!is.na(plaque_features_filtered$m.eccentricity)]
        if (length(ecc_stats) > 0) {
          message(paste("  Eccentricity statistics (after filtering):"))
          message(paste("    Min:", round(min(ecc_stats, na.rm = TRUE), 3),
                        "Max:", round(max(ecc_stats, na.rm = TRUE), 3),
                        "Mean:", round(mean(ecc_stats, na.rm = TRUE), 3)))
        }
      }
      if (num_filtered_by_ecc > 0) {
        message(paste("  Filtered out", num_filtered_by_ecc, "plaques with eccentricity >=", max_eccentricity))
      }
    }
  } else if (!is.null(max_eccentricity)) {
    message("  WARNING: m.eccentricity column not found. Skipping eccentricity filter.")
  }
  
  # Filter 3: Radius CV (coefficient of variation: sd/mean)
  if (!is.null(max_radius_cv) && 
      "s.radius.sd" %in% colnames(plaque_features_filtered) && 
      "s.radius.mean" %in% colnames(plaque_features_filtered)) {
    before_count <- nrow(plaque_features_filtered)
    # Calculate radius CV, avoiding division by zero
    valid_radius_mean <- !is.na(plaque_features_filtered$s.radius.mean) & 
                         plaque_features_filtered$s.radius.mean > 0
    valid_radius_sd <- !is.na(plaque_features_filtered$s.radius.sd)
    valid_radius_cv <- valid_radius_mean & valid_radius_sd
    
    if (sum(valid_radius_cv) > 0) {
      radius_cv <- plaque_features_filtered$s.radius.sd[valid_radius_cv] / 
                   plaque_features_filtered$s.radius.mean[valid_radius_cv]
      keep_mask <- rep(FALSE, nrow(plaque_features_filtered))
      keep_mask[valid_radius_cv] <- radius_cv < max_radius_cv
      plaque_features_filtered <- plaque_features_filtered[keep_mask, , drop = FALSE]
      
      message(paste("  After radius CV filtering (<", max_radius_cv, "):", 
                    nrow(plaque_features_filtered), "plaques"))
      num_filtered_by_cv <- before_count - nrow(plaque_features_filtered)
      if (num_filtered_by_cv > 0) {
        message(paste("  Filtered out", num_filtered_by_cv, "plaques with radius_cv >=", max_radius_cv))
      }
    } else {
      message("  WARNING: Cannot calculate radius CV (missing or invalid radius data). Skipping.")
    }
  }
  
  # Filter 4: Radius ratio (max/min)
  if (!is.null(max_radius_ratio) && 
      "s.radius.max" %in% colnames(plaque_features_filtered) && 
      "s.radius.min" %in% colnames(plaque_features_filtered)) {
    before_count <- nrow(plaque_features_filtered)
    # Calculate radius ratio, avoiding division by zero
    valid_radius_min <- !is.na(plaque_features_filtered$s.radius.min) & 
                        plaque_features_filtered$s.radius.min > 0
    valid_radius_max <- !is.na(plaque_features_filtered$s.radius.max)
    valid_radius_ratio <- valid_radius_min & valid_radius_max
    
    if (sum(valid_radius_ratio) > 0) {
      radius_ratio <- plaque_features_filtered$s.radius.max[valid_radius_ratio] / 
                      plaque_features_filtered$s.radius.min[valid_radius_ratio]
      keep_mask <- rep(FALSE, nrow(plaque_features_filtered))
      keep_mask[valid_radius_ratio] <- radius_ratio < max_radius_ratio
      plaque_features_filtered <- plaque_features_filtered[keep_mask, , drop = FALSE]
      
      message(paste("  After radius ratio filtering (<", max_radius_ratio, "):", 
                    nrow(plaque_features_filtered), "plaques"))
      num_filtered_by_ratio <- before_count - nrow(plaque_features_filtered)
      if (num_filtered_by_ratio > 0) {
        message(paste("  Filtered out", num_filtered_by_ratio, "plaques with radius_ratio >=", max_radius_ratio))
      }
    } else {
      message("  WARNING: Cannot calculate radius ratio (missing or invalid radius data). Skipping.")
    }
  }
  
  message(paste("Found", nrow(plaque_features_filtered), "plaques after all filtering"))
  
  # Create filtered labeled mask (only filtered plaques, with unique integer IDs starting from 1)
  # Also add label_id to plaque_features_filtered to match the mask labels
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Creating filtered labeled mask..."))
  if (nrow(plaque_features_filtered) > 0) {
    keep_ids <- plaque_features_filtered$plaque_id
    unique_old_ids <- sort(unique(keep_ids))
    
    # Create mapping from old_id to new_id (label_id) - optimized vectorized approach
    id_mapping <- data.frame(
      plaque_id = unique_old_ids,
      label_id = seq_along(unique_old_ids),
      stringsAsFactors = FALSE
    )
    
    # Use vectorized lookup table approach for much faster processing
    # Create lookup table: old_id -> new_id (0 for IDs not in keep_ids)
    max_old_id <- max(plaque_labels, na.rm = TRUE)
    lookup_table <- integer(max_old_id + 1)  # +1 because IDs start from 1
    lookup_table[unique_old_ids + 1] <- seq_along(unique_old_ids)  # +1 for 1-based indexing
    
    message("  Applying lookup table to remap labels (this may take a while for large images)...")
    # Get image data and apply lookup table directly (vectorized, much faster)
    labels_data <- imageData(plaque_labels)
    # Vectorized remapping using lookup table
    labels_data_remapped <- lookup_table[as.integer(labels_data) + 1]
    labels_data_remapped[is.na(labels_data_remapped)] <- 0
    
    # Create new Image object with remapped labels
    plaque_labels_filtered <- Image(labels_data_remapped, dim = dim(plaque_labels))
    # Add label_id column to plaque_features_filtered to match mask labels
    plaque_features_filtered <- merge(plaque_features_filtered, id_mapping, by = "plaque_id", all.x = TRUE)
    # Reorder columns to put label_id right after plaque_id
    other_cols <- colnames(plaque_features_filtered)[colnames(plaque_features_filtered) != "label_id" & colnames(plaque_features_filtered) != "plaque_id"]
    col_order <- c("plaque_id", "label_id", other_cols)
    plaque_features_filtered <- plaque_features_filtered[, col_order]
    message(paste("  Created filtered labeled mask with", length(unique_old_ids), "plaques (IDs: 1 to", length(unique_old_ids), ")"))
    message(paste("  Added label_id column to CSV to match mask labels (label_id corresponds to mask pixel values)"))
  } else {
    # No plaques after filtering, create empty mask
    plaque_labels_filtered <- plaque_labels * 0
    plaque_features_filtered$label_id <- integer(0)
    message("  No plaques after filtering, created empty labeled mask")
  }
  
  # Diagnostic: If no plaques after filtering, provide suggestions
  if (nrow(plaque_features_filtered) == 0 && nrow(plaque_features) > 0) {
    max_area <- max(plaque_features$s.area, na.rm = TRUE)
    message("  WARNING: All plaques were filtered out!")
    message(paste("  Largest plaque area:", round(max_area, 2), "pixels"))
    message(paste("  Current min_area threshold:", min_area, "pixels"))
    if ("m.eccentricity" %in% colnames(plaque_features)) {
      valid_ecc <- !is.na(plaque_features$m.eccentricity)
      if (sum(valid_ecc) > 0) {
        max_ecc <- max(plaque_features$m.eccentricity[valid_ecc], na.rm = TRUE)
        message(paste("  Maximum plaque eccentricity:", round(max_ecc, 3)))
        message(paste("  Current max_eccentricity threshold:", max_eccentricity))
      }
    }
    message("  Suggestions:")
    message("    1. Try reducing min_area parameter")
    message("    2. Try increasing max_eccentricity parameter")
  } else if (nrow(plaque_features) == 0) {
    message("  WARNING: No plaques detected. Possible issues:")
    message("    1. Threshold may be too high (try lower threshold)")
    message("    2. Threshold may be too low (try higher threshold)")
    message("    3. Image may need different preprocessing")
    message(paste("    Current threshold:", threshold, 
                  "| Image range:", round(img_min, 2), "-", round(img_max, 2)))
  }
  
  # Create visualization plots if requested
  if (show_plots) {
    message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Creating visualization plots..."))
    plot_start <- Sys.time()
    # Helper function to safely plot EBImage objects with automatic downsampling for large images
    safe_plot <- function(image_obj, main_title, downsample_for_display = NULL) {
      tryCatch({
        # Check if image object is valid
        if (is.null(image_obj)) {
          stop("Image object is NULL")
        }
        
        # Get dimensions and validate
        img_dims <- dim(image_obj)
        if (length(img_dims) < 2 || any(img_dims <= 0) || any(is.infinite(img_dims)) || any(is.na(img_dims))) {
          stop(paste("Invalid dimensions:", paste(img_dims, collapse = " x ")))
        }
        
        # Convert EBImage object to plain array/matrix for safer plotting
        if (is(image_obj, "Image")) {
          # Extract the image data as a plain array
          image_data <- imageData(image_obj)
          
          # Handle multi-dimensional images
          if (length(dim(image_data)) > 2) {
            # If multi-dimensional, take first frame
            image_data <- image_data[,,1]
          }
          
          # Ensure 2D matrix
          if (length(dim(image_data)) != 2) {
            stop("Image must be 2D for plotting")
          }
          
          # Downsample if needed and requested
          if (!is.null(downsample_for_display) && downsample_for_display > 1) {
            current_dims <- dim(image_data)
            new_width <- round(current_dims[1] / downsample_for_display)
            new_height <- round(current_dims[2] / downsample_for_display)
            # Use resize function from EBImage for downsampling
            image_data <- resize(image_data, w = new_width, h = new_height)
            # Update title to indicate downsampling
            main_title <- paste0(main_title, " (downsampled ", round(downsample_for_display, 1), "x)")
          }
          
          # For labeled images (integer values), convert to normalized float
          max_val <- max(image_data, na.rm = TRUE)
          if (is.finite(max_val) && max_val > 1) {
            image_data <- image_data / max_val
          }
          
          # Ensure values are in valid range [0, 1] for display
          image_data[image_data < 0] <- 0
          image_data[image_data > 1] <- 1
          image_data[is.na(image_data)] <- 0
          image_data[is.infinite(image_data)] <- 0
          
          # Convert to Image object for plotting (EBImage's plot method)
          image_obj <- Image(image_data)
        }
        
        # Plot with error handling
        plot(image_obj)
        title(main = main_title)
      }, error = function(e) {
        message(paste("Warning: Could not plot", main_title, "-", e$message))
        plot.new()
        title(main = paste(main_title, "(plot failed)"))
        mtext("Plotting error", col = "red")
      })
    }
    
    # Calculate downsampling factor if needed
    plot_downsample_factor <- NULL
    if (needs_downsample && auto_downsample) {
      plot_downsample_factor <- max(img_width, img_height) / max_plot_size
    }
    
    # Create a combined plot showing the segmentation pipeline
    # Increase margins to create more spacing between panels
    par(mfrow = c(2, 3), mar = c(5, 5, 5, 5), oma = c(2, 2, 2, 2))
    
    # Original image
    safe_plot(img, "Original Image", plot_downsample_factor)
    
    # Filtered image
    safe_plot(img_filtered, paste("Gaussian Filtered (sigma =", sigma, ")"), plot_downsample_factor)
    
    # Binary image
    safe_plot(img_binary, paste("Binary (threshold =", threshold, ")"), plot_downsample_factor)
    
    # Labeled image - safe_plot will handle normalization
    safe_plot(plaque_labels, "Labeled Plaques", plot_downsample_factor)
    
    # Filtered labeled image (only plaques > min_area)
    if (nrow(plaque_features_filtered) > 0) {
      tryCatch({
        # Create a new labeled image with only filtered plaques
        plaque_labels_filtered <- plaque_labels
        keep_ids <- plaque_features_filtered$plaque_id
        plaque_labels_filtered[!plaque_labels %in% keep_ids] <- 0
        safe_plot(plaque_labels_filtered, paste("Filtered Plaques (>", min_area, "px)"), plot_downsample_factor)
      }, error = function(e) {
        message(paste("Warning: Could not plot filtered labels -", e$message))
        plot.new()
        title(main = paste("Filtered Plaques (>", min_area, "px)"))
      })
    } else {
      safe_plot(img, paste("No plaques >", min_area, "px"), plot_downsample_factor)
    }
    
    # Centers overlay (replacing area distribution in the grid)
    safe_plot(img, "Filtered Plaque Centers", plot_downsample_factor)
    if (nrow(plaque_features_filtered) > 0) {
      tryCatch({
        cx <- plaque_features_filtered$m.cx
        cy <- plaque_features_filtered$m.cy
        
        # Adjust coordinates if image was downsampled for display
        if (!is.null(plot_downsample_factor) && plot_downsample_factor > 1) {
          cx <- cx / plot_downsample_factor
          cy <- cy / plot_downsample_factor
        }
        
        cross_half_len <- 6
        # Get displayed image dimensions (after potential downsampling)
        if (!is.null(plot_downsample_factor) && plot_downsample_factor > 1) {
          display_width <- round(img_width / plot_downsample_factor)
          display_height <- round(img_height / plot_downsample_factor)
        } else {
          display_width <- img_width
          display_height <- img_height
        }
        
        # Ensure crosses stay within image bounds
        x_left  <- pmax(1, cx - cross_half_len)
        x_right <- pmin(display_width, cx + cross_half_len)
        y_down  <- pmax(1, cy - cross_half_len)
        y_up    <- pmin(display_height, cy + cross_half_len)
        op_xpd <- par(xpd = NA)
        segments(x_left,  cy, x_right, cy, col = "red", lwd = 3)
        segments(cx, y_down, cx, y_up,   col = "red", lwd = 3)
        par(op_xpd)
      }, error = function(e) {
        message(paste("Warning: Could not draw centers -", e$message))
      })
    } else {
      mtext("No plaques after filtering", col = "red")
    }
    
    par(mfrow = c(1, 1))  # Reset plot layout

    # Extra page: Plaque area distribution (moved here)
    if (nrow(plaque_features) > 0) {
      # Make the histogram appear smaller by increasing margins and reducing cex
      op_hist <- par(mar = c(9, 9, 7, 7))
      hist(plaque_features$s.area,
           main = "Plaque Area Distribution",
           xlab = "Area (pixels)",
           col = "lightblue",
           breaks = 20,
           cex.main = 0.85,
           cex.lab = 0.85,
           cex.axis = 0.75)
      abline(v = min_area, col = "red", lwd = 2, lty = 2)
      legend("topright", paste("Min area =", min_area), col = "red", lty = 2, lwd = 2, cex = 0.85)
      par(op_hist)
    } else {
      plot.new(); title(main = "Plaque Area Distribution")
      mtext("No plaques detected", col = "red")
    }
    plot_time <- difftime(Sys.time(), plot_start, units = "secs")
    message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Plotting completed in", round(plot_time, 2), "seconds"))
  }
  
  # Return results
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Segmentation completed! Total time:", round(total_time, 2), "seconds"))
  results <- list(
    original_image = img,
    filtered_image = img_filtered,
    binary_image = img_binary,
    labeled_image = plaque_labels,
    labeled_image_filtered = plaque_labels_filtered,
    plaque_features = plaque_features,
    plaque_features_filtered = plaque_features_filtered,
    segmentation_params = list(
      sigma = sigma,
      threshold = threshold,
      min_area = min_area,
      max_eccentricity = max_eccentricity,
      max_radius_cv = max_radius_cv,
      max_radius_ratio = max_radius_ratio
    )
  )
  
  return(results)
}



# Helper function to test different thresholds and find optimal settings
test_thresholds <- function(image_path, 
                            sigma = 2,
                            thresholds = seq(10, 100, by = 10),
                            min_area = 400) {
  message("Testing different thresholds to find optimal settings...")
  message(paste(rep("=", 60), collapse = ""))
  
  # Load image
  if (is.character(image_path)) {
    img <- readImage(image_path)
  } else {
    img <- image_path
  }
  
  if (colorMode(img) == Color) {
    img <- channel(img, "gray")
  }
  
  # Apply Gaussian filter
  img_filtered <- gblur(img, sigma = sigma)
  img_max <- max(img_filtered, na.rm = TRUE)
  img_min <- min(img_filtered, na.rm = TRUE)
  img_mean <- mean(img_filtered, na.rm = TRUE)
  
  message(paste("Image statistics: Min =", round(img_min, 2), 
                "Max =", round(img_max, 2), 
                "Mean =", round(img_mean, 2)))
  message("")
  
  results_summary <- data.frame(
    threshold = numeric(),
    foreground_pct = numeric(),
    num_components = numeric(),
    num_plaques = numeric(),
    max_area = numeric(),
    mean_area = numeric()
  )
  
  for (thresh in thresholds) {
    # Binarize
    if (img_max <= 1) {
      img_binary <- img_filtered > (thresh / 255)
    } else {
      img_binary <- img_filtered > thresh
    }
    
    # Count foreground
    foreground_pct <- (sum(img_binary, na.rm = TRUE) / length(img_binary)) * 100
    
    # Segment
    plaque_labels <- bwlabel(img_binary)
    num_components <- max(plaque_labels, na.rm = TRUE)
    
    if (num_components > 0) {
      plaque_shapes <- computeFeatures.shape(plaque_labels)
      if (nrow(plaque_shapes) > 0) {
        num_plaques <- sum(plaque_shapes$s.area > min_area, na.rm = TRUE)
        max_area <- max(plaque_shapes$s.area, na.rm = TRUE)
        mean_area <- mean(plaque_shapes$s.area[plaque_shapes$s.area > min_area], na.rm = TRUE)
        if (is.na(mean_area)) mean_area <- 0
      } else {
        num_plaques <- 0
        max_area <- 0
        mean_area <- 0
      }
    } else {
      num_plaques <- 0
      max_area <- 0
      mean_area <- 0
    }
    
    results_summary <- rbind(results_summary, data.frame(
      threshold = thresh,
      foreground_pct = round(foreground_pct, 2),
      num_components = num_components,
      num_plaques = num_plaques,
      max_area = round(max_area, 2),
      mean_area = round(mean_area, 2)
    ))
  }
  
  # Print results
  message("Threshold Testing Results:")
  message("")
  print(results_summary)
  message("")
  message("Recommendation: Choose threshold that gives reasonable number of plaques")
  message("  (not too many small noise components, not too few)")
  
  return(results_summary)
}

# Function to save segmentation results
# save_images: if TRUE, save intermediate TIF images; if FALSE, only save CSV files
save_plaque_results <- function(results, output_dir = "plaque_segmentation_results", save_images = FALSE) {
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] Saving results to:", output_dir))
  save_start <- Sys.time()
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save images if requested
  if (save_images) {
    message("  Saving images...")
    writeImage(results$original_image, file.path(output_dir, "01_original.tif"))
    message("    - Original image saved")
    writeImage(results$filtered_image, file.path(output_dir, "02_filtered.tif"))
    message("    - Filtered image saved")
    writeImage(results$binary_image, file.path(output_dir, "03_binary.tif"))
    message("    - Binary image saved")
    writeImage(results$labeled_image, file.path(output_dir, "04_labeled.tif"))
    message("    - Labeled image saved")
    
    # Save an overlay image with red markers at plaque centers (filtered)
    if (!is.null(results$plaque_features_filtered) && nrow(results$plaque_features_filtered) > 0) {
      message("  Creating overlay image with plaque centers...")
      tryCatch({
        img_rgb <- toRGB(results$original_image)
        dims <- dim(img_rgb)
        img_w <- dims[1]
        img_h <- dims[2]
        
        # Check if image is too large for overlay (avoid memory issues)
        total_pixels <- img_w * img_h
        if (total_pixels > 1e9) {  # > 1 billion pixels
          message("    - WARNING: Image too large for overlay (>1B pixels). Skipping overlay image.")
          message("      Use 04_labeled.tif to visualize plaque locations.")
        } else {
          mark_radius <- 8
          centers_x <- round(results$plaque_features_filtered$m.cx)
          centers_y <- round(results$plaque_features_filtered$m.cy)
          
          # Draw red crosses at each center; clip within bounds to avoid incomplete lines
          for (i in seq_along(centers_x)) {
            cx <- centers_x[i]
            cy <- centers_y[i]
            if (is.na(cx) || is.na(cy)) next
            # Horizontal line
            x_range <- max(1, cx - mark_radius):min(img_w, cx + mark_radius)
            for (x in x_range) {
              img_rgb[x, cy, 1] <- 1
              img_rgb[x, cy, 2] <- 0
              img_rgb[x, cy, 3] <- 0
              # Thicken line by coloring neighboring pixels if within bounds
              if (cy - 1 >= 1) { img_rgb[x, cy - 1, 1] <- 1; img_rgb[x, cy - 1, 2] <- 0; img_rgb[x, cy - 1, 3] <- 0 }
              if (cy + 1 <= img_h) { img_rgb[x, cy + 1, 1] <- 1; img_rgb[x, cy + 1, 2] <- 0; img_rgb[x, cy + 1, 3] <- 0 }
            }
            # Vertical line
            y_range <- max(1, cy - mark_radius):min(img_h, cy + mark_radius)
            for (y in y_range) {
              img_rgb[cx, y, 1] <- 1
              img_rgb[cx, y, 2] <- 0
              img_rgb[cx, y, 3] <- 0
              # Thicken line by coloring neighboring pixels if within bounds
              if (cx - 1 >= 1) { img_rgb[cx - 1, y, 1] <- 1; img_rgb[cx - 1, y, 2] <- 0; img_rgb[cx - 1, y, 3] <- 0 }
              if (cx + 1 <= img_w) { img_rgb[cx + 1, y, 1] <- 1; img_rgb[cx + 1, y, 2] <- 0; img_rgb[cx + 1, y, 3] <- 0 }
            }
          }
          writeImage(img_rgb, file.path(output_dir, "06_centers_overlay.tif"))
          message("    - Overlay image saved (06_centers_overlay.tif)")
        }
      }, error = function(e) {
        message("    - WARNING: Could not save overlay image: ", e$message)
        message("      Use 04_labeled.tif to visualize plaque locations.")
      })
    }
  }
  
  # Save binary image of filtered plaques (always save, as it's a core output)
  if (!is.null(results$labeled_image_filtered)) {
    message("  Saving binary image of filtered plaques...")
    tryCatch({
      # Create binary mask: non-zero pixels (filtered plaques) = 1, zero = 0
      binary_filtered <- results$labeled_image_filtered > 0
      writeImage(binary_filtered, file.path(output_dir, "05_binary_filtered.tif"))
      message("    - Binary image of filtered plaques saved (05_binary_filtered.tif)")
    }, error = function(e) {
      message("    - WARNING: Could not save filtered binary image: ", e$message)
    })
  }
  
  # Save features (CSV files)
  message("  Saving CSV files...")
  write.csv(results$plaque_features, file.path(output_dir, "plaque_features_all.csv"), row.names = FALSE)
  message("    - plaque_features_all.csv saved")
  write.csv(results$plaque_features_filtered, file.path(output_dir, "plaque_features_filtered.csv"), row.names = FALSE)
  message("    - plaque_features_filtered.csv saved")
  
  # Save parameters
  # Convert NULL values to NA before converting to data.frame
  params <- results$segmentation_params
  params[sapply(params, is.null)] <- NA
  write.csv(data.frame(params, stringsAsFactors = FALSE), file.path(output_dir, "segmentation_parameters.csv"), row.names = FALSE)
  message("    - segmentation_parameters.csv saved")
  
  save_time <- difftime(Sys.time(), save_start, units = "secs")
  message(paste("[", format(Sys.time(), "%H:%M:%S"), "] All results saved in", round(save_time, 2), "seconds"))
  message(paste("Results saved to:", output_dir))
}
