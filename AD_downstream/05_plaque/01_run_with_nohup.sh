#!/bin/bash
# Script to run plaque segmentation with nohup
# Usage: bash run_with_nohup.sh
#
# Configuration parameters (modify these as needed):
# ===================================================

# Input image file path (relative to WORK_DIR or absolute path)
IMAGE_PATH="14mADabeta.tif"

# Output directory name (will be created in WORK_DIR)
OUTPUT_DIR="14mADabeta"

# Segmentation parameters
SIGMA=3          # Gaussian filter sigma
THRESHOLD=100     # Intensity threshold (0-255 range)
MIN_AREA=800     # Minimum plaque area in pixels
MAX_ECCENTRICITY=0.95  # Maximum eccentricity (0-1, where 0 is circle, 1 is line)

# Optional radius-based filters (set to empty string "" to disable)
# MAX_RADIUS_CV: Maximum radius coefficient of variation (s.radius.sd / s.radius.mean)
#                Lower values mean more uniform radius distribution
MAX_RADIUS_CV=""  # Set to "" to disable this filter

# MAX_RADIUS_RATIO: Maximum radius ratio (s.radius.max / s.radius.min)
#                   Lower values mean more circular shapes
MAX_RADIUS_RATIO=""  # Set to "" to disable this filter

# Whether to save intermediate TIF images (true/false)
# Set to "true" to save: 01_original.tif, 02_filtered.tif, 03_binary.tif, 
#                        04_labeled.tif, 05_centers_overlay.tif
# Set to "false" to only save CSV files
SAVE_IMAGES="false"

# Working directory (where the R scripts are located)
WORK_DIR="/storage/lingyuan2/STATES_analysis/AD_downstream/05_plaque"

# Log file name
LOG_FILE="plaque_segmentation_14mAD.log"

# ===================================================
# End of configuration
# ===================================================

# Set working directory
cd "$WORK_DIR"

# Check if required R packages are available
echo "Checking R environment..."
Rscript -e "if (!require('EBImage', quietly=TRUE)) { stop('EBImage package not found. Please install it first.') } else { cat('EBImage: OK\n') }"
if [ $? -ne 0 ]; then
    echo "Error: Required R packages not found!"
    echo "Please install packages in RStudio or run: Rscript -e \"install.packages(c('EBImage', 'ggplot2'), repos='https://cloud.r-project.org')\""
    exit 1
fi

Rscript -e "if (!require('ggplot2', quietly=TRUE)) { stop('ggplot2 package not found.') } else { cat('ggplot2: OK\n') }"
if [ $? -ne 0 ]; then
    echo "Error: ggplot2 package not found!"
    exit 1
fi

echo "R environment check passed!"
echo ""

# Display configuration
echo "========================================"
echo "Configuration:"
echo "  Image file:        $IMAGE_PATH"
echo "  Output dir:       $OUTPUT_DIR"
echo "  Sigma:            $SIGMA"
echo "  Threshold:        $THRESHOLD"
echo "  Min area:         $MIN_AREA"
echo "  Max eccentricity: $MAX_ECCENTRICITY"
if [ -n "$MAX_RADIUS_CV" ] && [ "$MAX_RADIUS_CV" != "" ]; then
    echo "  Max radius CV:    $MAX_RADIUS_CV"
else
    echo "  Max radius CV:    (disabled)"
fi
if [ -n "$MAX_RADIUS_RATIO" ] && [ "$MAX_RADIUS_RATIO" != "" ]; then
    echo "  Max radius ratio: $MAX_RADIUS_RATIO"
else
    echo "  Max radius ratio: (disabled)"
fi
echo "  Save images:      $SAVE_IMAGES"
echo "  Working dir:      $WORK_DIR"
echo "  Log file:         $LOG_FILE"
echo "========================================"
echo ""

# Run R script with nohup and redirect output to log file
# Pass parameters as command line arguments
# Use empty string or "NULL" for optional parameters if not set
echo "Starting plaque segmentation in background..."
nohup Rscript run_plaque_segmentation.R \
    "$IMAGE_PATH" \
    "$OUTPUT_DIR" \
    "$SIGMA" \
    "$THRESHOLD" \
    "$MIN_AREA" \
    "$MAX_ECCENTRICITY" \
    "${MAX_RADIUS_CV:-NULL}" \
    "${MAX_RADIUS_RATIO:-NULL}" \
    "$SAVE_IMAGES" \
    "$WORK_DIR" \
    > "$LOG_FILE" 2>&1 &

# Print the process ID
PID=$!
echo "========================================"
echo "Job started in background with PID: $PID"
echo "Log file: $LOG_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "========================================"
echo ""
echo "Useful commands:"
echo "  Monitor progress:    tail -f $LOG_FILE"
echo "  Check if running:    ps -p $PID"
echo "  Stop the job:        kill $PID"
echo ""

