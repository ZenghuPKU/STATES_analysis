## ===============================
## 0. env
## ===============================
rm(list = ls())

library(Seurat)
library(Matrix)
library(dplyr)

setwd("/storage/lingyuan2/STATES_analysis/AD_downstream/12_SDEG")

states <- readRDS("states_with_plaque_info.rds")

rings_to_test <- c(
  "ring1_0_10um",
  "ring2_10_20um",
  "ring3_20_30um",
  "ring4_30_40um",
  "ring5_40_50um"
)

outer_ring <- "outer_50um_plus"

## ===============================
## main loop (only analyze 14mAD)
## ===============================

tp <- "14mAD"
obj <- subset(states, subset = type == tp)

message("Processing type: ", tp)

# get totalRNA assay counts and data
totalRNA_counts_mat <- GetAssayData(obj, assay = "totalRNA", slot = "counts")
totalRNA_data_mat   <- GetAssayData(obj, assay = "totalRNA", slot = "data")

totalRNA_counts_mat <- as(totalRNA_counts_mat, "dgCMatrix")
totalRNA_data_mat   <- as.matrix(totalRNA_data_mat) 

meta <- obj@meta.data

## outer cells
outer_cells <- rownames(meta)[meta$ring == outer_ring]

## ===============================
## each ring vs outer
## ===============================
for (rg in rings_to_test) {

  message("  Ring vs outer: ", rg)

  ring_cells <- rownames(meta)[meta$ring == rg]

  ## check the number of cells in the ring / outer
  if (length(ring_cells) < 10 | length(outer_cells) < 10) {
    message("  Skip ", rg, " (not enough cells overall)")
    next
  }

  ## --------------------------------
  ## 2.1 gene filtering (based on expression in the ring)
  ## --------------------------------
  total_ring_counts <- totalRNA_counts_mat[, ring_cells, drop = FALSE]

  valid_gene <- Matrix::rowSums(total_ring_counts >= 1) >= 10
  genes_use <- rownames(totalRNA_counts_mat)[valid_gene]

  message("    Genes tested (after ring filter): ", length(genes_use))

  ## --------------------------------
  ## 2.2 Wilcoxon test (gene by gene, using totalRNA data slot)
  ## --------------------------------
  res_list <- lapply(genes_use, function(g) {

    n_ring  <- length(ring_cells)
    n_outer <- length(outer_cells)

    if (n_ring < 10 | n_outer < 10) return(NULL)

    expr_ring  <- as.numeric(totalRNA_data_mat[g, ring_cells])
    expr_outer <- as.numeric(totalRNA_data_mat[g, outer_cells])

    wt <- wilcox.test(expr_ring, expr_outer, exact = FALSE)

    mean_ring <- mean(expr_ring)
    mean_outer <- mean(expr_outer)

    pseudo <- 1e-6
    log2FC <- log2((mean_ring + pseudo) / (mean_outer + pseudo))

    data.frame(
      gene       = g,
      ring       = rg,
      type       = tp,
      n_ring     = n_ring,
      n_outer    = n_outer,
      mean_ring  = mean_ring,
      mean_outer = mean_outer,
      delta_mean = mean_ring - mean_outer,
      log2FC     = log2FC,
      p_value    = wt$p.value
    )
  })

  res_df <- bind_rows(res_list)

  if (nrow(res_df) == 0) {
    message("    No genes left after per-gene filter.")
    next
  }

  ## --------------------------------
  ## 2.3 multiple correction
  ## --------------------------------
  res_df$padj <- p.adjust(res_df$p_value, method = "BH")

  ## --------------------------------
  ## 2.4 save results
  ## --------------------------------
  out_file <- paste0(
    "totalRNA_Wilcox_dataSlot_countsFilter_",
    tp, "_",
    rg, "_vs_outer.csv"
  )

  write.csv(res_df, out_file, row.names = FALSE)

  message("    Saved: ", out_file)
}

message("Finished type: ", tp)
