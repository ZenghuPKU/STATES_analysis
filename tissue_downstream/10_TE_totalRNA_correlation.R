# 10_TE_totalRNA_correlation
# Load libraries and set environment
library(Matrix)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(readr)
library(tidyr)
library(anndata)
library(reticulate)
rm(list = ls())
setwd("~/tissue_downstream/")
use_condaenv("scanpy_env", required = TRUE)
py_config()

point_col <- "#2166AC"

# Load annotated Seurat object and corresponding h5ad file
load("states_celltypes_identification.RData")
ad <- read_h5ad("mousebrain_harmony.h5ad")

# Extract raw count matrices (gene × cell)
totalRNA_raw_mat <- t(ad$layers['totalRNA_raw'])       
rbRNA_raw_mat    <- t(ad$layers['rbRNA'])   
meta_ad          <- ad$obs

# Convert into Seurat objects for consistency
s_total <- CreateSeuratObject(counts = totalRNA_raw_mat, meta.data = meta_ad)
s_ribo  <- CreateSeuratObject(counts = rbRNA_raw_mat,    meta.data = meta_ad)

# Extract normalized totalRNA matrix from Seurat object
norm_total_mat <- as.matrix(GetAssayData(states, assay = "totalRNA", slot = "data"))

# Ensure consistent set of cells across modalities and metadata
keep_cells <- Reduce(intersect, list(colnames(s_total), colnames(s_ribo),
                                     colnames(norm_total_mat), rownames(states@meta.data)))
stopifnot(length(keep_cells) > 0)

total_raw <- as.matrix(GetAssayData(s_total, slot = "counts"))[, keep_cells, drop = FALSE]
ribo_raw  <- as.matrix(GetAssayData(s_ribo,  slot = "counts"))[, keep_cells, drop = FALSE]
norm_tot  <- norm_total_mat[, keep_cells, drop = FALSE]
meta_df   <- states@meta.data[keep_cells, , drop = FALSE]

# Harmonize the set of genes across matrices
common_genes <- Reduce(intersect, list(rownames(total_raw), rownames(ribo_raw), rownames(norm_tot)))
stopifnot(length(common_genes) > 0)
total_raw <- total_raw[common_genes, , drop = FALSE]
ribo_raw  <- ribo_raw[common_genes,  , drop = FALSE]
norm_tot  <- norm_tot[common_genes,  , drop = FALSE]

# Define color schemes, cell type order, and marker sets
label3celltype_colors <- c(
  "AC1"="#ccba33","AC2"="#ffbe85","AC3"="#e3782b","CHOR"="#7f52a9","EPEN"="#c4b0d4",
  "CHO/PEP"="#97f4f7","INH_Sst"="#96abeb","INH_Pvalb"="#96cad4","INH_Cnr1_Vip"="#a8e1eb",
  "MLG"="#8597c6","OPC"="#667872","OLG1"="#e4f768","OLG2"="#e6db17","VLMC"="#1f76b3",
  "VSMC"="#00aeef","Peri/VEC"="#d3a59c","DE/MEN"="#b274e8","MSN"="#D96DA1","DGGRC"="#a6e8a6",
  "TEGLU CA1"="#77ed8f","TEGLU CA2"="#82ad2d","TEGLU CA3"="#28330b","TEGLU L2/3"="#cbfc60",
  "TEGLU L4"="#96db00","TEGLU L5"="#04b361","TEGLU L5/6"="#40d102","TEGLU L6"="#32a630",
  "TEGLU L6b"="#406e27","TEGLU Mix"="#c5fcc5"
)
label2celltype_colors <- c(
  "OPC"="#667872","DE/MEN"="#b274e8","CHO/PEP"="#97f4f7","OLG"="#e6db17","CHOR/EPEN"="#F1B2E9",
  "TEPN"="#6C9F35","MLG"="#8597c6","VAS"="#17BED0","INH"="#3182BD","AC"="#FDC06F"
)

order_l3 <- c("TEGLU CA1","TEGLU CA2","TEGLU CA3","TEGLU L2/3","TEGLU L4","TEGLU L5",
              "TEGLU L5/6","TEGLU L6","TEGLU L6b","DGGRC","MSN","INH_Pvalb","INH_Sst",
              "INH_Cnr1_Vip","DE/MEN","CHO/PEP","AC1","AC2","AC3","OLG1","OLG2",
              "OPC","MLG","CHOR","EPEN","Peri/VEC","VLMC","VSMC")
order_l2 <- c("TEPN","INH","DE/MEN","CHO/PEP","AC","OLG","OPC","MLG","CHOR/EPEN","VAS")

markers_l3 <- c("Myh11","Myl9","Ptgds","Gjb2","Rgs5","Ly6c1","Vim","Rarres2","Ttr","Enpp2",
                "Hexb","C1qa","Pdgfra","Cacng4","Mbp","Mobp","Mag","Mog","Gfap","Mt2",
                "Apoe","Clu","Slc7a10","Mfge8","Gng8","Cadps2","Prkcd","Ntng1","Cnr1","Vip",
                "Sst","Gad2","Gad1","Pvalb","Ppp1r1b","Adora2a","C1ql2","Cplx3","Rprm",
                "Sncb","Cplx1","Fezf2","Pde1a","Rgs4","Lamp5","Cux2","Chgb","Rgs14","Ppp3r1","Wfs1")
markers_l2 <- c("Flt1","Ly6c1","Ly6e","Ptgds","Rgs5","Rarres2","Clic6","Folr1","Enpp2","Ttr",
                "Ctss","C1qc","Csf1r","C1qa","Hexb","Cacng4","Olig2","Sox10","Pdgfra","Ptprz1",
                "Trf","Cnp","Fth1","Mobp","Plp1","Clu","Mt2","Atp1b2","Gja1","Aldoc",
                "Resp18","Hap1","Cadps2","Gng8","Tac2","Tubb5","Ntng1","Prkcd","Sparc","Sparcl1",
                "Npy","Sst","Pvalb","Gad2","Gad1","Hpca","Nell2","Mapk1","Ppp3r1","Nrgn")

# Construct long-format dataset
build_long_norel <- function(label_col, genes_keep = NULL,
                             exclude_levels = c("Mix"), min_cells_per_ct = 1) {
  stopifnot(label_col %in% colnames(meta_df))
  labs <- as.character(meta_df[[label_col]])
  names(labs) <- colnames(norm_tot)
  
  ct_tab <- sort(table(labs), decreasing = TRUE)
  ct_use <- setdiff(names(ct_tab), exclude_levels)
  ct_use <- ct_use[ct_tab[ct_use] >= min_cells_per_ct]
  if (!length(ct_use)) stop("No CT left after excluding: ", paste(exclude_levels, collapse=", "))
  
  genes <- if (is.null(genes_keep)) rownames(norm_tot) else intersect(genes_keep, rownames(norm_tot))
  if (!length(genes)) stop("No genes to keep in matrices.")
  
  out <- lapply(ct_use, function(ct) {
    cells <- names(labs)[labs == ct]
    if (!length(cells)) return(NULL)
    
    mean_norm_ct <- rowMeans(norm_tot[genes, cells, drop = FALSE], na.rm = TRUE)
    
    sum_tot_ct <- rowSums(total_raw[genes, cells, drop = FALSE], na.rm = TRUE)
    sum_rib_ct <- rowSums(ribo_raw[genes,  cells, drop = FALSE], na.rm = TRUE)
    TE_ct      <- sum_rib_ct / sum_tot_ct
    TE_ct[!is.finite(TE_ct)] <- NA_real_
    
    data.frame(
      Gene = genes,
      CellType = ct,
      mean_norm_totalRNA = as.numeric(mean_norm_ct),
      TE_raw = as.numeric(TE_ct),
      stringsAsFactors = FALSE
    )
  })
  
  L <- dplyr::bind_rows(out)
  L <- L %>% dplyr::filter(is.finite(mean_norm_totalRNA), is.finite(TE_raw))
  L
}

# All target genes correlation across cell types
per_gene_corr_linear <- function(L) {
  if (!nrow(L)) {
    return(data.frame(Gene = character(0), N_Celltypes = integer(0),
                      Pearson_r = numeric(0), Pvalue = numeric(0), stringsAsFactors = FALSE))
  }
  genes <- unique(L$Gene)
  res <- lapply(genes, function(g) {
    df <- L[L$Gene == g, , drop = FALSE]
    x <- df$mean_norm_totalRNA
    y <- df$TE_raw
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    if (n_ok >= 3) {
      ct <- tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = "pearson")),
                     error = function(e) NULL)
      r <- if (is.null(ct)) NA_real_ else as.numeric(unname(ct$estimate))
      p <- if (is.null(ct)) NA_real_ else as.numeric(ct$p.value)
    } else { r <- NA_real_; p <- NA_real_ }
    data.frame(Gene = g, N_Celltypes = n_ok, Pearson_r = r, Pvalue = p, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

# Marker-gene scatter plots across cell types
plot_marker_scatter_linear <- function(L, label_tag, order_vec, color_vec,
                                       point_size = 2.2, label_size = 2.6,
                                       min_ct_points = 3) {
  out_dir <- file.path(paste0("out_", label_tag)); if (!dir.exists(out_dir)) dir.create(out_dir, FALSE)
  L$CellType <- factor(L$CellType, levels = intersect(order_vec, unique(L$CellType)))
  colors_use <- color_vec[levels(L$CellType)]
  
  pdf_file <- file.path(out_dir, paste0(label_tag, "_MARKERS_scatter_totalNorm_vs_TE_BY_GENE.pdf"))
  pdf(pdf_file, width = 6.8, height = 6.2)
  for (g in unique(L$Gene)) {
    df <- L %>% dplyr::filter(Gene == g) %>% dplyr::filter(is.finite(mean_norm_totalRNA), is.finite(TE_raw))
    r_val <- NA_real_; p_val <- NA_real_
    if (nrow(df) >= min_ct_points) {
      ct <- suppressWarnings(cor.test(df$mean_norm_totalRNA, df$TE_raw, method = "pearson"))
      r_val <- unname(ct$estimate); p_val <- ct$p.value
    }
    p <- ggplot(df, aes(x = mean_norm_totalRNA, y = TE_raw, color = CellType, label = CellType)) +
      geom_point(size = point_size, alpha = 0.9) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8) +
      ggrepel::geom_text_repel(size = label_size, max.overlaps = Inf, segment.alpha = 0.5) +
      scale_color_manual(values = colors_use, drop = FALSE) +
      labs(
        title = paste0(label_tag, " | ", g, "  (per gene across cell types)"),
        subtitle = paste0("Pearson r = ", ifelse(is.finite(r_val), sprintf('%.3f', r_val), "NA"),
                          if (is.finite(p_val)) paste0(" | p = ", signif(p_val, 3)) else "",
                          " | n CT = ", nrow(df)),
        x = "totalRNA",
        y = "TE"
      ) +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
    print(p)
  }
  dev.off()
  message("Saved marker scatters: ", pdf_file)
}

# Distributions of correlation coefficients
plot_marker_r_density <- function(per_gene_df, label_tag, markers_vec, fill = "#3182BD") {
  out_dir <- file.path(paste0("out_", label_tag)); if (!dir.exists(out_dir)) dir.create(out_dir, FALSE)
  df <- per_gene_df %>% dplyr::filter(is.finite(Pearson_r), Gene %in% markers_vec)
  if (!nrow(df)) { message("No markers with finite r for ", label_tag); return(invisible(NULL)) }
  df <- df %>% dplyr::arrange(Pearson_r)
  readr::write_csv(df, file.path(out_dir, paste0(label_tag, "_MARKERS_per_gene_correlation.csv")))
  med_r <- median(df$Pearson_r)
  p <- ggplot(df, aes(x = Pearson_r)) +
    geom_density(fill = fill, alpha = 0.35, color = fill, linewidth = 1) +
    geom_rug(aes(x = Pearson_r), sides = "b", color = fill, alpha = 0.85, size = 0.55) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = med_r, linetype = "solid") +
    labs(
      title = paste0(label_tag, " | Marker r density"),
      subtitle = paste0("n markers = ", nrow(df), "  |  median r = ", sprintf('%.3f', med_r)),
      x = "Pearson r",
      y = "Density"
    ) + theme_minimal(base_size = 12)
  ggsave(filename = file.path(out_dir, paste0(label_tag, "_MARKERS_r_density.pdf")),
         plot = p, width = 6.8, height = 4.8, dpi = 300)
  ggsave(filename = file.path(out_dir, paste0(label_tag, "_MARKERS_r_density.png")),
         plot = p, width = 6.8, height = 4.8, dpi = 300)
}

plot_allgenes_r_density <- function(per_gene_df, label_tag, fill = "#2166AC", add_rug = FALSE, rug_max_n = 1000) {
  out_dir <- file.path(paste0("out_", label_tag)); if (!dir.exists(out_dir)) dir.create(out_dir, FALSE)
  df <- per_gene_df[is.finite(per_gene_df$Pearson_r), , drop = FALSE]
  if (!nrow(df)) { message("No finite r for all genes under ", label_tag); return(invisible(NULL)) }
  n        <- nrow(df)
  median_r <- median(df$Pearson_r)
  q25      <- as.numeric(quantile(df$Pearson_r, 0.25, names = FALSE))
  q75      <- as.numeric(quantile(df$Pearson_r, 0.75, names = FALSE))
  readr::write_csv(df, file.path(out_dir, paste0(label_tag, "_ALLGENES_per_gene_correlation.csv")))
  p <- ggplot(df, aes(x = Pearson_r)) +
    geom_density(fill = fill, alpha = 0.30, color = fill, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = median_r, linetype = "solid") +
    labs(
      title = paste0(label_tag, " | All target genes r density"),
      subtitle = paste0("n = ", n, "  |  median r = ", sprintf('%.3f', median_r),
                        "  |  IQR [", sprintf('%.3f', q25), ", ", sprintf('%.3f', q75), "]"),
      x = "Pearson r",
      y = "Density"
    ) + theme_minimal(base_size = 12)
  if (isTRUE(add_rug)) {
    set.seed(1L)
    idx <- if (nrow(df) > rug_max_n) sample.int(nrow(df), rug_max_n) else seq_len(nrow(df))
    p <- p + geom_rug(data = df[idx, , drop = FALSE], aes(x = Pearson_r), sides = "b",
                      color = fill, alpha = 0.5, size = 0.4)
  }
  ggsave(filename = file.path(out_dir, paste0(label_tag, "_ALLGENES_r_density.pdf")),
         plot = p, width = 6.8, height = 4.8, dpi = 300)
  ggsave(filename = file.path(out_dir, paste0(label_tag, "_ALLGENES_r_density.png")),
         plot = p, width = 6.8, height = 4.8, dpi = 300)
}

# Unified pipeline runner (label2 / label3)
run_norel_pipeline <- function(label_col, label_tag, order_vec, color_vec, markers_vec,
                               exclude_levels = c("Mix")) {
  message("== Running NOREL pipeline: ", label_tag, " ==")
  L_all <- build_long_norel(label_col, genes_keep = NULL, exclude_levels = exclude_levels)
  per_gene_all <- per_gene_corr_linear(L_all)
  plot_allgenes_r_density(per_gene_all, label_tag, fill = "#2166AC", add_rug = FALSE)
  mk <- intersect(markers_vec, unique(L_all$Gene))
  if (length(mk)) {
    L_mk <- L_all %>% dplyr::filter(Gene %in% mk)
    per_gene_mk <- per_gene_all %>% dplyr::filter(Gene %in% mk)
    plot_marker_scatter_linear(L_mk, label_tag, order_vec, color_vec)
    plot_marker_r_density(per_gene_mk, label_tag, mk, fill = "#3182BD")
  } else {
    message("No markers found in dataset for ", label_tag)
  }
  invisible(list(L_all = L_all, per_gene_all = per_gene_all))
}

# Execution of analyses
run_norel_pipeline("states_nn_alg1_label2", "label2", order_l2, label2celltype_colors, markers_l2)
run_norel_pipeline("states_nn_alg1_label3", "label3", order_l3, label3celltype_colors, markers_l3)

