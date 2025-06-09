# README

## Overview
This project consists of analyses in STATES (Spatially-resolved transcriptomics and translatomics employed simultaneously)

## Cellline Analysis

### Upstream
The upstream analysis for cell lines is located in the `cellline_upstream` folder. It includes the following steps: registration, spotfinding, decoding, segmentation, stitching, assignment.

Key resources:
- The `genes.csv` file contains the target gene tri-barcode information.
- The core MATLAB script `core_matlab.m` and related MATLAB functions are organized in `STATE-matlab-cellline`.

#### Step-by-Step Process
- **01_registration_spotfinding**:
  - Perform global registration across all amplification rounds.
  - Split and conduct local registration.
  - Stitch back together.
  - Register the staining rounds accordingly.

- **02_IF_point_stitch**:
  - Stitch all FOVs for DAPI.
  - Apply the stitched transformation to other stainings and amplification rounds.

- **03_IF_mask**:
  - Use Cellpose to perform 2D segmentation on DAPI and Flamingo staining.
  - Obtain a 3D mask by element-wise multiplication of the mask and the original image.

- **04_assign**:
  - Assign amplification spots to the 3D mask.
  - Generate the final output.

### Downstream 
### Downstream  
The downstream analysis for spatial samples is organized into four main folders. It includes data preprocessing, stress trajectory modeling, TE (translation efficiency) and total RNA comparison, and subcellular localization analysis.

#### Step-by-Step Process  
- **01_preprocessing**:  
  - `01_combineh5ad.py`: Merge multiple `.zarr` files into a single integrated AnnData object.  
  - `02_filter.ipynb`: Filter out low-quality cells and genes based on standard QC metrics.

- **02_stress**:  
  - `01_meanTEcell.R`: Calculate mean TE per cell in each stress time.  
  - `02_cellcycle.ipynb`: Assign cell cycle phase labels using canonical marker genes.  
  - `03_wnnuamp_pseudotime.R`: Perform joint profiling of transcription and translation to generate a UMAP embedding, and infer pseudotime trajectories based on the integrated space.

- **03_TE_totalRNA**:  
  - `01_wilcoxon_totalRNA.ipynb`: Perform differential expression analysis on total RNA.  
  - `02_wilcoxon_TE.ipynb`: Identify genes with significant TE differences.  
  - `03_geneclustering.R`: Cluster genes based on TE and totalRNA.  
  - `04_overlap.R`: Draw a Venn diagram to visualize the overlap between two sets of differentially expressed genes, and perform a hypergeometric test to assess the statistical significance of the intersection.  
  - `05_geneinpseudotime.ipynb`: Visualize gene-level TE and totalRNA expression along pseudotime.  
  - `06_GO.R`: Conduct Gene Ontology enrichment for significant gene clusters.

- **04_subcellular**:  
  - `01_2dassign.py`: Assign transcripts to nucleus or cytoplasma compartments in 2D space.  
  - `02_2d_DR_zarr.py`: Calculate distance ratio for each point.  
  - `03_3bincellanalysis.ipynb`: Divide each cell into inner, middle, peripheral bins.  
  - `04_violinplot.R`: Generate violin plots comparing TE across subcellular compartments.  
  - `05_3bingeneanalysis.ipynb`: Analyze TE in each bins for individual genes.  
  - `06_diff_draw.R`: Plot subcellular difference of TE for individual genes.  
  - `07_genecurve.R`: Plot selected genes' TE across compartments.  
  - `08_DRspatialplot.ipynb`: Visualize spatial pattern of specific gene groups.  
  - `09_corr_distribution.ipynb`: Compute correlation distributions of TE/DR across genes.  
  - `10_geneinpseudotimeTEDR.ipynb`: Project gene-level TE and DR changes along pseudotime.

## Tissue Analysis

### Upstream
The upstream analysis for tissue samples is located in the `tissue_upstream` folder and follows a similar workflow as the cell line analysis.

### Downstream
The downstream analysis for tissue samples is located in the tissue_downstream folder. It includes normalization, cell type identification, TE (translation efficiency) calculation, spatial visualization, and functional enrichment analysis.

#### Step-by-Step Process
- **01_preprocessing**:
  - Load two raw `.h5ad` files and the `outer_counts` table from `remain_readsouter_extracted.csv`.  
  - Distinguish total RNA and ribosomal RNA reads, and computes the TE layer.

- **02_preprocessing_h5ad**:
  - Load and merge two replicates, add replicate identifiers and syncing any auxiliary read-info tables.  
  - Perform QC filtering on cells and genes, plot a replicate-to-replicate correlation scatter and max-count histograms, and filter low-expressed genes.  
  - Normalize, log-transform, scale, and regress out total counts; create cytoplasmic and ribosomal count layers; run PCA and Harmony batch correction; generate all diagnostic plots; and write out the final integrated `.h5ad` file.

- **03_Normalization_Visualization**:
  - Perform normalization of the integrated data.
  - Generate visualizations of global expression distributions after normalization.

- **04_Mixfind**:
  - Prior to cell type annotation, cells with abnormal distributions in the UMAP embedding—potentially representing mixtures of multiple cell types—were identified and excluded.
 
- **05_celltypes_identification**:
  - `01_celltypes_identification.R`: Annotate major cell types based on marker gene expression and clustering.
  - `02_Annotated_celltypes_UMAP_visualization.R`: Output UMAP coordinates with overlaid cell type labels. Additionally, generate separate UMAP plots for each annotated cell type label.
  - `03_Annotated_celltypes_spatial_visualization.R`: Map cell type annotations back to spatial coordinates.
  - `04_Marker_gene_dot_plot.R`: Visualize the expression of key marker genes in total RNA, rbRNA, and TE levels using dot plots.
  
- **06_TE_Distribution**:
  - `01_TE_Distribution_label2.R`: Calculate and visualize TE (translation efficiency) distributions for major cell types.
  - `02_TE_Distribution_label3.R`: Calculate and visualize TE distributions for sub cell types.
  - `03_TE_Spatial_map.R`: Map per-cell TE values to their original spatial locations.

- **07_Gene_clustering_based_on_TE_and_GO**:
  - Cluster genes into 10 gene modules based on their TE profiles.
  - Perform Gene Ontology (GO) enrichment analysis on each module.

- **08_TE_logFC_major_celltypes**:
  - Perform differential TE analysis between major cell types.

- **09_TE_logFC_HIP_vs_Cortex_and_GO**:
  - Conduct TE differential analysis specifically between hippocampal and cortical neuron subtypes.
  - Perform GO enrichment analysis on differentially expressed genes.

---
For more details, refer to the specific subdirectories and corresponding scripts in each folder.


