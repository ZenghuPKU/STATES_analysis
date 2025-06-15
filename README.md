# README

## Overview
This project consists of analyses in STATES (Spatially-resolved transcriptomics and translatomics employed simultaneously).

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
The downstream analysis for cell lines is organized into four main modules. It includes data preprocessing, stress trajectory modeling, TE (translation efficiency) and total RNA comparison, and subcellular localization analysis.

## Tissue Analysis

### Upstream
The upstream analysis for tissue samples is located in the `tissue_upstream` folder and follows a similar workflow as the cell line analysis.

### Downstream
The downstream analysis for tissue samples is located in the tissue_downstream folder. It includes normalization, cell type identification, TE calculation, spatial visualization, functional enrichment analysis, and subcellular translation heterogeneity analysis.

---
For more details, refer to the specific subdirectories and corresponding scripts in each folder.


