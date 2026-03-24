# README

## Overview
This project consists of analyses in STATES (Spatially-resolved transcriptomics and translatomics employed simultaneously).

## 📦 Requirements

All environment files are stored in the [`requirements`](./requirements) folder.

- **MATLAB**: R2023b was used for upstream image processing (no additional installation required beyond MATLAB itself).

- **Python**: `STATES_analysis_PYenv.yml`  
  Typical installation time: ~20–40 min (depends on network and hardware). 
  ```bash
  conda env create -f requirements/STATES_analysis_PYenv.yml
  conda activate STATES_analysis
  ```
- **R**:  
  STATES_cellline_sessionInfo.txt — cell line downstream analysis  
  STATES_tissue_sessionInfo.txt — tissue downstream analysis
  STATES_AD_sessionInfo.txt - AD downstream analysis  
  Typical package installation time: ~20–40 min (longer if Bioconductor dependencies are compiled).

Most computations were performed on a local server (Intel Xeon Platinum 8368 CPU, 152 cores, 251 GiB RAM) and on the High-Performance Computing Platform of the Center for Life Sciences (Peking University).

## **🔬** Cellline Analysis

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

## 🧠 Tissue Analysis

### Upstream
The upstream analysis for tissue samples is located in the `tissue_upstream` folder and follows a similar workflow as the cell line analysis.

### Downstream
The downstream analysis for tissue samples is located in the `tissue_downstream` folder. It includes normalization, cell type identification, TE calculation, spatial visualization, functional enrichment analysis, and subcellular translation heterogeneity analysis.

## 🧬 AD Analysis

### Upstream
The upstream analysis for the 5xFAD samples followed a workflow similar to that used in the `tissue_upstream` folder.

### Downstream
The downstream analysis for the 5xFAD samples is located in the `AD_downstream` folder. It includes normalization, cell type identification, TE calculation, functional enrichment analysis, trajectory analysis, and spatial analysis.

---
For more details, refer to the specific subdirectories and corresponding scripts in each folder.

## **📂** Demo

The `demo` folder provides example scripts illustrating how to process multi-channel, multi-round imaging data into gene identities with corresponding 3D coordinates.  

### Sample Data

The sample data can be downloaded directly from Google Drive:  

👉 [Download link](https://drive.google.com/drive/folders/1LJ7fk0XnUjWcBQzZzLkGZOmIRMB6YeYq)

The folder contains following files and folders:

- `genes.csv`: The gene barcode file used for decoding.
- `round001/` to `round011/`: Contain the image data for eleven imaging rounds.
- Each `round` folder includes images for multiple channels.

> Note: The demo scripts rely on additional MATLAB scripts located outside this folder:  
> - `tissue_upstream/core_matlab.m`  
> - `tissue_upstream/STATES-matlab-tissue/`

### Slurm Submission System Example

This demo was benchmarked on a high-memory HPC node:
- Resources: 4 CPU cores, ~168 GB RAM
- Runtime: ~15 minutes for 11 imaging rounds (3072 × 3072 × 42 pixels each)

On a normal desktop computer, memory may not be sufficient to process the full dataset directly. In such cases, we recommend block-wise processing or downsampling of the image stacks, which will increase the runtime accordingly.

**Global Registration**:
  ```bash
  # globally align multi-rounds
  sbatch global.srp
  ```
**Local Registration, Spotfinding & Decoding**: 
  ```bash
  # calls scripts in the local_position_subtile_tasks directory
  sbatch dynamics_submit.sh
  ```
**Stitching**:
  ```bash  
  # stitch the local decoding results into a single, complete file 
  sbatch stitch.sh
  ```

### Output File

After the process is complete, the final decoding results will be saved in the `goodPoints_max3d.csv` file.

This file contains detailed information for each decoded spot, including:

- `x, y, z` coordinates: The position of the spot in 3D space.
- `gene`: The name of the decoded gene.

## **📑 Data Availability**

The complete STATES sequencing dataset will be made publicly available on Zenodo upon publication.  

👉 https://doi.org/10.5281/zenodo.15682594
