# README

## Overview
This project consists of analyses in STATES (Spatially-resolved transcriptomics and translatomics employed simultaneously).

## 📦 Requirements

All environment files are stored in the [`requirements`](./requirements) folder.

- **Python**: `STATES_analysis_PYenv.yml`  

  ```bash
  conda env create -f requirements/STATES_analysis_PYenv.yml
  conda activate STATES_analysis
  ```

- **R**:
  - STATES_cellline_sessionInfo.txt — cell line downstream analysis
  - STATES_tissue_sessionInfo.txt — tissue downstream analysis

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
The downstream analysis for tissue samples is located in the tissue_downstream folder. It includes normalization, cell type identification, TE calculation, spatial visualization, functional enrichment analysis, and subcellular translation heterogeneity analysis.

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
- Each `round` folder includes images for multiple channels, with a resolution of 3072x3072x42 pixels.

> Note: The demo scripts rely on additional MATLAB scripts located outside this folder:  
> - `tissue_upstream/core_matlab.m`  
> - `tissue_upstream/STATES-matlab-tissue/`

### Slurm Submission System Example

This example demonstrates the data decoding process using a Slurm system.

**Global Registration**: Use the `global.srp` file to globally align the ten imaging rounds.

**Local Registration, Spotfinding & Decoding**: Run the `dynamics_submit.sh` script, which calls scripts in the `local_position_subtile_tasks` directory to perform local registration, detect molecular spots, and decode the data.

**Stitching**: Use the `stitch.sh` script to stitch the local decoding results into a single, complete file.

### Output File

After the process is complete, the final decoding results will be saved in the `goodPoints_max3d.csv` file.

This file contains detailed information for each decoded spot, including:

- `x, y, z` coordinates: The position of the spot in 3D space.
- `gene`: The name of the decoded gene.

## **📑 Data Availability**

The complete STATES sequencing dataset will be made publicly available on Zenodo upon publication.  

👉 https://doi.org/10.5281/zenodo.15682594
