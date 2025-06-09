#!/bin/bash
#SBATCH -p gpu_4l
#SBATCH -o zarr_%j.out        
#SBATCH -e zarr_%j.err        
#SBATCH -N 1 
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengg4c
#SBATCH -c 5

python mask2zarr_2_5_msk_yly8.py --cell_2d maskcell/maskcell_segmentation.tif --nuc_2d masknucleitorch/masknucleitorch_segmentation.tif --cell_3d 3dmask/cell.tif --nuc_3d 3dmask/nuclei.tif --molecular_file merged_spots/merged_goodPoints_max3d.csv --prefix cellline0101_C5Tg4h_final_6_tstgpu --json_cell maskcell/maskcell_geodataframe.json --json_nuc masknuclei/masknuclei_geodataframe.json --output_path .  --npy_dir .
