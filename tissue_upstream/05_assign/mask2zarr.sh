#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o zarr_%j.out
#SBATCH -e zarr_%j.err
#SBATCH -N 1 
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 5

python mask2zarr.py -c /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/06_mask_v13/clustermap_mask.json -n /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/masknucleitorch/masknucleitorch_geodataframe.json -nt /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/masknucleitorch/masknucleitorch_segmentation.tif -p B4_WT_mousebrain2 -o 07_zarr -t 5 -m /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/remain_readsv13.csv

