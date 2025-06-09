#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o point_%j.out
#SBATCH -e point_%j.err
#SBATCH -N 1 
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 4



/appsnew/home/huzeng_pkuhpc/mambaforge/envs/starfinder/bin/python new_stitch0220_v13.py \
        -IXY 3072 \
        -ID /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2 \
        -IO /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/orderlist \
        -IS /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/03_segmentation/clustermap \
        -OD /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch \
        -ITR  /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/stitchlinks/TileConfiguration.registered.txt \
        -OC /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/cell_center.csv \
        -OR /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/remain_reads.csv\
        -IT /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/stitchlinks/TileConfiguration.txt \
        -SM clustermap
