#!/bin/bash -l
# Generate stitchlinks for DAPI stitching 

user_dir="/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata"
sample="B4_WT_mousebrain2"

getOutpath() {
    echo "${user_dir}/${sample}/$1"
}

orderlist=$(getOutpath "04_stitch")"/orderlist"
segmethod="clustermap"
outpath=$(getOutpath "04_stitch")"/stitchlinks"
segoutpath=$(getOutpath "03_segmentation")"/"$segmethod


mkdir -p $outpath 

i=1
for j in `cat $orderlist`
do 
  j=`printf "%03d" $j`
  if [ "$j" = "0" ] || [ ! -s $segoutpath/Position$j/max_rotated_dapi.tif ] 
  then
    ln -sf ../blank.tif $outpath/tile_$i.tif
  else
    ln -sf ../../03_segmentation/$segmethod/Position$j/max_rotated_dapi.tif $outpath/tile_$i.tif 
  fi
  i=$((i+1))
done
