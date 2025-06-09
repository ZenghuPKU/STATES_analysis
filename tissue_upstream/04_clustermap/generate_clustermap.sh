#!/bin/bash
mkdir -p clustermap_tasks

for pos in $(seq -f "Position%03g" 1 126)
do
    script="clustermap_tasks/clustermap_task_${pos}.sh"
    cat << EOF > "$script"
#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o clustermap_${pos}_%j.out
#SBATCH -e clustermap_${pos}_%j.err
#SBATCH -N 1 
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 4

/appsnew/home/huzeng_pkuhpc/mambaforge/envs/starfinder/bin/python new_clustermaptest.py \\
         -IP ${pos} \\
         -IZ 42 -IEP F \\
         -IDir /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2 \\
         -Igenecsv /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/01_data/genes.csv \\
         -IDapi_path /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/01_data/round001/${pos} \\
         -IXY 3072 -IC 0.1 -ID 3 -IDR 1 -IPF 0.01 -ICR 35,10 \\
         -OP /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/03_segmentation/clustermap/${pos} \\
         -Igood_points_max3d /home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/02_registration/${pos}/goodPoints_max3d.csv -IDN 4
EOF
    chmod +x "$script"
done
echo "Generated 126 task scripts in 'clustermap_tasks' directory." 
