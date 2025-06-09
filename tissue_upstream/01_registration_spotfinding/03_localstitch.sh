#!/bin/bash

mkdir -p stitch_tasks

for pos in $(seq -f "Position%03g" 1 128)
do
    script="stitch_tasks/stitch_task_${pos}.sh"
    cat << EOF > "$script"
#!/bin/bash
#SBATCH -p cn-long
#SBATCH -o stitch_${pos}_%j.out
#SBATCH -e stitch_${pos}_%j.err
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengcnl
#SBATCH -c 4

matlab -batch "core_matlab('B4_WT_mousebrain2', 'stitch', '${pos}', 3072, 42, 1, 3, 11, '/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata', '01_data', '02_registration', 'log', 'spotfinding_method', 'max3d', 'sqrt_pieces', 4)"
EOF
    chmod +x "$script"
done
