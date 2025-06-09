#!/bin/bash

mkdir -p stitch_tasks

for pos in $(seq -f "Position%03g" 1 49)
do
    script="stitch_tasks/stitch_task_${pos}.sh"
    cat << EOF > "$script"
#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o stitch_${pos}_%j.out
#SBATCH -e stitch_${pos}_%j.err
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 4

matlab -batch "core_matlab('C5Tg4h', 'stitch', '${pos}', 3072, 34, 1, 3, 10, '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata', '01_data', '02_registration', 'log', 'spotfinding_method', 'max3d', 'sqrt_pieces', 4)"
EOF
    chmod +x "$script"
done
echo "Generated 49 task scripts in 'stitch_tasks' directory."
