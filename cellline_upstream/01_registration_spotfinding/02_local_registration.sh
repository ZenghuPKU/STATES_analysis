#!/bin/bash

output_dir="local_position_subtile_tasks2"
mkdir -p "$output_dir"

base_command="core_matlab('C5Tg4h', 'local_registration', '{position}', 3072, 34, 1, 3, 10, '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata', '01_data', '02_registration', 'log', 'spotfinding_method', 'max3d', 'sqrt_pieces', 4, 'subtile', {subtile}, 'voxel_size', [1, 1, 1], 'end_bases', ['G','G','A','A','A'], 'barcode_mode', 'tri','split_loc', [5,11],'intensity_threshold', 0.2)"

for position in $(seq -f "Position%03g" 10 49)  
do
    for subtile in $(seq 1 16) 
    do
        script_file="${output_dir}/local_task_${position}_subtile${subtile}.sh"

        command="${base_command//\{position\}/$position}"
        command="${command//\{subtile\}/$subtile}"

        cat <<EOF > "$script_file"
#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o ${position}_subtile${subtile}_%j.out
#SBATCH -e ${position}_subtile${subtile}_%j.err
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 4

matlab -batch "$command"
EOF

        echo "Generated $script_file"
    done
done
