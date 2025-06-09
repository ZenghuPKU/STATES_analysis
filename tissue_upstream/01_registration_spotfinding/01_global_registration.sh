#!/bin/bash

output_dir="global_position_tasks"
mkdir -p "$output_dir"

for i in $(seq -f "Position%03g" 1 128)
do
    script_name="${output_dir}/global_task_${i}.sh"

    cat <<EOF > $script_name
#!/bin/bash
#SBATCH -p fat4way               
#SBATCH -o global_${i}_%j.out    
#SBATCH -e global_${i}_%j.err    
#SBATCH -N 1                     
#SBATCH -c 4                     
#SBATCH -A huzeng_g1             
#SBATCH --qos=huzengf4w          


matlab -batch "core_matlab('B4_WT_mousebrain2', 'global_registration', '${i}', 3072, 42, 1, 3, 11, '/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata', '01_data', '02_registration', 'log', 'sqrt_pieces', 4)"
EOF

    echo "Generated $script_name"
done
