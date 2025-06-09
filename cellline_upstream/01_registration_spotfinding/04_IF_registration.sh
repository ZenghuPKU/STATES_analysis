#!/bin/bash

output_dir="protein_position_tasks"
mkdir -p "$output_dir"

for i in $(seq -f "Position%03g" 1 49)
do
    script_name="${output_dir}/protein_task_${i}.sh"

    cat <<EOF > $script_name
#!/bin/bash
#SBATCH -p fat4way               
#SBATCH -o protein_${i}_%j.out    
#SBATCH -e protein_${i}_%j.err    
#SBATCH -N 1                     
#SBATCH -c 4                     
#SBATCH -A huzeng_g1             
#SBATCH --qos=huzengf4w          

matlab -batch "core_matlab('C5Tg4h', 'nuclei_protein_registration', '${i}', 3072, 34, 1, 3, 10, '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata', '01_data', '02_registration', 'log', 'protein_round', 'IF', 'protein_stains', {'dapi', 'cona', 'flamingo'})"
EOF

    echo "Generated $script_name"
done
