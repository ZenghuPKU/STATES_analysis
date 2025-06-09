#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o stitchdapi_%j.out       
#SBATCH -e stitchdapi_%j.err       
#SBATCH -N 1 
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 10

/lustre2/huzeng_pkuhpc/yly/Fiji.app/ImageJ-linux64 --headless --console -macro fiji2.ijm > log2.txt
