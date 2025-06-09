#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o 3dmask_%j.out        
#SBATCH -e 3dmask_%j.err        
#SBATCH -N 1 
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 10


module load matlab/2021a
matlab -nodisplay -nosplash -r "run('IFmask3D_new015.m'); exit;"
