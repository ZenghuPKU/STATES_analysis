#!/bin/bash
#SBATCH -p fat4way
#SBATCH -o Position001_subtile5_%j.out
#SBATCH -e Position001_subtile5_%j.err
#SBATCH -N 1
#SBATCH --no-requeue
#SBATCH -A huzeng_g1
#SBATCH --qos=huzengf4w
#SBATCH -c 4

matlab -batch "core_matlab('p1test', 'local_registration', 'Position001', 3072, 42, 1, 3, 11, '/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain', '01_data', '02_registration', 'log', 'spotfinding_method', 'max3d', 'sqrt_pieces', 4, 'subtile', 5, 'voxel_size', [1, 1, 1], 'end_bases', ['G','G','A','A','A'], 'barcode_mode', 'tri','split_loc', [6,12],'intensity_threshold', 0.2)"
