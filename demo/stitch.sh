for pos in $(seq -f "Position%03g" 1 1)
do
    matlab -batch "core_matlab('p1test', 'stitch', '$pos', 3072, 42, 1, 3, 11, '/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain', '01_data', '02_registration', 'log', 'spotfinding_method', 'max3d', 'sqrt_pieces', 4)"
done

