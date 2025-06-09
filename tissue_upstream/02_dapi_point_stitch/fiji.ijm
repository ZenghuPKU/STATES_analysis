run("Grid/Collection stitching", "type=[Grid: column-by-column] order=[Down & Right                ] grid_size_x=12 grid_size_y=11 tile_overlap=10 first_file_index_i=1 directory=/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/stitchlinks file_names=tile_{i}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
print("Finished Grid/Collection stitching.");

saveAs("Tiff", "/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata/B4_WT_mousebrain2/04_stitch/dapibig2d.tif");

close();
