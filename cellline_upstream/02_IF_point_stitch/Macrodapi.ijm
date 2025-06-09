print("Starting Grid/Collection stitching...");
run("Grid/Collection stitching", "type=[Grid: row-by-row] order=[Right & Up] grid_size_x=7 grid_size_y=7 tile_overlap=10 first_file_index_i=1 directory=/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/dapi file_names=Position{iii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
print("Finished Grid/Collection stitching.");

saveAs("Tiff", "/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/dapi/dapibig3dnew.tif");

run("Z Project...", "projection=[Max Intensity]");
saveAs("Tiff", "/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/dapi/dapibig2dnew.tif");

close();
close();
