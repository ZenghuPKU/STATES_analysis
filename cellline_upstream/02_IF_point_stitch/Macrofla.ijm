print("Starting Grid/Collection stitching...");

run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory=/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/flamingo layout_file=TileConfiguration.registered.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");

print("Finished Grid/Collection stitching.");

saveAs("Tiff", "/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/flamingo/flamingobig3dnew.tif");

run("Z Project...", "projection=[Max Intensity]");
saveAs("Tiff", "/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/flamingo/flamingobig2dnew.tif");

close();
close();
