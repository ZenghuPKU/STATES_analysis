
% 2D segmentation reference images
ref_cell_seg_path = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/maskcell/maskcell_segmentation.tif'; 
ref_dapi_seg_path = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/masknucleitorch/masknucleitorch_segmentation.tif';  

% 3D images (stitched images)
overlay_img_path  = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/flamingo/flamingobig3dnew.tif';          
dapi_img_path     = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/dapi/dapibig3dnew.tif';               

cell_output_path       = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/3dmask/cell_new_f015_fgray.tif';
dapi_output_path       = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/3dmask/nuclei_new_f015_fgray.tif';

options.compress  = 'lzw';
options.overwrite = true;
options.big       = true;


% Load the cell segmentation reference image
ref_cell_seg = imread(ref_cell_seg_path);

% IMPORTANT: Ensure that ref_cell_seg has the same dimensions as each slice of the overlay image.
% Load the stitched overlay image (automatically reads image dimensions)
overlay_img  = imread_big(overlay_img_path);

% Apply median filtering to each slice of the overlay image
cell_filter = zeros(size(overlay_img), 'uint8');
for z = 1:size(overlay_img, 3)
    cell_filter(:,:,z) = medfilt2(overlay_img(:,:,z), [10 10]);
end

% Threshold segmentation (using a fixed threshold)
cell_threshold = 0.015;
cell_bw = imbinarize(cell_filter, cell_threshold);

% Generate the cell segmentation image
cell_seg = zeros(size(overlay_img), 'uint16');
se       = strel('disk', 10);
for z = 1:size(overlay_img, 3)
    % Fill holes in the binary mask
    current_mask = imfill(cell_bw(:,:,z), 'holes');
    % Remove small objects (area less than 200 pixels)
    current_mask = bwareaopen(current_mask, 200);
    % Dilate the mask with a disk-shaped structuring element
    current_mask = imdilate(current_mask, se);
    % Fill holes again after dilation
    current_mask = imfill(current_mask, 'holes');
    
    % Multiply the processed mask with the reference segmentation.
    % NOTE: This requires that ref_cell_seg has the same spatial dimensions as overlay_img(:,:,z).
    cell_seg(:,:,z) = uint16(current_mask) .* ref_cell_seg;
end

% Save the cell segmentation result as a TIFF file
saveastiff(cell_seg, cell_output_path, options);

% Load the DAPI segmentation reference image
ref_dapi_seg = imread(ref_dapi_seg_path);

% Load the stitched DAPI image (automatically reads image dimensions)
dapi_img = imread_big(dapi_img_path);

% Apply median filtering to each slice of the DAPI image
dapi_filter = zeros(size(dapi_img), 'uint8');
for z = 1:size(dapi_img, 3)
    dapi_filter(:,:,z) = medfilt2(dapi_img(:,:,z), [10 10]);
    % Alternatively, you can use Gaussian filtering:
    % dapi_filter(:,:,z) = imgaussfilt(dapi_img(:,:,z));
end

% Automatically compute threshold and binarize the DAPI image
dapi_threshold = graythresh(dapi_filter);
dapi_bw = imbinarize(dapi_filter, dapi_threshold);

% Generate the DAPI segmentation image
dapi_seg = zeros(size(dapi_img), 'uint16');
se = strel('disk', 5);
for z = 1:size(dapi_img, 3)
    current_mask = imfill(dapi_bw(:,:,z), 'holes');
    current_mask = bwareaopen(current_mask, 10);
    current_mask = imdilate(current_mask, se);
    % Multiply the processed mask with the DAPI reference segmentation.
    % Ensure that ref_dapi_seg has the same dimensions as dapi_img(:,:,z).
    dapi_seg(:,:,z) = uint16(current_mask) .* ref_dapi_seg;
end

% Save the DAPI segmentation result.
% First, save using saveastiff, then delete the file and write each slice to create a multi-page TIFF.
saveastiff(dapi_seg, dapi_output_path, options);
if exist(dapi_output_path, 'file') == 2
    delete(dapi_output_path);
end
for j = 1:size(dapi_seg, 3)
    imwrite(dapi_seg(:,:,j), dapi_output_path, 'writemode', 'append');
end

