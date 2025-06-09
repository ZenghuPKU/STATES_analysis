% merge dots 
stitch_file = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/dapi/TileConfiguration.registered.txt';
%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [5, Inf];
opts.Delimiter = ["(", ")", ",", ";"];

% Specify column names and types
opts.VariableNames = ["Definetheimagecoordinates", "Var2", "Var3", "VarName4", "VarName5", "VarName6"];
opts.SelectedVariableNames = ["Definetheimagecoordinates", "VarName4", "VarName5", "VarName6"];
opts.VariableTypes = ["double", "string", "string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var3"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Definetheimagecoordinates", "TrimNonNumeric", true);
opts = setvaropts(opts, "Definetheimagecoordinates", "ThousandsSeparator", ",");

% Import the data
TileConfiguration = readtable(stitch_file, opts);
TileConfiguration.Properties.VariableNames = {'tile', 'x', 'y', 'z'};
TileConfiguration.x = int32(fix(TileConfiguration.x));
TileConfiguration.y = int32(fix(TileConfiguration.y));
TileConfiguration.z = int32(fix(TileConfiguration.z));

head(TileConfiguration, 5)

%% Clear temporary variables
clear opts

% Global offsets 
offset_x = abs(min(TileConfiguration.x));
offset_y = abs(min(TileConfiguration.y));
offset_z = abs(min(TileConfiguration.z));

TileConfiguration.x = TileConfiguration.x + offset_x + 1;
TileConfiguration.y = TileConfiguration.y + offset_y + 1;
TileConfiguration.z = TileConfiguration.z + offset_z + 1;
%% Merge dots 


% load stitched dapi image
dapi_file = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration/IF/dapi/dapibig2dnew.tif';
dapi_2d = imread(dapi_file);

if ndims(dapi_2d) == 3
    dapi_max = max(dapi_2d, [], 3);
else
    dapi_max = dapi_2d;
end

% create empty holders 
merged_points = int32([]);
merged_genes = {};
merged_region = [];
%%

for r=1:size(TileConfiguration,1)
    
    % get each row of tileconfig
    curr_row = table2cell(TileConfiguration(r, :));
    [tile, x, y, z] = curr_row{:};
    
    % load dots of each tile
    curr_position_dir = sprintf("Position%03d", tile);
    curr_dot_file = fullfile('/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/02_registration', curr_position_dir, "goodPoints_max3d.csv");
    tile_data = readtable(curr_dot_file); % Use readtable to read CSV file
    
    if ~isempty(tile_data)
        tile_goodSpots = tile_data{:, 1:3}; % Extract x, y, z columns
        tile_genes = tile_data{:, 4}; % Extract the gene information
        

        
        tile_goodSpots = int32(tile_goodSpots) + int32([x y z]);
     
        % construct dots region
        current_min = min(tile_goodSpots, [], 1);
        current_max = max(tile_goodSpots, [], 1);
        if current_max(1) > size(dapi_max, 2)
            current_max(1) = size(dapi_max, 2);
            toKeep = tile_goodSpots(:, 1) <= current_max(1);
            tile_goodSpots = tile_goodSpots(toKeep, :);
            tile_genes = tile_genes(toKeep);

        elseif current_max(2) > size(dapi_max, 1)
            current_max(2) = size(dapi_max, 1);   
            toKeep = tile_goodSpots(:, 2) <= current_max(2);
            tile_goodSpots = tile_goodSpots(toKeep, :);
            tile_genes = tile_genes(toKeep);

        end
        
        current_region = zeros(size(dapi_max));
        current_region(current_min(2):current_max(2), current_min(1):current_max(1)) = 1;

        % merge dots 
        if isempty(merged_region)
            merged_region = current_region;
        else
            current_overlap = merged_region & current_region;
            merged_region = merged_region | current_region;
            current_region = current_region - current_overlap; 
            
            temp_cell = num2cell(tile_goodSpots, 2); 
            current_lindex = cellfun(@(x) sub2ind([size(dapi_2d)], x(2), x(1)), temp_cell);
            current_logical = logical(current_region(current_lindex));
            tile_goodSpots = tile_goodSpots(current_logical, :);
            tile_genes = tile_genes(current_logical);
        end
        
        % save current 
        merged_points = [merged_points; tile_goodSpots];
        merged_genes = [merged_genes; cellstr(tile_genes)]; % Ensure gene information is a cell array of strings
    end
end
%%
merged_points(:,1) = merged_points(:,1) - 1;
merged_points(:,2) = merged_points(:,2) - 1;
merged_points(:,3) = merged_points(:,3) - 1;



%%
spot_out_path = '/home/huzeng_pkuhpc/gpfs3/yly/20250101rj/alldata/C5Tg4h/merged_spots';

if ~exist(spot_out_path, 'dir')
    mkdir(spot_out_path)
end

% Convert merged_points to a table
merged_points_table = array2table(merged_points);

% Flatten merged_genes
%flattened_genes = vertcat(merged_genes{:});

% Ensure flattened_genes is a column cell array
%flattened_genes = reshape(flattened_genes, [], 1);

% Convert merged_genes to a table
merged_genes_table = cell2table(merged_genes, 'VariableNames', {'Gene'});



% Ensure the number of rows in both tables are the same
if height(merged_points_table) ~= height(merged_genes_table)
    error('The number of rows in merged_points_table and merged_genes_table must be the same.');
end

% Concatenate the tables
merged_table = [merged_points_table, merged_genes_table];

% Define column names for the table
merged_table.Properties.VariableNames = {'x', 'y', 'z', 'Gene'};

% Write to CSV
writetable(merged_table, fullfile(spot_out_path, 'merged_goodPoints_max3d.csv'));
