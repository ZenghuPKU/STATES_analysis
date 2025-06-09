import os, sys, copy, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance
from tqdm import tqdm
from glob import glob
import matplotlib
from itertools import chain
import warnings
import argparse
warnings.filterwarnings("ignore")

############################# FUNCTIONS #############################

# Read in FIJI Tile Coordinates (TileConfiguration.registered.txt and TileConfiguration.txt)
def get_coords(path, grid=False):
      f = open(path) 
      line = f.readline()
      list = []
      while line:
            if line.startswith('tile'):
                  a = np.array(line.replace('tile_','').replace('.tif; ; (',',').replace(', ',',').replace(')\n','').split(','))
                  if not grid:
                        a = [float(x) for x in a] # Remove rounding for raw coordinates
                  else:
                        a = a.astype(float)
                        tile_num = int(a[0])
                        if tile_num <= 11:  # 第一列
                            grid_x = 0
                            grid_y = tile_num - 1
                        else:
                            grid_x = (tile_num - 1) // 11  
                            grid_y = (tile_num - 1) % 11  
                        a = [a[0], grid_x, grid_y]
                  list.append(a)
            line = f.readline()
      coords_df = np.array(list)
      f.close
      return coords_df

# Find closest cell center for multiassigned reads
def closest_node(node, nodes):
      closest_index = distance.cdist([node], nodes).argmin()
      return nodes[closest_index], closest_index

def idxs_within_coords(df, top_edge, bottom_edge, right_edge, left_edge, include=True):
      if include:
            idxs = df[(df['column'] <= right_edge) & (df['column'] >= left_edge) & (df['row'] >= top_edge) & (df['row'] <= bottom_edge)]['cell_barcode'].tolist()
      else:
            idxs = df[(df['column'] < right_edge) & (df['column'] > left_edge) & (df['row'] > top_edge) & (df['row'] < bottom_edge)]['cell_barcode'].tolist()
      return idxs

def filter_multi_assign(grouping):
    
      (id, df) = grouping
      xyz = id.split('-')[0:3] # xyz coord of multi-assigned reads
      # identify cells that reads are assigned to
      repeat_reads_cell_index = cell_center['cell_barcode'].isin(df['cell_barcode']) 
      repeat_reads_cell = cell_center.loc[repeat_reads_cell_index,:]
      # calculate closest cell (according to cell center)
      closest_index = closest_node(xyz,np.array(repeat_reads_cell.loc[:,['column', 'row', 'z_axis']]).tolist())[1] 
      selected_cell = repeat_reads_cell.iloc[closest_index,0] # barcode of closest cell
      # get indices of other reads to be filtered from farther cells
      filtered_read_idxs = df.index[np.logical_not(remain_reads.loc[df.index,'cell_barcode'] == selected_cell)].tolist() 
      return({id:filtered_read_idxs})

# Get figsize according to image size
def get_figsize(coords_df_tuned, scale=5):
      max_col_coord = coords_df_tuned['column_coord_obs'].max()
      max_row_coord = coords_df_tuned['row_coord_obs'].max()
      return([max_col_coord/100/scale * 2, max_row_coord/100/scale * 2])


def command_args():
    desc = "stitch reads"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-IXY', '--input_xy', type=str, required=True,help='XY')
    parser.add_argument('-ID', '--input_dir', type=str, required=True, help='the input_dir')
    parser.add_argument('-IO', '--input_orderlist', type=str, required=True, help='input_orderlist')
    parser.add_argument('-IS', '--input_segmentation', type=str, required=True, help='input_segmentation')
    parser.add_argument('-OC', '--output_cell_center', type=str, required=True, help='output_cell_center')
    parser.add_argument('-OD', '--output_dir', type=str, required=True, help='output_cell_center')
    parser.add_argument('-OR', '--output_remain_reads', type=str, required=True, help='output_remain_reads')
    parser.add_argument('-ITR', '--Tile_registered', type=str, required=True, help='Tile_registered')
    parser.add_argument('-IT', '--Tile', type=str, required=True, help='TileConfiguration.txt')
    parser.add_argument('-SM', '--seg_method', type=str, required=True, help='clustermap or watershed')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
      args = command_args()
      alignment_thresh = 0.5
      img_c, img_r = [int(args.input_xy), int(args.input_xy)]
      data_dir = args.input_dir
      outputpath = args.output_dir
      os.makedirs(outputpath, exist_ok=True)
      orderlist_path = args.input_orderlist
      readspath = args.input_segmentation
      cell_center_path = args.output_cell_center
      remain_reads_path = args.output_remain_reads
      


      cell_barcode_min = 0
      middle_edge = 0
      remain_reads = pd.DataFrame({'gene_name':[],'spot_location_1':[],'spot_location_2':[],'spot_location_3':[],'gene':[],'is_noise':[],'cell_barcode':[],'raw_cell_barcode':[],'gridc_gridr_tilenum':[]})
      cell_center = pd.DataFrame({'cell_barcode':[],'column':[],'row':[],'z_axis':[],'raw_cell_barcode':[],'gridc_gridr_tilenum':[]})

      #### Read in tile coordinates and orderlist
      # Read coordinates
      obs_coords = get_coords(args.Tile_registered)
      exp_coords = get_coords(args.Tile)
      grid_order = get_coords(args.Tile, grid=True)

      # Format into dataframes
      coords_df = pd.DataFrame(obs_coords, columns=['tile','column_coord_obs','row_coord_obs'])
      coords_exp_df = pd.DataFrame(exp_coords, columns=['tile','column_coord_exp','row_coord_exp'])
      grid_df = pd.DataFrame(grid_order, columns=['tile', 'column_count', 'row_count'])

      # Convert tile numbers to integers
      coords_df['tile'] = coords_df['tile'].astype(int)
      coords_exp_df['tile'] = coords_exp_df['tile'].astype(int)
      grid_df['tile'] = grid_df['tile'].astype(int)

      # Merge coordinates data - keep the tile column for merging
      coords_df = coords_df.merge(
            coords_exp_df,
            on='tile'
      )
      coords_df = coords_df.merge(
            grid_df,
            on='tile'
      )
      
      # Save coordinates
      coords_df.to_csv(os.path.join(outputpath,'coords.csv'))

      # Zero-center (tune) registered coordinates
      print("Tuning coordinates...")
      coords_df_without_blank = coords_df.loc[coords_df['tile'] > 0,:]
      min_column, min_row = [np.min(coords_df_without_blank['column_coord_obs']), np.min(coords_df_without_blank['row_coord_obs'])]
      max_column, max_row = [np.max(coords_df_without_blank['column_coord_obs']), np.max(coords_df_without_blank['row_coord_obs'])]
      shape_column, shape_row = [max_column - min_column + img_c, max_row - min_row + img_r]
      
      coords_df_tuned = copy.deepcopy(coords_df)
      coords_df_tuned['column_coord_obs'] = coords_df['column_coord_obs'] - min_column
      coords_df_tuned['row_coord_obs'] = coords_df['row_coord_obs'] - min_row

      # Calculate grid dimensions based on actual coordinates
      grid_c = len(np.unique(coords_df_without_blank['column_coord_exp'] // (img_c * 0.9)))
      grid_r = len(np.unique(coords_df_without_blank['row_coord_exp'] // (img_r * 0.9)))

      # Save tuned coordinates
      coords_df_tuned.to_csv(os.path.join(outputpath,'tuned_coords.csv'))

      ### STITCH
      print("============STITCHING============")
      cell_barcode_min = 0
      remain_reads = pd.DataFrame({'gene_name':[],'spot_location_1':[],'spot_location_2':[],'spot_location_3':[],'gene':[],'is_noise':[],'cell_barcode':[],'raw_cell_barcode':[],'gridc_gridr_tilenum':[]})
      cell_center = pd.DataFrame({'cell_barcode':[],'column':[],'row':[],'z_axis':[],'raw_cell_barcode':[],'gridc_gridr_tilenum':[]})

      # Create mapping of tile numbers to grid positions
      tile_to_grid = {}
      for _, row in coords_df_tuned.iterrows():
            if row['tile'] > 0:  # Skip blank tiles
                  tile_num = int(row['tile'])
                  if tile_num <= 11:  
                        grid_x = 0
                        grid_y = tile_num - 1
                  else:
                        grid_x = (tile_num - 1) // 11  
                        grid_y = (tile_num - 1) % 11   
                  tile_to_grid[tile_num] = (grid_x, grid_y)

      sorted_tiles = sorted(coords_df_tuned['tile'].unique())
      sorted_tiles = [t for t in sorted_tiles if t > 0]

      print("\nAll available tiles:", sorted_tiles)
      print("\nGrid mapping:")
      for tile, (grid_x, grid_y) in tile_to_grid.items():
            print(f"Tile {tile}: grid position ({grid_x}, {grid_y})")

      processed_tiles = set()
      for tilenum in sorted_tiles:
            if tilenum <= 0:
                  continue
            
            if tilenum in processed_tiles:
                  print(f"Skipping already processed tile {tilenum}")
                  continue
            processed_tiles.add(tilenum)
            
            t_grid_c, t_grid_r = tile_to_grid[tilenum]
            print(f"\nProcessing tile {tilenum} at grid position ({t_grid_c}, {t_grid_r})")
            
            # Get median column coordinate for approx alignment
            median_col = coords_df_tuned[(coords_df_tuned.column_count == t_grid_c) & (coords_df_tuned.tile != 0)]['column_coord_obs']
            median_row = coords_df_tuned[(coords_df_tuned.row_count == t_grid_r) & (coords_df_tuned.tile != 0)]['row_coord_obs']
            
            if len(median_col) == 0 or len(median_row) == 0:
                  print(f"Warning: No median coordinates found for tile {tilenum}")
                  continue
            
            median_col_coord = np.median(median_col)
            median_row_coord = np.median(median_row)
            
            # Get tile coordinate and grid order information
            try:
                  order = coords_df_tuned.index[(coords_df_tuned['column_count']==t_grid_c) & (coords_df_tuned['row_count']==t_grid_r)][0]
            except IndexError:
                  print(f"Error: Could not find grid position for tile {tilenum}")
                  continue
            
            # Print alignment check information
            upper_left = coords_df_tuned.loc[order, ['column_coord_obs', 'row_coord_obs']]
            left_thresh = median_col_coord - (1+alignment_thresh)*img_c
            right_thresh = median_col_coord + (1+alignment_thresh)*img_c
            upper_thresh = median_row_coord - (1+alignment_thresh)*img_r
            lower_thresh = median_row_coord + (1+alignment_thresh)*img_r
            
            print(f"Alignment check for tile {tilenum}:")
            print(f"Position: {upper_left[0]:.2f}, {upper_left[1]:.2f}")
            print(f"Thresholds: {left_thresh:.2f} < x < {right_thresh:.2f}, {upper_thresh:.2f} < y < {lower_thresh:.2f}")
            
            tile_order = coords_df_tuned['tile'][order] # Tile number (blanks = 0)
            if pd.isna(tile_order) or tile_order == 0: # skip blanks
                  continue
            tile_order = int(tile_order)  # Convert to integer

            upper_left_new = copy.deepcopy(upper_left)

            # If either column or row alignment shifted more than 1.5x tile width/height away from median coordinate, throw out
            if upper_left[0] >= right_thresh or upper_left[0] <= left_thresh or upper_left[1] >= lower_thresh or upper_left[1] <= upper_thresh:
                  print(f"- Tile {tile_order} is aligned too far away from its expected position.")
                  print(f"\tTile coord: [{upper_left[0]}, {upper_left[1]}]. Median coord: [{median_col_coord}, {median_row_coord}]")
                  continue

            # Get tile read and cell center information
            dfpath = readspath + f'/Position{tile_order:03d}'
            print(f"\nProcessing tile {tile_order}:")
            print(f"Reading from path: {dfpath}")
            
            if not os.path.exists(os.path.join(dfpath, 'remain_reads.csv')):
                  print(f"- Tile {tile_order}: Reads file does not exist. [{t_grid_c},{t_grid_r}]")
                  continue
                  
            remain_reads_t = pd.read_csv(os.path.join(dfpath,'remain_reads.csv'),index_col=0)
            cell_center_t = pd.read_csv(os.path.join(dfpath,'cell_center.csv'),index_col=0)
            
            print(f"Initial reads count: {remain_reads_t.shape[0]}")
            print(f"Initial cells count: {cell_center_t.shape[0]}")

            if args.seg_method == 'clustermap':
                  remain_reads_t.rename(columns = {'clustermap':'cell_barcode'}, inplace = True)
            elif args.seg_method == 'watershed':
                  remain_reads_t.rename(
                  columns = {
                    'Gene':'gene', 
                    'x':'spot_location_1',
                    'y':'spot_location_2',
                    'z':'spot_location_3'
                  }, 
                  inplace=True)
                  cell_center_t.rename(
                  columns = {
                    'x':'column',
                    'y':'row'
                  },
                  inplace=True)

            # Format reads data and add in grid/order information   
            remain_reads_t['gridc_gridr_tilenum'] = str(t_grid_c)+","+str(t_grid_r)+","+str(tile_order)
            cell_center_t['gridc_gridr_tilenum'] = str(t_grid_c)+","+str(t_grid_r)+","+str(tile_order)

            # Adjust read coordinates and cell center coordinates using observed upper left tile coordinate
            remain_reads_t['spot_location_1'] = remain_reads_t['spot_location_1'] + upper_left[0]# col
            cell_center_t['column'] = cell_center_t['column']  + upper_left[0]
            remain_reads_t['spot_location_2'] = remain_reads_t['spot_location_2'] + upper_left[1]# col
            cell_center_t['row'] = cell_center_t['row']  + upper_left[1]

            # Keep tile-by-tile barcodes as raw barcodes, cell_barcode is cumulative
            remain_reads_t['raw_cell_barcode'] = remain_reads_t['cell_barcode']
            cell_center_t['raw_cell_barcode'] = cell_center_t['cell_barcode']
            remain_reads_t['cell_barcode'] = remain_reads_t['cell_barcode'] + cell_barcode_min + 1
            cell_center_t['cell_barcode'] = cell_center_t['cell_barcode'] + cell_barcode_min + 1

            top_middle_edge = None
            left_middle_edge = None
            cell_center_drop_barcodes = []
            cell_center_t_drop_barcodes = []

            # Evaluate previous neighbors (i.e. tiles directly left and up from current tile)
            t_grid_c_previous = t_grid_c - 1
            t_grid_r_previous = t_grid_r - 1

            # Get middle edge with left tile if it isn't blank/doesn't exist
            if t_grid_c_previous >= 0:  # current tile not on left edge
                left_indices = coords_df_tuned.index[
                    (coords_df_tuned['column_count']==t_grid_c_previous) & 
                    (coords_df_tuned['row_count']==t_grid_r)
                ]
                if len(left_indices) > 0:
                    order_left = left_indices[0]
                    if coords_df_tuned.loc[order_left,'tile'] != 0:
                        upper_left_left = coords_df_tuned.loc[order_left, ['column_coord_obs', 'row_coord_obs']]
                        
                        left_overlap_left = upper_left[0] 
                        left_overlap_right = upper_left_left[0] + img_c 
                        left_middle_edge = (left_overlap_left + left_overlap_right) / 2
                        
                        if left_middle_edge >= left_overlap_left:
                            cell_center_drop_barcodes += idxs_within_coords(
                                cell_center, 
                                upper_left[1],  # top
                                upper_left[1] + img_r,  # bottom
                                left_overlap_right,  # right
                                left_middle_edge,  # left
                                include=True
                            )
                            
                            cell_center_t_drop_barcodes += idxs_within_coords(
                                cell_center_t,
                                upper_left[1],  # top
                                upper_left[1] + img_r,  # bottom
                                left_middle_edge,  # right
                                left_overlap_left,  # left
                                include=True
                            )

            # Get middle edge with top tile if it isn't blank/doesn't exist
            if t_grid_r_previous >= 0:  # current tile not on top edge (first row)
                top_indices = coords_df_tuned.index[
                    (coords_df_tuned['column_count']==t_grid_c) & 
                    (coords_df_tuned['row_count']==t_grid_r_previous)
                ]
                if len(top_indices) > 0:
                    order_top = top_indices[0]
                    if coords_df_tuned.loc[order_top,'tile'] != 0:  # check top tile not blank
                        # Calculate overlap region between new tile and top tile
                        upper_left_top = coords_df_tuned.loc[order_top, ['column_coord_obs', 'row_coord_obs']]
                        new_right = upper_left[0] > upper_left_top[0]
                        top_overlap_top = upper_left[1]  # row coord
                        top_overlap_bot = upper_left_top[1] + img_r  # row coord
                        top_overlap_left = upper_left[0]  # col coord
                        top_overlap_right = upper_left_top[0] + img_c if new_right else upper_left[0] + img_c  # col coord
                        top_middle_edge = top_overlap_top + (top_overlap_bot - top_overlap_top) / 2

                        overlap_cells_old = cell_center[
                            (cell_center['row'] >= top_overlap_top) & 
                            (cell_center['row'] <= top_overlap_bot) &
                            (cell_center['column'] >= top_overlap_left) &
                            (cell_center['column'] <= top_overlap_right)
                        ]
                        
                        overlap_cells_new = cell_center_t[
                            (cell_center_t['row'] >= top_overlap_top) & 
                            (cell_center_t['row'] <= top_overlap_bot) &
                            (cell_center_t['column'] >= top_overlap_left) &
                            (cell_center_t['column'] <= top_overlap_right)
                        ]

                        for _, old_cell in overlap_cells_old.iterrows():
                            if len(overlap_cells_new) > 0:
                                distances = np.sqrt(
                                    (overlap_cells_new['column'] - old_cell['column'])**2 + 
                                    (overlap_cells_new['row'] - old_cell['row'])**2
                                )
                                min_dist_idx = distances.idxmin()
                                min_dist = distances.min()
                                
                                if min_dist < img_c * 0.1:  
                                    new_cell = overlap_cells_new.loc[min_dist_idx]
                                    old_dist = abs(old_cell['row'] - top_overlap_top)
                                    new_dist = abs(new_cell['row'] - top_overlap_bot)
                                    
                                    if old_dist > new_dist:
                                        cell_center_drop_barcodes.append(old_cell['cell_barcode'])
                                    else:
                                        cell_center_t_drop_barcodes.append(new_cell['cell_barcode'])

            # Add debug information
            print(f"\nProcessing overlaps for tile {tilenum}:")
            print(f"Left edge exists: {left_middle_edge is not None}")
            print(f"Top edge exists: {top_middle_edge is not None}")

            rescued_barcodes = []
            for i in np.unique(cell_center_drop_barcodes):
                row_coord, col_coord = cell_center.loc[cell_center['cell_barcode']==i, ['row', 'column']].values[0]
                if top_middle_edge is not None: 
                    if row_coord < top_middle_edge:
                        rescued_barcodes.append(i)
                if left_middle_edge is not None: 
                    if col_coord < left_middle_edge:
                        rescued_barcodes.append(i)
            cell_center_drop_barcodes = [idx for idx in cell_center_drop_barcodes if idx not in rescued_barcodes]
        
            # Drop cells and associated reads
            cell_center_drop_barcodes = np.unique(cell_center_drop_barcodes)
            cell_center_t_drop_barcodes = np.unique(cell_center_t_drop_barcodes)
            previous_cell_count = len(cell_center)
            matched_cell_count = len(cell_center[cell_center['cell_barcode'].isin(cell_center_drop_barcodes)].index)
            #print(f'Number of unique barcodes to drop: {len(cell_center_drop_barcodes)} | Number of barcodes that match in cell_center: {matched_cell_count}')
            if len(cell_center_drop_barcodes) != matched_cell_count:
                  print(f"================== WARNING! Number of unique barcodes to drop: {len(cell_center_drop_barcodes)} != Number of barcodes that match in cell_center: {matched_cell_count}")
            #cell_center = cell_center[~cell_center['cell_barcode'].isin(cell_center_drop_barcodes)] # using this method because indices are not unique
            cell_center.drop(cell_center[cell_center['cell_barcode'].isin(cell_center_drop_barcodes)].index, inplace=True)
            if len(cell_center_drop_barcodes) != previous_cell_count - len(cell_center):
                  print(f'{len(cell_center_drop_barcodes)} cells should have been dropped from old region (reality: {previous_cell_count - len(cell_center)})')
            remain_reads = remain_reads.loc[remain_reads['cell_barcode'].isin(cell_center['cell_barcode']),:]
            cell_center_t.drop(cell_center_t[cell_center_t['cell_barcode'].isin(cell_center_t_drop_barcodes)].index, inplace=True)
            remain_reads_t = remain_reads_t.loc[remain_reads_t['cell_barcode'].isin(cell_center_t['cell_barcode']),:]
            #print(f'{len(cell_center_drop_barcodes)} cells dropped from old region')
            #print(f'{len(cell_center_t_drop_barcodes)} cells dropped from new tile')

            ## append
            cell_center = cell_center.append(cell_center_t)
            cell_center.reset_index(inplace=True, drop=True)
            remain_reads = remain_reads.append(remain_reads_t)
            remain_reads.reset_index(inplace=True, drop=True)
            if cell_center_t.shape[0] > 0:
                  cell_barcode_min = np.max(cell_center_t['cell_barcode']) + 1

            print(f"After overlap processing:")
            print(f"Remaining reads count: {remain_reads_t.shape[0]}")
            print(f"Remaining cells count: {cell_center_t.shape[0]}")
            print(f"Cells dropped from old region: {len(cell_center_drop_barcodes)}")
            print(f"Cells dropped from new tile: {len(cell_center_t_drop_barcodes)}")

            print(f"\nAlignment check for tile {tile_order}:")
            print(f"Position: {upper_left[0]:.2f}, {upper_left[1]:.2f}")
            print(f"Median coord: [{median_col_coord}, {median_row_coord}]")
            print(f"Thresholds: left={left_thresh:.2f}, right={right_thresh:.2f}")
            print(f"           upper={upper_thresh:.2f}, lower={lower_thresh:.2f}")

            print(f"\nFinal merge for tile {tile_order}:")
            print(f"Current total reads: {remain_reads.shape[0]}")
            print(f"Current total cells: {cell_center.shape[0]}")
            print(f"Adding reads: {remain_reads_t.shape[0]}")
            print(f"Adding cells: {cell_center_t.shape[0]}")

      print(f"\tTotal Reads: {remain_reads.shape[0]} | Total Cells: {cell_center.shape[0]}")

      ################################ polish after stitch ################################
      print("============FILTERING============")
      print("Finding multi-assigned reads...")
      # filter the repeated reads
      remain_reads = remain_reads.drop_duplicates(subset = None, keep = 'first')
      cell_center.reset_index(inplace = True,drop = True)
      remain_reads.reset_index(inplace = True,drop = True)
      remain_reads.rename(columns = {'spot_location_1':'column', 'spot_location_2':'row','spot_location_3':'z'}, inplace = True)
      
      cell_center.rename(columns = {'z_axis':'z'}, inplace = True)

      remain_reads = remain_reads.astype(
      {'column':'int', 'row':'int', 'z':'int', 'cell_barcode':'int', 'raw_cell_barcode':'int'}
            )
      cell_center = cell_center.astype(
      {'cell_barcode':'int', 'column':'int', 'row':'int', 'z':'int', 'raw_cell_barcode':'int'}
            )
      #remain_reads = remain_reads.drop(['gene_name', 'is_noise'], axis=1)
      
      

      print("Read counts after filtering multi-assigned reads: " + str(remain_reads.shape[0]))
      cell_center.to_csv(cell_center_path)
      remain_reads.to_csv(remain_reads_path)

      plt.style.use('dark_background')

      plt.figure(figsize=get_figsize(coords_df_tuned, scale=5))
      plt.subplot(2,2,1)
      plt.title('Reads with cell centers')
      #plt.scatter(remain_reads.loc[:,'column'],shape_row - remain_reads.loc[:,'row'],s = 0.1,alpha = 0.2,color=remain_reads['cell_barcode'])
      plt.scatter(remain_reads.loc[:,'column'],shape_row - remain_reads.loc[:,'row'],s = 0.1,alpha = 0.8,c=pd.Categorical(np.array(remain_reads['raw_cell_barcode'])).codes, cmap= matplotlib.colors.ListedColormap ( np.random.rand ( 256,3)))
      plt.scatter(cell_center.loc[:,'column'],shape_row - cell_center.loc[:,'row'],s = 1,c='red',alpha = 1)
      plt.axis('off')

      plt.subplot(2,2,2)
      plt.title('Reads with cell centers and tile order')
      #plt.scatter(remain_reads.loc[:,'column'],shape_row - remain_reads.loc[:,'row'],s = 0.1,alpha = 0.2)
      plt.scatter(remain_reads.loc[:,'column'],shape_row - remain_reads.loc[:,'row'],s = 0.1,alpha = 0.8,c=pd.Categorical(np.array(remain_reads['raw_cell_barcode'])).codes, cmap= matplotlib.colors.ListedColormap ( np.random.rand ( 256,3)))
      plt.scatter(cell_center.loc[:,'column'],shape_row - cell_center.loc[:,'row'],s = 1,c='red',alpha = 1)
      plt.axis('off')
      y_reverse = True
      list_t = ['column_coord_obs','row_coord_obs','tile']
      coords_df = coords_df_tuned.copy()
      idx_t = coords_df[list_t[2]] != 0
      coords0 = np.array(coords_df.loc[idx_t,list_t])
      if y_reverse:
            coords0[:,1] = coords0[:,1].max() - coords0[:,1] 
      plt.scatter(x=coords0[:,0],y=coords0[:,1],c=coords0[:,2])
      for i in range(coords0.shape[0]):
            plt.text(x=coords0[i,0],y=coords0[i,1],s=coords0[i,2],fontdict=dict(fontsize=20))

      plt.subplot(2,2,3)
      plt.title('Cell centers')
      plt.scatter(cell_center.loc[:,'column'],shape_row - cell_center.loc[:,'row'],s = 10,c='red',alpha = 0.8)
      plt.axis('off')

      plt.subplot(2,2,4)
      plt.title('Cell centers and tile order')
      plt.scatter(cell_center.loc[:,'column'],shape_row - cell_center.loc[:,'row'],s = 10,c='red',alpha = 0.8)
      plt.scatter(x=coords0[:,0],y=coords0[:,1],c=coords0[:,2])
      for i in range(coords0.shape[0]):
            plt.text(x=coords0[i,0],y=coords0[i,1],s=coords0[i,2],fontdict=dict(fontsize=20))
      plt.axis('off')

      plt.tight_layout()

      plt.savefig(os.path.join(outputpath,'cell_reads_profile.png'))

def process_overlap_region(cell_center_old, cell_center_new, overlap_start, overlap_end, axis='column'):
    overlap_middle = (overlap_start + overlap_end) / 2
    cells_to_drop_old = []
    cells_to_drop_new = []
    
    old_cells = cell_center_old[(cell_center_old[axis] >= overlap_start) & 
                               (cell_center_old[axis] <= overlap_end)]
    new_cells = cell_center_new[(cell_center_new[axis] >= overlap_start) & 
                               (cell_center_new[axis] <= overlap_end)]
    
    for _, old_cell in old_cells.iterrows():
        for _, new_cell in new_cells.iterrows():
            dist = np.sqrt((old_cell['column'] - new_cell['column'])**2 + 
                         (old_cell['row'] - new_cell['row'])**2)
            
            if dist < img_c * 0.1: 
                old_dist = abs(old_cell[axis] - overlap_middle)
                new_dist = abs(new_cell[axis] - overlap_middle)
                
                if old_dist > new_dist:
                    cells_to_drop_old.append(old_cell['cell_barcode'])
                else:
                    cells_to_drop_new.append(new_cell['cell_barcode'])
    
    return cells_to_drop_old, cells_to_drop_new

if t_grid_c_previous >= 0:  
    cells_to_drop_old, cells_to_drop_new = process_overlap_region(
        cell_center, 
        cell_center_t,
        left_overlap_left,
        left_overlap_right,
        'column'
    )
    cell_center_drop_barcodes.extend(cells_to_drop_old)
    cell_center_t_drop_barcodes.extend(cells_to_drop_new)

if t_grid_r_previous >= 0:  
    cells_to_drop_old, cells_to_drop_new = process_overlap_region(
        cell_center, 
        cell_center_t,
        top_overlap_top,
        top_overlap_bot,
        'row'
    )
    cell_center_drop_barcodes.extend(cells_to_drop_old)
    cell_center_t_drop_barcodes.extend(cells_to_drop_new)
