import numpy as np
import pandas as pd
import scanpy as sc
from tifffile import imread
from scipy.ndimage import distance_transform_edt, label
from tqdm import tqdm

adata = sc.read("../01_preprocessing/8mAD_mousebrain.h5ad")
plaque_mask = imread("8mADabeta/05_binary_filtered.tif")
plaque_features = pd.read_csv("8mADabeta/plaque_features_filtered.csv")

label_to_plaque_id = dict(zip(plaque_features['label_id'], plaque_features['plaque_id']))

labeled_plaques, _ = label(plaque_mask > 0)

print("Distance transforming...")
dist_pix, indices = distance_transform_edt(plaque_mask == 0, return_indices=True)

PIXEL_SIZE_UM = 94.6 / 1000 
dist_um = dist_pix * PIXEL_SIZE_UM
pixel_plaque_labels = labeled_plaques[indices[0], indices[1]]

bins = [0, 10, 20, 30, 40, 50, np.inf]
ring_labels = ["ring1_0_10um", "ring2_10_20um", "ring3_20_30um", "ring4_30_40um", "ring5_40_50um", "outer_50um_plus"]

print("Updating adata.obs ...")
x = adata.obs["column"].astype(int).values
y = adata.obs["row"].astype(int).values

adata.obs["assigned_plaque_id"] = pd.Series(pixel_plaque_labels[y, x]).map(label_to_plaque_id).fillna(0).values
adata.obs["min_border_dist_um"] = dist_um[y, x]
ring_idx = np.digitize(adata.obs["min_border_dist_um"], bins) - 1
adata.obs["ring_assignment"] = [ring_labels[i] for i in ring_idx]

print("Calculating...")
valid_pixels = pixel_plaque_labels > 0
df_pixels = pd.DataFrame({
    'label_id': pixel_plaque_labels[valid_pixels].ravel(),
    'dist': dist_um[valid_pixels].ravel()
})
df_pixels['plaque_id'] = df_pixels['label_id'].map(label_to_plaque_id)
df_pixels['ring_idx'] = np.digitize(df_pixels['dist'], bins) - 1

PIXEL_AREA_UM2 = PIXEL_SIZE_UM ** 2
area_stats = df_pixels.groupby(['plaque_id', 'ring_idx']).size().reset_index(name='pixel_count')
area_stats['area_um2'] = area_stats['pixel_count'] * PIXEL_AREA_UM2

print("Summarizing...")
cell_stats = adata.obs[adata.obs['assigned_plaque_id'] > 0].copy()
cell_summary = cell_stats.groupby(['assigned_plaque_id', 'ring_assignment']).apply(
    lambda g: pd.Series({
        'cell_count': len(g),
        'cell_ids': ";".join(g.index.astype(str))
    })
).reset_index()
cell_summary.rename(columns={'assigned_plaque_id': 'plaque_id', 'ring_assignment': 'ring'}, inplace=True)

print("Combining and saving...")
area_stats['ring'] = area_stats['ring_idx'].map(lambda i: ring_labels[i])
final_csv = pd.merge(
    area_stats[['plaque_id', 'ring', 'area_um2']], 
    cell_summary, 
    on=['plaque_id', 'ring'], 
    how='left'
)

final_csv['cell_count'] = final_csv['cell_count'].fillna(0).astype(int)
final_csv['cell_ids'] = final_csv['cell_ids'].fillna("")

final_csv = final_df = final_csv[['plaque_id', 'ring', 'area_um2', 'cell_count', 'cell_ids']]
final_csv.to_csv("8mADabeta/8mAD_plaque_ring_area_and_cells.csv", index=False)

adata.write("8mAD_mousebrain_with_plaque_info.h5ad")

print("Done!")