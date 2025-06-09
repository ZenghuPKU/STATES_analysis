# conda: bigfish_env
import numpy as np
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import spatialdata as sd
import scanpy as sc

c3data = sd.read_zarr('cellline0101_C3control.zarr')
b4data = sd.read_zarr('cellline0101_B4Tg15min.zarr')
b5data = sd.read_zarr('cellline0101_B5Tg30min.zarr')
b6data = sd.read_zarr('cellline0101_B6Tg1h.zarr')
c4data = sd.read_zarr('cellline0101_C4Tg2h.zarr')
c5data = sd.read_zarr('cellline0101_C5Tg4h.zarr')

datac3 = c3data.table
datab4 = b4data.table
datab5 = b5data.table
datab6 = b6data.table
datac4 = c4data.table
datac5 = c5data.table

c3cell_mask = c3data.labels["3D_cell_mask"]
c3cell_mask_np = c3cell_mask.data.compute() 
c3cell_areas = np.bincount(c3cell_mask_np.ravel())[1:]  
datac3.obs["area"] = c3cell_areas
print('c3 done')

b4cell_mask = b4data.labels["3D_cell_mask"]
b4cell_mask_np = b4cell_mask.data.compute()
b4cell_areas = np.bincount(b4cell_mask_np.ravel())[1:]  
datab4.obs["area"] = b4cell_areas
print('b4 done')

b5cell_mask = b5data.labels["3D_cell_mask"]
b5cell_mask_np = b5cell_mask.data.compute()
b5cell_areas = np.bincount(b5cell_mask_np.ravel())[1:] 
datab5.obs["area"] = b5cell_areas
print('b5 done')

b6cell_mask = b6data.labels["3D_cell_mask"]
b6cell_mask_np = b6cell_mask.data.compute()
b6cell_areas = np.bincount(b6cell_mask_np.ravel())[1:]  
datab6.obs["area"] = b6cell_areas
print('b6 done')

c4cell_mask = c4data.labels["3D_cell_mask"]
c4cell_mask_np = c4cell_mask.data.compute()
c4cell_areas = np.bincount(c4cell_mask_np.ravel())[1:] 
datac4.obs["area"] = c4cell_areas
print('c4 done')

c5cell_mask = c5data.labels["3D_cell_mask"]
c5cell_mask_np = c5cell_mask.data.compute()
c5cell_areas = np.bincount(c5cell_mask_np.ravel())[1:] 
datac5.obs["area"] = c5cell_areas
print('c5 done')

print('area done')
datac3.obs['condition'] = 'C3control'
datab4.obs['condition'] = 'B4Tg15min'
datab5.obs['condition'] = 'B5Tg30min'
datab6.obs['condition'] = 'B6Tg1h'
datac4.obs['condition'] = 'C4Tg2h'
datac5.obs['condition'] = 'C5Tg4h'

datac3.obs['sample'] = 'C3control'
datab4.obs['sample'] = 'B4Tg15min'
datab5.obs['sample'] = 'B5Tg30min'
datab6.obs['sample'] = 'B6Tg1h'
datac4.obs['sample'] = 'C4Tg2h'
datac5.obs['sample'] = 'C5Tg4h'

adata_combined = datac3.concatenate(datab4, datab5, datab6, datac4, datac5, batch_key="condition")


adata_combined_ntRNA = adata_combined[:, adata_combined.var['feature_name'].str.endswith('_ntRNA')]
adata_combined_rbRNA = adata_combined[:, adata_combined.var['feature_name'].str.endswith('_rbRNA')]

new_feature_names = adata_combined.var['feature_name'].str.replace(r'(_rbRNA|_ntRNA)', '', regex=True).drop_duplicates()


new_X = adata_combined_ntRNA.X + adata_combined_rbRNA.X


new_data = sc.AnnData(X=new_X, var=pd.DataFrame(index=new_feature_names),obs=adata_combined.obs.copy())

new_data.layers['ntRNA'] = adata_combined_ntRNA.X
new_data.layers['rbRNA'] = adata_combined_rbRNA.X

new_data.uns = adata_combined.uns.copy()

new_data.write_h5ad('combinedall.h5ad')