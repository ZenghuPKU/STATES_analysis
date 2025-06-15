import json
import pandas as pd
import numpy as np
import spatialdata as sd
from scipy.spatial import cKDTree

def process_2d_assignment(
    input_zarr_path: str,
    output_zarr_path: str,
    table_key: str = "table",
    nucleus_mask_key: str = '2D_nucleus_mask'
) -> sd.SpatialData:

    sdata = sd.read_zarr(input_zarr_path)
    adata = sdata.tables[table_key]
    points_json = adata.uns["points"]
    df = pd.DataFrame(json.loads(points_json))

    nucleus_mask = sdata.labels[nucleus_mask_key].values
    nucleus_coords = np.argwhere(nucleus_mask)

    cell_mask = df['mask_label'] == 'Cell'
    cell_coords = df.loc[cell_mask, ['y', 'x']].values
    tree = cKDTree(nucleus_coords)
    distances, _ = tree.query(cell_coords, k=1, distance_upper_bound=0.5)
    is_in_nucleus = (distances == 0)
    df.loc[cell_mask, 'mask_label'] = np.where(is_in_nucleus, 'Nuclei', 'Cell')

    nuclei_mask = df['mask_label'] == 'Nuclei'
    df.loc[nuclei_mask, 'nuclei'] = df.loc[nuclei_mask, 'cell']

    updated_points_json = df.to_json(orient='records')
    adata.uns["points"] = updated_points_json
    sdata.tables[table_key] = adata
    sdata.write(output_zarr_path, overwrite=True)

    return sdata

if __name__ == "__main__":
    input_path = "/storage/lingyuan2/STATES_data/cellline0101_C3control"
    output_path = "/storage/lingyuan2/STATES_data/cellline0101_C3control_2d.zarr"
    process_2d_assignment(input_path, output_path)

    input_path = "/storage/lingyuan2/STATES_data/cellline0101_B4Tg15min"
    output_path = "/storage/lingyuan2/STATES_data/cellline0101_B4Tg15min_2d.zarr"
    process_2d_assignment(input_path, output_path)

    input_path = "/storage/lingyuan2/STATES_data/cellline0101_B5Tg30min"
    output_path = "/storage/lingyuan2/STATES_data/cellline0101_B5Tg30min_2d.zarr"
    process_2d_assignment(input_path, output_path)

    input_path = "/storage/lingyuan2/STATES_data/cellline0101_B6Tg1h"
    output_path = "/storage/lingyuan2/STATES_data/cellline0101_B6Tg1h_2d.zarr"
    process_2d_assignment(input_path, output_path)

    input_path = "/storage/lingyuan2/STATES_data/cellline0101_C4Tg2h"
    output_path = "/storage/lingyuan2/STATES_data/cellline0101_C4Tg2h_2d.zarr"
    process_2d_assignment(input_path, output_path)

    input_path = "/storage/lingyuan2/STATES_data/cellline0101_C5Tg4h"
    output_path = "/storage/lingyuan2/STATES_data/cellline0101_C5Tg4h_2d.zarr"
    process_2d_assignment(input_path, output_path)
