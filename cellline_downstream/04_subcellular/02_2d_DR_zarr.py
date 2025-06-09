import json
import numpy as np
import pandas as pd
import spatialdata as sd
from skimage.measure import regionprops
from scipy.ndimage import distance_transform_edt as dt

def add_DR_2d(
    original_zarr_path: str,
    new_zarr_path: str,
    table_key: str = "table",
    sampling=(1,1)
):
    print(f"Read SpatialData from: {original_zarr_path}")
    sdata = sd.read_zarr(original_zarr_path)

    cell_mask_da = sdata.labels["2D_cell_mask"]
    nuclei_mask_da = sdata.labels["2D_nucleus_mask"]
    cell_mask_2d = cell_mask_da.values
    nuclei_mask_2d = nuclei_mask_da.values

    print("Calculating regionprops for 2D_cell_mask...")
    rp = regionprops(cell_mask_2d)

    adata = sdata.tables[table_key]
    if "points" not in adata.uns:
        raise KeyError("Key 'points' not found in adata.uns. Please make sure point information is stored in adata.uns['points'].")

    points_json = adata.uns["points"]
    df = pd.DataFrame(json.loads(points_json))

    required_columns = ["x", "y", "cell", "nuclei"]
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' does not exist in points DataFrame, please check.")

    df["DR"] = np.nan
    df["out_of_bounds"] = False

    unique_cells = sorted(df["cell"].unique())
    unique_cells = [c for c in unique_cells if c > 0]
    print(f"Processing {len(unique_cells)} cells...")

    for curr_cell_label in unique_cells:
        idx_rp = curr_cell_label - 1
        if idx_rp < 0 or idx_rp >= len(rp):
            print(f"Skipping cell label {curr_cell_label}: out of regionprops range.")
            continue

        cb = rp[idx_rp].bbox
        cell_mask_local = rp[idx_rp].image

        nuc_mask_local = nuclei_mask_2d[cb[0]:cb[2], cb[1]:cb[3]].copy()
        nuc_mask_local[nuc_mask_local != curr_cell_label] = 0

        cdt = dt(cell_mask_local, sampling=sampling)
        ndt = dt(np.logical_not(nuc_mask_local), sampling=sampling)
        area = np.logical_xor(cell_mask_local, nuc_mask_local)
        rdt = np.zeros_like(cdt, dtype=np.float32)
        with np.errstate(divide='ignore', invalid='ignore'):
            rdt[area] = np.where(cdt[area] != 0, ndt[area] / cdt[area], 0)

        mask_cytosol = (df["cell"] == curr_cell_label) & (df["nuclei"] == 0)
        idxs_this_cell = df.index[mask_cytosol]
        if len(idxs_this_cell) == 0:
            continue

        points_cytosol = df.loc[idxs_this_cell, ["y", "x"]].values

        loc_y = points_cytosol[:, 0] + 1 - cb[0]
        loc_x = points_cytosol[:, 1] + 1 - cb[1]

        within_bounds = (
            (loc_y >= 0) & (loc_y < rdt.shape[0]) &
            (loc_x >= 0) & (loc_x < rdt.shape[1])
        )

        valid_y = loc_y[within_bounds].astype(int)
        valid_x = loc_x[within_bounds].astype(int)

        r_vals = rdt[valid_y, valid_x]
        dr_vals = r_vals / (r_vals + 1)
        dr_vals = np.where((r_vals + 1) != 0, dr_vals, np.nan)

        valid_idxs = idxs_this_cell[within_bounds]
        df.loc[valid_idxs, "DR"] = dr_vals

        out_of_bounds_idxs = idxs_this_cell[~within_bounds]
        df.loc[out_of_bounds_idxs, "out_of_bounds"] = True
        num_out_of_bounds = len(out_of_bounds_idxs)
        if num_out_of_bounds > 0:
            print(f"Cell {curr_cell_label}: {num_out_of_bounds} points are out of bounds and have been marked.")

    updated_points_list = df.to_dict(orient="records")
    adata.uns["points"] = json.dumps(updated_points_list)

    print(f"Writing updated SpatialData to: {new_zarr_path}")
    sdata.write(new_zarr_path, overwrite=True)
    print("Done! The new Zarr file now contains point information with DR column and out-of-bounds flag.")

add_DR_2d(
    original_zarr_path="cellline0101_C3control_2d.zarr",
    new_zarr_path="cellline0101_C3control_2d_dr.zarr",
    table_key="table",
    sampling=(1,1)
)
add_DR_2d(
    original_zarr_path="cellline0101_B4Tg15min_2d.zarr",
    new_zarr_path="cellline0101_B4Tg15min_2d_dr.zarr",
    table_key="table",
    sampling=(1,1)
)
add_DR_2d(
    original_zarr_path="cellline0101_B5Tg30min_2d.zarr",
    new_zarr_path="cellline0101_B5Tg30min_2d_dr.zarr",
    table_key="table",
    sampling=(1,1)
)
add_DR_2d(
    original_zarr_path="cellline0101_B6Tg1h_2d.zarr",
    new_zarr_path="cellline0101_B6Tg1h_2d_dr.zarr",
    table_key="table",
    sampling=(1,1)
)
add_DR_2d(
    original_zarr_path="cellline0101_C4Tg2h_2d.zarr",
    new_zarr_path="cellline0101_C4Tg2h_2d_dr.zarr",
    table_key="table",
    sampling=(1,1)
)
add_DR_2d(
    original_zarr_path="cellline0101_C5Tg4h_2d.zarr",
    new_zarr_path="cellline0101_C5Tg4h_2d_dr.zarr",
    table_key="table",
    sampling=(1,1)
)