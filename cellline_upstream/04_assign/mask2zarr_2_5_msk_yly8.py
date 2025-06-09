#!/usr/bin/env python3
import pandas as pd
import os
import pyproj
if not os.environ.get('PROJ_LIB'):
    proj_data_dir = "/appsnew/home/huzeng_pkuhpc/mambaforge/envs/sopa/share/proj"
    os.environ['PROJ_LIB'] = proj_data_dir
    pyproj.datadir.set_data_dir(proj_data_dir)
import argparse
import shutil
import json
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, Point
from skimage.io import imread
from skimage.measure import label, find_contours, regionprops
import matplotlib.pyplot as plt
import seaborn as sns
import spatialdata as sd
import anndata as ad
from spatialdata.models import PointsModel, ShapesModel, Labels2DModel, Labels3DModel
from rasterio.features import rasterize
from functools import partial, lru_cache
from multiprocessing import Pool
import time
import multiprocessing
import psutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc
import logging
from collections import defaultdict
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


class Labels3DModel(Labels2DModel):
    class dims:
        dims = ("z", "y", "x")

    @classmethod
    def parse(cls, data: np.ndarray):
        if data.ndim != 3:
            raise ValueError(f"Expected 3D array, got shape={data.shape}")
        return super().parse(data)


def preprocess_mask(input_file):
    """
    Preprocess the input mask file (either an image or a GeoJSON file) and convert it to polygons.
    """
    if input_file.endswith('.tif') or input_file.endswith('.png'):
        print(f"\nProcessing mask file: {input_file}")
        cell_mask = imread(input_file)
        print(f"Found {cell_mask.max()} objects")
        
        contours = []
        for i in tqdm(range(1, cell_mask.max()+1), desc="Finding contours"):
            contours.append(find_contours(cell_mask == i, 0.5)[0])

        print("Converting contours to polygons...")
        polygons = [Polygon(p[:, [1, 0]]) for p in tqdm(contours, desc="Creating polygons")]
        polygons_gdf = gpd.GeoDataFrame(geometry=polygons)
        return polygons_gdf

    elif input_file.endswith('.json') or input_file.endswith('.geojson'):
        polygons_gdf = gpd.read_file(input_file)
        if polygons_gdf.empty:
            raise ValueError("No valid polygons found in the GeoJSON file.")
        return polygons_gdf

    else:
        raise ValueError(f"Unsupported file format: {input_file}. Please provide a .tif, .png, or .json/.geojson file.")


def mask_to_polygons_2d(mask_2d):
    """
    Convert a 2D mask to polygons (boundaries).
    Returns a GeoDataFrame containing 'geometry' and 'label' columns.
    """
    geoms = []
    labels_ = []
    unique_lbs = np.unique(mask_2d[mask_2d > 0])
    for lb in unique_lbs:
        c_list = find_contours(mask_2d == lb, 0.5)
        if not c_list:
            continue
        c_big = max(c_list, key=lambda x: len(x))
        coords_xy = [(pt[1], pt[0]) for pt in c_big]
        poly = Polygon(coords_xy)
        if not poly.is_valid:
            poly = poly.buffer(0)
            if not poly.is_valid:
                continue
        geoms.append(poly)
        labels_.append(lb)
    gdf = gpd.GeoDataFrame({'geometry': geoms, 'label': labels_})
    return gdf

def preprocess_and_match_cells_nuclei(cell_gdf: gpd.GeoDataFrame,
                                      nuc_gdf: gpd.GeoDataFrame,
                                      overlap_threshold=0.2):
    """
    Preprocess and match cell and nucleus polygon data.

    Matching logic:
    1. If cell or nucleus data is empty, return empty result
    2. For each cell:
       a. Iterate through all nuclei:
          - If cell completely contains nucleus, add to candidate list
          - Otherwise calculate intersection area and overlap ratio 
            (intersection area / cell area)
          - If overlap ratio > threshold, add to candidate list
       b. Merge candidate lists and select largest nucleus as best match
       c. Update nucleus geometry to keep only the overlapping region

    Parameters:
    - cell_gdf: GeoDataFrame containing cell polygons
    - nuc_gdf: GeoDataFrame containing nucleus polygons  
    - overlap_threshold: Minimum overlap ratio required for matching

    Returns:
    - matched_cells_gdf: GeoDataFrame of matched cells
    - matched_nuclei_gdf: GeoDataFrame of matched nuclei with updated geometries
    """
    import geopandas as gpd
    from tqdm import tqdm

    # Return empty result if input data is empty
    if cell_gdf.empty or nuc_gdf.empty:
        print("Cell or nucleus data is empty, returning empty result.")
        return cell_gdf, nuc_gdf

    print("\nStarting cell-nucleus matching...")
    print(f"Initial counts - Cells: {len(cell_gdf)}, Nuclei: {len(nuc_gdf)}")

    matched_cells = []
    matched_nuclei = []

    # Iterate through each cell
    for idx, cell in tqdm(cell_gdf.iterrows(), total=len(cell_gdf), desc="Matching cells"):
        cell_geom = cell.geometry
        cell_area = cell_geom.area

        candidate_contained = []  # Fully contained nuclei
        candidate_overlap = []    # Nuclei meeting overlap threshold

        # Iterate through all nuclei
        for jdx, nucleus in nuc_gdf.iterrows():
            nuc_geom = nucleus.geometry
            nuc_area = nuc_geom.area
            
            # Skip if nucleus area is larger than cell area
            if nuc_area >= cell_area:
                continue

            # Check if cell completely contains nucleus
            if cell_geom.contains(nuc_geom):
                candidate_contained.append(nucleus)
            else:
                # Calculate intersection and overlap ratio
                intersection = cell_geom.intersection(nuc_geom)
                # Ensure intersection is a polygon
                if not intersection.is_empty:
                    if intersection.geom_type == 'GeometryCollection':
                        # Extract largest polygon from GeometryCollection
                        polygons = [geom for geom in intersection.geoms if geom.geom_type == 'Polygon']
                        if polygons:
                            intersection = max(polygons, key=lambda x: x.area)
                        else:
                            continue
                    elif intersection.geom_type != 'Polygon':
                        continue
                    
                    intersection_area = intersection.area
                    overlap_ratio = intersection_area / cell_area if cell_area > 0 else 0
                    if overlap_ratio > 0:
                        print(f"Cell {idx} and nucleus {jdx} have an overlap ratio of {overlap_ratio:.3f}")
                    if overlap_ratio >= overlap_threshold:
                        candidate_overlap.append(nucleus)

        # Merge candidate lists
        candidates = candidate_contained + candidate_overlap
        candidates = [nucleus for nucleus in candidates if nucleus.geometry.area >= 100]

        best_nucleus = None
        if candidates:
            best_nucleus = max(candidates, key=lambda n: n.geometry.area)

        # Update nucleus geometry with intersection
        if best_nucleus is not None:
            new_nucleus_geom = cell_geom.intersection(best_nucleus.geometry)
            # Ensure new geometry is a polygon
            if not new_nucleus_geom.is_empty:
                if new_nucleus_geom.geom_type == 'GeometryCollection':
                    polygons = [geom for geom in new_nucleus_geom.geoms if geom.geom_type == 'Polygon']
                    if polygons:
                        new_nucleus_geom = max(polygons, key=lambda x: x.area)
                    else:
                        continue
                elif new_nucleus_geom.geom_type != 'Polygon':
                    continue
                
                best_nucleus = best_nucleus.copy()
                best_nucleus.geometry = new_nucleus_geom
                matched_cells.append(cell)
                matched_nuclei.append(best_nucleus)

    # Create GeoDataFrames from matches
    matched_cells_gdf = gpd.GeoDataFrame(matched_cells, columns=cell_gdf.columns)
    matched_nuclei_gdf = gpd.GeoDataFrame(matched_nuclei, columns=nuc_gdf.columns)

    print(f"Final match count: {len(matched_cells_gdf)}")
    
    print("\nVerifying matches for mismatches...")
    mismatches = 0
    for i in tqdm(range(len(matched_cells_gdf)), desc="Verifying"):
        if not matched_cells_gdf['geometry'].iloc[i].contains(matched_nuclei_gdf['geometry'].iloc[i]):
            mismatches += 1

    if mismatches > 0:
        print(f"Warning: Found {mismatches} mismatches")
    else:
        print("All matches verified successfully")

    print(f"Final counts - Matched pairs: {len(matched_cells_gdf)}")
    
    return matched_cells_gdf.reset_index(drop=True), matched_nuclei_gdf.reset_index(drop=True)

def preprocess_two_polygons(cell_polygons, nuclei_polygons):
    """
    Preprocess the input two polygon mask files and convert them to polygons.
    Returns filtered and reindexed cell and nucleus GeoDataFrames.
    """
    cell_polygons = cell_polygons[cell_polygons.geometry.apply(
        lambda cell: any(cell.contains(nucleus) for nucleus in nuclei_polygons.geometry)
    )]
    nuclei_polygons = nuclei_polygons[nuclei_polygons.geometry.apply(
        lambda nucleus: any(cell.contains(nucleus) for cell in cell_polygons.geometry)
    )]
    cell_polygons = cell_polygons.reset_index(drop=True)
    nucleus_ordered = []
    for cell in cell_polygons['geometry']:
        for idx, nucleus in nuclei_polygons['geometry'].items():
            if cell.contains(nucleus):
                nucleus_ordered.append(nucleus)
                break
    nucleus_ordered_gdf = gpd.GeoDataFrame(geometry=nucleus_ordered)
    is_match = True
    for i in range(len(cell_polygons)):
        if not cell_polygons['geometry'].iloc[i].contains(nucleus_ordered_gdf['geometry'].iloc[i]):
            is_match = False
            print(f"Cell boundary at index {i} does not contain the corresponding nucleus boundary.")
            break
    if not is_match:
        print("There is a mismatch between cell and nucleus boundaries.")
    cell_polygons.index = range(1, len(cell_polygons) + 1)
    nucleus_ordered_gdf.index = range(1, len(nucleus_ordered_gdf) + 1)
    return cell_polygons, nucleus_ordered_gdf


def chunk_array(array, chunk_size=6000):
    Z, H, W = array.shape
    total_pixels = H * W
    optimal_chunks = int(np.sqrt(total_pixels / 1e6))
    chunk_size = max(min(chunk_size, H // optimal_chunks, W // optimal_chunks), 1000)
    print(f"Chunk size: {chunk_size}x{chunk_size}")
    total_chunks = ((H + chunk_size - 1) // chunk_size) * ((W + chunk_size - 1) // chunk_size)
    print(f"Total chunks per slice: {total_chunks}")
    for z in range(Z):
        for i in range(0, H, chunk_size):
            for j in range(0, W, chunk_size):
                h_end = min(i + chunk_size, H)
                w_end = min(j + chunk_size, W)
                yield z, (i, h_end), (j, w_end), array[z:z+1, i:h_end, j:w_end].copy()


@lru_cache(maxsize=1024)
def get_cell_at_point(x, y, cell_bounds):
    point = (x, y)
    for cell_id, bounds in cell_bounds.items():
        if bounds[0] <= x <= bounds[2] and bounds[1] <= y <= bounds[3]:
            return cell_id
    return None


def preprocess_cell_boundaries(cell_polygons):
    cell_bounds = {}
    cell_pixels = defaultdict(set)
    for idx, row in cell_polygons.iterrows():
        bounds = row['geometry'].bounds
        cell_bounds[idx] = bounds
        boundary = list(row['geometry'].exterior.coords)
        for x, y in boundary:
            cell_pixels[idx].add((int(x), int(y)))
    return cell_bounds, cell_pixels


def optimize_chunk_processing(chunk, y_start, x_start, cell_bounds, cell_pixels):
    result = np.zeros_like(chunk, dtype=np.uint32)
    mask = chunk > 0
    if not np.any(mask):
        return result
    y_coords, x_coords = np.where(mask)
    if len(y_coords) == 0:
        return result
    global_y = y_coords + y_start
    global_x = x_coords + x_start
    for i in range(len(y_coords)):
        x, y = global_x[i], global_y[i]
        cell_id = get_cell_at_point(x, y, tuple(cell_bounds.items()))
        if cell_id is not None:
            result[y_coords[i], x_coords[i]] = cell_id
    return result


def process_chunk(args):
    z, (y_start, y_end), (x_start, x_end), chunk, cell_polygons, cell_bounds, cell_pixels = args
    if chunk.ndim == 3:
        chunk = chunk[0]
    result = optimize_chunk_processing(chunk, y_start, x_start, cell_bounds, cell_pixels)
    return result, z, (y_start, y_end), (x_start, x_end)


def process_chunk_mask(args):
    z, (y_start, y_end), (x_start, x_end), chunk, ref_mask = args
    if chunk.ndim == 3:
        chunk = chunk[0]
    result = np.zeros_like(chunk, dtype=np.uint16)
    mask = chunk > 0
    if np.any(mask):
        result[mask] = ref_mask[y_start:y_end, x_start:x_end][mask]
    return result, z, (y_start, y_end), (x_start, x_end)


###############################################
###############################################
def process_3d_mask(mask_3d_img, ref_mask, chunk_size, max_workers, total_memory, checkpoint_prefix, npy_dir=None):
    start_time = time.time()
    Z, H, W = mask_3d_img.shape
    print(f"\nProcessing {checkpoint_prefix}:")
    print(f"Shape: {mask_3d_img.shape}")
    print(f"Input dtype: {mask_3d_img.dtype}")
    print(f"Memory usage: {mask_3d_img.nbytes / 1024**3:.2f} GB")
    mask_3d_final = np.zeros((Z, H, W), dtype=np.uint16)
    if npy_dir is not None:
        print(f"Loading preprocessed npy files from {npy_dir} for {checkpoint_prefix}")
        for z in range(Z):
            file_path = os.path.join(npy_dir, f"{checkpoint_prefix}_final_slice_{z}.npy")
            if os.path.exists(file_path):
                mask_3d_final[z] = np.load(file_path)
                print(f"Loaded slice {z} from {file_path}")
            else:
                print(f"Warning: {file_path} not found. Processing slice {z} normally.")
                chunk_args = []
                for _, (y_start, y_end), (x_start, x_end), chunk in chunk_array(mask_3d_img[z:z+1], chunk_size):
                    chunk_args.append((z, (y_start, y_end), (x_start, x_end), chunk, ref_mask))
                for args_item in chunk_args:
                    result, z_idx, (y_start, y_end), (x_start, x_end) = process_chunk_mask(args_item)
                    mask_3d_final[z_idx, y_start:y_end, x_start:x_end] = result
                np.save(file_path, mask_3d_final[z])
        elapsed = time.time() - start_time
        speed = Z / elapsed
        print(f"\n{checkpoint_prefix} loading completed:")
        print(f"Total time: {elapsed:.2f}s")
        print(f"Speed: {speed:.2f} slices/second")
        print(f"Final shape: {mask_3d_final.shape}")
        return mask_3d_final

    checkpoint_file = f'{checkpoint_prefix}_checkpoint.txt'
    start_z = 0
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as f:
            start_z = int(f.read().strip()) + 1
            print(f"Resuming {checkpoint_prefix} from slice {start_z}/{Z}")
    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            with tqdm(total=Z-start_z, desc=f"Processing {checkpoint_prefix}") as pbar:
                for z in range(start_z, Z):
                    mem = psutil.virtual_memory()
                    if mem.available < total_memory * 0.2:
                        print(f"\nLow memory warning: {mem.available/1024**3:.2f}GB available")
                        print("Waiting for memory cleanup...")
                        gc.collect()
                        time.sleep(30)
                    chunk_args = []
                    chunks_count = 0
                    for _, (y_start, y_end), (x_start, x_end), chunk in chunk_array(mask_3d_img[z:z+1], chunk_size):
                        chunk_args.append((z, (y_start, y_end), (x_start, x_end), chunk, ref_mask))
                        chunks_count += 1
                    print(f"\nProcessing slice {z+1}/{Z} with {chunks_count} chunks")
                    futures = [executor.submit(process_chunk_mask, args_item) for args_item in chunk_args]
                    with tqdm(total=len(futures), desc=f"Chunks in slice {z+1}") as chunk_pbar:
                        for future in as_completed(futures):
                            result, z_idx, (y_start, y_end), (x_start, x_end) = future.result()
                            mask_3d_final[z_idx, y_start:y_end, x_start:x_end] = result
                            chunk_pbar.update(1)
                    with open(checkpoint_file, 'w') as f:
                        f.write(str(z))
                    np.save(f'{checkpoint_prefix}_final_slice_{z}.npy', mask_3d_final[z])
                    pbar.update(1)
                    mem = psutil.virtual_memory()
                    print(f"Memory usage: {(mem.total-mem.available)/1024**3:.2f}GB / {mem.total/1024**3:.2f}GB")
    except Exception as e:
        print(f"\nError in {checkpoint_prefix} processing: {str(e)}")
        np.save(f'{checkpoint_prefix}_final_partial.npy', mask_3d_final)
        raise e
    finally:
        if os.path.exists(checkpoint_file):
            os.remove(checkpoint_file)
    elapsed = time.time() - start_time
    speed = Z / elapsed
    print(f"\n{checkpoint_prefix} processing completed:")
    print(f"Total time: {elapsed:.2f}s")
    print(f"Speed: {speed:.2f} slices/second")
    print(f"Final shape: {mask_3d_final.shape}")
    return mask_3d_final


def label_3d_masks_per_slice(cell_3d_img, nuc_3d_img, label_ref_gdf, npy_dir=None):
    total_memory = psutil.virtual_memory().total
    available_memory = psutil.virtual_memory().available
    chunk_size = min(8000, int(np.sqrt((available_memory * 0.3) / (2 * 2))))
    max_workers = min(
        16,
        int(available_memory / (2 * 1024 * 1024 * 1024))
    )
    print(f"System memory status:")
    print(f"Total memory: {total_memory / (1024**3):.1f} GB")
    print(f"Available memory: {available_memory / (1024**3):.1f} GB")
    print(f"Using chunk_size: {chunk_size}")
    print(f"Using {max_workers} workers")
    print(f"Input dtypes: cell={cell_3d_img.dtype}, nuc={nuc_3d_img.dtype}")

    input_shape = cell_3d_img.shape[1:]
    if 'label_ref' in label_ref_gdf.columns:
        cell_polygons = label_ref_gdf[['geometry', 'label_ref']].copy()
        cell_polygons = cell_polygons.rename(columns={'label_ref': 'cell_label'})
    elif 'label' in label_ref_gdf.columns:
        cell_polygons = label_ref_gdf[['geometry', 'label']].copy()
        cell_polygons = cell_polygons.rename(columns={'label': 'cell_label'})
    else:
        label_ref_gdf = label_ref_gdf.copy()
        label_ref_gdf['label_ref'] = label_ref_gdf.index
        cell_polygons = label_ref_gdf[['geometry', 'label_ref']].copy()
        cell_polygons = cell_polygons.rename(columns={'label_ref': 'cell_label'})
    ref_mask = PolygonToMask(input_shape, cell_polygons, dtype=np.uint16)
    cell_3d_final = process_3d_mask(cell_3d_img, ref_mask, chunk_size, max_workers, total_memory, 'cell_3d', npy_dir)
    nuc_3d_final  = process_3d_mask(nuc_3d_img, ref_mask, chunk_size, max_workers, total_memory, 'nuc_3d', npy_dir)
    return cell_3d_final, nuc_3d_final


def assign_labels_to_points(cell_mask, nuclei_mask, mol):
    """
    Assign labels to molecular points based on the provided masks.
    """
    if cell_mask.ndim == 2:
        H, W = cell_mask.shape
        x_int = np.clip(mol['x'].values - 1, 0, W - 1).astype(int)
        y_int = np.clip(mol['y'].values - 1, 0, H - 1).astype(int)
        c_vals = cell_mask[y_int, x_int]
        n_vals = nuclei_mask[y_int, x_int]
    else:
        Z, H, W = cell_mask.shape
        x_int = np.clip(mol['x'].values - 1, 0, W - 1).astype(int)
        y_int = np.clip(mol['y'].values - 1, 0, H - 1).astype(int)
        z_int = np.clip(mol['z'].values - 1, 0, Z - 1).astype(int)
        c_vals = cell_mask[z_int, y_int, x_int]
        n_vals = nuclei_mask[z_int, y_int, x_int]
    mol['cell'] = c_vals
    mol['nuclei'] = n_vals
    conditions = [
        mol['nuclei'] > 0,
        mol['cell'] > 0
    ]
    choices = ['Nuclei', 'Cell']
    mol['mask_label'] = np.select(conditions, choices, default='None')
    return mol


def multiprocess_mask_to_polygon(cell_mask, core_number):
    """
    Convert mask to polygons using multiprocessing.
    """
    start_time = time.time()
    with Pool(processes=core_number) as pool:
        func = partial(process_cell, cell_mask)
        cell_polygons = pool.map(func, range(1, cell_mask.max() + 1))
    cell_polygons = [poly for poly in cell_polygons if poly is not None]
    polygons_gdf = gpd.GeoDataFrame({'geometry': cell_polygons})
    end_time = time.time()
    return polygons_gdf


def preprocess_mask_multi(input_file, ncores):
    """
    Preprocess the input mask file (supports multiprocessing) and convert it to polygons.
    """
    if input_file.endswith('.tif') or input_file.endswith('.png'):
        cell_mask = imread(input_file)
        polygons_gdf = multiprocess_mask_to_polygon(cell_mask, ncores)
        return polygons_gdf
    elif input_file.endswith('.json') or input_file.endswith('.geojson'):
        polygons_gdf = gpd.read_file(input_file)
        if polygons_gdf.empty:
            raise ValueError("No valid polygons found in the GeoJSON file.")
        return polygons_gdf
    else:
        raise ValueError(f"Unsupported file format: {input_file}. Please provide a .tif, .png, or .json/.geojson file.")


def PolygonToMask(shape, polygon, dtype=np.uint16):
    """
    Convert polygons to mask.
    """
    shapes_and_values_boundary = [(geom, idx) for geom, idx in zip(polygon.geometry, polygon.index)]
    img = rasterize(shapes_and_values_boundary, out_shape=shape, dtype=dtype)
    return img


def process_cell(cell_mask, i):
    """
    Process a single cell mask and convert it to a polygon.
    """
    contour_list = find_contours(cell_mask == i, 0.5)
    if not contour_list:
        return None
    contour = max(contour_list, key=lambda x: len(x))
    if len(contour) < 3:
        return None
    polygon = Polygon(contour[:, [1, 0]])
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    if polygon.is_valid:
        return polygon
    return None


def main():
    parser = argparse.ArgumentParser(description="""
    Logic:
      (1) 2D cell/nuclei => polygons => match => label_ref => visualize
      (2) 3D cell/nuclei => per-slice labeling based on 2D overlap (max intersection) => visualize => final
      (3) molecules => .zarr
    """)
    parser.add_argument("--cell_2d", required=True, help="2D cell mask tif")
    parser.add_argument("--nuc_2d", required=True, help="2D nuclei mask tif")
    parser.add_argument("--cell_3d", required=True, help="2.5D or 3D cell mask tif")
    parser.add_argument("--nuc_3d", required=True, help="2.5D or 3D nuclei mask tif")
    parser.add_argument("--molecular_file", required=False, help="CSV with x,y,z,Gene")
    parser.add_argument("--prefix", default="sample")
    parser.add_argument("--output_path", default=".")
    parser.add_argument("--threads", default=1, type=int, help="Number of threads for processing")
    parser.add_argument("--json_cell", default=None, help="Path to the precomputed cell polygon JSON file (default: None)")
    parser.add_argument("--json_nuc", default=None, help="Path to the precomputed nuclei polygon JSON file (default: None)")
    parser.add_argument("--npy_dir", default=None, help="Directory containing preprocessed npy files for 3D masks (if available)")
    args = parser.parse_args()

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path, exist_ok=True)

    print("\n=== Starting Processing ===")
    
    print("\n1. Loading precomputed 2D polygon JSON files...")
    if args.threads > 1:
        print(f"Using {args.threads} threads for parallel processing")
        with tqdm(total=2, desc="Loading 2D polygon JSON") as pbar:
            cell_gdf = gpd.read_file(args.json_cell)
            pbar.update(1)
            nuc_gdf  = gpd.read_file(args.json_nuc)
            pbar.update(1)
    else:
        print("Using single thread processing")
        with tqdm(total=2, desc="Loading 2D polygon JSON") as pbar:
            cell_gdf = gpd.read_file(args.json_cell)
            pbar.update(1)
            nuc_gdf  = gpd.read_file(args.json_nuc)
            pbar.update(1)

    print("\n2. Matching cells and nuclei...")
    #matched_cell_gdf = cell_gdf.copy()
    #matched_nuc_gdf = nuc_gdf.copy()
    matched_cell_gdf, matched_nuc_gdf = preprocess_and_match_cells_nuclei(cell_gdf, nuc_gdf)

    print("\n3. Processing 3D data...")
    print("Reading 3D images...")
    with tqdm(total=2, desc="Loading 3D images") as pbar:
        cell_3d_img = imread(args.cell_3d)
        pbar.update(1)
        nuc_3d_img = imread(args.nuc_3d)
        pbar.update(1)
    
    print("\n4. Processing molecular data...")
    if args.molecular_file and os.path.exists(args.molecular_file):
        mol_df = pd.read_csv(args.molecular_file)
        print(f"Loaded {len(mol_df)} molecular points")
        if "z" not in mol_df.columns:
            mol_df["z"] = 1
        if "Gene" in mol_df.columns:
            mol_df.rename(columns={"Gene": "feature_name"}, inplace=True)
        mol_present = True
    else:
        print("No molecular data provided or file not found")
        mol_df = pd.DataFrame()
        mol_present = False
    
    print("\n5. Creating SpatialData object...")
    with tqdm(total=5, desc="Building final dataset") as pbar:
        label_ref_gdf = matched_cell_gdf.copy()
        print("\nSystem memory status before processing:")
        print(f"Total memory: {psutil.virtual_memory().total / (1024**3):.1f} GB")
        print(f"Available memory: {psutil.virtual_memory().available / (1024**3):.1f} GB")
        
        cell_3d_final, nuc_3d_final = label_3d_masks_per_slice(
            cell_3d_img, nuc_3d_img, label_ref_gdf, npy_dir=args.npy_dir
        )
        
        if not mol_present:
            mol_df = pd.DataFrame()
        else:
            mol_df = assign_labels_to_points(cell_3d_final, nuc_3d_final, mol_df)
        
        input_shape = imread(args.cell_2d).shape
        mask1 = PolygonToMask(input_shape, matched_cell_gdf)
        mask2 = PolygonToMask(input_shape, matched_nuc_gdf)
        label_2d_cell_mask = Labels2DModel.parse(mask1)
        label_2d_nucleus_mask = Labels2DModel.parse(mask2)
        labels_dict = {}
        labels_dict["2D_cell_mask"] = label_2d_cell_mask
        labels_dict["2D_nucleus_mask"] = label_2d_nucleus_mask

        if cell_3d_final.ndim == 3:
            label_model_3d = Labels3DModel
        else:
            label_model_3d = Labels2DModel

        labels_dict["3D_cell_mask"] = label_model_3d.parse(cell_3d_final)
        labels_dict["3D_nuclei_mask"] = label_model_3d.parse(nuc_3d_final)
        
        points_dict = {}
        if not mol_df.empty:
            points_dict["transcripts"] = PointsModel.parse(
                mol_df,
                coordinates={"x": "x", "y": "y"},
                feature_key="feature_name" if "feature_name" in mol_df.columns else "Gene"
            )

        sdata = sd.SpatialData(labels=labels_dict, points=points_dict)
        """
        if not cell_gdf.empty and not nuc_gdf.empty:
            print("\nCell GDF columns:", cell_gdf.columns.tolist())
            print("Nuclei GDF columns:", nuc_gdf.columns.tolist())
            
            if not any(col in cell_gdf.columns for col in ['label', 'label_ref']):
                print("Creating sequential index for cell_gdf")
                cell_gdf = cell_gdf.copy()
                cell_gdf.index = range(1, len(cell_gdf) + 1)
                
            if not any(col in nuc_gdf.columns for col in ['label', 'label_ref']):
                print("Creating sequential index for nuc_gdf")
                nuc_gdf = nuc_gdf.copy()
                nuc_gdf.index = range(1, len(nuc_gdf) + 1)
            
            shapes = {
                'cell_boundaries': ShapesModel.parse(cell_gdf),
                'nucleus_boundaries': ShapesModel.parse(nuc_gdf)
            }
            sdata.shapes = shapes
        else:
            print("[WARNING] No cell or nucleus boundaries to add to shapes.")
        """
        if not matched_cell_gdf.empty and not matched_nuc_gdf.empty:
            if "label" not in matched_cell_gdf.columns:
                matched_cell_gdf["label"] = matched_cell_gdf.index
            if "label" not in matched_nuc_gdf.columns:
                matched_nuc_gdf["label"] = matched_nuc_gdf.index

            matched_cell_gdf = matched_cell_gdf.set_index('label')
            matched_nuc_gdf = matched_nuc_gdf.set_index('label')
            shapes = {
                'cell_boundaries': ShapesModel.parse(matched_cell_gdf),
                'nucleus_boundaries': ShapesModel.parse(matched_nuc_gdf)
            }
            sdata.shapes = shapes
        else:
            print("[WARNING] No cell or nucleus boundaries to add to shapes.")

        if not mol_df.empty and "cell" in mol_df.columns:
            mol_f = mol_df[mol_df["cell"] > 0]
            if not mol_f.empty and "feature_name" in mol_df.columns:
                count_df = mol_f.groupby(["cell", "feature_name"]).size().unstack(fill_value=0)
                adata_table = ad.AnnData(X=count_df.values)
                adata_table.obs["cell_id"] = count_df.index.astype(str)
                adata_table.var["feature_name"] = count_df.columns.astype(str)
            else:
                adata_table = ad.AnnData()
        else:
            adata_table = ad.AnnData()

        if not adata_table.obs.empty:
            adata_table.obs["region"] = "cell_mask"
            adata_table.obs["region"] = adata_table.obs["region"].astype('category')

        if not mol_df.empty:
            mol_list = mol_df.to_dict(orient="records")
            mol_json = json.dumps(mol_list, ensure_ascii=False)
            adata_table.uns["points"] = mol_json

        adata_table.uns["spatialdata_attrs"] = {
            'region': 'cell_mask',
            'region_key': 'region',
            'instance_key': 'cell_id'
        }
        sdata.tables["table"] = adata_table
        out_zarr = os.path.join(args.output_path, f"{args.prefix}_processed.zarr")
        if os.path.exists(out_zarr):
            shutil.rmtree(out_zarr)
        sdata.write(out_zarr, overwrite=True)
        print("[INFO] Wrote sdata =>", out_zarr)
        pbar.update(1)
    print("\n=== Processing Complete ===")


if __name__ == "__main__":
    main()
