import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import spatialdata as sd
from cellpose import models
from shapely.geometry import Polygon
from skimage.io import imread
from skimage.measure import find_contours
from spatialdata.models import PointsModel, ShapesModel
from spatialdata.datasets import blobs
import argparse
import os
import json
import argparse
import os
from rasterio.features import rasterize
import anndata as ad
from spatialdata.datasets import Labels2DModel
from multiprocessing import Pool, cpu_count
from functools import partial
import time 
import dask.array as da

def preprocess_mask(input_file):
    """
    Preprocess the input mask file (either an image (.tif, .png) or a GeoJSON file) and convert it to polygons.

    Parameters:
    input_file (str): Path to the input file, either an image (.tif, .png) or a GeoJSON file.

    Returns:
    gpd.GeoDataFrame: A GeoDataFrame containing the polygons.
    """

    # Check if the input is an image file (.tif, .png)
    if input_file.endswith('.tif') or input_file.endswith('.png'):
        # Load the image (assuming it's a 2D mask image)
        cell_mask = imread(input_file)
        contours = []
        for i in range(1, cell_mask.max()+1):
            contours.append(find_contours(cell_mask == i, 0.5)[0])

        polygons_gdf = np.array([Polygon(p[:, [1, 0]]) for p in contours])
        polygons_gdf = gpd.GeoDataFrame(geometry=polygons_gdf)
        return polygons_gdf

    # Check if the input is a GeoJSON file
    elif input_file.endswith('.json') or input_file.endswith('.geojson'):
        # Use geopandas to read the GeoJSON file
        polygons_gdf = gpd.read_file(input_file)
        
        # Ensure the file contains valid geometries
        if polygons_gdf.empty:
            raise ValueError("No valid polygons found in the GeoJSON file.")
        
        return polygons_gdf

    else:
        raise ValueError(f"Unsupported file format: {input_file}. Please provide a .tif, .png, or .json/.geojson file.")
def process_cell(cell_mask, i):
    contour = find_contours(cell_mask == i, 0.5)[0]
    polygon = Polygon(contour[:, [1, 0]])
    return polygon

def multiprocess_mask_to_polygon(cell_mask, core_number):
    start_time = time.time()
    with Pool(processes=core_number) as pool:
        func = partial(process_cell, cell_mask)
        cell_polygons = pool.map(func, range(1, cell_mask.max() + 1))
    cell_polygons = gpd.GeoDataFrame(geometry=cell_polygons)
    end_time = time.time()
    print(f"mask to polygon execution time: {end_time - start_time} seconds")
    return cell_polygons

def preprocess_mask_multi(input_file,ncores):
    """
    Preprocess the input mask file (either an image (.tif, .png) or a GeoJSON file) and convert it to polygons.

    Parameters:
    input_file (str): Path to the input file, either an image (.tif, .png) or a GeoJSON file.

    
    Returns:
    gpd.GeoDataFrame: A GeoDataFrame containing the polygons.
    
    """
    if input_file.endswith('.tif') or input_file.endswith('.png'):
        # Load the image (assuming it's a 2D mask image)
        
        print("mask_loaded")
        polygons_gdf = multiprocess_mask_to_polygon(cell_mask,ncores)
        return polygons_gdf
        
        # Check if the input is a GeoJSON file
    elif input_file.endswith('.json') or input_file.endswith('.geojson'):
        # Use geopandas to read the GeoJSON file
        polygons_gdf = gpd.read_file(input_file)
        
        polygons_gdf = polygons_gdf[polygons_gdf['geometry'] != None]

        print("mask_loaded")
        # Ensure the file contains valid geometries
        if polygons_gdf.empty:
            raise ValueError("No valid polygons found in the GeoJSON file.")
        
        return polygons_gdf

    else:
        raise ValueError(f"Unsupported file format: {input_file}. Please provide a .tif, .png, or .json/.geojson file.")
        
    

def preprocess_two_ploygons(cell_polygons,nuclei_polygons):
    """
    Preprocess the input 2 polygon mask file (either an image (.tif, .png) or a GeoJSON file) and convert it to polygons.

    Parameters:
    cell_polygons,large polygon preprocessed by previous function 

    nuclei_polygons,small polygon preprocessed by previous function 

    Returns:
    cell_boundary,nuclei_boundary,reindex 1 to 1 match boundary ploygon sets , index start from 1
    """
    # 1.1 filter out no  nuclei boundary
    cell_polygons = cell_polygons[cell_polygons.geometry.apply(
        lambda cell: any(cell.contains(nucleus) for nucleus in nuclei_polygons.geometry)
    )]

    # 1.2 filter out no cell nuclei 
    nuclei_polygons = nuclei_polygons[nuclei_polygons.geometry.apply(
        lambda nucleus: any(cell.contains(nucleus) for cell in cell_polygons.geometry)
    )]
    # reset index 
    cell_polygons = cell_polygons.reset_index(drop=True)


    # 2. reindex nucleus_boundaries
    nucleus_ordered = []
    
    for cell in cell_polygons['geometry']:
        for idx, nucleus in nuclei_polygons['geometry'].items():
            if cell.contains(nucleus):
                nucleus_ordered.append(nucleus)
                break  # assume one cell to one nuclei
    
    # create new GeoDataFrame
    nucleus_ordered_gdf = gpd.GeoDataFrame(geometry=nucleus_ordered)

    # 3. check order
    is_match = True
    
    for i in range(len(cell_polygons)):
        cell_polygon = cell_polygons['geometry'].iloc[i]
        nucleus_polygon = nucleus_ordered_gdf['geometry'].iloc[i]
        
        if not cell_polygon.contains(nucleus_polygon):
            is_match = False
            print(f"Cell boundary at index {i} does not contain the corresponding nucleus boundary.")
            break
    
    if is_match:
        print("All cell boundaries contain the corresponding nucleus boundaries in order.")
    else:
        print("There is a mismatch between cell and nucleus boundaries.")
    # reset index start from 1 to generate mask 
    cell_polygons.index = range(1, len(cell_polygons) + 1)
    nucleus_ordered_gdf.index = range(1, len(nucleus_ordered_gdf) + 1)
    return cell_polygons,nucleus_ordered_gdf

def PolygonToMask(shape,polygon):
    shapes_and_values_boundrary = [(geom, idx) for geom, idx in zip(polygon.geometry,polygon.index)]
    img1=rasterize(shapes_and_values_boundrary,out_shape=shape,dtype=np.uint16)
    return img1
def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Script to process cell and nuclei mask files")

    # Add command-line arguments
    parser.add_argument('-c','--mask_cell', type=str, help="Cell mask file, supports .tif or .json format")
    parser.add_argument('-n','--mask_nuclei', type=str, help="Nuclei mask file, supports .tif or .json format")
    parser.add_argument('-m','--molecular_file', type=str, help="molecular file, supports .csvformat, and must have x,y,Gene column")
    parser.add_argument('-p','--prefix', type=str, help="Prefix for the output file")
    parser.add_argument('-o','--output_path', type=str, help="Path to save the output file")
    parser.add_argument('-t','--threads',default=1,type=int, help="n of threads to process data")
    parser.add_argument('-nt','--mask_nuclei_tif',type=str, help="Nuclei mask tif file,aim to crete shape")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Print parsed arguments
    print(f"Cell mask file: {args.mask_cell}")
    print(f"Nuclei mask file: {args.mask_nuclei}")
    print(f"Prefix: {args.prefix}")
    print(f"Output path: {args.output_path}")
    print(f"thread: {args.threads}")

    # Check if the input files exist
    if not os.path.exists(args.mask_cell):
        raise FileNotFoundError(f"Cell mask file {args.mask_cell} does not exist")

    if not os.path.exists(args.mask_nuclei):
        raise FileNotFoundError(f"Nuclei mask file {args.mask_nuclei} does not exist")

    # Check if the output path exists, create it if it doesn't
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
        print(f"Output path {args.output_path} does not exist, created it")

    #1.parse paramter


    mol =pd.read_csv(args.molecular_file)
    print("molecular_loaded")
    if args.threads > 1:
        print("use multi threads")
        cell=preprocess_mask_multi(args.mask_cell, args.threads)
        nuclei=preprocess_mask_multi(args.mask_nuclei, args.threads)
    else:
        print("use single threads")
        cell = preprocess_mask(args.mask_cell)
        nuclei = preprocess_mask(args.mask_nuclei)


    
    input_shape=imread(args.mask_nuclei_tif).shape

    #rename gene to feature_name
    mol = mol.rename(columns={'gene': 'feature_name'})
    mol = mol.rename({'column':'x','row':'y'},axis=1) 
    #2.preprocessing input
    cell_polygons,nuclei_polygons=preprocess_two_ploygons(cell,nuclei)
    mask1=PolygonToMask(input_shape,cell_polygons)
    
    mask2=PolygonToMask(input_shape,nuclei_polygons)
    
    shapes = dict(cell_boundaries=ShapesModel.parse(cell_polygons),nuclei_boundaries=ShapesModel.parse(nuclei_polygons))
    points = dict(transcripts=PointsModel.parse(mol, coordinates={"x": "x", "y": "y"}, feature_key="feature_name"))
    #3.create sdata raw
    sdata = sd.SpatialData(shapes=shapes,points=points)

    label_1 = Labels2DModel.parse(mask1)
    label_2 = Labels2DModel.parse(mask2)
    sdata.labels["cell_mask"]=label_1
    sdata.labels["nucleus_mask"]=label_2
    x_coords = mol['x'].values
    y_coords = mol['y'].values
    #4.extract cell and nuclei location
    
    x_clipped = np.clip(x_coords, 0, input_shape[0] - 1)
    y_clipped = np.clip(y_coords, 0, input_shape[1] - 1)

    cell_values = sdata.labels['cell_mask'].values[y_clipped-1, x_clipped-1]
    nuclei_values =sdata.labels['nucleus_mask'].values[y_clipped-1, x_clipped-1]
    
    #5.add cell and nuclei info to dataframe
    mol['cell'] = cell_values
    mol['nuclei'] = nuclei_values
    mol_filtered = mol[mol['cell'] != 0]
    count_df = mol_filtered.groupby(['cell', 'feature_name']).size().unstack(fill_value=0)
    
    adata = ad.AnnData(X=count_df.values)
    
    #6.set name and basic of sdata
    adata.obs['cell_id'] = count_df.index.astype(str) 
    adata.var['feature_name'] = count_df.columns.astype(str) 
    adata.obs["region"]="cell_mask"
    adata.uns["points"]=mol
    adata.obs["region"]=adata.obs["region"].astype('category')
    adata.uns["spatialdata_attrs"]={'region': 'cell_mask',
     'region_key': 'region',
     'instance_key': 'cell_id'}
    
    sdata.table=adata
    # Define the output file path
    output_file = os.path.join(args.output_path, f'{args.prefix}_processed.zarr')

    # Write the SpatialData to a Zarr file
    sdata.write(output_file, overwrite=True)

    # Print the result
    print(f"Processing result saved to {output_file}")

if __name__ == "__main__":
    main()
