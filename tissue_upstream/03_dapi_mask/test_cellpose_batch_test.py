import numpy as np
import rasterio
from rasterio.features import rasterize
import geopandas as gpd
from shapely.geometry import Polygon
from cellpose import models
import cellpose
import matplotlib.pyplot as plt
import time
import rasterio.features
import torch
from skimage.measure import find_contours
import argparse
import json
import os
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
import time

def percennorm(img,miper=0,maper=100):
    img = np.array(img, dtype='float32')
    datamin =  np.percentile(img, miper, interpolation='midpoint')
    datamax = np.percentile(img, maper, interpolation='midpoint')
    output = (  img - datamin) / (  datamax - datamin)
    output[output > 1] = 1
    output[output < 0] = 0
    return output

def run_cellpose2D(image,diam, model_path,model_type='cyto3', gpudevice='cuda:0', Flowthreshold=0.3,type='cell'):
    
    model = models.CellposeModel(pretrained_model=model_path + '/' +model_type, gpu=True)
    start = time.time()
    image=percennorm(image,0,99.975)
    

    if type == 'cell':
        masks, flows, styles = model.eval(image, diameter=diam, flow_threshold=Flowthreshold)
        end = time.time()
        elapsed_time = end - start
        print('Elapsed time is %f seconds.' % elapsed_time)
    
    elif type == 'nuclei':
        masks, flows, styles = model.eval(image, diameter=diam, flow_threshold=Flowthreshold,resample=False,channels=[[0, 0]])
        end = time.time()
        elapsed_time = end - start
        print('Elapsed time is %f seconds.' % elapsed_time)
    
    return model, masks, flows, styles



def process_cell(cell_mask, i):
    contour = find_contours(cell_mask == i, 0.5)[0]
    polygon = Polygon(contour[:, [1, 0]])
    return polygon


def multiprocess_mask_to_polygon(cell_mask, core_number):
    start_time = time.time()
    with Pool(processes=core_number) as pool:
        # 使用 partial 将 cell_mask 作为参数传递给 process_cell
        func = partial(process_cell, cell_mask)
        cell_polygons = pool.map(func, range(1, cell_mask.max() + 1))
    cell_polygons = gpd.GeoDataFrame(geometry=cell_polygons)
    end_time = time.time()
    print(f"mask to polygon execution time: {end_time - start_time} seconds")
    return cell_polygons


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process images with Cellpose and extract polygons.')
    parser.add_argument('-i', '--image', type=str, required=True, help='Path to the input image file.')
    parser.add_argument('-m', '--methods', type=str, default='2d', help='Method to use, currently supports only "2d".')
    parser.add_argument('-M', '--model', type=str, default='cyto3', help='Model type for Cellpose.')
    parser.add_argument('-d', '--device', type=str, default='cuda:0', help='Device to use for computation.')
    parser.add_argument('-f', '--flow_threshold', type=float, default=0.3, help='Flow threshold for cell segmentation.')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix for output filenames.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Directory to store output files.')
    parser.add_argument('-MD', '--model_path', type=str, required=True, help='Cellpose model path.')
    parser.add_argument('-D', '--diameter', type=float, default='None', help='estimater Diameter of segmented components ,when you do ER or membrane segmentation ,please increase it like 100-200')
    parser.add_argument('-T', '--type', type=str, default='cell', help='input nuclei or cell determines resample')
    parser.add_argument('-R', '--rotate_num', type=int, help='the input tif need to be rotated by K 90 degrees',choices = [0,1,2,3,4])
    parser.add_argument('-C', '--core_num', type=int, help='the core num',required=True)
    args = parser.parse_args()

    print('ready to run')

    img =  cellpose.io.imread(args.image)
    
    img = np.rot90(img,k=args.rotate_num)

    print(img.dtype)
    print(img.shape)
    
    if len(img.shape) > 2:
        img = np.max(img,axis=0)
    else:
        pass

    model, masks, flows, styles = run_cellpose2D(img,diam=args.diameter, model_type=args.model, gpudevice=args.device, Flowthreshold=args.flow_threshold,type=args.type,model_path=args.model_path) 
    

    plt.imsave(os.path.join(args.output, f'{args.prefix}_segmentation.png'), masks)
    cellpose.io.imsave(os.path.join(args.output, f'{args.prefix}_segmentation.tif'), masks)

    gdf = multiprocess_mask_to_polygon(masks,args.core_num)
    gdf.to_file(os.path.join(args.output, f'{args.prefix}_geodataframe.json'), driver="GeoJSON")
    
    # Plot and save max projection with segmentation overlay
    fig, ax = plt.subplots()
    ax.imshow(masks, cmap='gray')
    ax.imshow(masks, alpha=0.5, cmap='jet')
    plt.savefig(os.path.join(args.output, f'{args.prefix}_maxprojection.png'))

    print('read ok')
