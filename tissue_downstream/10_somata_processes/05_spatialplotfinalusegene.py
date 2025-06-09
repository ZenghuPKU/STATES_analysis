import numpy as np
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import spatialdata as sd
import scanpy as sc
from shapely.geometry import Polygon, box

adata = sc.read_h5ad('mousebrain_0503.h5ad')
data = adata.uns['remain_reads_info_new']
data['region'] = np.where(data['nuclei'] == -1, 'processes', 'somata')

import geopandas as gpd
import os


genes = ['Cnr1', 'Kcna6', 'Mrpl43', 'Olfm1','Smoc1', 'Spink8', 'Tcf20']

# Define your coordinate filter range
xmin, xmax = 0, 18000
ymin, ymax = 2000, 14000
filter_box = box(xmin, ymin, xmax, ymax)

# Ensure the output directory exists
output_dir = "finalusehippo"
os.makedirs(output_dir, exist_ok=True)

def plot_gene_spatial_region_hippo_filtered(gene_name, ax1, ax2):
    for ax in [ax1, ax2]:
        ax.grid(False)
        ax.axis('off')
        ax.invert_xaxis()

    # Load and filter polygons
    gdf = gpd.read_file('C1.json')
    gdf = gdf[
        (gdf.geometry.is_valid) &
        (~gdf.geometry.is_empty) &
        (gdf.geometry.type == 'Polygon')
    ]
    gdf = gdf[gdf.geometry.apply(lambda poly: len(poly.exterior.coords) >= 3)]
    
    # Clip polygons to the bounding box
    gdf.geometry = gdf.geometry.intersection(filter_box)
    gdf = gdf[~gdf.geometry.is_empty]
    gdf["color"] = "#F0F0F0"

    gdf.plot(ax=ax1, color=gdf["color"], edgecolor="gray", linewidth=0, aspect=1)
    gdf.plot(ax=ax2, color=gdf["color"], edgecolor="gray", linewidth=0, aspect=1)

    # Filter gene data to the bounding box
    gene_data = data[(data['gene'] == gene_name) & (data['replicate'] == 'C1')]
    gene_data = gene_data[
        (gene_data['column'] >= xmin) & (gene_data['column'] <= xmax) &
        (gene_data['row'] >= ymin) & (gene_data['row'] <= ymax)
    ]

    # Plot somata
    somata = gene_data[gene_data['region'] == 'somata']
    rbRNA_somata = somata[somata['type'] == 'rbRNA']
    ntRNA_somata = somata[somata['type'] == 'ntRNA']

    ax1.scatter(rbRNA_somata['column'], rbRNA_somata['row'], c='#09c609', label='rbRNA', alpha=1, s=1)
    ax1.scatter(ntRNA_somata['column'], ntRNA_somata['row'], c='#b3005f', label='ntRNA', alpha=1, s=1)
    ax1.set_title(f'{gene_name} Somata', fontsize=20)
    #ax1.legend(loc='lower right')

    # Plot processes
    processes = gene_data[gene_data['region'] == 'processes']
    rbRNA_processes = processes[processes['type'] == 'rbRNA']
    ntRNA_processes = processes[processes['type'] == 'ntRNA']

    ax2.scatter(rbRNA_processes['column'], rbRNA_processes['row'], c='#09c609', label='rbRNA', alpha=1, s=1)
    ax2.scatter(ntRNA_processes['column'], ntRNA_processes['row'], c='#b3005f', label='ntRNA', alpha=1, s=1)
    ax2.set_title(f'{gene_name} Processes', fontsize=20)
    #ax2.legend(loc='lower right')

for gene in genes:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    plot_gene_spatial_region_hippo_filtered(gene, ax1, ax2)
    pdf_path = os.path.join(output_dir, f"{gene}_spatialplotcolor.pdf")
    fig.savefig(pdf_path, bbox_inches='tight', dpi=300)  # Added dpi for higher quality
    plt.close(fig)
