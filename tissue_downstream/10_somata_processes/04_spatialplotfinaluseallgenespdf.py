import numpy as np
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import spatialdata as sd
import scanpy as sc
from shapely.geometry import box

adata = sc.read_h5ad('mousebrain_0503.h5ad')
data = adata.uns['remain_reads_info_new']
data['region'] = np.where(data['nuclei'] == -1, 'processes', 'somata')

import geopandas as gpd
import os

# Define your coordinate filter range
xmin, xmax = 8000, 12500
ymin, ymax = 11500, 14000
filter_box = box(xmin, ymin, xmax, ymax)

# Ensure the output directory exists
output_dir = "finaluse"
os.makedirs(output_dir, exist_ok=True)

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

# Filter all gene data to the bounding box and replicate C1
all_gene_data = data[(data['replicate'] == 'C1')]
all_gene_data = all_gene_data[
    (all_gene_data['column'] >= xmin) & (all_gene_data['column'] <= xmax) &
    (all_gene_data['row'] >= ymin) & (all_gene_data['row'] <= ymax)
]

# Plot somata
fig_somata, ax_somata = plt.subplots(figsize=(12, 12))
ax_somata.grid(False)
ax_somata.axis('off')
ax_somata.invert_xaxis()
gdf.plot(ax=ax_somata, facecolor="#F0F0F0", edgecolor="black", linewidth=0.5, aspect=1)

somata = all_gene_data[all_gene_data['region'] == 'somata']
rbRNA_somata = somata[somata['type'] == 'rbRNA']
ntRNA_somata = somata[somata['type'] == 'ntRNA']

ax_somata.scatter(rbRNA_somata['column'], rbRNA_somata['row'], c='#00913a', alpha=1, s=1, marker='o', edgecolors='none')
ax_somata.scatter(ntRNA_somata['column'], ntRNA_somata['row'], c='#e83828', alpha=1, s=1, marker='o', edgecolors='none')

somata_pdf_path = os.path.join(output_dir, "all_genes_somata_color.pdf")
fig_somata.savefig(somata_pdf_path, bbox_inches='tight', dpi=300)
plt.close(fig_somata)

# Plot processes
fig_processes, ax_processes = plt.subplots(figsize=(12, 12))
ax_processes.grid(False)
ax_processes.axis('off')
ax_processes.invert_xaxis()
gdf.plot(ax=ax_processes, facecolor="#F0F0F0", edgecolor="black", linewidth=0.5, aspect=1)

processes = all_gene_data[all_gene_data['region'] == 'processes']
rbRNA_processes = processes[processes['type'] == 'rbRNA']
ntRNA_processes = processes[processes['type'] == 'ntRNA']

ax_processes.scatter(rbRNA_processes['column'], rbRNA_processes['row'], c='#00913a', alpha=1, s=0.5, marker='o', edgecolors='none')
ax_processes.scatter(ntRNA_processes['column'], ntRNA_processes['row'], c='#e83828', alpha=1, s=0.5, marker='o', edgecolors='none')

processes_pdf_path = os.path.join(output_dir, "all_genes_processes_color.pdf")
fig_processes.savefig(processes_pdf_path, bbox_inches='tight', dpi=300)
plt.close(fig_processes)
