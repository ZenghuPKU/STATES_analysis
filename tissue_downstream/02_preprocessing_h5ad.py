#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running example:
nohup python 02_preprocessing_h5ad.py \
    --input_path /h5ad \
    --output_path /output\
    --rep1_file B4_WT_mousebrain.h5ad \
    --rep2_file C1_WT_mousebrain.h5ad \
    --qc_n 3 \
    --gene_threshold 3 \
    --min_cells 10 \
    --min_genes 10\
    --final_fn mousebrain_harmony.h5ad \
    > preprocessing_h5ad.out 2>&1 &

"""

import os
# Set thread environment variables to limit multi-threading for reproducibility
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import scanpy.external as sce
import harmonypy
import warnings
import matplotlib

matplotlib.use('Agg')
warnings.filterwarnings('ignore')


def load_h5ad_files(input_path, rep1_file, rep2_file):
    """
    Load two h5ad files from the given directory, preserve the original barcodes,
    and add replicate information to the 'remain_reads_info' if available.
    """
    adata1 = sc.read_h5ad(os.path.join(input_path, rep1_file))
    adata2 = sc.read_h5ad(os.path.join(input_path, rep2_file))
    
    # Save the original cell barcode for future reference
    adata1.obs["original_barcode"] = adata1.obs.index.astype(str).copy()
    adata2.obs["original_barcode"] = adata2.obs.index.astype(str).copy()
    adata1.obs["replicate-cellbarcode"] = "B4-" + adata1.obs["original_barcode"].astype(str)
    adata2.obs["replicate-cellbarcode"] = "C1-" + adata2.obs["original_barcode"].astype(str)

    # Add replicate identifiers and add a new column "replicate-cellbarcode" to remain_reads_info if available
    if 'remain_reads_info' in adata1.uns:
        adata1.uns['remain_reads_info']['replicate'] = 'B4'
        adata1.uns['remain_reads_info']['replicate-cellbarcode'] = "B4-" + adata1.uns['remain_reads_info']['cell_barcode'].astype(str)
    if 'remain_reads_info' in adata2.uns:
        adata2.uns['remain_reads_info']['replicate'] = 'C1'
        adata2.uns['remain_reads_info']['replicate-cellbarcode'] = "C1-" + adata2.uns['remain_reads_info']['cell_barcode'].astype(str)
    return adata1, adata2


def perform_qc(adata, qc_n):
    """
    Compute quality control metrics and determine lower and upper bounds for filtering 
    based on the median absolute deviation (MAD). 
    """
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    dense_matrix = adata.X.toarray()
    # Store the maximum counts for each gene in adata.var
    adata.var['max_counts_sample'] = dense_matrix.max(axis=0)
    mad = stats.median_abs_deviation(adata.obs['log1p_total_counts'], scale=1)
    median_val = adata.obs['log1p_total_counts'].median()
    lower_bd = median_val - qc_n * mad
    upper_bd = median_val + qc_n * mad
    return lower_bd, upper_bd


def filter_cells(adata, lower_bd, upper_bd, min_genes):
    """
    Filter cells based on:
      - A minimum number of genes expressed per cell.
      - Minimum and maximum total counts.
    The function prints the number of cells and genes after each filtering step.
    """
    print("Filter Cells: {} cells, {} genes initially".format(adata.n_obs, adata.n_vars))
    
    # Filter cells with fewer than the specified minimum number of genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print("After min_genes>={}: {} cells, {} genes".format(min_genes, adata.n_obs, adata.n_vars))
    
    # Filter cells with total counts below the lower bound
    sc.pp.filter_cells(adata, min_counts=np.expm1(lower_bd))
    print("After min_counts>= {}: {} cells, {} genes".format(np.expm1(lower_bd), adata.n_obs, adata.n_vars))
    
    # Filter cells with total counts above the upper bound
    sc.pp.filter_cells(adata, max_counts=np.expm1(upper_bd))
    print("After max_counts<= {}: {} cells, {} genes".format(np.expm1(upper_bd), adata.n_obs, adata.n_vars))
    
    # Save a copy of the raw count matrix in a new layer
    adata.layers['totalRNA_raw'] = adata.X.copy()


def filter_genes(adata, min_cells):
    """
    Filter genes based on being detected in a minimum number of cells.
    Prints the number of genes remaining after filtering and displays the first 50 passing genes.
    """
    print("Filter Genes: {} genes initially".format(adata.n_vars))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print("After min_cells>={}: {} genes remain".format(min_cells, adata.n_vars))
    passed_genes = adata.var.index.tolist()
    print("First 50 genes passing filter: {}".format(", ".join(passed_genes[:50])))


def sync_remain_reads_info(adata):
    """
    Sync the 'remain_reads_info' DataFrame with the filtered cells using the original barcodes.
    """
    if 'remain_reads_info' not in adata.uns:
        print("No remain_reads_info; skipping sync.")
        return
    df = adata.uns['remain_reads_info'].copy()
    before_count = df.shape[0]
    #df['cell_barcode'] = df['cell_barcode'].astype(str)
    #valid_barcodes = adata.obs["original_barcode"].astype(str).tolist()
    #df = df[(df['cell_barcode'].isin(valid_barcodes)) | (df['cell_barcode'] == '-1')]
    #print("remain_reads_info: {} -> {} after syncing".format(before_count, df.shape[0]))
    #adata.uns['remain_reads_info'] = df
    df['replicate-cellbarcode'] = df['replicate-cellbarcode'].astype(str)
    valid_barcodes = adata.obs["replicate-cellbarcode"].astype(str).tolist()
    df = df[(df['replicate-cellbarcode'].isin(valid_barcodes)) | (df['replicate-cellbarcode'].str.endswith('-1'))]
    print("remain_reads_info: {} -> {} records after syncing".format(before_count, df.shape[0]))
    adata.uns['remain_reads_info'] = df


def plot_correlation(adata1, adata2, output_path):
    """
    Plot a scatter plot of log2(total counts) between two replicates and calculate Pearson's correlation.
    Save the plot as 'correlation_plot.png'.
    """
    vec1 = np.log2(np.array(adata1.X.sum(axis=0)).flatten() + 1)
    vec2 = np.log2(np.array(adata2.X.sum(axis=0)).flatten() + 1)
    p_corr = stats.pearsonr(vec1, vec2)
    df = pd.DataFrame({'B4_WT': vec1, 'C1_WT': vec2})
    g = sns.lmplot(x='B4_WT', y='C1_WT', data=df, scatter_kws={'s': 1}, line_kws={'color': 'r'})
    g.set_axis_labels('B4_WT - log2(total counts)', 'C1_WT - log2(total counts)')
    plt.title(f"Pearson's correlation coefficient: {round(p_corr[0], 3)}")
    plt.savefig(os.path.join(output_path, 'correlation_plot.png'))
    plt.close()


def combine_replicates(adata1, adata2):
    """
    Combine two replicates into one AnnData object.
    Update metadata and add a combined 'protocol-replicate' field.
    """
    adata1.obs['protocol'] = 'STATES'
    adata1.obs['replicate'] = 'B4'
    adata2.obs['protocol'] = 'STATES'
    adata2.obs['replicate'] = 'C1'
    adata = ad.concat([adata1, adata2])
    print("After combining replicates: {} cells, {} genes".format(adata.n_obs, adata.n_vars))
    
    # Ensure the original barcode remains as string
    adata.obs["original_barcode"] = adata.obs["original_barcode"].astype(str)
    # Store max counts from each replicate in the variable annotations
    adata.var['max_counts_B4'] = adata1.var.get('max_counts_sample', np.nan).values
    adata.var['max_counts_C1'] = adata2.var.get('max_counts_sample', np.nan).values
    # Create a new column for combined protocol and replicate info
    adata.obs['protocol-replicate'] = adata.obs['protocol'].astype(str) + '-' + adata.obs['replicate'].astype(str)
    adata.obs['protocol-replicate'] = adata.obs['protocol-replicate'].astype('category')
    return adata


def merge_remain_reads_info(adata1, adata2, adata):
    """
    Merge the 'remain_reads_info' from both replicates into the combined AnnData object.
    """
    df1 = adata1.uns['remain_reads_info'] if 'remain_reads_info' in adata1.uns else None
    df2 = adata2.uns['remain_reads_info'] if 'remain_reads_info' in adata2.uns else None
    if df1 is not None and df2 is not None:
        merged = pd.concat([df1, df2], axis=0)
        adata.uns['remain_reads_info'] = merged
        print("Merged remain_reads_info shape:", merged.shape)
    elif df1 is not None:
        adata.uns['remain_reads_info'] = df1
    elif df2 is not None:
        adata.uns['remain_reads_info'] = df2
    else:
        print("No remain_reads_info found.")


def filter_low_expressed_genes(adata, gene_threshold):
    """
    Mark genes as detected if their max counts in both replicates exceed the specified threshold.
    """
    passed = ((adata.var['max_counts_B4'] > gene_threshold) &
              (adata.var['max_counts_C1'] > gene_threshold))
    print("Low expressed gene filter (threshold {}): {} genes pass".format(gene_threshold, passed.sum()))
    adata.var['detected'] = passed
    adata.var['highly_variable'] = passed


def plot_max_counts_histogram(adata, output_path, gene_threshold):
    """
    Plot histograms of maximum counts per gene for each replicate.
    Save the plots as PNG files.
    """
    for label in ['B4', 'C1']:
        counts = adata.var[f'max_counts_{label}']
        counts = counts[~np.isnan(counts)]
        bins = np.arange(0, 61, 1)
        plt.figure(figsize=(8, 5))
        sns.histplot(counts, bins=bins, kde=False)
        plt.axvline(gene_threshold, color='red', linestyle='--', label=f'Threshold = {gene_threshold}')
        plt.xlabel(f'Max Counts per Gene ({label})')
        plt.ylabel('Frequency')
        plt.title(f'Distribution of Max Counts per Gene ({label})')
        plt.xlim(0, 100)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f'max_counts_per_gene_histogram_{label}.png'))
        plt.close()


def normalize_and_scale(adata):
    """
    Normalize total counts, log-transform the data, scale the data, and regress out total counts.
    Save the normalized and scaled matrix in new layers.
    """
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    #adata.raw = adata 
    sc.pp.scale(adata)
    adata.layers['totalRNA_scaled'] = adata.X.copy()
    sc.pp.regress_out(adata, ['total_counts'])
    adata.layers['totalRNA_corrected'] = adata.X.copy()
    print("After normalization: {} cells, {} genes".format(adata.n_obs, adata.n_vars))


def map_nuclei(adata):
    """
    Map nuclei information from 'remain_reads_info' to the cell metadata based on the original barcode.
    """
    if 'remain_reads_info' not in adata.uns:
        print("No remain_reads_info; skipping nuclei mapping.")
        return
    df = adata.uns['remain_reads_info'].drop_duplicates('cell_barcode').set_index('cell_barcode')
    if 'nuclei' not in df.columns:
        print("No 'nuclei' column; skipping nuclei mapping.")
        return
    adata.obs['nuclei'] = adata.obs["original_barcode"].map(df['nuclei'])
    print("Nuclei mapping complete:\n", adata.obs['nuclei'].value_counts(dropna=False))


def create_cyto_layers(adata):
    """
    Create new layers for cytoplasmic raw counts and ribosomal counts.
    Gene names are cleaned by removing specific suffixes.
    """
    if 'remain_reads_info' not in adata.uns:
        print("No remain_reads_info; cannot create cyto layers.")
        return
    df = adata.uns['remain_reads_info'].copy()
    df['cell_barcode'] = df['cell_barcode'].astype(str)
    # Clean gene names by removing suffixes _R, _S, and _TE
    df['gene_name_clean'] = df['gene_name'].str.replace(r'(_R|_S|_TE)$', '', regex=True)
    df_cyto = df[df['nuclei'] == 0].copy()
    totalRNA_cyto = df_cyto.groupby(['cell_barcode', 'gene_name_clean']).size().unstack(fill_value=0)
    
    # For ribosomal reads, only include those whose original gene name ends with '_S'
    df_ntRNA_cyto = df_cyto[df_cyto['gene_name'].str.endswith('_S')]
    df_ntRNA_cyto['gene_name_clean'] = df_ntRNA_cyto['gene_name'].str.replace(r'(_R|_S|_TE)$', '', regex=True)
    ntRNA_cyto = df_ntRNA_cyto.groupby(['cell_barcode', 'gene_name_clean']).size().unstack(fill_value=0)
    
    # Align matrices to the same shape as the main data
    cell_list = adata.obs["original_barcode"].astype(str)
    gene_list = adata.var.index
    totalRNA_cyto = totalRNA_cyto.reindex(index=cell_list, columns=gene_list, fill_value=0)
    ntRNA_cyto = ntRNA_cyto.reindex(index=cell_list, columns=gene_list, fill_value=0)
    adata.layers['totalRNA_cyto'] = totalRNA_cyto.to_numpy(dtype=np.float32)
    adata.layers['ntRNA_cyto'] = ntRNA_cyto.to_numpy(dtype=np.float32)
    print("Cyto layers created.")


def perform_pca(adata, output_path):
    """
    Perform PCA on the processed data and generate several visualizations.
    Save the plots to the specified output directory.
    """
    sc.tl.pca(adata, svd_solver='full', use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, log=False)
    plt.savefig(os.path.join(output_path, 'pca_variance_ratio.png'))
    plt.close()
    sc.pl.pca(adata, color='total_counts')
    plt.savefig(os.path.join(output_path, 'pca_total_counts.png'))
    plt.close()
    sc.pl.pca(adata, color='n_genes')
    plt.savefig(os.path.join(output_path, 'pca_n_genes.png'))
    plt.close()
    g = sns.jointplot(x=adata.obsm['X_pca'][:, 0],
                      y=adata.obsm['X_pca'][:, 1],
                      hue=adata.obs['protocol-replicate'],
                      s=1)
    g.set_axis_labels('PC1', 'PC2')
    g.savefig(os.path.join(output_path, 'pca_jointplot.png'))
    plt.close()


def harmony_integration(adata, output_path):
    """
    Perform batch correction using Harmony and plot the corrected PCA.
    """
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    sce.pp.harmony_integrate(adata, 'protocol-replicate')
    g = sns.jointplot(x=adata.obsm['X_pca_harmony'][:, 0],
                      y=adata.obsm['X_pca_harmony'][:, 1],
                      hue=adata.obs['protocol-replicate'],
                      s=1)
    g.set_axis_labels('PC1', 'PC2')
    g.savefig(os.path.join(output_path, 'harmony_integration_jointplot.png'))
    plt.close()


def save_data(adata, output_path, final_fn):
    """
    Save the final AnnData object to an h5ad file.
    """
    final_file = os.path.join(output_path, final_fn)
    adata.write_h5ad(final_file)
    print("Data saved to:", final_file)

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess h5ad data with QC, filtering, normalization, PCA, and Harmony integration."
    )
    parser.add_argument("--input_path", type=str, required=True, help="Directory of input h5ad files")
    parser.add_argument("--output_path", type=str, default=None,
                        help="Output directory (default: input_path/250228_output)")
    parser.add_argument("--rep1_file", type=str, default="B4_WT_mousebrain2_0328.h5ad", help="Filename for replicate 1")
    parser.add_argument("--rep2_file", type=str, default="C1_WT_mousebrain1_0328.h5ad", help="Filename for replicate 2")
    parser.add_argument("--qc_n", type=float, default=3, help="Multiplier for MAD in QC (default: 3)")
    parser.add_argument("--gene_threshold", type=float, default=3, help="Threshold for low expressed genes (default: 3)")
    parser.add_argument("--min_cells", type=int, default=10, help="Minimum cells for gene filtering (default: 10)")
    parser.add_argument("--min_genes", type=int, default=10, help="Minimum genes per cell (default: 10)")
    parser.add_argument("--final_fn", type=str, default="0326-tissue-harmony.h5ad", help="Final output filename")
    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path if args.output_path else os.path.join(input_path, '250228_output')
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    sc.settings.figdir = output_path

    # Load replicate data
    adata1, adata2 = load_h5ad_files(input_path, args.rep1_file, args.rep2_file)
    
    # Perform QC on each replicate and get filtering bounds
    lower_bd1, upper_bd1 = perform_qc(adata1, args.qc_n)
    lower_bd2, upper_bd2 = perform_qc(adata2, args.qc_n)

    # Filter cells and genes for each replicate
    filter_cells(adata1, lower_bd1, upper_bd1, args.min_genes)
    filter_cells(adata2, lower_bd2, upper_bd2, args.min_genes)
    filter_genes(adata1, min_cells=args.min_cells)
    filter_genes(adata2, min_cells=args.min_cells)

    # Sync and merge additional read info
    sync_remain_reads_info(adata1)
    sync_remain_reads_info(adata2)
    plot_correlation(adata1, adata2, output_path)
    adata = combine_replicates(adata1, adata2)
    merge_remain_reads_info(adata1, adata2, adata)
    
    # Filter low expressed genes and plot histograms
    filter_low_expressed_genes(adata, args.gene_threshold)
    plot_max_counts_histogram(adata, output_path, args.gene_threshold)
    
    # Normalize, scale, and regress out total counts
    normalize_and_scale(adata)
    
    # Map nuclei info and create additional layers
    map_nuclei(adata)
    create_cyto_layers(adata)
    
    # Perform PCA and Harmony integration, then save the final object
    perform_pca(adata, output_path)
    harmony_integration(adata, output_path)
    save_data(adata, output_path, args.final_fn)


if __name__ == '__main__':
    main()