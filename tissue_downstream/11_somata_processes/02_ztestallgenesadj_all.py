# conda: bigfish_env
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import fdrcorrection
from joblib import Parallel, delayed
from tqdm import tqdm

def main():
    # Load data
    adata = sc.read_h5ad('/storage/lingyuan2/STATES_data/mousebrain_harmony_20250614.h5ad')
    adata = adata[:, adata.var['detected'] == True]

    # Read processes TE data and create a dictionary
    print("Loading processes TE data...")
    df = pd.read_csv('/storage/lingyuan2/STATES_data/soma_process_eachgene_all_0503.csv')
    processes_te_dict = df.set_index('gene')['processes_TE'].to_dict()

    # Prepare for parallel processing
    print("Preparing data for parallel processing...")
    raw_values = adata.layers['totalRNA_raw'].toarray() if not isinstance(adata.layers['totalRNA_raw'], np.ndarray) else adata.layers['totalRNA_raw']
    te_values = adata.layers['TE'].toarray() if not isinstance(adata.layers['TE'], np.ndarray) else adata.layers['TE']
    
    # Find common genes and build index mapping
    common_genes = list(set(adata.var_names) & set(processes_te_dict.keys()))
    gene_to_idx = {gene: idx for idx, gene in enumerate(adata.var_names) if gene in common_genes}

    # Define function to process a single gene
    def process_gene(gene):
        idx = gene_to_idx[gene]
        raw = raw_values[:, idx]
        te = te_values[:, idx]
        
        # Vectorized filtering
        mask = (raw > 1) & (~np.isnan(te))
        filtered_te = te[mask]
        
        if len(filtered_te) > 1:
            stat, pval = ttest_1samp(filtered_te, processes_te_dict[gene])
            return {
                'gene': gene,
                'somata_mean_TE': filtered_te.mean(),
                'processes_TE': processes_te_dict[gene],
                't_statistic': stat,
                'p_value': pval,
                'n_cells': len(filtered_te)
            }
        return None

    # Process all genes in parallel
    print(f"Processing {len(common_genes)} genes in parallel...")
    results = Parallel(n_jobs=-1, verbose=10)(
        delayed(process_gene)(gene) for gene in tqdm(common_genes)
    )

    # Convert to DataFrame and filter out None results
    results_df = pd.DataFrame([r for r in results if r is not None])
    
    # Multiple testing correction
    print("Performing FDR correction...")
    results_df['adj_pvalue'] = fdrcorrection(results_df['p_value'])[1]
    results_df['significant'] = results_df['adj_pvalue'] < 0.05
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('adj_pvalue')

    # Save results
    print("Saving results...")
    results_df.to_csv('/storage/lingyuan2/STATES_data/gene_TE_comparison_results_with_FDR_all0503.csv', index=False)
    
    # Print top significant results
    print("\nTop significant results:")
    print(results_df.head(10))
    
    # Results summary
    n_sig = sum(results_df['significant'])
    print(f"\nFound {n_sig} significant genes at FDR < 0.05 (out of {len(results_df)} tested)")

if __name__ == '__main__':
    main()