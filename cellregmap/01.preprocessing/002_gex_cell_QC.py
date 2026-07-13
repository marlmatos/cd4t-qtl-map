#-------------------------------------------------------
# Create metacells from sce experimen for CD4 T+ Cells
# Adapted by: Marliette Matos from CellRegMap analysis
# Github: https://github.com/annacuomo/CellRegMap_analyses/blob/main/neuroseq/preprocessing/create_metacells.py
#-------------------------------------------------------


import sys
import os
import re
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.offsetbox import AnchoredText
import seaborn as sns
from scipy import sparse


def plot_QC(adata, title, fig_home_dir, filtered=False):
    sc.settings.figdir = os.path.join(fig_home_dir, "figures/QC_figures_filtered/" if filtered else "figures/QC_figures/")
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    adata.var["MITO"] = adata.var_names.str.startswith("MT-") 
    adata.var["RIBO"] = adata.var_names.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["MITO", "RIBO"], percent_top=None, log1p=True, inplace=True)
   
    sc.pl.violin(adata, keys=['n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_MITO','pct_counts_RIBO'], 
                 use_raw=False, scale="count", jitter=False, groupby="cd4_subtypes_l1", rotation=45, 
                 xlabel=title, save=f"_QC_plots_{title.replace(' ', '-')}.png")

    sc.pl.scatter(adata, x='total_counts', y='total_counts_MITO', color="cd4_subtypes_l1", 
                  size=6, title=title, save=f"_plot_MITO_{title.replace(' ', '-')}.png")
    sc.pl.scatter(adata, x='total_counts', y='total_counts_RIBO', color="cd4_subtypes_l1", 
                  size=6, title=title, save=f"_plot_RIBO_{title.replace(' ', '-')}.png")

    for metric, label in [('log1p_total_counts', 'log1p'), ('total_counts', 'allcells'), ('n_genes_by_counts', 'genes')]:
        tmp = sns.distplot(adata.obs[metric], kde=False).set_title(title)
        tmp.get_figure().savefig(os.path.join(sc.settings.figdir, f"QC_{label}_per_barcode_{title.replace(' ', '-')}.png"))
        tmp.get_figure().clf()

    if not filtered:
        tmp = sns.distplot(adata.obs['total_counts'][adata.obs['total_counts']<10000], kde=False).set_title(title)
        tmp.get_figure().savefig(os.path.join(sc.settings.figdir, f"QC_n_counts_per_barcode_lower_tail_{title.replace(' ', '-')}.png"))
        tmp.get_figure().clf()
        
        for threshold in [2000, 2000]:
            compare = 'lower' if threshold == 2000 else 'upper'
            tmp = sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"].apply(lambda x: x < threshold if compare == 'lower' else x > threshold)], kde=False).set_title(title)
            tmp.get_figure().savefig(os.path.join(sc.settings.figdir, f"QC_n_genes_per_barcode_{compare}_tail_{title.replace(' ', '-')}.png"))
            tmp.get_figure().clf()

    sc.pl.scatter(adata, x="log1p_total_counts", y="n_genes_by_counts", size=6, color="pct_counts_MITO", 
                  title=title, save=f"_QC_combined_{title.replace(' ', '-')}log1p.png")
    
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", size=6, color="pct_counts_MITO", 
                  title=title, save=f"_QC_combined_{title.replace(' ', '-')}.png")
    
    sc.pl.scatter(adata[adata.obs.loc[(adata.obs.total_counts<10000) & (adata.obs.n_genes_by_counts <3500)].index], 
                  x="total_counts", y="n_genes_by_counts", size=6, color="pct_counts_MITO", 
                  title=title, save=f"_QC_combined_zoom-in_{title.replace(' ', '-')}.png")
    
    
def filter_scanpy_object(adata, obs_key='cd4_subtypes_l1'):
    adata_filtered = adata.copy()
    cell_types = adata_filtered.obs[obs_key].unique()
    stats_list = []
    combined_mask = np.zeros(adata_filtered.n_obs, dtype=bool)
    
    for cell_type in cell_types:
        mask = adata_filtered.obs[obs_key] == cell_type
        original_count = np.sum(mask)
        
        q1 = np.percentile(adata_filtered[mask].obs['n_genes_by_counts'], 5)
        q10 = np.percentile(adata_filtered[mask].obs['n_genes_by_counts'], 95)
        q99 = np.percentile(adata_filtered[mask].obs['total_counts'], 99)
        q1_log = np.percentile(adata_filtered[mask].obs['log1p_total_counts'], 1)
        q5_log = np.percentile(adata_filtered[mask].obs['log1p_total_counts'], 5)
        q10_ribo = np.percentile(adata_filtered[mask].obs['pct_counts_RIBO'], 98)
        
        # Use q5_log if cell_type is "cytotoxic", otherwise use q1_log
        log_threshold = q5_log if cell_type == "Cytotoxic" else q1_log
        
        final_mask = (
            (adata_filtered.obs[obs_key] == cell_type) & 
            (adata_filtered.obs['n_genes_by_counts'] > q1) & 
            (adata_filtered.obs['n_genes_by_counts'] <= q10) &
            (adata_filtered.obs['total_counts'] <= min(10000, q99)) &
            (adata_filtered.obs['log1p_total_counts'] > log_threshold) &  # Conditional threshold applied here 
            (adata_filtered.obs['pct_counts_RIBO'] <= q10_ribo)
        )
        
        filtered_count = np.sum(final_mask)
        percentage_remaining = (filtered_count / original_count) * 100
        
        stats_list.append({
            'CellType': cell_type,
            'OriginalCount': original_count,
            'FilteredCount': filtered_count,
            'PercentageRemaining': percentage_remaining
        })
        
        # Update the combined mask for this cell type
        combined_mask |= final_mask
    
    # Apply the combined mask to filter the adata object
    adata_filtered = adata_filtered[combined_mask]
    
    # Create a DataFrame with the statistics
    stats = pd.DataFrame(stats_list)
    return adata_filtered, stats



def main():
    # Set up environment variables
    sc.settings.verbosity = 3 
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    # Define working directory and read data
    wd = "/gpfs/commons/home/mmatos/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector"
    ct = 'allcells'
    filename = 'allcells_cd4_aging_filt_sce.h5ad'
    full_path = os.path.join(wd, filename)
    
    ct_all = sc.read_h5ad(full_path)
    print(ct_all)

    # Data preprocessing
    ct_all.obs['cd4_subtypes_l1'] = ct_all.obs['cd4_subtypes_l1'].astype(str).str.replace(' ', '_')
    print(f"The gene expression data has the following structure: {ct_all.shape}")
    print(f"The available reductions are the following: {ct_all.obsm}")

    # Perform initial QC
    plot_QC(ct_all, str(ct), wd)

    # Filter data
    ct_all_filtered, stats = filter_scanpy_object(ct_all)

    # Save filtering stats
    stats.to_csv(os.path.join(wd, 'preprocessing_filtering_stats.csv'), index=False)

    # Perform post-filtering QC
    plot_QC(ct_all_filtered, f"{ct}_filtered", wd, filtered=True)
    

    print("Script finished!")
    print("Goodbye!")

if __name__ == "__main__":
    main()

    
print("Script finished!")
