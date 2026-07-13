#-------------------------------------------------------
# Create metacells from sce experiment for CD4 T+ Cells
# Adapted from CellRegMap analysis with optimized plotting
# Original: https://github.com/annacuomo/CellRegMap_analyses/blob/main/neuroseq/preprocessing/create_metacells.py
#-------------------------------------------------------

#------------------libraries------------------------------
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
from datetime import datetime

def log_message(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def create_additional_plots(adata, con, wd):
    """Create additional analysis plots for a condition"""
    log_message(f"Creating additional plots for condition: {con}")
    
    # Calculate cells per donor
    cells_donors = adata.obs.groupby("scRNA_Sample_Name", observed=True).size()
    log_message(f"Calculated cells per donor. Found {len(cells_donors)} donors")
    
    # Calculate pseudocells for each resolution
    resolutions = {'34': '3.4', '4': '4.0', '45': '4.5', '5': '5.0'}
    pseudocells_data = {}
    avg_cells_data = {}
    
    for res_key in resolutions.keys():
        col_name = f"leiden_res_{res_key}_euclidean"
        
        # Calculate number of pseudocells
        pseudocells_data[res_key] = adata.obs.groupby(
            "scRNA_Sample_Name", observed=True
        ).apply(lambda x: x[col_name].nunique())
        
        # Calculate average cells per pseudocell
        avg_cells_data[res_key] = adata.obs.groupby(
            ["scRNA_Sample_Name", col_name], observed=True
        ).size().groupby("scRNA_Sample_Name").mean()
    
    log_message("Creating cells vs pseudocells plot")
    # Plot 1: Cells vs Pseudocells
    plt.figure(figsize=(10, 8))
    for res_key, res_val in resolutions.items():
        plt.scatter(cells_donors, pseudocells_data[res_key], label=res_val)
    
    plt.xlabel('Cells per donor', size=16)
    plt.ylabel('Pseudocells per donor', size=17)
    plt.title(f"Cells vs Pseudocells - {con}", size=20)
    plt.legend(title="Leiden res", loc='center left', 
              bbox_to_anchor=(1.03, 0.5), frameon=False, 
              fontsize=14, title_fontsize=16)
    plt.tight_layout()
    
    output_path = os.path.join(wd, "figures/figures_per_donor", str(con),
                              f"all_cells_vs_Pseudocells_per_donor_Leiden_{con}.png")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    log_message(f"Saved cells vs pseudocells plot to: {output_path}")
    
    log_message("Creating cells vs average cells plot")
    # Plot 2: Cells vs Average Cells per Pseudocell
    plt.figure(figsize=(10, 8))
    for res_key, res_val in resolutions.items():
        plt.scatter(cells_donors, avg_cells_data[res_key], label=res_val)
    
    plt.xlabel('Cells per donor', size=16)
    plt.ylabel('Average cells per pseudocell per donor', size=17)
    plt.title(f"Cells vs Average Cells per Pseudocell - {con}", size=20)
    plt.legend(title="Leiden res", loc='center left', 
              bbox_to_anchor=(1.03, 0.5), frameon=False, 
              fontsize=14, title_fontsize=16)
    plt.tight_layout()
    
    output_path = os.path.join(wd, "figures/figures_per_donor", str(con),
                              f"all_cells_vs_Avg-Cells-Pseudocell_per_donor_Leiden_{con}.png")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    log_message(f"Saved average cells plot to: {output_path}")

def filter_scanpy_object_by_celltypes(adata, obs_key='cd4_subtypes_l1'):
    """Filter Scanpy object based on cell types and quality metrics"""
    adata_filtered = adata.copy()
    filtered_data = []
    
    # Get unique cell types from the obs_key column
    cell_types = adata_filtered.obs[obs_key].unique()
    
    for cell_type in cell_types:
        # Filter based on the current cell type
        mask = adata_filtered.obs[obs_key] == cell_type
        original_count = np.sum(mask)
        
        if original_count == 0:
            print(f"No cells found for cell type: {cell_type}")
            continue
        
        # Calculate percentiles
        q1 = np.percentile(adata_filtered[mask].obs['n_genes_by_counts'], 5)
        q10 = np.percentile(adata_filtered[mask].obs['n_genes_by_counts'], 95)
        q99 = np.percentile(adata_filtered[mask].obs['total_counts'], 99)
        q1_log = np.percentile(adata_filtered[mask].obs['log1p_total_counts'], 1)
        q5_log = np.percentile(adata_filtered[mask].obs['log1p_total_counts'], 5)
        q10_ribo = np.percentile(adata_filtered[mask].obs['pct_counts_RIBO'], 98)
        
        # Use q5_log if the cell type is "Cytotoxic", otherwise use q1_log
        log_threshold = q5_log if cell_type.lower() == "Cytotoxic" else q1_log
        
        # Create the final mask for filtering
        final_mask = (
            (adata_filtered.obs[obs_key] == cell_type) &
            (adata_filtered.obs['n_genes_by_counts'] > q1) &
            (adata_filtered.obs['n_genes_by_counts'] <= q10) &
            (adata_filtered.obs['total_counts'] <= min(10000, q99)) &
            (adata_filtered.obs['log1p_total_counts'] > log_threshold) &
            (adata_filtered.obs['pct_counts_RIBO'] <= q10_ribo)
        )
        
        # Apply the mask to filter the AnnData object
        filtered_subset = adata_filtered[final_mask]
        filtered_data.append(filtered_subset)
        
        # Calculate the number of remaining cells after filtering
        filtered_count = np.sum(final_mask)
        percentage_remaining = (filtered_count / original_count) * 100
        
        # Print some statistics about the filtering process
        print(f"Cell type: {cell_type}")
        print(f"Original count: {original_count}")
        print(f"Filtered count: {filtered_count}")
        print(f"Percentage remaining: {percentage_remaining:.2f}%")
        print("---")
    
    # Concatenate all filtered subsets
    if filtered_data:
        final_adata = sc.concat(filtered_data)
        print(f"Total cells after filtering: {final_adata.n_obs}")
        return final_adata
    else:
        print("No cells remained after filtering.")
        return None

def main():
    log_message("Starting CD4 T-Cell analysis")
    
    inputArgs = sys.argv[1:]
    #ct = inputArgs[0] if inputArgs else "allcells"
    
    # Set working directory
    wd = "/gpfs/commons/home/mmatos/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector"
    log_message(f"Working directory: {wd}")
    
    #----------------- read single cell experiment ------------------------
    ct = 'allcells'
    filename = 'allcells_cd4_aging_filt_sce.h5ad'
    # Construct the full path
    full_path = os.path.join(wd, filename)
    
    # Read the .h5ad file
    log_message(f"Reading data file: {full_path}")
    batch_all = sc.read_h5ad(full_path)
    print(batch_all)
    
    batch_all.var["MITO"] = batch_all.var_names.str.startswith("MT-") 
    batch_all.var["RIBO"] = batch_all.var_names.str.startswith(("RPS","RPL"))
    
    sc.pp.calculate_qc_metrics(batch_all, qc_vars=["MITO", "RIBO"], percent_top=None, log1p=True, inplace=True)
    
    #--------------filter low quality cells for scQTL standards-------------------------
    log_message("Filtering low quality cells")
    filtered_batch_all = filter_scanpy_object_by_celltypes(batch_all)
    
    # Create an empty list to store processed AnnData objects
    processed_adatas = []
    
    # Process each batch separately
    for batch in filtered_batch_all.obs['scRNA_batch_Dave.x'].unique():
        log_message(f"Processing batch: {batch}")
        batch_i = filtered_batch_all[filtered_batch_all.obs['scRNA_batch_Dave.x'] == batch]  
    
        # Create figures directory if it doesn't exist
        if os.path.isdir(os.path.join(wd, "figures/preprocessing_per_batch")):
            sc.settings.figdir = wd+"/figures/preprocessing_per_batch"
        else:
            os.makedirs(os.path.join(wd, "figures/preprocessing_per_batch"))
            sc.settings.figdir = wd+"/figures/preprocessing_per_batch"
    
        # Normalize data
        log_message(f"Normalizing batch {batch}")
        sc.pp.normalize_total(batch_i, target_sum=1e6, exclude_highly_expressed=True, max_fraction=0.10, key_added="Norm_factors", inplace=True)
    
        # Log2 transform
        sc.pp.log1p(batch_i, base=2)
    
        # Regress out confounding factors
        sc.pp.regress_out(batch_i, keys=['nFeature_RNA', 'nCount_RNA'])
    
        # Perform PCA
        sc.tl.pca(batch_i, n_comps=50, zero_center=True, svd_solver='arpack', random_state=32, use_highly_variable=False)
    
        # Create PCA plots
        log_message(f"Creating PCA plots for batch {batch}")
        fig = sc.pl.pca(batch_i, 
            color="cd4_subtypes_l1",
            annotate_var_explained=True,
            ncols=3,
            components=['1,2', '2,3', '3,4'],
            return_fig=True)
    
        fig.tight_layout()
        plt.savefig(wd+"/figures/preprocessing_per_batch/"+str(batch) + "_PCA_projection_celltypes_l1.png", bbox_inches='tight')
    
        fig = sc.pl.pca(batch_i, 
            color="predicted.celltype.l3",
            annotate_var_explained=True,
            ncols=3,
            components=['1,2', '2,3', '3,4'],
            return_fig=True)
    
        fig.tight_layout()
        plt.savefig(wd+"/figures/preprocessing_per_batch/"+str(batch) + "_PCA_projection_celltypes_seuratl3.png", bbox_inches='tight')
    
        fig = sc.pl.pca(batch_i, color=["nCount_RNA", "nFeature_RNA", "Sex", "Group", "scRNA_Sample_Name"], 
                        annotate_var_explained=True, dimensions=[(0,1)])
        plt.savefig(wd+"/figures/preprocessing_per_batch/" +str(batch) + "_PCA_projection_other_covariates.png", bbox_inches='tight')
        
        # Store processed batch
        processed_adatas.append(batch_i)
    
    # Get the list of unique batches to use as batch categories
    log_message("Concatenating processed data")
    batches = list(filtered_batch_all.obs['scRNA_batch_Dave.x'].unique())
    
    # Concatenate all processed data
    adata = processed_adatas[0].concatenate(
        processed_adatas[1:],
        join='inner',
        batch_key='scRNA_batch_Dave.x',
        batch_categories=batches,
        index_unique="_"
    )
    
    # Save normalized data
    log_message("Saving normalized data")
    os.makedirs(os.path.join(wd, "Data"), exist_ok=True)
    adata.write(os.path.join(wd, "Data", "allcells_cd4_aging_filt_sce_normalized.h5ad"))
    
    del processed_adatas  # Delete to improve memory usage
    

    # Harmony integration
    log_message("Performing Harmony integration")
    sce.pp.harmony_integrate(adata, key="scRNAseq_batch_combo", basis='X_pca', adjusted_basis='X_pca_harmony')
    
    # Save the Harmony integrated data
    log_message("Saving Harmony integrated data")
    adata.write(os.path.join(wd, "Data", "allcells_cd4_aging_filt_sce_norm_harmony_integrated.h5ad"))
    
    log_message("Re-reading normalized data after code failure")
    adata = sc.read_h5ad(os.path.join(wd, "Data", "allcells_cd4_aging_filt_sce_norm_harmony_integrated.h5ad"))
    # Create PCA plot with Harmony correction
    log_message("Creating Harmony-corrected PCA projection plot")
    fig = sc.pl.embedding(adata, 
                    basis='X_pca_harmony',
                    color=["cd4_subtypes_l1", 'scRNA_batch_Dave.x', "Sex", "Group"],
                    components=['1,2'],
                    title='Harmony-corrected PCA')
    plt.savefig(wd+"/figures/PCA__harmony_projection_color_covariates.png", bbox_inches='tight')
    
    #----------------- Donor representation per condition ------------------------
    log_message("Creating donor representation plots") 
    nd = adata.obs.scRNA_Sample_Name.nunique()
    
    # Create figures directory for donor plots
    if os.path.isdir(os.path.join(wd, "figures/figures_per_donor")):
        sc.settings.figdir = wd+"/figures/figures_per_donor"
    else:
        os.makedirs(os.path.join(wd, "figures/figures_per_donor"))
        sc.settings.figdir = wd+"/figures/figures_per_donor"
    
    # Donor representation per condition
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(20,20), constrained_layout=True)
    tmp = pd.crosstab(adata.obs['Group'], adata.obs['scRNA_Sample_Name'], normalize='index')
    tmp.plot.bar(stacked=True, ax=axs)
    axs.set_title('Donor representation per condition', fontsize=24)
    axs.set_xlabel('Group', fontsize=24)
    axs.set_ylabel('scRNA_Sample_Name', fontsize=24)
    plt.savefig(wd + "/figures/figures_per_donor/Donor_representation_per_AgeGroup.png", bbox_inches='tight')
    
    # Donor count per batch pie chart
    lib_pools = pd.DataFrame.from_dict(
        adata.obs.groupby(by='scRNA_batch_Dave.x').apply(lambda x: x.scRNA_Sample_Name.nunique()).to_dict(), 
        orient="index"
    ).rename(columns={0:"counts"})
    
    fig, ax = plt.subplots(ncols=1, figsize=(10,10), constrained_layout=True, subplot_kw=dict(aspect="equal"))
    ax.pie(
        x=[lib_pools.loc[lib_pools.counts == c].shape[0] for c in lib_pools.counts.unique()], 
        labels=lib_pools.counts.unique(), 
        autopct='%1.1f%%', 
        colors=["mediumblue", "r", "y"], 
        wedgeprops={'linewidth': 5}, 
        labeldistance=0.35, 
        radius=1.7, 
        textprops=dict(color="w", weight="bold", va="top", ma="center")
    )
    fig.suptitle("Number of donor per batch")
    plt.savefig(wd + "/figures/figures_per_donor/Number_of_donors_per_batch.png", bbox_inches='tight')
    
    # Process each condition
    for con in adata.obs.Group.unique():
        log_message(f"Processing Group: {con}")
        ct_con = adata[adata.obs.Group == con]
        
        # Create directory for condition plots
        if os.path.isdir(os.path.join(wd, "figures/figures_per_donor", str(con))):
            sc.settings.figdir = os.path.join(wd, "figures/figures_per_donor", str(con))
        else:
            os.makedirs(os.path.join(wd, "figures/figures_per_donor", str(con)))
            sc.settings.figdir = os.path.join(wd, "figures/figures_per_donor", str(con))
    
        #----------------- Rank donors according to number of cells ------------------------ 
        donors_cells = {}
        for d in ct_con.obs.scRNA_Sample_Name.unique():
            donors_cells[d] = ct_con.obs.loc[ct_con.obs.scRNA_Sample_Name == d].shape[0]
    
        # Histogram of cells per donor
        fig, axs = plt.subplots(ncols=1, figsize=(15,10), constrained_layout=True)
        sns.distplot(list(donors_cells.values()), bins='auto', 
                    hist_kws={'alpha':0.8, 'edgecolor':"darkblue", "facecolor":"mediumblue"}, 
                    kde=False, ax=axs)
        plt.grid(axis='both', alpha=0.5)
        plt.xlabel('Cells')
        plt.ylabel('Counts')
        
        if bool(re.search("52", str(con))):
            plt.title("Number of cells per donor in "+str(con).replace("_", " ")+"eated")
        else:
            plt.title("Number of cells per donor in the "+str(con)+" Group")
            
        axs.add_artist(AnchoredText(
            f'$median={int(np.median(list(donors_cells.values())))}$\n$max={max(list(donors_cells.values()))}$\n$min={min(list(donors_cells.values()))}$', 
            loc="upper right", frameon=False, 
            prop={"color":"slategray", "fontfamily":"fantasy", "fontsize":"small", "fontstretch":"expanded", "ma":"right"}
        ))
        plt.savefig(wd+"/figures/figures_per_donor/"+str(con)+"/Hist_number_of_all-cells_per_donor_"+str(con)+".png")
    
        #----------------- Remove donors with less than 20 cells in the given condition ------------------------ 
        low_cell_count_donors = [(k, v) for k, v in donors_cells.items() if v < 20]
        log_message(f"Found {len(low_cell_count_donors)} donors with less than 20 cells in condition {con}")
        
        for k, v in low_cell_count_donors:
            ct_con = ct_con[ct_con.obs.loc[ct_con.obs.scRNA_Sample_Name != k].index]
    
        # Process each donor
        i = 1
        all_donors_ct = None
        for d in ct_con.obs.scRNA_Sample_Name.unique():
            log_message(f"Processing donor {i} out of {ct_con.obs.scRNA_Sample_Name.nunique()} for condition {con}")
            ct_donor = ct_con[ct_con.obs.scRNA_Sample_Name == d]
            assert ct_donor.obs.scRNA_Sample_Name.nunique() == 1
    
            # Create directory for donor plots
            donor_dir = os.path.join(wd, "figures/figures_per_donor", str(con), str(d)+"_figures")
            if not os.path.isdir(donor_dir):
                os.makedirs(donor_dir)
            sc.settings.figdir = donor_dir
            
            # Compute neighbors and UMAP
            sc.pp.neighbors(ct_donor, n_neighbors=10, use_rep="X_pca_harmony", n_pcs=50, 
                           knn=True, method='umap', metric='euclidean', key_added="neighbors_euclidean")
            sc.tl.umap(ct_donor, min_dist=0.1, spread=1.0, n_components=2, maxiter=None, 
                      alpha=1.0, gamma=1.0, negative_sample_rate=5, init_pos='spectral', 
                      random_state=0, a=None, b=None, copy=False, method='umap', 
                      neighbors_key="neighbors_euclidean")
            
            # Generate Leiden clusters at different resolutions
            res = [3.4, 4, 4.5, 5]
            for r in res:
                rr = str(r).replace(".", "")
                sc.tl.leiden(ct_donor, resolution=r, neighbors_key="neighbors_euclidean", 
                            directed=False, key_added="leiden_res_"+rr+"_euclidean")
                sc.pl.umap(ct_donor, color="leiden_res_"+rr+"_euclidean", size=14, 
                          save="_Harmony_all-cells_"+str(d)+"_Leiden-res"+rr+"-Euclidean.png") 
            
            # Create histogram of cells per pseudocell for different resolutions
            fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(40, 20), constrained_layout=True)
    
            for r_idx, r in enumerate(res[:3]):  # Plot only the first 3 resolutions to match subplot count
                rr = str(r).replace(".", "")
                clusters_cells = ct_donor.obs.groupby(by="leiden_res_" + rr + "_euclidean").size()
    
                # Plotting histogram
                sns.histplot(clusters_cells, bins='auto', kde=False, alpha=0.8, 
                            edgecolor="darkblue", facecolor="mediumblue", ax=axs[r_idx])
                axs[r_idx].grid(axis='both', alpha=0.9)
                axs[r_idx].set_xlabel('Number of cells per pseudocell')
                axs[r_idx].set_ylabel('Counts')
                axs[r_idx].set_title(f'Leiden res={r}')
    
                # Adding anchored text annotation
                median_value = int(np.median(clusters_cells))
                max_value = max(clusters_cells)
                min_value = min(clusters_cells)
                annotation_text = f"$median={median_value}$\n$max={max_value}$\n$min={min_value}$"
                anchored_text = AnchoredText(annotation_text, loc="upper right", frameon=False, 
                                            prop={"color": "slategray", "fontfamily": "fantasy", 
                                                 "fontsize": "small", "fontstretch": "expanded"})
                axs[r_idx].add_artist(anchored_text)
    
            # Adding 4th resolution plot
            r = res[3]
            rr = str(r).replace(".", "")
            clusters_cells = ct_donor.obs.groupby(by="leiden_res_" + rr + "_euclidean").size()
            sns.histplot(clusters_cells, bins='auto', kde=False, alpha=0.8, 
                        edgecolor="darkblue", facecolor="mediumblue", ax=axs[3])
            axs[3].grid(axis='both', alpha=0.9)
            axs[3].set_xlabel('Number of cells per pseudocell')
            axs[3].set_ylabel('Counts')
            axs[3].set_title(f'Leiden res={r}')
            
            # Add annotation for the 4th plot
            median_value = int(np.median(clusters_cells))
            max_value = max(clusters_cells)
            min_value = min(clusters_cells)
            annotation_text = f"$median={median_value}$\n$max={max_value}$\n$min={min_value}$"
            anchored_text = AnchoredText(annotation_text, loc="upper right", frameon=False, 
                                        prop={"color": "slategray", "fontfamily": "fantasy", 
                                             "fontsize": "small", "fontstretch": "expanded"})
            axs[3].add_artist(anchored_text)
    
            # Add figure title and save
            fig.suptitle(f'{d} - {con}')
            plt.savefig(os.path.join(donor_dir, f"Number_of_cells_per_pseudocell_{con}_{d}_Euclidean_Leiden.png"))
            plt.close(fig)
    
            # Concatenate donor data
            if i == 1:
                all_donors_ct = ct_donor
            else:
                all_donors_ct = all_donors_ct.concatenate(ct_donor, join="inner", batch_key="DONOR")
    
            log_message(f"Current concatenated shape: {all_donors_ct.shape}")
            i += 1
    
        #----------------- Save metacells per cells per condition ------------------------ 
        log_message(f"Saving metacells for condition: {con}")
        dir_path = os.path.join(wd, "Data/Metacells_per_donor", str(con))
    
        # Create directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
            log_message(f"Directory {dir_path} created.")
        else:
            log_message(f"Directory {dir_path} already exists.")
    
        # Rename the '_index' column if it exists in the var DataFrame
        if '_index' in all_donors_ct.var.columns:
            all_donors_ct.var.rename(columns={'_index': 'index'}, inplace=True)
    
        # Check and rename '_index' in raw.var if it exists
        if all_donors_ct.raw is not None and '_index' in all_donors_ct.raw.var.columns:
            all_donors_ct.raw.var.rename(columns={'_index': 'index'}, inplace=True)
    
        # Save condition data
        all_donors_ct.write(os.path.join(dir_path, f"allcells_{con}_donors.h5"))
        
        # Create additional plots using the optimized function
        create_additional_plots(all_donors_ct, con, wd)
    
    log_message("Script finished!")
    log_message("Goodbye!")

if __name__ == "__main__":
    main()
