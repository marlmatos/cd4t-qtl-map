#!/bin/bash/python

#-------------------------------------------------------
#### Prepare pseudobulk metacells for MOFA modeling
#Author: Marliette Matos
#Date: 07/22/2024
#-------------------------------------------------------


#------------Import Modules --------------------------
import os
import scanpy as sc
import anndata
from mofapy2.run.entry_point import mofa
import pandas as pd
import numpy as np

wd="~/cd4_CellRegMap/001_preprocessing/results_01312025"

#set output directories
# Set output directories
os.makedirs(os.path.join(wd, "mofa_factors/mofa_trained/5factors"), exist_ok=True)
os.makedirs(os.path.join(wd, "mofa_factors/inputs/5factors"), exist_ok=True)
    
# Loop through the different resolution values
resolutions = [34, 4, 45, 5]

for res in resolutions: 
  result=f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_res{res}.hdf5"
  #-----------Read Data-------------------
  metacells_file=f'~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res{res}.h5'
  sce = sc.read_h5ad(metacells_file) #after pseudobulk variable genes and pca need to be recalculated 
  
  sce.var["MITO"] = sce.var_names.str.startswith("MT-")
  
  sce.var["RIBO"] = sce.var_names.str.startswith(("RPS", "RPL"))
  
  sc.pp.calculate_qc_metrics(sce, qc_vars=["MITO", "RIBO"], percent_top=None, log1p=True, inplace=True)
  
  #sc.pp.regress_out(sce, keys=['n_genes_by_counts'], n_jobs=None, copy=False)
  
  
  #-----------get top 500 variable features ------
  
  sc.pp.highly_variable_genes(sce, n_top_genes=500)
  
  #-----------subset for just varfeatures--------
  sce = sce[:, sce.var.highly_variable]
  
  
  # Run PCA
  sc.pp.scale(sce, max_value=10)
  sc.tl.pca(sce, svd_solver='arpack')
  
  loadings = sce.varm['PCs']
  pc_scores = sce.obsm['X_pca']
  # Create a DataFrame with PC scores
  pc_scores_df = pd.DataFrame(pc_scores, index=sce.obs_names, columns=[f'PC{i+1}' for i in range(50)])
  
  # Display the first few rows of the DataFrame
  print(pc_scores_df.head())
  
  # Save PCs
  file_path = f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/inputs/PCA_all_cells_all_conditions_Leiden_res{res}.csv"
  pc_scores_df.to_csv(file_path)
  
  # Save the AnnData object to an .h5ad file
  sce_path = f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/inputs/topHVF_Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res{res}.h5ad"
  
  sce.write(sce_path)
  
  #train mofa model
  m = mofa(sce, 
           expectations=["W","Z","AlphaW","AlphaZ"],
           n_factors=5,
           convergence_mode='medium',
           #groups_label="Group",
           outfile=result,
           seed=42,
           quiet=False)
  
  print("End of Script")
  print("Good-bye!")


