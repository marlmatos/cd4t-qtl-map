#################################################################
#### Script to extract MOFA modeling results for CD4T Cell RNAseq 
#### Author: Marliette Matos                                   
#### Project: CD4 Aging: CellRegMap of CD4 T cells             
#### Date: 0722/2024
#################################################################

#--------------Import Modules---------
import mofax as mfx
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad


#---------------Interate through all resolutions
res = "34"
#for res in resolutions:
print(f"Processing resolution: {res}")
m_results=f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_res{res}.hdf5"
adata=f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/inputs/topHVF_Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res{res}.h5ad"
adata = sc.read_h5ad(adata)

m = mfx.mofa_model(m_results)
m
df = pd.DataFrame(m.get_factors())
# cells x MOFA factors
df.shape
df.head()

# Name columns as MOFA1-MOFA25
df.columns = [f"MOFA{i}" for i in range(1, 6)]
print("Head of df with renamed columns:")
print(df.head())

# Name rows as cell names of the pseudocells
df.index = adata.obs_names
print("Head of df with renamed rows:")
print(df.head())
df.to_csv(f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_model_factors_res{res}.csv")


df_weights = pd.DataFrame(data = m.get_weights(), index=m.features['rna'])

# highly variable genes x MOFA factors
df_weights.shape

# Name columns as MOFA1-MOFA25
df_weights.columns = [f"MOFA{i}" for i in range(1, 6)]
print("Head of df with renamed columns:")
print(df_weights.head())

df_weights.head()
df_weights.to_csv(f"~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_model_weights_res{res}.csv")

