#-------------------------------------------------------
#### Prepare eQTL variant-gene pair file 
#Author: Marliette Matos
#Date: 07/24/2024


#------------Import Modules --------------------------
import os
import pandas as pd
import anndata
import scanpy as sc
from sklearn.preprocessing import StandardScaler, MinMaxScaler


missingness_threshold = 0.90  # Adjust threshold as needed

wd="~/cd4_CellRegMap/001_preprocessing/results_01312025"
os.makedirs(os.path.join(wd, "genotypes"), exist_ok=True)

resolutions=["34", "4"]
results="~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/"

for res in resolutions:
    #------------ SET UP EGENE LIST TO TEST --------------------------
    # Load and preprocess the cis-eQTL file
    cis_eqtl_file = "~/cd4_QTL_analysis/03_Run_cisQTL_perchr/analysis/cis_eQTLs/summary_independent_allcells_independent_cis_qtl_pairs.csv"
    cd4_eqtl = pd.read_csv(cis_eqtl_file)

    # Extract chromosome information and select relevant columns
    cd4_eqtl["chrom"] = cd4_eqtl["variant_id"].str.split(":").str[0].astype(int)
    cd4_eqtl = cd4_eqtl[["variant_id", "phenotype_id", "Cell.type", "chrom", "pval_nominal"]]  # Include pval_nominal

    # Sort by phenotype_id and p-value (ascending order)
    cd4_eqtl = cd4_eqtl.sort_values(by=['phenotype_id', 'pval_nominal'])

    # Drop duplicates based on phenotype_id, keeping the first occurrence
    # (keeping the one with the smaller p-value due to sorting)
    cd4_eqtl = cd4_eqtl.drop_duplicates(subset='phenotype_id')
    gene_list = cd4_eqtl[["phenotype_id"]]
    gene_list.shape

    # Output the number of filtered eGenes
    print(f"The number of eGenes to test in CellRegMap after filtering for sparsity at pseudobulk res {res} is:", len(gene_list))

    # Save the filtered gene list
    #gene_list_filtered = gene_list[gene_list['phenotype_id'].isin(egenes_filtered)]

    list_path = f"~/cd4_CellRegMap/001_preprocessing/results_01312025/genotypes/egene_list_res{res}.txt"
    gene_list.to_csv(list_path, index=False, header=False, quoting=2)


    # Filter the original eGene DataFrame to contain only the eGenes that pass the threshold
    #cd4_eqtl_filtered = cd4_eqtl[cd4_eqtl['phenotype_id'].isin(egenes_filtered)]

    # Save the filtered eGene file
    fvf_filename = f"~/cd4_CellRegMap/001_preprocessing/results_01312025/genotypes/cis_eQTLS_allcells_res{res}.csv"
    cd4_eqtl.to_csv(fvf_filename, index=False, quoting=2)

    
    #-------------SET UP MAPPING FILE ------------------

    #------------- Read Metacells Pseudobulk for each resolution ------------------
    # Load AnnData object from HDF5 file

    # Load gene expression matrix
    sample_mapping_file = results +'Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res'+res+'.h5'
    adata = anndata.read_h5ad(sample_mapping_file)
    print(adata)


    ## this file will map cells to donors 
    # Create a table of cell IDs and donor IDs
    ## "WGS_sampleID" are donor IDs, as found in the genotype matrix (G) and GRM covariance (K)
    ## "scRNA_Sample_Name" are cell IDs, as found in the scRNA-seq phenotype vector (y) and cell context covariance (C)
    ##  scRNA-seq phenotype vector (y)
    ## pseudocell_ID == contains the cellular ids for all the cells the make up the metacell pseudobulk 

    sample_mapping = pd.DataFrame({
        'Cell_ID': adata.obs.index.tolist(),
        'pseudocell_ID': adata.obs['cell_ids'].tolist(),
        'scRNA_Sample_Name': adata.obs['scRNA_Sample_Name'].tolist(),  # IDs for single RNA-seq experiment
        'WGS_sampleID': adata.obs['WGS_sampleID'].tolist()  # IDs for genotype WGS experiment
    })

    sample_mapping['WGS_sampleID'] = sample_mapping['WGS_sampleID'].apply(lambda x: str(int(x)) if isinstance(x, float) else str(x))
    sample_mapping['WGS_sampleID'] = 'g' + sample_mapping['WGS_sampleID']
    sample_mapping.head()
    donors = sample_mapping['WGS_sampleID'].unique()
    donors.sort()
    print("Number of unique donors: {}".format(len(donors)))
    sample_mapping.to_csv("~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_sample_mapping_res" + str(res) + ".csv", index=False,quoting=2)


    #-------------SET UP SAMPLE COVARIATE FILE -------------------

    # Calculate the percentage of RIBO counts per cell
    adata.var["MITO"] = adata.var_names.str.startswith("MT-") 
    adata.var["RIBO"] = adata.var_names.str.startswith(("RPS","RPL")) 
    
    adata.obs['total_counts'] = adata.X.sum(axis=1)
    ribo_genes = adata.var_names[adata.var['RIBO']]
    adata.obs['ribo_counts'] = adata[:, ribo_genes].X.sum(axis=1)
    adata.obs['percent_ribo'] = (adata.obs['ribo_counts'] / adata.obs['total_counts']) * 100
    print(adata.obs[['total_counts', 'ribo_counts', 'percent_ribo']].head())

    # Calculate the percentage of MITO counts per cell
    adata.obs['total_counts'] = adata.X.sum(axis=1)
    mito_genes = adata.var_names[adata.var['MITO']]
    adata.obs['mito_counts'] = adata[:, mito_genes].X.sum(axis=1)
    adata.obs['percent_mito'] = (adata.obs['mito_counts'] / adata.obs['total_counts']) * 100
    print(adata.obs[['total_counts', 'mito_counts', 'percent_mito']].head())

    covariates = pd.DataFrame({
        'Cell_ID': adata.obs.index.tolist(),
        'scRNA_Sample_Name': adata.obs['scRNA_Sample_Name'].tolist(),
        'Age_Group': adata.obs['Group'].tolist(),
        'Sex': adata.obs['Sex'].tolist(),
        'Batch': adata.obs['scRNA_batch_Dave.x'].tolist(),  # IDs for single RNA-seq experiment
        'WGS_sampleID': adata.obs['WGS_sampleID'].tolist(),  # IDs for genotype WGS experiment
        'percent_ribo': adata.obs['percent_ribo'].tolist(),
        'percent_mito': adata.obs['percent_mito'].tolist(),
        'cell_subtype': adata.obs['cd4_subtypes_l1'].tolist(),
    })


    #binarise all covariates

    covariates['Sex'] = pd.factorize(covariates['Sex'])[0] + 1
    covariates['Age_Group'] = pd.factorize(covariates['Age_Group'])[0] + 1

    # Perform one-hot encoding on the 'cell_subtype' column
    one_hot_encoded = pd.get_dummies(covariates['cell_subtype'], prefix='cs')
    one_hot_encoded = one_hot_encoded.applymap(lambda x: 1 if x else 0)
    print(one_hot_encoded.head())
    covariates = pd.concat([covariates, one_hot_encoded], axis=1)
    covariates = covariates.drop(columns=['cell_subtype'])
    print(covariates.head())

    # Perform one-hot encoding on the 'Batch' column
    one_hot_encoded = pd.get_dummies(covariates['Batch'], prefix='Batch')
    one_hot_encoded = one_hot_encoded.applymap(lambda x: 1 if x else 0)
    print(one_hot_encoded.head())
    covariates = pd.concat([covariates, one_hot_encoded], axis=1)
    covariates = covariates.drop(columns=['Batch'])
    print(covariates.head())


    #read age of individuals to add to covariates

    age_meta="~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/aging_cohort_age.csv"
    age_meta = pd.read_csv(age_meta)
    age_meta[["scRNA_Sample_Name", "age"]]

    covariates = pd.merge(covariates, age_meta[["scRNA_Sample_Name", "age"]], on='scRNA_Sample_Name', how='left')
    print(covariates.head())
    print(covariates.shape)
    covariates.to_csv("~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_covariates_res" + str(res) + ".csv", index=False,quoting=2)
    
print("CellRegMap Inputs prepared sucessfully!")
print("GoodBye!")

